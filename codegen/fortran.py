from __future__ import annotations

import dataclasses
import textwrap
import typing
from abc import ABC, abstractmethod
from itertools import combinations

from .context import ctx, get_params
from .structs import TypeInformation
from .util import wrap_line

if typing.TYPE_CHECKING:
    from .routines import FortranRoutine, RoutineArg


def change_attributes(ftype: str, add: list[str] | None = None, remove: list[str] | None = None) -> str:
    exclude = frozenset(remove or {})
    if "intent" in exclude:
        exclude = frozenset(exclude | {"intent(in)", "intent(out)", "intent(inout)"})

    existing = {part.lower().strip() for part in ftype.split(",")}
    parts = ftype.split(",") + [item for item in add or [] if item.lower() not in existing]
    return ", ".join(part.strip() for part in parts if part.strip().lower() not in exclude)


@dataclasses.dataclass
class FortranWrapperArgument(ABC):
    arg: RoutineArg
    routine: FortranRoutine
    lines: list[str]
    have_err_flag: bool = False
    custom_call_arg_name: str | None = None

    @staticmethod
    def from_arg(
        arg: RoutineArg,
        routine: FortranRoutine,
        lines: list[str],
        have_err_flag: bool = True,
    ) -> FortranWrapperArgument:
        """Factory function to create the appropriate argument handler."""
        assert arg.member is not None
        arg_type = arg.member.type.lower()

        if arg_type == "character":
            return FortranWrapperStringArgument(arg, routine, lines, have_err_flag)
        if arg.array:
            if arg_type == "type":
                if arg.is_dynamic_array:
                    return FortranWrapperTypeArrayAllocatableArgument(arg, routine, lines, have_err_flag)
                return FortranWrapperTypeArrayArgument(arg, routine, lines, have_err_flag)
            if arg.is_dynamic_array:
                return FortranWrapperGeneralArrayAllocatableArgument(arg, routine, lines, have_err_flag)
            return FortranWrapperGeneralArrayArgument(arg, routine, lines, have_err_flag)

        if arg_type == "type":
            return FortranWrapperTypeArgument(arg, routine, lines, have_err_flag)
        return FortranWrapperGeneralArgument(arg, routine, lines, have_err_flag)

    @property
    def intent(self) -> str:
        """Get the intent of this argument, accounting for function results."""
        return "out" if self.arg.c_name == self.routine.result_name else self.arg.intent

    @property
    def is_function_result(self) -> bool:
        return self.arg.c_name == self.routine.result_name

    @property
    def f_name(self) -> str:
        return self.arg.f_name

    @property
    def c_name(self) -> str:
        return self.arg.c_name

    @staticmethod
    def if_block(cond: str, then: str) -> list[str]:
        """Add an if block to the lines."""
        if "\n" in then:
            then = textwrap.indent(then, " " * 4)
            return [
                f"  if ({cond}) then",
                then,
                "  endif",
            ]
        return [f"  if ({cond}) {then}"]

    @staticmethod
    def if_then_else_block(cond: str, then: str, else_: str) -> list[str]:
        """Add an if-then-else block to the lines."""
        then = textwrap.indent(then, " " * 4)
        else_ = textwrap.indent(else_, " " * 4)
        return [
            f"  if ({cond}) then",
            then,
            "  else",
            else_,
            "  endif",
        ]

    def error_if(self, cond: str, else_: str = "") -> list[str]:
        """Add error handling code to the lines."""
        else_lines = []
        if else_:
            else_lines = ["else", textwrap.indent(else_, "  ")]

        if "err_flag" in self.routine.declarations and self.have_err_flag:
            return [
                f"  if ({cond}) then",
                "    call c_f_pointer(err_flag, f_err_flag_ptr)",
                "    f_err_flag_ptr = .true.",
                "    return",
                *else_lines,
                "  endif",
            ]

        return [f"  if ({cond}) return"]

    def native_fortran_argument_type(self) -> str:
        """Extract Fortran type from definition and clean it up."""
        typ = TypeInformation.from_line(self.arg.member.definition)

        name = typ.type
        if typ.kind:
            kind = typ.kind.replace("fgsl_double", "c_double")
            return f"{name}({kind})"
        return name

    @abstractmethod
    def get_declarations(self) -> list[str]:
        """Get the declarations for this argument."""

    @abstractmethod
    def get_input_conversion(self) -> list[str]:
        """Get the code to convert this argument from C to Fortran."""

    @abstractmethod
    def get_output_conversion(self) -> list[str]:
        """Get the code to convert this argument from Fortran to C."""

    def get_call_arg_name(self) -> str:
        """Get the name of this argument for use in the function call."""
        if self.custom_call_arg_name:
            return self.custom_call_arg_name
        if self.arg.member.type_info.pointer and self.arg.full_type.type != "type":
            return f"{self.f_name}_ptr"
        return self.f_name


@dataclasses.dataclass
class FortranWrapperGeneralArrayAllocatableArgument(FortranWrapperArgument):
    """Handler for allocatable arrays of standard numerical/logical arguments."""

    @property
    def container_type(self) -> str:
        return f"{self.arg.type}_container_alloc"

    def get_declarations(self) -> list[str]:
        arg_ctype = "type(c_ptr), intent(in), value"
        arg_ftype = f"type({self.container_type}), pointer"

        return [
            f"  {arg_ctype} :: {self.c_name}",
            f"  {arg_ftype} :: {self.f_name}",
        ]

    def get_input_conversion(self) -> list[str]:
        result = [f"  !! container general array ({self.arg.full_type})"]
        code = f"  call c_f_pointer({self.c_name}, {self.f_name})"
        result.extend(self.if_block(f"assc({self.c_name})", code))
        self.custom_call_arg_name = f"{self.f_name}%data"
        return result

    def get_output_conversion(self) -> list[str]:
        return []


@dataclasses.dataclass
class FortranWrapperGeneralArrayArgument(FortranWrapperArgument):
    """Handler for arrays of standard numerical/logical arguments."""

    def get_declarations(self) -> list[str]:
        arg_ftype = self.native_fortran_argument_type()
        arg_ctype = "type(c_ptr), intent(in), value"
        tf_type = change_attributes(self.arg.transform.fortran_type, add=["pointer"], remove=["value"])
        arg_ftype = change_attributes(arg_ftype, remove=["intent"])
        dims = ",".join(self.arg.array)

        return [
            f"  {arg_ctype} :: {self.c_name}",
            f"  {arg_ftype} :: {self.f_name}({dims})",
            f"  {tf_type} :: {self.f_name}_ptr(:)",
        ]

    def get_input_conversion(self) -> list[str]:
        if self.intent == "out" or self.is_function_result:
            return []

        result = [f"  !! general array ({self.arg.full_type})"]
        arr = self.arg.array
        dimensions = "*".join(str(dim) for dim in self.arg.f_dims)

        if len(arr) == 1:
            code = textwrap.dedent(f"""\
                call c_f_pointer({self.c_name}, {self.f_name}_ptr, [{dimensions}])
                {self.f_name} = {self.f_name}_ptr(:)""")
        elif len(arr) == 2:
            code = textwrap.dedent(f"""\
                call c_f_pointer({self.c_name}, {self.f_name}_ptr, [{dimensions}])
                call vec2mat({self.f_name}_ptr, {self.f_name})""")
        elif len(arr) == 3:
            code = textwrap.dedent(f"""\
                call c_f_pointer({self.c_name}, {self.f_name}_ptr, [{dimensions}])
                call vec2tensor({self.f_name}_ptr, {self.f_name})""")
        else:
            raise NotImplementedError(len(arr))

        result.extend(self.if_block(f"assc({self.c_name})", code))
        return result

    def get_output_conversion(self) -> list[str]:
        if self.intent != "out" or not self.arg.array:
            return []

        result = [f"  ! {self.intent}: {self.f_name} {self.arg.full_type}"]
        arr = self.arg.array

        if len(arr) == 1:
            result.append(f"  if (assc({self.c_name})) {self.f_name}_ptr = {self.f_name}")
        elif len(arr) in {2, 3}:
            result.append(f"! TODO: output {len(arr)}D array: {self.c_name}")
        else:
            raise NotImplementedError(len(arr))
        return result


@dataclasses.dataclass
class FortranWrapperGeneralArgument(FortranWrapperArgument):
    """Handler for standard numerical/logical arguments."""

    native_type_conversion: bool = False

    def get_declarations(self) -> list[str]:
        arg_ftype = self.native_fortran_argument_type()
        arg_ctype = "type(c_ptr), intent(in), value"
        tf_type = self.arg.transform.fortran_type
        pointer_declarations = []

        assert self.arg.member is not None
        pointer_type = change_attributes(tf_type, add=["pointer"], remove=["value", "target"])

        if self.intent == "out":
            arg_ftype = change_attributes(arg_ftype, add=[], remove=["value", "target", "pointer"])
            pointer_declarations.append(f"  {pointer_type} :: {self.f_name}_ptr")
        elif self.intent == "inout" or self.arg.member.type_info.optional:
            if self.arg.type in {"logical", "real16"}:
                self.native_type_conversion = True
                native_type = self.arg.transform.fortran_native_type
                arg_ftype = f"{self.arg.transform.fortran_type}, pointer"
                pointer_declarations.append(f"  {native_type} :: {self.f_name}_native")
                self.custom_call_arg_name = f"{self.f_name}_native"
            else:
                arg_ftype = tf_type
            pointer_declarations.append(f"  {pointer_type} :: {self.f_name}_ptr")
        elif self.arg.member.type_info.pointer:
            self.custom_call_arg_name = self.f_name
            arg_ftype = change_attributes(arg_ftype, add=["pointer"], remove=["value", "target"])
        else:
            arg_ctype = self.arg.transform.fortran_type

        arg_ftype = change_attributes(arg_ftype, add=[], remove=["intent"])

        return [
            f"  {arg_ctype} :: {self.c_name}  ! {self.arg.full_type}",
            f"  {arg_ftype} :: {self.f_name}",
            *pointer_declarations,
        ]

    def get_input_conversion(self) -> list[str]:
        if self.is_function_result or self.intent == "out":
            return []

        result = [f"  ! {self.intent}: {self.f_name} {self.arg.full_type}"]

        if self.native_type_conversion:
            code = textwrap.dedent(f"""\
                call c_f_pointer({self.c_name}, {self.f_name}_ptr)
                {self.f_name}_native = {self.f_name}_ptr""")
            result.extend(self.if_then_else_block(f"assc({self.c_name})", code, f"! {self.c_name} unset"))

        elif self.arg.member.type_info.optional or self.intent == "inout":
            code = textwrap.dedent(f"""\
                call c_f_pointer({self.c_name}, {self.f_name}_ptr)
                {self.f_name} = {self.f_name}_ptr""")
            result.extend(self.if_then_else_block(f"assc({self.c_name})", code, f"! {self.c_name} unset"))

        elif self.arg.member.type_info.pointer:
            code = f"call c_f_pointer({self.c_name}, {self.f_name}) ! new case"
            result.extend(self.if_then_else_block(f"assc({self.c_name})", code, f"! {self.c_name} unset"))

        else:
            result.append(f"  {self.f_name} = {self.c_name}")

        return result

    def get_output_conversion(self) -> list[str]:
        if self.intent == "in":
            return []

        result = [f"  ! {self.intent}: {self.f_name} {self.arg.full_type}"]

        if self.native_type_conversion:
            code = textwrap.dedent(f"""\
                call c_f_pointer({self.c_name}, {self.f_name}_ptr)
                {self.f_name}_ptr = {self.f_name}_native""")
            result.extend(self.if_then_else_block(f"assc({self.c_name})", code, f"! {self.f_name} unset"))

        elif (
            self.arg.member.type_info.optional or self.intent == "inout" or self.arg.member.type_info.pointer
        ):
            code = textwrap.dedent(f"""\
                call c_f_pointer({self.c_name}, {self.f_name}_ptr)
                {self.f_name}_ptr = {self.f_name}""")
            result.extend(self.if_block(f"assc({self.c_name})", code))

        else:
            result.append(f"  call c_f_pointer({self.c_name}, {self.f_name}_ptr)")
            result.append(f"  {self.f_name}_ptr = {self.f_name}")

        return result


@dataclasses.dataclass
class FortranWrapperStringArgument(FortranWrapperArgument):
    """Handler for character/string arguments."""

    def get_declarations(self) -> list[str]:
        return [
            f"  type(c_ptr), intent(in), value :: {self.c_name}",
            f"  character(len=4096) :: {self.f_name}",
            f"  character(kind=c_char), pointer :: {self.f_name}_ptr(:)",
        ]

    def get_input_conversion(self) -> list[str]:
        if self.intent == "out" or self.is_function_result:
            return []

        result = [f"  ! {self.intent}: {self.f_name} {self.arg.full_type}"]

        code = textwrap.dedent(f"""\
            call c_f_pointer({self.c_name}, {self.f_name}_ptr, [huge(0)])
            call to_f_str({self.f_name}_ptr, {self.f_name})""")

        if self.arg.member.type_info.optional:
            result.extend(self.if_block(f"assc({self.c_name})", code))
        else:
            result.extend(self.error_if(f".not. assc({self.c_name})"))
            result.append(f"  call c_f_pointer({self.c_name}, {self.f_name}_ptr, [huge(0)])")
            result.append(f"  call to_f_str({self.f_name}_ptr, {self.f_name})")

        return result

    def get_output_conversion(self) -> list[str]:
        if self.intent == "in":
            return []

        result = [f"  ! {self.intent}: {self.f_name} {self.arg.full_type}"]
        if self.intent == "inout":
            result.append("  ! TODO i/o string (max length issue; buffer overflow...)")
        else:
            result.append(
                f"  call c_f_pointer({self.c_name}, {self.f_name}_ptr, [len_trim({self.f_name}) + 1]) ! output-only string"
            )
            result.append(f"  call to_c_str({self.f_name}, {self.f_name}_ptr)")

        return result


@dataclasses.dataclass
class FortranWrapperTypeArgument(FortranWrapperArgument):
    """Handler for TYPE arguments (derived types)."""

    @property
    def short_name(self) -> str:
        """Get the shortened type name without '_struct'."""
        assert self.arg.member.kind is not None
        return self.arg.member.kind.removesuffix("_struct")

    def get_declarations(self) -> list[str]:
        arg_ctype = self.arg.transform.fortran_type
        arg_ftype = change_attributes(self.native_fortran_argument_type(), add=["pointer"])

        return [
            f"  {arg_ctype} :: {self.c_name}  ! {self.arg.full_type}",
            f"  {arg_ftype} :: {self.f_name}",
        ]

    def get_input_conversion(self) -> list[str]:
        if self.is_function_result:
            return []

        result = [f"  ! {self.intent}: {self.f_name} {self.arg.full_type}"]
        code = f"  call c_f_pointer({self.c_name}, {self.f_name})"

        if self.arg.is_optional:
            result.extend(self.if_block(f"assc({self.c_name})", code))
        else:
            result.extend(self.error_if(f".not. assc({self.c_name})"))
            result.append(code)

        return result

    def get_output_conversion(self) -> list[str]:
        if self.intent != "out":
            return []
        return [
            f"  ! {self.intent}: {self.f_name} {self.arg.full_type}",
            f"  ! TODO may require output conversion? {self.arg.full_type}",
        ]


@dataclasses.dataclass
class FortranWrapperTypeArrayAllocatableArgument(FortranWrapperTypeArgument):
    @property
    def container_type(self) -> str:
        return f"{self.arg.member.type_info.kind}_container_alloc"

    def get_declarations(self) -> list[str]:
        arg_ctype = "type(c_ptr), intent(in), value"
        arg_ftype = f"type({self.container_type}), pointer"

        return [
            f"  {arg_ctype} :: {self.c_name}",
            f"  {arg_ftype} :: {self.f_name}",
        ]

    def get_input_conversion(self) -> list[str]:
        result = [f"  !! container type array ({self.arg.full_type})"]
        code = f"  call c_f_pointer({self.c_name}, {self.f_name})"
        result.extend(self.if_block(f"assc({self.c_name})", code))
        self.custom_call_arg_name = f"{self.f_name}%data"
        return result

    def get_output_conversion(self) -> list[str]:
        return []


@dataclasses.dataclass
class FortranWrapperTypeArrayArgument(FortranWrapperTypeArgument):
    def get_declarations(self) -> list[str]:
        arg_ctype = "type(c_ptr), intent(in), value"
        arg_ftype = f"type({self.arg.member.type_info.kind}), pointer"

        return [
            f"  {arg_ctype} :: {self.c_name}",
            f"  {arg_ftype} :: {self.f_name}(:)",
        ]

    def get_input_conversion(self) -> list[str]:
        result = [f"  !! type array ({self.arg.full_type})"]
        dimensions = "*".join(str(dim) for dim in self.arg.f_dims)
        result.append(f"  call c_f_pointer({self.c_name}, {self.f_name}, [{dimensions}])")
        return result

    def get_output_conversion(self) -> list[str]:
        if self.intent != "out":
            return []
        return [
            f"  ! {self.intent}: {self.f_name} {self.arg.full_type}",
            f"  ! TODO may require output conversion? {self.arg.full_type}",
        ]


def generate_fortran_routine_with_c_binding(routine: FortranRoutine) -> str:
    assert routine.docstring is not None
    lines = []
    args = routine.args
    structs = ctx().codegen_structs_by_name

    imports = {}
    for arg in args:
        if arg.type == "type":
            module = structs[arg.kind].module
            imports.setdefault(module, set())
            imports[module].add(arg.member.kind)

    have_err_flag = any(arg.c_name == "err_flag" and arg.intent == "out" for arg in args)
    arg_names = [arg.c_name for arg in routine.args]
    routine_and_args = f"fortran_{routine.name} ({', '.join(arg_names)})"

    lines.append(wrap_line(f"subroutine {routine_and_args} bind(c)"))

    for module, imps in imports.items():
        lines.append(f"  use {module}, only: " + ", ".join(sorted(imps)))

    handlers = [FortranWrapperArgument.from_arg(arg, routine, lines, have_err_flag) for arg in args]

    lines.append("  implicit none")

    by_intent = {"in": [], "out": [], "inout": []}
    for handler in handlers:
        by_intent[handler.intent].extend(handler.get_declarations())

    for intent, defns in by_intent.items():
        if defns:
            lines.extend([f"  ! ** {intent.capitalize()} parameters **", *defns])

    lines.append("  ! ** End of parameters **")

    for handler in handlers:
        lines.extend(handler.get_input_conversion())

    required_args = [handler for handler in handlers if not handler.arg.is_optional]
    optional_args = [handler for handler in handlers if handler.arg.is_optional]
    optional_combinations = []
    for i in range(len(optional_args), -1, -1):
        optional_combinations.extend([list(combo) for combo in combinations(optional_args, i)])

    handlers_by_name = {handler.arg.c_name: handler for handler in handlers}
    bmad_arg_names = {
        arg.c_name: decl for decl, arg in zip(routine.arg_names_with_result, routine.args, strict=True)
    }

    def get_f_args(arg_handlers: list[FortranWrapperArgument]):
        for handler in arg_handlers:
            name = handler.c_name
            if name == routine.result_name:
                continue
            decl = bmad_arg_names[name]
            yield f"{decl}=" + handler.get_call_arg_name()

    def add_call(f_args: str, indent: str = "  "):
        if routine.result_name:
            f_to_c_name = get_params().c_side_name_translation
            res_name = f_to_c_name.get(routine.result_name, routine.result_name)
            res = handlers_by_name[res_name].custom_call_arg_name or f"f_{routine.result_name}"
            lines.append(wrap_line(f"{res} = {routine.name}({f_args})", indent=indent))
        else:
            lines.append(wrap_line(f"call {routine.name}({f_args})", indent=indent))

    if len(optional_combinations) == 1 and not optional_combinations[0]:
        f_args = ", ".join(get_f_args(required_args))
        add_call(f_args, indent="  ")
    else:
        for idx, optional_combo in enumerate(optional_combinations):
            assoc = " .and. ".join(f"assc({arg.c_name})" for arg in optional_combo)
            call_args = get_f_args([*required_args, *optional_combo])

            if_ = "if" if idx == 0 else "elseif"
            if not assoc:
                lines.append("  else")
            else:
                lines.append(wrap_line(f"{if_} ({assoc}) then", indent="  "))

            f_args = ", ".join(call_args)
            add_call(f_args, indent="    ")

        lines.append("  endif")

    for handler in handlers:
        lines.extend(handler.get_output_conversion())

    lines.append("end subroutine")
    return "\n".join(lines)
