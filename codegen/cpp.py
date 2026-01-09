from __future__ import annotations

import dataclasses
import re
import string
import textwrap
import typing
from abc import ABC, abstractmethod

from .arg import Argument, CodegenStructure
from .context import get_params
from .proxy import templates as proxy_templates
from .types import STANDARD_TYPES, native_type_containers
from .util import snake_to_camel, sorted_routines

if typing.TYPE_CHECKING:
    from .routines import FortranRoutine, RoutineArg, RoutineSettings


@dataclasses.dataclass
class CppWrapperArgument(ABC):
    arg: RoutineArg
    decl_arg_name: str = dataclasses.field(init=False)
    fortran_call_arg_name: str = dataclasses.field(init=False)

    def __post_init__(self):
        self.fortran_call_arg_name = f"_{self.arg.c_name}"
        f_to_c_name = get_params().c_side_name_translation
        self.decl_arg_name = f_to_c_name.get(self.arg.c_name, self.arg.c_name)

    @classmethod
    def from_arg(cls, arg: RoutineArg):
        """Factory function to create the appropriate argument type."""
        assert arg.member is not None
        if len(arg.array):
            if arg.member.type.lower() == "type":
                if arg.member.type_info.allocatable or arg.is_dynamic_array:
                    return CppWrapperTypeArgumentAllocArray(arg)
                return CppWrapperTypeArgumentArray(arg)
            if arg.type == "character":
                return CppWrapperStringArgumentArray(arg)
            if arg.is_dynamic_array:
                return CppWrapperGeneralArgumentAllocArray(arg)
            return CppWrapperGeneralArgumentArray(arg)
        if arg.member.type.lower() == "type":
            return CppWrapperTypeArgument(arg)
        if arg.type == "character":
            return CppWrapperStringArgument(arg)
        return CppWrapperGeneralArgument(arg)

    def unwrap_optional(self, if_set: str, type_: str = "auto*", comment: str = "") -> str:
        """Helper function for optional arguments."""
        comment = f"// {comment}" if comment else ""
        return f"{type_} {self.fortran_call_arg_name} = {self.arg.c_name}.has_value() ? {if_set} : nullptr;{comment}"

    @property
    def cpp_decl(self) -> str:
        return f"{self.arg.transform.cpp_type} {self.arg.c_name}"

    @abstractmethod
    def pre_call_lines(self) -> list[str]:
        """Lines of code to execute before the Fortran call."""
        return []

    @abstractmethod
    def call_argument(self) -> str:
        """The argument string to pass to the Fortran call."""

    def post_call_lines(self) -> list[str]:
        """Lines of code to execute after the Fortran call."""
        return []

    def return_value(self) -> str | None:
        """The value to return from the wrapper function, if any."""
        if self.arg.intent == "out":
            return self.output_arg_name
        return None

    @property
    def output_arg_name(self) -> str | None:
        """Name of the variable holding the output value."""
        if self.arg.intent == "out":
            return self.fortran_call_arg_name
        return None

    def struct_decl(self, ignore_intent: bool = False) -> tuple[str, str] | None:  # noqa: ARG002
        """Return (type, name) for struct declaration if this is an output argument."""
        return None


@dataclasses.dataclass
class CppWrapperTypeArgumentAllocArray(CppWrapperArgument):
    def pre_call_lines(self) -> list[str]:
        lines = [f"// intent={self.arg.intent} allocatable type array"]
        argname = self.arg.c_name
        clsname = self.arg.kind_as_cpp_class
        container_cls = f"{clsname}Alloc1D"

        if self.arg.intent in {"in", "inout"}:
            if self.arg.is_optional:
                lines.append(
                    self.unwrap_optional(
                        f"{argname}->get().get_fortran_ptr()",
                        type_="auto*",
                        comment="input, optional",
                    )
                )
        elif self.arg.intent == "out":
            lines.append(f"auto {argname} {{ {container_cls}() }};")

        return lines

    def call_argument(self) -> str:
        argname = self.arg.c_name
        if self.arg.intent in {"in", "inout"}:
            if self.arg.is_optional:
                return self.fortran_call_arg_name
            return f"{argname}.get_fortran_ptr()"
        if self.arg.intent == "out":
            return f"{argname}.get_fortran_ptr()"
        return ""

    @property
    def output_arg_name(self) -> str | None:
        if self.arg.intent == "out":
            return f"std::move({self.arg.c_name})"
        return super().output_arg_name

    def struct_decl(self, ignore_intent: bool = False) -> tuple[str, str] | None:
        if self.arg.intent == "out" or ignore_intent:
            clsname = self.arg.kind_as_cpp_class
            container_cls = f"{clsname}Alloc1D"
            return container_cls, self.arg.c_name
        return None


@dataclasses.dataclass
class CppWrapperTypeArgumentArray(CppWrapperArgument):
    def pre_call_lines(self) -> list[str]:
        lines = []
        argname = self.arg.c_name
        clsname = self.arg.kind_as_cpp_class

        if self.arg.intent in {"in", "inout"}:
            if self.arg.is_optional:
                lines.append(
                    self.unwrap_optional(
                        f"{argname}->data()",
                        type_="auto*",
                        comment="input, optional",
                    )
                )
        elif self.arg.intent == "out":
            lines.append(f"""\
        // Output-only type array
        auto {argname} = {clsname}Array1D::allocate({self.arg.c_dim1}, 1);
        """)
        return lines

    def call_argument(self) -> str:
        argname = self.arg.c_name
        if self.arg.intent in {"in", "inout"}:
            if self.arg.is_optional:
                return self.fortran_call_arg_name
            return f"{argname}.data()"
        if self.arg.intent == "out":
            return f"{argname}.get_fortran_ptr()"
        return ""

    @property
    def output_arg_name(self) -> str | None:
        if self.arg.intent == "out":
            return f"std::move({self.arg.c_name})"
        return super().output_arg_name

    def struct_decl(self, ignore_intent: bool = False) -> tuple[str, str] | None:
        if self.arg.intent == "out" or ignore_intent:
            clsname = self.arg.kind_as_cpp_class
            return f"{clsname}Array1D", self.arg.c_name
        return None


@dataclasses.dataclass
class CppWrapperTypeArgument(CppWrapperArgument):
    def pre_call_lines(self) -> list[str]:
        lines = []
        arg = self.arg
        intent = arg.intent

        if intent in {"in", "inout"} and not arg.is_optional:
            if self.arg.member.type_info.pointer:  # pointer
                lines.append(
                    f"auto {self.fortran_call_arg_name} = &{self.arg.c_name}; // input, required, pointer"
                )
        elif intent in {"in", "inout"} and arg.is_optional:
            lines.append(
                self.unwrap_optional(
                    f"{self.arg.c_name}->get().get_fortran_ptr()",
                    type_="auto*",
                    comment="input, optional",
                )
            )
        elif intent == "out":
            if self.arg.member.type_info.pointer:
                lines.append(f"void *{self.fortran_call_arg_name};")
            else:
                lines.append(f"{self.arg.c_class} {self.fortran_call_arg_name};")

        return lines

    def call_argument(self) -> str:
        intent = self.arg.intent

        if intent in {"in", "inout"} and not self.arg.is_optional:
            if self.arg.member.type_info.pointer:
                return "&" + self.arg.c_name
            return f"{self.arg.c_name}.get_fortran_ptr()"
        if intent in {"in", "inout"} and self.arg.is_optional:
            return self.fortran_call_arg_name
        if intent == "out":
            if self.arg.member.type_info.pointer:
                return f"&{self.fortran_call_arg_name}"
            return f"{self.fortran_call_arg_name}.get_fortran_ptr()"
        return ""

    @property
    def output_arg_name(self) -> str | None:
        if self.arg.intent == "in":
            return self.arg.c_name
        if self.arg.intent == "out" and self.arg.member.type_info.pointer:
            return f"{self.arg.c_class}({self.fortran_call_arg_name})"
        return self.fortran_call_arg_name

    def struct_decl(self, ignore_intent: bool = False) -> tuple[str, str] | None:
        if self.arg.intent == "out" or ignore_intent:
            return self.arg.c_class, self.arg.c_name
        return None


@dataclasses.dataclass
class CppWrapperStringArgumentArray(CppWrapperArgument):
    def pre_call_lines(self) -> list[str]:
        lines = []
        argname = self.arg.c_name
        callname = self.fortran_call_arg_name
        clsname = "char"
        pbtype = self.arg.transform.cpp_declare_type

        if self.arg.intent in {"out", "inout"}:
            if self.arg.intent == "out":
                lines.append(
                    textwrap.dedent(f"""\
                    // Output-only string array
                    {pbtype} {argname};
                    std::vector<const {clsname}*> {callname}_ptr{{ {self.arg.c_dim1} }};
                    for (size_t i{{0}}; i < {argname}.size(); i++) {{
                        {callname}_ptr.push_back(&{argname}.data()[i]);
                    }}
                    """)
                )
        elif self.arg.is_optional:
            lines.append(
                textwrap.dedent(f"""\
                // Optional string array
                std::vector<const {clsname}*> {callname}{{ {argname}.has_value() ? {argname}->size() : 0 }};
                for (size_t i{{0}}; i < {callname}.size(); i++) {{
                    {callname}.push_back({argname}->data()[i].data());
                }}
                """)
            )
        else:
            lines.append(
                textwrap.dedent(f"""\
                // String array
                std::vector<const {clsname}*> {callname}{{ {argname}.size() }};
                for (size_t i{{0}}; i < {argname}.size(); i++) {{
                    {callname}.push_back({argname}[i].data());
                }}
                """)
            )
        return lines

    def call_argument(self) -> str:
        callname = self.fortran_call_arg_name
        if self.arg.intent == "out":
            return f"{callname}_ptr.data()"
        return f"{callname}.size() ? {callname}.data() : nullptr"

    @property
    def output_arg_name(self) -> str | None:
        if self.arg.intent == "out":
            return self.arg.c_name
        return super().output_arg_name

    def struct_decl(self, ignore_intent: bool = False) -> tuple[str, str] | None:
        if self.arg.intent == "out" or ignore_intent:
            return self.arg.transform.cpp_declare_type, self.arg.c_name
        return None


@dataclasses.dataclass
class CppWrapperStringArgument(CppWrapperArgument):
    def pre_call_lines(self) -> list[str]:
        lines = []
        as_pointer = self.arg.is_optional or self.arg.intent in {"inout", "out"}

        if not as_pointer:
            lines.append(f"auto {self.fortran_call_arg_name} = {self.arg.c_name}.c_str();")
        elif self.arg.intent == "in":
            if self.arg.is_optional:
                lines.append(self.unwrap_optional(f"{self.arg.c_name}->c_str()", type_="const char *"))
            else:
                lines.append(f"auto {self.fortran_call_arg_name} = {self.arg.c_name}.c_str();")
        elif self.arg.intent == "inout":
            if self.arg.is_optional:
                lines.append(self.unwrap_optional(f"{self.arg.c_name}->get().c_str()", type_="const char *"))
            else:
                lines.append(
                    f"auto {self.fortran_call_arg_name} = {self.arg.c_name}.c_str();// ptr, inout, required"
                )
        elif self.arg.intent == "out":
            lines.append(f"char {self.fortran_call_arg_name}[4096];")  # Fixed size buffer
        else:
            raise NotImplementedError("unknown string handling")
        return lines

    def call_argument(self) -> str:
        return self.fortran_call_arg_name

    def struct_decl(self, ignore_intent: bool = False) -> tuple[str, str] | None:
        if self.arg.intent == "out" or ignore_intent:
            return "std::string", self.arg.c_name
        return None


_alloc_array_type_map = {
    ctr.name: ctr.cpp_container_name
    for ctr in native_type_containers
    # "real": "RealAlloc1D",
    # "real8": "Real8Alloc1D",
    # "integer": "IntAlloc1D",
    # "integer8": "Int8Alloc1D",
    # "logical": "BoolAlloc1D",
    # "complex": "ComplexAlloc1D",
}


@dataclasses.dataclass
class CppWrapperGeneralArgumentAllocArray(CppWrapperArgument):
    def pre_call_lines(self) -> list[str]:
        lines = [f"// intent={self.arg.intent} allocatable general array"]
        argname = self.arg.c_name

        try:
            container_cls = _alloc_array_type_map[self.arg.type]
        except KeyError as ex:
            raise NotImplementedError(f"Allocatable array for type {self.arg.type} not implemented") from ex

        if self.arg.intent in {"in", "inout"}:
            if self.arg.is_optional:
                lines.append(
                    self.unwrap_optional(
                        f"{argname}->get().get_fortran_ptr()",
                        type_="auto*",
                        comment="input, optional",
                    )
                )
        elif self.arg.intent == "out":
            lines.append(f"auto {argname} {{ {container_cls}() }};")

        return lines

    def call_argument(self) -> str:
        argname = self.arg.c_name
        if self.arg.intent in {"in", "inout"}:
            if self.arg.is_optional:
                return self.fortran_call_arg_name
            return f"{argname}.get_fortran_ptr()"
        if self.arg.intent == "out":
            return f"{argname}.get_fortran_ptr()"
        return ""

    @property
    def output_arg_name(self) -> str | None:
        if self.arg.intent == "out":
            return f"std::move({self.arg.c_name})"
        return super().output_arg_name

    def struct_decl(self, ignore_intent: bool = False) -> tuple[str, str] | None:
        if self.arg.intent == "out" or ignore_intent:
            container_cls = _alloc_array_type_map[self.arg.type]
            return container_cls, self.arg.c_name
        return None


@dataclasses.dataclass
class CppWrapperGeneralArgumentArray(CppWrapperArgument):
    def pre_call_lines(self) -> list[str]:
        lines = []
        arr = self.arg.array
        opt_val = f"{self.arg.c_name}.value().data()"

        if len(arr) == 1:
            if self.arg.intent == "out":
                pbtype = self.arg.transform.cpp_declare_type
                lines.append(f"{pbtype} {self.fortran_call_arg_name};")
            elif self.arg.is_optional:
                lines.append(self.unwrap_optional(opt_val, type_=self.arg.transform.cpp_call_fortran_type))
            else:
                lines.append(
                    f"auto *{self.fortran_call_arg_name} = {self.arg.c_name}.data(); // CppWrapperGeneralArgument"
                )

        elif len(arr) in {2, 3}:
            ctype = STANDARD_TYPES[self.arg.type].c_type
            if self.arg.intent == "out":
                lines.append(f"{self.arg.transform.cpp_declare_type} {self.arg.c_name};")

            dims = "*".join(str(dim) for dim in self.arg.c_dims)
            vec_name = f"{self.fortran_call_arg_name}_vec"
            lines.append(f"{ctype} {vec_name}[{dims}];")

            prefix = "matrix" if len(arr) == 2 else "tensor"

            if self.arg.is_optional and self.arg.intent != "out":
                lines.append(f"const {ctype} *{self.fortran_call_arg_name} = nullptr;")
                lines.append(f"if ({self.arg.c_name}.has_value()) {{")
                lines.append(f"  {prefix}_to_vec({self.arg.c_name}.value(), {vec_name});")
                lines.append(f"  {self.fortran_call_arg_name} = {vec_name};")
                lines.append("}")
            elif self.arg.intent != "out":
                if self.arg.is_optional:
                    lines.append(f"if ({self.arg.c_name})")
                lines.append(f"{prefix}_to_vec({self.arg.c_name}, {vec_name});")
        else:
            raise NotImplementedError(len(arr))
        return lines

    def call_argument(self) -> str:
        arr = self.arg.array
        if len(arr) == 1:
            if self.arg.intent == "out":
                return f"{self.fortran_call_arg_name}.data()"
            return self.fortran_call_arg_name
        if len(arr) in {2, 3}:
            return f"{self.fortran_call_arg_name}_vec"
        return ""

    def post_call_lines(self) -> list[str]:
        lines = []
        arr = self.arg.array
        if len(arr) in {2, 3} and self.arg.intent in {"out", "inout"}:
            prefix = "matrix" if len(arr) == 2 else "tensor"
            vec_name = f"{self.fortran_call_arg_name}_vec"
            if self.arg.is_optional and self.arg.intent == "inout":
                lines.append(f"if ({self.arg.c_name}.has_value())")
                lines.append(f"  vec_to_{prefix}({vec_name}, {self.arg.c_name}.value());")
            else:
                lines.append(f"vec_to_{prefix}({vec_name}, {self.arg.c_name});")
        return lines

    @property
    def output_arg_name(self) -> str | None:
        if len(self.arg.array) == 1 and self.arg.intent != "out":
            return self.arg.c_name
        if len(self.arg.array) in {2, 3}:
            return self.arg.c_name
        return super().output_arg_name

    def struct_decl(self, ignore_intent: bool = False) -> tuple[str, str] | None:
        if self.arg.intent == "out" or ignore_intent:
            if len(self.arg.array) == 1:
                return self.arg.transform.cpp_declare_type, self.arg.c_name
            if len(self.arg.array) in {2, 3}:
                return self.arg.transform.cpp_type, self.arg.c_name
        return None


@dataclasses.dataclass
class CppWrapperGeneralArgument(CppWrapperArgument):
    def pre_call_lines(self) -> list[str]:
        lines = []
        arg = self.arg
        as_pointer = arg.is_optional or arg.intent in {"out", "inout"} or arg.member.type_info.pointer

        if not as_pointer:
            return []

        if arg.is_optional and arg.intent == "in":
            lvalue = f"{self.arg.c_name}_lvalue"
            lines.append(f"{STANDARD_TYPES[arg.type].c_type} {lvalue};")
            lines.append(f"auto *{self.fortran_call_arg_name} {{ &{lvalue} }};")
            lines.append(f"if ({self.arg.c_name}.has_value()) {{")
            lines.append(f"  {lvalue} = {self.arg.c_name}.value();")
            lines.append("} else {")
            lines.append(f"  {self.fortran_call_arg_name} = nullptr;")
            lines.append("}")

        elif arg.is_optional and arg.intent == "inout":
            if self.arg.transform.is_optional_ref:
                unwrapped = f"&{self.arg.c_name}->get()"
            else:
                unwrapped = f"&{self.arg.c_name}.value()"
            lines.append(
                self.unwrap_optional(
                    unwrapped,
                    type_="auto*",
                    comment="inout, optional",
                )
            )

        elif arg.intent == "out":
            lines.append(f"{STANDARD_TYPES[arg.type].c_type} {self.fortran_call_arg_name}{{}};")

        return lines

    def call_argument(self) -> str:
        arg = self.arg
        as_pointer = arg.is_optional or arg.intent in {"out", "inout"} or arg.member.type_info.pointer

        if not as_pointer:
            return arg.c_name

        if arg.is_optional and arg.intent in {"in", "inout"}:
            return self.fortran_call_arg_name
        if arg.intent in {"in", "inout"}:
            if arg.member.type_info.pointer:
                return "&" + arg.c_name
            return arg.c_name
        if arg.intent == "out":
            if self.arg.transform.is_reference:
                return self.fortran_call_arg_name
            return "&" + self.fortran_call_arg_name
        return arg.c_name

    @property
    def output_arg_name(self) -> str | None:
        if self.arg.intent in {"in", "inout"}:
            if self.arg.member.type_info.pointer:
                return "&" + self.arg.c_name
            return self.arg.c_name
        return super().output_arg_name

    def struct_decl(self, ignore_intent: bool = False) -> tuple[str, str] | None:
        if self.arg.intent == "out" or ignore_intent:
            return STANDARD_TYPES[self.arg.type].c_type, self.arg.c_name
        return None


def generate_routine_return_value_struct(routine: FortranRoutine) -> list[str]:
    cpp_wrapper_args = [CppWrapperArgument.from_arg(arg) for arg in routine.args]
    outputs = [arg for arg in cpp_wrapper_args if arg.arg.intent == "out"]
    if not routine.usable or len(outputs) <= 1:
        return []

    name = snake_to_camel(routine.name)
    lines = [f"struct {name} {{"]

    for output in outputs:
        decl = output.struct_decl()
        if decl:
            type_, name = decl
            lines.append(f"  {type_} {name};")

    lines.append("};")
    return lines


def generate_routine_cpp_wrapper(routine: FortranRoutine) -> list[str]:
    assert routine.docstring is not None
    lines = []

    decl = routine.get_cpp_decl(defaults=False, namespace=True)
    lines.append(decl + " {")

    cpp_wrapper_args = [CppWrapperArgument.from_arg(arg) for arg in routine.args]

    for arg in cpp_wrapper_args:
        lines.extend(arg.pre_call_lines())

    # Call the routine
    call_args = [
        f"/* {arg.arg.transform.cpp_call_fortran_type.strip()} */ {arg.call_argument()}"
        for arg in cpp_wrapper_args
    ]
    lines.append(f"  fortran_{routine.name}(" + ", ".join(call_args) + ");")

    for arg in cpp_wrapper_args:
        lines.extend(arg.post_call_lines())

    outputs = [arg for arg in cpp_wrapper_args if arg.arg.intent == "out"]
    if not outputs:
        pass

    elif len(outputs) == 1:
        (output,) = outputs
        return_expr = output.output_arg_name
        if return_expr and isinstance(
            output, (CppWrapperTypeArgument, CppWrapperTypeArgumentArray, CppWrapperStringArgumentArray)
        ):
            return_expr = f"std::move({return_expr})"

        if return_expr:
            lines.append(f"return {return_expr};")
    else:
        output_vals = []
        for arg in outputs:
            val = arg.output_arg_name
            if not val:
                continue

            if isinstance(
                arg, (CppWrapperTypeArgument, CppWrapperTypeArgumentArray, CppWrapperStringArgumentArray)
            ):
                val = f"std::move({val})"

            output_vals.append(val)

        output_struct_name = snake_to_camel(routine.name)
        lines.append(f"return {output_struct_name}{{" + ", ".join(output_vals) + "};")

    lines.append("}")
    return lines


def generate_routines_header(
    template: str,
    routines: dict[str, FortranRoutine],
    settings: RoutineSettings,
) -> str:
    forward_decls = []
    for routine in sorted_routines(routines):
        if not routine.usable:
            forward_decls.append(f"\n// Skipped unusable routine {routine.name}:")
            for reason in "\n".join(routine.unusable_reason).splitlines():
                forward_decls.append(f"// - {reason}")
        else:
            forward_decls.append(routine.fortran_forward_declaration)
            forward_decls.extend(generate_routine_return_value_struct(routine))

            decl = routine.get_cpp_decl(defaults=True, namespace=False)
            forward_decls.append(decl + ";")

    tpl = string.Template(re.sub("^// ", "", template, flags=re.MULTILINE))
    return tpl.substitute(forward_declarations="\n".join(forward_decls), namespace=settings.cpp_namespace)


def generate_to_string_header(
    template: str, structs: list[CodegenStructure], routines: dict[str, FortranRoutine]
) -> str:
    decls = [f"  std::string to_string(const {struct.cpp_class}& self);" for struct in structs]
    for routine in sorted_routines(routines):
        if routine.usable and len(routine.outputs) > 1:
            decls.append(f"  std::string to_string(const {routine.cpp_return_type}& self);")

    tpl = string.Template(template.replace("// ${", "${"))
    return tpl.substitute(decls="\n".join(decls))


repr_kwargs_by_struct = {
    "lat_struct": {
        "use_name": "self.use_name()",
        "#branch": "to_string(self.branch().size())",
    },
    "ele_struct": {
        "name": "self.name()",
        "ix_branch": "to_string(self.ix_branch())",
        "ix_ele": "to_string(self.ix_ele())",
    },
}


def arg_to_cpp_string(arg: Argument, use_call: bool = True) -> str | None:
    if arg.full_type not in proxy_templates:
        # print("TODO:", arg.f_name, arg.full_type)
        return None

    if "parent" in arg.comment.lower():
        return '"..."'

    if use_call:
        arg_name = f"self.{arg.c_name}()"
    else:
        arg_name = f"self.{arg.c_name}"

    # if arg.pointer_type != "NOT":
    #     return f"to_string(self.{arg.c_name}())"
    # elif arg.array:
    #     return '"[...]"'
    if arg.full_type.type == "type":
        if len(arg.array) == 0:
            return f"to_string({arg_name})"
        return '"[...]"'
    if arg.full_type.type == "character":
        if arg.array:
            return f"to_string({arg_name})"
        return arg_name
    # elif arg.full_type.type == "complex":
    #     return '"cplxTODO"'

    return f"to_string({arg_name})"


def generate_cpp_to_string(ptr: str, cpp_class_name: str, kwarg_pairs: dict[str, str | None]) -> list[str]:
    formatted_kwargs = ", ".join(
        f'std::pair{{"{key}", {value}}}' for key, value in kwarg_pairs.items() if value is not None
    )

    return [
        f"std::string to_string(const {cpp_class_name}& self) {{",
        f'    return repr({ptr}, "{cpp_class_name}", {{ {formatted_kwargs} }});',
        "}",
    ]


def generate_to_string_code(
    template: str, structs: list[CodegenStructure], routines: dict[str, FortranRoutine]
) -> str:
    code_lines = []
    for struct in structs:
        if struct.f_name not in repr_kwargs_by_struct:
            kw = {}
            repr_kwargs_by_struct[struct.f_name] = kw
            for arg in struct.arg:
                if not arg.is_component:
                    continue
                as_str = arg_to_cpp_string(arg)
                if as_str:
                    kw[arg.c_name] = as_str

        code_lines.extend(
            generate_cpp_to_string(
                ptr="self.get_fortran_ptr()",
                cpp_class_name=struct.cpp_class,
                kwarg_pairs=repr_kwargs_by_struct.get(struct.f_name, {}),
            )
        )

    for routine in sorted_routines(routines):
        if routine.usable and len(routine.outputs) > 1:
            code_lines.extend(
                generate_cpp_to_string(
                    ptr="&self",
                    cpp_class_name=routine.cpp_return_type,
                    kwarg_pairs={
                        output.c_name: arg_to_cpp_string(output, use_call=False) for output in routine.outputs
                    },
                )
            )

    tpl = string.Template(template.replace("// ${", "${"))
    return tpl.substitute(decls="\n".join(code_lines))
