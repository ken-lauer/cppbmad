from __future__ import annotations

import dataclasses
import logging
import os
import pathlib
import re
import string
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Literal

from .arg import Argument as InterfaceArgument
from .context import CodegenConfig, RoutineSettings, config_context, get_params
from .cpp import CppWrapperArgument, generate_routine_cpp_wrapper, generate_routines_header
from .docstring import DocstringParameter, RoutineDocstring, parse_routine_comment_block
from .enums import EnumValue, replace_enums_with_cpp
from .exceptions import RoutineNotFoundError, UnsupportedTypeError
from .fortran import generate_fortran_routine_with_c_binding
from .paths import CODEGEN_ROOT, CPPBMAD_ROOT
from .structs import (
    FileLine,
    ParsedStructure,
    StructureMember,
    TypeInformation,
    join_ampersand_lines,
    parse_declaration,
    split_comment,
)
from .types import ArgumentType, Intent, RoutineType, get_type_transform
from .util import (
    snake_to_camel,
    sorted_routines,
    struct_to_proxy_class_name,
    wrap_line,
)

logger = logging.getLogger(__name__)

GLOBAL_INTERFACE = "(global)"
TEST_BUILD = int(os.environ.get("PYBMAD_TEST_BUILD", "0"))
TEST_ROUTINES = {"bmad_parser"}
docstring_hotfixes = {}

docstring_hotfixes["init_coord1"] = [
    DocstringParameter(
        name="orb",
        description="Input orbit",
        data_type="coord_struct",
        is_input=True,
        is_output=True,
    )
]
docstring_hotfixes["init_coord2"] = [
    DocstringParameter(
        name="orb_in",
        description="Input orbit",
        data_type="coord_struct",
        is_input=True,
        is_output=True,
    ),
    DocstringParameter(
        name="orb_out",
        description="Initialized coordinate",
        data_type="coord_struct",
        is_input=False,
        is_output=True,
    ),
]
docstring_hotfixes["init_coord3"] = [
    DocstringParameter(
        name="orb",
        description="Input orbit",
        data_type="coord_struct",
        is_input=True,
        is_output=True,
    ),
]
docstring_hotfixes["bbu_track_all"] = [
    DocstringParameter(
        name="hom_voltage_normalized",
        description="HOM voltage normalized",
        data_type="real(rp)",
        is_input=False,
        is_output=True,
    ),
    DocstringParameter(
        name="growth_rate",
        description="Growth rate",
        data_type="real(rp)",
        is_input=False,
        is_output=True,
    ),
    DocstringParameter(
        name="lost",
        description="Lost",
        data_type="logical",
        is_input=False,
        is_output=True,
    ),
    DocstringParameter(
        name="irep",
        description="irep",
        data_type="integer",
        is_input=False,
        is_output=True,
    ),
]


def normalize_intent(
    typ: TypeInformation,
    arg_type: ArgumentType,
    doc_intent: Intent | None = None,
    ndim: int = 0,
    is_dynamic_array: bool = False,
    guessed: bool = False,
) -> Intent:
    if doc_intent:
        intent = doc_intent
    elif not typ.intent:
        intent = "inout"
    else:
        intent = typ.intent.replace(" ", "").lower()

    if intent not in ("in", "inout", "out"):
        raise ValueError(f"Unsupported intent: {intent}")

    if is_dynamic_array and intent == "out" and not typ.allocatable:
        # The user will need to preallocate this for us - allow them to do
        # it by making it 'inout'
        return "inout"
        # member.type_info = member.type_info.replace(allocatable=True)

    if (
        guessed  # 'inout' is default guess
        and intent == "inout"
        and not typ.intent
        and ndim == 0
        and arg_type
        in (
            "real",
            "real16",
            "complex",
            "integer",
            "integer8",
            "logical",
            "character",
            # "type",
        )
    ):
        # For scalar values, assume it's not by reference
        intent = "in"

    return intent


@dataclasses.dataclass
class _InterfaceOverload:
    name: str
    docstring: list[str]


@dataclasses.dataclass
class RoutineArg(InterfaceArgument):
    intent: Intent = "inout"
    description: str = ""
    doc_data_type: str | None = None
    doc_is_optional: bool = False

    @property
    def is_optional(self):
        return self.doc_is_optional or self.member.type_info.optional

    @property
    def is_python_immutable(self):
        if self.intent != "inout":
            return False

        return len(self.array) == 0 and self.type != "type"

    @property
    def is_input(self):
        return self.intent in {"in", "inout"}

    @property
    def is_output(self):
        return self.intent == "out"

    @property
    def kind_as_cpp_class(self) -> str:
        kind = self.kind
        return struct_to_proxy_class_name(kind)

    @property
    def c_to_fortran_decl(self) -> str:
        return f"{self.transform.cpp_call_fortran_type} {self.c_name}"

    @property
    def cpp(self):
        return CppWrapperArgument.from_arg(self)

    def _check_enums_if_dynamic(self, enums: dict[str, EnumValue]):
        def is_integer(s: str) -> bool:
            try:
                int(s)
                return True
            except ValueError:
                return False

        for dim in self.array:
            if ":" in dim:
                parts = dim.split(":")
                for part in parts:
                    if not part:
                        continue
                    if is_integer(part):
                        continue
                    if part.upper().rstrip("$") in enums:
                        continue
                    return True
            elif dim == "*":
                return True
            else:
                part = dim
                if not is_integer(part) and part.upper().rstrip("$") not in enums:
                    return True

        return False

    def maybe_make_dynamic_array(self):
        # if self._check_enums_if_dynamic(ctx().enums_by_name):
        #     self.orig_array = self.array
        #     self.array = [":" for _ in self.array]
        pass

    @classmethod
    def from_routine(
        cls,
        routine: FortranRoutine,
        member: StructureMember,
        doc: DocstringParameter,
        params: CodegenConfig | None = None,
    ):
        if params is None:
            params = get_params()

        arg = cls.from_fstruct(routine, member, params)

        # TODO some confusion:
        #   c_name is the argument name C++ callers see on the Fortran side
        #   f_name is the internal name used when converting C->F
        #     f_name then becomes what we use to call the original function
        arg.f_name = f"f_{member.name}"

        f_to_c_name = params.c_side_name_translation
        arg.c_name = f_to_c_name.get(member.name, member.name)

        arg.description = doc.description
        arg.doc_data_type = doc.data_type
        arg.doc_is_optional = doc.is_optional

        if arg.type in {"real16", "logical"} and len(arg.array):
            # Sorry, we need to use allocatable for this as the Fortran
            # representation doesn't match up with the C++ representation
            # logical: bitset/1 byte vs 4 bytes
            # real16:  ...
            #
            # This lets us get around copy conversions, with the downside of
            # making API usage more inconvenient
            arg.pointer_type = "ALLOC"
            arg.member.type_info = arg.member.type_info.replace(allocatable=True)

        arg.intent = normalize_intent(
            arg.member.type_info,
            doc_intent=None if doc.guessed else doc.intent,
            ndim=len(arg.array),
            is_dynamic_array=arg.is_dynamic_array,
            arg_type=arg.type,
            guessed=doc.guessed,
        )

        if arg.intent == "out" and arg.member.type_info.pointer:
            # if we get back a null ptr, we need to be able to give std::nullopt
            arg.member.type_info = arg.member.type_info.replace(optional=True)

        return arg

    @property
    def transform(self):
        tf = get_type_transform(
            self.full_type,
            intent=self.intent,
            is_optional=self.is_optional,
            kind=self.kind,
            is_dynamic_array=self.is_dynamic_array,
        )
        return self._replace_transform_placeholders(tf)

    @property
    def c_class(self) -> str:
        if self.type != "type":
            raise ValueError(f"Not a type: {self.type}")
        return struct_to_proxy_class_name(self.member.kind)


ProcedureType = Literal[
    "subroutine",
    "function",
    "elemental_function",
    "pure_function",
    "recursive_function",
    "module_function",
    "module_subroutine",
    "pure_subroutine",
    "elemental_subroutine",
]


def _get_docstring_arg(lower_doc: dict[str, DocstringParameter], arg_name: str):
    return lower_doc.get(
        arg_name.lower(),
        DocstringParameter(
            name=arg_name,
            description="Undocumented",
            is_input=True,
            is_output=True,
            is_optional=False,
            guessed=True,
        ),
    )


@dataclass
class FortranRoutine:
    """Represents a Fortran procedure (subroutine or function)."""

    filename: pathlib.Path
    name: str
    proc_type: ProcedureType
    start_line: FileLine
    end_line: FileLine | None = None
    doc_comment: list[str] = field(default_factory=list)
    declared_argument_list: list[str] = field(default_factory=list)
    declarations: dict[str, StructureMember] = field(default_factory=dict)
    body: list[FileLine] = field(default_factory=list)
    result_name: str | None = None
    docstring: RoutineDocstring | None = None
    interface: str | None = None
    module: str | None = None
    private: bool = False
    cpp_namespace: str = ""
    overload_disabled: bool = False
    # Arguments - in order and including result value (if applicable)
    args: list[RoutineArg] = field(default_factory=list)

    @property
    def is_function(self) -> bool:
        """True for any function variant (plain, elemental, pure, recursive, module)."""
        return self.proc_type.endswith("function")

    @property
    def overloaded_name(self) -> str:
        """Some routines have interfaces which result in overloaded functions in C++"""
        if not self.interface or self.interface == GLOBAL_INTERFACE:
            return self.name
        if self.overload_disabled:
            return self.name
        return self.interface

    @property
    def arg_names_with_result(self) -> list[str]:
        if self.result_name:
            result_arg = [self.result_name.lower()]
        else:
            result_arg = []

        return [*self.declared_argument_list, *result_arg]

    @property
    def wrapper_args(self):
        required = []
        optional = []
        for arg in self.args:
            if arg.is_optional:
                optional.append(arg)
            else:
                required.append(arg)
        return [*required, *optional]

    def parse(self, *, config: CodegenConfig) -> None:
        """
        Parses the docstring comment, matches up arguments, and creates all translated args.
        """
        self.docstring = self._parse_comment()
        if self.is_function and not self.result_name:
            self.result_name = self.name.lower()
        self.declarations = self._parse_argument_types()
        if not self.docstring:
            self.docstring = RoutineDocstring(
                name=self.name,
                filename=self.filename,
                lineno=-1,
                routine_type=RoutineType.UNKNOWN,
                description=[],
                params=[
                    DocstringParameter(
                        name=arg_name, is_input=True, is_output=True, is_optional=False, guessed=True
                    )
                    for arg_name in self.arg_names_with_result
                ],
                result_variable=None,
                is_overloaded=False,
                overloaded_versions=[],
                related_routines=[],
                notes=[],
            )
        for arg in docstring_hotfixes.get(self.name.lower(), []):
            self.docstring.update_parameter(arg)

        if self.is_function and self.result_name == self.name.lower():
            assert self.result_name is not None
            retval = self.declarations.pop(self.result_name)
            self.declarations["func_retval__"] = retval
            self.result_name = "func_retval__"
            retval.name = "func_retval__"

        self.args = self.translate_args(config)

    def _parse_argument_types(self) -> dict[str, StructureMember]:
        skips = {
            "private",
            "sequence",  # ?
            "contains",  # ?
            "procedure next_in_branch",  # ?
            "call",
        }
        # Brute force it!
        arguments = {}
        expected_args = {arg.lower() for arg in self.arg_names_with_result}
        for body_line in self.body[1:]:
            line, comment = body_line.split_comment("!")
            if not line or line.lower() in skips:
                continue

            start = line.split()[0].lower()
            if start in skips:
                continue
            try:
                for decl in parse_declaration(line):
                    arg = decl.name.lower()
                    if arg not in expected_args or arg in arguments:
                        continue
                    arguments[arg] = StructureMember(
                        name=decl.name,
                        type_info=decl.type,
                        line=body_line.lineno,
                        definition=line,
                        comment=comment,
                        default=decl.default,
                    )
            except Exception:
                # logger.warning("Line not a member? %s -> %s", line, ex)
                continue

        for arg in self.arg_names_with_result:
            assert arg.lower() in arguments, arg
        return arguments

    def translate_args(self, config: CodegenConfig):
        assert self.docstring is not None
        lower_members = {name.lower(): member for name, member in self.declarations.items()}
        lower_doc = {name.lower(): doc for name, doc in self.docstring.arguments_by_name.items()}

        args = [
            RoutineArg.from_routine(
                self,
                lower_members[arg_name.lower()],
                _get_docstring_arg(lower_doc, arg_name),
                params=config,
            )
            for arg_name in self.arg_names_with_result
        ]
        for arg in args:
            # smoke check that we can transform all args
            _ = arg.transform

        if self.is_function:
            args_by_c_name = {arg.c_name: arg for arg in args}
            if self.result_name in args_by_c_name:
                args_by_c_name[self.result_name].intent = "out"
            else:
                logger.error(f"Function result not in arg list? {list(args_by_c_name)}")

        return args

    def _parse_comment(self) -> RoutineDocstring | None:
        return parse_routine_comment_block(
            self.filename,
            "\n".join(self.doc_comment),
            starting_lineno=self.start_line.lineno,
        )

    def __str__(self) -> str:
        doc = "\n".join(self.doc_comment) if self.doc_comment else "No documentation"
        return f"{self.proc_type.upper()} {self.name}\n{doc}\nArguments: {', '.join(self.declared_argument_list)}"

    def _get_cpp_decl_spec(self, allow_defaults: bool) -> list[str]:
        # only a continuous block at the end can have defaults
        specs = []

        for arg in reversed(self.wrapper_args):
            if not arg.is_input:
                continue
            if allow_defaults and arg.transform.cpp_default:
                default_str = arg.transform.cpp_default
            else:
                default_str = ""
                allow_defaults = False

            spec = f"{arg.cpp.cpp_decl}{default_str}"
            specs.append(spec)

        return list(reversed(specs))

    @property
    def cpp_return_type(self) -> str:
        outputs = self.outputs
        if len(outputs) == 1:
            (output,) = outputs
            # return_type = output.transform.cpp_return_type
            return_type, _ = output.cpp.struct_decl()
        elif len(outputs) == 0:
            return_type = "void"
        else:
            return_type = snake_to_camel(self.name)
            if self.cpp_namespace:
                return "::".join((self.cpp_namespace, return_type))
        return return_type

    @property
    def python_class_return_type(self) -> tuple[str, str]:
        immut_args = [arg for arg in self.args if arg.is_python_immutable]
        outputs = [*self.outputs, *immut_args]

        if len(outputs) <= 1 and not immut_args:
            # TODO: immut args even for 1 -> struct
            raise ValueError("Returns native type")

        clsname = snake_to_camel(self.name)
        if immut_args:
            full_clsname = f"Py{clsname}"
        else:
            full_clsname = self.cpp_return_type
            clsname = full_clsname.split("::")[1]
        return clsname, full_clsname

    def get_cpp_decl(self, defaults: bool, namespace: bool) -> str:
        decl_args = ", ".join(self._get_cpp_decl_spec(defaults))

        # name = self.overloaded_name if not self.overload_disabled else self.name
        routine_and_args = f"{self.overloaded_name}({decl_args})"

        if self.cpp_namespace and namespace:
            routine_and_args = "::".join((self.cpp_namespace, routine_and_args))

        return f"{self.cpp_return_type} {routine_and_args}"

    @property
    def outputs(self) -> list[RoutineArg]:
        return [arg for arg in self.args if arg.intent == "out"]

    @property
    def args_by_c_name(self):
        return {arg.c_name: arg for arg in self.args}

    @property
    def args_by_c_name_lower(self):
        return {arg.c_name.lower(): arg for arg in self.args}

    @property
    def usable(self) -> bool:
        """Usable routine for pybmad?"""
        return not self.unusable_reason

    @property
    def unusable_reason(self) -> list[str]:
        # if self.interface:
        #     return ["Interface definition"]
        if self.private:
            return ["Marked as private"]
        if self.module is None:
            return ["Module name unset"]
        if self.docstring is None:
            return ["No matching docstring"]
        if TEST_BUILD and self.name not in TEST_ROUTINES:
            return [f"PYBMAD_TEST_BUILD enabled; {self.name} is skipped"]
        ctx = config_context.get()

        if ctx.params.should_skip(self.name.lower()):
            return ["Routine in configuration skip list"]
        if ctx.params.should_skip(self.module.lower()):
            return ["Routine module in configuration skip list"]
        if self.module.lower() in ctx.params.skips:
            return [f"Routine module ({self.module}) in configuration skip list"]
        structs_by_name = ctx.codegen_structs_by_name

        assert self.docstring is not None
        lower_members = {name.lower(): member for name, member in self.declarations.items()}
        lower_doc = {name.lower(): doc for name, doc in self.docstring.arguments_by_name.items()}

        enums_by_name = ctx.enums_by_name
        reasons = []

        for arg_name in self.arg_names_with_result:
            try:
                member = lower_members[arg_name.lower()]
            except KeyError:
                reasons.append(f"Argument not defined: {arg_name} (have: {list(lower_members)})")
                continue

            doc_arg = _get_docstring_arg(lower_doc, arg_name)

            try:
                arg = RoutineArg.from_routine(self, member, doc_arg)
            except Exception as ex:
                reasons.append(f"Exception for argument {member.name}: {ex.__class__.__name__} {ex}")
                if "debug" in str(ex):
                    raise
                logger.warning("RoutineArg.from_routine failed", exc_info=True)
            else:
                if arg.type == "type" and arg.kind.lower() not in structs_by_name:
                    reasons.append(f"Untranslated type: {arg.kind.lower()} ({len(arg.array)}D)")
                if arg.type == "type" and len(arg.array) > 1:
                    reasons.append(f"TODO type arrays: {arg.c_class} {arg.intent=} {arg.full_type}")
                if len(arg.array) > 1 and arg.type in {"character", "logical"}:
                    arr = ",".join(arg.array)
                    reasons.append(
                        f"{len(arg.array)}D array handling not supported for {arg.type}: {arg.c_name}({arr}) {arg.full_type}"
                    )
                if ":" in arg.array or "0:" in arg.array or "*" in arg.array:
                    arr = ",".join(arg.array)
                    if len(arg.array) > 1:
                        if arg.member.type_info.pointer:
                            reasons.append(f"Pointer to variable {arg.intent} sized array: {arg.full_type}")
                        if arg.type != "type":
                            reasons.append(f"Variable {arg.intent} sized array: {arg.full_type}")

                try:
                    replace_enums_with_cpp(arg.lbound, enums_by_name)
                    replace_enums_with_cpp(arg.ubound, enums_by_name)
                except Exception as ex:
                    reasons.append(f"Array bounds handling: {ex}")

        if len(self.arg_names_with_result) != len(self.args):
            reasons.append("Translated arg count mismatch (unsupported?)")

        return reasons

    @property
    def needs_python_wrapper(self) -> bool:
        return any(arg.is_python_immutable for arg in self.args)

    @property
    def fortran_param_spec(self) -> list[str]:
        return [f"{arg.c_to_fortran_decl} /* {arg.full_type} {arg.intent} */" for arg in self.args]

    @property
    def fortran_forward_declaration(self):
        return_type = "bool" if self.result_name else "void"
        param_spec = ", ".join(self.fortran_param_spec)
        return f'extern "C" {return_type} fortran_{self.name}({param_spec});'


class FortranParser:
    """Parser for Fortran source files."""

    # Regular expressions for identifying Fortran constructs
    COMMENT_PATTERN = re.compile(r"^\s*[!c*]", re.IGNORECASE)

    # Match different procedure types
    PROCEDURE_START_PATTERN = re.compile(
        r"^\s*(?:(recursive|pure|elemental|module)\s+)?"
        + r"((?:subroutine|function))\s+([a-z0-9_]+)(?:\s*\((.*?)\))?\s*(?:result\s*\(\s*([a-z0-9_]+)\s*\))?.*$",
        re.IGNORECASE,
    )

    # Match 'end subroutine' or 'end function' statements
    PROCEDURE_END_PATTERN = re.compile(
        r"^\s*end\s+(subroutine|function)(?:\s+([a-z0-9_]+))?\s*$", re.IGNORECASE
    )

    def __init__(self):
        """Initialize the Fortran parser."""
        self.procedures: list[FortranRoutine] = []

    def _is_comment(self, line: str) -> bool:
        """Check if a line is a comment."""
        return bool(self.COMMENT_PATTERN.match(line))

    def _extract_comment_content(self, line: str) -> str:
        """Extract the content of a comment line."""
        match = self.COMMENT_PATTERN.match(line)
        if match:
            # Remove the comment character and any leading/trailing whitespace
            return line[match.end() :].strip()
        return ""

    def _extract_procedure_info(
        self, line: FileLine
    ) -> tuple[str, ProcedureType, list[str], str | None] | None:
        """Extract procedure name, type and arguments from a procedure declaration line."""
        match = self.PROCEDURE_START_PATTERN.match(line.line)
        if not match:
            return None

        modifier, proc_type, name, args_str, result_name = match.groups()

        # Determine the full procedure type
        if modifier:
            modifier = modifier.lower()
            proc_type = proc_type.lower()
            full_type: ProcedureType = f"{modifier}_{proc_type}"
        else:
            full_type = proc_type.lower()

        # Parse the argument list
        args = []
        if args_str:
            args = [arg.strip() for arg in args_str.split(",")]

        return name, full_type, args, result_name

    def parse_file(self, file_path: pathlib.Path) -> list[FortranRoutine]:
        """Parse a Fortran file and extract procedures."""
        lines = FileLine.from_file(file_path)
        return self.parse_lines(lines)

    def parse_lines(self, lines: Sequence[FileLine]) -> list[FortranRoutine]:
        """Parse a sequence of FileLine objects and extract procedures."""
        self.procedures = []
        pending_comments: list[str] = []
        current_procedure: FortranRoutine | None = None
        in_interface: str | None = None
        in_module: str | None = None
        interface_doc: list[str] = []
        interface_overloads = {}

        # Track nesting level and procedure stack
        nesting_level = 0
        procedure_stack: list[tuple[FortranRoutine, bool]] = []  # (procedure, in_contains)
        marked_private = set()

        for line in lines:
            stripped_line = line.line.strip()

            # Skip completely empty lines when not in procedure
            if not stripped_line and not current_procedure:
                continue

            # Check if this is a comment
            if self._is_comment(stripped_line):
                if current_procedure is None:
                    # This is a comment before a procedure, save it for later
                    pending_comments.append(stripped_line)
                else:
                    # This is a comment inside a procedure, add it to the body
                    current_procedure.body.append(line)
                continue

            if stripped_line.lower().startswith("interface"):
                parts = stripped_line.lower().split()
                if len(parts) > 1 and parts[1] not in {"assignment", "operator"}:
                    in_interface = parts[1]
                else:
                    in_interface = GLOBAL_INTERFACE
                interface_doc = list(pending_comments)
                pending_comments = []
            elif stripped_line.lower().startswith("end interface"):
                in_interface = None
                interface_doc = []

            pre_comment, _comment = split_comment(stripped_line.lower())
            pre_comment_parts = pre_comment.split()

            if pre_comment_parts:
                if pre_comment_parts[0] == "module":
                    if len(pre_comment_parts) == 2:
                        in_module = pre_comment_parts[1]
                elif pre_comment_parts == ["end", "module"]:
                    in_module = None

                if pre_comment_parts[0] == "private":
                    for name in pre_comment_parts[1:]:
                        marked_private.add(name.strip(", "))

                if in_interface and pre_comment_parts[:2] == ["module", "procedure"]:
                    interface_overloads[pre_comment_parts[2].lower()] = _InterfaceOverload(
                        name=in_interface, docstring=list(interface_doc)
                    )

            # Handle the "contains" keyword
            if stripped_line.lower() == "contains" and current_procedure:
                # Mark that we're now in the contains section of this procedure
                current_procedure.body.append(line)

                # Update the current procedure's in_contains status in the stack
                if procedure_stack and procedure_stack[-1][0] == current_procedure:
                    procedure_stack[-1] = (current_procedure, True)
                continue

            # Check if this is the start of a procedure
            proc_info = self._extract_procedure_info(line)
            if proc_info:
                name, proc_type, args, result_name = proc_info

                # Check if we're inside another procedure's contains section
                if current_procedure and procedure_stack and procedure_stack[-1][1]:
                    # This is a nested procedure, just add it to parent's body
                    current_procedure.body.append(line)
                    # Increment nesting level
                    nesting_level += 1
                    continue

                if in_interface and len(interface_doc) > len(pending_comments):
                    pending_comments = list(interface_doc)

                # Start a new top-level procedure
                current_procedure = FortranRoutine(
                    filename=line.filename,
                    name=name.lower(),
                    proc_type=proc_type,
                    start_line=line,
                    doc_comment=list(pending_comments),
                    declared_argument_list=args,
                    body=[line],
                    result_name=result_name,
                    interface=in_interface,
                    module=in_module,
                    private=name.lower() in marked_private or len(procedure_stack) > 0,
                )

                # Add to procedures list and put on stack with in_contains=False
                self.procedures.append(current_procedure)
                procedure_stack.append((current_procedure, False))
                pending_comments = []  # Clear pending comments
                continue

            # Check if this is the end of a procedure
            end_match = self.PROCEDURE_END_PATTERN.match(stripped_line)
            if end_match and current_procedure:
                # Always add the line to the current procedure's body
                current_procedure.body.append(line)

                if nesting_level > 0:
                    # This is closing a nested procedure, not the main one
                    nesting_level -= 1
                    continue

                # This is closing the current top-level procedure
                _end_proc_type, end_proc_name = end_match.groups()

                # End current procedure if no name specified or names match
                if not end_proc_name or end_proc_name.lower() == current_procedure.name.lower():
                    current_procedure.end_line = line
                    procedure_stack.pop()  # Remove from stack

                    # Set current procedure to previous one in the stack if available
                    if procedure_stack:
                        current_procedure = procedure_stack[-1][0]
                    else:
                        current_procedure = None
                continue

            # If we're inside a procedure, add this line to its body
            if current_procedure:
                current_procedure.body.append(line)
                continue

            # If we reached here, this line is neither a comment nor part of a procedure
            # Clear any pending comments if we encounter non-comment code
            pending_comments = []

        for proc in self.procedures:
            if proc.name.lower() in interface_overloads:
                ovd = interface_overloads[proc.name.lower()]
                proc.interface = ovd.name
                proc.private = ovd.name in marked_private  # is this a thing?
                # if len(ovd.docstring) >= len(proc.doc_comment):
                proc.doc_comment = list(ovd.docstring)
            elif proc.name.lower() in marked_private:
                proc.private = True

        return self.procedures


def parse_fortran_file_for_routines(file_path: pathlib.Path) -> list[FortranRoutine]:
    """Convenience function to parse a Fortran file."""
    parser = FortranParser()
    return parser.parse_file(file_path)


def procedures_from_source_directory(
    source_path: pathlib.Path, skip_files: set[str], skip_procedures: set[str]
) -> list[FortranRoutine]:
    procedures: list[FortranRoutine] = []
    for source_fn in source_path.rglob("*.f90", case_sensitive=False):
        if source_fn.name.lower() in skip_files:
            continue
        for proc in parse_fortran_file_for_routines(source_fn):
            if proc.name.lower() in skip_procedures:
                logger.debug("User-skipped procedure: %s from %s", proc.name, source_fn)
            else:
                procedures.append(proc)
    return procedures


def parse_bmad_routines(settings: RoutineSettings, config: CodegenConfig):
    procedures: list[FortranRoutine] = []
    for source_path in settings.source_paths:
        for proc in procedures_from_source_directory(
            source_path,
            skip_files=settings.skip_files,
            skip_procedures=settings.skip_procedures,
        ):
            proc.cpp_namespace = settings.cpp_namespace
            proc.overload_disabled = bool(
                proc.interface and proc.interface.lower() in settings.do_not_overload
            )
            procedures.append(proc)

    for proc in procedures:
        try:
            proc.parse(config=config)
        except UnsupportedTypeError as ex:
            logger.warning("Parsing failed for routine %r: unsupported type (%s)", proc.name, ex)
        except Exception as ex:
            logger.warning("Parsing failed for routine %r: %s", proc.name, ex)
            logger.debug("Parsing failed for routine %r", proc.name, exc_info=True)

    return procedures


def prune_routines(procedures: list[FortranRoutine], config: CodegenConfig):
    unusable_procs = [proc for proc in procedures if not proc.usable]
    usable_procs = [proc for proc in procedures if proc.usable]
    usable_by_name: dict[str, FortranRoutine] = {proc.name: proc for proc in usable_procs}

    for unusable in unusable_procs:
        others = [other for other in unusable_procs if other.name == unusable.name and other is not unusable]
        usable = [other for other in usable_procs if other.name == unusable.name]
        logger.debug(
            "%d instance(s) of procedure %s found. %d unusable and %d usable",
            len(others) + len(usable) + 1,
            unusable.name,
            len(others) + 1,
            len(usable),
        )
        for idx, option in enumerate([unusable, *others, *usable], start=1):
            logger.debug(
                "%s option #%d %s:%d: %s\n%s",
                option.name,
                idx,
                option.filename.name,
                option.start_line.lineno,
                "usable" if option.usable else "unusable",
                "\n   ".join(option.unusable_reason),
            )

    for unusable in unusable_procs:
        options = [other for other in usable_procs if other.name == unusable.name]
        usable = usable_by_name.get(unusable.name, None)
        if usable is None:
            continue
        for option in options:
            option.module = option.module or unusable.module

    def usability_metric(proc: FortranRoutine):
        # if proc.interface is not None:
        #     return -3
        if proc.docstring is None:
            return -2
        if not proc.usable:
            return -1
        return len("\n".join(proc.doc_comment))

    def get_interface_defn() -> FortranRoutine | None:
        for opt in options:
            if opt.interface is not None:
                return opt
        return None

    proc_by_name: dict[str, FortranRoutine] = {}
    for usable in procedures:
        if usable.name in proc_by_name:
            continue

        options = [other for other in procedures if other.name == usable.name]

        best_option = sorted(options, key=usability_metric)[-1]
        if len(options) > 1:
            logger.debug(
                "Choosing %s option %s:%d with usability metric %d",
                best_option.name,
                best_option.filename.name,
                best_option.start_line.lineno,
                usability_metric(best_option),
            )

        if best_option.docstring:
            missing_arg_names = [
                name
                for name in best_option.arg_names_with_result
                if name.lower() not in best_option.docstring.arguments_by_name
                or best_option.docstring.arguments_by_name[name.lower()].guessed
            ]

            for missing_name in missing_arg_names:
                for other in options:
                    if other is best_option:
                        continue
                    if not other.docstring:
                        continue
                    try:
                        other_doc = other.docstring.arguments_by_name[missing_name.lower()]
                    except KeyError:
                        continue
                    if other_doc.guessed:
                        continue

                    best_option.docstring.update_parameter(other_doc)

            if missing_arg_names:
                try:
                    best_option.args = best_option.translate_args(config)
                except Exception as ex:
                    if "debug" in str(ex):
                        raise
                    logger.debug(
                        f"Reparse failed after docstring inclusion: {best_option.name} {missing_arg_names} {ex}"
                    )

        intf = get_interface_defn()
        if intf:
            best_option.module = intf.module

        proc_by_name[usable.name] = best_option
        if best_option.docstring:
            fix_routine(best_option, best_option.docstring)

    return proc_by_name


def parse_bmad_routine_file(
    fortran_code: list[FileLine],
) -> dict[tuple[str, tuple[str, ...]], list[FileLine]]:
    """
    Parse a Fortran file containing subroutines and return a dictionary
    mapping subroutine names to their content.
    """
    routines = {}

    lines = join_ampersand_lines(fortran_code)

    current_subroutine = None
    current_content = []
    have_contains = False

    for line in lines:
        lower = line.line.lower()
        if lower.startswith(
            (
                "subroutine ",
                "recursive subroutine ",
                "function ",
                "recursive function",
            )
        ):
            if lower.startswith("recursive "):
                lower = lower.removeprefix("recursive ")
            # If we were already collecting a subroutine, save it before starting a new one
            if current_subroutine:
                routines[current_subroutine] = current_content

            subroutine_name = lower.split()[1].split("(")[0]
            arguments = tuple(
                arg.strip() for arg in lower.split("(")[1].split(")")[0].split(",") if arg.strip()
            )

            current_subroutine = (subroutine_name, arguments)
            current_content = [line]
        elif lower.startswith(("end subroutine", "end function")):
            # current_content.append(line)
            if not have_contains:
                assert current_subroutine is not None

            if current_subroutine is not None:
                routines[current_subroutine] = current_content
                current_subroutine = None
                current_content = []
        elif lower == "contains" or lower.startswith("contains "):
            have_contains = True
            if current_subroutine is not None:
                routines[current_subroutine] = current_content
                current_subroutine = None
                current_content = []

            break  # TODO: this may be only for bmad_routine_interface.f90?
        elif current_subroutine:
            current_content.append(line)

    # In case there's a final subroutine without an explicit end
    if current_subroutine:
        routines[current_subroutine] = current_content

    return routines


def find_renames(filename: pathlib.Path, routine: FortranRoutine):
    lines = FileLine.from_file(filename)
    for name, args in parse_bmad_routine_file(lines):
        if routine.name == name:
            return args, routine.result_name
    raise RoutineNotFoundError(routine.name)


def fix_routine(routine: FortranRoutine, docstring: RoutineDocstring):
    # logger.debug("Checking routine arguments against docstring: %s", routine.name)
    #
    # have_all_args = all(arg.lower() in docstring.arguments_by_name for arg in routine.declared_argument_list)
    # if not have_all_args:
    #     # replace arguments with those from the actual subroutine; the interface file was wrong?
    #     renames, result_name = find_renames(docstring.filename, routine)
    #
    #     if len(routine.declared_argument_list) != len(renames):
    #         raise RenameError(f"{routine.declared_argument_list} -> {renames} length mismatch")
    #
    #     all_renames = list(zip(routine.declared_argument_list, renames, strict=False))
    #
    #     if result_name and routine.result_name and result_name != routine.result_name:
    #         all_renames.append((routine.result_name, result_name))
    #         logger.warning(f"Function return name changed: {routine.result_name} -> {result_name}")
    #
    #     for interface_arg_name, source_arg_name in all_renames:
    #         if interface_arg_name != source_arg_name:
    #             arg = routine.declarations.pop(interface_arg_name)
    #             routine.declarations[source_arg_name] = arg
    #             arg.name = source_arg_name
    #             logger.warning(
    #                 f"Routine: {routine.name} arg rename {interface_arg_name} -> {source_arg_name}"
    #             )
    #
    #     routine.declared_argument_list = list(renames)
    #
    # extra_args = {arg for arg in routine.declarations if arg.lower() not in routine.declared_argument_list}
    # for arg in extra_args:
    #     logger.warning(
    #         f"Routine: {routine.name} interface has extra arg defined {arg!r} ({routine.filename}:{routine.declarations[arg].line})"
    #     )
    #     routine.declarations.pop(arg)
    #
    doc_args = docstring.arguments_by_name
    for arg in routine.args:
        doc_arg = doc_args.get(arg.c_name.lower(), None)
        if doc_arg is None:
            # reparse failure?
            continue
        doc_arg.arg_type = arg.member.type_info
        if arg.is_optional:
            doc_arg.is_optional = True
        if not doc_arg.data_type:
            # this is not exactly accurate for array data
            doc_arg.data_type = arg.type


def generate_cpp_routine_code(template: str, routines: dict[str, FortranRoutine], settings: RoutineSettings):
    lines = []
    for routine in sorted_routines(routines):
        if not routine.usable:
            continue
        cpp_wrapper = generate_routine_cpp_wrapper(routine)
        lines.extend(cpp_wrapper)

    tpl = string.Template(template.replace("// ${", "${"))
    code = tpl.substitute(header_filename=settings.cpp_header_filename)
    return code + "\n".join(lines)


def generate_fortran_routine_code(
    template: str, routines: dict[str, FortranRoutine], settings: RoutineSettings
):
    lines = []
    using = {}

    tpl = string.Template(template.replace("! ${", "${"))
    for routine in sorted_routines(routines):
        if not routine.usable:
            continue

        using.setdefault(routine.module, set())
        using[routine.module].add(routine.overloaded_name)

        lines.append(generate_fortran_routine_with_c_binding(routine))

    use_lines = [
        wrap_line(f"use {module}, only: {', '.join(sorted(names))}") for module, names in using.items()
    ]
    return tpl.substitute(
        module_name=settings.fortran_module_name,
        use_lines="\n".join(use_lines),
        routines="\n".join(lines),
    )


SettingsAndRoutines = list[tuple[RoutineSettings, list[FortranRoutine]]]


def load_routines(
    params: CodegenConfig,
    parsed_structs: list[ParsedStructure],
    config: CodegenConfig,
) -> tuple[SettingsAndRoutines, list[str]]:
    structs_seen = set()
    settings_and_routines = []
    for settings in config.routines:
        logger.info("Parsing routines for namespace: %s", settings.cpp_namespace)
        routines = [routine for routine in parse_bmad_routines(settings, params) if not routine.private]
        settings_and_routines.append((settings, routines))
        for routine in routines:
            if routine.private:
                continue

            structs_seen |= {arg.kind.lower() for arg in routine.args if arg.type == "type"}

    parsed_structs_by_name = {struct.name.lower(): struct for struct in parsed_structs}

    structs_to_use = set()
    structs_to_skip = set()
    check = list(structs_seen)
    while check:
        typ = check.pop(0)
        if typ not in parsed_structs_by_name:
            logger.warning("Skipping %s as it's not in the bmad parsed structure list", typ)
            structs_to_skip.add(typ)
            continue

        parsed_st = parsed_structs_by_name[typ]
        if (
            config.should_skip(parsed_st.module)
            or typ.lower() in config.skips
            or typ.lower() in config.component_no_translate_list  # legacy - remove?
        ):
            structs_to_skip.add(typ)
            continue

        structs_to_use.add(typ)
        for member in parsed_st.members.values():
            if member.type.lower() not in {"type", "class"}:
                continue
            assert member.kind
            kind = member.kind.lower()
            if kind not in structs_to_skip and kind not in structs_to_use and kind not in check:
                check.append(kind)

    return settings_and_routines, sorted(structs_to_use)


def generate_routines(params: CodegenConfig, settings_and_routines: SettingsAndRoutines):
    all_routines = []
    all_routines_by_name = {}
    to_write = {}
    for settings, routines in settings_and_routines:
        all_routines.extend(routines)
        logger.info(f"Pruning routines ({settings.fortran_output_filename})")
        routines_by_name = prune_routines(routines, params)
        all_routines_by_name.update(routines_by_name)
        logger.info(f"Generating code: routines Fortran side ({settings.fortran_module_name}")

        routines_header = generate_routines_header(
            template=(CODEGEN_ROOT / "routines.tpl.hpp").read_text(),
            routines=routines_by_name,
            settings=settings,
        )
        logger.info(f"Generating code: routines C++ side ({settings.cpp_namespace})")
        cpp_routine_code = generate_cpp_routine_code(
            template=(CODEGEN_ROOT / "routines.tpl.cpp").read_text(),
            routines=routines_by_name,
            settings=settings,
        )
        fortran_routine_code = generate_fortran_routine_code(
            template=(CODEGEN_ROOT / "routines.tpl.f90").read_text(),
            routines=routines_by_name,
            settings=settings,
        )

        generated = CPPBMAD_ROOT / "src" / "generated"
        header_path = CPPBMAD_ROOT / "include" / "bmad" / "generated" / settings.cpp_header_filename
        to_write[header_path] = routines_header
        to_write[generated / settings.fortran_output_filename] = fortran_routine_code
        to_write[generated / settings.cpp_output_filename] = cpp_routine_code

    unique_routines = {rt.name: rt for rt in all_routines}
    usable = [rt for rt in unique_routines.values() if rt.usable]

    # enums_by_name = ctx().enums_by_name
    # for rt in usable:
    #     for arg in rt.args:
    #         arg.check_enums(enums_by_name)
    logger.info("Procedures: %d usable / %d total unique", len(usable), len(all_routines_by_name))
    return all_routines, all_routines_by_name, to_write
