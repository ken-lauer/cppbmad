from __future__ import annotations

import logging
import pathlib
import re
import textwrap
from dataclasses import dataclass, field

import pydantic
import pydantic.alias_generators

from .structs import (
    TypeInformation,
    get_in_parenthesis,
)
from .types import Intent, RoutineType

logger = logging.getLogger(__name__)


base_fortran_to_python_type = {
    "real": "float",
    "real16": "float",
    "complex": "complex",
    "integer": "int",
    "integer8": "int",
    "logical": "bool",
    "character(*)": "str",
}


@dataclass
class DocstringParameter:
    """Class for representing a parameter/argument in a Fortran routine."""

    name: str = ""
    description: str = ""
    data_type: str | None = None
    is_optional: bool = False
    is_input: bool = False
    is_output: bool = False
    guessed: bool = False
    arg_type: TypeInformation | None = None

    @property
    def is_input_only(self):
        return self.is_input and not self.is_output

    @property
    def is_output_only(self):
        return not self.is_input and self.is_output

    @property
    def is_input_output(self):
        return self.is_input and self.is_output

    @property
    def intent(self) -> Intent:
        if self.is_input and self.is_output:
            return "inout"
        if self.is_input:
            return "in"
        if self.is_output:
            return "out"
        raise NotImplementedError(f"No intent ({self})")

    def fix(self) -> None:
        if not self.description or ":" not in self.description:
            return

        type_info, desc = self.description.split(":", 1)
        if type_info.endswith("("):
            # e.g., real64(:,:,:), optional:
            in_paren = get_in_parenthesis(self.description)
            start = self.description.index(in_paren) + len(in_paren) + 1
            end = self.description.index(":", start)
            type_info = self.description[:end]
            self.description = self.description[len(type_info) :]
        else:
            self.description = desc.strip()
            type_info = [part.strip().lower() for part in type_info.split(",")]
            self.is_optional = "optional" in type_info

        self.description = self.description.replace("%", ".")
        self.data_type = type_info[0].strip()

        # try:
        #     typ = TypeInformation.from_line(f"{data_type} :: {self.name}")
        #     self.data_type = type_information_to_python_type(typ)
        # except Exception:
        #     logger.warning(f"TODO: arg parsing {self.description}")

        # if "(" in self.name or "%" in self.name:
        #     self.description = f"({self.name}) {self.description}"
        #     self.name = self.name.split("(")[0]
        #     self.name = self.name.split("%")[0]


def _wrap_docstring_lines(text: str, width: int = 110, indent: str = ""):
    """Wrap numpy docstring lines to specified width, preserving indentation."""

    wrapped = []
    for line in text.splitlines():
        # Determine the line's indentation
        leading_space = len(line) - len(line.lstrip())
        if leading_space > 0:
            line_indent = line[:leading_space]
            content = line[leading_space:]
        else:
            line_indent = indent
            content = line

        # If line is just whitespace, preserve it
        if not content:
            wrapped.append(line.strip())
            continue

        # Wrap the content, preserving the indentation
        wrapped_content = textwrap.fill(
            content,
            width=width - len(line_indent),
            initial_indent="",
            subsequent_indent="",
        )

        # Add the indentation back to each line
        for wl in wrapped_content.splitlines():
            wrapped.append(f"{line_indent}{wl}")

    return "\n".join(wrapped)


@dataclass
class RoutineDocstring:
    """Class for representing a Fortran routine (function or subroutine)."""

    name: str
    filename: pathlib.Path
    lineno: int
    routine_type: RoutineType
    description: list[str] = field(default_factory=list)
    params: list[DocstringParameter] = field(default_factory=list)
    result_variable: str | None = None
    is_overloaded: bool = False
    overloaded_versions: list[str] = field(default_factory=list)
    related_routines: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def arguments_by_name(self) -> dict[str, DocstringParameter]:
        return {arg.name.lower(): arg for arg in self.params}

    def update_parameter(self, param: DocstringParameter):
        name = param.name.lower()
        if name not in self.arguments_by_name:
            self.params.append(param)
            return None

        if param.guessed:
            return None

        old = self.arguments_by_name[name.lower()]
        old.name = param.name
        old.description = param.description
        old.data_type = param.data_type
        old.is_optional = param.is_optional
        old.is_input = param.is_input
        old.is_output = param.is_output
        old.guessed = param.guessed
        return old

    @property
    def inputs(self) -> list[DocstringParameter]:
        return [param for param in self.params if param.is_input]

    @property
    def outputs(self) -> list[DocstringParameter]:
        return [param for param in self.params if param.is_output]

    def to_description(self) -> str:
        """Bmad-style string representation of the routine."""
        result = f"{self.routine_type.value}: {self.name}"
        if self.result_variable:
            result += f" result({self.result_variable})"
        result += "\n"

        if self.description:
            desc = textwrap.dedent("\n".join(self.description))
            for line in desc.splitlines():
                result += f"  {line}\n"

        if self.notes:
            result += "\nNotes:\n"
            for note in self.notes:
                result += f"  {note}\n"

        if self.is_overloaded:
            result += "\nOverloaded versions:\n"
            for version in self.overloaded_versions:
                result += f"  {version}\n"

        if self.related_routines:
            result += "\nRelated routines:\n"
            for routine in self.related_routines:
                result += f"  {routine}\n"

        def add_param(param: DocstringParameter):
            nonlocal result

            optional_str = " (Optional)" if param.is_optional else ""
            type_str = f": {param.data_type}" if param.data_type else ""
            result += f"  {param.name}{type_str}{optional_str} -- {param.description}\n"

        if self.inputs:
            result += "\nInputs:\n"
            for input_param in self.inputs:
                add_param(input_param)

        if self.outputs:
            result += "\nOutputs:\n"
            for param in self.outputs:
                add_param(param)

        if self.result_variable:
            result += f"\nResult: {self.result_variable}\n"

        return result

    def to_numpy_docstring(self) -> str:
        """Create a NumPy-style docstring from the routine information."""
        lines = []

        if self.description:
            desc = textwrap.dedent("\n".join(self.description))
            lines.extend(desc.splitlines())
        else:
            lines.append(f"Wrapper for Fortran routine {self.name}")

        lines.append("")

        def add_param(param: DocstringParameter, is_last: bool):
            if param.name.startswith("%"):
                return

            if param.arg_type:
                # we have parsed code argument info, use that
                data_type = type_information_to_python_type(param.arg_type)
            else:
                # otherwise, use the data type from the docstring
                data_type = param.data_type

            if not data_type and not param.description.strip():
                logger.warning("Unknown parameter in docstring? %r: %s", self.name, param)
                return

            type_str = str(data_type)
            if param.is_optional:
                type_str += ", optional"
            lines.append(f"{param.name} : {type_str}")
            lines.extend(_wrap_docstring_lines(param.description.lstrip(), indent="    ").splitlines())
            if not is_last and (type_str or param.description.lstrip()):
                lines.append("")

        if self.inputs:
            lines.extend(["Parameters", "----------"])
            for i, param in enumerate(self.inputs):
                add_param(param, is_last=(i == len(self.inputs) - 1))

        has_returns = bool(self.result_variable or self.outputs)
        if has_returns:
            lines.extend(["", "Returns", "-------"])

            if self.result_variable and not self.outputs:
                lines.append(f"{self.result_variable}")

            for i, param in enumerate(self.outputs):
                add_param(param, is_last=(i == len(self.outputs) - 1))

        if self.notes or self.related_routines or self.is_overloaded:
            lines.extend(["", "Notes", "-----"])
            for note in self.notes:
                lines.append(_wrap_docstring_lines(note))

            if self.related_routines:
                lines.append("Related routines:")
                lines.append(_wrap_docstring_lines(" ".join(self.related_routines)))

            if self.is_overloaded:
                lines.append(
                    _wrap_docstring_lines("Overloaded versions: " + ", ".join(self.overloaded_versions))
                )

        lines.append("")

        def remove_multiple_blanks(lines: list[str]):
            last = None
            for line in lines:
                if last == line == "":
                    pass
                else:
                    yield line.rstrip()
                last = line

        return "\n".join(remove_multiple_blanks(lines))


def type_information_to_python_type(dt: TypeInformation) -> str:
    if dt.type.lower() == "type":
        assert dt.kind is not None
        type_name = pydantic.alias_generators.to_pascal(dt.kind)
    elif dt.type.lower().endswith("_struct"):
        type_name = pydantic.alias_generators.to_pascal(dt.type.split()[0])
    else:
        try:
            type_name = base_fortran_to_python_type[dt.type.lower()]
        except Exception:
            type_name = dt.type.lower()

    if not dt.dimension:
        return type_name

    dim = dt.dimension.count(",") + 1
    parts = [part.strip() for part in dt.dimension.split(",") if part.strip() not in ("0:", ":")]
    detailed = ""
    if parts:
        detailed = f" (shape: {dt.dimension})"
    return f"{dim}D array of {type_name}{detailed}"


separator_regex = re.compile(r"^!\-+$")


def _extract_related_routines(line: str) -> list[str]:
    """Extract related routine names from a line."""

    if ":" in line:
        line = line.split(":", 1)[1]

    related = [name.strip() for name in re.split(r"[,\s]+", line) if name.strip()]
    return [name for name in related if name and name.lower() not in {"and", "or", "see", "also"}]


def split_param_names(name: str) -> list[str]:
    """
    Split a string of parameter names, optionally containing array dimensions,
    into a list of pure parameter names.
    """
    if not name:
        return []

    chunks: list[str] = []
    current_chunk: list[str] = []
    paren_depth = 0

    for char in name:
        if char == "(":
            paren_depth += 1
            current_chunk.append(char)
        elif char == ")":
            paren_depth -= 1
            current_chunk.append(char)
        elif char == "," and paren_depth == 0:
            chunks.append("".join(current_chunk).strip())
            current_chunk = []
        else:
            current_chunk.append(char)

    if current_chunk:
        chunks.append("".join(current_chunk).strip())

    cleaned_names: list[str] = []

    for chunk in chunks:
        base_name = chunk.split("(", 1)[0].strip()
        if base_name:
            cleaned_names.append(base_name)

    return cleaned_names


def parse_routine_comment_block(
    filename: pathlib.Path, block: str, starting_lineno: int
) -> RoutineDocstring | None:
    """Parse a comment block that documents a routine."""
    lines = [line.strip() for line in block.splitlines()]

    definition_line = None
    definition_lineno = starting_lineno
    routine_type = RoutineType.UNKNOWN
    name = ""
    result_var = None

    for lineno, line in enumerate(lines, start=starting_lineno):
        if "function" in line.lower() and "result" in line.lower():
            match = re.match(
                r"!\+?\s*Function\s+(\w+)\s*\((.*?)\)\s*result\s*\(\s*(.*?)\s*\)",
                line,
                re.IGNORECASE,
            )
            if match:
                routine_type = RoutineType.FUNCTION
                name = match.group(1)
                result_var = match.group(3)
                definition_line = line
                definition_lineno = lineno
                break
        elif "function" in line.lower():
            match = re.match(r"!\+?\s*Function\s+(\w+)\s*\((.*?)\)", line, re.IGNORECASE)
            if match:
                routine_type = RoutineType.FUNCTION
                name = match.group(1)
                result_var = None
                definition_line = line
                definition_lineno = lineno
                break

            match = re.match(r"!\+?\s*Function\s+(\w+)\s*", line, re.IGNORECASE)
            if match:
                routine_type = RoutineType.FUNCTION
                name = match.group(1)
                result_var = None
                definition_line = line
                definition_lineno = lineno
                break
        elif "subroutine" in line.lower() and "(" in line:
            match = re.match(r"!\+?\s*Subroutine\s+(\w+)\s*\((.*?)\)", line, re.IGNORECASE)
            if match:
                routine_type = RoutineType.SUBROUTINE
                name = match.group(1)
                result_var = None
                definition_line = line
                definition_lineno = lineno
                break

            match = re.match(r"!\+?\s*Subroutine\s+(\w+)\s*", line, re.IGNORECASE)
            if match:
                routine_type = RoutineType.SUBROUTINE
                name = match.group(1)
                result_var = None
                definition_line = line
                definition_lineno = lineno
                break

    if definition_line is None:
        return None

    def clean_line(line: str) -> str:
        line = line.strip()
        m = re.match(r"^![-+]*(.*)$", line)
        if m:
            return m.groups()[0]
        return ""

    description = []
    for line in lines:
        line = clean_line(line)
        if line.lower().strip() in (
            "input:",
            "output:",
            "also see:",
        ):
            break
        if line or description:
            description.append(line)

    docstring = RoutineDocstring(
        filename=filename,
        lineno=definition_lineno,
        name=name,
        routine_type=routine_type,
        result_variable=result_var,
        description=description or ["No docstring available."],
    )

    current_section = "description"
    current_params = []

    param_name = ""
    param_desc = ""
    i = lines.index(definition_line) + 1

    while i < len(lines):
        line = lines[i].strip()
        next_line = lines[i + 1].strip() if i < len(lines) - 1 else ""

        if not line or re.match(r"^![-+]*$", line):
            i += 1
            continue

        if line.startswith("!"):
            line = line.removeprefix("!").strip()

        if line.lower().rstrip(": ") == "input":
            current_section = "inputs"
            i += 1
            continue
        if line.lower().rstrip(": ") == "output":
            current_section = "outputs"
            i += 1
            continue
        if "also see:" in line.lower():
            current_section = "related"
            docstring.related_routines.extend(_extract_related_routines(line))
            i += 1
            continue

        if current_section == "description":
            if "is an overloaded name for" in line.lower():
                docstring.is_overloaded = True
                i += 1

                while i < len(lines):
                    next_line = lines[i].strip()
                    if next_line.startswith("!"):
                        next_line = next_line[1:].strip()

                    if re.match(r"^\s*Function\s+\w+", next_line):
                        docstring.overloaded_versions.append(next_line)
                        i += 1
                    else:
                        break

            else:
                if re.match(r"^remember(?:\s*\:|\s+)", line.lower()):
                    docstring.notes.append(line)
                i += 1

        elif current_section in {"inputs", "outputs"}:
            param_name = ""
            param_desc = ""
            if "--" in line:
                param_name, param_desc = [part.strip() for part in line.split("--", 1)]
            elif "--" in next_line:
                # long lines can be in the form:
                #     really_long_name
                #           -- type info
                next_line_has_name = bool(next_line.lstrip("! ").split("--", 1)[0].strip())
                if not next_line_has_name:
                    param_name = line.strip()

            if param_name:
                data_type = None
                is_optional = "optional" in param_desc.lower()

                current_params = []
                for name in split_param_names(param_name):
                    param = DocstringParameter(
                        name=name,
                        description=param_desc,
                        data_type=data_type,
                        is_optional=is_optional,
                        is_input=current_section == "inputs",
                        is_output=current_section == "outputs",
                    )

                    existing_arg = docstring.arguments_by_name.get(param.name, None)
                    if existing_arg:
                        if current_section == "inputs":
                            existing_arg.is_input = True
                        elif current_section == "outputs":
                            existing_arg.is_output = True
                        desc = param.description.strip()
                        if ":" in desc:
                            desc = desc.split(":")[1].strip()
                        if desc not in existing_arg.description:
                            if existing_arg.is_input and existing_arg.is_output:
                                existing_arg.description += (
                                    "\nThis parameter is an input/output and is modified in-place."
                                )
                            existing_arg.description += f"\nAs an output, {existing_arg.name}: {desc}"
                    else:
                        docstring.params.append(param)

                    current_params.append(param)

            elif current_params:
                for current_param in current_params:
                    if line not in current_param.description:
                        current_param.description += " " + line

            i += 1
        elif current_section == "related":
            docstring.related_routines.extend(_extract_related_routines(line))
            i += 1
        else:
            i += 1

    for arg in docstring.params:
        arg.fix()

    mistaken_optional_outputs = [arg for arg in docstring.outputs if arg.is_optional]
    if mistaken_optional_outputs:
        for arg in mistaken_optional_outputs:
            # In cppbmad/pybmad, outputs are always there no matter what;
            # Is there actual meaning in the fortran version for this?
            arg.is_optional = False

        logger.debug(
            f"Optional output annotations for: {docstring.name} ({docstring.filename}:{docstring.lineno})"
        )
    return docstring
