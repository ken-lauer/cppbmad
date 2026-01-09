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
        data_type = type_info[0]

        try:
            typ = TypeInformation.from_line(f"{data_type} :: {self.name}")
            self.data_type = type_information_to_python_type(typ)
        except Exception:
            logger.warning(f"TODO: arg parsing {self.description}")

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
    inputs: list[DocstringParameter] = field(default_factory=list)
    outputs: list[DocstringParameter] = field(default_factory=list)
    result_variable: str | None = None
    is_overloaded: bool = False
    overloaded_versions: list[str] = field(default_factory=list)
    related_routines: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def arguments_by_name(self) -> dict[str, DocstringParameter]:
        return {arg.name.lower(): arg for arg in self.inputs + self.outputs}

    def update_parameter(self, param: DocstringParameter):
        name = param.name.lower()
        if name not in self.arguments_by_name:
            if param.is_input:
                self.inputs.append(param)
            else:
                self.outputs.append(param)
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

    def to_description(self) -> str:
        """String representation of the routine."""
        result = f"{self.routine_type.value}: {self.name}"
        if self.result_variable:
            result += f" result({self.result_variable})"
        result += "\n"

        if self.description:
            result += "\nDescription:\n"
            for line in self.description:
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

        if self.inputs:
            result += "\nInputs:\n"
            for input_param in self.inputs:
                optional_str = " (Optional)" if input_param.is_optional else ""
                type_str = f": {input_param.data_type}" if input_param.data_type else ""
                result += f"  {input_param.name}{type_str}{optional_str} -- {input_param.description}\n"

        if self.outputs:
            result += "\nOutputs:\n"
            for output_param in self.outputs:
                optional_str = " (Optional)" if output_param.is_optional else ""
                type_str = f": {output_param.data_type}" if output_param.data_type else ""
                result += f"  {output_param.name}{type_str}{optional_str} -- {output_param.description}\n"

        if self.result_variable:
            result += f"\nResult: {self.result_variable}\n"

        return result

    def to_numpy_docstring(self) -> str:
        """Create a NumPy-style docstring from the routine information."""
        lines = []

        if self.description:
            lines.append(self.description[0])
            lines.append("")
            if len(self.description) > 1:
                lines.extend(self.description[1:])
                lines.append("")

        if self.inputs:
            lines.extend(["Parameters", "----------"])
            for p in self.inputs:
                type_str = f"{p.data_type}" if p.data_type else ""
                if p.is_optional:
                    type_str += ", optional"
                lines.append(f"{p.name} : {type_str}")
                lines.extend(_wrap_docstring_lines(p.description.lstrip(), indent="    ").splitlines())

        has_returns = bool(self.result_variable or self.outputs)
        if has_returns:
            lines.extend(["", "Returns", "-------"])

            if self.result_variable and not self.outputs:
                lines.append(f"{self.result_variable}")

            for p in self.outputs:
                type_str = f"{p.data_type}" if p.data_type else ""
                lines.append(f"{p.name} : {type_str}")
                lines.extend(_wrap_docstring_lines(p.description, indent="    ").splitlines())

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
        return "\n".join(lines)


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
            type_name = "unknown"

    if not dt.dimension:
        return type_name

    return f"array of {type_name} ({dt.dimension})"


separator_regex = re.compile(r"^!\-+$")


def _split_comment_blocks(comments: str) -> list[tuple[int, str]]:
    """Split the comments into separate blocks based on separator lines or empty lines.
    Returns a list of tuples containing the starting line number and the block text."""

    lines = comments.split("\n")
    blocks = []
    current_block = []
    block_start_line = 0  # Track the starting line number of the current block

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if separator_regex.match(line):
            if current_block:
                blocks.append((block_start_line, "\n".join(current_block)))
                current_block = []

            while i < len(lines) and separator_regex.match(lines[i].strip()):
                i += 1

            # Update block_start_line for the next block
            block_start_line = i + 1
            continue

        if not line and current_block and i + 1 < len(lines):
            next_line = lines[i + 1].strip()
            if (
                not next_line
                or next_line.startswith("!+")
                or next_line.lower().startswith("function")
                or next_line.lower().startswith("subroutine")
            ):
                blocks.append((block_start_line, "\n".join(current_block)))
                current_block = []
                block_start_line = i + 2  # Update start line for next block

        if not current_block and len(current_block) == 0:
            # If this is the first line of a new block, update the start line
            block_start_line = i + 1

        current_block.append(lines[i])
        i += 1

    if current_block:
        blocks.append((block_start_line, "\n".join(current_block)))

    return blocks


def _extract_related_routines(line: str) -> list[str]:
    """Extract related routine names from a line."""

    if ":" in line:
        line = line.split(":", 1)[1]

    related = [name.strip() for name in re.split(r"[,\s]+", line) if name.strip()]
    return [name for name in related if name and name.lower() not in {"and", "or", "see", "also"}]


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

    docstring = RoutineDocstring(
        filename=filename,
        lineno=definition_lineno,
        name=name,
        routine_type=routine_type,
        result_variable=result_var,
    )

    current_section = "description"
    current_param = None
    in_note = False
    current_note = ""

    i = lines.index(definition_line) + 1

    while i < len(lines):
        line = lines[i].strip()
        next_line = lines[i + 1].strip() if i < len(lines) - 1 else ""

        if not line or re.match(r"^![-+]*$", line):
            i += 1
            continue

        if line.startswith("!"):
            line = line.removeprefix("!").strip()

        if line.lower() == "input:":
            current_section = "inputs"
            in_note = False
            i += 1
            continue
        if line.lower() == "output:":
            current_section = "outputs"
            in_note = False
            i += 1
            continue
        if "also see:" in line.lower():
            current_section = "related"
            in_note = False
            docstring.related_routines.extend(_extract_related_routines(line))
            i += 1
            continue

        if current_section == "description":
            if "is an overloaded name for" in line.lower():
                docstring.is_overloaded = True
                in_note = False
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

            elif re.match(r"^note(?:\s*\:|\s+that\b)", line.lower()):
                in_note = True
                current_note = line
                i += 1

            elif in_note:
                current_note += " " + line
                i += 1
            else:
                if re.match(r"^remember(?:\s*\:|\s+)", line.lower()):
                    docstring.notes.append(line)
                else:
                    docstring.description.append(line)
                i += 1

            if in_note and (i >= len(lines) or lines[i].strip() == "" or ":" in lines[i]):
                docstring.notes.append(current_note)
                in_note = False
        elif current_section in {"inputs", "outputs"}:
            param_match = re.match(r"^(\w+)(?:\s*\(.*?\))?\s*--\s*(.*?)$", line)
            param_name = None
            param_desc = ""
            if param_match:
                param_name = param_match.group(1)
                param_desc = param_match.group(2)
            elif "--" in next_line:
                if m := re.match(r"^!\s{3,4}(\S+)$", lines[i]):
                    param_name = m.groups()[0]
                    param_desc = next_line

            if param_name:
                data_type = None
                if "--" in param_desc:
                    parts = param_desc.split("--", 1)
                    data_type = parts[0].strip()
                    param_desc = parts[1].strip()

                is_optional = "optional" in param_desc.lower()

                param = DocstringParameter(
                    name=param_name,
                    description=param_desc,
                    data_type=data_type,
                    is_optional=is_optional,
                    is_input=current_section == "inputs",
                    is_output=current_section == "outputs",
                )

                if current_section == "inputs":
                    docstring.inputs.append(param)
                else:
                    existing_arg = docstring.arguments_by_name.get(param_name, None)
                    if existing_arg:
                        existing_arg.is_output = True
                        desc = param.description.strip()
                        if ":" in desc:
                            desc = desc.split(":")[1].strip()
                        if desc not in existing_arg.description:
                            existing_arg.description += f"\nThis parameter is an input/output and is modified in-place. As an output: {desc}"
                    else:
                        docstring.outputs.append(param)

                current_param = param
            elif current_param:
                if line not in current_param.description:
                    current_param.description += " " + line

            i += 1
        elif current_section == "related":
            docstring.related_routines.extend(_extract_related_routines(line))
            i += 1
        else:
            i += 1

    if in_note and current_note:
        docstring.notes.append(current_note)

    for arg in docstring.inputs + docstring.outputs:
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
