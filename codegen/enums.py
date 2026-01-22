#!/usr/bin/env python3

# Note:
#   Run this script in the cpp_bmad_interface directory.
#   This script searches a set of Bmad files and generates corresponding constants for use with C++ code.
#   The constants file is: bmad/generated/enums.hpp
#   For example, the proton$ parameter on the Fortran side is translated to PROTON on the C++ side.
#
from __future__ import annotations

import pathlib
import re
from dataclasses import dataclass

from .paths import CPPBMAD_INCLUDE
from .structs import FileLine

ENUM_FILENAME = CPPBMAD_INCLUDE / "bmad" / "generated" / "enums.hpp"

re_int = re.compile("INTEGER, *PARAMETER *:: *")
re_real = re.compile(r"REAL\(RP\), *PARAMETER *:: *")
re_d_exp = re.compile(r"\dD[+-]?\d")
re_equal = re.compile(r"\=.*\dD[+-]?\d")


@dataclass
class EnumValue:
    name: str
    value: str
    type: str
    comment: str


def _clean_line(line: str) -> str:
    # Strip datatype information
    line = re_int.sub("", line)
    line = re_real.sub("", line)

    if "Z'" in line:
        line = line.replace("Z'", "0x").replace("'", "")
    if "z'" in line:
        line = line.replace("z'", "0x").replace("'", "")

    # Handle Fortran-style exponents, replace `D` with `E`
    if re_equal.search(line):
        sub = re_d_exp.search(line)
        if sub:
            d_exp_value = sub.group(0).replace("D", "E")  # Replace "3D6" with "3E6"
            line = re_d_exp.sub(d_exp_value, line)
    return line


name_fixes = {
    "TRUE": "TRUE_",
    "FALSE": "FALSE_",
}


def parse_fortran_enums(fn: pathlib.Path) -> dict[str, EnumValue]:
    """
    Parse Fortran-like constants into a dictionary of {enum_name: enum_value}.
    """
    enum_dict = {}

    for file_line in FileLine.from_file(fn):
        line, comment = file_line.split_comment()
        line = line.upper()

        # Skip parameter arrays or irrelevant lines
        if "[" in line:
            continue

        if re_int.match(line):
            type_ = "int"
        elif re_real.match(line):
            type_ = "double"
        else:
            continue

        line = _clean_line(line)
        # Skip lines with arrays
        if "(" in line:
            continue

        for idx, part in enumerate(line.split(",")):
            part = part.strip()

            name, value = part.split("=", 1)
            name = name.strip().rstrip("$")
            value = value.strip().replace("$", "").replace("_RP", "")
            enum_dict[name] = EnumValue(
                name=name_fixes.get(name, name),
                value=value,
                type=type_,
                comment=comment if idx == 0 else "",
            )

    return enum_dict


def get_enums_in_range(enums: dict[str, EnumValue], start_key: str, num_key: str) -> list[EnumValue]:
    # Length will always be the first attribute.
    # We go from L and search until we exceed NUM_ELE_ATTRIB to find attributes.
    first_attr = enums[start_key]
    num_enums = int(enums[num_key].value)
    attrs = list(enums.values())
    attrs = attrs[attrs.index(first_attr) :]
    result = []
    for attr in attrs:
        try:
            cur_value = int(attr.value)
            if cur_value > num_enums:
                break
        except TypeError:
            break

        if result:
            last_value = int(result[-1].value)
            if last_value == num_enums and cur_value < num_enums:
                # Wrapping around to 1 for the next set of enums
                break
        result.append(attr)

    return result


def get_ele_attributes(enums: dict[str, EnumValue]):
    return get_enums_in_range(enums, start_key="L", num_key="NUM_ELE_ATTRIB")


def get_ele_keys(enums: dict[str, EnumValue]):
    return get_enums_in_range(enums, start_key="DRIFT", num_key="N_KEY")


def get_class_code(clsname: str, enums: list[EnumValue], offset: int = 0) -> str:
    code = [f"enum class {clsname} : size_t {{"]
    if offset > 0:
        offset_code = f" + {offset}" if offset else ""
    else:
        offset_code = f" - {abs(offset)}" if offset else ""
    for attr in enums:
        if attr.comment:
            code.append(f"// {attr.comment}")
        code.append(f"  {attr.name} = {attr.value}{offset_code},")
    code.append(f"}}; // enum class {clsname}")
    return "\n".join(code)


def parse_all_enums(enum_filenames: list[pathlib.Path]):
    return {fn.name: parse_fortran_enums(fn) for fn in enum_filenames}


def is_integer(s: str) -> bool:
    if not s:
        return False
    if s[0] == "-":
        return s[1:].isdigit()
    return s.isdigit()


def get_enums_from_bounds(bounds: list[str]) -> set[str]:
    enums = set()
    for bound in bounds:
        bound = bound.strip()

        if not bound or is_integer(bound):
            continue

        if bound.endswith("$"):
            bound = bound.removesuffix("$")

        enums.add(bound.lower())

    return enums


def to_cpp_enum(enum_str: str, namespace: str = "Bmad") -> str:
    clean_str = enum_str.strip()
    if clean_str.endswith("$"):
        clean_str = clean_str[:-1]

    return f"{namespace}::{clean_str.upper()}"


def replace_enums_with_cpp(bounds: list[str], enums: dict[str, EnumValue]) -> list[str]:
    result = []
    for b in bounds:
        b = b.strip()

        if is_integer(b):
            result.append(b)
            continue
        key = b.removesuffix("$").upper()

        if "(" in key:
            raise ValueError("Calls in array bounds are not supported")

        if not key:
            result.append(key)
        elif key in enums:
            result.append(to_cpp_enum(key))
        else:
            raise KeyError(f"Enum '{key}' found in bounds '{b}' but not in provided map.")

    return result


def get_enum_code(enums_by_filename: dict[str, dict[str, EnumValue]]):
    result = [
        """
//+
// C++ constants equivalent to Bmad parameters.
//
// This file is generated as part of cppbmad.
// The code generation files can be found in cppbmad/codegen.
//
// DO NOT EDIT THIS FILE DIRECTLY! 
//-

#pragma once

#include <cstddef>

namespace Bmad {
"""
    ]

    for fn, enums in enums_by_filename.items():
        result.append("")
        result.append(f"// Enums from {fn}")
        for enum in enums.values():
            if enum.comment:
                result.append(f"// {enum.comment}")

            result.append(f"const {enum.type} {enum.name} = {enum.value};")

        if fn == "bmad_struct.f90":
            result.append(
                get_class_code(
                    "EleAttribute",
                    get_ele_attributes(enums),
                    offset=-1,
                )
            )
            result.append(
                get_class_code(
                    "EleAttributeFortran",
                    get_ele_attributes(enums),
                )
            )
            result.append(get_class_code("EleKey", get_ele_keys(enums)))

    result.append("""
}

""")
    return "\n".join(result)
