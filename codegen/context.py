from __future__ import annotations

import contextvars
import pathlib
import re
import typing
from dataclasses import dataclass, field
from typing import Annotated

import pydantic
import tomllib

from . import paths
from .paths import CODEGEN_ROOT
from .proxy import struct_to_proxy_class_name

if typing.TYPE_CHECKING:
    from .arg import CodegenStructure
    from .enums import EnumValue
    from .routines import FortranRoutine
    from .structs import ParsedStructure


NormalizedPath = Annotated[pathlib.Path, pydantic.BeforeValidator(paths.normalize)]


class RoutineSettings(pydantic.BaseModel):
    fortran_output_filename: str
    cpp_output_filename: str
    cpp_namespace: str
    # interface_path: pathlib.Path | None
    source_paths: list[NormalizedPath]
    skip_files: set[str]
    skip_procedures: set[str]
    do_not_overload: set[str] = set()

    @property
    def fortran_module_name(self):
        stem = pathlib.Path(self.fortran_output_filename).stem
        return f"cppbmad_{stem}"

    @property
    def cpp_header_filename(self):
        return pathlib.Path(self.cpp_output_filename).with_suffix(".hpp").name


class CodegenConfig(pydantic.BaseModel):
    struct_list: list[str]
    component_no_translate_list: list[str]
    interface_ignore_list: list[str]
    structs_defined_externally: list[str]
    include_header_files: list[str]
    c_side_name_translation: dict[str, str]
    c_to_python_name_translation: dict[str, str]
    skips: list[str]
    routines: list[RoutineSettings] = []
    enum_filenames: list[NormalizedPath] = []

    @classmethod
    def from_file(cls, filename: pathlib.Path) -> CodegenConfig:
        with filename.open("rb") as fp:
            return CodegenConfig.model_validate(tomllib.load(fp))

    def should_skip_routine(self, name: str) -> bool:
        if name in self.skips:
            return True

        patterns = [skip for skip in self.skips if "." in skip or "*" in skip]
        return any(re.match(pat, name) for pat in patterns)


SUPPORTED_ARRAY_DIMS = (1, 2, 3)


@dataclass
class ConfigContext:
    params: CodegenConfig
    parsed_structs: list[ParsedStructure] = field(default_factory=list)
    codegen_structs: list[CodegenStructure] = field(default_factory=list)
    enums: dict[str, dict[str, EnumValue]] = field(default_factory=dict)
    missing_struct_attrs: dict[str, dict[str, str]] = field(default_factory=dict)
    routines: list[FortranRoutine] = field(default_factory=list)
    routines_by_name: dict[str, FortranRoutine] = field(default_factory=dict)
    pybmad_files: dict[pathlib.Path, str] = field(default_factory=dict)
    routine_files: dict[pathlib.Path, str] = field(default_factory=dict)

    @property
    def report_html(self) -> str:
        from .coverage import generate_coverage_report

        return generate_coverage_report(self.routines, self.codegen_structs, self.missing_struct_attrs)

    @property
    def cpp_to_string_header(self) -> str:
        from .cpp import generate_to_string_header

        return generate_to_string_header(
            template=(CODEGEN_ROOT / "to_string.tpl.hpp").read_text(),
            structs=self.codegen_structs,
            routines=self.routines_by_name,
        )

    @property
    def cpp_to_string_code(self) -> str:
        from .cpp import generate_to_string_code

        return generate_to_string_code(
            template=(CODEGEN_ROOT / "to_string.tpl.cpp").read_text(),
            structs=self.codegen_structs,
            routines=self.routines_by_name,
        )

    @property
    def enums_by_name(self) -> dict[str, EnumValue]:
        all_enums = {}
        for enums in self.enums.values():
            for name, enum in enums.items():
                all_enums[name.upper()] = enum
        return all_enums

    @property
    def parser_structs_by_cpp_class(self) -> dict[str, ParsedStructure]:
        assert len(self.parsed_structs)
        return {struct_to_proxy_class_name(struct.name): struct for struct in self.parsed_structs}

    @property
    def codegen_structs_by_cpp_class(self) -> dict[str, CodegenStructure]:
        assert len(self.codegen_structs)
        return {struct.cpp_class: struct for struct in self.codegen_structs}

    @property
    def codegen_structs_by_name(self) -> dict[str, CodegenStructure]:
        assert len(self.codegen_structs)
        return {struct.f_name: struct for struct in self.codegen_structs}


config_context: contextvars.ContextVar[ConfigContext] = contextvars.ContextVar("config_context")


def get_params() -> CodegenConfig:
    return config_context.get().params


def ctx() -> ConfigContext:
    return config_context.get()
