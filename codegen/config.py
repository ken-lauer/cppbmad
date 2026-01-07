from __future__ import annotations

import pathlib
from typing import Annotated

import pydantic
import tomllib

from . import paths

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

    @classmethod
    def from_file(cls, filename: pathlib.Path) -> CodegenConfig:
        with filename.open("rb") as fp:
            return CodegenConfig.model_validate(tomllib.load(fp))


SUPPORTED_ARRAY_DIMS = (1, 2, 3)
