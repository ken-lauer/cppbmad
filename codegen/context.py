from __future__ import annotations

import contextvars
import logging
import pathlib
import re
import subprocess
import typing
from dataclasses import dataclass, field
from typing import Annotated

import pydantic
import tomllib

from . import paths
from .paths import CODEGEN_ROOT
from .proxy import struct_to_proxy_class_name

logger = logging.getLogger(__name__)

if typing.TYPE_CHECKING:
    from .arg import CodegenStructure
    from .enums import EnumValue
    from .routines import FortranRoutine
    from .structs import ParsedStructure


NormalizedPath = Annotated[pathlib.Path, pydantic.BeforeValidator(paths.normalize)]


class ProjectSettings(pydantic.BaseModel):
    name: str
    title: str
    description: str
    cpp_namespace: str = ""


class RoutineSettings(pydantic.BaseModel):
    fortran_output_filename: str
    cpp_output_filename: str
    cpp_namespace: str
    python_submodule: str = ""
    python_top_level: bool = True
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


class HookSpec(pydantic.BaseModel):
    """
    A single Bmad/Tao callback hook whose C ABI contract is generated.
    """

    # "bmad" or "tao" -- selects the C++ namespace
    project: str
    # full hook name, e.g. "track1_wake"
    name: str
    # C++ std::function alias; default snake_to_camel(name) + "Hook"
    type_name: str = ""
    # parsed abstract-interface name; default below
    def_name: str = ""
    # real/integer args the hook writes (cross by reference)
    inout_scalars: list[str] = []
    # args the C++ callback returns by value (default: none -> void, out-params by ref)
    returns: list[str] = []

    @property
    def cpp_namespace(self) -> str:
        return {"bmad": "Bmad", "tao": "Tao"}[self.project]

    @property
    def register_c_name(self) -> str:
        return f"{self.project}_hook_register_{self.name}"

    @property
    def def_routine_name(self) -> str:
        if self.def_name:
            return self.def_name
        return f"tao_hook_{self.name}_def" if self.project == "tao" else f"{self.name}_def"


@dataclass
class UpstreamInfo:
    """Git version info for the upstream Bmad source."""

    tag: str = ""
    commit_hash: str = ""

    @classmethod
    def from_source_dir(cls, source_dir: pathlib.Path) -> UpstreamInfo:
        """Read git tag and commit hash from the upstream source directory."""
        try:
            tag = subprocess.check_output(
                ["git", "describe", "--tags", "--always"],
                cwd=source_dir,
                text=True,
            ).strip()
            commit_hash = subprocess.check_output(
                ["git", "rev-parse", "HEAD"],
                cwd=source_dir,
                text=True,
            ).strip()
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("Could not read git info from %s", source_dir)
            return cls()
        return cls(tag=tag, commit_hash=commit_hash)


class CodegenConfig(pydantic.BaseModel):
    struct_list: list[str]
    component_no_translate_list: list[str]
    interface_ignore_list: list[str]
    structs_defined_externally: list[str]
    include_header_files: list[str]
    c_side_name_translation: dict[str, str]
    c_to_python_name_translation: dict[str, str]
    skips: list[str]
    source_dir: NormalizedPath = pathlib.Path()
    upstream_source_url: str = ""
    projects: list[ProjectSettings] = []
    routines: list[RoutineSettings] = []
    hooks: list[HookSpec] = []
    enum_filenames: list[NormalizedPath] = []
    python_imports: dict[str, list[str]] = {}
    python_module_name: str = "_pybmad"

    @classmethod
    def from_file(cls, filename: pathlib.Path) -> CodegenConfig:
        with filename.open("rb") as fp:
            return CodegenConfig.model_validate(tomllib.load(fp))

    def should_skip(self, name: str) -> bool:
        if name in self.skips:
            return True

        patterns = [skip for skip in self.skips if "." in skip or "*" in skip or "$" in skip]
        return any(re.match(pat, name) for pat in patterns)


SUPPORTED_ARRAY_DIMS = (1, 2, 3)


@dataclass
class ConfigContext:
    params: CodegenConfig
    upstream: UpstreamInfo = field(default_factory=UpstreamInfo)
    parsed_structs: list[ParsedStructure] = field(default_factory=list)
    codegen_structs: list[CodegenStructure] = field(default_factory=list)
    enums: dict[str, dict[str, EnumValue]] = field(default_factory=dict)
    missing_struct_attrs: dict[str, dict[str, str]] = field(default_factory=dict)
    routines: list[FortranRoutine] = field(default_factory=list)
    routines_by_name: dict[str, FortranRoutine] = field(default_factory=dict)
    pybmad_files: dict[pathlib.Path, str] = field(default_factory=dict)
    proxy_files: dict[pathlib.Path, str] = field(default_factory=dict)
    routine_files: dict[pathlib.Path, str] = field(default_factory=dict)
    docs_files: dict[pathlib.Path, str] = field(default_factory=dict)

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
        return {struct.f_name.lower(): struct for struct in self.codegen_structs}


config_context: contextvars.ContextVar[ConfigContext] = contextvars.ContextVar("config_context")


def get_params() -> CodegenConfig:
    return config_context.get().params


def ctx() -> ConfigContext:
    return config_context.get()
