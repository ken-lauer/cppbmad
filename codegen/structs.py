from __future__ import annotations

import contextlib
import importlib
import logging
import pathlib
import sys

from .paths import ACC_ROOT_DIR, CPPBMAD_ROOT


def find_bmad_struct_parser(
    bmad_root: pathlib.Path | None = None,
    subdirectory: pathlib.Path = pathlib.Path("structs"),
):
    with contextlib.suppress(ImportError):
        importlib.import_module("bmad_struct_parser")

    if bmad_root is None:
        if ACC_ROOT_DIR is not None:
            bmad_root = ACC_ROOT_DIR
        else:
            raise ImportError(
                "Unable to import bmad_struct_parser: ACC_ROOT_DIR is unset and no Bmad path was given"
            )

    structs_path = str(bmad_root / subdirectory)
    logging.info(f"Adding {structs_path} to sys.path...")
    sys.path.insert(0, structs_path)


find_bmad_struct_parser()


import bmad_struct_parser

json_path = pathlib.Path(bmad_struct_parser.__file__).resolve().absolute().parents[1] / "json"


from bmad_struct_parser import Structure as ParsedStructure
from bmad_struct_parser import TypeInformation
from bmad_struct_parser import (
    load_configured_structures as _load_bmad_parser_structures,
)
from bmad_struct_parser.config import DEFAULT_CONFIG_FILE, SourceConfig
from bmad_struct_parser.parser import (
    ParserConfig,
    StructureMember,
    find_structs_in_file,
    get_in_parenthesis,
    parse_declaration,
)
from bmad_struct_parser.util import FileLine, join_ampersand_lines, split_comment


def load_bmad_parser_structures(
    config_file: pathlib.Path = DEFAULT_CONFIG_FILE, output_subpath: str = "structs"
):
    structs = _load_bmad_parser_structures(config_file, output_subpath=output_subpath)
    parser_config = ParserConfig.from_file(config_file)
    test_f90 = CPPBMAD_ROOT / "src" / "test_struct_defs.f90"
    source_config = SourceConfig(
        source_dir=test_f90.parent,
        fortran_filename=test_f90,
        json_filename="tests.json",
        function_prefix="",
        precision_module="",
        skip_includes=(),
        # json_config: JsonConfig = dataclasses.field(default_factory=JsonConfig)
        # skip_structs: tuple[str, ...] = dataclasses.field(default_factory=tuple)
        # include_dirs: tuple[pathlib.Path, ...] = dataclasses.field(default_factory=tuple)
    )
    for struct in find_structs_in_file(parser_config, source_config, test_f90):
        struct.parse()
        structs.append(struct)
    return structs


__all__ = [
    "FileLine",
    "ParsedStructure",
    "StructureMember",
    "TypeInformation",
    "get_in_parenthesis",
    "join_ampersand_lines",
    "load_bmad_parser_structures",
    "parse_declaration",
    "split_comment",
]
