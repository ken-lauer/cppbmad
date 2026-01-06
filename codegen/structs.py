from __future__ import annotations

import pathlib

from .paths import CPPBMAD_ROOT, DEFAULT_STRUCT_CONFIG_FILE
from .struct_parser import Structure as ParsedStructure
from .struct_parser import TypeInformation
from .struct_parser.config import SourceConfig
from .struct_parser.parser import (
    ParserConfig,
    StructureMember,
    find_structs_in_file,
    get_in_parenthesis,
    parse_all_structures,
    parse_declaration,
)
from .struct_parser.util import FileLine, join_ampersand_lines, split_comment


def load_bmad_parser_structures(config_file: pathlib.Path = DEFAULT_STRUCT_CONFIG_FILE):
    parser_conf = ParserConfig.from_file(config_file)
    # structs = _load_bmad_parser_structures(config_file, output_subpath=output_subpath)
    structs = parse_all_structures(parser_conf)

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
