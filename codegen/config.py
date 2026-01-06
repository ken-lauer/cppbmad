from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class CodegenConfig:
    struct_list: list[str]
    component_no_translate_list: list[str]
    interface_ignore_list: list[str]
    structs_defined_externally: list[str]
    include_header_files: list[str]
    c_side_name_translation: dict[str, str]
    c_to_python_name_translation: dict[str, str]
    skips: list[str]

    equality_use_statements: list[str] = field(default_factory=list)  # TODO remove
    test_use_statements: list[str] = field(default_factory=list)  # TODO remove
