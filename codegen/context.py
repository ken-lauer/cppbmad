from __future__ import annotations

import contextvars
import dataclasses
import typing

from .proxy import struct_to_proxy_class_name

if typing.TYPE_CHECKING:
    from .arg import CodegenStructure
    from .config import CodegenConfig
    from .enums import EnumValue
    from .structs import ParsedStructure


@dataclasses.dataclass
class ConfigContext:
    params: CodegenConfig
    parsed_structs: list[ParsedStructure] = dataclasses.field(default_factory=list)
    codegen_structs: list[CodegenStructure] = dataclasses.field(default_factory=list)
    enums: dict[str, dict[str, EnumValue]] = dataclasses.field(default_factory=dict)

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
