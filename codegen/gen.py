"""
Script to create:

* C++ classes for Fortran structures
* Translator between Fortran structure and C++ class
* Routines to check for equality between instances of a given Fortran structure.
* Routines to check for equality between instances of a given C++ class
* Program to check the Fortran / C++ translator

Note: The corresponding C++ class component for a pointer or allocatable Fortran
scalar struct component is an array whose length is zero if the Fortran component
is nullified and whose length is 1 otherwise.
"""

from __future__ import annotations

import argparse
import logging
import pathlib
import typing

from .arg import Argument, CodegenStructure
from .config import CodegenConfig
from .context import ConfigContext, config_context
from .coverage import generate_coverage_report
from .cpp import generate_to_string_code, generate_to_string_header
from .enums import ENUM_FILENAME, write_enums
from .paths import CODEGEN_ROOT, CPPBMAD_INCLUDE, CPPBMAD_ROOT, CPPBMAD_SRC, REPO_ROOT
from .proxy import create_cpp_proxy_header, create_cpp_proxy_impl, create_fortran_proxy_code
from .py import generate_pybmad
from .routines import generate_routines, parse_bmad_routines
from .structs import ParsedStructure, load_bmad_parser_structures
from .util import write_contents_if_differs, write_if_differs

if typing.TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)

DEBUG = False  # Change to True to enable more verbose printout


def match_structure_definition(
    parsed_structures: list[ParsedStructure],
    struct: CodegenStructure,
    params: CodegenConfig,
):
    for fstruct in parsed_structures:
        if struct.f_name.lower() == fstruct.name.lower():
            break
    else:
        raise RuntimeError(f"Structure not found: {struct.f_name!r}")

    struct.f_name = fstruct.name
    struct.arg = [
        Argument.from_fstruct(fstruct, member, params=params) for member in fstruct.members.values()
    ]
    struct.module = fstruct.module
    struct.parsed = fstruct


def write_parsed_structures(structs, fn):
    """
    Write parsed structure definitions to a file.
    """
    with pathlib.Path(fn).open("w") as f_out:
        for struct in structs:
            f_out.write("******************************************\n")
            f_out.write(f"{struct.f_name}    {len(struct.arg)}\n")
            for arg in struct.arg:
                f_out.write(f"    {arg.original_repr()}\n")


def check_missing(structs: list[CodegenStructure]):
    struct_names = {struct.f_name.lower() for struct in structs}

    missing_by_struct = {}
    for struct in structs:
        to_remove = []
        for idx, fld in enumerate(list(struct.arg)):
            if fld.type != "type" or fld.kind in params.structs_defined_externally:
                continue

            if fld.kind.lower() not in struct_names:
                logger.warning(
                    f"Missing definition for struct '{fld.kind}' which is used in '{struct.f_name}'"
                )
                missing_by_struct.setdefault(struct.f_name, {})
                missing_by_struct[struct.f_name][fld.f_name] = fld.kind

                to_remove.append(idx)
        for remove in sorted(to_remove, reverse=True):
            struct.arg.pop(remove)

    return missing_by_struct


def filter_structs(
    structs: list[CodegenStructure], allowed_modules: set[str], ignored_modules: set[str]
) -> list[CodegenStructure]:
    return [
        struct
        for struct in structs
        if struct.module in allowed_modules or struct.module not in ignored_modules
    ]


# def get_to_json_source(struct: CodegenStructure) -> list[str]:
#     args = [arg for arg in struct.args_to_convert if arg.is_component and arg.member is not None]
#
#     name_to_value = {arg.c_name: f"obj.{arg.c_name}" for arg in args}
#     # if struct.cpp_class == "CPP_ele":
#     #     name_to_value.pop("lord")
#     #     fixup_lines = [
#     #         "if (obj.lord.has_value()) {",
#     #         '    j["lord"] = json{*obj.lord.value()};',
#     #         "}",
#     #     ]
#     # else:
#     fixup_lines = []
#
#     members = ", ".join("{" + f'"{name}", {value}' + "}" for name, value in name_to_value.items())
#
#     return [
#         f"void to_json(json &j, const {struct.cpp_class} &obj) {{",
#         f"j = json {{ {members} }};",
#         *fixup_lines,
#         "}",
#         f"""
#         ostream &operator<<(ostream &os, const {struct.cpp_class} &obj) {{
#           json j;
#           to_json(j, obj);
#           std::string str = nlohmann::to_string(j);
#           os << str;
#           return os;
#         }}
#         """,
#     ]
#
#
# def write_cpp_json_source(file, structs: list[CodegenStructure]) -> None:
#     """Write C++ JSON serialization code."""
#     header_template = string.Template(
#         textwrap.dedent(
#             """\
#             //+
#             // C++ JSON helpers for Bmad / C++ structure interface.
#             //
#             // This file is generated as part of the Bmad/C++ interface code generation.
#             // The code generation files can be found in cpp_bmad_interface.
#             //
#             // DO NOT EDIT THIS FILE DIRECTLY!
#             //-
#
#             #include <iostream>
#             #include <memory>
#             #include <optional>
#
#             #include "bmad/convert.h"
#             #include "json.hpp"
#             ${include_headers}
#
#             using namespace Bmad;
#             using std::ostream;
#             using std::size_t;
#             using json = nlohmann::json;
#
#             namespace std {
#             template<typename T>
#             void to_json(json& j, const complex<T>& d) {
#                 j = {d.real(), d.imag()};
#             }
#             void from_json(const json& j, Complex &d) {
#                 d.real(j.at(0).get<double>());
#                 d.imag(j.at(1).get<double>());
#             }
#             } // namespace: std
#
#             namespace Bmad {
#             //--------------------------------------------------------------------
#             ${json_helpers}
#             //--------------------------------------------------------------------
#             } // namespace Bmad
#             """
#         )
#     )
#
#     include_headers = "\n".join(params.include_header_files)
#     json_helpers = "\n".join("\n".join(get_to_json_source(struct)) for struct in structs)
#     file.write(header_template.substitute(include_headers=include_headers, json_helpers=json_helpers))


def get_structure_definitions(
    params: CodegenConfig,
    parsed_structures: list[ParsedStructure],
) -> list[CodegenStructure]:
    structs: list[CodegenStructure] = []

    for name in params.struct_list:
        struct = CodegenStructure(name)
        structs.append(struct)
        match_structure_definition(parsed_structures, struct, params)
        struct.arg = [arg for arg in struct.arg if arg.should_translate(struct.f_name, params)]
    return structs


def write_proxy_classes(params: CodegenConfig, structs: list[CodegenStructure]) -> None:
    if DEBUG:
        write_parsed_structures(structs, "f_structs.parsed")

    generated = CPPBMAD_ROOT / "src" / "generated"

    write_if_differs(
        create_fortran_proxy_code,
        generated / "proxy_mod.f90",
        params,
        structs,
    )
    cpp_proxy_header_template = (CODEGEN_ROOT / "proxy.tpl.hpp").read_text()
    cpp_proxy_cpp_template = (CODEGEN_ROOT / "proxy.tpl.cpp").read_text()
    write_if_differs(
        create_cpp_proxy_header,
        CPPBMAD_ROOT / "include" / "bmad" / "generated" / "proxy.hpp",
        params,
        cpp_proxy_header_template,
        cpp_proxy_cpp_template,
        structs,
    )
    write_if_differs(
        create_cpp_proxy_impl,
        generated / "proxy.cpp",
        params,
        cpp_proxy_header_template,
        cpp_proxy_cpp_template,
        structs,
    )


def load_routines(parsed_structs: list[ParsedStructure], config: CodegenConfig):
    structs_seen = set()
    settings_and_routines = []
    for settings in config.routines:
        routines = parse_bmad_routines(settings, config)
        settings_and_routines.append((settings, routines))
        for routine in routines:
            if routine.private:
                continue

            structs_seen |= {arg.kind.lower() for arg in routine.args if arg.type == "type"}

    parsed_structs_by_name = {struct.name.lower(): struct for struct in parsed_structs}

    structs_to_use = set()
    structs_to_skip = set()
    check = list(structs_seen)
    while check:
        typ = check.pop(0)
        if typ not in parsed_structs_by_name:
            logger.warning("Skipping %s as it's not in the bmad parsed struture list", typ)
            structs_to_skip.add(typ)
            continue

        parsed_st = parsed_structs_by_name[typ]
        if parsed_st.module in params.skips or typ in params.skips:
            structs_to_skip.add(typ)
            continue

        structs_to_use.add(typ)
        for member in parsed_st.members.values():
            if member.type != "type":
                continue
            assert member.kind
            kind = member.kind.lower()
            if kind not in structs_to_skip and kind not in structs_to_use and kind not in check:
                check.append(kind)

    return settings_and_routines, structs_to_use


def generate(
    config_file: pathlib.Path = CODEGEN_ROOT / "default.toml",
    pybmad: bool = True,
):
    # TODO refactor globals
    global params  # noqa: PLW0603

    logger.info(f"Config file: {config_file}")
    params = CodegenConfig.from_file(config_file)

    parsed_structs = load_bmad_parser_structures()
    _settings_and_routines, _routine_structs = load_routines(parsed_structs, params)
    # TODO
    # params.struct_list = sorted(_routine_structs)

    structs = get_structure_definitions(params, parsed_structs)

    config_context.set(
        ConfigContext(
            params=params,
            codegen_structs=structs,
            parsed_structs=parsed_structs,
        )
    )
    assert len(config_context.get().codegen_structs)

    # Print diagnostics
    logger.info(f"Number of structs in input list: {len(structs)}")
    logger.info(f"Number of structs found:         {len(parsed_structs)}")

    missing_struct_attrs = check_missing(structs)
    write_proxy_classes(params, structs)
    write_if_differs(write_enums, ENUM_FILENAME)

    routines, routines_by_name = generate_routines(params)

    report_html = generate_coverage_report(routines, structs, missing_struct_attrs)
    write_contents_if_differs(REPO_ROOT / "coverage.html", report_html)

    to_string_header = generate_to_string_header(
        template=(CODEGEN_ROOT / "to_string.tpl.hpp").read_text(),
        structs=structs,
        routines=routines_by_name,
    )
    write_contents_if_differs(
        CPPBMAD_INCLUDE / "bmad" / "generated" / "to_string.hpp", contents=to_string_header
    )

    to_string_code = generate_to_string_code(
        template=(CODEGEN_ROOT / "to_string.tpl.cpp").read_text(),
        structs=structs,
        routines=routines_by_name,
    )
    write_contents_if_differs(CPPBMAD_SRC / "generated" / "to_string.cpp", contents=to_string_code)

    if pybmad:
        pybmad_files = generate_pybmad(structs, routines_by_name)
        # existing_files = set(OUTPUT_PATH.glob("*.cpp")) | set(OUTPUT_PATH.glob("*.hpp"))
        # for fn in existing_files:
        #     if fn.name not in module_file_to_source:
        #         logger.warning(f"Removing stale file from previous generation: {fn}")
        #         fn.unlink()
        for fn, source in pybmad_files.items():
            write_contents_if_differs(target_path=fn, contents=source)

    return structs, parsed_structs, routines, routines_by_name


def main():
    parser = argparse.ArgumentParser(description="Run the cppbmad code generator.")
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level.",
    )
    parser.add_argument(
        "--config-file",
        type=pathlib.Path,
        default=CODEGEN_ROOT / "default.toml",
        help="Path to the configuration file.",
    )
    parser.add_argument(
        "--no-pybmad",
        action="store_true",
        help="Do not generate pybmad pybind11 bindings",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug mode",
    )
    args = parser.parse_args()

    logging.basicConfig(level=args.log_level)
    logging.getLogger("codegen").setLevel(args.log_level)

    return generate(config_file=args.config_file, pybmad=not args.no_pybmad)


if __name__ == "__main__":
    structs, parsed_structs, routines, routines_by_name = main()
