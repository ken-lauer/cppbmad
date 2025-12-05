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
import sys
import tomllib
import typing

from .arg import Argument, CodegenStructure
from .config import CodegenConfig
from .context import ConfigContext, config_context
from .cpp import generate_routines_header
from .enums import ENUM_FILENAME, write_enums
from .paths import CODEGEN_ROOT, CPPBMAD_ROOT
from .proxy import (
    create_cpp_proxy_header,
    create_cpp_proxy_impl,
    create_fortran_proxy_code,
    struct_to_proxy_class_name,
)
from .routines import (
    generate_cpp_routine_code,
    generate_fortran_routine_code,
    parse_bmad_routines,
    prune_routines,
    routine_settings,
)
from .structs import ParsedStructure, load_bmad_parser_structures
from .types import (
    STRUCT,
)
from .util import write_contents_if_differs, write_if_differs

if typing.TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)

DEBUG = False  # Change to True to enable more verbose printout


##################################################################################
##################################################################################
def match_structure_definition(
    parsed_structures: list[ParsedStructure],
    struct: CodegenStructure,
):
    for fstruct in parsed_structures:
        if struct.f_name == fstruct.name:
            break
    else:
        raise RuntimeError(f"Structure not found: {struct.f_name!r}")

    struct.f_name = fstruct.name
    struct.short_name = fstruct.name.removesuffix("_struct")
    struct.cpp_class = struct_to_proxy_class_name(fstruct.name)
    struct.arg = [Argument.from_fstruct(fstruct, member) for member in fstruct.members.values()]
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
    # Report any structs not found
    missing_structs = [struct.f_name for struct in structs if struct.short_name == ""]
    for name in missing_structs:
        logger.error(f"NOT FOUND: {name}")

    # Exit if any structs are missing
    if missing_structs:
        sys.exit("COULD NOT FIND ALL THE STRUCTS! STOPPING HERE!")

    # Create set of defined struct names
    defined_struct_names = {struct.f_name for struct in structs}
    # Track all missing struct definitions
    missing_struct_definitions = []

    # Check that all referenced struct types have definitions
    for parent_struct in structs:
        for fld in parent_struct.arg:
            # Skip non-struct fields and externally defined structs
            if fld.type != STRUCT or fld.kind in params.structs_defined_externally:
                continue

            # Check if the struct type is defined
            if fld.kind not in defined_struct_names:
                missing_struct_definitions.append(
                    f"Missing definition for struct '{fld.kind}' which is used in '{parent_struct.short_name}'"
                )

    # Exit with error if any struct definitions are missing
    if missing_struct_definitions:
        for error_message in missing_struct_definitions:
            logger.error(error_message)
        sys.exit(1)


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
    params: CodegenConfig, parsed_structures: list[ParsedStructure]
) -> list[CodegenStructure]:
    structs: list[CodegenStructure] = []

    for name in params.struct_list:
        struct = CodegenStructure(name)
        match_structure_definition(parsed_structures, struct)
        struct.arg = [arg for arg in struct.arg if arg.should_translate(struct.f_name)]
        logger.debug(f"Struct: {struct}")
        structs.append(struct)
    return structs


def generate_routines():
    logger.info("Parsing routines")

    all_routines = []
    all_routines_by_name = {}
    to_write = {}
    for settings in routine_settings:
        # Filter out private routines to start with:
        routines = [routine for routine in parse_bmad_routines(settings) if not routine.private]
        all_routines.extend(routines)
        logger.info(f"Pruning routines ({settings.fortran_output_filename})")
        routines_by_name = prune_routines(routines)
        all_routines_by_name.update(routines_by_name)
        logger.info("Generating code: routines Fortran side")

        routines_header = generate_routines_header(
            template=(CODEGEN_ROOT / "routines.tpl.hpp").read_text(),
            routines=routines_by_name,
            settings=settings,
        )
        logger.info("Generating code: routines C++ side")
        cpp_routine_code = generate_cpp_routine_code(
            template=(CODEGEN_ROOT / "routines.tpl.cpp").read_text(),
            routines=routines_by_name,
            settings=settings,
        )
        fortran_routine_code = generate_fortran_routine_code(
            template=(CODEGEN_ROOT / "routines.tpl.f90").read_text(),
            routines=routines_by_name,
            settings=settings,
        )

        generated = CPPBMAD_ROOT / "src" / "generated"
        header_path = CPPBMAD_ROOT / "include" / "bmad" / "generated" / settings.cpp_header_filename
        to_write[header_path] = routines_header
        to_write[generated / settings.fortran_output_filename] = fortran_routine_code
        to_write[generated / settings.cpp_output_filename] = cpp_routine_code

    for fn, contents in to_write.items():
        write_contents_if_differs(
            target_path=fn,
            contents=contents,
        )

    unique_routines = {rt.name: rt for rt in all_routines}
    usable = [rt for rt in unique_routines.values() if rt.usable]
    logger.info("Usable procedures: %d / %d total", len(usable), len(all_routines))
    logger.info("Done")
    return all_routines, all_routines_by_name


def generate(config_file: pathlib.Path = CODEGEN_ROOT / "default.toml"):
    # TODO refactor globals
    global params  # noqa: PLW0603

    with config_file.open("rb") as fp:
        params = CodegenConfig(**tomllib.load(fp))

    logger.info(f"Config file: {config_file}")

    parsed_structs = load_bmad_parser_structures()
    config_context.set(
        ConfigContext(
            params=params,
            codegen_structs=[],
            parsed_structs=parsed_structs,
        )
    )

    # TODO this needs context.params for some reason
    structs = get_structure_definitions(params, parsed_structs)

    config_context.set(
        ConfigContext(
            params=params,
            codegen_structs=structs,
            parsed_structs=parsed_structs,
        )
    )
    assert len(config_context.get().codegen_structs)

    n_found = sum(1 for struct in structs if struct.short_name)

    # Print diagnostics
    logger.info(f"Number of structs in input list: {len(structs)}")
    logger.info(f"Number of structs found:         {n_found}")

    check_missing(structs)
    write_output(params, structs)
    write_if_differs(write_enums, ENUM_FILENAME)

    routines, routines_by_name = generate_routines()
    return structs, parsed_structs, routines, routines_by_name


def write_output(params: CodegenConfig, structs: list[CodegenStructure]) -> None:
    if DEBUG:
        write_parsed_structures(structs, "f_structs.parsed")

    generated = CPPBMAD_ROOT / "src" / "generated"

    # We aren't generating equality these now that it's moved out of bmad-ecosystem
    # bmad_structs = [struct for struct in structs if struct.parsed.filename.parts[1].lower() not in {"tao"}]
    # tao_structs = [struct for struct in structs if struct.parsed.filename.parts[1].lower() == "tao"]
    # write_if_differs(
    #     create_fortran_equality_check_code,
    #     ACC_ROOT_DIR / "bmad" / "modules" / "equality_mod.f90",
    #     bmad_structs,
    #     module_name="equality_mod",
    # )
    # write_if_differs(
    #     create_fortran_equality_check_code,
    #     # ACC_ROOT_DIR / "tao" / "src" / "tao_equality_mod.f90",
    #     generated / "tao_equality_mod.f90",
    #     tao_structs,
    #     module_name="tao_equality_mod",
    # )

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
        "--debug",
        action="store_true",
        help="Enable debug mode",
    )
    args = parser.parse_args()

    logging.basicConfig(level=args.log_level)
    logging.getLogger("codegen").setLevel(args.log_level)

    return generate(config_file=args.config_file)


if __name__ == "__main__":
    structs, parsed_structs, routines, routines_by_name = main()
