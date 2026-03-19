from __future__ import annotations

import argparse
import logging
import pathlib
import sys

import pydantic

from .cpp import generate_routine_cpp_wrapper, generate_routine_return_value_struct
from .fortran import generate_fortran_routine_with_c_binding
from .gen import load_context
from .paths import CODEGEN_ROOT
from .routines import FortranRoutine

logger = logging.getLogger(__name__)


class GeneratedRoutine(pydantic.BaseModel):
    routine: FortranRoutine
    cpp: list[str]
    cpp_header: list[str]
    fortran: list[str]


def main():
    logging.basicConfig(level=logging.CRITICAL, stream=sys.stderr)

    parser = argparse.ArgumentParser(description="Inspect a Bmad routine and output JSON.")
    parser.add_argument("routine_name", help="The name of the routine to inspect.", nargs="+")
    parser.add_argument(
        "--config-file",
        type=pathlib.Path,
        default=CODEGEN_ROOT / "default.toml",
        help="Path to the configuration file.",
    )
    parser.add_argument(
        "--info",
        choices=(
            "json",
            "docstring",
        ),
        default="json",
        help="Information to inspect",
    )
    args = parser.parse_args()

    ctx = load_context(args.config_file)

    for name in args.routine_name:
        if name not in ctx.routines_by_name:
            # Check if it exists in all routines but unusable
            found = [r for r in ctx.routines if r.name == name]
            if found:
                print(f"Routine '{name}' found but marked unusable:", file=sys.stderr)
                for r in found:
                    print(f" - {r.filename}:{r.start_line.lineno}: {r.unusable_reason}", file=sys.stderr)
            else:
                print(f"Routine '{name}' not found in any form.", file=sys.stderr)
            continue
        rt = ctx.routines_by_name[name]
        cpp_header = []
        cpp_header.extend(generate_routine_return_value_struct(rt))

        decl = rt.get_cpp_decl(defaults=True, namespace=False)
        cpp_header.append(decl + ";")
        data = GeneratedRoutine(
            routine=rt,
            cpp=generate_routine_cpp_wrapper(rt),
            fortran=generate_fortran_routine_with_c_binding(rt).splitlines(),
            cpp_header=cpp_header,
        )

        if args.info == "json":
            print(data.model_dump_json(indent=2))
        elif args.info == "docstring":
            if rt.docstring is None:
                print("No associated docstring")
            else:
                print(rt.docstring.to_numpy_docstring())


if __name__ == "__main__":
    main()
