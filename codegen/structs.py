from __future__ import annotations

import contextlib
import importlib
import logging
import pathlib
import sys

from .paths import ACC_ROOT_DIR


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


import bmad_struct_parser  # noqa: E402

json_path = pathlib.Path(bmad_struct_parser.__file__).resolve().absolute().parents[1] / "json"


from bmad_struct_parser import Structure as ParsedStructure  # noqa: E402
from bmad_struct_parser import load_structures_by_filename  # noqa: E402
from bmad_struct_parser.parser import StructureMember  # noqa: E402
from bmad_struct_parser.util import FileLine  # noqa: E402

__all__ = [
    "FileLine",
    "ParsedStructure",
    "StructureMember",
    "load_structures_by_filename",
]
