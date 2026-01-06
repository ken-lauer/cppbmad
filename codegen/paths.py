from __future__ import annotations

import os
import pathlib
import shutil


def normalize(path: str) -> pathlib.Path:
    path = os.path.expandvars(path)
    return pathlib.Path(path).expanduser().resolve().absolute()


if "ACC_ROOT_DIR" not in os.environ:
    raise RuntimeError("Environment variable ACC_ROOT_DIR is unset")

ACC_ROOT_DIR = normalize(os.path.expandvars(os.environ["ACC_ROOT_DIR"]))
CLANG_FORMAT_PATH = os.environ.get("CLANG_FORMAT_PATH", shutil.which("clang-format"))
CODEGEN_ROOT = pathlib.Path(__file__).resolve().absolute().parent
REPO_ROOT = CODEGEN_ROOT.parent
CPPBMAD_ROOT = REPO_ROOT
CPPBMAD_SRC = REPO_ROOT / "src"
CPPBMAD_INCLUDE = REPO_ROOT / "include"
PYBMAD_ROOT = REPO_ROOT / "python"
PYBMAD_INCLUDE = REPO_ROOT / "python" / "include"
PYBMAD_SRC = REPO_ROOT / "python" / "src"

DEFAULT_STRUCT_CONFIG_FILE = CODEGEN_ROOT / "struct_config.toml"
