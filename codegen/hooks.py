"""
Generate the C ABI contract for the hand-written Bmad/Tao callback hooks.

Hooks are the reverse of routines:
Fortran calls *into* C++/Python via procedure pointers.

Each hook's argument list is read from the parsed ``<def_routine_name>``
abstract interface, so it tracks upstream Bmad rather than being hand-copied.
From that it emits:

* ``src/generated/hook_interfaces.f90`` -- module ``hook_c_interfaces`` with one
  ``ci_<name>`` abstract interface per hook, and
* ``include/bmad/generated/hook_abi.hpp`` -- the matching ``extern "C"`` register
  prototypes plus the ``std::function`` typedefs and set/clear declarations.

The trampolines (``src/hooks/*.f90`` / ``*.cpp``) and the nanobind marshalling
(``python/src/*_hooks.cpp``) are hand-written. If the definitions of the hooks
change, it will be a compilation error.
"""

from __future__ import annotations

import pathlib
import typing
from dataclasses import dataclass, field

from .context import CodegenConfig, HookSpec
from .paths import CPPBMAD_INCLUDE, CPPBMAD_SRC
from .util import snake_to_camel, struct_to_proxy_class_name

if typing.TYPE_CHECKING:
    from .routines import FortranRoutine, RoutineArg

FORTRAN_MODULE = "hook_c_interfaces"
FORTRAN_PATH = CPPBMAD_SRC / "generated" / "hook_interfaces.f90"
HPP_PATH = CPPBMAD_INCLUDE / "bmad" / "generated" / "hook_abi.hpp"


@dataclass
class ResolvedArg:
    """One interface argument resolved to its three per-layer renderings.

    ``fortran`` / ``c_abi`` may expand to several entries (a derived array
    crosses as base pointer + bounds + element size).
    """

    name: str  # argument name (matches the parsed interface)
    fortran: list[tuple[str, str]]  # (declaration_attr, name) e.g. ("type(c_ptr), value", "bunch")
    c_abi: list[str]  # extern "C" parameter types, e.g. ["void *"]
    cpp: str  # std::function parameter type as a by-reference out-param, e.g. "bool &"
    value_cpp: str = ""  # by-value C++ type when this arg is returned instead (e.g. "bool")
    imports: set[str] = field(default_factory=set)  # iso_c_binding names used


def _proxy(struct: str, namespace: str) -> str:
    """C++ proxy class name for a Fortran struct, qualified when outside Bmad."""
    name = struct_to_proxy_class_name(struct)
    return name if namespace == "Bmad" else f"Bmad::{name}"


def resolve_arg(arg: RoutineArg, namespace: str, inout_scalars: set[str]) -> ResolvedArg:
    """Map one parsed interface argument to its Fortran / C ABI / C++ renderings."""
    name = arg.c_name
    is_array = bool(getattr(arg, "array", None))

    if arg.type == "type":  # derived type -> opaque pointer / proxy
        proxy = _proxy(arg.kind, namespace)
        if is_array:
            # Assumed-shape derived array: base pointer + inclusive bounds + element size.
            return ResolvedArg(
                name=name,
                fortran=[
                    ("type(c_ptr), value", f"{name}_data"),
                    ("integer(c_int), value", f"{name}_lb"),
                    ("integer(c_int), value", f"{name}_ub"),
                    ("integer(c_size_t), value", f"{name}_esize"),
                ],
                c_abi=["void *", "int", "int", "std::size_t"],
                cpp=f"{proxy}Array1D {'*' if arg.is_optional else '&'}",
                imports={"c_ptr", "c_int", "c_size_t"},
            )
        return ResolvedArg(
            name=name,
            fortran=[("type(c_ptr), value", name)],
            c_abi=["void *"],
            cpp=f"{proxy} {'*' if arg.is_optional else '&'}",
            imports={"c_ptr"},
        )

    if arg.is_optional:
        # Optional scalar int/logical crosses as a nullable c_ptr (target temp).
        return ResolvedArg(
            name=name,
            fortran=[("type(c_ptr), value", name)],
            c_abi=["void *"],
            cpp="bool *" if arg.type == "logical" else "int *",
            imports={"c_ptr"},
        )

    if arg.type == "logical":
        # Logicals cross as integer(c_int) by reference: a by-ref out-param, or the
        # return value if named in the hook's `returns`.
        return ResolvedArg(
            name=name,
            fortran=[("integer(c_int)", name)],
            c_abi=["int *"],
            cpp="bool &",
            value_cpp="bool",
            imports={"c_int"},
        )

    by_ref = name in inout_scalars

    if arg.type == "real":
        return ResolvedArg(
            name=name,
            fortran=[("real(c_double)" if by_ref else "real(c_double), value", name)],
            c_abi=["double *" if by_ref else "double"],
            cpp="double &" if by_ref else "double",
            value_cpp="double",
            imports={"c_double"},
        )

    if arg.type in ("integer", "integer8"):
        return ResolvedArg(
            name=name,
            fortran=[("integer(c_int)" if by_ref else "integer(c_int), value", name)],
            c_abi=["int *" if by_ref else "int"],
            cpp="int &" if by_ref else "int",
            value_cpp="int",
            imports={"c_int"},
        )

    raise ValueError(f"Unsupported hook argument type {arg.type!r} for {name!r}")


def _resolved_args(hook: HookSpec, routines_by_name: dict[str, FortranRoutine]) -> list[ResolvedArg]:
    routine = routines_by_name.get(hook.def_routine_name)
    if routine is None:
        raise KeyError(
            f"Hook {hook.name!r}: parsed interface {hook.def_routine_name!r} not found. "
            f"Check `def_name`, and that its source is in a [[routines]] source_paths."
        )
    inout = set(hook.inout_scalars)
    return [resolve_arg(a, hook.cpp_namespace, inout) for a in routine.args]


def hook_type_name(hook: HookSpec) -> str:
    return hook.type_name or (snake_to_camel(hook.name) + "Hook")


# ---------------------------------------------------------------------------
# Fortran abstract interfaces
# ---------------------------------------------------------------------------


def _fortran_interface(hook: HookSpec, resolved: list[ResolvedArg]) -> str:
    params = [name for r in resolved for _, name in r.fortran]
    imports = sorted({imp for r in resolved for imp in r.imports})
    lines = [f"    subroutine ci_{hook.name}({', '.join(params)}) bind(c)"]
    if imports:
        lines.append(f"      import :: {', '.join(imports)}")
    for r in resolved:
        for attr, name in r.fortran:
            lines.append(f"      {attr} :: {name}")
    lines.append("    end subroutine")
    return "\n".join(lines)


def generate_fortran_interfaces(hooks: list[HookSpec], routines_by_name: dict[str, FortranRoutine]) -> str:
    bodies = "\n\n".join(_fortran_interface(h, _resolved_args(h, routines_by_name)) for h in hooks)
    return "\n".join(
        [
            "! Autogenerated; do not edit.",
            "",
            f"module {FORTRAN_MODULE}",
            "",
            "  use, intrinsic :: iso_c_binding",
            "  implicit none",
            "",
            "  abstract interface",
            "",
            bodies,
            "",
            "  end interface",
            "",
            f"end module {FORTRAN_MODULE}",
            "",
        ]
    )


# ---------------------------------------------------------------------------
# C++ header: extern "C" prototypes + std::function typedefs + set/clear decls
# ---------------------------------------------------------------------------


def _extern_c_prototype(hook: HookSpec, resolved: list[ResolvedArg]) -> str:
    c_types = [t for r in resolved for t in r.c_abi]
    return f"void {hook.register_c_name}(void (*)({', '.join(c_types)}));"


def _typedef_signature(hook: HookSpec, resolved: list[ResolvedArg]) -> str:
    """The ``std::function<...>`` type for one hook.

    Args named in ``hook.returns`` are passed to the callback by value and become
    the return value (a tuple if there is more than one). Every other out-param is
    written through a by-reference argument, so with no ``returns`` the callback
    returns ``void``.
    """
    returns = set(hook.returns)
    params: list[str] = []
    ret_types: list[str] = []
    for r in resolved:
        if r.name in returns:
            params.append(r.value_cpp)
            ret_types.append(r.value_cpp)
        else:
            params.append(r.cpp)

    if not ret_types:
        ret = "void"
    elif len(ret_types) == 1:
        ret = ret_types[0]
    else:
        ret = f"std::tuple<{', '.join(ret_types)}>"
    return f"std::function<{ret}({', '.join(params)})>"


def generate_hook_abi_header(hooks: list[HookSpec], routines_by_name: dict[str, FortranRoutine]) -> str:
    resolved_by_hook = {h.name: _resolved_args(h, routines_by_name) for h in hooks}
    externs = "\n".join(_extern_c_prototype(h, resolved_by_hook[h.name]) for h in hooks)

    namespaces: dict[str, list[HookSpec]] = {}
    for h in hooks:
        namespaces.setdefault(h.cpp_namespace, []).append(h)

    ns_blocks = []
    for ns, ns_hooks in namespaces.items():
        entries = []
        for h in ns_hooks:
            tn = hook_type_name(h)
            entries.append(
                f"using {tn} = {_typedef_signature(h, resolved_by_hook[h.name])};\n"
                f"void set_{h.name}_hook({tn} fn);\n"
                f"void clear_{h.name}_hook();"
            )
        ns_blocks.append(f"namespace {ns} {{\n\n" + "\n\n".join(entries) + f"\n\n}} // namespace {ns}")

    return "\n".join(
        [
            "#pragma once",
            "",
            "// Autogenerated; do not edit.",
            "",
            "#include <cstddef>",
            "#include <functional>",
            "",
            '#include "bmad/generated/proxy.hpp"',
            "",
            'extern "C" {',
            externs,
            "}",
            "",
            "\n\n".join(ns_blocks),
            "",
        ]
    )


def generate_hooks(
    config: CodegenConfig, routines_by_name: dict[str, FortranRoutine]
) -> dict[pathlib.Path, str]:
    """Return {path: contents} for the generated hook ABI files."""
    if not config.hooks:
        return {}
    return {
        FORTRAN_PATH: generate_fortran_interfaces(config.hooks, routines_by_name),
        HPP_PATH: generate_hook_abi_header(config.hooks, routines_by_name),
    }
