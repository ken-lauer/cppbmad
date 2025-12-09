from __future__ import annotations

from dataclasses import dataclass, fields

from .types import STANDARD_TYPES, FullType, Intent
from .util import struct_to_proxy_class_name


@dataclass
class Transform:
    # The C++ type for arguments
    cpp_type: str = ""
    # The C++ type when declared
    cpp_declare_type: str = ""
    # The C++ type used when calling the exported Fortran function
    cpp_call_fortran_type: str = ""
    # The native Fortran type
    fortran_type: str = ""
    # The default value for the C++ argument
    cpp_default: str = ""
    # The return type for the C++ function
    cpp_return_type: str = ""
    # The native Fortran type (if different from fortran_type)
    fortran_native_type: str = ""

    def is_reference(self):
        return self.is_optional_ref or "&" in self.cpp_type

    @property
    def is_std_optional(self):
        return self.cpp_type.startswith("std::optional<")

    @property
    def is_optional_ref(self):
        return self.cpp_type.startswith("optional_ref<")

    def replace_all(self, old: str, new: str) -> None:
        for fld in fields(self):
            value = getattr(self, fld.name)

            if isinstance(value, str):
                setattr(self, fld.name, value.replace(old, new))
            else:
                setattr(self, fld.name, [v.replace(old, new) for v in value])


def get_cpp_container(cpp_base: str, dim: int, ptr: str) -> str:
    """Determines standard Local C++ container syntax."""
    if dim == 0:
        return cpp_base if ptr == "NOT" else f"std::optional<{cpp_base}>"

    if cpp_base == "PROXYCLS" and ptr != "NOT":
        return f"PROXYCLSAllocatable{dim}D&"

    if ptr == "NOT":
        dims = ", ".join(f"DIM{i + 1}" for i in range(dim))
        return f"FixedArray{dim}D<{cpp_base}, {dims}>"

    return f"VariableArray{dim}D<{cpp_base}>"


class UnsupportedTypeError(Exception):
    pass


def get_type_transform(
    ft: FullType,
    intent: Intent = "in",
    is_optional: bool = False,
    kind: str = "",
    is_dynamic_array: bool = False,
) -> Transform:
    """
    Generates a specific Transform object for a given type configuration.

    Parameters
    ----------
    ft : FullType
        The full type description (type name, dimension, pointer status).
    intent : Intent
        The argument intent (in, out, inout).
    is_optional : bool
        Whether the argument is optional.
    kind : str
        The kind of the type (e.g. struct name), if applicable.

    Returns
    -------
    Transform
    """

    info = STANDARD_TYPES[ft.type]
    proxy_cls = None

    # Logic:
    # - Scalars (Fixed/NOT) are passed by Reference (&).
    # - Everything else (Arrays/Allocatable/Pointers) uses Array Typedefs.

    if ft.type == "type" or is_dynamic_array:
        cpp_call_type = "void*"
    else:
        # cpp_call_type = info.cpp_call_fortran_type or info.cpp_call_ref
        if ft.dim == 0 and ft.ptr != "PTR":
            cpp_call_type = info.cpp_call_ref
        else:
            cpp_call_type = info.cpp_call_arr

        if intent in {"in", "inout"} and is_optional:
            cpp_call_type = cpp_call_type.replace("&", "*")

    # --- C++ Type Logic ---
    if ft.type == "type":
        proxy_cls = struct_to_proxy_class_name(kind) if kind else "PROXYCLS"
        cpp_base = f"{proxy_cls}&"
    elif ft.type == "character":
        cpp_base = "std::string"
    else:
        cpp_base = info.c_type

    if ft.dim > 0:
        # Array handling
        if ft.dim == 1 and ft.ptr != "NOT" and info.allocatable_container:
            # Native allocatable array
            alloc_container = info.allocatable_container
            if is_optional:
                cpp_type = f"optional_ref<{alloc_container}>"
            else:
                cpp_type = f"{alloc_container}&"
        elif ft.type == "type" and ft.ptr != "NOT":
            alloc_container = f"{proxy_cls}Alloc{ft.dim}D"
            if is_optional:
                cpp_type = f"optional_ref<{alloc_container}>"
            else:
                cpp_type = f"{alloc_container}&"
        elif is_optional:
            # Optional array
            container = get_cpp_container(info.cpp_bmad_type, ft.dim, ft.ptr)
            cpp_type = f"std::optional<{container}>"
        else:
            # Standard array
            cpp_type = get_cpp_container(info.cpp_bmad_type, ft.dim, ft.ptr)
    elif ft.type == "type":
        # Scalar type
        cpp_type = cpp_base
        if is_optional:
            cpp_type = cpp_type.replace("&", "")
            cpp_type = f"optional_ref<{cpp_type}>"
    elif is_optional:
        # Optional scalar
        if intent == "inout":
            cpp_type = f"optional_ref<{cpp_base}>"
        else:
            cpp_type = f"std::optional<{cpp_base}>"
    else:
        # Standard scalar
        cpp_type = cpp_base

    # Clean up reference in optional
    if "std::optional" in cpp_type:
        cpp_type = cpp_type.replace("&", "")

    # --- C++ Default Value Logic ---
    cpp_default = ""
    if is_optional:
        if cpp_type.startswith(("optional_ref", "std::optional")):
            cpp_default = " = std::nullopt"
        else:
            cpp_default = " = nullptr"

    # --- C++ Return Type Logic ---
    cpp_return_type = cpp_base
    if ft.type == "type":
        assert proxy_cls is not None
        if intent in {"in", "inout"} or ft.dim == 0:
            cpp_return_type = proxy_cls
        elif ft.ptr == "NOT":
            cpp_return_type = f"{proxy_cls}Array{ft.dim}D"
        else:
            cpp_return_type = f"{proxy_cls}Alloc{ft.dim}D"
    elif ft.dim > 0:
        if is_dynamic_array:
            cpp_return_type = info.allocatable_container
            assert cpp_return_type
        else:
            cpp_return_type = get_cpp_container(info.cpp_bmad_type, ft.dim, ft.ptr)

    # --- Fortran Type Logic ---
    fortran_type = info.fortran_type
    # Note: Attributes like pointer, allocatable, value, intent are handled in fortran.py
    # based on the argument properties, but we can store the base type here.

    fortran_native_type = info.fortran_native_type or ""
    if ft.type == "type" and kind:
        fortran_native_type = f"type({kind})"

    def remove_optional(typ: str) -> str:
        if typ.startswith("optional_ref<"):
            return typ.removeprefix("optional_ref<")[:-1]
        if typ.startswith("std::optional"):
            return typ.removeprefix("std::optional<")[:-1]
        return typ

    return Transform(
        fortran_type=fortran_type,
        cpp_type=cpp_type,
        cpp_declare_type=remove_optional(cpp_type),
        cpp_call_fortran_type=cpp_call_type,
        cpp_default=cpp_default,
        cpp_return_type=cpp_return_type,
        fortran_native_type=fortran_native_type,
    )
