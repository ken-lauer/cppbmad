from __future__ import annotations

from dataclasses import dataclass, fields
from enum import Enum
from typing import Literal, NamedTuple

from .exceptions import UnsupportedTypeError
from .util import struct_to_proxy_class_name


class RoutineType(Enum):
    FUNCTION = "Function"
    SUBROUTINE = "Subroutine"
    UNKNOWN = "Unknown"


REAL = "real"
REAL16 = "real16"  # quad precision
CMPLX = "complex"
INT = "integer"
INT8 = "integer8"
LOGIC = "logical"
CHAR = "character"
STRUCT = "type"
ArgumentType = Literal["real", "real16", "complex", "integer", "integer8", "logical", "character", "type"]

NOT = "NOT"
PTR = "PTR"
ALLOC = "ALLOC"
PointerType = Literal["NOT", "PTR", "ALLOC"]

Intent = Literal["in", "out", "inout"]


class FullType(NamedTuple):
    type: ArgumentType
    dim: int
    ptr: PointerType

    def sort_key(self):
        return (self.dim, self.ptr, self.type)

    def __str__(self) -> str:
        return f"{self.dim}D_{self.ptr}_{self.type}"

    @staticmethod
    def from_template(type: str) -> FullType:
        dim, ptr, type_name = type.split("_")

        try:
            dim = int(dim.lower().rstrip("d"))
        except TypeError:
            raise ValueError(f"Dimension of type from template is not integer: {type=} {dim=}") from None

        if type_name not in (
            "real",
            "real16",
            "complex",
            "integer",
            "integer8",
            "logical",
            "character",
            "type",
        ):
            raise ValueError(f"Unexpected type: {type_name}")
        if ptr not in ("NOT", "PTR", "ALLOC"):
            raise ValueError(f"Invalid pointer type: {type=} {ptr=}")
        return FullType(type_name, dim, ptr)


@dataclass
class TypeInfo:
    """
    Information about a specific data type used in the code generation process.

    Parameters
    ----------
    name : ArgumentType
        The internal name of the type (e.g., 'real', 'integer').
    fortran_type : str
        The Fortran type declaration used in the generated Fortran wrapper (e.g., 'real(c_double)').
        This is used in the `bind(c)` interface.
    c_type : str
        The corresponding C/C++ type (e.g., 'double', 'int').
        This is used in the generated C++ header and implementation.
    cpp_bmad_type : str
        The C++ type alias or class name used in the Bmad C++ library (e.g., 'Real', 'Int').
        This is often a typedef or a class wrapper.
    cpp_call_ref : str
        The C++ type signature used when passing this type by reference to the Fortran wrapper.
        (e.g., 'c_Real&').
    cpp_call_arr : str
        The C++ type signature used when passing an array of this type to the Fortran wrapper.
        (e.g., 'c_RealArr').
    fortran_native_type : str | None, optional
        The native Fortran type, if different from `fortran_type`.
        Used when converting between C-compatible types and native Fortran types (e.g., 'real(rp)').
    cpp_call_fortran_type : str | None, optional
        The C++ type used specifically in the `extern "C"` declaration for the Fortran function,
        if different from the standard mapping.
    allocatable_container : str | None, optional
        The C++ container type used for allocatable arrays of this type (e.g., 'RealAlloc1D').
    """

    name: ArgumentType
    fortran_type: str
    c_type: str
    cpp_bmad_type: str
    cpp_call_ref: str
    cpp_call_arr: str
    fortran_native_type: str | None = None
    cpp_call_fortran_type: str | None = None
    allocatable_container: str | None = None


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


def remove_optional(typ: str) -> str:
    typ = typ.strip()
    if typ.startswith("optional_ref<"):
        return typ.removeprefix("optional_ref<").removesuffix(">")
    if typ.startswith("std::optional"):
        return typ.removeprefix("std::optional<").removesuffix(">")
    return typ


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
    # Standard scalar
    elif intent == "inout":
        cpp_type = f"{cpp_base}&"
    else:
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
            if not cpp_return_type:
                raise UnsupportedTypeError(f"{ft} {ft.dim} {is_dynamic_array=} allocatable container")
        else:
            cpp_return_type = get_cpp_container(info.cpp_bmad_type, ft.dim, ft.ptr)

    # --- Fortran Type Logic ---
    fortran_type = info.fortran_type
    # Note: Attributes like pointer, allocatable, value, intent are handled in fortran.py
    # based on the argument properties, but we can store the base type here.
    #

    fortran_native_type = info.fortran_native_type or ""
    if ft.type == "type" and kind:
        fortran_native_type = f"type({kind})"

    if intent == "out" and ft.ptr == PTR:
        # need to pass back a pointer to C++
        fortran_type = fortran_type.replace(", value", "")

    return Transform(
        fortran_type=fortran_type,
        cpp_type=cpp_type,
        cpp_declare_type=remove_optional(cpp_type),
        cpp_call_fortran_type=cpp_call_type,
        cpp_default=cpp_default,
        cpp_return_type=cpp_return_type,
        fortran_native_type=fortran_native_type,
    )


STANDARD_TYPES: dict[ArgumentType, TypeInfo] = {
    "integer": TypeInfo(
        name="integer",
        fortran_type="integer(c_int)",
        c_type="int",
        cpp_bmad_type="Int",
        cpp_call_ref="int&",
        cpp_call_arr="int*",
        fortran_native_type="integer",
        allocatable_container="IntAlloc1D",
    ),
    "integer8": TypeInfo(
        name="integer8",
        fortran_type="integer(c_int64_t)",
        c_type="int64_t",
        cpp_bmad_type="Int8",
        cpp_call_ref="int64_t&",
        cpp_call_arr="int64_t*",
        fortran_native_type="integer(8)",
        allocatable_container="Int8Alloc1D",
    ),
    "real": TypeInfo(
        name="real",
        fortran_type="real(c_double)",
        c_type="double",
        cpp_bmad_type="Real",
        cpp_call_ref="double&",
        cpp_call_arr="double*",
        fortran_native_type="real",
        allocatable_container="RealAlloc1D",
    ),
    "real16": TypeInfo(
        name="real16",
        fortran_type="real(c_long_double)",
        c_type="long double",  # std::float128_t in C++23
        cpp_bmad_type="long double",
        cpp_call_ref="long double&",
        cpp_call_arr="long double*",
        fortran_native_type="real(16)",
        allocatable_container="Real16Alloc1D",
    ),
    "complex": TypeInfo(
        name="complex",
        fortran_type="complex(c_double_complex)",
        c_type="std::complex<double>",
        cpp_bmad_type="Complex",
        cpp_call_ref="std::complex<double>&",
        cpp_call_arr="std::complex<double>*",
        fortran_native_type="complex",
        allocatable_container="ComplexAlloc1D",
    ),
    "logical": TypeInfo(
        name="logical",
        fortran_type="logical(c_bool)",
        c_type="bool",
        cpp_bmad_type="Bool",
        cpp_call_ref="bool&",
        cpp_call_arr="bool*",
        fortran_native_type="logical",
        allocatable_container="BoolAlloc1D",
    ),
    "character": TypeInfo(
        name="character",
        fortran_type="character(c_char)",
        c_type="std::string",
        cpp_bmad_type="string",
        cpp_call_ref="const char*",
        cpp_call_arr="const char**",
        fortran_native_type="character",
        cpp_call_fortran_type="char*",
    ),
    "type": TypeInfo(
        name="type",
        fortran_type="type(c_ptr), value",
        c_type="PROXYCLS",
        cpp_bmad_type="PROXYCLS",
        cpp_call_ref="PROXYCLS&",
        cpp_call_arr="void*",
        fortran_native_type="type",
        cpp_call_fortran_type="void*",
    ),
}


class ContainerTypeInfo(NamedTuple):
    name: str
    cpp_type: str
    cpp_container_name: str
    fortran_type: str
    fortran_container_struct: str


native_type_containers = [
    ContainerTypeInfo(
        name="real",
        cpp_type="double",
        cpp_container_name="RealAlloc1D",
        fortran_type="real(rp)",
        fortran_container_struct="real_container_alloc",
    ),
    ContainerTypeInfo(
        name="real16",
        cpp_type="long double",
        cpp_container_name="Real16Alloc1D",
        fortran_type="real(qp)",
        fortran_container_struct="real16_container_alloc",
    ),  # std::float128_t in C++23
    ContainerTypeInfo(
        name="integer",
        cpp_type="int",
        cpp_container_name="IntAlloc1D",
        fortran_type="integer",
        fortran_container_struct="integer_container_alloc",
    ),
    ContainerTypeInfo(
        name="integer8",
        cpp_type="int64_t",
        cpp_container_name="Int8Alloc1D",
        fortran_type="integer(8)",
        fortran_container_struct="integer8_container_alloc",
    ),
    ContainerTypeInfo(
        name="logical",
        cpp_type="bool",
        cpp_container_name="BoolAlloc1D",
        fortran_type="logical",
        fortran_container_struct="logical_container_alloc",
    ),
    ContainerTypeInfo(
        name="complex",
        cpp_type="std::complex<double>",
        cpp_container_name="ComplexAlloc1D",
        fortran_type="complex(rp)",
        fortran_container_struct="complex_container_alloc",
    ),
]
