from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Literal, NamedTuple


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
SIZE = "size"
ArgumentType = Literal[
    "real", "real16", "complex", "integer", "integer8", "logical", "character", "type", "size"
]

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
            "size",
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


STANDARD_TYPES: dict[ArgumentType, TypeInfo] = {
    "integer": TypeInfo(
        name="integer",
        fortran_type="integer(c_int)",
        c_type="int",
        cpp_bmad_type="Int",
        cpp_call_ref="c_Int&",
        cpp_call_arr="c_IntArr",
        fortran_native_type="integer",
        allocatable_container="IntAlloc1D",
    ),
    "integer8": TypeInfo(
        name="integer8",
        fortran_type="integer(c_int64_t)",
        c_type="int64_t",
        cpp_bmad_type="Int8",
        cpp_call_ref="c_Int8&",
        cpp_call_arr="c_Int8Arr",
        fortran_native_type="integer(8)",
        allocatable_container="Int8Alloc1D",
    ),
    "real": TypeInfo(
        name="real",
        fortran_type="real(c_double)",
        c_type="double",
        cpp_bmad_type="Real",
        cpp_call_ref="c_Real&",
        cpp_call_arr="c_RealArr",
        fortran_native_type="real",
        allocatable_container="RealAlloc1D",
    ),
    "real16": TypeInfo(
        name="real16",
        fortran_type="real(c_long_double)",
        c_type="long double",  # std::float128_t in C++23
        cpp_bmad_type="long double",
        cpp_call_ref="const long double&",
        cpp_call_arr="const long double*",
        fortran_native_type="real(16)",
    ),
    "complex": TypeInfo(
        name="complex",
        fortran_type="complex(c_double_complex)",
        c_type="std::complex<double>",
        cpp_bmad_type="Complex",
        cpp_call_ref="c_Complex&",
        cpp_call_arr="c_ComplexArr",
        fortran_native_type="complex",
        allocatable_container="ComplexAlloc1D",
    ),
    "logical": TypeInfo(
        name="logical",
        fortran_type="logical(c_bool)",
        c_type="bool",
        cpp_bmad_type="Bool",
        cpp_call_ref="c_Bool&",
        cpp_call_arr="c_BoolArr",
        fortran_native_type="logical",
        allocatable_container="BoolAlloc1D",
    ),
    "character": TypeInfo(
        name="character",
        fortran_type="character(c_char)",
        c_type="std::string",
        cpp_bmad_type="string",
        cpp_call_ref="c_Char",
        cpp_call_arr="c_Char*",
        fortran_native_type="character",
        cpp_call_fortran_type="c_Char",
    ),
    "type": TypeInfo(
        name="type",
        fortran_type="type(c_ptr), value",
        c_type="PROXYCLS",
        cpp_bmad_type="PROXYCLS",
        cpp_call_ref="PROXYCLS&",
        cpp_call_arr="void*",
        fortran_native_type="character",
        cpp_call_fortran_type="void*",
    ),
}
