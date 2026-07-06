"""pybmad"""

from collections.abc import Iterator, Sequence
from typing import Annotated, overload

import numpy
from numpy.typing import NDArray

from . import (
    bmad as bmad,
    bsim as bsim,
    extra as extra,
    simutils as simutils,
    tao as tao,
    test as test
)


class RealArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: RealAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    @property
    def lower_bound(self) -> int: ...

    @property
    def upper_bound(self) -> int: ...

    @overload
    def __getitem__(self, arg: int, /) -> float: ...

    @overload
    def __getitem__(self, arg: slice, /) -> Annotated[NDArray[numpy.float64], dict(shape=(None,))]: ...

    @overload
    def __setitem__(self, arg0: int, arg1: float, /) -> None: ...

    @overload
    def __setitem__(self, arg0: slice, arg1: Sequence, /) -> None: ...

    def to_list(self) -> list[float]: ...

    def __iter__(self) -> Iterator[float]: ...

    def __array__(self, **kwargs) -> Annotated[NDArray[numpy.float64], dict(shape=(None,))]: ...

class RealAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: RealArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[float], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> RealArray1D: ...

    def __iter__(self) -> Iterator[float]: ...

    def __array__(self, **kwargs) -> Annotated[NDArray[numpy.float64], dict(shape=(None,))]: ...

    @overload
    def __getitem__(self, arg: int, /) -> float: ...

    @overload
    def __getitem__(self, arg: slice, /) -> Annotated[NDArray[numpy.float64], dict(shape=(None,))]: ...

    @overload
    def __setitem__(self, arg0: int, arg1: float, /) -> None: ...

    @overload
    def __setitem__(self, arg0: slice, arg1: Sequence, /) -> None: ...

class IntArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: IntAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    @property
    def lower_bound(self) -> int: ...

    @property
    def upper_bound(self) -> int: ...

    @overload
    def __getitem__(self, arg: int, /) -> int: ...

    @overload
    def __getitem__(self, arg: slice, /) -> Annotated[NDArray[numpy.int32], dict(shape=(None,))]: ...

    @overload
    def __setitem__(self, arg0: int, arg1: int, /) -> None: ...

    @overload
    def __setitem__(self, arg0: slice, arg1: Sequence, /) -> None: ...

    def to_list(self) -> list[int]: ...

    def __iter__(self) -> Iterator[int]: ...

    def __array__(self, **kwargs) -> Annotated[NDArray[numpy.int32], dict(shape=(None,))]: ...

class IntAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: IntArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[int], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> IntArray1D: ...

    def __iter__(self) -> Iterator[int]: ...

    def __array__(self, **kwargs) -> Annotated[NDArray[numpy.int32], dict(shape=(None,))]: ...

    @overload
    def __getitem__(self, arg: int, /) -> int: ...

    @overload
    def __getitem__(self, arg: slice, /) -> Annotated[NDArray[numpy.int32], dict(shape=(None,))]: ...

    @overload
    def __setitem__(self, arg0: int, arg1: int, /) -> None: ...

    @overload
    def __setitem__(self, arg0: slice, arg1: Sequence, /) -> None: ...

class Int8Array1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: Int8Alloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    @property
    def lower_bound(self) -> int: ...

    @property
    def upper_bound(self) -> int: ...

    @overload
    def __getitem__(self, arg: int, /) -> int: ...

    @overload
    def __getitem__(self, arg: slice, /) -> Annotated[NDArray[numpy.int64], dict(shape=(None,))]: ...

    @overload
    def __setitem__(self, arg0: int, arg1: int, /) -> None: ...

    @overload
    def __setitem__(self, arg0: slice, arg1: Sequence, /) -> None: ...

    def to_list(self) -> list[int]: ...

    def __iter__(self) -> Iterator[int]: ...

    def __array__(self, **kwargs) -> Annotated[NDArray[numpy.int64], dict(shape=(None,))]: ...

class Int8Alloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: Int8Array1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[int], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> Int8Array1D: ...

    def __iter__(self) -> Iterator[int]: ...

    def __array__(self, **kwargs) -> Annotated[NDArray[numpy.int64], dict(shape=(None,))]: ...

    @overload
    def __getitem__(self, arg: int, /) -> int: ...

    @overload
    def __getitem__(self, arg: slice, /) -> Annotated[NDArray[numpy.int64], dict(shape=(None,))]: ...

    @overload
    def __setitem__(self, arg0: int, arg1: int, /) -> None: ...

    @overload
    def __setitem__(self, arg0: slice, arg1: Sequence, /) -> None: ...

class BoolArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: BoolAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    @property
    def lower_bound(self) -> int: ...

    @property
    def upper_bound(self) -> int: ...

    @overload
    def __getitem__(self, arg: int, /) -> bool: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    @overload
    def __setitem__(self, arg0: int, arg1: bool, /) -> None: ...

    @overload
    def __setitem__(self, arg0: slice, arg1: Sequence, /) -> None: ...

    def to_list(self) -> list[bool]: ...

    def __iter__(self) -> object: ...

class BoolAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: BoolArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[bool], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> BoolArray1D: ...

    def __iter__(self) -> object: ...

    @overload
    def __getitem__(self, arg: int, /) -> bool: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    @overload
    def __setitem__(self, arg0: int, arg1: bool, /) -> None: ...

    @overload
    def __setitem__(self, arg0: slice, arg1: Sequence, /) -> None: ...

class ComplexArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: ComplexAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    @property
    def lower_bound(self) -> int: ...

    @property
    def upper_bound(self) -> int: ...

    @overload
    def __getitem__(self, arg: int, /) -> complex: ...

    @overload
    def __getitem__(self, arg: slice, /) -> Annotated[NDArray[numpy.complex128], dict(shape=(None,))]: ...

    @overload
    def __setitem__(self, arg0: int, arg1: complex, /) -> None: ...

    @overload
    def __setitem__(self, arg0: slice, arg1: Sequence, /) -> None: ...

    def to_list(self) -> list[complex]: ...

    def __iter__(self) -> Iterator[complex]: ...

    def __array__(self, **kwargs) -> Annotated[NDArray[numpy.complex128], dict(shape=(None,))]: ...

class ComplexAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: ComplexArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[complex], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> ComplexArray1D: ...

    def __iter__(self) -> Iterator[complex]: ...

    def __array__(self, **kwargs) -> Annotated[NDArray[numpy.complex128], dict(shape=(None,))]: ...

    @overload
    def __getitem__(self, arg: int, /) -> complex: ...

    @overload
    def __getitem__(self, arg: slice, /) -> Annotated[NDArray[numpy.complex128], dict(shape=(None,))]: ...

    @overload
    def __setitem__(self, arg0: int, arg1: complex, /) -> None: ...

    @overload
    def __setitem__(self, arg0: slice, arg1: Sequence, /) -> None: ...

class RealArray2D:
    def __init__(self) -> None: ...

    def is_valid(self) -> bool: ...

    def __len__(self) -> int: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    def to_list(self) -> list[float]: ...

    @property
    def total_size(self) -> int: ...

    def __array__(self, **kwargs) -> Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='F')]: ...

class IntArray2D:
    def __init__(self) -> None: ...

    def is_valid(self) -> bool: ...

    def __len__(self) -> int: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    def to_list(self) -> list[int]: ...

    @property
    def total_size(self) -> int: ...

    def __array__(self, **kwargs) -> Annotated[NDArray[numpy.int32], dict(shape=(None, None), order='F')]: ...

class Int8Array2D:
    def __init__(self) -> None: ...

    def is_valid(self) -> bool: ...

    def __len__(self) -> int: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    def to_list(self) -> list[int]: ...

    @property
    def total_size(self) -> int: ...

    def __array__(self, **kwargs) -> Annotated[NDArray[numpy.int64], dict(shape=(None, None), order='F')]: ...

class BoolArray2D:
    def __init__(self) -> None: ...

    def is_valid(self) -> bool: ...

    def __len__(self) -> int: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    def to_list(self) -> list[bool]: ...

    @property
    def total_size(self) -> int: ...

class ComplexArray2D:
    def __init__(self) -> None: ...

    def is_valid(self) -> bool: ...

    def __len__(self) -> int: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    def to_list(self) -> list[complex]: ...

    @property
    def total_size(self) -> int: ...

    def __array__(self, **kwargs) -> Annotated[NDArray[numpy.complex128], dict(shape=(None, None), order='F')]: ...

class RealArray3D:
    def __init__(self) -> None: ...

    def is_valid(self) -> bool: ...

    def __len__(self) -> int: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    def to_list(self) -> list[float]: ...

    @property
    def total_size(self) -> int: ...

    def __array__(self, **kwargs) -> Annotated[NDArray[numpy.float64], dict(shape=(None, None, None), order='F')]: ...

class IntArray3D:
    def __init__(self) -> None: ...

    def is_valid(self) -> bool: ...

    def __len__(self) -> int: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    def to_list(self) -> list[int]: ...

    @property
    def total_size(self) -> int: ...

    def __array__(self, **kwargs) -> Annotated[NDArray[numpy.int32], dict(shape=(None, None, None), order='F')]: ...

class Int8Array3D:
    def __init__(self) -> None: ...

    def is_valid(self) -> bool: ...

    def __len__(self) -> int: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    def to_list(self) -> list[int]: ...

    @property
    def total_size(self) -> int: ...

    def __array__(self, **kwargs) -> Annotated[NDArray[numpy.int64], dict(shape=(None, None, None), order='F')]: ...

class BoolArray3D:
    def __init__(self) -> None: ...

    def is_valid(self) -> bool: ...

    def __len__(self) -> int: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    def to_list(self) -> list[bool]: ...

    @property
    def total_size(self) -> int: ...

class ComplexArray3D:
    def __init__(self) -> None: ...

    def is_valid(self) -> bool: ...

    def __len__(self) -> int: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    def to_list(self) -> list[complex]: ...

    @property
    def total_size(self) -> int: ...

    def __array__(self, **kwargs) -> Annotated[NDArray[numpy.complex128], dict(shape=(None, None, None), order='F')]: ...

class FCharArray1D:
    def __init__(self) -> None: ...

    def is_valid(self) -> bool: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> "Bmad::FCharArray1D::StringProxy": ...

    def __setitem__(self, arg0: int, arg1: str, /) -> None: ...

    def __iter__(self) -> Iterator["Bmad::FCharArray1D::StringProxy"]: ...

    def to_list(self) -> list[str]: ...

    def __str__(self) -> str: ...

class CharacterAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int, str_len: int = 200) -> None: ...

    def resize(self, n: int, str_len: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def string_length(self) -> int: ...

    def to_list(self) -> list[str]: ...

    def __getitem__(self, arg: int, /) -> str: ...

    def __setitem__(self, arg0: int, arg1: str, /) -> None: ...

    def __iter__(self) -> object: ...

    def view(self) -> FCharArray1D: ...

    def __repr__(self) -> str: ...

    def __str__(self) -> str: ...

class AcKickerFreqStruct:
    """Fortran struct: ac_kicker_freq_struct"""

    def __init__(self, f: float | None = None, amp: float | None = None, phi: float | None = None) -> None: ...

    @property
    def f(self) -> float: ...

    @f.setter
    def f(self, arg: float, /) -> None: ...

    @property
    def amp(self) -> float: ...

    @amp.setter
    def amp(self, arg: float, /) -> None: ...

    @property
    def phi(self) -> float: ...

    @phi.setter
    def phi(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> AcKickerFreqStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> AcKickerFreqStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> AcKickerFreqStruct: ...

    def __deepcopy__(self, arg: dict, /) -> AcKickerFreqStruct: ...

    def __eq__(self, arg: AcKickerFreqStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class AcKickerStruct:
    """Fortran struct: ac_kicker_struct"""

    def __init__(self) -> None: ...

    @property
    def amp_vs_time(self) -> AcKickerTimeStructAlloc1D: ...

    @property
    def frequency(self) -> AcKickerFreqStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> AcKickerStruct: ...

    def __deepcopy__(self, arg: dict, /) -> AcKickerStruct: ...

    def __eq__(self, arg: AcKickerStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class AcKickerTimeStruct:
    """Fortran struct: ac_kicker_time_struct"""

    def __init__(self, amp: float | None = None, time: float | None = None, spline: SplineStruct | None = None) -> None: ...

    @property
    def amp(self) -> float: ...

    @amp.setter
    def amp(self, arg: float, /) -> None: ...

    @property
    def time(self) -> float: ...

    @time.setter
    def time(self, arg: float, /) -> None: ...

    @property
    def spline(self) -> SplineStruct: ...

    @spline.setter
    def spline(self, arg: SplineStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> AcKickerTimeStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> AcKickerTimeStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> AcKickerTimeStruct: ...

    def __deepcopy__(self, arg: dict, /) -> AcKickerTimeStruct: ...

    def __eq__(self, arg: AcKickerTimeStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class AllEncompassingStruct:
    """Fortran struct: all_encompassing_struct"""

    def __init__(self, real_rp_0d: float | None = None, real_rp_1d: Sequence[float] | None = None, real_rp_2d: Sequence[Sequence[float]] | None = None, real_rp_3d: Sequence[Sequence[Sequence[float]]] | None = None, real_rp_0d_ptr: float | None = None, real_rp_1d_ptr: Sequence[float] | None = None, real_rp_2d_ptr: Sequence[Sequence[float]] | None = None, real_rp_3d_ptr: Sequence[Sequence[Sequence[float]]] | None = None, real_rp_1d_alloc: Sequence[float] | None = None, real_rp_2d_alloc: Sequence[Sequence[float]] | None = None, real_rp_3d_alloc: Sequence[Sequence[Sequence[float]]] | None = None, real_dp_0d: float | None = None, real_dp_1d: Sequence[float] | None = None, real_dp_2d: Sequence[Sequence[float]] | None = None, real_dp_3d: Sequence[Sequence[Sequence[float]]] | None = None, real_dp_0d_ptr: float | None = None, real_dp_1d_ptr: Sequence[float] | None = None, real_dp_2d_ptr: Sequence[Sequence[float]] | None = None, real_dp_3d_ptr: Sequence[Sequence[Sequence[float]]] | None = None, real_dp_1d_alloc: Sequence[float] | None = None, real_dp_2d_alloc: Sequence[Sequence[float]] | None = None, real_dp_3d_alloc: Sequence[Sequence[Sequence[float]]] | None = None, complex_dp_0d: complex | None = None, complex_dp_1d: Sequence[complex] | None = None, complex_dp_2d: Sequence[Sequence[complex]] | None = None, complex_dp_3d: Sequence[Sequence[Sequence[complex]]] | None = None, complex_dp_0d_ptr: complex | None = None, complex_dp_1d_ptr: Sequence[complex] | None = None, complex_dp_2d_ptr: Sequence[Sequence[complex]] | None = None, complex_dp_3d_ptr: Sequence[Sequence[Sequence[complex]]] | None = None, complex_dp_1d_alloc: Sequence[complex] | None = None, complex_dp_2d_alloc: Sequence[Sequence[complex]] | None = None, complex_dp_3d_alloc: Sequence[Sequence[Sequence[complex]]] | None = None, int_0d: int | None = None, int_1d: Sequence[int] | None = None, int_2d: Sequence[Sequence[int]] | None = None, int_3d: Sequence[Sequence[Sequence[int]]] | None = None, int_0d_ptr: int | None = None, int_1d_ptr: Sequence[int] | None = None, int_2d_ptr: Sequence[Sequence[int]] | None = None, int_3d_ptr: Sequence[Sequence[Sequence[int]]] | None = None, int_1d_alloc: Sequence[int] | None = None, int_2d_alloc: Sequence[Sequence[int]] | None = None, int_3d_alloc: Sequence[Sequence[Sequence[int]]] | None = None, int8_0d: int | None = None, int8_1d: Sequence[int] | None = None, int8_2d: Sequence[Sequence[int]] | None = None, int8_3d: Sequence[Sequence[Sequence[int]]] | None = None, int8_0d_ptr: int | None = None, int8_1d_ptr: Sequence[int] | None = None, int8_2d_ptr: Sequence[Sequence[int]] | None = None, int8_3d_ptr: Sequence[Sequence[Sequence[int]]] | None = None, int8_1d_alloc: Sequence[int] | None = None, int8_2d_alloc: Sequence[Sequence[int]] | None = None, int8_3d_alloc: Sequence[Sequence[Sequence[int]]] | None = None, logical_0d: bool | None = None, logical_1d: Sequence[bool] | None = None, logical_2d: Sequence[Sequence[bool]] | None = None, logical_3d: Sequence[Sequence[Sequence[bool]]] | None = None, logical_0d_ptr: bool | None = None, logical_1d_ptr: Sequence[bool] | None = None, type_0d: TestSubStruct | None = None, type_0d_ptr: TestSubStruct | None = None) -> None: ...

    @property
    def real_rp_0d(self) -> float: ...

    @real_rp_0d.setter
    def real_rp_0d(self, arg: float, /) -> None: ...

    @property
    def real_rp_1d(self) -> RealArray1D: ...

    @real_rp_1d.setter
    def real_rp_1d(self, arg: Sequence[float], /) -> None: ...

    @property
    def real_rp_2d(self) -> RealArray2D: ...

    @real_rp_2d.setter
    def real_rp_2d(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def real_rp_3d(self) -> RealArray3D: ...

    @real_rp_3d.setter
    def real_rp_3d(self, arg: Sequence[Sequence[Sequence[float]]], /) -> None: ...

    @property
    def real_rp_0d_ptr(self) -> float | None: ...

    @real_rp_0d_ptr.setter
    def real_rp_0d_ptr(self, arg: float, /) -> None: ...

    @property
    def real_rp_1d_ptr(self) -> RealArray1D: ...

    @real_rp_1d_ptr.setter
    def real_rp_1d_ptr(self, arg: Sequence[float], /) -> None: ...

    @property
    def real_rp_2d_ptr(self) -> RealArray2D: ...

    @real_rp_2d_ptr.setter
    def real_rp_2d_ptr(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def real_rp_3d_ptr(self) -> RealArray3D: ...

    @real_rp_3d_ptr.setter
    def real_rp_3d_ptr(self, arg: Sequence[Sequence[Sequence[float]]], /) -> None: ...

    @property
    def real_rp_1d_alloc(self) -> RealAlloc1D: ...

    @real_rp_1d_alloc.setter
    def real_rp_1d_alloc(self, arg: Sequence[float], /) -> None: ...

    @property
    def real_rp_2d_alloc(self) -> RealArray2D: ...

    @real_rp_2d_alloc.setter
    def real_rp_2d_alloc(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def real_rp_3d_alloc(self) -> RealArray3D:
        """Real(dp)"""

    @real_rp_3d_alloc.setter
    def real_rp_3d_alloc(self, arg: Sequence[Sequence[Sequence[float]]], /) -> None: ...

    @property
    def real_dp_0d(self) -> float: ...

    @real_dp_0d.setter
    def real_dp_0d(self, arg: float, /) -> None: ...

    @property
    def real_dp_1d(self) -> RealArray1D: ...

    @real_dp_1d.setter
    def real_dp_1d(self, arg: Sequence[float], /) -> None: ...

    @property
    def real_dp_2d(self) -> RealArray2D: ...

    @real_dp_2d.setter
    def real_dp_2d(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def real_dp_3d(self) -> RealArray3D: ...

    @real_dp_3d.setter
    def real_dp_3d(self, arg: Sequence[Sequence[Sequence[float]]], /) -> None: ...

    @property
    def real_dp_0d_ptr(self) -> float | None: ...

    @real_dp_0d_ptr.setter
    def real_dp_0d_ptr(self, arg: float, /) -> None: ...

    @property
    def real_dp_1d_ptr(self) -> RealArray1D: ...

    @real_dp_1d_ptr.setter
    def real_dp_1d_ptr(self, arg: Sequence[float], /) -> None: ...

    @property
    def real_dp_2d_ptr(self) -> RealArray2D: ...

    @real_dp_2d_ptr.setter
    def real_dp_2d_ptr(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def real_dp_3d_ptr(self) -> RealArray3D: ...

    @real_dp_3d_ptr.setter
    def real_dp_3d_ptr(self, arg: Sequence[Sequence[Sequence[float]]], /) -> None: ...

    @property
    def real_dp_1d_alloc(self) -> RealAlloc1D: ...

    @real_dp_1d_alloc.setter
    def real_dp_1d_alloc(self, arg: Sequence[float], /) -> None: ...

    @property
    def real_dp_2d_alloc(self) -> RealArray2D: ...

    @real_dp_2d_alloc.setter
    def real_dp_2d_alloc(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def real_dp_3d_alloc(self) -> RealArray3D:
        """complex(dp)"""

    @real_dp_3d_alloc.setter
    def real_dp_3d_alloc(self, arg: Sequence[Sequence[Sequence[float]]], /) -> None: ...

    @property
    def complex_dp_0d(self) -> complex: ...

    @complex_dp_0d.setter
    def complex_dp_0d(self, arg: complex, /) -> None: ...

    @property
    def complex_dp_1d(self) -> ComplexArray1D: ...

    @complex_dp_1d.setter
    def complex_dp_1d(self, arg: Sequence[complex], /) -> None: ...

    @property
    def complex_dp_2d(self) -> ComplexArray2D: ...

    @complex_dp_2d.setter
    def complex_dp_2d(self, arg: Sequence[Sequence[complex]], /) -> None: ...

    @property
    def complex_dp_3d(self) -> ComplexArray3D: ...

    @complex_dp_3d.setter
    def complex_dp_3d(self, arg: Sequence[Sequence[Sequence[complex]]], /) -> None: ...

    @property
    def complex_dp_0d_ptr(self) -> complex | None: ...

    @complex_dp_0d_ptr.setter
    def complex_dp_0d_ptr(self, arg: complex, /) -> None: ...

    @property
    def complex_dp_1d_ptr(self) -> ComplexArray1D: ...

    @complex_dp_1d_ptr.setter
    def complex_dp_1d_ptr(self, arg: Sequence[complex], /) -> None: ...

    @property
    def complex_dp_2d_ptr(self) -> ComplexArray2D: ...

    @complex_dp_2d_ptr.setter
    def complex_dp_2d_ptr(self, arg: Sequence[Sequence[complex]], /) -> None: ...

    @property
    def complex_dp_3d_ptr(self) -> ComplexArray3D: ...

    @complex_dp_3d_ptr.setter
    def complex_dp_3d_ptr(self, arg: Sequence[Sequence[Sequence[complex]]], /) -> None: ...

    @property
    def complex_dp_1d_alloc(self) -> ComplexAlloc1D: ...

    @complex_dp_1d_alloc.setter
    def complex_dp_1d_alloc(self, arg: Sequence[complex], /) -> None: ...

    @property
    def complex_dp_2d_alloc(self) -> ComplexArray2D: ...

    @complex_dp_2d_alloc.setter
    def complex_dp_2d_alloc(self, arg: Sequence[Sequence[complex]], /) -> None: ...

    @property
    def complex_dp_3d_alloc(self) -> ComplexArray3D:
        """Integer"""

    @complex_dp_3d_alloc.setter
    def complex_dp_3d_alloc(self, arg: Sequence[Sequence[Sequence[complex]]], /) -> None: ...

    @property
    def int_0d(self) -> int: ...

    @int_0d.setter
    def int_0d(self, arg: int, /) -> None: ...

    @property
    def int_1d(self) -> IntArray1D: ...

    @int_1d.setter
    def int_1d(self, arg: Sequence[int], /) -> None: ...

    @property
    def int_2d(self) -> IntArray2D: ...

    @int_2d.setter
    def int_2d(self, arg: Sequence[Sequence[int]], /) -> None: ...

    @property
    def int_3d(self) -> IntArray3D: ...

    @int_3d.setter
    def int_3d(self, arg: Sequence[Sequence[Sequence[int]]], /) -> None: ...

    @property
    def int_0d_ptr(self) -> int | None: ...

    @int_0d_ptr.setter
    def int_0d_ptr(self, arg: int, /) -> None: ...

    @property
    def int_1d_ptr(self) -> IntArray1D: ...

    @int_1d_ptr.setter
    def int_1d_ptr(self, arg: Sequence[int], /) -> None: ...

    @property
    def int_2d_ptr(self) -> IntArray2D: ...

    @int_2d_ptr.setter
    def int_2d_ptr(self, arg: Sequence[Sequence[int]], /) -> None: ...

    @property
    def int_3d_ptr(self) -> IntArray3D: ...

    @int_3d_ptr.setter
    def int_3d_ptr(self, arg: Sequence[Sequence[Sequence[int]]], /) -> None: ...

    @property
    def int_1d_alloc(self) -> IntAlloc1D: ...

    @int_1d_alloc.setter
    def int_1d_alloc(self, arg: Sequence[int], /) -> None: ...

    @property
    def int_2d_alloc(self) -> IntArray2D: ...

    @int_2d_alloc.setter
    def int_2d_alloc(self, arg: Sequence[Sequence[int]], /) -> None: ...

    @property
    def int_3d_alloc(self) -> IntArray3D:
        """Integer8"""

    @int_3d_alloc.setter
    def int_3d_alloc(self, arg: Sequence[Sequence[Sequence[int]]], /) -> None: ...

    @property
    def int8_0d(self) -> int: ...

    @int8_0d.setter
    def int8_0d(self, arg: int, /) -> None: ...

    @property
    def int8_1d(self) -> Int8Array1D: ...

    @int8_1d.setter
    def int8_1d(self, arg: Sequence[int], /) -> None: ...

    @property
    def int8_2d(self) -> Int8Array2D: ...

    @int8_2d.setter
    def int8_2d(self, arg: Sequence[Sequence[int]], /) -> None: ...

    @property
    def int8_3d(self) -> Int8Array3D: ...

    @int8_3d.setter
    def int8_3d(self, arg: Sequence[Sequence[Sequence[int]]], /) -> None: ...

    @property
    def int8_0d_ptr(self) -> int | None: ...

    @int8_0d_ptr.setter
    def int8_0d_ptr(self, arg: int, /) -> None: ...

    @property
    def int8_1d_ptr(self) -> Int8Array1D: ...

    @int8_1d_ptr.setter
    def int8_1d_ptr(self, arg: Sequence[int], /) -> None: ...

    @property
    def int8_2d_ptr(self) -> Int8Array2D: ...

    @int8_2d_ptr.setter
    def int8_2d_ptr(self, arg: Sequence[Sequence[int]], /) -> None: ...

    @property
    def int8_3d_ptr(self) -> Int8Array3D: ...

    @int8_3d_ptr.setter
    def int8_3d_ptr(self, arg: Sequence[Sequence[Sequence[int]]], /) -> None: ...

    @property
    def int8_1d_alloc(self) -> Int8Alloc1D: ...

    @int8_1d_alloc.setter
    def int8_1d_alloc(self, arg: Sequence[int], /) -> None: ...

    @property
    def int8_2d_alloc(self) -> Int8Array2D: ...

    @int8_2d_alloc.setter
    def int8_2d_alloc(self, arg: Sequence[Sequence[int]], /) -> None: ...

    @property
    def int8_3d_alloc(self) -> Int8Array3D:
        """logical"""

    @int8_3d_alloc.setter
    def int8_3d_alloc(self, arg: Sequence[Sequence[Sequence[int]]], /) -> None: ...

    @property
    def logical_0d(self) -> bool: ...

    @logical_0d.setter
    def logical_0d(self, arg: bool, /) -> None: ...

    @property
    def logical_1d(self) -> BoolArray1D: ...

    @logical_1d.setter
    def logical_1d(self, arg: Sequence[bool], /) -> None: ...

    @property
    def logical_2d(self) -> BoolArray2D: ...

    @logical_2d.setter
    def logical_2d(self, arg: Sequence[Sequence[bool]], /) -> None: ...

    @property
    def logical_3d(self) -> BoolArray3D: ...

    @logical_3d.setter
    def logical_3d(self, arg: Sequence[Sequence[Sequence[bool]]], /) -> None: ...

    @property
    def logical_0d_ptr(self) -> bool | None: ...

    @logical_0d_ptr.setter
    def logical_0d_ptr(self, arg: bool, /) -> None: ...

    @property
    def logical_1d_ptr(self) -> BoolArray1D:
        """
        logical, pointer :: logical_2d_ptr(:,:) logical, pointer :: logical_3d_ptr(:,:,:) logical, allocatable :: logical_1d_alloc(:) logical, allocatable :: logical_2d_alloc(:,:) logical, allocatable :: logical_3d_alloc(:,:,:) type
        """

    @logical_1d_ptr.setter
    def logical_1d_ptr(self, arg: Sequence[bool], /) -> None: ...

    @property
    def type_0d(self) -> TestSubStruct: ...

    @type_0d.setter
    def type_0d(self, arg: TestSubStruct, /) -> None: ...

    @property
    def type_1d(self) -> TestSubStructArray1D: ...

    @property
    def type_2d(self) -> TestSubStructArray2D: ...

    @property
    def type_3d(self) -> TestSubStructArray3D: ...

    @property
    def type_0d_ptr(self) -> TestSubStruct | None: ...

    @type_0d_ptr.setter
    def type_0d_ptr(self, arg: TestSubStruct, /) -> None: ...

    @property
    def type_1d_ptr(self) -> TestSubStructArray1D: ...

    @property
    def type_2d_ptr(self) -> TestSubStructArray2D: ...

    @property
    def type_3d_ptr(self) -> TestSubStructArray3D: ...

    @property
    def type_1d_alloc(self) -> TestSubStructAlloc1D: ...

    @property
    def type_2d_alloc(self) -> TestSubStructArray2D: ...

    @property
    def type_3d_alloc(self) -> TestSubStructArray3D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> AllEncompassingStruct: ...

    def __deepcopy__(self, arg: dict, /) -> AllEncompassingStruct: ...

    def __eq__(self, arg: AllEncompassingStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class AllPointerStruct:
    """Fortran struct: all_pointer_struct"""

    def __init__(self, r: float | None = None, q: float | None = None, i: int | None = None, l: bool | None = None, r1: Sequence[float] | None = None, i1: Sequence[int] | None = None) -> None: ...

    @property
    def r(self) -> float | None: ...

    @r.setter
    def r(self, arg: float, /) -> None: ...

    @property
    def q(self) -> float | None: ...

    @q.setter
    def q(self, arg: float, /) -> None: ...

    @property
    def i(self) -> int | None: ...

    @i.setter
    def i(self, arg: int, /) -> None: ...

    @property
    def l(self) -> bool | None: ...

    @l.setter
    def l(self, arg: bool, /) -> None: ...

    @property
    def r1(self) -> RealArray1D: ...

    @r1.setter
    def r1(self, arg: Sequence[float], /) -> None: ...

    @property
    def i1(self) -> IntArray1D: ...

    @i1.setter
    def i1(self, arg: Sequence[int], /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> AllPointerStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> AllPointerStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> AllPointerStruct: ...

    def __deepcopy__(self, arg: dict, /) -> AllPointerStruct: ...

    def __eq__(self, arg: AllPointerStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class AnormalModeStruct:
    """Fortran struct: anormal_mode_struct"""

    def __init__(self, emittance: float | None = None, emittance_no_vert: float | None = None, synch_int: Sequence[float] | None = None, j_damp: float | None = None, alpha_damp: float | None = None, chrom: float | None = None, tune: float | None = None) -> None: ...

    @property
    def emittance(self) -> float:
        """Beam emittance (unnormalized). Includes vertical photon opening angle."""

    @emittance.setter
    def emittance(self, arg: float, /) -> None: ...

    @property
    def emittance_no_vert(self) -> float:
        """
        Unnormalized beam emittance without the vertical photon opening angle taken into account.
        """

    @emittance_no_vert.setter
    def emittance_no_vert(self, arg: float, /) -> None: ...

    @property
    def synch_int(self) -> RealArray1D:
        """Synchrotron integrals"""

    @synch_int.setter
    def synch_int(self, arg: Sequence[float], /) -> None: ...

    @property
    def j_damp(self) -> float:
        """damping partition number"""

    @j_damp.setter
    def j_damp(self, arg: float, /) -> None: ...

    @property
    def alpha_damp(self) -> float:
        """damping per turn"""

    @alpha_damp.setter
    def alpha_damp(self, arg: float, /) -> None: ...

    @property
    def chrom(self) -> float:
        """Chromaticity"""

    @chrom.setter
    def chrom(self, arg: float, /) -> None: ...

    @property
    def tune(self) -> float:
        """'Fractional' tune in radians"""

    @tune.setter
    def tune(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> AnormalModeStruct: ...

    def __deepcopy__(self, arg: dict, /) -> AnormalModeStruct: ...

    def __eq__(self, arg: AnormalModeStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ApertureParamStruct:
    """Fortran struct: aperture_param_struct"""

    def __init__(self, min_angle: float | None = None, max_angle: float | None = None, n_angle: int | None = None, n_turn: int | None = None, x_init: float | None = None, y_init: float | None = None, rel_accuracy: float | None = None, abs_accuracy: float | None = None, start_ele: str | None = None) -> None: ...

    @property
    def min_angle(self) -> float: ...

    @min_angle.setter
    def min_angle(self, arg: float, /) -> None: ...

    @property
    def max_angle(self) -> float: ...

    @max_angle.setter
    def max_angle(self, arg: float, /) -> None: ...

    @property
    def n_angle(self) -> int: ...

    @n_angle.setter
    def n_angle(self, arg: int, /) -> None: ...

    @property
    def n_turn(self) -> int:
        """Number of turns a particle must survive."""

    @n_turn.setter
    def n_turn(self, arg: int, /) -> None: ...

    @property
    def x_init(self) -> float:
        """Initial x coordinate to start with for theta_xy = 0."""

    @x_init.setter
    def x_init(self, arg: float, /) -> None: ...

    @property
    def y_init(self) -> float:
        """Initial y coordinate to start with for theta_xy = pi/2."""

    @y_init.setter
    def y_init(self, arg: float, /) -> None: ...

    @property
    def rel_accuracy(self) -> float:
        """Relative resolution of bracketed aperture."""

    @rel_accuracy.setter
    def rel_accuracy(self, arg: float, /) -> None: ...

    @property
    def abs_accuracy(self) -> float:
        """Absolute resolution of bracketed aperture (meters)."""

    @abs_accuracy.setter
    def abs_accuracy(self, arg: float, /) -> None: ...

    @property
    def start_ele(self) -> str:
        """Element to start tracking at."""

    @start_ele.setter
    def start_ele(self, arg: str, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ApertureParamStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ApertureParamStruct: ...

    def __eq__(self, arg: ApertureParamStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class AperturePointStruct:
    """Fortran struct: aperture_point_struct"""

    def __init__(self, x: float | None = None, y: float | None = None, plane: int | None = None, ix_ele: int | None = None, i_turn: int | None = None) -> None: ...

    @property
    def x(self) -> float:
        """(x,y) aperture point with respect to the reference orbit."""

    @x.setter
    def x(self, arg: float, /) -> None: ...

    @property
    def y(self) -> float:
        """(x,y) aperture point with respect to the reference orbit."""

    @y.setter
    def y(self, arg: float, /) -> None: ...

    @property
    def plane(self) -> int:
        """plane determining loss"""

    @plane.setter
    def plane(self, arg: int, /) -> None: ...

    @property
    def ix_ele(self) -> int:
        """ele index particle lost at"""

    @ix_ele.setter
    def ix_ele(self, arg: int, /) -> None: ...

    @property
    def i_turn(self) -> int:
        """turn particle lost at"""

    @i_turn.setter
    def i_turn(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> AperturePointStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> AperturePointStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> AperturePointStruct: ...

    def __deepcopy__(self, arg: dict, /) -> AperturePointStruct: ...

    def __eq__(self, arg: AperturePointStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ApertureScanStruct:
    """Fortran struct: aperture_scan_struct"""

    def __init__(self, ref_orb: CoordStruct | None = None, pz_start: float | None = None) -> None: ...

    @property
    def point(self) -> AperturePointStructAlloc1D:
        """Set of aperture points at different angles."""

    @property
    def ref_orb(self) -> CoordStruct:
        """Ref orbit around which the scan is made."""

    @ref_orb.setter
    def ref_orb(self, arg: CoordStruct, /) -> None: ...

    @property
    def pz_start(self) -> float:
        """Starting pz."""

    @pz_start.setter
    def pz_start(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> ApertureScanStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> ApertureScanStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ApertureScanStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ApertureScanStruct: ...

    def __eq__(self, arg: ApertureScanStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class AstraLatticeParamStruct:
    """Fortran struct: astra_lattice_param_struct"""

    def __init__(self, fieldmap_dimension: int | None = None) -> None: ...

    @property
    def fieldmap_dimension(self) -> int:
        """Dimensions for field map. 1 or 3"""

    @fieldmap_dimension.setter
    def fieldmap_dimension(self, arg: int, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> AstraLatticeParamStruct: ...

    def __deepcopy__(self, arg: dict, /) -> AstraLatticeParamStruct: ...

    def __eq__(self, arg: AstraLatticeParamStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class BaseLineEleStruct:
    """Fortran struct: base_line_ele_struct"""

    def __init__(self, name: str | None = None, tag: str | None = None, ix_multi: int | None = None, orientation: int | None = None, ix_ele_in_in_lat: int | None = None, ele_order_reflect: bool | None = None) -> None: ...

    @property
    def name(self) -> str:
        """Name of sequence or element"""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def tag(self) -> str:
        """Tag name."""

    @tag.setter
    def tag(self, arg: str, /) -> None: ...

    @property
    def ix_multi(self) -> int:
        """Multipass indentifier"""

    @ix_multi.setter
    def ix_multi(self, arg: int, /) -> None: ...

    @property
    def orientation(self) -> int:
        """Element reversed?"""

    @orientation.setter
    def orientation(self, arg: int, /) -> None: ...

    @property
    def ix_ele_in_in_lat(self) -> int: ...

    @ix_ele_in_in_lat.setter
    def ix_ele_in_in_lat(self, arg: int, /) -> None: ...

    @property
    def ele_order_reflect(self) -> bool:
        """Part of reflection or reversed line?"""

    @ele_order_reflect.setter
    def ele_order_reflect(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> BaseLineEleStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> BaseLineEleStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> BaseLineEleStruct: ...

    def __deepcopy__(self, arg: dict, /) -> BaseLineEleStruct: ...

    def __eq__(self, arg: BaseLineEleStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class BbuBeamStruct:
    """Fortran struct: bbu_beam_struct"""

    def __init__(self, ix_ele_bunch: Sequence[int] | None = None, ix_bunch_head: int | None = None, ix_bunch_end: int | None = None, n_bunch_in_lat: int | None = None, ix_stage_voltage_max: int | None = None, hom_voltage_max: float | None = None, time_now: float | None = None, one_turn_time: float | None = None, rf_wavelength_max: float | None = None) -> None: ...

    @property
    def bunch(self) -> BunchStructAlloc1D:
        """Bunches in the lattice"""

    @property
    def stage(self) -> BbuStageStructAlloc1D: ...

    @property
    def ix_ele_bunch(self) -> IntAlloc1D:
        """element where bunch is"""

    @ix_ele_bunch.setter
    def ix_ele_bunch(self, arg: Sequence[int], /) -> None: ...

    @property
    def ix_bunch_head(self) -> int:
        """Index to head bunch(:)"""

    @ix_bunch_head.setter
    def ix_bunch_head(self, arg: int, /) -> None: ...

    @property
    def ix_bunch_end(self) -> int:
        """Index of the end bunch(:). -1 -> no bunches."""

    @ix_bunch_end.setter
    def ix_bunch_end(self, arg: int, /) -> None: ...

    @property
    def n_bunch_in_lat(self) -> int:
        """Number of bunches transversing the lattice."""

    @n_bunch_in_lat.setter
    def n_bunch_in_lat(self, arg: int, /) -> None: ...

    @property
    def ix_stage_voltage_max(self) -> int: ...

    @ix_stage_voltage_max.setter
    def ix_stage_voltage_max(self, arg: int, /) -> None: ...

    @property
    def hom_voltage_max(self) -> float: ...

    @hom_voltage_max.setter
    def hom_voltage_max(self, arg: float, /) -> None: ...

    @property
    def time_now(self) -> float: ...

    @time_now.setter
    def time_now(self, arg: float, /) -> None: ...

    @property
    def one_turn_time(self) -> float: ...

    @one_turn_time.setter
    def one_turn_time(self, arg: float, /) -> None: ...

    @property
    def rf_wavelength_max(self) -> float: ...

    @rf_wavelength_max.setter
    def rf_wavelength_max(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> BbuBeamStruct: ...

    def __deepcopy__(self, arg: dict, /) -> BbuBeamStruct: ...

    def __eq__(self, arg: BbuBeamStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class BbuParamStruct:
    """Fortran struct: bbu_param_struct"""

    def __init__(self, lat_filename: str | None = None, lat2_filename: str | None = None, bunch_by_bunch_info_file: str | None = None, hybridize: bool | None = None, write_digested_hybrid_lat: bool | None = None, write_voltage_vs_time_dat: bool | None = None, keep_overlays_and_groups: bool | None = None, keep_all_lcavities: bool | None = None, use_taylor_for_hybrids: bool | None = None, stable_orbit_anal: bool | None = None, limit_factor: float | None = None, simulation_turns_max: float | None = None, bunch_freq: float | None = None, init_particle_offset: float | None = None, current: float | None = None, rel_tol: float | None = None, drscan: bool | None = None, use_interpolated_threshold: bool | None = None, write_hom_info: bool | None = None, elindex: int | None = None, elname: str | None = None, nstep: int | None = None, begdr: float | None = None, enddr: float | None = None, nrep: int | None = None, ran_seed: int | None = None, hom_order_cutoff: int | None = None, ran_gauss_sigma_cut: float | None = None, ele_track_end: str | None = None, ix_ele_track_end: int | None = None, regression: bool | None = None, normalize_z_to_rf: bool | None = None, ramp_on: bool | None = None, ramp_pattern: Sequence[float] | None = None, ramp_n_start: int | None = None, n_ramp_pattern: int | None = None) -> None: ...

    @property
    def lat_filename(self) -> str:
        """Bmad lattice file name"""

    @lat_filename.setter
    def lat_filename(self, arg: str, /) -> None: ...

    @property
    def lat2_filename(self) -> str:
        """Bmad lattice2 file name for secondary parser"""

    @lat2_filename.setter
    def lat2_filename(self, arg: str, /) -> None: ...

    @property
    def bunch_by_bunch_info_file(self) -> str:
        """For outputting bunch-by-bunch info."""

    @bunch_by_bunch_info_file.setter
    def bunch_by_bunch_info_file(self, arg: str, /) -> None: ...

    @property
    def hybridize(self) -> bool:
        """Combine non-hom elements to speed up simulation?"""

    @hybridize.setter
    def hybridize(self, arg: bool, /) -> None: ...

    @property
    def write_digested_hybrid_lat(self) -> bool:
        """For debugging purposes."""

    @write_digested_hybrid_lat.setter
    def write_digested_hybrid_lat(self, arg: bool, /) -> None: ...

    @property
    def write_voltage_vs_time_dat(self) -> bool:
        """For debugging purposes."""

    @write_voltage_vs_time_dat.setter
    def write_voltage_vs_time_dat(self, arg: bool, /) -> None: ...

    @property
    def keep_overlays_and_groups(self) -> bool:
        """Keep when hybridizing?"""

    @keep_overlays_and_groups.setter
    def keep_overlays_and_groups(self, arg: bool, /) -> None: ...

    @property
    def keep_all_lcavities(self) -> bool:
        """Keep when hybridizing?"""

    @keep_all_lcavities.setter
    def keep_all_lcavities(self, arg: bool, /) -> None: ...

    @property
    def use_taylor_for_hybrids(self) -> bool:
        """
        Use taylor map for hybrids when true. Otherwise tracking method is linear.
        """

    @use_taylor_for_hybrids.setter
    def use_taylor_for_hybrids(self, arg: bool, /) -> None: ...

    @property
    def stable_orbit_anal(self) -> bool:
        """Write stable_orbit.out and hom_voltage.out?"""

    @stable_orbit_anal.setter
    def stable_orbit_anal(self, arg: bool, /) -> None: ...

    @property
    def limit_factor(self) -> float:
        """Init_hom_amp * limit_factor = simulation unstable limit"""

    @limit_factor.setter
    def limit_factor(self, arg: float, /) -> None: ...

    @property
    def simulation_turns_max(self) -> float:
        """Sets the duration of the simulation."""

    @simulation_turns_max.setter
    def simulation_turns_max(self, arg: float, /) -> None: ...

    @property
    def bunch_freq(self) -> float:
        """Freq in Hz."""

    @bunch_freq.setter
    def bunch_freq(self, arg: float, /) -> None: ...

    @property
    def init_particle_offset(self) -> float:
        """Initial particle offset for particles born in the first turn period."""

    @init_particle_offset.setter
    def init_particle_offset(self, arg: float, /) -> None: ...

    @property
    def current(self) -> float:
        """Starting current (amps)"""

    @current.setter
    def current(self, arg: float, /) -> None: ...

    @property
    def rel_tol(self) -> float:
        """Final threshold current accuracy."""

    @rel_tol.setter
    def rel_tol(self, arg: float, /) -> None: ...

    @property
    def drscan(self) -> bool:
        """If true, scan DR variable as in PRSTAB 7 (2004) Fig. 3."""

    @drscan.setter
    def drscan(self, arg: bool, /) -> None: ...

    @property
    def use_interpolated_threshold(self) -> bool: ...

    @use_interpolated_threshold.setter
    def use_interpolated_threshold(self, arg: bool, /) -> None: ...

    @property
    def write_hom_info(self) -> bool:
        """Write HOM parameters to main output file?"""

    @write_hom_info.setter
    def write_hom_info(self, arg: bool, /) -> None: ...

    @property
    def elindex(self) -> int: ...

    @elindex.setter
    def elindex(self, arg: int, /) -> None: ...

    @property
    def elname(self) -> str:
        """Element to step length for DRSCAN"""

    @elname.setter
    def elname(self, arg: str, /) -> None: ...

    @property
    def nstep(self) -> int:
        """Number of steps for DRSCAN."""

    @nstep.setter
    def nstep(self, arg: int, /) -> None: ...

    @property
    def begdr(self) -> float:
        """Beginning DR value for DRSCAN."""

    @begdr.setter
    def begdr(self, arg: float, /) -> None: ...

    @property
    def enddr(self) -> float:
        """End DR value for DRSCAN."""

    @enddr.setter
    def enddr(self, arg: float, /) -> None: ...

    @property
    def nrep(self) -> int:
        """Number of times to repeat threshold calculation"""

    @nrep.setter
    def nrep(self, arg: int, /) -> None: ...

    @property
    def ran_seed(self) -> int:
        """If set to 0, the output results will vary from run to run."""

    @ran_seed.setter
    def ran_seed(self, arg: int, /) -> None: ...

    @property
    def hom_order_cutoff(self) -> int:
        """If positive -> ignore HOM's with order greater than this."""

    @hom_order_cutoff.setter
    def hom_order_cutoff(self, arg: int, /) -> None: ...

    @property
    def ran_gauss_sigma_cut(self) -> float: ...

    @ran_gauss_sigma_cut.setter
    def ran_gauss_sigma_cut(self, arg: float, /) -> None: ...

    @property
    def ele_track_end(self) -> str: ...

    @ele_track_end.setter
    def ele_track_end(self, arg: str, /) -> None: ...

    @property
    def ix_ele_track_end(self) -> int:
        """Default: set to last element with a wake"""

    @ix_ele_track_end.setter
    def ix_ele_track_end(self, arg: int, /) -> None: ...

    @property
    def regression(self) -> bool:
        """Do regression test?"""

    @regression.setter
    def regression(self, arg: bool, /) -> None: ...

    @property
    def normalize_z_to_rf(self) -> bool:
        """make starting z = mod(z, rf_wavelength)? Ramp parameters"""

    @normalize_z_to_rf.setter
    def normalize_z_to_rf(self, arg: bool, /) -> None: ...

    @property
    def ramp_on(self) -> bool: ...

    @ramp_on.setter
    def ramp_on(self, arg: bool, /) -> None: ...

    @property
    def ramp_pattern(self) -> RealArray1D: ...

    @ramp_pattern.setter
    def ramp_pattern(self, arg: Sequence[float], /) -> None: ...

    @property
    def ramp_n_start(self) -> int:
        """Index of start of ramp Internal parameters"""

    @ramp_n_start.setter
    def ramp_n_start(self, arg: int, /) -> None: ...

    @property
    def n_ramp_pattern(self) -> int:
        """Number of valid ramp_pattern"""

    @n_ramp_pattern.setter
    def n_ramp_pattern(self, arg: int, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> BbuParamStruct: ...

    def __deepcopy__(self, arg: dict, /) -> BbuParamStruct: ...

    def __eq__(self, arg: BbuParamStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class BbuStageStruct:
    """Fortran struct: bbu_stage_struct"""

    def __init__(self, ix_ele_lr_wake: int | None = None, ix_ele_stage_end: int | None = None, ix_pass: int | None = None, ix_stage_pass1: int | None = None, ix_head_bunch: int | None = None, ix_hom_max: int | None = None, hom_voltage_max: float | None = None, time_at_wake_ele: float | None = None, ave_orb: Sequence[float] | None = None, rms_orb: Sequence[float] | None = None, min_orb: Sequence[float] | None = None, max_orb: Sequence[float] | None = None, n_orb: int | None = None) -> None: ...

    @property
    def ix_ele_lr_wake(self) -> int:
        """Element index of element with the wake"""

    @ix_ele_lr_wake.setter
    def ix_ele_lr_wake(self, arg: int, /) -> None: ...

    @property
    def ix_ele_stage_end(self) -> int:
        """Element at end of stage."""

    @ix_ele_stage_end.setter
    def ix_ele_stage_end(self, arg: int, /) -> None: ...

    @property
    def ix_pass(self) -> int:
        """Pass index when in multipass section"""

    @ix_pass.setter
    def ix_pass(self, arg: int, /) -> None: ...

    @property
    def ix_stage_pass1(self) -> int:
        """Index of corresponding stage on first pass"""

    @ix_stage_pass1.setter
    def ix_stage_pass1(self, arg: int, /) -> None: ...

    @property
    def ix_head_bunch(self) -> int: ...

    @ix_head_bunch.setter
    def ix_head_bunch(self, arg: int, /) -> None: ...

    @property
    def ix_hom_max(self) -> int: ...

    @ix_hom_max.setter
    def ix_hom_max(self, arg: int, /) -> None: ...

    @property
    def hom_voltage_max(self) -> float: ...

    @hom_voltage_max.setter
    def hom_voltage_max(self, arg: float, /) -> None: ...

    @property
    def time_at_wake_ele(self) -> float: ...

    @time_at_wake_ele.setter
    def time_at_wake_ele(self, arg: float, /) -> None: ...

    @property
    def ave_orb(self) -> RealArray1D: ...

    @ave_orb.setter
    def ave_orb(self, arg: Sequence[float], /) -> None: ...

    @property
    def rms_orb(self) -> RealArray1D: ...

    @rms_orb.setter
    def rms_orb(self, arg: Sequence[float], /) -> None: ...

    @property
    def min_orb(self) -> RealArray1D: ...

    @min_orb.setter
    def min_orb(self, arg: Sequence[float], /) -> None: ...

    @property
    def max_orb(self) -> RealArray1D: ...

    @max_orb.setter
    def max_orb(self, arg: Sequence[float], /) -> None: ...

    @property
    def n_orb(self) -> int: ...

    @n_orb.setter
    def n_orb(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> BbuStageStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> BbuStageStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> BbuStageStruct: ...

    def __deepcopy__(self, arg: dict, /) -> BbuStageStruct: ...

    def __eq__(self, arg: BbuStageStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class BeamInitStruct:
    """Fortran struct: beam_init_struct"""

    def __init__(self, position_file: str | None = None, spin: Sequence[float] | None = None, KV: KvBeamInitStruct | None = None, center_jitter: Sequence[float] | None = None, emit_jitter: Sequence[float] | None = None, sig_z_jitter: float | None = None, sig_pz_jitter: float | None = None, n_particle: int | None = None, renorm_center: bool | None = None, renorm_sigma: bool | None = None, random_engine: str | None = None, random_gauss_converter: str | None = None, random_sigma_cutoff: float | None = None, a_norm_emit: float | None = None, b_norm_emit: float | None = None, a_emit: float | None = None, b_emit: float | None = None, dPz_dz: float | None = None, center: Sequence[float] | None = None, t_offset: float | None = None, dt_bunch: float | None = None, sig_z: float | None = None, sig_pz: float | None = None, bunch_charge: float | None = None, n_bunch: int | None = None, ix_turn: int | None = None, species: str | None = None, full_6D_coupling_calc: bool | None = None, use_particle_start: bool | None = None, use_t_coords: bool | None = None, use_z_as_t: bool | None = None, file_name: str | None = None) -> None: ...

    @property
    def position_file(self) -> str:
        """File with particle positions."""

    @position_file.setter
    def position_file(self, arg: str, /) -> None: ...

    @property
    def distribution_type(self) -> FCharArray1D:
        r"""
        distribution type (in x-px, y-py, and z-pz planes) 'ELLIPSE', 'KV', 'GRID', 'FILE', 'RAN_GAUSS' or \'\' = 'RAN_GAUSS'
        """

    @property
    def spin(self) -> RealArray1D:
        """Spin (x, y, z)"""

    @spin.setter
    def spin(self, arg: Sequence[float], /) -> None: ...

    @property
    def ellipse(self) -> EllipseBeamInitStructArray1D:
        """Ellipse beam distribution"""

    @property
    def KV(self) -> KvBeamInitStruct:
        """KV beam distribution"""

    @KV.setter
    def KV(self, arg: KvBeamInitStruct, /) -> None: ...

    @property
    def grid(self) -> GridBeamInitStructArray1D:
        """Grid beam distribution"""

    @property
    def center_jitter(self) -> RealArray1D:
        """Bunch center rms jitter"""

    @center_jitter.setter
    def center_jitter(self, arg: Sequence[float], /) -> None: ...

    @property
    def emit_jitter(self) -> RealArray1D:
        """a and b bunch emittance rms jitter normalized to emittance"""

    @emit_jitter.setter
    def emit_jitter(self, arg: Sequence[float], /) -> None: ...

    @property
    def sig_z_jitter(self) -> float:
        """bunch length RMS jitter"""

    @sig_z_jitter.setter
    def sig_z_jitter(self, arg: float, /) -> None: ...

    @property
    def sig_pz_jitter(self) -> float:
        """RMS pz spread jitter"""

    @sig_pz_jitter.setter
    def sig_pz_jitter(self, arg: float, /) -> None: ...

    @property
    def n_particle(self) -> int:
        """Number of particles per bunch."""

    @n_particle.setter
    def n_particle(self, arg: int, /) -> None: ...

    @property
    def renorm_center(self) -> bool:
        """Renormalize centroid?"""

    @renorm_center.setter
    def renorm_center(self, arg: bool, /) -> None: ...

    @property
    def renorm_sigma(self) -> bool:
        """Renormalize sigma?"""

    @renorm_sigma.setter
    def renorm_sigma(self, arg: bool, /) -> None: ...

    @property
    def random_engine(self) -> str:
        """Or 'quasi'. Random number engine to use."""

    @random_engine.setter
    def random_engine(self, arg: str, /) -> None: ...

    @property
    def random_gauss_converter(self) -> str:
        """Or 'quick' or 'exact'. Uniform to gauss conversion method."""

    @random_gauss_converter.setter
    def random_gauss_converter(self, arg: str, /) -> None: ...

    @property
    def random_sigma_cutoff(self) -> float:
        """Cut-off in sigmas."""

    @random_sigma_cutoff.setter
    def random_sigma_cutoff(self, arg: float, /) -> None: ...

    @property
    def a_norm_emit(self) -> float:
        """a-mode normalized emittance (emit * beta * gamma)"""

    @a_norm_emit.setter
    def a_norm_emit(self, arg: float, /) -> None: ...

    @property
    def b_norm_emit(self) -> float:
        """b-mode normalized emittance (emit * beta * gamma)"""

    @b_norm_emit.setter
    def b_norm_emit(self, arg: float, /) -> None: ...

    @property
    def a_emit(self) -> float:
        """a-mode emittance"""

    @a_emit.setter
    def a_emit(self, arg: float, /) -> None: ...

    @property
    def b_emit(self) -> float:
        """b-mode emittance"""

    @b_emit.setter
    def b_emit(self, arg: float, /) -> None: ...

    @property
    def dPz_dz(self) -> float:
        """Correlation of Pz with long position."""

    @dPz_dz.setter
    def dPz_dz(self, arg: float, /) -> None: ...

    @property
    def center(self) -> RealArray1D:
        """Bench phase space center offset relative to reference."""

    @center.setter
    def center(self, arg: Sequence[float], /) -> None: ...

    @property
    def t_offset(self) -> float:
        """Time center offset"""

    @t_offset.setter
    def t_offset(self, arg: float, /) -> None: ...

    @property
    def dt_bunch(self) -> float:
        """Time between bunches."""

    @dt_bunch.setter
    def dt_bunch(self, arg: float, /) -> None: ...

    @property
    def sig_z(self) -> float:
        """Z sigma in m."""

    @sig_z.setter
    def sig_z(self, arg: float, /) -> None: ...

    @property
    def sig_pz(self) -> float:
        """pz sigma"""

    @sig_pz.setter
    def sig_pz(self, arg: float, /) -> None: ...

    @property
    def bunch_charge(self) -> float:
        """charge (Coul) in a bunch."""

    @bunch_charge.setter
    def bunch_charge(self, arg: float, /) -> None: ...

    @property
    def n_bunch(self) -> int:
        """Number of bunches."""

    @n_bunch.setter
    def n_bunch(self, arg: int, /) -> None: ...

    @property
    def ix_turn(self) -> int:
        """Turn index used to adjust particles time if needed."""

    @ix_turn.setter
    def ix_turn(self, arg: int, /) -> None: ...

    @property
    def species(self) -> str:
        r"""'positron', etc. \'\' => use referece particle."""

    @species.setter
    def species(self, arg: str, /) -> None: ...

    @property
    def full_6D_coupling_calc(self) -> bool:
        """
        Use V from 6x6 1-turn mat to match distribution? Else use 4x4 1-turn mat used.
        """

    @full_6D_coupling_calc.setter
    def full_6D_coupling_calc(self, arg: bool, /) -> None: ...

    @property
    def use_particle_start(self) -> bool:
        """
        Use lat%particle_start instead of beam_init%center, %t_offset, and %spin?
        """

    @use_particle_start.setter
    def use_particle_start(self, arg: bool, /) -> None: ...

    @property
    def use_t_coords(self) -> bool:
        """If true, the distributions will be taken as in t-coordinates"""

    @use_t_coords.setter
    def use_t_coords(self, arg: bool, /) -> None: ...

    @property
    def use_z_as_t(self) -> bool:
        """
        Only used if  use_t_coords = .true. If true,  z describes the t distribution If false, z describes the s distribution
        """

    @use_z_as_t.setter
    def use_z_as_t(self, arg: bool, /) -> None: ...

    @property
    def file_name(self) -> str:
        """OLD!! DO NOT USE!!"""

    @file_name.setter
    def file_name(self, arg: str, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> BeamInitStruct: ...

    def __deepcopy__(self, arg: dict, /) -> BeamInitStruct: ...

    def __eq__(self, arg: BeamInitStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class BeamStruct:
    """Fortran struct: beam_struct"""

    def __init__(self) -> None: ...

    @property
    def bunch(self) -> BunchStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> BeamStruct: ...

    def __deepcopy__(self, arg: dict, /) -> BeamStruct: ...

    def __eq__(self, arg: BeamStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class BicubicCmplxCoefStruct:
    """Fortran struct: bicubic_cmplx_coef_struct"""

    def __init__(self, coef: Sequence[Sequence[complex]] | None = None, i_box: Sequence[int] | None = None) -> None: ...

    @property
    def coef(self) -> ComplexArray2D:
        """Coefs"""

    @coef.setter
    def coef(self, arg: Sequence[Sequence[complex]], /) -> None: ...

    @property
    def i_box(self) -> IntArray1D:
        """index at lower box corner."""

    @i_box.setter
    def i_box(self, arg: Sequence[int], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> BicubicCmplxCoefStruct: ...

    def __deepcopy__(self, arg: dict, /) -> BicubicCmplxCoefStruct: ...

    def __eq__(self, arg: BicubicCmplxCoefStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class BicubicCoefStruct:
    """Fortran struct: bicubic_coef_struct"""

    def __init__(self, coef: Sequence[Sequence[float]] | None = None, i_box: Sequence[int] | None = None) -> None: ...

    @property
    def coef(self) -> RealArray2D:
        """Coefs"""

    @coef.setter
    def coef(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def i_box(self) -> IntArray1D:
        """index at lower box corner."""

    @i_box.setter
    def i_box(self, arg: Sequence[int], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> BicubicCoefStruct: ...

    def __deepcopy__(self, arg: dict, /) -> BicubicCoefStruct: ...

    def __eq__(self, arg: BicubicCoefStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class BinStruct:
    """Fortran struct: bin_struct"""

    def __init__(self, count: Sequence[float] | None = None, min: float | None = None, max: float | None = None, delta: float | None = None, n: int | None = None) -> None: ...

    @property
    def count(self) -> RealAlloc1D:
        """Counts (or weight) in each bin"""

    @count.setter
    def count(self, arg: Sequence[float], /) -> None: ...

    @property
    def min(self) -> float:
        """Bounds for the bins"""

    @min.setter
    def min(self, arg: float, /) -> None: ...

    @property
    def max(self) -> float: ...

    @max.setter
    def max(self, arg: float, /) -> None: ...

    @property
    def delta(self) -> float:
        """Size of a bin"""

    @delta.setter
    def delta(self, arg: float, /) -> None: ...

    @property
    def n(self) -> int:
        """Number of bins"""

    @n.setter
    def n(self, arg: int, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> BinStruct: ...

    def __deepcopy__(self, arg: dict, /) -> BinStruct: ...

    def __eq__(self, arg: BinStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class BmadCommonStruct:
    """Fortran struct: bmad_common_struct"""

    def __init__(self, max_aperture_limit: float | None = None, d_orb: Sequence[float] | None = None, default_ds_step: float | None = None, significant_length: float | None = None, rel_tol_tracking: float | None = None, abs_tol_tracking: float | None = None, rel_tol_adaptive_tracking: float | None = None, abs_tol_adaptive_tracking: float | None = None, init_ds_adaptive_tracking: float | None = None, min_ds_adaptive_tracking: float | None = None, fatal_ds_adaptive_tracking: float | None = None, autoscale_amp_abs_tol: float | None = None, autoscale_amp_rel_tol: float | None = None, autoscale_phase_tol: float | None = None, electric_dipole_moment: float | None = None, synch_rad_scale: float | None = None, sad_eps_scale: float | None = None, sad_amp_max: float | None = None, sad_n_div_max: int | None = None, taylor_order: int | None = None, runge_kutta_order: int | None = None, default_integ_order: int | None = None, max_num_runge_kutta_step: int | None = None, rf_phase_below_transition_ref: bool | None = None, sr_wakes_on: bool | None = None, lr_wakes_on: bool | None = None, auto_bookkeeper: bool | None = None, high_energy_space_charge_on: bool | None = None, high_energy_space_charge_linear: bool | None = None, csr_and_space_charge_on: bool | None = None, spin_tracking_on: bool | None = None, spin_sokolov_ternov_flipping_on: bool | None = None, radiation_damping_on: bool | None = None, radiation_zero_average: bool | None = None, radiation_fluctuations_on: bool | None = None, conserve_taylor_maps: bool | None = None, absolute_time_tracking: bool | None = None, absolute_time_ref_shift: bool | None = None, convert_to_kinetic_momentum: bool | None = None, normalize_twiss: bool | None = None, aperture_limit_on: bool | None = None, spin_n0_direction_user_set: bool | None = None, debug: bool | None = None) -> None: ...

    @property
    def max_aperture_limit(self) -> float:
        """Max Aperture."""

    @max_aperture_limit.setter
    def max_aperture_limit(self, arg: float, /) -> None: ...

    @property
    def d_orb(self) -> RealArray1D:
        """Orbit deltas for the mat6 via tracking calc."""

    @d_orb.setter
    def d_orb(self, arg: Sequence[float], /) -> None: ...

    @property
    def default_ds_step(self) -> float:
        """Default integration step for eles without an explicit step calc."""

    @default_ds_step.setter
    def default_ds_step(self, arg: float, /) -> None: ...

    @property
    def significant_length(self) -> float:
        """meter"""

    @significant_length.setter
    def significant_length(self, arg: float, /) -> None: ...

    @property
    def rel_tol_tracking(self) -> float:
        """Closed orbit relative tolerance."""

    @rel_tol_tracking.setter
    def rel_tol_tracking(self, arg: float, /) -> None: ...

    @property
    def abs_tol_tracking(self) -> float:
        """Closed orbit absolute tolerance."""

    @abs_tol_tracking.setter
    def abs_tol_tracking(self, arg: float, /) -> None: ...

    @property
    def rel_tol_adaptive_tracking(self) -> float:
        """Runge-Kutta tracking relative tolerance."""

    @rel_tol_adaptive_tracking.setter
    def rel_tol_adaptive_tracking(self, arg: float, /) -> None: ...

    @property
    def abs_tol_adaptive_tracking(self) -> float:
        """Runge-Kutta tracking absolute tolerance."""

    @abs_tol_adaptive_tracking.setter
    def abs_tol_adaptive_tracking(self, arg: float, /) -> None: ...

    @property
    def init_ds_adaptive_tracking(self) -> float:
        """Initial step size"""

    @init_ds_adaptive_tracking.setter
    def init_ds_adaptive_tracking(self, arg: float, /) -> None: ...

    @property
    def min_ds_adaptive_tracking(self) -> float:
        """Min step size to take."""

    @min_ds_adaptive_tracking.setter
    def min_ds_adaptive_tracking(self, arg: float, /) -> None: ...

    @property
    def fatal_ds_adaptive_tracking(self) -> float:
        """If actual step size is below this particle is lost."""

    @fatal_ds_adaptive_tracking.setter
    def fatal_ds_adaptive_tracking(self, arg: float, /) -> None: ...

    @property
    def autoscale_amp_abs_tol(self) -> float:
        """Autoscale absolute amplitude tolerance (eV)."""

    @autoscale_amp_abs_tol.setter
    def autoscale_amp_abs_tol(self, arg: float, /) -> None: ...

    @property
    def autoscale_amp_rel_tol(self) -> float:
        """Autoscale relative amplitude tolerance"""

    @autoscale_amp_rel_tol.setter
    def autoscale_amp_rel_tol(self, arg: float, /) -> None: ...

    @property
    def autoscale_phase_tol(self) -> float:
        """Autoscale phase tolerance."""

    @autoscale_phase_tol.setter
    def autoscale_phase_tol(self, arg: float, /) -> None: ...

    @property
    def electric_dipole_moment(self) -> float:
        """Particle's EDM. Call set_ptc to transfer value to PTC."""

    @electric_dipole_moment.setter
    def electric_dipole_moment(self, arg: float, /) -> None: ...

    @property
    def synch_rad_scale(self) -> float:
        """Synch radiation kick scale. 1 => normal, 0 => no kicks."""

    @synch_rad_scale.setter
    def synch_rad_scale(self, arg: float, /) -> None: ...

    @property
    def sad_eps_scale(self) -> float:
        """Used in sad_mult step length calc."""

    @sad_eps_scale.setter
    def sad_eps_scale(self, arg: float, /) -> None: ...

    @property
    def sad_amp_max(self) -> float:
        """Used in sad_mult step length calc."""

    @sad_amp_max.setter
    def sad_amp_max(self, arg: float, /) -> None: ...

    @property
    def sad_n_div_max(self) -> int:
        """Used in sad_mult step length calc."""

    @sad_n_div_max.setter
    def sad_n_div_max(self, arg: int, /) -> None: ...

    @property
    def taylor_order(self) -> int:
        """Taylor order to use. 0 -> default = ptc_private%taylor_order_saved."""

    @taylor_order.setter
    def taylor_order(self, arg: int, /) -> None: ...

    @property
    def runge_kutta_order(self) -> int:
        """Runge Kutta order."""

    @runge_kutta_order.setter
    def runge_kutta_order(self, arg: int, /) -> None: ...

    @property
    def default_integ_order(self) -> int:
        """PTC integration order."""

    @default_integ_order.setter
    def default_integ_order(self, arg: int, /) -> None: ...

    @property
    def max_num_runge_kutta_step(self) -> int:
        """Maximum number of RK steps before particle is considered lost."""

    @max_num_runge_kutta_step.setter
    def max_num_runge_kutta_step(self, arg: int, /) -> None: ...

    @property
    def rf_phase_below_transition_ref(self) -> bool:
        """Autoscale uses below transition stable point for RFCavities?"""

    @rf_phase_below_transition_ref.setter
    def rf_phase_below_transition_ref(self, arg: bool, /) -> None: ...

    @property
    def sr_wakes_on(self) -> bool:
        """Short range wakefields?"""

    @sr_wakes_on.setter
    def sr_wakes_on(self, arg: bool, /) -> None: ...

    @property
    def lr_wakes_on(self) -> bool:
        """Long range wakefields"""

    @lr_wakes_on.setter
    def lr_wakes_on(self, arg: bool, /) -> None: ...

    @property
    def auto_bookkeeper(self) -> bool:
        """Deprecated and no longer used."""

    @auto_bookkeeper.setter
    def auto_bookkeeper(self, arg: bool, /) -> None: ...

    @property
    def high_energy_space_charge_on(self) -> bool:
        """High energy space charge effect switch."""

    @high_energy_space_charge_on.setter
    def high_energy_space_charge_on(self, arg: bool, /) -> None: ...

    @property
    def high_energy_space_charge_linear(self) -> bool:
        """High energy space charge effect switch."""

    @high_energy_space_charge_linear.setter
    def high_energy_space_charge_linear(self, arg: bool, /) -> None: ...

    @property
    def csr_and_space_charge_on(self) -> bool:
        """Space charge switch."""

    @csr_and_space_charge_on.setter
    def csr_and_space_charge_on(self, arg: bool, /) -> None: ...

    @property
    def spin_tracking_on(self) -> bool:
        """spin tracking?"""

    @spin_tracking_on.setter
    def spin_tracking_on(self, arg: bool, /) -> None: ...

    @property
    def spin_sokolov_ternov_flipping_on(self) -> bool:
        """Spin flipping during synchrotron radiation emission?"""

    @spin_sokolov_ternov_flipping_on.setter
    def spin_sokolov_ternov_flipping_on(self, arg: bool, /) -> None: ...

    @property
    def radiation_damping_on(self) -> bool:
        """Radiation damping toggle."""

    @radiation_damping_on.setter
    def radiation_damping_on(self, arg: bool, /) -> None: ...

    @property
    def radiation_zero_average(self) -> bool:
        """Shift damping to be zero on the zero orbit to get rid of sawtooth?"""

    @radiation_zero_average.setter
    def radiation_zero_average(self, arg: bool, /) -> None: ...

    @property
    def radiation_fluctuations_on(self) -> bool:
        """Radiation fluctuations toggle."""

    @radiation_fluctuations_on.setter
    def radiation_fluctuations_on(self, arg: bool, /) -> None: ...

    @property
    def conserve_taylor_maps(self) -> bool:
        """Enable bookkeeper to set ele%taylor_map_includes_offsets = F?"""

    @conserve_taylor_maps.setter
    def conserve_taylor_maps(self, arg: bool, /) -> None: ...

    @property
    def absolute_time_tracking(self) -> bool:
        """Absolute or relative time tracking?"""

    @absolute_time_tracking.setter
    def absolute_time_tracking(self, arg: bool, /) -> None: ...

    @property
    def absolute_time_ref_shift(self) -> bool:
        """Apply reference time shift when using absolute time tracking?"""

    @absolute_time_ref_shift.setter
    def absolute_time_ref_shift(self, arg: bool, /) -> None: ...

    @property
    def convert_to_kinetic_momentum(self) -> bool:
        """
        Cancel kicks due to finite vector potential when doing symplectic tracking? Set to True to test symp_lie_bmad against runge_kutta.
        """

    @convert_to_kinetic_momentum.setter
    def convert_to_kinetic_momentum(self, arg: bool, /) -> None: ...

    @property
    def normalize_twiss(self) -> bool:
        """Normalize matrix when computing Twiss for off-energy ref?"""

    @normalize_twiss.setter
    def normalize_twiss(self, arg: bool, /) -> None: ...

    @property
    def aperture_limit_on(self) -> bool:
        """Use apertures in tracking?"""

    @aperture_limit_on.setter
    def aperture_limit_on(self, arg: bool, /) -> None: ...

    @property
    def spin_n0_direction_user_set(self) -> bool:
        """User sets direction of n0 for closed geometry branches?"""

    @spin_n0_direction_user_set.setter
    def spin_n0_direction_user_set(self, arg: bool, /) -> None: ...

    @property
    def debug(self) -> bool:
        """Used for code debugging."""

    @debug.setter
    def debug(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> BmadCommonStruct: ...

    def __deepcopy__(self, arg: dict, /) -> BmadCommonStruct: ...

    def __eq__(self, arg: BmadCommonStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class BmadNormalFormStruct:
    """Fortran struct: bmad_normal_form_struct"""

    def __init__(self, ele_origin: EleStruct | None = None) -> None: ...

    @property
    def ele_origin(self) -> EleStruct | None:
        """Element at which the on-turn map was created."""

    @ele_origin.setter
    def ele_origin(self, arg: EleStruct, /) -> None: ...

    @property
    def M(self) -> TaylorStructArray1D:
        """One-turn taylor map: M = A o N o A_inv, N = exp(:h:)"""

    @property
    def A(self) -> TaylorStructArray1D:
        """Map from Floquet -> Lab coordinates"""

    @property
    def A_inv(self) -> TaylorStructArray1D:
        """Map from Lab -> Floquet coordinates"""

    @property
    def dhdj(self) -> TaylorStructArray1D:
        """Nonlinear tune function operating on Floquet coordinates"""

    @property
    def F(self) -> ComplexTaylorStructArray1D:
        """Vector field factorization in phasor basis:"""

    @property
    def L(self) -> ComplexTaylorStructArray1D:
        """L component"""

    @property
    def h(self) -> ResonanceHStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> BmadNormalFormStruct: ...

    def __deepcopy__(self, arg: dict, /) -> BmadNormalFormStruct: ...

    def __eq__(self, arg: BmadNormalFormStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class BookkeepingStateStruct:
    """Fortran struct: bookkeeping_state_struct"""

    def __init__(self, attributes: int | None = None, control: int | None = None, floor_position: int | None = None, s_position: int | None = None, ref_energy: int | None = None, mat6: int | None = None, rad_int: int | None = None, ptc: int | None = None, has_misalign: bool | None = None) -> None: ...

    @property
    def attributes(self) -> int:
        """Element dependent attributes: super_ok$, ok$ or stale$"""

    @attributes.setter
    def attributes(self, arg: int, /) -> None: ...

    @property
    def control(self) -> int:
        """Lord/slave bookkeeping status: super_ok$, ok$ or stale$"""

    @control.setter
    def control(self, arg: int, /) -> None: ...

    @property
    def floor_position(self) -> int:
        """Global (floor) geometry: super_ok$, ok$ or stale$"""

    @floor_position.setter
    def floor_position(self, arg: int, /) -> None: ...

    @property
    def s_position(self) -> int:
        """Longitudinal position & element length: super_ok$, ok$ or stale$"""

    @s_position.setter
    def s_position(self, arg: int, /) -> None: ...

    @property
    def ref_energy(self) -> int:
        """Reference energy and ref time: super_ok$, ok$ or stale$"""

    @ref_energy.setter
    def ref_energy(self, arg: int, /) -> None: ...

    @property
    def mat6(self) -> int:
        """Linear transfer map status: super_ok$, ok$ or stale$"""

    @mat6.setter
    def mat6(self, arg: int, /) -> None: ...

    @property
    def rad_int(self) -> int:
        """Radiation integrals cache status"""

    @rad_int.setter
    def rad_int(self, arg: int, /) -> None: ...

    @property
    def ptc(self) -> int:
        """Associated PTC fibre (or layout) status."""

    @ptc.setter
    def ptc(self, arg: int, /) -> None: ...

    @property
    def has_misalign(self) -> bool:
        """Used to avoid unnecessary calls to offset_particle."""

    @has_misalign.setter
    def has_misalign(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> BookkeepingStateStruct: ...

    def __deepcopy__(self, arg: dict, /) -> BookkeepingStateStruct: ...

    def __eq__(self, arg: BookkeepingStateStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class BpmPhaseCouplingStruct:
    """Fortran struct: bpm_phase_coupling_struct"""

    def __init__(self, K_22a: float | None = None, K_12a: float | None = None, K_11b: float | None = None, K_12b: float | None = None, Cbar22_a: float | None = None, Cbar12_a: float | None = None, Cbar11_b: float | None = None, Cbar12_b: float | None = None, phi_a: float | None = None, phi_b: float | None = None) -> None: ...

    @property
    def K_22a(self) -> float:
        """In-phase y/x for a-mode oscillations."""

    @K_22a.setter
    def K_22a(self, arg: float, /) -> None: ...

    @property
    def K_12a(self) -> float:
        """Out-of-phase y/x for a-mode oscillations."""

    @K_12a.setter
    def K_12a(self, arg: float, /) -> None: ...

    @property
    def K_11b(self) -> float:
        """In-phase x/y for b-mode oscillations."""

    @K_11b.setter
    def K_11b(self, arg: float, /) -> None: ...

    @property
    def K_12b(self) -> float:
        """Out-of-phase x/y for b-mode oscillations."""

    @K_12b.setter
    def K_12b(self, arg: float, /) -> None: ...

    @property
    def Cbar22_a(self) -> float:
        """Cbar22 as calculated from K_22a."""

    @Cbar22_a.setter
    def Cbar22_a(self, arg: float, /) -> None: ...

    @property
    def Cbar12_a(self) -> float:
        """Cbar12 as calculated from K_12a."""

    @Cbar12_a.setter
    def Cbar12_a(self, arg: float, /) -> None: ...

    @property
    def Cbar11_b(self) -> float:
        """Cbar11 as calculated from K_11b."""

    @Cbar11_b.setter
    def Cbar11_b(self, arg: float, /) -> None: ...

    @property
    def Cbar12_b(self) -> float:
        """Cbar12 as calculated from K_12b."""

    @Cbar12_b.setter
    def Cbar12_b(self, arg: float, /) -> None: ...

    @property
    def phi_a(self) -> float:
        """a-mode betatron phase."""

    @phi_a.setter
    def phi_a(self, arg: float, /) -> None: ...

    @property
    def phi_b(self) -> float:
        """b-mode betatron phase."""

    @phi_b.setter
    def phi_b(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> BpmPhaseCouplingStruct: ...

    def __deepcopy__(self, arg: dict, /) -> BpmPhaseCouplingStruct: ...

    def __eq__(self, arg: BpmPhaseCouplingStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class BranchPointerStruct:
    """Fortran struct: branch_pointer_struct"""

    def __init__(self, branch: BranchStruct | None = None) -> None: ...

    @property
    def branch(self) -> BranchStruct | None: ...

    @branch.setter
    def branch(self, arg: BranchStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> BranchPointerStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> BranchPointerStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> BranchPointerStruct: ...

    def __deepcopy__(self, arg: dict, /) -> BranchPointerStruct: ...

    def __eq__(self, arg: BranchPointerStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class BranchStruct:
    """Fortran struct: branch_struct"""

    def __init__(self, name: str | None = None, ix_branch: int | None = None, ix_from_branch: int | None = None, ix_from_ele: int | None = None, ix_to_ele: int | None = None, ix_fixer: int | None = None, n_ele_track: int | None = None, n_ele_max: int | None = None, lat: LatStruct | None = None, a: ModeInfoStruct | None = None, b: ModeInfoStruct | None = None, z: ModeInfoStruct | None = None, param: LatParamStruct | None = None, particle_start: CoordStruct | None = None, ptc: PtcBranch1Struct | None = None) -> None: ...

    @property
    def name(self) -> str:
        """Name of line that defines the branch."""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def ix_branch(self) -> int:
        """Index of this branch. 0 => Main branch"""

    @ix_branch.setter
    def ix_branch(self, arg: int, /) -> None: ...

    @property
    def ix_from_branch(self) -> int:
        """-1 => No creating fork element to this branch."""

    @ix_from_branch.setter
    def ix_from_branch(self, arg: int, /) -> None: ...

    @property
    def ix_from_ele(self) -> int:
        """Index of creating fork element which forks to this branch."""

    @ix_from_ele.setter
    def ix_from_ele(self, arg: int, /) -> None: ...

    @property
    def ix_to_ele(self) -> int:
        """Index of element in this branch that creating fork element forks to."""

    @ix_to_ele.setter
    def ix_to_ele(self, arg: int, /) -> None: ...

    @property
    def ix_fixer(self) -> int:
        """Index of active fixer or beginning_ele element."""

    @ix_fixer.setter
    def ix_fixer(self, arg: int, /) -> None: ...

    @property
    def n_ele_track(self) -> int: ...

    @n_ele_track.setter
    def n_ele_track(self, arg: int, /) -> None: ...

    @property
    def n_ele_max(self) -> int: ...

    @n_ele_max.setter
    def n_ele_max(self, arg: int, /) -> None: ...

    @property
    def lat(self) -> LatStruct | None: ...

    @lat.setter
    def lat(self, arg: LatStruct, /) -> None: ...

    @property
    def a(self) -> ModeInfoStruct:
        """Note: Tunes are the fractional part."""

    @a.setter
    def a(self, arg: ModeInfoStruct, /) -> None: ...

    @property
    def b(self) -> ModeInfoStruct:
        """Note: Tunes are the fractional part."""

    @b.setter
    def b(self, arg: ModeInfoStruct, /) -> None: ...

    @property
    def z(self) -> ModeInfoStruct:
        """Note: Tunes are the fractional part."""

    @z.setter
    def z(self, arg: ModeInfoStruct, /) -> None: ...

    @property
    def ele(self) -> EleStructArray1D: ...

    @property
    def param(self) -> LatParamStruct: ...

    @param.setter
    def param(self, arg: LatParamStruct, /) -> None: ...

    @property
    def particle_start(self) -> CoordStruct: ...

    @particle_start.setter
    def particle_start(self, arg: CoordStruct, /) -> None: ...

    @property
    def wall3d(self) -> Wall3DStructArray1D: ...

    @property
    def ptc(self) -> PtcBranch1Struct:
        """
        Pointer to layout. Note: ptc info not transferred with 'branch1 = branch2' set.
        """

    @ptc.setter
    def ptc(self, arg: PtcBranch1Struct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> BranchStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> BranchStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> BranchStruct: ...

    def __deepcopy__(self, arg: dict, /) -> BranchStruct: ...

    def __eq__(self, arg: BranchStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class BunchParamsStruct:
    """Fortran struct: bunch_params_struct"""

    def __init__(self, centroid: CoordStruct | None = None, x: TwissStruct | None = None, y: TwissStruct | None = None, z: TwissStruct | None = None, a: TwissStruct | None = None, b: TwissStruct | None = None, c: TwissStruct | None = None, sigma: Sequence[Sequence[float]] | None = None, rel_max: Sequence[float] | None = None, rel_min: Sequence[float] | None = None, s: float | None = None, t: float | None = None, sigma_t: float | None = None, charge_live: float | None = None, charge_tot: float | None = None, n_particle_tot: int | None = None, n_particle_live: int | None = None, n_particle_lost_in_ele: int | None = None, n_good_steps: int | None = None, n_bad_steps: int | None = None, ix_ele: int | None = None, location: int | None = None, twiss_valid: bool | None = None) -> None: ...

    @property
    def centroid(self) -> CoordStruct:
        """Lab frame"""

    @centroid.setter
    def centroid(self, arg: CoordStruct, /) -> None: ...

    @property
    def x(self) -> TwissStruct:
        """Projected Twiss parameters"""

    @x.setter
    def x(self, arg: TwissStruct, /) -> None: ...

    @property
    def y(self) -> TwissStruct:
        """Projected Twiss parameters"""

    @y.setter
    def y(self, arg: TwissStruct, /) -> None: ...

    @property
    def z(self) -> TwissStruct:
        """Projected Twiss parameters"""

    @z.setter
    def z(self, arg: TwissStruct, /) -> None: ...

    @property
    def a(self) -> TwissStruct:
        """Normal mode twiss parameters"""

    @a.setter
    def a(self, arg: TwissStruct, /) -> None: ...

    @property
    def b(self) -> TwissStruct:
        """Normal mode twiss parameters"""

    @b.setter
    def b(self, arg: TwissStruct, /) -> None: ...

    @property
    def c(self) -> TwissStruct:
        """Normal mode twiss parameters"""

    @c.setter
    def c(self, arg: TwissStruct, /) -> None: ...

    @property
    def sigma(self) -> RealArray2D:
        """beam size matrix"""

    @sigma.setter
    def sigma(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def rel_max(self) -> RealArray1D:
        """Max orbit relative to centroid. 7 -> time."""

    @rel_max.setter
    def rel_max(self, arg: Sequence[float], /) -> None: ...

    @property
    def rel_min(self) -> RealArray1D:
        """Min orbit relative to_centroid. 7 -> time."""

    @rel_min.setter
    def rel_min(self, arg: Sequence[float], /) -> None: ...

    @property
    def s(self) -> float:
        """Longitudinal position."""

    @s.setter
    def s(self, arg: float, /) -> None: ...

    @property
    def t(self) -> float:
        """Time."""

    @t.setter
    def t(self, arg: float, /) -> None: ...

    @property
    def sigma_t(self) -> float:
        """RMS of time spread."""

    @sigma_t.setter
    def sigma_t(self, arg: float, /) -> None: ...

    @property
    def charge_live(self) -> float:
        """Charge of all non-lost particle"""

    @charge_live.setter
    def charge_live(self, arg: float, /) -> None: ...

    @property
    def charge_tot(self) -> float:
        """Charge of all particles."""

    @charge_tot.setter
    def charge_tot(self, arg: float, /) -> None: ...

    @property
    def n_particle_tot(self) -> int:
        """Total number of particles"""

    @n_particle_tot.setter
    def n_particle_tot(self, arg: int, /) -> None: ...

    @property
    def n_particle_live(self) -> int:
        """Number of non-lost particles"""

    @n_particle_live.setter
    def n_particle_live(self, arg: int, /) -> None: ...

    @property
    def n_particle_lost_in_ele(self) -> int:
        """Number lost in element (not calculated by Bmad)"""

    @n_particle_lost_in_ele.setter
    def n_particle_lost_in_ele(self, arg: int, /) -> None: ...

    @property
    def n_good_steps(self) -> int:
        """Number of good steps (set when tracking with space charge)"""

    @n_good_steps.setter
    def n_good_steps(self, arg: int, /) -> None: ...

    @property
    def n_bad_steps(self) -> int:
        """Number of bad steps (set when tracking with space charge)"""

    @n_bad_steps.setter
    def n_bad_steps(self, arg: int, /) -> None: ...

    @property
    def ix_ele(self) -> int:
        """Lattice element where params evaluated at."""

    @ix_ele.setter
    def ix_ele(self, arg: int, /) -> None: ...

    @property
    def location(self) -> int:
        """Location in element: upstream_end$, inside$, or downstream_end$"""

    @location.setter
    def location(self, arg: int, /) -> None: ...

    @property
    def twiss_valid(self) -> bool:
        """
        Is the data here valid? Note: IF there is no energy variation (RF off) twiss_valid may be true but in this case the z-twiss will not be valid.
        """

    @twiss_valid.setter
    def twiss_valid(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> BunchParamsStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> BunchParamsStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> BunchParamsStruct: ...

    def __deepcopy__(self, arg: dict, /) -> BunchParamsStruct: ...

    def __eq__(self, arg: BunchParamsStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class BunchStruct:
    """Fortran struct: bunch_struct"""

    def __init__(self, ix_z: Sequence[int] | None = None, charge_tot: float | None = None, charge_live: float | None = None, z_center: float | None = None, t_center: float | None = None, t0: float | None = None, drift_between_t_and_s: bool | None = None, ix_ele: int | None = None, ix_bunch: int | None = None, ix_turn: int | None = None, n_live: int | None = None, n_good: int | None = None, n_bad: int | None = None) -> None: ...

    @property
    def particle(self) -> CoordStructAlloc1D: ...

    @property
    def ix_z(self) -> IntAlloc1D:
        """bunch%ix_z(1) is index of head particle, etc."""

    @ix_z.setter
    def ix_z(self, arg: Sequence[int], /) -> None: ...

    @property
    def charge_tot(self) -> float:
        """Total charge in a bunch (Coul)."""

    @charge_tot.setter
    def charge_tot(self, arg: float, /) -> None: ...

    @property
    def charge_live(self) -> float:
        """Charge of live particles (Coul)."""

    @charge_live.setter
    def charge_live(self, arg: float, /) -> None: ...

    @property
    def z_center(self) -> float:
        """
        Longitudinal center of bunch at creation time. Note: Generally, z_center of bunch #1 is 0 and z_center of the other bunches is negative.
        """

    @z_center.setter
    def z_center(self, arg: float, /) -> None: ...

    @property
    def t_center(self) -> float:
        """Center of bunch at creation time relative to head bunch."""

    @t_center.setter
    def t_center(self, arg: float, /) -> None: ...

    @property
    def t0(self) -> float:
        """
        Used by track1_bunch_space_charge for tracking so particles have constant t.
        """

    @t0.setter
    def t0(self, arg: float, /) -> None: ...

    @property
    def drift_between_t_and_s(self) -> bool:
        """
        Drift (ignore any fields) instead of tracking to speed up the calculation? This can only be done under certain circumstances.
        """

    @drift_between_t_and_s.setter
    def drift_between_t_and_s(self, arg: bool, /) -> None: ...

    @property
    def ix_ele(self) -> int:
        """
        Nominal element bunch is at. But, EG, dead particles can be someplace else.
        """

    @ix_ele.setter
    def ix_ele(self, arg: int, /) -> None: ...

    @property
    def ix_bunch(self) -> int:
        """Bunch index. Head bunch = 1, etc."""

    @ix_bunch.setter
    def ix_bunch(self, arg: int, /) -> None: ...

    @property
    def ix_turn(self) -> int:
        """
        Turn index for long term tracking. ix_turn = 0 before end of first turn, etc.
        """

    @ix_turn.setter
    def ix_turn(self, arg: int, /) -> None: ...

    @property
    def n_live(self) -> int: ...

    @n_live.setter
    def n_live(self, arg: int, /) -> None: ...

    @property
    def n_good(self) -> int:
        """Number of accepted steps when using adaptive step size control."""

    @n_good.setter
    def n_good(self, arg: int, /) -> None: ...

    @property
    def n_bad(self) -> int:
        """Number of rejected steps when using adaptive step size control."""

    @n_bad.setter
    def n_bad(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> BunchStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> BunchStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> BunchStruct: ...

    def __deepcopy__(self, arg: dict, /) -> BunchStruct: ...

    def __eq__(self, arg: BunchStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class BunchTrackStruct:
    """Fortran struct: bunch_track_struct"""

    def __init__(self, ds_save: float | None = None, n_pt: int | None = None) -> None: ...

    @property
    def pt(self) -> BunchParamsStructAlloc1D:
        """Array indexed from 0"""

    @property
    def ds_save(self) -> float:
        """Min distance between points."""

    @ds_save.setter
    def ds_save(self, arg: float, /) -> None: ...

    @property
    def n_pt(self) -> int:
        """Track upper bound"""

    @n_pt.setter
    def n_pt(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> BunchTrackStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> BunchTrackStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> BunchTrackStruct: ...

    def __deepcopy__(self, arg: dict, /) -> BunchTrackStruct: ...

    def __eq__(self, arg: BunchTrackStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class CartesianMapStruct:
    """Fortran struct: cartesian_map_struct"""

    def __init__(self, field_scale: float | None = None, r0: Sequence[float] | None = None, master_parameter: int | None = None, ele_anchor_pt: int | None = None, field_type: int | None = None, ptr: CartesianMapTermStruct | None = None) -> None: ...

    @property
    def field_scale(self) -> float:
        """Factor to scale the fields by"""

    @field_scale.setter
    def field_scale(self, arg: float, /) -> None: ...

    @property
    def r0(self) -> RealArray1D:
        """Field origin offset."""

    @r0.setter
    def r0(self, arg: Sequence[float], /) -> None: ...

    @property
    def master_parameter(self) -> int:
        """Master parameter in ele%value(:) array to use for scaling the field."""

    @master_parameter.setter
    def master_parameter(self, arg: int, /) -> None: ...

    @property
    def ele_anchor_pt(self) -> int:
        """anchor_beginning$, anchor_center$, or anchor_end$"""

    @ele_anchor_pt.setter
    def ele_anchor_pt(self, arg: int, /) -> None: ...

    @property
    def field_type(self) -> int:
        """or electric$"""

    @field_type.setter
    def field_type(self, arg: int, /) -> None: ...

    @property
    def ptr(self) -> CartesianMapTermStruct | None: ...

    @ptr.setter
    def ptr(self, arg: CartesianMapTermStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> CartesianMapStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> CartesianMapStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> CartesianMapStruct: ...

    def __deepcopy__(self, arg: dict, /) -> CartesianMapStruct: ...

    def __eq__(self, arg: CartesianMapStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class CartesianMapTerm1Struct:
    """Fortran struct: cartesian_map_term1_struct"""

    def __init__(self, coef: float | None = None, kx: float | None = None, ky: float | None = None, kz: float | None = None, x0: float | None = None, y0: float | None = None, phi_z: float | None = None, family: int | None = None, form: int | None = None) -> None: ...

    @property
    def coef(self) -> float: ...

    @coef.setter
    def coef(self, arg: float, /) -> None: ...

    @property
    def kx(self) -> float: ...

    @kx.setter
    def kx(self, arg: float, /) -> None: ...

    @property
    def ky(self) -> float: ...

    @ky.setter
    def ky(self, arg: float, /) -> None: ...

    @property
    def kz(self) -> float: ...

    @kz.setter
    def kz(self, arg: float, /) -> None: ...

    @property
    def x0(self) -> float: ...

    @x0.setter
    def x0(self, arg: float, /) -> None: ...

    @property
    def y0(self) -> float: ...

    @y0.setter
    def y0(self, arg: float, /) -> None: ...

    @property
    def phi_z(self) -> float: ...

    @phi_z.setter
    def phi_z(self, arg: float, /) -> None: ...

    @property
    def family(self) -> int:
        """family_x$, etc."""

    @family.setter
    def family(self, arg: int, /) -> None: ...

    @property
    def form(self) -> int:
        """hyper_y$, etc."""

    @form.setter
    def form(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> CartesianMapTerm1StructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> CartesianMapTerm1StructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> CartesianMapTerm1Struct: ...

    def __deepcopy__(self, arg: dict, /) -> CartesianMapTerm1Struct: ...

    def __eq__(self, arg: CartesianMapTerm1Struct, /) -> bool: ...

    def __hash__(self) -> int: ...

class CartesianMapTermStruct:
    """Fortran struct: cartesian_map_term_struct"""

    def __init__(self, file: str | None = None, n_link: int | None = None) -> None: ...

    @property
    def file(self) -> str:
        """Input file name. Used also as ID for instances."""

    @file.setter
    def file(self, arg: str, /) -> None: ...

    @property
    def n_link(self) -> int:
        """For memory management of %term"""

    @n_link.setter
    def n_link(self, arg: int, /) -> None: ...

    @property
    def term(self) -> CartesianMapTerm1StructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> CartesianMapTermStruct: ...

    def __deepcopy__(self, arg: dict, /) -> CartesianMapTermStruct: ...

    def __eq__(self, arg: CartesianMapTermStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class CmplxField1At2DPtStruct:
    """Fortran struct: cmplx_field1_at_2D_pt_struct"""

    def __init__(self, f: complex | None = None, df_dx: complex | None = None, df_dy: complex | None = None, d2f_dxdy: complex | None = None) -> None: ...

    @property
    def f(self) -> complex:
        """Field"""

    @f.setter
    def f(self, arg: complex, /) -> None: ...

    @property
    def df_dx(self) -> complex:
        """Normalized field 1st derivatives"""

    @df_dx.setter
    def df_dx(self, arg: complex, /) -> None: ...

    @property
    def df_dy(self) -> complex:
        """Normalized field 1st derivatives"""

    @df_dy.setter
    def df_dy(self, arg: complex, /) -> None: ...

    @property
    def d2f_dxdy(self) -> complex:
        """Normalized field 2nd derivative"""

    @d2f_dxdy.setter
    def d2f_dxdy(self, arg: complex, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> CmplxField1At2DPtStruct: ...

    def __deepcopy__(self, arg: dict, /) -> CmplxField1At2DPtStruct: ...

    def __eq__(self, arg: CmplxField1At2DPtStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class CmplxField1At3DPtStruct:
    """Fortran struct: cmplx_field1_at_3D_pt_struct"""

    def __init__(self, f: complex | None = None, df_dx: complex | None = None, df_dy: complex | None = None, df_dz: complex | None = None, d2f_dxdy: complex | None = None, d2f_dxdz: complex | None = None, d2f_dydz: complex | None = None, d3f_dxdydz: complex | None = None) -> None: ...

    @property
    def f(self) -> complex:
        """Field"""

    @f.setter
    def f(self, arg: complex, /) -> None: ...

    @property
    def df_dx(self) -> complex:
        """Normalized field 1st derivatives"""

    @df_dx.setter
    def df_dx(self, arg: complex, /) -> None: ...

    @property
    def df_dy(self) -> complex:
        """Normalized field 1st derivatives"""

    @df_dy.setter
    def df_dy(self, arg: complex, /) -> None: ...

    @property
    def df_dz(self) -> complex:
        """Normalized field 1st derivatives"""

    @df_dz.setter
    def df_dz(self, arg: complex, /) -> None: ...

    @property
    def d2f_dxdy(self) -> complex:
        """Normalized field 2nd derivatives"""

    @d2f_dxdy.setter
    def d2f_dxdy(self, arg: complex, /) -> None: ...

    @property
    def d2f_dxdz(self) -> complex:
        """Normalized field 2nd derivatives"""

    @d2f_dxdz.setter
    def d2f_dxdz(self, arg: complex, /) -> None: ...

    @property
    def d2f_dydz(self) -> complex:
        """Normalized field 2nd derivatives"""

    @d2f_dydz.setter
    def d2f_dydz(self, arg: complex, /) -> None: ...

    @property
    def d3f_dxdydz(self) -> complex:
        """Normalized field 3rd derivative"""

    @d3f_dxdydz.setter
    def d3f_dxdydz(self, arg: complex, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> CmplxField1At3DPtStruct: ...

    def __deepcopy__(self, arg: dict, /) -> CmplxField1At3DPtStruct: ...

    def __eq__(self, arg: CmplxField1At3DPtStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class CmplxFieldAt2DBoxStruct:
    """Fortran struct: cmplx_field_at_2D_box_struct"""

    def __init__(self, i_box: Sequence[int] | None = None) -> None: ...

    @property
    def pt(self) -> CmplxField1At2DPtStructArray2D: ...

    @property
    def i_box(self) -> IntArray1D:
        """index at lower box corner."""

    @i_box.setter
    def i_box(self, arg: Sequence[int], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> CmplxFieldAt2DBoxStruct: ...

    def __deepcopy__(self, arg: dict, /) -> CmplxFieldAt2DBoxStruct: ...

    def __eq__(self, arg: CmplxFieldAt2DBoxStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class CmplxFieldAt3DBoxStruct:
    """Fortran struct: cmplx_field_at_3D_box_struct"""

    def __init__(self, i_box: Sequence[int] | None = None) -> None: ...

    @property
    def pt(self) -> CmplxField1At3DPtStructArray3D: ...

    @property
    def i_box(self) -> IntArray1D:
        """index at lower box corner."""

    @i_box.setter
    def i_box(self, arg: Sequence[int], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> CmplxFieldAt3DBoxStruct: ...

    def __deepcopy__(self, arg: dict, /) -> CmplxFieldAt3DBoxStruct: ...

    def __eq__(self, arg: CmplxFieldAt3DBoxStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ComplexTaylorStruct:
    """Fortran struct: complex_taylor_struct"""

    def __init__(self, ref: complex | None = None) -> None: ...

    @property
    def ref(self) -> complex: ...

    @ref.setter
    def ref(self, arg: complex, /) -> None: ...

    @property
    def term(self) -> ComplexTaylorTermStructArray1D: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> ComplexTaylorStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> ComplexTaylorStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ComplexTaylorStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ComplexTaylorStruct: ...

    def __eq__(self, arg: ComplexTaylorStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ComplexTaylorTermStruct:
    """Fortran struct: complex_taylor_term_struct"""

    def __init__(self, coef: complex | None = None, expn: Sequence[int] | None = None) -> None: ...

    @property
    def coef(self) -> complex: ...

    @coef.setter
    def coef(self, arg: complex, /) -> None: ...

    @property
    def expn(self) -> IntArray1D: ...

    @expn.setter
    def expn(self, arg: Sequence[int], /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> ComplexTaylorTermStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> ComplexTaylorTermStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ComplexTaylorTermStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ComplexTaylorTermStruct: ...

    def __eq__(self, arg: ComplexTaylorTermStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ControlRamp1Struct:
    """Fortran struct: control_ramp1_struct"""

    def __init__(self, y_knot: Sequence[float] | None = None, attribute: str | None = None, slave_name: str | None = None, is_controller: bool | None = None) -> None: ...

    @property
    def y_knot(self) -> RealAlloc1D: ...

    @y_knot.setter
    def y_knot(self, arg: Sequence[float], /) -> None: ...

    @property
    def stack(self) -> ExpressionAtomStructAlloc1D:
        """Evaluation stack"""

    @property
    def attribute(self) -> str:
        """
        Name of attribute controlled. Set to 'FIELD_OVERLAPS' for field overlaps.
        """

    @attribute.setter
    def attribute(self, arg: str, /) -> None: ...

    @property
    def slave_name(self) -> str:
        """Name of slave."""

    @slave_name.setter
    def slave_name(self, arg: str, /) -> None: ...

    @property
    def is_controller(self) -> bool:
        """Is the slave a controller? If so bookkeeping is different."""

    @is_controller.setter
    def is_controller(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> ControlRamp1StructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> ControlRamp1StructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ControlRamp1Struct: ...

    def __deepcopy__(self, arg: dict, /) -> ControlRamp1Struct: ...

    def __eq__(self, arg: ControlRamp1Struct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ControlStruct:
    """Fortran struct: control_struct"""

    def __init__(self, value: float | None = None, y_knot: Sequence[float] | None = None, slave: LatEleLocStruct | None = None, lord: LatEleLocStruct | None = None, slave_name: str | None = None, attribute: str | None = None, ix_attrib: int | None = None) -> None: ...

    @property
    def value(self) -> float:
        """Used by group, and overlay elements."""

    @value.setter
    def value(self, arg: float, /) -> None: ...

    @property
    def y_knot(self) -> RealAlloc1D: ...

    @y_knot.setter
    def y_knot(self, arg: Sequence[float], /) -> None: ...

    @property
    def stack(self) -> ExpressionAtomStructAlloc1D:
        """Evaluation stack"""

    @property
    def slave(self) -> LatEleLocStruct: ...

    @slave.setter
    def slave(self, arg: LatEleLocStruct, /) -> None: ...

    @property
    def lord(self) -> LatEleLocStruct: ...

    @lord.setter
    def lord(self, arg: LatEleLocStruct, /) -> None: ...

    @property
    def slave_name(self) -> str:
        """Name of slave."""

    @slave_name.setter
    def slave_name(self, arg: str, /) -> None: ...

    @property
    def attribute(self) -> str:
        """
        Name of attribute controlled. Set to 'FIELD_OVERLAPS' for field overlaps. Set to 'INPUT' or 'OUTPUT' for feedback slaves.
        """

    @attribute.setter
    def attribute(self, arg: str, /) -> None: ...

    @property
    def ix_attrib(self) -> int:
        """Index of attribute controlled. See note above!"""

    @ix_attrib.setter
    def ix_attrib(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> ControlStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> ControlStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ControlStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ControlStruct: ...

    def __eq__(self, arg: ControlStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ControlVar1Struct:
    """Fortran struct: control_var1_struct"""

    def __init__(self, name: str | None = None, value: float | None = None, old_value: float | None = None) -> None: ...

    @property
    def name(self) -> str: ...

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def value(self) -> float: ...

    @value.setter
    def value(self, arg: float, /) -> None: ...

    @property
    def old_value(self) -> float: ...

    @old_value.setter
    def old_value(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> ControlVar1StructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> ControlVar1StructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ControlVar1Struct: ...

    def __deepcopy__(self, arg: dict, /) -> ControlVar1Struct: ...

    def __eq__(self, arg: ControlVar1Struct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ControllerStruct:
    """Fortran struct: controller_struct"""

    def __init__(self, x_knot: Sequence[float] | None = None) -> None: ...

    @property
    def var(self) -> ControlVar1StructAlloc1D: ...

    @property
    def ramp(self) -> ControlRamp1StructAlloc1D:
        """For ramper lord elements"""

    @property
    def ramper_lord(self) -> RamperLordStructAlloc1D:
        """Ramper lord info for this slave"""

    @property
    def x_knot(self) -> RealAlloc1D: ...

    @x_knot.setter
    def x_knot(self, arg: Sequence[float], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ControllerStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ControllerStruct: ...

    def __eq__(self, arg: ControllerStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ConverterDir1DStruct:
    """Fortran struct: converter_dir_1D_struct"""

    def __init__(self, pc_out: float | None = None, poly: Sequence[float] | None = None) -> None: ...

    @property
    def pc_out(self) -> float:
        """pc_out value at fit"""

    @pc_out.setter
    def pc_out(self, arg: float, /) -> None: ...

    @property
    def poly(self) -> RealArray1D:
        """param(r) = Sum: poly(i) * r^i"""

    @poly.setter
    def poly(self, arg: Sequence[float], /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> ConverterDir1DStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> ConverterDir1DStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ConverterDir1DStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ConverterDir1DStruct: ...

    def __eq__(self, arg: ConverterDir1DStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ConverterDir2DStruct:
    """Fortran struct: converter_dir_2D_struct"""

    def __init__(self, k: float | None = None, poly: Sequence[float] | None = None) -> None: ...

    @property
    def k(self) -> float: ...

    @k.setter
    def k(self, arg: float, /) -> None: ...

    @property
    def poly(self) -> RealArray1D: ...

    @poly.setter
    def poly(self, arg: Sequence[float], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ConverterDir2DStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ConverterDir2DStruct: ...

    def __eq__(self, arg: ConverterDir2DStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ConverterDirCoefStruct:
    """Fortran struct: converter_dir_coef_struct"""

    def __init__(self, fit_2d_r: ConverterDir2DStruct | None = None, fit_2d_pc: ConverterDir2DStruct | None = None, c0: float | None = None) -> None: ...

    @property
    def fit_1d_r(self) -> ConverterDir1DStructAlloc1D: ...

    @property
    def fit_2d_r(self) -> ConverterDir2DStruct: ...

    @fit_2d_r.setter
    def fit_2d_r(self, arg: ConverterDir2DStruct, /) -> None: ...

    @property
    def fit_2d_pc(self) -> ConverterDir2DStruct: ...

    @fit_2d_pc.setter
    def fit_2d_pc(self, arg: ConverterDir2DStruct, /) -> None: ...

    @property
    def c0(self) -> float: ...

    @c0.setter
    def c0(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ConverterDirCoefStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ConverterDirCoefStruct: ...

    def __eq__(self, arg: ConverterDirCoefStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ConverterDirectionOutStruct:
    """Fortran struct: converter_direction_out_struct"""

    def __init__(self, beta: ConverterDirCoefStruct | None = None, alpha_x: ConverterDirCoefStruct | None = None, alpha_y: ConverterDirCoefStruct | None = None, dxds_min: ConverterDirCoefStruct | None = None, dxds_max: ConverterDirCoefStruct | None = None, dyds_max: ConverterDirCoefStruct | None = None, c_x: ConverterDirCoefStruct | None = None) -> None: ...

    @property
    def beta(self) -> ConverterDirCoefStruct: ...

    @beta.setter
    def beta(self, arg: ConverterDirCoefStruct, /) -> None: ...

    @property
    def alpha_x(self) -> ConverterDirCoefStruct: ...

    @alpha_x.setter
    def alpha_x(self, arg: ConverterDirCoefStruct, /) -> None: ...

    @property
    def alpha_y(self) -> ConverterDirCoefStruct: ...

    @alpha_y.setter
    def alpha_y(self, arg: ConverterDirCoefStruct, /) -> None: ...

    @property
    def dxds_min(self) -> ConverterDirCoefStruct: ...

    @dxds_min.setter
    def dxds_min(self, arg: ConverterDirCoefStruct, /) -> None: ...

    @property
    def dxds_max(self) -> ConverterDirCoefStruct: ...

    @dxds_max.setter
    def dxds_max(self, arg: ConverterDirCoefStruct, /) -> None: ...

    @property
    def dyds_max(self) -> ConverterDirCoefStruct: ...

    @dyds_max.setter
    def dyds_max(self, arg: ConverterDirCoefStruct, /) -> None: ...

    @property
    def c_x(self) -> ConverterDirCoefStruct: ...

    @c_x.setter
    def c_x(self, arg: ConverterDirCoefStruct, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ConverterDirectionOutStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ConverterDirectionOutStruct: ...

    def __eq__(self, arg: ConverterDirectionOutStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ConverterDistributionStruct:
    """Fortran struct: converter_distribution_struct"""

    def __init__(self, thickness: float | None = None) -> None: ...

    @property
    def thickness(self) -> float: ...

    @thickness.setter
    def thickness(self, arg: float, /) -> None: ...

    @property
    def sub_dist(self) -> ConverterSubDistributionStructAlloc1D:
        """Distribution at various pc_in values."""

    @staticmethod
    def new_array1d(sz: int = 0) -> ConverterDistributionStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> ConverterDistributionStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ConverterDistributionStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ConverterDistributionStruct: ...

    def __eq__(self, arg: ConverterDistributionStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ConverterProbPcRStruct:
    """Fortran struct: converter_prob_pc_r_struct"""

    def __init__(self, pc_out: Sequence[float] | None = None, r: Sequence[float] | None = None, prob: Sequence[Sequence[float]] | None = None, spin_z: Sequence[Sequence[float]] | None = None, pc_out_min: float | None = None, pc_out_max: float | None = None, integrated_prob: float | None = None, p_norm: Sequence[Sequence[float]] | None = None, integ_pc_out: Sequence[float] | None = None, integ_r: Sequence[Sequence[float]] | None = None, integ_r_ave: Sequence[float] | None = None) -> None: ...

    @property
    def pc_out(self) -> RealAlloc1D:
        """Grid pc_out values."""

    @pc_out.setter
    def pc_out(self, arg: Sequence[float], /) -> None: ...

    @property
    def r(self) -> RealAlloc1D:
        """Grid r_out values."""

    @r.setter
    def r(self, arg: Sequence[float], /) -> None: ...

    @property
    def prob(self) -> RealArray2D:
        """Probability grid."""

    @prob.setter
    def prob(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def spin_z(self) -> RealArray2D:
        """
        Z polarization grid. Stuff below is calculated rather than read in from the lattice file.
        """

    @spin_z.setter
    def spin_z(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def pc_out_min(self) -> float: ...

    @pc_out_min.setter
    def pc_out_min(self, arg: float, /) -> None: ...

    @property
    def pc_out_max(self) -> float: ...

    @pc_out_max.setter
    def pc_out_max(self, arg: float, /) -> None: ...

    @property
    def integrated_prob(self) -> float:
        """Integrated probability over (pc_out, r) with restrictions factered in."""

    @integrated_prob.setter
    def integrated_prob(self, arg: float, /) -> None: ...

    @property
    def p_norm(self) -> RealArray2D:
        """
        Normalized probability taking into account. angle_out_max, pc_out_min, and pc_out_max restrictions.
        """

    @p_norm.setter
    def p_norm(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def integ_pc_out(self) -> RealAlloc1D:
        """Normalized probability integrated from min pc_out up."""

    @integ_pc_out.setter
    def integ_pc_out(self, arg: Sequence[float], /) -> None: ...

    @property
    def integ_r(self) -> RealArray2D: ...

    @integ_r.setter
    def integ_r(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def integ_r_ave(self) -> RealAlloc1D: ...

    @integ_r_ave.setter
    def integ_r_ave(self, arg: Sequence[float], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ConverterProbPcRStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ConverterProbPcRStruct: ...

    def __eq__(self, arg: ConverterProbPcRStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ConverterStruct:
    """Fortran struct: converter_struct"""

    def __init__(self, species_out: int | None = None, material_type: str | None = None) -> None: ...

    @property
    def species_out(self) -> int:
        """Output species"""

    @species_out.setter
    def species_out(self, arg: int, /) -> None: ...

    @property
    def material_type(self) -> str: ...

    @material_type.setter
    def material_type(self, arg: str, /) -> None: ...

    @property
    def dist(self) -> ConverterDistributionStructAlloc1D:
        """Distribution at various thicknesses"""

    def __repr__(self) -> str: ...

    def __copy__(self) -> ConverterStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ConverterStruct: ...

    def __eq__(self, arg: ConverterStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ConverterSubDistributionStruct:
    """Fortran struct: converter_sub_distribution_struct"""

    def __init__(self, pc_in: float | None = None, spin_in: Sequence[float] | None = None, prob_pc_r: ConverterProbPcRStruct | None = None, dir_out: ConverterDirectionOutStruct | None = None) -> None: ...

    @property
    def pc_in(self) -> float: ...

    @pc_in.setter
    def pc_in(self, arg: float, /) -> None: ...

    @property
    def spin_in(self) -> RealArray1D: ...

    @spin_in.setter
    def spin_in(self, arg: Sequence[float], /) -> None: ...

    @property
    def prob_pc_r(self) -> ConverterProbPcRStruct: ...

    @prob_pc_r.setter
    def prob_pc_r(self, arg: ConverterProbPcRStruct, /) -> None: ...

    @property
    def dir_out(self) -> ConverterDirectionOutStruct: ...

    @dir_out.setter
    def dir_out(self, arg: ConverterDirectionOutStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> ConverterSubDistributionStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> ConverterSubDistributionStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ConverterSubDistributionStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ConverterSubDistributionStruct: ...

    def __eq__(self, arg: ConverterSubDistributionStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class CoordArrayStruct:
    """Fortran struct: coord_array_struct"""

    def __init__(self) -> None: ...

    @property
    def orbit(self) -> CoordStructAlloc1D: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> CoordArrayStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> CoordArrayStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> CoordArrayStruct: ...

    def __deepcopy__(self, arg: dict, /) -> CoordArrayStruct: ...

    def __eq__(self, arg: CoordArrayStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class CoordStruct:
    """Fortran struct: coord_struct"""

    def __init__(self, vec: Sequence[float] | None = None, s: float | None = None, t: float | None = None, spin: Sequence[float] | None = None, field: Sequence[float] | None = None, phase: Sequence[float] | None = None, charge: float | None = None, dt_ref: float | None = None, r: float | None = None, p0c: float | None = None, E_potential: float | None = None, beta: float | None = None, ix_ele: int | None = None, ix_branch: int | None = None, ix_turn: int | None = None, ix_user: int | None = None, state: int | None = None, direction: int | None = None, time_dir: int | None = None, species: int | None = None, location: int | None = None) -> None: ...

    @property
    def vec(self) -> RealArray1D:
        """
        (x, px, y, py, z, pz). Generally phase space for charged particles. See Bmad manual.
        """

    @vec.setter
    def vec(self, arg: Sequence[float], /) -> None: ...

    @property
    def s(self) -> float:
        """Longitudinal position"""

    @s.setter
    def s(self, arg: float, /) -> None: ...

    @property
    def t(self) -> float:
        """Absolute time (not relative to reference). Note: Quad precision!"""

    @t.setter
    def t(self, arg: float, /) -> None: ...

    @property
    def spin(self) -> RealArray1D:
        """Spin."""

    @spin.setter
    def spin(self, arg: Sequence[float], /) -> None: ...

    @property
    def field(self) -> RealArray1D:
        """Photon E-field intensity (x,y)."""

    @field.setter
    def field(self, arg: Sequence[float], /) -> None: ...

    @property
    def phase(self) -> RealArray1D:
        """
        Photon E-field phase (x,y). For charged particles, phase(1) is RF phase.
        """

    @phase.setter
    def phase(self, arg: Sequence[float], /) -> None: ...

    @property
    def charge(self) -> float:
        """
        Macroparticle weight (which is different from particle species charge). For some space charge calcs the weight is in Coulombs.
        """

    @charge.setter
    def charge(self, arg: float, /) -> None: ...

    @property
    def dt_ref(self) -> float:
        """
        Used in: * time tracking for computing z. * by coherent photons = path_length/c_light.
        """

    @dt_ref.setter
    def dt_ref(self, arg: float, /) -> None: ...

    @property
    def r(self) -> float:
        """For general use. Not used by Bmad."""

    @r.setter
    def r(self, arg: float, /) -> None: ...

    @property
    def p0c(self) -> float:
        """
        For non-photons: Reference momentum. For photons: Photon momentum (not reference).
        """

    @p0c.setter
    def p0c(self, arg: float, /) -> None: ...

    @property
    def E_potential(self) -> float:
        """Potential energy."""

    @E_potential.setter
    def E_potential(self, arg: float, /) -> None: ...

    @property
    def beta(self) -> float:
        """Velocity / c_light."""

    @beta.setter
    def beta(self, arg: float, /) -> None: ...

    @property
    def ix_ele(self) -> int:
        """
        Index of the lattice element the particle is in. May be -1 if element is not associated with a lattice.
        """

    @ix_ele.setter
    def ix_ele(self, arg: int, /) -> None: ...

    @property
    def ix_branch(self) -> int:
        """Index of the lattice branch the particle is in."""

    @ix_branch.setter
    def ix_branch(self, arg: int, /) -> None: ...

    @property
    def ix_turn(self) -> int:
        """Turn index for multiturn tracking."""

    @ix_turn.setter
    def ix_turn(self, arg: int, /) -> None: ...

    @property
    def ix_user(self) -> int:
        """For general use, not used by Bmad."""

    @ix_user.setter
    def ix_user(self, arg: int, /) -> None: ...

    @property
    def state(self) -> int:
        """alive$, lost$, lost_neg_x_aperture$, lost_pz$, etc."""

    @state.setter
    def state(self, arg: int, /) -> None: ...

    @property
    def direction(self) -> int:
        """
        +1 or -1. Sign of longitudinal direction of motion (ds/dt). This is independent of the element orientation.
        """

    @direction.setter
    def direction(self, arg: int, /) -> None: ...

    @property
    def time_dir(self) -> int:
        """+1 or -1. Time direction. -1 => Traveling backwards in time."""

    @time_dir.setter
    def time_dir(self, arg: int, /) -> None: ...

    @property
    def species(self) -> int:
        """positron$, proton$, etc."""

    @species.setter
    def species(self, arg: int, /) -> None: ...

    @property
    def location(self) -> int:
        """upstream_end$, inside$, or downstream_end$"""

    @location.setter
    def location(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> CoordStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> CoordStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> CoordStruct: ...

    def __deepcopy__(self, arg: dict, /) -> CoordStruct: ...

    def __eq__(self, arg: CoordStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class CrystalParamStruct:
    """Fortran struct: crystal_param_struct"""

    def __init__(self, cap_gamma: float | None = None, dtheta_sin_2theta: float | None = None, b_eff: float | None = None, wavelength: float | None = None, old_vvec: Sequence[float] | None = None, new_vvec: Sequence[float] | None = None) -> None: ...

    @property
    def cap_gamma(self) -> float: ...

    @cap_gamma.setter
    def cap_gamma(self, arg: float, /) -> None: ...

    @property
    def dtheta_sin_2theta(self) -> float: ...

    @dtheta_sin_2theta.setter
    def dtheta_sin_2theta(self, arg: float, /) -> None: ...

    @property
    def b_eff(self) -> float: ...

    @b_eff.setter
    def b_eff(self, arg: float, /) -> None: ...

    @property
    def wavelength(self) -> float: ...

    @wavelength.setter
    def wavelength(self, arg: float, /) -> None: ...

    @property
    def old_vvec(self) -> RealArray1D: ...

    @old_vvec.setter
    def old_vvec(self, arg: Sequence[float], /) -> None: ...

    @property
    def new_vvec(self) -> RealArray1D: ...

    @new_vvec.setter
    def new_vvec(self, arg: Sequence[float], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> CrystalParamStruct: ...

    def __deepcopy__(self, arg: dict, /) -> CrystalParamStruct: ...

    def __eq__(self, arg: CrystalParamStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class CsrBunchSliceStruct:
    """Fortran struct: csr_bunch_slice_struct"""

    def __init__(self, x0: float | None = None, y0: float | None = None, z0_edge: float | None = None, z1_edge: float | None = None, z_center: float | None = None, sig_x: float | None = None, sig_y: float | None = None, charge: float | None = None, dcharge_density_dz: float | None = None, edge_dcharge_density_dz: float | None = None, kick_csr: float | None = None, coef_lsc_plus: Sequence[Sequence[float]] | None = None, coef_lsc_minus: Sequence[Sequence[float]] | None = None, kick_lsc: float | None = None, n_particle: float | None = None) -> None: ...

    @property
    def x0(self) -> float:
        """Transverse center of the particle distrubution"""

    @x0.setter
    def x0(self, arg: float, /) -> None: ...

    @property
    def y0(self) -> float:
        """Transverse center of the particle distrubution"""

    @y0.setter
    def y0(self, arg: float, /) -> None: ...

    @property
    def z0_edge(self) -> float:
        """Left (min z) edge of bin"""

    @z0_edge.setter
    def z0_edge(self, arg: float, /) -> None: ...

    @property
    def z1_edge(self) -> float:
        """Right (max z) edge of bin"""

    @z1_edge.setter
    def z1_edge(self, arg: float, /) -> None: ...

    @property
    def z_center(self) -> float:
        """z at center of bin."""

    @z_center.setter
    def z_center(self, arg: float, /) -> None: ...

    @property
    def sig_x(self) -> float:
        """particle's RMS width"""

    @sig_x.setter
    def sig_x(self, arg: float, /) -> None: ...

    @property
    def sig_y(self) -> float:
        """particle's RMS width"""

    @sig_y.setter
    def sig_y(self, arg: float, /) -> None: ...

    @property
    def charge(self) -> float:
        """charge of the particles"""

    @charge.setter
    def charge(self, arg: float, /) -> None: ...

    @property
    def dcharge_density_dz(self) -> float:
        """Charge density gradient"""

    @dcharge_density_dz.setter
    def dcharge_density_dz(self, arg: float, /) -> None: ...

    @property
    def edge_dcharge_density_dz(self) -> float:
        """gradient between this and preceeding bin. [Evaluated at bin edge.]"""

    @edge_dcharge_density_dz.setter
    def edge_dcharge_density_dz(self, arg: float, /) -> None: ...

    @property
    def kick_csr(self) -> float:
        """CSR kick"""

    @kick_csr.setter
    def kick_csr(self, arg: float, /) -> None: ...

    @property
    def coef_lsc_plus(self) -> RealArray2D:
        """LSC Kick coefs."""

    @coef_lsc_plus.setter
    def coef_lsc_plus(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def coef_lsc_minus(self) -> RealArray2D:
        """LSC Kick coefs."""

    @coef_lsc_minus.setter
    def coef_lsc_minus(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def kick_lsc(self) -> float: ...

    @kick_lsc.setter
    def kick_lsc(self, arg: float, /) -> None: ...

    @property
    def n_particle(self) -> float:
        """
        Number of particles in slice can be a fraction since particles span multiple bins.
        """

    @n_particle.setter
    def n_particle(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> CsrBunchSliceStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> CsrBunchSliceStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> CsrBunchSliceStruct: ...

    def __deepcopy__(self, arg: dict, /) -> CsrBunchSliceStruct: ...

    def __eq__(self, arg: CsrBunchSliceStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class CsrEleInfoStruct:
    """Fortran struct: csr_ele_info_struct"""

    def __init__(self, ele: EleStruct | None = None, orbit0: CoordStruct | None = None, orbit1: CoordStruct | None = None, floor0: FloorPositionStruct | None = None, floor1: FloorPositionStruct | None = None, ref_floor0: FloorPositionStruct | None = None, ref_floor1: FloorPositionStruct | None = None, spline: SplineStruct | None = None, theta_chord: float | None = None, L_chord: float | None = None, dL_s: float | None = None) -> None: ...

    @property
    def ele(self) -> EleStruct | None:
        """lattice element"""

    @ele.setter
    def ele(self, arg: EleStruct, /) -> None: ...

    @property
    def orbit0(self) -> CoordStruct:
        """centroid orbit at entrance/exit ends"""

    @orbit0.setter
    def orbit0(self, arg: CoordStruct, /) -> None: ...

    @property
    def orbit1(self) -> CoordStruct:
        """centroid orbit at entrance/exit ends"""

    @orbit1.setter
    def orbit1(self, arg: CoordStruct, /) -> None: ...

    @property
    def floor0(self) -> FloorPositionStruct:
        """Floor position of centroid at entrance/exit ends"""

    @floor0.setter
    def floor0(self, arg: FloorPositionStruct, /) -> None: ...

    @property
    def floor1(self) -> FloorPositionStruct:
        """Floor position of centroid at entrance/exit ends"""

    @floor1.setter
    def floor1(self, arg: FloorPositionStruct, /) -> None: ...

    @property
    def ref_floor0(self) -> FloorPositionStruct:
        """Floor position of element ref coords at entrance/exit ends"""

    @ref_floor0.setter
    def ref_floor0(self, arg: FloorPositionStruct, /) -> None: ...

    @property
    def ref_floor1(self) -> FloorPositionStruct:
        """Floor position of element ref coords at entrance/exit ends"""

    @ref_floor1.setter
    def ref_floor1(self, arg: FloorPositionStruct, /) -> None: ...

    @property
    def spline(self) -> SplineStruct:
        """
        Spline for centroid orbit. spline%x = distance along chord. The spline is zero at the ends by construction.
        """

    @spline.setter
    def spline(self, arg: SplineStruct, /) -> None: ...

    @property
    def theta_chord(self) -> float:
        """Reference angle of chord in z-x plane"""

    @theta_chord.setter
    def theta_chord(self, arg: float, /) -> None: ...

    @property
    def L_chord(self) -> float:
        """Chord Length. Negative if bunch moves backwards in element."""

    @L_chord.setter
    def L_chord(self, arg: float, /) -> None: ...

    @property
    def dL_s(self) -> float:
        """L_s(of element) - L_chord"""

    @dL_s.setter
    def dL_s(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> CsrEleInfoStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> CsrEleInfoStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> CsrEleInfoStruct: ...

    def __deepcopy__(self, arg: dict, /) -> CsrEleInfoStruct: ...

    def __eq__(self, arg: CsrEleInfoStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class CsrKick1Struct:
    """Fortran struct: csr_kick1_struct"""

    def __init__(self, I_csr: float | None = None, I_int_csr: float | None = None, image_kick_csr: float | None = None, L_vec: Sequence[float] | None = None, L: float | None = None, dL: float | None = None, dz_particles: float | None = None, s_chord_source: float | None = None, theta_L: float | None = None, theta_sl: float | None = None, theta_lk: float | None = None, ix_ele_source: int | None = None, floor_s: FloorPositionStruct | None = None) -> None: ...

    @property
    def I_csr(self) -> float:
        """Kick integral."""

    @I_csr.setter
    def I_csr(self, arg: float, /) -> None: ...

    @property
    def I_int_csr(self) -> float:
        """Integrated Kick integral."""

    @I_int_csr.setter
    def I_int_csr(self, arg: float, /) -> None: ...

    @property
    def image_kick_csr(self) -> float:
        """kick."""

    @image_kick_csr.setter
    def image_kick_csr(self, arg: float, /) -> None: ...

    @property
    def L_vec(self) -> RealArray1D:
        """L vector in global coordinates."""

    @L_vec.setter
    def L_vec(self, arg: Sequence[float], /) -> None: ...

    @property
    def L(self) -> float:
        """Distance between source and kick points."""

    @L.setter
    def L(self, arg: float, /) -> None: ...

    @property
    def dL(self) -> float:
        """= epsilon_L = Ls - L"""

    @dL.setter
    def dL(self, arg: float, /) -> None: ...

    @property
    def dz_particles(self) -> float:
        """Kicked particle - source particle position at constant time."""

    @dz_particles.setter
    def dz_particles(self, arg: float, /) -> None: ...

    @property
    def s_chord_source(self) -> float:
        """Source point coordinate along chord."""

    @s_chord_source.setter
    def s_chord_source(self, arg: float, /) -> None: ...

    @property
    def theta_L(self) -> float:
        """Angle of L vector"""

    @theta_L.setter
    def theta_L(self, arg: float, /) -> None: ...

    @property
    def theta_sl(self) -> float:
        """Angle between velocity of particle at source pt and L"""

    @theta_sl.setter
    def theta_sl(self, arg: float, /) -> None: ...

    @property
    def theta_lk(self) -> float:
        """Angle between L and velocity of kicked particle"""

    @theta_lk.setter
    def theta_lk(self, arg: float, /) -> None: ...

    @property
    def ix_ele_source(self) -> int:
        """Source element index."""

    @ix_ele_source.setter
    def ix_ele_source(self, arg: int, /) -> None: ...

    @property
    def floor_s(self) -> FloorPositionStruct:
        """Floor position of source pt"""

    @floor_s.setter
    def floor_s(self, arg: FloorPositionStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> CsrKick1StructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> CsrKick1StructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> CsrKick1Struct: ...

    def __deepcopy__(self, arg: dict, /) -> CsrKick1Struct: ...

    def __eq__(self, arg: CsrKick1Struct, /) -> bool: ...

    def __hash__(self) -> int: ...

class CsrParticlePositionStruct:
    """Fortran struct: csr_particle_position_struct"""

    def __init__(self, r: Sequence[float] | None = None, charge: float | None = None) -> None: ...

    @property
    def r(self) -> RealArray1D:
        """particle position"""

    @r.setter
    def r(self, arg: Sequence[float], /) -> None: ...

    @property
    def charge(self) -> float:
        """particle charge"""

    @charge.setter
    def charge(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> CsrParticlePositionStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> CsrParticlePositionStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> CsrParticlePositionStruct: ...

    def __deepcopy__(self, arg: dict, /) -> CsrParticlePositionStruct: ...

    def __eq__(self, arg: CsrParticlePositionStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class CsrStruct:
    """Fortran struct: csr_struct"""

    def __init__(self, gamma: float | None = None, gamma2: float | None = None, rel_mass: float | None = None, beta: float | None = None, dz_slice: float | None = None, ds_track_step: float | None = None, s_kick: float | None = None, s_chord_kick: float | None = None, y_source: float | None = None, kick_factor: float | None = None, actual_track_step: float | None = None, x0_bunch: float | None = None, y0_bunch: float | None = None, floor_k: FloorPositionStruct | None = None, species: int | None = None, ix_ele_kick: int | None = None, kick_ele: EleStruct | None = None, mesh3d: Mesh3DStruct | None = None) -> None: ...

    @property
    def gamma(self) -> float:
        """Relativistic gamma factor."""

    @gamma.setter
    def gamma(self, arg: float, /) -> None: ...

    @property
    def gamma2(self) -> float:
        """Relativistic gamma factor."""

    @gamma2.setter
    def gamma2(self, arg: float, /) -> None: ...

    @property
    def rel_mass(self) -> float:
        """m_particle / m_electron"""

    @rel_mass.setter
    def rel_mass(self, arg: float, /) -> None: ...

    @property
    def beta(self) -> float:
        """Relativistic beta factor."""

    @beta.setter
    def beta(self, arg: float, /) -> None: ...

    @property
    def dz_slice(self) -> float:
        """Bin width"""

    @dz_slice.setter
    def dz_slice(self, arg: float, /) -> None: ...

    @property
    def ds_track_step(self) -> float:
        """True step size"""

    @ds_track_step.setter
    def ds_track_step(self, arg: float, /) -> None: ...

    @property
    def s_kick(self) -> float:
        """
        Kick point longitudinal location (element ref coords) from entrance end
        """

    @s_kick.setter
    def s_kick(self, arg: float, /) -> None: ...

    @property
    def s_chord_kick(self) -> float:
        """Kick point along beam centroid line"""

    @s_chord_kick.setter
    def s_chord_kick(self, arg: float, /) -> None: ...

    @property
    def y_source(self) -> float:
        """Height of source particle."""

    @y_source.setter
    def y_source(self, arg: float, /) -> None: ...

    @property
    def kick_factor(self) -> float:
        """Coefficient to scale the kick"""

    @kick_factor.setter
    def kick_factor(self, arg: float, /) -> None: ...

    @property
    def actual_track_step(self) -> float:
        """ds_track_step scalled by Length_centroid_chord / Length_element ratio"""

    @actual_track_step.setter
    def actual_track_step(self, arg: float, /) -> None: ...

    @property
    def x0_bunch(self) -> float:
        """Bunch centroid"""

    @x0_bunch.setter
    def x0_bunch(self, arg: float, /) -> None: ...

    @property
    def y0_bunch(self) -> float:
        """Bunch centroid"""

    @y0_bunch.setter
    def y0_bunch(self, arg: float, /) -> None: ...

    @property
    def floor_k(self) -> FloorPositionStruct:
        """Floor coords at kick point"""

    @floor_k.setter
    def floor_k(self, arg: FloorPositionStruct, /) -> None: ...

    @property
    def species(self) -> int:
        """Particle type"""

    @species.setter
    def species(self, arg: int, /) -> None: ...

    @property
    def ix_ele_kick(self) -> int:
        """Same as element being tracked through."""

    @ix_ele_kick.setter
    def ix_ele_kick(self, arg: int, /) -> None: ...

    @property
    def slice(self) -> CsrBunchSliceStructAlloc1D:
        """slice(i) refers to the i^th bunch slice."""

    @property
    def kick1(self) -> CsrKick1StructAlloc1D:
        """kick1(i) referes to the kick between two slices i bins apart."""

    @property
    def eleinfo(self) -> CsrEleInfoStructAlloc1D:
        """Element-by-element information."""

    @property
    def kick_ele(self) -> EleStruct | None:
        """Element where the kick pt is == ele tracked through."""

    @kick_ele.setter
    def kick_ele(self, arg: EleStruct, /) -> None: ...

    @property
    def mesh3d(self) -> Mesh3DStruct: ...

    @mesh3d.setter
    def mesh3d(self, arg: Mesh3DStruct, /) -> None: ...

    @property
    def position(self) -> CsrParticlePositionStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> CsrStruct: ...

    def __deepcopy__(self, arg: dict, /) -> CsrStruct: ...

    def __eq__(self, arg: CsrStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class CylindricalMapStruct:
    """Fortran struct: cylindrical_map_struct"""

    def __init__(self, m: int | None = None, harmonic: int | None = None, phi0_fieldmap: float | None = None, theta0_azimuth: float | None = None, field_scale: float | None = None, master_parameter: int | None = None, ele_anchor_pt: int | None = None, dz: float | None = None, r0: Sequence[float] | None = None, ptr: CylindricalMapTermStruct | None = None) -> None: ...

    @property
    def m(self) -> int:
        """Azimuthal Mode: varies as cos(m*phi - theta0_azimuth)"""

    @m.setter
    def m(self, arg: int, /) -> None: ...

    @property
    def harmonic(self) -> int:
        """Harmonic of fundamental"""

    @harmonic.setter
    def harmonic(self, arg: int, /) -> None: ...

    @property
    def phi0_fieldmap(self) -> float:
        """Mode oscillates as: twopi * (f * t + phi0_fieldmap)"""

    @phi0_fieldmap.setter
    def phi0_fieldmap(self, arg: float, /) -> None: ...

    @property
    def theta0_azimuth(self) -> float:
        """Azimuthal ((x, y) plane) orientation of mode."""

    @theta0_azimuth.setter
    def theta0_azimuth(self, arg: float, /) -> None: ...

    @property
    def field_scale(self) -> float:
        """Factor to scale the fields by"""

    @field_scale.setter
    def field_scale(self, arg: float, /) -> None: ...

    @property
    def master_parameter(self) -> int:
        """Master parameter in ele%value(:) array to use for scaling the field."""

    @master_parameter.setter
    def master_parameter(self, arg: int, /) -> None: ...

    @property
    def ele_anchor_pt(self) -> int:
        """anchor_beginning$, anchor_center$, or anchor_end$"""

    @ele_anchor_pt.setter
    def ele_anchor_pt(self, arg: int, /) -> None: ...

    @property
    def dz(self) -> float:
        """Distance between sampled field points."""

    @dz.setter
    def dz(self, arg: float, /) -> None: ...

    @property
    def r0(self) -> RealArray1D:
        """Field origin offset."""

    @r0.setter
    def r0(self, arg: Sequence[float], /) -> None: ...

    @property
    def ptr(self) -> CylindricalMapTermStruct | None: ...

    @ptr.setter
    def ptr(self, arg: CylindricalMapTermStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> CylindricalMapStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> CylindricalMapStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> CylindricalMapStruct: ...

    def __deepcopy__(self, arg: dict, /) -> CylindricalMapStruct: ...

    def __eq__(self, arg: CylindricalMapStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class CylindricalMapTerm1Struct:
    """Fortran struct: cylindrical_map_term1_struct"""

    def __init__(self, e_coef: complex | None = None, b_coef: complex | None = None) -> None: ...

    @property
    def e_coef(self) -> complex: ...

    @e_coef.setter
    def e_coef(self, arg: complex, /) -> None: ...

    @property
    def b_coef(self) -> complex: ...

    @b_coef.setter
    def b_coef(self, arg: complex, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> CylindricalMapTerm1StructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> CylindricalMapTerm1StructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> CylindricalMapTerm1Struct: ...

    def __deepcopy__(self, arg: dict, /) -> CylindricalMapTerm1Struct: ...

    def __eq__(self, arg: CylindricalMapTerm1Struct, /) -> bool: ...

    def __hash__(self) -> int: ...

class CylindricalMapTermStruct:
    """Fortran struct: cylindrical_map_term_struct"""

    def __init__(self, file: str | None = None, n_link: int | None = None) -> None: ...

    @property
    def file(self) -> str:
        """Input file name. Used also as ID for instances."""

    @file.setter
    def file(self, arg: str, /) -> None: ...

    @property
    def n_link(self) -> int:
        """For memory management of this structure"""

    @n_link.setter
    def n_link(self, arg: int, /) -> None: ...

    @property
    def term(self) -> CylindricalMapTerm1StructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> CylindricalMapTermStruct: ...

    def __deepcopy__(self, arg: dict, /) -> CylindricalMapTermStruct: ...

    def __eq__(self, arg: CylindricalMapTermStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class DiffuseParamStruct:
    """Fortran struct: diffuse_param_struct"""

    def __init__(self, x: float | None = None, y: float | None = None, lambda_: float | None = None, c_norm: float | None = None, chx_norm: float | None = None, n_pt_spline: int | None = None) -> None: ...

    @property
    def x(self) -> float: ...

    @x.setter
    def x(self, arg: float, /) -> None: ...

    @property
    def y(self) -> float: ...

    @y.setter
    def y(self, arg: float, /) -> None: ...

    @property
    def c_norm(self) -> float: ...

    @c_norm.setter
    def c_norm(self, arg: float, /) -> None: ...

    @property
    def chx_norm(self) -> float: ...

    @chx_norm.setter
    def chx_norm(self, arg: float, /) -> None: ...

    @property
    def prob_spline(self) -> SplineStructArray1D: ...

    @property
    def n_pt_spline(self) -> int: ...

    @n_pt_spline.setter
    def n_pt_spline(self, arg: int, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> DiffuseParamStruct: ...

    def __deepcopy__(self, arg: dict, /) -> DiffuseParamStruct: ...

    def __eq__(self, arg: DiffuseParamStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class DoLoopStruct:
    """Fortran struct: do_loop_struct"""

    def __init__(self, name: str | None = None, index: int | None = None, start: int | None = None, end: int | None = None, step: int | None = None, n_line_start: int | None = None, n_line_end: int | None = None, value: int | None = None) -> None: ...

    @property
    def name(self) -> str:
        """do loop index name"""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def index(self) -> int:
        """for do loops"""

    @index.setter
    def index(self, arg: int, /) -> None: ...

    @property
    def start(self) -> int:
        """for do loops"""

    @start.setter
    def start(self, arg: int, /) -> None: ...

    @property
    def end(self) -> int:
        """for do loops"""

    @end.setter
    def end(self, arg: int, /) -> None: ...

    @property
    def step(self) -> int:
        """for do loops"""

    @step.setter
    def step(self, arg: int, /) -> None: ...

    @property
    def n_line_start(self) -> int:
        """lines in each nested loop"""

    @n_line_start.setter
    def n_line_start(self, arg: int, /) -> None: ...

    @property
    def n_line_end(self) -> int:
        """lines in each nested loop"""

    @n_line_end.setter
    def n_line_end(self, arg: int, /) -> None: ...

    @property
    def value(self) -> int: ...

    @value.setter
    def value(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> DoLoopStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> DoLoopStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> DoLoopStruct: ...

    def __deepcopy__(self, arg: dict, /) -> DoLoopStruct: ...

    def __eq__(self, arg: DoLoopStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class EleAttributeStruct:
    """Fortran struct: ele_attribute_struct"""

    def __init__(self, name: str | None = None, state: int | None = None, kind: int | None = None, units: str | None = None, ix_attrib: int | None = None, value: float | None = None) -> None: ...

    @property
    def name(self) -> str: ...

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def state(self) -> int:
        """See above."""

    @state.setter
    def state(self, arg: int, /) -> None: ...

    @property
    def kind(self) -> int:
        """Is_switch$, is_real$, etc. See attribute_type routine."""

    @kind.setter
    def kind(self, arg: int, /) -> None: ...

    @property
    def units(self) -> str:
        """EG: 'T*m'."""

    @units.setter
    def units(self, arg: str, /) -> None: ...

    @property
    def ix_attrib(self) -> int:
        """
        Attribute index. Frequently will be where in the ele%value(:) array the attribute is.
        """

    @ix_attrib.setter
    def ix_attrib(self, arg: int, /) -> None: ...

    @property
    def value(self) -> float:
        """Used by type_ele."""

    @value.setter
    def value(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> EleAttributeStruct: ...

    def __deepcopy__(self, arg: dict, /) -> EleAttributeStruct: ...

    def __eq__(self, arg: EleAttributeStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ElePointerStruct:
    """Fortran struct: ele_pointer_struct"""

    def __init__(self, ele: EleStruct | None = None, loc: LatEleLocStruct | None = None, id: int | None = None) -> None: ...

    @property
    def ele(self) -> EleStruct | None: ...

    @ele.setter
    def ele(self, arg: EleStruct, /) -> None: ...

    @property
    def loc(self) -> LatEleLocStruct: ...

    @loc.setter
    def loc(self, arg: LatEleLocStruct, /) -> None: ...

    @property
    def id(self) -> int:
        """
        For general use. Not used by Bmad. In particular, used by Tao to designate universe ele is in.
        """

    @id.setter
    def id(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> ElePointerStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> ElePointerStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ElePointerStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ElePointerStruct: ...

    def __eq__(self, arg: ElePointerStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class EleStruct:
    """Fortran struct: ele_struct"""

    def __init__(self, name: str | None = None, type: str | None = None, alias: str | None = None, component_name: str | None = None, descrip: str | None = None, a: TwissStruct | None = None, b: TwissStruct | None = None, z: TwissStruct | None = None, x: XyDispStruct | None = None, y: XyDispStruct | None = None, ac_kick: AcKickerStruct | None = None, bookkeeping_state: BookkeepingStateStruct | None = None, branch: BranchStruct | None = None, control: ControllerStruct | None = None, converter: ConverterStruct | None = None, rf: RfEleStruct | None = None, foil: FoilStruct | None = None, lord: EleStruct | None = None, floor: FloorPositionStruct | None = None, high_energy_space_charge: HighEnergySpaceChargeStruct | None = None, mode3: Mode3Struct | None = None, photon: PhotonElementStruct | None = None, rad_map: RadMapEleStruct | None = None, spin_taylor_ref_orb_in: Sequence[float] | None = None, wake: WakeStruct | None = None, map_ref_orb_in: CoordStruct | None = None, map_ref_orb_out: CoordStruct | None = None, time_ref_orb_in: CoordStruct | None = None, time_ref_orb_out: CoordStruct | None = None, value: Sequence[float] | None = None, old_value: Sequence[float] | None = None, spin_q: Sequence[Sequence[float]] | None = None, vec0: Sequence[float] | None = None, mat6: Sequence[Sequence[float]] | None = None, c_mat: Sequence[Sequence[float]] | None = None, dc_mat_dpz: Sequence[Sequence[float]] | None = None, gamma_c: float | None = None, s_start: float | None = None, s: float | None = None, ref_time: float | None = None, a_pole: Sequence[float] | None = None, b_pole: Sequence[float] | None = None, a_pole_elec: Sequence[float] | None = None, b_pole_elec: Sequence[float] | None = None, custom: Sequence[float] | None = None, r: Sequence[Sequence[Sequence[float]]] | None = None, key: int | None = None, sub_key: int | None = None, ix_ele: int | None = None, ix_branch: int | None = None, lord_status: int | None = None, n_slave: int | None = None, n_slave_field: int | None = None, ix1_slave: int | None = None, slave_status: int | None = None, n_lord: int | None = None, n_lord_field: int | None = None, n_lord_ramper: int | None = None, ic1_lord: int | None = None, ix_pointer: int | None = None, ixx: int | None = None, iyy: int | None = None, izz: int | None = None, mat6_calc_method: int | None = None, tracking_method: int | None = None, spin_tracking_method: int | None = None, csr_method: int | None = None, space_charge_method: int | None = None, ptc_integration_type: int | None = None, field_calc: int | None = None, aperture_at: int | None = None, aperture_type: int | None = None, ref_species: int | None = None, orientation: int | None = None, symplectify: bool | None = None, mode_flip: bool | None = None, multipoles_on: bool | None = None, scale_multipoles: bool | None = None, taylor_map_includes_offsets: bool | None = None, field_master: bool | None = None, is_on: bool | None = None, logic: bool | None = None, bmad_logic: bool | None = None, select: bool | None = None, offset_moves_aperture: bool | None = None) -> None: ...

    @property
    def name(self) -> str:
        """name of element."""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def type(self) -> str:
        """type name."""

    @type.setter
    def type(self, arg: str, /) -> None: ...

    @property
    def alias(self) -> str:
        """Another name."""

    @alias.setter
    def alias(self, arg: str, /) -> None: ...

    @property
    def component_name(self) -> str:
        """Used by overlays, multipass patch, etc."""

    @component_name.setter
    def component_name(self, arg: str, /) -> None: ...

    @property
    def descrip(self) -> str:
        """Description string."""

    @descrip.setter
    def descrip(self, arg: str, /) -> None: ...

    @property
    def a(self) -> TwissStruct:
        """Twiss parameters at end of element"""

    @a.setter
    def a(self, arg: TwissStruct, /) -> None: ...

    @property
    def b(self) -> TwissStruct:
        """Twiss parameters at end of element"""

    @b.setter
    def b(self, arg: TwissStruct, /) -> None: ...

    @property
    def z(self) -> TwissStruct:
        """Twiss parameters at end of element"""

    @z.setter
    def z(self, arg: TwissStruct, /) -> None: ...

    @property
    def x(self) -> XyDispStruct:
        """Projected dispersions."""

    @x.setter
    def x(self, arg: XyDispStruct, /) -> None: ...

    @property
    def y(self) -> XyDispStruct:
        """Projected dispersions."""

    @y.setter
    def y(self, arg: XyDispStruct, /) -> None: ...

    @property
    def ac_kick(self) -> AcKickerStruct | None:
        """ac_kicker element parameters."""

    @ac_kick.setter
    def ac_kick(self, arg: AcKickerStruct, /) -> None: ...

    @property
    def bookkeeping_state(self) -> BookkeepingStateStruct:
        """Attribute bookkeeping"""

    @bookkeeping_state.setter
    def bookkeeping_state(self, arg: BookkeepingStateStruct, /) -> None: ...

    @property
    def branch(self) -> BranchStruct | None:
        """Pointer to branch containing element."""

    @branch.setter
    def branch(self, arg: BranchStruct, /) -> None: ...

    @property
    def control(self) -> ControllerStruct | None:
        """group & overlay variables."""

    @control.setter
    def control(self, arg: ControllerStruct, /) -> None: ...

    @property
    def converter(self) -> ConverterStruct | None:
        """EG: Positron converter in linac."""

    @converter.setter
    def converter(self, arg: ConverterStruct, /) -> None: ...

    @property
    def rf(self) -> RfEleStruct | None:
        """RF parameters."""

    @rf.setter
    def rf(self, arg: RfEleStruct, /) -> None: ...

    @property
    def foil(self) -> FoilStruct | None: ...

    @foil.setter
    def foil(self, arg: FoilStruct, /) -> None: ...

    @property
    def lord(self) -> EleStruct | None:
        """Pointer to a slice lord."""

    @lord.setter
    def lord(self, arg: EleStruct, /) -> None: ...

    @property
    def floor(self) -> FloorPositionStruct: ...

    @floor.setter
    def floor(self, arg: FloorPositionStruct, /) -> None: ...

    @property
    def high_energy_space_charge(self) -> HighEnergySpaceChargeStruct | None: ...

    @high_energy_space_charge.setter
    def high_energy_space_charge(self, arg: HighEnergySpaceChargeStruct, /) -> None: ...

    @property
    def mode3(self) -> Mode3Struct | None:
        """6D normal mode structure."""

    @mode3.setter
    def mode3(self, arg: Mode3Struct, /) -> None: ...

    @property
    def photon(self) -> PhotonElementStruct | None: ...

    @photon.setter
    def photon(self, arg: PhotonElementStruct, /) -> None: ...

    @property
    def rad_map(self) -> RadMapEleStruct | None:
        """
        Radiation kick parameters Note: The reference orbits for spin and orbit Taylor maps are not necessarily the same. For example, Sprint spin Taylor maps can be with respect to the zero orbit independent of the orbital map.
        """

    @rad_map.setter
    def rad_map(self, arg: RadMapEleStruct, /) -> None: ...

    @property
    def taylor(self) -> TaylorStructArray1D:
        """Phase space Taylor map."""

    @property
    def spin_taylor_ref_orb_in(self) -> RealArray1D: ...

    @spin_taylor_ref_orb_in.setter
    def spin_taylor_ref_orb_in(self, arg: Sequence[float], /) -> None: ...

    @property
    def spin_taylor(self) -> TaylorStructArray1D:
        """Quaternion Spin Taylor map."""

    @property
    def wake(self) -> WakeStruct | None:
        """Wakes"""

    @wake.setter
    def wake(self, arg: WakeStruct, /) -> None: ...

    @property
    def wall3d(self) -> Wall3DStructArray1D:
        """Chamber or capillary wall E/M field structs."""

    @property
    def cartesian_map(self) -> CartesianMapStructArray1D:
        """Used to define E/M fields"""

    @property
    def cylindrical_map(self) -> CylindricalMapStructArray1D:
        """Used to define E/M fields"""

    @property
    def gen_grad_map(self) -> GenGradMapStructArray1D:
        """Used to define E/M fields."""

    @property
    def grid_field(self) -> GridFieldStructArray1D:
        """
        Used to define E/M fields. The difference between map_ref_orb and time_ref_orb is that map_ref_orb is the reference orbit for the 1st order spin/orbit map which, in general, is non-zero while time_ref_orb follows the reference particle which is generally the zero orbit (non-zero, for example, in the second slice of a sliced wiggler).
        """

    @property
    def map_ref_orb_in(self) -> CoordStruct:
        """Entrance end transfer map ref orbit"""

    @map_ref_orb_in.setter
    def map_ref_orb_in(self, arg: CoordStruct, /) -> None: ...

    @property
    def map_ref_orb_out(self) -> CoordStruct:
        """Exit end transfer map ref orbit"""

    @map_ref_orb_out.setter
    def map_ref_orb_out(self, arg: CoordStruct, /) -> None: ...

    @property
    def time_ref_orb_in(self) -> CoordStruct:
        """Reference orbit at entrance end for ref_time calc."""

    @time_ref_orb_in.setter
    def time_ref_orb_in(self, arg: CoordStruct, /) -> None: ...

    @property
    def time_ref_orb_out(self) -> CoordStruct:
        """Reference orbit at exit end for ref_time calc."""

    @time_ref_orb_out.setter
    def time_ref_orb_out(self, arg: CoordStruct, /) -> None: ...

    @property
    def value(self) -> RealArray1D:
        """attribute values."""

    @value.setter
    def value(self, arg: Sequence[float], /) -> None: ...

    @property
    def old_value(self) -> RealArray1D:
        """
        Used to see if %value(:) array has changed. Note: The reference orbit for spin/orbit matrices is %map_ref_orb_in/out
        """

    @old_value.setter
    def old_value(self, arg: Sequence[float], /) -> None: ...

    @property
    def spin_q(self) -> RealArray2D:
        """0th and 1st order Spin transport quaternion."""

    @spin_q.setter
    def spin_q(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def vec0(self) -> RealArray1D:
        """0th order transport vector."""

    @vec0.setter
    def vec0(self, arg: Sequence[float], /) -> None: ...

    @property
    def mat6(self) -> RealArray2D:
        """1st order transport matrix."""

    @mat6.setter
    def mat6(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def c_mat(self) -> RealArray2D:
        """2x2 C coupling matrix"""

    @c_mat.setter
    def c_mat(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def dc_mat_dpz(self) -> RealArray2D:
        """d(c_mat)/dpz variation."""

    @dc_mat_dpz.setter
    def dc_mat_dpz(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def gamma_c(self) -> float:
        """gamma associated with C matrix"""

    @gamma_c.setter
    def gamma_c(self, arg: float, /) -> None: ...

    @property
    def s_start(self) -> float:
        """longitudinal ref position at entrance_end"""

    @s_start.setter
    def s_start(self, arg: float, /) -> None: ...

    @property
    def s(self) -> float:
        """longitudinal ref position at the exit end."""

    @s.setter
    def s(self, arg: float, /) -> None: ...

    @property
    def ref_time(self) -> float:
        """Time ref particle passes exit end."""

    @ref_time.setter
    def ref_time(self, arg: float, /) -> None: ...

    @property
    def a_pole(self) -> RealArray1D:
        """knl for multipole elements."""

    @a_pole.setter
    def a_pole(self, arg: Sequence[float], /) -> None: ...

    @property
    def b_pole(self) -> RealArray1D:
        """tilt for multipole elements."""

    @b_pole.setter
    def b_pole(self, arg: Sequence[float], /) -> None: ...

    @property
    def a_pole_elec(self) -> RealArray1D:
        """Electrostatic multipoles. ksnl for multipole elements."""

    @a_pole_elec.setter
    def a_pole_elec(self, arg: Sequence[float], /) -> None: ...

    @property
    def b_pole_elec(self) -> RealArray1D:
        """Electrostatic multipoles."""

    @b_pole_elec.setter
    def b_pole_elec(self, arg: Sequence[float], /) -> None: ...

    @property
    def custom(self) -> RealArray1D:
        """Custom attributes."""

    @custom.setter
    def custom(self, arg: Sequence[float], /) -> None: ...

    @property
    def r(self) -> RealArray3D:
        """For general use. Not used by Bmad."""

    @r.setter
    def r(self, arg: Sequence[Sequence[Sequence[float]]], /) -> None: ...

    @property
    def key(self) -> int:
        """Element class (quadrupole, etc.)."""

    @key.setter
    def key(self, arg: int, /) -> None: ...

    @property
    def sub_key(self) -> int:
        """Records bend input type."""

    @sub_key.setter
    def sub_key(self, arg: int, /) -> None: ...

    @property
    def ix_ele(self) -> int:
        """
        Index in branch ele(0:) array. Set to ix_slice_slave$ = -2 for slice_slave$ elements.
        """

    @ix_ele.setter
    def ix_ele(self, arg: int, /) -> None: ...

    @property
    def ix_branch(self) -> int:
        """Index in lat%branch(:) array. Note: lat%ele => lat%branch(0)."""

    @ix_branch.setter
    def ix_branch(self, arg: int, /) -> None: ...

    @property
    def lord_status(self) -> int:
        """Type of lord element this is. overlay_lord$, etc."""

    @lord_status.setter
    def lord_status(self, arg: int, /) -> None: ...

    @property
    def n_slave(self) -> int:
        """Number of slaves (except field overlap slaves) of this element."""

    @n_slave.setter
    def n_slave(self, arg: int, /) -> None: ...

    @property
    def n_slave_field(self) -> int:
        """Number of field slaves of this element."""

    @n_slave_field.setter
    def n_slave_field(self, arg: int, /) -> None: ...

    @property
    def ix1_slave(self) -> int:
        """Pointer index to this element's slaves."""

    @ix1_slave.setter
    def ix1_slave(self, arg: int, /) -> None: ...

    @property
    def slave_status(self) -> int:
        """Type of slave element this is. multipass_slave$, slice_slave$, etc."""

    @slave_status.setter
    def slave_status(self, arg: int, /) -> None: ...

    @property
    def n_lord(self) -> int:
        """Number of lords (except field overlap and ramper lords)."""

    @n_lord.setter
    def n_lord(self, arg: int, /) -> None: ...

    @property
    def n_lord_field(self) -> int:
        """Number of field lords of this element."""

    @n_lord_field.setter
    def n_lord_field(self, arg: int, /) -> None: ...

    @property
    def n_lord_ramper(self) -> int:
        """Number of ramper lords."""

    @n_lord_ramper.setter
    def n_lord_ramper(self, arg: int, /) -> None: ...

    @property
    def ic1_lord(self) -> int:
        """Pointer index to this element's lords."""

    @ic1_lord.setter
    def ic1_lord(self, arg: int, /) -> None: ...

    @property
    def ix_pointer(self) -> int:
        """For general use. Not used by Bmad."""

    @ix_pointer.setter
    def ix_pointer(self, arg: int, /) -> None: ...

    @property
    def ixx(self) -> int:
        """Index for Bmad internal use."""

    @ixx.setter
    def ixx(self, arg: int, /) -> None: ...

    @property
    def iyy(self) -> int:
        """Index for Bmad internal use."""

    @iyy.setter
    def iyy(self, arg: int, /) -> None: ...

    @property
    def izz(self) -> int:
        """Index for Bmad internal use."""

    @izz.setter
    def izz(self, arg: int, /) -> None: ...

    @property
    def mat6_calc_method(self) -> int:
        """taylor$, symp_lie_ptc$, etc."""

    @mat6_calc_method.setter
    def mat6_calc_method(self, arg: int, /) -> None: ...

    @property
    def tracking_method(self) -> int:
        """taylor$, linear$, etc."""

    @tracking_method.setter
    def tracking_method(self, arg: int, /) -> None: ...

    @property
    def spin_tracking_method(self) -> int:
        """symp_lie_ptc$, etc."""

    @spin_tracking_method.setter
    def spin_tracking_method(self, arg: int, /) -> None: ...

    @property
    def csr_method(self) -> int:
        """or one_dim$ ('1_dim'), steady_state_3d$"""

    @csr_method.setter
    def csr_method(self, arg: int, /) -> None: ...

    @property
    def space_charge_method(self) -> int:
        """
        slice$, slice_longitudinal$, slice_transverse$, fft_3D$, cathode_fft_3d$
        """

    @space_charge_method.setter
    def space_charge_method(self, arg: int, /) -> None: ...

    @property
    def ptc_integration_type(self) -> int:
        """drift_kick$, matrix_kick$, or ripken_kick$"""

    @ptc_integration_type.setter
    def ptc_integration_type(self, arg: int, /) -> None: ...

    @property
    def field_calc(self) -> int:
        """no_field$, fieldmap$, refer_to_lords$, or custom$"""

    @field_calc.setter
    def field_calc(self, arg: int, /) -> None: ...

    @property
    def aperture_at(self) -> int:
        """Aperture location: entrance_end$, ..."""

    @aperture_at.setter
    def aperture_at(self, arg: int, /) -> None: ...

    @property
    def aperture_type(self) -> int:
        """rectangular$, elliptical$, auto_aperture$, ..."""

    @aperture_type.setter
    def aperture_type(self, arg: int, /) -> None: ...

    @property
    def ref_species(self) -> int:
        """Reference species"""

    @ref_species.setter
    def ref_species(self, arg: int, /) -> None: ...

    @property
    def orientation(self) -> int:
        """-1 -> Element is longitudinally reversed. +1 -> Normal."""

    @orientation.setter
    def orientation(self, arg: int, /) -> None: ...

    @property
    def symplectify(self) -> bool:
        """Symplectify mat6 matrices."""

    @symplectify.setter
    def symplectify(self, arg: bool, /) -> None: ...

    @property
    def mode_flip(self) -> bool:
        """Have the normal modes traded places?"""

    @mode_flip.setter
    def mode_flip(self, arg: bool, /) -> None: ...

    @property
    def multipoles_on(self) -> bool:
        """For turning multipoles on/off"""

    @multipoles_on.setter
    def multipoles_on(self, arg: bool, /) -> None: ...

    @property
    def scale_multipoles(self) -> bool:
        """
        Are ab_multipoles within other elements (EG: quads, etc.) scaled by the strength of the element?
        """

    @scale_multipoles.setter
    def scale_multipoles(self, arg: bool, /) -> None: ...

    @property
    def taylor_map_includes_offsets(self) -> bool:
        """Taylor map calculated with element misalignments?"""

    @taylor_map_includes_offsets.setter
    def taylor_map_includes_offsets(self, arg: bool, /) -> None: ...

    @property
    def field_master(self) -> bool:
        """Calculate strength from the field value?"""

    @field_master.setter
    def field_master(self, arg: bool, /) -> None: ...

    @property
    def is_on(self) -> bool:
        """For turning element on/off."""

    @is_on.setter
    def is_on(self, arg: bool, /) -> None: ...

    @property
    def logic(self) -> bool:
        """For general use. Not used by Bmad (except during lattice parsing)."""

    @logic.setter
    def logic(self, arg: bool, /) -> None: ...

    @property
    def bmad_logic(self) -> bool:
        """For Bmad internal use only."""

    @bmad_logic.setter
    def bmad_logic(self, arg: bool, /) -> None: ...

    @property
    def select(self) -> bool:
        """For Bmad internal use only."""

    @select.setter
    def select(self, arg: bool, /) -> None: ...

    @property
    def offset_moves_aperture(self) -> bool:
        """element offsets affects aperture? ! final :: ele_finalizer"""

    @offset_moves_aperture.setter
    def offset_moves_aperture(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> EleStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> EleStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> EleStruct: ...

    def __deepcopy__(self, arg: dict, /) -> EleStruct: ...

    def __eq__(self, arg: EleStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class EllipseBeamInitStruct:
    """Fortran struct: ellipse_beam_init_struct"""

    def __init__(self, part_per_ellipse: int | None = None, n_ellipse: int | None = None, sigma_cutoff: float | None = None) -> None: ...

    @property
    def part_per_ellipse(self) -> int:
        """number of particles per ellipse"""

    @part_per_ellipse.setter
    def part_per_ellipse(self, arg: int, /) -> None: ...

    @property
    def n_ellipse(self) -> int:
        """number of ellipses (>= 1)"""

    @n_ellipse.setter
    def n_ellipse(self, arg: int, /) -> None: ...

    @property
    def sigma_cutoff(self) -> float:
        """sigma cutoff of the representation"""

    @sigma_cutoff.setter
    def sigma_cutoff(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> EllipseBeamInitStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> EllipseBeamInitStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> EllipseBeamInitStruct: ...

    def __deepcopy__(self, arg: dict, /) -> EllipseBeamInitStruct: ...

    def __eq__(self, arg: EllipseBeamInitStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class EmFieldStruct:
    """Fortran struct: em_field_struct"""

    def __init__(self, E: Sequence[float] | None = None, B: Sequence[float] | None = None, dE: Sequence[Sequence[float]] | None = None, dB: Sequence[Sequence[float]] | None = None, phi: float | None = None, phi_B: float | None = None, A: Sequence[float] | None = None) -> None: ...

    @property
    def E(self) -> RealArray1D:
        """electric field."""

    @E.setter
    def E(self, arg: Sequence[float], /) -> None: ...

    @property
    def B(self) -> RealArray1D:
        """magnetic field."""

    @B.setter
    def B(self, arg: Sequence[float], /) -> None: ...

    @property
    def dE(self) -> RealArray2D:
        """electric field gradient."""

    @dE.setter
    def dE(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def dB(self) -> RealArray2D:
        """magnetic field gradient."""

    @dB.setter
    def dB(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def phi(self) -> float:
        """Electric scalar potential."""

    @phi.setter
    def phi(self, arg: float, /) -> None: ...

    @property
    def phi_B(self) -> float:
        """Magnetic scalar potential."""

    @phi_B.setter
    def phi_B(self, arg: float, /) -> None: ...

    @property
    def A(self) -> RealArray1D:
        """Magnetic vector potential."""

    @A.setter
    def A(self, arg: Sequence[float], /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> EmFieldStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> EmFieldStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> EmFieldStruct: ...

    def __deepcopy__(self, arg: dict, /) -> EmFieldStruct: ...

    def __eq__(self, arg: EmFieldStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ExpressionAtomStruct:
    """Fortran struct: expression_atom_struct"""

    def __init__(self, name: str | None = None, type: int | None = None, value: float | None = None) -> None: ...

    @property
    def name(self) -> str: ...

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def type(self) -> int:
        """
        plus$, minum$, sin$, cos$, etc. To convert to string use: expression_op_name
        """

    @type.setter
    def type(self, arg: int, /) -> None: ...

    @property
    def value(self) -> float: ...

    @value.setter
    def value(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> ExpressionAtomStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> ExpressionAtomStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ExpressionAtomStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ExpressionAtomStruct: ...

    def __eq__(self, arg: ExpressionAtomStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ExpressionTreeStruct:
    """Fortran struct: expression_tree_struct"""

    def __init__(self, name: str | None = None, type: int | None = None, value: float | None = None) -> None: ...

    @property
    def name(self) -> str: ...

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def type(self) -> int:
        """plus$, minum$, sin$, cos$, etc."""

    @type.setter
    def type(self, arg: int, /) -> None: ...

    @property
    def value(self) -> float: ...

    @value.setter
    def value(self, arg: float, /) -> None: ...

    @property
    def node(self) -> ExpressionTreeStructArray1D:
        """
        Child nodes. Note: Pointer used here since Ifort does not support allocatable.
        """

    @staticmethod
    def new_array1d(sz: int = 0) -> ExpressionTreeStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> ExpressionTreeStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ExpressionTreeStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ExpressionTreeStruct: ...

    def __eq__(self, arg: ExpressionTreeStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ExtraParsingInfoStruct:
    """Fortran struct: extra_parsing_info_struct"""

    def __init__(self, ran_state: RandomStateStruct | None = None, ran_seed: int | None = None, undeterministic_ran_function_called: bool | None = None, d_orb_set: bool | None = None, max_aperture_limit_set: bool | None = None, default_ds_step_set: bool | None = None, significant_length_set: bool | None = None, rel_tol_tracking_set: bool | None = None, abs_tol_tracking_set: bool | None = None, rel_tol_adaptive_tracking_set: bool | None = None, abs_tol_adaptive_tracking_set: bool | None = None, init_ds_adaptive_tracking_set: bool | None = None, min_ds_adaptive_tracking_set: bool | None = None, fatal_ds_adaptive_tracking_set: bool | None = None, synch_rad_scale_set: bool | None = None, autoscale_amp_abs_tol_set: bool | None = None, autoscale_amp_rel_tol_set: bool | None = None, autoscale_phase_tol_set: bool | None = None, rf_phase_below_transition_ref_set: bool | None = None, electric_dipole_moment_set: bool | None = None, taylor_order_set: bool | None = None, runge_kutta_order_set: bool | None = None, default_integ_order_set: bool | None = None, sr_wakes_on_set: bool | None = None, lr_wakes_on_set: bool | None = None, high_energy_space_charge_on_set: bool | None = None, high_energy_space_charge_linear_set: bool | None = None, csr_and_space_charge_on_set: bool | None = None, spin_tracking_on_set: bool | None = None, spin_sokolov_ternov_flipping_on_set: bool | None = None, radiation_damping_on_set: bool | None = None, radiation_zero_average_set: bool | None = None, radiation_fluctuations_on_set: bool | None = None, conserve_taylor_maps_set: bool | None = None, absolute_time_tracking_set: bool | None = None, absolute_time_ref_shift_set: bool | None = None, convert_to_kinetic_momentum_set: bool | None = None, aperture_limit_on_set: bool | None = None, normalize_twiss_set: bool | None = None, sad_eps_scale_set: bool | None = None, sad_amp_max_set: bool | None = None, sad_n_div_max_set: bool | None = None, max_num_runge_kutta_step_set: bool | None = None, spin_n0_direction_user_set_set: bool | None = None, debug_set: bool | None = None, ds_track_step_set: bool | None = None, dt_track_step_set: bool | None = None, cathode_strength_cutoff_set: bool | None = None, sc_rel_tol_tracking_set: bool | None = None, sc_abs_tol_tracking_set: bool | None = None, beam_chamber_height_set: bool | None = None, lsc_sigma_cutoff_set: bool | None = None, particle_sigma_cutoff_set: bool | None = None, space_charge_mesh_size_set: bool | None = None, csr3d_mesh_size_set: bool | None = None, n_bin_set: bool | None = None, particle_bin_span_set: bool | None = None, n_shield_images_set: bool | None = None, sc_min_in_bin_set: bool | None = None, lsc_kick_transverse_dependence_set: bool | None = None, sc_debug_set: bool | None = None, diagnostic_output_file_set: bool | None = None, old_integrator_set: bool | None = None, use_orientation_patches_set: bool | None = None, print_info_messages_set: bool | None = None, max_fringe_order_set: bool | None = None, exact_model_set: bool | None = None, exact_misalign_set: bool | None = None, vertical_kick_set: bool | None = None, cut_factor_set: bool | None = None, translate_patch_drift_time_set: bool | None = None) -> None: ...

    @property
    def ran_state(self) -> RandomStateStruct: ...

    @ran_state.setter
    def ran_state(self, arg: RandomStateStruct, /) -> None: ...

    @property
    def ran_seed(self) -> int: ...

    @ran_seed.setter
    def ran_seed(self, arg: int, /) -> None: ...

    @property
    def undeterministic_ran_function_called(self) -> bool:
        """Used with bmad_com"""

    @undeterministic_ran_function_called.setter
    def undeterministic_ran_function_called(self, arg: bool, /) -> None: ...

    @property
    def d_orb_set(self) -> bool: ...

    @d_orb_set.setter
    def d_orb_set(self, arg: bool, /) -> None: ...

    @property
    def max_aperture_limit_set(self) -> bool: ...

    @max_aperture_limit_set.setter
    def max_aperture_limit_set(self, arg: bool, /) -> None: ...

    @property
    def default_ds_step_set(self) -> bool: ...

    @default_ds_step_set.setter
    def default_ds_step_set(self, arg: bool, /) -> None: ...

    @property
    def significant_length_set(self) -> bool: ...

    @significant_length_set.setter
    def significant_length_set(self, arg: bool, /) -> None: ...

    @property
    def rel_tol_tracking_set(self) -> bool: ...

    @rel_tol_tracking_set.setter
    def rel_tol_tracking_set(self, arg: bool, /) -> None: ...

    @property
    def abs_tol_tracking_set(self) -> bool: ...

    @abs_tol_tracking_set.setter
    def abs_tol_tracking_set(self, arg: bool, /) -> None: ...

    @property
    def rel_tol_adaptive_tracking_set(self) -> bool: ...

    @rel_tol_adaptive_tracking_set.setter
    def rel_tol_adaptive_tracking_set(self, arg: bool, /) -> None: ...

    @property
    def abs_tol_adaptive_tracking_set(self) -> bool: ...

    @abs_tol_adaptive_tracking_set.setter
    def abs_tol_adaptive_tracking_set(self, arg: bool, /) -> None: ...

    @property
    def init_ds_adaptive_tracking_set(self) -> bool: ...

    @init_ds_adaptive_tracking_set.setter
    def init_ds_adaptive_tracking_set(self, arg: bool, /) -> None: ...

    @property
    def min_ds_adaptive_tracking_set(self) -> bool: ...

    @min_ds_adaptive_tracking_set.setter
    def min_ds_adaptive_tracking_set(self, arg: bool, /) -> None: ...

    @property
    def fatal_ds_adaptive_tracking_set(self) -> bool: ...

    @fatal_ds_adaptive_tracking_set.setter
    def fatal_ds_adaptive_tracking_set(self, arg: bool, /) -> None: ...

    @property
    def synch_rad_scale_set(self) -> bool: ...

    @synch_rad_scale_set.setter
    def synch_rad_scale_set(self, arg: bool, /) -> None: ...

    @property
    def autoscale_amp_abs_tol_set(self) -> bool: ...

    @autoscale_amp_abs_tol_set.setter
    def autoscale_amp_abs_tol_set(self, arg: bool, /) -> None: ...

    @property
    def autoscale_amp_rel_tol_set(self) -> bool: ...

    @autoscale_amp_rel_tol_set.setter
    def autoscale_amp_rel_tol_set(self, arg: bool, /) -> None: ...

    @property
    def autoscale_phase_tol_set(self) -> bool: ...

    @autoscale_phase_tol_set.setter
    def autoscale_phase_tol_set(self, arg: bool, /) -> None: ...

    @property
    def rf_phase_below_transition_ref_set(self) -> bool: ...

    @rf_phase_below_transition_ref_set.setter
    def rf_phase_below_transition_ref_set(self, arg: bool, /) -> None: ...

    @property
    def electric_dipole_moment_set(self) -> bool: ...

    @electric_dipole_moment_set.setter
    def electric_dipole_moment_set(self, arg: bool, /) -> None: ...

    @property
    def taylor_order_set(self) -> bool: ...

    @taylor_order_set.setter
    def taylor_order_set(self, arg: bool, /) -> None: ...

    @property
    def runge_kutta_order_set(self) -> bool: ...

    @runge_kutta_order_set.setter
    def runge_kutta_order_set(self, arg: bool, /) -> None: ...

    @property
    def default_integ_order_set(self) -> bool: ...

    @default_integ_order_set.setter
    def default_integ_order_set(self, arg: bool, /) -> None: ...

    @property
    def sr_wakes_on_set(self) -> bool: ...

    @sr_wakes_on_set.setter
    def sr_wakes_on_set(self, arg: bool, /) -> None: ...

    @property
    def lr_wakes_on_set(self) -> bool: ...

    @lr_wakes_on_set.setter
    def lr_wakes_on_set(self, arg: bool, /) -> None: ...

    @property
    def high_energy_space_charge_on_set(self) -> bool: ...

    @high_energy_space_charge_on_set.setter
    def high_energy_space_charge_on_set(self, arg: bool, /) -> None: ...

    @property
    def high_energy_space_charge_linear_set(self) -> bool: ...

    @high_energy_space_charge_linear_set.setter
    def high_energy_space_charge_linear_set(self, arg: bool, /) -> None: ...

    @property
    def csr_and_space_charge_on_set(self) -> bool: ...

    @csr_and_space_charge_on_set.setter
    def csr_and_space_charge_on_set(self, arg: bool, /) -> None: ...

    @property
    def spin_tracking_on_set(self) -> bool: ...

    @spin_tracking_on_set.setter
    def spin_tracking_on_set(self, arg: bool, /) -> None: ...

    @property
    def spin_sokolov_ternov_flipping_on_set(self) -> bool: ...

    @spin_sokolov_ternov_flipping_on_set.setter
    def spin_sokolov_ternov_flipping_on_set(self, arg: bool, /) -> None: ...

    @property
    def radiation_damping_on_set(self) -> bool: ...

    @radiation_damping_on_set.setter
    def radiation_damping_on_set(self, arg: bool, /) -> None: ...

    @property
    def radiation_zero_average_set(self) -> bool: ...

    @radiation_zero_average_set.setter
    def radiation_zero_average_set(self, arg: bool, /) -> None: ...

    @property
    def radiation_fluctuations_on_set(self) -> bool: ...

    @radiation_fluctuations_on_set.setter
    def radiation_fluctuations_on_set(self, arg: bool, /) -> None: ...

    @property
    def conserve_taylor_maps_set(self) -> bool: ...

    @conserve_taylor_maps_set.setter
    def conserve_taylor_maps_set(self, arg: bool, /) -> None: ...

    @property
    def absolute_time_tracking_set(self) -> bool: ...

    @absolute_time_tracking_set.setter
    def absolute_time_tracking_set(self, arg: bool, /) -> None: ...

    @property
    def absolute_time_ref_shift_set(self) -> bool: ...

    @absolute_time_ref_shift_set.setter
    def absolute_time_ref_shift_set(self, arg: bool, /) -> None: ...

    @property
    def convert_to_kinetic_momentum_set(self) -> bool: ...

    @convert_to_kinetic_momentum_set.setter
    def convert_to_kinetic_momentum_set(self, arg: bool, /) -> None: ...

    @property
    def aperture_limit_on_set(self) -> bool: ...

    @aperture_limit_on_set.setter
    def aperture_limit_on_set(self, arg: bool, /) -> None: ...

    @property
    def normalize_twiss_set(self) -> bool: ...

    @normalize_twiss_set.setter
    def normalize_twiss_set(self, arg: bool, /) -> None: ...

    @property
    def sad_eps_scale_set(self) -> bool: ...

    @sad_eps_scale_set.setter
    def sad_eps_scale_set(self, arg: bool, /) -> None: ...

    @property
    def sad_amp_max_set(self) -> bool: ...

    @sad_amp_max_set.setter
    def sad_amp_max_set(self, arg: bool, /) -> None: ...

    @property
    def sad_n_div_max_set(self) -> bool: ...

    @sad_n_div_max_set.setter
    def sad_n_div_max_set(self, arg: bool, /) -> None: ...

    @property
    def max_num_runge_kutta_step_set(self) -> bool: ...

    @max_num_runge_kutta_step_set.setter
    def max_num_runge_kutta_step_set(self, arg: bool, /) -> None: ...

    @property
    def spin_n0_direction_user_set_set(self) -> bool: ...

    @spin_n0_direction_user_set_set.setter
    def spin_n0_direction_user_set_set(self, arg: bool, /) -> None: ...

    @property
    def debug_set(self) -> bool:
        """Used with space_charge_com"""

    @debug_set.setter
    def debug_set(self, arg: bool, /) -> None: ...

    @property
    def ds_track_step_set(self) -> bool: ...

    @ds_track_step_set.setter
    def ds_track_step_set(self, arg: bool, /) -> None: ...

    @property
    def dt_track_step_set(self) -> bool: ...

    @dt_track_step_set.setter
    def dt_track_step_set(self, arg: bool, /) -> None: ...

    @property
    def cathode_strength_cutoff_set(self) -> bool: ...

    @cathode_strength_cutoff_set.setter
    def cathode_strength_cutoff_set(self, arg: bool, /) -> None: ...

    @property
    def sc_rel_tol_tracking_set(self) -> bool:
        """For: space_charge_com%rel_tol_tracking"""

    @sc_rel_tol_tracking_set.setter
    def sc_rel_tol_tracking_set(self, arg: bool, /) -> None: ...

    @property
    def sc_abs_tol_tracking_set(self) -> bool:
        """For: space_charge_com%abs_tol_tracking"""

    @sc_abs_tol_tracking_set.setter
    def sc_abs_tol_tracking_set(self, arg: bool, /) -> None: ...

    @property
    def beam_chamber_height_set(self) -> bool: ...

    @beam_chamber_height_set.setter
    def beam_chamber_height_set(self, arg: bool, /) -> None: ...

    @property
    def lsc_sigma_cutoff_set(self) -> bool: ...

    @lsc_sigma_cutoff_set.setter
    def lsc_sigma_cutoff_set(self, arg: bool, /) -> None: ...

    @property
    def particle_sigma_cutoff_set(self) -> bool: ...

    @particle_sigma_cutoff_set.setter
    def particle_sigma_cutoff_set(self, arg: bool, /) -> None: ...

    @property
    def space_charge_mesh_size_set(self) -> bool: ...

    @space_charge_mesh_size_set.setter
    def space_charge_mesh_size_set(self, arg: bool, /) -> None: ...

    @property
    def csr3d_mesh_size_set(self) -> bool: ...

    @csr3d_mesh_size_set.setter
    def csr3d_mesh_size_set(self, arg: bool, /) -> None: ...

    @property
    def n_bin_set(self) -> bool: ...

    @n_bin_set.setter
    def n_bin_set(self, arg: bool, /) -> None: ...

    @property
    def particle_bin_span_set(self) -> bool: ...

    @particle_bin_span_set.setter
    def particle_bin_span_set(self, arg: bool, /) -> None: ...

    @property
    def n_shield_images_set(self) -> bool: ...

    @n_shield_images_set.setter
    def n_shield_images_set(self, arg: bool, /) -> None: ...

    @property
    def sc_min_in_bin_set(self) -> bool: ...

    @sc_min_in_bin_set.setter
    def sc_min_in_bin_set(self, arg: bool, /) -> None: ...

    @property
    def lsc_kick_transverse_dependence_set(self) -> bool: ...

    @lsc_kick_transverse_dependence_set.setter
    def lsc_kick_transverse_dependence_set(self, arg: bool, /) -> None: ...

    @property
    def sc_debug_set(self) -> bool: ...

    @sc_debug_set.setter
    def sc_debug_set(self, arg: bool, /) -> None: ...

    @property
    def diagnostic_output_file_set(self) -> bool:
        """Used with ptc_com"""

    @diagnostic_output_file_set.setter
    def diagnostic_output_file_set(self, arg: bool, /) -> None: ...

    @property
    def old_integrator_set(self) -> bool: ...

    @old_integrator_set.setter
    def old_integrator_set(self, arg: bool, /) -> None: ...

    @property
    def use_orientation_patches_set(self) -> bool: ...

    @use_orientation_patches_set.setter
    def use_orientation_patches_set(self, arg: bool, /) -> None: ...

    @property
    def print_info_messages_set(self) -> bool: ...

    @print_info_messages_set.setter
    def print_info_messages_set(self, arg: bool, /) -> None: ...

    @property
    def max_fringe_order_set(self) -> bool: ...

    @max_fringe_order_set.setter
    def max_fringe_order_set(self, arg: bool, /) -> None: ...

    @property
    def exact_model_set(self) -> bool: ...

    @exact_model_set.setter
    def exact_model_set(self, arg: bool, /) -> None: ...

    @property
    def exact_misalign_set(self) -> bool: ...

    @exact_misalign_set.setter
    def exact_misalign_set(self, arg: bool, /) -> None: ...

    @property
    def vertical_kick_set(self) -> bool: ...

    @vertical_kick_set.setter
    def vertical_kick_set(self, arg: bool, /) -> None: ...

    @property
    def cut_factor_set(self) -> bool: ...

    @cut_factor_set.setter
    def cut_factor_set(self, arg: bool, /) -> None: ...

    @property
    def translate_patch_drift_time_set(self) -> bool: ...

    @translate_patch_drift_time_set.setter
    def translate_patch_drift_time_set(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ExtraParsingInfoStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ExtraParsingInfoStruct: ...

    def __eq__(self, arg: ExtraParsingInfoStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class Fibre:
    """Fortran struct: fibre"""

    def __init__(self, DIR: int | None = None, pos: int | None = None, BETA0: float | None = None, GAMMA0I: float | None = None, GAMBET: float | None = None, MASS: float | None = None, CHARGE: float | None = None, AG: float | None = None, loc: int | None = None) -> None: ...

    @property
    def DIR(self) -> int | None: ...

    @DIR.setter
    def DIR(self, arg: int, /) -> None: ...

    @property
    def pos(self) -> int | None:
        """POSITION IN LAYOUT NEW STUFF...."""

    @pos.setter
    def pos(self, arg: int, /) -> None: ...

    @property
    def BETA0(self) -> float | None:
        """,P0C"""

    @BETA0.setter
    def BETA0(self, arg: float, /) -> None: ...

    @property
    def GAMMA0I(self) -> float | None:
        """,P0C"""

    @GAMMA0I.setter
    def GAMMA0I(self, arg: float, /) -> None: ...

    @property
    def GAMBET(self) -> float | None:
        """,P0C"""

    @GAMBET.setter
    def GAMBET(self, arg: float, /) -> None: ...

    @property
    def MASS(self) -> float | None:
        """,P0C"""

    @MASS.setter
    def MASS(self, arg: float, /) -> None: ...

    @property
    def CHARGE(self) -> float | None: ...

    @CHARGE.setter
    def CHARGE(self, arg: float, /) -> None: ...

    @property
    def AG(self) -> float | None:
        """spin g-2 TO TIE LAYOUTS"""

    @AG.setter
    def AG(self, arg: float, /) -> None: ...

    @property
    def loc(self) -> int | None: ...

    @loc.setter
    def loc(self, arg: int, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> Fibre: ...

    def __deepcopy__(self, arg: dict, /) -> Fibre: ...

    def __eq__(self, arg: Fibre, /) -> bool: ...

    def __hash__(self) -> int: ...

class Field1At2DPtStruct:
    """Fortran struct: field1_at_2D_pt_struct"""

    def __init__(self, f: float | None = None, df_dx: float | None = None, df_dy: float | None = None, d2f_dxdy: float | None = None) -> None: ...

    @property
    def f(self) -> float:
        """Field"""

    @f.setter
    def f(self, arg: float, /) -> None: ...

    @property
    def df_dx(self) -> float:
        """Normalized field 1st derivatives"""

    @df_dx.setter
    def df_dx(self, arg: float, /) -> None: ...

    @property
    def df_dy(self) -> float:
        """Normalized field 1st derivatives"""

    @df_dy.setter
    def df_dy(self, arg: float, /) -> None: ...

    @property
    def d2f_dxdy(self) -> float:
        """Normalized field 2nd derivative"""

    @d2f_dxdy.setter
    def d2f_dxdy(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> Field1At2DPtStruct: ...

    def __deepcopy__(self, arg: dict, /) -> Field1At2DPtStruct: ...

    def __eq__(self, arg: Field1At2DPtStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class Field1At3DPtStruct:
    """Fortran struct: field1_at_3D_pt_struct"""

    def __init__(self, f: float | None = None, df_dx: float | None = None, df_dy: float | None = None, df_dz: float | None = None, d2f_dxdy: float | None = None, d2f_dxdz: float | None = None, d2f_dydz: float | None = None, d3f_dxdydz: float | None = None) -> None: ...

    @property
    def f(self) -> float:
        """Field"""

    @f.setter
    def f(self, arg: float, /) -> None: ...

    @property
    def df_dx(self) -> float:
        """Normalized field 1st derivatives"""

    @df_dx.setter
    def df_dx(self, arg: float, /) -> None: ...

    @property
    def df_dy(self) -> float:
        """Normalized field 1st derivatives"""

    @df_dy.setter
    def df_dy(self, arg: float, /) -> None: ...

    @property
    def df_dz(self) -> float:
        """Normalized field 1st derivatives"""

    @df_dz.setter
    def df_dz(self, arg: float, /) -> None: ...

    @property
    def d2f_dxdy(self) -> float:
        """Normalized field 2nd derivatives"""

    @d2f_dxdy.setter
    def d2f_dxdy(self, arg: float, /) -> None: ...

    @property
    def d2f_dxdz(self) -> float:
        """Normalized field 2nd derivatives"""

    @d2f_dxdz.setter
    def d2f_dxdz(self, arg: float, /) -> None: ...

    @property
    def d2f_dydz(self) -> float:
        """Normalized field 2nd derivatives"""

    @d2f_dydz.setter
    def d2f_dydz(self, arg: float, /) -> None: ...

    @property
    def d3f_dxdydz(self) -> float:
        """Normalized field 3rd derivative"""

    @d3f_dxdydz.setter
    def d3f_dxdydz(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> Field1At3DPtStruct: ...

    def __deepcopy__(self, arg: dict, /) -> Field1At3DPtStruct: ...

    def __eq__(self, arg: Field1At3DPtStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class FieldAt2DBoxStruct:
    """Fortran struct: field_at_2D_box_struct"""

    def __init__(self, i_box: Sequence[int] | None = None) -> None: ...

    @property
    def pt(self) -> Field1At2DPtStructArray2D: ...

    @property
    def i_box(self) -> IntArray1D:
        """index at lower box corner."""

    @i_box.setter
    def i_box(self, arg: Sequence[int], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> FieldAt2DBoxStruct: ...

    def __deepcopy__(self, arg: dict, /) -> FieldAt2DBoxStruct: ...

    def __eq__(self, arg: FieldAt2DBoxStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class FieldAt3DBoxStruct:
    """Fortran struct: field_at_3D_box_struct"""

    def __init__(self, i_box: Sequence[int] | None = None) -> None: ...

    @property
    def pt(self) -> Field1At3DPtStructArray3D: ...

    @property
    def i_box(self) -> IntArray1D:
        """index at lower box corner."""

    @i_box.setter
    def i_box(self, arg: Sequence[int], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> FieldAt3DBoxStruct: ...

    def __deepcopy__(self, arg: dict, /) -> FieldAt3DBoxStruct: ...

    def __eq__(self, arg: FieldAt3DBoxStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class FloorPositionStruct:
    """Fortran struct: floor_position_struct"""

    def __init__(self, r: Sequence[float] | None = None, w: Sequence[Sequence[float]] | None = None, theta: float | None = None, phi: float | None = None, psi: float | None = None) -> None: ...

    @property
    def r(self) -> RealArray1D:
        """(x, y, z) offset from origin"""

    @r.setter
    def r(self, arg: Sequence[float], /) -> None: ...

    @property
    def w(self) -> RealArray2D:
        """W matrix. Columns are unit vectors of the frame axes."""

    @w.setter
    def w(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def theta(self) -> float:
        """angular orientation consistent with W matrix"""

    @theta.setter
    def theta(self, arg: float, /) -> None: ...

    @property
    def phi(self) -> float:
        """angular orientation consistent with W matrix"""

    @phi.setter
    def phi(self, arg: float, /) -> None: ...

    @property
    def psi(self) -> float:
        """angular orientation consistent with W matrix"""

    @psi.setter
    def psi(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> FloorPositionStruct: ...

    def __deepcopy__(self, arg: dict, /) -> FloorPositionStruct: ...

    def __eq__(self, arg: FloorPositionStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class FoilStruct:
    """Fortran struct: foil_struct"""

    def __init__(self) -> None: ...

    @property
    def material(self) -> MaterialStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> FoilStruct: ...

    def __deepcopy__(self, arg: dict, /) -> FoilStruct: ...

    def __eq__(self, arg: FoilStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class FringeFieldInfoStruct:
    """Fortran struct: fringe_field_info_struct"""

    def __init__(self, hard_ele: EleStruct | None = None, s_edge_hard: float | None = None, ds_edge: float | None = None, particle_at: int | None = None, hard_location: int | None = None, location: Sequence[int] | None = None, has_fringe: bool | None = None) -> None: ...

    @property
    def hard_ele(self) -> EleStruct | None: ...

    @hard_ele.setter
    def hard_ele(self, arg: EleStruct, /) -> None: ...

    @property
    def s_edge_hard(self) -> float: ...

    @s_edge_hard.setter
    def s_edge_hard(self, arg: float, /) -> None: ...

    @property
    def ds_edge(self) -> float:
        """Distance from particle to edge in hard_ele frame."""

    @ds_edge.setter
    def ds_edge(self, arg: float, /) -> None: ...

    @property
    def particle_at(self) -> int:
        """first_track_edge$, second_track_edge$, or none$"""

    @particle_at.setter
    def particle_at(self, arg: int, /) -> None: ...

    @property
    def hard_location(self) -> int | None:
        """Particle location wrt hard_ele. Points to element in location(:)."""

    @hard_location.setter
    def hard_location(self, arg: int, /) -> None: ...

    @property
    def location(self) -> IntAlloc1D:
        """
        Particle location in an element. entrance_end$, inside$, or exit_end$ Elements in list are the tracking element or its lords.
        """

    @location.setter
    def location(self, arg: Sequence[int], /) -> None: ...

    @property
    def has_fringe(self) -> bool:
        """Has a fringe to worry about?"""

    @has_fringe.setter
    def has_fringe(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> FringeFieldInfoStruct: ...

    def __deepcopy__(self, arg: dict, /) -> FringeFieldInfoStruct: ...

    def __eq__(self, arg: FringeFieldInfoStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class GenGrad1Struct:
    """Fortran struct: gen_grad1_struct"""

    def __init__(self, m: int | None = None, sincos: int | None = None, n_deriv_max: int | None = None, deriv: Sequence[Sequence[float]] | None = None) -> None: ...

    @property
    def m(self) -> int:
        """Azimuthal index"""

    @m.setter
    def m(self, arg: int, /) -> None: ...

    @property
    def sincos(self) -> int:
        """sin$ or cos$"""

    @sincos.setter
    def sincos(self, arg: int, /) -> None: ...

    @property
    def n_deriv_max(self) -> int:
        """
        Max GG derivative The derivative matrix is extended to include the interpolating spline polynomial.
        """

    @n_deriv_max.setter
    def n_deriv_max(self, arg: int, /) -> None: ...

    @property
    def deriv(self) -> RealArray2D:
        """Range: (iz0:iz1, 0:2*n_deriv_max+1)"""

    @deriv.setter
    def deriv(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> GenGrad1StructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> GenGrad1StructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> GenGrad1Struct: ...

    def __deepcopy__(self, arg: dict, /) -> GenGrad1Struct: ...

    def __eq__(self, arg: GenGrad1Struct, /) -> bool: ...

    def __hash__(self) -> int: ...

class GenGradMapStruct:
    """Fortran struct: gen_grad_map_struct"""

    def __init__(self, file: str | None = None, ele_anchor_pt: int | None = None, field_type: int | None = None, iz0: int | None = None, iz1: int | None = None, dz: float | None = None, r0: Sequence[float] | None = None, field_scale: float | None = None, master_parameter: int | None = None, curved_ref_frame: bool | None = None) -> None: ...

    @property
    def file(self) -> str:
        """Input file name. Used also as ID for instances."""

    @file.setter
    def file(self, arg: str, /) -> None: ...

    @property
    def gg(self) -> GenGrad1StructAlloc1D: ...

    @property
    def ele_anchor_pt(self) -> int:
        """anchor_beginning$, anchor_center$, or anchor_end$"""

    @ele_anchor_pt.setter
    def ele_anchor_pt(self, arg: int, /) -> None: ...

    @property
    def field_type(self) -> int:
        """or electric$"""

    @field_type.setter
    def field_type(self, arg: int, /) -> None: ...

    @property
    def iz0(self) -> int:
        """gg%deriv(iz0:iz1, :) lower bound."""

    @iz0.setter
    def iz0(self, arg: int, /) -> None: ...

    @property
    def iz1(self) -> int:
        """gg%deriv(iz0:iz1, :) upper bound."""

    @iz1.setter
    def iz1(self, arg: int, /) -> None: ...

    @property
    def dz(self) -> float:
        """Point spacing."""

    @dz.setter
    def dz(self, arg: float, /) -> None: ...

    @property
    def r0(self) -> RealArray1D:
        """field origin relative to ele_anchor_pt."""

    @r0.setter
    def r0(self, arg: Sequence[float], /) -> None: ...

    @property
    def field_scale(self) -> float:
        """Factor to scale the fields by"""

    @field_scale.setter
    def field_scale(self, arg: float, /) -> None: ...

    @property
    def master_parameter(self) -> int:
        """Master parameter in ele%value(:) array to use for scaling the field."""

    @master_parameter.setter
    def master_parameter(self, arg: int, /) -> None: ...

    @property
    def curved_ref_frame(self) -> bool: ...

    @curved_ref_frame.setter
    def curved_ref_frame(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> GenGradMapStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> GenGradMapStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> GenGradMapStruct: ...

    def __deepcopy__(self, arg: dict, /) -> GenGradMapStruct: ...

    def __eq__(self, arg: GenGradMapStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class GeneralBinStruct:
    """Fortran struct: general_bin_struct"""

    def __init__(self, count: Sequence[float] | None = None, min: Sequence[float] | None = None, max: Sequence[float] | None = None, delta: Sequence[float] | None = None, dim: int | None = None, n: Sequence[int] | None = None) -> None: ...

    @property
    def count(self) -> RealAlloc1D:
        """Counts (or weight) in each bin"""

    @count.setter
    def count(self, arg: Sequence[float], /) -> None: ...

    @property
    def min(self) -> RealArray1D:
        """Bounds for the bins"""

    @min.setter
    def min(self, arg: Sequence[float], /) -> None: ...

    @property
    def max(self) -> RealArray1D: ...

    @max.setter
    def max(self, arg: Sequence[float], /) -> None: ...

    @property
    def delta(self) -> RealArray1D:
        """Size of a bin"""

    @delta.setter
    def delta(self, arg: Sequence[float], /) -> None: ...

    @property
    def dim(self) -> int:
        """Number of dimensions"""

    @dim.setter
    def dim(self, arg: int, /) -> None: ...

    @property
    def n(self) -> IntArray1D:
        """number of bins in each dimension"""

    @n.setter
    def n(self, arg: Sequence[int], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> GeneralBinStruct: ...

    def __deepcopy__(self, arg: dict, /) -> GeneralBinStruct: ...

    def __eq__(self, arg: GeneralBinStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class GgTaylorStruct:
    """Fortran struct: gg_taylor_struct"""

    def __init__(self, ref: float | None = None) -> None: ...

    @property
    def ref(self) -> float: ...

    @ref.setter
    def ref(self, arg: float, /) -> None: ...

    @property
    def term(self) -> GgTaylorTermStructAlloc1D: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> GgTaylorStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> GgTaylorStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> GgTaylorStruct: ...

    def __deepcopy__(self, arg: dict, /) -> GgTaylorStruct: ...

    def __eq__(self, arg: GgTaylorStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class GgTaylorTermStruct:
    """Fortran struct: gg_taylor_term_struct"""

    def __init__(self, coef: float | None = None, expn: Sequence[int] | None = None) -> None: ...

    @property
    def coef(self) -> float: ...

    @coef.setter
    def coef(self, arg: float, /) -> None: ...

    @property
    def expn(self) -> IntArray1D: ...

    @expn.setter
    def expn(self, arg: Sequence[int], /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> GgTaylorTermStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> GgTaylorTermStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> GgTaylorTermStruct: ...

    def __deepcopy__(self, arg: dict, /) -> GgTaylorTermStruct: ...

    def __eq__(self, arg: GgTaylorTermStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class GptLatParamStruct:
    """Fortran struct: gpt_lat_param_struct"""

    def __init__(self, fieldmap_dimension: int | None = None, only_write_autophase_parameters: bool | None = None, gpt_filename: str | None = None, header_file_name: str | None = None, tracking_end_element: str | None = None) -> None: ...

    @property
    def fieldmap_dimension(self) -> int:
        """Dimensions for field map. 1 or 3"""

    @fieldmap_dimension.setter
    def fieldmap_dimension(self, arg: int, /) -> None: ...

    @property
    def only_write_autophase_parameters(self) -> bool:
        """Option to only write phasing info"""

    @only_write_autophase_parameters.setter
    def only_write_autophase_parameters(self, arg: bool, /) -> None: ...

    @property
    def gpt_filename(self) -> str:
        """Blank => Append '.gpt' to Bmad lattice file name."""

    @gpt_filename.setter
    def gpt_filename(self, arg: str, /) -> None: ...

    @property
    def header_file_name(self) -> str:
        """Header file to include in gpt file."""

    @header_file_name.setter
    def header_file_name(self, arg: str, /) -> None: ...

    @property
    def tracking_end_element(self) -> str:
        """Bmad lattice element name or index."""

    @tracking_end_element.setter
    def tracking_end_element(self, arg: str, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> GptLatParamStruct: ...

    def __deepcopy__(self, arg: dict, /) -> GptLatParamStruct: ...

    def __eq__(self, arg: GptLatParamStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class GridBeamInitStruct:
    """Fortran struct: grid_beam_init_struct"""

    def __init__(self, n_x: int | None = None, n_px: int | None = None, x_min: float | None = None, x_max: float | None = None, px_min: float | None = None, px_max: float | None = None) -> None: ...

    @property
    def n_x(self) -> int:
        """Number of columns."""

    @n_x.setter
    def n_x(self, arg: int, /) -> None: ...

    @property
    def n_px(self) -> int:
        """Number of rows."""

    @n_px.setter
    def n_px(self, arg: int, /) -> None: ...

    @property
    def x_min(self) -> float:
        """Lower x limit."""

    @x_min.setter
    def x_min(self, arg: float, /) -> None: ...

    @property
    def x_max(self) -> float:
        """Upper x limit."""

    @x_max.setter
    def x_max(self, arg: float, /) -> None: ...

    @property
    def px_min(self) -> float:
        """Lower px limit."""

    @px_min.setter
    def px_min(self, arg: float, /) -> None: ...

    @property
    def px_max(self) -> float:
        """Upper px limit."""

    @px_max.setter
    def px_max(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> GridBeamInitStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> GridBeamInitStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> GridBeamInitStruct: ...

    def __deepcopy__(self, arg: dict, /) -> GridBeamInitStruct: ...

    def __eq__(self, arg: GridBeamInitStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class GridFieldPt1Struct:
    """Fortran struct: grid_field_pt1_struct"""

    def __init__(self, E: Sequence[complex] | None = None, B: Sequence[complex] | None = None) -> None: ...

    @property
    def E(self) -> ComplexArray1D: ...

    @E.setter
    def E(self, arg: Sequence[complex], /) -> None: ...

    @property
    def B(self) -> ComplexArray1D: ...

    @B.setter
    def B(self, arg: Sequence[complex], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> GridFieldPt1Struct: ...

    def __deepcopy__(self, arg: dict, /) -> GridFieldPt1Struct: ...

    def __eq__(self, arg: GridFieldPt1Struct, /) -> bool: ...

    def __hash__(self) -> int: ...

class GridFieldPtStruct:
    """Fortran struct: grid_field_pt_struct"""

    def __init__(self, file: str | None = None, n_link: int | None = None) -> None: ...

    @property
    def file(self) -> str:
        """Input file name. Used also as ID for instances."""

    @file.setter
    def file(self, arg: str, /) -> None: ...

    @property
    def n_link(self) -> int:
        """For memory management of this structure"""

    @n_link.setter
    def n_link(self, arg: int, /) -> None: ...

    @property
    def pt(self) -> GridFieldPt1StructArray3D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> GridFieldPtStruct: ...

    def __deepcopy__(self, arg: dict, /) -> GridFieldPtStruct: ...

    def __eq__(self, arg: GridFieldPtStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class GridFieldStruct:
    """Fortran struct: grid_field_struct"""

    def __init__(self, geometry: int | None = None, harmonic: int | None = None, phi0_fieldmap: float | None = None, field_scale: float | None = None, field_type: int | None = None, master_parameter: int | None = None, ele_anchor_pt: int | None = None, interpolation_order: int | None = None, dr: Sequence[float] | None = None, r0: Sequence[float] | None = None, curved_ref_frame: bool | None = None, ptr: GridFieldPtStruct | None = None) -> None: ...

    @property
    def geometry(self) -> int:
        """Type of grid: xyz$, or rotationally_symmetric_rz$"""

    @geometry.setter
    def geometry(self, arg: int, /) -> None: ...

    @property
    def harmonic(self) -> int:
        """Harmonic of fundamental for AC fields."""

    @harmonic.setter
    def harmonic(self, arg: int, /) -> None: ...

    @property
    def phi0_fieldmap(self) -> float:
        """Mode oscillates as: twopi * (f * t + phi0_fieldmap)"""

    @phi0_fieldmap.setter
    def phi0_fieldmap(self, arg: float, /) -> None: ...

    @property
    def field_scale(self) -> float:
        """Factor to scale the fields by"""

    @field_scale.setter
    def field_scale(self, arg: float, /) -> None: ...

    @property
    def field_type(self) -> int:
        """or magnetic$ or electric$"""

    @field_type.setter
    def field_type(self, arg: int, /) -> None: ...

    @property
    def master_parameter(self) -> int:
        """Master parameter in ele%value(:) array to use for scaling the field."""

    @master_parameter.setter
    def master_parameter(self, arg: int, /) -> None: ...

    @property
    def ele_anchor_pt(self) -> int:
        """anchor_beginning$, anchor_center$, or anchor_end$"""

    @ele_anchor_pt.setter
    def ele_anchor_pt(self, arg: int, /) -> None: ...

    @property
    def interpolation_order(self) -> int:
        """Possibilities are 1 or 3."""

    @interpolation_order.setter
    def interpolation_order(self, arg: int, /) -> None: ...

    @property
    def dr(self) -> RealArray1D:
        """Grid spacing."""

    @dr.setter
    def dr(self, arg: Sequence[float], /) -> None: ...

    @property
    def r0(self) -> RealArray1D:
        """Field origin relative to ele_anchor_pt."""

    @r0.setter
    def r0(self, arg: Sequence[float], /) -> None: ...

    @property
    def curved_ref_frame(self) -> bool: ...

    @curved_ref_frame.setter
    def curved_ref_frame(self, arg: bool, /) -> None: ...

    @property
    def ptr(self) -> GridFieldPtStruct | None: ...

    @ptr.setter
    def ptr(self, arg: GridFieldPtStruct, /) -> None: ...

    @property
    def bi_coef(self) -> BicubicCmplxCoefStructArray3D:
        """Save computed coefs for faster tracking"""

    @property
    def tri_coef(self) -> TricubicCmplxCoefStructArray3D:
        """Save computed coefs for faster tracking"""

    @staticmethod
    def new_array1d(sz: int = 0) -> GridFieldStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> GridFieldStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> GridFieldStruct: ...

    def __deepcopy__(self, arg: dict, /) -> GridFieldStruct: ...

    def __eq__(self, arg: GridFieldStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class HighEnergySpaceChargeStruct:
    """Fortran struct: high_energy_space_charge_struct"""

    def __init__(self, closed_orb: CoordStruct | None = None, kick_const: float | None = None, sig_x: float | None = None, sig_y: float | None = None, phi: float | None = None, sin_phi: float | None = None, cos_phi: float | None = None, sig_z: float | None = None) -> None: ...

    @property
    def closed_orb(self) -> CoordStruct:
        """beam orbit"""

    @closed_orb.setter
    def closed_orb(self, arg: CoordStruct, /) -> None: ...

    @property
    def kick_const(self) -> float: ...

    @kick_const.setter
    def kick_const(self, arg: float, /) -> None: ...

    @property
    def sig_x(self) -> float: ...

    @sig_x.setter
    def sig_x(self, arg: float, /) -> None: ...

    @property
    def sig_y(self) -> float: ...

    @sig_y.setter
    def sig_y(self, arg: float, /) -> None: ...

    @property
    def phi(self) -> float:
        """Rotation angle to go from lab frame to rotated frame."""

    @phi.setter
    def phi(self, arg: float, /) -> None: ...

    @property
    def sin_phi(self) -> float: ...

    @sin_phi.setter
    def sin_phi(self, arg: float, /) -> None: ...

    @property
    def cos_phi(self) -> float: ...

    @cos_phi.setter
    def cos_phi(self, arg: float, /) -> None: ...

    @property
    def sig_z(self) -> float: ...

    @sig_z.setter
    def sig_z(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> HighEnergySpaceChargeStruct: ...

    def __deepcopy__(self, arg: dict, /) -> HighEnergySpaceChargeStruct: ...

    def __eq__(self, arg: HighEnergySpaceChargeStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class IbsLifetimeStruct:
    """Fortran struct: ibs_lifetime_struct"""

    def __init__(self, Tlx: float | None = None, Tly: float | None = None, Tlp: float | None = None) -> None: ...

    @property
    def Tlx(self) -> float: ...

    @Tlx.setter
    def Tlx(self, arg: float, /) -> None: ...

    @property
    def Tly(self) -> float: ...

    @Tly.setter
    def Tly(self, arg: float, /) -> None: ...

    @property
    def Tlp(self) -> float: ...

    @Tlp.setter
    def Tlp(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> IbsLifetimeStruct: ...

    def __deepcopy__(self, arg: dict, /) -> IbsLifetimeStruct: ...

    def __eq__(self, arg: IbsLifetimeStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class IbsMaxratioStruct:
    """Fortran struct: ibs_maxratio_struct"""

    def __init__(self, rx: float | None = None, ry: float | None = None, r_p: float | None = None) -> None: ...

    @property
    def rx(self) -> float: ...

    @rx.setter
    def rx(self, arg: float, /) -> None: ...

    @property
    def ry(self) -> float: ...

    @ry.setter
    def ry(self, arg: float, /) -> None: ...

    @property
    def r_p(self) -> float: ...

    @r_p.setter
    def r_p(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> IbsMaxratioStruct: ...

    def __deepcopy__(self, arg: dict, /) -> IbsMaxratioStruct: ...

    def __eq__(self, arg: IbsMaxratioStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class IbsSimParamStruct:
    """Fortran struct: ibs_sim_param_struct"""

    def __init__(self, tau_a: float | None = None, clog_to_use: int | None = None, set_dispersion: bool | None = None, eta_set: float | None = None, etap_set: float | None = None, do_pwd: bool | None = None, inductance: float | None = None, formula: str | None = None) -> None: ...

    @property
    def tau_a(self) -> float:
        """horizontal damping rate (needed for coulomb log tail cut)"""

    @tau_a.setter
    def tau_a(self, arg: float, /) -> None: ...

    @property
    def clog_to_use(self) -> int:
        """
        see multi_coulomb_log subroutine for valid settings.  Set to 1 to disable tail-cut.  Set to 1 for linacs.
        """

    @clog_to_use.setter
    def clog_to_use(self, arg: int, /) -> None: ...

    @property
    def set_dispersion(self) -> bool:
        """
        True: add vertical dispersion to transfer matrix.  Valid for kubo method.
        """

    @set_dispersion.setter
    def set_dispersion(self, arg: bool, /) -> None: ...

    @property
    def eta_set(self) -> float:
        """
        If set_dispersion, then this value is used to add y-z coupling to the transfer matrix.
        """

    @eta_set.setter
    def eta_set(self, arg: float, /) -> None: ...

    @property
    def etap_set(self) -> float:
        """
        If set_dispersion, then this value is used to add y-z coupling to the transfer matrix.
        """

    @etap_set.setter
    def etap_set(self, arg: float, /) -> None: ...

    @property
    def do_pwd(self) -> bool:
        """
        If true, then use potential well distortion to calculate bunch lengths.  If false, bunch length is proportional to energy spread.
        """

    @do_pwd.setter
    def do_pwd(self, arg: bool, /) -> None: ...

    @property
    def inductance(self) -> float:
        """Inductive part of impedance for pwd calc."""

    @inductance.setter
    def inductance(self, arg: float, /) -> None: ...

    @property
    def formula(self) -> str:
        """
        Which IBS formulation to use.  See subroutine ibs1 for a list. real(rp) :: fake_3HC = -1   ! If greater than zero, divide growth rates by this factor.
        """

    @formula.setter
    def formula(self, arg: str, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> IbsSimParamStruct: ...

    def __deepcopy__(self, arg: dict, /) -> IbsSimParamStruct: ...

    def __eq__(self, arg: IbsSimParamStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class IbsStruct:
    """Fortran struct: ibs_struct"""

    def __init__(self, inv_Ta: float | None = None, inv_Tb: float | None = None, inv_Tz: float | None = None) -> None: ...

    @property
    def inv_Ta(self) -> float: ...

    @inv_Ta.setter
    def inv_Ta(self, arg: float, /) -> None: ...

    @property
    def inv_Tb(self) -> float: ...

    @inv_Tb.setter
    def inv_Tb(self, arg: float, /) -> None: ...

    @property
    def inv_Tz(self) -> float: ...

    @inv_Tz.setter
    def inv_Tz(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> IbsStruct: ...

    def __deepcopy__(self, arg: dict, /) -> IbsStruct: ...

    def __eq__(self, arg: IbsStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class Interval1CoefStruct:
    """Fortran struct: interval1_coef_struct"""

    def __init__(self, c0: float | None = None, c1: float | None = None, n_exp: float | None = None) -> None: ...

    @property
    def c0(self) -> float: ...

    @c0.setter
    def c0(self, arg: float, /) -> None: ...

    @property
    def c1(self) -> float: ...

    @c1.setter
    def c1(self, arg: float, /) -> None: ...

    @property
    def n_exp(self) -> float: ...

    @n_exp.setter
    def n_exp(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> Interval1CoefStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> Interval1CoefStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> Interval1CoefStruct: ...

    def __deepcopy__(self, arg: dict, /) -> Interval1CoefStruct: ...

    def __eq__(self, arg: Interval1CoefStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class KvBeamInitStruct:
    """Fortran struct: kv_beam_init_struct"""

    def __init__(self, part_per_phi: Sequence[int] | None = None, n_I2: int | None = None, A: float | None = None) -> None: ...

    @property
    def part_per_phi(self) -> IntArray1D:
        """number of particles per angle variable."""

    @part_per_phi.setter
    def part_per_phi(self, arg: Sequence[int], /) -> None: ...

    @property
    def n_I2(self) -> int:
        """number of I2"""

    @n_I2.setter
    def n_I2(self, arg: int, /) -> None: ...

    @property
    def A(self) -> float:
        """A = I1/e"""

    @A.setter
    def A(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> KvBeamInitStruct: ...

    def __deepcopy__(self, arg: dict, /) -> KvBeamInitStruct: ...

    def __eq__(self, arg: KvBeamInitStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class LatEleLocStruct:
    """Fortran struct: lat_ele_loc_struct"""

    def __init__(self, ix_ele: int | None = None, ix_branch: int | None = None) -> None: ...

    @property
    def ix_ele(self) -> int: ...

    @ix_ele.setter
    def ix_ele(self, arg: int, /) -> None: ...

    @property
    def ix_branch(self) -> int: ...

    @ix_branch.setter
    def ix_branch(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> LatEleLocStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> LatEleLocStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> LatEleLocStruct: ...

    def __deepcopy__(self, arg: dict, /) -> LatEleLocStruct: ...

    def __eq__(self, arg: LatEleLocStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class LatEleOrder1Struct:
    """Fortran struct: lat_ele_order1_struct"""

    def __init__(self, ix_branch: int | None = None, ix_order: int | None = None) -> None: ...

    @property
    def ix_branch(self) -> int:
        """Branch index"""

    @ix_branch.setter
    def ix_branch(self, arg: int, /) -> None: ...

    @property
    def ix_order(self) -> int:
        """Order index. -1 -> Unique in lattice, 0 -> unique in branch."""

    @ix_order.setter
    def ix_order(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> LatEleOrder1StructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> LatEleOrder1StructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> LatEleOrder1Struct: ...

    def __deepcopy__(self, arg: dict, /) -> LatEleOrder1Struct: ...

    def __eq__(self, arg: LatEleOrder1Struct, /) -> bool: ...

    def __hash__(self) -> int: ...

class LatEleOrderArrayStruct:
    """Fortran struct: lat_ele_order_array_struct"""

    def __init__(self) -> None: ...

    @property
    def ele(self) -> LatEleOrder1StructAlloc1D: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> LatEleOrderArrayStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> LatEleOrderArrayStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> LatEleOrderArrayStruct: ...

    def __deepcopy__(self, arg: dict, /) -> LatEleOrderArrayStruct: ...

    def __eq__(self, arg: LatEleOrderArrayStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class LatEleOrderStruct:
    """Fortran struct: lat_ele_order_struct"""

    def __init__(self) -> None: ...

    @property
    def branch(self) -> LatEleOrderArrayStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> LatEleOrderStruct: ...

    def __deepcopy__(self, arg: dict, /) -> LatEleOrderStruct: ...

    def __eq__(self, arg: LatEleOrderStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class LatParamStruct:
    """Fortran struct: lat_param_struct"""

    def __init__(self, n_part: float | None = None, total_length: float | None = None, unstable_factor: float | None = None, t1_with_RF: Sequence[Sequence[float]] | None = None, t1_no_RF: Sequence[Sequence[float]] | None = None, spin_tune: float | None = None, particle: int | None = None, default_tracking_species: int | None = None, geometry: int | None = None, ixx: int | None = None, stable: bool | None = None, live_branch: bool | None = None, g1_integral: float | None = None, g2_integral: float | None = None, g3_integral: float | None = None, bookkeeping_state: BookkeepingStateStruct | None = None, beam_init: BeamInitStruct | None = None) -> None: ...

    @property
    def n_part(self) -> float:
        """Particles/bunch (for BeamBeam elements)."""

    @n_part.setter
    def n_part(self, arg: float, /) -> None: ...

    @property
    def total_length(self) -> float:
        """total_length of branch. Warning: branch may not start at s = 0."""

    @total_length.setter
    def total_length(self, arg: float, /) -> None: ...

    @property
    def unstable_factor(self) -> float:
        """
        If positive: Growth rate/turn if unstable in closed branches or |orbit-aperture|/aperture if particle hits wall. Zero otherwise.
        """

    @unstable_factor.setter
    def unstable_factor(self, arg: float, /) -> None: ...

    @property
    def t1_with_RF(self) -> RealArray2D:
        """Full 1-turn matrix with RF on."""

    @t1_with_RF.setter
    def t1_with_RF(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def t1_no_RF(self) -> RealArray2D:
        """Full 1-turn matrix with RF off."""

    @t1_no_RF.setter
    def t1_no_RF(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def spin_tune(self) -> float:
        """Closed orbit spin tune."""

    @spin_tune.setter
    def spin_tune(self, arg: float, /) -> None: ...

    @property
    def particle(self) -> int:
        """
        Reference particle: positron$, electron$, etc. Call lattice_bookkeeper if this is changed.
        """

    @particle.setter
    def particle(self, arg: int, /) -> None: ...

    @property
    def default_tracking_species(self) -> int:
        """Default particle type to use in tracking."""

    @default_tracking_species.setter
    def default_tracking_species(self, arg: int, /) -> None: ...

    @property
    def geometry(self) -> int:
        """open$ or closed$"""

    @geometry.setter
    def geometry(self, arg: int, /) -> None: ...

    @property
    def ixx(self) -> int:
        """Integer for general use"""

    @ixx.setter
    def ixx(self, arg: int, /) -> None: ...

    @property
    def stable(self) -> bool:
        """is closed lat stable?"""

    @stable.setter
    def stable(self, arg: bool, /) -> None: ...

    @property
    def live_branch(self) -> bool:
        """Should tracking be done on the branch?"""

    @live_branch.setter
    def live_branch(self, arg: bool, /) -> None: ...

    @property
    def g1_integral(self) -> float:
        """Approximate |g| (bending strength) integral of branch."""

    @g1_integral.setter
    def g1_integral(self, arg: float, /) -> None: ...

    @property
    def g2_integral(self) -> float:
        """Approximate g^2 integral of branch."""

    @g2_integral.setter
    def g2_integral(self, arg: float, /) -> None: ...

    @property
    def g3_integral(self) -> float:
        """Approximate g^2 integral of branch."""

    @g3_integral.setter
    def g3_integral(self, arg: float, /) -> None: ...

    @property
    def bookkeeping_state(self) -> BookkeepingStateStruct:
        """Overall status for the branch."""

    @bookkeeping_state.setter
    def bookkeeping_state(self, arg: BookkeepingStateStruct, /) -> None: ...

    @property
    def beam_init(self) -> BeamInitStruct:
        """For beam initialization."""

    @beam_init.setter
    def beam_init(self, arg: BeamInitStruct, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> LatParamStruct: ...

    def __deepcopy__(self, arg: dict, /) -> LatParamStruct: ...

    def __eq__(self, arg: LatParamStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class LatPointerStruct:
    """Fortran struct: lat_pointer_struct"""

    def __init__(self, lat: LatStruct | None = None) -> None: ...

    @property
    def lat(self) -> LatStruct | None: ...

    @lat.setter
    def lat(self, arg: LatStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> LatPointerStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> LatPointerStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> LatPointerStruct: ...

    def __deepcopy__(self, arg: dict, /) -> LatPointerStruct: ...

    def __eq__(self, arg: LatPointerStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class LatStruct:
    """Fortran struct: lat_struct"""

    def __init__(self, use_name: str | None = None, lattice: str | None = None, machine: str | None = None, input_file_name: str | None = None, title: str | None = None, a: ModeInfoStruct | None = None, b: ModeInfoStruct | None = None, z: ModeInfoStruct | None = None, param: LatParamStruct | None = None, lord_state: BookkeepingStateStruct | None = None, ele_init: EleStruct | None = None, particle_start: CoordStruct | None = None, beam_init: BeamInitStruct | None = None, pre_tracker: PreTrackerStruct | None = None, nametable: NametableStruct | None = None, custom: Sequence[float] | None = None, version: int | None = None, n_ele_track: int | None = None, n_ele_max: int | None = None, n_control_max: int | None = None, n_ic_max: int | None = None, input_taylor_order: int | None = None, ic: Sequence[int] | None = None, photon_type: int | None = None, creation_hash: int | None = None, ramper_slave_bookkeeping: int | None = None, parser_make_xfer_mats: bool | None = None) -> None: ...

    @property
    def use_name(self) -> str:
        """Name of lat given by USE statement"""

    @use_name.setter
    def use_name(self, arg: str, /) -> None: ...

    @property
    def lattice(self) -> str:
        """Lattice"""

    @lattice.setter
    def lattice(self, arg: str, /) -> None: ...

    @property
    def machine(self) -> str:
        """Name of the machine the lattice is for ('LHC', etc)."""

    @machine.setter
    def machine(self, arg: str, /) -> None: ...

    @property
    def input_file_name(self) -> str:
        """Name of the lattice input file"""

    @input_file_name.setter
    def input_file_name(self, arg: str, /) -> None: ...

    @property
    def title(self) -> str:
        """General title"""

    @title.setter
    def title(self, arg: str, /) -> None: ...

    @property
    def print_str(self) -> FCharArray1D:
        """Saved print statements."""

    @property
    def constant(self) -> ExpressionAtomStructAlloc1D:
        """Constants defined in the lattice"""

    @property
    def a(self) -> ModeInfoStruct | None:
        """Tunes (fractional part), etc."""

    @a.setter
    def a(self, arg: ModeInfoStruct, /) -> None: ...

    @property
    def b(self) -> ModeInfoStruct | None:
        """Tunes (fractional part), etc."""

    @b.setter
    def b(self, arg: ModeInfoStruct, /) -> None: ...

    @property
    def z(self) -> ModeInfoStruct | None:
        """Tunes (fractional part), etc."""

    @z.setter
    def z(self, arg: ModeInfoStruct, /) -> None: ...

    @property
    def param(self) -> LatParamStruct | None:
        """Parameters"""

    @param.setter
    def param(self, arg: LatParamStruct, /) -> None: ...

    @property
    def lord_state(self) -> BookkeepingStateStruct:
        """lord bookkeeping status."""

    @lord_state.setter
    def lord_state(self, arg: BookkeepingStateStruct, /) -> None: ...

    @property
    def ele_init(self) -> EleStruct:
        """For use by any program"""

    @ele_init.setter
    def ele_init(self, arg: EleStruct, /) -> None: ...

    @property
    def ele(self) -> EleStructArray1D:
        """Array of elements [=> branch(0)]."""

    @property
    def branch(self) -> BranchStructAlloc1D:
        """Branch(0:) array"""

    @property
    def control(self) -> ControlStructAlloc1D:
        """Control list"""

    @property
    def particle_start(self) -> CoordStruct | None:
        """Starting particle_coords."""

    @particle_start.setter
    def particle_start(self, arg: CoordStruct, /) -> None: ...

    @property
    def beam_init(self) -> BeamInitStruct:
        """Beam initialization."""

    @beam_init.setter
    def beam_init(self, arg: BeamInitStruct, /) -> None: ...

    @property
    def pre_tracker(self) -> PreTrackerStruct:
        """For OPAL/IMPACT-T"""

    @pre_tracker.setter
    def pre_tracker(self, arg: PreTrackerStruct, /) -> None: ...

    @property
    def nametable(self) -> NametableStruct:
        """For quick searching by element name."""

    @nametable.setter
    def nametable(self, arg: NametableStruct, /) -> None: ...

    @property
    def custom(self) -> RealAlloc1D:
        """Custom attributes."""

    @custom.setter
    def custom(self, arg: Sequence[float], /) -> None: ...

    @property
    def version(self) -> int:
        """Version number"""

    @version.setter
    def version(self, arg: int, /) -> None: ...

    @property
    def n_ele_track(self) -> int | None:
        """Number of lat elements to track through."""

    @n_ele_track.setter
    def n_ele_track(self, arg: int, /) -> None: ...

    @property
    def n_ele_max(self) -> int | None:
        """Index of last valid element in %ele(:) array"""

    @n_ele_max.setter
    def n_ele_max(self, arg: int, /) -> None: ...

    @property
    def n_control_max(self) -> int:
        """Last index used in control_array"""

    @n_control_max.setter
    def n_control_max(self, arg: int, /) -> None: ...

    @property
    def n_ic_max(self) -> int:
        """Last index used in ic_array"""

    @n_ic_max.setter
    def n_ic_max(self, arg: int, /) -> None: ...

    @property
    def input_taylor_order(self) -> int:
        """As set in the input file"""

    @input_taylor_order.setter
    def input_taylor_order(self, arg: int, /) -> None: ...

    @property
    def ic(self) -> IntAlloc1D:
        """Index to %control(:) from slaves."""

    @ic.setter
    def ic(self, arg: Sequence[int], /) -> None: ...

    @property
    def photon_type(self) -> int:
        """Or coherent$. For X-ray simulations."""

    @photon_type.setter
    def photon_type(self, arg: int, /) -> None: ...

    @property
    def creation_hash(self) -> int:
        """
        Set by bmad_parser. creation_hash will vary if any of the lattice files are modified.
        """

    @creation_hash.setter
    def creation_hash(self, arg: int, /) -> None: ...

    @property
    def ramper_slave_bookkeeping(self) -> int: ...

    @ramper_slave_bookkeeping.setter
    def ramper_slave_bookkeeping(self, arg: int, /) -> None: ...

    @property
    def parser_make_xfer_mats(self) -> bool:
        """Is Bmad parser to make element transfer matrices?"""

    @parser_make_xfer_mats.setter
    def parser_make_xfer_mats(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> LatStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> LatStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> LatStruct: ...

    def __deepcopy__(self, arg: dict, /) -> LatStruct: ...

    def __eq__(self, arg: LatStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class Layout:
    """Fortran struct: layout"""

    def __init__(self, NAME: str | None = None, INDEX: int | None = None, HARMONIC_NUMBER: float | None = None, CLOSED: bool | None = None, N: int | None = None, NTHIN: int | None = None, THIN: float | None = None, LASTPOS: int | None = None) -> None: ...

    @property
    def NAME(self) -> str:
        """IDENTIFICATION"""

    @NAME.setter
    def NAME(self, arg: str, /) -> None: ...

    @property
    def INDEX(self) -> int | None:
        """IDENTIFICATION, CHARGE SIGN"""

    @INDEX.setter
    def INDEX(self, arg: int, /) -> None: ...

    @property
    def HARMONIC_NUMBER(self) -> float | None: ...

    @HARMONIC_NUMBER.setter
    def HARMONIC_NUMBER(self, arg: float, /) -> None: ...

    @property
    def CLOSED(self) -> bool | None: ...

    @CLOSED.setter
    def CLOSED(self, arg: bool, /) -> None: ...

    @property
    def N(self) -> int | None:
        """TOTAL ELEMENT IN THE CHAIN"""

    @N.setter
    def N(self, arg: int, /) -> None: ...

    @property
    def NTHIN(self) -> int | None:
        """NUMBER IF THIN LENSES IN COLLECTION  (FOR SPEED ESTIMATES)"""

    @NTHIN.setter
    def NTHIN(self, arg: int, /) -> None: ...

    @property
    def THIN(self) -> float | None:
        """
        PARAMETER USED FOR AUTOMATIC CUTTING INTO THIN LENS POINTERS OF LINK LAYOUT
        """

    @THIN.setter
    def THIN(self, arg: float, /) -> None: ...

    @property
    def LASTPOS(self) -> int | None:
        """POSITION OF LAST VISITED"""

    @LASTPOS.setter
    def LASTPOS(self, arg: int, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> Layout: ...

    def __deepcopy__(self, arg: dict, /) -> Layout: ...

    def __eq__(self, arg: Layout, /) -> bool: ...

    def __hash__(self) -> int: ...

class LinacNormalModeStruct:
    """Fortran struct: linac_normal_mode_struct"""

    def __init__(self, i2_E4: float | None = None, i3_E7: float | None = None, i5a_E6: float | None = None, i5b_E6: float | None = None, sig_E1: float | None = None, a_emittance_end: float | None = None, b_emittance_end: float | None = None) -> None: ...

    @property
    def i2_E4(self) -> float:
        """Integral: g^2 * gamma^4"""

    @i2_E4.setter
    def i2_E4(self, arg: float, /) -> None: ...

    @property
    def i3_E7(self) -> float:
        """Integral: g^3 * gamma^7"""

    @i3_E7.setter
    def i3_E7(self, arg: float, /) -> None: ...

    @property
    def i5a_E6(self) -> float:
        """Integral: (g^3 * H_a) * gamma^6"""

    @i5a_E6.setter
    def i5a_E6(self, arg: float, /) -> None: ...

    @property
    def i5b_E6(self) -> float:
        """Integral: (g^3 * H_b) * gamma^6"""

    @i5b_E6.setter
    def i5b_E6(self, arg: float, /) -> None: ...

    @property
    def sig_E1(self) -> float:
        """Energy spread after 1 pass (eV)"""

    @sig_E1.setter
    def sig_E1(self, arg: float, /) -> None: ...

    @property
    def a_emittance_end(self) -> float:
        """a mode emittance at end of linac"""

    @a_emittance_end.setter
    def a_emittance_end(self, arg: float, /) -> None: ...

    @property
    def b_emittance_end(self) -> float:
        """b mode emittance at end of linac"""

    @b_emittance_end.setter
    def b_emittance_end(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> LinacNormalModeStruct: ...

    def __deepcopy__(self, arg: dict, /) -> LinacNormalModeStruct: ...

    def __eq__(self, arg: LinacNormalModeStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class LinearEleIsfStruct:
    """Fortran struct: linear_ele_isf_struct"""

    def __init__(self) -> None: ...

    @property
    def node(self) -> LinearIsf1StructAlloc1D:
        """Array per PTC integration node."""

    @staticmethod
    def new_array1d(sz: int = 0) -> LinearEleIsfStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> LinearEleIsfStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> LinearEleIsfStruct: ...

    def __deepcopy__(self, arg: dict, /) -> LinearEleIsfStruct: ...

    def __eq__(self, arg: LinearEleIsfStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class LinearIsf1Struct:
    """Fortran struct: linear_isf1_struct"""

    def __init__(self, orb0: Sequence[float] | None = None, isf: Sequence[Sequence[float]] | None = None, s: float | None = None) -> None: ...

    @property
    def orb0(self) -> RealArray1D:
        """Closed orbit."""

    @orb0.setter
    def orb0(self, arg: Sequence[float], /) -> None: ...

    @property
    def isf(self) -> RealArray2D:
        """Linear ISF map at a given point."""

    @isf.setter
    def isf(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def s(self) -> float:
        """
        Offset from beginning of element. !! real(rp) :: m_1turn(6,6) = 0   ! Orbital 1-turn matrix.
        """

    @s.setter
    def s(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> LinearIsf1StructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> LinearIsf1StructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> LinearIsf1Struct: ...

    def __deepcopy__(self, arg: dict, /) -> LinearIsf1Struct: ...

    def __eq__(self, arg: LinearIsf1Struct, /) -> bool: ...

    def __hash__(self) -> int: ...

class MadEnergyStruct:
    """Fortran struct: mad_energy_struct"""

    def __init__(self, total: float | None = None, beta: float | None = None, gamma: float | None = None, kinetic: float | None = None, p0c: float | None = None, particle: int | None = None) -> None: ...

    @property
    def total(self) -> float: ...

    @total.setter
    def total(self, arg: float, /) -> None: ...

    @property
    def beta(self) -> float:
        """normalized velocity: v/c"""

    @beta.setter
    def beta(self, arg: float, /) -> None: ...

    @property
    def gamma(self) -> float:
        """relativistic factor: 1/sqrt(1-beta^2)"""

    @gamma.setter
    def gamma(self, arg: float, /) -> None: ...

    @property
    def kinetic(self) -> float:
        """kinetic energy"""

    @kinetic.setter
    def kinetic(self, arg: float, /) -> None: ...

    @property
    def p0c(self) -> float:
        """particle momentum"""

    @p0c.setter
    def p0c(self, arg: float, /) -> None: ...

    @property
    def particle(self) -> int:
        """particle species"""

    @particle.setter
    def particle(self, arg: int, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> MadEnergyStruct: ...

    def __deepcopy__(self, arg: dict, /) -> MadEnergyStruct: ...

    def __eq__(self, arg: MadEnergyStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class MadMapStruct:
    """Fortran struct: mad_map_struct"""

    def __init__(self, k: Sequence[float] | None = None, r: Sequence[Sequence[float]] | None = None, t: Sequence[Sequence[Sequence[float]]] | None = None) -> None: ...

    @property
    def k(self) -> RealArray1D:
        """0th order map."""

    @k.setter
    def k(self, arg: Sequence[float], /) -> None: ...

    @property
    def r(self) -> RealArray2D:
        """1st order map."""

    @r.setter
    def r(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def t(self) -> RealArray3D:
        """2nd order map."""

    @t.setter
    def t(self, arg: Sequence[Sequence[Sequence[float]]], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> MadMapStruct: ...

    def __deepcopy__(self, arg: dict, /) -> MadMapStruct: ...

    def __eq__(self, arg: MadMapStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class MaterialStruct:
    """Fortran struct: material_struct"""

    def __init__(self, species: int | None = None, number: int | None = None, density: float | None = None, density_used: float | None = None, area_density: float | None = None, area_density_used: float | None = None, radiation_length: float | None = None, radiation_length_used: float | None = None) -> None: ...

    @property
    def species(self) -> int: ...

    @species.setter
    def species(self, arg: int, /) -> None: ...

    @property
    def number(self) -> int:
        """Relative number"""

    @number.setter
    def number(self, arg: int, /) -> None: ...

    @property
    def density(self) -> float: ...

    @density.setter
    def density(self, arg: float, /) -> None: ...

    @property
    def density_used(self) -> float: ...

    @density_used.setter
    def density_used(self, arg: float, /) -> None: ...

    @property
    def area_density(self) -> float: ...

    @area_density.setter
    def area_density(self, arg: float, /) -> None: ...

    @property
    def area_density_used(self) -> float: ...

    @area_density_used.setter
    def area_density_used(self, arg: float, /) -> None: ...

    @property
    def radiation_length(self) -> float: ...

    @radiation_length.setter
    def radiation_length(self, arg: float, /) -> None: ...

    @property
    def radiation_length_used(self) -> float: ...

    @radiation_length_used.setter
    def radiation_length_used(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> MaterialStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> MaterialStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> MaterialStruct: ...

    def __deepcopy__(self, arg: dict, /) -> MaterialStruct: ...

    def __eq__(self, arg: MaterialStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class Mesh3DStruct:
    """Fortran struct: mesh3d_struct"""

    def __init__(self, nlo: Sequence[int] | None = None, nhi: Sequence[int] | None = None, npad: Sequence[int] | None = None, min: Sequence[float] | None = None, max: Sequence[float] | None = None, delta: Sequence[float] | None = None, gamma: float | None = None, charge: float | None = None, rho: Sequence[Sequence[Sequence[float]]] | None = None, phi: Sequence[Sequence[Sequence[float]]] | None = None) -> None: ...

    @property
    def nlo(self) -> IntArray1D:
        """
        Lowest  grid index in x, y, z (m) of rho and the quantity being computed (phi or E)
        """

    @nlo.setter
    def nlo(self, arg: Sequence[int], /) -> None: ...

    @property
    def nhi(self) -> IntArray1D:
        """
        Highest grid index in x, y, z (m) of rho and the quantity being computed (phi or E)
        """

    @nhi.setter
    def nhi(self, arg: Sequence[int], /) -> None: ...

    @property
    def npad(self) -> IntArray1D:
        """Array padding for cyclic convolution"""

    @npad.setter
    def npad(self, arg: Sequence[int], /) -> None: ...

    @property
    def min(self) -> RealArray1D:
        """Minimim in each dimension"""

    @min.setter
    def min(self, arg: Sequence[float], /) -> None: ...

    @property
    def max(self) -> RealArray1D:
        """Maximum in each dimension"""

    @max.setter
    def max(self, arg: Sequence[float], /) -> None: ...

    @property
    def delta(self) -> RealArray1D:
        """Grid spacing"""

    @delta.setter
    def delta(self, arg: Sequence[float], /) -> None: ...

    @property
    def gamma(self) -> float:
        """Relativistic gamma"""

    @gamma.setter
    def gamma(self, arg: float, /) -> None: ...

    @property
    def charge(self) -> float:
        """Total charge on mesh"""

    @charge.setter
    def charge(self, arg: float, /) -> None: ...

    @property
    def rho(self) -> RealArray3D:
        """Charge density grid"""

    @rho.setter
    def rho(self, arg: Sequence[Sequence[Sequence[float]]], /) -> None: ...

    @property
    def phi(self) -> RealArray3D:
        """electric potential grid"""

    @phi.setter
    def phi(self, arg: Sequence[Sequence[Sequence[float]]], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> Mesh3DStruct: ...

    def __deepcopy__(self, arg: dict, /) -> Mesh3DStruct: ...

    def __eq__(self, arg: Mesh3DStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class Mode3Struct:
    """Fortran struct: mode3_struct"""

    def __init__(self, v: Sequence[Sequence[float]] | None = None, a: TwissStruct | None = None, b: TwissStruct | None = None, c: TwissStruct | None = None, x: TwissStruct | None = None, y: TwissStruct | None = None) -> None: ...

    @property
    def v(self) -> RealArray2D: ...

    @v.setter
    def v(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def a(self) -> TwissStruct: ...

    @a.setter
    def a(self, arg: TwissStruct, /) -> None: ...

    @property
    def b(self) -> TwissStruct: ...

    @b.setter
    def b(self, arg: TwissStruct, /) -> None: ...

    @property
    def c(self) -> TwissStruct: ...

    @c.setter
    def c(self, arg: TwissStruct, /) -> None: ...

    @property
    def x(self) -> TwissStruct: ...

    @x.setter
    def x(self, arg: TwissStruct, /) -> None: ...

    @property
    def y(self) -> TwissStruct: ...

    @y.setter
    def y(self, arg: TwissStruct, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> Mode3Struct: ...

    def __deepcopy__(self, arg: dict, /) -> Mode3Struct: ...

    def __eq__(self, arg: Mode3Struct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ModeInfoStruct:
    """Fortran struct: mode_info_struct"""

    def __init__(self, stable: bool | None = None, tune: float | None = None, emit: float | None = None, chrom: float | None = None, sigma: float | None = None, sigmap: float | None = None) -> None: ...

    @property
    def stable(self) -> bool:
        """Is the mode stable?"""

    @stable.setter
    def stable(self, arg: bool, /) -> None: ...

    @property
    def tune(self) -> float:
        """'fractional' tune in radians"""

    @tune.setter
    def tune(self, arg: float, /) -> None: ...

    @property
    def emit(self) -> float:
        """Emittance (unnormalized)."""

    @emit.setter
    def emit(self, arg: float, /) -> None: ...

    @property
    def chrom(self) -> float:
        """Chromaticity."""

    @chrom.setter
    def chrom(self, arg: float, /) -> None: ...

    @property
    def sigma(self) -> float:
        """Beam size."""

    @sigma.setter
    def sigma(self, arg: float, /) -> None: ...

    @property
    def sigmap(self) -> float:
        """Beam divergence."""

    @sigmap.setter
    def sigmap(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ModeInfoStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ModeInfoStruct: ...

    def __eq__(self, arg: ModeInfoStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class MolecularComponentStruct:
    """Fortran struct: molecular_component_struct"""

    def __init__(self, atom: str | None = None, number: int | None = None) -> None: ...

    @property
    def atom(self) -> str: ...

    @atom.setter
    def atom(self, arg: str, /) -> None: ...

    @property
    def number(self) -> int: ...

    @number.setter
    def number(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> MolecularComponentStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> MolecularComponentStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> MolecularComponentStruct: ...

    def __deepcopy__(self, arg: dict, /) -> MolecularComponentStruct: ...

    def __eq__(self, arg: MolecularComponentStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class MomentumApertureStruct:
    """Fortran struct: momentum_aperture_struct"""

    def __init__(self, s: float | None = None, pos: float | None = None, neg: float | None = None) -> None: ...

    @property
    def s(self) -> float: ...

    @s.setter
    def s(self, arg: float, /) -> None: ...

    @property
    def pos(self) -> float: ...

    @pos.setter
    def pos(self, arg: float, /) -> None: ...

    @property
    def neg(self) -> float: ...

    @neg.setter
    def neg(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> MomentumApertureStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> MomentumApertureStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> MomentumApertureStruct: ...

    def __deepcopy__(self, arg: dict, /) -> MomentumApertureStruct: ...

    def __eq__(self, arg: MomentumApertureStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class MultipassAllInfoStruct:
    """Fortran struct: multipass_all_info_struct"""

    def __init__(self) -> None: ...

    @property
    def lord(self) -> MultipassLordInfoStructAlloc1D:
        """Array of lords"""

    @property
    def branch(self) -> MultipassBranchInfoStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> MultipassAllInfoStruct: ...

    def __deepcopy__(self, arg: dict, /) -> MultipassAllInfoStruct: ...

    def __eq__(self, arg: MultipassAllInfoStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class MultipassBranchInfoStruct:
    """Fortran struct: multipass_branch_info_struct"""

    def __init__(self) -> None: ...

    @property
    def ele(self) -> MultipassEleInfoStructAlloc1D: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> MultipassBranchInfoStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> MultipassBranchInfoStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> MultipassBranchInfoStruct: ...

    def __deepcopy__(self, arg: dict, /) -> MultipassBranchInfoStruct: ...

    def __eq__(self, arg: MultipassBranchInfoStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class MultipassEleInfoStruct:
    """Fortran struct: multipass_ele_info_struct"""

    def __init__(self, multipass: bool | None = None, ix_pass: int | None = None, ix_lord: Sequence[int] | None = None, ix_super: Sequence[int] | None = None) -> None: ...

    @property
    def multipass(self) -> bool:
        """True if involved in multipass. False otherwise"""

    @multipass.setter
    def multipass(self, arg: bool, /) -> None: ...

    @property
    def ix_pass(self) -> int:
        """Pass number"""

    @ix_pass.setter
    def ix_pass(self, arg: int, /) -> None: ...

    @property
    def ix_lord(self) -> IntAlloc1D:
        """Pointers to lord(:) array"""

    @ix_lord.setter
    def ix_lord(self, arg: Sequence[int], /) -> None: ...

    @property
    def ix_super(self) -> IntAlloc1D:
        """Indexes to slave(ix_pass, super_slave%ix_ele) matrix"""

    @ix_super.setter
    def ix_super(self, arg: Sequence[int], /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> MultipassEleInfoStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> MultipassEleInfoStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> MultipassEleInfoStruct: ...

    def __deepcopy__(self, arg: dict, /) -> MultipassEleInfoStruct: ...

    def __eq__(self, arg: MultipassEleInfoStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class MultipassLordInfoStruct:
    """Fortran struct: multipass_lord_info_struct"""

    def __init__(self, lord: EleStruct | None = None, n_pass: int | None = None, n_super_slave: int | None = None) -> None: ...

    @property
    def lord(self) -> EleStruct | None:
        """Lord element"""

    @lord.setter
    def lord(self, arg: EleStruct, /) -> None: ...

    @property
    def n_pass(self) -> int:
        """Number of passes (= number of slaves)"""

    @n_pass.setter
    def n_pass(self, arg: int, /) -> None: ...

    @property
    def n_super_slave(self) -> int:
        """Number of super_slaves per super_lord."""

    @n_super_slave.setter
    def n_super_slave(self, arg: int, /) -> None: ...

    @property
    def super_lord(self) -> ElePointerStructAlloc1D:
        """Super_lord list if they exist."""

    @property
    def slave(self) -> ElePointerStructArray2D:
        """Slaves list in tracking part."""

    @staticmethod
    def new_array1d(sz: int = 0) -> MultipassLordInfoStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> MultipassLordInfoStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> MultipassLordInfoStruct: ...

    def __deepcopy__(self, arg: dict, /) -> MultipassLordInfoStruct: ...

    def __eq__(self, arg: MultipassLordInfoStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class MultipassRegionBranchStruct:
    """Fortran struct: multipass_region_branch_struct"""

    def __init__(self) -> None: ...

    @property
    def ele(self) -> MultipassRegionEleStructAlloc1D: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> MultipassRegionBranchStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> MultipassRegionBranchStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> MultipassRegionBranchStruct: ...

    def __deepcopy__(self, arg: dict, /) -> MultipassRegionBranchStruct: ...

    def __eq__(self, arg: MultipassRegionBranchStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class MultipassRegionEleStruct:
    """Fortran struct: multipass_region_ele_struct"""

    def __init__(self, ix_region: int | None = None, region_start_pt: bool | None = None, region_stop_pt: bool | None = None) -> None: ...

    @property
    def ix_region(self) -> int: ...

    @ix_region.setter
    def ix_region(self, arg: int, /) -> None: ...

    @property
    def region_start_pt(self) -> bool: ...

    @region_start_pt.setter
    def region_start_pt(self, arg: bool, /) -> None: ...

    @property
    def region_stop_pt(self) -> bool: ...

    @region_stop_pt.setter
    def region_stop_pt(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> MultipassRegionEleStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> MultipassRegionEleStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> MultipassRegionEleStruct: ...

    def __deepcopy__(self, arg: dict, /) -> MultipassRegionEleStruct: ...

    def __eq__(self, arg: MultipassRegionEleStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class MultipassRegionLatStruct:
    """Fortran struct: multipass_region_lat_struct"""

    def __init__(self) -> None: ...

    @property
    def branch(self) -> MultipassRegionBranchStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> MultipassRegionLatStruct: ...

    def __deepcopy__(self, arg: dict, /) -> MultipassRegionLatStruct: ...

    def __eq__(self, arg: MultipassRegionLatStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class MultipoleCacheStruct:
    """Fortran struct: multipole_cache_struct"""

    def __init__(self, a_pole_mag: Sequence[float] | None = None, b_pole_mag: Sequence[float] | None = None, a_kick_mag: Sequence[float] | None = None, b_kick_mag: Sequence[float] | None = None, ix_pole_mag_max: int | None = None, ix_kick_mag_max: int | None = None, mag_valid: bool | None = None, a_pole_elec: Sequence[float] | None = None, b_pole_elec: Sequence[float] | None = None, a_kick_elec: Sequence[float] | None = None, b_kick_elec: Sequence[float] | None = None, ix_pole_elec_max: int | None = None, ix_kick_elec_max: int | None = None, elec_valid: bool | None = None) -> None: ...

    @property
    def a_pole_mag(self) -> RealAlloc1D: ...

    @a_pole_mag.setter
    def a_pole_mag(self, arg: Sequence[float], /) -> None: ...

    @property
    def b_pole_mag(self) -> RealAlloc1D: ...

    @b_pole_mag.setter
    def b_pole_mag(self, arg: Sequence[float], /) -> None: ...

    @property
    def a_kick_mag(self) -> RealAlloc1D: ...

    @a_kick_mag.setter
    def a_kick_mag(self, arg: Sequence[float], /) -> None: ...

    @property
    def b_kick_mag(self) -> RealAlloc1D: ...

    @b_kick_mag.setter
    def b_kick_mag(self, arg: Sequence[float], /) -> None: ...

    @property
    def ix_pole_mag_max(self) -> int: ...

    @ix_pole_mag_max.setter
    def ix_pole_mag_max(self, arg: int, /) -> None: ...

    @property
    def ix_kick_mag_max(self) -> int: ...

    @ix_kick_mag_max.setter
    def ix_kick_mag_max(self, arg: int, /) -> None: ...

    @property
    def mag_valid(self) -> bool:
        """From elseparator hkick and vkick."""

    @mag_valid.setter
    def mag_valid(self, arg: bool, /) -> None: ...

    @property
    def a_pole_elec(self) -> RealAlloc1D: ...

    @a_pole_elec.setter
    def a_pole_elec(self, arg: Sequence[float], /) -> None: ...

    @property
    def b_pole_elec(self) -> RealAlloc1D: ...

    @b_pole_elec.setter
    def b_pole_elec(self, arg: Sequence[float], /) -> None: ...

    @property
    def a_kick_elec(self) -> RealAlloc1D: ...

    @a_kick_elec.setter
    def a_kick_elec(self, arg: Sequence[float], /) -> None: ...

    @property
    def b_kick_elec(self) -> RealAlloc1D: ...

    @b_kick_elec.setter
    def b_kick_elec(self, arg: Sequence[float], /) -> None: ...

    @property
    def ix_pole_elec_max(self) -> int: ...

    @ix_pole_elec_max.setter
    def ix_pole_elec_max(self, arg: int, /) -> None: ...

    @property
    def ix_kick_elec_max(self) -> int: ...

    @ix_kick_elec_max.setter
    def ix_kick_elec_max(self, arg: int, /) -> None: ...

    @property
    def elec_valid(self) -> bool: ...

    @elec_valid.setter
    def elec_valid(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> MultipoleCacheStruct: ...

    def __deepcopy__(self, arg: dict, /) -> MultipoleCacheStruct: ...

    def __eq__(self, arg: MultipoleCacheStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class NamedNumberStruct:
    """Fortran struct: named_number_struct"""

    def __init__(self, name: str | None = None, value: float | None = None) -> None: ...

    @property
    def name(self) -> str: ...

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def value(self) -> float: ...

    @value.setter
    def value(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> NamedNumberStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> NamedNumberStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> NamedNumberStruct: ...

    def __deepcopy__(self, arg: dict, /) -> NamedNumberStruct: ...

    def __eq__(self, arg: NamedNumberStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class NametableStruct:
    """Fortran struct: nametable_struct"""

    def __init__(self, index: Sequence[int] | None = None, n_min: int | None = None, n_max: int | None = None) -> None: ...

    @property
    def name(self) -> FCharArray1D:
        """Array of names."""

    @property
    def index(self) -> IntAlloc1D:
        """
        Sorted index for names(:) array. names(an_index(i)) is in alphabetical order.
        """

    @index.setter
    def index(self, arg: Sequence[int], /) -> None: ...

    @property
    def n_min(self) -> int:
        """Set to 0 for use in a lattice."""

    @n_min.setter
    def n_min(self, arg: int, /) -> None: ...

    @property
    def n_max(self) -> int:
        """Use only names(n_min:n_max) part of array."""

    @n_max.setter
    def n_max(self, arg: int, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> NametableStruct: ...

    def __deepcopy__(self, arg: dict, /) -> NametableStruct: ...

    def __eq__(self, arg: NametableStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class NormalModesStruct:
    """Fortran struct: normal_modes_struct"""

    def __init__(self, synch_int: Sequence[float] | None = None, sigE_E: float | None = None, sig_z: float | None = None, e_loss: float | None = None, rf_voltage: float | None = None, pz_aperture: float | None = None, pz_average: float | None = None, momentum_compaction: float | None = None, dpz_damp: float | None = None, a: AnormalModeStruct | None = None, b: AnormalModeStruct | None = None, z: AnormalModeStruct | None = None, lin: LinacNormalModeStruct | None = None) -> None: ...

    @property
    def synch_int(self) -> RealArray1D:
        """Synchrotron integrals I0, I1, I2, and I3"""

    @synch_int.setter
    def synch_int(self, arg: Sequence[float], /) -> None: ...

    @property
    def sigE_E(self) -> float:
        """SigmaE/E"""

    @sigE_E.setter
    def sigE_E(self, arg: float, /) -> None: ...

    @property
    def sig_z(self) -> float:
        """Sigma_Z"""

    @sig_z.setter
    def sig_z(self, arg: float, /) -> None: ...

    @property
    def e_loss(self) -> float:
        """Energy loss / turn (eV)"""

    @e_loss.setter
    def e_loss(self, arg: float, /) -> None: ...

    @property
    def rf_voltage(self) -> float:
        """Total rfcavity voltage (eV)"""

    @rf_voltage.setter
    def rf_voltage(self, arg: float, /) -> None: ...

    @property
    def pz_aperture(self) -> float:
        """pz aperture limit. Used with Touschek calculations."""

    @pz_aperture.setter
    def pz_aperture(self, arg: float, /) -> None: ...

    @property
    def pz_average(self) -> float:
        """Average over branch due to damping."""

    @pz_average.setter
    def pz_average(self, arg: float, /) -> None: ...

    @property
    def momentum_compaction(self) -> float: ...

    @momentum_compaction.setter
    def momentum_compaction(self, arg: float, /) -> None: ...

    @property
    def dpz_damp(self) -> float:
        """Change in pz without RF"""

    @dpz_damp.setter
    def dpz_damp(self, arg: float, /) -> None: ...

    @property
    def a(self) -> AnormalModeStruct: ...

    @a.setter
    def a(self, arg: AnormalModeStruct, /) -> None: ...

    @property
    def b(self) -> AnormalModeStruct: ...

    @b.setter
    def b(self, arg: AnormalModeStruct, /) -> None: ...

    @property
    def z(self) -> AnormalModeStruct: ...

    @z.setter
    def z(self, arg: AnormalModeStruct, /) -> None: ...

    @property
    def lin(self) -> LinacNormalModeStruct: ...

    @lin.setter
    def lin(self, arg: LinacNormalModeStruct, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> NormalModesStruct: ...

    def __deepcopy__(self, arg: dict, /) -> NormalModesStruct: ...

    def __eq__(self, arg: NormalModesStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class OutIoOutputDirectStruct:
    """Fortran struct: out_io_output_direct_struct"""

    def __init__(self, print_and_capture: Sequence[bool] | None = None, file_unit: Sequence[int] | None = None) -> None: ...

    @property
    def print_and_capture(self) -> BoolArray1D: ...

    @print_and_capture.setter
    def print_and_capture(self, arg: Sequence[bool], /) -> None: ...

    @property
    def file_unit(self) -> IntArray1D: ...

    @file_unit.setter
    def file_unit(self, arg: Sequence[int], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> OutIoOutputDirectStruct: ...

    def __deepcopy__(self, arg: dict, /) -> OutIoOutputDirectStruct: ...

    def __eq__(self, arg: OutIoOutputDirectStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ParserControllerStruct:
    """Fortran struct: parser_controller_struct"""

    def __init__(self, name: str | None = None, attrib_name: str | None = None, y_knot: Sequence[float] | None = None, n_stk: int | None = None) -> None: ...

    @property
    def name(self) -> str: ...

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def attrib_name(self) -> str: ...

    @attrib_name.setter
    def attrib_name(self, arg: str, /) -> None: ...

    @property
    def stack(self) -> ExpressionAtomStructAlloc1D:
        """Arithmetic expression stack"""

    @property
    def y_knot(self) -> RealAlloc1D: ...

    @y_knot.setter
    def y_knot(self, arg: Sequence[float], /) -> None: ...

    @property
    def n_stk(self) -> int: ...

    @n_stk.setter
    def n_stk(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> ParserControllerStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> ParserControllerStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ParserControllerStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ParserControllerStruct: ...

    def __eq__(self, arg: ParserControllerStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ParserEleStruct:
    """Fortran struct: parser_ele_struct"""

    def __init__(self, ref_name: str | None = None, ix_super_ref_multipass: int | None = None, ele_name: str | None = None, lat_file: str | None = None, offset: float | None = None, ix_line_in_file: int | None = None, ix_count: int | None = None, ele_pt: int | None = None, ref_pt: int | None = None, index: int | None = None, superposition_command_here: bool | None = None, superposition_has_been_set: bool | None = None, wrap_superimpose: bool | None = None, create_jumbo_slave: bool | None = None, is_range: bool | None = None, default_attrib: str | None = None) -> None: ...

    @property
    def control(self) -> ParserControllerStructAlloc1D: ...

    @property
    def field_overlaps(self) -> FCharArray1D: ...

    @property
    def ref_name(self) -> str: ...

    @ref_name.setter
    def ref_name(self, arg: str, /) -> None: ...

    @property
    def ix_super_ref_multipass(self) -> int:
        """Multipass index for superimpose reference element."""

    @ix_super_ref_multipass.setter
    def ix_super_ref_multipass(self, arg: int, /) -> None: ...

    @property
    def ele_name(self) -> str:
        """For fork element or superimpose statement."""

    @ele_name.setter
    def ele_name(self, arg: str, /) -> None: ...

    @property
    def names1(self) -> FCharArray1D:
        """Currently just used by feedback element."""

    @property
    def names2(self) -> FCharArray1D:
        """Currently just used by feedback element."""

    @property
    def lat_file(self) -> str:
        """File where element was defined."""

    @lat_file.setter
    def lat_file(self, arg: str, /) -> None: ...

    @property
    def offset(self) -> float: ...

    @offset.setter
    def offset(self, arg: float, /) -> None: ...

    @property
    def ix_line_in_file(self) -> int:
        """Line in file where element was defined."""

    @ix_line_in_file.setter
    def ix_line_in_file(self, arg: int, /) -> None: ...

    @property
    def ix_count(self) -> int: ...

    @ix_count.setter
    def ix_count(self, arg: int, /) -> None: ...

    @property
    def ele_pt(self) -> int: ...

    @ele_pt.setter
    def ele_pt(self, arg: int, /) -> None: ...

    @property
    def ref_pt(self) -> int: ...

    @ref_pt.setter
    def ref_pt(self, arg: int, /) -> None: ...

    @property
    def index(self) -> int: ...

    @index.setter
    def index(self, arg: int, /) -> None: ...

    @property
    def superposition_command_here(self) -> bool: ...

    @superposition_command_here.setter
    def superposition_command_here(self, arg: bool, /) -> None: ...

    @property
    def superposition_has_been_set(self) -> bool: ...

    @superposition_has_been_set.setter
    def superposition_has_been_set(self, arg: bool, /) -> None: ...

    @property
    def wrap_superimpose(self) -> bool: ...

    @wrap_superimpose.setter
    def wrap_superimpose(self, arg: bool, /) -> None: ...

    @property
    def create_jumbo_slave(self) -> bool: ...

    @create_jumbo_slave.setter
    def create_jumbo_slave(self, arg: bool, /) -> None: ...

    @property
    def is_range(self) -> bool:
        """For girders"""

    @is_range.setter
    def is_range(self, arg: bool, /) -> None: ...

    @property
    def default_attrib(self) -> str:
        """For group/overlay elements: slave attribute"""

    @default_attrib.setter
    def default_attrib(self, arg: str, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> ParserEleStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> ParserEleStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ParserEleStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ParserEleStruct: ...

    def __eq__(self, arg: ParserEleStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ParserLatStruct:
    """Fortran struct: parser_lat_struct"""

    def __init__(self) -> None: ...

    @property
    def ele(self) -> ParserEleStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ParserLatStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ParserLatStruct: ...

    def __eq__(self, arg: ParserLatStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class PhotonCoordStruct:
    """Fortran struct: photon_coord_struct"""

    def __init__(self, orb: CoordStruct | None = None, track_len: float | None = None, ix_section: int | None = None) -> None: ...

    @property
    def orb(self) -> CoordStruct:
        """Phase space: orb%vec = (x, vx/c, y, vy/c, s, vs/c)"""

    @orb.setter
    def orb(self, arg: CoordStruct, /) -> None: ...

    @property
    def track_len(self) -> float:
        """Total track length from the start of the element."""

    @track_len.setter
    def track_len(self, arg: float, /) -> None: ...

    @property
    def ix_section(self) -> int:
        """Cross section index"""

    @ix_section.setter
    def ix_section(self, arg: int, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> PhotonCoordStruct: ...

    def __deepcopy__(self, arg: dict, /) -> PhotonCoordStruct: ...

    def __eq__(self, arg: PhotonCoordStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class PhotonElementStruct:
    """Fortran struct: photon_element_struct"""

    def __init__(self, curvature: SurfaceCurvatureStruct | None = None, target: PhotonTargetStruct | None = None, material: PhotonMaterialStruct | None = None, segmented: SurfaceSegmentedStruct | None = None, h_misalign: SurfaceHMisalignStruct | None = None, displacement: SurfaceDisplacementStruct | None = None, pixel: PixelDetecStruct | None = None, reflectivity_table_type: int | None = None, reflectivity_table_sigma: PhotonReflectTableStruct | None = None, reflectivity_table_pi: PhotonReflectTableStruct | None = None, integrated_init_energy_prob: Sequence[float] | None = None) -> None: ...

    @property
    def curvature(self) -> SurfaceCurvatureStruct: ...

    @curvature.setter
    def curvature(self, arg: SurfaceCurvatureStruct, /) -> None: ...

    @property
    def target(self) -> PhotonTargetStruct: ...

    @target.setter
    def target(self, arg: PhotonTargetStruct, /) -> None: ...

    @property
    def material(self) -> PhotonMaterialStruct: ...

    @material.setter
    def material(self, arg: PhotonMaterialStruct, /) -> None: ...

    @property
    def segmented(self) -> SurfaceSegmentedStruct: ...

    @segmented.setter
    def segmented(self, arg: SurfaceSegmentedStruct, /) -> None: ...

    @property
    def h_misalign(self) -> SurfaceHMisalignStruct: ...

    @h_misalign.setter
    def h_misalign(self, arg: SurfaceHMisalignStruct, /) -> None: ...

    @property
    def displacement(self) -> SurfaceDisplacementStruct: ...

    @displacement.setter
    def displacement(self, arg: SurfaceDisplacementStruct, /) -> None: ...

    @property
    def pixel(self) -> PixelDetecStruct: ...

    @pixel.setter
    def pixel(self, arg: PixelDetecStruct, /) -> None: ...

    @property
    def reflectivity_table_type(self) -> int: ...

    @reflectivity_table_type.setter
    def reflectivity_table_type(self, arg: int, /) -> None: ...

    @property
    def reflectivity_table_sigma(self) -> PhotonReflectTableStruct:
        """If polarization is ignored use sigma table."""

    @reflectivity_table_sigma.setter
    def reflectivity_table_sigma(self, arg: PhotonReflectTableStruct, /) -> None: ...

    @property
    def reflectivity_table_pi(self) -> PhotonReflectTableStruct: ...

    @reflectivity_table_pi.setter
    def reflectivity_table_pi(self, arg: PhotonReflectTableStruct, /) -> None: ...

    @property
    def init_energy_prob(self) -> SplineStructAlloc1D:
        """Initial energy probability density"""

    @property
    def integrated_init_energy_prob(self) -> RealAlloc1D: ...

    @integrated_init_energy_prob.setter
    def integrated_init_energy_prob(self, arg: Sequence[float], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> PhotonElementStruct: ...

    def __deepcopy__(self, arg: dict, /) -> PhotonElementStruct: ...

    def __eq__(self, arg: PhotonElementStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class PhotonInitSplinesStruct:
    """Fortran struct: photon_init_splines_struct"""

    def __init__(self, source_type: str | None = None, spline_space_dimensions: int | None = None) -> None: ...

    @property
    def source_type(self) -> str:
        """'bend', 'wiggler', 'undulator'"""

    @source_type.setter
    def source_type(self, arg: str, /) -> None: ...

    @property
    def spline_space_dimensions(self) -> int:
        """Dimensions: [energy, y_angle, x_angle, x, y]"""

    @spline_space_dimensions.setter
    def spline_space_dimensions(self, arg: int, /) -> None: ...

    @property
    def energy_prob(self) -> SplineStructAlloc1D: ...

    @property
    def y_angle(self) -> PhotonInitYAngleSplineStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> PhotonInitSplinesStruct: ...

    def __deepcopy__(self, arg: dict, /) -> PhotonInitSplinesStruct: ...

    def __eq__(self, arg: PhotonInitSplinesStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class PhotonInitXAngleSplineStruct:
    """Fortran struct: photon_init_x_angle_spline_struct"""

    def __init__(self) -> None: ...

    @property
    def prob(self) -> SplineStructAlloc1D: ...

    @property
    def pl(self) -> SplineStructAlloc1D: ...

    @property
    def pc(self) -> SplineStructAlloc1D: ...

    @property
    def pl45(self) -> SplineStructAlloc1D: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> PhotonInitXAngleSplineStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> PhotonInitXAngleSplineStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> PhotonInitXAngleSplineStruct: ...

    def __deepcopy__(self, arg: dict, /) -> PhotonInitXAngleSplineStruct: ...

    def __eq__(self, arg: PhotonInitXAngleSplineStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class PhotonInitYAngleSplineStruct:
    """Fortran struct: photon_init_y_angle_spline_struct"""

    def __init__(self) -> None: ...

    @property
    def prob(self) -> SplineStructAlloc1D: ...

    @property
    def pl(self) -> SplineStructAlloc1D: ...

    @property
    def pc(self) -> SplineStructAlloc1D: ...

    @property
    def pl45(self) -> SplineStructAlloc1D: ...

    @property
    def x_angle(self) -> PhotonInitXAngleSplineStructAlloc1D: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> PhotonInitYAngleSplineStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> PhotonInitYAngleSplineStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> PhotonInitYAngleSplineStruct: ...

    def __deepcopy__(self, arg: dict, /) -> PhotonInitYAngleSplineStruct: ...

    def __eq__(self, arg: PhotonInitYAngleSplineStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class PhotonMaterialStruct:
    """Fortran struct: photon_material_struct"""

    def __init__(self, f0_m1: complex | None = None, f0_m2: complex | None = None, f_0: complex | None = None, f_h: complex | None = None, f_hbar: complex | None = None, f_hkl: complex | None = None, h_norm: Sequence[float] | None = None, l_ref: Sequence[float] | None = None) -> None: ...

    @property
    def f0_m1(self) -> complex:
        """For multilayer_mirror only."""

    @f0_m1.setter
    def f0_m1(self, arg: complex, /) -> None: ...

    @property
    def f0_m2(self) -> complex:
        """For multilayer_mirror only."""

    @f0_m2.setter
    def f0_m2(self, arg: complex, /) -> None: ...

    @property
    def f_0(self) -> complex: ...

    @f_0.setter
    def f_0(self, arg: complex, /) -> None: ...

    @property
    def f_h(self) -> complex:
        """Structure factor for H direction."""

    @f_h.setter
    def f_h(self, arg: complex, /) -> None: ...

    @property
    def f_hbar(self) -> complex:
        """Structure factor for -H direction."""

    @f_hbar.setter
    def f_hbar(self, arg: complex, /) -> None: ...

    @property
    def f_hkl(self) -> complex:
        """= sqrt(f_h * f_hbar)"""

    @f_hkl.setter
    def f_hkl(self, arg: complex, /) -> None: ...

    @property
    def h_norm(self) -> RealArray1D:
        """Normalized H vector for crystals."""

    @h_norm.setter
    def h_norm(self, arg: Sequence[float], /) -> None: ...

    @property
    def l_ref(self) -> RealArray1D:
        """Crystal reference orbit displacement vector in element coords."""

    @l_ref.setter
    def l_ref(self, arg: Sequence[float], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> PhotonMaterialStruct: ...

    def __deepcopy__(self, arg: dict, /) -> PhotonMaterialStruct: ...

    def __eq__(self, arg: PhotonMaterialStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class PhotonReflectSurfaceStruct:
    """Fortran struct: photon_reflect_surface_struct"""

    def __init__(self, name: str | None = None, description: str | None = None, reflectivity_file: str | None = None, surface_roughness_rms: float | None = None, roughness_correlation_len: float | None = None, ix_surface: int | None = None) -> None: ...

    @property
    def name(self) -> str: ...

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def description(self) -> str:
        """Descriptive name"""

    @description.setter
    def description(self, arg: str, /) -> None: ...

    @property
    def reflectivity_file(self) -> str: ...

    @reflectivity_file.setter
    def reflectivity_file(self, arg: str, /) -> None: ...

    @property
    def table(self) -> PhotonReflectTableStructAlloc1D: ...

    @property
    def surface_roughness_rms(self) -> float:
        """sigma in Dugan's notation"""

    @surface_roughness_rms.setter
    def surface_roughness_rms(self, arg: float, /) -> None: ...

    @property
    def roughness_correlation_len(self) -> float:
        """T in Dugan's notation"""

    @roughness_correlation_len.setter
    def roughness_correlation_len(self, arg: float, /) -> None: ...

    @property
    def ix_surface(self) -> int: ...

    @ix_surface.setter
    def ix_surface(self, arg: int, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> PhotonReflectSurfaceStruct: ...

    def __deepcopy__(self, arg: dict, /) -> PhotonReflectSurfaceStruct: ...

    def __eq__(self, arg: PhotonReflectSurfaceStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class PhotonReflectTableStruct:
    """Fortran struct: photon_reflect_table_struct"""

    def __init__(self, angle: Sequence[float] | None = None, energy: Sequence[float] | None = None, p_reflect: Sequence[Sequence[float]] | None = None, max_energy: float | None = None, p_reflect_scratch: Sequence[float] | None = None, bragg_angle: Sequence[float] | None = None) -> None: ...

    @property
    def angle(self) -> RealAlloc1D:
        """Vector of angle values for %p_reflect"""

    @angle.setter
    def angle(self, arg: Sequence[float], /) -> None: ...

    @property
    def energy(self) -> RealAlloc1D:
        """Vector of energy values for %p_reflect"""

    @energy.setter
    def energy(self, arg: Sequence[float], /) -> None: ...

    @property
    def int1(self) -> Interval1CoefStructAlloc1D: ...

    @property
    def p_reflect(self) -> RealArray2D:
        """(angle, ev) probability. Log used for smooth surface reflection"""

    @p_reflect.setter
    def p_reflect(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def max_energy(self) -> float:
        """maximum energy for this table"""

    @max_energy.setter
    def max_energy(self, arg: float, /) -> None: ...

    @property
    def p_reflect_scratch(self) -> RealAlloc1D:
        """Scratch space"""

    @p_reflect_scratch.setter
    def p_reflect_scratch(self, arg: Sequence[float], /) -> None: ...

    @property
    def bragg_angle(self) -> RealAlloc1D:
        """Bragg angle at energy values."""

    @bragg_angle.setter
    def bragg_angle(self, arg: Sequence[float], /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> PhotonReflectTableStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> PhotonReflectTableStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> PhotonReflectTableStruct: ...

    def __deepcopy__(self, arg: dict, /) -> PhotonReflectTableStruct: ...

    def __eq__(self, arg: PhotonReflectTableStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class PhotonTargetStruct:
    """Fortran struct: photon_target_struct"""

    def __init__(self, type: int | None = None, n_corner: int | None = None, ele_loc: LatEleLocStruct | None = None, center: TargetPointStruct | None = None) -> None: ...

    @property
    def type(self) -> int:
        """or rectangular$"""

    @type.setter
    def type(self, arg: int, /) -> None: ...

    @property
    def n_corner(self) -> int: ...

    @n_corner.setter
    def n_corner(self, arg: int, /) -> None: ...

    @property
    def ele_loc(self) -> LatEleLocStruct: ...

    @ele_loc.setter
    def ele_loc(self, arg: LatEleLocStruct, /) -> None: ...

    @property
    def corner(self) -> TargetPointStructArray1D: ...

    @property
    def center(self) -> TargetPointStruct: ...

    @center.setter
    def center(self, arg: TargetPointStruct, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> PhotonTargetStruct: ...

    def __deepcopy__(self, arg: dict, /) -> PhotonTargetStruct: ...

    def __eq__(self, arg: PhotonTargetStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class PhotonTrackStruct:
    """Fortran struct: photon_track_struct"""

    def __init__(self, old: PhotonCoordStruct | None = None, now: PhotonCoordStruct | None = None) -> None: ...

    @property
    def old(self) -> PhotonCoordStruct: ...

    @old.setter
    def old(self, arg: PhotonCoordStruct, /) -> None: ...

    @property
    def now(self) -> PhotonCoordStruct: ...

    @now.setter
    def now(self, arg: PhotonCoordStruct, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> PhotonTrackStruct: ...

    def __deepcopy__(self, arg: dict, /) -> PhotonTrackStruct: ...

    def __eq__(self, arg: PhotonTrackStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class PixelDetecStruct:
    """Fortran struct: pixel_detec_struct"""

    def __init__(self, dr: Sequence[float] | None = None, r0: Sequence[float] | None = None, n_track_tot: int | None = None, n_hit_detec: int | None = None, n_hit_pixel: int | None = None) -> None: ...

    @property
    def dr(self) -> RealArray1D: ...

    @dr.setter
    def dr(self, arg: Sequence[float], /) -> None: ...

    @property
    def r0(self) -> RealArray1D: ...

    @r0.setter
    def r0(self, arg: Sequence[float], /) -> None: ...

    @property
    def n_track_tot(self) -> int:
        """How many photons were launched from source element."""

    @n_track_tot.setter
    def n_track_tot(self, arg: int, /) -> None: ...

    @property
    def n_hit_detec(self) -> int:
        """How many photons hit the detector."""

    @n_hit_detec.setter
    def n_hit_detec(self, arg: int, /) -> None: ...

    @property
    def n_hit_pixel(self) -> int:
        """How many photons hit the pixel grid of the detector."""

    @n_hit_pixel.setter
    def n_hit_pixel(self, arg: int, /) -> None: ...

    @property
    def pt(self) -> PixelPtStructArray2D:
        """Grid of pixels"""

    def __repr__(self) -> str: ...

    def __copy__(self) -> PixelDetecStruct: ...

    def __deepcopy__(self, arg: dict, /) -> PixelDetecStruct: ...

    def __eq__(self, arg: PixelDetecStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class PixelPtStruct:
    """Fortran struct: pixel_pt_struct"""

    def __init__(self, n_photon: int | None = None, E_x: complex | None = None, E_y: complex | None = None, intensity_x: float | None = None, intensity_y: float | None = None, intensity: float | None = None, orbit: Sequence[float] | None = None, orbit_rms: Sequence[float] | None = None, init_orbit: Sequence[float] | None = None, init_orbit_rms: Sequence[float] | None = None) -> None: ...

    @property
    def n_photon(self) -> int: ...

    @n_photon.setter
    def n_photon(self, arg: int, /) -> None: ...

    @property
    def E_x(self) -> complex: ...

    @E_x.setter
    def E_x(self, arg: complex, /) -> None: ...

    @property
    def E_y(self) -> complex: ...

    @E_y.setter
    def E_y(self, arg: complex, /) -> None: ...

    @property
    def intensity_x(self) -> float: ...

    @intensity_x.setter
    def intensity_x(self, arg: float, /) -> None: ...

    @property
    def intensity_y(self) -> float: ...

    @intensity_y.setter
    def intensity_y(self, arg: float, /) -> None: ...

    @property
    def intensity(self) -> float: ...

    @intensity.setter
    def intensity(self, arg: float, /) -> None: ...

    @property
    def orbit(self) -> RealArray1D:
        """x, Vx/c, y, Vy/c, dummy, E - E_ref."""

    @orbit.setter
    def orbit(self, arg: Sequence[float], /) -> None: ...

    @property
    def orbit_rms(self) -> RealArray1D:
        """RMS statistics."""

    @orbit_rms.setter
    def orbit_rms(self, arg: Sequence[float], /) -> None: ...

    @property
    def init_orbit(self) -> RealArray1D:
        """Initial orbit at start of lattice statistics."""

    @init_orbit.setter
    def init_orbit(self, arg: Sequence[float], /) -> None: ...

    @property
    def init_orbit_rms(self) -> RealArray1D:
        """Initial orbit at start of lattice RMS statistics."""

    @init_orbit_rms.setter
    def init_orbit_rms(self, arg: Sequence[float], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> PixelPtStruct: ...

    def __deepcopy__(self, arg: dict, /) -> PixelPtStruct: ...

    def __eq__(self, arg: PixelPtStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class PmdHeaderStruct:
    """Fortran struct: pmd_header_struct"""

    def __init__(self, openPMD: str | None = None, openPMDextension: str | None = None, basePath: str | None = None, particlesPath: str | None = None, meshesPath: str | None = None, author: str | None = None, software: str | None = None, softwareVersion: str | None = None, date: str | None = None, latticeFile: str | None = None, latticeName: str | None = None) -> None: ...

    @property
    def openPMD(self) -> str: ...

    @openPMD.setter
    def openPMD(self, arg: str, /) -> None: ...

    @property
    def openPMDextension(self) -> str: ...

    @openPMDextension.setter
    def openPMDextension(self, arg: str, /) -> None: ...

    @property
    def basePath(self) -> str: ...

    @basePath.setter
    def basePath(self, arg: str, /) -> None: ...

    @property
    def particlesPath(self) -> str: ...

    @particlesPath.setter
    def particlesPath(self, arg: str, /) -> None: ...

    @property
    def meshesPath(self) -> str: ...

    @meshesPath.setter
    def meshesPath(self, arg: str, /) -> None: ...

    @property
    def author(self) -> str: ...

    @author.setter
    def author(self, arg: str, /) -> None: ...

    @property
    def software(self) -> str: ...

    @software.setter
    def software(self, arg: str, /) -> None: ...

    @property
    def softwareVersion(self) -> str: ...

    @softwareVersion.setter
    def softwareVersion(self, arg: str, /) -> None: ...

    @property
    def date(self) -> str: ...

    @date.setter
    def date(self, arg: str, /) -> None: ...

    @property
    def latticeFile(self) -> str: ...

    @latticeFile.setter
    def latticeFile(self, arg: str, /) -> None: ...

    @property
    def latticeName(self) -> str: ...

    @latticeName.setter
    def latticeName(self, arg: str, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> PmdHeaderStruct: ...

    def __deepcopy__(self, arg: dict, /) -> PmdHeaderStruct: ...

    def __eq__(self, arg: PmdHeaderStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class PreTrackerStruct:
    """Fortran struct: pre_tracker_struct"""

    def __init__(self, who: int | None = None, ix_ele_start: int | None = None, ix_ele_end: int | None = None, input_file: str | None = None) -> None: ...

    @property
    def who(self) -> int:
        """Can be opal$, or impactt$"""

    @who.setter
    def who(self, arg: int, /) -> None: ...

    @property
    def ix_ele_start(self) -> int: ...

    @ix_ele_start.setter
    def ix_ele_start(self, arg: int, /) -> None: ...

    @property
    def ix_ele_end(self) -> int: ...

    @ix_ele_end.setter
    def ix_ele_end(self, arg: int, /) -> None: ...

    @property
    def input_file(self) -> str: ...

    @input_file.setter
    def input_file(self, arg: str, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> PreTrackerStruct: ...

    def __deepcopy__(self, arg: dict, /) -> PreTrackerStruct: ...

    def __eq__(self, arg: PreTrackerStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class PtcBranch1Struct:
    """Fortran struct: ptc_branch1_struct"""

    def __init__(self) -> None: ...

    @property
    def m_u_layout(self) -> PtcLayoutPointerStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> PtcBranch1Struct: ...

    def __deepcopy__(self, arg: dict, /) -> PtcBranch1Struct: ...

    def __eq__(self, arg: PtcBranch1Struct, /) -> bool: ...

    def __hash__(self) -> int: ...

class PtcLayoutPointerStruct:
    """Fortran struct: ptc_layout_pointer_struct"""

    def __init__(self) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> PtcLayoutPointerStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> PtcLayoutPointerStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> PtcLayoutPointerStruct: ...

    def __deepcopy__(self, arg: dict, /) -> PtcLayoutPointerStruct: ...

    def __eq__(self, arg: PtcLayoutPointerStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class PtcNormalFormStruct:
    """Fortran struct: ptc_normal_form_struct"""

    def __init__(self, ele_origin: EleStruct | None = None, orb0: Sequence[float] | None = None, valid_map: bool | None = None) -> None: ...

    @property
    def ele_origin(self) -> EleStruct | None:
        """Element at which the on-turn map was created."""

    @ele_origin.setter
    def ele_origin(self, arg: EleStruct, /) -> None: ...

    @property
    def orb0(self) -> RealArray1D:
        """Closed orbit at element."""

    @orb0.setter
    def orb0(self, arg: Sequence[float], /) -> None: ...

    @property
    def valid_map(self) -> bool: ...

    @valid_map.setter
    def valid_map(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> PtcNormalFormStruct: ...

    def __deepcopy__(self, arg: dict, /) -> PtcNormalFormStruct: ...

    def __eq__(self, arg: PtcNormalFormStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class PtcRadMapStruct:
    """Fortran struct: ptc_rad_map_struct"""

    def __init__(self, lattice_file: str | None = None, dref_time: float | None = None, p0c_start: float | None = None, p0c_end: float | None = None, s_end: float | None = None, map_order: int | None = None, radiation_damping_on: bool | None = None, ix_branch: int | None = None, ix_ele_start: int | None = None, ix_ele_end: int | None = None, nodamp_mat: Sequence[Sequence[float]] | None = None, damp_mat: Sequence[Sequence[float]] | None = None, stoc_mat: Sequence[Sequence[float]] | None = None, ref0: Sequence[float] | None = None, ref1: Sequence[float] | None = None) -> None: ...

    @property
    def lattice_file(self) -> str:
        """Name of the lattice file"""

    @lattice_file.setter
    def lattice_file(self, arg: str, /) -> None: ...

    @property
    def dref_time(self) -> float:
        """Time ref particle takes."""

    @dref_time.setter
    def dref_time(self, arg: float, /) -> None: ...

    @property
    def p0c_start(self) -> float:
        """ref momentum at start"""

    @p0c_start.setter
    def p0c_start(self, arg: float, /) -> None: ...

    @property
    def p0c_end(self) -> float:
        """ref momentum at end"""

    @p0c_end.setter
    def p0c_end(self, arg: float, /) -> None: ...

    @property
    def s_end(self) -> float:
        """Ending s-position"""

    @s_end.setter
    def s_end(self, arg: float, /) -> None: ...

    @property
    def map_order(self) -> int: ...

    @map_order.setter
    def map_order(self, arg: int, /) -> None: ...

    @property
    def radiation_damping_on(self) -> bool: ...

    @radiation_damping_on.setter
    def radiation_damping_on(self, arg: bool, /) -> None: ...

    @property
    def ix_branch(self) -> int: ...

    @ix_branch.setter
    def ix_branch(self, arg: int, /) -> None: ...

    @property
    def ix_ele_start(self) -> int:
        """Start point for making the map"""

    @ix_ele_start.setter
    def ix_ele_start(self, arg: int, /) -> None: ...

    @property
    def ix_ele_end(self) -> int:
        """End point for making the map"""

    @ix_ele_end.setter
    def ix_ele_end(self, arg: int, /) -> None: ...

    @property
    def nodamp_mat(self) -> RealArray2D:
        """Nondamped orbital matrix. M_orbit = M_damp * M_nodamp"""

    @nodamp_mat.setter
    def nodamp_mat(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def damp_mat(self) -> RealArray2D:
        """
        Damping 'correction' to orbital matrix. Stoc_mat is referenced to the start of the map. That is, it is applied before the transport matrix.
        """

    @damp_mat.setter
    def damp_mat(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def stoc_mat(self) -> RealArray2D:
        """Stochatic matrix for the orbit."""

    @stoc_mat.setter
    def stoc_mat(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def ref0(self) -> RealArray1D:
        """Reference orbit at start."""

    @ref0.setter
    def ref0(self, arg: Sequence[float], /) -> None: ...

    @property
    def ref1(self) -> RealArray1D:
        """Reference orbit at end."""

    @ref1.setter
    def ref1(self, arg: Sequence[float], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> PtcRadMapStruct: ...

    def __deepcopy__(self, arg: dict, /) -> PtcRadMapStruct: ...

    def __eq__(self, arg: PtcRadMapStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class QpAxisStruct:
    """Fortran struct: qp_axis_struct"""

    def __init__(self, label: str | None = None, min: float | None = None, max: float | None = None, tick_min: float | None = None, tick_max: float | None = None, eval_min: float | None = None, eval_max: float | None = None, dtick: float | None = None, number_offset: float | None = None, label_offset: float | None = None, major_tick_len: float | None = None, minor_tick_len: float | None = None, label_color: str | None = None, major_div: int | None = None, major_div_nominal: int | None = None, minor_div: int | None = None, minor_div_max: int | None = None, places: int | None = None, type: str | None = None, bounds: str | None = None, tick_side: int | None = None, number_side: int | None = None, draw_label: bool | None = None, draw_numbers: bool | None = None) -> None: ...

    @property
    def label(self) -> str: ...

    @label.setter
    def label(self, arg: str, /) -> None: ...

    @property
    def min(self) -> float:
        """Axis min/max in data units."""

    @min.setter
    def min(self, arg: float, /) -> None: ...

    @property
    def max(self) -> float:
        """Axis min/max in data units."""

    @max.setter
    def max(self, arg: float, /) -> None: ...

    @property
    def tick_min(self) -> float:
        """Min tick location along axis in data units."""

    @tick_min.setter
    def tick_min(self, arg: float, /) -> None: ...

    @property
    def tick_max(self) -> float:
        """Max tick location along axis in data units."""

    @tick_max.setter
    def tick_max(self, arg: float, /) -> None: ...

    @property
    def eval_min(self) -> float:
        """For general use. Not set by quick_plot."""

    @eval_min.setter
    def eval_min(self, arg: float, /) -> None: ...

    @property
    def eval_max(self) -> float:
        """For general use. Not set by quick_plot."""

    @eval_max.setter
    def eval_max(self, arg: float, /) -> None: ...

    @property
    def dtick(self) -> float:
        """
        Distance between ticks. In data units. Ticks will be drawn between %min and %max.
        """

    @dtick.setter
    def dtick(self, arg: float, /) -> None: ...

    @property
    def number_offset(self) -> float:
        """Offset from axis line in inches."""

    @number_offset.setter
    def number_offset(self, arg: float, /) -> None: ...

    @property
    def label_offset(self) -> float:
        """Offset from numbers in inches."""

    @label_offset.setter
    def label_offset(self, arg: float, /) -> None: ...

    @property
    def major_tick_len(self) -> float:
        """In inches."""

    @major_tick_len.setter
    def major_tick_len(self, arg: float, /) -> None: ...

    @property
    def minor_tick_len(self) -> float:
        """In inches."""

    @minor_tick_len.setter
    def minor_tick_len(self, arg: float, /) -> None: ...

    @property
    def label_color(self) -> str:
        """Color of the label."""

    @label_color.setter
    def label_color(self, arg: str, /) -> None: ...

    @property
    def major_div(self) -> int:
        """Actual numbrer of major divisions"""

    @major_div.setter
    def major_div(self, arg: int, /) -> None: ...

    @property
    def major_div_nominal(self) -> int:
        """Nominal value."""

    @major_div_nominal.setter
    def major_div_nominal(self, arg: int, /) -> None: ...

    @property
    def minor_div(self) -> int:
        """0 = auto choose."""

    @minor_div.setter
    def minor_div(self, arg: int, /) -> None: ...

    @property
    def minor_div_max(self) -> int:
        """Max number for auto choose."""

    @minor_div_max.setter
    def minor_div_max(self, arg: int, /) -> None: ...

    @property
    def places(self) -> int:
        """Number of places after the decimal point to print."""

    @places.setter
    def places(self, arg: int, /) -> None: ...

    @property
    def type(self) -> str:
        """Or 'LOG', or 'CUSTOM'"""

    @type.setter
    def type(self, arg: str, /) -> None: ...

    @property
    def bounds(self) -> str:
        """Or 'ZERO_AT_END' or 'ZERO_SYMMETRIC'"""

    @bounds.setter
    def bounds(self, arg: str, /) -> None: ...

    @property
    def tick_side(self) -> int:
        """
        +1 = Draw on the side inside the graph, 0 = both (longer tick), -1 = outside.
        """

    @tick_side.setter
    def tick_side(self, arg: int, /) -> None: ...

    @property
    def number_side(self) -> int:
        """+1 = Draw to the side inside the graph, -1 = outside."""

    @number_side.setter
    def number_side(self, arg: int, /) -> None: ...

    @property
    def draw_label(self) -> bool: ...

    @draw_label.setter
    def draw_label(self, arg: bool, /) -> None: ...

    @property
    def draw_numbers(self) -> bool: ...

    @draw_numbers.setter
    def draw_numbers(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> QpAxisStruct: ...

    def __deepcopy__(self, arg: dict, /) -> QpAxisStruct: ...

    def __eq__(self, arg: QpAxisStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class QpLegendStruct:
    """Fortran struct: qp_legend_struct"""

    def __init__(self, row_spacing: float | None = None, line_length: float | None = None, text_offset: float | None = None, draw_line: bool | None = None, draw_symbol: bool | None = None, draw_text: bool | None = None) -> None: ...

    @property
    def row_spacing(self) -> float:
        """Spacing between rows."""

    @row_spacing.setter
    def row_spacing(self, arg: float, /) -> None: ...

    @property
    def line_length(self) -> float:
        """Length of the line in points."""

    @line_length.setter
    def line_length(self, arg: float, /) -> None: ...

    @property
    def text_offset(self) -> float:
        """Horizontal offset in points between the line and the text."""

    @text_offset.setter
    def text_offset(self, arg: float, /) -> None: ...

    @property
    def draw_line(self) -> bool:
        """Draw lines?"""

    @draw_line.setter
    def draw_line(self, arg: bool, /) -> None: ...

    @property
    def draw_symbol(self) -> bool:
        """Draw symbols?"""

    @draw_symbol.setter
    def draw_symbol(self, arg: bool, /) -> None: ...

    @property
    def draw_text(self) -> bool:
        """Draw text?"""

    @draw_text.setter
    def draw_text(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> QpLegendStruct: ...

    def __deepcopy__(self, arg: dict, /) -> QpLegendStruct: ...

    def __eq__(self, arg: QpLegendStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class QpLineStruct:
    """Fortran struct: qp_line_struct"""

    def __init__(self, width: int | None = None, color: str | None = None, pattern: str | None = None) -> None: ...

    @property
    def width(self) -> int: ...

    @width.setter
    def width(self, arg: int, /) -> None: ...

    @property
    def color(self) -> str: ...

    @color.setter
    def color(self, arg: str, /) -> None: ...

    @property
    def pattern(self) -> str: ...

    @pattern.setter
    def pattern(self, arg: str, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> QpLineStruct: ...

    def __deepcopy__(self, arg: dict, /) -> QpLineStruct: ...

    def __eq__(self, arg: QpLineStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class QpPointStruct:
    """Fortran struct: qp_point_struct"""

    def __init__(self, x: float | None = None, y: float | None = None, units: str | None = None) -> None: ...

    @property
    def x(self) -> float: ...

    @x.setter
    def x(self, arg: float, /) -> None: ...

    @property
    def y(self) -> float: ...

    @y.setter
    def y(self, arg: float, /) -> None: ...

    @property
    def units(self) -> str: ...

    @units.setter
    def units(self, arg: str, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> QpPointStruct: ...

    def __deepcopy__(self, arg: dict, /) -> QpPointStruct: ...

    def __eq__(self, arg: QpPointStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class QpRectStruct:
    """Fortran struct: qp_rect_struct"""

    def __init__(self, x1: float | None = None, x2: float | None = None, y1: float | None = None, y2: float | None = None, units: str | None = None) -> None: ...

    @property
    def x1(self) -> float: ...

    @x1.setter
    def x1(self, arg: float, /) -> None: ...

    @property
    def x2(self) -> float: ...

    @x2.setter
    def x2(self, arg: float, /) -> None: ...

    @property
    def y1(self) -> float: ...

    @y1.setter
    def y1(self, arg: float, /) -> None: ...

    @property
    def y2(self) -> float: ...

    @y2.setter
    def y2(self, arg: float, /) -> None: ...

    @property
    def units(self) -> str: ...

    @units.setter
    def units(self, arg: str, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> QpRectStruct: ...

    def __deepcopy__(self, arg: dict, /) -> QpRectStruct: ...

    def __eq__(self, arg: QpRectStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class QpSymbolStruct:
    """Fortran struct: qp_symbol_struct"""

    def __init__(self, type: str | None = None, height: float | None = None, color: str | None = None, fill_pattern: str | None = None, line_width: int | None = None) -> None: ...

    @property
    def type(self) -> str: ...

    @type.setter
    def type(self, arg: str, /) -> None: ...

    @property
    def height(self) -> float:
        """in points (same as text height)"""

    @height.setter
    def height(self, arg: float, /) -> None: ...

    @property
    def color(self) -> str: ...

    @color.setter
    def color(self, arg: str, /) -> None: ...

    @property
    def fill_pattern(self) -> str: ...

    @fill_pattern.setter
    def fill_pattern(self, arg: str, /) -> None: ...

    @property
    def line_width(self) -> int: ...

    @line_width.setter
    def line_width(self, arg: int, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> QpSymbolStruct: ...

    def __deepcopy__(self, arg: dict, /) -> QpSymbolStruct: ...

    def __eq__(self, arg: QpSymbolStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class RadInt1Struct:
    """Fortran struct: rad_int1_struct"""

    def __init__(self, i0: float | None = None, i1: float | None = None, i2: float | None = None, i3: float | None = None, i4a: float | None = None, i4b: float | None = None, i4z: float | None = None, i5a: float | None = None, i5b: float | None = None, i6b: float | None = None, lin_i2_E4: float | None = None, lin_i3_E7: float | None = None, lin_i5a_E6: float | None = None, lin_i5b_E6: float | None = None, lin_norm_emit_a: float | None = None, lin_norm_emit_b: float | None = None, lin_sig_E: float | None = None, n_steps: float | None = None) -> None: ...

    @property
    def i0(self) -> float: ...

    @i0.setter
    def i0(self, arg: float, /) -> None: ...

    @property
    def i1(self) -> float: ...

    @i1.setter
    def i1(self, arg: float, /) -> None: ...

    @property
    def i2(self) -> float: ...

    @i2.setter
    def i2(self, arg: float, /) -> None: ...

    @property
    def i3(self) -> float: ...

    @i3.setter
    def i3(self, arg: float, /) -> None: ...

    @property
    def i4a(self) -> float: ...

    @i4a.setter
    def i4a(self, arg: float, /) -> None: ...

    @property
    def i4b(self) -> float: ...

    @i4b.setter
    def i4b(self, arg: float, /) -> None: ...

    @property
    def i4z(self) -> float: ...

    @i4z.setter
    def i4z(self, arg: float, /) -> None: ...

    @property
    def i5a(self) -> float: ...

    @i5a.setter
    def i5a(self, arg: float, /) -> None: ...

    @property
    def i5b(self) -> float: ...

    @i5b.setter
    def i5b(self, arg: float, /) -> None: ...

    @property
    def i6b(self) -> float: ...

    @i6b.setter
    def i6b(self, arg: float, /) -> None: ...

    @property
    def lin_i2_E4(self) -> float: ...

    @lin_i2_E4.setter
    def lin_i2_E4(self, arg: float, /) -> None: ...

    @property
    def lin_i3_E7(self) -> float: ...

    @lin_i3_E7.setter
    def lin_i3_E7(self, arg: float, /) -> None: ...

    @property
    def lin_i5a_E6(self) -> float: ...

    @lin_i5a_E6.setter
    def lin_i5a_E6(self, arg: float, /) -> None: ...

    @property
    def lin_i5b_E6(self) -> float: ...

    @lin_i5b_E6.setter
    def lin_i5b_E6(self, arg: float, /) -> None: ...

    @property
    def lin_norm_emit_a(self) -> float:
        """Running sum"""

    @lin_norm_emit_a.setter
    def lin_norm_emit_a(self, arg: float, /) -> None: ...

    @property
    def lin_norm_emit_b(self) -> float:
        """Running sum"""

    @lin_norm_emit_b.setter
    def lin_norm_emit_b(self, arg: float, /) -> None: ...

    @property
    def lin_sig_E(self) -> float:
        """Running sum"""

    @lin_sig_E.setter
    def lin_sig_E(self, arg: float, /) -> None: ...

    @property
    def n_steps(self) -> float:
        """number of qromb steps needed"""

    @n_steps.setter
    def n_steps(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> RadInt1StructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> RadInt1StructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> RadInt1Struct: ...

    def __deepcopy__(self, arg: dict, /) -> RadInt1Struct: ...

    def __eq__(self, arg: RadInt1Struct, /) -> bool: ...

    def __hash__(self) -> int: ...

class RadIntAllEleStruct:
    """Fortran struct: rad_int_all_ele_struct"""

    def __init__(self) -> None: ...

    @property
    def branch(self) -> RadIntBranchStructAlloc1D:
        """Array is indexed from 0"""

    def __repr__(self) -> str: ...

    def __copy__(self) -> RadIntAllEleStruct: ...

    def __deepcopy__(self, arg: dict, /) -> RadIntAllEleStruct: ...

    def __eq__(self, arg: RadIntAllEleStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class RadIntBranchStruct:
    """Fortran struct: rad_int_branch_struct"""

    def __init__(self) -> None: ...

    @property
    def ele(self) -> RadInt1StructAlloc1D:
        """Array is indexed from 0"""

    @staticmethod
    def new_array1d(sz: int = 0) -> RadIntBranchStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> RadIntBranchStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> RadIntBranchStruct: ...

    def __deepcopy__(self, arg: dict, /) -> RadIntBranchStruct: ...

    def __eq__(self, arg: RadIntBranchStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class RadIntCache1Struct:
    """Fortran struct: rad_int_cache1_struct"""

    def __init__(self, n_pt: int | None = None, cache_type: int | None = None) -> None: ...

    @property
    def pt(self) -> RadIntTrackPointStructAlloc1D:
        """pt(0:n_pt)"""

    @property
    def n_pt(self) -> int:
        """Upper bound of pt(0:n_pt)"""

    @n_pt.setter
    def n_pt(self, arg: int, /) -> None: ...

    @property
    def cache_type(self) -> int: ...

    @cache_type.setter
    def cache_type(self, arg: int, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> RadIntCache1Struct: ...

    def __deepcopy__(self, arg: dict, /) -> RadIntCache1Struct: ...

    def __eq__(self, arg: RadIntCache1Struct, /) -> bool: ...

    def __hash__(self) -> int: ...

class RadIntInfoStruct:
    """Fortran struct: rad_int_info_struct"""

    def __init__(self, branch: BranchStruct | None = None, ele: EleStruct | None = None, a: TwissStruct | None = None, b: TwissStruct | None = None, cache_ele: RadIntCache1Struct | None = None, eta_a: Sequence[float] | None = None, eta_b: Sequence[float] | None = None, g: float | None = None, g2: float | None = None, g_x: float | None = None, g_y: float | None = None, dg2_x: float | None = None, dg2_y: float | None = None) -> None: ...

    @property
    def branch(self) -> BranchStruct | None: ...

    @branch.setter
    def branch(self, arg: BranchStruct, /) -> None: ...

    @property
    def ele(self) -> EleStruct | None: ...

    @ele.setter
    def ele(self, arg: EleStruct, /) -> None: ...

    @property
    def orbit(self) -> CoordStructArray1D: ...

    @property
    def a(self) -> TwissStruct: ...

    @a.setter
    def a(self, arg: TwissStruct, /) -> None: ...

    @property
    def b(self) -> TwissStruct: ...

    @b.setter
    def b(self, arg: TwissStruct, /) -> None: ...

    @property
    def cache_ele(self) -> RadIntCache1Struct | None:
        """pointer to cache in use"""

    @cache_ele.setter
    def cache_ele(self, arg: RadIntCache1Struct, /) -> None: ...

    @property
    def eta_a(self) -> RealArray1D: ...

    @eta_a.setter
    def eta_a(self, arg: Sequence[float], /) -> None: ...

    @property
    def eta_b(self) -> RealArray1D: ...

    @eta_b.setter
    def eta_b(self, arg: Sequence[float], /) -> None: ...

    @property
    def g(self) -> float:
        """bending strength (1/bending_radius)"""

    @g.setter
    def g(self, arg: float, /) -> None: ...

    @property
    def g2(self) -> float:
        """bending strength (1/bending_radius)"""

    @g2.setter
    def g2(self, arg: float, /) -> None: ...

    @property
    def g_x(self) -> float:
        """components in x-y plane"""

    @g_x.setter
    def g_x(self, arg: float, /) -> None: ...

    @property
    def g_y(self) -> float:
        """components in x-y plane"""

    @g_y.setter
    def g_y(self, arg: float, /) -> None: ...

    @property
    def dg2_x(self) -> float: ...

    @dg2_x.setter
    def dg2_x(self, arg: float, /) -> None: ...

    @property
    def dg2_y(self) -> float: ...

    @dg2_y.setter
    def dg2_y(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> RadIntInfoStruct: ...

    def __deepcopy__(self, arg: dict, /) -> RadIntInfoStruct: ...

    def __eq__(self, arg: RadIntInfoStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class RadIntTrackPointStruct:
    """Fortran struct: rad_int_track_point_struct"""

    def __init__(self, s_body: float | None = None, mat6: Sequence[Sequence[float]] | None = None, vec0: Sequence[float] | None = None, ref_orb_in: CoordStruct | None = None, ref_orb_out: CoordStruct | None = None, g_x0: float | None = None, g_y0: float | None = None, dgx_dx: float | None = None, dgx_dy: float | None = None, dgy_dx: float | None = None, dgy_dy: float | None = None) -> None: ...

    @property
    def s_body(self) -> float: ...

    @s_body.setter
    def s_body(self, arg: float, /) -> None: ...

    @property
    def mat6(self) -> RealArray2D: ...

    @mat6.setter
    def mat6(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def vec0(self) -> RealArray1D: ...

    @vec0.setter
    def vec0(self, arg: Sequence[float], /) -> None: ...

    @property
    def ref_orb_in(self) -> CoordStruct: ...

    @ref_orb_in.setter
    def ref_orb_in(self, arg: CoordStruct, /) -> None: ...

    @property
    def ref_orb_out(self) -> CoordStruct: ...

    @ref_orb_out.setter
    def ref_orb_out(self, arg: CoordStruct, /) -> None: ...

    @property
    def g_x0(self) -> float:
        """Additional g factors for bends."""

    @g_x0.setter
    def g_x0(self, arg: float, /) -> None: ...

    @property
    def g_y0(self) -> float:
        """Additional g factors for bends."""

    @g_y0.setter
    def g_y0(self, arg: float, /) -> None: ...

    @property
    def dgx_dx(self) -> float:
        """bending strength gradient"""

    @dgx_dx.setter
    def dgx_dx(self, arg: float, /) -> None: ...

    @property
    def dgx_dy(self) -> float:
        """bending strength gradient"""

    @dgx_dy.setter
    def dgx_dy(self, arg: float, /) -> None: ...

    @property
    def dgy_dx(self) -> float:
        """bending strength gradient"""

    @dgy_dx.setter
    def dgy_dx(self, arg: float, /) -> None: ...

    @property
    def dgy_dy(self) -> float:
        """bending strength gradient"""

    @dgy_dy.setter
    def dgy_dy(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> RadIntTrackPointStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> RadIntTrackPointStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> RadIntTrackPointStruct: ...

    def __deepcopy__(self, arg: dict, /) -> RadIntTrackPointStruct: ...

    def __eq__(self, arg: RadIntTrackPointStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class RadMapEleStruct:
    """Fortran struct: rad_map_ele_struct"""

    def __init__(self, rm0: RadMapStruct | None = None, rm1: RadMapStruct | None = None, stale: bool | None = None) -> None: ...

    @property
    def rm0(self) -> RadMapStruct:
        """Upstream half and downstream half matrices for an element."""

    @rm0.setter
    def rm0(self, arg: RadMapStruct, /) -> None: ...

    @property
    def rm1(self) -> RadMapStruct:
        """Upstream half and downstream half matrices for an element."""

    @rm1.setter
    def rm1(self, arg: RadMapStruct, /) -> None: ...

    @property
    def stale(self) -> bool: ...

    @stale.setter
    def stale(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> RadMapEleStruct: ...

    def __deepcopy__(self, arg: dict, /) -> RadMapEleStruct: ...

    def __eq__(self, arg: RadMapEleStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class RadMapStruct:
    """Fortran struct: rad_map_struct"""

    def __init__(self, ref_orb: Sequence[float] | None = None, damp_dmat: Sequence[Sequence[float]] | None = None, xfer_damp_vec: Sequence[float] | None = None, xfer_damp_mat: Sequence[Sequence[float]] | None = None, stoc_mat: Sequence[Sequence[float]] | None = None) -> None: ...

    @property
    def ref_orb(self) -> RealArray1D:
        """Reference point around which damp_mat is calculated."""

    @ref_orb.setter
    def ref_orb(self, arg: Sequence[float], /) -> None: ...

    @property
    def damp_dmat(self) -> RealArray2D:
        """damp_correction = xfer_mat_with_damping - xfer_mat_without_damping."""

    @damp_dmat.setter
    def damp_dmat(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def xfer_damp_vec(self) -> RealArray1D:
        """Transfer map with damping 0th order vector."""

    @xfer_damp_vec.setter
    def xfer_damp_vec(self, arg: Sequence[float], /) -> None: ...

    @property
    def xfer_damp_mat(self) -> RealArray2D:
        """1st order matrix: xfer_no_damp_mat + xfer_damp_correction."""

    @xfer_damp_mat.setter
    def xfer_damp_mat(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def stoc_mat(self) -> RealArray2D:
        """Stochastic variance or 'kick' (Cholesky decomposed) matrix."""

    @stoc_mat.setter
    def stoc_mat(self, arg: Sequence[Sequence[float]], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> RadMapStruct: ...

    def __deepcopy__(self, arg: dict, /) -> RadMapStruct: ...

    def __eq__(self, arg: RadMapStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class RamperLordStruct:
    """Fortran struct: ramper_lord_struct"""

    def __init__(self, ix_ele: int | None = None, ix_con: int | None = None, attrib_ptr: float | None = None) -> None: ...

    @property
    def ix_ele(self) -> int:
        """Lord index"""

    @ix_ele.setter
    def ix_ele(self, arg: int, /) -> None: ...

    @property
    def ix_con(self) -> int:
        """Index in lord%control%ramp(:) array"""

    @ix_con.setter
    def ix_con(self, arg: int, /) -> None: ...

    @property
    def attrib_ptr(self) -> float | None:
        """Pointer to attribute in this element."""

    @attrib_ptr.setter
    def attrib_ptr(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> RamperLordStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> RamperLordStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> RamperLordStruct: ...

    def __deepcopy__(self, arg: dict, /) -> RamperLordStruct: ...

    def __eq__(self, arg: RamperLordStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class RandomStateStruct:
    """Fortran struct: random_state_struct"""

    def __init__(self, ix: int | None = None, iy: int | None = None, number_stored: bool | None = None, h_saved: float | None = None, engine: int | None = None, seed: int | None = None, am: float | None = None, gauss_converter: int | None = None, gauss_sigma_cut: float | None = None, in_sobseq: int | None = None, ix_sobseq: Sequence[int] | None = None, x_sobseq: Sequence[float] | None = None) -> None: ...

    @property
    def ix(self) -> int: ...

    @ix.setter
    def ix(self, arg: int, /) -> None: ...

    @property
    def iy(self) -> int: ...

    @iy.setter
    def iy(self, arg: int, /) -> None: ...

    @property
    def number_stored(self) -> bool: ...

    @number_stored.setter
    def number_stored(self, arg: bool, /) -> None: ...

    @property
    def h_saved(self) -> float: ...

    @h_saved.setter
    def h_saved(self, arg: float, /) -> None: ...

    @property
    def engine(self) -> int:
        """Params"""

    @engine.setter
    def engine(self, arg: int, /) -> None: ...

    @property
    def seed(self) -> int: ...

    @seed.setter
    def seed(self, arg: int, /) -> None: ...

    @property
    def am(self) -> float: ...

    @am.setter
    def am(self, arg: float, /) -> None: ...

    @property
    def gauss_converter(self) -> int: ...

    @gauss_converter.setter
    def gauss_converter(self, arg: int, /) -> None: ...

    @property
    def gauss_sigma_cut(self) -> float:
        """Only used if positive."""

    @gauss_sigma_cut.setter
    def gauss_sigma_cut(self, arg: float, /) -> None: ...

    @property
    def in_sobseq(self) -> int: ...

    @in_sobseq.setter
    def in_sobseq(self, arg: int, /) -> None: ...

    @property
    def ix_sobseq(self) -> Int8Array1D: ...

    @ix_sobseq.setter
    def ix_sobseq(self, arg: Sequence[int], /) -> None: ...

    @property
    def x_sobseq(self) -> RealArray1D: ...

    @x_sobseq.setter
    def x_sobseq(self, arg: Sequence[float], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> RandomStateStruct: ...

    def __deepcopy__(self, arg: dict, /) -> RandomStateStruct: ...

    def __eq__(self, arg: RandomStateStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class ResonanceHStruct:
    """Fortran struct: resonance_h_struct"""

    def __init__(self, id: str | None = None, c_val: complex | None = None) -> None: ...

    @property
    def id(self) -> str:
        """6 digit ID. EG: '003100'"""

    @id.setter
    def id(self, arg: str, /) -> None: ...

    @property
    def c_val(self) -> complex:
        """Resonance value"""

    @c_val.setter
    def c_val(self, arg: complex, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> ResonanceHStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> ResonanceHStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> ResonanceHStruct: ...

    def __deepcopy__(self, arg: dict, /) -> ResonanceHStruct: ...

    def __eq__(self, arg: ResonanceHStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class RfEleStruct:
    """Fortran struct: rf_ele_struct"""

    def __init__(self, ds_step: float | None = None) -> None: ...

    @property
    def steps(self) -> RfStairStepStructAlloc1D:
        """Energy stair step array indexed from zero."""

    @property
    def ds_step(self) -> float:
        """length of a stair step."""

    @ds_step.setter
    def ds_step(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> RfEleStruct: ...

    def __deepcopy__(self, arg: dict, /) -> RfEleStruct: ...

    def __eq__(self, arg: RfEleStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class RfStairStepStruct:
    """Fortran struct: rf_stair_step_struct"""

    def __init__(self, E_tot0: float | None = None, E_tot1: float | None = None, p0c: float | None = None, p1c: float | None = None, scale: float | None = None, time: float | None = None, s0: float | None = None, s: float | None = None, ix_step: int | None = None) -> None: ...

    @property
    def E_tot0(self) -> float:
        """Reference energy in the drift region (before the kick point)."""

    @E_tot0.setter
    def E_tot0(self, arg: float, /) -> None: ...

    @property
    def E_tot1(self) -> float:
        """Reference energy after the kick point."""

    @E_tot1.setter
    def E_tot1(self, arg: float, /) -> None: ...

    @property
    def p0c(self) -> float:
        """Reference momentum in the drift region (before the kick point)."""

    @p0c.setter
    def p0c(self, arg: float, /) -> None: ...

    @property
    def p1c(self) -> float:
        """Reference momentum after the kick point."""

    @p1c.setter
    def p1c(self, arg: float, /) -> None: ...

    @property
    def scale(self) -> float:
        """
        Scale for multipole kick at the kick point. Sum over all steps will be 1.
        """

    @scale.setter
    def scale(self, arg: float, /) -> None: ...

    @property
    def time(self) -> float:
        """
        Reference particle time at the kick point with respect to beginning of element.
        """

    @time.setter
    def time(self, arg: float, /) -> None: ...

    @property
    def s0(self) -> float:
        """
        S-position at beginning of drift region relative to the beginning of the element.
        """

    @s0.setter
    def s0(self, arg: float, /) -> None: ...

    @property
    def s(self) -> float:
        """S-position at the kick point relative to the beginning of the element."""

    @s.setter
    def s(self, arg: float, /) -> None: ...

    @property
    def ix_step(self) -> int:
        """Step index in ele%rf%steps(:) array"""

    @ix_step.setter
    def ix_step(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> RfStairStepStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> RfStairStepStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> RfStairStepStruct: ...

    def __deepcopy__(self, arg: dict, /) -> RfStairStepStruct: ...

    def __eq__(self, arg: RfStairStepStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class SeqEleStruct:
    """Fortran struct: seq_ele_struct"""

    def __init__(self, name: str | None = None, tag: str | None = None, slice_start: str | None = None, slice_end: str | None = None, type: int | None = None, ix_ele: int | None = None, ix_arg: int | None = None, rep_count: int | None = None, ele_order_reflect: bool | None = None, ele_orientation: int | None = None) -> None: ...

    @property
    def name(self) -> str:
        """name of element, subline, or sublist"""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def actual_arg(self) -> FCharArray1D: ...

    @property
    def tag(self) -> str:
        """tag name."""

    @tag.setter
    def tag(self, arg: str, /) -> None: ...

    @property
    def slice_start(self) -> str:
        """For 'my_line[start:end]' slice constructs."""

    @slice_start.setter
    def slice_start(self, arg: str, /) -> None: ...

    @property
    def slice_end(self) -> str:
        """For 'my_line[start:end]' slice constructs."""

    @slice_end.setter
    def slice_end(self, arg: str, /) -> None: ...

    @property
    def type(self) -> int:
        """LINE$, REPLACEMENT_LINE$, LIST$, ELEMENT$"""

    @type.setter
    def type(self, arg: int, /) -> None: ...

    @property
    def ix_ele(self) -> int:
        """
        if an element: pointer to ELE array if a line or list: pointer to SEQ array
        """

    @ix_ele.setter
    def ix_ele(self, arg: int, /) -> None: ...

    @property
    def ix_arg(self) -> int:
        """index in arg list (for replacement lines)"""

    @ix_arg.setter
    def ix_arg(self, arg: int, /) -> None: ...

    @property
    def rep_count(self) -> int:
        """how many copies of an element"""

    @rep_count.setter
    def rep_count(self, arg: int, /) -> None: ...

    @property
    def ele_order_reflect(self) -> bool:
        """Travel through ele sequence in reverse order"""

    @ele_order_reflect.setter
    def ele_order_reflect(self, arg: bool, /) -> None: ...

    @property
    def ele_orientation(self) -> int:
        """element has reverse orientation."""

    @ele_orientation.setter
    def ele_orientation(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> SeqEleStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> SeqEleStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> SeqEleStruct: ...

    def __deepcopy__(self, arg: dict, /) -> SeqEleStruct: ...

    def __eq__(self, arg: SeqEleStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class SeqStruct:
    """Fortran struct: seq_struct"""

    def __init__(self, name: str | None = None, type: int | None = None, ix_list: int | None = None, list_upcount: int | None = None, index: int | None = None, file_name: str | None = None, ix_file_line: int | None = None, multipass: bool | None = None, ptc_layout: bool | None = None, active: bool | None = None) -> None: ...

    @property
    def name(self) -> str:
        """name of sequence"""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def ele(self) -> SeqEleStructAlloc1D:
        """Elements in the sequence"""

    @property
    def dummy_arg(self) -> FCharArray1D: ...

    @property
    def corresponding_actual_arg(self) -> FCharArray1D: ...

    @property
    def type(self) -> int:
        """LINE$, REPLACEMENT_LINE$ or LIST$"""

    @type.setter
    def type(self, arg: int, /) -> None: ...

    @property
    def ix_list(self) -> int:
        """Current index for lists"""

    @ix_list.setter
    def ix_list(self, arg: int, /) -> None: ...

    @property
    def list_upcount(self) -> int: ...

    @list_upcount.setter
    def list_upcount(self, arg: int, /) -> None: ...

    @property
    def index(self) -> int:
        """Alphabetical order sorted index"""

    @index.setter
    def index(self, arg: int, /) -> None: ...

    @property
    def file_name(self) -> str:
        """File where sequence is defined"""

    @file_name.setter
    def file_name(self, arg: str, /) -> None: ...

    @property
    def ix_file_line(self) -> int:
        """Line number in file where sequence is defined"""

    @ix_file_line.setter
    def ix_file_line(self, arg: int, /) -> None: ...

    @property
    def multipass(self) -> bool: ...

    @multipass.setter
    def multipass(self, arg: bool, /) -> None: ...

    @property
    def ptc_layout(self) -> bool:
        """Put in separate PTC layout"""

    @ptc_layout.setter
    def ptc_layout(self, arg: bool, /) -> None: ...

    @property
    def active(self) -> bool:
        """Used to prevent infinite loops."""

    @active.setter
    def active(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> SeqStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> SeqStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> SeqStruct: ...

    def __deepcopy__(self, arg: dict, /) -> SeqStruct: ...

    def __eq__(self, arg: SeqStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class SpaceChargeCommonStruct:
    """Fortran struct: space_charge_common_struct"""

    def __init__(self, ds_track_step: float | None = None, dt_track_step: float | None = None, cathode_strength_cutoff: float | None = None, rel_tol_tracking: float | None = None, abs_tol_tracking: float | None = None, beam_chamber_height: float | None = None, lsc_sigma_cutoff: float | None = None, particle_sigma_cutoff: float | None = None, mesh_growth_factor: float | None = None, mesh_shrink_factor: float | None = None, space_charge_mesh_size: Sequence[int] | None = None, csr3d_mesh_size: Sequence[int] | None = None, n_bin: int | None = None, particle_bin_span: int | None = None, n_shield_images: int | None = None, sc_min_in_bin: int | None = None, lsc_kick_transverse_dependence: bool | None = None, debug: bool | None = None, diagnostic_output_file: str | None = None) -> None: ...

    @property
    def ds_track_step(self) -> float:
        """CSR tracking step size"""

    @ds_track_step.setter
    def ds_track_step(self, arg: float, /) -> None: ...

    @property
    def dt_track_step(self) -> float:
        """Time Runge kutta initial step."""

    @dt_track_step.setter
    def dt_track_step(self, arg: float, /) -> None: ...

    @property
    def cathode_strength_cutoff(self) -> float:
        """Cutoff for the cathode field calc."""

    @cathode_strength_cutoff.setter
    def cathode_strength_cutoff(self, arg: float, /) -> None: ...

    @property
    def rel_tol_tracking(self) -> float:
        """Relative tolerance for tracking."""

    @rel_tol_tracking.setter
    def rel_tol_tracking(self, arg: float, /) -> None: ...

    @property
    def abs_tol_tracking(self) -> float:
        """Absolute tolerance for tracking."""

    @abs_tol_tracking.setter
    def abs_tol_tracking(self, arg: float, /) -> None: ...

    @property
    def beam_chamber_height(self) -> float:
        """Used in shielding calculation."""

    @beam_chamber_height.setter
    def beam_chamber_height(self, arg: float, /) -> None: ...

    @property
    def lsc_sigma_cutoff(self) -> float:
        """
        Cutoff for the 1-dim longitudinal SC calc. If a bin sigma is < cutoff * sigma_ave then ignore.
        """

    @lsc_sigma_cutoff.setter
    def lsc_sigma_cutoff(self, arg: float, /) -> None: ...

    @property
    def particle_sigma_cutoff(self) -> float:
        """
        3D SC calc cutoff for particles with (x,y,z) position far from the center. Negative or zero means ignore.
        """

    @particle_sigma_cutoff.setter
    def particle_sigma_cutoff(self, arg: float, /) -> None: ...

    @property
    def mesh_growth_factor(self) -> float:
        """
        Fractional padding when growing SC mesh (default: 10%). Set to 0 for tight-fit (no caching speedup).
        """

    @mesh_growth_factor.setter
    def mesh_growth_factor(self, arg: float, /) -> None: ...

    @property
    def mesh_shrink_factor(self) -> float:
        """
        Fractional threshold for shrinking SC mesh (default: 10%). Mesh shrinks when bunch fills < (1-this) of the mesh range.
        """

    @mesh_shrink_factor.setter
    def mesh_shrink_factor(self, arg: float, /) -> None: ...

    @property
    def space_charge_mesh_size(self) -> IntArray1D:
        """Gird size for fft_3d space charge calc."""

    @space_charge_mesh_size.setter
    def space_charge_mesh_size(self, arg: Sequence[int], /) -> None: ...

    @property
    def csr3d_mesh_size(self) -> IntArray1D:
        """Gird size for CSR."""

    @csr3d_mesh_size.setter
    def csr3d_mesh_size(self, arg: Sequence[int], /) -> None: ...

    @property
    def n_bin(self) -> int:
        """Number of bins used"""

    @n_bin.setter
    def n_bin(self, arg: int, /) -> None: ...

    @property
    def particle_bin_span(self) -> int:
        """Longitudinal particle length / dz_bin"""

    @particle_bin_span.setter
    def particle_bin_span(self, arg: int, /) -> None: ...

    @property
    def n_shield_images(self) -> int:
        """Chamber wall shielding. 0 = no shielding."""

    @n_shield_images.setter
    def n_shield_images(self, arg: int, /) -> None: ...

    @property
    def sc_min_in_bin(self) -> int:
        """Minimum number of particles in a bin for sigmas to be valid."""

    @sc_min_in_bin.setter
    def sc_min_in_bin(self, arg: int, /) -> None: ...

    @property
    def lsc_kick_transverse_dependence(self) -> bool: ...

    @lsc_kick_transverse_dependence.setter
    def lsc_kick_transverse_dependence(self, arg: bool, /) -> None: ...

    @property
    def debug(self) -> bool: ...

    @debug.setter
    def debug(self, arg: bool, /) -> None: ...

    @property
    def diagnostic_output_file(self) -> str:
        """If non-blank write a diagnostic (EG wake) file"""

    @diagnostic_output_file.setter
    def diagnostic_output_file(self, arg: str, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> SpaceChargeCommonStruct: ...

    def __deepcopy__(self, arg: dict, /) -> SpaceChargeCommonStruct: ...

    def __eq__(self, arg: SpaceChargeCommonStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class SpinAxisStruct:
    """Fortran struct: spin_axis_struct"""

    def __init__(self, l: Sequence[float] | None = None, n0: Sequence[float] | None = None, m: Sequence[float] | None = None) -> None: ...

    @property
    def l(self) -> RealArray1D:
        """Transverse axis."""

    @l.setter
    def l(self, arg: Sequence[float], /) -> None: ...

    @property
    def n0(self) -> RealArray1D:
        """Invariant spin axis on closed orbit."""

    @n0.setter
    def n0(self, arg: Sequence[float], /) -> None: ...

    @property
    def m(self) -> RealArray1D:
        """Transverse axis."""

    @m.setter
    def m(self, arg: Sequence[float], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> SpinAxisStruct: ...

    def __deepcopy__(self, arg: dict, /) -> SpinAxisStruct: ...

    def __eq__(self, arg: SpinAxisStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class SpinEigenStruct:
    """Fortran struct: spin_eigen_struct"""

    def __init__(self, vec: Sequence[complex] | None = None, val: complex | None = None) -> None: ...

    @property
    def vec(self) -> ComplexArray1D: ...

    @vec.setter
    def vec(self, arg: Sequence[complex], /) -> None: ...

    @property
    def val(self) -> complex: ...

    @val.setter
    def val(self, arg: complex, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> SpinEigenStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> SpinEigenStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> SpinEigenStruct: ...

    def __deepcopy__(self, arg: dict, /) -> SpinEigenStruct: ...

    def __eq__(self, arg: SpinEigenStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class SpinMatchingStruct:
    """Fortran struct: spin_matching_struct"""

    def __init__(self, axis: SpinAxisStruct | None = None, dn_dpz: Sequence[float] | None = None, alpha: Sequence[float] | None = None, beta: Sequence[float] | None = None, orb0: Sequence[float] | None = None, M_1turn: Sequence[Sequence[float]] | None = None, M_ele: Sequence[Sequence[float]] | None = None, sq_ele: Sequence[float] | None = None, sq_1turn: Sequence[float] | None = None, valid: bool | None = None) -> None: ...

    @property
    def axis(self) -> SpinAxisStruct: ...

    @axis.setter
    def axis(self, arg: SpinAxisStruct, /) -> None: ...

    @property
    def eigen(self) -> SpinEigenStructArray1D: ...

    @property
    def dn_dpz(self) -> RealArray1D:
        """Invariant spin derivative"""

    @dn_dpz.setter
    def dn_dpz(self, arg: Sequence[float], /) -> None: ...

    @property
    def alpha(self) -> RealArray1D:
        """Alpha vector"""

    @alpha.setter
    def alpha(self, arg: Sequence[float], /) -> None: ...

    @property
    def beta(self) -> RealArray1D:
        """Beta vector"""

    @beta.setter
    def beta(self, arg: Sequence[float], /) -> None: ...

    @property
    def orb0(self) -> RealArray1D:
        """Closed orbit"""

    @orb0.setter
    def orb0(self, arg: Sequence[float], /) -> None: ...

    @property
    def M_1turn(self) -> RealArray2D:
        """1-turn matrix"""

    @M_1turn.setter
    def M_1turn(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def M_ele(self) -> RealArray2D:
        """Transfer matrix through element."""

    @M_ele.setter
    def M_ele(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def sq_ele(self) -> RealArray1D: ...

    @sq_ele.setter
    def sq_ele(self, arg: Sequence[float], /) -> None: ...

    @property
    def sq_1turn(self) -> RealArray1D: ...

    @sq_1turn.setter
    def sq_1turn(self, arg: Sequence[float], /) -> None: ...

    @property
    def valid(self) -> bool: ...

    @valid.setter
    def valid(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> SpinMatchingStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> SpinMatchingStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> SpinMatchingStruct: ...

    def __deepcopy__(self, arg: dict, /) -> SpinMatchingStruct: ...

    def __eq__(self, arg: SpinMatchingStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class SpinOrbitMap1Struct:
    """Fortran struct: spin_orbit_map1_struct"""

    def __init__(self, orb_mat: Sequence[Sequence[float]] | None = None, vec0: Sequence[float] | None = None, spin_q: Sequence[Sequence[float]] | None = None) -> None: ...

    @property
    def orb_mat(self) -> RealArray2D:
        """Orbital matrix"""

    @orb_mat.setter
    def orb_mat(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def vec0(self) -> RealArray1D:
        """Orbital 0th order map: r_out = mat6 * r_in + vec0"""

    @vec0.setter
    def vec0(self, arg: Sequence[float], /) -> None: ...

    @property
    def spin_q(self) -> RealArray2D:
        """0th and 1st order quaternion spin map"""

    @spin_q.setter
    def spin_q(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> SpinOrbitMap1StructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> SpinOrbitMap1StructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> SpinOrbitMap1Struct: ...

    def __deepcopy__(self, arg: dict, /) -> SpinOrbitMap1Struct: ...

    def __eq__(self, arg: SpinOrbitMap1Struct, /) -> bool: ...

    def __hash__(self) -> int: ...

class SpinPolarStruct:
    """Fortran struct: spin_polar_struct"""

    def __init__(self, polarization: float | None = None, theta: float | None = None, phi: float | None = None, xi: float | None = None) -> None: ...

    @property
    def polarization(self) -> float: ...

    @polarization.setter
    def polarization(self, arg: float, /) -> None: ...

    @property
    def theta(self) -> float:
        """Spherical coords: Angle from z-axis."""

    @theta.setter
    def theta(self, arg: float, /) -> None: ...

    @property
    def phi(self) -> float:
        """Spherical coords: Angle in (x,y) plane."""

    @phi.setter
    def phi(self, arg: float, /) -> None: ...

    @property
    def xi(self) -> float:
        """Spinor phase angle (See Bmad manual)."""

    @xi.setter
    def xi(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> SpinPolarStruct: ...

    def __deepcopy__(self, arg: dict, /) -> SpinPolarStruct: ...

    def __eq__(self, arg: SpinPolarStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class SplineStruct:
    """Fortran struct: spline_struct"""

    def __init__(self, x0: float | None = None, y0: float | None = None, x1: float | None = None, coef: Sequence[float] | None = None) -> None: ...

    @property
    def x0(self) -> float:
        """Point at start of spline"""

    @x0.setter
    def x0(self, arg: float, /) -> None: ...

    @property
    def y0(self) -> float:
        """Point at start of spline"""

    @y0.setter
    def y0(self, arg: float, /) -> None: ...

    @property
    def x1(self) -> float:
        """Point at end of spline"""

    @x1.setter
    def x1(self, arg: float, /) -> None: ...

    @property
    def coef(self) -> RealArray1D:
        """coefficients for cubic spline"""

    @coef.setter
    def coef(self, arg: Sequence[float], /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> SplineStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> SplineStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> SplineStruct: ...

    def __deepcopy__(self, arg: dict, /) -> SplineStruct: ...

    def __eq__(self, arg: SplineStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class StrIndexStruct:
    """Fortran struct: str_index_struct"""

    def __init__(self, index: Sequence[int] | None = None, n_min: int | None = None, n_max: int | None = None) -> None: ...

    @property
    def name(self) -> VarLengthStringStructAlloc1D:
        """Array of names."""

    @property
    def index(self) -> IntAlloc1D:
        """
        Sorted index for names(:) array. names(an_index(i)) is in alphabetical order.
        """

    @index.setter
    def index(self, arg: Sequence[int], /) -> None: ...

    @property
    def n_min(self) -> int: ...

    @n_min.setter
    def n_min(self, arg: int, /) -> None: ...

    @property
    def n_max(self) -> int:
        """Use only names(n_min:n_max) part of array."""

    @n_max.setter
    def n_max(self, arg: int, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> StrIndexStruct: ...

    def __deepcopy__(self, arg: dict, /) -> StrIndexStruct: ...

    def __eq__(self, arg: StrIndexStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class StrongBeamStruct:
    """Fortran struct: strong_beam_struct"""

    def __init__(self, ix_slice: int | None = None, x_center: float | None = None, y_center: float | None = None, x_sigma: float | None = None, y_sigma: float | None = None, dx: float | None = None, dy: float | None = None) -> None: ...

    @property
    def ix_slice(self) -> int:
        """0 -> at element center and not at slice."""

    @ix_slice.setter
    def ix_slice(self, arg: int, /) -> None: ...

    @property
    def x_center(self) -> float:
        """Strong beam slice center."""

    @x_center.setter
    def x_center(self, arg: float, /) -> None: ...

    @property
    def y_center(self) -> float:
        """Strong beam slice center."""

    @y_center.setter
    def y_center(self, arg: float, /) -> None: ...

    @property
    def x_sigma(self) -> float:
        """Strong beam slice sigma."""

    @x_sigma.setter
    def x_sigma(self, arg: float, /) -> None: ...

    @property
    def y_sigma(self) -> float:
        """Strong beam slice sigma."""

    @y_sigma.setter
    def y_sigma(self, arg: float, /) -> None: ...

    @property
    def dx(self) -> float:
        """Particle - beam slice distance."""

    @dx.setter
    def dx(self, arg: float, /) -> None: ...

    @property
    def dy(self) -> float:
        """Particle - beam slice distance."""

    @dy.setter
    def dy(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> StrongBeamStruct: ...

    def __deepcopy__(self, arg: dict, /) -> StrongBeamStruct: ...

    def __eq__(self, arg: StrongBeamStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class SummationRdtStruct:
    """Fortran struct: summation_rdt_struct"""

    def __init__(self, h11001: complex | None = None, h00111: complex | None = None, h20001: complex | None = None, h00201: complex | None = None, h10002: complex | None = None, h21000: complex | None = None, h30000: complex | None = None, h10110: complex | None = None, h10020: complex | None = None, h10200: complex | None = None, h31000: complex | None = None, h40000: complex | None = None, h20110: complex | None = None, h11200: complex | None = None, h20020: complex | None = None, h20200: complex | None = None, h00310: complex | None = None, h00400: complex | None = None, h22000: complex | None = None, h00220: complex | None = None, h11110: complex | None = None) -> None: ...

    @property
    def h11001(self) -> complex: ...

    @h11001.setter
    def h11001(self, arg: complex, /) -> None: ...

    @property
    def h00111(self) -> complex: ...

    @h00111.setter
    def h00111(self, arg: complex, /) -> None: ...

    @property
    def h20001(self) -> complex: ...

    @h20001.setter
    def h20001(self, arg: complex, /) -> None: ...

    @property
    def h00201(self) -> complex: ...

    @h00201.setter
    def h00201(self, arg: complex, /) -> None: ...

    @property
    def h10002(self) -> complex: ...

    @h10002.setter
    def h10002(self, arg: complex, /) -> None: ...

    @property
    def h21000(self) -> complex: ...

    @h21000.setter
    def h21000(self, arg: complex, /) -> None: ...

    @property
    def h30000(self) -> complex: ...

    @h30000.setter
    def h30000(self, arg: complex, /) -> None: ...

    @property
    def h10110(self) -> complex: ...

    @h10110.setter
    def h10110(self, arg: complex, /) -> None: ...

    @property
    def h10020(self) -> complex: ...

    @h10020.setter
    def h10020(self, arg: complex, /) -> None: ...

    @property
    def h10200(self) -> complex:
        """2nd order in K2 moments"""

    @h10200.setter
    def h10200(self, arg: complex, /) -> None: ...

    @property
    def h31000(self) -> complex: ...

    @h31000.setter
    def h31000(self, arg: complex, /) -> None: ...

    @property
    def h40000(self) -> complex: ...

    @h40000.setter
    def h40000(self, arg: complex, /) -> None: ...

    @property
    def h20110(self) -> complex: ...

    @h20110.setter
    def h20110(self, arg: complex, /) -> None: ...

    @property
    def h11200(self) -> complex: ...

    @h11200.setter
    def h11200(self, arg: complex, /) -> None: ...

    @property
    def h20020(self) -> complex: ...

    @h20020.setter
    def h20020(self, arg: complex, /) -> None: ...

    @property
    def h20200(self) -> complex: ...

    @h20200.setter
    def h20200(self, arg: complex, /) -> None: ...

    @property
    def h00310(self) -> complex: ...

    @h00310.setter
    def h00310(self, arg: complex, /) -> None: ...

    @property
    def h00400(self) -> complex: ...

    @h00400.setter
    def h00400(self, arg: complex, /) -> None: ...

    @property
    def h22000(self) -> complex: ...

    @h22000.setter
    def h22000(self, arg: complex, /) -> None: ...

    @property
    def h00220(self) -> complex: ...

    @h00220.setter
    def h00220(self, arg: complex, /) -> None: ...

    @property
    def h11110(self) -> complex: ...

    @h11110.setter
    def h11110(self, arg: complex, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> SummationRdtStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> SummationRdtStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> SummationRdtStruct: ...

    def __deepcopy__(self, arg: dict, /) -> SummationRdtStruct: ...

    def __eq__(self, arg: SummationRdtStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class SurfaceCurvatureStruct:
    """Fortran struct: surface_curvature_struct"""

    def __init__(self, xy: Sequence[Sequence[float]] | None = None, spherical: float | None = None, elliptical: Sequence[float] | None = None, has_curvature: bool | None = None) -> None: ...

    @property
    def xy(self) -> RealArray2D: ...

    @xy.setter
    def xy(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def spherical(self) -> float: ...

    @spherical.setter
    def spherical(self, arg: float, /) -> None: ...

    @property
    def elliptical(self) -> RealArray1D:
        """Total curvature = elliptical + spherical"""

    @elliptical.setter
    def elliptical(self, arg: Sequence[float], /) -> None: ...

    @property
    def has_curvature(self) -> bool:
        """Dependent var. Will be set by Bmad"""

    @has_curvature.setter
    def has_curvature(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> SurfaceCurvatureStruct: ...

    def __deepcopy__(self, arg: dict, /) -> SurfaceCurvatureStruct: ...

    def __eq__(self, arg: SurfaceCurvatureStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class SurfaceDisplacementPtStruct:
    """Fortran struct: surface_displacement_pt_struct"""

    def __init__(self, x0: float | None = None, y0: float | None = None, z0: float | None = None, dz_dx: float | None = None, dz_dy: float | None = None, d2z_dxdy: float | None = None) -> None: ...

    @property
    def x0(self) -> float:
        """Position at center"""

    @x0.setter
    def x0(self, arg: float, /) -> None: ...

    @property
    def y0(self) -> float:
        """Position at center"""

    @y0.setter
    def y0(self, arg: float, /) -> None: ...

    @property
    def z0(self) -> float: ...

    @z0.setter
    def z0(self, arg: float, /) -> None: ...

    @property
    def dz_dx(self) -> float: ...

    @dz_dx.setter
    def dz_dx(self, arg: float, /) -> None: ...

    @property
    def dz_dy(self) -> float: ...

    @dz_dy.setter
    def dz_dy(self, arg: float, /) -> None: ...

    @property
    def d2z_dxdy(self) -> float: ...

    @d2z_dxdy.setter
    def d2z_dxdy(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> SurfaceDisplacementPtStruct: ...

    def __deepcopy__(self, arg: dict, /) -> SurfaceDisplacementPtStruct: ...

    def __eq__(self, arg: SurfaceDisplacementPtStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class SurfaceDisplacementStruct:
    """Fortran struct: surface_displacement_struct"""

    def __init__(self, active: bool | None = None, dr: Sequence[float] | None = None, r0: Sequence[float] | None = None) -> None: ...

    @property
    def active(self) -> bool: ...

    @active.setter
    def active(self, arg: bool, /) -> None: ...

    @property
    def dr(self) -> RealArray1D: ...

    @dr.setter
    def dr(self, arg: Sequence[float], /) -> None: ...

    @property
    def r0(self) -> RealArray1D: ...

    @r0.setter
    def r0(self, arg: Sequence[float], /) -> None: ...

    @property
    def pt(self) -> SurfaceDisplacementPtStructArray2D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> SurfaceDisplacementStruct: ...

    def __deepcopy__(self, arg: dict, /) -> SurfaceDisplacementStruct: ...

    def __eq__(self, arg: SurfaceDisplacementStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class SurfaceHMisalignPtStruct:
    """Fortran struct: surface_h_misalign_pt_struct"""

    def __init__(self, x0: float | None = None, y0: float | None = None, rot_y: float | None = None, rot_t: float | None = None, rot_y_rms: float | None = None, rot_t_rms: float | None = None) -> None: ...

    @property
    def x0(self) -> float:
        """Position at center"""

    @x0.setter
    def x0(self, arg: float, /) -> None: ...

    @property
    def y0(self) -> float:
        """Position at center"""

    @y0.setter
    def y0(self, arg: float, /) -> None: ...

    @property
    def rot_y(self) -> float:
        """rot_t = x-rotation for Bragg and z-rotation for Laue."""

    @rot_y.setter
    def rot_y(self, arg: float, /) -> None: ...

    @property
    def rot_t(self) -> float:
        """rot_t = x-rotation for Bragg and z-rotation for Laue."""

    @rot_t.setter
    def rot_t(self, arg: float, /) -> None: ...

    @property
    def rot_y_rms(self) -> float:
        """rot_t = x-rotation for Bragg and z-rotation for Laue."""

    @rot_y_rms.setter
    def rot_y_rms(self, arg: float, /) -> None: ...

    @property
    def rot_t_rms(self) -> float:
        """rot_t = x-rotation for Bragg and z-rotation for Laue."""

    @rot_t_rms.setter
    def rot_t_rms(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> SurfaceHMisalignPtStruct: ...

    def __deepcopy__(self, arg: dict, /) -> SurfaceHMisalignPtStruct: ...

    def __eq__(self, arg: SurfaceHMisalignPtStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class SurfaceHMisalignStruct:
    """Fortran struct: surface_h_misalign_struct"""

    def __init__(self, active: bool | None = None, dr: Sequence[float] | None = None, r0: Sequence[float] | None = None) -> None: ...

    @property
    def active(self) -> bool: ...

    @active.setter
    def active(self, arg: bool, /) -> None: ...

    @property
    def dr(self) -> RealArray1D: ...

    @dr.setter
    def dr(self, arg: Sequence[float], /) -> None: ...

    @property
    def r0(self) -> RealArray1D: ...

    @r0.setter
    def r0(self, arg: Sequence[float], /) -> None: ...

    @property
    def pt(self) -> SurfaceHMisalignPtStructArray2D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> SurfaceHMisalignStruct: ...

    def __deepcopy__(self, arg: dict, /) -> SurfaceHMisalignStruct: ...

    def __eq__(self, arg: SurfaceHMisalignStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class SurfaceSegmentedPtStruct:
    """Fortran struct: surface_segmented_pt_struct"""

    def __init__(self, x0: float | None = None, y0: float | None = None, z0: float | None = None, dz_dx: float | None = None, dz_dy: float | None = None) -> None: ...

    @property
    def x0(self) -> float:
        """Position at center"""

    @x0.setter
    def x0(self, arg: float, /) -> None: ...

    @property
    def y0(self) -> float:
        """Position at center"""

    @y0.setter
    def y0(self, arg: float, /) -> None: ...

    @property
    def z0(self) -> float:
        """Position at center"""

    @z0.setter
    def z0(self, arg: float, /) -> None: ...

    @property
    def dz_dx(self) -> float:
        """Slope at center"""

    @dz_dx.setter
    def dz_dx(self, arg: float, /) -> None: ...

    @property
    def dz_dy(self) -> float:
        """Slope at center"""

    @dz_dy.setter
    def dz_dy(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> SurfaceSegmentedPtStruct: ...

    def __deepcopy__(self, arg: dict, /) -> SurfaceSegmentedPtStruct: ...

    def __eq__(self, arg: SurfaceSegmentedPtStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class SurfaceSegmentedStruct:
    """Fortran struct: surface_segmented_struct"""

    def __init__(self, active: bool | None = None, dr: Sequence[float] | None = None, r0: Sequence[float] | None = None) -> None: ...

    @property
    def active(self) -> bool: ...

    @active.setter
    def active(self, arg: bool, /) -> None: ...

    @property
    def dr(self) -> RealArray1D: ...

    @dr.setter
    def dr(self, arg: Sequence[float], /) -> None: ...

    @property
    def r0(self) -> RealArray1D: ...

    @r0.setter
    def r0(self, arg: Sequence[float], /) -> None: ...

    @property
    def pt(self) -> SurfaceSegmentedPtStructArray2D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> SurfaceSegmentedStruct: ...

    def __deepcopy__(self, arg: dict, /) -> SurfaceSegmentedStruct: ...

    def __eq__(self, arg: SurfaceSegmentedStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoAliasStruct:
    """Fortran struct: tao_alias_struct"""

    def __init__(self, name: str | None = None, expanded_str: str | None = None) -> None: ...

    @property
    def name(self) -> str: ...

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def expanded_str(self) -> str: ...

    @expanded_str.setter
    def expanded_str(self, arg: str, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoAliasStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoAliasStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoAliasStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoAliasStruct: ...

    def __eq__(self, arg: TaoAliasStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoBeamBranchStruct:
    """Fortran struct: tao_beam_branch_struct"""

    def __init__(self, beam_at_start: BeamStruct | None = None, beam_init: BeamInitStruct | None = None, beam_init_used: BeamInitStruct | None = None, init_starting_distribution: bool | None = None, track_start: str | None = None, track_end: str | None = None, ix_branch: int | None = None, ix_track_start: int | None = None, ix_track_end: int | None = None) -> None: ...

    @property
    def beam_at_start(self) -> BeamStruct:
        """Initial beam"""

    @beam_at_start.setter
    def beam_at_start(self, arg: BeamStruct, /) -> None: ...

    @property
    def beam_init(self) -> BeamInitStruct:
        """User set beam distrubution at track start."""

    @beam_init.setter
    def beam_init(self, arg: BeamInitStruct, /) -> None: ...

    @property
    def beam_init_used(self) -> BeamInitStruct:
        """beam distribution with emit values set."""

    @beam_init_used.setter
    def beam_init_used(self, arg: BeamInitStruct, /) -> None: ...

    @property
    def init_starting_distribution(self) -> bool:
        """Init beam"""

    @init_starting_distribution.setter
    def init_starting_distribution(self, arg: bool, /) -> None: ...

    @property
    def track_start(self) -> str:
        """Tracking start element."""

    @track_start.setter
    def track_start(self, arg: str, /) -> None: ...

    @property
    def track_end(self) -> str: ...

    @track_end.setter
    def track_end(self, arg: str, /) -> None: ...

    @property
    def ix_branch(self) -> int:
        """
        Branch tracked. If track_start or track_end is a lord, ix_track_start/end index will be a index of slave.
        """

    @ix_branch.setter
    def ix_branch(self, arg: int, /) -> None: ...

    @property
    def ix_track_start(self) -> int:
        """Element track start index."""

    @ix_track_start.setter
    def ix_track_start(self, arg: int, /) -> None: ...

    @property
    def ix_track_end(self) -> int:
        """Element track end index"""

    @ix_track_end.setter
    def ix_track_end(self, arg: int, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoBeamBranchStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoBeamBranchStruct: ...

    def __eq__(self, arg: TaoBeamBranchStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoBeamUniStruct:
    """Fortran struct: tao_beam_uni_struct"""

    def __init__(self, saved_at: str | None = None, dump_file: str | None = None, dump_at: str | None = None, track_beam_in_universe: bool | None = None, always_reinit: bool | None = None) -> None: ...

    @property
    def saved_at(self) -> str: ...

    @saved_at.setter
    def saved_at(self, arg: str, /) -> None: ...

    @property
    def dump_file(self) -> str: ...

    @dump_file.setter
    def dump_file(self, arg: str, /) -> None: ...

    @property
    def dump_at(self) -> str: ...

    @dump_at.setter
    def dump_at(self, arg: str, /) -> None: ...

    @property
    def track_beam_in_universe(self) -> bool:
        """Beam tracking enabled in this universe?"""

    @track_beam_in_universe.setter
    def track_beam_in_universe(self, arg: bool, /) -> None: ...

    @property
    def always_reinit(self) -> bool: ...

    @always_reinit.setter
    def always_reinit(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoBeamUniStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoBeamUniStruct: ...

    def __eq__(self, arg: TaoBeamUniStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoBuildingWallOrientationStruct:
    """Fortran struct: tao_building_wall_orientation_struct"""

    def __init__(self, theta: float | None = None, x_offset: float | None = None, z_offset: float | None = None) -> None: ...

    @property
    def theta(self) -> float: ...

    @theta.setter
    def theta(self, arg: float, /) -> None: ...

    @property
    def x_offset(self) -> float: ...

    @x_offset.setter
    def x_offset(self, arg: float, /) -> None: ...

    @property
    def z_offset(self) -> float: ...

    @z_offset.setter
    def z_offset(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoBuildingWallOrientationStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoBuildingWallOrientationStruct: ...

    def __eq__(self, arg: TaoBuildingWallOrientationStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoBuildingWallPointStruct:
    """Fortran struct: tao_building_wall_point_struct"""

    def __init__(self, z: float | None = None, x: float | None = None, radius: float | None = None, z_center: float | None = None, x_center: float | None = None) -> None: ...

    @property
    def z(self) -> float:
        """Global floor position"""

    @z.setter
    def z(self, arg: float, /) -> None: ...

    @property
    def x(self) -> float:
        """Global floor position"""

    @x.setter
    def x(self, arg: float, /) -> None: ...

    @property
    def radius(self) -> float:
        """Arc radius. +r -> CW rotation, same as bends."""

    @radius.setter
    def radius(self, arg: float, /) -> None: ...

    @property
    def z_center(self) -> float:
        """Arc center."""

    @z_center.setter
    def z_center(self, arg: float, /) -> None: ...

    @property
    def x_center(self) -> float:
        """Arc center."""

    @x_center.setter
    def x_center(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoBuildingWallPointStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoBuildingWallPointStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoBuildingWallPointStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoBuildingWallPointStruct: ...

    def __eq__(self, arg: TaoBuildingWallPointStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoBuildingWallSectionStruct:
    """Fortran struct: tao_building_wall_section_struct"""

    def __init__(self, name: str | None = None, constraint: str | None = None) -> None: ...

    @property
    def name(self) -> str: ...

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def constraint(self) -> str:
        """'left_side' or 'right_side' constraint."""

    @constraint.setter
    def constraint(self, arg: str, /) -> None: ...

    @property
    def point(self) -> TaoBuildingWallPointStructAlloc1D: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoBuildingWallSectionStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoBuildingWallSectionStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoBuildingWallSectionStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoBuildingWallSectionStruct: ...

    def __eq__(self, arg: TaoBuildingWallSectionStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoBuildingWallStruct:
    """Fortran struct: tao_building_wall_struct"""

    def __init__(self, orientation: TaoBuildingWallOrientationStruct | None = None) -> None: ...

    @property
    def orientation(self) -> TaoBuildingWallOrientationStruct: ...

    @orientation.setter
    def orientation(self, arg: TaoBuildingWallOrientationStruct, /) -> None: ...

    @property
    def section(self) -> TaoBuildingWallSectionStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoBuildingWallStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoBuildingWallStruct: ...

    def __eq__(self, arg: TaoBuildingWallStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoCmdHistoryStruct:
    """Fortran struct: tao_cmd_history_struct"""

    def __init__(self, cmd: str | None = None, ix: int | None = None) -> None: ...

    @property
    def cmd(self) -> str:
        """The command"""

    @cmd.setter
    def cmd(self, arg: str, /) -> None: ...

    @property
    def ix(self) -> int:
        """
        Command index (1st command has ix = 1, etc.) Note: Commands from command files will be assigned an index.
        """

    @ix.setter
    def ix(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoCmdHistoryStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoCmdHistoryStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoCmdHistoryStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoCmdHistoryStruct: ...

    def __eq__(self, arg: TaoCmdHistoryStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoCommandFileStruct:
    """Fortran struct: tao_command_file_struct"""

    def __init__(self, full_name: str | None = None, dir: str | None = None, ix_unit: int | None = None, quiet: str | None = None, paused: bool | None = None, n_line: int | None = None, reset_at_end: bool | None = None, lattice_calc_save: bool | None = None, plot_save: bool | None = None, multi_cmd: str | None = None) -> None: ...

    @property
    def full_name(self) -> str: ...

    @full_name.setter
    def full_name(self, arg: str, /) -> None: ...

    @property
    def dir(self) -> str: ...

    @dir.setter
    def dir(self, arg: str, /) -> None: ...

    @property
    def ix_unit(self) -> int: ...

    @ix_unit.setter
    def ix_unit(self, arg: int, /) -> None: ...

    @property
    def cmd_arg(self) -> FCharArray1D:
        """Command file arguments."""

    @property
    def quiet(self) -> str: ...

    @quiet.setter
    def quiet(self, arg: str, /) -> None: ...

    @property
    def paused(self) -> bool:
        """Is the command file paused?"""

    @paused.setter
    def paused(self, arg: bool, /) -> None: ...

    @property
    def n_line(self) -> int:
        """Current line number"""

    @n_line.setter
    def n_line(self, arg: int, /) -> None: ...

    @property
    def reset_at_end(self) -> bool:
        """Reset lattice_calc_on and plot_on at end of file?"""

    @reset_at_end.setter
    def reset_at_end(self, arg: bool, /) -> None: ...

    @property
    def lattice_calc_save(self) -> bool: ...

    @lattice_calc_save.setter
    def lattice_calc_save(self, arg: bool, /) -> None: ...

    @property
    def plot_save(self) -> bool: ...

    @plot_save.setter
    def plot_save(self, arg: bool, /) -> None: ...

    @property
    def multi_cmd(self) -> str:
        """Commands not yet executed when there are mulitple commands on a line"""

    @multi_cmd.setter
    def multi_cmd(self, arg: str, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoCommandFileStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoCommandFileStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoCommandFileStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoCommandFileStruct: ...

    def __eq__(self, arg: TaoCommandFileStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoCommonStruct:
    """Fortran struct: tao_common_struct"""

    def __init__(self, covar: Sequence[Sequence[float]] | None = None, alpha: Sequence[Sequence[float]] | None = None, dummy_target: float | None = None, n_alias: int | None = None, cmd_file_level: int | None = None, ix_key_bank: int | None = None, ix_history: int | None = None, n_history: int | None = None, lev_loop: int | None = None, n_err_messages_printed: int | None = None, n_universes: int | None = None, ix_beam_track_active_element: int | None = None, cmd_file_paused: bool | None = None, use_cmd_here: bool | None = None, cmd_from_cmd_file: bool | None = None, use_saved_beam_in_tracking: bool | None = None, single_mode: bool | None = None, combine_consecutive_elements_of_like_name: bool | None = None, have_tracked_beam: bool | None = None, init_plot_needed: bool | None = None, init_beam: bool | None = None, init_var: bool | None = None, init_read_lat_info: bool | None = None, optimizer_running: bool | None = None, have_datums_using_expressions: bool | None = None, print_to_terminal: bool | None = None, lattice_calc_done: bool | None = None, add_measurement_noise: bool | None = None, is_err_message_printed: Sequence[bool] | None = None, command_arg_has_been_executed: bool | None = None, all_merit_weights_positive: bool | None = None, multi_turn_orbit_is_plotted: bool | None = None, force_chrom_calc: bool | None = None, force_rad_int_calc: bool | None = None, rad_int_ri_calc_on: bool | None = None, rad_int_6d_calc_on: bool | None = None, single_mode_buffer: str | None = None, cmd: str | None = None) -> None: ...

    @property
    def alias(self) -> TaoAliasStructArray1D: ...

    @property
    def key(self) -> TaoAliasStructArray1D: ...

    @property
    def cmd_file(self) -> TaoCommandFileStructAlloc1D: ...

    @property
    def symbolic_num(self) -> NamedNumberStructAlloc1D:
        """Named numbers"""

    @property
    def plot_place_buffer(self) -> TaoPlotRegionStructAlloc1D:
        """Used when %external_plotting is on."""

    @property
    def do_loop(self) -> DoLoopStructAlloc1D: ...

    @property
    def covar(self) -> RealArray2D: ...

    @covar.setter
    def covar(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def alpha(self) -> RealArray2D: ...

    @alpha.setter
    def alpha(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def dummy_target(self) -> float:
        """Dummy varaible"""

    @dummy_target.setter
    def dummy_target(self, arg: float, /) -> None: ...

    @property
    def n_alias(self) -> int: ...

    @n_alias.setter
    def n_alias(self, arg: int, /) -> None: ...

    @property
    def cmd_file_level(self) -> int:
        """For nested command files. 0 -> no command file."""

    @cmd_file_level.setter
    def cmd_file_level(self, arg: int, /) -> None: ...

    @property
    def ix_key_bank(self) -> int:
        """For single mode."""

    @ix_key_bank.setter
    def ix_key_bank(self, arg: int, /) -> None: ...

    @property
    def ix_history(self) -> int:
        """Index to latest command in the history circular buffer."""

    @ix_history.setter
    def ix_history(self, arg: int, /) -> None: ...

    @property
    def n_history(self) -> int:
        """Number of commands issued from beginning of starting Tao."""

    @n_history.setter
    def n_history(self, arg: int, /) -> None: ...

    @property
    def lev_loop(self) -> int:
        """in do loop nest level"""

    @lev_loop.setter
    def lev_loop(self, arg: int, /) -> None: ...

    @property
    def n_err_messages_printed(self) -> int:
        """Used by tao_set_invalid to limit number of messages."""

    @n_err_messages_printed.setter
    def n_err_messages_printed(self, arg: int, /) -> None: ...

    @property
    def n_universes(self) -> int: ...

    @n_universes.setter
    def n_universes(self, arg: int, /) -> None: ...

    @property
    def ix_beam_track_active_element(self) -> int:
        """Element being tracked through `tao_beam_track`."""

    @ix_beam_track_active_element.setter
    def ix_beam_track_active_element(self, arg: int, /) -> None: ...

    @property
    def cmd_file_paused(self) -> bool: ...

    @cmd_file_paused.setter
    def cmd_file_paused(self, arg: bool, /) -> None: ...

    @property
    def use_cmd_here(self) -> bool:
        """Used for commands recalled from the cmd history stack"""

    @use_cmd_here.setter
    def use_cmd_here(self, arg: bool, /) -> None: ...

    @property
    def cmd_from_cmd_file(self) -> bool:
        """was command from a command file?"""

    @cmd_from_cmd_file.setter
    def cmd_from_cmd_file(self, arg: bool, /) -> None: ...

    @property
    def use_saved_beam_in_tracking(self) -> bool: ...

    @use_saved_beam_in_tracking.setter
    def use_saved_beam_in_tracking(self, arg: bool, /) -> None: ...

    @property
    def single_mode(self) -> bool: ...

    @single_mode.setter
    def single_mode(self, arg: bool, /) -> None: ...

    @property
    def combine_consecutive_elements_of_like_name(self) -> bool: ...

    @combine_consecutive_elements_of_like_name.setter
    def combine_consecutive_elements_of_like_name(self, arg: bool, /) -> None: ...

    @property
    def have_tracked_beam(self) -> bool:
        """Used to catch error when beam plotting without having tracked a beam."""

    @have_tracked_beam.setter
    def have_tracked_beam(self, arg: bool, /) -> None: ...

    @property
    def init_plot_needed(self) -> bool:
        """reinitialize plotting?"""

    @init_plot_needed.setter
    def init_plot_needed(self, arg: bool, /) -> None: ...

    @property
    def init_beam(self) -> bool:
        """Used by custom programs to control Tao init"""

    @init_beam.setter
    def init_beam(self, arg: bool, /) -> None: ...

    @property
    def init_var(self) -> bool:
        """Used by custom programs to control Tao init"""

    @init_var.setter
    def init_var(self, arg: bool, /) -> None: ...

    @property
    def init_read_lat_info(self) -> bool:
        """Used by custom programs to control Tao init"""

    @init_read_lat_info.setter
    def init_read_lat_info(self, arg: bool, /) -> None: ...

    @property
    def optimizer_running(self) -> bool: ...

    @optimizer_running.setter
    def optimizer_running(self, arg: bool, /) -> None: ...

    @property
    def have_datums_using_expressions(self) -> bool: ...

    @have_datums_using_expressions.setter
    def have_datums_using_expressions(self, arg: bool, /) -> None: ...

    @property
    def print_to_terminal(self) -> bool:
        """Print command prompt to the terminal? For use with GUIs."""

    @print_to_terminal.setter
    def print_to_terminal(self, arg: bool, /) -> None: ...

    @property
    def lattice_calc_done(self) -> bool:
        """Used by GUI for deciding when to refresh."""

    @lattice_calc_done.setter
    def lattice_calc_done(self, arg: bool, /) -> None: ...

    @property
    def add_measurement_noise(self) -> bool:
        """Turn off to take data derivatives."""

    @add_measurement_noise.setter
    def add_measurement_noise(self, arg: bool, /) -> None: ...

    @property
    def is_err_message_printed(self) -> BoolArray1D:
        """Used by tao_set_invalid"""

    @is_err_message_printed.setter
    def is_err_message_printed(self, arg: Sequence[bool], /) -> None: ...

    @property
    def command_arg_has_been_executed(self) -> bool:
        """Has the -command command line argument been executed?"""

    @command_arg_has_been_executed.setter
    def command_arg_has_been_executed(self, arg: bool, /) -> None: ...

    @property
    def all_merit_weights_positive(self) -> bool: ...

    @all_merit_weights_positive.setter
    def all_merit_weights_positive(self, arg: bool, /) -> None: ...

    @property
    def multi_turn_orbit_is_plotted(self) -> bool:
        """Is a multi_turn_orbit being plotted?"""

    @multi_turn_orbit_is_plotted.setter
    def multi_turn_orbit_is_plotted(self, arg: bool, /) -> None: ...

    @property
    def force_chrom_calc(self) -> bool:
        """Used by a routine to force a single chromaticity calculation."""

    @force_chrom_calc.setter
    def force_chrom_calc(self, arg: bool, /) -> None: ...

    @property
    def force_rad_int_calc(self) -> bool:
        """Used by a routine to force a single radiation integrals calculation"""

    @force_rad_int_calc.setter
    def force_rad_int_calc(self, arg: bool, /) -> None: ...

    @property
    def rad_int_ri_calc_on(self) -> bool:
        """'Classical' radiation integrals calculation on/off."""

    @rad_int_ri_calc_on.setter
    def rad_int_ri_calc_on(self, arg: bool, /) -> None: ...

    @property
    def rad_int_6d_calc_on(self) -> bool:
        """6D Radiation integrals calculation on/off."""

    @rad_int_6d_calc_on.setter
    def rad_int_6d_calc_on(self, arg: bool, /) -> None: ...

    @property
    def valid_plot_who(self) -> FCharArray1D:
        """model, base, ref etc..."""

    @property
    def single_mode_buffer(self) -> str: ...

    @single_mode_buffer.setter
    def single_mode_buffer(self, arg: str, /) -> None: ...

    @property
    def cmd(self) -> str:
        """Used for the cmd history"""

    @cmd.setter
    def cmd(self, arg: str, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoCommonStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoCommonStruct: ...

    def __eq__(self, arg: TaoCommonStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoCurveArrayStruct:
    """Fortran struct: tao_curve_array_struct"""

    def __init__(self, c: TaoCurveStruct | None = None) -> None: ...

    @property
    def c(self) -> TaoCurveStruct | None: ...

    @c.setter
    def c(self, arg: TaoCurveStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoCurveArrayStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoCurveArrayStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoCurveArrayStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoCurveArrayStruct: ...

    def __eq__(self, arg: TaoCurveArrayStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoCurveColorStruct:
    """Fortran struct: tao_curve_color_struct"""

    def __init__(self, data_type: str | None = None, is_on: bool | None = None, min: float | None = None, max: float | None = None, autoscale: bool | None = None) -> None: ...

    @property
    def data_type(self) -> str:
        """Datum type to use for z-axis."""

    @data_type.setter
    def data_type(self, arg: str, /) -> None: ...

    @property
    def is_on(self) -> bool:
        """On/Off"""

    @is_on.setter
    def is_on(self, arg: bool, /) -> None: ...

    @property
    def min(self) -> float:
        """Min and max values for mapping z-axis to color."""

    @min.setter
    def min(self, arg: float, /) -> None: ...

    @property
    def max(self) -> float:
        """Min and max values for mapping z-axis to color."""

    @max.setter
    def max(self, arg: float, /) -> None: ...

    @property
    def autoscale(self) -> bool:
        """Set %min, %max automatically to the limits of %data_type"""

    @autoscale.setter
    def autoscale(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoCurveColorStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoCurveColorStruct: ...

    def __eq__(self, arg: TaoCurveColorStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoCurveOrbitStruct:
    """Fortran struct: tao_curve_orbit_struct"""

    def __init__(self, x: float | None = None, y: float | None = None, t: float | None = None) -> None: ...

    @property
    def x(self) -> float:
        """Transverse offset"""

    @x.setter
    def x(self, arg: float, /) -> None: ...

    @property
    def y(self) -> float:
        """Transverse offset"""

    @y.setter
    def y(self, arg: float, /) -> None: ...

    @property
    def t(self) -> float:
        """Time"""

    @t.setter
    def t(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoCurveOrbitStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoCurveOrbitStruct: ...

    def __eq__(self, arg: TaoCurveOrbitStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoCurveStruct:
    """Fortran struct: tao_curve_struct"""

    def __init__(self, name: str | None = None, data_source: str | None = None, data_index: str | None = None, data_type_x: str | None = None, data_type: str | None = None, ele_ref_name: str | None = None, legend_text: str | None = None, message_text: str | None = None, component: str | None = None, why_invalid: str | None = None, g: TaoGraphStruct | None = None, hist: TaoHistogramStruct | None = None, z_color: TaoCurveColorStruct | None = None, x_line: Sequence[float] | None = None, y_line: Sequence[float] | None = None, y2_line: Sequence[float] | None = None, ix_line: Sequence[int] | None = None, x_symb: Sequence[float] | None = None, y_symb: Sequence[float] | None = None, z_symb: Sequence[float] | None = None, err_symb: Sequence[float] | None = None, symb_size: Sequence[float] | None = None, ix_symb: Sequence[int] | None = None, y_axis_scale_factor: float | None = None, line: QpLineStruct | None = None, symbol: QpSymbolStruct | None = None, orbit: TaoCurveOrbitStruct | None = None, ix_universe: int | None = None, symbol_every: int | None = None, ix_branch: int | None = None, ix_bunch: int | None = None, n_turn: int | None = None, use_y2: bool | None = None, draw_line: bool | None = None, draw_symbols: bool | None = None, draw_symbol_index: bool | None = None, draw_error_bars: bool | None = None, smooth_line_calc: bool | None = None, valid: bool | None = None) -> None: ...

    @property
    def name(self) -> str:
        """Name identifying the curve."""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def data_source(self) -> str:
        """'lat', 'beam', 'data' (deprecated: 'dat'), 'var', 'multi_turn_orbit'"""

    @data_source.setter
    def data_source(self, arg: str, /) -> None: ...

    @property
    def data_index(self) -> str:
        """Used for calculating %ix_symb(:)."""

    @data_index.setter
    def data_index(self, arg: str, /) -> None: ...

    @property
    def data_type_x(self) -> str:
        """Used for data slices and phase space plots."""

    @data_type_x.setter
    def data_type_x(self, arg: str, /) -> None: ...

    @property
    def data_type(self) -> str:
        """'orbit.x', etc."""

    @data_type.setter
    def data_type(self, arg: str, /) -> None: ...

    @property
    def ele_ref_name(self) -> str:
        """Reference element."""

    @ele_ref_name.setter
    def ele_ref_name(self, arg: str, /) -> None: ...

    @property
    def legend_text(self) -> str:
        """String to draw in a curve legend."""

    @legend_text.setter
    def legend_text(self, arg: str, /) -> None: ...

    @property
    def message_text(self) -> str:
        """Informational message to draw with graph."""

    @message_text.setter
    def message_text(self, arg: str, /) -> None: ...

    @property
    def component(self) -> str:
        """Who to plot. Eg: 'meas - design'"""

    @component.setter
    def component(self, arg: str, /) -> None: ...

    @property
    def why_invalid(self) -> str:
        """Informative string to print."""

    @why_invalid.setter
    def why_invalid(self, arg: str, /) -> None: ...

    @property
    def g(self) -> TaoGraphStruct | None:
        """pointer to parent graph"""

    @g.setter
    def g(self, arg: TaoGraphStruct, /) -> None: ...

    @property
    def hist(self) -> TaoHistogramStruct: ...

    @hist.setter
    def hist(self, arg: TaoHistogramStruct, /) -> None: ...

    @property
    def z_color(self) -> TaoCurveColorStruct: ...

    @z_color.setter
    def z_color(self, arg: TaoCurveColorStruct, /) -> None: ...

    @property
    def x_line(self) -> RealAlloc1D:
        """Coords for drawing a curve"""

    @x_line.setter
    def x_line(self, arg: Sequence[float], /) -> None: ...

    @property
    def y_line(self) -> RealAlloc1D: ...

    @y_line.setter
    def y_line(self, arg: Sequence[float], /) -> None: ...

    @property
    def y2_line(self) -> RealAlloc1D:
        """Second array needed for beam chamber curve."""

    @y2_line.setter
    def y2_line(self, arg: Sequence[float], /) -> None: ...

    @property
    def ix_line(self) -> IntAlloc1D:
        """Used by wave and aperture curves."""

    @ix_line.setter
    def ix_line(self, arg: Sequence[int], /) -> None: ...

    @property
    def x_symb(self) -> RealAlloc1D:
        """Coords for drawing the symbols"""

    @x_symb.setter
    def x_symb(self, arg: Sequence[float], /) -> None: ...

    @property
    def y_symb(self) -> RealAlloc1D: ...

    @y_symb.setter
    def y_symb(self, arg: Sequence[float], /) -> None: ...

    @property
    def z_symb(self) -> RealAlloc1D:
        """Symbol color"""

    @z_symb.setter
    def z_symb(self, arg: Sequence[float], /) -> None: ...

    @property
    def err_symb(self) -> RealAlloc1D:
        """Error bars"""

    @err_symb.setter
    def err_symb(self, arg: Sequence[float], /) -> None: ...

    @property
    def symb_size(self) -> RealAlloc1D:
        """Symbol size. Used with symbol_size_scale."""

    @symb_size.setter
    def symb_size(self, arg: Sequence[float], /) -> None: ...

    @property
    def ix_symb(self) -> IntAlloc1D:
        """Corresponding index in d1_data%d(:) array."""

    @ix_symb.setter
    def ix_symb(self, arg: Sequence[int], /) -> None: ...

    @property
    def y_axis_scale_factor(self) -> float:
        """y-axis conversion from internal to plotting units."""

    @y_axis_scale_factor.setter
    def y_axis_scale_factor(self, arg: float, /) -> None: ...

    @property
    def line(self) -> QpLineStruct:
        """Line attributes"""

    @line.setter
    def line(self, arg: QpLineStruct, /) -> None: ...

    @property
    def symbol(self) -> QpSymbolStruct:
        """Symbol attributes"""

    @symbol.setter
    def symbol(self, arg: QpSymbolStruct, /) -> None: ...

    @property
    def orbit(self) -> TaoCurveOrbitStruct:
        """Used for E/B field plotting."""

    @orbit.setter
    def orbit(self, arg: TaoCurveOrbitStruct, /) -> None: ...

    @property
    def ix_universe(self) -> int:
        """Universe where data is. -1 => use s%global%default_universe"""

    @ix_universe.setter
    def ix_universe(self, arg: int, /) -> None: ...

    @property
    def symbol_every(self) -> int:
        """Symbol every how many points."""

    @symbol_every.setter
    def symbol_every(self, arg: int, /) -> None: ...

    @property
    def ix_branch(self) -> int: ...

    @ix_branch.setter
    def ix_branch(self, arg: int, /) -> None: ...

    @property
    def ix_bunch(self) -> int:
        """Bunch to plot."""

    @ix_bunch.setter
    def ix_bunch(self, arg: int, /) -> None: ...

    @property
    def n_turn(self) -> int:
        """Used for multi_turn_orbit plotting"""

    @n_turn.setter
    def n_turn(self, arg: int, /) -> None: ...

    @property
    def use_y2(self) -> bool:
        """Use y2 axis?"""

    @use_y2.setter
    def use_y2(self, arg: bool, /) -> None: ...

    @property
    def draw_line(self) -> bool:
        """Draw a line through the data points?"""

    @draw_line.setter
    def draw_line(self, arg: bool, /) -> None: ...

    @property
    def draw_symbols(self) -> bool:
        """Draw a symbol at the data points?"""

    @draw_symbols.setter
    def draw_symbols(self, arg: bool, /) -> None: ...

    @property
    def draw_symbol_index(self) -> bool:
        """Draw the symbol index number curve%ix_symb?"""

    @draw_symbol_index.setter
    def draw_symbol_index(self, arg: bool, /) -> None: ...

    @property
    def draw_error_bars(self) -> bool:
        """
        Draw error bars based upon data%error_rms if drawing data? !! logical :: draw_rms = .false.          ! Show mean and RMS values with legend?
        """

    @draw_error_bars.setter
    def draw_error_bars(self, arg: bool, /) -> None: ...

    @property
    def smooth_line_calc(self) -> bool:
        """Calculate data between element edge points?"""

    @smooth_line_calc.setter
    def smooth_line_calc(self, arg: bool, /) -> None: ...

    @property
    def valid(self) -> bool:
        """valid data?"""

    @valid.setter
    def valid(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoCurveStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoCurveStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoCurveStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoCurveStruct: ...

    def __eq__(self, arg: TaoCurveStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoD1DataArrayStruct:
    """Fortran struct: tao_d1_data_array_struct"""

    def __init__(self, d1: TaoD1DataStruct | None = None) -> None: ...

    @property
    def d1(self) -> TaoD1DataStruct | None: ...

    @d1.setter
    def d1(self, arg: TaoD1DataStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoD1DataArrayStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoD1DataArrayStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoD1DataArrayStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoD1DataArrayStruct: ...

    def __eq__(self, arg: TaoD1DataArrayStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoD1DataStruct:
    """Fortran struct: tao_d1_data_struct"""

    def __init__(self, name: str | None = None, d2: TaoD2DataStruct | None = None) -> None: ...

    @property
    def name(self) -> str:
        """Eg: 'x', etc."""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def d2(self) -> TaoD2DataStruct | None:
        """ptr to parent d2_data"""

    @d2.setter
    def d2(self, arg: TaoD2DataStruct, /) -> None: ...

    @property
    def d(self) -> TaoDataStructArray1D:
        """Pointer to the appropriate section in u%data"""

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoD1DataStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoD1DataStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoD1DataStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoD1DataStruct: ...

    def __eq__(self, arg: TaoD1DataStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoD2DataArrayStruct:
    """Fortran struct: tao_d2_data_array_struct"""

    def __init__(self, d2: TaoD2DataStruct | None = None) -> None: ...

    @property
    def d2(self) -> TaoD2DataStruct | None: ...

    @d2.setter
    def d2(self, arg: TaoD2DataStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoD2DataArrayStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoD2DataArrayStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoD2DataArrayStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoD2DataArrayStruct: ...

    def __eq__(self, arg: TaoD2DataArrayStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoD2DataStruct:
    """Fortran struct: tao_d2_data_struct"""

    def __init__(self, name: str | None = None, data_file_name: str | None = None, ref_file_name: str | None = None, data_date: str | None = None, ref_date: str | None = None, ix_universe: int | None = None, ix_d2_data: int | None = None, ix_ref: int | None = None, data_read_in: bool | None = None, ref_read_in: bool | None = None) -> None: ...

    @property
    def name(self) -> str:
        """Name to be used with commands."""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def data_file_name(self) -> str:
        """Data file name ."""

    @data_file_name.setter
    def data_file_name(self, arg: str, /) -> None: ...

    @property
    def ref_file_name(self) -> str:
        """Reference file name."""

    @ref_file_name.setter
    def ref_file_name(self, arg: str, /) -> None: ...

    @property
    def data_date(self) -> str:
        """Data measurement date."""

    @data_date.setter
    def data_date(self, arg: str, /) -> None: ...

    @property
    def ref_date(self) -> str:
        """Reference data measurement date."""

    @ref_date.setter
    def ref_date(self, arg: str, /) -> None: ...

    @property
    def descrip(self) -> FCharArray1D:
        """Array for descriptive information."""

    @property
    def d1(self) -> TaoD1DataStructAlloc1D:
        """Points to children"""

    @property
    def ix_universe(self) -> int:
        """Index of universe this is in."""

    @ix_universe.setter
    def ix_universe(self, arg: int, /) -> None: ...

    @property
    def ix_d2_data(self) -> int:
        """Index in u%d2_data(:) array."""

    @ix_d2_data.setter
    def ix_d2_data(self, arg: int, /) -> None: ...

    @property
    def ix_ref(self) -> int:
        """Index of the reference data set."""

    @ix_ref.setter
    def ix_ref(self, arg: int, /) -> None: ...

    @property
    def data_read_in(self) -> bool:
        """A data set has been read in?"""

    @data_read_in.setter
    def data_read_in(self, arg: bool, /) -> None: ...

    @property
    def ref_read_in(self) -> bool:
        """A reference data set has been read in?"""

    @ref_read_in.setter
    def ref_read_in(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoD2DataStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoD2DataStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoD2DataStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoD2DataStruct: ...

    def __eq__(self, arg: TaoD2DataStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoDataArrayStruct:
    """Fortran struct: tao_data_array_struct"""

    def __init__(self, d: TaoDataStruct | None = None) -> None: ...

    @property
    def d(self) -> TaoDataStruct | None: ...

    @d.setter
    def d(self, arg: TaoDataStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoDataArrayStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoDataArrayStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoDataArrayStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoDataArrayStruct: ...

    def __eq__(self, arg: TaoDataArrayStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoDataStruct:
    """Fortran struct: tao_data_struct"""

    def __init__(self, ele_name: str | None = None, ele_start_name: str | None = None, ele_ref_name: str | None = None, data_type: str | None = None, merit_type: str | None = None, id: str | None = None, data_source: str | None = None, why_invalid: str | None = None, ix_uni: int | None = None, ix_bunch: int | None = None, ix_branch: int | None = None, ix_ele: int | None = None, ix_ele_start: int | None = None, ix_ele_ref: int | None = None, ix_ele_merit: int | None = None, ix_d1: int | None = None, ix_data: int | None = None, ix_dModel: int | None = None, eval_point: int | None = None, meas_value: float | None = None, ref_value: float | None = None, model_value: float | None = None, design_value: float | None = None, old_value: float | None = None, base_value: float | None = None, error_rms: float | None = None, delta_merit: float | None = None, weight: float | None = None, invalid_value: float | None = None, merit: float | None = None, s: float | None = None, s_offset: float | None = None, ref_s_offset: float | None = None, err_message_printed: bool | None = None, exists: bool | None = None, good_model: bool | None = None, good_base: bool | None = None, good_design: bool | None = None, good_meas: bool | None = None, good_ref: bool | None = None, good_user: bool | None = None, good_opt: bool | None = None, good_plot: bool | None = None, useit_plot: bool | None = None, useit_opt: bool | None = None, spin_map: TaoSpinMapStruct | None = None, d1: TaoD1DataStruct | None = None) -> None: ...

    @property
    def ele_name(self) -> str:
        """Name of the lattice element where datum is evaluated."""

    @ele_name.setter
    def ele_name(self, arg: str, /) -> None: ...

    @property
    def ele_start_name(self) -> str:
        """Name of starting lattice element when there is a range"""

    @ele_start_name.setter
    def ele_start_name(self, arg: str, /) -> None: ...

    @property
    def ele_ref_name(self) -> str:
        """Name of reference lattice element"""

    @ele_ref_name.setter
    def ele_ref_name(self, arg: str, /) -> None: ...

    @property
    def data_type(self) -> str:
        """Type of data: 'orbit.x', etc."""

    @data_type.setter
    def data_type(self, arg: str, /) -> None: ...

    @property
    def merit_type(self) -> str:
        """Type of constraint: 'target', 'max', 'min', etc."""

    @merit_type.setter
    def merit_type(self, arg: str, /) -> None: ...

    @property
    def id(self) -> str:
        """Used by Tao extension code. Not used by Tao directly."""

    @id.setter
    def id(self, arg: str, /) -> None: ...

    @property
    def data_source(self) -> str:
        """'lat', 'beam', 'data' or 'var'. Last two used for expressions."""

    @data_source.setter
    def data_source(self, arg: str, /) -> None: ...

    @property
    def why_invalid(self) -> str:
        """Informational string if there is a problem."""

    @why_invalid.setter
    def why_invalid(self, arg: str, /) -> None: ...

    @property
    def ix_uni(self) -> int:
        """Universe index of datum."""

    @ix_uni.setter
    def ix_uni(self, arg: int, /) -> None: ...

    @property
    def ix_bunch(self) -> int:
        """Bunch number to get the data from."""

    @ix_bunch.setter
    def ix_bunch(self, arg: int, /) -> None: ...

    @property
    def ix_branch(self) -> int:
        """Index of the associated lattice branch."""

    @ix_branch.setter
    def ix_branch(self, arg: int, /) -> None: ...

    @property
    def ix_ele(self) -> int:
        """Index of the lattice element corresponding to ele_name"""

    @ix_ele.setter
    def ix_ele(self, arg: int, /) -> None: ...

    @property
    def ix_ele_start(self) -> int:
        """Index of lattice elment when there is a range"""

    @ix_ele_start.setter
    def ix_ele_start(self, arg: int, /) -> None: ...

    @property
    def ix_ele_ref(self) -> int:
        """Index of lattice elment when there is a reference."""

    @ix_ele_ref.setter
    def ix_ele_ref(self, arg: int, /) -> None: ...

    @property
    def ix_ele_merit(self) -> int:
        """Index of lattice elment where merit is evaluated."""

    @ix_ele_merit.setter
    def ix_ele_merit(self, arg: int, /) -> None: ...

    @property
    def ix_d1(self) -> int:
        """Index number in u%d2_data(i)%d1_data(j)%d(:) array."""

    @ix_d1.setter
    def ix_d1(self, arg: int, /) -> None: ...

    @property
    def ix_data(self) -> int:
        """Index of this datum in the u%data(:) array of data_structs."""

    @ix_data.setter
    def ix_data(self, arg: int, /) -> None: ...

    @property
    def ix_dModel(self) -> int:
        """Row number in the dModel_dVar derivative matrix."""

    @ix_dModel.setter
    def ix_dModel(self, arg: int, /) -> None: ...

    @property
    def eval_point(self) -> int:
        """
        or anchor_center$, anchor_beginning$. Where to evaluate data relative to the element.
        """

    @eval_point.setter
    def eval_point(self, arg: int, /) -> None: ...

    @property
    def meas_value(self) -> float:
        """Measured datum value."""

    @meas_value.setter
    def meas_value(self, arg: float, /) -> None: ...

    @property
    def ref_value(self) -> float:
        """Measured datum value from the reference data set."""

    @ref_value.setter
    def ref_value(self, arg: float, /) -> None: ...

    @property
    def model_value(self) -> float:
        """Datum value as calculated from the model."""

    @model_value.setter
    def model_value(self, arg: float, /) -> None: ...

    @property
    def design_value(self) -> float:
        """What the datum value is in the design lattice."""

    @design_value.setter
    def design_value(self, arg: float, /) -> None: ...

    @property
    def old_value(self) -> float:
        """The model_value at some previous time."""

    @old_value.setter
    def old_value(self, arg: float, /) -> None: ...

    @property
    def base_value(self) -> float:
        """The value as calculated from the base model."""

    @base_value.setter
    def base_value(self, arg: float, /) -> None: ...

    @property
    def error_rms(self) -> float:
        """Measurement error RMS. Used in plotting."""

    @error_rms.setter
    def error_rms(self, arg: float, /) -> None: ...

    @property
    def delta_merit(self) -> float:
        """Diff used to calculate the merit function term."""

    @delta_merit.setter
    def delta_merit(self, arg: float, /) -> None: ...

    @property
    def weight(self) -> float:
        """Weight for the merit function term."""

    @weight.setter
    def weight(self, arg: float, /) -> None: ...

    @property
    def invalid_value(self) -> float:
        """
        Value used in merit calc if good_model = F (or possibly good_design & good_base).
        """

    @invalid_value.setter
    def invalid_value(self, arg: float, /) -> None: ...

    @property
    def merit(self) -> float:
        """Merit function term value: weight * delta_merit^2"""

    @merit.setter
    def merit(self, arg: float, /) -> None: ...

    @property
    def s(self) -> float:
        """longitudinal position of ele."""

    @s.setter
    def s(self, arg: float, /) -> None: ...

    @property
    def s_offset(self) -> float:
        """Offset of the evaluation point."""

    @s_offset.setter
    def s_offset(self, arg: float, /) -> None: ...

    @property
    def ref_s_offset(self) -> float:
        """Offset of the reference point. In development."""

    @ref_s_offset.setter
    def ref_s_offset(self, arg: float, /) -> None: ...

    @property
    def err_message_printed(self) -> bool:
        """Used to prevent zillions of error messages being generated"""

    @err_message_printed.setter
    def err_message_printed(self, arg: bool, /) -> None: ...

    @property
    def exists(self) -> bool:
        """See above"""

    @exists.setter
    def exists(self, arg: bool, /) -> None: ...

    @property
    def good_model(self) -> bool:
        """See above"""

    @good_model.setter
    def good_model(self, arg: bool, /) -> None: ...

    @property
    def good_base(self) -> bool:
        """See above"""

    @good_base.setter
    def good_base(self, arg: bool, /) -> None: ...

    @property
    def good_design(self) -> bool:
        """See above"""

    @good_design.setter
    def good_design(self, arg: bool, /) -> None: ...

    @property
    def good_meas(self) -> bool:
        """See above"""

    @good_meas.setter
    def good_meas(self, arg: bool, /) -> None: ...

    @property
    def good_ref(self) -> bool:
        """See above"""

    @good_ref.setter
    def good_ref(self, arg: bool, /) -> None: ...

    @property
    def good_user(self) -> bool:
        """See above"""

    @good_user.setter
    def good_user(self, arg: bool, /) -> None: ...

    @property
    def good_opt(self) -> bool:
        """See above"""

    @good_opt.setter
    def good_opt(self, arg: bool, /) -> None: ...

    @property
    def good_plot(self) -> bool:
        """See above"""

    @good_plot.setter
    def good_plot(self, arg: bool, /) -> None: ...

    @property
    def useit_plot(self) -> bool:
        """See above"""

    @useit_plot.setter
    def useit_plot(self, arg: bool, /) -> None: ...

    @property
    def useit_opt(self) -> bool:
        """See above"""

    @useit_opt.setter
    def useit_opt(self, arg: bool, /) -> None: ...

    @property
    def spin_map(self) -> TaoSpinMapStruct: ...

    @spin_map.setter
    def spin_map(self, arg: TaoSpinMapStruct, /) -> None: ...

    @property
    def d1(self) -> TaoD1DataStruct | None:
        """Pointer to the parent d1_data_struct"""

    @d1.setter
    def d1(self, arg: TaoD1DataStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoDataStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoDataStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoDataStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoDataStruct: ...

    def __eq__(self, arg: TaoDataStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoDataVarComponentStruct:
    """Fortran struct: tao_data_var_component_struct"""

    def __init__(self, name: str | None = None, sign: float | None = None) -> None: ...

    @property
    def name(self) -> str:
        """Eg: 'meas', 'ref', 'model', etc."""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def sign(self) -> float:
        """+1 or -1"""

    @sign.setter
    def sign(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoDataVarComponentStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoDataVarComponentStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoDataVarComponentStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoDataVarComponentStruct: ...

    def __eq__(self, arg: TaoDataVarComponentStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoDrawingStruct:
    """Fortran struct: tao_drawing_struct"""

    def __init__(self) -> None: ...

    @property
    def ele_shape(self) -> TaoEleShapeStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoDrawingStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoDrawingStruct: ...

    def __eq__(self, arg: TaoDrawingStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoDynamicApertureStruct:
    """Fortran struct: tao_dynamic_aperture_struct"""

    def __init__(self, param: ApertureParamStruct | None = None, pz: Sequence[float] | None = None, ellipse_scale: float | None = None, a_emit: float | None = None, b_emit: float | None = None) -> None: ...

    @property
    def param(self) -> ApertureParamStruct: ...

    @param.setter
    def param(self, arg: ApertureParamStruct, /) -> None: ...

    @property
    def scan(self) -> ApertureScanStructAlloc1D:
        """One scan for each pz."""

    @property
    def pz(self) -> RealAlloc1D: ...

    @pz.setter
    def pz(self, arg: Sequence[float], /) -> None: ...

    @property
    def ellipse_scale(self) -> float: ...

    @ellipse_scale.setter
    def ellipse_scale(self, arg: float, /) -> None: ...

    @property
    def a_emit(self) -> float: ...

    @a_emit.setter
    def a_emit(self, arg: float, /) -> None: ...

    @property
    def b_emit(self) -> float: ...

    @b_emit.setter
    def b_emit(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoDynamicApertureStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoDynamicApertureStruct: ...

    def __eq__(self, arg: TaoDynamicApertureStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoElePointerStruct:
    """Fortran struct: tao_ele_pointer_struct"""

    def __init__(self, n_loc: int | None = None) -> None: ...

    @property
    def eles(self) -> ElePointerStructAlloc1D: ...

    @property
    def n_loc(self) -> int: ...

    @n_loc.setter
    def n_loc(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoElePointerStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoElePointerStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoElePointerStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoElePointerStruct: ...

    def __eq__(self, arg: TaoElePointerStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoEleShapeInput:
    """Fortran struct: tao_ele_shape_input"""

    def __init__(self, ele_id: str | None = None, shape: str | None = None, color: str | None = None, size: float | None = None, label: str | None = None, draw: bool | None = None, multi: bool | None = None, line_width: int | None = None, offset: float | None = None) -> None: ...

    @property
    def ele_id(self) -> str:
        """element 'key::name' to match to."""

    @ele_id.setter
    def ele_id(self, arg: str, /) -> None: ...

    @property
    def shape(self) -> str:
        """Shape to draw"""

    @shape.setter
    def shape(self, arg: str, /) -> None: ...

    @property
    def color(self) -> str:
        """Color of shape"""

    @color.setter
    def color(self, arg: str, /) -> None: ...

    @property
    def size(self) -> float:
        """plot vertical height"""

    @size.setter
    def size(self, arg: float, /) -> None: ...

    @property
    def label(self) -> str:
        """Can be: 'name', 's', 'none'"""

    @label.setter
    def label(self, arg: str, /) -> None: ...

    @property
    def draw(self) -> bool:
        """Draw the shape?"""

    @draw.setter
    def draw(self, arg: bool, /) -> None: ...

    @property
    def multi(self) -> bool:
        """Can be part of a multi-shape."""

    @multi.setter
    def multi(self, arg: bool, /) -> None: ...

    @property
    def line_width(self) -> int:
        """Width of lines used to draw the shape."""

    @line_width.setter
    def line_width(self, arg: int, /) -> None: ...

    @property
    def offset(self) -> float:
        """Vertical offset."""

    @offset.setter
    def offset(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoEleShapeInput: ...

    def __deepcopy__(self, arg: dict, /) -> TaoEleShapeInput: ...

    def __eq__(self, arg: TaoEleShapeInput, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoEleShapeStruct:
    """Fortran struct: tao_ele_shape_struct"""

    def __init__(self, ele_id: str | None = None, shape: str | None = None, color: str | None = None, size: float | None = None, label: str | None = None, draw: bool | None = None, multi: bool | None = None, line_width: int | None = None, offset: float | None = None, ix_key: int | None = None, name_ele: str | None = None) -> None: ...

    @property
    def ele_id(self) -> str:
        """element 'key::name' to match to."""

    @ele_id.setter
    def ele_id(self, arg: str, /) -> None: ...

    @property
    def shape(self) -> str:
        """Shape to draw"""

    @shape.setter
    def shape(self, arg: str, /) -> None: ...

    @property
    def color(self) -> str:
        """Color of shape"""

    @color.setter
    def color(self, arg: str, /) -> None: ...

    @property
    def size(self) -> float:
        """plot vertical height"""

    @size.setter
    def size(self, arg: float, /) -> None: ...

    @property
    def label(self) -> str:
        """Can be: 'name', 's', 'none'"""

    @label.setter
    def label(self, arg: str, /) -> None: ...

    @property
    def draw(self) -> bool:
        """Draw the shape?"""

    @draw.setter
    def draw(self, arg: bool, /) -> None: ...

    @property
    def multi(self) -> bool:
        """Can be part of a multi-shape."""

    @multi.setter
    def multi(self, arg: bool, /) -> None: ...

    @property
    def line_width(self) -> int:
        """Width of lines used to draw the shape."""

    @line_width.setter
    def line_width(self, arg: int, /) -> None: ...

    @property
    def offset(self) -> float:
        """Vertical offset."""

    @offset.setter
    def offset(self, arg: float, /) -> None: ...

    @property
    def ix_key(self) -> int:
        """Extracted from ele_id. 0 => all classes (quadrupole, etc.)"""

    @ix_key.setter
    def ix_key(self, arg: int, /) -> None: ...

    @property
    def name_ele(self) -> str:
        """Name of element."""

    @name_ele.setter
    def name_ele(self, arg: str, /) -> None: ...

    @property
    def uni(self) -> TaoElePointerStructAlloc1D: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoEleShapeStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoEleShapeStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoEleShapeStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoEleShapeStruct: ...

    def __eq__(self, arg: TaoEleShapeStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoEvalNodeStruct:
    """Fortran struct: tao_eval_node_struct"""

    def __init__(self, type: int | None = None, name: str | None = None, scale: float | None = None, value: Sequence[float] | None = None) -> None: ...

    @property
    def type(self) -> int: ...

    @type.setter
    def type(self, arg: int, /) -> None: ...

    @property
    def name(self) -> str: ...

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def scale(self) -> float:
        """Scale factor for ping data"""

    @scale.setter
    def scale(self, arg: float, /) -> None: ...

    @property
    def value(self) -> RealAlloc1D: ...

    @value.setter
    def value(self, arg: Sequence[float], /) -> None: ...

    @property
    def info(self) -> TaoExpressionInfoStructAlloc1D: ...

    @property
    def value_ptr(self) -> TaoRealPointerStructAlloc1D:
        """Used to point to data, lattice parameters, etc"""

    @property
    def node(self) -> TaoEvalNodeStructArray1D:
        """Child nodes for tree construction."""

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoEvalNodeStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoEvalNodeStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoEvalNodeStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoEvalNodeStruct: ...

    def __eq__(self, arg: TaoEvalNodeStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoExpressionInfoStruct:
    """Fortran struct: tao_expression_info_struct"""

    def __init__(self, good: bool | None = None, ele: EleStruct | None = None, s: float | None = None) -> None: ...

    @property
    def good(self) -> bool:
        """Expression is valid."""

    @good.setter
    def good(self, arg: bool, /) -> None: ...

    @property
    def ele(self) -> EleStruct | None:
        """Associated ele if it exists"""

    @ele.setter
    def ele(self, arg: EleStruct, /) -> None: ...

    @property
    def s(self) -> float:
        """Longitudinal position of expression."""

    @s.setter
    def s(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoExpressionInfoStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoExpressionInfoStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoExpressionInfoStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoExpressionInfoStruct: ...

    def __eq__(self, arg: TaoExpressionInfoStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoFloorPlanStruct:
    """Fortran struct: tao_floor_plan_struct"""

    def __init__(self, view: str | None = None, rotation: float | None = None, correct_distortion: bool | None = None, flip_label_side: bool | None = None, size_is_absolute: bool | None = None, draw_only_first_pass: bool | None = None, draw_building_wall: bool | None = None, orbit_scale: float | None = None, orbit_color: str | None = None, orbit_pattern: str | None = None, orbit_lattice: str | None = None, orbit_width: int | None = None) -> None: ...

    @property
    def view(self) -> str:
        """or 'xz'."""

    @view.setter
    def view(self, arg: str, /) -> None: ...

    @property
    def rotation(self) -> float:
        """Rotation of floor plan plot: 1.0 -> 360^deg"""

    @rotation.setter
    def rotation(self, arg: float, /) -> None: ...

    @property
    def correct_distortion(self) -> bool:
        """T -> Shrink one axis so x-scale = y-scale."""

    @correct_distortion.setter
    def correct_distortion(self, arg: bool, /) -> None: ...

    @property
    def flip_label_side(self) -> bool:
        """Draw element label on other side of element?"""

    @flip_label_side.setter
    def flip_label_side(self, arg: bool, /) -> None: ...

    @property
    def size_is_absolute(self) -> bool:
        """Are shape sizes in meters or window pixels?"""

    @size_is_absolute.setter
    def size_is_absolute(self, arg: bool, /) -> None: ...

    @property
    def draw_only_first_pass(self) -> bool:
        """Draw only first pass with multipass elements?"""

    @draw_only_first_pass.setter
    def draw_only_first_pass(self, arg: bool, /) -> None: ...

    @property
    def draw_building_wall(self) -> bool:
        """Draw the building wall?"""

    @draw_building_wall.setter
    def draw_building_wall(self, arg: bool, /) -> None: ...

    @property
    def orbit_scale(self) -> float:
        """Scale factor for drawing orbits. 0 -> Do not draw."""

    @orbit_scale.setter
    def orbit_scale(self, arg: float, /) -> None: ...

    @property
    def orbit_color(self) -> str: ...

    @orbit_color.setter
    def orbit_color(self, arg: str, /) -> None: ...

    @property
    def orbit_pattern(self) -> str: ...

    @orbit_pattern.setter
    def orbit_pattern(self, arg: str, /) -> None: ...

    @property
    def orbit_lattice(self) -> str:
        """Or 'design' or 'base'"""

    @orbit_lattice.setter
    def orbit_lattice(self, arg: str, /) -> None: ...

    @property
    def orbit_width(self) -> int: ...

    @orbit_width.setter
    def orbit_width(self, arg: int, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoFloorPlanStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoFloorPlanStruct: ...

    def __eq__(self, arg: TaoFloorPlanStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoGlobalStruct:
    """Fortran struct: tao_global_struct"""

    def __init__(self, beam_dead_cutoff: float | None = None, lm_opt_deriv_reinit: float | None = None, de_lm_step_ratio: float | None = None, de_var_to_population_factor: float | None = None, lmdif_eps: float | None = None, lmdif_negligible_merit: float | None = None, svd_cutoff: float | None = None, unstable_penalty: float | None = None, merit_stop_value: float | None = None, dmerit_stop_value: float | None = None, random_sigma_cutoff: float | None = None, delta_e_chrom: float | None = None, max_plot_time: float | None = None, default_universe: int | None = None, default_branch: int | None = None, n_opti_cycles: int | None = None, n_opti_loops: int | None = None, n_threads: int | None = None, phase_units: int | None = None, bunch_to_plot: int | None = None, random_seed: int | None = None, n_top10_merit: int | None = None, srdt_gen_n_slices: int | None = None, datum_err_messages_max: int | None = None, srdt_sxt_n_slices: int | None = None, srdt_use_cache: bool | None = None, quiet: str | None = None, random_engine: str | None = None, random_gauss_converter: str | None = None, track_type: str | None = None, lat_sigma_calc_uses_emit_from: str | None = None, prompt_string: str | None = None, prompt_color: str | None = None, optimizer: str | None = None, print_command: str | None = None, var_out_file: str | None = None, history_file: str | None = None, beam_timer_on: bool | None = None, box_plots: bool | None = None, blank_line_between_commands: bool | None = None, cmd_file_abort_on_error: bool | None = None, concatenate_maps: bool | None = None, derivative_recalc: bool | None = None, derivative_uses_design: bool | None = None, disable_smooth_line_calc: bool | None = None, draw_curve_off_scale_warn: bool | None = None, external_plotting: bool | None = None, label_lattice_elements: bool | None = None, label_keys: bool | None = None, lattice_calc_on: bool | None = None, only_limit_opt_vars: bool | None = None, opt_with_ref: bool | None = None, opt_with_base: bool | None = None, opt_match_auto_recalc: bool | None = None, opti_write_var_file: bool | None = None, optimizer_allow_user_abort: bool | None = None, optimizer_var_limit_warn: bool | None = None, plot_on: bool | None = None, rad_int_user_calc_on: bool | None = None, rf_on: bool | None = None, single_step: bool | None = None, stop_on_error: bool | None = None, svd_retreat_on_merit_increase: bool | None = None, var_limits_on: bool | None = None, wait_for_CR_in_single_mode: bool | None = None, symbol_import: bool | None = None, debug_on: bool | None = None, expression_tree_on: bool | None = None, verbose_on: bool | None = None) -> None: ...

    @property
    def beam_dead_cutoff(self) -> float:
        """Percentage of dead particles at which beam tracking is stopped."""

    @beam_dead_cutoff.setter
    def beam_dead_cutoff(self, arg: float, /) -> None: ...

    @property
    def lm_opt_deriv_reinit(self) -> float:
        """Reinit derivative matrix cutoff"""

    @lm_opt_deriv_reinit.setter
    def lm_opt_deriv_reinit(self, arg: float, /) -> None: ...

    @property
    def de_lm_step_ratio(self) -> float:
        """Scaling for step sizes between DE and LM optimizers."""

    @de_lm_step_ratio.setter
    def de_lm_step_ratio(self, arg: float, /) -> None: ...

    @property
    def de_var_to_population_factor(self) -> float:
        """DE population = max(n_var*factor, 20)"""

    @de_var_to_population_factor.setter
    def de_var_to_population_factor(self, arg: float, /) -> None: ...

    @property
    def lmdif_eps(self) -> float:
        """Tollerance for lmdif optimizer."""

    @lmdif_eps.setter
    def lmdif_eps(self, arg: float, /) -> None: ...

    @property
    def lmdif_negligible_merit(self) -> float: ...

    @lmdif_negligible_merit.setter
    def lmdif_negligible_merit(self, arg: float, /) -> None: ...

    @property
    def svd_cutoff(self) -> float:
        """SVD singular value cutoff."""

    @svd_cutoff.setter
    def svd_cutoff(self, arg: float, /) -> None: ...

    @property
    def unstable_penalty(self) -> float:
        """Used in unstable_ring datum merit calculation."""

    @unstable_penalty.setter
    def unstable_penalty(self, arg: float, /) -> None: ...

    @property
    def merit_stop_value(self) -> float:
        """Merit value below which an optimizer will stop."""

    @merit_stop_value.setter
    def merit_stop_value(self, arg: float, /) -> None: ...

    @property
    def dmerit_stop_value(self) -> float:
        """Fractional Merit change below which an optimizer will stop."""

    @dmerit_stop_value.setter
    def dmerit_stop_value(self, arg: float, /) -> None: ...

    @property
    def random_sigma_cutoff(self) -> float:
        """Cut-off in sigmas."""

    @random_sigma_cutoff.setter
    def random_sigma_cutoff(self, arg: float, /) -> None: ...

    @property
    def delta_e_chrom(self) -> float:
        """Delta E used from chrom calc."""

    @delta_e_chrom.setter
    def delta_e_chrom(self, arg: float, /) -> None: ...

    @property
    def max_plot_time(self) -> float:
        """If plotting time (seconds) exceeds this than a message is generated."""

    @max_plot_time.setter
    def max_plot_time(self, arg: float, /) -> None: ...

    @property
    def default_universe(self) -> int:
        """Default universe to work with."""

    @default_universe.setter
    def default_universe(self, arg: int, /) -> None: ...

    @property
    def default_branch(self) -> int:
        """Default lattice branch to work with."""

    @default_branch.setter
    def default_branch(self, arg: int, /) -> None: ...

    @property
    def n_opti_cycles(self) -> int:
        """Number of optimization cycles"""

    @n_opti_cycles.setter
    def n_opti_cycles(self, arg: int, /) -> None: ...

    @property
    def n_opti_loops(self) -> int:
        """Number of optimization loops"""

    @n_opti_loops.setter
    def n_opti_loops(self, arg: int, /) -> None: ...

    @property
    def n_threads(self) -> int:
        """Number of OpenMP threads for parallel calculations."""

    @n_threads.setter
    def n_threads(self, arg: int, /) -> None: ...

    @property
    def phase_units(self) -> int:
        """Phase units on output."""

    @phase_units.setter
    def phase_units(self, arg: int, /) -> None: ...

    @property
    def bunch_to_plot(self) -> int:
        """Which bunch to plot"""

    @bunch_to_plot.setter
    def bunch_to_plot(self, arg: int, /) -> None: ...

    @property
    def random_seed(self) -> int:
        """Use system clock by default"""

    @random_seed.setter
    def random_seed(self, arg: int, /) -> None: ...

    @property
    def n_top10_merit(self) -> int:
        """Number of top merit constraints to print."""

    @n_top10_merit.setter
    def n_top10_merit(self, arg: int, /) -> None: ...

    @property
    def srdt_gen_n_slices(self) -> int:
        """Number times to slice elements for summation RDT calculation"""

    @srdt_gen_n_slices.setter
    def srdt_gen_n_slices(self, arg: int, /) -> None: ...

    @property
    def datum_err_messages_max(self) -> int:
        """Maximum number of error messages per call to lattice_calc."""

    @datum_err_messages_max.setter
    def datum_err_messages_max(self, arg: int, /) -> None: ...

    @property
    def srdt_sxt_n_slices(self) -> int:
        """Number times to slice sextupoles for summation RDT calculation"""

    @srdt_sxt_n_slices.setter
    def srdt_sxt_n_slices(self, arg: int, /) -> None: ...

    @property
    def srdt_use_cache(self) -> bool:
        """
        Create cache for SRDT calculations.  Can use lots of memory if srdt_*_n_slices large.
        """

    @srdt_use_cache.setter
    def srdt_use_cache(self, arg: bool, /) -> None: ...

    @property
    def quiet(self) -> str:
        """Print I/O when running a command file?"""

    @quiet.setter
    def quiet(self, arg: str, /) -> None: ...

    @property
    def random_engine(self) -> str:
        """Non-beam random number engine"""

    @random_engine.setter
    def random_engine(self, arg: str, /) -> None: ...

    @property
    def random_gauss_converter(self) -> str:
        """Non-beam"""

    @random_gauss_converter.setter
    def random_gauss_converter(self, arg: str, /) -> None: ...

    @property
    def track_type(self) -> str:
        """or 'beam'"""

    @track_type.setter
    def track_type(self, arg: str, /) -> None: ...

    @property
    def lat_sigma_calc_uses_emit_from(self) -> str:
        """
        Lattice derived sigma matrix uses emit values from where? Other possibilities: 'beam', 'beam_init'.
        """

    @lat_sigma_calc_uses_emit_from.setter
    def lat_sigma_calc_uses_emit_from(self, arg: str, /) -> None: ...

    @property
    def prompt_string(self) -> str: ...

    @prompt_string.setter
    def prompt_string(self, arg: str, /) -> None: ...

    @property
    def prompt_color(self) -> str:
        """See read_a_line routine for possible settings."""

    @prompt_color.setter
    def prompt_color(self, arg: str, /) -> None: ...

    @property
    def optimizer(self) -> str:
        """optimizer to use."""

    @optimizer.setter
    def optimizer(self, arg: str, /) -> None: ...

    @property
    def print_command(self) -> str: ...

    @print_command.setter
    def print_command(self, arg: str, /) -> None: ...

    @property
    def var_out_file(self) -> str: ...

    @var_out_file.setter
    def var_out_file(self, arg: str, /) -> None: ...

    @property
    def history_file(self) -> str: ...

    @history_file.setter
    def history_file(self, arg: str, /) -> None: ...

    @property
    def beam_timer_on(self) -> bool:
        """For timing the beam tracking calculation."""

    @beam_timer_on.setter
    def beam_timer_on(self, arg: bool, /) -> None: ...

    @property
    def box_plots(self) -> bool:
        """For debugging plot layout issues."""

    @box_plots.setter
    def box_plots(self, arg: bool, /) -> None: ...

    @property
    def blank_line_between_commands(self) -> bool:
        """Add a blank line between command output?"""

    @blank_line_between_commands.setter
    def blank_line_between_commands(self, arg: bool, /) -> None: ...

    @property
    def cmd_file_abort_on_error(self) -> bool:
        """Abort open command files if there is an error?"""

    @cmd_file_abort_on_error.setter
    def cmd_file_abort_on_error(self, arg: bool, /) -> None: ...

    @property
    def concatenate_maps(self) -> bool:
        """False => tracking using DA."""

    @concatenate_maps.setter
    def concatenate_maps(self, arg: bool, /) -> None: ...

    @property
    def derivative_recalc(self) -> bool:
        """Recalc before each optimizer run?"""

    @derivative_recalc.setter
    def derivative_recalc(self, arg: bool, /) -> None: ...

    @property
    def derivative_uses_design(self) -> bool:
        """Derivative calc uses design lattice instead of model?"""

    @derivative_uses_design.setter
    def derivative_uses_design(self, arg: bool, /) -> None: ...

    @property
    def disable_smooth_line_calc(self) -> bool:
        """Global disable of the smooth line calculation."""

    @disable_smooth_line_calc.setter
    def disable_smooth_line_calc(self, arg: bool, /) -> None: ...

    @property
    def draw_curve_off_scale_warn(self) -> bool:
        """Display warning on graphs?"""

    @draw_curve_off_scale_warn.setter
    def draw_curve_off_scale_warn(self, arg: bool, /) -> None: ...

    @property
    def external_plotting(self) -> bool:
        """Used with matplotlib and gui."""

    @external_plotting.setter
    def external_plotting(self, arg: bool, /) -> None: ...

    @property
    def label_lattice_elements(self) -> bool:
        """For lat_layout plots"""

    @label_lattice_elements.setter
    def label_lattice_elements(self, arg: bool, /) -> None: ...

    @property
    def label_keys(self) -> bool:
        """For lat_layout plots"""

    @label_keys.setter
    def label_keys(self, arg: bool, /) -> None: ...

    @property
    def lattice_calc_on(self) -> bool:
        """Turn on/off beam and single particle calculations."""

    @lattice_calc_on.setter
    def lattice_calc_on(self, arg: bool, /) -> None: ...

    @property
    def only_limit_opt_vars(self) -> bool:
        """Only apply limits to variables used in optimization."""

    @only_limit_opt_vars.setter
    def only_limit_opt_vars(self, arg: bool, /) -> None: ...

    @property
    def opt_with_ref(self) -> bool:
        """Use reference data in optimization?"""

    @opt_with_ref.setter
    def opt_with_ref(self, arg: bool, /) -> None: ...

    @property
    def opt_with_base(self) -> bool:
        """Use base data in optimization?"""

    @opt_with_base.setter
    def opt_with_base(self, arg: bool, /) -> None: ...

    @property
    def opt_match_auto_recalc(self) -> bool:
        """Set recalc = True for match elements before each cycle?"""

    @opt_match_auto_recalc.setter
    def opt_match_auto_recalc(self, arg: bool, /) -> None: ...

    @property
    def opti_write_var_file(self) -> bool:
        """'run' command writes var_out_file"""

    @opti_write_var_file.setter
    def opti_write_var_file(self, arg: bool, /) -> None: ...

    @property
    def optimizer_allow_user_abort(self) -> bool:
        """See Tao manual for more details."""

    @optimizer_allow_user_abort.setter
    def optimizer_allow_user_abort(self, arg: bool, /) -> None: ...

    @property
    def optimizer_var_limit_warn(self) -> bool:
        """Warn when vars reach a limit with optimization."""

    @optimizer_var_limit_warn.setter
    def optimizer_var_limit_warn(self, arg: bool, /) -> None: ...

    @property
    def plot_on(self) -> bool:
        """Do plotting?"""

    @plot_on.setter
    def plot_on(self, arg: bool, /) -> None: ...

    @property
    def rad_int_user_calc_on(self) -> bool:
        """User set radiation integrals calculation on/off."""

    @rad_int_user_calc_on.setter
    def rad_int_user_calc_on(self, arg: bool, /) -> None: ...

    @property
    def rf_on(self) -> bool:
        """RFcavities on or off? Does not affect lcavities."""

    @rf_on.setter
    def rf_on(self, arg: bool, /) -> None: ...

    @property
    def single_step(self) -> bool:
        """For debugging and demonstrations: Single step through a command file?"""

    @single_step.setter
    def single_step(self, arg: bool, /) -> None: ...

    @property
    def stop_on_error(self) -> bool:
        """For debugging: False prevents tao from exiting on an error."""

    @stop_on_error.setter
    def stop_on_error(self, arg: bool, /) -> None: ...

    @property
    def svd_retreat_on_merit_increase(self) -> bool: ...

    @svd_retreat_on_merit_increase.setter
    def svd_retreat_on_merit_increase(self, arg: bool, /) -> None: ...

    @property
    def var_limits_on(self) -> bool:
        """Respect the variable limits?"""

    @var_limits_on.setter
    def var_limits_on(self, arg: bool, /) -> None: ...

    @property
    def wait_for_CR_in_single_mode(self) -> bool:
        """For use with a python GUI."""

    @wait_for_CR_in_single_mode.setter
    def wait_for_CR_in_single_mode(self, arg: bool, /) -> None: ...

    @property
    def symbol_import(self) -> bool:
        """Import symbols from lattice file(s)? Internal stuff"""

    @symbol_import.setter
    def symbol_import(self, arg: bool, /) -> None: ...

    @property
    def debug_on(self) -> bool:
        """For debugging."""

    @debug_on.setter
    def debug_on(self, arg: bool, /) -> None: ...

    @property
    def expression_tree_on(self) -> bool:
        """Use an expression tree instead of a stack?"""

    @expression_tree_on.setter
    def expression_tree_on(self, arg: bool, /) -> None: ...

    @property
    def verbose_on(self) -> bool:
        """For verbose output. Used with debugging."""

    @verbose_on.setter
    def verbose_on(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoGlobalStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoGlobalStruct: ...

    def __eq__(self, arg: TaoGlobalStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoGraphArrayStruct:
    """Fortran struct: tao_graph_array_struct"""

    def __init__(self, g: TaoGraphStruct | None = None) -> None: ...

    @property
    def g(self) -> TaoGraphStruct | None: ...

    @g.setter
    def g(self, arg: TaoGraphStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoGraphArrayStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoGraphArrayStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoGraphArrayStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoGraphArrayStruct: ...

    def __eq__(self, arg: TaoGraphArrayStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoGraphStruct:
    """Fortran struct: tao_graph_struct"""

    def __init__(self, name: str | None = None, type: str | None = None, title: str | None = None, title_suffix: str | None = None, why_invalid: str | None = None, p: TaoPlotStruct | None = None, floor_plan: TaoFloorPlanStruct | None = None, text_legend_origin: QpPointStruct | None = None, curve_legend_origin: QpPointStruct | None = None, curve_legend: QpLegendStruct | None = None, x: QpAxisStruct | None = None, y: QpAxisStruct | None = None, x2: QpAxisStruct | None = None, y2: QpAxisStruct | None = None, margin: QpRectStruct | None = None, scale_margin: QpRectStruct | None = None, x_axis_scale_factor: float | None = None, symbol_size_scale: float | None = None, box: Sequence[int] | None = None, ix_branch: int | None = None, ix_universe: int | None = None, clip: bool | None = None, y2_mirrors_y: bool | None = None, limited: bool | None = None, draw_axes: bool | None = None, draw_curve_legend: bool | None = None, draw_grid: bool | None = None, draw_title: bool | None = None, draw_only_good_user_data_or_vars: bool | None = None, allow_wrap_around: bool | None = None, is_valid: bool | None = None) -> None: ...

    @property
    def name(self) -> str:
        """Name identifying the graph"""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def type(self) -> str:
        """'data', 'lat_layout', 'phase_space', 'histogram', 'dynamic_aperture'"""

    @type.setter
    def type(self, arg: str, /) -> None: ...

    @property
    def title(self) -> str: ...

    @title.setter
    def title(self, arg: str, /) -> None: ...

    @property
    def title_suffix(self) -> str: ...

    @title_suffix.setter
    def title_suffix(self, arg: str, /) -> None: ...

    @property
    def text_legend(self) -> FCharArray1D:
        """Array for holding descriptive info."""

    @property
    def text_legend_out(self) -> FCharArray1D:
        """Array for holding descriptive info."""

    @property
    def why_invalid(self) -> str:
        """Informative string to print."""

    @why_invalid.setter
    def why_invalid(self, arg: str, /) -> None: ...

    @property
    def curve(self) -> TaoCurveStructAlloc1D: ...

    @property
    def p(self) -> TaoPlotStruct | None:
        """pointer to parent plot"""

    @p.setter
    def p(self, arg: TaoPlotStruct, /) -> None: ...

    @property
    def floor_plan(self) -> TaoFloorPlanStruct: ...

    @floor_plan.setter
    def floor_plan(self, arg: TaoFloorPlanStruct, /) -> None: ...

    @property
    def text_legend_origin(self) -> QpPointStruct: ...

    @text_legend_origin.setter
    def text_legend_origin(self, arg: QpPointStruct, /) -> None: ...

    @property
    def curve_legend_origin(self) -> QpPointStruct: ...

    @curve_legend_origin.setter
    def curve_legend_origin(self, arg: QpPointStruct, /) -> None: ...

    @property
    def curve_legend(self) -> QpLegendStruct: ...

    @curve_legend.setter
    def curve_legend(self, arg: QpLegendStruct, /) -> None: ...

    @property
    def x(self) -> QpAxisStruct:
        """X-axis parameters."""

    @x.setter
    def x(self, arg: QpAxisStruct, /) -> None: ...

    @property
    def y(self) -> QpAxisStruct:
        """Y-axis attributes."""

    @y.setter
    def y(self, arg: QpAxisStruct, /) -> None: ...

    @property
    def x2(self) -> QpAxisStruct:
        """X2-axis attributes (Not currently used)."""

    @x2.setter
    def x2(self, arg: QpAxisStruct, /) -> None: ...

    @property
    def y2(self) -> QpAxisStruct:
        """Y2-axis attributes."""

    @y2.setter
    def y2(self, arg: QpAxisStruct, /) -> None: ...

    @property
    def margin(self) -> QpRectStruct:
        """Margin around the graph."""

    @margin.setter
    def margin(self, arg: QpRectStruct, /) -> None: ...

    @property
    def scale_margin(self) -> QpRectStruct:
        """Margin for scaling"""

    @scale_margin.setter
    def scale_margin(self, arg: QpRectStruct, /) -> None: ...

    @property
    def x_axis_scale_factor(self) -> float:
        """x-axis conversion from internal to plotting units."""

    @x_axis_scale_factor.setter
    def x_axis_scale_factor(self, arg: float, /) -> None: ...

    @property
    def symbol_size_scale(self) -> float:
        """Symbol size scale factor for phase_space plots."""

    @symbol_size_scale.setter
    def symbol_size_scale(self, arg: float, /) -> None: ...

    @property
    def box(self) -> IntArray1D:
        """Defines which box the plot is put in."""

    @box.setter
    def box(self, arg: Sequence[int], /) -> None: ...

    @property
    def ix_branch(self) -> int:
        """Branch in lattice. Used when there are no associated curves."""

    @ix_branch.setter
    def ix_branch(self, arg: int, /) -> None: ...

    @property
    def ix_universe(self) -> int:
        """Used for lat_layout plots."""

    @ix_universe.setter
    def ix_universe(self, arg: int, /) -> None: ...

    @property
    def clip(self) -> bool:
        """Clip plot at graph boundary."""

    @clip.setter
    def clip(self, arg: bool, /) -> None: ...

    @property
    def y2_mirrors_y(self) -> bool:
        """Y2-axis same as Y-axis?"""

    @y2_mirrors_y.setter
    def y2_mirrors_y(self, arg: bool, /) -> None: ...

    @property
    def limited(self) -> bool:
        """True if at least one data point past graph bounds."""

    @limited.setter
    def limited(self, arg: bool, /) -> None: ...

    @property
    def draw_axes(self) -> bool:
        """Draw axes, labels, etc?"""

    @draw_axes.setter
    def draw_axes(self, arg: bool, /) -> None: ...

    @property
    def draw_curve_legend(self) -> bool:
        """Legend for displaying curve info."""

    @draw_curve_legend.setter
    def draw_curve_legend(self, arg: bool, /) -> None: ...

    @property
    def draw_grid(self) -> bool:
        """Draw a grid?"""

    @draw_grid.setter
    def draw_grid(self, arg: bool, /) -> None: ...

    @property
    def draw_title(self) -> bool: ...

    @draw_title.setter
    def draw_title(self, arg: bool, /) -> None: ...

    @property
    def draw_only_good_user_data_or_vars(self) -> bool: ...

    @draw_only_good_user_data_or_vars.setter
    def draw_only_good_user_data_or_vars(self, arg: bool, /) -> None: ...

    @property
    def allow_wrap_around(self) -> bool:
        """'Wrap' curves to extend past lattice boundaries?"""

    @allow_wrap_around.setter
    def allow_wrap_around(self, arg: bool, /) -> None: ...

    @property
    def is_valid(self) -> bool:
        """EG: Bad x_axis_type."""

    @is_valid.setter
    def is_valid(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoGraphStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoGraphStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoGraphStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoGraphStruct: ...

    def __eq__(self, arg: TaoGraphStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoHistogramStruct:
    """Fortran struct: tao_histogram_struct"""

    def __init__(self, density_normalized: bool | None = None, weight_by_charge: bool | None = None, minimum: float | None = None, maximum: float | None = None, width: float | None = None, center: float | None = None, number: int | None = None) -> None: ...

    @property
    def density_normalized(self) -> bool: ...

    @density_normalized.setter
    def density_normalized(self, arg: bool, /) -> None: ...

    @property
    def weight_by_charge(self) -> bool: ...

    @weight_by_charge.setter
    def weight_by_charge(self, arg: bool, /) -> None: ...

    @property
    def minimum(self) -> float:
        """Computed by Tao. Not User settable."""

    @minimum.setter
    def minimum(self, arg: float, /) -> None: ...

    @property
    def maximum(self) -> float:
        """Computed by Tao. Not User settable."""

    @maximum.setter
    def maximum(self, arg: float, /) -> None: ...

    @property
    def width(self) -> float: ...

    @width.setter
    def width(self, arg: float, /) -> None: ...

    @property
    def center(self) -> float: ...

    @center.setter
    def center(self, arg: float, /) -> None: ...

    @property
    def number(self) -> int: ...

    @number.setter
    def number(self, arg: int, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoHistogramStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoHistogramStruct: ...

    def __eq__(self, arg: TaoHistogramStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoInitStruct:
    """Fortran struct: tao_init_struct"""

    def __init__(self, parse_cmd_args: bool | None = None, debug_switch: bool | None = None, external_plotting_switch: bool | None = None, init_name: str | None = None, hook_init_file: str | None = None, hook_lat_file: str | None = None, hook_beam_file: str | None = None, hook_data_file: str | None = None, hook_plot_file: str | None = None, hook_startup_file: str | None = None, hook_var_file: str | None = None, hook_building_wall_file: str | None = None, init_file_arg_path: str | None = None, lattice_file_arg: str | None = None, hook_init_file_arg: str | None = None, init_file_arg: str | None = None, beam_file_arg: str | None = None, beam_init_position_file_arg: str | None = None, command_arg: str | None = None, data_file_arg: str | None = None, plot_file_arg: str | None = None, startup_file_arg: str | None = None, var_file_arg: str | None = None, building_wall_file_arg: str | None = None, geometry_arg: str | None = None, slice_lattice_arg: str | None = None, start_branch_at_arg: str | None = None, log_startup_arg: str | None = None, no_stopping_arg: str | None = None, noplot_arg: str | None = None, no_rad_int_arg: str | None = None, reverse_arg: str | None = None, debug_arg: str | None = None, disable_smooth_line_calc_arg: str | None = None, rf_on_arg: str | None = None, prompt_color_arg: str | None = None, quiet_arg: str | None = None, noinit_arg: str | None = None, nostartup_arg: str | None = None, symbol_import_arg: str | None = None, unique_name_suffix: str | None = None) -> None: ...

    @property
    def parse_cmd_args(self) -> bool:
        """Used by custom programs to control Tao init"""

    @parse_cmd_args.setter
    def parse_cmd_args(self, arg: bool, /) -> None: ...

    @property
    def debug_switch(self) -> bool:
        """Is the '-debug' switch present?"""

    @debug_switch.setter
    def debug_switch(self, arg: bool, /) -> None: ...

    @property
    def external_plotting_switch(self) -> bool:
        """Is '-external_plotting' switch present?"""

    @external_plotting_switch.setter
    def external_plotting_switch(self, arg: bool, /) -> None: ...

    @property
    def init_name(self) -> str:
        """label for initialization"""

    @init_name.setter
    def init_name(self, arg: str, /) -> None: ...

    @property
    def hook_init_file(self) -> str: ...

    @hook_init_file.setter
    def hook_init_file(self, arg: str, /) -> None: ...

    @property
    def hook_lat_file(self) -> str:
        """To be set by tao_hook_parse_command_args"""

    @hook_lat_file.setter
    def hook_lat_file(self, arg: str, /) -> None: ...

    @property
    def hook_beam_file(self) -> str:
        """To be set by tao_hook_parse_command_args"""

    @hook_beam_file.setter
    def hook_beam_file(self, arg: str, /) -> None: ...

    @property
    def hook_data_file(self) -> str:
        """To be set by tao_hook_parse_command_args"""

    @hook_data_file.setter
    def hook_data_file(self, arg: str, /) -> None: ...

    @property
    def hook_plot_file(self) -> str:
        """To be set by tao_hook_parse_command_args"""

    @hook_plot_file.setter
    def hook_plot_file(self, arg: str, /) -> None: ...

    @property
    def hook_startup_file(self) -> str:
        """To be set by tao_hook_parse_command_args"""

    @hook_startup_file.setter
    def hook_startup_file(self, arg: str, /) -> None: ...

    @property
    def hook_var_file(self) -> str:
        """To be set by tao_hook_parse_command_args"""

    @hook_var_file.setter
    def hook_var_file(self, arg: str, /) -> None: ...

    @property
    def hook_building_wall_file(self) -> str:
        """To be set by tao_hook_parse_command_args"""

    @hook_building_wall_file.setter
    def hook_building_wall_file(self, arg: str, /) -> None: ...

    @property
    def init_file_arg_path(self) -> str:
        """Path part of init_tao_file"""

    @init_file_arg_path.setter
    def init_file_arg_path(self, arg: str, /) -> None: ...

    @property
    def lattice_file_arg(self) -> str:
        """-lattice_file        command line argument."""

    @lattice_file_arg.setter
    def lattice_file_arg(self, arg: str, /) -> None: ...

    @property
    def hook_init_file_arg(self) -> str:
        """-hook_init_file      command line argument"""

    @hook_init_file_arg.setter
    def hook_init_file_arg(self, arg: str, /) -> None: ...

    @property
    def init_file_arg(self) -> str:
        """-init_file           command line argument."""

    @init_file_arg.setter
    def init_file_arg(self, arg: str, /) -> None: ...

    @property
    def beam_file_arg(self) -> str:
        """-beam_file           command line argument."""

    @beam_file_arg.setter
    def beam_file_arg(self, arg: str, /) -> None: ...

    @property
    def beam_init_position_file_arg(self) -> str:
        """-beam_init_position_file command line argument."""

    @beam_init_position_file_arg.setter
    def beam_init_position_file_arg(self, arg: str, /) -> None: ...

    @property
    def command_arg(self) -> str:
        """-command             command line argument."""

    @command_arg.setter
    def command_arg(self, arg: str, /) -> None: ...

    @property
    def data_file_arg(self) -> str:
        """-data_file           command line argument."""

    @data_file_arg.setter
    def data_file_arg(self, arg: str, /) -> None: ...

    @property
    def plot_file_arg(self) -> str:
        """-plot_file           command line argument."""

    @plot_file_arg.setter
    def plot_file_arg(self, arg: str, /) -> None: ...

    @property
    def startup_file_arg(self) -> str:
        """-startup_file        command line argument."""

    @startup_file_arg.setter
    def startup_file_arg(self, arg: str, /) -> None: ...

    @property
    def var_file_arg(self) -> str:
        """-var_file            command line argument."""

    @var_file_arg.setter
    def var_file_arg(self, arg: str, /) -> None: ...

    @property
    def building_wall_file_arg(self) -> str:
        """-building_wall_file  command line argument."""

    @building_wall_file_arg.setter
    def building_wall_file_arg(self, arg: str, /) -> None: ...

    @property
    def geometry_arg(self) -> str:
        """-geometry            command line argument."""

    @geometry_arg.setter
    def geometry_arg(self, arg: str, /) -> None: ...

    @property
    def slice_lattice_arg(self) -> str:
        """-slice_lattice       command line argument."""

    @slice_lattice_arg.setter
    def slice_lattice_arg(self, arg: str, /) -> None: ...

    @property
    def start_branch_at_arg(self) -> str:
        """-start_branch_at     command line argument."""

    @start_branch_at_arg.setter
    def start_branch_at_arg(self, arg: str, /) -> None: ...

    @property
    def log_startup_arg(self) -> str:
        """-log_startup         command line argument"""

    @log_startup_arg.setter
    def log_startup_arg(self, arg: str, /) -> None: ...

    @property
    def no_stopping_arg(self) -> str:
        """-no_stopping         command line argument"""

    @no_stopping_arg.setter
    def no_stopping_arg(self, arg: str, /) -> None: ...

    @property
    def noplot_arg(self) -> str:
        """-noplot              command line argument"""

    @noplot_arg.setter
    def noplot_arg(self, arg: str, /) -> None: ...

    @property
    def no_rad_int_arg(self) -> str:
        """-no_rad_int          command line argument"""

    @no_rad_int_arg.setter
    def no_rad_int_arg(self, arg: str, /) -> None: ...

    @property
    def reverse_arg(self) -> str:
        """-reverse             command line argument"""

    @reverse_arg.setter
    def reverse_arg(self, arg: str, /) -> None: ...

    @property
    def debug_arg(self) -> str:
        """-debug               command line argument"""

    @debug_arg.setter
    def debug_arg(self, arg: str, /) -> None: ...

    @property
    def disable_smooth_line_calc_arg(self) -> str:
        """-disable_smooth_line_calc"""

    @disable_smooth_line_calc_arg.setter
    def disable_smooth_line_calc_arg(self, arg: str, /) -> None: ...

    @property
    def rf_on_arg(self) -> str:
        """-rf_on               command line argument"""

    @rf_on_arg.setter
    def rf_on_arg(self, arg: str, /) -> None: ...

    @property
    def prompt_color_arg(self) -> str:
        """-prompt_color        command line argument"""

    @prompt_color_arg.setter
    def prompt_color_arg(self, arg: str, /) -> None: ...

    @property
    def quiet_arg(self) -> str:
        """-quiet               command line argument"""

    @quiet_arg.setter
    def quiet_arg(self, arg: str, /) -> None: ...

    @property
    def noinit_arg(self) -> str:
        """-noinit              command line argument"""

    @noinit_arg.setter
    def noinit_arg(self, arg: str, /) -> None: ...

    @property
    def nostartup_arg(self) -> str:
        """-nostartup           command line argument"""

    @nostartup_arg.setter
    def nostartup_arg(self, arg: str, /) -> None: ...

    @property
    def symbol_import_arg(self) -> str:
        """-symbol_import       command line argument"""

    @symbol_import_arg.setter
    def symbol_import_arg(self, arg: str, /) -> None: ...

    @property
    def unique_name_suffix(self) -> str: ...

    @unique_name_suffix.setter
    def unique_name_suffix(self, arg: str, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoInitStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoInitStruct: ...

    def __eq__(self, arg: TaoInitStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoIntegerArrayStruct:
    """Fortran struct: tao_integer_array_struct"""

    def __init__(self, i: int | None = None) -> None: ...

    @property
    def i(self) -> int | None: ...

    @i.setter
    def i(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoIntegerArrayStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoIntegerArrayStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoIntegerArrayStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoIntegerArrayStruct: ...

    def __eq__(self, arg: TaoIntegerArrayStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoLatSigmaStruct:
    """Fortran struct: tao_lat_sigma_struct"""

    def __init__(self, mat: Sequence[Sequence[float]] | None = None) -> None: ...

    @property
    def mat(self) -> RealArray2D: ...

    @mat.setter
    def mat(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoLatSigmaStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoLatSigmaStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoLatSigmaStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoLatSigmaStruct: ...

    def __eq__(self, arg: TaoLatSigmaStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoLatticeBranchStruct:
    """Fortran struct: tao_lattice_branch_struct"""

    def __init__(self, tao_lat: TaoLatticeStruct | None = None, spin: TaoSpinPolarizationStruct | None = None, srdt: SummationRdtStruct | None = None, orb0: CoordStruct | None = None, modes_ri: NormalModesStruct | None = None, modes_6d: NormalModesStruct | None = None, ptc_normal_form: PtcNormalFormStruct | None = None, bmad_normal_form: BmadNormalFormStruct | None = None, cache_x_min: float | None = None, cache_x_max: float | None = None, comb_ds_save: float | None = None, ix_ref_taylor: int | None = None, ix_ele_taylor: int | None = None, track_state: int | None = None, cache_n_pts: int | None = None, ix_rad_int_cache: int | None = None, has_open_match_element: bool | None = None, plot_cache_valid: bool | None = None, spin_map_valid: bool | None = None, twiss_valid: bool | None = None, mode_flip_here: bool | None = None, chrom_calc_ok: bool | None = None, rad_int_calc_ok: bool | None = None, emit_6d_calc_ok: bool | None = None, sigma_track_ok: bool | None = None) -> None: ...

    @property
    def tao_lat(self) -> TaoLatticeStruct | None:
        """Parent tao_lat"""

    @tao_lat.setter
    def tao_lat(self, arg: TaoLatticeStruct, /) -> None: ...

    @property
    def lat_sigma(self) -> TaoLatSigmaStructAlloc1D:
        """Sigma matrix derived from lattice (not beam)."""

    @property
    def spin_ele(self) -> TaoSpinEleStructAlloc1D:
        """Spin stuff"""

    @property
    def bunch_params(self) -> BunchParamsStructAlloc1D:
        """Per element"""

    @property
    def bunch_params_comb(self) -> BunchTrackStructAlloc1D:
        """A comb for each bunch in beam."""

    @property
    def orbit(self) -> CoordStructAlloc1D: ...

    @property
    def plot_cache(self) -> TaoPlotCacheStructAlloc1D:
        """Plotting data cache"""

    @property
    def spin(self) -> TaoSpinPolarizationStruct: ...

    @spin.setter
    def spin(self, arg: TaoSpinPolarizationStruct, /) -> None: ...

    @property
    def srdt(self) -> SummationRdtStruct: ...

    @srdt.setter
    def srdt(self, arg: SummationRdtStruct, /) -> None: ...

    @property
    def orb0(self) -> CoordStruct:
        """
        For saving beginning orbit in closed geometry branches. orb0 can then be used as an initial guess when closed_orbit is called again.
        """

    @orb0.setter
    def orb0(self, arg: CoordStruct, /) -> None: ...

    @property
    def modes_ri(self) -> NormalModesStruct:
        """Synchrotron integrals stuff"""

    @modes_ri.setter
    def modes_ri(self, arg: NormalModesStruct, /) -> None: ...

    @property
    def modes_6d(self) -> NormalModesStruct:
        """6D radiation matrices."""

    @modes_6d.setter
    def modes_6d(self, arg: NormalModesStruct, /) -> None: ...

    @property
    def ptc_normal_form(self) -> PtcNormalFormStruct:
        """Collection of normal form structures defined in PTC"""

    @ptc_normal_form.setter
    def ptc_normal_form(self, arg: PtcNormalFormStruct, /) -> None: ...

    @property
    def bmad_normal_form(self) -> BmadNormalFormStruct:
        """Collection of normal form structures defined in Bmad"""

    @bmad_normal_form.setter
    def bmad_normal_form(self, arg: BmadNormalFormStruct, /) -> None: ...

    @property
    def high_E_orb(self) -> CoordStructAlloc1D: ...

    @property
    def low_E_orb(self) -> CoordStructAlloc1D: ...

    @property
    def taylor_save(self) -> TaylorStructArray1D:
        """Save to reduce computation time."""

    @property
    def cache_x_min(self) -> float: ...

    @cache_x_min.setter
    def cache_x_min(self, arg: float, /) -> None: ...

    @property
    def cache_x_max(self) -> float: ...

    @cache_x_max.setter
    def cache_x_max(self, arg: float, /) -> None: ...

    @property
    def comb_ds_save(self) -> float:
        """Master parameter for %bunch_params_comb(:)%ds_save"""

    @comb_ds_save.setter
    def comb_ds_save(self, arg: float, /) -> None: ...

    @property
    def ix_ref_taylor(self) -> int: ...

    @ix_ref_taylor.setter
    def ix_ref_taylor(self, arg: int, /) -> None: ...

    @property
    def ix_ele_taylor(self) -> int: ...

    @ix_ele_taylor.setter
    def ix_ele_taylor(self, arg: int, /) -> None: ...

    @property
    def track_state(self) -> int: ...

    @track_state.setter
    def track_state(self, arg: int, /) -> None: ...

    @property
    def cache_n_pts(self) -> int: ...

    @cache_n_pts.setter
    def cache_n_pts(self, arg: int, /) -> None: ...

    @property
    def ix_rad_int_cache(self) -> int:
        """Radiation integrals cache index."""

    @ix_rad_int_cache.setter
    def ix_rad_int_cache(self, arg: int, /) -> None: ...

    @property
    def has_open_match_element(self) -> bool: ...

    @has_open_match_element.setter
    def has_open_match_element(self, arg: bool, /) -> None: ...

    @property
    def plot_cache_valid(self) -> bool:
        """Valid plotting data cache?"""

    @plot_cache_valid.setter
    def plot_cache_valid(self, arg: bool, /) -> None: ...

    @property
    def spin_map_valid(self) -> bool: ...

    @spin_map_valid.setter
    def spin_map_valid(self, arg: bool, /) -> None: ...

    @property
    def twiss_valid(self) -> bool:
        """
        Invalid EG with unstable 1-turn matrix with a closed branch. With open branch: twiss_valid = T even if some Twiss (and orbit) is invalid.
        """

    @twiss_valid.setter
    def twiss_valid(self, arg: bool, /) -> None: ...

    @property
    def mode_flip_here(self) -> bool:
        """Twiss parameter mode flip seen?"""

    @mode_flip_here.setter
    def mode_flip_here(self, arg: bool, /) -> None: ...

    @property
    def chrom_calc_ok(self) -> bool: ...

    @chrom_calc_ok.setter
    def chrom_calc_ok(self, arg: bool, /) -> None: ...

    @property
    def rad_int_calc_ok(self) -> bool: ...

    @rad_int_calc_ok.setter
    def rad_int_calc_ok(self, arg: bool, /) -> None: ...

    @property
    def emit_6d_calc_ok(self) -> bool: ...

    @emit_6d_calc_ok.setter
    def emit_6d_calc_ok(self, arg: bool, /) -> None: ...

    @property
    def sigma_track_ok(self) -> bool: ...

    @sigma_track_ok.setter
    def sigma_track_ok(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoLatticeBranchStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoLatticeBranchStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoLatticeBranchStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoLatticeBranchStruct: ...

    def __eq__(self, arg: TaoLatticeBranchStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoLatticeStruct:
    """Fortran struct: tao_lattice_struct"""

    def __init__(self, name: str | None = None, lat: LatStruct | None = None, high_E_lat: LatStruct | None = None, low_E_lat: LatStruct | None = None, u: TaoUniverseStruct | None = None, rad_int_by_ele_ri: RadIntAllEleStruct | None = None, rad_int_by_ele_6d: RadIntAllEleStruct | None = None) -> None: ...

    @property
    def name(self) -> str:
        """'model', 'base', or 'design'."""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def lat(self) -> LatStruct:
        """lattice structures"""

    @lat.setter
    def lat(self, arg: LatStruct, /) -> None: ...

    @property
    def high_E_lat(self) -> LatStruct:
        """For chrom calc."""

    @high_E_lat.setter
    def high_E_lat(self, arg: LatStruct, /) -> None: ...

    @property
    def low_E_lat(self) -> LatStruct:
        """For chrom calc."""

    @low_E_lat.setter
    def low_E_lat(self, arg: LatStruct, /) -> None: ...

    @property
    def u(self) -> TaoUniverseStruct | None:
        """Parent universe"""

    @u.setter
    def u(self, arg: TaoUniverseStruct, /) -> None: ...

    @property
    def rad_int_by_ele_ri(self) -> RadIntAllEleStruct: ...

    @rad_int_by_ele_ri.setter
    def rad_int_by_ele_ri(self, arg: RadIntAllEleStruct, /) -> None: ...

    @property
    def rad_int_by_ele_6d(self) -> RadIntAllEleStruct: ...

    @rad_int_by_ele_6d.setter
    def rad_int_by_ele_6d(self, arg: RadIntAllEleStruct, /) -> None: ...

    @property
    def tao_branch(self) -> TaoLatticeBranchStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoLatticeStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoLatticeStruct: ...

    def __eq__(self, arg: TaoLatticeStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoLogicalArrayStruct:
    """Fortran struct: tao_logical_array_struct"""

    def __init__(self, l: bool | None = None) -> None: ...

    @property
    def l(self) -> bool | None: ...

    @l.setter
    def l(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoLogicalArrayStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoLogicalArrayStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoLogicalArrayStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoLogicalArrayStruct: ...

    def __eq__(self, arg: TaoLogicalArrayStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoModelBranchStruct:
    """Fortran struct: tao_model_branch_struct"""

    def __init__(self, beam: TaoBeamBranchStruct | None = None) -> None: ...

    @property
    def ele(self) -> TaoModelElementStructAlloc1D:
        """Per element information"""

    @property
    def beam(self) -> TaoBeamBranchStruct: ...

    @beam.setter
    def beam(self, arg: TaoBeamBranchStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoModelBranchStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoModelBranchStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoModelBranchStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoModelBranchStruct: ...

    def __eq__(self, arg: TaoModelBranchStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoModelElementStruct:
    """Fortran struct: tao_model_element_struct"""

    def __init__(self, beam: BeamStruct | None = None, save_beam_internally: bool | None = None, save_beam_to_file: bool | None = None) -> None: ...

    @property
    def beam(self) -> BeamStruct:
        """Beam distribution at element."""

    @beam.setter
    def beam(self, arg: BeamStruct, /) -> None: ...

    @property
    def save_beam_internally(self) -> bool:
        """Save beam here? Beam also saved at fork elements and at track ends."""

    @save_beam_internally.setter
    def save_beam_internally(self, arg: bool, /) -> None: ...

    @property
    def save_beam_to_file(self) -> bool:
        """
        Save beam to a file? Beam also saved at fork elements and at track ends.
        """

    @save_beam_to_file.setter
    def save_beam_to_file(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoModelElementStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoModelElementStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoModelElementStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoModelElementStruct: ...

    def __eq__(self, arg: TaoModelElementStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoPingScaleStruct:
    """Fortran struct: tao_ping_scale_struct"""

    def __init__(self, a_mode_meas: float | None = None, a_mode_ref: float | None = None, b_mode_meas: float | None = None, b_mode_ref: float | None = None) -> None: ...

    @property
    def a_mode_meas(self) -> float: ...

    @a_mode_meas.setter
    def a_mode_meas(self, arg: float, /) -> None: ...

    @property
    def a_mode_ref(self) -> float: ...

    @a_mode_ref.setter
    def a_mode_ref(self, arg: float, /) -> None: ...

    @property
    def b_mode_meas(self) -> float: ...

    @b_mode_meas.setter
    def b_mode_meas(self, arg: float, /) -> None: ...

    @property
    def b_mode_ref(self) -> float: ...

    @b_mode_ref.setter
    def b_mode_ref(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoPingScaleStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoPingScaleStruct: ...

    def __eq__(self, arg: TaoPingScaleStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoPlotArrayStruct:
    """Fortran struct: tao_plot_array_struct"""

    def __init__(self, p: TaoPlotStruct | None = None) -> None: ...

    @property
    def p(self) -> TaoPlotStruct | None: ...

    @p.setter
    def p(self, arg: TaoPlotStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoPlotArrayStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoPlotArrayStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoPlotArrayStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoPlotArrayStruct: ...

    def __eq__(self, arg: TaoPlotArrayStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoPlotCacheStruct:
    """Fortran struct: tao_plot_cache_struct"""

    def __init__(self, ele_to_s: EleStruct | None = None, orbit: CoordStruct | None = None, err: bool | None = None) -> None: ...

    @property
    def ele_to_s(self) -> EleStruct:
        """
        Integrated element from branch beginning. Will be marked as a hybrid element.
        """

    @ele_to_s.setter
    def ele_to_s(self, arg: EleStruct, /) -> None: ...

    @property
    def orbit(self) -> CoordStruct: ...

    @orbit.setter
    def orbit(self, arg: CoordStruct, /) -> None: ...

    @property
    def err(self) -> bool: ...

    @err.setter
    def err(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoPlotCacheStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoPlotCacheStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoPlotCacheStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoPlotCacheStruct: ...

    def __eq__(self, arg: TaoPlotCacheStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoPlotPageInput:
    """Fortran struct: tao_plot_page_input"""

    def __init__(self, title: TaoTitleStruct | None = None, subtitle: TaoTitleStruct | None = None, border: QpRectStruct | None = None, plot_display_type: str | None = None, size: Sequence[float] | None = None, text_height: float | None = None, main_title_text_scale: float | None = None, graph_title_text_scale: float | None = None, axis_number_text_scale: float | None = None, axis_label_text_scale: float | None = None, legend_text_scale: float | None = None, key_table_text_scale: float | None = None, floor_plan_shape_scale: float | None = None, floor_plan_text_scale: float | None = None, lat_layout_shape_scale: float | None = None, lat_layout_text_scale: float | None = None, curve_legend_line_len: float | None = None, curve_legend_text_offset: float | None = None, n_curve_pts: int | None = None, delete_overlapping_plots: bool | None = None, draw_graph_title_suffix: bool | None = None) -> None: ...

    @property
    def title(self) -> TaoTitleStruct:
        """Title  at top of page."""

    @title.setter
    def title(self, arg: TaoTitleStruct, /) -> None: ...

    @property
    def subtitle(self) -> TaoTitleStruct:
        """Subtitle at top of page."""

    @subtitle.setter
    def subtitle(self, arg: TaoTitleStruct, /) -> None: ...

    @property
    def border(self) -> QpRectStruct:
        """Border around plots edge of page."""

    @border.setter
    def border(self, arg: QpRectStruct, /) -> None: ...

    @property
    def plot_display_type(self) -> str: ...

    @plot_display_type.setter
    def plot_display_type(self, arg: str, /) -> None: ...

    @property
    def size(self) -> RealArray1D:
        """width and height of window in pixels."""

    @size.setter
    def size(self, arg: Sequence[float], /) -> None: ...

    @property
    def text_height(self) -> float:
        """In points. Scales the height of all text"""

    @text_height.setter
    def text_height(self, arg: float, /) -> None: ...

    @property
    def main_title_text_scale(self) -> float:
        """Relative to text_height"""

    @main_title_text_scale.setter
    def main_title_text_scale(self, arg: float, /) -> None: ...

    @property
    def graph_title_text_scale(self) -> float:
        """Relative to text_height"""

    @graph_title_text_scale.setter
    def graph_title_text_scale(self, arg: float, /) -> None: ...

    @property
    def axis_number_text_scale(self) -> float:
        """Relative to text_height"""

    @axis_number_text_scale.setter
    def axis_number_text_scale(self, arg: float, /) -> None: ...

    @property
    def axis_label_text_scale(self) -> float:
        """Relative to text_height"""

    @axis_label_text_scale.setter
    def axis_label_text_scale(self, arg: float, /) -> None: ...

    @property
    def legend_text_scale(self) -> float:
        """Relative to text_height"""

    @legend_text_scale.setter
    def legend_text_scale(self, arg: float, /) -> None: ...

    @property
    def key_table_text_scale(self) -> float:
        """Relative to text_height"""

    @key_table_text_scale.setter
    def key_table_text_scale(self, arg: float, /) -> None: ...

    @property
    def floor_plan_shape_scale(self) -> float: ...

    @floor_plan_shape_scale.setter
    def floor_plan_shape_scale(self, arg: float, /) -> None: ...

    @property
    def floor_plan_text_scale(self) -> float:
        """Scale used = floor_plan_text_scale * legend_text_scale"""

    @floor_plan_text_scale.setter
    def floor_plan_text_scale(self, arg: float, /) -> None: ...

    @property
    def lat_layout_shape_scale(self) -> float: ...

    @lat_layout_shape_scale.setter
    def lat_layout_shape_scale(self, arg: float, /) -> None: ...

    @property
    def lat_layout_text_scale(self) -> float:
        """Scale used = lat_layout_text_scale * legend_text_scale"""

    @lat_layout_text_scale.setter
    def lat_layout_text_scale(self, arg: float, /) -> None: ...

    @property
    def curve_legend_line_len(self) -> float:
        """OLD STYLE. Points."""

    @curve_legend_line_len.setter
    def curve_legend_line_len(self, arg: float, /) -> None: ...

    @property
    def curve_legend_text_offset(self) -> float:
        """OLD STYLE. Points."""

    @curve_legend_text_offset.setter
    def curve_legend_text_offset(self, arg: float, /) -> None: ...

    @property
    def n_curve_pts(self) -> int:
        """Number of points for plotting a smooth curve"""

    @n_curve_pts.setter
    def n_curve_pts(self, arg: int, /) -> None: ...

    @property
    def delete_overlapping_plots(self) -> bool:
        """Delete overlapping plots when a plot is placed?"""

    @delete_overlapping_plots.setter
    def delete_overlapping_plots(self, arg: bool, /) -> None: ...

    @property
    def draw_graph_title_suffix(self) -> bool: ...

    @draw_graph_title_suffix.setter
    def draw_graph_title_suffix(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoPlotPageInput: ...

    def __deepcopy__(self, arg: dict, /) -> TaoPlotPageInput: ...

    def __eq__(self, arg: TaoPlotPageInput, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoPlotPageStruct:
    """Fortran struct: tao_plot_page_struct"""

    def __init__(self, title: TaoTitleStruct | None = None, subtitle: TaoTitleStruct | None = None, border: QpRectStruct | None = None, floor_plan: TaoDrawingStruct | None = None, lat_layout: TaoDrawingStruct | None = None, plot_display_type: str | None = None, size: Sequence[float] | None = None, text_height: float | None = None, main_title_text_scale: float | None = None, graph_title_text_scale: float | None = None, axis_number_text_scale: float | None = None, axis_label_text_scale: float | None = None, legend_text_scale: float | None = None, key_table_text_scale: float | None = None, floor_plan_shape_scale: float | None = None, floor_plan_text_scale: float | None = None, lat_layout_shape_scale: float | None = None, lat_layout_text_scale: float | None = None, n_curve_pts: int | None = None, id_window: int | None = None, delete_overlapping_plots: bool | None = None, draw_graph_title_suffix: bool | None = None) -> None: ...

    @property
    def title(self) -> TaoTitleStruct:
        """Title  at top of page."""

    @title.setter
    def title(self, arg: TaoTitleStruct, /) -> None: ...

    @property
    def subtitle(self) -> TaoTitleStruct:
        """Subtitle below title at top of page."""

    @subtitle.setter
    def subtitle(self, arg: TaoTitleStruct, /) -> None: ...

    @property
    def border(self) -> QpRectStruct:
        """Border around plots edge of page."""

    @border.setter
    def border(self, arg: QpRectStruct, /) -> None: ...

    @property
    def floor_plan(self) -> TaoDrawingStruct: ...

    @floor_plan.setter
    def floor_plan(self, arg: TaoDrawingStruct, /) -> None: ...

    @property
    def lat_layout(self) -> TaoDrawingStruct: ...

    @lat_layout.setter
    def lat_layout(self, arg: TaoDrawingStruct, /) -> None: ...

    @property
    def pattern(self) -> TaoShapePatternStructAlloc1D: ...

    @property
    def region(self) -> TaoPlotRegionStructAlloc1D: ...

    @property
    def plot_display_type(self) -> str:
        """'X' or 'TK'"""

    @plot_display_type.setter
    def plot_display_type(self, arg: str, /) -> None: ...

    @property
    def size(self) -> RealArray1D:
        """width and height of plot window in pixels."""

    @size.setter
    def size(self, arg: Sequence[float], /) -> None: ...

    @property
    def text_height(self) -> float:
        """In points. Scales the height of all text"""

    @text_height.setter
    def text_height(self, arg: float, /) -> None: ...

    @property
    def main_title_text_scale(self) -> float:
        """Relative to text_height"""

    @main_title_text_scale.setter
    def main_title_text_scale(self, arg: float, /) -> None: ...

    @property
    def graph_title_text_scale(self) -> float:
        """Relative to text_height"""

    @graph_title_text_scale.setter
    def graph_title_text_scale(self, arg: float, /) -> None: ...

    @property
    def axis_number_text_scale(self) -> float:
        """Relative to text_height"""

    @axis_number_text_scale.setter
    def axis_number_text_scale(self, arg: float, /) -> None: ...

    @property
    def axis_label_text_scale(self) -> float:
        """Relative to text_height"""

    @axis_label_text_scale.setter
    def axis_label_text_scale(self, arg: float, /) -> None: ...

    @property
    def legend_text_scale(self) -> float:
        """Relative to text_height. For legends, plot_page, and lat_layout"""

    @legend_text_scale.setter
    def legend_text_scale(self, arg: float, /) -> None: ...

    @property
    def key_table_text_scale(self) -> float:
        """Relative to text_height"""

    @key_table_text_scale.setter
    def key_table_text_scale(self, arg: float, /) -> None: ...

    @property
    def floor_plan_shape_scale(self) -> float: ...

    @floor_plan_shape_scale.setter
    def floor_plan_shape_scale(self, arg: float, /) -> None: ...

    @property
    def floor_plan_text_scale(self) -> float:
        """Scale used = floor_plan_text_scale * legend_text_scale"""

    @floor_plan_text_scale.setter
    def floor_plan_text_scale(self, arg: float, /) -> None: ...

    @property
    def lat_layout_shape_scale(self) -> float: ...

    @lat_layout_shape_scale.setter
    def lat_layout_shape_scale(self, arg: float, /) -> None: ...

    @property
    def lat_layout_text_scale(self) -> float:
        """Scale used = lat_layout_text_scale * legend_text_scale"""

    @lat_layout_text_scale.setter
    def lat_layout_text_scale(self, arg: float, /) -> None: ...

    @property
    def n_curve_pts(self) -> int:
        """Default number of points for plotting a smooth curve."""

    @n_curve_pts.setter
    def n_curve_pts(self, arg: int, /) -> None: ...

    @property
    def id_window(self) -> int:
        """X window id number."""

    @id_window.setter
    def id_window(self, arg: int, /) -> None: ...

    @property
    def delete_overlapping_plots(self) -> bool:
        """Delete overlapping plots when a plot is placed?"""

    @delete_overlapping_plots.setter
    def delete_overlapping_plots(self, arg: bool, /) -> None: ...

    @property
    def draw_graph_title_suffix(self) -> bool:
        """Draw the graph title suffix?"""

    @draw_graph_title_suffix.setter
    def draw_graph_title_suffix(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoPlotPageStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoPlotPageStruct: ...

    def __eq__(self, arg: TaoPlotPageStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoPlotRegionStruct:
    """Fortran struct: tao_plot_region_struct"""

    def __init__(self, name: str | None = None, plot: TaoPlotStruct | None = None, location: Sequence[float] | None = None, visible: bool | None = None, list_with_show_plot_command: bool | None = None, setup_done: bool | None = None) -> None: ...

    @property
    def name(self) -> str:
        """Region name. Eg: 'r13', etc."""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def plot(self) -> TaoPlotStruct:
        """Plot associated with this region"""

    @plot.setter
    def plot(self, arg: TaoPlotStruct, /) -> None: ...

    @property
    def location(self) -> RealArray1D:
        """[x1, x2, y1, y2] location on page."""

    @location.setter
    def location(self, arg: Sequence[float], /) -> None: ...

    @property
    def visible(self) -> bool:
        """To draw or not to draw."""

    @visible.setter
    def visible(self, arg: bool, /) -> None: ...

    @property
    def list_with_show_plot_command(self) -> bool:
        """False used for default plots to shorten the output of 'show plot'"""

    @list_with_show_plot_command.setter
    def list_with_show_plot_command(self, arg: bool, /) -> None: ...

    @property
    def setup_done(self) -> bool:
        """Used for plot bookkeeping."""

    @setup_done.setter
    def setup_done(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoPlotRegionStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoPlotRegionStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoPlotRegionStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoPlotRegionStruct: ...

    def __eq__(self, arg: TaoPlotRegionStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoPlotStruct:
    """Fortran struct: tao_plot_struct"""

    def __init__(self, name: str | None = None, description: str | None = None, r: TaoPlotRegionStruct | None = None, ix_plot: int | None = None, n_curve_pts: int | None = None, type: str | None = None, x_axis_type: str | None = None, autoscale_x: bool | None = None, autoscale_y: bool | None = None, autoscale_gang_x: bool | None = None, autoscale_gang_y: bool | None = None, list_with_show_plot_command: bool | None = None, phantom: bool | None = None, default_plot: bool | None = None) -> None: ...

    @property
    def name(self) -> str:
        """Identifying name. Rule: If name is blank, plot is not valid."""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def description(self) -> str:
        """Descriptive string."""

    @description.setter
    def description(self, arg: str, /) -> None: ...

    @property
    def graph(self) -> TaoGraphStructAlloc1D:
        """individual graphs of a plot"""

    @property
    def r(self) -> TaoPlotRegionStruct | None:
        """pointer to parent."""

    @r.setter
    def r(self, arg: TaoPlotRegionStruct, /) -> None: ...

    @property
    def ix_plot(self) -> int:
        """Index in s%plot_page%template(:) or %region(:) arrays."""

    @ix_plot.setter
    def ix_plot(self, arg: int, /) -> None: ...

    @property
    def n_curve_pts(self) -> int:
        """Overrides s%plot_page%n_curve_pts."""

    @n_curve_pts.setter
    def n_curve_pts(self, arg: int, /) -> None: ...

    @property
    def type(self) -> str:
        """or 'wave'"""

    @type.setter
    def type(self, arg: str, /) -> None: ...

    @property
    def x_axis_type(self) -> str:
        """'index', 'ele_index', 's', 'none', 'floor', 'phase_space', etc."""

    @x_axis_type.setter
    def x_axis_type(self, arg: str, /) -> None: ...

    @property
    def autoscale_x(self) -> bool:
        """Horizontal autoscale."""

    @autoscale_x.setter
    def autoscale_x(self, arg: bool, /) -> None: ...

    @property
    def autoscale_y(self) -> bool:
        """Vertical autoscale."""

    @autoscale_y.setter
    def autoscale_y(self, arg: bool, /) -> None: ...

    @property
    def autoscale_gang_x(self) -> bool:
        """scale cmd scales graphs together?"""

    @autoscale_gang_x.setter
    def autoscale_gang_x(self, arg: bool, /) -> None: ...

    @property
    def autoscale_gang_y(self) -> bool:
        """scale cmd scales graphs together?"""

    @autoscale_gang_y.setter
    def autoscale_gang_y(self, arg: bool, /) -> None: ...

    @property
    def list_with_show_plot_command(self) -> bool:
        """False used for default plots to shorten the output of 'show plot'"""

    @list_with_show_plot_command.setter
    def list_with_show_plot_command(self, arg: bool, /) -> None: ...

    @property
    def phantom(self) -> bool:
        """Used by tao_plot_init to add info lines to 'show plot -templates'"""

    @phantom.setter
    def phantom(self, arg: bool, /) -> None: ...

    @property
    def default_plot(self) -> bool:
        """One of Tao's default plots?"""

    @default_plot.setter
    def default_plot(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoPlotStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoPlotStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoPlotStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoPlotStruct: ...

    def __eq__(self, arg: TaoPlotStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoRealPointerStruct:
    """Fortran struct: tao_real_pointer_struct"""

    def __init__(self, r: float | None = None, good_value: bool | None = None, good_user: bool | None = None) -> None: ...

    @property
    def r(self) -> float | None: ...

    @r.setter
    def r(self, arg: float, /) -> None: ...

    @property
    def good_value(self) -> bool | None: ...

    @good_value.setter
    def good_value(self, arg: bool, /) -> None: ...

    @property
    def good_user(self) -> bool | None: ...

    @good_user.setter
    def good_user(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoRealPointerStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoRealPointerStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoRealPointerStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoRealPointerStruct: ...

    def __eq__(self, arg: TaoRealPointerStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoShapePatternPointStruct:
    """Fortran struct: tao_shape_pattern_point_struct"""

    def __init__(self, s: float | None = None, y: float | None = None, radius: float | None = None) -> None: ...

    @property
    def s(self) -> float: ...

    @s.setter
    def s(self, arg: float, /) -> None: ...

    @property
    def y(self) -> float: ...

    @y.setter
    def y(self, arg: float, /) -> None: ...

    @property
    def radius(self) -> float: ...

    @radius.setter
    def radius(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoShapePatternPointStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoShapePatternPointStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoShapePatternPointStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoShapePatternPointStruct: ...

    def __eq__(self, arg: TaoShapePatternPointStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoShapePatternStruct:
    """Fortran struct: tao_shape_pattern_struct"""

    def __init__(self, name: str | None = None, line: QpLineStruct | None = None) -> None: ...

    @property
    def name(self) -> str: ...

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def line(self) -> QpLineStruct:
        """Line color and pattern set by shape using this pattern."""

    @line.setter
    def line(self, arg: QpLineStruct, /) -> None: ...

    @property
    def pt(self) -> TaoShapePatternPointStructAlloc1D: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoShapePatternStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoShapePatternStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoShapePatternStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoShapePatternStruct: ...

    def __eq__(self, arg: TaoShapePatternStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoSpinDnDpzStruct:
    """Fortran struct: tao_spin_dn_dpz_struct"""

    def __init__(self, vec: Sequence[float] | None = None, partial: Sequence[Sequence[float]] | None = None, partial2: Sequence[Sequence[float]] | None = None) -> None: ...

    @property
    def vec(self) -> RealArray1D:
        """n0 derivative wrt pz."""

    @vec.setter
    def vec(self, arg: Sequence[float], /) -> None: ...

    @property
    def partial(self) -> RealArray2D:
        """
        partial(i:) is spin n0 derivative wrt pz for i^th oscillation mode (1 => a-mode, etc.)
        """

    @partial.setter
    def partial(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def partial2(self) -> RealArray2D:
        """
        partial(i:) is spin n0 derivative wrt pz with i^th oscillation mode missing (1 => a-mode, etc.)
        """

    @partial2.setter
    def partial2(self, arg: Sequence[Sequence[float]], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoSpinDnDpzStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoSpinDnDpzStruct: ...

    def __eq__(self, arg: TaoSpinDnDpzStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoSpinEleStruct:
    """Fortran struct: tao_spin_ele_struct"""

    def __init__(self, dn_dpz: TaoSpinDnDpzStruct | None = None, orb_eigen_val: Sequence[float] | None = None, orb_eigen_vec: Sequence[Sequence[float]] | None = None, spin_eigen_vec: Sequence[Sequence[float]] | None = None, valid: bool | None = None) -> None: ...

    @property
    def dn_dpz(self) -> TaoSpinDnDpzStruct: ...

    @dn_dpz.setter
    def dn_dpz(self, arg: TaoSpinDnDpzStruct, /) -> None: ...

    @property
    def orb_eigen_val(self) -> RealArray1D: ...

    @orb_eigen_val.setter
    def orb_eigen_val(self, arg: Sequence[float], /) -> None: ...

    @property
    def orb_eigen_vec(self) -> RealArray2D:
        """(j,:) is j^th vector"""

    @orb_eigen_vec.setter
    def orb_eigen_vec(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def spin_eigen_vec(self) -> RealArray2D:
        """(j,:) is j^th vector"""

    @spin_eigen_vec.setter
    def spin_eigen_vec(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def valid(self) -> bool: ...

    @valid.setter
    def valid(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoSpinEleStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoSpinEleStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoSpinEleStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoSpinEleStruct: ...

    def __eq__(self, arg: TaoSpinEleStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoSpinMapStruct:
    """Fortran struct: tao_spin_map_struct"""

    def __init__(self, valid: bool | None = None, map1: SpinOrbitMap1Struct | None = None, axis_input: SpinAxisStruct | None = None, axis0: SpinAxisStruct | None = None, axis1: SpinAxisStruct | None = None, ix_ele: int | None = None, ix_ref: int | None = None, ix_uni: int | None = None, ix_branch: int | None = None, mat8: Sequence[Sequence[float]] | None = None) -> None: ...

    @property
    def valid(self) -> bool: ...

    @valid.setter
    def valid(self, arg: bool, /) -> None: ...

    @property
    def map1(self) -> SpinOrbitMap1Struct: ...

    @map1.setter
    def map1(self, arg: SpinOrbitMap1Struct, /) -> None: ...

    @property
    def axis_input(self) -> SpinAxisStruct:
        """Input axes."""

    @axis_input.setter
    def axis_input(self, arg: SpinAxisStruct, /) -> None: ...

    @property
    def axis0(self) -> SpinAxisStruct:
        """Initial axes."""

    @axis0.setter
    def axis0(self, arg: SpinAxisStruct, /) -> None: ...

    @property
    def axis1(self) -> SpinAxisStruct:
        """Final axes."""

    @axis1.setter
    def axis1(self, arg: SpinAxisStruct, /) -> None: ...

    @property
    def ix_ele(self) -> int: ...

    @ix_ele.setter
    def ix_ele(self, arg: int, /) -> None: ...

    @property
    def ix_ref(self) -> int: ...

    @ix_ref.setter
    def ix_ref(self, arg: int, /) -> None: ...

    @property
    def ix_uni(self) -> int: ...

    @ix_uni.setter
    def ix_uni(self, arg: int, /) -> None: ...

    @property
    def ix_branch(self) -> int: ...

    @ix_branch.setter
    def ix_branch(self, arg: int, /) -> None: ...

    @property
    def mat8(self) -> RealArray2D: ...

    @mat8.setter
    def mat8(self, arg: Sequence[Sequence[float]], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoSpinMapStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoSpinMapStruct: ...

    def __eq__(self, arg: TaoSpinMapStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoSpinPolarizationStruct:
    """Fortran struct: tao_spin_polarization_struct"""

    def __init__(self, tune: float | None = None, pol_limit_st: float | None = None, pol_limit_dk: float | None = None, pol_limit_dk_partial: Sequence[float] | None = None, pol_limit_dk_partial2: Sequence[float] | None = None, pol_rate_bks: float | None = None, depol_rate: float | None = None, depol_rate_partial: Sequence[float] | None = None, depol_rate_partial2: Sequence[float] | None = None, integral_bn: float | None = None, integral_bdn: float | None = None, integral_1ns: float | None = None, integral_dn2: float | None = None, valid: bool | None = None, q_1turn: SpinOrbitMap1Struct | None = None) -> None: ...

    @property
    def tune(self) -> float: ...

    @tune.setter
    def tune(self, arg: float, /) -> None: ...

    @property
    def pol_limit_st(self) -> float:
        """Polarization calculated using Sokolov-Ternov formula."""

    @pol_limit_st.setter
    def pol_limit_st(self, arg: float, /) -> None: ...

    @property
    def pol_limit_dk(self) -> float:
        """
        Equalibrium Polarization calculated via the Derbenev-Kondratenko-Mane formula.
        """

    @pol_limit_dk.setter
    def pol_limit_dk(self, arg: float, /) -> None: ...

    @property
    def pol_limit_dk_partial(self) -> RealArray1D:
        """Limit using only single mode to calc dn_dpz"""

    @pol_limit_dk_partial.setter
    def pol_limit_dk_partial(self, arg: Sequence[float], /) -> None: ...

    @property
    def pol_limit_dk_partial2(self) -> RealArray1D:
        """Limit using only single mode to calc dn_dpz"""

    @pol_limit_dk_partial2.setter
    def pol_limit_dk_partial2(self, arg: Sequence[float], /) -> None: ...

    @property
    def pol_rate_bks(self) -> float:
        """BKS Polarization rate (1/sec)."""

    @pol_rate_bks.setter
    def pol_rate_bks(self, arg: float, /) -> None: ...

    @property
    def depol_rate(self) -> float:
        """Depolarization rate (1/sec)."""

    @depol_rate.setter
    def depol_rate(self, arg: float, /) -> None: ...

    @property
    def depol_rate_partial(self) -> RealArray1D:
        """Depolarization rate (1/sec) using only single mode to calc dn_dpz."""

    @depol_rate_partial.setter
    def depol_rate_partial(self, arg: Sequence[float], /) -> None: ...

    @property
    def depol_rate_partial2(self) -> RealArray1D:
        """Depolarization rate (1/sec) using only two modes to calc dn_dpz."""

    @depol_rate_partial2.setter
    def depol_rate_partial2(self, arg: Sequence[float], /) -> None: ...

    @property
    def integral_bn(self) -> float:
        """Integral of g^3 * b_hat * n_0"""

    @integral_bn.setter
    def integral_bn(self, arg: float, /) -> None: ...

    @property
    def integral_bdn(self) -> float:
        """Integral of g^3 * b_hat * dn/ddelta"""

    @integral_bdn.setter
    def integral_bdn(self, arg: float, /) -> None: ...

    @property
    def integral_1ns(self) -> float:
        """Integral of g^3 (1 - 2(n * s_hat)/9)"""

    @integral_1ns.setter
    def integral_1ns(self, arg: float, /) -> None: ...

    @property
    def integral_dn2(self) -> float:
        """Integral of g^3 * 11 (dn/ddelta)^2 / 9"""

    @integral_dn2.setter
    def integral_dn2(self, arg: float, /) -> None: ...

    @property
    def valid(self) -> bool: ...

    @valid.setter
    def valid(self, arg: bool, /) -> None: ...

    @property
    def q_1turn(self) -> SpinOrbitMap1Struct:
        """Save results from spin_concat_linear_maps in tao_spin_polarization."""

    @q_1turn.setter
    def q_1turn(self, arg: SpinOrbitMap1Struct, /) -> None: ...

    @property
    def q_ele(self) -> SpinOrbitMap1StructAlloc1D:
        """Save results from spin_concat_linear_maps in tao_spin_polarization."""

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoSpinPolarizationStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoSpinPolarizationStruct: ...

    def __eq__(self, arg: TaoSpinPolarizationStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoStringArrayStruct:
    """Fortran struct: tao_string_array_struct"""

    def __init__(self, s: str | None = None) -> None: ...

    @property
    def s(self) -> str: ...

    @s.setter
    def s(self, arg: str, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoStringArrayStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoStringArrayStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoStringArrayStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoStringArrayStruct: ...

    def __eq__(self, arg: TaoStringArrayStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoSuperUniverseStruct:
    """Fortran struct: tao_super_universe_struct"""

    def __init__(self, global_: TaoGlobalStruct | None = None, init: TaoInitStruct | None = None, com: TaoCommonStruct | None = None, plot_page: TaoPlotPageStruct | None = None, key: Sequence[int] | None = None, building_wall: TaoBuildingWallStruct | None = None, wave: TaoWaveStruct | None = None, n_var_used: int | None = None, n_v1_var_used: int | None = None, initialized: bool | None = None) -> None: ...

    @property
    def init(self) -> TaoInitStruct:
        """Initialization parameters"""

    @init.setter
    def init(self, arg: TaoInitStruct, /) -> None: ...

    @property
    def com(self) -> TaoCommonStruct:
        """Non-initialization common parameters"""

    @com.setter
    def com(self, arg: TaoCommonStruct, /) -> None: ...

    @property
    def plot_page(self) -> TaoPlotPageStruct:
        """Defines the plot window."""

    @plot_page.setter
    def plot_page(self, arg: TaoPlotPageStruct, /) -> None: ...

    @property
    def v1_var(self) -> TaoV1VarStructAlloc1D:
        """The variable types"""

    @property
    def var(self) -> TaoVarStructAlloc1D:
        """array of all variables."""

    @property
    def u(self) -> TaoUniverseStructAlloc1D:
        """array of universes."""

    @property
    def key(self) -> IntAlloc1D: ...

    @key.setter
    def key(self, arg: Sequence[int], /) -> None: ...

    @property
    def building_wall(self) -> TaoBuildingWallStruct: ...

    @building_wall.setter
    def building_wall(self, arg: TaoBuildingWallStruct, /) -> None: ...

    @property
    def wave(self) -> TaoWaveStruct: ...

    @wave.setter
    def wave(self, arg: TaoWaveStruct, /) -> None: ...

    @property
    def n_var_used(self) -> int: ...

    @n_var_used.setter
    def n_var_used(self, arg: int, /) -> None: ...

    @property
    def n_v1_var_used(self) -> int: ...

    @n_v1_var_used.setter
    def n_v1_var_used(self, arg: int, /) -> None: ...

    @property
    def history(self) -> TaoCmdHistoryStructArray1D:
        """command history"""

    @property
    def initialized(self) -> bool:
        """Does tao_init() need to be called?"""

    @initialized.setter
    def initialized(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoSuperUniverseStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoSuperUniverseStruct: ...

    def __eq__(self, arg: TaoSuperUniverseStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoTitleStruct:
    """Fortran struct: tao_title_struct"""

    def __init__(self, string: str | None = None, x: float | None = None, y: float | None = None, units: str | None = None, justify: str | None = None, draw_it: bool | None = None) -> None: ...

    @property
    def string(self) -> str:
        """title character string."""

    @string.setter
    def string(self, arg: str, /) -> None: ...

    @property
    def x(self) -> float:
        """x, y rwt lower left corner"""

    @x.setter
    def x(self, arg: float, /) -> None: ...

    @property
    def y(self) -> float:
        """x, y rwt lower left corner"""

    @y.setter
    def y(self, arg: float, /) -> None: ...

    @property
    def units(self) -> str:
        """%BOX, POINTS, etc..."""

    @units.setter
    def units(self, arg: str, /) -> None: ...

    @property
    def justify(self) -> str:
        """Left, Center, or Right justification."""

    @justify.setter
    def justify(self, arg: str, /) -> None: ...

    @property
    def draw_it(self) -> bool:
        """draw the title?"""

    @draw_it.setter
    def draw_it(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoTitleStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoTitleStruct: ...

    def __eq__(self, arg: TaoTitleStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoTop10Struct:
    """Fortran struct: tao_top10_struct"""

    def __init__(self, name: str | None = None, value: float | None = None, index: int | None = None, valid: bool | None = None) -> None: ...

    @property
    def name(self) -> str:
        """name of contributor"""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def value(self) -> float:
        """contribution to the merit function"""

    @value.setter
    def value(self, arg: float, /) -> None: ...

    @property
    def index(self) -> int:
        """index of contributor."""

    @index.setter
    def index(self, arg: int, /) -> None: ...

    @property
    def valid(self) -> bool:
        """valid entry?"""

    @valid.setter
    def valid(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoTop10StructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoTop10StructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoTop10Struct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoTop10Struct: ...

    def __eq__(self, arg: TaoTop10Struct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoUniverseCalcStruct:
    """Fortran struct: tao_universe_calc_struct"""

    def __init__(self, srdt_for_data: int | None = None, rad_int_for_data: bool | None = None, rad_int_for_plotting: bool | None = None, chrom_for_data: bool | None = None, chrom_for_plotting: bool | None = None, lat_sigma_for_data: bool | None = None, lat_sigma_for_plotting: bool | None = None, dynamic_aperture: bool | None = None, one_turn_map: bool | None = None, lattice: bool | None = None, twiss: bool | None = None, track: bool | None = None, spin_matrices: bool | None = None) -> None: ...

    @property
    def srdt_for_data(self) -> int:
        """0 = false, 1 = 1st order, 2 = 1st & 2nd order"""

    @srdt_for_data.setter
    def srdt_for_data(self, arg: int, /) -> None: ...

    @property
    def rad_int_for_data(self) -> bool:
        """Do the radiation integrals need to be computed for"""

    @rad_int_for_data.setter
    def rad_int_for_data(self, arg: bool, /) -> None: ...

    @property
    def rad_int_for_plotting(self) -> bool:
        """data or plotting?"""

    @rad_int_for_plotting.setter
    def rad_int_for_plotting(self, arg: bool, /) -> None: ...

    @property
    def chrom_for_data(self) -> bool:
        """Does the chromaticity need to be computed for"""

    @chrom_for_data.setter
    def chrom_for_data(self, arg: bool, /) -> None: ...

    @property
    def chrom_for_plotting(self) -> bool:
        """data or plotting?"""

    @chrom_for_plotting.setter
    def chrom_for_plotting(self, arg: bool, /) -> None: ...

    @property
    def lat_sigma_for_data(self) -> bool:
        """Do the beam sigmas need to be computed for"""

    @lat_sigma_for_data.setter
    def lat_sigma_for_data(self, arg: bool, /) -> None: ...

    @property
    def lat_sigma_for_plotting(self) -> bool:
        """data or plotting?"""

    @lat_sigma_for_plotting.setter
    def lat_sigma_for_plotting(self, arg: bool, /) -> None: ...

    @property
    def dynamic_aperture(self) -> bool:
        """Do the dynamic_aperture calc?"""

    @dynamic_aperture.setter
    def dynamic_aperture(self, arg: bool, /) -> None: ...

    @property
    def one_turn_map(self) -> bool:
        """Compute the one turn map?"""

    @one_turn_map.setter
    def one_turn_map(self, arg: bool, /) -> None: ...

    @property
    def lattice(self) -> bool:
        """Used to indicate which lattices need tracking done."""

    @lattice.setter
    def lattice(self, arg: bool, /) -> None: ...

    @property
    def twiss(self) -> bool:
        """calc linear transfer matrix?"""

    @twiss.setter
    def twiss(self, arg: bool, /) -> None: ...

    @property
    def track(self) -> bool:
        """tracking needs to be done?"""

    @track.setter
    def track(self, arg: bool, /) -> None: ...

    @property
    def spin_matrices(self) -> bool:
        """Calculate G and D spin matrices?"""

    @spin_matrices.setter
    def spin_matrices(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoUniverseCalcStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoUniverseCalcStruct: ...

    def __eq__(self, arg: TaoUniverseCalcStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoUniversePointerStruct:
    """Fortran struct: tao_universe_pointer_struct"""

    def __init__(self, u: TaoUniverseStruct | None = None) -> None: ...

    @property
    def u(self) -> TaoUniverseStruct | None: ...

    @u.setter
    def u(self, arg: TaoUniverseStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoUniversePointerStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoUniversePointerStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoUniversePointerStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoUniversePointerStruct: ...

    def __eq__(self, arg: TaoUniversePointerStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoUniverseStruct:
    """Fortran struct: tao_universe_struct"""

    def __init__(self, model: TaoLatticeStruct | None = None, design: TaoLatticeStruct | None = None, base: TaoLatticeStruct | None = None, beam: TaoBeamUniStruct | None = None, dynamic_aperture: TaoDynamicApertureStruct | None = None, ping_scale: TaoPingScaleStruct | None = None, scratch_lat: LatStruct | None = None, calc: TaoUniverseCalcStruct | None = None, ele_order: LatEleOrderStruct | None = None, spin_map: TaoSpinMapStruct | None = None, dModel_dVar: Sequence[Sequence[float]] | None = None, ix_uni: int | None = None, n_d2_data_used: int | None = None, n_data_used: int | None = None, is_on: bool | None = None, design_same_as_previous: bool | None = None, picked_uni: bool | None = None) -> None: ...

    @property
    def model(self) -> TaoLatticeStruct | None: ...

    @model.setter
    def model(self, arg: TaoLatticeStruct, /) -> None: ...

    @property
    def design(self) -> TaoLatticeStruct | None: ...

    @design.setter
    def design(self, arg: TaoLatticeStruct, /) -> None: ...

    @property
    def base(self) -> TaoLatticeStruct | None: ...

    @base.setter
    def base(self, arg: TaoLatticeStruct, /) -> None: ...

    @property
    def beam(self) -> TaoBeamUniStruct: ...

    @beam.setter
    def beam(self, arg: TaoBeamUniStruct, /) -> None: ...

    @property
    def dynamic_aperture(self) -> TaoDynamicApertureStruct: ...

    @dynamic_aperture.setter
    def dynamic_aperture(self, arg: TaoDynamicApertureStruct, /) -> None: ...

    @property
    def model_branch(self) -> TaoModelBranchStructArray1D:
        """model specific information"""

    @property
    def d2_data(self) -> TaoD2DataStructAlloc1D:
        """The data types"""

    @property
    def data(self) -> TaoDataStructAlloc1D:
        """Array of all data."""

    @property
    def ping_scale(self) -> TaoPingScaleStruct: ...

    @ping_scale.setter
    def ping_scale(self, arg: TaoPingScaleStruct, /) -> None: ...

    @property
    def scratch_lat(self) -> LatStruct:
        """Scratch area."""

    @scratch_lat.setter
    def scratch_lat(self, arg: LatStruct, /) -> None: ...

    @property
    def calc(self) -> TaoUniverseCalcStruct:
        """What needs to be calculated?"""

    @calc.setter
    def calc(self, arg: TaoUniverseCalcStruct, /) -> None: ...

    @property
    def ele_order(self) -> LatEleOrderStruct:
        """Order of elements with same name."""

    @ele_order.setter
    def ele_order(self, arg: LatEleOrderStruct, /) -> None: ...

    @property
    def spin_map(self) -> TaoSpinMapStruct: ...

    @spin_map.setter
    def spin_map(self, arg: TaoSpinMapStruct, /) -> None: ...

    @property
    def dModel_dVar(self) -> RealArray2D:
        """Derivative matrix."""

    @dModel_dVar.setter
    def dModel_dVar(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @property
    def ix_uni(self) -> int:
        """Universe index."""

    @ix_uni.setter
    def ix_uni(self, arg: int, /) -> None: ...

    @property
    def n_d2_data_used(self) -> int:
        """Number of used %d2_data(:) components."""

    @n_d2_data_used.setter
    def n_d2_data_used(self, arg: int, /) -> None: ...

    @property
    def n_data_used(self) -> int:
        """Number of used %data(:) components."""

    @n_data_used.setter
    def n_data_used(self, arg: int, /) -> None: ...

    @property
    def is_on(self) -> bool:
        """universe turned on"""

    @is_on.setter
    def is_on(self, arg: bool, /) -> None: ...

    @property
    def design_same_as_previous(self) -> bool:
        """Design lat same as the previous uni?"""

    @design_same_as_previous.setter
    def design_same_as_previous(self, arg: bool, /) -> None: ...

    @property
    def picked_uni(self) -> bool:
        """Scratch logical."""

    @picked_uni.setter
    def picked_uni(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoUniverseStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoUniverseStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoUniverseStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoUniverseStruct: ...

    def __eq__(self, arg: TaoUniverseStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoV1VarArrayStruct:
    """Fortran struct: tao_v1_var_array_struct"""

    def __init__(self, v1: TaoV1VarStruct | None = None) -> None: ...

    @property
    def v1(self) -> TaoV1VarStruct | None: ...

    @v1.setter
    def v1(self, arg: TaoV1VarStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoV1VarArrayStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoV1VarArrayStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoV1VarArrayStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoV1VarArrayStruct: ...

    def __eq__(self, arg: TaoV1VarArrayStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoV1VarStruct:
    """Fortran struct: tao_v1_var_struct"""

    def __init__(self, name: str | None = None, ix_v1_var: int | None = None) -> None: ...

    @property
    def name(self) -> str:
        """V1 variable name. Eg: 'quad_k1'."""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def ix_v1_var(self) -> int:
        """Index to s%v1_var(:) array"""

    @ix_v1_var.setter
    def ix_v1_var(self, arg: int, /) -> None: ...

    @property
    def v(self) -> TaoVarStructArray1D:
        """Pointer to the appropriate section in s%var."""

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoV1VarStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoV1VarStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoV1VarStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoV1VarStruct: ...

    def __eq__(self, arg: TaoV1VarStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoVarArrayStruct:
    """Fortran struct: tao_var_array_struct"""

    def __init__(self, v: TaoVarStruct | None = None) -> None: ...

    @property
    def v(self) -> TaoVarStruct | None: ...

    @v.setter
    def v(self, arg: TaoVarStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoVarArrayStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoVarArrayStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoVarArrayStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoVarArrayStruct: ...

    def __eq__(self, arg: TaoVarArrayStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoVarSlaveStruct:
    """Fortran struct: tao_var_slave_struct"""

    def __init__(self, ix_uni: int | None = None, ix_branch: int | None = None, ix_ele: int | None = None, model_value: float | None = None, base_value: float | None = None) -> None: ...

    @property
    def ix_uni(self) -> int:
        """universe index."""

    @ix_uni.setter
    def ix_uni(self, arg: int, /) -> None: ...

    @property
    def ix_branch(self) -> int: ...

    @ix_branch.setter
    def ix_branch(self, arg: int, /) -> None: ...

    @property
    def ix_ele(self) -> int:
        """Index of element in the u%lattice%ele(:) array."""

    @ix_ele.setter
    def ix_ele(self, arg: int, /) -> None: ...

    @property
    def model_value(self) -> float | None:
        """Pointer to the variable in the model lat."""

    @model_value.setter
    def model_value(self, arg: float, /) -> None: ...

    @property
    def base_value(self) -> float | None:
        """Pointer to the variable in the base lat."""

    @base_value.setter
    def base_value(self, arg: float, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoVarSlaveStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoVarSlaveStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoVarSlaveStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoVarSlaveStruct: ...

    def __eq__(self, arg: TaoVarSlaveStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoVarStruct:
    """Fortran struct: tao_var_struct"""

    def __init__(self, ele_name: str | None = None, attrib_name: str | None = None, id: str | None = None, ix_v1: int | None = None, ix_var: int | None = None, ix_dvar: int | None = None, ix_attrib: int | None = None, ix_key_table: int | None = None, model_value: float | None = None, base_value: float | None = None, design_value: float | None = None, scratch_value: float | None = None, old_value: float | None = None, meas_value: float | None = None, ref_value: float | None = None, correction_value: float | None = None, high_lim: float | None = None, low_lim: float | None = None, step: float | None = None, weight: float | None = None, delta_merit: float | None = None, merit: float | None = None, dMerit_dVar: float | None = None, key_val0: float | None = None, key_delta: float | None = None, s: float | None = None, extend_val: float | None = None, merit_type: str | None = None, exists: bool | None = None, good_var: bool | None = None, good_user: bool | None = None, good_opt: bool | None = None, good_plot: bool | None = None, useit_opt: bool | None = None, useit_plot: bool | None = None, key_bound: bool | None = None, v1: TaoV1VarStruct | None = None) -> None: ...

    @property
    def ele_name(self) -> str:
        """Associated lattice element name."""

    @ele_name.setter
    def ele_name(self, arg: str, /) -> None: ...

    @property
    def attrib_name(self) -> str:
        """Name of the attribute to vary."""

    @attrib_name.setter
    def attrib_name(self, arg: str, /) -> None: ...

    @property
    def id(self) -> str:
        """Used by Tao extension code. Not used by Tao directly."""

    @id.setter
    def id(self, arg: str, /) -> None: ...

    @property
    def slave(self) -> TaoVarSlaveStructAlloc1D: ...

    @property
    def ix_v1(self) -> int:
        """Index of this var in the s%v1_var(i)%v(:) array."""

    @ix_v1.setter
    def ix_v1(self, arg: int, /) -> None: ...

    @property
    def ix_var(self) -> int:
        """Index number of this var in the s%var(:) array."""

    @ix_var.setter
    def ix_var(self, arg: int, /) -> None: ...

    @property
    def ix_dvar(self) -> int:
        """Column in the dData_dVar derivative matrix."""

    @ix_dvar.setter
    def ix_dvar(self, arg: int, /) -> None: ...

    @property
    def ix_attrib(self) -> int:
        """Index in ele%value(:) array if appropriate."""

    @ix_attrib.setter
    def ix_attrib(self, arg: int, /) -> None: ...

    @property
    def ix_key_table(self) -> int:
        """Has a key binding?"""

    @ix_key_table.setter
    def ix_key_table(self, arg: int, /) -> None: ...

    @property
    def model_value(self) -> float | None:
        """Model value."""

    @model_value.setter
    def model_value(self, arg: float, /) -> None: ...

    @property
    def base_value(self) -> float | None:
        """Base value."""

    @base_value.setter
    def base_value(self, arg: float, /) -> None: ...

    @property
    def design_value(self) -> float:
        """Design value from the design lattice."""

    @design_value.setter
    def design_value(self, arg: float, /) -> None: ...

    @property
    def scratch_value(self) -> float:
        """Scratch space used by Tao."""

    @scratch_value.setter
    def scratch_value(self, arg: float, /) -> None: ...

    @property
    def old_value(self) -> float:
        """Scratch space used by Tao."""

    @old_value.setter
    def old_value(self, arg: float, /) -> None: ...

    @property
    def meas_value(self) -> float:
        """The value when the data measurement was taken."""

    @meas_value.setter
    def meas_value(self, arg: float, /) -> None: ...

    @property
    def ref_value(self) -> float:
        """Value when the reference measurement was taken."""

    @ref_value.setter
    def ref_value(self, arg: float, /) -> None: ...

    @property
    def correction_value(self) -> float:
        """Value determined by a fit to correct the lattice."""

    @correction_value.setter
    def correction_value(self, arg: float, /) -> None: ...

    @property
    def high_lim(self) -> float:
        """High limit for the model_value."""

    @high_lim.setter
    def high_lim(self, arg: float, /) -> None: ...

    @property
    def low_lim(self) -> float:
        """Low limit for the model_value."""

    @low_lim.setter
    def low_lim(self, arg: float, /) -> None: ...

    @property
    def step(self) -> float:
        """Sets what is a small step for varying this var."""

    @step.setter
    def step(self, arg: float, /) -> None: ...

    @property
    def weight(self) -> float:
        """Weight for the merit function term."""

    @weight.setter
    def weight(self, arg: float, /) -> None: ...

    @property
    def delta_merit(self) -> float:
        """Diff used to calculate the merit function term."""

    @delta_merit.setter
    def delta_merit(self, arg: float, /) -> None: ...

    @property
    def merit(self) -> float:
        """merit_term = weight * delta^2."""

    @merit.setter
    def merit(self, arg: float, /) -> None: ...

    @property
    def dMerit_dVar(self) -> float:
        """Merit derivative."""

    @dMerit_dVar.setter
    def dMerit_dVar(self, arg: float, /) -> None: ...

    @property
    def key_val0(self) -> float:
        """Key base value"""

    @key_val0.setter
    def key_val0(self, arg: float, /) -> None: ...

    @property
    def key_delta(self) -> float:
        """Change in value when a key is pressed."""

    @key_delta.setter
    def key_delta(self, arg: float, /) -> None: ...

    @property
    def s(self) -> float:
        """longitudinal position of ele."""

    @s.setter
    def s(self, arg: float, /) -> None: ...

    @property
    def extend_val(self) -> float:
        """For extension code. Not used by Tao."""

    @extend_val.setter
    def extend_val(self, arg: float, /) -> None: ...

    @property
    def merit_type(self) -> str:
        """'target' or 'limit'"""

    @merit_type.setter
    def merit_type(self, arg: str, /) -> None: ...

    @property
    def exists(self) -> bool:
        """See above"""

    @exists.setter
    def exists(self, arg: bool, /) -> None: ...

    @property
    def good_var(self) -> bool:
        """See above"""

    @good_var.setter
    def good_var(self, arg: bool, /) -> None: ...

    @property
    def good_user(self) -> bool:
        """See above"""

    @good_user.setter
    def good_user(self, arg: bool, /) -> None: ...

    @property
    def good_opt(self) -> bool:
        """See above"""

    @good_opt.setter
    def good_opt(self, arg: bool, /) -> None: ...

    @property
    def good_plot(self) -> bool:
        """See above"""

    @good_plot.setter
    def good_plot(self, arg: bool, /) -> None: ...

    @property
    def useit_opt(self) -> bool:
        """See above"""

    @useit_opt.setter
    def useit_opt(self, arg: bool, /) -> None: ...

    @property
    def useit_plot(self) -> bool:
        """See above"""

    @useit_plot.setter
    def useit_plot(self, arg: bool, /) -> None: ...

    @property
    def key_bound(self) -> bool:
        """Variable bound to keyboard key?"""

    @key_bound.setter
    def key_bound(self, arg: bool, /) -> None: ...

    @property
    def v1(self) -> TaoV1VarStruct | None:
        """Pointer to the parent."""

    @v1.setter
    def v1(self, arg: TaoV1VarStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoVarStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoVarStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoVarStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoVarStruct: ...

    def __eq__(self, arg: TaoVarStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoWaveKickPtStruct:
    """Fortran struct: tao_wave_kick_pt_struct"""

    def __init__(self, phi_s: float | None = None, phi_r: float | None = None, phi: float | None = None, amp: float | None = None, s: float | None = None, ix_dat_before_kick: int | None = None, ele: EleStruct | None = None) -> None: ...

    @property
    def phi_s(self) -> float: ...

    @phi_s.setter
    def phi_s(self, arg: float, /) -> None: ...

    @property
    def phi_r(self) -> float: ...

    @phi_r.setter
    def phi_r(self, arg: float, /) -> None: ...

    @property
    def phi(self) -> float: ...

    @phi.setter
    def phi(self, arg: float, /) -> None: ...

    @property
    def amp(self) -> float: ...

    @amp.setter
    def amp(self, arg: float, /) -> None: ...

    @property
    def s(self) -> float:
        """s-position of kick"""

    @s.setter
    def s(self, arg: float, /) -> None: ...

    @property
    def ix_dat_before_kick(self) -> int:
        """Index of datum in data array just before the kick."""

    @ix_dat_before_kick.setter
    def ix_dat_before_kick(self, arg: int, /) -> None: ...

    @property
    def ele(self) -> EleStruct | None:
        """lattice element at position of kick."""

    @ele.setter
    def ele(self, arg: EleStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaoWaveKickPtStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaoWaveKickPtStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoWaveKickPtStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoWaveKickPtStruct: ...

    def __eq__(self, arg: TaoWaveKickPtStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaoWaveStruct:
    """Fortran struct: tao_wave_struct"""

    def __init__(self, data_type: str | None = None, rms_rel_a: float | None = None, rms_rel_b: float | None = None, rms_rel_as: float | None = None, rms_rel_bs: float | None = None, rms_rel_ar: float | None = None, rms_rel_br: float | None = None, rms_rel_k: float | None = None, rms_rel_ks: float | None = None, rms_rel_kr: float | None = None, rms_phi: float | None = None, rms_phi_s: float | None = None, rms_phi_r: float | None = None, amp_ba_s: float | None = None, amp_ba_r: float | None = None, chi_a: float | None = None, chi_c: float | None = None, chi_ba: float | None = None, amp_a: Sequence[float] | None = None, amp_b: Sequence[float] | None = None, amp_ba: Sequence[float] | None = None, coef_a: Sequence[float] | None = None, coef_b: Sequence[float] | None = None, coef_ba: Sequence[float] | None = None, n_func: int | None = None, ix_a1: int | None = None, ix_a2: int | None = None, ix_b1: int | None = None, ix_b2: int | None = None, i_a1: int | None = None, i_a2: int | None = None, i_b1: int | None = None, i_b2: int | None = None, n_a: int | None = None, n_b: int | None = None, i_curve_wrap_pt: int | None = None, ix_data: Sequence[int] | None = None, n_kick: int | None = None, base_graph: TaoGraphStruct | None = None, region: TaoPlotRegionStruct | None = None, d1_dat: TaoD1DataStruct | None = None) -> None: ...

    @property
    def data_type(self) -> str: ...

    @data_type.setter
    def data_type(self, arg: str, /) -> None: ...

    @property
    def rms_rel_a(self) -> float: ...

    @rms_rel_a.setter
    def rms_rel_a(self, arg: float, /) -> None: ...

    @property
    def rms_rel_b(self) -> float: ...

    @rms_rel_b.setter
    def rms_rel_b(self, arg: float, /) -> None: ...

    @property
    def rms_rel_as(self) -> float: ...

    @rms_rel_as.setter
    def rms_rel_as(self, arg: float, /) -> None: ...

    @property
    def rms_rel_bs(self) -> float: ...

    @rms_rel_bs.setter
    def rms_rel_bs(self, arg: float, /) -> None: ...

    @property
    def rms_rel_ar(self) -> float: ...

    @rms_rel_ar.setter
    def rms_rel_ar(self, arg: float, /) -> None: ...

    @property
    def rms_rel_br(self) -> float: ...

    @rms_rel_br.setter
    def rms_rel_br(self, arg: float, /) -> None: ...

    @property
    def rms_rel_k(self) -> float: ...

    @rms_rel_k.setter
    def rms_rel_k(self, arg: float, /) -> None: ...

    @property
    def rms_rel_ks(self) -> float: ...

    @rms_rel_ks.setter
    def rms_rel_ks(self, arg: float, /) -> None: ...

    @property
    def rms_rel_kr(self) -> float: ...

    @rms_rel_kr.setter
    def rms_rel_kr(self, arg: float, /) -> None: ...

    @property
    def rms_phi(self) -> float: ...

    @rms_phi.setter
    def rms_phi(self, arg: float, /) -> None: ...

    @property
    def rms_phi_s(self) -> float: ...

    @rms_phi_s.setter
    def rms_phi_s(self, arg: float, /) -> None: ...

    @property
    def rms_phi_r(self) -> float: ...

    @rms_phi_r.setter
    def rms_phi_r(self, arg: float, /) -> None: ...

    @property
    def amp_ba_s(self) -> float: ...

    @amp_ba_s.setter
    def amp_ba_s(self, arg: float, /) -> None: ...

    @property
    def amp_ba_r(self) -> float: ...

    @amp_ba_r.setter
    def amp_ba_r(self, arg: float, /) -> None: ...

    @property
    def chi_a(self) -> float: ...

    @chi_a.setter
    def chi_a(self, arg: float, /) -> None: ...

    @property
    def chi_c(self) -> float: ...

    @chi_c.setter
    def chi_c(self, arg: float, /) -> None: ...

    @property
    def chi_ba(self) -> float: ...

    @chi_ba.setter
    def chi_ba(self, arg: float, /) -> None: ...

    @property
    def amp_a(self) -> RealArray1D: ...

    @amp_a.setter
    def amp_a(self, arg: Sequence[float], /) -> None: ...

    @property
    def amp_b(self) -> RealArray1D: ...

    @amp_b.setter
    def amp_b(self, arg: Sequence[float], /) -> None: ...

    @property
    def amp_ba(self) -> RealArray1D: ...

    @amp_ba.setter
    def amp_ba(self, arg: Sequence[float], /) -> None: ...

    @property
    def coef_a(self) -> RealArray1D: ...

    @coef_a.setter
    def coef_a(self, arg: Sequence[float], /) -> None: ...

    @property
    def coef_b(self) -> RealArray1D: ...

    @coef_b.setter
    def coef_b(self, arg: Sequence[float], /) -> None: ...

    @property
    def coef_ba(self) -> RealArray1D: ...

    @coef_ba.setter
    def coef_ba(self, arg: Sequence[float], /) -> None: ...

    @property
    def n_func(self) -> int:
        """Number of functions used in the fit."""

    @n_func.setter
    def n_func(self, arg: int, /) -> None: ...

    @property
    def ix_a1(self) -> int: ...

    @ix_a1.setter
    def ix_a1(self, arg: int, /) -> None: ...

    @property
    def ix_a2(self) -> int: ...

    @ix_a2.setter
    def ix_a2(self, arg: int, /) -> None: ...

    @property
    def ix_b1(self) -> int: ...

    @ix_b1.setter
    def ix_b1(self, arg: int, /) -> None: ...

    @property
    def ix_b2(self) -> int: ...

    @ix_b2.setter
    def ix_b2(self, arg: int, /) -> None: ...

    @property
    def i_a1(self) -> int: ...

    @i_a1.setter
    def i_a1(self, arg: int, /) -> None: ...

    @property
    def i_a2(self) -> int: ...

    @i_a2.setter
    def i_a2(self, arg: int, /) -> None: ...

    @property
    def i_b1(self) -> int: ...

    @i_b1.setter
    def i_b1(self, arg: int, /) -> None: ...

    @property
    def i_b2(self) -> int: ...

    @i_b2.setter
    def i_b2(self, arg: int, /) -> None: ...

    @property
    def n_a(self) -> int: ...

    @n_a.setter
    def n_a(self, arg: int, /) -> None: ...

    @property
    def n_b(self) -> int: ...

    @n_b.setter
    def n_b(self, arg: int, /) -> None: ...

    @property
    def i_curve_wrap_pt(self) -> int:
        """Index of last point before wrap in curve array."""

    @i_curve_wrap_pt.setter
    def i_curve_wrap_pt(self, arg: int, /) -> None: ...

    @property
    def ix_data(self) -> IntAlloc1D:
        """Translates from plot point to datum index"""

    @ix_data.setter
    def ix_data(self, arg: Sequence[int], /) -> None: ...

    @property
    def n_kick(self) -> int: ...

    @n_kick.setter
    def n_kick(self, arg: int, /) -> None: ...

    @property
    def kick(self) -> TaoWaveKickPtStructAlloc1D: ...

    @property
    def base_graph(self) -> TaoGraphStruct:
        """Graph before curves extended to 1.5 periods."""

    @base_graph.setter
    def base_graph(self, arg: TaoGraphStruct, /) -> None: ...

    @property
    def region(self) -> TaoPlotRegionStruct | None:
        """Where the wave plot is"""

    @region.setter
    def region(self, arg: TaoPlotRegionStruct, /) -> None: ...

    @property
    def d1_dat(self) -> TaoD1DataStruct | None:
        """D1 data for analysis"""

    @d1_dat.setter
    def d1_dat(self, arg: TaoD1DataStruct, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaoWaveStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaoWaveStruct: ...

    def __eq__(self, arg: TaoWaveStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TargetPointStruct:
    """Fortran struct: target_point_struct"""

    def __init__(self, r: Sequence[float] | None = None) -> None: ...

    @property
    def r(self) -> RealArray1D:
        """(x, y, z)"""

    @r.setter
    def r(self, arg: Sequence[float], /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TargetPointStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TargetPointStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TargetPointStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TargetPointStruct: ...

    def __eq__(self, arg: TargetPointStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaylorStruct:
    """Fortran struct: taylor_struct"""

    def __init__(self, ref: float | None = None) -> None: ...

    @property
    def ref(self) -> float: ...

    @ref.setter
    def ref(self, arg: float, /) -> None: ...

    @property
    def term(self) -> TaylorTermStructArray1D: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaylorStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaylorStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaylorStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaylorStruct: ...

    def __eq__(self, arg: TaylorStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TaylorTermStruct:
    """Fortran struct: taylor_term_struct"""

    def __init__(self, coef: float | None = None, expn: Sequence[int] | None = None) -> None: ...

    @property
    def coef(self) -> float: ...

    @coef.setter
    def coef(self, arg: float, /) -> None: ...

    @property
    def expn(self) -> IntArray1D: ...

    @expn.setter
    def expn(self, arg: Sequence[int], /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TaylorTermStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TaylorTermStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TaylorTermStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TaylorTermStruct: ...

    def __eq__(self, arg: TaylorTermStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TestSubStruct:
    """Fortran struct: test_sub_struct"""

    def __init__(self, sr: TestSubSubStruct | None = None) -> None: ...

    @property
    def sr(self) -> TestSubSubStruct: ...

    @sr.setter
    def sr(self, arg: TestSubSubStruct, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TestSubStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TestSubStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TestSubStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TestSubStruct: ...

    def __eq__(self, arg: TestSubStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TestSubSubStruct:
    """Fortran struct: test_sub_sub_struct"""

    def __init__(self, aaa: int | None = None, bbb: int | None = None, file: str | None = None, t_ref: float | None = None, freq_spread: float | None = None) -> None: ...

    @property
    def aaa(self) -> int: ...

    @aaa.setter
    def aaa(self, arg: int, /) -> None: ...

    @property
    def bbb(self) -> int: ...

    @bbb.setter
    def bbb(self, arg: int, /) -> None: ...

    @property
    def file(self) -> str: ...

    @file.setter
    def file(self, arg: str, /) -> None: ...

    @property
    def t_ref(self) -> float:
        """
        time reference value for computing the wake amplitude. This is used to prevent value overflow with long trains.
        """

    @t_ref.setter
    def t_ref(self, arg: float, /) -> None: ...

    @property
    def freq_spread(self) -> float:
        """Random frequency spread of long range modes."""

    @freq_spread.setter
    def freq_spread(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TestSubSubStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TestSubSubStruct: ...

    def __eq__(self, arg: TestSubSubStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TrackPointStruct:
    """Fortran struct: track_point_struct"""

    def __init__(self, s_lab: float | None = None, s_body: float | None = None, orb: CoordStruct | None = None, field: EmFieldStruct | None = None, strong_beam: StrongBeamStruct | None = None, vec0: Sequence[float] | None = None, mat6: Sequence[Sequence[float]] | None = None) -> None: ...

    @property
    def s_lab(self) -> float:
        """Longitudinal lab coord with respect to the upstream end."""

    @s_lab.setter
    def s_lab(self, arg: float, /) -> None: ...

    @property
    def s_body(self) -> float:
        """Longitudinal body coord with respect to the entrance end."""

    @s_body.setter
    def s_body(self, arg: float, /) -> None: ...

    @property
    def orb(self) -> CoordStruct:
        """Particle position in lab coords."""

    @orb.setter
    def orb(self, arg: CoordStruct, /) -> None: ...

    @property
    def field(self) -> EmFieldStruct:
        """E&M fields in lab coordinates."""

    @field.setter
    def field(self, arg: EmFieldStruct, /) -> None: ...

    @property
    def strong_beam(self) -> StrongBeamStruct:
        """Strong beam info for beambeam element."""

    @strong_beam.setter
    def strong_beam(self, arg: StrongBeamStruct, /) -> None: ...

    @property
    def vec0(self) -> RealArray1D:
        """0th order part of xfer map from the beginning."""

    @vec0.setter
    def vec0(self, arg: Sequence[float], /) -> None: ...

    @property
    def mat6(self) -> RealArray2D:
        """1st order part of xfer map (transfer matrix)."""

    @mat6.setter
    def mat6(self, arg: Sequence[Sequence[float]], /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> TrackPointStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> TrackPointStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TrackPointStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TrackPointStruct: ...

    def __eq__(self, arg: TrackPointStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TrackStruct:
    """Fortran struct: track_struct"""

    def __init__(self, ds_save: float | None = None, n_pt: int | None = None, n_bad: int | None = None, n_ok: int | None = None) -> None: ...

    @property
    def pt(self) -> TrackPointStructAlloc1D:
        """Array of track points indexed from 0."""

    @property
    def ds_save(self) -> float:
        """Min distance between points. Not positive => Save at all points."""

    @ds_save.setter
    def ds_save(self, arg: float, /) -> None: ...

    @property
    def n_pt(self) -> int:
        """
        Track upper bound for %pt(0:) array. n_bad and n_ok are used by adaptive trackers to record the number of times the step length had to be shortened.
        """

    @n_pt.setter
    def n_pt(self, arg: int, /) -> None: ...

    @property
    def n_bad(self) -> int:
        """Number of 'bad' steps where the step length was shortened."""

    @n_bad.setter
    def n_bad(self, arg: int, /) -> None: ...

    @property
    def n_ok(self) -> int:
        """Number of 'good' steps where the step length was not shortened."""

    @n_ok.setter
    def n_ok(self, arg: int, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TrackStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TrackStruct: ...

    def __eq__(self, arg: TrackStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TricubicCmplxCoefStruct:
    """Fortran struct: tricubic_cmplx_coef_struct"""

    def __init__(self, coef: Sequence[Sequence[Sequence[complex]]] | None = None, i_box: Sequence[int] | None = None) -> None: ...

    @property
    def coef(self) -> ComplexArray3D:
        """Coefs"""

    @coef.setter
    def coef(self, arg: Sequence[Sequence[Sequence[complex]]], /) -> None: ...

    @property
    def i_box(self) -> IntArray1D:
        """index at lower box corner."""

    @i_box.setter
    def i_box(self, arg: Sequence[int], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TricubicCmplxCoefStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TricubicCmplxCoefStruct: ...

    def __eq__(self, arg: TricubicCmplxCoefStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TricubicCoefStruct:
    """Fortran struct: tricubic_coef_struct"""

    def __init__(self, coef: Sequence[Sequence[Sequence[float]]] | None = None, i_box: Sequence[int] | None = None) -> None: ...

    @property
    def coef(self) -> RealArray3D:
        """Coefs"""

    @coef.setter
    def coef(self, arg: Sequence[Sequence[Sequence[float]]], /) -> None: ...

    @property
    def i_box(self) -> IntArray1D:
        """index at lower box corner."""

    @i_box.setter
    def i_box(self, arg: Sequence[int], /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TricubicCoefStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TricubicCoefStruct: ...

    def __eq__(self, arg: TricubicCoefStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class TwissStruct:
    """Fortran struct: twiss_struct"""

    def __init__(self, beta: float | None = None, alpha: float | None = None, gamma: float | None = None, phi: float | None = None, eta: float | None = None, etap: float | None = None, deta_ds: float | None = None, sigma: float | None = None, sigma_p: float | None = None, emit: float | None = None, norm_emit: float | None = None, chrom: float | None = None, dbeta_dpz: float | None = None, dalpha_dpz: float | None = None, deta_dpz: float | None = None, detap_dpz: float | None = None) -> None: ...

    @property
    def beta(self) -> float: ...

    @beta.setter
    def beta(self, arg: float, /) -> None: ...

    @property
    def alpha(self) -> float: ...

    @alpha.setter
    def alpha(self, arg: float, /) -> None: ...

    @property
    def gamma(self) -> float: ...

    @gamma.setter
    def gamma(self, arg: float, /) -> None: ...

    @property
    def phi(self) -> float: ...

    @phi.setter
    def phi(self, arg: float, /) -> None: ...

    @property
    def eta(self) -> float: ...

    @eta.setter
    def eta(self, arg: float, /) -> None: ...

    @property
    def etap(self) -> float: ...

    @etap.setter
    def etap(self, arg: float, /) -> None: ...

    @property
    def deta_ds(self) -> float: ...

    @deta_ds.setter
    def deta_ds(self, arg: float, /) -> None: ...

    @property
    def sigma(self) -> float: ...

    @sigma.setter
    def sigma(self, arg: float, /) -> None: ...

    @property
    def sigma_p(self) -> float: ...

    @sigma_p.setter
    def sigma_p(self, arg: float, /) -> None: ...

    @property
    def emit(self) -> float: ...

    @emit.setter
    def emit(self, arg: float, /) -> None: ...

    @property
    def norm_emit(self) -> float: ...

    @norm_emit.setter
    def norm_emit(self, arg: float, /) -> None: ...

    @property
    def chrom(self) -> float: ...

    @chrom.setter
    def chrom(self, arg: float, /) -> None: ...

    @property
    def dbeta_dpz(self) -> float: ...

    @dbeta_dpz.setter
    def dbeta_dpz(self, arg: float, /) -> None: ...

    @property
    def dalpha_dpz(self) -> float: ...

    @dalpha_dpz.setter
    def dalpha_dpz(self, arg: float, /) -> None: ...

    @property
    def deta_dpz(self) -> float: ...

    @deta_dpz.setter
    def deta_dpz(self, arg: float, /) -> None: ...

    @property
    def detap_dpz(self) -> float: ...

    @detap_dpz.setter
    def detap_dpz(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> TwissStruct: ...

    def __deepcopy__(self, arg: dict, /) -> TwissStruct: ...

    def __eq__(self, arg: TwissStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class VarLengthStringStruct:
    """Fortran struct: var_length_string_struct"""

    def __init__(self, str: str | None = None) -> None: ...

    @property
    def str(self) -> str: ...

    @str.setter
    def str(self, arg: str, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> VarLengthStringStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> VarLengthStringStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> VarLengthStringStruct: ...

    def __deepcopy__(self, arg: dict, /) -> VarLengthStringStruct: ...

    def __eq__(self, arg: VarLengthStringStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class WakeLrModeStruct:
    """Fortran struct: wake_lr_mode_struct"""

    def __init__(self, freq: float | None = None, freq_in: float | None = None, R_over_Q: float | None = None, Q: float | None = None, damp: float | None = None, phi: float | None = None, angle: float | None = None, b_sin: float | None = None, b_cos: float | None = None, a_sin: float | None = None, a_cos: float | None = None, m: int | None = None, polarized: bool | None = None) -> None: ...

    @property
    def freq(self) -> float:
        """Actual Frequency in Hz."""

    @freq.setter
    def freq(self, arg: float, /) -> None: ...

    @property
    def freq_in(self) -> float:
        """Input frequency in Hz."""

    @freq_in.setter
    def freq_in(self, arg: float, /) -> None: ...

    @property
    def R_over_Q(self) -> float:
        """Strength in V/C/m^(2*m_mode)."""

    @R_over_Q.setter
    def R_over_Q(self, arg: float, /) -> None: ...

    @property
    def Q(self) -> float:
        """Used for backwards compatability."""

    @Q.setter
    def Q(self, arg: float, /) -> None: ...

    @property
    def damp(self) -> float:
        """Damping factor = omega / 2 * Q = pi * freq / Q"""

    @damp.setter
    def damp(self, arg: float, /) -> None: ...

    @property
    def phi(self) -> float:
        """Phase in radians/2pi."""

    @phi.setter
    def phi(self, arg: float, /) -> None: ...

    @property
    def angle(self) -> float:
        """polarization angle (radians/2pi)."""

    @angle.setter
    def angle(self, arg: float, /) -> None: ...

    @property
    def b_sin(self) -> float:
        """non-skew sin-like component of the wake."""

    @b_sin.setter
    def b_sin(self, arg: float, /) -> None: ...

    @property
    def b_cos(self) -> float:
        """non-skew cos-like component of the wake."""

    @b_cos.setter
    def b_cos(self, arg: float, /) -> None: ...

    @property
    def a_sin(self) -> float:
        """skew sin-like component of the wake."""

    @a_sin.setter
    def a_sin(self, arg: float, /) -> None: ...

    @property
    def a_cos(self) -> float:
        """skew cos-like component of the wake."""

    @a_cos.setter
    def a_cos(self, arg: float, /) -> None: ...

    @property
    def m(self) -> int:
        """Mode order (1 = dipole, 2 = quad, etc.)"""

    @m.setter
    def m(self, arg: int, /) -> None: ...

    @property
    def polarized(self) -> bool:
        """Polaraized mode?"""

    @polarized.setter
    def polarized(self, arg: bool, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> WakeLrModeStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> WakeLrModeStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> WakeLrModeStruct: ...

    def __deepcopy__(self, arg: dict, /) -> WakeLrModeStruct: ...

    def __eq__(self, arg: WakeLrModeStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class WakeLrStruct:
    """Fortran struct: wake_lr_struct"""

    def __init__(self, file: str | None = None, t_ref: float | None = None, freq_spread: float | None = None, amp_scale: float | None = None, time_scale: float | None = None, self_wake_on: bool | None = None) -> None: ...

    @property
    def file(self) -> str: ...

    @file.setter
    def file(self, arg: str, /) -> None: ...

    @property
    def mode(self) -> WakeLrModeStructAlloc1D: ...

    @property
    def t_ref(self) -> float:
        """
        time reference value for computing the wake amplitude. This is used to prevent value overflow with long trains.
        """

    @t_ref.setter
    def t_ref(self, arg: float, /) -> None: ...

    @property
    def freq_spread(self) -> float:
        """Random frequency spread of long range modes."""

    @freq_spread.setter
    def freq_spread(self, arg: float, /) -> None: ...

    @property
    def amp_scale(self) -> float:
        """Wake amplitude scale factor."""

    @amp_scale.setter
    def amp_scale(self, arg: float, /) -> None: ...

    @property
    def time_scale(self) -> float:
        """time scale factor."""

    @time_scale.setter
    def time_scale(self, arg: float, /) -> None: ...

    @property
    def self_wake_on(self) -> bool:
        """Long range self-wake used in tracking?"""

    @self_wake_on.setter
    def self_wake_on(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> WakeLrStruct: ...

    def __deepcopy__(self, arg: dict, /) -> WakeLrStruct: ...

    def __eq__(self, arg: WakeLrStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class WakeSrModeStruct:
    """Fortran struct: wake_sr_mode_struct"""

    def __init__(self, amp: float | None = None, damp: float | None = None, k: float | None = None, phi: float | None = None, b_sin: float | None = None, b_cos: float | None = None, a_sin: float | None = None, a_cos: float | None = None, polarization: int | None = None, position_dependence: int | None = None) -> None: ...

    @property
    def amp(self) -> float:
        """Amplitude"""

    @amp.setter
    def amp(self, arg: float, /) -> None: ...

    @property
    def damp(self) -> float:
        """Dampling factor."""

    @damp.setter
    def damp(self, arg: float, /) -> None: ...

    @property
    def k(self) -> float:
        """k factor"""

    @k.setter
    def k(self, arg: float, /) -> None: ...

    @property
    def phi(self) -> float:
        """Phase in radians/2pi"""

    @phi.setter
    def phi(self, arg: float, /) -> None: ...

    @property
    def b_sin(self) -> float:
        """non-skew (x) sin-like component of the wake"""

    @b_sin.setter
    def b_sin(self, arg: float, /) -> None: ...

    @property
    def b_cos(self) -> float:
        """non-skew (x) cos-like component of the wake"""

    @b_cos.setter
    def b_cos(self, arg: float, /) -> None: ...

    @property
    def a_sin(self) -> float:
        """skew (y) sin-like component of the wake"""

    @a_sin.setter
    def a_sin(self, arg: float, /) -> None: ...

    @property
    def a_cos(self) -> float:
        """skew (y) cos-like component of the wake"""

    @a_cos.setter
    def a_cos(self, arg: float, /) -> None: ...

    @property
    def polarization(self) -> int:
        """Transverse: none$, x_axis$, y_axis$. Not used for longitudinal."""

    @polarization.setter
    def polarization(self, arg: int, /) -> None: ...

    @property
    def position_dependence(self) -> int:
        """
        Transverse: leading$, trailing$, none$ Longitudinal: x_leading$, ..., y_trailing$, none$
        """

    @position_dependence.setter
    def position_dependence(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> WakeSrModeStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> WakeSrModeStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> WakeSrModeStruct: ...

    def __deepcopy__(self, arg: dict, /) -> WakeSrModeStruct: ...

    def __eq__(self, arg: WakeSrModeStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class WakeSrStruct:
    """Fortran struct: wake_sr_struct"""

    def __init__(self, file: str | None = None, z_long: WakeSrZLongStruct | None = None, z_ref_long: float | None = None, z_ref_trans: float | None = None, z_max: float | None = None, amp_scale: float | None = None, z_scale: float | None = None, scale_with_length: bool | None = None) -> None: ...

    @property
    def file(self) -> str: ...

    @file.setter
    def file(self, arg: str, /) -> None: ...

    @property
    def z_long(self) -> WakeSrZLongStruct: ...

    @z_long.setter
    def z_long(self, arg: WakeSrZLongStruct, /) -> None: ...

    @property
    def long_wake(self) -> WakeSrModeStructAlloc1D: ...

    @property
    def trans_wake(self) -> WakeSrModeStructAlloc1D: ...

    @property
    def z_ref_long(self) -> float:
        """z reference value for computing the wake amplitude."""

    @z_ref_long.setter
    def z_ref_long(self, arg: float, /) -> None: ...

    @property
    def z_ref_trans(self) -> float:
        """This is used to prevent value overflow with long bunches."""

    @z_ref_trans.setter
    def z_ref_trans(self, arg: float, /) -> None: ...

    @property
    def z_max(self) -> float:
        """Max allowable z value. 0-> ignore"""

    @z_max.setter
    def z_max(self, arg: float, /) -> None: ...

    @property
    def amp_scale(self) -> float:
        """Wake amplitude scale factor."""

    @amp_scale.setter
    def amp_scale(self, arg: float, /) -> None: ...

    @property
    def z_scale(self) -> float:
        """z-distance scale factor."""

    @z_scale.setter
    def z_scale(self, arg: float, /) -> None: ...

    @property
    def scale_with_length(self) -> bool:
        """Scale wake with element length?"""

    @scale_with_length.setter
    def scale_with_length(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> WakeSrStruct: ...

    def __deepcopy__(self, arg: dict, /) -> WakeSrStruct: ...

    def __eq__(self, arg: WakeSrStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class WakeSrZLongStruct:
    """Fortran struct: wake_sr_z_long_struct"""

    def __init__(self, w: Sequence[float] | None = None, fw: Sequence[complex] | None = None, fbunch: Sequence[complex] | None = None, w_out: Sequence[complex] | None = None, dz: float | None = None, z0: float | None = None, smoothing_sigma: float | None = None, position_dependence: int | None = None, time_based: bool | None = None) -> None: ...

    @property
    def w(self) -> RealAlloc1D:
        """Input single particle Wake. Indexed from 1."""

    @w.setter
    def w(self, arg: Sequence[float], /) -> None: ...

    @property
    def fw(self) -> ComplexAlloc1D:
        """Fourier transform of w."""

    @fw.setter
    def fw(self, arg: Sequence[complex], /) -> None: ...

    @property
    def fbunch(self) -> ComplexAlloc1D:
        """Scratch space."""

    @fbunch.setter
    def fbunch(self, arg: Sequence[complex], /) -> None: ...

    @property
    def w_out(self) -> ComplexAlloc1D:
        """Scratch space."""

    @w_out.setter
    def w_out(self, arg: Sequence[complex], /) -> None: ...

    @property
    def dz(self) -> float:
        """Distance between points. If zero there is no wake."""

    @dz.setter
    def dz(self, arg: float, /) -> None: ...

    @property
    def z0(self) -> float:
        """Wake extent is [-z0, z0]."""

    @z0.setter
    def z0(self, arg: float, /) -> None: ...

    @property
    def smoothing_sigma(self) -> float:
        """0 => No smoothing."""

    @smoothing_sigma.setter
    def smoothing_sigma(self, arg: float, /) -> None: ...

    @property
    def position_dependence(self) -> int:
        """
        Transverse: leading$, trailing$, none$ Longitudinal: x_leading$, ..., y_trailing$, none$
        """

    @position_dependence.setter
    def position_dependence(self, arg: int, /) -> None: ...

    @property
    def time_based(self) -> bool:
        """Was input time based?"""

    @time_based.setter
    def time_based(self, arg: bool, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> WakeSrZLongStruct: ...

    def __deepcopy__(self, arg: dict, /) -> WakeSrZLongStruct: ...

    def __eq__(self, arg: WakeSrZLongStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class WakeStruct:
    """Fortran struct: wake_struct"""

    def __init__(self, sr: WakeSrStruct | None = None, lr: WakeLrStruct | None = None) -> None: ...

    @property
    def sr(self) -> WakeSrStruct:
        """Short-range wake"""

    @sr.setter
    def sr(self, arg: WakeSrStruct, /) -> None: ...

    @property
    def lr(self) -> WakeLrStruct:
        """Long-range wake"""

    @lr.setter
    def lr(self, arg: WakeLrStruct, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> WakeStruct: ...

    def __deepcopy__(self, arg: dict, /) -> WakeStruct: ...

    def __eq__(self, arg: WakeStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class Wall3DSectionStruct:
    """Fortran struct: wall3d_section_struct"""

    def __init__(self, name: str | None = None, material: str | None = None, surface: PhotonReflectSurfaceStruct | None = None, type: int | None = None, n_vertex_input: int | None = None, ix_ele: int | None = None, ix_branch: int | None = None, vertices_state: int | None = None, patch_in_region: bool | None = None, thickness: float | None = None, s: float | None = None, r0: Sequence[float] | None = None, dx0_ds: float | None = None, dy0_ds: float | None = None, x0_coef: Sequence[float] | None = None, y0_coef: Sequence[float] | None = None, dr_ds: float | None = None, p1_coef: Sequence[float] | None = None, p2_coef: Sequence[float] | None = None) -> None: ...

    @property
    def name(self) -> str:
        """Identifying name"""

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def material(self) -> str:
        """Material."""

    @material.setter
    def material(self, arg: str, /) -> None: ...

    @property
    def v(self) -> Wall3DVertexStructAlloc1D:
        """Array of vertices. Always stored relative."""

    @property
    def surface(self) -> PhotonReflectSurfaceStruct | None:
        """Surface reflectivity tables."""

    @surface.setter
    def surface(self, arg: PhotonReflectSurfaceStruct, /) -> None: ...

    @property
    def type(self) -> int:
        """normal$, clear$, opaque$, wall_start$, wall_end$"""

    @type.setter
    def type(self, arg: int, /) -> None: ...

    @property
    def n_vertex_input(self) -> int:
        """Number of vertices specified by the user."""

    @n_vertex_input.setter
    def n_vertex_input(self, arg: int, /) -> None: ...

    @property
    def ix_ele(self) -> int:
        """index of lattice element containing section"""

    @ix_ele.setter
    def ix_ele(self, arg: int, /) -> None: ...

    @property
    def ix_branch(self) -> int:
        """Index of branch lattice element is in."""

    @ix_branch.setter
    def ix_branch(self, arg: int, /) -> None: ...

    @property
    def vertices_state(self) -> int:
        """
        absolute$, or shifted_to_relative$. If set to absolute$ on input, will be changed to shifted_to_relative$ by section initalizer.
        """

    @vertices_state.setter
    def vertices_state(self, arg: int, /) -> None: ...

    @property
    def patch_in_region(self) -> bool:
        """Patch element exists between this section and previous one?"""

    @patch_in_region.setter
    def patch_in_region(self, arg: bool, /) -> None: ...

    @property
    def thickness(self) -> float:
        """Material thickness."""

    @thickness.setter
    def thickness(self, arg: float, /) -> None: ...

    @property
    def s(self) -> float:
        """Longitudinal position"""

    @s.setter
    def s(self, arg: float, /) -> None: ...

    @property
    def r0(self) -> RealArray1D:
        """
        Center of section Section-to-section spline interpolation of the center of the section
        """

    @r0.setter
    def r0(self, arg: Sequence[float], /) -> None: ...

    @property
    def dx0_ds(self) -> float:
        """Center of wall derivative"""

    @dx0_ds.setter
    def dx0_ds(self, arg: float, /) -> None: ...

    @property
    def dy0_ds(self) -> float:
        """Center of wall derivative"""

    @dy0_ds.setter
    def dy0_ds(self, arg: float, /) -> None: ...

    @property
    def x0_coef(self) -> RealArray1D:
        """Spline coefs for x-center"""

    @x0_coef.setter
    def x0_coef(self, arg: Sequence[float], /) -> None: ...

    @property
    def y0_coef(self) -> RealArray1D:
        """
        Spline coefs for y-center Section-to_section spline interpolation of the wall.
        """

    @y0_coef.setter
    def y0_coef(self, arg: Sequence[float], /) -> None: ...

    @property
    def dr_ds(self) -> float:
        """derivative of wall radius"""

    @dr_ds.setter
    def dr_ds(self, arg: float, /) -> None: ...

    @property
    def p1_coef(self) -> RealArray1D:
        """Spline coefs for p0 function"""

    @p1_coef.setter
    def p1_coef(self, arg: Sequence[float], /) -> None: ...

    @property
    def p2_coef(self) -> RealArray1D:
        """Spline coefs for p1 function"""

    @p2_coef.setter
    def p2_coef(self, arg: Sequence[float], /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> Wall3DSectionStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> Wall3DSectionStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> Wall3DSectionStruct: ...

    def __deepcopy__(self, arg: dict, /) -> Wall3DSectionStruct: ...

    def __eq__(self, arg: Wall3DSectionStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class Wall3DStruct:
    """Fortran struct: wall3d_struct"""

    def __init__(self, name: str | None = None, type: int | None = None, ix_wall3d: int | None = None, n_link: int | None = None, thickness: float | None = None, clear_material: str | None = None, opaque_material: str | None = None, superimpose: bool | None = None, ele_anchor_pt: int | None = None) -> None: ...

    @property
    def name(self) -> str: ...

    @name.setter
    def name(self, arg: str, /) -> None: ...

    @property
    def type(self) -> int:
        """or mask_plate$"""

    @type.setter
    def type(self, arg: int, /) -> None: ...

    @property
    def ix_wall3d(self) -> int:
        """Index in branch%wall3d(:) array."""

    @ix_wall3d.setter
    def ix_wall3d(self, arg: int, /) -> None: ...

    @property
    def n_link(self) -> int:
        """For memory management of ele%wall3d"""

    @n_link.setter
    def n_link(self, arg: int, /) -> None: ...

    @property
    def thickness(self) -> float:
        """For diffraction_plate elements"""

    @thickness.setter
    def thickness(self, arg: float, /) -> None: ...

    @property
    def clear_material(self) -> str: ...

    @clear_material.setter
    def clear_material(self, arg: str, /) -> None: ...

    @property
    def opaque_material(self) -> str: ...

    @opaque_material.setter
    def opaque_material(self, arg: str, /) -> None: ...

    @property
    def superimpose(self) -> bool:
        """Can overlap another wall"""

    @superimpose.setter
    def superimpose(self, arg: bool, /) -> None: ...

    @property
    def ele_anchor_pt(self) -> int:
        """anchor_beginning$, anchor_center$, or anchor_end$"""

    @ele_anchor_pt.setter
    def ele_anchor_pt(self, arg: int, /) -> None: ...

    @property
    def section(self) -> Wall3DSectionStructAlloc1D:
        """Indexed from 1."""

    @staticmethod
    def new_array1d(sz: int = 0) -> Wall3DStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> Wall3DStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> Wall3DStruct: ...

    def __deepcopy__(self, arg: dict, /) -> Wall3DStruct: ...

    def __eq__(self, arg: Wall3DStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class Wall3DVertexStruct:
    """Fortran struct: wall3d_vertex_struct"""

    def __init__(self, x: float | None = None, y: float | None = None, radius_x: float | None = None, radius_y: float | None = None, tilt: float | None = None, angle: float | None = None, x0: float | None = None, y0: float | None = None, type: int | None = None) -> None: ...

    @property
    def x(self) -> float:
        """Coordinates of the vertex."""

    @x.setter
    def x(self, arg: float, /) -> None: ...

    @property
    def y(self) -> float:
        """Coordinates of the vertex."""

    @y.setter
    def y(self, arg: float, /) -> None: ...

    @property
    def radius_x(self) -> float:
        """Radius of arc or ellipse x-axis half width. 0 => Straight line."""

    @radius_x.setter
    def radius_x(self, arg: float, /) -> None: ...

    @property
    def radius_y(self) -> float:
        """Ellipse y-axis half height."""

    @radius_y.setter
    def radius_y(self, arg: float, /) -> None: ...

    @property
    def tilt(self) -> float:
        """Tilt of ellipse"""

    @tilt.setter
    def tilt(self, arg: float, /) -> None: ...

    @property
    def angle(self) -> float:
        """Angle of (x, y) point."""

    @angle.setter
    def angle(self, arg: float, /) -> None: ...

    @property
    def x0(self) -> float:
        """Center of ellipse"""

    @x0.setter
    def x0(self, arg: float, /) -> None: ...

    @property
    def y0(self) -> float:
        """Center of ellipse"""

    @y0.setter
    def y0(self, arg: float, /) -> None: ...

    @property
    def type(self) -> int:
        """No longer used."""

    @type.setter
    def type(self, arg: int, /) -> None: ...

    @staticmethod
    def new_array1d(sz: int = 0) -> Wall3DVertexStructAlloc1D: ...

    @staticmethod
    def new_array1d_bounds(lbound: int, ubound: int) -> Wall3DVertexStructAlloc1D: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> Wall3DVertexStruct: ...

    def __deepcopy__(self, arg: dict, /) -> Wall3DVertexStruct: ...

    def __eq__(self, arg: Wall3DVertexStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class XyDispStruct:
    """Fortran struct: xy_disp_struct"""

    def __init__(self, eta: float | None = None, etap: float | None = None, deta_ds: float | None = None, sigma: float | None = None, deta_dpz: float | None = None, detap_dpz: float | None = None) -> None: ...

    @property
    def eta(self) -> float: ...

    @eta.setter
    def eta(self, arg: float, /) -> None: ...

    @property
    def etap(self) -> float: ...

    @etap.setter
    def etap(self, arg: float, /) -> None: ...

    @property
    def deta_ds(self) -> float: ...

    @deta_ds.setter
    def deta_ds(self, arg: float, /) -> None: ...

    @property
    def sigma(self) -> float: ...

    @sigma.setter
    def sigma(self, arg: float, /) -> None: ...

    @property
    def deta_dpz(self) -> float: ...

    @deta_dpz.setter
    def deta_dpz(self, arg: float, /) -> None: ...

    @property
    def detap_dpz(self) -> float: ...

    @detap_dpz.setter
    def detap_dpz(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

    def __copy__(self) -> XyDispStruct: ...

    def __deepcopy__(self, arg: dict, /) -> XyDispStruct: ...

    def __eq__(self, arg: XyDispStruct, /) -> bool: ...

    def __hash__(self) -> int: ...

class AcKickerFreqStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: AcKickerFreqStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> AcKickerFreqStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: AcKickerFreqStruct, /) -> None: ...

    def __iter__(self) -> Iterator[AcKickerFreqStruct]: ...

    def to_list(self) -> list[AcKickerFreqStruct]: ...

class AcKickerFreqStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: AcKickerFreqStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[AcKickerFreqStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> AcKickerFreqStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> AcKickerFreqStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: AcKickerFreqStruct, /) -> None: ...

    def __iter__(self) -> Iterator[AcKickerFreqStruct]: ...

class AcKickerTimeStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: AcKickerTimeStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> AcKickerTimeStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: AcKickerTimeStruct, /) -> None: ...

    def __iter__(self) -> Iterator[AcKickerTimeStruct]: ...

    def to_list(self) -> list[AcKickerTimeStruct]: ...

class AcKickerTimeStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: AcKickerTimeStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[AcKickerTimeStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> AcKickerTimeStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> AcKickerTimeStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: AcKickerTimeStruct, /) -> None: ...

    def __iter__(self) -> Iterator[AcKickerTimeStruct]: ...

class AllPointerStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: AllPointerStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> AllPointerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: AllPointerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[AllPointerStruct]: ...

    def to_list(self) -> list[AllPointerStruct]: ...

class AllPointerStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: AllPointerStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[AllPointerStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> AllPointerStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> AllPointerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: AllPointerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[AllPointerStruct]: ...

class AperturePointStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: AperturePointStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> AperturePointStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: AperturePointStruct, /) -> None: ...

    def __iter__(self) -> Iterator[AperturePointStruct]: ...

    def to_list(self) -> list[AperturePointStruct]: ...

class AperturePointStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: AperturePointStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[AperturePointStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> AperturePointStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> AperturePointStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: AperturePointStruct, /) -> None: ...

    def __iter__(self) -> Iterator[AperturePointStruct]: ...

class ApertureScanStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: ApertureScanStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> ApertureScanStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ApertureScanStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ApertureScanStruct]: ...

    def to_list(self) -> list[ApertureScanStruct]: ...

class ApertureScanStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: ApertureScanStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[ApertureScanStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> ApertureScanStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> ApertureScanStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ApertureScanStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ApertureScanStruct]: ...

class BaseLineEleStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: BaseLineEleStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> BaseLineEleStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: BaseLineEleStruct, /) -> None: ...

    def __iter__(self) -> Iterator[BaseLineEleStruct]: ...

    def to_list(self) -> list[BaseLineEleStruct]: ...

class BaseLineEleStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: BaseLineEleStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[BaseLineEleStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> BaseLineEleStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> BaseLineEleStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: BaseLineEleStruct, /) -> None: ...

    def __iter__(self) -> Iterator[BaseLineEleStruct]: ...

class BbuStageStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: BbuStageStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> BbuStageStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: BbuStageStruct, /) -> None: ...

    def __iter__(self) -> Iterator[BbuStageStruct]: ...

    def to_list(self) -> list[BbuStageStruct]: ...

class BbuStageStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: BbuStageStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[BbuStageStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> BbuStageStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> BbuStageStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: BbuStageStruct, /) -> None: ...

    def __iter__(self) -> Iterator[BbuStageStruct]: ...

class BicubicCmplxCoefStructArray3D:
    def __init__(self) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    @overload
    def __getitem__(self, arg: int, /) -> BicubicCmplxCoefStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: BicubicCmplxCoefStruct, /) -> None: ...

    def __iter__(self) -> Iterator[BicubicCmplxCoefStruct]: ...

    def __str__(self) -> str: ...

class BranchPointerStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: BranchPointerStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> BranchPointerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: BranchPointerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[BranchPointerStruct]: ...

    def to_list(self) -> list[BranchPointerStruct]: ...

class BranchPointerStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: BranchPointerStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[BranchPointerStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> BranchPointerStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> BranchPointerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: BranchPointerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[BranchPointerStruct]: ...

class BranchStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: BranchStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> BranchStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: BranchStruct, /) -> None: ...

    def __iter__(self) -> Iterator[BranchStruct]: ...

    def to_list(self) -> list[BranchStruct]: ...

class BranchStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: BranchStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[BranchStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> BranchStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> BranchStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: BranchStruct, /) -> None: ...

    def __iter__(self) -> Iterator[BranchStruct]: ...

class BunchParamsStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: BunchParamsStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> BunchParamsStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: BunchParamsStruct, /) -> None: ...

    def __iter__(self) -> Iterator[BunchParamsStruct]: ...

    def to_list(self) -> list[BunchParamsStruct]: ...

class BunchParamsStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: BunchParamsStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[BunchParamsStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> BunchParamsStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> BunchParamsStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: BunchParamsStruct, /) -> None: ...

    def __iter__(self) -> Iterator[BunchParamsStruct]: ...

class BunchStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: BunchStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> BunchStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: BunchStruct, /) -> None: ...

    def __iter__(self) -> Iterator[BunchStruct]: ...

    def to_list(self) -> list[BunchStruct]: ...

class BunchStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: BunchStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[BunchStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> BunchStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> BunchStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: BunchStruct, /) -> None: ...

    def __iter__(self) -> Iterator[BunchStruct]: ...

class BunchTrackStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: BunchTrackStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> BunchTrackStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: BunchTrackStruct, /) -> None: ...

    def __iter__(self) -> Iterator[BunchTrackStruct]: ...

    def to_list(self) -> list[BunchTrackStruct]: ...

class BunchTrackStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: BunchTrackStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[BunchTrackStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> BunchTrackStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> BunchTrackStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: BunchTrackStruct, /) -> None: ...

    def __iter__(self) -> Iterator[BunchTrackStruct]: ...

class CartesianMapStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: CartesianMapStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> CartesianMapStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CartesianMapStruct, /) -> None: ...

    def __iter__(self) -> Iterator[CartesianMapStruct]: ...

    def to_list(self) -> list[CartesianMapStruct]: ...

class CartesianMapStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: CartesianMapStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[CartesianMapStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> CartesianMapStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> CartesianMapStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CartesianMapStruct, /) -> None: ...

    def __iter__(self) -> Iterator[CartesianMapStruct]: ...

class CartesianMapTerm1StructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: CartesianMapTerm1StructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> CartesianMapTerm1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CartesianMapTerm1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[CartesianMapTerm1Struct]: ...

    def to_list(self) -> list[CartesianMapTerm1Struct]: ...

class CartesianMapTerm1StructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: CartesianMapTerm1StructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[CartesianMapTerm1Struct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> CartesianMapTerm1StructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> CartesianMapTerm1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CartesianMapTerm1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[CartesianMapTerm1Struct]: ...

class CmplxField1At2DPtStructArray2D:
    def __init__(self) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    @overload
    def __getitem__(self, arg: int, /) -> CmplxField1At2DPtStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CmplxField1At2DPtStruct, /) -> None: ...

    def __iter__(self) -> Iterator[CmplxField1At2DPtStruct]: ...

    def __str__(self) -> str: ...

class CmplxField1At3DPtStructArray3D:
    def __init__(self) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    @overload
    def __getitem__(self, arg: int, /) -> CmplxField1At3DPtStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CmplxField1At3DPtStruct, /) -> None: ...

    def __iter__(self) -> Iterator[CmplxField1At3DPtStruct]: ...

    def __str__(self) -> str: ...

class ComplexTaylorStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: ComplexTaylorStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> ComplexTaylorStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ComplexTaylorStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ComplexTaylorStruct]: ...

    def to_list(self) -> list[ComplexTaylorStruct]: ...

class ComplexTaylorStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: ComplexTaylorStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[ComplexTaylorStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> ComplexTaylorStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> ComplexTaylorStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ComplexTaylorStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ComplexTaylorStruct]: ...

class ComplexTaylorTermStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: ComplexTaylorTermStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> ComplexTaylorTermStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ComplexTaylorTermStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ComplexTaylorTermStruct]: ...

    def to_list(self) -> list[ComplexTaylorTermStruct]: ...

class ComplexTaylorTermStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: ComplexTaylorTermStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[ComplexTaylorTermStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> ComplexTaylorTermStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> ComplexTaylorTermStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ComplexTaylorTermStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ComplexTaylorTermStruct]: ...

class ControlRamp1StructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: ControlRamp1StructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> ControlRamp1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ControlRamp1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[ControlRamp1Struct]: ...

    def to_list(self) -> list[ControlRamp1Struct]: ...

class ControlRamp1StructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: ControlRamp1StructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[ControlRamp1Struct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> ControlRamp1StructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> ControlRamp1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ControlRamp1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[ControlRamp1Struct]: ...

class ControlStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: ControlStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> ControlStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ControlStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ControlStruct]: ...

    def to_list(self) -> list[ControlStruct]: ...

class ControlStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: ControlStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[ControlStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> ControlStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> ControlStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ControlStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ControlStruct]: ...

class ControlVar1StructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: ControlVar1StructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> ControlVar1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ControlVar1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[ControlVar1Struct]: ...

    def to_list(self) -> list[ControlVar1Struct]: ...

class ControlVar1StructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: ControlVar1StructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[ControlVar1Struct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> ControlVar1StructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> ControlVar1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ControlVar1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[ControlVar1Struct]: ...

class ConverterDir1DStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: ConverterDir1DStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> ConverterDir1DStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ConverterDir1DStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ConverterDir1DStruct]: ...

    def to_list(self) -> list[ConverterDir1DStruct]: ...

class ConverterDir1DStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: ConverterDir1DStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[ConverterDir1DStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> ConverterDir1DStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> ConverterDir1DStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ConverterDir1DStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ConverterDir1DStruct]: ...

class ConverterDistributionStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: ConverterDistributionStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> ConverterDistributionStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ConverterDistributionStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ConverterDistributionStruct]: ...

    def to_list(self) -> list[ConverterDistributionStruct]: ...

class ConverterDistributionStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: ConverterDistributionStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[ConverterDistributionStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> ConverterDistributionStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> ConverterDistributionStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ConverterDistributionStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ConverterDistributionStruct]: ...

class ConverterSubDistributionStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: ConverterSubDistributionStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> ConverterSubDistributionStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ConverterSubDistributionStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ConverterSubDistributionStruct]: ...

    def to_list(self) -> list[ConverterSubDistributionStruct]: ...

class ConverterSubDistributionStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: ConverterSubDistributionStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[ConverterSubDistributionStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> ConverterSubDistributionStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> ConverterSubDistributionStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ConverterSubDistributionStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ConverterSubDistributionStruct]: ...

class CoordArrayStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: CoordArrayStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> CoordArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CoordArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[CoordArrayStruct]: ...

    def to_list(self) -> list[CoordArrayStruct]: ...

class CoordArrayStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: CoordArrayStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[CoordArrayStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> CoordArrayStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> CoordArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CoordArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[CoordArrayStruct]: ...

class CoordStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: CoordStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> CoordStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CoordStruct, /) -> None: ...

    def __iter__(self) -> Iterator[CoordStruct]: ...

    def to_list(self) -> list[CoordStruct]: ...

class CoordStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: CoordStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[CoordStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> CoordStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> CoordStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CoordStruct, /) -> None: ...

    def __iter__(self) -> Iterator[CoordStruct]: ...

class CsrBunchSliceStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: CsrBunchSliceStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> CsrBunchSliceStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CsrBunchSliceStruct, /) -> None: ...

    def __iter__(self) -> Iterator[CsrBunchSliceStruct]: ...

    def to_list(self) -> list[CsrBunchSliceStruct]: ...

class CsrBunchSliceStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: CsrBunchSliceStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[CsrBunchSliceStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> CsrBunchSliceStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> CsrBunchSliceStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CsrBunchSliceStruct, /) -> None: ...

    def __iter__(self) -> Iterator[CsrBunchSliceStruct]: ...

class CsrEleInfoStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: CsrEleInfoStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> CsrEleInfoStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CsrEleInfoStruct, /) -> None: ...

    def __iter__(self) -> Iterator[CsrEleInfoStruct]: ...

    def to_list(self) -> list[CsrEleInfoStruct]: ...

class CsrEleInfoStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: CsrEleInfoStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[CsrEleInfoStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> CsrEleInfoStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> CsrEleInfoStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CsrEleInfoStruct, /) -> None: ...

    def __iter__(self) -> Iterator[CsrEleInfoStruct]: ...

class CsrKick1StructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: CsrKick1StructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> CsrKick1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CsrKick1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[CsrKick1Struct]: ...

    def to_list(self) -> list[CsrKick1Struct]: ...

class CsrKick1StructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: CsrKick1StructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[CsrKick1Struct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> CsrKick1StructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> CsrKick1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CsrKick1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[CsrKick1Struct]: ...

class CsrParticlePositionStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: CsrParticlePositionStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> CsrParticlePositionStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CsrParticlePositionStruct, /) -> None: ...

    def __iter__(self) -> Iterator[CsrParticlePositionStruct]: ...

    def to_list(self) -> list[CsrParticlePositionStruct]: ...

class CsrParticlePositionStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: CsrParticlePositionStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[CsrParticlePositionStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> CsrParticlePositionStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> CsrParticlePositionStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CsrParticlePositionStruct, /) -> None: ...

    def __iter__(self) -> Iterator[CsrParticlePositionStruct]: ...

class CylindricalMapStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: CylindricalMapStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> CylindricalMapStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CylindricalMapStruct, /) -> None: ...

    def __iter__(self) -> Iterator[CylindricalMapStruct]: ...

    def to_list(self) -> list[CylindricalMapStruct]: ...

class CylindricalMapStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: CylindricalMapStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[CylindricalMapStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> CylindricalMapStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> CylindricalMapStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CylindricalMapStruct, /) -> None: ...

    def __iter__(self) -> Iterator[CylindricalMapStruct]: ...

class CylindricalMapTerm1StructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: CylindricalMapTerm1StructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> CylindricalMapTerm1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CylindricalMapTerm1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[CylindricalMapTerm1Struct]: ...

    def to_list(self) -> list[CylindricalMapTerm1Struct]: ...

class CylindricalMapTerm1StructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: CylindricalMapTerm1StructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[CylindricalMapTerm1Struct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> CylindricalMapTerm1StructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> CylindricalMapTerm1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: CylindricalMapTerm1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[CylindricalMapTerm1Struct]: ...

class DoLoopStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: DoLoopStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> DoLoopStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: DoLoopStruct, /) -> None: ...

    def __iter__(self) -> Iterator[DoLoopStruct]: ...

    def to_list(self) -> list[DoLoopStruct]: ...

class DoLoopStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: DoLoopStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[DoLoopStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> DoLoopStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> DoLoopStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: DoLoopStruct, /) -> None: ...

    def __iter__(self) -> Iterator[DoLoopStruct]: ...

class ElePointerStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: ElePointerStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> ElePointerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ElePointerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ElePointerStruct]: ...

    def to_list(self) -> list[ElePointerStruct]: ...

class ElePointerStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: ElePointerStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[ElePointerStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> ElePointerStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> ElePointerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ElePointerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ElePointerStruct]: ...

class ElePointerStructArray2D:
    def __init__(self) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    @overload
    def __getitem__(self, arg: int, /) -> ElePointerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ElePointerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ElePointerStruct]: ...

    def __str__(self) -> str: ...

class EleStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: EleStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> EleStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: EleStruct, /) -> None: ...

    def __iter__(self) -> Iterator[EleStruct]: ...

    def to_list(self) -> list[EleStruct]: ...

class EleStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: EleStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[EleStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> EleStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> EleStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: EleStruct, /) -> None: ...

    def __iter__(self) -> Iterator[EleStruct]: ...

class EllipseBeamInitStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: EllipseBeamInitStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> EllipseBeamInitStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: EllipseBeamInitStruct, /) -> None: ...

    def __iter__(self) -> Iterator[EllipseBeamInitStruct]: ...

    def to_list(self) -> list[EllipseBeamInitStruct]: ...

class EllipseBeamInitStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: EllipseBeamInitStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[EllipseBeamInitStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> EllipseBeamInitStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> EllipseBeamInitStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: EllipseBeamInitStruct, /) -> None: ...

    def __iter__(self) -> Iterator[EllipseBeamInitStruct]: ...

class EmFieldStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: EmFieldStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> EmFieldStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: EmFieldStruct, /) -> None: ...

    def __iter__(self) -> Iterator[EmFieldStruct]: ...

    def to_list(self) -> list[EmFieldStruct]: ...

class EmFieldStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: EmFieldStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[EmFieldStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> EmFieldStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> EmFieldStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: EmFieldStruct, /) -> None: ...

    def __iter__(self) -> Iterator[EmFieldStruct]: ...

class ExpressionAtomStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: ExpressionAtomStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> ExpressionAtomStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ExpressionAtomStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ExpressionAtomStruct]: ...

    def to_list(self) -> list[ExpressionAtomStruct]: ...

class ExpressionAtomStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: ExpressionAtomStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[ExpressionAtomStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> ExpressionAtomStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> ExpressionAtomStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ExpressionAtomStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ExpressionAtomStruct]: ...

class ExpressionTreeStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: ExpressionTreeStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> ExpressionTreeStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ExpressionTreeStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ExpressionTreeStruct]: ...

    def to_list(self) -> list[ExpressionTreeStruct]: ...

class ExpressionTreeStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: ExpressionTreeStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[ExpressionTreeStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> ExpressionTreeStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> ExpressionTreeStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ExpressionTreeStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ExpressionTreeStruct]: ...

class Field1At2DPtStructArray2D:
    def __init__(self) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    @overload
    def __getitem__(self, arg: int, /) -> Field1At2DPtStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: Field1At2DPtStruct, /) -> None: ...

    def __iter__(self) -> Iterator[Field1At2DPtStruct]: ...

    def __str__(self) -> str: ...

class Field1At3DPtStructArray3D:
    def __init__(self) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    @overload
    def __getitem__(self, arg: int, /) -> Field1At3DPtStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: Field1At3DPtStruct, /) -> None: ...

    def __iter__(self) -> Iterator[Field1At3DPtStruct]: ...

    def __str__(self) -> str: ...

class GenGrad1StructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: GenGrad1StructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> GenGrad1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: GenGrad1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[GenGrad1Struct]: ...

    def to_list(self) -> list[GenGrad1Struct]: ...

class GenGrad1StructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: GenGrad1StructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[GenGrad1Struct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> GenGrad1StructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> GenGrad1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: GenGrad1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[GenGrad1Struct]: ...

class GenGradMapStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: GenGradMapStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> GenGradMapStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: GenGradMapStruct, /) -> None: ...

    def __iter__(self) -> Iterator[GenGradMapStruct]: ...

    def to_list(self) -> list[GenGradMapStruct]: ...

class GenGradMapStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: GenGradMapStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[GenGradMapStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> GenGradMapStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> GenGradMapStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: GenGradMapStruct, /) -> None: ...

    def __iter__(self) -> Iterator[GenGradMapStruct]: ...

class GgTaylorStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: GgTaylorStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> GgTaylorStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: GgTaylorStruct, /) -> None: ...

    def __iter__(self) -> Iterator[GgTaylorStruct]: ...

    def to_list(self) -> list[GgTaylorStruct]: ...

class GgTaylorStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: GgTaylorStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[GgTaylorStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> GgTaylorStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> GgTaylorStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: GgTaylorStruct, /) -> None: ...

    def __iter__(self) -> Iterator[GgTaylorStruct]: ...

class GgTaylorTermStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: GgTaylorTermStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> GgTaylorTermStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: GgTaylorTermStruct, /) -> None: ...

    def __iter__(self) -> Iterator[GgTaylorTermStruct]: ...

    def to_list(self) -> list[GgTaylorTermStruct]: ...

class GgTaylorTermStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: GgTaylorTermStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[GgTaylorTermStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> GgTaylorTermStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> GgTaylorTermStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: GgTaylorTermStruct, /) -> None: ...

    def __iter__(self) -> Iterator[GgTaylorTermStruct]: ...

class GridBeamInitStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: GridBeamInitStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> GridBeamInitStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: GridBeamInitStruct, /) -> None: ...

    def __iter__(self) -> Iterator[GridBeamInitStruct]: ...

    def to_list(self) -> list[GridBeamInitStruct]: ...

class GridBeamInitStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: GridBeamInitStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[GridBeamInitStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> GridBeamInitStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> GridBeamInitStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: GridBeamInitStruct, /) -> None: ...

    def __iter__(self) -> Iterator[GridBeamInitStruct]: ...

class GridFieldPt1StructArray3D:
    def __init__(self) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    @overload
    def __getitem__(self, arg: int, /) -> GridFieldPt1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: GridFieldPt1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[GridFieldPt1Struct]: ...

    def __str__(self) -> str: ...

class GridFieldStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: GridFieldStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> GridFieldStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: GridFieldStruct, /) -> None: ...

    def __iter__(self) -> Iterator[GridFieldStruct]: ...

    def to_list(self) -> list[GridFieldStruct]: ...

class GridFieldStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: GridFieldStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[GridFieldStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> GridFieldStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> GridFieldStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: GridFieldStruct, /) -> None: ...

    def __iter__(self) -> Iterator[GridFieldStruct]: ...

class Interval1CoefStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: Interval1CoefStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> Interval1CoefStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: Interval1CoefStruct, /) -> None: ...

    def __iter__(self) -> Iterator[Interval1CoefStruct]: ...

    def to_list(self) -> list[Interval1CoefStruct]: ...

class Interval1CoefStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: Interval1CoefStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[Interval1CoefStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> Interval1CoefStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> Interval1CoefStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: Interval1CoefStruct, /) -> None: ...

    def __iter__(self) -> Iterator[Interval1CoefStruct]: ...

class LatEleLocStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: LatEleLocStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> LatEleLocStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: LatEleLocStruct, /) -> None: ...

    def __iter__(self) -> Iterator[LatEleLocStruct]: ...

    def to_list(self) -> list[LatEleLocStruct]: ...

class LatEleLocStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: LatEleLocStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[LatEleLocStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> LatEleLocStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> LatEleLocStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: LatEleLocStruct, /) -> None: ...

    def __iter__(self) -> Iterator[LatEleLocStruct]: ...

class LatEleOrder1StructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: LatEleOrder1StructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> LatEleOrder1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: LatEleOrder1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[LatEleOrder1Struct]: ...

    def to_list(self) -> list[LatEleOrder1Struct]: ...

class LatEleOrder1StructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: LatEleOrder1StructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[LatEleOrder1Struct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> LatEleOrder1StructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> LatEleOrder1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: LatEleOrder1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[LatEleOrder1Struct]: ...

class LatEleOrderArrayStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: LatEleOrderArrayStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> LatEleOrderArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: LatEleOrderArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[LatEleOrderArrayStruct]: ...

    def to_list(self) -> list[LatEleOrderArrayStruct]: ...

class LatEleOrderArrayStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: LatEleOrderArrayStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[LatEleOrderArrayStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> LatEleOrderArrayStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> LatEleOrderArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: LatEleOrderArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[LatEleOrderArrayStruct]: ...

class LatPointerStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: LatPointerStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> LatPointerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: LatPointerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[LatPointerStruct]: ...

    def to_list(self) -> list[LatPointerStruct]: ...

class LatPointerStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: LatPointerStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[LatPointerStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> LatPointerStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> LatPointerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: LatPointerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[LatPointerStruct]: ...

class LatStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: LatStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> LatStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: LatStruct, /) -> None: ...

    def __iter__(self) -> Iterator[LatStruct]: ...

    def to_list(self) -> list[LatStruct]: ...

class LatStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: LatStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[LatStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> LatStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> LatStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: LatStruct, /) -> None: ...

    def __iter__(self) -> Iterator[LatStruct]: ...

class LinearEleIsfStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: LinearEleIsfStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> LinearEleIsfStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: LinearEleIsfStruct, /) -> None: ...

    def __iter__(self) -> Iterator[LinearEleIsfStruct]: ...

    def to_list(self) -> list[LinearEleIsfStruct]: ...

class LinearEleIsfStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: LinearEleIsfStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[LinearEleIsfStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> LinearEleIsfStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> LinearEleIsfStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: LinearEleIsfStruct, /) -> None: ...

    def __iter__(self) -> Iterator[LinearEleIsfStruct]: ...

class LinearIsf1StructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: LinearIsf1StructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> LinearIsf1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: LinearIsf1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[LinearIsf1Struct]: ...

    def to_list(self) -> list[LinearIsf1Struct]: ...

class LinearIsf1StructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: LinearIsf1StructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[LinearIsf1Struct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> LinearIsf1StructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> LinearIsf1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: LinearIsf1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[LinearIsf1Struct]: ...

class MaterialStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: MaterialStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> MaterialStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: MaterialStruct, /) -> None: ...

    def __iter__(self) -> Iterator[MaterialStruct]: ...

    def to_list(self) -> list[MaterialStruct]: ...

class MaterialStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: MaterialStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[MaterialStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> MaterialStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> MaterialStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: MaterialStruct, /) -> None: ...

    def __iter__(self) -> Iterator[MaterialStruct]: ...

class MolecularComponentStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: MolecularComponentStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> MolecularComponentStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: MolecularComponentStruct, /) -> None: ...

    def __iter__(self) -> Iterator[MolecularComponentStruct]: ...

    def to_list(self) -> list[MolecularComponentStruct]: ...

class MolecularComponentStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: MolecularComponentStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[MolecularComponentStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> MolecularComponentStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> MolecularComponentStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: MolecularComponentStruct, /) -> None: ...

    def __iter__(self) -> Iterator[MolecularComponentStruct]: ...

class MomentumApertureStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: MomentumApertureStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> MomentumApertureStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: MomentumApertureStruct, /) -> None: ...

    def __iter__(self) -> Iterator[MomentumApertureStruct]: ...

    def to_list(self) -> list[MomentumApertureStruct]: ...

class MomentumApertureStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: MomentumApertureStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[MomentumApertureStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> MomentumApertureStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> MomentumApertureStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: MomentumApertureStruct, /) -> None: ...

    def __iter__(self) -> Iterator[MomentumApertureStruct]: ...

class MultipassBranchInfoStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: MultipassBranchInfoStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> MultipassBranchInfoStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: MultipassBranchInfoStruct, /) -> None: ...

    def __iter__(self) -> Iterator[MultipassBranchInfoStruct]: ...

    def to_list(self) -> list[MultipassBranchInfoStruct]: ...

class MultipassBranchInfoStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: MultipassBranchInfoStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[MultipassBranchInfoStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> MultipassBranchInfoStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> MultipassBranchInfoStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: MultipassBranchInfoStruct, /) -> None: ...

    def __iter__(self) -> Iterator[MultipassBranchInfoStruct]: ...

class MultipassEleInfoStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: MultipassEleInfoStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> MultipassEleInfoStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: MultipassEleInfoStruct, /) -> None: ...

    def __iter__(self) -> Iterator[MultipassEleInfoStruct]: ...

    def to_list(self) -> list[MultipassEleInfoStruct]: ...

class MultipassEleInfoStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: MultipassEleInfoStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[MultipassEleInfoStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> MultipassEleInfoStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> MultipassEleInfoStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: MultipassEleInfoStruct, /) -> None: ...

    def __iter__(self) -> Iterator[MultipassEleInfoStruct]: ...

class MultipassLordInfoStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: MultipassLordInfoStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> MultipassLordInfoStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: MultipassLordInfoStruct, /) -> None: ...

    def __iter__(self) -> Iterator[MultipassLordInfoStruct]: ...

    def to_list(self) -> list[MultipassLordInfoStruct]: ...

class MultipassLordInfoStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: MultipassLordInfoStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[MultipassLordInfoStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> MultipassLordInfoStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> MultipassLordInfoStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: MultipassLordInfoStruct, /) -> None: ...

    def __iter__(self) -> Iterator[MultipassLordInfoStruct]: ...

class MultipassRegionBranchStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: MultipassRegionBranchStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> MultipassRegionBranchStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: MultipassRegionBranchStruct, /) -> None: ...

    def __iter__(self) -> Iterator[MultipassRegionBranchStruct]: ...

    def to_list(self) -> list[MultipassRegionBranchStruct]: ...

class MultipassRegionBranchStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: MultipassRegionBranchStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[MultipassRegionBranchStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> MultipassRegionBranchStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> MultipassRegionBranchStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: MultipassRegionBranchStruct, /) -> None: ...

    def __iter__(self) -> Iterator[MultipassRegionBranchStruct]: ...

class MultipassRegionEleStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: MultipassRegionEleStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> MultipassRegionEleStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: MultipassRegionEleStruct, /) -> None: ...

    def __iter__(self) -> Iterator[MultipassRegionEleStruct]: ...

    def to_list(self) -> list[MultipassRegionEleStruct]: ...

class MultipassRegionEleStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: MultipassRegionEleStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[MultipassRegionEleStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> MultipassRegionEleStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> MultipassRegionEleStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: MultipassRegionEleStruct, /) -> None: ...

    def __iter__(self) -> Iterator[MultipassRegionEleStruct]: ...

class NamedNumberStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: NamedNumberStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> NamedNumberStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: NamedNumberStruct, /) -> None: ...

    def __iter__(self) -> Iterator[NamedNumberStruct]: ...

    def to_list(self) -> list[NamedNumberStruct]: ...

class NamedNumberStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: NamedNumberStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[NamedNumberStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> NamedNumberStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> NamedNumberStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: NamedNumberStruct, /) -> None: ...

    def __iter__(self) -> Iterator[NamedNumberStruct]: ...

class ParserControllerStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: ParserControllerStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> ParserControllerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ParserControllerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ParserControllerStruct]: ...

    def to_list(self) -> list[ParserControllerStruct]: ...

class ParserControllerStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: ParserControllerStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[ParserControllerStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> ParserControllerStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> ParserControllerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ParserControllerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ParserControllerStruct]: ...

class ParserEleStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: ParserEleStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> ParserEleStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ParserEleStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ParserEleStruct]: ...

    def to_list(self) -> list[ParserEleStruct]: ...

class ParserEleStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: ParserEleStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[ParserEleStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> ParserEleStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> ParserEleStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ParserEleStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ParserEleStruct]: ...

class PhotonInitXAngleSplineStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: PhotonInitXAngleSplineStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> PhotonInitXAngleSplineStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: PhotonInitXAngleSplineStruct, /) -> None: ...

    def __iter__(self) -> Iterator[PhotonInitXAngleSplineStruct]: ...

    def to_list(self) -> list[PhotonInitXAngleSplineStruct]: ...

class PhotonInitXAngleSplineStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: PhotonInitXAngleSplineStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[PhotonInitXAngleSplineStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> PhotonInitXAngleSplineStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> PhotonInitXAngleSplineStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: PhotonInitXAngleSplineStruct, /) -> None: ...

    def __iter__(self) -> Iterator[PhotonInitXAngleSplineStruct]: ...

class PhotonInitYAngleSplineStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: PhotonInitYAngleSplineStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> PhotonInitYAngleSplineStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: PhotonInitYAngleSplineStruct, /) -> None: ...

    def __iter__(self) -> Iterator[PhotonInitYAngleSplineStruct]: ...

    def to_list(self) -> list[PhotonInitYAngleSplineStruct]: ...

class PhotonInitYAngleSplineStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: PhotonInitYAngleSplineStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[PhotonInitYAngleSplineStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> PhotonInitYAngleSplineStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> PhotonInitYAngleSplineStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: PhotonInitYAngleSplineStruct, /) -> None: ...

    def __iter__(self) -> Iterator[PhotonInitYAngleSplineStruct]: ...

class PhotonReflectTableStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: PhotonReflectTableStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> PhotonReflectTableStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: PhotonReflectTableStruct, /) -> None: ...

    def __iter__(self) -> Iterator[PhotonReflectTableStruct]: ...

    def to_list(self) -> list[PhotonReflectTableStruct]: ...

class PhotonReflectTableStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: PhotonReflectTableStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[PhotonReflectTableStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> PhotonReflectTableStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> PhotonReflectTableStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: PhotonReflectTableStruct, /) -> None: ...

    def __iter__(self) -> Iterator[PhotonReflectTableStruct]: ...

class PixelPtStructArray2D:
    def __init__(self) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    @overload
    def __getitem__(self, arg: int, /) -> PixelPtStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: PixelPtStruct, /) -> None: ...

    def __iter__(self) -> Iterator[PixelPtStruct]: ...

    def __str__(self) -> str: ...

class PtcLayoutPointerStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: PtcLayoutPointerStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> PtcLayoutPointerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: PtcLayoutPointerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[PtcLayoutPointerStruct]: ...

    def to_list(self) -> list[PtcLayoutPointerStruct]: ...

class PtcLayoutPointerStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: PtcLayoutPointerStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[PtcLayoutPointerStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> PtcLayoutPointerStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> PtcLayoutPointerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: PtcLayoutPointerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[PtcLayoutPointerStruct]: ...

class RadInt1StructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: RadInt1StructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> RadInt1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: RadInt1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[RadInt1Struct]: ...

    def to_list(self) -> list[RadInt1Struct]: ...

class RadInt1StructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: RadInt1StructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[RadInt1Struct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> RadInt1StructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> RadInt1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: RadInt1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[RadInt1Struct]: ...

class RadIntBranchStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: RadIntBranchStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> RadIntBranchStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: RadIntBranchStruct, /) -> None: ...

    def __iter__(self) -> Iterator[RadIntBranchStruct]: ...

    def to_list(self) -> list[RadIntBranchStruct]: ...

class RadIntBranchStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: RadIntBranchStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[RadIntBranchStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> RadIntBranchStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> RadIntBranchStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: RadIntBranchStruct, /) -> None: ...

    def __iter__(self) -> Iterator[RadIntBranchStruct]: ...

class RadIntTrackPointStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: RadIntTrackPointStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> RadIntTrackPointStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: RadIntTrackPointStruct, /) -> None: ...

    def __iter__(self) -> Iterator[RadIntTrackPointStruct]: ...

    def to_list(self) -> list[RadIntTrackPointStruct]: ...

class RadIntTrackPointStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: RadIntTrackPointStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[RadIntTrackPointStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> RadIntTrackPointStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> RadIntTrackPointStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: RadIntTrackPointStruct, /) -> None: ...

    def __iter__(self) -> Iterator[RadIntTrackPointStruct]: ...

class RamperLordStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: RamperLordStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> RamperLordStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: RamperLordStruct, /) -> None: ...

    def __iter__(self) -> Iterator[RamperLordStruct]: ...

    def to_list(self) -> list[RamperLordStruct]: ...

class RamperLordStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: RamperLordStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[RamperLordStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> RamperLordStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> RamperLordStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: RamperLordStruct, /) -> None: ...

    def __iter__(self) -> Iterator[RamperLordStruct]: ...

class ResonanceHStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: ResonanceHStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> ResonanceHStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ResonanceHStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ResonanceHStruct]: ...

    def to_list(self) -> list[ResonanceHStruct]: ...

class ResonanceHStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: ResonanceHStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[ResonanceHStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> ResonanceHStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> ResonanceHStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: ResonanceHStruct, /) -> None: ...

    def __iter__(self) -> Iterator[ResonanceHStruct]: ...

class RfStairStepStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: RfStairStepStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> RfStairStepStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: RfStairStepStruct, /) -> None: ...

    def __iter__(self) -> Iterator[RfStairStepStruct]: ...

    def to_list(self) -> list[RfStairStepStruct]: ...

class RfStairStepStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: RfStairStepStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[RfStairStepStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> RfStairStepStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> RfStairStepStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: RfStairStepStruct, /) -> None: ...

    def __iter__(self) -> Iterator[RfStairStepStruct]: ...

class SeqEleStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: SeqEleStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> SeqEleStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: SeqEleStruct, /) -> None: ...

    def __iter__(self) -> Iterator[SeqEleStruct]: ...

    def to_list(self) -> list[SeqEleStruct]: ...

class SeqEleStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: SeqEleStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[SeqEleStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> SeqEleStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> SeqEleStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: SeqEleStruct, /) -> None: ...

    def __iter__(self) -> Iterator[SeqEleStruct]: ...

class SeqStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: SeqStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> SeqStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: SeqStruct, /) -> None: ...

    def __iter__(self) -> Iterator[SeqStruct]: ...

    def to_list(self) -> list[SeqStruct]: ...

class SeqStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: SeqStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[SeqStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> SeqStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> SeqStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: SeqStruct, /) -> None: ...

    def __iter__(self) -> Iterator[SeqStruct]: ...

class SpinEigenStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: SpinEigenStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> SpinEigenStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: SpinEigenStruct, /) -> None: ...

    def __iter__(self) -> Iterator[SpinEigenStruct]: ...

    def to_list(self) -> list[SpinEigenStruct]: ...

class SpinEigenStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: SpinEigenStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[SpinEigenStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> SpinEigenStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> SpinEigenStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: SpinEigenStruct, /) -> None: ...

    def __iter__(self) -> Iterator[SpinEigenStruct]: ...

class SpinMatchingStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: SpinMatchingStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> SpinMatchingStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: SpinMatchingStruct, /) -> None: ...

    def __iter__(self) -> Iterator[SpinMatchingStruct]: ...

    def to_list(self) -> list[SpinMatchingStruct]: ...

class SpinMatchingStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: SpinMatchingStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[SpinMatchingStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> SpinMatchingStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> SpinMatchingStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: SpinMatchingStruct, /) -> None: ...

    def __iter__(self) -> Iterator[SpinMatchingStruct]: ...

class SpinOrbitMap1StructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: SpinOrbitMap1StructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> SpinOrbitMap1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: SpinOrbitMap1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[SpinOrbitMap1Struct]: ...

    def to_list(self) -> list[SpinOrbitMap1Struct]: ...

class SpinOrbitMap1StructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: SpinOrbitMap1StructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[SpinOrbitMap1Struct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> SpinOrbitMap1StructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> SpinOrbitMap1Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: SpinOrbitMap1Struct, /) -> None: ...

    def __iter__(self) -> Iterator[SpinOrbitMap1Struct]: ...

class SplineStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: SplineStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> SplineStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: SplineStruct, /) -> None: ...

    def __iter__(self) -> Iterator[SplineStruct]: ...

    def to_list(self) -> list[SplineStruct]: ...

class SplineStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: SplineStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[SplineStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> SplineStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> SplineStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: SplineStruct, /) -> None: ...

    def __iter__(self) -> Iterator[SplineStruct]: ...

class SummationRdtStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: SummationRdtStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> SummationRdtStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: SummationRdtStruct, /) -> None: ...

    def __iter__(self) -> Iterator[SummationRdtStruct]: ...

    def to_list(self) -> list[SummationRdtStruct]: ...

class SummationRdtStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: SummationRdtStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[SummationRdtStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> SummationRdtStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> SummationRdtStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: SummationRdtStruct, /) -> None: ...

    def __iter__(self) -> Iterator[SummationRdtStruct]: ...

class SurfaceDisplacementPtStructArray2D:
    def __init__(self) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    @overload
    def __getitem__(self, arg: int, /) -> SurfaceDisplacementPtStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: SurfaceDisplacementPtStruct, /) -> None: ...

    def __iter__(self) -> Iterator[SurfaceDisplacementPtStruct]: ...

    def __str__(self) -> str: ...

class SurfaceHMisalignPtStructArray2D:
    def __init__(self) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    @overload
    def __getitem__(self, arg: int, /) -> SurfaceHMisalignPtStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: SurfaceHMisalignPtStruct, /) -> None: ...

    def __iter__(self) -> Iterator[SurfaceHMisalignPtStruct]: ...

    def __str__(self) -> str: ...

class SurfaceSegmentedPtStructArray2D:
    def __init__(self) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    @overload
    def __getitem__(self, arg: int, /) -> SurfaceSegmentedPtStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: SurfaceSegmentedPtStruct, /) -> None: ...

    def __iter__(self) -> Iterator[SurfaceSegmentedPtStruct]: ...

    def __str__(self) -> str: ...

class TaoAliasStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoAliasStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoAliasStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoAliasStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoAliasStruct]: ...

    def to_list(self) -> list[TaoAliasStruct]: ...

class TaoAliasStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoAliasStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoAliasStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoAliasStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoAliasStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoAliasStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoAliasStruct]: ...

class TaoBuildingWallPointStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoBuildingWallPointStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoBuildingWallPointStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoBuildingWallPointStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoBuildingWallPointStruct]: ...

    def to_list(self) -> list[TaoBuildingWallPointStruct]: ...

class TaoBuildingWallPointStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoBuildingWallPointStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoBuildingWallPointStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoBuildingWallPointStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoBuildingWallPointStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoBuildingWallPointStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoBuildingWallPointStruct]: ...

class TaoBuildingWallSectionStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoBuildingWallSectionStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoBuildingWallSectionStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoBuildingWallSectionStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoBuildingWallSectionStruct]: ...

    def to_list(self) -> list[TaoBuildingWallSectionStruct]: ...

class TaoBuildingWallSectionStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoBuildingWallSectionStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoBuildingWallSectionStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoBuildingWallSectionStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoBuildingWallSectionStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoBuildingWallSectionStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoBuildingWallSectionStruct]: ...

class TaoCmdHistoryStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoCmdHistoryStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoCmdHistoryStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoCmdHistoryStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoCmdHistoryStruct]: ...

    def to_list(self) -> list[TaoCmdHistoryStruct]: ...

class TaoCmdHistoryStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoCmdHistoryStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoCmdHistoryStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoCmdHistoryStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoCmdHistoryStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoCmdHistoryStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoCmdHistoryStruct]: ...

class TaoCommandFileStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoCommandFileStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoCommandFileStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoCommandFileStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoCommandFileStruct]: ...

    def to_list(self) -> list[TaoCommandFileStruct]: ...

class TaoCommandFileStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoCommandFileStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoCommandFileStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoCommandFileStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoCommandFileStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoCommandFileStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoCommandFileStruct]: ...

class TaoCurveArrayStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoCurveArrayStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoCurveArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoCurveArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoCurveArrayStruct]: ...

    def to_list(self) -> list[TaoCurveArrayStruct]: ...

class TaoCurveArrayStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoCurveArrayStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoCurveArrayStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoCurveArrayStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoCurveArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoCurveArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoCurveArrayStruct]: ...

class TaoCurveStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoCurveStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoCurveStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoCurveStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoCurveStruct]: ...

    def to_list(self) -> list[TaoCurveStruct]: ...

class TaoCurveStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoCurveStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoCurveStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoCurveStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoCurveStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoCurveStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoCurveStruct]: ...

class TaoD1DataArrayStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoD1DataArrayStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoD1DataArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoD1DataArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoD1DataArrayStruct]: ...

    def to_list(self) -> list[TaoD1DataArrayStruct]: ...

class TaoD1DataArrayStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoD1DataArrayStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoD1DataArrayStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoD1DataArrayStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoD1DataArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoD1DataArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoD1DataArrayStruct]: ...

class TaoD1DataStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoD1DataStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoD1DataStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoD1DataStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoD1DataStruct]: ...

    def to_list(self) -> list[TaoD1DataStruct]: ...

class TaoD1DataStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoD1DataStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoD1DataStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoD1DataStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoD1DataStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoD1DataStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoD1DataStruct]: ...

class TaoD2DataArrayStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoD2DataArrayStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoD2DataArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoD2DataArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoD2DataArrayStruct]: ...

    def to_list(self) -> list[TaoD2DataArrayStruct]: ...

class TaoD2DataArrayStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoD2DataArrayStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoD2DataArrayStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoD2DataArrayStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoD2DataArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoD2DataArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoD2DataArrayStruct]: ...

class TaoD2DataStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoD2DataStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoD2DataStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoD2DataStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoD2DataStruct]: ...

    def to_list(self) -> list[TaoD2DataStruct]: ...

class TaoD2DataStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoD2DataStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoD2DataStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoD2DataStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoD2DataStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoD2DataStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoD2DataStruct]: ...

class TaoDataArrayStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoDataArrayStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoDataArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoDataArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoDataArrayStruct]: ...

    def to_list(self) -> list[TaoDataArrayStruct]: ...

class TaoDataArrayStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoDataArrayStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoDataArrayStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoDataArrayStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoDataArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoDataArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoDataArrayStruct]: ...

class TaoDataStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoDataStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoDataStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoDataStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoDataStruct]: ...

    def to_list(self) -> list[TaoDataStruct]: ...

class TaoDataStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoDataStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoDataStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoDataStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoDataStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoDataStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoDataStruct]: ...

class TaoDataVarComponentStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoDataVarComponentStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoDataVarComponentStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoDataVarComponentStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoDataVarComponentStruct]: ...

    def to_list(self) -> list[TaoDataVarComponentStruct]: ...

class TaoDataVarComponentStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoDataVarComponentStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoDataVarComponentStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoDataVarComponentStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoDataVarComponentStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoDataVarComponentStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoDataVarComponentStruct]: ...

class TaoElePointerStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoElePointerStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoElePointerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoElePointerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoElePointerStruct]: ...

    def to_list(self) -> list[TaoElePointerStruct]: ...

class TaoElePointerStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoElePointerStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoElePointerStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoElePointerStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoElePointerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoElePointerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoElePointerStruct]: ...

class TaoEleShapeStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoEleShapeStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoEleShapeStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoEleShapeStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoEleShapeStruct]: ...

    def to_list(self) -> list[TaoEleShapeStruct]: ...

class TaoEleShapeStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoEleShapeStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoEleShapeStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoEleShapeStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoEleShapeStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoEleShapeStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoEleShapeStruct]: ...

class TaoEvalNodeStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoEvalNodeStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoEvalNodeStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoEvalNodeStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoEvalNodeStruct]: ...

    def to_list(self) -> list[TaoEvalNodeStruct]: ...

class TaoEvalNodeStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoEvalNodeStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoEvalNodeStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoEvalNodeStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoEvalNodeStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoEvalNodeStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoEvalNodeStruct]: ...

class TaoExpressionInfoStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoExpressionInfoStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoExpressionInfoStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoExpressionInfoStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoExpressionInfoStruct]: ...

    def to_list(self) -> list[TaoExpressionInfoStruct]: ...

class TaoExpressionInfoStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoExpressionInfoStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoExpressionInfoStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoExpressionInfoStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoExpressionInfoStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoExpressionInfoStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoExpressionInfoStruct]: ...

class TaoGraphArrayStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoGraphArrayStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoGraphArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoGraphArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoGraphArrayStruct]: ...

    def to_list(self) -> list[TaoGraphArrayStruct]: ...

class TaoGraphArrayStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoGraphArrayStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoGraphArrayStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoGraphArrayStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoGraphArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoGraphArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoGraphArrayStruct]: ...

class TaoGraphStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoGraphStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoGraphStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoGraphStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoGraphStruct]: ...

    def to_list(self) -> list[TaoGraphStruct]: ...

class TaoGraphStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoGraphStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoGraphStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoGraphStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoGraphStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoGraphStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoGraphStruct]: ...

class TaoIntegerArrayStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoIntegerArrayStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoIntegerArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoIntegerArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoIntegerArrayStruct]: ...

    def to_list(self) -> list[TaoIntegerArrayStruct]: ...

class TaoIntegerArrayStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoIntegerArrayStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoIntegerArrayStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoIntegerArrayStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoIntegerArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoIntegerArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoIntegerArrayStruct]: ...

class TaoLatSigmaStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoLatSigmaStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoLatSigmaStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoLatSigmaStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoLatSigmaStruct]: ...

    def to_list(self) -> list[TaoLatSigmaStruct]: ...

class TaoLatSigmaStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoLatSigmaStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoLatSigmaStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoLatSigmaStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoLatSigmaStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoLatSigmaStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoLatSigmaStruct]: ...

class TaoLatticeBranchStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoLatticeBranchStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoLatticeBranchStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoLatticeBranchStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoLatticeBranchStruct]: ...

    def to_list(self) -> list[TaoLatticeBranchStruct]: ...

class TaoLatticeBranchStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoLatticeBranchStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoLatticeBranchStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoLatticeBranchStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoLatticeBranchStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoLatticeBranchStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoLatticeBranchStruct]: ...

class TaoLogicalArrayStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoLogicalArrayStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoLogicalArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoLogicalArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoLogicalArrayStruct]: ...

    def to_list(self) -> list[TaoLogicalArrayStruct]: ...

class TaoLogicalArrayStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoLogicalArrayStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoLogicalArrayStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoLogicalArrayStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoLogicalArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoLogicalArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoLogicalArrayStruct]: ...

class TaoModelBranchStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoModelBranchStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoModelBranchStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoModelBranchStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoModelBranchStruct]: ...

    def to_list(self) -> list[TaoModelBranchStruct]: ...

class TaoModelBranchStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoModelBranchStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoModelBranchStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoModelBranchStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoModelBranchStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoModelBranchStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoModelBranchStruct]: ...

class TaoModelElementStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoModelElementStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoModelElementStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoModelElementStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoModelElementStruct]: ...

    def to_list(self) -> list[TaoModelElementStruct]: ...

class TaoModelElementStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoModelElementStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoModelElementStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoModelElementStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoModelElementStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoModelElementStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoModelElementStruct]: ...

class TaoPlotArrayStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoPlotArrayStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoPlotArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoPlotArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoPlotArrayStruct]: ...

    def to_list(self) -> list[TaoPlotArrayStruct]: ...

class TaoPlotArrayStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoPlotArrayStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoPlotArrayStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoPlotArrayStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoPlotArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoPlotArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoPlotArrayStruct]: ...

class TaoPlotCacheStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoPlotCacheStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoPlotCacheStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoPlotCacheStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoPlotCacheStruct]: ...

    def to_list(self) -> list[TaoPlotCacheStruct]: ...

class TaoPlotCacheStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoPlotCacheStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoPlotCacheStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoPlotCacheStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoPlotCacheStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoPlotCacheStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoPlotCacheStruct]: ...

class TaoPlotRegionStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoPlotRegionStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoPlotRegionStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoPlotRegionStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoPlotRegionStruct]: ...

    def to_list(self) -> list[TaoPlotRegionStruct]: ...

class TaoPlotRegionStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoPlotRegionStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoPlotRegionStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoPlotRegionStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoPlotRegionStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoPlotRegionStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoPlotRegionStruct]: ...

class TaoPlotStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoPlotStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoPlotStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoPlotStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoPlotStruct]: ...

    def to_list(self) -> list[TaoPlotStruct]: ...

class TaoPlotStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoPlotStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoPlotStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoPlotStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoPlotStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoPlotStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoPlotStruct]: ...

class TaoRealPointerStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoRealPointerStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoRealPointerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoRealPointerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoRealPointerStruct]: ...

    def to_list(self) -> list[TaoRealPointerStruct]: ...

class TaoRealPointerStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoRealPointerStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoRealPointerStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoRealPointerStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoRealPointerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoRealPointerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoRealPointerStruct]: ...

class TaoShapePatternPointStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoShapePatternPointStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoShapePatternPointStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoShapePatternPointStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoShapePatternPointStruct]: ...

    def to_list(self) -> list[TaoShapePatternPointStruct]: ...

class TaoShapePatternPointStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoShapePatternPointStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoShapePatternPointStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoShapePatternPointStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoShapePatternPointStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoShapePatternPointStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoShapePatternPointStruct]: ...

class TaoShapePatternStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoShapePatternStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoShapePatternStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoShapePatternStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoShapePatternStruct]: ...

    def to_list(self) -> list[TaoShapePatternStruct]: ...

class TaoShapePatternStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoShapePatternStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoShapePatternStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoShapePatternStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoShapePatternStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoShapePatternStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoShapePatternStruct]: ...

class TaoSpinEleStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoSpinEleStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoSpinEleStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoSpinEleStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoSpinEleStruct]: ...

    def to_list(self) -> list[TaoSpinEleStruct]: ...

class TaoSpinEleStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoSpinEleStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoSpinEleStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoSpinEleStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoSpinEleStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoSpinEleStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoSpinEleStruct]: ...

class TaoStringArrayStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoStringArrayStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoStringArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoStringArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoStringArrayStruct]: ...

    def to_list(self) -> list[TaoStringArrayStruct]: ...

class TaoStringArrayStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoStringArrayStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoStringArrayStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoStringArrayStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoStringArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoStringArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoStringArrayStruct]: ...

class TaoTop10StructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoTop10StructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoTop10Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoTop10Struct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoTop10Struct]: ...

    def to_list(self) -> list[TaoTop10Struct]: ...

class TaoTop10StructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoTop10StructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoTop10Struct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoTop10StructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoTop10Struct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoTop10Struct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoTop10Struct]: ...

class TaoUniversePointerStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoUniversePointerStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoUniversePointerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoUniversePointerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoUniversePointerStruct]: ...

    def to_list(self) -> list[TaoUniversePointerStruct]: ...

class TaoUniversePointerStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoUniversePointerStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoUniversePointerStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoUniversePointerStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoUniversePointerStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoUniversePointerStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoUniversePointerStruct]: ...

class TaoUniverseStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoUniverseStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoUniverseStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoUniverseStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoUniverseStruct]: ...

    def to_list(self) -> list[TaoUniverseStruct]: ...

class TaoUniverseStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoUniverseStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoUniverseStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoUniverseStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoUniverseStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoUniverseStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoUniverseStruct]: ...

class TaoV1VarArrayStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoV1VarArrayStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoV1VarArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoV1VarArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoV1VarArrayStruct]: ...

    def to_list(self) -> list[TaoV1VarArrayStruct]: ...

class TaoV1VarArrayStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoV1VarArrayStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoV1VarArrayStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoV1VarArrayStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoV1VarArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoV1VarArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoV1VarArrayStruct]: ...

class TaoV1VarStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoV1VarStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoV1VarStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoV1VarStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoV1VarStruct]: ...

    def to_list(self) -> list[TaoV1VarStruct]: ...

class TaoV1VarStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoV1VarStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoV1VarStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoV1VarStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoV1VarStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoV1VarStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoV1VarStruct]: ...

class TaoVarArrayStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoVarArrayStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoVarArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoVarArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoVarArrayStruct]: ...

    def to_list(self) -> list[TaoVarArrayStruct]: ...

class TaoVarArrayStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoVarArrayStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoVarArrayStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoVarArrayStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoVarArrayStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoVarArrayStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoVarArrayStruct]: ...

class TaoVarSlaveStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoVarSlaveStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoVarSlaveStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoVarSlaveStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoVarSlaveStruct]: ...

    def to_list(self) -> list[TaoVarSlaveStruct]: ...

class TaoVarSlaveStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoVarSlaveStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoVarSlaveStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoVarSlaveStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoVarSlaveStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoVarSlaveStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoVarSlaveStruct]: ...

class TaoVarStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoVarStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoVarStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoVarStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoVarStruct]: ...

    def to_list(self) -> list[TaoVarStruct]: ...

class TaoVarStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoVarStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoVarStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoVarStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoVarStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoVarStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoVarStruct]: ...

class TaoWaveKickPtStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaoWaveKickPtStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoWaveKickPtStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoWaveKickPtStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoWaveKickPtStruct]: ...

    def to_list(self) -> list[TaoWaveKickPtStruct]: ...

class TaoWaveKickPtStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaoWaveKickPtStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaoWaveKickPtStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaoWaveKickPtStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaoWaveKickPtStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaoWaveKickPtStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaoWaveKickPtStruct]: ...

class TargetPointStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TargetPointStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TargetPointStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TargetPointStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TargetPointStruct]: ...

    def to_list(self) -> list[TargetPointStruct]: ...

class TargetPointStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TargetPointStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TargetPointStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TargetPointStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TargetPointStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TargetPointStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TargetPointStruct]: ...

class TaylorStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaylorStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaylorStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaylorStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaylorStruct]: ...

    def to_list(self) -> list[TaylorStruct]: ...

class TaylorStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaylorStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaylorStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaylorStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaylorStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaylorStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaylorStruct]: ...

class TaylorTermStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TaylorTermStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaylorTermStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaylorTermStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaylorTermStruct]: ...

    def to_list(self) -> list[TaylorTermStruct]: ...

class TaylorTermStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TaylorTermStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TaylorTermStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TaylorTermStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TaylorTermStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TaylorTermStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TaylorTermStruct]: ...

class TestSubStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TestSubStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TestSubStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TestSubStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TestSubStruct]: ...

    def to_list(self) -> list[TestSubStruct]: ...

class TestSubStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TestSubStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TestSubStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TestSubStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TestSubStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TestSubStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TestSubStruct]: ...

class TestSubStructArray2D:
    def __init__(self) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    @overload
    def __getitem__(self, arg: int, /) -> TestSubStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TestSubStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TestSubStruct]: ...

    def __str__(self) -> str: ...

class TestSubStructArray3D:
    def __init__(self) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    @overload
    def __getitem__(self, arg: int, /) -> TestSubStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TestSubStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TestSubStruct]: ...

    def __str__(self) -> str: ...

class TrackPointStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: TrackPointStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> TrackPointStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TrackPointStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TrackPointStruct]: ...

    def to_list(self) -> list[TrackPointStruct]: ...

class TrackPointStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: TrackPointStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[TrackPointStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> TrackPointStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> TrackPointStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TrackPointStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TrackPointStruct]: ...

class TricubicCmplxCoefStructArray3D:
    def __init__(self) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    @overload
    def __getitem__(self, arg: int, /) -> TricubicCmplxCoefStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: TricubicCmplxCoefStruct, /) -> None: ...

    def __iter__(self) -> Iterator[TricubicCmplxCoefStruct]: ...

    def __str__(self) -> str: ...

class VarLengthStringStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: VarLengthStringStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> VarLengthStringStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: VarLengthStringStruct, /) -> None: ...

    def __iter__(self) -> Iterator[VarLengthStringStruct]: ...

    def to_list(self) -> list[VarLengthStringStruct]: ...

class VarLengthStringStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: VarLengthStringStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[VarLengthStringStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> VarLengthStringStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> VarLengthStringStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: VarLengthStringStruct, /) -> None: ...

    def __iter__(self) -> Iterator[VarLengthStringStruct]: ...

class WakeLrModeStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: WakeLrModeStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> WakeLrModeStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: WakeLrModeStruct, /) -> None: ...

    def __iter__(self) -> Iterator[WakeLrModeStruct]: ...

    def to_list(self) -> list[WakeLrModeStruct]: ...

class WakeLrModeStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: WakeLrModeStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[WakeLrModeStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> WakeLrModeStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> WakeLrModeStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: WakeLrModeStruct, /) -> None: ...

    def __iter__(self) -> Iterator[WakeLrModeStruct]: ...

class WakeSrModeStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: WakeSrModeStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> WakeSrModeStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: WakeSrModeStruct, /) -> None: ...

    def __iter__(self) -> Iterator[WakeSrModeStruct]: ...

    def to_list(self) -> list[WakeSrModeStruct]: ...

class WakeSrModeStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: WakeSrModeStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[WakeSrModeStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> WakeSrModeStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> WakeSrModeStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: WakeSrModeStruct, /) -> None: ...

    def __iter__(self) -> Iterator[WakeSrModeStruct]: ...

class Wall3DSectionStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: Wall3DSectionStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> Wall3DSectionStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: Wall3DSectionStruct, /) -> None: ...

    def __iter__(self) -> Iterator[Wall3DSectionStruct]: ...

    def to_list(self) -> list[Wall3DSectionStruct]: ...

class Wall3DSectionStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: Wall3DSectionStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[Wall3DSectionStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> Wall3DSectionStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> Wall3DSectionStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: Wall3DSectionStruct, /) -> None: ...

    def __iter__(self) -> Iterator[Wall3DSectionStruct]: ...

class Wall3DStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: Wall3DStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> Wall3DStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: Wall3DStruct, /) -> None: ...

    def __iter__(self) -> Iterator[Wall3DStruct]: ...

    def to_list(self) -> list[Wall3DStruct]: ...

class Wall3DStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: Wall3DStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[Wall3DStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> Wall3DStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> Wall3DStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: Wall3DStruct, /) -> None: ...

    def __iter__(self) -> Iterator[Wall3DStruct]: ...

class Wall3DVertexStructArray1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: Wall3DVertexStructAlloc1D, /) -> None: ...

    def __len__(self) -> int: ...

    def is_valid(self) -> bool: ...

    def __str__(self) -> str: ...

    @overload
    def __getitem__(self, arg: int, /) -> Wall3DVertexStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: Wall3DVertexStruct, /) -> None: ...

    def __iter__(self) -> Iterator[Wall3DVertexStruct]: ...

    def to_list(self) -> list[Wall3DVertexStruct]: ...

class Wall3DVertexStructAlloc1D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, n: int) -> None: ...

    @overload
    def __init__(self, arg: Wall3DVertexStructArray1D, /) -> None: ...

    @overload
    def __init__(self, arg: Sequence[Wall3DVertexStruct], /) -> None: ...

    def resize(self, n: int) -> None: ...

    def resize_bounds(self, lbound: int, ubound: int) -> None: ...

    def clear(self) -> None: ...

    def __len__(self) -> int: ...

    def view(self) -> Wall3DVertexStructArray1D: ...

    @overload
    def __getitem__(self, arg: int, /) -> Wall3DVertexStruct: ...

    @overload
    def __getitem__(self, arg: slice, /) -> list: ...

    def __setitem__(self, arg0: int, arg1: Wall3DVertexStruct, /) -> None: ...

    def __iter__(self) -> Iterator[Wall3DVertexStruct]: ...

def get_bmad_com() -> BmadCommonStruct:
    """Get the shared BmadCommon structure"""

def get_space_charge_com() -> SpaceChargeCommonStruct:
    """Get the shared SpaceChargeCommon structure"""

def get_super_universe() -> TaoSuperUniverseStruct:
    """Get the shared TaoSuperUniverse structure"""
