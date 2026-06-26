"""CppBmadTest routines"""

import _pybmad


class TestBunchStructArray:
    """test_bunch_struct_array return type"""

    @property
    def arr_out(self) -> _pybmad.BunchStructAlloc1D: ...

    @property
    def opt_status(self) -> list[int]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def test_bunch_struct_array(arr_in: _pybmad.BunchStructArray1D, arr_inout: _pybmad.BunchStructArray1D, arr_in_opt: _pybmad.BunchStructArray1D | None = None, arr_inout_opt: _pybmad.BunchStructArray1D | None = None) -> TestBunchStructArray:
    """
    Wrapper for Fortran routine test_bunch_struct_array

    Parameters
    ----------
    arr_in : 1D array of BunchStruct

    arr_inout : 1D array of BunchStruct

    arr_in_opt : 1D array of BunchStruct, optional

    arr_inout_opt : 1D array of BunchStruct, optional

    Returns
    -------
    arr_out : 1D array of BunchStruct

    opt_status : 1D array of int (shape: 2)
    """

class TestBunchStructScalar:
    """test_bunch_struct_scalar return type"""

    @property
    def val_out(self) -> _pybmad.BunchStruct: ...

    @property
    def opt_status(self) -> list[int]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def test_bunch_struct_scalar(val_in: _pybmad.BunchStruct, val_inout: _pybmad.BunchStruct, val_in_opt: _pybmad.BunchStruct | None = None, val_inout_opt: _pybmad.BunchStruct | None = None) -> TestBunchStructScalar:
    """
    Wrapper for Fortran routine test_bunch_struct_scalar

    Parameters
    ----------
    val_in : BunchStruct

    val_inout : BunchStruct

    val_in_opt : BunchStruct, optional

    val_inout_opt : BunchStruct, optional

    Returns
    -------
    val_out : BunchStruct

    opt_status : 1D array of int (shape: 2)
    """

class TestCharacterArray:
    """test_character_array return type"""

    @property
    def arr_out(self) -> _pybmad.CharacterAlloc1D: ...

    @property
    def opt_status(self) -> list[int]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def test_character_array(arr_in: _pybmad.CharacterAlloc1D, arr_inout: _pybmad.CharacterAlloc1D, arr_in_opt: _pybmad.CharacterAlloc1D | None = None, arr_inout_opt: _pybmad.CharacterAlloc1D | None = None) -> TestCharacterArray:
    """
    Wrapper for Fortran routine test_character_array

    Parameters
    ----------
    arr_in : 1D array of str

    arr_inout : 1D array of str

    arr_in_opt : 1D array of str, optional

    arr_inout_opt : 1D array of str, optional

    Returns
    -------
    arr_out : 1D array of str

    opt_status : 1D array of int (shape: 2)
    """

class TestCharacterScalar:
    """test_character_scalar return type"""

    @property
    def val_out(self) -> str: ...

    @property
    def opt_status(self) -> list[int]: ...

    @property
    def val_inout(self) -> str: ...

    @property
    def val_inout_opt(self) -> str | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def test_character_scalar(val_in: str, val_inout: str, val_in_opt: str | None = None, val_inout_opt: str | None = None) -> TestCharacterScalar:
    """
    Wrapper for Fortran routine test_character_scalar

    Parameters
    ----------
    val_in : str

    val_inout : str

    val_in_opt : str, optional

    val_inout_opt : str, optional

    Returns
    -------
    val_inout : str

    val_out : str

    opt_status : 1D array of int (shape: 2)

    val_inout_opt : str, optional
    """

class TestComplexArray:
    """test_complex_array return type"""

    @property
    def arr_out(self) -> _pybmad.ComplexAlloc1D: ...

    @property
    def opt_status(self) -> list[int]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def test_complex_array(arr_in: _pybmad.ComplexArray1D, arr_inout: _pybmad.ComplexArray1D, arr_in_opt: _pybmad.ComplexArray1D | None = None, arr_inout_opt: _pybmad.ComplexArray1D | None = None) -> TestComplexArray:
    """
    Wrapper for Fortran routine test_complex_array

    Parameters
    ----------
    arr_in : 1D array of complex

    arr_inout : 1D array of complex

    arr_in_opt : 1D array of complex, optional

    arr_inout_opt : 1D array of complex, optional

    Returns
    -------
    arr_out : 1D array of complex

    opt_status : 1D array of int (shape: 2)
    """

class TestComplexScalar:
    """test_complex_scalar return type"""

    @property
    def val_out(self) -> complex: ...

    @property
    def opt_status(self) -> list[int]: ...

    @property
    def val_inout(self) -> complex: ...

    @property
    def val_inout_opt(self) -> complex | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def test_complex_scalar(val_in: complex, val_inout: complex, val_in_opt: complex | None = None, val_inout_opt: complex | None = None) -> TestComplexScalar:
    """
    Wrapper for Fortran routine test_complex_scalar

    Parameters
    ----------
    val_in : complex

    val_inout : complex

    val_in_opt : complex, optional

    val_inout_opt : complex, optional

    Returns
    -------
    val_inout : complex

    val_out : complex

    opt_status : 1D array of int (shape: 2)

    val_inout_opt : complex, optional
    """

class TestInteger8Array:
    """test_integer8_array return type"""

    @property
    def arr_out(self) -> _pybmad.Int8Alloc1D: ...

    @property
    def opt_status(self) -> list[int]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def test_integer8_array(arr_in: _pybmad.Int8Array1D, arr_inout: _pybmad.Int8Array1D, arr_in_opt: _pybmad.Int8Array1D | None = None, arr_inout_opt: _pybmad.Int8Array1D | None = None) -> TestInteger8Array:
    """
    Wrapper for Fortran routine test_integer8_array

    Parameters
    ----------
    arr_in : 1D array of int

    arr_inout : 1D array of int

    arr_in_opt : 1D array of int, optional

    arr_inout_opt : 1D array of int, optional

    Returns
    -------
    arr_out : 1D array of int

    opt_status : 1D array of int (shape: 2)
    """

class TestInteger8Scalar:
    """test_integer8_scalar return type"""

    @property
    def val_out(self) -> int: ...

    @property
    def opt_status(self) -> list[int]: ...

    @property
    def val_inout(self) -> int: ...

    @property
    def val_inout_opt(self) -> int | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def test_integer8_scalar(val_in: int, val_inout: int, val_in_opt: int | None = None, val_inout_opt: int | None = None) -> TestInteger8Scalar:
    """
    Wrapper for Fortran routine test_integer8_scalar

    Parameters
    ----------
    val_in : int

    val_inout : int

    val_in_opt : int, optional

    val_inout_opt : int, optional

    Returns
    -------
    val_inout : int

    val_out : int

    opt_status : 1D array of int (shape: 2)

    val_inout_opt : int, optional
    """

class TestIntegerArray:
    """test_integer_array return type"""

    @property
    def arr_out(self) -> _pybmad.IntAlloc1D: ...

    @property
    def opt_status(self) -> list[int]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def test_integer_array(arr_in: _pybmad.IntArray1D, arr_inout: _pybmad.IntArray1D, arr_in_opt: _pybmad.IntArray1D | None = None, arr_inout_opt: _pybmad.IntArray1D | None = None) -> TestIntegerArray:
    """
    Wrapper for Fortran routine test_integer_array

    Parameters
    ----------
    arr_in : 1D array of int

    arr_inout : 1D array of int

    arr_in_opt : 1D array of int, optional

    arr_inout_opt : 1D array of int, optional

    Returns
    -------
    arr_out : 1D array of int

    opt_status : 1D array of int (shape: 2)
    """

class TestIntegerScalar:
    """test_integer_scalar return type"""

    @property
    def val_out(self) -> int: ...

    @property
    def opt_status(self) -> list[int]: ...

    @property
    def val_inout(self) -> int: ...

    @property
    def val_inout_opt(self) -> int | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def test_integer_scalar(val_in: int, val_inout: int, val_in_opt: int | None = None, val_inout_opt: int | None = None) -> TestIntegerScalar:
    """
    Wrapper for Fortran routine test_integer_scalar

    Parameters
    ----------
    val_in : int

    val_inout : int

    val_in_opt : int, optional

    val_inout_opt : int, optional

    Returns
    -------
    val_inout : int

    val_out : int

    opt_status : 1D array of int (shape: 2)

    val_inout_opt : int, optional
    """

class TestLogicalArray:
    """test_logical_array return type"""

    @property
    def arr_out(self) -> _pybmad.BoolAlloc1D: ...

    @property
    def opt_status(self) -> list[int]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def test_logical_array(arr_in: _pybmad.BoolAlloc1D, arr_inout: _pybmad.BoolAlloc1D, arr_in_opt: _pybmad.BoolAlloc1D | None = None, arr_inout_opt: _pybmad.BoolAlloc1D | None = None) -> TestLogicalArray:
    """
    Wrapper for Fortran routine test_logical_array

    Parameters
    ----------
    arr_in : 1D array of bool

    arr_inout : 1D array of bool

    arr_in_opt : 1D array of bool, optional

    arr_inout_opt : 1D array of bool, optional

    Returns
    -------
    arr_out : 1D array of bool

    opt_status : 1D array of int (shape: 2)
    """

class TestLogicalScalar:
    """test_logical_scalar return type"""

    @property
    def val_out(self) -> bool: ...

    @property
    def opt_status(self) -> list[int]: ...

    @property
    def val_inout(self) -> bool: ...

    @property
    def val_inout_opt(self) -> bool | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def test_logical_scalar(val_in: bool, val_inout: bool, val_in_opt: bool | None = None, val_inout_opt: bool | None = None) -> TestLogicalScalar:
    """
    Wrapper for Fortran routine test_logical_scalar

    Parameters
    ----------
    val_in : bool

    val_inout : bool

    val_in_opt : bool, optional

    val_inout_opt : bool, optional

    Returns
    -------
    val_inout : bool

    val_out : bool

    opt_status : 1D array of int (shape: 2)

    val_inout_opt : bool, optional
    """

class TestReal16Array:
    """test_real16_array return type"""

    @property
    def arr_out(self) -> _pybmad.Real16Alloc1D: ...

    @property
    def opt_status(self) -> list[int]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def test_real16_array(arr_in: _pybmad.Real16Alloc1D, arr_inout: _pybmad.Real16Alloc1D, arr_in_opt: _pybmad.Real16Alloc1D | None = None, arr_inout_opt: _pybmad.Real16Alloc1D | None = None) -> TestReal16Array:
    """
    Wrapper for Fortran routine test_real16_array

    Parameters
    ----------
    arr_in : 1D array of float

    arr_inout : 1D array of float

    arr_in_opt : 1D array of float, optional

    arr_inout_opt : 1D array of float, optional

    Returns
    -------
    arr_out : 1D array of float

    opt_status : 1D array of int (shape: 2)
    """

class TestReal16Scalar:
    """test_real16_scalar return type"""

    @property
    def val_out(self) -> float: ...

    @property
    def opt_status(self) -> list[int]: ...

    @property
    def val_inout(self) -> float: ...

    @property
    def val_inout_opt(self) -> float | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def test_real16_scalar(val_in: float, val_inout: float, val_in_opt: float | None = None, val_inout_opt: float | None = None) -> TestReal16Scalar:
    """
    Wrapper for Fortran routine test_real16_scalar

    Parameters
    ----------
    val_in : float

    val_inout : float

    val_in_opt : float, optional

    val_inout_opt : float, optional

    Returns
    -------
    val_inout : float

    val_out : float

    opt_status : 1D array of int (shape: 2)

    val_inout_opt : float, optional
    """

class TestRealArray:
    """test_real_array return type"""

    @property
    def arr_out(self) -> _pybmad.RealAlloc1D: ...

    @property
    def opt_status(self) -> list[int]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def test_real_array(arr_in: _pybmad.RealArray1D, arr_inout: _pybmad.RealArray1D, arr_in_opt: _pybmad.RealArray1D | None = None, arr_inout_opt: _pybmad.RealArray1D | None = None) -> TestRealArray:
    """
    Wrapper for Fortran routine test_real_array

    Parameters
    ----------
    arr_in : 1D array of float

    arr_inout : 1D array of float

    arr_in_opt : 1D array of float, optional

    arr_inout_opt : 1D array of float, optional

    Returns
    -------
    arr_out : 1D array of float

    opt_status : 1D array of int (shape: 2)
    """

class TestRealScalar:
    """test_real_scalar return type"""

    @property
    def val_out(self) -> float: ...

    @property
    def opt_status(self) -> list[int]: ...

    @property
    def val_inout(self) -> float: ...

    @property
    def val_inout_opt(self) -> float | None: ...

    def __len__(self) -> int: ...

    def __getitem__(self, arg: int, /) -> object: ...

def test_real_scalar(val_in: float, val_inout: float, val_in_opt: float | None = None, val_inout_opt: float | None = None) -> TestRealScalar:
    """
    Wrapper for Fortran routine test_real_scalar

    Parameters
    ----------
    val_in : float

    val_inout : float

    val_in_opt : float, optional

    val_inout_opt : float, optional

    Returns
    -------
    val_inout : float

    val_out : float

    opt_status : 1D array of int (shape: 2)

    val_inout_opt : float, optional
    """
