#include "pybmad/generated/SimUtils_routines_r.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_r(py::module &m) {
  m.def(
      "ran_default_state",
      &SimUtils::ran_default_state,
      py::arg("set_state") = py::none(),
      R"""(Subroutine ran_default_state (set_state, get_state)

Routine to set or get the state of the default random number generator.
See the ran_seed_put documentation for more details

Parameters
----------
set_state : RandomStateStruct, optional
    State to set the default generator to.

Returns
-------
get_state : RandomStateStruct
    Returns the state of the default generator.
)"""
  );
  m.def(
      "ran_engine",
      &SimUtils::ran_engine,
      py::arg("set") = py::none(),
      py::arg("get") = py::none(),
      py::arg("ran_state") = py::none(),
      R"""(Subroutine ran_engine (set, get, ran_state)

Routine to set what random number generator algorithm is used.
If this routine is never called then pseudo_random$ is used.
With sobseq quasi-random numbers the maximum dimension is 6.

Parameters
----------
set : unknown, optional
    Set the random number engine. Possibilities are: 'pseudo' -> Uses ran from Numerical Recipies (F90).
    'quasi'  -> Uses sobseq from Numerical Recipes.
''       -> Do nothing. : 
get : unknown, optional
    Get the current (before any set) random number engine.
ran_state : RandomStateStruct, optional
    Internal state. See the ran_seed_put documentation for more details.
)"""
  );
  py::class_<SimUtils::RanGaussConverter, std::unique_ptr<SimUtils::RanGaussConverter>>(
      m,
      "RanGaussConverter",
      "ran_gauss_converter return type"
  )
      .def_readonly("get", &SimUtils::RanGaussConverter::get)
      .def_readonly("get_sigma_cut", &SimUtils::RanGaussConverter::get_sigma_cut)
      .def("__len__", [](const SimUtils::RanGaussConverter &) { return 2; })
      .def("__getitem__", [](const SimUtils::RanGaussConverter &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.get);
        if (i == 1)
          return py::cast(s.get_sigma_cut);
        throw py::index_error();
      });
  m.def(
      "ran_gauss_converter",
      &SimUtils::ran_gauss_converter,
      py::arg("set") = py::none(),
      py::arg("set_sigma_cut") = py::none(),
      py::arg("ran_state") = py::none(),
      R"""(Subroutine ran_gauss_converter (set, set_sigma_cut, get, get_sigma_cut, ran_state)

Routine to set what conversion routine is used for converting
uniformly distributed random numbers to Gaussian distributed random numbers.

If this routine is not called then exact_gaussian$ is used.

exact_gaussian$ is a straight forward converter as explained in Numerical recipes.

quick_gaussian$ is a quick a dirty approximation with a cutoff so that no
numbers will be generated beyound what is set for sigma_cut.

A negative sigma_cut means that the exact_gaussian$ will not be limited
and the quick_gaussian$ will use a default of 10.0

Note: Because of technical issues, when using the quasi_random$ number generator
(see the ran_engine routine), the quick_gaussian$ method will automatically be
used independent of what was set with this routine.

Parameters
----------
set : unknown, optional
    Set the random number engine. Possibilities are: 'exact' 'quick'  ! Old deprecated: 'limited'
''       ! Do nothing : 
set_sigma_cut : float, optional
    Sigma cutoff. Initially: sigma_cut = -1.
ran_state : RandomStateStruct, optional
    Internal state. See the ran_seed_put documentation for more details.

Returns
-------
get : unknown
    Get the current (before any set) gaussian converter.
get_sigma_cut : float
    Get the current (before any set) sigma cutoff.
)"""
  );
  m.def(
      "ran_gauss_scalar",
      &SimUtils::ran_gauss_scalar,
      py::arg("ran_state") = py::none(),
      py::arg("sigma_cut") = py::none(),
      py::arg("index_quasi") = py::none(),
      R"""(Subroutine ran_gauss (harvest, ran_state, sigma_cut)

Routine to return a gaussian distributed random number with unit sigma.
This routine uses the same algorithm as gasdev from Numerical Recipes.

Note: ran_gauss is an overloaded name for:
    ran_gauss_scalar   ! harvest is a scalar
    ran_gauss_vector   ! harvest is a 1-D array.

Note: Use ran_seed_put for initialization.
Note: Use ran_engine to set which random number generator to use.
Note: Use ran_gauss_converter to set which conversion routine to use.

Parameters
----------
ran_state : RandomStateStruct, optional
    Internal state.
See the ran_seed_put documentation for more details. : 
sigma_cut : float, optional
    If present and positive will override setting of ran_state.gauss_sigma_cut.

Returns
-------
harvest : float
    Random number.
    This parameter is an input/output and is modified in-place. As an output: Random number array.
Or : 

Notes
-----
Overloaded versions:
)"""
  );
  m.def(
      "ran_gauss_vector",
      &SimUtils::ran_gauss_vector,
      py::arg("harvest"),
      py::arg("ran_state") = py::none(),
      py::arg("sigma_cut") = py::none(),
      R"""(Subroutine ran_gauss (harvest, ran_state, sigma_cut)

Routine to return a gaussian distributed random number with unit sigma.
This routine uses the same algorithm as gasdev from Numerical Recipes.

Note: ran_gauss is an overloaded name for:
    ran_gauss_scalar   ! harvest is a scalar
    ran_gauss_vector   ! harvest is a 1-D array.

Note: Use ran_seed_put for initialization.
Note: Use ran_engine to set which random number generator to use.
Note: Use ran_gauss_converter to set which conversion routine to use.

Parameters
----------
ran_state : RandomStateStruct, optional
    Internal state.
See the ran_seed_put documentation for more details. : 
sigma_cut : float, optional
    If present and positive will override setting of ran_state.gauss_sigma_cut.

Returns
-------
harvest : float
    Random number.
    This parameter is an input/output and is modified in-place. As an output: Random number array.
Or : 

Notes
-----
Overloaded versions:
)"""
  );
  m.def(
      "ran_seed_get",
      &SimUtils::ran_seed_get,
      R"""(Subroutine ran_seed_get (seed)

Routine to return the seed used for the random number generator.

Parameters
----------
ran_state : RandomStateStruct, optional
    Internal state. See the ran_seed_put documentation for more details.

Returns
-------
seed : int
    Random number seed used.
)"""
  );
  m.def(
      "ran_seed_put",
      &SimUtils::ran_seed_put,
      py::arg("seed"),
      py::arg("mpi_offset") = py::none(),
      R"""(Subroutine ran_seed_put (seed, mpi_offset)

Routine to seed a random number generator.

If a program never calls ran_seed_put, or ran_seed_put is called with seed = 0,
the system clock will be used to generate the seed.

Note: The seed is only used with the pseudo_random$ engine.
Note: Use the subroutine ran_seed_get(seed) to get the seed used.
Note: Use pointer_to_ran_state() to access the ran state directly.

Parameters
----------
seed : int
    Seed number. If seed = 0 then a
seed will be choosen based upon the system clock. : 
mpi_offset : int, optional
    Offset added to seed. Default is zero. Used with MPI processes ensure different threads use different
    random numbers.
)"""
  );
  m.def(
      "ran_uniform",
      py::overload_cast<optional_ref<RandomStateStruct>, std::optional<int>>(&SimUtils::ran_uniform
      ),
      py::arg("ran_state") = py::none(),
      py::arg("index_quasi") = py::none(),
      R"""(Subroutine ran_uniform (harvest, ran_state)

Routine to return a random number uniformly distributed in the
interval [0, 1]. This routine uses the same algorithm as ran or sobseq
from Numberical Recipes in Fortran90.
See ran_engine.

Note: ran_uniform is an overloaded name for:
    ran_uniform_scalar   ! harvest is a scalar
    ran_uniform_vector   ! harvest is a 1-D array.

Note: Use ran_seed_put for initialization.
Note: Use ran_engine to set which random number generator to use.

Parameters
----------
ran_state : RandomStateStruct, optional
    Internal state. See the ran_seed_put documentation for more details.

Returns
-------
harvest : float
    Random number.
    This parameter is an input/output and is modified in-place. As an output: Random number array.
Or : 

Notes
-----
Overloaded versions:
)"""
  );
  m.def(
      "ran_uniform",
      py::overload_cast<FArray1D<Real> &, optional_ref<RandomStateStruct>>(&SimUtils::ran_uniform),
      py::arg("harvest"),
      py::arg("ran_state") = py::none(),
      R"""(Subroutine ran_uniform (harvest, ran_state)

Routine to return a random number uniformly distributed in the
interval [0, 1]. This routine uses the same algorithm as ran or sobseq
from Numberical Recipes in Fortran90.
See ran_engine.

Note: ran_uniform is an overloaded name for:
    ran_uniform_scalar   ! harvest is a scalar
    ran_uniform_vector   ! harvest is a 1-D array.

Note: Use ran_seed_put for initialization.
Note: Use ran_engine to set which random number generator to use.

Parameters
----------
ran_state : RandomStateStruct, optional
    Internal state. See the ran_seed_put documentation for more details.

Returns
-------
harvest : float
    Random number.
    This parameter is an input/output and is modified in-place. As an output: Random number array.
Or : 

Notes
-----
Overloaded versions:
)"""
  );
  m.def(
      "real_num_fortran_format",
      &SimUtils::real_num_fortran_format,
      py::arg("number"),
      py::arg("width"),
      py::arg("n_blanks") = py::none(),
      py::arg("fmt_str"),
      R"""(Parameters
----------
number : 
width : 
n_blanks : 
fmt_str : 
)"""
  );
  m.def(
      "real_path",
      &SimUtils::real_path,
      py::arg("path_in"),
      py::arg("path_out"),
      py::arg("is_ok"),
      R"""(Parameters
----------
path_in : 
path_out : 
is_ok : 
)"""
  );
  m.def(
      "real_str",
      &SimUtils::real_str,
      py::arg("r_num"),
      py::arg("n_signif") = py::none(),
      py::arg("n_decimal") = py::none(),
      py::arg("str"),
      R"""(Parameters
----------
r_num : 
n_signif : 
n_decimal : 
str : 
)"""
  );
  m.def(
      "real_to_string",
      &SimUtils::real_to_string,
      py::arg("real_num"),
      py::arg("width"),
      py::arg("n_signif") = py::none(),
      py::arg("n_decimal") = py::none(),
      py::arg("str"),
      R"""(Parameters
----------
real_num : 
width : 
n_signif : 
n_decimal : 
str : 
)"""
  );
  m.def(
      "reallocate_spline",
      &SimUtils::reallocate_spline,
      py::arg("spline"),
      py::arg("n"),
      py::arg("n_min") = py::none(),
      py::arg("exact") = py::none(),
      R"""(Subroutine reallocate_spline (spline, n, n_min, exact)

Subroutine to allocate an allocatable spline_struct array.
The data of the array is preserved but data at the end of the
array will be lost if n is less than the original size of the array


Parameters
----------
spline : SplineStruct
    Spline to reallocate.
    This parameter is an input/output and is modified in-place. As an output: Allocated spline.
n : int
    Upper bound needed for 1-dimensional arrays.
n_min : int, optional
    Lower bound of spline array. Default is 1.
exact : bool, optional
    If present and False then the size of the output array is permitted to be larger than n. Default is True.
)"""
  );
  py::class_<SimUtils::RmsValue, std::unique_ptr<SimUtils::RmsValue>>(
      m,
      "RmsValue",
      "rms_value return type"
  )
      .def_readonly("ave_val", &SimUtils::RmsValue::ave_val)
      .def_readonly("rms_val", &SimUtils::RmsValue::rms_val)
      .def("__len__", [](const SimUtils::RmsValue &) { return 2; })
      .def("__getitem__", [](const SimUtils::RmsValue &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.ave_val);
        if (i == 1)
          return py::cast(s.rms_val);
        throw py::index_error();
      });
  m.def(
      "rms_value",
      &SimUtils::rms_value,
      py::arg("val_arr"),
      py::arg("good_val") = py::none(),
      R"""(Parameters
----------
val_arr : float
    Array of reals.
good_val : bool, optional
    If present, only calculate RMS where good_val(i) = True.
ave_val : float
    average value.
rms_val : float
    RMS value. Set to real_garbage$ if there is a problem.
)"""
  );
  m.def(
      "rot_2d",
      &SimUtils::rot_2d,
      py::arg("vec_in"),
      py::arg("angle"),
      R"""(Parameters
----------
vec_in : float
    Init vec
angle : float
    angle in radians.
vec_out : float
    Rotated vec.
)"""
  );
  m.def(
      "rotate_vec",
      &SimUtils::rotate_vec,
      py::arg("vec"),
      py::arg("axis"),
      py::arg("angle"),
      R"""(Subroutine rotate_vec (vec, axis, angle)

Basic routine to rotate vector components around the x, y, or z axis.

Parameters
----------
vec : float
    vector
    This parameter is an input/output and is modified in-place. As an output: Rotated vector.
axis : int
    x_axis$, y_axis$, or z_axis$
angle : float
    angle to rotate.
)"""
  );
  m.def(
      "rotate_vec_given_axis_angle",
      &SimUtils::rotate_vec_given_axis_angle,
      py::arg("vec_in"),
      py::arg("axis"),
      py::arg("angle"),
      R"""(Function rotate_vec_given_axis_angle (vec_in, axis, angle) result (vec_out)

Routine to rotate a vector.

Parameters
----------
vec_in : float
    Initial vector.
axis : float
    Axis of rotation. Must be normalized to 1.
angle : float
    Angle to rotate by

Returns
-------
vec_out : float
    Final vector.
)"""
  );
  m.def(
      "rp8",
      &SimUtils::rp8,
      py::arg("int_in"),
      R"""(Function rp8(int_in) result (re_out)

Routine to convert from integer to real of type rp.
This routine is used to avoid the implicit integer to single precision that happens when
multiplying int*real(rp).

Parameters
----------
int_in : int
    Input integer.

Returns
-------
re_out : float
    Equiv real.
)"""
  );
  m.def(
      "run_timer",
      &SimUtils::run_timer,
      py::arg("command"),
      py::arg("time") = py::none(),
      py::arg("time0") = py::none(),
      R"""(Parameters
----------
command : 
time : 
time0 : 
)"""
  );
}
