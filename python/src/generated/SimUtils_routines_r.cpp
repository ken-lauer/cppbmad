#include "pybmad/generated/SimUtils_routines_r.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_r(nb::module_ &m) {
  m.def(
      "ran_default_state",
      [](RandomStateStruct *set_state) {
        auto fn = static_cast<RandomStateStruct (*)(optional_ref<RandomStateStruct>)>(
            &SimUtils::ran_default_state
        );
        return fn(ptr_to_opt_ref(set_state));
      },
      nb::arg("set_state") = nb::none(),
      R"""(Routine to set or get the state of the default random number generator.
See the ran_seed_put documentation for more details

Parameters
----------
set_state : RandomStateStruct, optional
    State to set the default generator to.

Returns
-------
get_state : RandomStateStruct, optional
    Returns the state of the default generator.
)"""
  );
  m.def(
      "ran_engine",
      [](std::optional<std::string> set,
         std::optional<std::string> get,
         RandomStateStruct *ran_state) {
        auto fn = static_cast<
            void (*)(std::optional<std::string>, std::optional<std::string>, optional_ref<RandomStateStruct>)>(
            &SimUtils::ran_engine
        );
        return fn(set, get, ptr_to_opt_ref(ran_state));
      },
      nb::arg("set") = nb::none(),
      nb::arg("get") = nb::none(),
      nb::arg("ran_state") = nb::none(),
      R"""(Routine to set what random number generator algorithm is used.
If this routine is never called then pseudo_random$ is used.
With sobseq quasi-random numbers the maximum dimension is 6.

Parameters
----------
set : str, optional
    Set the random number engine. Possibilities are: 'pseudo' -> Uses ran from Numerical Recipies (F90).
    'quasi'  -> Uses sobseq from Numerical Recipes. ''       -> Do nothing.

get : str, optional
    Get the current (before any set) random number engine.

ran_state : RandomStateStruct, optional
    Internal state. See the ran_seed_put documentation for more details.
)"""
  );
  nb::class_<SimUtils::RanGaussConverter>(m, "RanGaussConverter", "ran_gauss_converter return type")
      .def_ro("get", &SimUtils::RanGaussConverter::get)
      .def_ro("get_sigma_cut", &SimUtils::RanGaussConverter::get_sigma_cut)
      .def("__len__", [](const SimUtils::RanGaussConverter &) { return 2; })
      .def("__getitem__", [](const SimUtils::RanGaussConverter &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.get);
        if (i == 1)
          return nb::cast(s.get_sigma_cut);
        throw nb::index_error();
      });
  m.def(
      "ran_gauss_converter",
      [](std::optional<std::string> set,
         std::optional<double> set_sigma_cut,
         RandomStateStruct *ran_state) {
        auto fn = static_cast<
            SimUtils::
                RanGaussConverter (*)(std::optional<std::string>, std::optional<double>, optional_ref<RandomStateStruct>)>(
            &SimUtils::ran_gauss_converter
        );
        return fn(set, set_sigma_cut, ptr_to_opt_ref(ran_state));
      },
      nb::arg("set") = nb::none(),
      nb::arg("set_sigma_cut") = nb::none(),
      nb::arg("ran_state") = nb::none(),
      R"""(Routine to set what conversion routine is used for converting
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
set : str, optional
    Set the random number engine. Possibilities are: 'exact' 'quick'  ! Old deprecated: 'limited' 'ziggurat'
    ''       ! Do nothing

set_sigma_cut : float, optional
    Sigma cutoff. Initially: sigma_cut = -1.

ran_state : RandomStateStruct, optional
    Internal state. See the ran_seed_put documentation for more details.

Returns
-------
get : str, optional
    Get the current (before any set) gaussian converter.

get_sigma_cut : float, optional
    Get the current (before any set) sigma cutoff.
)"""
  );
  m.def(
      "ran_gauss_scalar",
      [](RandomStateStruct *ran_state,
         std::optional<double> sigma_cut,
         std::optional<int> index_quasi) {
        auto fn = static_cast<
            double (*)(optional_ref<RandomStateStruct>, std::optional<double>, std::optional<int>)>(
            &SimUtils::ran_gauss_scalar
        );
        return fn(ptr_to_opt_ref(ran_state), sigma_cut, index_quasi);
      },
      nb::arg("ran_state") = nb::none(),
      nb::arg("sigma_cut") = nb::none(),
      nb::arg("index_quasi") = nb::none(),
      R"""(Routine to return a gaussian distributed random number with unit sigma.
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
    Internal state. See the ran_seed_put documentation for more details.

sigma_cut : float, optional
    If present and positive will override setting of ran_state.gauss_sigma_cut.

Returns
-------
harvest : float
    Random number. Or
    As an output, harvest: Random number array.
)"""
  );
  m.def(
      "ran_gauss_vector",
      [](FArray1D<Real> &harvest, RandomStateStruct *ran_state, std::optional<double> sigma_cut) {
        auto fn = static_cast<
            void (*)(FArray1D<Real> &, optional_ref<RandomStateStruct>, std::optional<double>)>(
            &SimUtils::ran_gauss_vector
        );
        return fn(harvest, ptr_to_opt_ref(ran_state), sigma_cut);
      },
      nb::arg("harvest"),
      nb::arg("ran_state") = nb::none(),
      nb::arg("sigma_cut") = nb::none(),
      R"""(Routine to return a gaussian distributed random number with unit sigma.
This routine uses the same algorithm as gasdev from Numerical Recipes.

Note: ran_gauss is an overloaded name for:
    ran_gauss_scalar   ! harvest is a scalar
    ran_gauss_vector   ! harvest is a 1-D array.

Note: Use ran_seed_put for initialization.
Note: Use ran_engine to set which random number generator to use.
Note: Use ran_gauss_converter to set which conversion routine to use.

Parameters
----------
harvest : 1D array of float
    Random number. Or
    As an output, harvest: Random number array.

ran_state : RandomStateStruct, optional
    Internal state. See the ran_seed_put documentation for more details.

sigma_cut : float, optional
    If present and positive will override setting of ran_state.gauss_sigma_cut.
)"""
  );
  m.def(
      "ran_seed_get",
      &SimUtils::ran_seed_get,
      R"""(Routine to return the seed used for the random number generator.

Parameters
----------

Returns
-------
seed : int
    Random number seed used.
)"""
  );
  m.def(
      "ran_seed_put",
      &SimUtils::ran_seed_put,
      nb::arg("seed"),
      nb::arg("mpi_offset") = nb::none(),
      R"""(Routine to seed a random number generator.

If a program never calls ran_seed_put, or ran_seed_put is called with seed = 0,
the system clock will be used to generate the seed.

Note: The seed is only used with the pseudo_random$ engine.
Note: Use the subroutine ran_seed_get(seed) to get the seed used.
Note: Use pointer_to_ran_state() to access the ran state directly.

Parameters
----------
seed : int
    Seed number. If seed = 0 then a seed will be choosen based upon the system clock.

mpi_offset : int, optional
    Offset added to seed. Default is zero. Used with MPI processes ensure different threads use different
    random numbers.
)"""
  );
  m.def(
      "ran_uniform",
      [](RandomStateStruct *ran_state, std::optional<int> index_quasi) {
        auto fn = static_cast<double (*)(optional_ref<RandomStateStruct>, std::optional<int>)>(
            &SimUtils::ran_uniform
        );
        return fn(ptr_to_opt_ref(ran_state), index_quasi);
      },
      nb::arg("ran_state") = nb::none(),
      nb::arg("index_quasi") = nb::none(),
      R"""(Routine to return a random number uniformly distributed in the
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
    Random number. Or
    As an output, harvest: Random number array.
)"""
  );
  m.def(
      "ran_uniform",
      [](FArray1D<Real> &harvest, RandomStateStruct *ran_state) {
        auto fn = static_cast<void (*)(FArray1D<Real> &, optional_ref<RandomStateStruct>)>(
            &SimUtils::ran_uniform
        );
        return fn(harvest, ptr_to_opt_ref(ran_state));
      },
      nb::arg("harvest"),
      nb::arg("ran_state") = nb::none(),
      R"""(Routine to return a random number uniformly distributed in the
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
harvest : 1D array of float
    Random number. Or
    As an output, harvest: Random number array.

ran_state : RandomStateStruct, optional
    Internal state. See the ran_seed_put documentation for more details.
)"""
  );
  m.def(
      "rcelbd",
      &SimUtils::rcelbd,
      nb::arg("mc"),
      nb::arg("elb"),
      nb::arg("eld"),
      R"""(Wrapper for Fortran routine rcelbd

Parameters
----------
mc : float

elb : float

eld : float
)"""
  );
  m.def(
      "read_a_line",
      &SimUtils::read_a_line,
      nb::arg("prompt"),
      nb::arg("trim_prompt") = nb::none(),
      nb::arg("prompt_color") = nb::none(),
      nb::arg("prompt_bold") = nb::none(),
      nb::arg("history_file") = nb::none(),
      R"""(Subroutine to read a line of input from the terminal.
The line is also add to the history buffer so that the up-arrow
and down-arrow keys can be used to recall past commands.

Also see:
  readline_read_history
  readline_write_history

System Libraries that need to be linked to:
  readline curses

Parameters
----------
prompt : str
    Prompt string to use.

trim_prompt : bool, optional
    If present and True then trim the prompt string and add a single blank before printing the prompt string.
    Default is True.

prompt_color : str, optional
    Color of the prompt. Possibilities are: 'BLACK', 'RED', 'GREEN', 'YELLOW', 'BLUE', 'MAGENTA', 'CYAN',
    'GRAY', 'DEFAULT'. The 'DEFAULT' setting (the default) does not set the prompt color.

prompt_bold : bool, optional
    If present and True then the prompt will be printed in bold.

history_file : str, optional
    If present, add line_out to a file whose name is given by history_file. History files are useful for
    saving the command history in between when a program is run multiple times.

Returns
-------
line_out : str
    Line typed by the user. Note: If cntl-D is pressed, line_out = achar(24).
)"""
  );
  m.def(
      "readline_read_history",
      &SimUtils::readline_read_history,
      nb::arg("history_file"),
      R"""(Routine to add the contents of a file to the readline history list.
Use this routine with the read_a_line routine.

Parameters
----------
history_file : str
    Name of the history file. EG: '~/.my_history'

Returns
-------
status : int
    0 = Success, otherwise failure.
)"""
  );
  m.def(
      "readline_write_history",
      &SimUtils::readline_write_history,
      nb::arg("history_file"),
      R"""(Routine to write the contents of the readline history list to a file.
Use this routine with the read_a_line routine.

Parameters
----------
history_file : str
    Name of the history file. EG: '~/.my_history'

Returns
-------
status : int
    0 = Success, otherwise failure.
)"""
  );
  m.def(
      "real_num_fortran_format",
      &SimUtils::real_num_fortran_format,
      nb::arg("number"),
      nb::arg("width"),
      nb::arg("n_blanks") = nb::none(),
      R"""(Wrapper for Fortran routine real_num_fortran_format

Parameters
----------
number : float

width : int

n_blanks : int, optional

Returns
-------
fmt_str : str
)"""
  );
  m.def(
      "real_path",
      &SimUtils::real_path,
      nb::arg("path_in"),
      nb::arg("path_out"),
      R"""(Wrapper for Fortran routine real_path

Parameters
----------
path_in : str

path_out : str

Returns
-------
is_ok : bool
)"""
  );
  m.def(
      "real_str",
      &SimUtils::real_str,
      nb::arg("r_num"),
      nb::arg("n_signif") = nb::none(),
      nb::arg("n_decimal") = nb::none(),
      R"""(Wrapper for Fortran routine real_str

Parameters
----------
r_num : float

n_signif : int, optional

n_decimal : int, optional

Returns
-------
str : str
)"""
  );
  m.def(
      "real_to_string",
      &SimUtils::real_to_string,
      nb::arg("real_num"),
      nb::arg("width"),
      nb::arg("n_signif") = nb::none(),
      nb::arg("n_decimal") = nb::none(),
      R"""(Wrapper for Fortran routine real_to_string

Parameters
----------
real_num : float

width : int

n_signif : int, optional

n_decimal : int, optional

Returns
-------
str : str
)"""
  );
  m.def(
      "reallocate_spline",
      &SimUtils::reallocate_spline,
      nb::arg("spline"),
      nb::arg("n"),
      nb::arg("n_min") = nb::none(),
      nb::arg("exact") = nb::none(),
      R"""(Subroutine to allocate an allocatable spline_struct array.
The data of the array is preserved but data at the end of the
array will be lost if n is less than the original size of the array

Parameters
----------
spline : 1D array of SplineStruct
    Spline to reallocate.
    This parameter is an input/output and is modified in-place.
    As an output, spline: Allocated spline.

n : int
    Upper bound needed for 1-dimensional arrays.

n_min : int, optional
    Lower bound of spline array. Default is 1.

exact : bool, optional
    If present and False then the size of the output array is permitted to be larger than n. Default is True.
)"""
  );
  m.def(
      "relbd",
      &SimUtils::relbd,
      nb::arg("phi"),
      nb::arg("phic"),
      nb::arg("mc"),
      nb::arg("b"),
      nb::arg("d"),
      R"""(Wrapper for Fortran routine relbd

Parameters
----------
phi : float

phic : float

mc : float

b : float

d : float
)"""
  );
  m.def(
      "relcbd",
      &SimUtils::relcbd,
      nb::arg("c0"),
      nb::arg("mc"),
      nb::arg("b"),
      nb::arg("dx"),
      R"""(Wrapper for Fortran routine relcbd

Parameters
----------
c0 : float

mc : float

b : float

dx : float
)"""
  );
  m.def(
      "relsbd",
      &SimUtils::relsbd,
      nb::arg("s0"),
      nb::arg("mc"),
      nb::arg("b"),
      nb::arg("d"),
      R"""(Wrapper for Fortran routine relsbd

Parameters
----------
s0 : float

mc : float

b : float

d : float
)"""
  );
  m.def(
      "rgelbd",
      &SimUtils::rgelbd,
      nb::arg("phi"),
      nb::arg("mc"),
      nb::arg("elb"),
      nb::arg("eld"),
      R"""(Wrapper for Fortran routine rgelbd

Parameters
----------
phi : float

mc : float

elb : float

eld : float
)"""
  );
  nb::class_<SimUtils::RmsValue>(m, "RmsValue", "rms_value return type")
      .def_ro("ave_val", &SimUtils::RmsValue::ave_val)
      .def_ro("rms_val", &SimUtils::RmsValue::rms_val)
      .def("__len__", [](const SimUtils::RmsValue &) { return 2; })
      .def("__getitem__", [](const SimUtils::RmsValue &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.ave_val);
        if (i == 1)
          return nb::cast(s.rms_val);
        throw nb::index_error();
      });
  m.def(
      "rms_value",
      [](FArray1D<Real> &val_arr, BoolAlloc1D *good_val) {
        auto fn = static_cast<SimUtils::RmsValue (*)(FArray1D<Real> &, optional_ref<BoolAlloc1D>)>(
            &SimUtils::rms_value
        );
        return fn(val_arr, ptr_to_opt_ref(good_val));
      },
      nb::arg("val_arr"),
      nb::arg("good_val") = nb::none(),
      R"""(Wrapper for Fortran routine rms_value

Parameters
----------
val_arr : 1D array of float
    Array of reals.

good_val : 1D array of bool, optional
    If present, only calculate RMS where good_val(i) = True.

Returns
-------
rms_val : float
    RMS value. Set to real_garbage$ if there is a problem.

ave_val : float, optional
    average value.
)"""
  );
  m.def(
      "rot_2d",
      &SimUtils::rot_2d,
      nb::arg("vec_in"),
      nb::arg("angle"),
      R"""(Wrapper for Fortran routine rot_2d

Parameters
----------
vec_in : 1D array of float (shape: 2)
    Init vec

angle : float
    angle in radians.

Returns
-------
vec_out : 1D array of float (shape: 2)
    Rotated vec.
)"""
  );
  m.def(
      "rotate_vec",
      &SimUtils::rotate_vec,
      nb::arg("vec"),
      nb::arg("axis"),
      nb::arg("angle"),
      R"""(Basic routine to rotate vector components around the x, y, or z axis.

Parameters
----------
vec : 1D array of float
    vector
    This parameter is an input/output and is modified in-place.
    As an output, vec: Rotated vector.

axis : int
    x_axis$, y_axis$, or z_axis$

angle : float
    angle to rotate.
)"""
  );
  m.def(
      "rotate_vec_given_axis_angle",
      &SimUtils::rotate_vec_given_axis_angle,
      nb::arg("vec_in"),
      nb::arg("axis"),
      nb::arg("angle"),
      R"""(Routine to rotate a vector.

Parameters
----------
vec_in : 1D array of float (shape: 3)
    Initial vector.

axis : 1D array of float
    Axis of rotation. Must be normalized to 1.

angle : float
    Angle to rotate by

Returns
-------
vec_out : 1D array of float (shape: 3)
    Final vector.
)"""
  );
  m.def(
      "rp8",
      &SimUtils::rp8,
      nb::arg("int_in"),
      R"""(Routine to convert from integer to real of type rp.
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
      "rserbd",
      &SimUtils::rserbd,
      nb::arg("y"),
      nb::arg("m"),
      nb::arg("b"),
      nb::arg("d"),
      R"""(Wrapper for Fortran routine rserbd

Parameters
----------
y : float

m : float

b : float

d : float
)"""
  );
  m.def(
      "run_timer",
      &SimUtils::run_timer,
      nb::arg("command"),
      nb::arg("time") = nb::none(),
      nb::arg("time0") = nb::none(),
      R"""(Wrapper for Fortran routine run_timer

Parameters
----------
command : str

time : float, optional

time0 : float, optional
)"""
  );
}
