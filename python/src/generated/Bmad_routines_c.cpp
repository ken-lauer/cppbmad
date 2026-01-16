#include "pybmad/generated/Bmad_routines_c.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Bmad_routines_c(py::module &m) {
  m.def(
      "c_to_cbar",
      &Bmad::c_to_cbar,
      py::arg("ele"),
      R"""(Parameters
----------
ele : EleStruct
    Element with C matrix and Twiss parameters.
cbar_mat : float
    Cbar matrix.
)"""
  );
  m.def(
      "calc_bunch_params",
      &Bmad::calc_bunch_params,
      py::arg("bunch"),
      py::arg("bunch_params"),
      py::arg("error"),
      py::arg("print_err") = py::none(),
      py::arg("n_mat") = py::none(),
      py::arg("is_time_coords") = py::none(),
      py::arg("ele") = py::none(),
      R"""(Subroutine calc_bunch_params (bunch, bunch_params, error, print_err, n_mat, is_time_coords, ele)

Finds all bunch parameters defined in bunch_params_struct, both normal-mode
and projected. Projected parameters are found purely from the geometrical
distribution of the beam. Normal-Mode parameters are found using the method
developed in:
  "Alternate approach to general coupled linear optics"
   A. Wolski, PRST AB 9, 024001 (2006)

Note: If less than two particle remain then the various parameters will be
set to zero.

Parameters
----------
bunch : 
    Bunch_struct
print_err : bool, optional
    If present and False then suppress
"no eigen-system found" messages. : 
is_time_coords : bool, optional
    Are particle coords using time coords. Default is False.
ele : EleStruct, optional
    Element being tracked through. Must be present if is_time_coords = True.
Output : 
bunch_params : BunchParamsStruct
error : bool
    Set True if there is an error.
n_mat : float, optional
    N matrix defined in Wolski Eq 44 and used to convert from action-angle coords to lab coords (Wolski Eq
    51.).
)"""
  );
  m.def(
      "calc_bunch_params_slice",
      &Bmad::calc_bunch_params_slice,
      py::arg("bunch"),
      py::arg("bunch_params"),
      py::arg("plane"),
      py::arg("slice_center"),
      py::arg("slice_spread"),
      py::arg("err"),
      py::arg("print_err") = py::none(),
      py::arg("is_time_coords") = py::none(),
      py::arg("ele") = py::none(),
      R"""(subroutine calc_bunch_params_slice (bunch, bunch_params, plane, slice_center, slice_spread, err, print_err, is_time_coords, ele)

Finds bunch parameters for a slice of the beam.

Parameters
----------
bunch : 
    bunch_struct
plane : int
    plane to slice through (x$, px$, & etc...)
slice_center : float
    Center to take slice about
slice_spread : float
    +/- spread in slice about center.
print_err : bool, optional
    If present and False then suppress
"no eigen-system found" messages. : 
is_time_coords : bool, optional
    Default is False. If True, input bunch is using time coordinates in which
case there will be a conversion to s-coords before bunch_params are computed. : 
ele : EleStruct, optional
    Element being tracked through. Must be present if is_time_coords = True.
Output : 
params : BunchParamsStruct
err : bool
    Set True if there is an error.
)"""
  );
  m.def(
      "calc_bunch_params_z_slice",
      &Bmad::calc_bunch_params_z_slice,
      py::arg("bunch"),
      py::arg("bunch_params"),
      py::arg("slice_bounds"),
      py::arg("err"),
      py::arg("print_err") = py::none(),
      py::arg("is_time_coords") = py::none(),
      py::arg("ele") = py::none(),
      R"""(subroutine calc_bunch_params_z_slice (bunch, bunch_params, slice_bounds, err, print_err, is_time_coords, ele)

Finds bunch parameters for a slice of the beam.

The slice is specified in terms of percentage of particles ordered by z-position.
For example, slice_bounds = [0.0, 0.5] specifies the trailing half of the bunch

Parameters
----------
bunch : 
    bunch_struct
slice_bounds : float
    Slice bounds in percentage of particles ordered by z-position.
0.0 is the back of the bunch and 1.0 is the front of the bunch. : 
print_err : bool, optional
    If present and False then suppress
"no eigen-system found" messages. : 
is_time_coords : bool, optional
    Default is False. If True, input bunch is using time coordinates in which
case there will be a conversion to s-coords before bunch_params are computed. : 
ele : EleStruct, optional
    Element being tracked through. Must be present if is_time_coords = True.
Output : 
params : BunchParamsStruct
err : bool
    Set True if there is an error.
)"""
  );
  m.def(
      "calc_bunch_sigma_matrix_etc",
      &Bmad::calc_bunch_sigma_matrix_etc,
      py::arg("particle"),
      py::arg("charge"),
      py::arg("is_time_coords") = py::none(),
      py::arg("ele") = py::none(),
      R"""(Subroutine calc_bunch_sigma_matrix_etc (particle, charge, bunch_params, is_time_coords, ele)

Routine to find the sigma matrix elements of a particle distribution.

Parameters
----------
particle : CoordStruct
    Array of particles.
charge : float
    Particle charge or photon intensity.

Returns
-------
bunch_params : BunchParamsStruct
    Bunch parameters. .sigma(6,6) .centroid.vec(6) .centroid.p0c .rel_max(6) .rel_min(6)
)"""
  );
  py::class_<
      Bmad::CalcEmittancesAndTwissFromSigmaMatrix,
      std::unique_ptr<Bmad::CalcEmittancesAndTwissFromSigmaMatrix>>(
      m,
      "CalcEmittancesAndTwissFromSigmaMatrix",
      "calc_emittances_and_twiss_from_sigma_matrix return type"
  )
      .def_readonly("bunch_params", &Bmad::CalcEmittancesAndTwissFromSigmaMatrix::bunch_params)
      .def_readonly("error", &Bmad::CalcEmittancesAndTwissFromSigmaMatrix::error)
      .def_readonly("n_mat", &Bmad::CalcEmittancesAndTwissFromSigmaMatrix::n_mat)
      .def("__len__", [](const Bmad::CalcEmittancesAndTwissFromSigmaMatrix &) { return 3; })
      .def(
          "__getitem__",
          [](const Bmad::CalcEmittancesAndTwissFromSigmaMatrix &s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.bunch_params);
            if (i == 1)
              return py::cast(s.error);
            if (i == 2)
              return py::cast(s.n_mat);
            throw py::index_error();
          }
      );
  m.def(
      "calc_emittances_and_twiss_from_sigma_matrix",
      &Bmad::calc_emittances_and_twiss_from_sigma_matrix,
      py::arg("sigma_mat"),
      py::arg("print_err") = py::none(),
      R"""(Subroutine calc_emittances_and_twiss_from_sigma_matrix(sigma_mat, bunch_params, error, print_err, n_mat)

Routine to calc emittances and Twiss function from a beam sigma matrix.
See: Andy Wolski "Alternative approach to general coupled linear optics".

Parameters
----------
sigma_mat : float
    Sigma matrix.
print_err : bool, optional
    If present and False then suppress "no eigen-system found" messages.

Returns
-------
bunch_params : BunchParamsStruct
    Holds Twiss and emittance info.
error : bool
    Set True if there is an error. Can happen if the emittance of a mode is zero.
n_mat : float
    N matrix defined in Wolski Eq 44 and used to convert from action-angle coords to lab coords (Wolski Eq
    51.).
)"""
  );
  m.def(
      "calc_spin_params",
      &Bmad::calc_spin_params,
      py::arg("bunch"),
      R"""(Subroutine calc_spin_params (bunch, bunch_params)

Rotine to calculate spin averages

Parameters
----------
bunch : BunchStruct
    Bunch of spins

Returns
-------
bunch_params : BunchParamStruct
    Structure holding average
centroid%spin : 
    (x,y,z) polarization.
)"""
  );
  m.def(
      "calc_super_slave_key",
      &Bmad::calc_super_slave_key,
      py::arg("lord1"),
      py::arg("lord2"),
      py::arg("create_jumbo_slave") = py::none(),
      R"""(Parameters
----------
lord1 : EleStruct
    First slave. .key
lord2 : EleStruct
    Second slave. .key .sub_key
slave : EleStruct
    Super_slave element.
create_jumbo_slave : unknown, optional
    If True then slave.key will be set to em_field. Default is False.
)"""
  );
  py::class_<Bmad::CalcWallRadius, std::unique_ptr<Bmad::CalcWallRadius>>(
      m,
      "CalcWallRadius",
      "calc_wall_radius return type"
  )
      .def_readonly("r_wall", &Bmad::CalcWallRadius::r_wall)
      .def_readonly("dr_dtheta", &Bmad::CalcWallRadius::dr_dtheta)
      .def_readonly("ix_vertex", &Bmad::CalcWallRadius::ix_vertex)
      .def("__len__", [](const Bmad::CalcWallRadius &) { return 3; })
      .def("__getitem__", [](const Bmad::CalcWallRadius &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.r_wall);
        if (i == 1)
          return py::cast(s.dr_dtheta);
        if (i == 2)
          return py::cast(s.ix_vertex);
        throw py::index_error();
      });
  m.def(
      "calc_wall_radius",
      &Bmad::calc_wall_radius,
      py::arg("v"),
      py::arg("cos_ang"),
      py::arg("sin_ang"),
      R"""(Subroutine calc_wall_radius (v, cos_ang, sin_ang, r_wall, dr_dtheta, ix_vertex)

Routine to calculate the wall radius at a given angle for a given cross-section
Additionally, the transverse directional derivative is calculated.

Module needed:
  use wall3d_mod

Parameters
----------
v : Wall3DVertexStruct
    Array of vertices that make up the cross-section.
cos_ang : float
    cosine of the transverse photon position.
sin_ang : float
    sine of the transverse photon position.

Returns
-------
r_wall : float
    Wall radius at given angle.
dr_dtheta : float
    derivative of r_wall.
ix_vertex : int
    Wall at given angle is between v(ix_vertex-1) and v(ix_vertex). If ix_vertex = 1 then Wall at given angle
    is between v(N) and v(1) where N = size(v).
)"""
  );
  m.def(
      "calc_z_tune",
      &Bmad::calc_z_tune,
      py::arg("branch"),
      R"""(Parameters
----------
branch : BranchStruct
    Lattice branch
    This parameter is an input/output and is modified in-place. As an output: Synchrotron tune (radians). If
    unstable tune = 0.
    This parameter is an input/output and is modified in-place. As an output: Is the mode stable? If no rf
    then tune is zero but is stable.
    This parameter is an input/output and is modified in-place. As an output: 6x6 1-turn matrix.
)"""
  );
  m.def(
      "canonical_to_angle_coords",
      &Bmad::canonical_to_angle_coords,
      py::arg("orbit"),
      py::arg("coord_type") = py::none(),
      R"""(Parameters
----------
orbit : CoordStruct
    Orbit in canonical coordinates.
    This parameter is an input/output and is modified in-place. As an output: Orbit in angular coordinates.
coord_type : unknown, optional
    Angular coordinates type '' (default): (x, x' = dx/ds, y, y' = dy/ds, z, pz) 'ZGOUBI':     (x, x' = dx/ds,
    y, y' = dy/ds, dt = -z / (beta * c), pz)
)"""
  );
  m.def(
      "cbar_to_c",
      &Bmad::cbar_to_c,
      py::arg("cbar_mat"),
      py::arg("a"),
      py::arg("b"),
      R"""(Parameters
----------
cbar_mat : float
    Cbar matrix.
a : TwissStruct
    a-mode Twiss parameters
b : TwissStruct
    b-mode Twiss parameters
c_mat : float
    C matrix.
)"""
  );
  m.def(
      "check_aperture_limit",
      &Bmad::check_aperture_limit,
      py::arg("orb"),
      py::arg("ele"),
      py::arg("particle_at"),
      py::arg("param"),
      py::arg("old_orb") = py::none(),
      py::arg("check_momentum") = py::none(),
      R"""(Parameters
----------
orb : CoordStruct
    coordinates of a particle.
ele : EleStruct
    Element holding the aperture
particle_at : int
    first_track_edge$, second_track_edge$, surface$, in_between$
param : LatParamStruct
    Lattice global parameter structure.
old_orb : CoordStruct, optional
    Old coordinates at last check. Needed if ele.aperture_at = wall_transition$.
check_momentum : bool, optional
    If present and false then checking of p_x and p_y will be disabled.
)"""
  );
  m.def(
      "check_controller_controls",
      &Bmad::check_controller_controls,
      py::arg("ele_key"),
      py::arg("contrl"),
      py::arg("name"),
      R"""(Parameters
----------
ele_key : int
    Element type. overlay$, etc.
contrl : ControlStruct
    control info. 1 element for each slave.
name : unknown
    Lord name. Used for error reporting.
err : bool
    Set true if there is a problem. False otherwise.
)"""
  );
  m.def(
      "check_for_superimpose_problem",
      &Bmad::check_for_superimpose_problem,
      py::arg("branch"),
      py::arg("super_ele"),
      py::arg("err_flag"),
      py::arg("ref_ele") = py::none(),
      py::arg("wrap"),
      R"""(Subroutine check_for_superimpose_problem (branch, super_ele, err_flag, ref_ele, wrap)

Subroutine to check if there is a problem superimposing an element when there is multipass.
In particular will check that:
  1) If the ref_ele is part of a multipass region then super_ele must be superimposed
     within the region.
Or:
  2) If the ref_ele is not part of a multipass region then super_ele must also not
     be part of a multipass region.

This subroutine is used by bmad_parser and bmad_parser2.
This subroutine is not intended for general use.

)"""
  );
  py::class_<Bmad::CheckIfSInBounds, std::unique_ptr<Bmad::CheckIfSInBounds>>(
      m,
      "CheckIfSInBounds",
      "check_if_s_in_bounds return type"
  )
      .def_readonly("err_flag", &Bmad::CheckIfSInBounds::err_flag)
      .def_readonly("translated_s", &Bmad::CheckIfSInBounds::translated_s)
      .def("__len__", [](const Bmad::CheckIfSInBounds &) { return 2; })
      .def("__getitem__", [](const Bmad::CheckIfSInBounds &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.translated_s);
        throw py::index_error();
      });
  m.def(
      "check_if_s_in_bounds",
      &Bmad::check_if_s_in_bounds,
      py::arg("branch"),
      py::arg("s"),
      py::arg("print_err") = py::none(),
      R"""(Parameters
----------
branch : BranchStruct
    Branch
s : float
    longitudinal position in the given branch.
err_flag : bool
    Set True if s position is out-of-bounds. False otherwise.
translated_s : float
    position translated to the range [0, branch_length]
print_err : bool, optional
    Print error message if there is an error? Default is True.
)"""
  );
  m.def(
      "choose_quads_for_set_tune",
      &Bmad::choose_quads_for_set_tune,
      py::arg("branch"),
      py::arg("dk1"),
      py::arg("eles"),
      py::arg("mask") = py::none(),
      R"""(Parameters
----------
branch : BranchStruct
    Lattice branch.
dk1 : float
    Weights for the quadrupoles. All values will be +1 or -1.
eles : ElePointerStruct
    eles(i).ele points to element with dk1(i) weight.
mask : unknown, optional
    If present, assign weight of zero for all quads that do not match. That is, no variation for matching
    quads.
err_flag : bool
    Set True if there is not one quad with positive dk1 and one quad with negative dk1.
)"""
  );
  py::class_<Bmad::ChromCalc, std::unique_ptr<Bmad::ChromCalc>>(
      m,
      "ChromCalc",
      "chrom_calc return type"
  )
      .def_readonly("chrom_a", &Bmad::ChromCalc::chrom_a)
      .def_readonly("chrom_b", &Bmad::ChromCalc::chrom_b)
      .def_readonly("err_flag", &Bmad::ChromCalc::err_flag)
      .def_readonly("low_E_lat", &Bmad::ChromCalc::low_E_lat)
      .def_readonly("high_E_lat", &Bmad::ChromCalc::high_E_lat)
      .def("__len__", [](const Bmad::ChromCalc &) { return 5; })
      .def("__getitem__", [](const Bmad::ChromCalc &s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.chrom_a);
        if (i == 1)
          return py::cast(s.chrom_b);
        if (i == 2)
          return py::cast(s.err_flag);
        if (i == 3)
          return py::cast(s.low_E_lat);
        if (i == 4)
          return py::cast(s.high_E_lat);
        throw py::index_error();
      });
  m.def(
      "chrom_calc",
      &Bmad::chrom_calc,
      py::arg("lat"),
      py::arg("delta_e"),
      py::arg("pz") = py::none(),
      py::arg("low_E_orb") = py::none(),
      py::arg("high_E_orb") = py::none(),
      py::arg("ix_branch") = py::none(),
      py::arg("orb0") = py::none(),
      R"""(Parameters
----------
lat : LatStruct
    Lat
delta_e : float
    +/- Delta energy used for the calculation. Notice that the energy difference
    This parameter is an input/output and is modified in-place. As an output: Set to 1.0d-4 if on input
    DELTA_E =< 0.
chrom_a : float
    a-mode chromaticity.
chrom_b : float
    b-mode chromaticity.
err_flag : bool
    Set true if there is an error. False otherwise.
pz : float, optional
    reference momentum about which to calculate. Default is 0.
low_E_lat : LatStruct
    Lattice with RF off and matrices computed at E_lat +pz - delta_e
high_E_lat : LatStruct
    Lattice with RF off and matrices computed at E_lat +pz + delta_e
low_E_orb : CoordStruct
    Orbit computed at E_lat + pz - delta_e.
high_E_orb : CoordStruct
    Orbit computed at E_lat + pz + delta_e.
ix_branch : int, optional
    Index of the lattice branch to use. Default is 0.
orb0 : CoordStruct, optional
    On-energy orbit at start (fixer point). Default is the branch.particle_start. Only needed if lattice
    branch has an open geometry.
)"""
  );
  m.def(
      "chrom_tune",
      &Bmad::chrom_tune,
      py::arg("lat"),
      py::arg("delta_e"),
      py::arg("target_x"),
      py::arg("target_y"),
      py::arg("err_tol"),
      R"""(Parameters
----------
lat : LatStruct
    Lat to use,
    This parameter is an input/output and is modified in-place. As an output: Lat with sextupole set
delta_e : float
    Delta energy used for the calculation.
    This parameter is an input/output and is modified in-place. As an output: Set to 1.0d-4 if on input
    DELTA_E =< 0.
target_x : float
    Target X Chromaticity
target_y : float
    Target Y Chromaticity
err_tol : float
    Max allowable Error: Error = | X_Target - X_Actual | + | Y_Target -Y_Actual | A good number is: err_tol =
    0.05_rp
err_flag : bool
    .false. if match successful, .true. if failed Fails if takes longer than 100 iterations. If it fails the
    sextupoles are set to the last value calculated. Note: This subroutine assumes the Twiss parameters have
    been computed.
)"""
  );
  m.def(
      "classical_radius",
      &Bmad::classical_radius,
      py::arg("species"),
      R"""(Parameters
----------
species : int
    Species of particle.
radius : float
    Classical radius.
)"""
  );
  m.def(
      "clear_lat_1turn_mats",
      &Bmad::clear_lat_1turn_mats,
      R"""(Parameters
----------
lat : LatStruct
    Lat with 1-turn matrices cleared.
)"""
  );
  m.def(
      "clear_taylor_maps_from_elements",
      &Bmad::clear_taylor_maps_from_elements,
      py::arg("lat"),
      R"""(Parameters
----------
lat : LatStruct
    Lattice
    This parameter is an input/output and is modified in-place. As an output: Lattice with all maps cleared
)"""
  );
  m.def(
      "closed_orbit_calc",
      &Bmad::closed_orbit_calc,
      py::arg("lat"),
      py::arg("closed_orb"),
      py::arg("i_dim") = py::none(),
      py::arg("direction") = py::none(),
      py::arg("ix_branch") = py::none(),
      py::arg("print_err") = py::none(),
      R"""(Parameters
----------
lat : LatStruct
    Lat to track through.
closed_orb : CoordStruct
    closed_orb(nt) is the initial guess where nt = 0 for direction = 1 and nt = lat.n_ele_track for direction
    = -1. Additionally, if i_dim = 4, then closed_orb(nt).vec(6) is used as the energy
    This parameter is an input/output and is modified in-place. As an output: Closed orbit. closed_orb(i)
i_dim : int, optional
    Phase space dimensions to use: = 4  Transverse closed orbit at constant energy (RF off). (dE/E =
    closed_orb(0).vec(6)) = 5 Transverse closed orbit at constant energy (RF off) with the energy adjusted so
    that vec(5) is the same at the beginning and at the end. = 6 True closed orbit.
direction : int, optional
    Direction of tracking.
ix_branch : int, optional
    Lattice branch to find the closed orbit of.
err_flag : bool
    Set true if there is an error. False otherwise.
print_err : bool, optional
    Print error message if calc does not converge? Default is True. Note: Condition messages like no RF
    voltage with i_dim = 6 will always be printed.
)"""
  );
  m.def(
      "closed_orbit_from_tracking",
      &Bmad::closed_orbit_from_tracking,
      py::arg("lat"),
      py::arg("closed_orb"),
      py::arg("i_dim"),
      py::arg("eps_rel") = py::none(),
      py::arg("eps_abs") = py::none(),
      py::arg("init_guess") = py::none(),
      R"""(Parameters
----------
lat : LatStruct
    Lat to track through.
closed_orb : CoordStruct
    closed orbit.
i_dim : int
    = 2,4  Transverse closed orbit at constant energy.
eps_rel : float, optional
    Relative allowed error.
eps_abs : float, optional
    Absolute allowed error.
init_guess : CoordStruct, optional
    Starting guess for the closed orbit at the start of the lattice. Set init_guess.vec(6) to the appropriate
    value of pz when calculating off-energy orbits. If not present then the origin will be used.
err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "cmplx_re_str",
      &Bmad::cmplx_re_str,
      py::arg("cmp"),
      py::arg("str_out"),
      R"""(Parameters
----------
cmp : 
str_out : 
)"""
  );
  m.def(
      "combine_consecutive_elements",
      &Bmad::combine_consecutive_elements,
      py::arg("lat"),
      R"""(Parameters
----------
lat : LatStruct
    Lattice.
    This parameter is an input/output and is modified in-place. As an output: Lattice with elements combined.
error : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "complex_taylor_clean",
      &Bmad::complex_taylor_clean,
      py::arg("complex_taylor"),
      R"""(Parameters
----------
complex_taylor : 
)"""
  );
  m.def(
      "complex_taylor_coef",
      py::overload_cast<ComplexTaylorStruct &, FArray1D<Int> &, std::complex<double>>(
          &Bmad::complex_taylor_coef
      ),
      py::arg("complex_taylor"),
      py::arg("exp"),
      py::arg("coef"),
      R"""(Function complex_taylor_coef (complex_taylor, exp)
Function complex_taylor_coef (complex_taylor, i1, i2, i3, i4, i5, i6, i7, i8, i9)

Function to return the coefficient for a particular complex_taylor term
from a complex_taylor Series.

Note: complex_taylor_coef is overloaded by:
  complex_taylor_coef1 (complex_taylor, exp)
  complex_taylor_coef2 (complex_taylor, i1, i2, i3, i4, i5, i6, i7, i8, i9)
Using the complex_taylor_coef2 form limits obtaining coefficients to 9th order
or less. Also: complex_taylor_coef2 does not check that all i1, ..., i9 are between
1 and 6.

For example: To get the 2nd order term corresponding to
  y(out) = Coef * p_z(in)^2
[This is somtimes refered to as the T_366 term]
The call would be:
  type (complex_taylor_struct) complex_taylor(6)      ! complex_taylor Map
  ...
  coef = complex_taylor_coef (complex_taylor(3), 6, 6)  ! 1st possibility or ...
  coef = complex_taylor_coef (complex_taylor(3), [0, 0, 0, 0, 0, 2 ])

Input (complex_taylor_coef1):
  complex_taylor -- complex_taylor_struct: complex_taylor series.
  exp(6)      -- Integer: Array of exponent indices.

Input (complex_taylor_coef2):
  complex_taylor -- complex_taylor_struct: complex_taylor series.
  i1, ..., i9 -- Integer, optional: indexes (each between 1 and 6).


Returns
-------
complex_taylor_coef : complex
    Coefficient.
)"""
  );
  m.def(
      "complex_taylor_coef",
      py::overload_cast<
          ComplexTaylorStruct &,
          std::optional<int>,
          std::optional<int>,
          std::optional<int>,
          std::optional<int>,
          std::optional<int>,
          std::optional<int>,
          std::optional<int>,
          std::optional<int>,
          std::optional<int>,
          std::complex<double>>(&Bmad::complex_taylor_coef),
      py::arg("complex_taylor"),
      py::arg("i1") = py::none(),
      py::arg("i2") = py::none(),
      py::arg("i3") = py::none(),
      py::arg("i4") = py::none(),
      py::arg("i5") = py::none(),
      py::arg("i6") = py::none(),
      py::arg("i7") = py::none(),
      py::arg("i8") = py::none(),
      py::arg("i9") = py::none(),
      py::arg("coef"),
      R"""(Function complex_taylor_coef (complex_taylor, exp)
Function complex_taylor_coef (complex_taylor, i1, i2, i3, i4, i5, i6, i7, i8, i9)

Function to return the coefficient for a particular complex_taylor term
from a complex_taylor Series.

Note: complex_taylor_coef is overloaded by:
  complex_taylor_coef1 (complex_taylor, exp)
  complex_taylor_coef2 (complex_taylor, i1, i2, i3, i4, i5, i6, i7, i8, i9)
Using the complex_taylor_coef2 form limits obtaining coefficients to 9th order
or less. Also: complex_taylor_coef2 does not check that all i1, ..., i9 are between
1 and 6.

For example: To get the 2nd order term corresponding to
  y(out) = Coef * p_z(in)^2
[This is somtimes refered to as the T_366 term]
The call would be:
  type (complex_taylor_struct) complex_taylor(6)      ! complex_taylor Map
  ...
  coef = complex_taylor_coef (complex_taylor(3), 6, 6)  ! 1st possibility or ...
  coef = complex_taylor_coef (complex_taylor(3), [0, 0, 0, 0, 0, 2 ])

Input (complex_taylor_coef1):
  complex_taylor -- complex_taylor_struct: complex_taylor series.
  exp(6)      -- Integer: Array of exponent indices.

Input (complex_taylor_coef2):
  complex_taylor -- complex_taylor_struct: complex_taylor series.
  i1, ..., i9 -- Integer, optional: indexes (each between 1 and 6).


Returns
-------
complex_taylor_coef : complex
    Coefficient.
)"""
  );
  m.def(
      "complex_taylor_equal_complex_taylor",
      &Bmad::complex_taylor_equal_complex_taylor,
      py::arg("complex_taylor1"),
      py::arg("complex_taylor2"),
      R"""(Parameters
----------
complex_taylor1 : 
complex_taylor2 : 
)"""
  );
  m.def(
      "complex_taylor_exponent_index",
      &Bmad::complex_taylor_exponent_index,
      py::arg("expn"),
      R"""(Function complex_taylor_exponent_index(expn) result(index)

Function to associate a unique number with a complex_taylor exponent.

The number associated with a complex_taylor_term that is used for the sort is:
    number = sum(exp(i))*10^6 + exp(6)*10^5 + ... + exp(1)*10^0
where exp(1) is the exponent for x, exp(2) is the exponent for P_x, etc.

Parameters
----------
expn : int
    complex_taylor exponent

Returns
-------
index : int
    Sorted complex_taylor series.
)"""
  );
  m.def(
      "complex_taylor_make_unit",
      &Bmad::complex_taylor_make_unit,
      py::arg("complex_taylor"),
      R"""(Subroutine complex_taylor_make_unit (complex_taylor)

Subroutine to make the unit complex_taylor map:
      r(out) = Map * r(in) = r(in)


Returns
-------
complex_taylor : ComplexTaylorStruct
    Unit complex_taylor map .
)"""
  );
  py::class_<Bmad::ComplexTaylorToMat6, std::unique_ptr<Bmad::ComplexTaylorToMat6>>(
      m,
      "ComplexTaylorToMat6",
      "complex_taylor_to_mat6 return type"
  )
      .def_readonly("vec0", &Bmad::ComplexTaylorToMat6::vec0)
      .def_readonly("mat6", &Bmad::ComplexTaylorToMat6::mat6)
      .def("__len__", [](const Bmad::ComplexTaylorToMat6 &) { return 2; })
      .def("__getitem__", [](const Bmad::ComplexTaylorToMat6 &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.vec0);
        if (i == 1)
          return py::cast(s.mat6);
        throw py::index_error();
      });
  m.def(
      "complex_taylor_to_mat6",
      &Bmad::complex_taylor_to_mat6,
      py::arg("a_complex_taylor"),
      py::arg("r_in"),
      py::arg("r_out") = py::none(),
      R"""(Subroutine complex_taylor_to_mat6 (a_complex_taylor, r_in, vec0, mat6, r_out)

Subroutine to calculate, from a complex_taylor map and about some trajectory:
  The 1st order (Jacobian) transfer matrix.

Parameters
----------
a_complex_taylor : ComplexTaylorStruct
    complex_taylor map.
r_in : complex
    Coordinates at the input.

Returns
-------
vec0 : complex
    0th order tranfsfer map
mat6 : complex
    1st order transfer map (6x6 matrix).
r_out : complex
    Coordinates at output.
)"""
  );
  m.def(
      "complex_taylors_equal_complex_taylors",
      &Bmad::complex_taylors_equal_complex_taylors,
      py::arg("complex_taylor1"),
      py::arg("complex_taylor2"),
      R"""(Parameters
----------
complex_taylor1 : 
complex_taylor2 : 
)"""
  );
  m.def(
      "compute_slave_coupler",
      &Bmad::compute_slave_coupler,
      py::arg("slave"),
      R"""(Subroutine compute_slave_coupler (slave)

This routine is not meant for general use.

)"""
  );
  m.def(
      "concat_ele_taylor",
      &Bmad::concat_ele_taylor,
      py::arg("orb_taylor"),
      py::arg("ele"),
      py::arg("err_flag"),
      py::arg("spin_taylor") = py::none(),
      R"""(Subroutine concat_ele_taylor (orb_taylor, ele, err_flag, spin_taylor)

Routine to concatinate an orbital taylor map and, optionally if present and
bmad_com%spin_tracking_on = T, a spin taylor map.

Transform:
  orb_taylor[x] -> ele_taylor(orb_taylor[x])
If ele%taylor_map_includes_offsets = True:  ele_taylor == ele%taylor
If ele%taylor_map_includes_offsets = False: ele_taylor == ele%taylor + offset corrections.

Also see: concat_taylor

Parameters
----------
orb_taylor : TaylorStruct
    Orbital Taylor map.
ele : EleStruct
    Element containing a Taylor map.
spin_taylor : TaylorStruct, optional
    Spin map to propagate
Output : 
orb_taylor : TaylorStruct
    Concatinated orbital map
err_flag : bool
    Set True if there is an error. False otherwise.
spin_taylor : TaylorStruct, optional
    Concatinated spin map.

Notes
-----
Related routines:
concat_taylor
)"""
  );
  m.def(
      "concat_taylor",
      &Bmad::concat_taylor,
      py::arg("taylor1"),
      py::arg("taylor2"),
      py::arg("taylor3"),
      R"""(Subroutine concat_taylor (taylor1, taylor2, taylor3)

Subroutine to concatinate two taylor maps:
  taylor3[x] = taylor2(taylor1[x])

Note: In general, if taylor2 is a component of an ele_struct, use
concat_ele_taylor instead.

Parameters
----------
taylor1 : TaylorStruct
    Taylor map.
taylor2 : TaylorStruct
    Taylor map.
Output : 
taylor3 : TaylorStruct
    Concatinated map
)"""
  );
  py::class_<Bmad::ConcatTransferMat, std::unique_ptr<Bmad::ConcatTransferMat>>(
      m,
      "ConcatTransferMat",
      "concat_transfer_mat return type"
  )
      .def_readonly("mat_out", &Bmad::ConcatTransferMat::mat_out)
      .def_readonly("vec_out", &Bmad::ConcatTransferMat::vec_out)
      .def("__len__", [](const Bmad::ConcatTransferMat &) { return 2; })
      .def("__getitem__", [](const Bmad::ConcatTransferMat &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.mat_out);
        if (i == 1)
          return py::cast(s.vec_out);
        throw py::index_error();
      });
  m.def(
      "concat_transfer_mat",
      &Bmad::concat_transfer_mat,
      py::arg("mat_1"),
      py::arg("vec_1"),
      py::arg("mat_0"),
      py::arg("vec_0"),
      R"""(Subroutine concat_transfer_mat (mat_1, vec_1, mat_0, vec_0, mat_out, vec_out)

Routine to concatinate two linear maps:
  mat_out = matmul(mat_1, mat_0)
  vec_out = matmul(mat_1, vec_0) + vec_1

Parameters
----------
mat_1 : float
    Map from s1 to s2
vec_1 : float
    Map from s1 to s2
mat_0 : float
    Map from s0 to s1
vec_0 : float
    Map from s0 to s1

Returns
-------
mat_out : float
    Map from s0 to s2
vec_out : float
    Map from s0 to s2
)"""
  );
  m.def(
      "control_bookkeeper",
      &Bmad::control_bookkeeper,
      py::arg("lat"),
      py::arg("ele") = py::none(),
      py::arg("err_flag") = py::none(),
      R"""(Parameters
----------
lat : LatStruct
    lattice to be used
ele : EleStruct, optional
    Element whose attribute values have been changed. If not present bookkeeping will be done
err_flag : bool, optional
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "convert_bend_exact_multipole",
      &Bmad::convert_bend_exact_multipole,
      py::arg("g"),
      py::arg("out_type"),
      py::arg("an"),
      py::arg("bn"),
      R"""(Parameters
----------
g : float
    1/rho bending strength.
out_type : int
    Output type: horizontally_pure$ or vertically_pure$.
an : float
    Skew multipoles.
    This parameter is an input/output and is modified in-place. As an output: Converted skew multipoles.
bn : float
    Non-skew multipoles.
    This parameter is an input/output and is modified in-place. As an output: Converted Non-skew multipoles.
)"""
  );
  py::class_<Bmad::ConvertCoords, std::unique_ptr<Bmad::ConvertCoords>>(
      m,
      "ConvertCoords",
      "convert_coords return type"
  )
      .def_readonly("out_type_str", &Bmad::ConvertCoords::out_type_str)
      .def_readonly("coord_out", &Bmad::ConvertCoords::coord_out)
      .def_readonly("err_flag", &Bmad::ConvertCoords::err_flag)
      .def("__len__", [](const Bmad::ConvertCoords &) { return 3; })
      .def("__getitem__", [](const Bmad::ConvertCoords &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.out_type_str);
        if (i == 1)
          return py::cast(s.coord_out);
        if (i == 2)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "convert_coords",
      &Bmad::convert_coords,
      py::arg("in_type_str"),
      py::arg("coord_in"),
      py::arg("ele"),
      R"""(Parameters
----------
in_type_str : unknown
    type of the input coords.
coord_in : CoordStruct
    Input coordinates.
ele : EleStruct
    Provides the Twiss parameters.
out_type_str : unknown
    type of the output coords.
coord_out : CoordStruct
    Output coordinates.
err_flag : bool
    Set True if there is an error. False otherwise. in_type_str and out_type_str can be: 'LAB'
    {x, x', y, y', z, z'} 'MODE'               {a, a', b, b', z, z'} 'NORMALIZED'         {a_bar, a'_bar,
    b_bar, b'_bar, z_bar, z'_bar} 'ACTION-ANGLE'       {j_a, phi_a, j_b, phi_b, j_z,  phi_z} x_vec = V_mat *
    (a_vec + eta_vec * z') a_bar  =  sqrt(2*j_a) * cos(phi_a) a'_bar = -sqrt(2*j_a) * sin(phi_a) Note: 1) If
    ELE.Z.BETA = 0 then ELE.Z.BETA is set to 1. 2) phases are in radians
)"""
  );
  m.def(
      "convert_field_ele_to_lab",
      &Bmad::convert_field_ele_to_lab,
      py::arg("ele"),
      py::arg("s_here"),
      py::arg("forward_transform"),
      py::arg("calc_dfield") = py::none(),
      py::arg("calc_potential") = py::none(),
      R"""(Subroutine convert_field_ele_to_lab (ele, s_here, forward_transform, field, calc_dfield, calc_potential)

Convert fields: ele to lab coords

Parameters
----------
ele : EleStruct
    Lattice element.
s_here : 
    real(rp) S-position.
forward_transform : 
    Transform foward (to lab) or reverse.
calc_dfield : bool, optional
    If present and True then calculate the field derivatives.
calc_potential : bool, optional
    Calc electric and magnetic potentials? Default is false. This is experimental and only implemented for
    wigglers at present.

Returns
-------
field : EmFieldStruct
    EM field.
)"""
  );
  m.def(
      "convert_local_cartesian_to_local_curvilinear",
      &Bmad::convert_local_cartesian_to_local_curvilinear,
      py::arg("x"),
      py::arg("z"),
      py::arg("g"),
      py::arg("xout"),
      py::arg("sout"),
      R"""(Parameters
----------
x : 
z : 
g : 
xout : 
sout : 
)"""
  );
  m.def(
      "convert_local_curvilinear_to_local_cartesian",
      &Bmad::convert_local_curvilinear_to_local_cartesian,
      py::arg("x"),
      py::arg("s"),
      py::arg("g"),
      py::arg("xout"),
      py::arg("zout"),
      R"""(Parameters
----------
x : 
s : 
g : 
xout : 
zout : 
)"""
  );
  m.def(
      "convert_particle_coordinates_s_to_t",
      &Bmad::convert_particle_coordinates_s_to_t,
      py::arg("particle"),
      py::arg("s_body"),
      py::arg("orientation"),
      R"""(Parameters
----------
particle : CoordStruct
    Particle with .vec(:) in s-coords.
s_body : float
    s-position in element body coords.
orientation : int
    ele.orientation for vec(6).
)"""
  );
  m.def(
      "convert_particle_coordinates_t_to_s",
      &Bmad::convert_particle_coordinates_t_to_s,
      py::arg("particle"),
      py::arg("ele"),
      py::arg("use_downstream_p0c") = py::none(),
      R"""(Parameters
----------
particle : CoordStruct
    Particle with .vec(:) in t-coords.
ele : EleStruct
    Element particle is going through.
s_body : float
    s-position in element body coords.
use_downstream_p0c : bool, optional
    If True (the default), use ele.value(p0c$) as the reference momentum. If False, use ele.value(p0c_start$)
    as the reference.
)"""
  );
  py::class_<Bmad::ConvertPcTo, std::unique_ptr<Bmad::ConvertPcTo>>(
      m,
      "ConvertPcTo",
      "convert_pc_to return type"
  )
      .def_readonly("E_tot", &Bmad::ConvertPcTo::E_tot)
      .def_readonly("gamma", &Bmad::ConvertPcTo::gamma)
      .def_readonly("kinetic", &Bmad::ConvertPcTo::kinetic)
      .def_readonly("beta", &Bmad::ConvertPcTo::beta)
      .def_readonly("brho", &Bmad::ConvertPcTo::brho)
      .def_readonly("beta1", &Bmad::ConvertPcTo::beta1)
      .def_readonly("err_flag", &Bmad::ConvertPcTo::err_flag)
      .def("__len__", [](const Bmad::ConvertPcTo &) { return 7; })
      .def("__getitem__", [](const Bmad::ConvertPcTo &s, int i) -> py::object {
        if (i < 0)
          i += 7;
        if (i == 0)
          return py::cast(s.E_tot);
        if (i == 1)
          return py::cast(s.gamma);
        if (i == 2)
          return py::cast(s.kinetic);
        if (i == 3)
          return py::cast(s.beta);
        if (i == 4)
          return py::cast(s.brho);
        if (i == 5)
          return py::cast(s.beta1);
        if (i == 6)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "convert_pc_to",
      &Bmad::convert_pc_to,
      py::arg("pc"),
      py::arg("particle"),
      R"""(Parameters
----------
pc : float
    Particle momentum
particle : int
    Type of particle. positron$, etc.
E_tot : float
    Total energy of the particle.
gamma : float
    Gamma factor.
kinetic : float
    Kinetic energy
beta : float
    velocity / c_light
brho : float
    Nominal B_field*rho_bend
beta1 : float
    1 - beta. Equal to 1/(2*gamma^2) in ultra-rel limit.
err_flag : bool
    Set true if there is an error. False otherwise.
)"""
  );
  py::class_<Bmad::ConvertTotalEnergyTo, std::unique_ptr<Bmad::ConvertTotalEnergyTo>>(
      m,
      "ConvertTotalEnergyTo",
      "convert_total_energy_to return type"
  )
      .def_readonly("gamma", &Bmad::ConvertTotalEnergyTo::gamma)
      .def_readonly("kinetic", &Bmad::ConvertTotalEnergyTo::kinetic)
      .def_readonly("beta", &Bmad::ConvertTotalEnergyTo::beta)
      .def_readonly("pc", &Bmad::ConvertTotalEnergyTo::pc)
      .def_readonly("brho", &Bmad::ConvertTotalEnergyTo::brho)
      .def_readonly("beta1", &Bmad::ConvertTotalEnergyTo::beta1)
      .def_readonly("err_flag", &Bmad::ConvertTotalEnergyTo::err_flag)
      .def("__len__", [](const Bmad::ConvertTotalEnergyTo &) { return 7; })
      .def("__getitem__", [](const Bmad::ConvertTotalEnergyTo &s, int i) -> py::object {
        if (i < 0)
          i += 7;
        if (i == 0)
          return py::cast(s.gamma);
        if (i == 1)
          return py::cast(s.kinetic);
        if (i == 2)
          return py::cast(s.beta);
        if (i == 3)
          return py::cast(s.pc);
        if (i == 4)
          return py::cast(s.brho);
        if (i == 5)
          return py::cast(s.beta1);
        if (i == 6)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "convert_total_energy_to",
      &Bmad::convert_total_energy_to,
      py::arg("E_tot"),
      py::arg("particle"),
      py::arg("print_err") = py::none(),
      R"""(Parameters
----------
E_tot : float
    Total energy of the particle.
particle : int
    Type of particle. positron$, etc.
gamma : float
    Gamma factor. Set to -1 for photons.
kinetic : float
    Kinetic energy
beta : float
    velocity / c_light
pc : float
    Particle momentum
brho : float
    Nominal B_field*rho_bend
beta1 : float
    1 - beta. Equal to 1/(2*gamma^2) in ultra-rel limit.
err_flag : bool
    Set true if there is an error. False otherwise.
print_err : bool, optional
    Print error message if E_tot < particle mass? Default is True.
)"""
  );
  py::class_<Bmad::ConverterDistributionParser, std::unique_ptr<Bmad::ConverterDistributionParser>>(
      m,
      "ConverterDistributionParser",
      "converter_distribution_parser return type"
  )
      .def_readonly("delim", &Bmad::ConverterDistributionParser::delim)
      .def_readonly("delim_found", &Bmad::ConverterDistributionParser::delim_found)
      .def_readonly("err_flag", &Bmad::ConverterDistributionParser::err_flag)
      .def("__len__", [](const Bmad::ConverterDistributionParser &) { return 3; })
      .def("__getitem__", [](const Bmad::ConverterDistributionParser &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.delim);
        if (i == 1)
          return py::cast(s.delim_found);
        if (i == 2)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "converter_distribution_parser",
      &Bmad::converter_distribution_parser,
      py::arg("ele"),
      R"""(Parameters
----------
ele : EleStruct
    Converter element.
    This parameter is an input/output and is modified in-place. As an output: Converter element with
    .converter field set.
delim : unknown
    Ending delimitor.
delim_found : bool
    Has a delimitor been found?
err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "coord_equal_coord",
      &Bmad::coord_equal_coord,
      py::arg("coord2"),
      R"""( Subroutine coord_equal_coord (coord1, coord2)

 Subroutine that is used to set one coord equal to another.

 Note: This subroutine is called by the overloaded equal sign:
		coord1 = coord2

Parameters
----------
coord2 : CoordStruct
    Input coord.

Returns
-------
coord1 : CoordStruct
    Output coord.
)"""
  );
  m.def(
      "coord_state_name",
      &Bmad::coord_state_name,
      py::arg("coord_state"),
      py::arg("one_word") = py::none(),
      R"""(Function coord_state_name (coord_state) result (state_str)

Routine to return the string representation of a coord%state state.

Parameters
----------
coord_state : int
    coord.state value

Returns
-------
state_str : unknown
    String representation.
)"""
  );
  m.def(
      "coords_body_to_local",
      &Bmad::coords_body_to_local,
      py::arg("body_position"),
      py::arg("ele"),
      py::arg("w_mat") = py::none(),
      py::arg("calculate_angles") = py::none(),
      py::arg("local_position"),
      R"""(Parameters
----------
body_position : FloorPositionStruct
    Element body frame coordinates.
ele : EleStruct
    element that local_position coordinates are relative to.
w_mat : float, optional
    W matrix at to transform vectors. v_local  = w_mat . v_body v_body   = transpose(w_mat) . v_local
calculate_angles : bool, optional
    calculate angles for local_position Default: True. False returns local_position angles (.theta, .phi,
    .psi) = 0.
local_position : FloorPositionStruct
    Local laboratory coordinates.
)"""
  );
  m.def(
      "coords_body_to_rel_exit",
      &Bmad::coords_body_to_rel_exit,
      py::arg("body_position"),
      py::arg("ele"),
      py::arg("w_mat") = py::none(),
      py::arg("calculate_angles") = py::none(),
      py::arg("rel_exit"),
      R"""(Parameters
----------
body_position : FloorPositionStruct
    Element body frame coordinates.
ele : EleStruct
    element that rel_exit coordinates are relative to.
w_mat : float, optional
    W matrix at to transform vectors. v_rel_exit = w_mat . v_body v_body     = transpose(w_mat) . v_rel_exit
calculate_angles : bool, optional
    calculate angles for rel_exit Default: True. False returns rel_exit angles (.theta, .phi, .psi) = 0.
rel_exit : FloorPositionStruct
    Cartesian coordinates relative to exit of the element.
)"""
  );
  py::class_<Bmad::CoordsCurvilinearToFloor, std::unique_ptr<Bmad::CoordsCurvilinearToFloor>>(
      m,
      "CoordsCurvilinearToFloor",
      "coords_curvilinear_to_floor return type"
  )
      .def_readonly("err_flag", &Bmad::CoordsCurvilinearToFloor::err_flag)
      .def_readonly("global_", &Bmad::CoordsCurvilinearToFloor::global)
      .def("__len__", [](const Bmad::CoordsCurvilinearToFloor &) { return 2; })
      .def("__getitem__", [](const Bmad::CoordsCurvilinearToFloor &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.err_flag);
        if (i == 1)
          return py::cast(s.global);
        throw py::index_error();
      });
  m.def(
      "coords_curvilinear_to_floor",
      &Bmad::coords_curvilinear_to_floor,
      py::arg("xys"),
      py::arg("branch"),
      R"""(Parameters
----------
xys : float
    (x, y, s) lab frame position vector.
branch : BranchStruct
    Lattice branch that defines the local reference coordinates.
err_flag : bool
    Set True if global floor position cannot be computed.
global : FloorPositionStruct
    Global floor position corresponding to (x, y, s) --    .w    -- W matrix to transform vectors: v_global =
    w_mat * v_local
)"""
  );
  py::class_<Bmad::CoordsFloorToCurvilinear, std::unique_ptr<Bmad::CoordsFloorToCurvilinear>>(
      m,
      "CoordsFloorToCurvilinear",
      "coords_floor_to_curvilinear return type"
  )
      .def_readonly("ele1", &Bmad::CoordsFloorToCurvilinear::ele1)
      .def_readonly("status", &Bmad::CoordsFloorToCurvilinear::status)
      .def_readonly("w_mat", &Bmad::CoordsFloorToCurvilinear::w_mat)
      .def_readonly("local_coords", &Bmad::CoordsFloorToCurvilinear::local_coords)
      .def("__len__", [](const Bmad::CoordsFloorToCurvilinear &) { return 4; })
      .def("__getitem__", [](const Bmad::CoordsFloorToCurvilinear &s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.ele1);
        if (i == 1)
          return py::cast(s.status);
        if (i == 2)
          return py::cast(s.w_mat);
        if (i == 3)
          return py::cast(s.local_coords);
        throw py::index_error();
      });
  m.def(
      "coords_floor_to_curvilinear",
      &Bmad::coords_floor_to_curvilinear,
      py::arg("floor_coords"),
      py::arg("ele0"),
      R"""(Parameters
----------
floor_coords : FloorPositionStruct
    .r = [X, Y, Z] position in global coordinates
ele0 : EleStruct
    Element to start the search at.
ele1 : EleStruct
    Element that local_coords is with respect to.
status : bool
    ok$             -> Local_coords found. patch_problem$  -> No solution due to a patch element.
w_mat : float
    W matrix at s, to transform vectors from floor to local. w_mat will only be well defined if status = ok$
local_coords : FloorPositionStruct
    .r = [x, y, s] position in curvilinear coordinates
)"""
  );
  py::class_<
      Bmad::CoordsFloorToLocalCurvilinear,
      std::unique_ptr<Bmad::CoordsFloorToLocalCurvilinear>>(
      m,
      "CoordsFloorToLocalCurvilinear",
      "coords_floor_to_local_curvilinear return type"
  )
      .def_readonly("status", &Bmad::CoordsFloorToLocalCurvilinear::status)
      .def_readonly("w_mat", &Bmad::CoordsFloorToLocalCurvilinear::w_mat)
      .def_readonly("local_position", &Bmad::CoordsFloorToLocalCurvilinear::local_position)
      .def("__len__", [](const Bmad::CoordsFloorToLocalCurvilinear &) { return 3; })
      .def("__getitem__", [](const Bmad::CoordsFloorToLocalCurvilinear &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.status);
        if (i == 1)
          return py::cast(s.w_mat);
        if (i == 2)
          return py::cast(s.local_position);
        throw py::index_error();
      });
  m.def(
      "coords_floor_to_local_curvilinear",
      &Bmad::coords_floor_to_local_curvilinear,
      py::arg("global_position"),
      py::arg("ele"),
      py::arg("relative_to") = py::none(),
      R"""(Parameters
----------
global_position : FloorPositionStruct
    .r = [X, Y, Z] position in global coordinates
ele : EleStruct
    element to find local coordinates of.
status : bool
    longitudinal position: inside$: Inside the element. upstream_end$: At upstream end of element or beyound.
w_mat : float
    W matrix at s, to transform vectors. v_global = w_mat.v_local v_local = transpose(w_mat).v_global
relative_to : int, optional
    not_set$ (default), upstream_end$, or downstream_end$. Force which end is used for z = 0. If
    upstream_end$, local_position.r(3) is relative to the upstream end which will not be the entrance end if
    ele.orientation = -1.
local_position : FloorPositionStruct
    .r = [x, y, z] position in local curvilinear coordinates.
)"""
  );
  m.def(
      "coords_floor_to_relative",
      &Bmad::coords_floor_to_relative,
      py::arg("floor0"),
      py::arg("global_position"),
      py::arg("calculate_angles") = py::none(),
      py::arg("is_delta_position") = py::none(),
      R"""(Parameters
----------
floor0 : FloorPositionStruct
    reference position
global_position : FloorPositionStruct
    global position
calculate_angles : bool, optional
    calculate angles for local_position Default: True.
is_delta_position : bool, optional
    If True then treat global_position.r as a difference position in global space and only rotate the position
    but not shift it. Default: False.
local_position : FloorPositionStruct
    position relative to floor0
)"""
  );
  m.def(
      "coords_local_curvilinear_to_body",
      &Bmad::coords_local_curvilinear_to_body,
      py::arg("local_position"),
      py::arg("ele"),
      py::arg("w_mat") = py::none(),
      py::arg("calculate_angles") = py::none(),
      py::arg("body_position"),
      R"""(Parameters
----------
local_position : FloorPositionStruct
    local coordinates.
ele : EleStruct
    element that coordinates are relative to.
w_mat : float, optional
    W matrix at to transform vectors. v_local  = w_mat . v_body v_body   = transpose(w_mat) . v_local
calculate_angles : bool, optional
    calculate angles for body_position Default: True. False returns body_position angles (.theta, .phi, .psi)
    = 0.
body_position : FloorPositionStruct
    Element coordinates relative to exit of the element.
)"""
  );
  py::class_<
      Bmad::CoordsLocalCurvilinearToFloor,
      std::unique_ptr<Bmad::CoordsLocalCurvilinearToFloor>>(
      m,
      "CoordsLocalCurvilinearToFloor",
      "coords_local_curvilinear_to_floor return type"
  )
      .def_readonly("w_mat", &Bmad::CoordsLocalCurvilinearToFloor::w_mat)
      .def_readonly("global_position", &Bmad::CoordsLocalCurvilinearToFloor::global_position)
      .def("__len__", [](const Bmad::CoordsLocalCurvilinearToFloor &) { return 2; })
      .def("__getitem__", [](const Bmad::CoordsLocalCurvilinearToFloor &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.w_mat);
        if (i == 1)
          return py::cast(s.global_position);
        throw py::index_error();
      });
  m.def(
      "coords_local_curvilinear_to_floor",
      &Bmad::coords_local_curvilinear_to_floor,
      py::arg("local_position"),
      py::arg("ele"),
      py::arg("in_body_frame") = py::none(),
      py::arg("calculate_angles") = py::none(),
      py::arg("relative_to") = py::none(),
      R"""(Parameters
----------
local_position : FloorPositionStruct
    Floor position in local curvilinear coordinates, with .r = [x, y, z_local] where z_local is wrt the
    entrance end of the element except if relative_to = downstream_end$. In this case, z_local is a distance
    -ele.value(l$)
ele : EleStruct
    element that local_position coordinates are relative to.
in_body_frame : bool, optional
    True => local_position is in ele body frame and includes misalignments.
w_mat : float
    W matrix at z, to transform vectors. v_global     = w_mat . v_local/body
calculate_angles : bool, optional
    calculate angles for global_position Default: True.
relative_to : int, optional
    not_set$ (default), upstream_end$, or downstream_end$. Force which end is used for z = 0. If
    upstream_end$, local_position.r(3) is relative to the upstream end which will not be the entrance end if
    ele.orientation = -1.
global_position : FloorPositionStruct
    Position in global coordinates.
)"""
  );
  m.def(
      "coords_relative_to_floor",
      &Bmad::coords_relative_to_floor,
      py::arg("floor0"),
      py::arg("dr"),
      py::arg("theta") = py::none(),
      py::arg("phi") = py::none(),
      py::arg("psi") = py::none(),
      R"""(Parameters
----------
floor0 : FloorPositionStruct
    Initial reference frame.
dr : float
    (x, y, z) positional shift of the reference frame.
theta : unknown, optional
    Angular shift of the reference frame. See the Bmad manual on the Global Coordinate system for more
    details. All angles must either be absent or present.
phi : unknown, optional
    Angular shift of the reference frame. See the Bmad manual on the Global Coordinate system for more
    details. All angles must either be absent or present.
psi : unknown, optional
    Angular shift of the reference frame. See the Bmad manual on the Global Coordinate system for more
    details. All angles must either be absent or present.
floor1 : FloorPositionStruct
    Shifted reference frame.
)"""
  );
  m.def(
      "coulombfun",
      &Bmad::coulombfun,
      py::arg("u"),
      py::arg("v"),
      py::arg("w"),
      py::arg("gam"),
      py::arg("res"),
      R"""(Parameters
----------
u : 
v : 
w : 
gam : 
res : 
)"""
  );
  m.def(
      "create_concatenated_wall3d",
      &Bmad::create_concatenated_wall3d,
      py::arg("lat"),
      py::arg("err"),
      R"""(Subroutine create_concatenated_wall3d (lat)

Routine to concatinate lat%branch(i)ele(:)%wall3d%section(:) arrays into
one lat%branch(i)%wall3d%section(:) array.

Exceptions: capillary and aperture elements do not have their walls included.

Module needed:
  use wall3d_mod

Parameters
----------
lat : LatStruct
    lattice
    This parameter is an input/output and is modified in-place. As an output: Lattice

Returns
-------
err_flag : bool
    Set True if there is an error, false otherwise.
)"""
  );
  py::class_<Bmad::CreateElementSlice, std::unique_ptr<Bmad::CreateElementSlice>>(
      m,
      "CreateElementSlice",
      "create_element_slice return type"
  )
      .def_readonly("sliced_ele", &Bmad::CreateElementSlice::sliced_ele)
      .def_readonly("err_flag", &Bmad::CreateElementSlice::err_flag)
      .def("__len__", [](const Bmad::CreateElementSlice &) { return 2; })
      .def("__getitem__", [](const Bmad::CreateElementSlice &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.sliced_ele);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "create_element_slice",
      &Bmad::create_element_slice,
      py::arg("ele_in"),
      py::arg("l_slice"),
      py::arg("offset"),
      py::arg("param"),
      py::arg("include_upstream_end"),
      py::arg("include_downstream_end"),
      py::arg("old_slice") = py::none(),
      py::arg("orb_in") = py::none(),
      R"""(Parameters
----------
sliced_ele : EleStruct
    Sliced_ele element with appropriate values set.
ele_in : EleStruct
    Original element to slice
l_slice : float
    Length of the slice
offset : float
    Offset of entrance end of sliced_ele from entrance end of ele_in.
param : LatParamStruct
    lattice paramters.
include_upstream_end : bool
    Sliced_ele contains the ele's entrance end?
include_downstream_end : bool
    Sliced_ele contains the ele's exit end?
err_flag : bool
    Set True if there is an error. False otherwise.
old_slice : EleStruct, optional
    Previous slice or, if offset = 0, the previous element. If present this saves computation time of the
    reference energy and time at the start of the present slice. Also makes the ref energy continuous (there
    can be some small differences when
orb_in : CoordStruct, optional
    Incoming orbit if calling routine is doing tracking through the slice. This is used when old_slice is not
    present and there may be an adjustment needed to the orbit ref energy (EG space charge tracking does not
    keep track of ref energy through an lcavity).
)"""
  );
  m.def(
      "create_field_overlap",
      &Bmad::create_field_overlap,
      py::arg("lat"),
      py::arg("lord_name"),
      py::arg("slave_name"),
      R"""(Parameters
----------
lat : LatStruct
    Lattice
lord_name : unknown
    Name of the element with a field extending beyound it's bounds.
slave_name : unknown
    Name of the element the lord's field overlaps.
err_flag : bool
    Set true if there is a problem (like no elements found).
)"""
  );
  m.def(
      "create_girder",
      &Bmad::create_girder,
      py::arg("lat"),
      py::arg("ix_girder"),
      py::arg("contrl"),
      py::arg("girder_info"),
      py::arg("err_flag"),
      R"""(Parameters
----------
lat : LatStruct
    Lat to modify.
    This parameter is an input/output and is modified in-place. As an output: Modified lattice.
ix_girder : int
    Index of girder element.
contrl : ControlStruct
    Array of elements that are supported by the girder.
girder_info : EleStruct
    Element containing attributes to be transfered to the Girder element: girder_info.name girder_info.alias
    girder_info.descrip girder_info.value(:)
err_flag : 
)"""
  );
  m.def(
      "create_group",
      &Bmad::create_group,
      py::arg("lord"),
      py::arg("contrl"),
      py::arg("err"),
      R"""(Parameters
----------
lord : EleStruct
    Group element.
    This parameter is an input/output and is modified in-place. As an output: Modified group elment
contrl : ControlStruct
    control info. 1 element for each slave.
err : bool
    Set True if an attribute is not free to be controlled.
)"""
  );
  m.def(
      "create_lat_ele_nametable",
      &Bmad::create_lat_ele_nametable,
      py::arg("lat"),
      py::arg("nametable"),
      R"""(Parameters
----------
lat : LatStruct
    Lattice.
nametable : NametableStruct
    Nametable of the elment names
)"""
  );
  m.def(
      "create_overlay",
      &Bmad::create_overlay,
      py::arg("lord"),
      py::arg("contrl"),
      py::arg("err"),
      R"""(Parameters
----------
lord : EleStruct
    Overlay element.
    This parameter is an input/output and is modified in-place. As an output: Modified overlay elment
contrl : ControlStruct
    control info. 1 element for each slave.
err : bool
    Set True if an attribute is not free to be controlled.
)"""
  );
  py::class_<Bmad::CreatePlanarWigglerModel, std::unique_ptr<Bmad::CreatePlanarWigglerModel>>(
      m,
      "CreatePlanarWigglerModel",
      "create_planar_wiggler_model return type"
  )
      .def_readonly("lat", &Bmad::CreatePlanarWigglerModel::lat)
      .def_readonly("err_flag", &Bmad::CreatePlanarWigglerModel::err_flag)
      .def("__len__", [](const Bmad::CreatePlanarWigglerModel &) { return 2; })
      .def("__getitem__", [](const Bmad::CreatePlanarWigglerModel &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.lat);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "create_planar_wiggler_model",
      &Bmad::create_planar_wiggler_model,
      py::arg("wiggler_in"),
      py::arg("print_err") = py::none(),
      R"""(Subroutine create_planar_wiggler_model (wiggler_in, lat, err_flag, print_err)

Routine to create series of bend and drift elements to serve as a replacement
model for a planar wiggler.

This routine is helpful for translating bmad lattices to a language that does not
implement the Bmad wiggler model.

This routine uses the mrqmin nonlinear optimizer to vary the parameters in the wiggler
model to match:
  Integral g^2 (I_2 radiation integral)
  Integral g^3 (I_3 radiation integral)
  Transfer matrix.
Also the endding horizontal transverse offset of the reference orbit (floor%r(1)) is
matched to zero.

Note: The resulting model does not have the vertical cubic nonlinearity that
the actual wiggler has.

Parameters
----------
wiggler : EleStruct
    Planar model wiggler to match to.
wig_model_com : WigglerModelingCommonStruct
    Global variable that can be used
to set weights and step sizes for the optimization. : 
print_err : bool, optional
    If True (default) print an error message if there is an error.

Returns
-------
lat : LatStruct
    Lattice containing the wiggler model
%ele : 
    Array of bends and drifts.
%n_ele_track : 
    Number of elements in the model.
err_flag : bool
    Set True if there is an error.
)"""
  );
  m.def(
      "create_ramper",
      &Bmad::create_ramper,
      py::arg("lord"),
      py::arg("contrl"),
      py::arg("err"),
      R"""(Parameters
----------
lord : EleStruct
    Ramper element.
    This parameter is an input/output and is modified in-place. As an output: Modified ramper elment
contrl : ControlStruct
    control info. 1 element for each slave.
err : bool
    Set True if an attribute is not free to be controlled.
)"""
  );
  m.def(
      "create_sol_quad_model",
      &Bmad::create_sol_quad_model,
      py::arg("sol_quad"),
      py::arg("lat"),
      R"""(Subroutine create_sol_quad_model (sol_quad, lat)

Routine to create series of solenoid and quadrupole elements to serve as a replacement
model for a sol_quad element.

This routine is helpful for translating bmad lattices to a language that does not
implement a combination solenoid/quadrupole.

Not yet implemented!

)"""
  );
  m.def(
      "create_unique_ele_names",
      &Bmad::create_unique_ele_names,
      py::arg("lat"),
      py::arg("key"),
      py::arg("suffix"),
      R"""(Parameters
----------
lat : LatStruct
    Lattice holding the elements.
    This parameter is an input/output and is modified in-place. As an output: Lattice with names made unique.
key : int
    Class key of elements to consider.
suffix : unknown
    Suffix string. Must have a single "?" character.
)"""
  );
  m.def(
      "create_wiggler_cartesian_map",
      &Bmad::create_wiggler_cartesian_map,
      py::arg("ele"),
      R"""(Parameters
----------
ele : EleStruct
    Wiggler or undulator element.
cart_map : CartesianMapStruct
    Cartesian map.
)"""
  );
  m.def(
      "crystal_attribute_bookkeeper",
      &Bmad::crystal_attribute_bookkeeper,
      py::arg("ele"),
      R"""(Parameters
----------
ele : EleStruct
    Crystal element. .value(bragg_angle_in$) .value(bragg_angle_out$) .value(tilt_corr$) ... etc.
)"""
  );
  m.def(
      "crystal_h_misalign",
      &Bmad::crystal_h_misalign,
      py::arg("ele"),
      py::arg("orbit"),
      py::arg("h_vec"),
      R"""(Subroutine crystal_h_misalign (ele, orbit, h_vec)

Routine reorient the crystal H vector due to local imperfections in the crystal lattice.

Parameters
----------
ele : EleStruct
    Crystal element
orbit : CoordStruct
    Photon position at crystal surface.
h_vec : float
    H vector before misalignment.
    This parameter is an input/output and is modified in-place. As an output: H vector after misalignment.
)"""
  );
  m.def(
      "crystal_type_to_crystal_params",
      &Bmad::crystal_type_to_crystal_params,
      py::arg("ele"),
      R"""(Subroutine crystal_type_to_crystal_params (ele, err_flag)

Routine to set the crystal parameters based upon the crystal type.

Crystal types are of the form:
  "ZZZ(ijk)"
Where "ZZZ" is the atomic formula of the crystal material and "ijk" is the reciprical lattice
vetor specifying the diffraction plans.

Parameters
----------
ele : EleStruct
    Crystal element.
    This parameter is an input/output and is modified in-place. As an output: Crystal element with computed
    parameter..
%component_name : unknown
    Crystal type name. Assumed upper case.
A blank name is not an error and results in nothing set. : 
%value : 
    Photon energy in eV.

Returns
-------
err_flag : bool
    Set True if crystal type is unrecognized. False otherwise.
)"""
  );
  m.def(
      "custom_attribute_ubound_index",
      &Bmad::custom_attribute_ubound_index,
      py::arg("ele_class"),
      R"""(Function custom_attribute_ubound_index(ele_class) result (ix_ubound)

Routine to return, for a given element class, the upper bound index for the ele%custom(:)
array which is needed to accomodate the registered custom attributes for that class.

Parameters
----------
ele_class : int
    Element class (key). EG: quadrupole$, etc.

Returns
-------
ix_ubound : int
    Maximum index needed.
)"""
  );
}
