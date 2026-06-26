#include "pybmad/generated/Bmad_routines_c.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

PyChromCalc python_chrom_calc(
    LatStruct &lat,
    double delta_e,
    std::optional<double> pz = std::nullopt,
    std::optional<int> ix_branch = std::nullopt,
    CoordStruct *orb0 = nullptr
) {
  auto _result = Bmad::chrom_calc(lat, delta_e, pz, ix_branch, ptr_to_opt_ref(orb0));
  auto py_result{PyChromCalc{_result, delta_e}};
  return py_result;
}
PyChromTune python_chrom_tune(
    LatStruct &lat,
    double delta_e,
    double target_x,
    double target_y,
    double err_tol
) {
  auto _result = Bmad::chrom_tune(lat, delta_e, target_x, target_y, err_tol);
  auto py_result{PyChromTune{_result, delta_e}};
  return py_result;
}

void init_Bmad_routines_c(nb::module_ &m) {
  m.def(
      "c_to_cbar",
      &Bmad::c_to_cbar,
      nb::arg("ele"),
      R"""(Wrapper for Fortran routine c_to_cbar

Parameters
----------
ele : EleStruct
    Element with C matrix and Twiss parameters.

Returns
-------
cbar_mat : 2D array of float (shape: 2,2)
    Cbar matrix.
)"""
  );
  nb::class_<Bmad::CalcBunchParams>(m, "CalcBunchParams", "calc_bunch_params return type")
      .def_ro("bunch_params", &Bmad::CalcBunchParams::bunch_params)
      .def_ro("error", &Bmad::CalcBunchParams::error)
      .def_ro("n_mat", &Bmad::CalcBunchParams::n_mat)
      .def("__len__", [](const Bmad::CalcBunchParams &) { return 3; })
      .def("__getitem__", [](const Bmad::CalcBunchParams &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.bunch_params);
        if (i == 1)
          return nb::cast(s.error);
        if (i == 2)
          return nb::cast(s.n_mat);
        throw nb::index_error();
      });
  m.def(
      "calc_bunch_params",
      [](BunchStruct &bunch,
         std::optional<bool> print_err,
         std::optional<bool> is_time_coords,
         EleStruct *ele) {
        auto fn = static_cast<
            Bmad::
                CalcBunchParams (*)(BunchStruct &, std::optional<bool>, std::optional<bool>, optional_ref<EleStruct>)>(
            &Bmad::calc_bunch_params
        );
        return fn(bunch, print_err, is_time_coords, ptr_to_opt_ref(ele));
      },
      nb::arg("bunch"),
      nb::arg("print_err") = nb::none(),
      nb::arg("is_time_coords") = nb::none(),
      nb::arg("ele") = nb::none(),
      R"""(Finds all bunch parameters defined in bunch_params_struct, both normal-mode
and projected. Projected parameters are found purely from the geometrical
distribution of the beam. Normal-Mode parameters are found using the method
developed in:
  "Alternate approach to general coupled linear optics"
   A. Wolski, PRST AB 9, 024001 (2006)

Note: If less than two particle remain then the various parameters will be
set to zero.

Parameters
----------
bunch : BunchStruct
    Bunch_struct

print_err : bool, optional
    If present and False then suppress "no eigen-system found" messages.

is_time_coords : bool, optional
    Are particle coords using time coords. Default is False.

ele : EleStruct, optional
    Element being tracked through. Must be present if is_time_coords = True.

Returns
-------
bunch_params : BunchParamsStruct

error : bool
    Set True if there is an error.

n_mat : 2D array of float (shape: 6,6), optional
    N matrix defined in Wolski Eq 44 and used to convert from action-angle coords to lab coords (Wolski Eq
    51.).
)"""
  );
  m.def(
      "calc_bunch_params_slice",
      [](BunchStruct &bunch,
         BunchParamsStruct &bunch_params,
         int plane,
         double slice_center,
         double slice_spread,
         std::optional<bool> print_err,
         std::optional<bool> is_time_coords,
         EleStruct *ele) {
        auto fn = static_cast<
            bool (*)(BunchStruct &, BunchParamsStruct &, int, double, double, std::optional<bool>, std::optional<bool>, optional_ref<EleStruct>)>(
            &Bmad::calc_bunch_params_slice
        );
        return fn(
            bunch,
            bunch_params,
            plane,
            slice_center,
            slice_spread,
            print_err,
            is_time_coords,
            ptr_to_opt_ref(ele)
        );
      },
      nb::arg("bunch"),
      nb::arg("bunch_params"),
      nb::arg("plane"),
      nb::arg("slice_center"),
      nb::arg("slice_spread"),
      nb::arg("print_err") = nb::none(),
      nb::arg("is_time_coords") = nb::none(),
      nb::arg("ele") = nb::none(),
      R"""(Finds bunch parameters for a slice of the beam.

Parameters
----------
bunch : BunchStruct
    bunch_struct

plane : int
    plane to slice through (x$, px$, & etc...)

slice_center : float
    Center to take slice about

slice_spread : float
    +/- spread in slice about center.

print_err : bool, optional
    If present and False then suppress "no eigen-system found" messages.

is_time_coords : bool, optional
    Default is False. If True, input bunch is using time coordinates in which case there will be a conversion
    to s-coords before bunch_params are computed.

ele : EleStruct, optional
    Element being tracked through. Must be present if is_time_coords = True.

Returns
-------
err : bool
    Set True if there is an error.
)"""
  );
  m.def(
      "calc_bunch_params_z_slice",
      [](BunchStruct &bunch,
         BunchParamsStruct &bunch_params,
         FixedArray1D<Real, 2> slice_bounds,
         std::optional<bool> print_err,
         std::optional<bool> is_time_coords,
         EleStruct *ele) {
        auto fn = static_cast<
            bool (*)(BunchStruct &, BunchParamsStruct &, FixedArray1D<Real, 2>, std::optional<bool>, std::optional<bool>, optional_ref<EleStruct>)>(
            &Bmad::calc_bunch_params_z_slice
        );
        return fn(
            bunch,
            bunch_params,
            slice_bounds,
            print_err,
            is_time_coords,
            ptr_to_opt_ref(ele)
        );
      },
      nb::arg("bunch"),
      nb::arg("bunch_params"),
      nb::arg("slice_bounds"),
      nb::arg("print_err") = nb::none(),
      nb::arg("is_time_coords") = nb::none(),
      nb::arg("ele") = nb::none(),
      R"""(Finds bunch parameters for a slice of the beam.

The slice is specified in terms of percentage of particles ordered by z-position.
For example, slice_bounds = [0.0, 0.5] specifies the trailing half of the bunch

Parameters
----------
bunch : BunchStruct
    bunch_struct

slice_bounds : 1D array of float (shape: 2)
    Slice bounds in percentage of particles ordered by z-position. 0.0 is the back of the bunch and 1.0 is the
    front of the bunch.

print_err : bool, optional
    If present and False then suppress "no eigen-system found" messages.

is_time_coords : bool, optional
    Default is False. If True, input bunch is using time coordinates in which case there will be a conversion
    to s-coords before bunch_params are computed.

ele : EleStruct, optional
    Element being tracked through. Must be present if is_time_coords = True.

Returns
-------
err : bool
    Set True if there is an error.
)"""
  );
  m.def(
      "calc_bunch_sigma_matrix_etc",
      [](CoordStructArray1D particle,
         FArray1D<Real> &charge,
         std::optional<bool> is_time_coords,
         EleStruct *ele) {
        auto fn = static_cast<
            BunchParamsStruct (*)(CoordStructArray1D, FArray1D<Real> &, std::optional<bool>, optional_ref<EleStruct>)>(
            &Bmad::calc_bunch_sigma_matrix_etc
        );
        return fn(particle, charge, is_time_coords, ptr_to_opt_ref(ele));
      },
      nb::arg("particle"),
      nb::arg("charge"),
      nb::arg("is_time_coords") = nb::none(),
      nb::arg("ele") = nb::none(),
      R"""(Routine to find the sigma matrix elements of a particle distribution.

Parameters
----------
particle : 1D array of CoordStruct
    Array of particles.

charge : 1D array of float
    Particle charge or photon intensity.

Returns
-------
bunch_params : BunchParamsStruct
    Bunch parameters. .sigma(6,6) .centroid.vec(6) .centroid.p0c .rel_max(6) .rel_min(6)
)"""
  );
  nb::class_<Bmad::CalcEmittancesAndTwissFromSigmaMatrix>(
      m,
      "CalcEmittancesAndTwissFromSigmaMatrix",
      "calc_emittances_and_twiss_from_sigma_matrix return type"
  )
      .def_ro("bunch_params", &Bmad::CalcEmittancesAndTwissFromSigmaMatrix::bunch_params)
      .def_ro("error", &Bmad::CalcEmittancesAndTwissFromSigmaMatrix::error)
      .def_ro("n_mat", &Bmad::CalcEmittancesAndTwissFromSigmaMatrix::n_mat)
      .def("__len__", [](const Bmad::CalcEmittancesAndTwissFromSigmaMatrix &) { return 3; })
      .def(
          "__getitem__",
          [](const Bmad::CalcEmittancesAndTwissFromSigmaMatrix &s, int i) -> nb::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return nb::cast(s.bunch_params);
            if (i == 1)
              return nb::cast(s.error);
            if (i == 2)
              return nb::cast(s.n_mat);
            throw nb::index_error();
          }
      );
  m.def(
      "calc_emittances_and_twiss_from_sigma_matrix",
      &Bmad::calc_emittances_and_twiss_from_sigma_matrix,
      nb::arg("sigma_mat"),
      nb::arg("print_err") = nb::none(),
      R"""(Routine to calc emittances and Twiss function from a beam sigma matrix.
See: Andy Wolski "Alternative approach to general coupled linear optics".

Parameters
----------
sigma_mat : 2D array of float (shape: 6,6)
    Sigma matrix.

print_err : bool, optional
    If present and False then suppress "no eigen-system found" messages.

Returns
-------
bunch_params : BunchParamsStruct
    Holds Twiss and emittance info.

error : bool
    Set True if there is an error. Can happen if the emittance of a mode is zero.

n_mat : 2D array of float (shape: 6,6), optional
    N matrix defined in Wolski Eq 44 and used to convert from action-angle coords to lab coords (Wolski Eq
    51.).
)"""
  );
  m.def(
      "calc_spin_params",
      &Bmad::calc_spin_params,
      nb::arg("bunch"),
      R"""(Rotine to calculate spin averages

Parameters
----------
bunch : BunchStruct
    Bunch of spins

Returns
-------
bunch_params : BunchParamsStruct
    Structure holding average
)"""
  );
  m.def(
      "calc_super_slave_key",
      &Bmad::calc_super_slave_key,
      nb::arg("lord1"),
      nb::arg("lord2"),
      nb::arg("create_jumbo_slave") = nb::none(),
      R"""(Wrapper for Fortran routine calc_super_slave_key

Parameters
----------
lord1 : EleStruct
    First slave. .key .sub_key

lord2 : EleStruct
    Second slave. .key .sub_key

create_jumbo_slave : bool, optional
    If True then slave.key will be set to em_field. Default is False.

Returns
-------
slave : EleStruct
    Super_slave element.
)"""
  );
  nb::class_<Bmad::CalcWallRadius>(m, "CalcWallRadius", "calc_wall_radius return type")
      .def_ro("r_wall", &Bmad::CalcWallRadius::r_wall)
      .def_ro("dr_dtheta", &Bmad::CalcWallRadius::dr_dtheta)
      .def_ro("ix_vertex", &Bmad::CalcWallRadius::ix_vertex)
      .def("__len__", [](const Bmad::CalcWallRadius &) { return 3; })
      .def("__getitem__", [](const Bmad::CalcWallRadius &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.r_wall);
        if (i == 1)
          return nb::cast(s.dr_dtheta);
        if (i == 2)
          return nb::cast(s.ix_vertex);
        throw nb::index_error();
      });
  m.def(
      "calc_wall_radius",
      &Bmad::calc_wall_radius,
      nb::arg("v"),
      nb::arg("cos_ang"),
      nb::arg("sin_ang"),
      R"""(Routine to calculate the wall radius at a given angle for a given cross-section
Additionally, the transverse directional derivative is calculated.

Module needed:
  use wall3d_mod

Parameters
----------
v : 1D array of Wall3dVertexStruct
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

ix_vertex : int, optional
    Wall at given angle is between v(ix_vertex-1) and v(ix_vertex). If ix_vertex = 1 then Wall at given angle
    is between v(N) and v(1) where N = size(v).
)"""
  );
  m.def(
      "calc_z_tune",
      &Bmad::calc_z_tune,
      nb::arg("branch"),
      R"""(Wrapper for Fortran routine calc_z_tune

Parameters
----------
branch : BranchStruct
    Lattice branch
    This parameter is an input/output and is modified in-place.
    As an output, branch: Synchrotron tune (radians). If unstable tune = 0.
    This parameter is an input/output and is modified in-place.
    As an output, branch: Is the mode stable? If no rf then tune is zero but is stable.
    This parameter is an input/output and is modified in-place.
    As an output, branch: 6x6 1-turn matrix.
)"""
  );
  m.def(
      "canonical_to_angle_coords",
      &Bmad::canonical_to_angle_coords,
      nb::arg("orbit"),
      nb::arg("coord_type") = nb::none(),
      R"""(Wrapper for Fortran routine canonical_to_angle_coords

Parameters
----------
orbit : CoordStruct
    Orbit in canonical coordinates.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Orbit in angular coordinates.

coord_type : str, optional
    Angular coordinates type '' (default): (x, x' = dx/ds, y, y' = dy/ds, z, pz) 'ZGOUBI':     (x, x' = dx/ds,
    y, y' = dy/ds, dt = -z / (beta * c), pz)
)"""
  );
  m.def(
      "cbar_to_c",
      &Bmad::cbar_to_c,
      nb::arg("cbar_mat"),
      nb::arg("a"),
      nb::arg("b"),
      R"""(Wrapper for Fortran routine cbar_to_c

Parameters
----------
cbar_mat : 2D array of float (shape: 2,2)
    Cbar matrix.

a : TwissStruct
    a-mode Twiss parameters

b : TwissStruct
    b-mode Twiss parameters

Returns
-------
c_mat : 2D array of float (shape: 2,2)
    C matrix.
)"""
  );
  m.def(
      "check_aperture_limit",
      [](CoordStruct &orb,
         EleStruct &ele,
         int particle_at,
         LatParamStruct &param,
         CoordStruct *old_orb,
         std::optional<bool> check_momentum) {
        auto fn = static_cast<
            void (*)(CoordStruct &, EleStruct &, int, LatParamStruct &, optional_ref<CoordStruct>, std::optional<bool>)>(
            &Bmad::check_aperture_limit
        );
        return fn(orb, ele, particle_at, param, ptr_to_opt_ref(old_orb), check_momentum);
      },
      nb::arg("orb"),
      nb::arg("ele"),
      nb::arg("particle_at"),
      nb::arg("param"),
      nb::arg("old_orb") = nb::none(),
      nb::arg("check_momentum") = nb::none(),
      R"""(Wrapper for Fortran routine check_aperture_limit

Parameters
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
    Old coordinates at last check. Needed if ele.aperture_at = wall_transition$. If not present then wall
    transitions will be ignored.

check_momentum : bool, optional
    If present and false then checking of p_x and p_y will be disabled.
)"""
  );
  m.def(
      "check_controller_controls",
      &Bmad::check_controller_controls,
      nb::arg("ele_key"),
      nb::arg("contrl"),
      nb::arg("name"),
      R"""(Wrapper for Fortran routine check_controller_controls

Parameters
----------
ele_key : int
    Element type. overlay$, etc.

contrl : 1D array of ControlStruct
    control info. 1 element for each slave.

name : str
    Lord name. Used for error reporting.

Returns
-------
err : bool
    Set true if there is a problem. False otherwise.
)"""
  );
  m.def(
      "check_for_superimpose_problem",
      [](BranchStruct &branch, EleStruct &super_ele, bool err_flag, bool wrap, EleStruct *ref_ele) {
        auto fn =
            static_cast<void (*)(BranchStruct &, EleStruct &, bool, bool, optional_ref<EleStruct>)>(
                &Bmad::check_for_superimpose_problem
            );
        return fn(branch, super_ele, err_flag, wrap, ptr_to_opt_ref(ref_ele));
      },
      nb::arg("branch"),
      nb::arg("super_ele"),
      nb::arg("err_flag"),
      nb::arg("wrap"),
      nb::arg("ref_ele") = nb::none(),
      R"""(Subroutine to check if there is a problem superimposing an element when there is multipass.
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
  nb::class_<Bmad::CheckIfSInBounds>(m, "CheckIfSInBounds", "check_if_s_in_bounds return type")
      .def_ro("err_flag", &Bmad::CheckIfSInBounds::err_flag)
      .def_ro("translated_s", &Bmad::CheckIfSInBounds::translated_s)
      .def("__len__", [](const Bmad::CheckIfSInBounds &) { return 2; })
      .def("__getitem__", [](const Bmad::CheckIfSInBounds &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.translated_s);
        throw nb::index_error();
      });
  m.def(
      "check_if_s_in_bounds",
      &Bmad::check_if_s_in_bounds,
      nb::arg("branch"),
      nb::arg("s"),
      nb::arg("print_err") = nb::none(),
      R"""(Wrapper for Fortran routine check_if_s_in_bounds

Parameters
----------
branch : BranchStruct
    Branch

s : float
    longitudinal position in the given branch.

print_err : bool, optional
    Print error message if there is an error? Default is True.

Returns
-------
err_flag : bool
    Set True if s position is out-of-bounds. False otherwise.

translated_s : float, optional
    position translated to the range [0, branch_length]
)"""
  );
  nb::class_<Bmad::ChooseQuadsForSetTune>(
      m,
      "ChooseQuadsForSetTune",
      "choose_quads_for_set_tune return type"
  )
      .def_ro("dk1", &Bmad::ChooseQuadsForSetTune::dk1)
      .def_ro("eles", &Bmad::ChooseQuadsForSetTune::eles)
      .def_ro("err_flag", &Bmad::ChooseQuadsForSetTune::err_flag)
      .def("__len__", [](const Bmad::ChooseQuadsForSetTune &) { return 3; })
      .def("__getitem__", [](const Bmad::ChooseQuadsForSetTune &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.dk1);
        if (i == 1)
          return nb::cast(s.eles);
        if (i == 2)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "choose_quads_for_set_tune",
      &Bmad::choose_quads_for_set_tune,
      nb::arg("branch"),
      nb::arg("mask") = nb::none(),
      R"""(Wrapper for Fortran routine choose_quads_for_set_tune

Parameters
----------
branch : BranchStruct
    Lattice branch.

mask : str, optional
    If present, assign weight of zero for all quads that do not match. That is, no variation for matching
    quads.

Returns
-------
dk1 : 1D array of float
    Weights for the quadrupoles. All values will be +1 or -1.

eles : 1D array of ElePointerStruct
    eles(i).ele points to element with dk1(i) weight.

err_flag : bool, optional
    Set True if there is not one quad with positive dk1 and one quad with negative dk1.
)"""
  );
  nb::class_<PyChromCalc>(m, "ChromCalc", "chrom_calc return type")
      .def_ro("chrom_a", &PyChromCalc::chrom_a)
      .def_ro("chrom_b", &PyChromCalc::chrom_b)
      .def_ro("err_flag", &PyChromCalc::err_flag)
      .def_ro("low_E_lat", &PyChromCalc::low_E_lat)
      .def_ro("high_E_lat", &PyChromCalc::high_E_lat)
      .def_ro("low_E_orb", &PyChromCalc::low_E_orb)
      .def_ro("high_E_orb", &PyChromCalc::high_E_orb)
      .def_ro("delta_e", &PyChromCalc::delta_e)
      .def("__len__", [](const PyChromCalc &) { return 8; })
      .def("__getitem__", [](const PyChromCalc &s, int i) -> nb::object {
        if (i < 0)
          i += 8;
        if (i == 0)
          return nb::cast(s.chrom_a);
        if (i == 1)
          return nb::cast(s.chrom_b);
        if (i == 2)
          return nb::cast(s.err_flag);
        if (i == 3)
          return nb::cast(s.low_E_lat);
        if (i == 4)
          return nb::cast(s.high_E_lat);
        if (i == 5)
          return nb::cast(s.low_E_orb);
        if (i == 6)
          return nb::cast(s.high_E_orb);
        if (i == 7)
          return nb::cast(s.delta_e);
        throw nb::index_error();
      });
  m.def(
      "chrom_calc",
      &python_chrom_calc,
      nb::arg("lat"),
      nb::arg("delta_e"),
      nb::arg("pz") = nb::none(),
      nb::arg("ix_branch") = nb::none(),
      nb::arg("orb0") = nb::none(),
      R"""(Wrapper for Fortran routine chrom_calc

Parameters
----------
lat : LatStruct
    Lat

delta_e : float
    +/- Delta energy used for the calculation. Notice that the energy difference between high and low is 2 *
    delta_e. If 0 then default of 1.0d-4 is used.
    This parameter is an input/output and is modified in-place.
    As an output, delta_e: Set to 1.0d-4 if on input DELTA_E =< 0.

pz : float, optional
    reference momentum about which to calculate. Default is 0.

ix_branch : int, optional
    Index of the lattice branch to use. Default is 0.

orb0 : CoordStruct, optional
    On-energy orbit at start (fixer point). Default is the branch.particle_start. Only needed if lattice
    branch has an open geometry.

Returns
-------
delta_e : float
    +/- Delta energy used for the calculation. Notice that the energy difference between high and low is 2 *
    delta_e. If 0 then default of 1.0d-4 is used.
    This parameter is an input/output and is modified in-place.
    As an output, delta_e: Set to 1.0d-4 if on input DELTA_E =< 0.

chrom_a : float
    a-mode chromaticity.

chrom_b : float
    b-mode chromaticity.

err_flag : bool, optional
    Set true if there is an error. False otherwise.

low_E_lat : LatStruct, optional
    Lattice with RF off and matrices computed at E_lat +pz - delta_e

high_E_lat : LatStruct, optional
    Lattice with RF off and matrices computed at E_lat +pz + delta_e

low_E_orb : 1D array of CoordStruct, optional
    Orbit computed at E_lat + pz - delta_e.

high_E_orb : 1D array of CoordStruct, optional
    Orbit computed at E_lat + pz + delta_e.
)"""
  );
  nb::class_<PyChromTune>(m, "ChromTune", "chrom_tune return type")
      .def_ro("err_flag", &PyChromTune::err_flag)
      .def_ro("delta_e", &PyChromTune::delta_e)
      .def("__len__", [](const PyChromTune &) { return 2; })
      .def("__getitem__", [](const PyChromTune &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.delta_e);
        throw nb::index_error();
      });
  m.def(
      "chrom_tune",
      &python_chrom_tune,
      nb::arg("lat"),
      nb::arg("delta_e"),
      nb::arg("target_x"),
      nb::arg("target_y"),
      nb::arg("err_tol"),
      R"""(Wrapper for Fortran routine chrom_tune

Parameters
----------
lat : LatStruct
    Lat to use,
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lat with sextupole set

delta_e : float
    Delta energy used for the calculation. If 0 then default of 1.0d-4 is used.
    This parameter is an input/output and is modified in-place.
    As an output, delta_e: Set to 1.0d-4 if on input DELTA_E =< 0.

target_x : float
    Target X Chromaticity

target_y : float
    Target Y Chromaticity

err_tol : float
    Max allowable Error: Error = | X_Target - X_Actual | + | Y_Target -Y_Actual | A good number is: err_tol =
    0.05_rp

Returns
-------
delta_e : float
    Delta energy used for the calculation. If 0 then default of 1.0d-4 is used.
    This parameter is an input/output and is modified in-place.
    As an output, delta_e: Set to 1.0d-4 if on input DELTA_E =< 0.

err_flag : bool
    .false. if match successful, .true. if failed Fails if takes longer than 100 iterations. If it fails the
    sextupoles are set to the last value calculated. Note: This subroutine assumes the Twiss parameters have
    been computed.
)"""
  );
  m.def(
      "classical_radius",
      &Bmad::classical_radius,
      nb::arg("species"),
      R"""(Wrapper for Fortran routine classical_radius

Parameters
----------
species : int
    Species of particle.

Returns
-------
radius : float
    Classical radius.
)"""
  );
  m.def(
      "clear_lat_1turn_mats",
      &Bmad::clear_lat_1turn_mats,
      R"""(Wrapper for Fortran routine clear_lat_1turn_mats

Returns
-------
lat : LatStruct
    Lat with 1-turn matrices cleared.
)"""
  );
  m.def(
      "clear_taylor_maps_from_elements",
      &Bmad::clear_taylor_maps_from_elements,
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine clear_taylor_maps_from_elements

Parameters
----------
lat : LatStruct
    Lattice
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lattice with all maps cleared
)"""
  );
  m.def(
      "closed_orbit_calc",
      &Bmad::closed_orbit_calc,
      nb::arg("lat"),
      nb::arg("closed_orb"),
      nb::arg("i_dim") = nb::none(),
      nb::arg("direction") = nb::none(),
      nb::arg("ix_branch") = nb::none(),
      nb::arg("print_err") = nb::none(),
      R"""(Wrapper for Fortran routine closed_orbit_calc

Parameters
----------
lat : LatStruct
    Lat to track through.

closed_orb : 1D array of CoordStruct
    closed_orb(nt) is the initial guess where nt = 0 for direction = 1 and nt = lat.n_ele_track for direction
    = -1. Additionally, if i_dim = 4, then closed_orb(nt).vec(6) is used as the energy around which the closed
    orbit is calculated.
    This parameter is an input/output and is modified in-place.
    As an output, closed_orb: Closed orbit. closed_orb(i)

i_dim : int, optional
    Phase space dimensions to use: = 4  Transverse closed orbit at constant energy (RF off). (dE/E =
    closed_orb(0).vec(6)) = 5 Transverse closed orbit at constant energy (RF off) with the energy adjusted so
    that vec(5) is the same at the beginning and at the end. = 6 True closed orbit. Default: 4 if RF is off, 6
    if RF is on.

direction : int, optional
    Direction of tracking.

ix_branch : int, optional
    Lattice branch to find the closed orbit of. Default is 0 (main branch).

print_err : bool, optional
    Print error message if calc does not converge? Default is True. Note: Condition messages like no RF
    voltage with i_dim = 6 will always be printed.

Returns
-------
err_flag : bool, optional
    Set true if there is an error. False otherwise.
)"""
  );
  nb::class_<Bmad::ClosedOrbitFromTracking>(
      m,
      "ClosedOrbitFromTracking",
      "closed_orbit_from_tracking return type"
  )
      .def_ro("closed_orb", &Bmad::ClosedOrbitFromTracking::closed_orb)
      .def_ro("err_flag", &Bmad::ClosedOrbitFromTracking::err_flag)
      .def("__len__", [](const Bmad::ClosedOrbitFromTracking &) { return 2; })
      .def("__getitem__", [](const Bmad::ClosedOrbitFromTracking &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.closed_orb);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "closed_orbit_from_tracking",
      [](LatStruct &lat,
         int i_dim,
         std::optional<FArray1D<Real>> eps_rel,
         std::optional<FArray1D<Real>> eps_abs,
         CoordStruct *init_guess) {
        auto fn = static_cast<
            Bmad::
                ClosedOrbitFromTracking (*)(LatStruct &, int, std::optional<FArray1D<Real>>, std::optional<FArray1D<Real>>, optional_ref<CoordStruct>)>(
            &Bmad::closed_orbit_from_tracking
        );
        return fn(lat, i_dim, eps_rel, eps_abs, ptr_to_opt_ref(init_guess));
      },
      nb::arg("lat"),
      nb::arg("i_dim"),
      nb::arg("eps_rel") = nb::none(),
      nb::arg("eps_abs") = nb::none(),
      nb::arg("init_guess") = nb::none(),
      R"""(Wrapper for Fortran routine closed_orbit_from_tracking

Parameters
----------
lat : LatStruct
    Lat to track through.

i_dim : int
    = 2,4  Transverse closed orbit at constant energy. = 6    Full closed orbit using the entire transfer 6x6
    matrix.

eps_rel : 1D array of float, optional
    Relative allowed error. Default is bmad_com.rel_tol_tracking

eps_abs : 1D array of float, optional
    Absolute allowed error. Default is bmad_com.abs_tol_tracking

init_guess : CoordStruct, optional
    Starting guess for the closed orbit at the start of the lattice. Set init_guess.vec(6) to the appropriate
    value of pz when calculating off-energy orbits. If not present then the origin will be used.

Returns
-------
closed_orb : 1D array of CoordStruct
    closed orbit. This routine will allocate this array for you.

err_flag : bool, optional
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "cmplx_re_str",
      &Bmad::cmplx_re_str,
      nb::arg("cmp"),
      R"""(Wrapper for Fortran routine cmplx_re_str

Parameters
----------
cmp : complex

Returns
-------
str_out : str
)"""
  );
  m.def(
      "combine_consecutive_elements",
      &Bmad::combine_consecutive_elements,
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine combine_consecutive_elements

Parameters
----------
lat : LatStruct
    Lattice.
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lattice with elements combined.

Returns
-------
error : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "complex_taylor_clean",
      &Bmad::complex_taylor_clean,
      nb::arg("complex_taylor"),
      R"""(Wrapper for Fortran routine complex_taylor_clean

Parameters
----------
complex_taylor : ComplexTaylorStruct
)"""
  );
  m.def(
      "complex_taylor_coef",
      nb::overload_cast<ComplexTaylorStruct &, FArray1D<Int> &>(&Bmad::complex_taylor_coef),
      nb::arg("complex_taylor"),
      nb::arg("exp"),
      R"""(Function to return the coefficient for a particular complex_taylor term
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
)"""
  );
  m.def(
      "complex_taylor_coef",
      nb::overload_cast<
          ComplexTaylorStruct &,
          std::optional<int>,
          std::optional<int>,
          std::optional<int>,
          std::optional<int>,
          std::optional<int>,
          std::optional<int>,
          std::optional<int>,
          std::optional<int>,
          std::optional<int>>(&Bmad::complex_taylor_coef),
      nb::arg("complex_taylor"),
      nb::arg("i1") = nb::none(),
      nb::arg("i2") = nb::none(),
      nb::arg("i3") = nb::none(),
      nb::arg("i4") = nb::none(),
      nb::arg("i5") = nb::none(),
      nb::arg("i6") = nb::none(),
      nb::arg("i7") = nb::none(),
      nb::arg("i8") = nb::none(),
      nb::arg("i9") = nb::none(),
      R"""(Function to return the coefficient for a particular complex_taylor term
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
)"""
  );
  m.def(
      "complex_taylor_equal_complex_taylor",
      &Bmad::complex_taylor_equal_complex_taylor,
      nb::arg("complex_taylor1"),
      nb::arg("complex_taylor2"),
      R"""(Wrapper for Fortran routine complex_taylor_equal_complex_taylor

Parameters
----------
complex_taylor1 : ComplexTaylorStruct

complex_taylor2 : ComplexTaylorStruct
)"""
  );
  m.def(
      "complex_taylor_exponent_index",
      &Bmad::complex_taylor_exponent_index,
      nb::arg("expn"),
      R"""(Function to associate a unique number with a complex_taylor exponent.

The number associated with a complex_taylor_term that is used for the sort is:
    number = sum(exp(i))*10^6 + exp(6)*10^5 + ... + exp(1)*10^0
where exp(1) is the exponent for x, exp(2) is the exponent for P_x, etc.

Parameters
----------
expn : 1D array of int (shape: 6)
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
      nb::arg("complex_taylor"),
      R"""(Subroutine to make the unit complex_taylor map:
      r(out) = Map * r(in) = r(in)
)"""
  );
  nb::class_<Bmad::ComplexTaylorToMat6>(
      m,
      "ComplexTaylorToMat6",
      "complex_taylor_to_mat6 return type"
  )
      .def_ro("vec0", &Bmad::ComplexTaylorToMat6::vec0)
      .def_ro("mat6", &Bmad::ComplexTaylorToMat6::mat6)
      .def("__len__", [](const Bmad::ComplexTaylorToMat6 &) { return 2; })
      .def("__getitem__", [](const Bmad::ComplexTaylorToMat6 &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.vec0);
        if (i == 1)
          return nb::cast(s.mat6);
        throw nb::index_error();
      });
  m.def(
      "complex_taylor_to_mat6",
      &Bmad::complex_taylor_to_mat6,
      nb::arg("a_complex_taylor"),
      nb::arg("r_in"),
      nb::arg("r_out") = nb::none(),
      R"""(Subroutine to calculate, from a complex_taylor map and about some trajectory:
  The 1st order (Jacobian) transfer matrix.

Parameters
----------
a_complex_taylor : 1D array of ComplexTaylorStruct (shape: 6)
    complex_taylor map.

r_in : 1D array of complex
    Coordinates at the input.

r_out : 1D array of complex, optional
    Coordinates at output.

Returns
-------
vec0 : 1D array of complex (shape: 6)
    0th order tranfsfer map

mat6 : 2D array of complex (shape: 6,6)
    1st order transfer map (6x6 matrix).
)"""
  );
  m.def(
      "complex_taylors_equal_complex_taylors",
      &Bmad::complex_taylors_equal_complex_taylors,
      nb::arg("complex_taylor1"),
      nb::arg("complex_taylor2"),
      R"""(Wrapper for Fortran routine complex_taylors_equal_complex_taylors

Parameters
----------
complex_taylor1 : 1D array of ComplexTaylorStruct

complex_taylor2 : 1D array of ComplexTaylorStruct
)"""
  );
  m.def(
      "compute_slave_coupler",
      &Bmad::compute_slave_coupler,
      nb::arg("slave"),
      R"""(This routine is not meant for general use.
)"""
  );
  m.def(
      "compute_super_lord_s",
      &Bmad::compute_super_lord_s,
      nb::arg("ref_ele"),
      nb::arg("super_ele"),
      nb::arg("pele"),
      nb::arg("ix_insert"),
      R"""(Wrapper for Fortran routine compute_super_lord_s

Parameters
----------
ref_ele : EleStruct

super_ele : EleStruct

pele : ParserEleStruct

ix_insert : int
)"""
  );
  m.def(
      "concat_ele_taylor",
      &Bmad::concat_ele_taylor,
      nb::arg("orb_taylor"),
      nb::arg("ele"),
      nb::arg("spin_taylor") = nb::none(),
      R"""(Routine to concatinate an orbital taylor map and, optionally if present and
bmad_com%spin_tracking_on = T, a spin taylor map.

Transform:
  orb_taylor[x] -> ele_taylor(orb_taylor[x])
If ele%taylor_map_includes_offsets = True:  ele_taylor == ele%taylor
If ele%taylor_map_includes_offsets = False: ele_taylor == ele%taylor + offset corrections.

Also see: concat_taylor

Parameters
----------
orb_taylor : 1D array of TaylorStruct
    Orbital Taylor map.
    This parameter is an input/output and is modified in-place.
    As an output, orb_taylor: Concatinated orbital map

ele : EleStruct
    Element containing a Taylor map.

spin_taylor : 1D array of TaylorStruct, optional
    Spin map to propagate
    This parameter is an input/output and is modified in-place.
    As an output, spin_taylor: Concatinated spin map.

Returns
-------
err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "concat_taylor",
      &Bmad::concat_taylor,
      nb::arg("taylor1"),
      nb::arg("taylor2"),
      nb::arg("taylor3"),
      R"""(Subroutine to concatinate two taylor maps:
  taylor3[x] = taylor2(taylor1[x])

Note: In general, if taylor2 is a component of an ele_struct, use
concat_ele_taylor instead.

Parameters
----------
taylor1 : 1D array of TaylorStruct
    Taylor map.

taylor2 : 1D array of TaylorStruct
    Taylor map.

taylor3 : 1D array of TaylorStruct
    Concatinated map
)"""
  );
  nb::class_<Bmad::ConcatTransferMat>(m, "ConcatTransferMat", "concat_transfer_mat return type")
      .def_ro("mat_out", &Bmad::ConcatTransferMat::mat_out)
      .def_ro("vec_out", &Bmad::ConcatTransferMat::vec_out)
      .def("__len__", [](const Bmad::ConcatTransferMat &) { return 2; })
      .def("__getitem__", [](const Bmad::ConcatTransferMat &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.mat_out);
        if (i == 1)
          return nb::cast(s.vec_out);
        throw nb::index_error();
      });
  m.def(
      "concat_transfer_mat",
      &Bmad::concat_transfer_mat,
      nb::arg("mat_1"),
      nb::arg("vec_1"),
      nb::arg("mat_0"),
      nb::arg("vec_0"),
      R"""(Routine to concatinate two linear maps:
  mat_out = matmul(mat_1, mat_0)
  vec_out = matmul(mat_1, vec_0) + vec_1

Parameters
----------
mat_1 : 2D array of float (shape: 6,6)
    Map from s1 to s2

vec_1 : 1D array of float (shape: 6)
    Map from s1 to s2

mat_0 : 2D array of float (shape: 6,6)
    Map from s0 to s1

vec_0 : 1D array of float (shape: 6)
    Map from s0 to s1

Returns
-------
mat_out : 2D array of float (shape: 6,6)
    Map from s0 to s2

vec_out : 1D array of float (shape: 6)
    Map from s0 to s2
)"""
  );
  m.def(
      "control_bookkeeper",
      [](LatStruct &lat, EleStruct *ele, std::optional<bool> err_flag) {
        auto fn = static_cast<void (*)(LatStruct &, optional_ref<EleStruct>, std::optional<bool>)>(
            &Bmad::control_bookkeeper
        );
        return fn(lat, ptr_to_opt_ref(ele), err_flag);
      },
      nb::arg("lat"),
      nb::arg("ele") = nb::none(),
      nb::arg("err_flag") = nb::none(),
      R"""(Wrapper for Fortran routine control_bookkeeper

Parameters
----------
lat : LatStruct
    lattice to be used

ele : EleStruct, optional
    Element whose attribute values have been changed. If not present bookkeeping will be done for all
    elements.

err_flag : bool, optional
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "convert_bend_exact_multipole",
      &Bmad::convert_bend_exact_multipole,
      nb::arg("g"),
      nb::arg("out_type"),
      nb::arg("an"),
      nb::arg("bn"),
      R"""(Wrapper for Fortran routine convert_bend_exact_multipole

Parameters
----------
g : float
    1/rho bending strength.

out_type : int
    Output type: horizontally_pure$ or vertically_pure$.

an : 1D array of float (shape: 0:n_pole_maxx)
    Skew multipoles.
    This parameter is an input/output and is modified in-place.
    As an output, an: Converted skew multipoles.

bn : 1D array of float (shape: 0:n_pole_maxx)
    Non-skew multipoles.
    This parameter is an input/output and is modified in-place.
    As an output, bn: Converted Non-skew multipoles.
)"""
  );
  nb::class_<Bmad::ConvertCoords>(m, "ConvertCoords", "convert_coords return type")
      .def_ro("out_type_str", &Bmad::ConvertCoords::out_type_str)
      .def_ro("coord_out", &Bmad::ConvertCoords::coord_out)
      .def_ro("err_flag", &Bmad::ConvertCoords::err_flag)
      .def("__len__", [](const Bmad::ConvertCoords &) { return 3; })
      .def("__getitem__", [](const Bmad::ConvertCoords &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.out_type_str);
        if (i == 1)
          return nb::cast(s.coord_out);
        if (i == 2)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "convert_coords",
      &Bmad::convert_coords,
      nb::arg("in_type_str"),
      nb::arg("coord_in"),
      nb::arg("ele"),
      R"""(Wrapper for Fortran routine convert_coords

Parameters
----------
in_type_str : str
    type of the input coords.

coord_in : CoordStruct
    Input coordinates.

ele : EleStruct
    Provides the Twiss parameters.

Returns
-------
out_type_str : str
    type of the output coords.

coord_out : CoordStruct
    Output coordinates.

err_flag : bool, optional
    Set True if there is an error. False otherwise. in_type_str and out_type_str can be: 'LAB'
    {x, x', y, y', z, z'} 'MODE'               {a, a', b, b', z, z'} 'NORMALIZED'         {a_bar, a'_bar,
    b_bar, b'_bar, z_bar, z'_bar} 'ACTION-ANGLE'       {j_a, phi_a, j_b, phi_b, j_z,  phi_z} x_vec = V_mat *
    (a_vec + eta_vec * z') a_bar  =  sqrt(2*j_a) * cos(phi_a) a'_bar = -sqrt(2*j_a) * sin(phi_a)
)"""
  );
  m.def(
      "convert_field_ele_to_lab",
      &Bmad::convert_field_ele_to_lab,
      nb::arg("ele"),
      nb::arg("s_here"),
      nb::arg("forward_transform"),
      nb::arg("calc_dfield") = nb::none(),
      nb::arg("calc_potential") = nb::none(),
      R"""(Convert fields: ele to lab coords

Parameters
----------
ele : EleStruct
    Lattice element.

s_here : float
    real(rp) S-position.

forward_transform : bool
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
      nb::arg("x"),
      nb::arg("z"),
      nb::arg("g"),
      nb::arg("xout"),
      nb::arg("sout"),
      R"""(Wrapper for Fortran routine convert_local_cartesian_to_local_curvilinear

Parameters
----------
x : float

z : float

g : float

xout : float

sout : float
)"""
  );
  m.def(
      "convert_local_curvilinear_to_local_cartesian",
      &Bmad::convert_local_curvilinear_to_local_cartesian,
      nb::arg("x"),
      nb::arg("s"),
      nb::arg("g"),
      nb::arg("xout"),
      nb::arg("zout"),
      R"""(Wrapper for Fortran routine convert_local_curvilinear_to_local_cartesian

Parameters
----------
x : float

s : float

g : float

xout : float

zout : float
)"""
  );
  m.def(
      "convert_particle_coordinates_s_to_t",
      &Bmad::convert_particle_coordinates_s_to_t,
      nb::arg("particle"),
      nb::arg("s_body"),
      nb::arg("orientation"),
      R"""(Wrapper for Fortran routine convert_particle_coordinates_s_to_t

Parameters
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
      nb::arg("particle"),
      nb::arg("ele"),
      nb::arg("use_downstream_p0c") = nb::none(),
      R"""(Wrapper for Fortran routine convert_particle_coordinates_t_to_s

Parameters
----------
particle : CoordStruct
    Particle with .vec(:) in t-coords.

ele : EleStruct
    Element particle is going through.

use_downstream_p0c : bool, optional
    If True (the default), use ele.value(p0c$) as the reference momentum. If False, use ele.value(p0c_start$)
    as the reference.

Returns
-------
s_body : float, optional
    s-position in element body coords.
)"""
  );
  nb::class_<Bmad::ConvertPcTo>(m, "ConvertPcTo", "convert_pc_to return type")
      .def_ro("E_tot", &Bmad::ConvertPcTo::E_tot)
      .def_ro("gamma", &Bmad::ConvertPcTo::gamma)
      .def_ro("kinetic", &Bmad::ConvertPcTo::kinetic)
      .def_ro("beta", &Bmad::ConvertPcTo::beta)
      .def_ro("brho", &Bmad::ConvertPcTo::brho)
      .def_ro("beta1", &Bmad::ConvertPcTo::beta1)
      .def_ro("err_flag", &Bmad::ConvertPcTo::err_flag)
      .def("__len__", [](const Bmad::ConvertPcTo &) { return 7; })
      .def("__getitem__", [](const Bmad::ConvertPcTo &s, int i) -> nb::object {
        if (i < 0)
          i += 7;
        if (i == 0)
          return nb::cast(s.E_tot);
        if (i == 1)
          return nb::cast(s.gamma);
        if (i == 2)
          return nb::cast(s.kinetic);
        if (i == 3)
          return nb::cast(s.beta);
        if (i == 4)
          return nb::cast(s.brho);
        if (i == 5)
          return nb::cast(s.beta1);
        if (i == 6)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "convert_pc_to",
      &Bmad::convert_pc_to,
      nb::arg("pc"),
      nb::arg("particle"),
      R"""(Wrapper for Fortran routine convert_pc_to

Parameters
----------
pc : float
    Particle momentum

particle : int
    Type of particle. positron$, etc.

Returns
-------
E_tot : float, optional
    Total energy of the particle.

gamma : float, optional
    Gamma factor.

kinetic : float, optional
    Kinetic energy

beta : float, optional
    velocity / c_light

brho : float, optional
    Nominal B_field*rho_bend

beta1 : float, optional
    1 - beta. Equal to 1/(2*gamma^2) in ultra-rel limit.

err_flag : bool, optional
    Set true if there is an error. False otherwise.
)"""
  );
  nb::class_<Bmad::ConvertTotalEnergyTo>(
      m,
      "ConvertTotalEnergyTo",
      "convert_total_energy_to return type"
  )
      .def_ro("gamma", &Bmad::ConvertTotalEnergyTo::gamma)
      .def_ro("kinetic", &Bmad::ConvertTotalEnergyTo::kinetic)
      .def_ro("beta", &Bmad::ConvertTotalEnergyTo::beta)
      .def_ro("pc", &Bmad::ConvertTotalEnergyTo::pc)
      .def_ro("brho", &Bmad::ConvertTotalEnergyTo::brho)
      .def_ro("beta1", &Bmad::ConvertTotalEnergyTo::beta1)
      .def_ro("err_flag", &Bmad::ConvertTotalEnergyTo::err_flag)
      .def("__len__", [](const Bmad::ConvertTotalEnergyTo &) { return 7; })
      .def("__getitem__", [](const Bmad::ConvertTotalEnergyTo &s, int i) -> nb::object {
        if (i < 0)
          i += 7;
        if (i == 0)
          return nb::cast(s.gamma);
        if (i == 1)
          return nb::cast(s.kinetic);
        if (i == 2)
          return nb::cast(s.beta);
        if (i == 3)
          return nb::cast(s.pc);
        if (i == 4)
          return nb::cast(s.brho);
        if (i == 5)
          return nb::cast(s.beta1);
        if (i == 6)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "convert_total_energy_to",
      &Bmad::convert_total_energy_to,
      nb::arg("E_tot"),
      nb::arg("particle"),
      nb::arg("print_err") = nb::none(),
      R"""(Wrapper for Fortran routine convert_total_energy_to

Parameters
----------
E_tot : float
    Total energy of the particle.

particle : int
    Type of particle. positron$, etc.

print_err : bool, optional
    Print error message if E_tot < particle mass? Default is True.

Returns
-------
gamma : float, optional
    Gamma factor. Set to -1 for photons.

kinetic : float, optional
    Kinetic energy

beta : float, optional
    velocity / c_light

pc : float, optional
    Particle momentum

brho : float, optional
    Nominal B_field*rho_bend

beta1 : float, optional
    1 - beta. Equal to 1/(2*gamma^2) in ultra-rel limit.

err_flag : bool, optional
    Set true if there is an error. False otherwise.
)"""
  );
  nb::class_<Bmad::ConverterDistributionParser>(
      m,
      "ConverterDistributionParser",
      "converter_distribution_parser return type"
  )
      .def_ro("delim", &Bmad::ConverterDistributionParser::delim)
      .def_ro("delim_found", &Bmad::ConverterDistributionParser::delim_found)
      .def_ro("err_flag", &Bmad::ConverterDistributionParser::err_flag)
      .def("__len__", [](const Bmad::ConverterDistributionParser &) { return 3; })
      .def("__getitem__", [](const Bmad::ConverterDistributionParser &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.delim);
        if (i == 1)
          return nb::cast(s.delim_found);
        if (i == 2)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "converter_distribution_parser",
      &Bmad::converter_distribution_parser,
      nb::arg("ele"),
      R"""(Wrapper for Fortran routine converter_distribution_parser

Parameters
----------
ele : EleStruct
    Converter element.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Converter element with .converter field set.

Returns
-------
delim : str
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
      nb::arg("coord2"),
      R"""( Subroutine that is used to set one coord equal to another.

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
      nb::arg("coord_state"),
      nb::arg("one_word") = nb::none(),
      R"""(Routine to return the string representation of a coord%state state.

Parameters
----------
coord_state : int
    coord.state value

Returns
-------
state_str : str
    String representation.
)"""
  );
  nb::class_<Bmad::CoordsBodyToLocal>(m, "CoordsBodyToLocal", "coords_body_to_local return type")
      .def_ro("w_mat", &Bmad::CoordsBodyToLocal::w_mat)
      .def_ro("local_position", &Bmad::CoordsBodyToLocal::local_position)
      .def("__len__", [](const Bmad::CoordsBodyToLocal &) { return 2; })
      .def("__getitem__", [](const Bmad::CoordsBodyToLocal &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.w_mat);
        if (i == 1)
          return nb::cast(s.local_position);
        throw nb::index_error();
      });
  m.def(
      "coords_body_to_local",
      &Bmad::coords_body_to_local,
      nb::arg("body_position"),
      nb::arg("ele"),
      nb::arg("calculate_angles") = nb::none(),
      R"""(Wrapper for Fortran routine coords_body_to_local

Parameters
----------
body_position : FloorPositionStruct
    Element body frame coordinates. .r(3)               [x, y, s] position with s = Position from entrance end
    of element.

ele : EleStruct
    element that local_position coordinates are relative to.

calculate_angles : bool, optional
    calculate angles for local_position Default: True. False returns local_position angles (.theta, .phi,
    .psi) = 0.

Returns
-------
local_position : FloorPositionStruct
    Local laboratory coordinates. .r(3)               [x, y, s] position with s = Position from entrance end
    of element.

w_mat : 2D array of float (shape: 3,3), optional
    W matrix at to transform vectors. v_local  = w_mat . v_body v_body   = transpose(w_mat) . v_local
)"""
  );
  nb::class_<Bmad::CoordsBodyToRelExit>(
      m,
      "CoordsBodyToRelExit",
      "coords_body_to_rel_exit return type"
  )
      .def_ro("w_mat", &Bmad::CoordsBodyToRelExit::w_mat)
      .def_ro("rel_exit", &Bmad::CoordsBodyToRelExit::rel_exit)
      .def("__len__", [](const Bmad::CoordsBodyToRelExit &) { return 2; })
      .def("__getitem__", [](const Bmad::CoordsBodyToRelExit &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.w_mat);
        if (i == 1)
          return nb::cast(s.rel_exit);
        throw nb::index_error();
      });
  m.def(
      "coords_body_to_rel_exit",
      &Bmad::coords_body_to_rel_exit,
      nb::arg("body_position"),
      nb::arg("ele"),
      nb::arg("calculate_angles") = nb::none(),
      R"""(Wrapper for Fortran routine coords_body_to_rel_exit

Parameters
----------
body_position : FloorPositionStruct
    Element body frame coordinates. .r                  [x, y, s] position with s = Position from entrance end
    of element .

ele : EleStruct
    element that rel_exit coordinates are relative to.

calculate_angles : bool, optional
    calculate angles for rel_exit Default: True. False returns rel_exit angles (.theta, .phi, .psi) = 0.

Returns
-------
rel_exit : FloorPositionStruct
    Cartesian coordinates relative to exit of the element.

w_mat : 2D array of float (shape: 3,3), optional
    W matrix at to transform vectors. v_rel_exit = w_mat . v_body v_body     = transpose(w_mat) . v_rel_exit
)"""
  );
  nb::class_<Bmad::CoordsCurvilinearToFloor>(
      m,
      "CoordsCurvilinearToFloor",
      "coords_curvilinear_to_floor return type"
  )
      .def_ro("err_flag", &Bmad::CoordsCurvilinearToFloor::err_flag)
      .def_ro("global_", &Bmad::CoordsCurvilinearToFloor::global)
      .def("__len__", [](const Bmad::CoordsCurvilinearToFloor &) { return 2; })
      .def("__getitem__", [](const Bmad::CoordsCurvilinearToFloor &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.global);
        throw nb::index_error();
      });
  m.def(
      "coords_curvilinear_to_floor",
      &Bmad::coords_curvilinear_to_floor,
      nb::arg("xys"),
      nb::arg("branch"),
      R"""(Wrapper for Fortran routine coords_curvilinear_to_floor

Parameters
----------
xys : 1D array of float (shape: 3)
    (x, y, s) lab frame position vector.

branch : BranchStruct
    Lattice branch that defines the local reference coordinates.

Returns
-------
err_flag : bool
    Set True if global floor position cannot be computed.

global : FloorPositionStruct
    Global floor position corresponding to (x, y, s) --    .w    -- W matrix to transform vectors: v_global =
    w_mat * v_local
)"""
  );
  nb::class_<Bmad::CoordsFloorToCurvilinear>(
      m,
      "CoordsFloorToCurvilinear",
      "coords_floor_to_curvilinear return type"
  )
      .def_ro("ele1", &Bmad::CoordsFloorToCurvilinear::ele1)
      .def_ro("status", &Bmad::CoordsFloorToCurvilinear::status)
      .def_ro("w_mat", &Bmad::CoordsFloorToCurvilinear::w_mat)
      .def_ro("local_coords", &Bmad::CoordsFloorToCurvilinear::local_coords)
      .def("__len__", [](const Bmad::CoordsFloorToCurvilinear &) { return 4; })
      .def("__getitem__", [](const Bmad::CoordsFloorToCurvilinear &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.ele1);
        if (i == 1)
          return nb::cast(s.status);
        if (i == 2)
          return nb::cast(s.w_mat);
        if (i == 3)
          return nb::cast(s.local_coords);
        throw nb::index_error();
      });
  m.def(
      "coords_floor_to_curvilinear",
      &Bmad::coords_floor_to_curvilinear,
      nb::arg("floor_coords"),
      nb::arg("ele0"),
      R"""(Wrapper for Fortran routine coords_floor_to_curvilinear

Parameters
----------
floor_coords : FloorPositionStruct
    .r = [X, Y, Z] position in global coordinates

ele0 : EleStruct
    Element to start the search at.

Returns
-------
status : int
    ok$             -> Local_coords found. patch_problem$  -> No solution due to a patch element. outside$
    -> Outside of lattice ends (for open lattices).

local_coords : FloorPositionStruct
    .r = [x, y, s] position in curvilinear coordinates with respect to ele1 with s relative to start the
    lattice branch.

ele1 : EleStruct, optional
    Element that local_coords is with respect to.

w_mat : 2D array of float (shape: 3,3), optional
    W matrix at s, to transform vectors from floor to local. w_mat will only be well defined if status = ok$
)"""
  );
  nb::class_<Bmad::CoordsFloorToLocalCurvilinear>(
      m,
      "CoordsFloorToLocalCurvilinear",
      "coords_floor_to_local_curvilinear return type"
  )
      .def_ro("status", &Bmad::CoordsFloorToLocalCurvilinear::status)
      .def_ro("w_mat", &Bmad::CoordsFloorToLocalCurvilinear::w_mat)
      .def_ro("local_position", &Bmad::CoordsFloorToLocalCurvilinear::local_position)
      .def("__len__", [](const Bmad::CoordsFloorToLocalCurvilinear &) { return 3; })
      .def("__getitem__", [](const Bmad::CoordsFloorToLocalCurvilinear &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.status);
        if (i == 1)
          return nb::cast(s.w_mat);
        if (i == 2)
          return nb::cast(s.local_position);
        throw nb::index_error();
      });
  m.def(
      "coords_floor_to_local_curvilinear",
      &Bmad::coords_floor_to_local_curvilinear,
      nb::arg("global_position"),
      nb::arg("ele"),
      nb::arg("relative_to") = nb::none(),
      R"""(Wrapper for Fortran routine coords_floor_to_local_curvilinear

Parameters
----------
global_position : FloorPositionStruct
    .r = [X, Y, Z] position in global coordinates

ele : EleStruct
    element to find local coordinates of.

relative_to : int, optional
    not_set$ (default), upstream_end$, or downstream_end$. Force which end is used for z = 0. If
    upstream_end$, local_position.r(3) is relative to the upstream end which will not be the entrance end if
    ele.orientation = -1.

Returns
-------
status : int
    longitudinal position: inside$: Inside the element. upstream_end$: At upstream end of element or beyound.
    downstream_end$: At downstream end of element or beyound.

local_position : FloorPositionStruct
    .r = [x, y, z] position in local curvilinear coordinates.

w_mat : 2D array of float (shape: 3,3), optional
    W matrix at s, to transform vectors. v_global = w_mat.v_local v_local = transpose(w_mat).v_global
)"""
  );
  m.def(
      "coords_floor_to_relative",
      &Bmad::coords_floor_to_relative,
      nb::arg("floor0"),
      nb::arg("global_position"),
      nb::arg("calculate_angles") = nb::none(),
      nb::arg("is_delta_position") = nb::none(),
      R"""(Wrapper for Fortran routine coords_floor_to_relative

Parameters
----------
floor0 : FloorPositionStruct
    reference position

global_position : FloorPositionStruct
    global position

calculate_angles : bool, optional
    calculate angles for local_position Default: True. False returns local_position angles (.theta, .phi,
    .psi) = 0.

is_delta_position : bool, optional
    If True then treat global_position.r as a difference position in global space and only rotate the position
    but not shift it. Default: False.

Returns
-------
local_position : FloorPositionStruct
    position relative to floor0
)"""
  );
  nb::class_<Bmad::CoordsLocalCurvilinearToBody>(
      m,
      "CoordsLocalCurvilinearToBody",
      "coords_local_curvilinear_to_body return type"
  )
      .def_ro("w_mat", &Bmad::CoordsLocalCurvilinearToBody::w_mat)
      .def_ro("body_position", &Bmad::CoordsLocalCurvilinearToBody::body_position)
      .def("__len__", [](const Bmad::CoordsLocalCurvilinearToBody &) { return 2; })
      .def("__getitem__", [](const Bmad::CoordsLocalCurvilinearToBody &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.w_mat);
        if (i == 1)
          return nb::cast(s.body_position);
        throw nb::index_error();
      });
  m.def(
      "coords_local_curvilinear_to_body",
      &Bmad::coords_local_curvilinear_to_body,
      nb::arg("local_position"),
      nb::arg("ele"),
      nb::arg("calculate_angles") = nb::none(),
      R"""(Wrapper for Fortran routine coords_local_curvilinear_to_body

Parameters
----------
local_position : FloorPositionStruct
    local coordinates. .r(3)               [x, y, s] position with s = Position from entrance end of element.

ele : EleStruct
    element that coordinates are relative to.

calculate_angles : bool, optional
    calculate angles for body_position Default: True. False returns body_position angles (.theta, .phi, .psi)
    = 0.

Returns
-------
body_position : FloorPositionStruct
    Element coordinates relative to exit of the element. .r(3)               [x, y, s] position with s =
    Position from entrance end of element.

w_mat : 2D array of float (shape: 3,3), optional
    W matrix at to transform vectors. v_local  = w_mat . v_body v_body   = transpose(w_mat) . v_local
)"""
  );
  nb::class_<Bmad::CoordsLocalCurvilinearToFloor>(
      m,
      "CoordsLocalCurvilinearToFloor",
      "coords_local_curvilinear_to_floor return type"
  )
      .def_ro("w_mat", &Bmad::CoordsLocalCurvilinearToFloor::w_mat)
      .def_ro("global_position", &Bmad::CoordsLocalCurvilinearToFloor::global_position)
      .def("__len__", [](const Bmad::CoordsLocalCurvilinearToFloor &) { return 2; })
      .def("__getitem__", [](const Bmad::CoordsLocalCurvilinearToFloor &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.w_mat);
        if (i == 1)
          return nb::cast(s.global_position);
        throw nb::index_error();
      });
  m.def(
      "coords_local_curvilinear_to_floor",
      &Bmad::coords_local_curvilinear_to_floor,
      nb::arg("local_position"),
      nb::arg("ele"),
      nb::arg("in_body_frame") = nb::none(),
      nb::arg("calculate_angles") = nb::none(),
      nb::arg("end_origin") = nb::none(),
      nb::arg("downstream_dir_ref") = nb::none(),
      R"""(Wrapper for Fortran routine coords_local_curvilinear_to_floor

Parameters
----------
local_position : FloorPositionStruct
    Floor position in local curvilinear coordinates, with .r = [x, y, z_local] where z_local is wrt the
    entrance end of the element except if end_origin = downstream_end$. In this case, z_local is a distance
    -ele.value(l$) from the exit end (important for patch elements).

ele : EleStruct
    element that local_position coordinates are relative to.

in_body_frame : bool, optional
    True => local_position is in ele body frame and includes misalignments. Ignored if element is a patch.
    Default: False.

calculate_angles : bool, optional
    calculate angles for global_position Default: True. False returns local_position angles (.theta, .phi,
    .psi) = 0.

end_origin : int, optional
    not_set$ (default), upstream_end$, or downstream_end$. Force which end is used for z = 0. If
    upstream_end$, local_position.r(3) is relative to the upstream end which will not be the entrance end if
    ele.orientation = -1.

downstream_dir_ref : bool, optional
    Default False. The output theta angle is calculated so that moduo 2pi this angle is near ele.floor.theta.
    If the element is reversed (ele.direction = -1), the element body coords point upstream which is not
    always wanted. If this arg is set True, ele.floor.theta+pi modulo to be in the range [-pi, pi] is the
    reference.

Returns
-------
global_position : FloorPositionStruct
    Position in global coordinates.

w_mat : 2D array of float (shape: 3,3), optional
    W matrix at z, to transform vectors. v_global     = w_mat . v_local/body v_local/body = transpose(w_mat) .
    v_global
)"""
  );
  m.def(
      "coords_relative_to_floor",
      &Bmad::coords_relative_to_floor,
      nb::arg("floor0"),
      nb::arg("dr"),
      nb::arg("theta") = nb::none(),
      nb::arg("phi") = nb::none(),
      nb::arg("psi") = nb::none(),
      R"""(Wrapper for Fortran routine coords_relative_to_floor

Parameters
----------
floor0 : FloorPositionStruct
    Initial reference frame.

dr : 1D array of float (shape: 3)
    (x, y, z) positional shift of the reference frame.

theta : float, optional
    Angular shift of the reference frame. See the Bmad manual on the Global Coordinate system for more
    details. All angles must either be absent or present.

phi : float, optional
    Angular shift of the reference frame. See the Bmad manual on the Global Coordinate system for more
    details. All angles must either be absent or present.

psi : float, optional
    Angular shift of the reference frame. See the Bmad manual on the Global Coordinate system for more
    details. All angles must either be absent or present.

Returns
-------
floor1 : FloorPositionStruct
    Shifted reference frame.
)"""
  );
  m.def(
      "coulombfun",
      &Bmad::coulombfun,
      nb::arg("u"),
      nb::arg("v"),
      nb::arg("w"),
      nb::arg("gam"),
      R"""(Wrapper for Fortran routine coulombfun

Parameters
----------
u : float

v : float

w : float

gam : float

Returns
-------
res : float
)"""
  );
  m.def(
      "create_concatenated_wall3d",
      &Bmad::create_concatenated_wall3d,
      nb::arg("lat"),
      nb::arg("err"),
      R"""(Routine to concatinate lat%branch(i)ele(:)%wall3d%section(:) arrays into
one lat%branch(i)%wall3d%section(:) array.

Exceptions: capillary and aperture elements do not have their walls included.

Module needed:
  use wall3d_mod

Parameters
----------
lat : LatStruct
    lattice
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lattice
)"""
  );
  nb::class_<Bmad::CreateElementSlice>(m, "CreateElementSlice", "create_element_slice return type")
      .def_ro("sliced_ele", &Bmad::CreateElementSlice::sliced_ele)
      .def_ro("err_flag", &Bmad::CreateElementSlice::err_flag)
      .def("__len__", [](const Bmad::CreateElementSlice &) { return 2; })
      .def("__getitem__", [](const Bmad::CreateElementSlice &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.sliced_ele);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "create_element_slice",
      [](EleStruct &ele_in,
         double l_slice,
         double offset,
         LatParamStruct &param,
         bool include_upstream_end,
         bool include_downstream_end,
         EleStruct *old_slice,
         CoordStruct *orb_in) {
        auto fn = static_cast<
            Bmad::
                CreateElementSlice (*)(EleStruct &, double, double, LatParamStruct &, bool, bool, optional_ref<EleStruct>, optional_ref<CoordStruct>)>(
            &Bmad::create_element_slice
        );
        return fn(
            ele_in,
            l_slice,
            offset,
            param,
            include_upstream_end,
            include_downstream_end,
            ptr_to_opt_ref(old_slice),
            ptr_to_opt_ref(orb_in)
        );
      },
      nb::arg("ele_in"),
      nb::arg("l_slice"),
      nb::arg("offset"),
      nb::arg("param"),
      nb::arg("include_upstream_end"),
      nb::arg("include_downstream_end"),
      nb::arg("old_slice") = nb::none(),
      nb::arg("orb_in") = nb::none(),
      R"""(Wrapper for Fortran routine create_element_slice

Parameters
----------
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

old_slice : EleStruct, optional
    Previous slice or, if offset = 0, the previous element. If present this saves computation time of the
    reference energy and time at the start of the present slice. Also makes the ref energy continuous (there
    can be some small differences when using, say, runge_kutta tracking due to tracking tolerances).

orb_in : CoordStruct, optional
    Incoming orbit if calling routine is doing tracking through the slice. This is used when old_slice is not
    present and there may be an adjustment needed to the orbit ref energy (EG space charge tracking does not
    keep track of ref energy through an lcavity).

Returns
-------
sliced_ele : EleStruct
    Sliced_ele element with appropriate values set.

err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "create_feedback",
      &Bmad::create_feedback,
      nb::arg("lord"),
      nb::arg("input"),
      nb::arg("output"),
      nb::arg("err_flag"),
      R"""(Wrapper for Fortran routine create_feedback

Parameters
----------
lord : EleStruct
    Feedback element.
    This parameter is an input/output and is modified in-place.
    As an output, lord: Modified feedback elment.

input : 1D array of str
    Names of input slaves.

output : 1D array of str
    Names of output slaves.

err_flag : bool
    Set True if there is a problem.
)"""
  );
  m.def(
      "create_field_overlap",
      &Bmad::create_field_overlap,
      nb::arg("lat"),
      nb::arg("lord_name"),
      nb::arg("slave_name"),
      R"""(Wrapper for Fortran routine create_field_overlap

Parameters
----------
lat : LatStruct
    Lattice

lord_name : str
    Name of the element with a field extending beyound it's bounds.

slave_name : str
    Name of the element the lord's field overlaps.

Returns
-------
err_flag : bool
    Set true if there is a problem (like no elements found).
)"""
  );
  m.def(
      "create_girder",
      &Bmad::create_girder,
      nb::arg("lat"),
      nb::arg("ix_girder"),
      nb::arg("contrl"),
      nb::arg("girder_info"),
      nb::arg("err_flag"),
      R"""(Wrapper for Fortran routine create_girder

Parameters
----------
lat : LatStruct
    Lat to modify.
    This parameter is an input/output and is modified in-place.
    As an output, lat: Modified lattice.

ix_girder : int
    Index of girder element.

contrl : 1D array of ControlStruct
    Array of elements that are supported by the girder.

girder_info : EleStruct
    Element containing attributes to be transfered to the Girder element: girder_info.name girder_info.alias
    girder_info.descrip girder_info.value(:)

err_flag : bool
)"""
  );
  m.def(
      "create_group",
      &Bmad::create_group,
      nb::arg("lord"),
      nb::arg("contrl"),
      nb::arg("err"),
      R"""(Wrapper for Fortran routine create_group

Parameters
----------
lord : EleStruct
    Group element. .control.type
    This parameter is an input/output and is modified in-place.
    As an output, lord: Modified group elment

contrl : 1D array of ControlStruct
    control info. 1 element for each slave.

err : bool
    Set True if an attribute is not free to be controlled.
)"""
  );
  m.def(
      "create_lat_ele_nametable",
      &Bmad::create_lat_ele_nametable,
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine create_lat_ele_nametable

Parameters
----------
lat : LatStruct
    Lattice.

Returns
-------
nametable : NametableStruct
    Nametable of the elment names
)"""
  );
  m.def(
      "create_overlay",
      &Bmad::create_overlay,
      nb::arg("lord"),
      nb::arg("contrl"),
      nb::arg("err"),
      R"""(Wrapper for Fortran routine create_overlay

Parameters
----------
lord : EleStruct
    Overlay element. .control.type
    This parameter is an input/output and is modified in-place.
    As an output, lord: Modified overlay elment

contrl : 1D array of ControlStruct
    control info. 1 element for each slave.

err : bool
    Set True if an attribute is not free to be controlled.
)"""
  );
  nb::class_<Bmad::CreatePlanarWigglerModel>(
      m,
      "CreatePlanarWigglerModel",
      "create_planar_wiggler_model return type"
  )
      .def_ro("lat", &Bmad::CreatePlanarWigglerModel::lat)
      .def_ro("err_flag", &Bmad::CreatePlanarWigglerModel::err_flag)
      .def("__len__", [](const Bmad::CreatePlanarWigglerModel &) { return 2; })
      .def("__getitem__", [](const Bmad::CreatePlanarWigglerModel &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.lat);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "create_planar_wiggler_model",
      &Bmad::create_planar_wiggler_model,
      nb::arg("wiggler_in"),
      nb::arg("print_err") = nb::none(),
      R"""(Routine to create series of bend and drift elements to serve as a replacement
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
print_err : bool, optional
    If True (default) print an error message if there is an error.

Returns
-------
lat : LatStruct
    Lattice containing the wiggler model

err_flag : bool, optional
    Set True if there is an error.
)"""
  );
  m.def(
      "create_ramper",
      &Bmad::create_ramper,
      nb::arg("lord"),
      nb::arg("contrl"),
      nb::arg("err"),
      R"""(Wrapper for Fortran routine create_ramper

Parameters
----------
lord : EleStruct
    Ramper element. .control.type
    This parameter is an input/output and is modified in-place.
    As an output, lord: Modified ramper elment

contrl : 1D array of ControlStruct
    control info. 1 element for each slave.

err : bool
    Set True if an attribute is not free to be controlled.
)"""
  );
  m.def(
      "create_sol_quad_model",
      &Bmad::create_sol_quad_model,
      nb::arg("sol_quad"),
      nb::arg("lat"),
      R"""(Routine to create series of solenoid and quadrupole elements to serve as a replacement
model for a sol_quad element.

This routine is helpful for translating bmad lattices to a language that does not
implement a combination solenoid/quadrupole.

Not yet implemented!
)"""
  );
  m.def(
      "create_unique_ele_names",
      &Bmad::create_unique_ele_names,
      nb::arg("lat"),
      nb::arg("key"),
      nb::arg("suffix"),
      R"""(Wrapper for Fortran routine create_unique_ele_names

Parameters
----------
lat : LatStruct
    Lattice holding the elements.
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lattice with names made unique.

key : int
    Class key of elements to consider.

suffix : str
    Suffix string. Must have a single "?" character.
)"""
  );
  m.def(
      "create_wiggler_cartesian_map",
      &Bmad::create_wiggler_cartesian_map,
      nb::arg("ele"),
      R"""(Wrapper for Fortran routine create_wiggler_cartesian_map

Parameters
----------
ele : EleStruct
    Wiggler or undulator element.

Returns
-------
cart_map : CartesianMapStruct
    Cartesian map.
)"""
  );
  m.def(
      "crystal_attribute_bookkeeper",
      &Bmad::crystal_attribute_bookkeeper,
      nb::arg("ele"),
      R"""(Wrapper for Fortran routine crystal_attribute_bookkeeper

Parameters
----------
ele : EleStruct
    Crystal element.
)"""
  );
  m.def(
      "crystal_h_misalign",
      &Bmad::crystal_h_misalign,
      nb::arg("ele"),
      nb::arg("orbit"),
      nb::arg("h_vec"),
      R"""(Routine reorient the crystal H vector due to local imperfections in the crystal lattice.

Parameters
----------
ele : EleStruct
    Crystal element

orbit : CoordStruct
    Photon position at crystal surface.

h_vec : 1D array of float (shape: 3)
    H vector before misalignment.
    This parameter is an input/output and is modified in-place.
    As an output, h_vec: H vector after misalignment.
)"""
  );
  m.def(
      "crystal_type_to_crystal_params",
      &Bmad::crystal_type_to_crystal_params,
      nb::arg("ele"),
      R"""(Routine to set the crystal parameters based upon the crystal type.

Crystal types are of the form:
  "ZZZ(ijk)"
Where "ZZZ" is the atomic formula of the crystal material and "ijk" is the reciprical lattice
vetor specifying the diffraction plans.

Parameters
----------
ele : EleStruct
    Crystal element.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Crystal element with computed parameter..

Returns
-------
err_flag : bool
    Set True if crystal type is unrecognized. False otherwise.
)"""
  );
  m.def(
      "custom_attribute_ubound_index",
      &Bmad::custom_attribute_ubound_index,
      nb::arg("ele_class"),
      R"""(Routine to return, for a given element class, the upper bound index for the ele%custom(:)
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
  nb::class_<Bmad::CustomEleAttribNameList>(
      m,
      "CustomEleAttribNameList",
      "custom_ele_attrib_name_list return type"
  )
      .def_ro("index_list", &Bmad::CustomEleAttribNameList::index_list)
      .def_ro("name_list", &Bmad::CustomEleAttribNameList::name_list)
      .def("__len__", [](const Bmad::CustomEleAttribNameList &) { return 2; })
      .def("__getitem__", [](const Bmad::CustomEleAttribNameList &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.index_list);
        if (i == 1)
          return nb::cast(s.name_list);
        throw nb::index_error();
      });
  m.def(
      "custom_ele_attrib_name_list",
      &Bmad::custom_ele_attrib_name_list,
      R"""(Routine to create an array (index_list(i), name_list(i)) of custom element attribute names and indexes.
Each name in the name_list is of the form:
  "{<class>::}<attribute_name>"
where:
  <class>:: is an optional class prefix.
  <attribute_name> is the name of the attribute.

Returns
-------
index_list : 1D array of int
    Index of the custom attribute.

name_list : 1D array of str
    List of custom attributes.
)"""
  );
}
