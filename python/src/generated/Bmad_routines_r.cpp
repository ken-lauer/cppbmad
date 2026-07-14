#include "pybmad/generated/Bmad_routines_r.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

PyRadiationIntegrals python_radiation_integrals(
    LatStruct &lat,
    CoordStructArray1D orb,
    std::optional<int> ix_cache = std::nullopt,
    std::optional<int> ix_branch = std::nullopt
) {
  auto _result = Bmad::radiation_integrals(lat, orb, make_opt_ref(ix_cache), ix_branch);
  auto py_result{PyRadiationIntegrals{_result, ix_cache}};
  return py_result;
}
PyReleaseRadIntCache python_release_rad_int_cache(int ix_cache) {
  Bmad::release_rad_int_cache(ix_cache);
  auto py_result{PyReleaseRadIntCache{ix_cache}};
  return py_result;
}

void init_Bmad_routines_r(nb::module_ &m) {
  nb::class_<Bmad::Rad1DampAndStocMats>(
      m,
      "Rad1DampAndStocMats",
      "rad1_damp_and_stoc_mats return type"
  )
      .def_ro("rad_map", &Bmad::Rad1DampAndStocMats::rad_map)
      .def_ro("err_flag", &Bmad::Rad1DampAndStocMats::err_flag)
      .def_ro("rad_int1", &Bmad::Rad1DampAndStocMats::rad_int1)
      .def("__len__", [](const Bmad::Rad1DampAndStocMats &) { return 3; })
      .def("__getitem__", [](const Bmad::Rad1DampAndStocMats &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.rad_map);
        if (i == 1)
          return nb::cast(s.err_flag);
        if (i == 2)
          return nb::cast(s.rad_int1);
        throw nb::index_error();
      });
  m.def(
      "rad1_damp_and_stoc_mats",
      [](EleStruct &ele,
         bool include_opening_angle,
         CoordStruct &orb_in,
         CoordStruct &orb_out,
         double g2_tol,
         double g3_tol,
         EleStruct *ele0) {
        auto fn = static_cast<Bmad::Rad1DampAndStocMats (*)(
            EleStruct &,
            bool,
            CoordStruct &,
            CoordStruct &,
            double,
            double,
            optional_ref<EleStruct>
        )>(&Bmad::rad1_damp_and_stoc_mats);
        return fn(
            ele,
            include_opening_angle,
            orb_in,
            orb_out,
            g2_tol,
            g3_tol,
            ptr_to_opt_ref(ele0)
        );
      },
      nb::arg("ele"),
      nb::arg("include_opening_angle"),
      nb::arg("orb_in"),
      nb::arg("orb_out"),
      nb::arg("g2_tol"),
      nb::arg("g3_tol"),
      nb::arg("ele0") = nb::none(),
      R"""(Routine to calculate the damping and stochastic matrices for a given lattice element.

Parameters
----------
ele : EleStruct
    Element under consideration.

include_opening_angle : bool
    If True include the effect of the vertical opening angle of emitted radiation. Generally use True unless
    comparing against other codes.

orb_in : CoordStruct
    Entrance orbit about which to compute the matrices.

orb_out : CoordStruct
    Exit orbit.

g2_tol : float
    Tollerance on g^2 per unit length (damping tolerance).

g3_tol : float
    Tollerance on g^3 per unit length (stocastic tolerance).

ele0 : EleStruct, optional
    Element before `ele`. Needed if and only if rad_int1 is present

Returns
-------
rad_map : RadMapStruct
    Damping and stochastic matrices.

err_flag : bool
    Set true if there is an error. False otherwise.

rad_int1 : RadInt1Struct, optional
    Radiation integrals
)"""
  );
  nb::class_<Bmad::RadDampAndStocMats>(
      m,
      "RadDampAndStocMats",
      "rad_damp_and_stoc_mats return type"
  )
      .def_ro("rmap", &Bmad::RadDampAndStocMats::rmap)
      .def_ro("mode", &Bmad::RadDampAndStocMats::mode)
      .def_ro("xfer_nodamp_mat", &Bmad::RadDampAndStocMats::xfer_nodamp_mat)
      .def_ro("err_flag", &Bmad::RadDampAndStocMats::err_flag)
      .def_ro("rad_int_branch", &Bmad::RadDampAndStocMats::rad_int_branch)
      .def("__len__", [](const Bmad::RadDampAndStocMats &) { return 5; })
      .def("__getitem__", [](const Bmad::RadDampAndStocMats &s, int i) -> nb::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return nb::cast(s.rmap);
        if (i == 1)
          return nb::cast(s.mode);
        if (i == 2)
          return nb::cast(s.xfer_nodamp_mat);
        if (i == 3)
          return nb::cast(s.err_flag);
        if (i == 4)
          return nb::cast(s.rad_int_branch);
        throw nb::index_error();
      });
  m.def(
      "rad_damp_and_stoc_mats",
      &Bmad::rad_damp_and_stoc_mats,
      nb::arg("ele1"),
      nb::arg("ele2"),
      nb::arg("include_opening_angle"),
      nb::arg("closed_orbit") = nb::none(),
      R"""(Routine to calculate the damping and stochastic variance matrices from exit end of ele1
to the exit end of ele2. Use ele1 = ele2 to get 1-turn matrices.

If ele2 is before ele1 the integration range if from ele1 to the branch end plus
from the beginning to ele2.

Note: The ele%mat6 matrices will be remade. By convention, these matrices
do not include damping.

Parameters
----------
ele1 : EleStruct
    Start element of integration range.

ele2 : EleStruct
    End element of integration range.

include_opening_angle : bool
    If True include the effect of the vertical opening angle of emitted radiation. Generally use True unless
    comparing against other codes.

closed_orbit : 1D array of CoordStruct, optional
    Closed orbit. If not present this routine will calculate it.

Returns
-------
rmap : RadMapStruct
    Damping and stochastic mats

mode : NormalModesStruct

xfer_nodamp_mat : 2D array of float (shape: 6,6)
    Transfer matrix without damping.

err_flag : bool
    Set true if there is a problem.

rad_int_branch : RadIntBranchStruct, optional
    Array of element-by-element radiation integrals.
)"""
  );
  nb::class_<Bmad::RadGIntegrals>(m, "RadGIntegrals", "rad_g_integrals return type")
      .def_ro("int_g", &Bmad::RadGIntegrals::int_g)
      .def_ro("int_g3", &Bmad::RadGIntegrals::int_g3)
      .def("__len__", [](const Bmad::RadGIntegrals &) { return 2; })
      .def("__getitem__", [](const Bmad::RadGIntegrals &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.int_g);
        if (i == 1)
          return nb::cast(s.int_g3);
        throw nb::index_error();
      });
  m.def(
      "rad_g_integrals",
      &Bmad::rad_g_integrals,
      nb::arg("ele"),
      nb::arg("where"),
      nb::arg("orb_in"),
      nb::arg("orb_out"),
      nb::arg("int_g2"),
      nb::arg("g_tol"),
      nb::arg("g2_tol"),
      nb::arg("g3_tol"),
      R"""(Routine to calculate bending strength integrals (g(s) = 1/trajectory_bending_radius(s)) in
laboratory coords.

Parameters
----------
ele : EleStruct
    Element under consideration.

where : int
    What part of ele to integrate over. upstream$ -> 1st half of element, downsteam$ -> 2nd half, all$ ->
    everything.

orb_in : CoordStruct
    Entrance orbit about which to compute the matrices.

orb_out : CoordStruct
    Exit orbit.

g_tol : float
    Tollerance on |g| per unit length.

g2_tol : float
    Tollerance on g^2 per unit length.

g3_tol : float
    Tollerance on g^3 per unit length.

Returns
-------
int_g : 1D array of float (shape: 2)
    Integrals of (gx,gy) vector.

int_g3 : float
    integrals of |g|^2 and |g|^3.
)"""
  );
  nb::class_<PyRadiationIntegrals>(m, "RadiationIntegrals", "radiation_integrals return type")
      .def_ro("mode", &PyRadiationIntegrals::mode)
      .def_ro("rad_int_by_ele", &PyRadiationIntegrals::rad_int_by_ele)
      .def_ro("ix_cache", &PyRadiationIntegrals::ix_cache)
      .def("__len__", [](const PyRadiationIntegrals &) { return 3; })
      .def("__getitem__", [](const PyRadiationIntegrals &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.mode);
        if (i == 1)
          return nb::cast(s.rad_int_by_ele);
        if (i == 2)
          return nb::cast(s.ix_cache);
        throw nb::index_error();
      });
  m.def(
      "radiation_integrals",
      &python_radiation_integrals,
      nb::arg("lat"),
      nb::arg("orb"),
      nb::arg("ix_cache") = nb::none(),
      nb::arg("ix_branch") = nb::none(),
      R"""(Wrapper for Fortran routine radiation_integrals

Parameters
----------
lat : LatStruct
    Lattice to use. The calculation assumes that the Twiss parameters have been calculated.

orb : 1D array of CoordStruct

ix_cache : int, optional
    Cache pointer.
    This parameter is an input/output and is modified in-place.
    As an output, ix_cache: Cache pointer. If ix_cache = 0 at input then

ix_branch : int, optional
    Lattice branch index. Default is 0.

Returns
-------
mode : NormalModesStruct
    Parameters for the ("horizontal like") a-mode, ("vertical like") b-mode, and the z-mode

ix_cache : int, optional
    Cache pointer.
    This parameter is an input/output and is modified in-place.
    As an output, ix_cache: Cache pointer. If ix_cache = 0 at input then

rad_int_by_ele : RadIntAllEleStruct, optional
    Radiation integrals element by element.
)"""
  );
  m.def(
      "radiation_map_setup",
      [](EleStruct &ele, CoordStruct *ref_orbit_in) {
        auto fn = static_cast<bool (*)(EleStruct &, optional_ref<CoordStruct>)>(
            &Bmad::radiation_map_setup
        );
        return fn(ele, ptr_to_opt_ref(ref_orbit_in));
      },
      nb::arg("ele"),
      nb::arg("ref_orbit_in") = nb::none(),
      R"""(Routine to calculate the radiation kick for a lattice element.

Parameters
----------
ele : EleStruct
    Element whose map is to be setup.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with map calculated.

Returns
-------
err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "ramper_slave_setup",
      &Bmad::ramper_slave_setup,
      nb::arg("lat"),
      nb::arg("do_setup") = nb::none(),
      R"""(Wrapper for Fortran routine ramper_slave_setup

Parameters
----------
lat : LatStruct
    Lattice to be setup.
    This parameter is an input/output and is modified in-place.
    As an output, lat: Lattice with ramper slaves setup.

do_setup : bool, optional
)"""
  );
  nb::class_<Bmad::RamperValue>(m, "RamperValue", "ramper_value return type")
      .def_ro("err_flag", &Bmad::RamperValue::err_flag)
      .def_ro("value", &Bmad::RamperValue::value)
      .def("__len__", [](const Bmad::RamperValue &) { return 2; })
      .def("__getitem__", [](const Bmad::RamperValue &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.err_flag);
        if (i == 1)
          return nb::cast(s.value);
        throw nb::index_error();
      });
  m.def(
      "ramper_value",
      &Bmad::ramper_value,
      nb::arg("ramper"),
      nb::arg("r1"),
      R"""(Wrapper for Fortran routine ramper_value

Parameters
----------
ramper : EleStruct
    Ramper lord.

r1 : ControlRamp1Struct
    Slave function.

Returns
-------
err_flag : bool
    Set True if there is an error, False otherwise.

value : float
    Value of the slave function.
)"""
  );
  m.def(
      "randomize_lr_wake_frequencies",
      &Bmad::randomize_lr_wake_frequencies,
      nb::arg("ele"),
      R"""(Routine to randomize the frequencies of the lr wake HOMs according to:
  freq = freq_in * (1 + lr_freq_spread) * rr)
where rr is a Gaussian distributed random number with unit variance.

Parameters
----------
ele : EleStruct
    Element with wake. If no wake then nothing is done.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Element with wake frequencies set.

Returns
-------
set_done : bool, optional
    Set True if there where lr wakes to be set. False otherwise.
)"""
  );
  m.def(
      "rchomp",
      &Bmad::rchomp,
      nb::arg("rel"),
      nb::arg("plc"),
      R"""(Wrapper for Fortran routine rchomp

Parameters
----------
rel : float

plc : int

Returns
-------
out : str
)"""
  );
  m.def(
      "re_allocate_eles",
      &Bmad::re_allocate_eles,
      nb::arg("eles"),
      nb::arg("n"),
      nb::arg("save_old") = nb::none(),
      nb::arg("exact") = nb::none(),
      R"""(Wrapper for Fortran routine re_allocate_eles

Parameters
----------
eles : 1D array of ElePointerStruct
    Array of element pointers with possible old data.
    This parameter is an input/output and is modified in-place.
    As an output, eles: Array of element pointers.

n : int
    Array size to set.

save_old : bool, optional
    If present and True then save the old data.

exact : bool, optional
    If present and True then eles will have size = n If False (default), reallcation will not be done if eles
    is already large enough
)"""
  );
  m.def(
      "re_allocate",
      nb::overload_cast<Wall3dSectionStructAlloc1D, int, std::optional<bool>>(&Bmad::re_allocate),
      nb::arg("section"),
      nb::arg("n"),
      nb::arg("exact") = nb::none(),
      R"""(Wrapper for Fortran routine re_allocate_wall3d_section_array

Parameters
----------
section : 1D array of Wall3dSectionStruct

n : int

exact : bool, optional
)"""
  );
  m.def(
      "re_allocate",
      nb::overload_cast<Wall3dVertexStructAlloc1D, int, std::optional<bool>>(&Bmad::re_allocate),
      nb::arg("v"),
      nb::arg("n"),
      nb::arg("exact") = nb::none(),
      R"""(Wrapper for Fortran routine re_allocate_wall3d_vertex_array

Parameters
----------
v : 1D array of Wall3dVertexStruct

n : int

exact : bool, optional
)"""
  );
  m.def(
      "re_associate_node_array",
      &Bmad::re_associate_node_array,
      nb::arg("tree"),
      nb::arg("n"),
      nb::arg("exact") = nb::none(),
      R"""(Routine to resize the tree%node(:) array.

Note: The data of the array is preserved but data at the end of the
array will be lost if n is less than the original size of the array

Parameters
----------
tree : ExpressionTreeStruct

n : int
    Size wanted.

exact : bool, optional
    Default is False. If False, the size of the output array is permitted to be larger than n.
)"""
  );
  m.def(
      "re_str",
      nb::overload_cast<long double>(&Bmad::re_str),
      nb::arg("rel"),
      R"""(Wrapper for Fortran routine re_str_qp

Parameters
----------
rel : float

Returns
-------
str_out : str
)"""
  );
  m.def(
      "re_str",
      nb::overload_cast<double>(&Bmad::re_str),
      nb::arg("rel"),
      R"""(Wrapper for Fortran routine re_str_rp

Parameters
----------
rel : float

Returns
-------
str_out : str
)"""
  );
  nb::class_<Bmad::ReadBeamAscii>(m, "ReadBeamAscii", "read_beam_ascii return type")
      .def_ro("beam", &Bmad::ReadBeamAscii::beam)
      .def_ro("err_flag", &Bmad::ReadBeamAscii::err_flag)
      .def("__len__", [](const Bmad::ReadBeamAscii &) { return 2; })
      .def("__getitem__", [](const Bmad::ReadBeamAscii &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.beam);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "read_beam_ascii",
      &Bmad::read_beam_ascii,
      nb::arg("file_name"),
      nb::arg("beam_init"),
      R"""(Subroutine to read in a beam definition file.
If non_zero, the following components of beam_init are used to rescale the beam:
    %n_bunch
    %n_particle
    %charge_tot

If the beam file has '.h5' or '.hdf5' suffix then the file is taken to be an HDF5 file.
Otherwise the file is assumed to be ASCII.

Parameters
----------
file_name : str
    Name of beam file.

beam_init : BeamInitStruct
    See above.

Returns
-------
beam : BeamStruct
    Structure holding the beam information.

err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  nb::class_<Bmad::ReadBeamFile>(m, "ReadBeamFile", "read_beam_file return type")
      .def_ro("beam", &Bmad::ReadBeamFile::beam)
      .def_ro("err_flag", &Bmad::ReadBeamFile::err_flag)
      .def("__len__", [](const Bmad::ReadBeamFile &) { return 2; })
      .def("__getitem__", [](const Bmad::ReadBeamFile &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.beam);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "read_beam_file",
      [](std::string file_name,
         BeamInitStruct &beam_init,
         EleStruct *ele,
         std::optional<bool> print_mom_shift_warning,
         std::optional<bool> conserve_momentum) {
        auto fn = static_cast<Bmad::ReadBeamFile (*)(
            std::string,
            BeamInitStruct &,
            optional_ref<EleStruct>,
            std::optional<bool>,
            std::optional<bool>
        )>(&Bmad::read_beam_file);
        return fn(
            file_name,
            beam_init,
            ptr_to_opt_ref(ele),
            print_mom_shift_warning,
            conserve_momentum
        );
      },
      nb::arg("file_name"),
      nb::arg("beam_init"),
      nb::arg("ele") = nb::none(),
      nb::arg("print_mom_shift_warning") = nb::none(),
      nb::arg("conserve_momentum") = nb::none(),
      R"""(Subroutine to read in a beam definition file.
If non_zero, the following components of beam_init are used to rescale the beam:
    %n_bunch
    %n_particle
    %bunch_charge -> charge_tot
    %species

If the beam file has '.h5' or '.hdf5' suffix then the file is taken to be an HDF5 file.
Otherwise the file is assumed to be ASCII.

Parameters
----------
file_name : str
    Name of beam file.

beam_init : BeamInitStruct
    See above.

ele : EleStruct, optional
    Element with reference energy, etc.

print_mom_shift_warning : bool, optional
    Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.

Returns
-------
beam : BeamStruct
    Structure holding the beam information.

err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "read_binary_cartesian_map",
      &Bmad::read_binary_cartesian_map,
      nb::arg("file_name"),
      nb::arg("ele"),
      nb::arg("cart_map"),
      nb::arg("err_flag"),
      R"""(Routine to read a binary cartesian_map structure.

Parameters
----------
file_name : str
    File to create.

ele : EleStruct
    Element associated with the map. Ouput:

cart_map : CartesianMapStruct
    cartesian_map_struct, cartesian map.

err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "read_binary_cylindrical_map",
      &Bmad::read_binary_cylindrical_map,
      nb::arg("file_name"),
      nb::arg("ele"),
      nb::arg("cl_map"),
      nb::arg("err_flag"),
      R"""(Routine to read a binary cylindrical_map structure.

Parameters
----------
file_name : str
    File to create.

ele : EleStruct
    Element associated with the map. Ouput:

cl_map : CylindricalMapStruct
    cylindrical_map_struct, cylindrical map.

err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "read_binary_grid_field",
      &Bmad::read_binary_grid_field,
      nb::arg("file_name"),
      nb::arg("ele"),
      nb::arg("g_field"),
      nb::arg("err_flag"),
      R"""(Routine to read a binary grid_field structure.

Parameters
----------
file_name : str
    File to create.

ele : EleStruct
    Element associated with the map. Ouput:

g_field : GridFieldStruct
    grid_field_struct, cylindrical map.

err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  nb::class_<Bmad::ReadDigestedBmadFile>(
      m,
      "ReadDigestedBmadFile",
      "read_digested_bmad_file return type"
  )
      .def_ro("lat", &Bmad::ReadDigestedBmadFile::lat)
      .def_ro("err_flag", &Bmad::ReadDigestedBmadFile::err_flag)
      .def_ro("parser_calling", &Bmad::ReadDigestedBmadFile::parser_calling)
      .def_ro("lat_files", &Bmad::ReadDigestedBmadFile::lat_files)
      .def("__len__", [](const Bmad::ReadDigestedBmadFile &) { return 4; })
      .def("__getitem__", [](const Bmad::ReadDigestedBmadFile &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.lat);
        if (i == 1)
          return nb::cast(s.err_flag);
        if (i == 2)
          return nb::cast(s.parser_calling);
        if (i == 3)
          return nb::cast(s.lat_files);
        throw nb::index_error();
      });
  m.def(
      "read_digested_bmad_file",
      &Bmad::read_digested_bmad_file,
      nb::arg("in_file_name"),
      nb::arg("version"),
      R"""(Wrapper for Fortran routine read_digested_bmad_file

Parameters
----------
in_file_name : str

version : int

Returns
-------
lat : LatStruct
    Output lattice structure

err_flag : bool, optional
    Set True if there is an error. False otherwise.

parser_calling : bool, optional
    Is this routine being called from a parser routine (like bmad_parser)? Default is False. This argument
    determines what are considered errors. For example, a moved digested file is considered an error if this
    routine is called from a parser but not otherwise. The reason for this dichotomy is that a parser is able
    to reread the original lattice file.

lat_files : 1D array of str, optional
    List of Bmad lattice files that defined this lattice.
)"""
  );
  m.def(
      "read_surface_reflection_file",
      &Bmad::read_surface_reflection_file,
      nb::arg("file_name"),
      R"""(Routine to read the reflection probability data for a given type of surface from a file.

Parameters
----------
file_name : str
    Name of the file.

Returns
-------
surface : PhotonReflectSurfaceStruct
    Surface info.
)"""
  );
  m.def(
      "reallocate_beam",
      &Bmad::reallocate_beam,
      nb::arg("beam"),
      nb::arg("n_bunch"),
      nb::arg("n_particle") = nb::none(),
      nb::arg("extend") = nb::none(),
      R"""(Wrapper for Fortran routine reallocate_beam

Parameters
----------
beam : BeamStruct
    Beam bunches are saved if save = True.
    This parameter is an input/output and is modified in-place.
    As an output, beam: Allocated beam_struct structure.

n_bunch : int
    Number of bunches.

n_particle : int, optional
    Number of particles. Must be non-negative. If save = True then the number of particles in existing bunches
    will not be touched. If not present, beam.bunch(i).particle(:) will be in an undefined state.

extend : bool, optional
)"""
  );
  m.def(
      "reallocate_bp_com_const",
      &Bmad::reallocate_bp_com_const,
      R"""(Wrapper for Fortran routine reallocate_bp_com_const
)"""
  );
  m.def(
      "reallocate_bunch",
      &Bmad::reallocate_bunch,
      nb::arg("n_particle"),
      nb::arg("save") = nb::none(),
      R"""(Wrapper for Fortran routine reallocate_bunch

Parameters
----------
n_particle : int
    Number of particles. Must be non-negative.

save : bool, optional
    If present and True then save the old bunch info.

Returns
-------
bunch : BunchStruct
    Allocated bunch_struct structure.
)"""
  );
  m.def(
      "reallocate_control",
      &Bmad::reallocate_control,
      nb::arg("lat"),
      nb::arg("n"),
      R"""(Wrapper for Fortran routine reallocate_control

Parameters
----------
lat : LatStruct
    Lattice.

n : int
    Array size for lat.control(:) and lat.ic(:).
)"""
  );
  m.def(
      "reallocate_coord",
      nb::overload_cast<CoordArrayStructAlloc1D, LatStruct &>(&Bmad::reallocate_coord),
      nb::arg("coord_array"),
      nb::arg("lat"),
      R"""(Routine to allocate or reallocate at allocatable coord_struct array.
reallocate_coord is an overloaded name for:
  reallocate_coord_n (coord, n_coord)
  reallocate_coord_lat (coord, lat, ix_branch)

Subroutine to allocate an allocatable coord_struct array to at least:
    coord(0:n_coord)                            if n_coord arg is used.
    coord(0:lat%branch(ix_branch)%n_ele_max)    if lat arg is used.

The old coordinates are saved
If, at input, coord(:) is not allocated, coord(0)%vec is set to zero.
In any case, coord(n)%vec for n > 0 is set to zero.

Parameters
----------
lat : LatStruct
    Lattice
)"""
  );
  m.def(
      "reallocate_coord",
      nb::overload_cast<CoordStructAlloc1D, LatStruct &, std::optional<int>>(
          &Bmad::reallocate_coord
      ),
      nb::arg("coord"),
      nb::arg("lat"),
      nb::arg("ix_branch") = nb::none(),
      R"""(Routine to allocate or reallocate at allocatable coord_struct array.
reallocate_coord is an overloaded name for:
  reallocate_coord_n (coord, n_coord)
  reallocate_coord_lat (coord, lat, ix_branch)

Subroutine to allocate an allocatable coord_struct array to at least:
    coord(0:n_coord)                            if n_coord arg is used.
    coord(0:lat%branch(ix_branch)%n_ele_max)    if lat arg is used.

The old coordinates are saved
If, at input, coord(:) is not allocated, coord(0)%vec is set to zero.
In any case, coord(n)%vec for n > 0 is set to zero.

Parameters
----------
coord : 1D array of CoordStruct
    Allocatable array.
    This parameter is an input/output and is modified in-place.
    As an output, coord: Allocated array.

lat : LatStruct
    Lattice

ix_branch : int, optional
    Branch to use. Default is 0 (main branch).
)"""
  );
  m.def(
      "reallocate_coord",
      nb::overload_cast<CoordStructAlloc1D, int>(&Bmad::reallocate_coord),
      nb::arg("coord"),
      nb::arg("n_coord"),
      R"""(Routine to allocate or reallocate at allocatable coord_struct array.
reallocate_coord is an overloaded name for:
  reallocate_coord_n (coord, n_coord)
  reallocate_coord_lat (coord, lat, ix_branch)

Subroutine to allocate an allocatable coord_struct array to at least:
    coord(0:n_coord)                            if n_coord arg is used.
    coord(0:lat%branch(ix_branch)%n_ele_max)    if lat arg is used.

The old coordinates are saved
If, at input, coord(:) is not allocated, coord(0)%vec is set to zero.
In any case, coord(n)%vec for n > 0 is set to zero.

Parameters
----------
coord : 1D array of CoordStruct
    Allocatable array.
    This parameter is an input/output and is modified in-place.
    As an output, coord: Allocated array.

n_coord : int
    Minimum array upper bound wanted.
)"""
  );
  m.def(
      "reallocate_expression_stack",
      &Bmad::reallocate_expression_stack,
      nb::arg("stack"),
      nb::arg("n"),
      nb::arg("exact") = nb::none(),
      R"""(Wrapper for Fortran routine reallocate_expression_stack

Parameters
----------
stack : 1D array of ExpressionAtomStruct
    Existing stack array.
    This parameter is an input/output and is modified in-place.
    As an output, stack: Resized stack.

n : int
    Array size needed.

exact : bool, optional
    If present and False then the size of the output array is permitted to be larger than n. Default is True.
)"""
  );
  m.def(
      "reallocate_sequence",
      &Bmad::reallocate_sequence,
      nb::arg("sequence"),
      nb::arg("n_seq"),
      R"""(No docstring available.
)"""
  );
  m.def(
      "rel_tracking_charge_to_mass",
      &Bmad::rel_tracking_charge_to_mass,
      nb::arg("orbit"),
      nb::arg("ref_species"),
      R"""(Wrapper for Fortran routine rel_tracking_charge_to_mass

Parameters
----------
orbit : CoordStruct
    Particle position structure.

ref_species : int
    Reference species

Returns
-------
rel_charge : float
    Relative charge/mass
)"""
  );
  m.def(
      "relative_mode_flip",
      &Bmad::relative_mode_flip,
      nb::arg("ele1"),
      nb::arg("ele2"),
      R"""(Wrapper for Fortran routine relative_mode_flip

Parameters
----------
ele1 : EleStruct
    Elements to compare.

ele2 : EleStruct
    Elements to compare.
)"""
  );
  nb::class_<PyReleaseRadIntCache>(m, "ReleaseRadIntCache", "release_rad_int_cache return type")
      .def_ro("ix_cache", &PyReleaseRadIntCache::ix_cache)
      .def("__len__", [](const PyReleaseRadIntCache &) { return 1; })
      .def("__getitem__", [](const PyReleaseRadIntCache &s, int i) -> nb::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return nb::cast(s.ix_cache);
        throw nb::index_error();
      });
  m.def(
      "release_rad_int_cache",
      &python_release_rad_int_cache,
      nb::arg("ix_cache"),
      R"""(Subroutine to release the memory associated with caching wiggler values.
See the radiation_integrals routine for further details.

Parameters
----------
ix_cache : int
    Cache number.
    This parameter is an input/output and is modified in-place.
    As an output, ix_cache: Cache number set to 0,

Returns
-------
ix_cache : int
    Cache number.
    This parameter is an input/output and is modified in-place.
    As an output, ix_cache: Cache number set to 0,
)"""
  );
  m.def(
      "remove_constant_taylor",
      &Bmad::remove_constant_taylor,
      nb::arg("taylor_in"),
      nb::arg("taylor_out"),
      nb::arg("c0"),
      nb::arg("remove_higher_order_terms"),
      R"""(Subroutine to remove the constant part of a taylor map.
Optionally terms that are higher order than bmad_com%taylor_order can
be removed.

Note: It is assumed that taylor_out has been deallocated before the call to
this routine. Calling this routine with the first two actual arguments the
same is prohibited.

Parameters
----------
taylor_in : 1D array of TaylorStruct
    Input taylor map.

taylor_out : 1D array of TaylorStruct
    Taylor with constant terms removed.

c0 : 1D array of float
    The constant part of the taylor map

remove_higher_order_terms : bool
    If True then terms that are higher order than bmad_com.taylor_order are removed.
)"""
  );
  m.def(
      "remove_dead_from_bunch",
      &Bmad::remove_dead_from_bunch,
      nb::arg("bunch_in"),
      R"""(Wrapper for Fortran routine remove_dead_from_bunch

Parameters
----------
bunch_in : BunchStruct
    Input bunch with alive and dead particles.

Returns
-------
bunch_out : BunchStruct
    Output bunch with only alive and pre_born particles. Note: bunch_out can be the same actual argument as
    bunch_in.
)"""
  );
  m.def(
      "remove_eles_from_lat",
      &Bmad::remove_eles_from_lat,
      nb::arg("lat"),
      nb::arg("check_sanity") = nb::none(),
      R"""(Wrapper for Fortran routine remove_eles_from_lat

Parameters
----------
lat : LatStruct
    Lattice to compress.
    This parameter is an input/output and is modified in-place.
    As an output, lat: Compressed lattice.

check_sanity : bool, optional
    If True (default) then call lat_sanity_check
)"""
  );
  m.def(
      "remove_lord_slave_link",
      &Bmad::remove_lord_slave_link,
      nb::arg("lord"),
      nb::arg("slave"),
      R"""(Wrapper for Fortran routine remove_lord_slave_link

Parameters
----------
lord : EleStruct
    Lord element
    This parameter is an input/output and is modified in-place.
    As an output, lord: Lord element with link info removed

slave : EleStruct
    Slave element
    This parameter is an input/output and is modified in-place.
    As an output, slave: Slave element with link info removed
)"""
  );
  m.def(
      "reverse_lat",
      &Bmad::reverse_lat,
      nb::arg("lat_in"),
      nb::arg("track_antiparticle") = nb::none(),
      R"""(Wrapper for Fortran routine reverse_lat

Parameters
----------
lat_in : LatStruct
    Input lattice to reverse.

track_antiparticle : bool, optional
    Set the particle species of the reversed lat to the anti-particle of lat_in? Default is True.

Returns
-------
lat_rev : LatStruct
    Reversed lattice.
)"""
  );
  m.def(
      "rf_clock_setup",
      &Bmad::rf_clock_setup,
      nb::arg("branch"),
      nb::arg("n_rf_included"),
      nb::arg("n_rf_excluded"),
      R"""(Wrapper for Fortran routine rf_clock_setup

Parameters
----------
branch : BranchStruct

n_rf_included : int

n_rf_excluded : int

Returns
-------
ok : bool
)"""
  );
  m.def(
      "rf_coupler_kick",
      &Bmad::rf_coupler_kick,
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("particle_at"),
      nb::arg("phase"),
      nb::arg("orbit"),
      nb::arg("mat6") = nb::none(),
      nb::arg("make_matrix") = nb::none(),
      R"""(Wrapper for Fortran routine rf_coupler_kick

Parameters
----------
ele : EleStruct
    Element being tracked through

param : LatParamStruct
    branch parameters.

particle_at : int
    first_track_edge$, or second_track_edge$.

phase : float
    phase of cavity

orbit : CoordStruct
    Position before kick.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Position after kick.

mat6 : 2D array of float (shape: 6,6), optional
    Transfer matrix before the element.
    This parameter is an input/output and is modified in-place.
    As an output, mat6: Transfer matrix through the element.

make_matrix : bool, optional
    Propagate the transfer matrix? Default is false.
)"""
  );
  m.def(
      "rf_is_on",
      &Bmad::rf_is_on,
      nb::arg("branch"),
      nb::arg("ix_ele1") = nb::none(),
      nb::arg("ix_ele2") = nb::none(),
      R"""(Wrapper for Fortran routine rf_is_on

Parameters
----------
branch : BranchStruct
    Lattice branch to check.

ix_ele1 : int, optional
    Start of range of elements to check. Default is 0.

ix_ele2 : int, optional
    End of range of elements to check. Default is branch.n_ele_track.

Returns
-------
is_on : bool
    True if any rfcavity is powered. False otherwise.
)"""
  );
  m.def(
      "rf_ref_time_offset",
      &Bmad::rf_ref_time_offset,
      nb::arg("ele"),
      nb::arg("ds") = nb::none(),
      R"""(Wrapper for Fortran routine rf_ref_time_offset

Parameters
----------
ele : EleStruct
    RF Element being tracked through.

ds : float, optional
    Distance of particle from start edge. Default is zero. Ouput:

Returns
-------
time : float
    Offset time.
)"""
  );
  m.def(
      "rfun",
      &Bmad::rfun,
      nb::arg("u"),
      nb::arg("v"),
      nb::arg("w"),
      nb::arg("gam"),
      nb::arg("a"),
      nb::arg("b"),
      nb::arg("hz"),
      nb::arg("i"),
      nb::arg("j"),
      R"""(Wrapper for Fortran routine rfun

Parameters
----------
u : float

v : float

w : float

gam : float

a : float

b : float

hz : float

i : int

j : int

Returns
-------
res : float
)"""
  );
  m.def(
      "rk_adaptive_time_step",
      [](EleStruct &ele,
         LatParamStruct &param,
         CoordStruct &orb,
         int t_dir,
         double rf_time,
         double dt_try,
         double dt_did,
         double dt_next,
         bool err_flag,
         EmFieldStruct *extra_field) {
        auto fn = static_cast<void (*)(
            EleStruct &,
            LatParamStruct &,
            CoordStruct &,
            int,
            double,
            double,
            double,
            double,
            bool,
            optional_ref<EmFieldStruct>
        )>(&Bmad::rk_adaptive_time_step);
        return fn(
            ele,
            param,
            orb,
            t_dir,
            rf_time,
            dt_try,
            dt_did,
            dt_next,
            err_flag,
            ptr_to_opt_ref(extra_field)
        );
      },
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("orb"),
      nb::arg("t_dir"),
      nb::arg("rf_time"),
      nb::arg("dt_try"),
      nb::arg("dt_did"),
      nb::arg("dt_next"),
      nb::arg("err_flag"),
      nb::arg("extra_field") = nb::none(),
      R"""(Wrapper for Fortran routine rk_adaptive_time_step

Parameters
----------
ele : EleStruct

param : LatParamStruct

orb : CoordStruct

t_dir : int

rf_time : float

dt_try : float

dt_did : float

dt_next : float

err_flag : bool

extra_field : EmFieldStruct, optional
)"""
  );
  m.def(
      "rk_time_step1",
      [](EleStruct &ele,
         LatParamStruct &param,
         double rf_time,
         CoordStruct &orb,
         double dt,
         CoordStruct &new_orb,
         bool err_flag,
         std::optional<FixedArray1D<Real, 10>> dr_dt,
         std::optional<bool> print_err,
         EmFieldStruct *extra_field) {
        auto fn = static_cast<FixedArray1D<Real, 10> (*)(
            EleStruct &,
            LatParamStruct &,
            double,
            CoordStruct &,
            double,
            CoordStruct &,
            bool,
            std::optional<FixedArray1D<Real, 10>>,
            std::optional<bool>,
            optional_ref<EmFieldStruct>
        )>(&Bmad::rk_time_step1);
        return fn(
            ele,
            param,
            rf_time,
            orb,
            dt,
            new_orb,
            err_flag,
            dr_dt,
            print_err,
            ptr_to_opt_ref(extra_field)
        );
      },
      nb::arg("ele"),
      nb::arg("param"),
      nb::arg("rf_time"),
      nb::arg("orb"),
      nb::arg("dt"),
      nb::arg("new_orb"),
      nb::arg("err_flag"),
      nb::arg("dr_dt") = nb::none(),
      nb::arg("print_err") = nb::none(),
      nb::arg("extra_field") = nb::none(),
      R"""(Wrapper for Fortran routine rk_time_step1

Parameters
----------
ele : EleStruct

param : LatParamStruct

rf_time : float

orb : CoordStruct

dt : float

new_orb : CoordStruct

err_flag : bool

dr_dt : 1D array of float (shape: 10), optional

print_err : bool, optional

extra_field : EmFieldStruct, optional

Returns
-------
r_err : 1D array of float (shape: 10)
)"""
  );
  m.def(
      "rotate3",
      &Bmad::rotate3,
      nb::arg("vec"),
      nb::arg("angle"),
      R"""(Wrapper for Fortran routine rotate3

Parameters
----------
vec : 1D array of float (shape: 3)

angle : float

Returns
-------
rvec : 1D array of float (shape: 3)
)"""
  );
  m.def(
      "rotate_em_field",
      &Bmad::rotate_em_field,
      nb::arg("field"),
      nb::arg("w_mat"),
      nb::arg("w_inv"),
      nb::arg("calc_dfield") = nb::none(),
      nb::arg("calc_potential") = nb::none(),
      R"""(Routine to transform the fields using the given rotation matrices.

Parameters
----------
field : EmFieldStruct
    E and B fields and derivatives.

w_mat : 2D array of float (shape: 3,3)
    rotation matrix.

w_inv : 2D array of float (shape: 3,3)
    rotation matrix inverse = transpose(w_mat)

calc_dfield : bool, optional
    If present and True then rotate the field derivatives.

calc_potential : bool, optional
    Rotate the magnetic vector potential? Default is false.
)"""
  );
  m.def(
      "rotate_field_zx",
      &Bmad::rotate_field_zx,
      nb::arg("field"),
      nb::arg("theta"),
      R"""(Wrapper for Fortran routine rotate_field_zx

Parameters
----------
field : EmFieldStruct

theta : float
)"""
  );
  m.def(
      "rotate_for_curved_surface",
      &Bmad::rotate_for_curved_surface,
      nb::arg("ele"),
      nb::arg("orbit"),
      nb::arg("set"),
      nb::arg("rot_mat"),
      R"""(Wrapper for Fortran routine rotate_for_curved_surface

Parameters
----------
ele : EleStruct
    reflecting element

orbit : CoordStruct
    Photon position.

set : bool
    True -> Transform body coords to local curved body coords. False -> Transform local curved body to body
    coords.

rot_mat : 2D array of float (shape: 3,3)
    When set = False, rotation matrix calculated from previous call with set = True.
    This parameter is an input/output and is modified in-place.
    As an output, rot_mat: When set = True, calculated rotation matrix.
)"""
  );
  m.def(
      "rotate_spin",
      &Bmad::rotate_spin,
      nb::arg("rot_vec"),
      nb::arg("spin"),
      R"""(Wrapper for Fortran routine rotate_spin

Parameters
----------
rot_vec : 1D array of float (shape: 3)
    Rotation axis. Magnitude of rot_vec is the rotation angle.

spin : 1D array of float (shape: 3)
    Initial coords.
    This parameter is an input/output and is modified in-place.
    As an output, spin: Final coords.

Returns
-------
qrot : 1D array of float (shape: 0:3), optional
    : rotation quaternion.
)"""
  );
  m.def(
      "rotate_spin_a_step",
      &Bmad::rotate_spin_a_step,
      nb::arg("orbit"),
      nb::arg("field"),
      nb::arg("ele"),
      nb::arg("ds"),
      R"""(Wrapper for Fortran routine rotate_spin_a_step

Parameters
----------
orbit : CoordStruct
    Initial orbit.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Orbit with rotated spin

field : EmFieldStruct
    EM Field

ele : EleStruct
    ele_struct, Element being tracked through.

ds : float
    Longitudinal step in element body frame.
)"""
  );
  m.def(
      "rotate_spin_given_field",
      &Bmad::rotate_spin_given_field,
      nb::arg("orbit"),
      nb::arg("sign_z_vel"),
      nb::arg("BL") = nb::none(),
      nb::arg("EL") = nb::none(),
      nb::arg("qrot") = nb::none(),
      R"""(Wrapper for Fortran routine rotate_spin_given_field

Parameters
----------
orbit : CoordStruct
    Initial orbit.
    This parameter is an input/output and is modified in-place.
    As an output, orbit: Orbit with rotated spin

sign_z_vel : int
    +/- 1. Sign of direction of travel relative to the element.

BL : 1D array of float (shape: 3), optional
    Integrated field strength. Assumed zero if not present.

EL : 1D array of float (shape: 3), optional
    Integrated field strength. Assumed zero if not present.

qrot : 1D array of float (shape: 0:3), optional
    Initial rotation quaternion.
    This parameter is an input/output and is modified in-place.
    As an output, qrot: Rotation quaternion with rotation due to the field added in.
)"""
  );
}
