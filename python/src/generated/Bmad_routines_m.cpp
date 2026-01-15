#include "pybmad/generated/Bmad_routines_m.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Bmad_routines_m(py::module& m) {
  m.def(
      "mad_add_offsets_and_multipoles",
      &Bmad::mad_add_offsets_and_multipoles,
      py::arg("ele"),
      R"""(Subroutine mad_add_offsets_and_multipoles (ele, map)

Subroutine to add in the effect of element offsets and/or multipoles
on the 2nd order transport map for the element.

Parameters
----------
ele : EleStruct
    Drift element.
energy : MadEnergyStruct
    particle energy structure.

Returns
-------
map : MadMapStruct
    Structure holding the transfer map.
)""");
  m.def(
      "mad_concat_map2",
      &Bmad::mad_concat_map2,
      py::arg("map1"),
      py::arg("map2"),
      R"""(Subroutine mad_concat_map2 (map1, map2, map3)

Subroutine to concatinate two 2nd order transport maps.
    map3 = map2(map1)
The equivalent MAD-8 routine is: TMCAT1

Parameters
----------
map1 : MadMapStruct
    First map in the beam line.
map2 : MadMapStruct
    Second map in the beam line.

Returns
-------
map3 : MadMapStruct
    Concatinated map.
)""");
  m.def(
      "mad_drift",
      &Bmad::mad_drift,
      py::arg("ele"),
      py::arg("energy"),
      R"""(Subroutine mad_drift (ele, energy, map)

Subroutine to make a transport map for a drift space.
The equivalent MAD-8 routine is: TMDRF

Parameters
----------
ele : EleStruct
    Drift element.
energy : MadEnergyStruct
    particle energy structure.

Returns
-------
map : MadMapStruct
    Structure holding the transfer map.
)""");
  m.def(
      "mad_elsep",
      &Bmad::mad_elsep,
      py::arg("ele"),
      py::arg("energy"),
      R"""(Subroutine mad_elsep (ele, energy, map)

Subroutine to make a transport map for an electric separator.
The equivalent MAD-8 routine is: TMSEP

Parameters
----------
ele : EleStruct
    Electric seperator element.
energy : MadEnergyStruct
    particle energy structure.

Returns
-------
map : MadMapStruct
    Structure holding the transfer map.
)""");
  m.def(
      "mad_map_to_taylor",
      &Bmad::mad_map_to_taylor,
      py::arg("map"),
      py::arg("energy"),
      py::arg("taylor"),
      R"""(Subroutine mad_map_to_taylor (map, energy, taylor)

Subroutine to convert a MAD order 2 map to a Bmad taylor map.
The conversion will also convert between MAD's (t, dE) and Bmad's (beta*t, dP) coords.

Parameters
----------
map : MadMapStruct
    Order 2 map.
energy : MadEnergyStruct
    Energy numbers.

Returns
-------
taylor : TaylorStruct
    Taylor map.
)""");
  m.def(
      "mad_quadrupole",
      &Bmad::mad_quadrupole,
      py::arg("ele"),
      py::arg("energy"),
      R"""(Subroutine mad_quadrupole (ele, energy, map)

Subroutine to make a transport map for an quadrupole element.
The equivalent MAD-8 routine is: TMSEXT

Parameters
----------
ele : EleStruct
    Quadrupole element.
energy : MadEnergyStruct
    particle energy structure.

Returns
-------
map : MadMapStruct
    Structure holding the transfer map.
)""");
  m.def(
      "mad_rfcavity",
      &Bmad::mad_rfcavity,
      py::arg("ele"),
      py::arg("energy"),
      R"""(Subroutine mad_rfcavity (ele, energy, map)

Subroutine to make a transport map for an rfcavity element.
The equivalent MAD-8 routine is: TMRF

Parameters
----------
ele : EleStruct
    Rfcavity element.
energy : MadEnergyStruct
    particle energy structure.

Returns
-------
map : MadMapStruct
    Structure holding the transfer map.
)""");
  m.def(
      "mad_sbend",
      &Bmad::mad_sbend,
      py::arg("ele"),
      py::arg("energy"),
      R"""(Subroutine mad_sbend (ele, energy, map)

Subroutine to make a transport map for a sector bend element.
The equivalent MAD-8 routine is: TMBEND

Parameters
----------
ele : EleStruct
    Sbend element.
energy : MadEnergyStruct
    particle energy structure.

Returns
-------
map : MadMapStruct
    Structure holding the transfer map.
)""");
  m.def(
      "mad_sbend_body",
      &Bmad::mad_sbend_body,
      py::arg("ele"),
      py::arg("energy"),
      R"""(Subroutine mad_sbend_body (ele, energy, map)

Subroutine to make a transport map for the body of a sector dipole.
The equivalent MAD-8 routine is: TMSECT

Parameters
----------
ele : EleStruct
    Solenoid element.
energy : MadEnergyStruct
    particle energy structure.
into : bool
    If True then map is for particle entering a dipole

Returns
-------
map : MadMapStruct
    Structure holding the transfer map.
)""");
  m.def(
      "mad_sbend_fringe",
      &Bmad::mad_sbend_fringe,
      py::arg("ele"),
      py::arg("energy"),
      py::arg("into"),
      R"""(Subroutine mad_sbend_fringe (ele, energy, into, map)

Subroutine to make a transport map for the fringe field of a dipole.
The equivalent MAD-8 routine is: TMFRNG

Parameters
----------
ele : EleStruct
    Solenoid element.
energy : MadEnergyStruct
    particle energy structure.
into : bool
    If True then map is for particle entering a dipole

Returns
-------
map : MadMapStruct
    Fringe dipole map. .k(6)     -- 0th order map. .r(6,6)   -- 1st order map. .t(6,6,6) -- 2nd order map.
)""");
  m.def(
      "mad_sextupole",
      &Bmad::mad_sextupole,
      py::arg("ele"),
      py::arg("energy"),
      R"""(Subroutine mad_sextupole (ele, energy, map)

Subroutine to make a transport map for an sextupole.
The equivalent MAD-8 routine is: TMSEXT

Parameters
----------
ele : EleStruct
    Sextupole element.
energy : MadEnergyStruct
    particle energy structure.

Returns
-------
map : MadMapStruct
    Structure holding the transfer map.
)""");
  m.def(
      "mad_solenoid",
      &Bmad::mad_solenoid,
      py::arg("ele"),
      py::arg("energy"),
      R"""(Subroutine mad_solenoid (ele, energy, map)

Subroutine to make a transport map for an solenoid.
The equivalent MAD-8 routine is: TMSEXT

Parameters
----------
ele : EleStruct
    Solenoid element.
energy : MadEnergyStruct
    particle energy structure.

Returns
-------
map : MadMapStruct
    Structure holding the transfer map.
)""");
  py::class_<Bmad::MadTmfoc, std::unique_ptr<Bmad::MadTmfoc>>(
      m, "MadTmfoc", "mad_tmfoc return type")
      .def_readonly("c", &Bmad::MadTmfoc::c)
      .def_readonly("s", &Bmad::MadTmfoc::s)
      .def_readonly("d", &Bmad::MadTmfoc::d)
      .def_readonly("f", &Bmad::MadTmfoc::f)
      .def("__len__", [](const Bmad::MadTmfoc&) { return 4; })
      .def("__getitem__", [](const Bmad::MadTmfoc& s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.c);
        if (i == 1)
          return py::cast(s.s);
        if (i == 2)
          return py::cast(s.d);
        if (i == 3)
          return py::cast(s.f);
        throw py::index_error();
      });
  m.def(
      "mad_tmfoc",
      &Bmad::mad_tmfoc,
      py::arg("el"),
      py::arg("sk1"),
      R"""(Subroutine mad_tmfoc (el, sk1, c, s, d, f)

Subroutine to compute the linear focussing functions.
The equivalent MAD-8 routine is: TMFOC

Parameters
----------
el : float
    Length.
sk1 : float
    Quadrupole strength.

Returns
-------
c : float
    Cosine-like function.             c(k,l)
s : float
    Sine-like function.               s(k,l)
d : float
    Dispersion function.              d(k,l)
f : float
    Integral of dispersion function.  f(k,l)
)""");
  m.def(
      "mad_tmsymm",
      &Bmad::mad_tmsymm,
      py::arg("te"),
      R"""(subroutine mad_tmsymm (te)

subroutine to symmertrize the 2nd order map t.
The equivalent MAD-8 routine is: tmsymm

Parameters
----------
te : float
    array to be symmertrized.
    This parameter is an input/output and is modified in-place. As an output: symmetrized array.
)""");
  m.def(
      "mad_tmtilt",
      &Bmad::mad_tmtilt,
      py::arg("map"),
      py::arg("tilt"),
      R"""(Subroutine mad_tmtilt (map, tilt)

Subroutine to apply a tilt to a transport map.
The equivalent MAD-8 routine is: TMTILT

Parameters
----------
map : MadMapStruct
    Unrotated transport map. .k(6)     -- 0th order map. .r(6,6)   -- 1st order map. .t(6,6,6) -- 2nd order
    map.
    This parameter is an input/output and is modified in-place. As an output: Rotated transport map.
tilt : float
    Tilt
)""");
  m.def(
      "mad_track1",
      &Bmad::mad_track1,
      py::arg("c0"),
      py::arg("map"),
      R"""(Subroutine mad_track1 (c0, map, c1)

Subroutine to track through a 2nd order transfer map.
The equivalent MAD-8 routine is: TMTRAK

Parameters
----------
c0 : CoordStruct
    Starting coords.
map : MadMapStruct
    2nd order map.

Returns
-------
c1 : CoordStruct
    Ending coords.
)""");
  m.def(
      "make_g2_mats",
      &Bmad::make_g2_mats,
      py::arg("twiss"),
      py::arg("g2_mat"),
      py::arg("g2_inv_mat"),
      R"""(Parameters
----------
twiss : TwissStruct
    Twiss parameters.
g2_mat : 
g2_inv_mat : 
)""");
  py::class_<Bmad::MakeGMats, std::unique_ptr<Bmad::MakeGMats>>(
      m, "MakeGMats", "make_g_mats return type")
      .def_readonly("g_mat", &Bmad::MakeGMats::g_mat)
      .def_readonly("g_inv_mat", &Bmad::MakeGMats::g_inv_mat)
      .def("__len__", [](const Bmad::MakeGMats&) { return 2; })
      .def("__getitem__", [](const Bmad::MakeGMats& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.g_mat);
        if (i == 1)
          return py::cast(s.g_inv_mat);
        throw py::index_error();
      });
  m.def(
      "make_g_mats",
      &Bmad::make_g_mats,
      py::arg("ele"),
      R"""(Parameters
----------
ele : EleStruct
    Element
g_mat : float
    Normal mode to betaless coords
g_inv_mat : float
    The inverse of G_MAT
)""");
  py::class_<Bmad::MakeHvbp, std::unique_ptr<Bmad::MakeHvbp>>(
      m, "MakeHvbp", "make_hvbp return type")
      .def_readonly("B", &Bmad::MakeHvbp::B)
      .def_readonly("V", &Bmad::MakeHvbp::V)
      .def_readonly("H", &Bmad::MakeHvbp::H)
      .def_readonly("Vbar", &Bmad::MakeHvbp::Vbar)
      .def_readonly("Hbar", &Bmad::MakeHvbp::Hbar)
      .def("__len__", [](const Bmad::MakeHvbp&) { return 5; })
      .def("__getitem__", [](const Bmad::MakeHvbp& s, int i) -> py::object {
        if (i < 0)
          i += 5;
        if (i == 0)
          return py::cast(s.B);
        if (i == 1)
          return py::cast(s.V);
        if (i == 2)
          return py::cast(s.H);
        if (i == 3)
          return py::cast(s.Vbar);
        if (i == 4)
          return py::cast(s.Hbar);
        throw py::index_error();
      });
  m.def(
      "make_hvbp",
      &Bmad::make_hvbp,
      py::arg("N"),
      R"""(Subroutine make_HVBP(N, B, V, H, Vbar, Hbar)

Parameterizes the eigen-decomposition of the 6x6 transfer matrix into HVBP as defined in:
"From the beam-envelop matrix to synchrotron-radiation integrals" by Ohmi, Hirata, and Oide.

This routine takes N, which is usually made from make_N (also in this module), and decomposes
it into H, V, B, and P.

N is defined by:
M = N.U.Inverse[N] where U is block diagonal and the blocks are 2x2 rotation matrices.
and it is decomposed by this subroutine as,
N = H.V.B.P
P has the same free parameters as B
B "Twiss matrix" has 6 free parameters (Twiss alphas and betas)
B blocks have the form /     sqrt(beta)         0       \
                       \ -alpha/sqrt(beta) 1/sqrt(beta) /
V "Teng matrix" has 4 free parameters (xy, xpy, ypx, and pxpy coupling)
H "Dispersion matrix" has 8 free parameters (xz, xpz, pxz, pxpz, yz, ypz, pyz, pypz coupling)


Parameters
----------
N : float
    Matrix of eigenvectors prepared by make_N

Returns
-------
B : float
    Block diagonal matrix of Twiss parameters
V : float
    horizontal-vertical coupling information
H : float
    horizontal-longitudinal and vertical-longitudinal coupling information
Vbar : float
    mat_symp_conj(B).V.B
Hbar : float
    mat_symp_conj(B).H.B
)""");
  m.def(
      "make_hybrid_lat",
      &Bmad::make_hybrid_lat,
      py::arg("lat_in"),
      py::arg("use_taylor") = py::none(),
      py::arg("orb0_arr") = py::none(),
      R"""(Parameters
----------
lat_in : LatStruct
    Input lattice. .branch(:).ele(:).select  -- Roughly: Set True to keep and False to hybridize. See above.
lat_out : LatStruct
    Lattice with hybrid elements. Note: Lat_out must not be the same actual argument as lat_in.
use_taylor : bool, optional
    If present and True then the hybrid elements will have a taylor series instead of a simple linear matrix.
    If an element to be concatenated has a taylor series then this taylor series will be concatenated with the
    other elements in the hybrid element.
orb0_arr : CoordArrayStruct, optional
    Central orbit for taylor stuff. Each orb0_arr(i).orbit(:) holds the orbit for the i^th lattice branch
)""");
  py::class_<Bmad::MakeMadMap, std::unique_ptr<Bmad::MakeMadMap>>(
      m, "MakeMadMap", "make_mad_map return type")
      .def_readonly("energy", &Bmad::MakeMadMap::energy)
      .def_readonly("map", &Bmad::MakeMadMap::map)
      .def("__len__", [](const Bmad::MakeMadMap&) { return 2; })
      .def("__getitem__", [](const Bmad::MakeMadMap& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.energy);
        if (i == 1)
          return py::cast(s.map);
        throw py::index_error();
      });
  m.def(
      "make_mad_map",
      &Bmad::make_mad_map,
      py::arg("ele"),
      py::arg("param"),
      R"""(Subroutine make_mad_map (ele, param, energy, map)

Subroutine to make a 2nd order transport map a la MAD.

Parameters
----------
ele : EleStruct
    Element
param : LatParamStruct
    particle id

Returns
-------
energy : MadEnergyStruct
    Energy of the particle
map : MadMapStruct
    Structure holding the transfer map.
)""");
  py::class_<Bmad::MakeMat6, std::unique_ptr<Bmad::MakeMat6>>(
      m, "MakeMat6", "make_mat6 return type")
      .def_readonly("end_orb", &Bmad::MakeMat6::end_orb)
      .def_readonly("err_flag", &Bmad::MakeMat6::err_flag)
      .def("__len__", [](const Bmad::MakeMat6&) { return 2; })
      .def("__getitem__", [](const Bmad::MakeMat6& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.end_orb);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "make_mat6",
      &Bmad::make_mat6,
      py::arg("ele"),
      py::arg("param"),
      py::arg("start_orb") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Element holding the transfer matrix.
param : LatParamStruct
    Lattice global parameters.
start_orb : CoordStruct, optional
    Reference coordinates at the beginning of element. If not present, default is to use the zero orbit.
end_orb : CoordStruct
    Reference coordinates at the end of element.
err_flag : bool
    Set True if there is an error. False otherwise.
)""");
  py::class_<Bmad::MakeMat6Bmad, std::unique_ptr<Bmad::MakeMat6Bmad>>(
      m, "MakeMat6Bmad", "make_mat6_bmad return type")
      .def_readonly("end_orb", &Bmad::MakeMat6Bmad::end_orb)
      .def_readonly("err", &Bmad::MakeMat6Bmad::err)
      .def("__len__", [](const Bmad::MakeMat6Bmad&) { return 2; })
      .def("__getitem__", [](const Bmad::MakeMat6Bmad& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.end_orb);
        if (i == 1)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "make_mat6_bmad",
      &Bmad::make_mat6_bmad,
      py::arg("ele"),
      py::arg("param"),
      py::arg("start_orb"),
      R"""(Parameters
----------
ele : EleStruct
    Element to track through.
    This parameter is an input/output and is modified in-place. As an output: Element with transfer matrix.
param : LatParamStruct
    Parameters are needed for some elements.
start_orb : CoordStruct
    Starting coords.
end_orb : CoordStruct
    Coordinates at the end of element.
err : bool
    Set True if there is an error. False otherwise.
)""");
  py::class_<
      Bmad::MakeMat6BmadPhoton,
      std::unique_ptr<Bmad::MakeMat6BmadPhoton>>(
      m, "MakeMat6BmadPhoton", "make_mat6_bmad_photon return type")
      .def_readonly("end_orb", &Bmad::MakeMat6BmadPhoton::end_orb)
      .def_readonly("err", &Bmad::MakeMat6BmadPhoton::err)
      .def("__len__", [](const Bmad::MakeMat6BmadPhoton&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::MakeMat6BmadPhoton& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.end_orb);
            if (i == 1)
              return py::cast(s.err);
            throw py::index_error();
          });
  m.def(
      "make_mat6_bmad_photon",
      &Bmad::make_mat6_bmad_photon,
      py::arg("ele"),
      py::arg("param"),
      py::arg("start_orb"),
      R"""(Parameters
----------
ele : EleStruct
    Element with transfer matrix
    This parameter is an input/output and is modified in-place. As an output: Element with transfer matrix.
param : LatParamStruct
    Parameters are needed for some elements.
start_orb : CoordStruct
    Coordinates at the beginning of element.
end_orb : CoordStruct
    Coordinates at the end of element.
err : bool
    Set True if there is an error. False otherwise.
)""");
  m.def(
      "make_mat6_high_energy_space_charge",
      &Bmad::make_mat6_high_energy_space_charge,
      py::arg("ele"),
      py::arg("param"),
      R"""(Subroutine make_mat6_high_energy_space_charge (ele, param)

Routine to add the ultra relativistic space charge kick to the element transfer matrix.
The routine setup_space_charge_calc must be called
initially before any tracking is done. This routine assumes a Gaussian
bunch and is only valid with relativistic particles where the effect
of the space charge is small.

Parameters
----------
ele : EleStruct
    Element tracked through.
param : LatParamStruct

Returns
-------
end : CoordStruct
    End position
)""");
  m.def(
      "make_mat6_mad",
      &Bmad::make_mat6_mad,
      py::arg("ele"),
      py::arg("param"),
      py::arg("c0"),
      R"""(Subroutine make_mat6_mad (ele, param, c0, c1)

Subroutine to make the 6x6 transfer matrix for an element from the
2nd order MAD transport map. The map is stored in ele%taylor.
If the map exists then it is simply used to calculate ele%mat6.
If ele%taylor doesn't exist then calculate it.

Parameters
----------
ele : EleStruct
    Element with transfer matrix.
param : LatParamStruct
    Lattice parameters.
map : MadMapStruct
    2nd order map.
c0 : CoordStruct
    Coordinates at the beginning of element.

Returns
-------
c1 : CoordStruct
    Coordinates at the end of element.
)""");
  m.def(
      "make_mat6_symp_lie_ptc",
      &Bmad::make_mat6_symp_lie_ptc,
      py::arg("ele"),
      py::arg("start_orb"),
      R"""(Parameters
----------
ele : EleStruct
    Element with transfer matrix
    This parameter is an input/output and is modified in-place. As an output: Element with transfer matrix.
start_orb : CoordStruct
    Coordinates at the beginning of element.
end_orb : CoordStruct
    Coordinates at end of element.
)""");
  m.def(
      "make_mat6_taylor",
      &Bmad::make_mat6_taylor,
      py::arg("ele"),
      py::arg("start_orb"),
      py::arg("err_flag") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Element to track through.
    This parameter is an input/output and is modified in-place. As an output: Element with transfer matrix.
start_orb : CoordStruct
    Starting coords.
end_orb : CoordStruct
    Coordinates at the end of element.
err_flag : 
)""");
  py::class_<Bmad::MakeMat6Tracking, std::unique_ptr<Bmad::MakeMat6Tracking>>(
      m, "MakeMat6Tracking", "make_mat6_tracking return type")
      .def_readonly("end_orb", &Bmad::MakeMat6Tracking::end_orb)
      .def_readonly("err_flag", &Bmad::MakeMat6Tracking::err_flag)
      .def("__len__", [](const Bmad::MakeMat6Tracking&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::MakeMat6Tracking& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.end_orb);
            if (i == 1)
              return py::cast(s.err_flag);
            throw py::index_error();
          });
  m.def(
      "make_mat6_tracking",
      &Bmad::make_mat6_tracking,
      py::arg("ele"),
      py::arg("param"),
      py::arg("start_orb"),
      py::arg("spin_only") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Element with transfer matrix
    This parameter is an input/output and is modified in-place. As an output: Element with transfer matrix.
param : LatParamStruct
    Parameters are needed for some elements.
start_orb : CoordStruct
    Coordinates at the beginning of element.
end_orb : CoordStruct
    Coordinates at the end of element.
err_flag : bool
    Set True if there is an error. False otherwise.
spin_only : bool, optional
    Default False. If True, only calculate ele.spin_taylor.
)""");
  py::class_<Bmad::MakeN, std::unique_ptr<Bmad::MakeN>>(
      m, "MakeN", "make_n return type")
      .def_readonly("N", &Bmad::MakeN::N)
      .def_readonly("err_flag", &Bmad::MakeN::err_flag)
      .def_readonly("tunes_out", &Bmad::MakeN::tunes_out)
      .def_readonly("U", &Bmad::MakeN::U)
      .def("__len__", [](const Bmad::MakeN&) { return 4; })
      .def("__getitem__", [](const Bmad::MakeN& s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.N);
        if (i == 1)
          return py::cast(s.err_flag);
        if (i == 2)
          return py::cast(s.tunes_out);
        if (i == 3)
          return py::cast(s.U);
        throw py::index_error();
      });
  m.def(
      "make_n",
      &Bmad::make_n,
      py::arg("t6"),
      py::arg("abz_tunes") = py::none(),
      R"""(Subroutine make_N(t6, N, err_flag, abz_tunes, tunes_out, U)

Given a 1-turn transfer matrix, this returns the matrix N.
N converts between normal invarients and phases and canonical coordinates:
X = N.J

N is obtained from the Eigen decomposition of the 1-turn transfer matrix.
It is obtained by applying certain normalizations to the matrix of Eigen vectors, then making
the result real using Q.

If abz_tunes is present, then the eigensystem is ordered by matching the tunes.
If abz_tunes is not present, then the eigensystem is ordered by plane dominance.

It is assumed that the synchrotron tune is less than pi.

Parameters
----------
t6 : float
    1-turn transfer matrix
abz_tunes : float, optional
    a-mode is abz_tunes(1), b-mode is abz_tunes(2), synch tune is abz_tunes(3)

Returns
-------
N : float
    X = N.J
err_flag : bool
    Set to true on error.  Often means Eigen decomposition failed.
tunes_out : float
    Fractional tune (in radians) of the 3 normal modes of t6.
U : float
    U = Inverse(N).t6.N.  Block diagonal matrix of 2x2 rotation matrices.
)""");
  py::class_<Bmad::MakePbrh, std::unique_ptr<Bmad::MakePbrh>>(
      m, "MakePbrh", "make_pbrh return type")
      .def_readonly("P", &Bmad::MakePbrh::P)
      .def_readonly("Bp", &Bmad::MakePbrh::Bp)
      .def_readonly("R", &Bmad::MakePbrh::R)
      .def_readonly("H", &Bmad::MakePbrh::H)
      .def("__len__", [](const Bmad::MakePbrh&) { return 4; })
      .def("__getitem__", [](const Bmad::MakePbrh& s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.P);
        if (i == 1)
          return py::cast(s.Bp);
        if (i == 2)
          return py::cast(s.R);
        if (i == 3)
          return py::cast(s.H);
        throw py::index_error();
      });
  m.def(
      "make_pbrh",
      &Bmad::make_pbrh,
      py::arg("M"),
      py::arg("abz_tunes"),
      R"""(subroutine make_PBRH(M, P, Bp, R, H, abz_tunes)

Decomposes the 1-turn transfer matrix into normal mode twiss-like parameters,
according to Sec. IIIB of Ohmi, Hirata, and Oide paper.

Note:  The Twiss parameters generated by this function are identical to those delivered
       by mode3_mod.

Parameters
----------
M : float
    1-turn transfer matrix
abz_tunes : float
    tunes for a,b, and c modes.  Used to identify which eigenvector is associated with which mode.

Returns
-------
P : complex
    Eqn. 97.  Phase advances.
Bp : complex
    Eqns. 89 & 101.  Beta functions.
R : complex
    Eqn. 99.  Transverse coupling.
H : complex
    Eqn. 100.  Longitudinal coupling.
)""");
  py::class_<Bmad::MakeSmatFromAbc, std::unique_ptr<Bmad::MakeSmatFromAbc>>(
      m, "MakeSmatFromAbc", "make_smat_from_abc return type")
      .def_readonly("sigma_mat", &Bmad::MakeSmatFromAbc::sigma_mat)
      .def_readonly("err_flag", &Bmad::MakeSmatFromAbc::err_flag)
      .def_readonly("Nout", &Bmad::MakeSmatFromAbc::Nout)
      .def("__len__", [](const Bmad::MakeSmatFromAbc&) { return 3; })
      .def(
          "__getitem__",
          [](const Bmad::MakeSmatFromAbc& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.sigma_mat);
            if (i == 1)
              return py::cast(s.err_flag);
            if (i == 2)
              return py::cast(s.Nout);
            throw py::index_error();
          });
  m.def(
      "make_smat_from_abc",
      &Bmad::make_smat_from_abc,
      py::arg("t6"),
      py::arg("mode"),
      R"""(Subroutine make_smat_from_abc(t6, mode, sigma_mat, err_flag, Nout)

Given the 1-turn transfer matrix and a normal_modes_struct containing the normal mode
emittances, this routine returns the beam envelop sigma matrix.

sigma_mat = N.D.transpose(N)
equivalent to: sigma_mat.S = N.D.mat_symp_conj(N)

One way to populate mode%a%tune and mode%b%tune:
  mode%a%tune = mod(lat%ele(lat%n_ele_track)%a%phi, twopi)
  mode%b%tune = mod(lat%ele(lat%n_ele_track)%b%phi, twopi)

Parameters
----------
t6 : float
    1-turn transfer matrix
mode : NormalModesStruct
    normal mode emittances .a.emittance -- real(rp): a-mode emittance .b.emittance -- real(rp): b-mode
    emittance .z.emittance -- real(rp): z-mode emittance .a.tune      -- real(rp): a-mode tune.  Used to
    associate emittances with the proper mode. .b.tune      -- real(rp): b-mode tune.  Used to associate
    emittances with the proper mode. .z.tune      -- real(rp): z-mode tune.  Used to associate emittances with
    the proper mode.

Returns
-------
sigma_mat : float
    beam envelop sigma matrix
err_flag : bool
    set to true if something goes wrong.  Usually means Eigen decomposition of the 1-turn matrix failed.
Nout : float
    Contains the normalized eigenvectors that were used to make the sigma matrix.
)""");
  m.def(
      "make_unit_mad_map",
      &Bmad::make_unit_mad_map,
      py::arg("map"),
      R"""(Subroutine make_unit_mad_map (map)

Subroutine to initialize a 2nd order transport map to unity.

Parameters
----------
map : MadMapStruct
    2nd order transport map.
    This parameter is an input/output and is modified in-place. As an output: Unity 2nd order map.
)""");
  m.def(
      "make_v",
      &Bmad::make_v,
      py::arg("M"),
      py::arg("V"),
      py::arg("abz_tunes"),
      R"""(subroutine make_V(M,V,abz_tunes)

For a one-turn transfer matrix M, this routine find the eigen matrix V.
V is ordered such that the per turn phase advance of its column pairs agree with abz_tunes.
It is normalized to be symplectic.

)""");
  py::class_<Bmad::MakeVMats, std::unique_ptr<Bmad::MakeVMats>>(
      m, "MakeVMats", "make_v_mats return type")
      .def_readonly("v_mat", &Bmad::MakeVMats::v_mat)
      .def_readonly("v_inv_mat", &Bmad::MakeVMats::v_inv_mat)
      .def("__len__", [](const Bmad::MakeVMats&) { return 2; })
      .def("__getitem__", [](const Bmad::MakeVMats& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.v_mat);
        if (i == 1)
          return py::cast(s.v_inv_mat);
        throw py::index_error();
      });
  m.def(
      "make_v_mats",
      &Bmad::make_v_mats,
      py::arg("ele"),
      R"""(Parameters
----------
ele : EleStruct
    Element
v_mat : float
    Normal mode to X-Y coords transformation
v_inv_mat : float
    X-Y coords to Normal mode transformation
)""");
  m.def(
      "makeup_control_slave",
      &Bmad::makeup_control_slave,
      py::arg("lat"),
      py::arg("slave"),
      py::arg("err_flag"),
      R"""(Subroutine makeup_control_slave (lat, slave, err_flag)

This routine is not meant for general use.

)""");
  m.def(
      "makeup_group_lord",
      &Bmad::makeup_group_lord,
      py::arg("lat"),
      py::arg("lord"),
      py::arg("err_flag"),
      R"""(Subroutine makeup_group_lord (lat, lord, err_flag)

Subroutine to calculate the attributes of group slave elements.
This routine is private to bookkeeper_mod.

)""");
  m.def(
      "makeup_multipass_slave",
      &Bmad::makeup_multipass_slave,
      py::arg("lat"),
      py::arg("slave"),
      py::arg("err_flag"),
      R"""(Subroutine makeup_multipass_slave (lat, slave, err_flag)

Subroutine to calcualte the attributes of multipass slave elements.
This routine is not meant for guse.

)""");
  m.def(
      "makeup_super_slave",
      &Bmad::makeup_super_slave,
      py::arg("lat"),
      py::arg("slave"),
      py::arg("err_flag"),
      R"""(Subroutine makeup_super_slave (lat, slave, err_flag)

Subroutine to calcualte the attributes of superposition slave elements.
This routine is not meant for general use.

)""");
  m.def(
      "makeup_super_slave1",
      &Bmad::makeup_super_slave1,
      py::arg("slave"),
      py::arg("lord"),
      py::arg("offset"),
      py::arg("param"),
      py::arg("include_upstream_end"),
      py::arg("include_downstream_end"),
      R"""(Subroutine makeup_super_slave1 (slave, lord, offset, param, include_upstream_end, include_downstream_end)

Routine to construct a super_slave from a super_lord when the slave has only one lord.
Note: Reference energy and times are not computed in this routine.

Parameters
----------
slave : EleStruct
    Slave element.
    This parameter is an input/output and is modified in-place. As an output: Slave element with appropriate
    values set.
lord : EleStruct
    Lord element.
offset : float
    offset of entrance end of slave from entrance end of the lord.
param : LatParamStruct
    lattice paramters.
include_upstream_end : bool
    Slave contains the lord's entrance end?
include_downstream_end : bool
    Slave contains the lord's exit end?

Returns
-------
err_flag : bool
    Set true if there is an error. False otherwise.
)""");
  m.def(
      "map1_inverse",
      &Bmad::map1_inverse,
      py::arg("map1"),
      py::arg("inv_map1"),
      R"""(Parameters
----------
map1 : SpinOrbitMap1Struct
    Input map.
inv_map1 : 
)""");
  m.def(
      "map1_make_unit",
      &Bmad::map1_make_unit,
      R"""(Parameters
----------
map1 : SpinOrbitMap1Struct
    Unit map.
)""");
  m.def(
      "map1_times_map1",
      &Bmad::map1_times_map1,
      py::arg("map2"),
      py::arg("map1"),
      py::arg("map_out"),
      R"""(Parameters
----------
map2 : 
map1 : 
map_out : 
)""");
  m.def(
      "map_to_angle_coords",
      &Bmad::map_to_angle_coords,
      py::arg("t_canon"),
      R"""(Parameters
----------
t_canon : TaylorStruct
    Taylor map in canonical coords.
t_angle : TaylorStruct
    Taylor map in angle coords.
)""");
  m.def(
      "mark_patch_regions",
      &Bmad::mark_patch_regions,
      py::arg("branch"),
      R"""(Subroutine mark_patch_regions (branch)

Routine to mark which regions in a wall3d structure contain patch elements.
This routine should be called by any routine that creates a beam chamber wall.

Parameters
----------
branch : BranchStruct
    Lattice branch with .wall3d beam chamber wall.
    This parameter is an input/output and is modified in-place. As an output: Lattice branch with
    .wall3d.section(i).patch_in_region marked.
)""");
  m.def(
      "master_parameter_value",
      &Bmad::master_parameter_value,
      py::arg("master_parameter"),
      py::arg("ele"),
      py::arg("value"),
      R"""(Parameters
----------
master_parameter : int
    Index of the master parameter.
ele : EleStruct
    Element containing the fieldmap.
value : 
)""");
  m.def(
      "mat4_multipole",
      &Bmad::mat4_multipole,
      py::arg("knl"),
      py::arg("tilt"),
      py::arg("n"),
      py::arg("orbit"),
      R"""(Parameters
----------
knl : float
    Strength of multipole
tilt : float
    Tilt of multipole
n : 
orbit : CoordStruct
    coordinates of particle
kick_mat : float
    Kick matrix (Jacobian) at orbit.
)""");
  m.def(
      "mat6_add_offsets",
      &Bmad::mat6_add_offsets,
      py::arg("ele"),
      py::arg("param"),
      R"""(Parameters
----------
ele : EleStruct
    Element with given orientation. .vec0(6)         -- 0th order part of the transfer map. .mat6(6,6)
    -- 1st order part of the transfer map (Jacobian). .map_ref_orb_in  -- Reference orbit at entrance end.
    .map_ref_orb_out -- Reference orbit at exit end. .value(x_offset$), .value(x_pitch$), .value(tilt$), etc.
    -- Offsets, tilts, and pitches
param : LatParamStruct
)""");
  m.def(
      "mat6_add_pitch",
      &Bmad::mat6_add_pitch,
      py::arg("x_pitch_tot"),
      py::arg("y_pitch_tot"),
      py::arg("orientation"),
      py::arg("mat6"),
      R"""(Parameters
----------
x_pitch_tot : float
    Horizontal pitch
y_pitch_tot : float
    Vertical pitch
orientation : int
    Element longitudinal orientation. +1 or -1.
mat6 : float
    1st order part of the transfer map (Jacobian).
    This parameter is an input/output and is modified in-place. As an output: 1st order xfer map with pitches.
)""");
  m.def(
      "mat6_to_complex_taylor",
      &Bmad::mat6_to_complex_taylor,
      py::arg("vec0"),
      py::arg("mat6"),
      R"""(Subroutine mat6_to_complex_taylor (vec0, mat6, complex_taylor)

Subroutine to form a first order complex_taylor map from the 6x6 transfer
matrix and the 0th order transfer vector.

Parameters
----------
vec0 : 
    0th order transfer vector.
mat6 : 
    6x6 transfer matrix.

Returns
-------
complex_taylor : ComplexTaylorStruct
    first order complex_taylor map.
)""");
  py::class_<Bmad::MatSympDecouple, std::unique_ptr<Bmad::MatSympDecouple>>(
      m, "MatSympDecouple", "mat_symp_decouple return type")
      .def_readonly("stat", &Bmad::MatSympDecouple::stat)
      .def_readonly("twiss1", &Bmad::MatSympDecouple::twiss1)
      .def_readonly("twiss2", &Bmad::MatSympDecouple::twiss2)
      .def_readonly("gamma", &Bmad::MatSympDecouple::gamma)
      .def("__len__", [](const Bmad::MatSympDecouple&) { return 4; })
      .def(
          "__getitem__",
          [](const Bmad::MatSympDecouple& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.stat);
            if (i == 1)
              return py::cast(s.twiss1);
            if (i == 2)
              return py::cast(s.twiss2);
            if (i == 3)
              return py::cast(s.gamma);
            throw py::index_error();
          });
  m.def(
      "mat_symp_decouple",
      &Bmad::mat_symp_decouple,
      py::arg("t0"),
      py::arg("U"),
      py::arg("V"),
      py::arg("Ubar"),
      py::arg("Vbar"),
      py::arg("G"),
      py::arg("type_out"),
      R"""(Parameters
----------
t0 : float
    Input matrix
stat : int
    status of results: ok$, in_stop_band$, or unstable$
U : 
V : 
Ubar : 
Vbar : 
G : 
twiss1 : TwissStruct
    Twiss params for the "upper left" mode. .phi     -- Rotation angle in radians, 0 < .PHI < twopi
twiss2 : TwissStruct
    Twiss params for the "lower right" mode. u(4,4), v(4,4), ubar(4,4), vbar(4,4), g(4,4) -- Real(rp): See MGB
    CBN 85-2 and PPB/DLR PAC89 papers for more info.
gamma : float
    gamma_c factor.
type_out : bool
    If .true. then an error message is typed out for a non ok$ STAT
)""");
  py::class_<Bmad::MatchEleToMat6, std::unique_ptr<Bmad::MatchEleToMat6>>(
      m, "MatchEleToMat6", "match_ele_to_mat6 return type")
      .def_readonly("mat6", &Bmad::MatchEleToMat6::mat6)
      .def_readonly("vec0", &Bmad::MatchEleToMat6::vec0)
      .def_readonly("err_flag", &Bmad::MatchEleToMat6::err_flag)
      .def("__len__", [](const Bmad::MatchEleToMat6&) { return 3; })
      .def(
          "__getitem__",
          [](const Bmad::MatchEleToMat6& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.mat6);
            if (i == 1)
              return py::cast(s.vec0);
            if (i == 2)
              return py::cast(s.err_flag);
            throw py::index_error();
          });
  m.def(
      "match_ele_to_mat6",
      &Bmad::match_ele_to_mat6,
      py::arg("ele"),
      py::arg("start_orb"),
      py::arg("include_delta_time") = py::none(),
      py::arg("set_trombone") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Match element.
start_orb : CoordStruct
    Starting orbit.
mat6 : float
    Transfer matrix (1st order part of xfer map).
vec0 : float
    0th order part of the transfer map.
err_flag : bool
    Set true if there is an error. False otherwise.
include_delta_time : bool, optional
    If False, ignore any finite ele.value(delta_time$). Default is True.
set_trombone : bool, optional
    Default is False. If True, set the beginning and ending Twiss values in the element to create a phase
    trombone.
)""");
  m.def(
      "mexp",
      &Bmad::mexp,
      py::arg("x"),
      py::arg("m"),
      py::arg("this_exp"),
      R"""(Parameters
----------
x : float
    Number.
m : int
    Exponent.
this_exp : 
)""");
  m.def(
      "mfft1",
      &Bmad::mfft1,
      py::arg("a"),
      py::arg("b"),
      py::arg("n"),
      py::arg("ndim"),
      py::arg("isn"),
      R"""(Parameters
----------
a : 
b : 
n : 
ndim : 
isn : 
ierr : 
)""");
  m.def(
      "misalign_ptc_fibre",
      &Bmad::misalign_ptc_fibre,
      py::arg("ele"),
      py::arg("use_offsets"),
      py::arg("for_layout"),
      R"""(Subroutine misalign_ptc_fibre (ele, use_offsets, ptc_fibre, for_layout)

Routine to misalign a fibre associated with a Bmad element.

Parameters
----------
ele : EleStruct
    Bmad element with misalignments.
use_offsets : bool
    Does ptc_fibre include element offsets, pitches and tilt? This argument is ignored if the element is a
    patch.
for_layout : bool
    If True then fibre is being created as part of a layout as opposed to a stand-alone fibre

Returns
-------
ptc_fibre : unknown
    PTC fibre element with misalignments.
)""");
  m.def(
      "momentum_compaction",
      &Bmad::momentum_compaction,
      py::arg("branch"),
      py::arg("mom_comp"),
      R"""(Parameters
----------
branch : BranchStruct
    Lattice branch to calculate on.
mom_comp : 
)""");
  py::class_<
      Bmad::MultiTurnTrackingAnalysis,
      std::unique_ptr<Bmad::MultiTurnTrackingAnalysis>>(
      m,
      "MultiTurnTrackingAnalysis",
      "multi_turn_tracking_analysis return type")
      .def_readonly("track0", &Bmad::MultiTurnTrackingAnalysis::track0)
      .def_readonly("ele", &Bmad::MultiTurnTrackingAnalysis::ele)
      .def_readonly("stable", &Bmad::MultiTurnTrackingAnalysis::stable)
      .def_readonly(
          "growth_rate", &Bmad::MultiTurnTrackingAnalysis::growth_rate)
      .def_readonly("chi", &Bmad::MultiTurnTrackingAnalysis::chi)
      .def_readonly("err_flag", &Bmad::MultiTurnTrackingAnalysis::err_flag)
      .def("__len__", [](const Bmad::MultiTurnTrackingAnalysis&) { return 6; })
      .def(
          "__getitem__",
          [](const Bmad::MultiTurnTrackingAnalysis& s, int i) -> py::object {
            if (i < 0)
              i += 6;
            if (i == 0)
              return py::cast(s.track0);
            if (i == 1)
              return py::cast(s.ele);
            if (i == 2)
              return py::cast(s.stable);
            if (i == 3)
              return py::cast(s.growth_rate);
            if (i == 4)
              return py::cast(s.chi);
            if (i == 5)
              return py::cast(s.err_flag);
            throw py::index_error();
          });
  m.def(
      "multi_turn_tracking_analysis",
      &Bmad::multi_turn_tracking_analysis,
      py::arg("track"),
      py::arg("i_dim"),
      R"""(Parameters
----------
track : CooordStruct
    multi-turn tracking data to analyze. track(i) is the particle position at a given point in the lat on the
    i^th turn.
i_dim : int
    number of dimensions used in the tracking: 2, or 4.
track0 : CoordStruct
    Closed orbit.
ele : EleStruct
    structure holding the 1-turn matrix and Twiss parameters. .mat6        -- Symplectified 1-turn matrix. If
    you want the true non-symplectified 1-turn matrix use the routine multi_turn_tracking_to_mat. .a.beta, etc
    -- a-mode beta,  etc. .a.phi       -- a-mode fractional tune in radians. .a.sigma     -- a-mode amplitude
    = sqrt(ele.a.beta * ele.a.sigma) .c_mat       -- c coupling matrix (only with i_dim = 4)
stable : bool
    Is motion stable?
growth_rate : float
    Unstable growth rate (= 0 if stable).
chi : float
    How symplectic the computed 1-turn matrix is. See mat_symp_check for more details.
err_flag : bool
    Set true if there is an error. False otherwise.
)""");
  m.def(
      "multilayer_type_to_multilayer_params",
      &Bmad::multilayer_type_to_multilayer_params,
      py::arg("ele"),
      R"""(Subroutine multilayer_type_to_multilayer_params (ele, err_flag)

Routine to set the multilayer parameters based upon the multilayer type.

Multilayer types are of the form:
  "AAA:BBB"
Where "AAA" is the atomic formula for the top layer crystal and "BBB" is the second layer atomic formula.

Parameters
----------
ele : EleStruct
    Multilayer element. .component_name -- Character: Multilayer type name. Assumed upper case. A blank name
    is not an error and results in nothing set. .value(e_tot$)  -- Photon energy in eV.

Returns
-------
err_flag : bool
    Set True if multilayer type is unrecognized. False otherwise.
)""");
  m.def(
      "multipass_chain",
      &Bmad::multipass_chain,
      py::arg("ele"),
      py::arg("ix_pass"),
      py::arg("n_links"),
      py::arg("chain_ele") = py::none(),
      py::arg("use_super_lord") = py::none(),
      R"""(Parameters
----------
ele : ElePointerStruct
    Element in a multipass chain.
ix_pass : int
    Multipass pass number of the input element. Set to -1 if input element is not in a multipass section.
n_links : int
    Number of times the physical element is passed through.
chain_ele : ElePointerStruct, optional
    pointers to the elements of the chain. Note: chain_ele(ix_pass).ele => ele
use_super_lord : bool, optional
    If present and True and if ele is a super_slave, construct the chain_ele(:) array using the corresponding
    super_lords. Output
)""");
  py::class_<Bmad::Multipole1AbToKt, std::unique_ptr<Bmad::Multipole1AbToKt>>(
      m, "Multipole1AbToKt", "multipole1_ab_to_kt return type")
      .def_readonly("knl", &Bmad::Multipole1AbToKt::knl)
      .def_readonly("tn", &Bmad::Multipole1AbToKt::tn)
      .def("__len__", [](const Bmad::Multipole1AbToKt&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::Multipole1AbToKt& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.knl);
            if (i == 1)
              return py::cast(s.tn);
            throw py::index_error();
          });
  m.def(
      "multipole1_ab_to_kt",
      &Bmad::multipole1_ab_to_kt,
      py::arg("an"),
      py::arg("bn"),
      py::arg("n"),
      R"""(Parameters
----------
an : float
    Skew multipole component.
bn : float
    Normal multipole component.
n : int
    Order of multipole.
knl : float
    Multitude magnatude.
tn : float
    Multipole angle.
)""");
  py::class_<Bmad::Multipole1KtToAb, std::unique_ptr<Bmad::Multipole1KtToAb>>(
      m, "Multipole1KtToAb", "multipole1_kt_to_ab return type")
      .def_readonly("an", &Bmad::Multipole1KtToAb::an)
      .def_readonly("bn", &Bmad::Multipole1KtToAb::bn)
      .def("__len__", [](const Bmad::Multipole1KtToAb&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::Multipole1KtToAb& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.an);
            if (i == 1)
              return py::cast(s.bn);
            throw py::index_error();
          });
  m.def(
      "multipole1_kt_to_ab",
      &Bmad::multipole1_kt_to_ab,
      py::arg("knl"),
      py::arg("knsl"),
      py::arg("tn"),
      py::arg("n"),
      R"""(Parameters
----------
knl : float
    Normal multitude component.
knsl : float
    Skew multitude component.
tn : float
    Multipole angle.
n : int
    Multipole order.
an : float
    Skew multipole component.
bn : float
    Normal multipole component.
)""");
  m.def(
      "multipole_ab_to_kt",
      &Bmad::multipole_ab_to_kt,
      py::arg("an"),
      py::arg("bn"),
      py::arg("knl"),
      py::arg("tn"),
      R"""(Parameters
----------
an : float
    Skew multipole component.
bn : float
    Normal multipole component.
knl : float
    Multitude magnatude.
tn : float
    Multipole angle.
)""");
  py::class_<Bmad::MultipoleEleToAb, std::unique_ptr<Bmad::MultipoleEleToAb>>(
      m, "MultipoleEleToAb", "multipole_ele_to_ab return type")
      .def_readonly("ix_pole_max", &Bmad::MultipoleEleToAb::ix_pole_max)
      .def_readonly("a", &Bmad::MultipoleEleToAb::a)
      .def_readonly("b", &Bmad::MultipoleEleToAb::b)
      .def_readonly("b1", &Bmad::MultipoleEleToAb::b1)
      .def("__len__", [](const Bmad::MultipoleEleToAb&) { return 4; })
      .def(
          "__getitem__",
          [](const Bmad::MultipoleEleToAb& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.ix_pole_max);
            if (i == 1)
              return py::cast(s.a);
            if (i == 2)
              return py::cast(s.b);
            if (i == 3)
              return py::cast(s.b1);
            throw py::index_error();
          });
  m.def(
      "multipole_ele_to_ab",
      &Bmad::multipole_ele_to_ab,
      py::arg("ele"),
      py::arg("use_ele_tilt"),
      py::arg("pole_type") = py::none(),
      py::arg("include_kicks") = py::none(),
      py::arg("original") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Element. .value()      -- ab_multipole values.
use_ele_tilt : bool
    If True then include ele.value(tilt_tot$) in calculations. use_ele_tilt is ignored in the case of
    multipole$ elements.
ix_pole_max : int
    Index of largest nonzero a(:) or b(:) pole. Set to -1 if all multipoles are zero. ix_pole_max is set
    independent of a nonzero b1 (if present).
a : float
    Array of multipole values.
b : float
    Array of multipole values.
pole_type : int, optional
    Type of multipole. magnetic$ (default) or electric$.
include_kicks : int, optional
    Ignored for for pole_type == electric$ for non-elseparator elements. Possibilities are: no$
    -- Default. Do not include any kick components in a and b multipoles. include_kicks$           -- Include
    hkick/vkick/dg in the n = 0 components. Also included are quad k1, sextupole k2 and octupole k3
    components.
b1 : float
    If present, b1 is set to the value of the b(1) component of the b(:) array and b(1) is set to zero. Also
    ix_pole_max is ajusted as needed. This is used by routines that want to handle b(1) in a special way in
    tracking.
original : bool, optional
    Default is false. If True, no scaling is applied.
)""");
  m.def(
      "multipole_ele_to_kt",
      &Bmad::multipole_ele_to_kt,
      py::arg("ele"),
      py::arg("use_ele_tilt"),
      py::arg("knl"),
      py::arg("tilt"),
      py::arg("pole_type") = py::none(),
      py::arg("include_kicks") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Lattice element.
use_ele_tilt : bool
    If True then include ele.value(tilt_tot$) in calculations. use_ele_tilt is ignored in the case of
    multipole$ elements.
ix_pole_max : int
    Index of largest nonzero pole.
knl : float
    Vector of strengths, MAD units.
tilt : float
    Vector of tilts.
pole_type : int, optional
    Type of multipole. magnetic$ (default) or electric$.
include_kicks : int, optional
    Possibilities are: no$                      -- Default. Do not include any kick components in a and b
    multipoles. include_kicks$           -- Include hkick/vkick/dg in the n = 0 components. Also included are
    quad k1, sextupole k2 and octupole k3 components.
)""");
  m.def(
      "multipole_init",
      &Bmad::multipole_init,
      py::arg("who"),
      py::arg("zero") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Element holding the multipoles. .a_pole(0:n_pole_maxx) -- Multipole An array .b_pole(0:n_pole_maxx) --
    Multipole Bn array
who : int
    electric$, magnetic$, or all$
zero : bool, optional
    If present and True then zero the arrays even if they already exist when this routine is called. Default
    is False which means that if the arrays already exist then this routine will do nothing.
)""");
  m.def(
      "multipole_kick",
      &Bmad::multipole_kick,
      py::arg("knl"),
      py::arg("tilt"),
      py::arg("n"),
      py::arg("ref_species"),
      py::arg("ele_orientation"),
      py::arg("coord"),
      py::arg("pole_type") = py::none(),
      py::arg("ref_orb_offset") = py::none(),
      R"""(Subroutine multipole_kick (knl, tilt, n, ref_species, ele_orientation, coord, pole_type, ref_orb_offset)

Subroutine to put in the kick due to a multipole.
Note: The kick for an electric multipole does not include any energy change.

Parameters
----------
knl : float
    Multipole integrated strength.
tilt : float
    Multipole tilt.
n : float
    Multipole order.
ref_species : int
    Reference species.
ele_orientation : int
    Element orientation +1 = normal, -1 = reversed.
coord : CoordStruct
    Particle position and direction of travel.
pole_type : int, optional
    Type of multipole. magnetic$ (default) or electric$.
ref_orb_offset : bool, optional
    If True and n = 0 then use the MAD convention and model the multipole as a zero length bend with bending
    angle knl. Default is False.
)""");
  m.def(
      "multipole_kick_mat",
      &Bmad::multipole_kick_mat,
      py::arg("knl"),
      py::arg("tilt"),
      py::arg("ref_species"),
      py::arg("ele"),
      py::arg("orbit"),
      py::arg("factor"),
      R"""(Parameters
----------
knl : float
    Strength of multipoles
tilt : float
    Tilt of multipoles
ref_species : int
    Reference species.
ele : EleStruct
    Lattice element containing multipoles.
orbit : CoordStruct
    coordinates of particle around which the multipole kick matrix is computed.
factor : float
    Factor to scale knl by.
mat6 : float
    matrix with kick values at mat6(2:4:2, 1:3:2). The rest of the matrix is untouched.
)""");
  m.def(
      "multipole_kicks",
      &Bmad::multipole_kicks,
      py::arg("knl"),
      py::arg("tilt"),
      py::arg("ele"),
      py::arg("orbit"),
      py::arg("pole_type") = py::none(),
      py::arg("ref_orb_offset") = py::none(),
      R"""(Subroutine multipole_kicks (knl, tilt, ele, orbit, pole_type, ref_orb_offset)

Subroutine to put in the kick due to a multipole element.
Also see the ab_multipole_kicks routine.

Parameters
----------
knl : float
    Multipole strengths.
tilt : float
    Multipole tilts.
ele : EleStruct
    Lattice element containing the multipoles.
orbit : CoordStruct
    Particle position.
    This parameter is an input/output and is modified in-place. As an output: Kicked particle.
pole_type : int, optional
    Type of multipole. magnetic$ (default) or electric$.
ref_orb_offset : bool, optional
    If present and n = 0 then the multipole simulates a zero length bend with bending angle knl.
)""");
  m.def(
      "multipole_kt_to_ab",
      &Bmad::multipole_kt_to_ab,
      py::arg("knl"),
      py::arg("knsl"),
      py::arg("tn"),
      py::arg("an"),
      py::arg("bn"),
      R"""(Parameters
----------
knl : float
    Normal multitude component.
knsl : float
    Skew multitude component.
tn : float
    Multipole angle.
an : float
    Skew multipole component.
bn : float
    Normal multipole component.
)""");
  m.def(
      "multipole_spin_tracking",
      &Bmad::multipole_spin_tracking,
      py::arg("ele"),
      py::arg("param"),
      py::arg("orbit"),
      R"""(Parameters
----------
ele : EleStruct
    Element
param : 
    Lat_param_struct
orbit : CoordStruct
    Particle coordinates.
)""");
  m.def(
      "mytan",
      &Bmad::mytan,
      py::arg("y"),
      py::arg("x"),
      py::arg("arg"),
      R"""(Parameters
----------
y : 
x : 
arg : 
)""");
}
