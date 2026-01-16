#include "pybmad/generated/Bmad_routines_n.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_Bmad_routines_n(py::module &m) {
  m.def(
      "n_attrib_string_max_len",
      &Bmad::n_attrib_string_max_len,
      R"""(Function n_attrib_string_max_len () result (max_len)

Routine to return the the maximum number of characters in any attribute
name known to bmad.


Returns
-------
max_len : int
    Maximum number of characters in any attribute name.
)"""
  );
  m.def(
      "new_control",
      &Bmad::new_control,
      py::arg("lat"),
      py::arg("ix_ele"),
      py::arg("ele_name") = py::none(),
      R"""(Parameters
----------
lat : LatStruct
    Lat used
ix_ele : int
    Index of the new control element
ele_name : unknown, optional
    Name of the new element. Output
)"""
  );
  m.def(
      "nint_chk",
      &Bmad::nint_chk,
      py::arg("re_val"),
      R"""(Function nint_chk (re_val) result (int_val)

Returns the nearest integer to re_val.
Also does out-of-bounds error checking.
Used with bmad parsing.

Parameters
----------
re_val : float
    Input real number.

Returns
-------
int_val : int
    Output nearest integer.
)"""
  );
  m.def(
      "normal_form_complex_taylors",
      &Bmad::normal_form_complex_taylors,
      py::arg("one_turn_taylor"),
      py::arg("rf_on"),
      py::arg("F") = py::none(),
      py::arg("L") = py::none(),
      py::arg("A") = py::none(),
      py::arg("A_inverse") = py::none(),
      py::arg("order") = py::none(),
      R"""(Parameters
----------
one_turn_taylor : 
rf_on : 
F : 
L : 
A : 
A_inverse : 
order : 
)"""
  );
  py::class_<Bmad::NormalFormTaylors, std::unique_ptr<Bmad::NormalFormTaylors>>(
      m,
      "NormalFormTaylors",
      "normal_form_taylors return type"
  )
      .def_readonly("dhdj", &Bmad::NormalFormTaylors::dhdj)
      .def_readonly("A", &Bmad::NormalFormTaylors::A)
      .def_readonly("A_inverse", &Bmad::NormalFormTaylors::A_inverse)
      .def("__len__", [](const Bmad::NormalFormTaylors &) { return 3; })
      .def("__getitem__", [](const Bmad::NormalFormTaylors &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.dhdj);
        if (i == 1)
          return py::cast(s.A);
        if (i == 2)
          return py::cast(s.A_inverse);
        throw py::index_error();
      });
  m.def(
      "normal_form_taylors",
      &Bmad::normal_form_taylors,
      py::arg("one_turn_taylor"),
      py::arg("rf_on"),
      R"""(Subroutine normal_form_taylors(one_turn_taylor, rf_on, dhdj, A, A_inverse)

Do a normal form decomposition on a one-turn taylor map M:
  M = A o R o A_inverse
where A maps Floquet (fully normalized) coordinates to lab coordinates.
In Floquet coordinates, the amplitudes are defined as J_i = (1/2) (x_i^2 + p_i^2).
The map R = exp(:h:) is a pure rotation with h = h(J) is a function of the amplitudes only.
The angles (phase advances) are given by phi_i = 2pi*dh/dJ_i.
The taylor terms of dhdj are therefore the tunes, chromaticities, amplitude dependent tune shifts, etc.

The mapping procedure for one turn is:
 z_Floquet_in = A_inverse o z_Lab_in
 [phi_a, phi_b, phi_c] = 2 pi * dhdj o z_Floquet_in
 z_Floquet_out = RotationMatrix(phi_a, phi_b, phi_c) . z_Floquet_in
 z_Lab_out = A o z_Floquet_out


Parameters
----------
one_turn_taylor : TaylorStruct
    one turn taylor map
rf_on : bool
    Was the map calculated with RF on?

Returns
-------
A : TaylorStruct
    Map from Floquet coordinates to Lab coordinates
A_inverse : TaylorStruct
    Map from Lab coordinates to Floquet coordinates
dhdj : TaylorStruct
    Map from Floquet coordinates to phase advances
)"""
  );
  py::class_<Bmad::NormalMode3Calc, std::unique_ptr<Bmad::NormalMode3Calc>>(
      m,
      "NormalMode3Calc",
      "normal_mode3_calc return type"
  )
      .def_readonly("tune", &Bmad::NormalMode3Calc::tune)
      .def_readonly("B", &Bmad::NormalMode3Calc::B)
      .def_readonly("HV", &Bmad::NormalMode3Calc::HV)
      .def("__len__", [](const Bmad::NormalMode3Calc &) { return 3; })
      .def("__getitem__", [](const Bmad::NormalMode3Calc &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.tune);
        if (i == 1)
          return py::cast(s.B);
        if (i == 2)
          return py::cast(s.HV);
        throw py::index_error();
      });
  m.def(
      "normal_mode3_calc",
      &Bmad::normal_mode3_calc,
      py::arg("t6"),
      py::arg("above_transition") = py::none(),
      py::arg("abz_tunes") = py::none(),
      R"""(Subroutine normal_mode3_calc (mat, tune, B, HV, above_transition)

Does an Eigen decomposition of the 1-turn transfer matrix (mat) and generates
B, V, H.

If the above_transition argument is present and false, then the 3rd (z) mode is assumed
to have a positive slip factor (z-mode rotates counter clockwise in phase space).
Default is True ==> z-mode has a negative slip factor so the mode rotates clock-wise in phase space.

Parameters
----------
mat : float
    1-turn transfer matrix
above_transition : bool, optional
    If present and false, then z-mode assumes positive slip factor. Else negative slip factor assumed.
abz_tunes : float, optional
    Tunes to order eigensystem by.

Returns
-------
tune : float
    Tunes of the 3 normal modes (radians)
B : float
    B is block diagonal and related to the normal mode Twiss parameters.
HV : float
    Transforms from normal mode coordinates to canonical coordinates: x = H.V.a
)"""
  );
  m.def(
      "normal_mode_dispersion",
      &Bmad::normal_mode_dispersion,
      py::arg("ele"),
      py::arg("reverse") = py::none(),
      R"""(Parameters
----------
ele : EleStruct
    Element whose dispersions are to be adjusted.
    This parameter is an input/output and is modified in-place. As an output: Element with adjusted
    dispersions.
reverse : bool, optional
    Default is False. If True, calculate the x,y dispersions from the normal mode ones.
)"""
  );
  m.def(
      "normalize_evecs",
      &Bmad::normalize_evecs,
      py::arg("evec"),
      R"""(Subroutine normalize_evecs(evec, err_flag)

Normalizes eigenvectors such that transpose(E).S.E = iS, where E = evec_r + i evec_i

Parameters
----------
evec : float
    complex eigenvectors arranged down columns.
    This parameter is an input/output and is modified in-place. As an output: Eigensystem normalized to be
    symplectic.

Returns
-------
err_flag : bool
    Set true of normalization is not possible due to amplitude is zero.
)"""
  );
  m.def(
      "num_field_eles",
      &Bmad::num_field_eles,
      py::arg("ele"),
      R"""(Parameters
----------
ele : EleStruct
    Element with sum number of associated field elements.
n_field_ele : int
    Number of associated field elements.
)"""
  );
  m.def(
      "num_lords",
      &Bmad::num_lords,
      py::arg("slave"),
      py::arg("lord_type"),
      R"""(Parameters
----------
slave : EleStruct
    Slave element.
lord_type : int
    Type of lord. super_lord$, multipass_lord$, girder_lord$, group_lord$, overlay_lord$, and governor$ (=
    group + overlay + control + girder)
num : int
    Number of lords of the given type.
)"""
  );
}
