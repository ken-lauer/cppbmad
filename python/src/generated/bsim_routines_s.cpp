#include "pybmad/generated/bsim_routines_s.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PySetTune3d python_set_tune_3d(
    BranchProxy& branch,
    FixedArray1D<Real, 3> target_tunes,
    std::optional<std::string> mask,
    std::optional<bool> use_phase_trombone,
    std::optional<bool> z_tune_set,
    std::optional<FixedArray1D<string, 2>> group_knobs,
    std::optional<bool> print_err,
    bool everything_ok) {
  bsim::set_tune_3d(
      branch,
      target_tunes,
      make_opt_ref(mask),
      use_phase_trombone,
      z_tune_set,
      group_knobs,
      print_err,
      everything_ok);
  auto py_result{PySetTune3d{mask, everything_ok}};
  return py_result;
}

void init_bsim_routines_s(py::module& m) {
  py::class_<PySetTune3d, std::unique_ptr<PySetTune3d>>(
      m, "SetTune3d", "Fortran routine set_tune_3d return value")
      .def_readonly("mask", &PySetTune3d::mask)
      .def_readonly("everything_ok", &PySetTune3d::everything_ok)
      .def("__len__", [](const PySetTune3d&) { return 2; })
      .def("__getitem__", [](const PySetTune3d& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.mask);
        if (i == 1)
          return py::cast(s.everything_ok);
        throw py::index_error();
      });
  m.def(
      "set_tune_3d",
      &python_set_tune_3d,
      py::arg("branch"),
      py::arg("target_tunes"),
      py::arg("mask") = py::none(),
      py::arg("use_phase_trombone") = py::none(),
      py::arg("z_tune_set") = py::none(),
      py::arg("group_knobs") = py::none(),
      py::arg("print_err") = py::none(),
      py::arg("everything_ok"),
      R"""(Parameters
  ----------
  branch : BranchStruct
      This parameter is an input/output and is modified in-place. As an output: with adjusted quads and RF to
      match desired tunes.
  target_tunes : float
      tunes for a, b, z modes (rad/2pi). Must include integer part.
  mask : 
  use_phase_trombone : bool, optional
      Default False. If true, use a match element in phase trombone mode to adjust the tunes. The match element
      must be the first element in the lattice. Use insert_phase_trombone to insert one.
  z_tune_set : bool, optional
      Default True. If false, do not try to set the synch tune.
  group_knobs : unknown, optional
      If set non-blank, use these group elements for tuning.
  print_err : bool, optional
      Print error message if there is a problem? Default is True.
  everything_ok : 
  )""");
}
