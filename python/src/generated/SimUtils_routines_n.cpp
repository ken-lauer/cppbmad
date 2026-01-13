#include "pybmad/generated/SimUtils_routines_n.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyNBinsAutomatic python_n_bins_automatic(int n_data, int n) {
  SimUtils::n_bins_automatic(n_data, n);
  auto py_result{PyNBinsAutomatic{n_data, n}};
  return py_result;
}
PyNChooseK python_n_choose_k(int n, int k, double nck) {
  SimUtils::n_choose_k(n, k, nck);
  auto py_result{PyNChooseK{n, k, nck}};
  return py_result;
}
PyNaff python_naff(
    ComplexAlloc1D& cdata,
    RealAlloc1D& freqs,
    ComplexAlloc1D& amps,
    std::optional<int> opt_dump_spectra = std::nullopt,
    std::optional<bool> opt_zero_first = std::nullopt) {
  SimUtils::naff(
      cdata,
      freqs,
      amps,
      make_opt_ref(opt_dump_spectra),
      make_opt_ref(opt_zero_first));
  auto py_result{PyNaff{opt_dump_spectra, opt_zero_first}};
  return py_result;
}
PyNametableAdd python_nametable_add(
    NametableProxy& nametable,
    std::string name,
    int ix_name) {
  SimUtils::nametable_add(nametable, name, ix_name);
  auto py_result{PyNametableAdd{name, ix_name}};
  return py_result;
}
PyNametableBracketIndexx python_nametable_bracket_indexx(
    NametableProxy& nametable,
    std::string name,
    std::optional<int> n_match,
    int ix_max) {
  SimUtils::nametable_bracket_indexx(
      nametable, name, make_opt_ref(n_match), ix_max);
  auto py_result{PyNametableBracketIndexx{name, n_match, ix_max}};
  return py_result;
}
PyNametableChange1 python_nametable_change1(
    NametableProxy& nametable,
    std::string name,
    int ix_name) {
  SimUtils::nametable_change1(nametable, name, ix_name);
  auto py_result{PyNametableChange1{name, ix_name}};
  return py_result;
}
PyNametableInit python_nametable_init(
    NametableProxy& nametable,
    std::optional<int> n_min = std::nullopt,
    std::optional<int> n_max = std::nullopt) {
  SimUtils::nametable_init(nametable, make_opt_ref(n_min), make_opt_ref(n_max));
  auto py_result{PyNametableInit{n_min, n_max}};
  return py_result;
}
PyNametableRemove python_nametable_remove(
    NametableProxy& nametable,
    int ix_name) {
  SimUtils::nametable_remove(nametable, ix_name);
  auto py_result{PyNametableRemove{ix_name}};
  return py_result;
}

void init_SimUtils_routines_n(py::module& m) {
  py::class_<PyNBinsAutomatic, std::unique_ptr<PyNBinsAutomatic>>(
      m, "NBinsAutomatic", "n_bins_automatic return type")
      .def_readonly("n_data", &PyNBinsAutomatic::n_data)
      .def_readonly("n", &PyNBinsAutomatic::n)
      .def("__len__", [](const PyNBinsAutomatic&) { return 2; })
      .def("__getitem__", [](const PyNBinsAutomatic& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.n_data);
        if (i == 1)
          return py::cast(s.n);
        throw py::index_error();
      });
  m.def(
      "n_bins_automatic",
      &python_n_bins_automatic,
      py::arg("n_data"),
      py::arg("n"),
      R"""(Function to automatically select the number of bins

  )""");
  py::class_<PyNChooseK, std::unique_ptr<PyNChooseK>>(
      m, "NChooseK", "n_choose_k return type")
      .def_readonly("n", &PyNChooseK::n)
      .def_readonly("k", &PyNChooseK::k)
      .def_readonly("nck", &PyNChooseK::nck)
      .def("__len__", [](const PyNChooseK&) { return 3; })
      .def("__getitem__", [](const PyNChooseK& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.n);
        if (i == 1)
          return py::cast(s.k);
        if (i == 2)
          return py::cast(s.nck);
        throw py::index_error();
      });
  m.def(
      "n_choose_k",
      &python_n_choose_k,
      py::arg("n"),
      py::arg("k"),
      py::arg("nck"),
      R"""(Parameters
  ----------
  n : 
  k : 
  nck : 
  )""");
  m.def(
      "n_spline_create",
      &SimUtils::n_spline_create,
      py::arg("deriv0"),
      py::arg("deriv1"),
      py::arg("x1"),
      R"""(Parameters
  ----------
  deriv0 : float
      Derivative vector from order 0 to some order n at x = 0.
  deriv1 : float
      Derivative vector from order 0 to some order n at x = x1.
  x1 : float
      Location where deriv1 derivatives have been evaluated.
  n_spline : 
      real(rp), Derivative vector from order 0 to order 2*n+1 of the interpolation spline.
  )""");
  py::class_<PyNaff, std::unique_ptr<PyNaff>>(m, "Naff", "naff return type")
      .def_readonly("opt_dump_spectra", &PyNaff::opt_dump_spectra)
      .def_readonly("opt_zero_first", &PyNaff::opt_zero_first)
      .def("__len__", [](const PyNaff&) { return 2; })
      .def("__getitem__", [](const PyNaff& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.opt_dump_spectra);
        if (i == 1)
          return py::cast(s.opt_zero_first);
        throw py::index_error();
      });
  m.def(
      "naff",
      &python_naff,
      py::arg("cdata"),
      py::arg("freqs"),
      py::arg("amps"),
      py::arg("opt_dump_spectra") = py::none(),
      py::arg("opt_zero_first") = py::none(),
      R"""(subroutine naff(cdata,freqs,amps,opt_dump_spectra,opt_zero_first)

  This subroutine implements the NAFF algorithm for calculating the spectra
  of periodic data.

  See naff_mod documentation for details.

  Frequencies returned are in units of 2pi. That is, freqs ranges from 0 to 1.

  freqs and amps must be allocated before hand.  This subroutine will repeat the
  decomposition loop until all elements of freqs and amps are populated.

  )""");
  py::class_<PyNametableAdd, std::unique_ptr<PyNametableAdd>>(
      m, "NametableAdd", "nametable_add return type")
      .def_readonly("name", &PyNametableAdd::name)
      .def_readonly("ix_name", &PyNametableAdd::ix_name)
      .def("__len__", [](const PyNametableAdd&) { return 2; })
      .def("__getitem__", [](const PyNametableAdd& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.name);
        if (i == 1)
          return py::cast(s.ix_name);
        throw py::index_error();
      });
  m.def(
      "nametable_add",
      &python_nametable_add,
      py::arg("nametable"),
      py::arg("name"),
      py::arg("ix_name"),
      R"""(Parameters
  ----------
  nametable : 
  name : 
  ix_name : 
  )""");
  py::class_<
      PyNametableBracketIndexx,
      std::unique_ptr<PyNametableBracketIndexx>>(
      m, "NametableBracketIndexx", "nametable_bracket_indexx return type")
      .def_readonly("name", &PyNametableBracketIndexx::name)
      .def_readonly("n_match", &PyNametableBracketIndexx::n_match)
      .def_readonly("ix_max", &PyNametableBracketIndexx::ix_max)
      .def("__len__", [](const PyNametableBracketIndexx&) { return 3; })
      .def(
          "__getitem__",
          [](const PyNametableBracketIndexx& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.name);
            if (i == 1)
              return py::cast(s.n_match);
            if (i == 2)
              return py::cast(s.ix_max);
            throw py::index_error();
          });
  m.def(
      "nametable_bracket_indexx",
      &python_nametable_bracket_indexx,
      py::arg("nametable"),
      py::arg("name"),
      py::arg("n_match") = py::none(),
      py::arg("ix_max"),
      R"""(Parameters
  ----------
  nametable : 
  name : 
  n_match : 
  ix_max : 
  )""");
  py::class_<PyNametableChange1, std::unique_ptr<PyNametableChange1>>(
      m, "NametableChange1", "nametable_change1 return type")
      .def_readonly("name", &PyNametableChange1::name)
      .def_readonly("ix_name", &PyNametableChange1::ix_name)
      .def("__len__", [](const PyNametableChange1&) { return 2; })
      .def("__getitem__", [](const PyNametableChange1& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.name);
        if (i == 1)
          return py::cast(s.ix_name);
        throw py::index_error();
      });
  m.def(
      "nametable_change1",
      &python_nametable_change1,
      py::arg("nametable"),
      py::arg("name"),
      py::arg("ix_name"),
      R"""(Parameters
  ----------
  nametable : 
  name : 
  ix_name : 
  )""");
  py::class_<PyNametableInit, std::unique_ptr<PyNametableInit>>(
      m, "NametableInit", "nametable_init return type")
      .def_readonly("n_min", &PyNametableInit::n_min)
      .def_readonly("n_max", &PyNametableInit::n_max)
      .def("__len__", [](const PyNametableInit&) { return 2; })
      .def("__getitem__", [](const PyNametableInit& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.n_min);
        if (i == 1)
          return py::cast(s.n_max);
        throw py::index_error();
      });
  m.def(
      "nametable_init",
      &python_nametable_init,
      py::arg("nametable"),
      py::arg("n_min") = py::none(),
      py::arg("n_max") = py::none(),
      R"""(Parameters
  ----------
  nametable : 
  n_min : 
  n_max : 
  )""");
  py::class_<PyNametableRemove, std::unique_ptr<PyNametableRemove>>(
      m, "NametableRemove", "nametable_remove return type")
      .def_readonly("ix_name", &PyNametableRemove::ix_name)
      .def("__len__", [](const PyNametableRemove&) { return 1; })
      .def("__getitem__", [](const PyNametableRemove& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.ix_name);
        throw py::index_error();
      });
  m.def(
      "nametable_remove",
      &python_nametable_remove,
      py::arg("nametable"),
      py::arg("ix_name"),
      R"""(Parameters
  ----------
  nametable : 
  ix_name : 
  )""");
}
