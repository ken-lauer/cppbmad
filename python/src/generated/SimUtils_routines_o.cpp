#include "pybmad/generated/SimUtils_routines_o.hpp"
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

// Wrappers
struct PyOrdinalStr {
  int n;
  std::string str;
};
PyOrdinalStr python_ordinal_str(int n, std::string str) {
  SimUtils::ordinal_str(n, str);
  auto py_result{PyOrdinalStr{n, str}};
  return py_result;
}

void init_SimUtils_routines_o(py::module& m) {
  m.def(
      "omega_to_quat",
      &SimUtils::omega_to_quat,
      py::arg("omega"),
      R"""(Function omega_to_quat (omega) result (quat)

  Routine to convert from omega + angle representation to a quaternion.

  Parameters
  ----------
  omega : float
      Axis of rotation + magnitude = rotation angle.

  Returns
  -------
  quat : float
      Rotation quaternion.
  )""");
  m.def(
      "openpmd_species_name",
      &SimUtils::openpmd_species_name,
      py::arg("species"),
      R"""(Function openpmd_species_name (species) result(pmd_name)

  Routine to return the openPMD name of a particle species given the Bmad species ID.
  Note: the pmd_name does not include the particle charge. For example, if species
  corresponds to He+ then the pmd_name will be "He".

  Parameters
  ----------
  species : int
      Bmad species ID number.

  Returns
  -------
  pmd_name : unknown
      Name of the species. Will return 'INVALID!' (= invalid_name) if index is not valid.
  )""");
  m.def(
      "ordinal_str",
      &python_ordinal_str,
      py::arg("n"),
      py::arg("str"),
      R"""(Parameters
  ----------
  n : 
  str : 
  )""");
  py::class_<PyOrdinalStr, std::unique_ptr<PyOrdinalStr>>(
      m, "OrdinalStr", "Fortran routine ordinal_str return value")
      .def_readonly("n", &PyOrdinalStr::n)
      .def_readonly("str", &PyOrdinalStr::str)
      .def("__len__", [](const PyOrdinalStr&) { return 2; })
      .def("__getitem__", [](const PyOrdinalStr& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.n);
        if (i == 1)
          return py::cast(s.str);
        throw py::index_error();
      });
}
