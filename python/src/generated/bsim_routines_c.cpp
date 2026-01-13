#include "pybmad/generated/bsim_routines_c.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyCheckRfFreq python_check_rf_freq(LatProxy& lat, double fb) {
  bsim::check_rf_freq(lat, fb);
  auto py_result{PyCheckRfFreq{fb}};
  return py_result;
}

void init_bsim_routines_c(py::module& m) {
  py::class_<PyCheckRfFreq, std::unique_ptr<PyCheckRfFreq>>(
      m, "CheckRfFreq", "check_rf_freq return type")
      .def_readonly("fb", &PyCheckRfFreq::fb)
      .def("__len__", [](const PyCheckRfFreq&) { return 1; })
      .def("__getitem__", [](const PyCheckRfFreq& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.fb);
        throw py::index_error();
      });
  m.def(
      "check_rf_freq",
      &python_check_rf_freq,
      py::arg("lat"),
      py::arg("fb"),
      R"""(Parameters
  ----------
  lat : 
  fb : 
  )""");
  m.def(
      "count_lines_in_file",
      &bsim::count_lines_in_file,
      py::arg("file_name"),
      R"""(Parameters
  ----------
  file_name : 
  lines : 
  )""");
}
