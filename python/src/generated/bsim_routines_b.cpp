#include "pybmad/generated/bsim_routines_b.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyBbuHomVoltageCalc python_bbu_hom_voltage_calc(
    LatProxy& lat,
    BbuBeamProxy& bbu_beam,
    int n_period,
    int ix_stage_last_tracked) {
  bsim::bbu_hom_voltage_calc(lat, bbu_beam, n_period, ix_stage_last_tracked);
  auto py_result{PyBbuHomVoltageCalc{n_period, ix_stage_last_tracked}};
  return py_result;
}
PyBbuSetup python_bbu_setup(
    LatProxy& lat,
    double dt_bunch,
    BbuParamProxy& bbu_param,
    BbuBeamProxy& bbu_beam) {
  bsim::bbu_setup(lat, dt_bunch, bbu_param, bbu_beam);
  auto py_result{PyBbuSetup{dt_bunch}};
  return py_result;
}
PyBbuTrackAStage python_bbu_track_a_stage(
    LatProxy& lat,
    BbuBeamProxy& bbu_beam,
    BbuParamProxy& bbu_param,
    bool lost,
    int ix_stage_tracked) {
  bsim::bbu_track_a_stage(lat, bbu_beam, bbu_param, lost, ix_stage_tracked);
  auto py_result{PyBbuTrackAStage{lost, ix_stage_tracked}};
  return py_result;
}
PyBbuTrackAll python_bbu_track_all(
    LatProxy& lat,
    BbuBeamProxy& bbu_beam,
    BbuParamProxy& bbu_param,
    BeamInitProxy& beam_init,
    double hom_voltage_normalized,
    double growth_rate,
    bool lost,
    int irep) {
  bsim::bbu_track_all(
      lat,
      bbu_beam,
      bbu_param,
      beam_init,
      hom_voltage_normalized,
      growth_rate,
      lost,
      irep);
  auto py_result{
      PyBbuTrackAll{hom_voltage_normalized, growth_rate, lost, irep}};
  return py_result;
}

void init_bsim_routines_b(py::module& m) {
  m.def(
      "bbu_add_a_bunch",
      &bsim::bbu_add_a_bunch,
      py::arg("lat"),
      py::arg("bbu_beam"),
      py::arg("bbu_param"),
      py::arg("beam_init"),
      R"""(Parameters
  ----------
  lat : 
  bbu_beam : 
  bbu_param : 
  beam_init : 
  )""");
  py::class_<PyBbuHomVoltageCalc, std::unique_ptr<PyBbuHomVoltageCalc>>(
      m, "BbuHomVoltageCalc", "bbu_hom_voltage_calc return type")
      .def_readonly("n_period", &PyBbuHomVoltageCalc::n_period)
      .def_readonly(
          "ix_stage_last_tracked", &PyBbuHomVoltageCalc::ix_stage_last_tracked)
      .def("__len__", [](const PyBbuHomVoltageCalc&) { return 2; })
      .def(
          "__getitem__", [](const PyBbuHomVoltageCalc& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.n_period);
            if (i == 1)
              return py::cast(s.ix_stage_last_tracked);
            throw py::index_error();
          });
  m.def(
      "bbu_hom_voltage_calc",
      &python_bbu_hom_voltage_calc,
      py::arg("lat"),
      py::arg("bbu_beam"),
      py::arg("n_period"),
      py::arg("ix_stage_last_tracked"),
      R"""(Parameters
  ----------
  lat : 
  bbu_beam : 
  n_period : 
  ix_stage_last_tracked : 
  )""");
  m.def(
      "bbu_remove_head_bunch",
      &bsim::bbu_remove_head_bunch,
      py::arg("bbu_beam"),
      R"""(Parameters
  ----------
  bbu_beam : 
  )""");
  py::class_<PyBbuSetup, std::unique_ptr<PyBbuSetup>>(
      m, "BbuSetup", "bbu_setup return type")
      .def_readonly("dt_bunch", &PyBbuSetup::dt_bunch)
      .def("__len__", [](const PyBbuSetup&) { return 1; })
      .def("__getitem__", [](const PyBbuSetup& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.dt_bunch);
        throw py::index_error();
      });
  m.def(
      "bbu_setup",
      &python_bbu_setup,
      py::arg("lat"),
      py::arg("dt_bunch"),
      py::arg("bbu_param"),
      py::arg("bbu_beam"),
      R"""(Parameters
  ----------
  lat : 
  dt_bunch : 
  bbu_param : 
  bbu_beam : 
  )""");
  py::class_<PyBbuTrackAStage, std::unique_ptr<PyBbuTrackAStage>>(
      m, "BbuTrackAStage", "bbu_track_a_stage return type")
      .def_readonly("lost", &PyBbuTrackAStage::lost)
      .def_readonly("ix_stage_tracked", &PyBbuTrackAStage::ix_stage_tracked)
      .def("__len__", [](const PyBbuTrackAStage&) { return 2; })
      .def("__getitem__", [](const PyBbuTrackAStage& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.lost);
        if (i == 1)
          return py::cast(s.ix_stage_tracked);
        throw py::index_error();
      });
  m.def(
      "bbu_track_a_stage",
      &python_bbu_track_a_stage,
      py::arg("lat"),
      py::arg("bbu_beam"),
      py::arg("bbu_param"),
      py::arg("lost"),
      py::arg("ix_stage_tracked"),
      R"""(Parameters
  ----------
  lat : 
  bbu_beam : 
  bbu_param : 
  lost : 
  ix_stage_tracked : 
  )""");
  py::class_<PyBbuTrackAll, std::unique_ptr<PyBbuTrackAll>>(
      m, "BbuTrackAll", "bbu_track_all return type")
      .def_readonly(
          "hom_voltage_normalized", &PyBbuTrackAll::hom_voltage_normalized)
      .def_readonly("growth_rate", &PyBbuTrackAll::growth_rate)
      .def_readonly("lost", &PyBbuTrackAll::lost)
      .def_readonly("irep", &PyBbuTrackAll::irep)
      .def("__len__", [](const PyBbuTrackAll&) { return 4; })
      .def("__getitem__", [](const PyBbuTrackAll& s, int i) -> py::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return py::cast(s.hom_voltage_normalized);
        if (i == 1)
          return py::cast(s.growth_rate);
        if (i == 2)
          return py::cast(s.lost);
        if (i == 3)
          return py::cast(s.irep);
        throw py::index_error();
      });
  m.def(
      "bbu_track_all",
      &python_bbu_track_all,
      py::arg("lat"),
      py::arg("bbu_beam"),
      py::arg("bbu_param"),
      py::arg("beam_init"),
      py::arg("hom_voltage_normalized"),
      py::arg("growth_rate"),
      py::arg("lost"),
      py::arg("irep"),
      R"""(Parameters
  ----------
  lat : 
  bbu_beam : 
  bbu_param : 
  beam_init : 
  hom_voltage_normalized : 
  growth_rate : 
  lost : 
  irep : 
  )""");
}
