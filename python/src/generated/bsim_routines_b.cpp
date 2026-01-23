#include "pybmad/generated/bsim_routines_b.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_bsim_routines_b(py::module &m) {
  m.def(
      "bbu_add_a_bunch",
      &bsim::bbu_add_a_bunch,
      py::arg("lat"),
      py::arg("bbu_beam"),
      py::arg("bbu_param"),
      py::arg("beam_init"),
      R"""(Wrapper for Fortran routine bbu_add_a_bunch

Parameters
----------
lat : 
bbu_beam : 
bbu_param : 
beam_init : 
)"""
  );
  m.def(
      "bbu_hom_voltage_calc",
      &bsim::bbu_hom_voltage_calc,
      py::arg("lat"),
      py::arg("bbu_beam"),
      py::arg("n_period"),
      py::arg("ix_stage_last_tracked"),
      R"""(Wrapper for Fortran routine bbu_hom_voltage_calc

Parameters
----------
lat : 
bbu_beam : 
n_period : 
ix_stage_last_tracked : 
)"""
  );
  m.def(
      "bbu_remove_head_bunch",
      &bsim::bbu_remove_head_bunch,
      py::arg("bbu_beam"),
      R"""(Wrapper for Fortran routine bbu_remove_head_bunch

Parameters
----------
bbu_beam : 
)"""
  );
  m.def(
      "bbu_setup",
      &bsim::bbu_setup,
      py::arg("lat"),
      py::arg("dt_bunch"),
      py::arg("bbu_param"),
      py::arg("bbu_beam"),
      R"""(Wrapper for Fortran routine bbu_setup

Parameters
----------
lat : 
dt_bunch : 
bbu_param : 
bbu_beam : 
)"""
  );
  m.def(
      "bbu_track_a_stage",
      &bsim::bbu_track_a_stage,
      py::arg("lat"),
      py::arg("bbu_beam"),
      py::arg("bbu_param"),
      py::arg("lost"),
      py::arg("ix_stage_tracked"),
      R"""(Wrapper for Fortran routine bbu_track_a_stage

Parameters
----------
lat : 
bbu_beam : 
bbu_param : 
lost : 
ix_stage_tracked : 
)"""
  );
  m.def(
      "bbu_track_all",
      &bsim::bbu_track_all,
      py::arg("lat"),
      py::arg("bbu_beam"),
      py::arg("bbu_param"),
      py::arg("beam_init"),
      py::arg("hom_voltage_normalized"),
      py::arg("growth_rate"),
      py::arg("lost"),
      py::arg("irep"),
      R"""(Wrapper for Fortran routine bbu_track_all

Parameters
----------
lat : 
bbu_beam : 
bbu_param : 
beam_init : 
hom_voltage_normalized : 
growth_rate : 
lost : 
irep : 
)"""
  );
}
