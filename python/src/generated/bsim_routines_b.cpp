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
lat : LatStruct

bbu_beam : BbuBeamStruct

bbu_param : BbuParamStruct

beam_init : BeamInitStruct

Returns
-------
lat : LatStruct

bbu_beam : BbuBeamStruct

bbu_param : BbuParamStruct

beam_init : BeamInitStruct
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
lat : LatStruct

bbu_beam : BbuBeamStruct

n_period : int

ix_stage_last_tracked : int

Returns
-------
lat : LatStruct

bbu_beam : BbuBeamStruct

n_period : int

ix_stage_last_tracked : int
)"""
  );
  m.def(
      "bbu_remove_head_bunch",
      &bsim::bbu_remove_head_bunch,
      py::arg("bbu_beam"),
      R"""(Wrapper for Fortran routine bbu_remove_head_bunch

Parameters
----------
bbu_beam : BbuBeamStruct

Returns
-------
bbu_beam : BbuBeamStruct
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
lat : LatStruct

dt_bunch : float

bbu_param : BbuParamStruct

bbu_beam : BbuBeamStruct

Returns
-------
lat : LatStruct

dt_bunch : float

bbu_param : BbuParamStruct

bbu_beam : BbuBeamStruct
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
lat : LatStruct

bbu_beam : BbuBeamStruct

bbu_param : BbuParamStruct

lost : bool

ix_stage_tracked : int

Returns
-------
lat : LatStruct

bbu_beam : BbuBeamStruct

bbu_param : BbuParamStruct

lost : bool

ix_stage_tracked : int
)"""
  );
  py::class_<bsim::BbuTrackAll, std::unique_ptr<bsim::BbuTrackAll>>(
      m,
      "BbuTrackAll",
      "bbu_track_all return type"
  )
      .def_readonly("hom_voltage_normalized", &bsim::BbuTrackAll::hom_voltage_normalized)
      .def_readonly("growth_rate", &bsim::BbuTrackAll::growth_rate)
      .def_readonly("lost", &bsim::BbuTrackAll::lost)
      .def_readonly("irep", &bsim::BbuTrackAll::irep)
      .def("__len__", [](const bsim::BbuTrackAll &) { return 4; })
      .def("__getitem__", [](const bsim::BbuTrackAll &s, int i) -> py::object {
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
      &bsim::bbu_track_all,
      py::arg("lat"),
      py::arg("bbu_beam"),
      py::arg("bbu_param"),
      py::arg("beam_init"),
      R"""(Wrapper for Fortran routine bbu_track_all

Parameters
----------
lat : LatStruct

bbu_beam : BbuBeamStruct

bbu_param : BbuParamStruct

beam_init : BeamInitStruct

Returns
-------
lat : LatStruct

bbu_beam : BbuBeamStruct

bbu_param : BbuParamStruct

beam_init : BeamInitStruct

hom_voltage_normalized : float
    HOM voltage normalized

growth_rate : float
    Growth rate

lost : bool
    Lost

irep : int
    irep
)"""
  );
}
