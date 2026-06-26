#include "pybmad/generated/bsim_routines_b.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_bsim_routines_b(nb::module_ &m) {
  m.def(
      "bbu_add_a_bunch",
      &bsim::bbu_add_a_bunch,
      nb::arg("lat"),
      nb::arg("bbu_beam"),
      nb::arg("bbu_param"),
      nb::arg("beam_init"),
      R"""(Wrapper for Fortran routine bbu_add_a_bunch

Parameters
----------
lat : LatStruct

bbu_beam : BbuBeamStruct

bbu_param : BbuParamStruct

beam_init : BeamInitStruct
)"""
  );
  m.def(
      "bbu_hom_voltage_calc",
      &bsim::bbu_hom_voltage_calc,
      nb::arg("lat"),
      nb::arg("bbu_beam"),
      nb::arg("n_period"),
      nb::arg("ix_stage_last_tracked"),
      R"""(Wrapper for Fortran routine bbu_hom_voltage_calc

Parameters
----------
lat : LatStruct

bbu_beam : BbuBeamStruct

n_period : int

ix_stage_last_tracked : int
)"""
  );
  m.def(
      "bbu_remove_head_bunch",
      &bsim::bbu_remove_head_bunch,
      nb::arg("bbu_beam"),
      R"""(Wrapper for Fortran routine bbu_remove_head_bunch

Parameters
----------
bbu_beam : BbuBeamStruct
)"""
  );
  m.def(
      "bbu_setup",
      &bsim::bbu_setup,
      nb::arg("lat"),
      nb::arg("dt_bunch"),
      nb::arg("bbu_param"),
      nb::arg("bbu_beam"),
      R"""(Wrapper for Fortran routine bbu_setup

Parameters
----------
lat : LatStruct

dt_bunch : float

bbu_param : BbuParamStruct

bbu_beam : BbuBeamStruct
)"""
  );
  m.def(
      "bbu_track_a_stage",
      &bsim::bbu_track_a_stage,
      nb::arg("lat"),
      nb::arg("bbu_beam"),
      nb::arg("bbu_param"),
      nb::arg("lost"),
      nb::arg("ix_stage_tracked"),
      R"""(Wrapper for Fortran routine bbu_track_a_stage

Parameters
----------
lat : LatStruct

bbu_beam : BbuBeamStruct

bbu_param : BbuParamStruct

lost : bool

ix_stage_tracked : int
)"""
  );
  nb::class_<bsim::BbuTrackAll>(m, "BbuTrackAll", "bbu_track_all return type")
      .def_ro("hom_voltage_normalized", &bsim::BbuTrackAll::hom_voltage_normalized)
      .def_ro("growth_rate", &bsim::BbuTrackAll::growth_rate)
      .def_ro("lost", &bsim::BbuTrackAll::lost)
      .def_ro("irep", &bsim::BbuTrackAll::irep)
      .def("__len__", [](const bsim::BbuTrackAll &) { return 4; })
      .def("__getitem__", [](const bsim::BbuTrackAll &s, int i) -> nb::object {
        if (i < 0)
          i += 4;
        if (i == 0)
          return nb::cast(s.hom_voltage_normalized);
        if (i == 1)
          return nb::cast(s.growth_rate);
        if (i == 2)
          return nb::cast(s.lost);
        if (i == 3)
          return nb::cast(s.irep);
        throw nb::index_error();
      });
  m.def(
      "bbu_track_all",
      &bsim::bbu_track_all,
      nb::arg("lat"),
      nb::arg("bbu_beam"),
      nb::arg("bbu_param"),
      nb::arg("beam_init"),
      R"""(Wrapper for Fortran routine bbu_track_all

Parameters
----------
lat : LatStruct

bbu_beam : BbuBeamStruct

bbu_param : BbuParamStruct

beam_init : BeamInitStruct

Returns
-------
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
