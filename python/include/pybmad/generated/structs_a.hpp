#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_astra_lattice_param_struct(nb::module_ &m, nb::class_<AstraLatticeParamStruct> &class_);
void init_ac_kicker_freq_struct(nb::module_ &m, nb::class_<AcKickerFreqStruct> &class_);
void init_ac_kicker_struct(nb::module_ &m, nb::class_<AcKickerStruct> &class_);
void init_ac_kicker_time_struct(nb::module_ &m, nb::class_<AcKickerTimeStruct> &class_);
void init_anormal_mode_struct(nb::module_ &m, nb::class_<AnormalModeStruct> &class_);
void init_aperture_param_struct(nb::module_ &m, nb::class_<ApertureParamStruct> &class_);
void init_aperture_point_struct(nb::module_ &m, nb::class_<AperturePointStruct> &class_);
void init_aperture_scan_struct(nb::module_ &m, nb::class_<ApertureScanStruct> &class_);
void init_all_pointer_struct(nb::module_ &m, nb::class_<AllPointerStruct> &class_);
void init_all_encompassing_struct(nb::module_ &m, nb::class_<AllEncompassingStruct> &class_);
