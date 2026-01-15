#pragma once
#include <pybind11/pybind11.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace py = pybind11;

using namespace Bmad;

// Per-struct init functions
void init_ac_kicker_freq_struct(py::module &m, py::class_<AcKickerFreqStruct> &class_);
void init_ac_kicker_struct(py::module &m, py::class_<AcKickerStruct> &class_);
void init_ac_kicker_time_struct(py::module &m, py::class_<AcKickerTimeStruct> &class_);
void init_anormal_mode_struct(py::module &m, py::class_<AnormalModeStruct> &class_);
void init_aperture_param_struct(py::module &m, py::class_<ApertureParamStruct> &class_);
void init_aperture_point_struct(py::module &m, py::class_<AperturePointStruct> &class_);
void init_aperture_scan_struct(py::module &m, py::class_<ApertureScanStruct> &class_);
void init_all_encompassing_struct(py::module &m, py::class_<AllEncompassingStruct> &class_);
