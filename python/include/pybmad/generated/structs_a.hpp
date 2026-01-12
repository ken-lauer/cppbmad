#pragma once
#include <pybind11/pybind11.h>
namespace py = pybind11;

void init_structs_a(py::module& m);

// Per-struct init functions
void init_ac_kicker_freq_struct(
    py::module& m,
    py::class_<AcKickerFreqProxy>& class_);
void init_ac_kicker_struct(py::module& m, py::class_<AcKickerProxy>& class_);
void init_ac_kicker_time_struct(
    py::module& m,
    py::class_<AcKickerTimeProxy>& class_);
void init_anormal_mode_struct(
    py::module& m,
    py::class_<AnormalModeProxy>& class_);
void init_aperture_param_struct(
    py::module& m,
    py::class_<ApertureParamProxy>& class_);
void init_aperture_point_struct(
    py::module& m,
    py::class_<AperturePointProxy>& class_);
void init_aperture_scan_struct(
    py::module& m,
    py::class_<ApertureScanProxy>& class_);
void init_all_encompassing_struct(
    py::module& m,
    py::class_<AllEncompassingProxy>& class_);
