#pragma once
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Per-struct init functions
void init_qp_axis_struct(py::module& m, py::class_<QpAxisProxy>& class_);
void init_qp_legend_struct(py::module& m, py::class_<QpLegendProxy>& class_);
void init_qp_line_struct(py::module& m, py::class_<QpLineProxy>& class_);
void init_qp_point_struct(py::module& m, py::class_<QpPointProxy>& class_);
void init_qp_rect_struct(py::module& m, py::class_<QpRectProxy>& class_);
void init_qp_symbol_struct(py::module& m, py::class_<QpSymbolProxy>& class_);
