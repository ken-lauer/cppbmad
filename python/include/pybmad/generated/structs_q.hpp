#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_qp_axis_struct(nb::module_ &m, nb::class_<QpAxisStruct> &class_);
void init_qp_legend_struct(nb::module_ &m, nb::class_<QpLegendStruct> &class_);
void init_qp_line_struct(nb::module_ &m, nb::class_<QpLineStruct> &class_);
void init_qp_point_struct(nb::module_ &m, nb::class_<QpPointStruct> &class_);
void init_qp_rect_struct(nb::module_ &m, nb::class_<QpRectStruct> &class_);
void init_qp_symbol_struct(nb::module_ &m, nb::class_<QpSymbolStruct> &class_);
