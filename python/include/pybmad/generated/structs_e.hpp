#pragma once
#include <pybind11/pybind11.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace py = pybind11;

using namespace Bmad;

// Per-struct init functions
void init_ele_pointer_struct(py::module &m, py::class_<ElePointerStruct> &class_);
void init_ele_struct(py::module &m, py::class_<EleStruct> &class_);
void init_ellipse_beam_init_struct(py::module &m, py::class_<EllipseBeamInitStruct> &class_);
void init_em_field_struct(py::module &m, py::class_<EmFieldStruct> &class_);
void init_expression_atom_struct(py::module &m, py::class_<ExpressionAtomStruct> &class_);
void init_expression_tree_struct(py::module &m, py::class_<ExpressionTreeStruct> &class_);
