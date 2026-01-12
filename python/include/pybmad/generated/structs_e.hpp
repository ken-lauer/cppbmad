#pragma once
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Per-struct init functions
void init_ele_pointer_struct(
    py::module& m,
    py::class_<ElePointerProxy>& class_);
void init_ele_struct(py::module& m, py::class_<EleProxy>& class_);
void init_ellipse_beam_init_struct(
    py::module& m,
    py::class_<EllipseBeamInitProxy>& class_);
void init_em_field_struct(py::module& m, py::class_<EmFieldProxy>& class_);
void init_em_taylor_struct(py::module& m, py::class_<EmTaylorProxy>& class_);
void init_em_taylor_term_struct(
    py::module& m,
    py::class_<EmTaylorTermProxy>& class_);
void init_expression_atom_struct(
    py::module& m,
    py::class_<ExpressionAtomProxy>& class_);
void init_expression_tree_struct(
    py::module& m,
    py::class_<ExpressionTreeProxy>& class_);
