#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_ele_attribute_struct(nb::module_ &m, nb::class_<EleAttributeStruct> &class_);
void init_ele_pointer_struct(nb::module_ &m, nb::class_<ElePointerStruct> &class_);
void init_ele_struct(nb::module_ &m, nb::class_<EleStruct> &class_);
void init_ellipse_beam_init_struct(nb::module_ &m, nb::class_<EllipseBeamInitStruct> &class_);
void init_em_field_struct(nb::module_ &m, nb::class_<EmFieldStruct> &class_);
void init_expression_atom_struct(nb::module_ &m, nb::class_<ExpressionAtomStruct> &class_);
void init_expression_tree_struct(nb::module_ &m, nb::class_<ExpressionTreeStruct> &class_);
void init_extra_parsing_info_struct(nb::module_ &m, nb::class_<ExtraParsingInfoStruct> &class_);
