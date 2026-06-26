#pragma once
#include <nanobind/nanobind.h>

#include "bmad/generated/proxy.hpp"
#include "pybmad/generated/structs.hpp"
namespace nb = nanobind;

using namespace Bmad;

// Per-struct init functions
void init_interval1_coef_struct(nb::module_ &m, nb::class_<Interval1CoefStruct> &class_);
void init_ibs_lifetime_struct(nb::module_ &m, nb::class_<IbsLifetimeStruct> &class_);
void init_ibs_maxratio_struct(nb::module_ &m, nb::class_<IbsMaxratioStruct> &class_);
void init_ibs_sim_param_struct(nb::module_ &m, nb::class_<IbsSimParamStruct> &class_);
void init_ibs_struct(nb::module_ &m, nb::class_<IbsStruct> &class_);
