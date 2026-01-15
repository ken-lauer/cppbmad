#include "pybmad/generated/structs_x.hpp"

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// xy_disp_struct
void init_xy_disp_struct(py::module &m, py::class_<XyDispStruct> &cls) {
  cls.def(py::init<>())
      // XyDispStruct.eta (0D_NOT_real -
      .def_property("eta", &XyDispStruct::eta, &XyDispStruct::set_eta)
      // XyDispStruct.etap (0D_NOT_real -
      .def_property("etap", &XyDispStruct::etap, &XyDispStruct::set_etap)
      // XyDispStruct.deta_ds (0D_NOT_real -
      .def_property("deta_ds", &XyDispStruct::deta_ds, &XyDispStruct::set_deta_ds)
      // XyDispStruct.sigma (0D_NOT_real -
      .def_property("sigma", &XyDispStruct::sigma, &XyDispStruct::set_sigma)
      // XyDispStruct.deta_dpz (0D_NOT_real -
      .def_property("deta_dpz", &XyDispStruct::deta_dpz, &XyDispStruct::set_deta_dpz)
      // XyDispStruct.detap_dpz (0D_NOT_real -
      .def_property("detap_dpz", &XyDispStruct::detap_dpz, &XyDispStruct::set_detap_dpz)

      .def("__repr__", [](const XyDispStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const XyDispStruct &self) {
            return XyDispStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const XyDispStruct &self, py::dict &memo) { return XyDispStruct(self); }
      )

      ;

  // 1D XyDispStruct arrays are not used in structs/routines
  // 2D XyDispStruct arrays are not used in structs/routines
  // 3D XyDispStruct arrays are not used in structs/routines
}