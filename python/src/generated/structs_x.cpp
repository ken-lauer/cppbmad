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
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("eta") = py::none(),
         py::arg("etap") = py::none(),
         py::arg("deta_ds") = py::none(),
         py::arg("sigma") = py::none(),
         py::arg("deta_dpz") = py::none(),
         py::arg("detap_dpz") = py::none()
  )
      .def_property("eta", &XyDispStruct::eta, &XyDispStruct::set_eta)
      .def_property("etap", &XyDispStruct::etap, &XyDispStruct::set_etap)
      .def_property("deta_ds", &XyDispStruct::deta_ds, &XyDispStruct::set_deta_ds)
      .def_property("sigma", &XyDispStruct::sigma, &XyDispStruct::set_sigma)
      .def_property("deta_dpz", &XyDispStruct::deta_dpz, &XyDispStruct::set_deta_dpz)
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