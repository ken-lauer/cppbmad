#include "pybmad/generated/structs_x.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// xy_disp_struct
void init_xy_disp_struct(py::module& m, py::class_<XyDispProxy>& cls) {
  cls.def(py::init<>())
      // XyDispProxy.eta (0D_NOT_real -
      .def_property("eta", &XyDispProxy::eta, &XyDispProxy::set_eta)
      // XyDispProxy.etap (0D_NOT_real -
      .def_property("etap", &XyDispProxy::etap, &XyDispProxy::set_etap)
      // XyDispProxy.deta_ds (0D_NOT_real -
      .def_property("deta_ds", &XyDispProxy::deta_ds, &XyDispProxy::set_deta_ds)
      // XyDispProxy.sigma (0D_NOT_real -
      .def_property("sigma", &XyDispProxy::sigma, &XyDispProxy::set_sigma)
      // XyDispProxy.deta_dpz (0D_NOT_real -
      .def_property(
          "deta_dpz", &XyDispProxy::deta_dpz, &XyDispProxy::set_deta_dpz)
      // XyDispProxy.detap_dpz (0D_NOT_real -
      .def_property(
          "detap_dpz", &XyDispProxy::detap_dpz, &XyDispProxy::set_detap_dpz)

      .def("__repr__", [](const XyDispProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const XyDispProxy& self) {
            return XyDispProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const XyDispProxy& self, py::dict& memo) {
            return XyDispProxy(self);
          })

      ;

  // 1D XyDispProxy arrays are not used in structs/routines
  // 2D XyDispProxy arrays are not used in structs/routines
  // 3D XyDispProxy arrays are not used in structs/routines
}