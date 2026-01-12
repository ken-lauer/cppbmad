#include "pybmad/generated/structs_f.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// floor_position_struct
void init_floor_position_struct(
    py::module& m,
    py::class_<FloorPositionProxy>& cls) {
  cls.def(py::init<>())
      // FloorPositionProxy.r (1D_NOT_real - (x, y, z) offset from origin
      .def_property_readonly("r", &FloorPositionProxy::r)
      // FloorPositionProxy.w (2D_NOT_real - W matrix. Columns are unit vectors of the frame axes.
      .def_property_readonly("w", &FloorPositionProxy::w)
      // FloorPositionProxy.theta (0D_NOT_real - angular orientation consistent with W matrix
      .def_property(
          "theta", &FloorPositionProxy::theta, &FloorPositionProxy::set_theta)
      // FloorPositionProxy.phi (0D_NOT_real - angular orientation consistent with W matrix
      .def_property(
          "phi", &FloorPositionProxy::phi, &FloorPositionProxy::set_phi)
      // FloorPositionProxy.psi (0D_NOT_real - angular orientation consistent with W matrix
      .def_property(
          "psi", &FloorPositionProxy::psi, &FloorPositionProxy::set_psi)

      .def(
          "__repr__",
          [](const FloorPositionProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const FloorPositionProxy& self) {
            return FloorPositionProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const FloorPositionProxy& self, py::dict& memo) {
            return FloorPositionProxy(self);
          })

      ;

  // 1D FloorPositionProxy arrays are not used in structs/routines
  // 2D FloorPositionProxy arrays are not used in structs/routines
  // 3D FloorPositionProxy arrays are not used in structs/routines
}