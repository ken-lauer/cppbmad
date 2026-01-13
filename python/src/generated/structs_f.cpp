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

// =============================================================================
// fibre
void init_fibre(py::module& m, py::class_<FibreRawStruct>& cls) {
  cls.def(py::init<>())
      // FibreRawStruct.DIR (0D_PTR_integer -
      .def_property("DIR", &FibreRawStruct::DIR, &FibreRawStruct::set_DIR)
      // FibreRawStruct.PREVIOUS (0D_PTR_type -
      .def_property(
          "PREVIOUS", &FibreRawStruct::PREVIOUS, &FibreRawStruct::set_PREVIOUS)
      // FibreRawStruct.NEXT (0D_PTR_type - POINTING TO PARENT LAYOUT AND PARENT FIBRE DATA
      .def_property("NEXT", &FibreRawStruct::NEXT, &FibreRawStruct::set_NEXT)
      // FibreRawStruct.PARENT_LAYOUT (0D_PTR_type -
      .def_property(
          "PARENT_LAYOUT",
          &FibreRawStruct::PARENT_LAYOUT,
          &FibreRawStruct::set_PARENT_LAYOUT)
      // FibreRawStruct.pos (0D_PTR_integer - POSITION IN LAYOUT NEW STUFF....
      .def_property("pos", &FibreRawStruct::pos, &FibreRawStruct::set_pos)
      // FibreRawStruct.BETA0 (0D_PTR_real - ,P0C
      .def_property("BETA0", &FibreRawStruct::BETA0, &FibreRawStruct::set_BETA0)
      // FibreRawStruct.GAMMA0I (0D_PTR_real - ,P0C
      .def_property(
          "GAMMA0I", &FibreRawStruct::GAMMA0I, &FibreRawStruct::set_GAMMA0I)
      // FibreRawStruct.GAMBET (0D_PTR_real - ,P0C
      .def_property(
          "GAMBET", &FibreRawStruct::GAMBET, &FibreRawStruct::set_GAMBET)
      // FibreRawStruct.MASS (0D_PTR_real - ,P0C
      .def_property("MASS", &FibreRawStruct::MASS, &FibreRawStruct::set_MASS)
      // FibreRawStruct.CHARGE (0D_PTR_real -
      .def_property(
          "CHARGE", &FibreRawStruct::CHARGE, &FibreRawStruct::set_CHARGE)
      // FibreRawStruct.AG (0D_PTR_real - spin g-2 TO TIE LAYOUTS
      .def_property("AG", &FibreRawStruct::AG, &FibreRawStruct::set_AG)
      // FibreRawStruct.P (0D_PTR_type - tying them in the so-called database universe M_u
      .def_property("P", &FibreRawStruct::P, &FibreRawStruct::set_P)
      // FibreRawStruct.N (0D_PTR_type -
      .def_property("N", &FibreRawStruct::N, &FibreRawStruct::set_N)
      // FibreRawStruct.loc (0D_PTR_integer -
      .def_property("loc", &FibreRawStruct::loc, &FibreRawStruct::set_loc)

      .def(
          "__repr__",
          [](const FibreRawStruct& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const FibreRawStruct& self) {
            return FibreRawStruct(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const FibreRawStruct& self, py::dict& memo) {
            return FibreRawStruct(self);
          })

      ;

  // 1D FibreRawStruct arrays are not used in structs/routines
  // 2D FibreRawStruct arrays are not used in structs/routines
  // 3D FibreRawStruct arrays are not used in structs/routines
}