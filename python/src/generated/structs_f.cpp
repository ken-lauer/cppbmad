#include "pybmad/generated/structs_f.hpp"

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// floor_position_struct
void init_floor_position_struct(py::module &m, py::class_<FloorPositionStruct> &cls) {
  cls.def(py::init<>())
      // FloorPositionStruct.r (1D_NOT_real - (x, y, z) offset from origin
      .def_property_readonly("r", &FloorPositionStruct::r)
      // FloorPositionStruct.w (2D_NOT_real - W matrix. Columns are unit vectors of the frame axes.
      .def_property_readonly("w", &FloorPositionStruct::w)
      // FloorPositionStruct.theta (0D_NOT_real - angular orientation consistent with W matrix
      .def_property("theta", &FloorPositionStruct::theta, &FloorPositionStruct::set_theta)
      // FloorPositionStruct.phi (0D_NOT_real - angular orientation consistent with W matrix
      .def_property("phi", &FloorPositionStruct::phi, &FloorPositionStruct::set_phi)
      // FloorPositionStruct.psi (0D_NOT_real - angular orientation consistent with W matrix
      .def_property("psi", &FloorPositionStruct::psi, &FloorPositionStruct::set_psi)

      .def("__repr__", [](const FloorPositionStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const FloorPositionStruct &self) {
            return FloorPositionStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const FloorPositionStruct &self, py::dict &memo) { return FloorPositionStruct(self); }
      )

      ;

  // 1D FloorPositionStruct arrays are not used in structs/routines
  // 2D FloorPositionStruct arrays are not used in structs/routines
  // 3D FloorPositionStruct arrays are not used in structs/routines
}

// =============================================================================
// fibre
void init_fibre(py::module &m, py::class_<Fibre> &cls) {
  cls.def(py::init<>())
      // Fibre.DIR (0D_PTR_integer -
      .def_property("DIR", &Fibre::DIR, &Fibre::set_DIR)
      // Fibre.PREVIOUS (0D_PTR_type -
      .def_property("PREVIOUS", &Fibre::PREVIOUS, &Fibre::set_PREVIOUS)
      // Fibre.NEXT (0D_PTR_type - POINTING TO PARENT LAYOUT AND PARENT FIBRE DATA
      .def_property("NEXT", &Fibre::NEXT, &Fibre::set_NEXT)
      // Fibre.PARENT_LAYOUT (0D_PTR_type -
      .def_property("PARENT_LAYOUT", &Fibre::PARENT_LAYOUT, &Fibre::set_PARENT_LAYOUT)
      // Fibre.pos (0D_PTR_integer - POSITION IN LAYOUT NEW STUFF....
      .def_property("pos", &Fibre::pos, &Fibre::set_pos)
      // Fibre.BETA0 (0D_PTR_real - ,P0C
      .def_property("BETA0", &Fibre::BETA0, &Fibre::set_BETA0)
      // Fibre.GAMMA0I (0D_PTR_real - ,P0C
      .def_property("GAMMA0I", &Fibre::GAMMA0I, &Fibre::set_GAMMA0I)
      // Fibre.GAMBET (0D_PTR_real - ,P0C
      .def_property("GAMBET", &Fibre::GAMBET, &Fibre::set_GAMBET)
      // Fibre.MASS (0D_PTR_real - ,P0C
      .def_property("MASS", &Fibre::MASS, &Fibre::set_MASS)
      // Fibre.CHARGE (0D_PTR_real -
      .def_property("CHARGE", &Fibre::CHARGE, &Fibre::set_CHARGE)
      // Fibre.AG (0D_PTR_real - spin g-2 TO TIE LAYOUTS
      .def_property("AG", &Fibre::AG, &Fibre::set_AG)
      // Fibre.P (0D_PTR_type - tying them in the so-called database universe M_u
      .def_property("P", &Fibre::P, &Fibre::set_P)
      // Fibre.N (0D_PTR_type -
      .def_property("N", &Fibre::N, &Fibre::set_N)
      // Fibre.loc (0D_PTR_integer -
      .def_property("loc", &Fibre::loc, &Fibre::set_loc)

      .def("__repr__", [](const Fibre &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Fibre &self) {
            return Fibre(self); // under-the-hood fortran copy
          }
      )
      .def("__deepcopy__", [](const Fibre &self, py::dict &memo) { return Fibre(self); })

      ;

  // 1D Fibre arrays are not used in structs/routines
  // 2D Fibre arrays are not used in structs/routines
  // 3D Fibre arrays are not used in structs/routines
}