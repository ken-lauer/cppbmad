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
  cls.def(
         py::init<
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<std::vector<double>>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("r") = py::none(),
         py::arg("w") = py::none(),
         py::arg("theta") = py::none(),
         py::arg("phi") = py::none(),
         py::arg("psi") = py::none()
  )
      .def_property(
          "r",
          &FloorPositionStruct::r,
          &FloorPositionStruct::set_r,
          "(x, y, z) offset from origin"
      )
      .def_property(
          "w",
          &FloorPositionStruct::w,
          &FloorPositionStruct::set_w,
          "W matrix. Columns are unit vectors of the frame axes."
      )
      .def_property(
          "theta",
          &FloorPositionStruct::theta,
          &FloorPositionStruct::set_theta,
          "angular orientation consistent with W matrix"
      )
      .def_property(
          "phi",
          &FloorPositionStruct::phi,
          &FloorPositionStruct::set_phi,
          "angular orientation consistent with W matrix"
      )
      .def_property(
          "psi",
          &FloorPositionStruct::psi,
          &FloorPositionStruct::set_psi,
          "angular orientation consistent with W matrix"
      )

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
  cls.def(
         py::init<
             std::optional<int>,
             optional_ref<const Fibre>,
             optional_ref<const Fibre>,
             optional_ref<const Layout>,
             std::optional<int>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             optional_ref<const Fibre>,
             optional_ref<const Fibre>,
             std::optional<int>>(),
         py::arg("DIR") = py::none(),
         py::arg("PREVIOUS") = py::none(),
         py::arg("NEXT") = py::none(),
         py::arg("PARENT_LAYOUT") = py::none(),
         py::arg("pos") = py::none(),
         py::arg("BETA0") = py::none(),
         py::arg("GAMMA0I") = py::none(),
         py::arg("GAMBET") = py::none(),
         py::arg("MASS") = py::none(),
         py::arg("CHARGE") = py::none(),
         py::arg("AG") = py::none(),
         py::arg("P") = py::none(),
         py::arg("N") = py::none(),
         py::arg("loc") = py::none()
  )
      .def_property("DIR", &Fibre::DIR, &Fibre::set_DIR)
      .def_property("PREVIOUS", &Fibre::PREVIOUS, &Fibre::set_PREVIOUS)
      .def_property(
          "NEXT",
          &Fibre::NEXT,
          &Fibre::set_NEXT,
          "POINTING TO PARENT LAYOUT AND PARENT FIBRE DATA"
      )
      .def_property("PARENT_LAYOUT", &Fibre::PARENT_LAYOUT, &Fibre::set_PARENT_LAYOUT)
      .def_property("pos", &Fibre::pos, &Fibre::set_pos, "POSITION IN LAYOUT NEW STUFF....")
      .def_property("BETA0", &Fibre::BETA0, &Fibre::set_BETA0, ",P0C")
      .def_property("GAMMA0I", &Fibre::GAMMA0I, &Fibre::set_GAMMA0I, ",P0C")
      .def_property("GAMBET", &Fibre::GAMBET, &Fibre::set_GAMBET, ",P0C")
      .def_property("MASS", &Fibre::MASS, &Fibre::set_MASS, ",P0C")
      .def_property("CHARGE", &Fibre::CHARGE, &Fibre::set_CHARGE)
      .def_property("AG", &Fibre::AG, &Fibre::set_AG, "spin g-2 TO TIE LAYOUTS")
      .def_property(
          "P",
          &Fibre::P,
          &Fibre::set_P,
          "tying them in the so-called database universe M_u"
      )
      .def_property("N", &Fibre::N, &Fibre::set_N)
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