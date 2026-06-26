#include "pybmad/generated/structs_f.hpp"

#include <cstdint>
#include <functional>

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

using namespace Pybmad;
namespace nb = nanobind;

// =============================================================================
// floor_position_struct
void init_floor_position_struct(nb::module_ &m, nb::class_<FloorPositionStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("r") = nb::none(),
         nb::arg("w") = nb::none(),
         nb::arg("theta") = nb::none(),
         nb::arg("phi") = nb::none(),
         nb::arg("psi") = nb::none()
  )
      .def_prop_rw(
          "r",
          &FloorPositionStruct::r,
          &FloorPositionStruct::set_r,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "(x, y, z) offset from origin"
      )
      .def_prop_rw(
          "w",
          &FloorPositionStruct::w,
          &FloorPositionStruct::set_w,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "W matrix. Columns are unit vectors of the frame axes."
      )
      .def_prop_rw(
          "theta",
          &FloorPositionStruct::theta,
          &FloorPositionStruct::set_theta,
          "angular orientation consistent with W matrix"
      )
      .def_prop_rw(
          "phi",
          &FloorPositionStruct::phi,
          &FloorPositionStruct::set_phi,
          "angular orientation consistent with W matrix"
      )
      .def_prop_rw(
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
          [](const FloorPositionStruct &self, nb::dict &memo) { return FloorPositionStruct(self); }
      )
      .def(
          "__eq__",
          [](const FloorPositionStruct &self, const FloorPositionStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const FloorPositionStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D FloorPositionStruct arrays are not used in structs/routines
  // 2D FloorPositionStruct arrays are not used in structs/routines
  // 3D FloorPositionStruct arrays are not used in structs/routines
}

// =============================================================================
// fibre
void init_fibre(nb::module_ &m, nb::class_<Fibre> &cls) {
  cls.def(
         "__init__",
         [](Fibre *self,
            std::optional<int> DIR,
            const Fibre *PREVIOUS,
            const Fibre *NEXT,
            const Layout *PARENT_LAYOUT,
            std::optional<int> pos,
            std::optional<double> BETA0,
            std::optional<double> GAMMA0I,
            std::optional<double> GAMBET,
            std::optional<double> MASS,
            std::optional<double> CHARGE,
            std::optional<double> AG,
            const Fibre *P,
            const Fibre *N,
            std::optional<int> loc) {
           new (self) Fibre(
               DIR,
               ptr_to_opt_ref(PREVIOUS),
               ptr_to_opt_ref(NEXT),
               ptr_to_opt_ref(PARENT_LAYOUT),
               pos,
               BETA0,
               GAMMA0I,
               GAMBET,
               MASS,
               CHARGE,
               AG,
               ptr_to_opt_ref(P),
               ptr_to_opt_ref(N),
               loc
           );
         },
         nb::arg("DIR") = nb::none(),
         nb::arg("PREVIOUS") = nb::none(),
         nb::arg("NEXT") = nb::none(),
         nb::arg("PARENT_LAYOUT") = nb::none(),
         nb::arg("pos") = nb::none(),
         nb::arg("BETA0") = nb::none(),
         nb::arg("GAMMA0I") = nb::none(),
         nb::arg("GAMBET") = nb::none(),
         nb::arg("MASS") = nb::none(),
         nb::arg("CHARGE") = nb::none(),
         nb::arg("AG") = nb::none(),
         nb::arg("P") = nb::none(),
         nb::arg("N") = nb::none(),
         nb::arg("loc") = nb::none()
  )
      .def_prop_rw("DIR", &Fibre::DIR, &Fibre::set_DIR)
      .def_prop_rw(
          "PREVIOUS",
          &Fibre::PREVIOUS,
          &Fibre::set_PREVIOUS,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "NEXT",
          &Fibre::NEXT,
          &Fibre::set_NEXT,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "POINTING TO PARENT LAYOUT AND PARENT FIBRE DATA"
      )
      .def_prop_rw(
          "PARENT_LAYOUT",
          &Fibre::PARENT_LAYOUT,
          &Fibre::set_PARENT_LAYOUT,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw("pos", &Fibre::pos, &Fibre::set_pos, "POSITION IN LAYOUT NEW STUFF....")
      .def_prop_rw("BETA0", &Fibre::BETA0, &Fibre::set_BETA0, ",P0C")
      .def_prop_rw("GAMMA0I", &Fibre::GAMMA0I, &Fibre::set_GAMMA0I, ",P0C")
      .def_prop_rw("GAMBET", &Fibre::GAMBET, &Fibre::set_GAMBET, ",P0C")
      .def_prop_rw("MASS", &Fibre::MASS, &Fibre::set_MASS, ",P0C")
      .def_prop_rw("CHARGE", &Fibre::CHARGE, &Fibre::set_CHARGE)
      .def_prop_rw("AG", &Fibre::AG, &Fibre::set_AG, "spin g-2 TO TIE LAYOUTS")
      .def_prop_rw(
          "P",
          &Fibre::P,
          &Fibre::set_P,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "tying them in the so-called database universe M_u"
      )
      .def_prop_rw("N", &Fibre::N, &Fibre::set_N, nb::for_getter(nb::keep_alive<0, 1>()))
      .def_prop_rw("loc", &Fibre::loc, &Fibre::set_loc)

      .def("__repr__", [](const Fibre &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Fibre &self) {
            return Fibre(self); // under-the-hood fortran copy
          }
      )
      .def("__deepcopy__", [](const Fibre &self, nb::dict &memo) { return Fibre(self); })
      .def(
          "__eq__",
          [](const Fibre &self, const Fibre &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const Fibre &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D Fibre arrays are not used in structs/routines
  // 2D Fibre arrays are not used in structs/routines
  // 3D Fibre arrays are not used in structs/routines
}