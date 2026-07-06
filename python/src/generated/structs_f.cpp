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
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D FloorPositionStruct arrays are not used in structs/routines
  // 2D FloorPositionStruct arrays are not used in structs/routines
  // 3D FloorPositionStruct arrays are not used in structs/routines
}

// =============================================================================
// foil_struct
void init_foil_struct(nb::module_ &m, nb::class_<FoilStruct> &cls) {
  cls.def(nb::init<>())
      .def_prop_ro("material", &FoilStruct::material, nb::keep_alive<0, 1>())

      .def("__repr__", [](const FoilStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const FoilStruct &self) {
            return FoilStruct(self); // under-the-hood fortran copy
          }
      )
      .def("__deepcopy__", [](const FoilStruct &self, nb::dict &memo) { return FoilStruct(self); })
      .def(
          "__eq__",
          [](const FoilStruct &self, const FoilStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const FoilStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D FoilStruct arrays are not used in structs/routines
  // 2D FoilStruct arrays are not used in structs/routines
  // 3D FoilStruct arrays are not used in structs/routines
}

// =============================================================================
// fringe_field_info_struct
void init_fringe_field_info_struct(nb::module_ &m, nb::class_<FringeFieldInfoStruct> &cls) {
  cls.def(
         "__init__",
         [](FringeFieldInfoStruct *self,
            const EleStruct *hard_ele,
            std::optional<double> s_edge_hard,
            std::optional<double> ds_edge,
            std::optional<int> particle_at,
            std::optional<int> hard_location,
            std::optional<std::vector<int>> location,
            std::optional<bool> has_fringe) {
           new (self) FringeFieldInfoStruct(
               ptr_to_opt_ref(hard_ele),
               s_edge_hard,
               ds_edge,
               particle_at,
               hard_location,
               location,
               has_fringe
           );
         },
         nb::arg("hard_ele") = nb::none(),
         nb::arg("s_edge_hard") = nb::none(),
         nb::arg("ds_edge") = nb::none(),
         nb::arg("particle_at") = nb::none(),
         nb::arg("hard_location") = nb::none(),
         nb::arg("location") = nb::none(),
         nb::arg("has_fringe") = nb::none()
  )
      .def_prop_rw(
          "hard_ele",
          &FringeFieldInfoStruct::hard_ele,
          &FringeFieldInfoStruct::set_hard_ele,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "s_edge_hard",
          &FringeFieldInfoStruct::s_edge_hard,
          &FringeFieldInfoStruct::set_s_edge_hard
      )
      .def_prop_rw(
          "ds_edge",
          &FringeFieldInfoStruct::ds_edge,
          &FringeFieldInfoStruct::set_ds_edge,
          "Distance from particle to edge in hard_ele frame."
      )
      .def_prop_rw(
          "particle_at",
          &FringeFieldInfoStruct::particle_at,
          &FringeFieldInfoStruct::set_particle_at,
          "first_track_edge$, second_track_edge$, or none$"
      )
      .def_prop_rw(
          "hard_location",
          &FringeFieldInfoStruct::hard_location,
          &FringeFieldInfoStruct::set_hard_location,
          "Particle location wrt hard_ele. Points to element in location(:)."
      )
      .def_prop_rw(
          "location",
          &FringeFieldInfoStruct::location,
          &FringeFieldInfoStruct::set_location,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Particle location in an element. entrance_end$, inside$, or exit_end$ Elements in list "
          "are the tracking element or its lords."
      )
      .def_prop_rw(
          "has_fringe",
          &FringeFieldInfoStruct::has_fringe,
          &FringeFieldInfoStruct::set_has_fringe,
          "Has a fringe to worry about?"
      )

      .def("__repr__", [](const FringeFieldInfoStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const FringeFieldInfoStruct &self) {
            return FringeFieldInfoStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const FringeFieldInfoStruct &self, nb::dict &memo) {
            return FringeFieldInfoStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const FringeFieldInfoStruct &self, const FringeFieldInfoStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const FringeFieldInfoStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D FringeFieldInfoStruct arrays are not used in structs/routines
  // 2D FringeFieldInfoStruct arrays are not used in structs/routines
  // 3D FringeFieldInfoStruct arrays are not used in structs/routines
}

// =============================================================================
// field1_at_2D_pt_struct
void init_field1_at_2D_pt_struct(nb::module_ &m, nb::class_<Field1At2dPtStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("f") = nb::none(),
         nb::arg("df_dx") = nb::none(),
         nb::arg("df_dy") = nb::none(),
         nb::arg("d2f_dxdy") = nb::none()
  )
      .def_prop_rw("f", &Field1At2dPtStruct::f, &Field1At2dPtStruct::set_f, "Field")
      .def_prop_rw(
          "df_dx",
          &Field1At2dPtStruct::df_dx,
          &Field1At2dPtStruct::set_df_dx,
          "Normalized field 1st derivatives"
      )
      .def_prop_rw(
          "df_dy",
          &Field1At2dPtStruct::df_dy,
          &Field1At2dPtStruct::set_df_dy,
          "Normalized field 1st derivatives"
      )
      .def_prop_rw(
          "d2f_dxdy",
          &Field1At2dPtStruct::d2f_dxdy,
          &Field1At2dPtStruct::set_d2f_dxdy,
          "Normalized field 2nd derivative"
      )

      .def("__repr__", [](const Field1At2dPtStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Field1At2dPtStruct &self) {
            return Field1At2dPtStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const Field1At2dPtStruct &self, nb::dict &memo) { return Field1At2dPtStruct(self); }
      )
      .def(
          "__eq__",
          [](const Field1At2dPtStruct &self, const Field1At2dPtStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const Field1At2dPtStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D Field1At2dPtStruct arrays are not used in structs/routines
  bind_FTypeArrayND<Field1At2dPtStructArray2D>(m, "Field1At2DPtStructArray2D");
  // 3D Field1At2dPtStruct arrays are not used in structs/routines
}

// =============================================================================
// field1_at_3D_pt_struct
void init_field1_at_3D_pt_struct(nb::module_ &m, nb::class_<Field1At3dPtStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("f") = nb::none(),
         nb::arg("df_dx") = nb::none(),
         nb::arg("df_dy") = nb::none(),
         nb::arg("df_dz") = nb::none(),
         nb::arg("d2f_dxdy") = nb::none(),
         nb::arg("d2f_dxdz") = nb::none(),
         nb::arg("d2f_dydz") = nb::none(),
         nb::arg("d3f_dxdydz") = nb::none()
  )
      .def_prop_rw("f", &Field1At3dPtStruct::f, &Field1At3dPtStruct::set_f, "Field")
      .def_prop_rw(
          "df_dx",
          &Field1At3dPtStruct::df_dx,
          &Field1At3dPtStruct::set_df_dx,
          "Normalized field 1st derivatives"
      )
      .def_prop_rw(
          "df_dy",
          &Field1At3dPtStruct::df_dy,
          &Field1At3dPtStruct::set_df_dy,
          "Normalized field 1st derivatives"
      )
      .def_prop_rw(
          "df_dz",
          &Field1At3dPtStruct::df_dz,
          &Field1At3dPtStruct::set_df_dz,
          "Normalized field 1st derivatives"
      )
      .def_prop_rw(
          "d2f_dxdy",
          &Field1At3dPtStruct::d2f_dxdy,
          &Field1At3dPtStruct::set_d2f_dxdy,
          "Normalized field 2nd derivatives"
      )
      .def_prop_rw(
          "d2f_dxdz",
          &Field1At3dPtStruct::d2f_dxdz,
          &Field1At3dPtStruct::set_d2f_dxdz,
          "Normalized field 2nd derivatives"
      )
      .def_prop_rw(
          "d2f_dydz",
          &Field1At3dPtStruct::d2f_dydz,
          &Field1At3dPtStruct::set_d2f_dydz,
          "Normalized field 2nd derivatives"
      )
      .def_prop_rw(
          "d3f_dxdydz",
          &Field1At3dPtStruct::d3f_dxdydz,
          &Field1At3dPtStruct::set_d3f_dxdydz,
          "Normalized field 3rd derivative"
      )

      .def("__repr__", [](const Field1At3dPtStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Field1At3dPtStruct &self) {
            return Field1At3dPtStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const Field1At3dPtStruct &self, nb::dict &memo) { return Field1At3dPtStruct(self); }
      )
      .def(
          "__eq__",
          [](const Field1At3dPtStruct &self, const Field1At3dPtStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const Field1At3dPtStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D Field1At3dPtStruct arrays are not used in structs/routines
  // 2D Field1At3dPtStruct arrays are not used in structs/routines
  bind_FTypeArrayND<Field1At3dPtStructArray3D>(m, "Field1At3DPtStructArray3D");
}

// =============================================================================
// field_at_2D_box_struct
void init_field_at_2D_box_struct(nb::module_ &m, nb::class_<FieldAt2dBoxStruct> &cls) {
  cls.def(nb::init<std::optional<std::vector<int>>>(), nb::arg("i_box") = nb::none())
      .def_prop_ro("pt", &FieldAt2dBoxStruct::pt, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "i_box",
          &FieldAt2dBoxStruct::i_box,
          &FieldAt2dBoxStruct::set_i_box,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "index at lower box corner."
      )

      .def("__repr__", [](const FieldAt2dBoxStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const FieldAt2dBoxStruct &self) {
            return FieldAt2dBoxStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const FieldAt2dBoxStruct &self, nb::dict &memo) { return FieldAt2dBoxStruct(self); }
      )
      .def(
          "__eq__",
          [](const FieldAt2dBoxStruct &self, const FieldAt2dBoxStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const FieldAt2dBoxStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D FieldAt2dBoxStruct arrays are not used in structs/routines
  // 2D FieldAt2dBoxStruct arrays are not used in structs/routines
  // 3D FieldAt2dBoxStruct arrays are not used in structs/routines
}

// =============================================================================
// field_at_3D_box_struct
void init_field_at_3D_box_struct(nb::module_ &m, nb::class_<FieldAt3dBoxStruct> &cls) {
  cls.def(nb::init<std::optional<std::vector<int>>>(), nb::arg("i_box") = nb::none())
      .def_prop_ro("pt", &FieldAt3dBoxStruct::pt, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "i_box",
          &FieldAt3dBoxStruct::i_box,
          &FieldAt3dBoxStruct::set_i_box,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "index at lower box corner."
      )

      .def("__repr__", [](const FieldAt3dBoxStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const FieldAt3dBoxStruct &self) {
            return FieldAt3dBoxStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const FieldAt3dBoxStruct &self, nb::dict &memo) { return FieldAt3dBoxStruct(self); }
      )
      .def(
          "__eq__",
          [](const FieldAt3dBoxStruct &self, const FieldAt3dBoxStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const FieldAt3dBoxStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D FieldAt3dBoxStruct arrays are not used in structs/routines
  // 2D FieldAt3dBoxStruct arrays are not used in structs/routines
  // 3D FieldAt3dBoxStruct arrays are not used in structs/routines
}

// =============================================================================
// fibre
void init_fibre(nb::module_ &m, nb::class_<Fibre> &cls) {
  cls.def(
         nb::init<
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>>(),
         nb::arg("DIR") = nb::none(),
         nb::arg("pos") = nb::none(),
         nb::arg("BETA0") = nb::none(),
         nb::arg("GAMMA0I") = nb::none(),
         nb::arg("GAMBET") = nb::none(),
         nb::arg("MASS") = nb::none(),
         nb::arg("CHARGE") = nb::none(),
         nb::arg("AG") = nb::none(),
         nb::arg("loc") = nb::none()
  )
      .def_prop_rw("DIR", &Fibre::DIR, &Fibre::set_DIR)
      .def_prop_rw("pos", &Fibre::pos, &Fibre::set_pos, "POSITION IN LAYOUT NEW STUFF....")
      .def_prop_rw("BETA0", &Fibre::BETA0, &Fibre::set_BETA0, ",P0C")
      .def_prop_rw("GAMMA0I", &Fibre::GAMMA0I, &Fibre::set_GAMMA0I, ",P0C")
      .def_prop_rw("GAMBET", &Fibre::GAMBET, &Fibre::set_GAMBET, ",P0C")
      .def_prop_rw("MASS", &Fibre::MASS, &Fibre::set_MASS, ",P0C")
      .def_prop_rw("CHARGE", &Fibre::CHARGE, &Fibre::set_CHARGE)
      .def_prop_rw("AG", &Fibre::AG, &Fibre::set_AG, "spin g-2 TO TIE LAYOUTS")
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
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D Fibre arrays are not used in structs/routines
  // 2D Fibre arrays are not used in structs/routines
  // 3D Fibre arrays are not used in structs/routines
}