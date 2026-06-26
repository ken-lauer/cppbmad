#include "pybmad/generated/structs_c.hpp"

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
// cartesian_map_struct
void init_cartesian_map_struct(nb::module_ &m, nb::class_<CartesianMapStruct> &cls) {
  cls.def(
         "__init__",
         [](CartesianMapStruct *self,
            std::optional<double> field_scale,
            std::optional<std::vector<double>> r0,
            std::optional<int> master_parameter,
            std::optional<int> ele_anchor_pt,
            std::optional<int> field_type,
            const CartesianMapTermStruct *ptr) {
           new (self) CartesianMapStruct(
               field_scale,
               r0,
               master_parameter,
               ele_anchor_pt,
               field_type,
               ptr_to_opt_ref(ptr)
           );
         },
         nb::arg("field_scale") = nb::none(),
         nb::arg("r0") = nb::none(),
         nb::arg("master_parameter") = nb::none(),
         nb::arg("ele_anchor_pt") = nb::none(),
         nb::arg("field_type") = nb::none(),
         nb::arg("ptr") = nb::none()
  )
      .def_prop_rw(
          "field_scale",
          &CartesianMapStruct::field_scale,
          &CartesianMapStruct::set_field_scale,
          "Factor to scale the fields by"
      )
      .def_prop_rw(
          "r0",
          &CartesianMapStruct::r0,
          &CartesianMapStruct::set_r0,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Field origin offset."
      )
      .def_prop_rw(
          "master_parameter",
          &CartesianMapStruct::master_parameter,
          &CartesianMapStruct::set_master_parameter,
          "Master parameter in ele%value(:) array to use for scaling the field."
      )
      .def_prop_rw(
          "ele_anchor_pt",
          &CartesianMapStruct::ele_anchor_pt,
          &CartesianMapStruct::set_ele_anchor_pt,
          "anchor_beginning$, anchor_center$, or anchor_end$"
      )
      .def_prop_rw(
          "field_type",
          &CartesianMapStruct::field_type,
          &CartesianMapStruct::set_field_type,
          "or electric$"
      )
      .def_prop_rw(
          "ptr",
          &CartesianMapStruct::ptr,
          &CartesianMapStruct::set_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return CartesianMapStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = CartesianMapStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const CartesianMapStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CartesianMapStruct &self) {
            return CartesianMapStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CartesianMapStruct &self, nb::dict &memo) { return CartesianMapStruct(self); }
      )
      .def(
          "__eq__",
          [](const CartesianMapStruct &self, const CartesianMapStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const CartesianMapStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<CartesianMapStructArray1D, CartesianMapStructAlloc1D>(
      m,
      "CartesianMapStructArray1D",
      "CartesianMapStructAlloc1D"
  );
  // 2D CartesianMapStruct arrays are not used in structs/routines
  // 3D CartesianMapStruct arrays are not used in structs/routines
}

// =============================================================================
// cartesian_map_term1_struct
void init_cartesian_map_term1_struct(nb::module_ &m, nb::class_<CartesianMapTerm1Struct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>,
             std::optional<int>>(),
         nb::arg("coef") = nb::none(),
         nb::arg("kx") = nb::none(),
         nb::arg("ky") = nb::none(),
         nb::arg("kz") = nb::none(),
         nb::arg("x0") = nb::none(),
         nb::arg("y0") = nb::none(),
         nb::arg("phi_z") = nb::none(),
         nb::arg("family") = nb::none(),
         nb::arg("form") = nb::none()
  )
      .def_prop_rw("coef", &CartesianMapTerm1Struct::coef, &CartesianMapTerm1Struct::set_coef)
      .def_prop_rw("kx", &CartesianMapTerm1Struct::kx, &CartesianMapTerm1Struct::set_kx)
      .def_prop_rw("ky", &CartesianMapTerm1Struct::ky, &CartesianMapTerm1Struct::set_ky)
      .def_prop_rw("kz", &CartesianMapTerm1Struct::kz, &CartesianMapTerm1Struct::set_kz)
      .def_prop_rw("x0", &CartesianMapTerm1Struct::x0, &CartesianMapTerm1Struct::set_x0)
      .def_prop_rw("y0", &CartesianMapTerm1Struct::y0, &CartesianMapTerm1Struct::set_y0)
      .def_prop_rw("phi_z", &CartesianMapTerm1Struct::phi_z, &CartesianMapTerm1Struct::set_phi_z)
      .def_prop_rw(
          "family",
          &CartesianMapTerm1Struct::family,
          &CartesianMapTerm1Struct::set_family,
          "family_x$, etc."
      )
      .def_prop_rw(
          "form",
          &CartesianMapTerm1Struct::form,
          &CartesianMapTerm1Struct::set_form,
          "hyper_y$, etc."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return CartesianMapTerm1StructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = CartesianMapTerm1StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const CartesianMapTerm1Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CartesianMapTerm1Struct &self) {
            return CartesianMapTerm1Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CartesianMapTerm1Struct &self, nb::dict &memo) {
            return CartesianMapTerm1Struct(self);
          }
      )
      .def(
          "__eq__",
          [](const CartesianMapTerm1Struct &self, const CartesianMapTerm1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const CartesianMapTerm1Struct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<CartesianMapTerm1StructArray1D, CartesianMapTerm1StructAlloc1D>(
      m,
      "CartesianMapTerm1StructArray1D",
      "CartesianMapTerm1StructAlloc1D"
  );
  // 2D CartesianMapTerm1Struct arrays are not used in structs/routines
  // 3D CartesianMapTerm1Struct arrays are not used in structs/routines
}

// =============================================================================
// cartesian_map_term_struct
void init_cartesian_map_term_struct(nb::module_ &m, nb::class_<CartesianMapTermStruct> &cls) {
  cls.def(
         nb::init<std::optional<std::string>, std::optional<int>>(),
         nb::arg("file") = nb::none(),
         nb::arg("n_link") = nb::none()
  )
      .def_prop_rw(
          "file",
          &CartesianMapTermStruct::file,
          &CartesianMapTermStruct::set_file,
          "Input file name. Used also as ID for instances."
      )
      .def_prop_rw(
          "n_link",
          &CartesianMapTermStruct::n_link,
          &CartesianMapTermStruct::set_n_link,
          "For memory management of %term"
      )
      .def_prop_ro("term", &CartesianMapTermStruct::term, nb::keep_alive<0, 1>())

      .def("__repr__", [](const CartesianMapTermStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CartesianMapTermStruct &self) {
            return CartesianMapTermStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CartesianMapTermStruct &self, nb::dict &memo) {
            return CartesianMapTermStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const CartesianMapTermStruct &self, const CartesianMapTermStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const CartesianMapTermStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D CartesianMapTermStruct arrays are not used in structs/routines
  // 2D CartesianMapTermStruct arrays are not used in structs/routines
  // 3D CartesianMapTermStruct arrays are not used in structs/routines
}

// =============================================================================
// complex_taylor_struct
void init_complex_taylor_struct(nb::module_ &m, nb::class_<ComplexTaylorStruct> &cls) {
  cls.def(nb::init<std::optional<std::complex<double>>>(), nb::arg("ref") = nb::none())
      .def_prop_rw("ref", &ComplexTaylorStruct::ref, &ComplexTaylorStruct::set_ref)
      .def_prop_ro("term", &ComplexTaylorStruct::term, nb::keep_alive<0, 1>())
      .def_static(
          "new_array1d",
          [](int sz) { return ComplexTaylorStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ComplexTaylorStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const ComplexTaylorStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ComplexTaylorStruct &self) {
            return ComplexTaylorStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ComplexTaylorStruct &self, nb::dict &memo) { return ComplexTaylorStruct(self); }
      )
      .def(
          "__eq__",
          [](const ComplexTaylorStruct &self, const ComplexTaylorStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ComplexTaylorStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<ComplexTaylorStructArray1D, ComplexTaylorStructAlloc1D>(
      m,
      "ComplexTaylorStructArray1D",
      "ComplexTaylorStructAlloc1D"
  );
  // 2D ComplexTaylorStruct arrays are not used in structs/routines
  // 3D ComplexTaylorStruct arrays are not used in structs/routines
}

// =============================================================================
// complex_taylor_term_struct
void init_complex_taylor_term_struct(nb::module_ &m, nb::class_<ComplexTaylorTermStruct> &cls) {
  cls.def(
         nb::init<std::optional<std::complex<double>>, std::optional<std::vector<int>>>(),
         nb::arg("coef") = nb::none(),
         nb::arg("expn") = nb::none()
  )
      .def_prop_rw("coef", &ComplexTaylorTermStruct::coef, &ComplexTaylorTermStruct::set_coef)
      .def_prop_rw(
          "expn",
          &ComplexTaylorTermStruct::expn,
          &ComplexTaylorTermStruct::set_expn,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return ComplexTaylorTermStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ComplexTaylorTermStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const ComplexTaylorTermStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ComplexTaylorTermStruct &self) {
            return ComplexTaylorTermStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ComplexTaylorTermStruct &self, nb::dict &memo) {
            return ComplexTaylorTermStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const ComplexTaylorTermStruct &self, const ComplexTaylorTermStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ComplexTaylorTermStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<ComplexTaylorTermStructArray1D, ComplexTaylorTermStructAlloc1D>(
      m,
      "ComplexTaylorTermStructArray1D",
      "ComplexTaylorTermStructAlloc1D"
  );
  // 2D ComplexTaylorTermStruct arrays are not used in structs/routines
  // 3D ComplexTaylorTermStruct arrays are not used in structs/routines
}

// =============================================================================
// control_ramp1_struct
void init_control_ramp1_struct(nb::module_ &m, nb::class_<ControlRamp1Struct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<double>>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<bool>>(),
         nb::arg("y_knot") = nb::none(),
         nb::arg("attribute") = nb::none(),
         nb::arg("slave_name") = nb::none(),
         nb::arg("is_controller") = nb::none()
  )
      .def_prop_rw(
          "y_knot",
          &ControlRamp1Struct::y_knot,
          &ControlRamp1Struct::set_y_knot,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro("stack", &ControlRamp1Struct::stack, nb::keep_alive<0, 1>(), "Evaluation stack")
      .def_prop_rw(
          "attribute",
          &ControlRamp1Struct::attribute,
          &ControlRamp1Struct::set_attribute,
          "Name of attribute controlled. Set to 'FIELD_OVERLAPS' for field overlaps."
      )
      .def_prop_rw(
          "slave_name",
          &ControlRamp1Struct::slave_name,
          &ControlRamp1Struct::set_slave_name,
          "Name of slave."
      )
      .def_prop_rw(
          "is_controller",
          &ControlRamp1Struct::is_controller,
          &ControlRamp1Struct::set_is_controller,
          "Is the slave a controller? If so bookkeeping is different."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return ControlRamp1StructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ControlRamp1StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const ControlRamp1Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ControlRamp1Struct &self) {
            return ControlRamp1Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ControlRamp1Struct &self, nb::dict &memo) { return ControlRamp1Struct(self); }
      )
      .def(
          "__eq__",
          [](const ControlRamp1Struct &self, const ControlRamp1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ControlRamp1Struct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<ControlRamp1StructArray1D, ControlRamp1StructAlloc1D>(
      m,
      "ControlRamp1StructArray1D",
      "ControlRamp1StructAlloc1D"
  );
  // 2D ControlRamp1Struct arrays are not used in structs/routines
  // 3D ControlRamp1Struct arrays are not used in structs/routines
}

// =============================================================================
// control_struct
void init_control_struct(nb::module_ &m, nb::class_<ControlStruct> &cls) {
  cls.def(
         "__init__",
         [](ControlStruct *self,
            std::optional<double> value,
            std::optional<std::vector<double>> y_knot,
            const LatEleLocStruct *slave,
            const LatEleLocStruct *lord,
            std::optional<std::string> slave_name,
            std::optional<std::string> attribute,
            std::optional<int> ix_attrib) {
           new (self) ControlStruct(
               value,
               y_knot,
               ptr_to_opt_ref(slave),
               ptr_to_opt_ref(lord),
               slave_name,
               attribute,
               ix_attrib
           );
         },
         nb::arg("value") = nb::none(),
         nb::arg("y_knot") = nb::none(),
         nb::arg("slave") = nb::none(),
         nb::arg("lord") = nb::none(),
         nb::arg("slave_name") = nb::none(),
         nb::arg("attribute") = nb::none(),
         nb::arg("ix_attrib") = nb::none()
  )
      .def_prop_rw(
          "value",
          &ControlStruct::value,
          &ControlStruct::set_value,
          "Used by group, and overlay elements."
      )
      .def_prop_rw(
          "y_knot",
          &ControlStruct::y_knot,
          &ControlStruct::set_y_knot,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro("stack", &ControlStruct::stack, nb::keep_alive<0, 1>(), "Evaluation stack")
      .def_prop_rw(
          "slave",
          &ControlStruct::slave,
          &ControlStruct::set_slave,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "lord",
          &ControlStruct::lord,
          &ControlStruct::set_lord,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "slave_name",
          &ControlStruct::slave_name,
          &ControlStruct::set_slave_name,
          "Name of slave."
      )
      .def_prop_rw(
          "attribute",
          &ControlStruct::attribute,
          &ControlStruct::set_attribute,
          "Name of attribute controlled. Set to 'FIELD_OVERLAPS' for field overlaps. Set to "
          "'INPUT' or 'OUTPUT' for feedback slaves."
      )
      .def_prop_rw(
          "ix_attrib",
          &ControlStruct::ix_attrib,
          &ControlStruct::set_ix_attrib,
          "Index of attribute controlled. See note above!"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return ControlStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ControlStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const ControlStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ControlStruct &self) {
            return ControlStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ControlStruct &self, nb::dict &memo) { return ControlStruct(self); }
      )
      .def(
          "__eq__",
          [](const ControlStruct &self, const ControlStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ControlStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<ControlStructArray1D, ControlStructAlloc1D>(
      m,
      "ControlStructArray1D",
      "ControlStructAlloc1D"
  );
  // 2D ControlStruct arrays are not used in structs/routines
  // 3D ControlStruct arrays are not used in structs/routines
}

// =============================================================================
// control_var1_struct
void init_control_var1_struct(nb::module_ &m, nb::class_<ControlVar1Struct> &cls) {
  cls.def(
         nb::init<std::optional<std::string>, std::optional<double>, std::optional<double>>(),
         nb::arg("name") = nb::none(),
         nb::arg("value") = nb::none(),
         nb::arg("old_value") = nb::none()
  )
      .def_prop_rw("name", &ControlVar1Struct::name, &ControlVar1Struct::set_name)
      .def_prop_rw("value", &ControlVar1Struct::value, &ControlVar1Struct::set_value)
      .def_prop_rw("old_value", &ControlVar1Struct::old_value, &ControlVar1Struct::set_old_value)
      .def_static(
          "new_array1d",
          [](int sz) { return ControlVar1StructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ControlVar1StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const ControlVar1Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ControlVar1Struct &self) {
            return ControlVar1Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ControlVar1Struct &self, nb::dict &memo) { return ControlVar1Struct(self); }
      )
      .def(
          "__eq__",
          [](const ControlVar1Struct &self, const ControlVar1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ControlVar1Struct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<ControlVar1StructArray1D, ControlVar1StructAlloc1D>(
      m,
      "ControlVar1StructArray1D",
      "ControlVar1StructAlloc1D"
  );
  // 2D ControlVar1Struct arrays are not used in structs/routines
  // 3D ControlVar1Struct arrays are not used in structs/routines
}

// =============================================================================
// controller_struct
void init_controller_struct(nb::module_ &m, nb::class_<ControllerStruct> &cls) {
  cls.def(nb::init<std::optional<std::vector<double>>>(), nb::arg("x_knot") = nb::none())
      .def_prop_ro("var", &ControllerStruct::var, nb::keep_alive<0, 1>())
      .def_prop_ro(
          "ramp",
          &ControllerStruct::ramp,
          nb::keep_alive<0, 1>(),
          "For ramper lord elements"
      )
      .def_prop_ro(
          "ramper_lord",
          &ControllerStruct::ramper_lord,
          nb::keep_alive<0, 1>(),
          "Ramper lord info for this slave"
      )
      .def_prop_rw(
          "x_knot",
          &ControllerStruct::x_knot,
          &ControllerStruct::set_x_knot,
          nb::for_getter(nb::keep_alive<0, 1>())
      )

      .def("__repr__", [](const ControllerStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ControllerStruct &self) {
            return ControllerStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ControllerStruct &self, nb::dict &memo) { return ControllerStruct(self); }
      )
      .def(
          "__eq__",
          [](const ControllerStruct &self, const ControllerStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ControllerStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D ControllerStruct arrays are not used in structs/routines
  // 2D ControllerStruct arrays are not used in structs/routines
  // 3D ControllerStruct arrays are not used in structs/routines
}

// =============================================================================
// coord_array_struct
void init_coord_array_struct(nb::module_ &m, nb::class_<CoordArrayStruct> &cls) {
  cls.def(nb::init<>())
      .def_prop_ro("orbit", &CoordArrayStruct::orbit, nb::keep_alive<0, 1>())
      .def_static(
          "new_array1d",
          [](int sz) { return CoordArrayStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = CoordArrayStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const CoordArrayStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CoordArrayStruct &self) {
            return CoordArrayStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CoordArrayStruct &self, nb::dict &memo) { return CoordArrayStruct(self); }
      )
      .def(
          "__eq__",
          [](const CoordArrayStruct &self, const CoordArrayStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const CoordArrayStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<CoordArrayStructArray1D, CoordArrayStructAlloc1D>(
      m,
      "CoordArrayStructArray1D",
      "CoordArrayStructAlloc1D"
  );
  // 2D CoordArrayStruct arrays are not used in structs/routines
  // 3D CoordArrayStruct arrays are not used in structs/routines
}

// =============================================================================
// coord_struct
void init_coord_struct(nb::module_ &m, nb::class_<CoordStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<double>>,
             std::optional<double>,
             std::optional<long double>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>>(),
         nb::arg("vec") = nb::none(),
         nb::arg("s") = nb::none(),
         nb::arg("t") = nb::none(),
         nb::arg("spin") = nb::none(),
         nb::arg("field") = nb::none(),
         nb::arg("phase") = nb::none(),
         nb::arg("charge") = nb::none(),
         nb::arg("dt_ref") = nb::none(),
         nb::arg("r") = nb::none(),
         nb::arg("p0c") = nb::none(),
         nb::arg("E_potential") = nb::none(),
         nb::arg("beta") = nb::none(),
         nb::arg("ix_ele") = nb::none(),
         nb::arg("ix_branch") = nb::none(),
         nb::arg("ix_turn") = nb::none(),
         nb::arg("ix_user") = nb::none(),
         nb::arg("state") = nb::none(),
         nb::arg("direction") = nb::none(),
         nb::arg("time_dir") = nb::none(),
         nb::arg("species") = nb::none(),
         nb::arg("location") = nb::none()
  )
      .def_prop_rw(
          "vec",
          &CoordStruct::vec,
          &CoordStruct::set_vec,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "(x, px, y, py, z, pz). Generally phase space for charged particles. See Bmad manual."
      )
      .def_prop_rw("s", &CoordStruct::s, &CoordStruct::set_s, "Longitudinal position")
      .def_prop_rw(
          "t",
          &CoordStruct::t,
          &CoordStruct::set_t,
          "Absolute time (not relative to reference). Note: Quad precision!"
      )
      .def_prop_rw(
          "spin",
          &CoordStruct::spin,
          &CoordStruct::set_spin,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Spin."
      )
      .def_prop_rw(
          "field",
          &CoordStruct::field,
          &CoordStruct::set_field,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Photon E-field intensity (x,y)."
      )
      .def_prop_rw(
          "phase",
          &CoordStruct::phase,
          &CoordStruct::set_phase,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Photon E-field phase (x,y). For charged particles, phase(1) is RF phase."
      )
      .def_prop_rw(
          "charge",
          &CoordStruct::charge,
          &CoordStruct::set_charge,
          "Macroparticle weight (which is different from particle species charge). For some space "
          "charge calcs the weight is in Coulombs."
      )
      .def_prop_rw(
          "dt_ref",
          &CoordStruct::dt_ref,
          &CoordStruct::set_dt_ref,
          "Used in: * time tracking for computing z. * by coherent photons = path_length/c_light."
      )
      .def_prop_rw("r", &CoordStruct::r, &CoordStruct::set_r, "For general use. Not used by Bmad.")
      .def_prop_rw(
          "p0c",
          &CoordStruct::p0c,
          &CoordStruct::set_p0c,
          "For non-photons: Reference momentum. For photons: Photon momentum (not reference)."
      )
      .def_prop_rw(
          "E_potential",
          &CoordStruct::E_potential,
          &CoordStruct::set_E_potential,
          "Potential energy."
      )
      .def_prop_rw("beta", &CoordStruct::beta, &CoordStruct::set_beta, "Velocity / c_light.")
      .def_prop_rw(
          "ix_ele",
          &CoordStruct::ix_ele,
          &CoordStruct::set_ix_ele,
          "Index of the lattice element the particle is in. May be -1 if element is not associated "
          "with a lattice."
      )
      .def_prop_rw(
          "ix_branch",
          &CoordStruct::ix_branch,
          &CoordStruct::set_ix_branch,
          "Index of the lattice branch the particle is in."
      )
      .def_prop_rw(
          "ix_turn",
          &CoordStruct::ix_turn,
          &CoordStruct::set_ix_turn,
          "Turn index for multiturn tracking."
      )
      .def_prop_rw(
          "ix_user",
          &CoordStruct::ix_user,
          &CoordStruct::set_ix_user,
          "For general use, not used by Bmad."
      )
      .def_prop_rw(
          "state",
          &CoordStruct::state,
          &CoordStruct::set_state,
          "alive$, lost$, lost_neg_x_aperture$, lost_pz$, etc."
      )
      .def_prop_rw(
          "direction",
          &CoordStruct::direction,
          &CoordStruct::set_direction,
          "+1 or -1. Sign of longitudinal direction of motion (ds/dt). This is independent of the "
          "element orientation."
      )
      .def_prop_rw(
          "time_dir",
          &CoordStruct::time_dir,
          &CoordStruct::set_time_dir,
          "+1 or -1. Time direction. -1 => Traveling backwards in time."
      )
      .def_prop_rw(
          "species",
          &CoordStruct::species,
          &CoordStruct::set_species,
          "positron$, proton$, etc."
      )
      .def_prop_rw(
          "location",
          &CoordStruct::location,
          &CoordStruct::set_location,
          "upstream_end$, inside$, or downstream_end$"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return CoordStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = CoordStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const CoordStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CoordStruct &self) {
            return CoordStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CoordStruct &self, nb::dict &memo) { return CoordStruct(self); }
      )
      .def(
          "__eq__",
          [](const CoordStruct &self, const CoordStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const CoordStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<CoordStructArray1D, CoordStructAlloc1D>(
      m,
      "CoordStructArray1D",
      "CoordStructAlloc1D"
  );
  // 2D CoordStruct arrays are not used in structs/routines
  // 3D CoordStruct arrays are not used in structs/routines
}

// =============================================================================
// cylindrical_map_struct
void init_cylindrical_map_struct(nb::module_ &m, nb::class_<CylindricalMapStruct> &cls) {
  cls.def(
         "__init__",
         [](CylindricalMapStruct *self,
            std::optional<int> m,
            std::optional<int> harmonic,
            std::optional<double> phi0_fieldmap,
            std::optional<double> theta0_azimuth,
            std::optional<double> field_scale,
            std::optional<int> master_parameter,
            std::optional<int> ele_anchor_pt,
            std::optional<double> dz,
            std::optional<std::vector<double>> r0,
            const CylindricalMapTermStruct *ptr) {
           new (self) CylindricalMapStruct(
               m,
               harmonic,
               phi0_fieldmap,
               theta0_azimuth,
               field_scale,
               master_parameter,
               ele_anchor_pt,
               dz,
               r0,
               ptr_to_opt_ref(ptr)
           );
         },
         nb::arg("m") = nb::none(),
         nb::arg("harmonic") = nb::none(),
         nb::arg("phi0_fieldmap") = nb::none(),
         nb::arg("theta0_azimuth") = nb::none(),
         nb::arg("field_scale") = nb::none(),
         nb::arg("master_parameter") = nb::none(),
         nb::arg("ele_anchor_pt") = nb::none(),
         nb::arg("dz") = nb::none(),
         nb::arg("r0") = nb::none(),
         nb::arg("ptr") = nb::none()
  )
      .def_prop_rw(
          "m",
          &CylindricalMapStruct::m,
          &CylindricalMapStruct::set_m,
          "Azimuthal Mode: varies as cos(m*phi - theta0_azimuth)"
      )
      .def_prop_rw(
          "harmonic",
          &CylindricalMapStruct::harmonic,
          &CylindricalMapStruct::set_harmonic,
          "Harmonic of fundamental"
      )
      .def_prop_rw(
          "phi0_fieldmap",
          &CylindricalMapStruct::phi0_fieldmap,
          &CylindricalMapStruct::set_phi0_fieldmap,
          "Mode oscillates as: twopi * (f * t + phi0_fieldmap)"
      )
      .def_prop_rw(
          "theta0_azimuth",
          &CylindricalMapStruct::theta0_azimuth,
          &CylindricalMapStruct::set_theta0_azimuth,
          "Azimuthal ((x, y) plane) orientation of mode."
      )
      .def_prop_rw(
          "field_scale",
          &CylindricalMapStruct::field_scale,
          &CylindricalMapStruct::set_field_scale,
          "Factor to scale the fields by"
      )
      .def_prop_rw(
          "master_parameter",
          &CylindricalMapStruct::master_parameter,
          &CylindricalMapStruct::set_master_parameter,
          "Master parameter in ele%value(:) array to use for scaling the field."
      )
      .def_prop_rw(
          "ele_anchor_pt",
          &CylindricalMapStruct::ele_anchor_pt,
          &CylindricalMapStruct::set_ele_anchor_pt,
          "anchor_beginning$, anchor_center$, or anchor_end$"
      )
      .def_prop_rw(
          "dz",
          &CylindricalMapStruct::dz,
          &CylindricalMapStruct::set_dz,
          "Distance between sampled field points."
      )
      .def_prop_rw(
          "r0",
          &CylindricalMapStruct::r0,
          &CylindricalMapStruct::set_r0,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Field origin offset."
      )
      .def_prop_rw(
          "ptr",
          &CylindricalMapStruct::ptr,
          &CylindricalMapStruct::set_ptr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return CylindricalMapStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = CylindricalMapStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const CylindricalMapStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CylindricalMapStruct &self) {
            return CylindricalMapStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CylindricalMapStruct &self, nb::dict &memo) {
            return CylindricalMapStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const CylindricalMapStruct &self, const CylindricalMapStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const CylindricalMapStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<CylindricalMapStructArray1D, CylindricalMapStructAlloc1D>(
      m,
      "CylindricalMapStructArray1D",
      "CylindricalMapStructAlloc1D"
  );
  // 2D CylindricalMapStruct arrays are not used in structs/routines
  // 3D CylindricalMapStruct arrays are not used in structs/routines
}

// =============================================================================
// cylindrical_map_term1_struct
void init_cylindrical_map_term1_struct(nb::module_ &m, nb::class_<CylindricalMapTerm1Struct> &cls) {
  cls.def(
         nb::init<std::optional<std::complex<double>>, std::optional<std::complex<double>>>(),
         nb::arg("e_coef") = nb::none(),
         nb::arg("b_coef") = nb::none()
  )
      .def_prop_rw(
          "e_coef",
          &CylindricalMapTerm1Struct::e_coef,
          &CylindricalMapTerm1Struct::set_e_coef
      )
      .def_prop_rw(
          "b_coef",
          &CylindricalMapTerm1Struct::b_coef,
          &CylindricalMapTerm1Struct::set_b_coef
      )
      .def_static(
          "new_array1d",
          [](int sz) { return CylindricalMapTerm1StructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = CylindricalMapTerm1StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const CylindricalMapTerm1Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CylindricalMapTerm1Struct &self) {
            return CylindricalMapTerm1Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CylindricalMapTerm1Struct &self, nb::dict &memo) {
            return CylindricalMapTerm1Struct(self);
          }
      )
      .def(
          "__eq__",
          [](const CylindricalMapTerm1Struct &self, const CylindricalMapTerm1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const CylindricalMapTerm1Struct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<CylindricalMapTerm1StructArray1D, CylindricalMapTerm1StructAlloc1D>(
      m,
      "CylindricalMapTerm1StructArray1D",
      "CylindricalMapTerm1StructAlloc1D"
  );
  // 2D CylindricalMapTerm1Struct arrays are not used in structs/routines
  // 3D CylindricalMapTerm1Struct arrays are not used in structs/routines
}

// =============================================================================
// cylindrical_map_term_struct
void init_cylindrical_map_term_struct(nb::module_ &m, nb::class_<CylindricalMapTermStruct> &cls) {
  cls.def(
         nb::init<std::optional<std::string>, std::optional<int>>(),
         nb::arg("file") = nb::none(),
         nb::arg("n_link") = nb::none()
  )
      .def_prop_rw(
          "file",
          &CylindricalMapTermStruct::file,
          &CylindricalMapTermStruct::set_file,
          "Input file name. Used also as ID for instances."
      )
      .def_prop_rw(
          "n_link",
          &CylindricalMapTermStruct::n_link,
          &CylindricalMapTermStruct::set_n_link,
          "For memory management of this structure"
      )
      .def_prop_ro("term", &CylindricalMapTermStruct::term, nb::keep_alive<0, 1>())

      .def("__repr__", [](const CylindricalMapTermStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CylindricalMapTermStruct &self) {
            return CylindricalMapTermStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CylindricalMapTermStruct &self, nb::dict &memo) {
            return CylindricalMapTermStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const CylindricalMapTermStruct &self, const CylindricalMapTermStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const CylindricalMapTermStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D CylindricalMapTermStruct arrays are not used in structs/routines
  // 2D CylindricalMapTermStruct arrays are not used in structs/routines
  // 3D CylindricalMapTermStruct arrays are not used in structs/routines
}