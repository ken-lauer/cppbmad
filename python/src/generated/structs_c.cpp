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
// converter_dir_1D_struct
void init_converter_dir_1D_struct(nb::module_ &m, nb::class_<ConverterDir1dStruct> &cls) {
  cls.def(
         nb::init<std::optional<double>, std::optional<std::vector<double>>>(),
         nb::arg("pc_out") = nb::none(),
         nb::arg("poly") = nb::none()
  )
      .def_prop_rw(
          "pc_out",
          &ConverterDir1dStruct::pc_out,
          &ConverterDir1dStruct::set_pc_out,
          "pc_out value at fit"
      )
      .def_prop_rw(
          "poly",
          &ConverterDir1dStruct::poly,
          &ConverterDir1dStruct::set_poly,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "param(r) = Sum: poly(i) * r^i"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return ConverterDir1dStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ConverterDir1dStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const ConverterDir1dStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ConverterDir1dStruct &self) {
            return ConverterDir1dStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ConverterDir1dStruct &self, nb::dict &memo) {
            return ConverterDir1dStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const ConverterDir1dStruct &self, const ConverterDir1dStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ConverterDir1dStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<ConverterDir1dStructArray1D, ConverterDir1dStructAlloc1D>(
      m,
      "ConverterDir1DStructArray1D",
      "ConverterDir1DStructAlloc1D"
  );
  // 2D ConverterDir1dStruct arrays are not used in structs/routines
  // 3D ConverterDir1dStruct arrays are not used in structs/routines
}

// =============================================================================
// converter_dir_2D_struct
void init_converter_dir_2D_struct(nb::module_ &m, nb::class_<ConverterDir2dStruct> &cls) {
  cls.def(
         nb::init<std::optional<double>, std::optional<std::vector<double>>>(),
         nb::arg("k") = nb::none(),
         nb::arg("poly") = nb::none()
  )
      .def_prop_rw("k", &ConverterDir2dStruct::k, &ConverterDir2dStruct::set_k)
      .def_prop_rw(
          "poly",
          &ConverterDir2dStruct::poly,
          &ConverterDir2dStruct::set_poly,
          nb::for_getter(nb::keep_alive<0, 1>())
      )

      .def("__repr__", [](const ConverterDir2dStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ConverterDir2dStruct &self) {
            return ConverterDir2dStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ConverterDir2dStruct &self, nb::dict &memo) {
            return ConverterDir2dStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const ConverterDir2dStruct &self, const ConverterDir2dStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ConverterDir2dStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D ConverterDir2dStruct arrays are not used in structs/routines
  // 2D ConverterDir2dStruct arrays are not used in structs/routines
  // 3D ConverterDir2dStruct arrays are not used in structs/routines
}

// =============================================================================
// converter_dir_coef_struct
void init_converter_dir_coef_struct(nb::module_ &m, nb::class_<ConverterDirCoefStruct> &cls) {
  cls.def(
         "__init__",
         [](ConverterDirCoefStruct *self,
            const ConverterDir2dStruct *fit_2d_r,
            const ConverterDir2dStruct *fit_2d_pc,
            std::optional<double> c0) {
           new (self)
               ConverterDirCoefStruct(ptr_to_opt_ref(fit_2d_r), ptr_to_opt_ref(fit_2d_pc), c0);
         },
         nb::arg("fit_2d_r") = nb::none(),
         nb::arg("fit_2d_pc") = nb::none(),
         nb::arg("c0") = nb::none()
  )
      .def_prop_ro("fit_1d_r", &ConverterDirCoefStruct::fit_1d_r, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "fit_2d_r",
          &ConverterDirCoefStruct::fit_2d_r,
          &ConverterDirCoefStruct::set_fit_2d_r,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "fit_2d_pc",
          &ConverterDirCoefStruct::fit_2d_pc,
          &ConverterDirCoefStruct::set_fit_2d_pc,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw("c0", &ConverterDirCoefStruct::c0, &ConverterDirCoefStruct::set_c0)

      .def("__repr__", [](const ConverterDirCoefStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ConverterDirCoefStruct &self) {
            return ConverterDirCoefStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ConverterDirCoefStruct &self, nb::dict &memo) {
            return ConverterDirCoefStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const ConverterDirCoefStruct &self, const ConverterDirCoefStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ConverterDirCoefStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D ConverterDirCoefStruct arrays are not used in structs/routines
  // 2D ConverterDirCoefStruct arrays are not used in structs/routines
  // 3D ConverterDirCoefStruct arrays are not used in structs/routines
}

// =============================================================================
// converter_direction_out_struct
void init_converter_direction_out_struct(
    nb::module_ &m,
    nb::class_<ConverterDirectionOutStruct> &cls
) {
  cls.def(
         "__init__",
         [](ConverterDirectionOutStruct *self,
            const ConverterDirCoefStruct *beta,
            const ConverterDirCoefStruct *alpha_x,
            const ConverterDirCoefStruct *alpha_y,
            const ConverterDirCoefStruct *dxds_min,
            const ConverterDirCoefStruct *dxds_max,
            const ConverterDirCoefStruct *dyds_max,
            const ConverterDirCoefStruct *c_x) {
           new (self) ConverterDirectionOutStruct(
               ptr_to_opt_ref(beta),
               ptr_to_opt_ref(alpha_x),
               ptr_to_opt_ref(alpha_y),
               ptr_to_opt_ref(dxds_min),
               ptr_to_opt_ref(dxds_max),
               ptr_to_opt_ref(dyds_max),
               ptr_to_opt_ref(c_x)
           );
         },
         nb::arg("beta") = nb::none(),
         nb::arg("alpha_x") = nb::none(),
         nb::arg("alpha_y") = nb::none(),
         nb::arg("dxds_min") = nb::none(),
         nb::arg("dxds_max") = nb::none(),
         nb::arg("dyds_max") = nb::none(),
         nb::arg("c_x") = nb::none()
  )
      .def_prop_rw(
          "beta",
          &ConverterDirectionOutStruct::beta,
          &ConverterDirectionOutStruct::set_beta,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "alpha_x",
          &ConverterDirectionOutStruct::alpha_x,
          &ConverterDirectionOutStruct::set_alpha_x,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "alpha_y",
          &ConverterDirectionOutStruct::alpha_y,
          &ConverterDirectionOutStruct::set_alpha_y,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "dxds_min",
          &ConverterDirectionOutStruct::dxds_min,
          &ConverterDirectionOutStruct::set_dxds_min,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "dxds_max",
          &ConverterDirectionOutStruct::dxds_max,
          &ConverterDirectionOutStruct::set_dxds_max,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "dyds_max",
          &ConverterDirectionOutStruct::dyds_max,
          &ConverterDirectionOutStruct::set_dyds_max,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "c_x",
          &ConverterDirectionOutStruct::c_x,
          &ConverterDirectionOutStruct::set_c_x,
          nb::for_getter(nb::keep_alive<0, 1>())
      )

      .def("__repr__", [](const ConverterDirectionOutStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ConverterDirectionOutStruct &self) {
            return ConverterDirectionOutStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ConverterDirectionOutStruct &self, nb::dict &memo) {
            return ConverterDirectionOutStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const ConverterDirectionOutStruct &self, const ConverterDirectionOutStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ConverterDirectionOutStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D ConverterDirectionOutStruct arrays are not used in structs/routines
  // 2D ConverterDirectionOutStruct arrays are not used in structs/routines
  // 3D ConverterDirectionOutStruct arrays are not used in structs/routines
}

// =============================================================================
// converter_distribution_struct
void init_converter_distribution_struct(
    nb::module_ &m,
    nb::class_<ConverterDistributionStruct> &cls
) {
  cls.def(nb::init<std::optional<double>>(), nb::arg("thickness") = nb::none())
      .def_prop_rw(
          "thickness",
          &ConverterDistributionStruct::thickness,
          &ConverterDistributionStruct::set_thickness
      )
      .def_prop_ro(
          "sub_dist",
          &ConverterDistributionStruct::sub_dist,
          nb::keep_alive<0, 1>(),
          "Distribution at various pc_in values."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return ConverterDistributionStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ConverterDistributionStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const ConverterDistributionStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ConverterDistributionStruct &self) {
            return ConverterDistributionStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ConverterDistributionStruct &self, nb::dict &memo) {
            return ConverterDistributionStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const ConverterDistributionStruct &self, const ConverterDistributionStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ConverterDistributionStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<ConverterDistributionStructArray1D, ConverterDistributionStructAlloc1D>(
      m,
      "ConverterDistributionStructArray1D",
      "ConverterDistributionStructAlloc1D"
  );
  // 2D ConverterDistributionStruct arrays are not used in structs/routines
  // 3D ConverterDistributionStruct arrays are not used in structs/routines
}

// =============================================================================
// converter_prob_pc_r_struct
void init_converter_prob_pc_r_struct(nb::module_ &m, nb::class_<ConverterProbPcRStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<double>>>(),
         nb::arg("pc_out") = nb::none(),
         nb::arg("r") = nb::none(),
         nb::arg("prob") = nb::none(),
         nb::arg("spin_z") = nb::none(),
         nb::arg("pc_out_min") = nb::none(),
         nb::arg("pc_out_max") = nb::none(),
         nb::arg("integrated_prob") = nb::none(),
         nb::arg("p_norm") = nb::none(),
         nb::arg("integ_pc_out") = nb::none(),
         nb::arg("integ_r") = nb::none(),
         nb::arg("integ_r_ave") = nb::none()
  )
      .def_prop_rw(
          "pc_out",
          &ConverterProbPcRStruct::pc_out,
          &ConverterProbPcRStruct::set_pc_out,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Grid pc_out values."
      )
      .def_prop_rw(
          "r",
          &ConverterProbPcRStruct::r,
          &ConverterProbPcRStruct::set_r,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Grid r_out values."
      )
      .def_prop_rw(
          "prob",
          &ConverterProbPcRStruct::prob,
          &ConverterProbPcRStruct::set_prob,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Probability grid."
      )
      .def_prop_rw(
          "spin_z",
          &ConverterProbPcRStruct::spin_z,
          &ConverterProbPcRStruct::set_spin_z,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Z polarization grid. Stuff below is calculated rather than read in from the lattice "
          "file."
      )
      .def_prop_rw(
          "pc_out_min",
          &ConverterProbPcRStruct::pc_out_min,
          &ConverterProbPcRStruct::set_pc_out_min
      )
      .def_prop_rw(
          "pc_out_max",
          &ConverterProbPcRStruct::pc_out_max,
          &ConverterProbPcRStruct::set_pc_out_max
      )
      .def_prop_rw(
          "integrated_prob",
          &ConverterProbPcRStruct::integrated_prob,
          &ConverterProbPcRStruct::set_integrated_prob,
          "Integrated probability over (pc_out, r) with restrictions factered in."
      )
      .def_prop_rw(
          "p_norm",
          &ConverterProbPcRStruct::p_norm,
          &ConverterProbPcRStruct::set_p_norm,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Normalized probability taking into account. angle_out_max, pc_out_min, and pc_out_max "
          "restrictions."
      )
      .def_prop_rw(
          "integ_pc_out",
          &ConverterProbPcRStruct::integ_pc_out,
          &ConverterProbPcRStruct::set_integ_pc_out,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Normalized probability integrated from min pc_out up."
      )
      .def_prop_rw(
          "integ_r",
          &ConverterProbPcRStruct::integ_r,
          &ConverterProbPcRStruct::set_integ_r,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "integ_r_ave",
          &ConverterProbPcRStruct::integ_r_ave,
          &ConverterProbPcRStruct::set_integ_r_ave,
          nb::for_getter(nb::keep_alive<0, 1>())
      )

      .def("__repr__", [](const ConverterProbPcRStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ConverterProbPcRStruct &self) {
            return ConverterProbPcRStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ConverterProbPcRStruct &self, nb::dict &memo) {
            return ConverterProbPcRStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const ConverterProbPcRStruct &self, const ConverterProbPcRStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ConverterProbPcRStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D ConverterProbPcRStruct arrays are not used in structs/routines
  // 2D ConverterProbPcRStruct arrays are not used in structs/routines
  // 3D ConverterProbPcRStruct arrays are not used in structs/routines
}

// =============================================================================
// converter_struct
void init_converter_struct(nb::module_ &m, nb::class_<ConverterStruct> &cls) {
  cls.def(
         nb::init<std::optional<int>, std::optional<std::string>>(),
         nb::arg("species_out") = nb::none(),
         nb::arg("material_type") = nb::none()
  )
      .def_prop_rw(
          "species_out",
          &ConverterStruct::species_out,
          &ConverterStruct::set_species_out,
          "Output species"
      )
      .def_prop_rw(
          "material_type",
          &ConverterStruct::material_type,
          &ConverterStruct::set_material_type
      )
      .def_prop_ro(
          "dist",
          &ConverterStruct::dist,
          nb::keep_alive<0, 1>(),
          "Distribution at various thicknesses"
      )

      .def("__repr__", [](const ConverterStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ConverterStruct &self) {
            return ConverterStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ConverterStruct &self, nb::dict &memo) { return ConverterStruct(self); }
      )
      .def(
          "__eq__",
          [](const ConverterStruct &self, const ConverterStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ConverterStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D ConverterStruct arrays are not used in structs/routines
  // 2D ConverterStruct arrays are not used in structs/routines
  // 3D ConverterStruct arrays are not used in structs/routines
}

// =============================================================================
// converter_sub_distribution_struct
void init_converter_sub_distribution_struct(
    nb::module_ &m,
    nb::class_<ConverterSubDistributionStruct> &cls
) {
  cls.def(
         "__init__",
         [](ConverterSubDistributionStruct *self,
            std::optional<double> pc_in,
            std::optional<std::vector<double>> spin_in,
            const ConverterProbPcRStruct *prob_pc_r,
            const ConverterDirectionOutStruct *dir_out) {
           new (self) ConverterSubDistributionStruct(
               pc_in,
               spin_in,
               ptr_to_opt_ref(prob_pc_r),
               ptr_to_opt_ref(dir_out)
           );
         },
         nb::arg("pc_in") = nb::none(),
         nb::arg("spin_in") = nb::none(),
         nb::arg("prob_pc_r") = nb::none(),
         nb::arg("dir_out") = nb::none()
  )
      .def_prop_rw(
          "pc_in",
          &ConverterSubDistributionStruct::pc_in,
          &ConverterSubDistributionStruct::set_pc_in
      )
      .def_prop_rw(
          "spin_in",
          &ConverterSubDistributionStruct::spin_in,
          &ConverterSubDistributionStruct::set_spin_in,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "prob_pc_r",
          &ConverterSubDistributionStruct::prob_pc_r,
          &ConverterSubDistributionStruct::set_prob_pc_r,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "dir_out",
          &ConverterSubDistributionStruct::dir_out,
          &ConverterSubDistributionStruct::set_dir_out,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return ConverterSubDistributionStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ConverterSubDistributionStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const ConverterSubDistributionStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ConverterSubDistributionStruct &self) {
            return ConverterSubDistributionStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ConverterSubDistributionStruct &self, nb::dict &memo) {
            return ConverterSubDistributionStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const ConverterSubDistributionStruct &self, const ConverterSubDistributionStruct &other
          ) { return self.get_fortran_ptr() == other.get_fortran_ptr(); },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ConverterSubDistributionStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<
      ConverterSubDistributionStructArray1D,
      ConverterSubDistributionStructAlloc1D>(
      m,
      "ConverterSubDistributionStructArray1D",
      "ConverterSubDistributionStructAlloc1D"
  );
  // 2D ConverterSubDistributionStruct arrays are not used in structs/routines
  // 3D ConverterSubDistributionStruct arrays are not used in structs/routines
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

// =============================================================================
// csr_bunch_slice_struct
void init_csr_bunch_slice_struct(nb::module_ &m, nb::class_<CsrBunchSliceStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("x0") = nb::none(),
         nb::arg("y0") = nb::none(),
         nb::arg("z0_edge") = nb::none(),
         nb::arg("z1_edge") = nb::none(),
         nb::arg("z_center") = nb::none(),
         nb::arg("sig_x") = nb::none(),
         nb::arg("sig_y") = nb::none(),
         nb::arg("charge") = nb::none(),
         nb::arg("dcharge_density_dz") = nb::none(),
         nb::arg("edge_dcharge_density_dz") = nb::none(),
         nb::arg("kick_csr") = nb::none(),
         nb::arg("coef_lsc_plus") = nb::none(),
         nb::arg("coef_lsc_minus") = nb::none(),
         nb::arg("kick_lsc") = nb::none(),
         nb::arg("n_particle") = nb::none()
  )
      .def_prop_rw(
          "x0",
          &CsrBunchSliceStruct::x0,
          &CsrBunchSliceStruct::set_x0,
          "Transverse center of the particle distrubution"
      )
      .def_prop_rw(
          "y0",
          &CsrBunchSliceStruct::y0,
          &CsrBunchSliceStruct::set_y0,
          "Transverse center of the particle distrubution"
      )
      .def_prop_rw(
          "z0_edge",
          &CsrBunchSliceStruct::z0_edge,
          &CsrBunchSliceStruct::set_z0_edge,
          "Left (min z) edge of bin"
      )
      .def_prop_rw(
          "z1_edge",
          &CsrBunchSliceStruct::z1_edge,
          &CsrBunchSliceStruct::set_z1_edge,
          "Right (max z) edge of bin"
      )
      .def_prop_rw(
          "z_center",
          &CsrBunchSliceStruct::z_center,
          &CsrBunchSliceStruct::set_z_center,
          "z at center of bin."
      )
      .def_prop_rw(
          "sig_x",
          &CsrBunchSliceStruct::sig_x,
          &CsrBunchSliceStruct::set_sig_x,
          "particle's RMS width"
      )
      .def_prop_rw(
          "sig_y",
          &CsrBunchSliceStruct::sig_y,
          &CsrBunchSliceStruct::set_sig_y,
          "particle's RMS width"
      )
      .def_prop_rw(
          "charge",
          &CsrBunchSliceStruct::charge,
          &CsrBunchSliceStruct::set_charge,
          "charge of the particles"
      )
      .def_prop_rw(
          "dcharge_density_dz",
          &CsrBunchSliceStruct::dcharge_density_dz,
          &CsrBunchSliceStruct::set_dcharge_density_dz,
          "Charge density gradient"
      )
      .def_prop_rw(
          "edge_dcharge_density_dz",
          &CsrBunchSliceStruct::edge_dcharge_density_dz,
          &CsrBunchSliceStruct::set_edge_dcharge_density_dz,
          "gradient between this and preceeding bin. [Evaluated at bin edge.]"
      )
      .def_prop_rw(
          "kick_csr",
          &CsrBunchSliceStruct::kick_csr,
          &CsrBunchSliceStruct::set_kick_csr,
          "CSR kick"
      )
      .def_prop_rw(
          "coef_lsc_plus",
          &CsrBunchSliceStruct::coef_lsc_plus,
          &CsrBunchSliceStruct::set_coef_lsc_plus,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "LSC Kick coefs."
      )
      .def_prop_rw(
          "coef_lsc_minus",
          &CsrBunchSliceStruct::coef_lsc_minus,
          &CsrBunchSliceStruct::set_coef_lsc_minus,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "LSC Kick coefs."
      )
      .def_prop_rw("kick_lsc", &CsrBunchSliceStruct::kick_lsc, &CsrBunchSliceStruct::set_kick_lsc)
      .def_prop_rw(
          "n_particle",
          &CsrBunchSliceStruct::n_particle,
          &CsrBunchSliceStruct::set_n_particle,
          "Number of particles in slice can be a fraction since particles span multiple bins."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return CsrBunchSliceStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = CsrBunchSliceStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const CsrBunchSliceStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CsrBunchSliceStruct &self) {
            return CsrBunchSliceStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CsrBunchSliceStruct &self, nb::dict &memo) { return CsrBunchSliceStruct(self); }
      )
      .def(
          "__eq__",
          [](const CsrBunchSliceStruct &self, const CsrBunchSliceStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const CsrBunchSliceStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<CsrBunchSliceStructArray1D, CsrBunchSliceStructAlloc1D>(
      m,
      "CsrBunchSliceStructArray1D",
      "CsrBunchSliceStructAlloc1D"
  );
  // 2D CsrBunchSliceStruct arrays are not used in structs/routines
  // 3D CsrBunchSliceStruct arrays are not used in structs/routines
}

// =============================================================================
// csr_ele_info_struct
void init_csr_ele_info_struct(nb::module_ &m, nb::class_<CsrEleInfoStruct> &cls) {
  cls.def(
         "__init__",
         [](CsrEleInfoStruct *self,
            const EleStruct *ele,
            const CoordStruct *orbit0,
            const CoordStruct *orbit1,
            const FloorPositionStruct *floor0,
            const FloorPositionStruct *floor1,
            const FloorPositionStruct *ref_floor0,
            const FloorPositionStruct *ref_floor1,
            const SplineStruct *spline,
            std::optional<double> theta_chord,
            std::optional<double> L_chord,
            std::optional<double> dL_s) {
           new (self) CsrEleInfoStruct(
               ptr_to_opt_ref(ele),
               ptr_to_opt_ref(orbit0),
               ptr_to_opt_ref(orbit1),
               ptr_to_opt_ref(floor0),
               ptr_to_opt_ref(floor1),
               ptr_to_opt_ref(ref_floor0),
               ptr_to_opt_ref(ref_floor1),
               ptr_to_opt_ref(spline),
               theta_chord,
               L_chord,
               dL_s
           );
         },
         nb::arg("ele") = nb::none(),
         nb::arg("orbit0") = nb::none(),
         nb::arg("orbit1") = nb::none(),
         nb::arg("floor0") = nb::none(),
         nb::arg("floor1") = nb::none(),
         nb::arg("ref_floor0") = nb::none(),
         nb::arg("ref_floor1") = nb::none(),
         nb::arg("spline") = nb::none(),
         nb::arg("theta_chord") = nb::none(),
         nb::arg("L_chord") = nb::none(),
         nb::arg("dL_s") = nb::none()
  )
      .def_prop_rw(
          "ele",
          &CsrEleInfoStruct::ele,
          &CsrEleInfoStruct::set_ele,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "lattice element"
      )
      .def_prop_rw(
          "orbit0",
          &CsrEleInfoStruct::orbit0,
          &CsrEleInfoStruct::set_orbit0,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "centroid orbit at entrance/exit ends"
      )
      .def_prop_rw(
          "orbit1",
          &CsrEleInfoStruct::orbit1,
          &CsrEleInfoStruct::set_orbit1,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "centroid orbit at entrance/exit ends"
      )
      .def_prop_rw(
          "floor0",
          &CsrEleInfoStruct::floor0,
          &CsrEleInfoStruct::set_floor0,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Floor position of centroid at entrance/exit ends"
      )
      .def_prop_rw(
          "floor1",
          &CsrEleInfoStruct::floor1,
          &CsrEleInfoStruct::set_floor1,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Floor position of centroid at entrance/exit ends"
      )
      .def_prop_rw(
          "ref_floor0",
          &CsrEleInfoStruct::ref_floor0,
          &CsrEleInfoStruct::set_ref_floor0,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Floor position of element ref coords at entrance/exit ends"
      )
      .def_prop_rw(
          "ref_floor1",
          &CsrEleInfoStruct::ref_floor1,
          &CsrEleInfoStruct::set_ref_floor1,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Floor position of element ref coords at entrance/exit ends"
      )
      .def_prop_rw(
          "spline",
          &CsrEleInfoStruct::spline,
          &CsrEleInfoStruct::set_spline,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Spline for centroid orbit. spline%x = distance along chord. The spline is zero at the "
          "ends by construction."
      )
      .def_prop_rw(
          "theta_chord",
          &CsrEleInfoStruct::theta_chord,
          &CsrEleInfoStruct::set_theta_chord,
          "Reference angle of chord in z-x plane"
      )
      .def_prop_rw(
          "L_chord",
          &CsrEleInfoStruct::L_chord,
          &CsrEleInfoStruct::set_L_chord,
          "Chord Length. Negative if bunch moves backwards in element."
      )
      .def_prop_rw(
          "dL_s",
          &CsrEleInfoStruct::dL_s,
          &CsrEleInfoStruct::set_dL_s,
          "L_s(of element) - L_chord"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return CsrEleInfoStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = CsrEleInfoStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const CsrEleInfoStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CsrEleInfoStruct &self) {
            return CsrEleInfoStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CsrEleInfoStruct &self, nb::dict &memo) { return CsrEleInfoStruct(self); }
      )
      .def(
          "__eq__",
          [](const CsrEleInfoStruct &self, const CsrEleInfoStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const CsrEleInfoStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<CsrEleInfoStructArray1D, CsrEleInfoStructAlloc1D>(
      m,
      "CsrEleInfoStructArray1D",
      "CsrEleInfoStructAlloc1D"
  );
  // 2D CsrEleInfoStruct arrays are not used in structs/routines
  // 3D CsrEleInfoStruct arrays are not used in structs/routines
}

// =============================================================================
// csr_kick1_struct
void init_csr_kick1_struct(nb::module_ &m, nb::class_<CsrKick1Struct> &cls) {
  cls.def(
         "__init__",
         [](CsrKick1Struct *self,
            std::optional<double> I_csr,
            std::optional<double> I_int_csr,
            std::optional<double> image_kick_csr,
            std::optional<std::vector<double>> L_vec,
            std::optional<double> L,
            std::optional<double> dL,
            std::optional<double> dz_particles,
            std::optional<double> s_chord_source,
            std::optional<double> theta_L,
            std::optional<double> theta_sl,
            std::optional<double> theta_lk,
            std::optional<int> ix_ele_source,
            const FloorPositionStruct *floor_s) {
           new (self) CsrKick1Struct(
               I_csr,
               I_int_csr,
               image_kick_csr,
               L_vec,
               L,
               dL,
               dz_particles,
               s_chord_source,
               theta_L,
               theta_sl,
               theta_lk,
               ix_ele_source,
               ptr_to_opt_ref(floor_s)
           );
         },
         nb::arg("I_csr") = nb::none(),
         nb::arg("I_int_csr") = nb::none(),
         nb::arg("image_kick_csr") = nb::none(),
         nb::arg("L_vec") = nb::none(),
         nb::arg("L") = nb::none(),
         nb::arg("dL") = nb::none(),
         nb::arg("dz_particles") = nb::none(),
         nb::arg("s_chord_source") = nb::none(),
         nb::arg("theta_L") = nb::none(),
         nb::arg("theta_sl") = nb::none(),
         nb::arg("theta_lk") = nb::none(),
         nb::arg("ix_ele_source") = nb::none(),
         nb::arg("floor_s") = nb::none()
  )
      .def_prop_rw("I_csr", &CsrKick1Struct::I_csr, &CsrKick1Struct::set_I_csr, "Kick integral.")
      .def_prop_rw(
          "I_int_csr",
          &CsrKick1Struct::I_int_csr,
          &CsrKick1Struct::set_I_int_csr,
          "Integrated Kick integral."
      )
      .def_prop_rw(
          "image_kick_csr",
          &CsrKick1Struct::image_kick_csr,
          &CsrKick1Struct::set_image_kick_csr,
          "kick."
      )
      .def_prop_rw(
          "L_vec",
          &CsrKick1Struct::L_vec,
          &CsrKick1Struct::set_L_vec,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "L vector in global coordinates."
      )
      .def_prop_rw(
          "L",
          &CsrKick1Struct::L,
          &CsrKick1Struct::set_L,
          "Distance between source and kick points."
      )
      .def_prop_rw("dL", &CsrKick1Struct::dL, &CsrKick1Struct::set_dL, "= epsilon_L = Ls - L")
      .def_prop_rw(
          "dz_particles",
          &CsrKick1Struct::dz_particles,
          &CsrKick1Struct::set_dz_particles,
          "Kicked particle - source particle position at constant time."
      )
      .def_prop_rw(
          "s_chord_source",
          &CsrKick1Struct::s_chord_source,
          &CsrKick1Struct::set_s_chord_source,
          "Source point coordinate along chord."
      )
      .def_prop_rw(
          "theta_L",
          &CsrKick1Struct::theta_L,
          &CsrKick1Struct::set_theta_L,
          "Angle of L vector"
      )
      .def_prop_rw(
          "theta_sl",
          &CsrKick1Struct::theta_sl,
          &CsrKick1Struct::set_theta_sl,
          "Angle between velocity of particle at source pt and L"
      )
      .def_prop_rw(
          "theta_lk",
          &CsrKick1Struct::theta_lk,
          &CsrKick1Struct::set_theta_lk,
          "Angle between L and velocity of kicked particle"
      )
      .def_prop_rw(
          "ix_ele_source",
          &CsrKick1Struct::ix_ele_source,
          &CsrKick1Struct::set_ix_ele_source,
          "Source element index."
      )
      .def_prop_rw(
          "floor_s",
          &CsrKick1Struct::floor_s,
          &CsrKick1Struct::set_floor_s,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Floor position of source pt"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return CsrKick1StructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = CsrKick1StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const CsrKick1Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CsrKick1Struct &self) {
            return CsrKick1Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CsrKick1Struct &self, nb::dict &memo) { return CsrKick1Struct(self); }
      )
      .def(
          "__eq__",
          [](const CsrKick1Struct &self, const CsrKick1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const CsrKick1Struct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<CsrKick1StructArray1D, CsrKick1StructAlloc1D>(
      m,
      "CsrKick1StructArray1D",
      "CsrKick1StructAlloc1D"
  );
  // 2D CsrKick1Struct arrays are not used in structs/routines
  // 3D CsrKick1Struct arrays are not used in structs/routines
}

// =============================================================================
// csr_particle_position_struct
void init_csr_particle_position_struct(nb::module_ &m, nb::class_<CsrParticlePositionStruct> &cls) {
  cls.def(
         nb::init<std::optional<std::vector<double>>, std::optional<double>>(),
         nb::arg("r") = nb::none(),
         nb::arg("charge") = nb::none()
  )
      .def_prop_rw(
          "r",
          &CsrParticlePositionStruct::r,
          &CsrParticlePositionStruct::set_r,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "particle position"
      )
      .def_prop_rw(
          "charge",
          &CsrParticlePositionStruct::charge,
          &CsrParticlePositionStruct::set_charge,
          "particle charge"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return CsrParticlePositionStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = CsrParticlePositionStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const CsrParticlePositionStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CsrParticlePositionStruct &self) {
            return CsrParticlePositionStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CsrParticlePositionStruct &self, nb::dict &memo) {
            return CsrParticlePositionStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const CsrParticlePositionStruct &self, const CsrParticlePositionStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const CsrParticlePositionStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<CsrParticlePositionStructArray1D, CsrParticlePositionStructAlloc1D>(
      m,
      "CsrParticlePositionStructArray1D",
      "CsrParticlePositionStructAlloc1D"
  );
  // 2D CsrParticlePositionStruct arrays are not used in structs/routines
  // 3D CsrParticlePositionStruct arrays are not used in structs/routines
}

// =============================================================================
// csr_struct
void init_csr_struct(nb::module_ &m, nb::class_<CsrStruct> &cls) {
  cls.def(
         "__init__",
         [](CsrStruct *self,
            std::optional<double> gamma,
            std::optional<double> gamma2,
            std::optional<double> rel_mass,
            std::optional<double> beta,
            std::optional<double> dz_slice,
            std::optional<double> ds_track_step,
            std::optional<double> s_kick,
            std::optional<double> s_chord_kick,
            std::optional<double> y_source,
            std::optional<double> kick_factor,
            std::optional<double> actual_track_step,
            std::optional<double> x0_bunch,
            std::optional<double> y0_bunch,
            const FloorPositionStruct *floor_k,
            std::optional<int> species,
            std::optional<int> ix_ele_kick,
            const EleStruct *kick_ele,
            const Mesh3dStruct *mesh3d) {
           new (self) CsrStruct(
               gamma,
               gamma2,
               rel_mass,
               beta,
               dz_slice,
               ds_track_step,
               s_kick,
               s_chord_kick,
               y_source,
               kick_factor,
               actual_track_step,
               x0_bunch,
               y0_bunch,
               ptr_to_opt_ref(floor_k),
               species,
               ix_ele_kick,
               ptr_to_opt_ref(kick_ele),
               ptr_to_opt_ref(mesh3d)
           );
         },
         nb::arg("gamma") = nb::none(),
         nb::arg("gamma2") = nb::none(),
         nb::arg("rel_mass") = nb::none(),
         nb::arg("beta") = nb::none(),
         nb::arg("dz_slice") = nb::none(),
         nb::arg("ds_track_step") = nb::none(),
         nb::arg("s_kick") = nb::none(),
         nb::arg("s_chord_kick") = nb::none(),
         nb::arg("y_source") = nb::none(),
         nb::arg("kick_factor") = nb::none(),
         nb::arg("actual_track_step") = nb::none(),
         nb::arg("x0_bunch") = nb::none(),
         nb::arg("y0_bunch") = nb::none(),
         nb::arg("floor_k") = nb::none(),
         nb::arg("species") = nb::none(),
         nb::arg("ix_ele_kick") = nb::none(),
         nb::arg("kick_ele") = nb::none(),
         nb::arg("mesh3d") = nb::none()
  )
      .def_prop_rw("gamma", &CsrStruct::gamma, &CsrStruct::set_gamma, "Relativistic gamma factor.")
      .def_prop_rw(
          "gamma2",
          &CsrStruct::gamma2,
          &CsrStruct::set_gamma2,
          "Relativistic gamma factor."
      )
      .def_prop_rw(
          "rel_mass",
          &CsrStruct::rel_mass,
          &CsrStruct::set_rel_mass,
          "m_particle / m_electron"
      )
      .def_prop_rw("beta", &CsrStruct::beta, &CsrStruct::set_beta, "Relativistic beta factor.")
      .def_prop_rw("dz_slice", &CsrStruct::dz_slice, &CsrStruct::set_dz_slice, "Bin width")
      .def_prop_rw(
          "ds_track_step",
          &CsrStruct::ds_track_step,
          &CsrStruct::set_ds_track_step,
          "True step size"
      )
      .def_prop_rw(
          "s_kick",
          &CsrStruct::s_kick,
          &CsrStruct::set_s_kick,
          "Kick point longitudinal location (element ref coords) from entrance end"
      )
      .def_prop_rw(
          "s_chord_kick",
          &CsrStruct::s_chord_kick,
          &CsrStruct::set_s_chord_kick,
          "Kick point along beam centroid line"
      )
      .def_prop_rw(
          "y_source",
          &CsrStruct::y_source,
          &CsrStruct::set_y_source,
          "Height of source particle."
      )
      .def_prop_rw(
          "kick_factor",
          &CsrStruct::kick_factor,
          &CsrStruct::set_kick_factor,
          "Coefficient to scale the kick"
      )
      .def_prop_rw(
          "actual_track_step",
          &CsrStruct::actual_track_step,
          &CsrStruct::set_actual_track_step,
          "ds_track_step scalled by Length_centroid_chord / Length_element ratio"
      )
      .def_prop_rw("x0_bunch", &CsrStruct::x0_bunch, &CsrStruct::set_x0_bunch, "Bunch centroid")
      .def_prop_rw("y0_bunch", &CsrStruct::y0_bunch, &CsrStruct::set_y0_bunch, "Bunch centroid")
      .def_prop_rw(
          "floor_k",
          &CsrStruct::floor_k,
          &CsrStruct::set_floor_k,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Floor coords at kick point"
      )
      .def_prop_rw("species", &CsrStruct::species, &CsrStruct::set_species, "Particle type")
      .def_prop_rw(
          "ix_ele_kick",
          &CsrStruct::ix_ele_kick,
          &CsrStruct::set_ix_ele_kick,
          "Same as element being tracked through."
      )
      .def_prop_ro(
          "slice",
          &CsrStruct::slice,
          nb::keep_alive<0, 1>(),
          "slice(i) refers to the i^th bunch slice."
      )
      .def_prop_ro(
          "kick1",
          &CsrStruct::kick1,
          nb::keep_alive<0, 1>(),
          "kick1(i) referes to the kick between two slices i bins apart."
      )
      .def_prop_ro(
          "eleinfo",
          &CsrStruct::eleinfo,
          nb::keep_alive<0, 1>(),
          "Element-by-element information."
      )
      .def_prop_rw(
          "kick_ele",
          &CsrStruct::kick_ele,
          &CsrStruct::set_kick_ele,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Element where the kick pt is == ele tracked through."
      )
      .def_prop_rw(
          "mesh3d",
          &CsrStruct::mesh3d,
          &CsrStruct::set_mesh3d,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro("position", &CsrStruct::position, nb::keep_alive<0, 1>())

      .def("__repr__", [](const CsrStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CsrStruct &self) {
            return CsrStruct(self); // under-the-hood fortran copy
          }
      )
      .def("__deepcopy__", [](const CsrStruct &self, nb::dict &memo) { return CsrStruct(self); })
      .def(
          "__eq__",
          [](const CsrStruct &self, const CsrStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const CsrStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D CsrStruct arrays are not used in structs/routines
  // 2D CsrStruct arrays are not used in structs/routines
  // 3D CsrStruct arrays are not used in structs/routines
}

// =============================================================================
// cmplx_field1_at_2D_pt_struct
void init_cmplx_field1_at_2D_pt_struct(nb::module_ &m, nb::class_<CmplxField1At2dPtStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>>(),
         nb::arg("f") = nb::none(),
         nb::arg("df_dx") = nb::none(),
         nb::arg("df_dy") = nb::none(),
         nb::arg("d2f_dxdy") = nb::none()
  )
      .def_prop_rw("f", &CmplxField1At2dPtStruct::f, &CmplxField1At2dPtStruct::set_f, "Field")
      .def_prop_rw(
          "df_dx",
          &CmplxField1At2dPtStruct::df_dx,
          &CmplxField1At2dPtStruct::set_df_dx,
          "Normalized field 1st derivatives"
      )
      .def_prop_rw(
          "df_dy",
          &CmplxField1At2dPtStruct::df_dy,
          &CmplxField1At2dPtStruct::set_df_dy,
          "Normalized field 1st derivatives"
      )
      .def_prop_rw(
          "d2f_dxdy",
          &CmplxField1At2dPtStruct::d2f_dxdy,
          &CmplxField1At2dPtStruct::set_d2f_dxdy,
          "Normalized field 2nd derivative"
      )

      .def("__repr__", [](const CmplxField1At2dPtStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CmplxField1At2dPtStruct &self) {
            return CmplxField1At2dPtStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CmplxField1At2dPtStruct &self, nb::dict &memo) {
            return CmplxField1At2dPtStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const CmplxField1At2dPtStruct &self, const CmplxField1At2dPtStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const CmplxField1At2dPtStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D CmplxField1At2dPtStruct arrays are not used in structs/routines
  bind_FTypeArrayND<CmplxField1At2dPtStructArray2D>(m, "CmplxField1At2DPtStructArray2D");
  // 3D CmplxField1At2dPtStruct arrays are not used in structs/routines
}

// =============================================================================
// cmplx_field1_at_3D_pt_struct
void init_cmplx_field1_at_3D_pt_struct(nb::module_ &m, nb::class_<CmplxField1At3dPtStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>>(),
         nb::arg("f") = nb::none(),
         nb::arg("df_dx") = nb::none(),
         nb::arg("df_dy") = nb::none(),
         nb::arg("df_dz") = nb::none(),
         nb::arg("d2f_dxdy") = nb::none(),
         nb::arg("d2f_dxdz") = nb::none(),
         nb::arg("d2f_dydz") = nb::none(),
         nb::arg("d3f_dxdydz") = nb::none()
  )
      .def_prop_rw("f", &CmplxField1At3dPtStruct::f, &CmplxField1At3dPtStruct::set_f, "Field")
      .def_prop_rw(
          "df_dx",
          &CmplxField1At3dPtStruct::df_dx,
          &CmplxField1At3dPtStruct::set_df_dx,
          "Normalized field 1st derivatives"
      )
      .def_prop_rw(
          "df_dy",
          &CmplxField1At3dPtStruct::df_dy,
          &CmplxField1At3dPtStruct::set_df_dy,
          "Normalized field 1st derivatives"
      )
      .def_prop_rw(
          "df_dz",
          &CmplxField1At3dPtStruct::df_dz,
          &CmplxField1At3dPtStruct::set_df_dz,
          "Normalized field 1st derivatives"
      )
      .def_prop_rw(
          "d2f_dxdy",
          &CmplxField1At3dPtStruct::d2f_dxdy,
          &CmplxField1At3dPtStruct::set_d2f_dxdy,
          "Normalized field 2nd derivatives"
      )
      .def_prop_rw(
          "d2f_dxdz",
          &CmplxField1At3dPtStruct::d2f_dxdz,
          &CmplxField1At3dPtStruct::set_d2f_dxdz,
          "Normalized field 2nd derivatives"
      )
      .def_prop_rw(
          "d2f_dydz",
          &CmplxField1At3dPtStruct::d2f_dydz,
          &CmplxField1At3dPtStruct::set_d2f_dydz,
          "Normalized field 2nd derivatives"
      )
      .def_prop_rw(
          "d3f_dxdydz",
          &CmplxField1At3dPtStruct::d3f_dxdydz,
          &CmplxField1At3dPtStruct::set_d3f_dxdydz,
          "Normalized field 3rd derivative"
      )

      .def("__repr__", [](const CmplxField1At3dPtStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CmplxField1At3dPtStruct &self) {
            return CmplxField1At3dPtStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CmplxField1At3dPtStruct &self, nb::dict &memo) {
            return CmplxField1At3dPtStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const CmplxField1At3dPtStruct &self, const CmplxField1At3dPtStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const CmplxField1At3dPtStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D CmplxField1At3dPtStruct arrays are not used in structs/routines
  // 2D CmplxField1At3dPtStruct arrays are not used in structs/routines
  bind_FTypeArrayND<CmplxField1At3dPtStructArray3D>(m, "CmplxField1At3DPtStructArray3D");
}

// =============================================================================
// cmplx_field_at_2D_box_struct
void init_cmplx_field_at_2D_box_struct(nb::module_ &m, nb::class_<CmplxFieldAt2dBoxStruct> &cls) {
  cls.def(nb::init<std::optional<std::vector<int>>>(), nb::arg("i_box") = nb::none())
      .def_prop_ro("pt", &CmplxFieldAt2dBoxStruct::pt, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "i_box",
          &CmplxFieldAt2dBoxStruct::i_box,
          &CmplxFieldAt2dBoxStruct::set_i_box,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "index at lower box corner."
      )

      .def("__repr__", [](const CmplxFieldAt2dBoxStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CmplxFieldAt2dBoxStruct &self) {
            return CmplxFieldAt2dBoxStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CmplxFieldAt2dBoxStruct &self, nb::dict &memo) {
            return CmplxFieldAt2dBoxStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const CmplxFieldAt2dBoxStruct &self, const CmplxFieldAt2dBoxStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const CmplxFieldAt2dBoxStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D CmplxFieldAt2dBoxStruct arrays are not used in structs/routines
  // 2D CmplxFieldAt2dBoxStruct arrays are not used in structs/routines
  // 3D CmplxFieldAt2dBoxStruct arrays are not used in structs/routines
}

// =============================================================================
// cmplx_field_at_3D_box_struct
void init_cmplx_field_at_3D_box_struct(nb::module_ &m, nb::class_<CmplxFieldAt3dBoxStruct> &cls) {
  cls.def(nb::init<std::optional<std::vector<int>>>(), nb::arg("i_box") = nb::none())
      .def_prop_ro("pt", &CmplxFieldAt3dBoxStruct::pt, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "i_box",
          &CmplxFieldAt3dBoxStruct::i_box,
          &CmplxFieldAt3dBoxStruct::set_i_box,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "index at lower box corner."
      )

      .def("__repr__", [](const CmplxFieldAt3dBoxStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CmplxFieldAt3dBoxStruct &self) {
            return CmplxFieldAt3dBoxStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CmplxFieldAt3dBoxStruct &self, nb::dict &memo) {
            return CmplxFieldAt3dBoxStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const CmplxFieldAt3dBoxStruct &self, const CmplxFieldAt3dBoxStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const CmplxFieldAt3dBoxStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D CmplxFieldAt3dBoxStruct arrays are not used in structs/routines
  // 2D CmplxFieldAt3dBoxStruct arrays are not used in structs/routines
  // 3D CmplxFieldAt3dBoxStruct arrays are not used in structs/routines
}

// =============================================================================
// crystal_param_struct
void init_crystal_param_struct(nb::module_ &m, nb::class_<CrystalParamStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>>(),
         nb::arg("cap_gamma") = nb::none(),
         nb::arg("dtheta_sin_2theta") = nb::none(),
         nb::arg("b_eff") = nb::none(),
         nb::arg("wavelength") = nb::none(),
         nb::arg("old_vvec") = nb::none(),
         nb::arg("new_vvec") = nb::none()
  )
      .def_prop_rw("cap_gamma", &CrystalParamStruct::cap_gamma, &CrystalParamStruct::set_cap_gamma)
      .def_prop_rw(
          "dtheta_sin_2theta",
          &CrystalParamStruct::dtheta_sin_2theta,
          &CrystalParamStruct::set_dtheta_sin_2theta
      )
      .def_prop_rw("b_eff", &CrystalParamStruct::b_eff, &CrystalParamStruct::set_b_eff)
      .def_prop_rw(
          "wavelength",
          &CrystalParamStruct::wavelength,
          &CrystalParamStruct::set_wavelength
      )
      .def_prop_rw(
          "old_vvec",
          &CrystalParamStruct::old_vvec,
          &CrystalParamStruct::set_old_vvec,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "new_vvec",
          &CrystalParamStruct::new_vvec,
          &CrystalParamStruct::set_new_vvec,
          nb::for_getter(nb::keep_alive<0, 1>())
      )

      .def("__repr__", [](const CrystalParamStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CrystalParamStruct &self) {
            return CrystalParamStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CrystalParamStruct &self, nb::dict &memo) { return CrystalParamStruct(self); }
      )
      .def(
          "__eq__",
          [](const CrystalParamStruct &self, const CrystalParamStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const CrystalParamStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D CrystalParamStruct arrays are not used in structs/routines
  // 2D CrystalParamStruct arrays are not used in structs/routines
  // 3D CrystalParamStruct arrays are not used in structs/routines
}