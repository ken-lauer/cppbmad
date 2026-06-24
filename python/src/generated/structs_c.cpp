#include "pybmad/generated/structs_c.hpp"

#include <cstdint>
#include <functional>

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// cartesian_map_struct
void init_cartesian_map_struct(py::module &m, py::class_<CartesianMapStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             optional_ref<const CartesianMapTermStruct>>(),
         py::arg("field_scale") = py::none(),
         py::arg("r0") = py::none(),
         py::arg("master_parameter") = py::none(),
         py::arg("ele_anchor_pt") = py::none(),
         py::arg("field_type") = py::none(),
         py::arg("ptr") = py::none()
  )
      .def_property(
          "field_scale",
          &CartesianMapStruct::field_scale,
          &CartesianMapStruct::set_field_scale,
          "Factor to scale the fields by"
      )
      .def_property(
          "r0",
          py::cpp_function(&CartesianMapStruct::r0, py::keep_alive<0, 1>()),
          &CartesianMapStruct::set_r0,
          "Field origin offset."
      )
      .def_property(
          "master_parameter",
          &CartesianMapStruct::master_parameter,
          &CartesianMapStruct::set_master_parameter,
          "Master parameter in ele%value(:) array to use for scaling the field."
      )
      .def_property(
          "ele_anchor_pt",
          &CartesianMapStruct::ele_anchor_pt,
          &CartesianMapStruct::set_ele_anchor_pt,
          "anchor_beginning$, anchor_center$, or anchor_end$"
      )
      .def_property(
          "field_type",
          &CartesianMapStruct::field_type,
          &CartesianMapStruct::set_field_type,
          "or electric$"
      )
      .def_property(
          "ptr",
          py::cpp_function(&CartesianMapStruct::ptr, py::keep_alive<0, 1>()),
          &CartesianMapStruct::set_ptr
      )
      .def_static(
          "new_array1d",
          [](int sz) { return CartesianMapStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = CartesianMapStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const CartesianMapStruct &self, py::dict &memo) { return CartesianMapStruct(self); }
      )
      .def(
          "__eq__",
          [](const CartesianMapStruct &self, const CartesianMapStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const CartesianMapStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
void init_cartesian_map_term1_struct(py::module &m, py::class_<CartesianMapTerm1Struct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>,
             std::optional<int>>(),
         py::arg("coef") = py::none(),
         py::arg("kx") = py::none(),
         py::arg("ky") = py::none(),
         py::arg("kz") = py::none(),
         py::arg("x0") = py::none(),
         py::arg("y0") = py::none(),
         py::arg("phi_z") = py::none(),
         py::arg("family") = py::none(),
         py::arg("form") = py::none()
  )
      .def_property("coef", &CartesianMapTerm1Struct::coef, &CartesianMapTerm1Struct::set_coef)
      .def_property("kx", &CartesianMapTerm1Struct::kx, &CartesianMapTerm1Struct::set_kx)
      .def_property("ky", &CartesianMapTerm1Struct::ky, &CartesianMapTerm1Struct::set_ky)
      .def_property("kz", &CartesianMapTerm1Struct::kz, &CartesianMapTerm1Struct::set_kz)
      .def_property("x0", &CartesianMapTerm1Struct::x0, &CartesianMapTerm1Struct::set_x0)
      .def_property("y0", &CartesianMapTerm1Struct::y0, &CartesianMapTerm1Struct::set_y0)
      .def_property("phi_z", &CartesianMapTerm1Struct::phi_z, &CartesianMapTerm1Struct::set_phi_z)
      .def_property(
          "family",
          &CartesianMapTerm1Struct::family,
          &CartesianMapTerm1Struct::set_family,
          "family_x$, etc."
      )
      .def_property(
          "form",
          &CartesianMapTerm1Struct::form,
          &CartesianMapTerm1Struct::set_form,
          "hyper_y$, etc."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return CartesianMapTerm1StructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = CartesianMapTerm1StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const CartesianMapTerm1Struct &self, py::dict &memo) {
            return CartesianMapTerm1Struct(self);
          }
      )
      .def(
          "__eq__",
          [](const CartesianMapTerm1Struct &self, const CartesianMapTerm1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const CartesianMapTerm1Struct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
void init_cartesian_map_term_struct(py::module &m, py::class_<CartesianMapTermStruct> &cls) {
  cls.def(
         py::init<std::optional<std::string>, std::optional<int>>(),
         py::arg("file") = py::none(),
         py::arg("n_link") = py::none()
  )
      .def_property(
          "file",
          &CartesianMapTermStruct::file,
          &CartesianMapTermStruct::set_file,
          "Input file name. Used also as ID for instances."
      )
      .def_property(
          "n_link",
          &CartesianMapTermStruct::n_link,
          &CartesianMapTermStruct::set_n_link,
          "For memory management of %term"
      )
      .def_property_readonly(
          "term",
          py::cpp_function(&CartesianMapTermStruct::term, py::keep_alive<0, 1>())
      )

      .def("__repr__", [](const CartesianMapTermStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CartesianMapTermStruct &self) {
            return CartesianMapTermStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CartesianMapTermStruct &self, py::dict &memo) {
            return CartesianMapTermStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const CartesianMapTermStruct &self, const CartesianMapTermStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const CartesianMapTermStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D CartesianMapTermStruct arrays are not used in structs/routines
  // 2D CartesianMapTermStruct arrays are not used in structs/routines
  // 3D CartesianMapTermStruct arrays are not used in structs/routines
}

// =============================================================================
// complex_taylor_struct
void init_complex_taylor_struct(py::module &m, py::class_<ComplexTaylorStruct> &cls) {
  cls.def(py::init<std::optional<std::complex<double>>>(), py::arg("ref") = py::none())
      .def_property("ref", &ComplexTaylorStruct::ref, &ComplexTaylorStruct::set_ref)
      .def_property_readonly(
          "term",
          py::cpp_function(&ComplexTaylorStruct::term, py::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return ComplexTaylorStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ComplexTaylorStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const ComplexTaylorStruct &self, py::dict &memo) { return ComplexTaylorStruct(self); }
      )
      .def(
          "__eq__",
          [](const ComplexTaylorStruct &self, const ComplexTaylorStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const ComplexTaylorStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
void init_complex_taylor_term_struct(py::module &m, py::class_<ComplexTaylorTermStruct> &cls) {
  cls.def(
         py::init<std::optional<std::complex<double>>, std::optional<std::vector<int>>>(),
         py::arg("coef") = py::none(),
         py::arg("expn") = py::none()
  )
      .def_property("coef", &ComplexTaylorTermStruct::coef, &ComplexTaylorTermStruct::set_coef)
      .def_property(
          "expn",
          py::cpp_function(&ComplexTaylorTermStruct::expn, py::keep_alive<0, 1>()),
          &ComplexTaylorTermStruct::set_expn
      )
      .def_static(
          "new_array1d",
          [](int sz) { return ComplexTaylorTermStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ComplexTaylorTermStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const ComplexTaylorTermStruct &self, py::dict &memo) {
            return ComplexTaylorTermStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const ComplexTaylorTermStruct &self, const ComplexTaylorTermStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const ComplexTaylorTermStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
void init_control_ramp1_struct(py::module &m, py::class_<ControlRamp1Struct> &cls) {
  cls.def(
         py::init<
             std::optional<std::vector<double>>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<bool>>(),
         py::arg("y_knot") = py::none(),
         py::arg("attribute") = py::none(),
         py::arg("slave_name") = py::none(),
         py::arg("is_controller") = py::none()
  )
      .def_property(
          "y_knot",
          py::cpp_function(&ControlRamp1Struct::y_knot, py::keep_alive<0, 1>()),
          &ControlRamp1Struct::set_y_knot
      )
      .def_property_readonly(
          "stack",
          py::cpp_function(&ControlRamp1Struct::stack, py::keep_alive<0, 1>()),
          "Evaluation stack"
      )
      .def_property(
          "attribute",
          &ControlRamp1Struct::attribute,
          &ControlRamp1Struct::set_attribute,
          "Name of attribute controlled. Set to 'FIELD_OVERLAPS' for field overlaps."
      )
      .def_property(
          "slave_name",
          &ControlRamp1Struct::slave_name,
          &ControlRamp1Struct::set_slave_name,
          "Name of slave."
      )
      .def_property(
          "is_controller",
          &ControlRamp1Struct::is_controller,
          &ControlRamp1Struct::set_is_controller,
          "Is the slave a controller? If so bookkeeping is different."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return ControlRamp1StructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ControlRamp1StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const ControlRamp1Struct &self, py::dict &memo) { return ControlRamp1Struct(self); }
      )
      .def(
          "__eq__",
          [](const ControlRamp1Struct &self, const ControlRamp1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const ControlRamp1Struct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
void init_control_struct(py::module &m, py::class_<ControlStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<std::vector<double>>,
             optional_ref<const LatEleLocStruct>,
             optional_ref<const LatEleLocStruct>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<int>>(),
         py::arg("value") = py::none(),
         py::arg("y_knot") = py::none(),
         py::arg("slave") = py::none(),
         py::arg("lord") = py::none(),
         py::arg("slave_name") = py::none(),
         py::arg("attribute") = py::none(),
         py::arg("ix_attrib") = py::none()
  )
      .def_property(
          "value",
          &ControlStruct::value,
          &ControlStruct::set_value,
          "Used by group, and overlay elements."
      )
      .def_property(
          "y_knot",
          py::cpp_function(&ControlStruct::y_knot, py::keep_alive<0, 1>()),
          &ControlStruct::set_y_knot
      )
      .def_property_readonly(
          "stack",
          py::cpp_function(&ControlStruct::stack, py::keep_alive<0, 1>()),
          "Evaluation stack"
      )
      .def_property(
          "slave",
          py::cpp_function(&ControlStruct::slave, py::keep_alive<0, 1>()),
          &ControlStruct::set_slave
      )
      .def_property(
          "lord",
          py::cpp_function(&ControlStruct::lord, py::keep_alive<0, 1>()),
          &ControlStruct::set_lord
      )
      .def_property(
          "slave_name",
          &ControlStruct::slave_name,
          &ControlStruct::set_slave_name,
          "Name of slave."
      )
      .def_property(
          "attribute",
          &ControlStruct::attribute,
          &ControlStruct::set_attribute,
          "Name of attribute controlled. Set to 'FIELD_OVERLAPS' for field overlaps. Set to "
          "'INPUT' or 'OUTPUT' for feedback slaves."
      )
      .def_property(
          "ix_attrib",
          &ControlStruct::ix_attrib,
          &ControlStruct::set_ix_attrib,
          "Index of attribute controlled. See note above!"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return ControlStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ControlStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const ControlStruct &self, py::dict &memo) { return ControlStruct(self); }
      )
      .def(
          "__eq__",
          [](const ControlStruct &self, const ControlStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const ControlStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
void init_control_var1_struct(py::module &m, py::class_<ControlVar1Struct> &cls) {
  cls.def(
         py::init<std::optional<std::string>, std::optional<double>, std::optional<double>>(),
         py::arg("name") = py::none(),
         py::arg("value") = py::none(),
         py::arg("old_value") = py::none()
  )
      .def_property("name", &ControlVar1Struct::name, &ControlVar1Struct::set_name)
      .def_property("value", &ControlVar1Struct::value, &ControlVar1Struct::set_value)
      .def_property("old_value", &ControlVar1Struct::old_value, &ControlVar1Struct::set_old_value)
      .def_static(
          "new_array1d",
          [](int sz) { return ControlVar1StructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ControlVar1StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const ControlVar1Struct &self, py::dict &memo) { return ControlVar1Struct(self); }
      )
      .def(
          "__eq__",
          [](const ControlVar1Struct &self, const ControlVar1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const ControlVar1Struct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
void init_controller_struct(py::module &m, py::class_<ControllerStruct> &cls) {
  cls.def(py::init<std::optional<std::vector<double>>>(), py::arg("x_knot") = py::none())
      .def_property_readonly(
          "var",
          py::cpp_function(&ControllerStruct::var, py::keep_alive<0, 1>())
      )
      .def_property_readonly(
          "ramp",
          py::cpp_function(&ControllerStruct::ramp, py::keep_alive<0, 1>()),
          "For ramper lord elements"
      )
      .def_property_readonly(
          "ramper_lord",
          py::cpp_function(&ControllerStruct::ramper_lord, py::keep_alive<0, 1>()),
          "Ramper lord info for this slave"
      )
      .def_property(
          "x_knot",
          py::cpp_function(&ControllerStruct::x_knot, py::keep_alive<0, 1>()),
          &ControllerStruct::set_x_knot
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
          [](const ControllerStruct &self, py::dict &memo) { return ControllerStruct(self); }
      )
      .def(
          "__eq__",
          [](const ControllerStruct &self, const ControllerStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const ControllerStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D ControllerStruct arrays are not used in structs/routines
  // 2D ControllerStruct arrays are not used in structs/routines
  // 3D ControllerStruct arrays are not used in structs/routines
}

// =============================================================================
// coord_array_struct
void init_coord_array_struct(py::module &m, py::class_<CoordArrayStruct> &cls) {
  cls.def(py::init<>())
      .def_property_readonly(
          "orbit",
          py::cpp_function(&CoordArrayStruct::orbit, py::keep_alive<0, 1>())
      )
      .def_static(
          "new_array1d",
          [](int sz) { return CoordArrayStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = CoordArrayStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const CoordArrayStruct &self, py::dict &memo) { return CoordArrayStruct(self); }
      )
      .def(
          "__eq__",
          [](const CoordArrayStruct &self, const CoordArrayStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const CoordArrayStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
void init_coord_struct(py::module &m, py::class_<CoordStruct> &cls) {
  cls.def(
         py::init<
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
         py::arg("vec") = py::none(),
         py::arg("s") = py::none(),
         py::arg("t") = py::none(),
         py::arg("spin") = py::none(),
         py::arg("field") = py::none(),
         py::arg("phase") = py::none(),
         py::arg("charge") = py::none(),
         py::arg("dt_ref") = py::none(),
         py::arg("r") = py::none(),
         py::arg("p0c") = py::none(),
         py::arg("E_potential") = py::none(),
         py::arg("beta") = py::none(),
         py::arg("ix_ele") = py::none(),
         py::arg("ix_branch") = py::none(),
         py::arg("ix_turn") = py::none(),
         py::arg("ix_user") = py::none(),
         py::arg("state") = py::none(),
         py::arg("direction") = py::none(),
         py::arg("time_dir") = py::none(),
         py::arg("species") = py::none(),
         py::arg("location") = py::none()
  )
      .def_property(
          "vec",
          py::cpp_function(&CoordStruct::vec, py::keep_alive<0, 1>()),
          &CoordStruct::set_vec,
          "(x, px, y, py, z, pz). Generally phase space for charged particles. See Bmad manual."
      )
      .def_property("s", &CoordStruct::s, &CoordStruct::set_s, "Longitudinal position")
      .def_property(
          "t",
          &CoordStruct::t,
          &CoordStruct::set_t,
          "Absolute time (not relative to reference). Note: Quad precision!"
      )
      .def_property(
          "spin",
          py::cpp_function(&CoordStruct::spin, py::keep_alive<0, 1>()),
          &CoordStruct::set_spin,
          "Spin."
      )
      .def_property(
          "field",
          py::cpp_function(&CoordStruct::field, py::keep_alive<0, 1>()),
          &CoordStruct::set_field,
          "Photon E-field intensity (x,y)."
      )
      .def_property(
          "phase",
          py::cpp_function(&CoordStruct::phase, py::keep_alive<0, 1>()),
          &CoordStruct::set_phase,
          "Photon E-field phase (x,y). For charged particles, phase(1) is RF phase."
      )
      .def_property(
          "charge",
          &CoordStruct::charge,
          &CoordStruct::set_charge,
          "Macroparticle weight (which is different from particle species charge). For some space "
          "charge calcs the weight is in Coulombs."
      )
      .def_property(
          "dt_ref",
          &CoordStruct::dt_ref,
          &CoordStruct::set_dt_ref,
          "Used in: * time tracking for computing z. * by coherent photons = path_length/c_light."
      )
      .def_property("r", &CoordStruct::r, &CoordStruct::set_r, "For general use. Not used by Bmad.")
      .def_property(
          "p0c",
          &CoordStruct::p0c,
          &CoordStruct::set_p0c,
          "For non-photons: Reference momentum. For photons: Photon momentum (not reference)."
      )
      .def_property(
          "E_potential",
          &CoordStruct::E_potential,
          &CoordStruct::set_E_potential,
          "Potential energy."
      )
      .def_property("beta", &CoordStruct::beta, &CoordStruct::set_beta, "Velocity / c_light.")
      .def_property(
          "ix_ele",
          &CoordStruct::ix_ele,
          &CoordStruct::set_ix_ele,
          "Index of the lattice element the particle is in. May be -1 if element is not associated "
          "with a lattice."
      )
      .def_property(
          "ix_branch",
          &CoordStruct::ix_branch,
          &CoordStruct::set_ix_branch,
          "Index of the lattice branch the particle is in."
      )
      .def_property(
          "ix_turn",
          &CoordStruct::ix_turn,
          &CoordStruct::set_ix_turn,
          "Turn index for multiturn tracking."
      )
      .def_property(
          "ix_user",
          &CoordStruct::ix_user,
          &CoordStruct::set_ix_user,
          "For general use, not used by Bmad."
      )
      .def_property(
          "state",
          &CoordStruct::state,
          &CoordStruct::set_state,
          "alive$, lost$, lost_neg_x_aperture$, lost_pz$, etc."
      )
      .def_property(
          "direction",
          &CoordStruct::direction,
          &CoordStruct::set_direction,
          "+1 or -1. Sign of longitudinal direction of motion (ds/dt). This is independent of the "
          "element orientation."
      )
      .def_property(
          "time_dir",
          &CoordStruct::time_dir,
          &CoordStruct::set_time_dir,
          "+1 or -1. Time direction. -1 => Traveling backwards in time."
      )
      .def_property(
          "species",
          &CoordStruct::species,
          &CoordStruct::set_species,
          "positron$, proton$, etc."
      )
      .def_property(
          "location",
          &CoordStruct::location,
          &CoordStruct::set_location,
          "upstream_end$, inside$, or downstream_end$"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return CoordStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = CoordStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const CoordStruct &self, py::dict &memo) { return CoordStruct(self); }
      )
      .def(
          "__eq__",
          [](const CoordStruct &self, const CoordStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const CoordStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
void init_cylindrical_map_struct(py::module &m, py::class_<CylindricalMapStruct> &cls) {
  cls.def(
         py::init<
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>,
             std::optional<int>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             optional_ref<const CylindricalMapTermStruct>>(),
         py::arg("m") = py::none(),
         py::arg("harmonic") = py::none(),
         py::arg("phi0_fieldmap") = py::none(),
         py::arg("theta0_azimuth") = py::none(),
         py::arg("field_scale") = py::none(),
         py::arg("master_parameter") = py::none(),
         py::arg("ele_anchor_pt") = py::none(),
         py::arg("dz") = py::none(),
         py::arg("r0") = py::none(),
         py::arg("ptr") = py::none()
  )
      .def_property(
          "m",
          &CylindricalMapStruct::m,
          &CylindricalMapStruct::set_m,
          "Azimuthal Mode: varies as cos(m*phi - theta0_azimuth)"
      )
      .def_property(
          "harmonic",
          &CylindricalMapStruct::harmonic,
          &CylindricalMapStruct::set_harmonic,
          "Harmonic of fundamental"
      )
      .def_property(
          "phi0_fieldmap",
          &CylindricalMapStruct::phi0_fieldmap,
          &CylindricalMapStruct::set_phi0_fieldmap,
          "Mode oscillates as: twopi * (f * t + phi0_fieldmap)"
      )
      .def_property(
          "theta0_azimuth",
          &CylindricalMapStruct::theta0_azimuth,
          &CylindricalMapStruct::set_theta0_azimuth,
          "Azimuthal ((x, y) plane) orientation of mode."
      )
      .def_property(
          "field_scale",
          &CylindricalMapStruct::field_scale,
          &CylindricalMapStruct::set_field_scale,
          "Factor to scale the fields by"
      )
      .def_property(
          "master_parameter",
          &CylindricalMapStruct::master_parameter,
          &CylindricalMapStruct::set_master_parameter,
          "Master parameter in ele%value(:) array to use for scaling the field."
      )
      .def_property(
          "ele_anchor_pt",
          &CylindricalMapStruct::ele_anchor_pt,
          &CylindricalMapStruct::set_ele_anchor_pt,
          "anchor_beginning$, anchor_center$, or anchor_end$"
      )
      .def_property(
          "dz",
          &CylindricalMapStruct::dz,
          &CylindricalMapStruct::set_dz,
          "Distance between sampled field points."
      )
      .def_property(
          "r0",
          py::cpp_function(&CylindricalMapStruct::r0, py::keep_alive<0, 1>()),
          &CylindricalMapStruct::set_r0,
          "Field origin offset."
      )
      .def_property(
          "ptr",
          py::cpp_function(&CylindricalMapStruct::ptr, py::keep_alive<0, 1>()),
          &CylindricalMapStruct::set_ptr
      )
      .def_static(
          "new_array1d",
          [](int sz) { return CylindricalMapStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = CylindricalMapStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const CylindricalMapStruct &self, py::dict &memo) {
            return CylindricalMapStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const CylindricalMapStruct &self, const CylindricalMapStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const CylindricalMapStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
void init_cylindrical_map_term1_struct(py::module &m, py::class_<CylindricalMapTerm1Struct> &cls) {
  cls.def(
         py::init<std::optional<std::complex<double>>, std::optional<std::complex<double>>>(),
         py::arg("e_coef") = py::none(),
         py::arg("b_coef") = py::none()
  )
      .def_property(
          "e_coef",
          &CylindricalMapTerm1Struct::e_coef,
          &CylindricalMapTerm1Struct::set_e_coef
      )
      .def_property(
          "b_coef",
          &CylindricalMapTerm1Struct::b_coef,
          &CylindricalMapTerm1Struct::set_b_coef
      )
      .def_static(
          "new_array1d",
          [](int sz) { return CylindricalMapTerm1StructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = CylindricalMapTerm1StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
          [](const CylindricalMapTerm1Struct &self, py::dict &memo) {
            return CylindricalMapTerm1Struct(self);
          }
      )
      .def(
          "__eq__",
          [](const CylindricalMapTerm1Struct &self, const CylindricalMapTerm1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const CylindricalMapTerm1Struct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
void init_cylindrical_map_term_struct(py::module &m, py::class_<CylindricalMapTermStruct> &cls) {
  cls.def(
         py::init<std::optional<std::string>, std::optional<int>>(),
         py::arg("file") = py::none(),
         py::arg("n_link") = py::none()
  )
      .def_property(
          "file",
          &CylindricalMapTermStruct::file,
          &CylindricalMapTermStruct::set_file,
          "Input file name. Used also as ID for instances."
      )
      .def_property(
          "n_link",
          &CylindricalMapTermStruct::n_link,
          &CylindricalMapTermStruct::set_n_link,
          "For memory management of this structure"
      )
      .def_property_readonly(
          "term",
          py::cpp_function(&CylindricalMapTermStruct::term, py::keep_alive<0, 1>())
      )

      .def("__repr__", [](const CylindricalMapTermStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const CylindricalMapTermStruct &self) {
            return CylindricalMapTermStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const CylindricalMapTermStruct &self, py::dict &memo) {
            return CylindricalMapTermStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const CylindricalMapTermStruct &self, const CylindricalMapTermStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const CylindricalMapTermStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D CylindricalMapTermStruct arrays are not used in structs/routines
  // 2D CylindricalMapTermStruct arrays are not used in structs/routines
  // 3D CylindricalMapTermStruct arrays are not used in structs/routines
}