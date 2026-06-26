#include "pybmad/generated/structs_d.hpp"

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
// diffuse_param_struct
void init_diffuse_param_struct(nb::module_ &m, nb::class_<DiffuseParamStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>>(),
         nb::arg("x") = nb::none(),
         nb::arg("y") = nb::none(),
         nb::arg("lambda_") = nb::none(),
         nb::arg("c_norm") = nb::none(),
         nb::arg("chx_norm") = nb::none(),
         nb::arg("n_pt_spline") = nb::none()
  )
      .def_prop_rw("x", &DiffuseParamStruct::x, &DiffuseParamStruct::set_x)
      .def_prop_rw("y", &DiffuseParamStruct::y, &DiffuseParamStruct::set_y)
      .def_prop_rw("lambda_", &DiffuseParamStruct::lambda, &DiffuseParamStruct::set_lambda)
      .def_prop_rw("c_norm", &DiffuseParamStruct::c_norm, &DiffuseParamStruct::set_c_norm)
      .def_prop_rw("chx_norm", &DiffuseParamStruct::chx_norm, &DiffuseParamStruct::set_chx_norm)
      .def_prop_ro("prob_spline", &DiffuseParamStruct::prob_spline, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "n_pt_spline",
          &DiffuseParamStruct::n_pt_spline,
          &DiffuseParamStruct::set_n_pt_spline
      )

      .def("__repr__", [](const DiffuseParamStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const DiffuseParamStruct &self) {
            return DiffuseParamStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const DiffuseParamStruct &self, nb::dict &memo) { return DiffuseParamStruct(self); }
      )
      .def(
          "__eq__",
          [](const DiffuseParamStruct &self, const DiffuseParamStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const DiffuseParamStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D DiffuseParamStruct arrays are not used in structs/routines
  // 2D DiffuseParamStruct arrays are not used in structs/routines
  // 3D DiffuseParamStruct arrays are not used in structs/routines
}

// =============================================================================
// do_loop_struct
void init_do_loop_struct(nb::module_ &m, nb::class_<DoLoopStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::string>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>>(),
         nb::arg("name") = nb::none(),
         nb::arg("index") = nb::none(),
         nb::arg("start") = nb::none(),
         nb::arg("end") = nb::none(),
         nb::arg("step") = nb::none(),
         nb::arg("n_line_start") = nb::none(),
         nb::arg("n_line_end") = nb::none(),
         nb::arg("value") = nb::none()
  )
      .def_prop_rw("name", &DoLoopStruct::name, &DoLoopStruct::set_name, "do loop index name")
      .def_prop_rw("index", &DoLoopStruct::index, &DoLoopStruct::set_index, "for do loops")
      .def_prop_rw("start", &DoLoopStruct::start, &DoLoopStruct::set_start, "for do loops")
      .def_prop_rw("end", &DoLoopStruct::end, &DoLoopStruct::set_end, "for do loops")
      .def_prop_rw("step", &DoLoopStruct::step, &DoLoopStruct::set_step, "for do loops")
      .def_prop_rw(
          "n_line_start",
          &DoLoopStruct::n_line_start,
          &DoLoopStruct::set_n_line_start,
          "lines in each nested loop"
      )
      .def_prop_rw(
          "n_line_end",
          &DoLoopStruct::n_line_end,
          &DoLoopStruct::set_n_line_end,
          "lines in each nested loop"
      )
      .def_prop_rw("value", &DoLoopStruct::value, &DoLoopStruct::set_value)
      .def_static(
          "new_array1d",
          [](int sz) { return DoLoopStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = DoLoopStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const DoLoopStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const DoLoopStruct &self) {
            return DoLoopStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const DoLoopStruct &self, nb::dict &memo) { return DoLoopStruct(self); }
      )
      .def(
          "__eq__",
          [](const DoLoopStruct &self, const DoLoopStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const DoLoopStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<DoLoopStructArray1D, DoLoopStructAlloc1D>(
      m,
      "DoLoopStructArray1D",
      "DoLoopStructAlloc1D"
  );
  // 2D DoLoopStruct arrays are not used in structs/routines
  // 3D DoLoopStruct arrays are not used in structs/routines
}