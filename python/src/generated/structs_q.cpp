#include "pybmad/generated/structs_q.hpp"

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
// qp_axis_struct
void init_qp_axis_struct(nb::module_ &m, nb::class_<QpAxisStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::string>,
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
             std::optional<std::string>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<bool>>(),
         nb::arg("label") = nb::none(),
         nb::arg("min") = nb::none(),
         nb::arg("max") = nb::none(),
         nb::arg("tick_min") = nb::none(),
         nb::arg("tick_max") = nb::none(),
         nb::arg("eval_min") = nb::none(),
         nb::arg("eval_max") = nb::none(),
         nb::arg("dtick") = nb::none(),
         nb::arg("number_offset") = nb::none(),
         nb::arg("label_offset") = nb::none(),
         nb::arg("major_tick_len") = nb::none(),
         nb::arg("minor_tick_len") = nb::none(),
         nb::arg("label_color") = nb::none(),
         nb::arg("major_div") = nb::none(),
         nb::arg("major_div_nominal") = nb::none(),
         nb::arg("minor_div") = nb::none(),
         nb::arg("minor_div_max") = nb::none(),
         nb::arg("places") = nb::none(),
         nb::arg("type") = nb::none(),
         nb::arg("bounds") = nb::none(),
         nb::arg("tick_side") = nb::none(),
         nb::arg("number_side") = nb::none(),
         nb::arg("draw_label") = nb::none(),
         nb::arg("draw_numbers") = nb::none()
  )
      .def_prop_rw("label", &QpAxisStruct::label, &QpAxisStruct::set_label)
      .def_prop_rw("min", &QpAxisStruct::min, &QpAxisStruct::set_min, "Axis min/max in data units.")
      .def_prop_rw("max", &QpAxisStruct::max, &QpAxisStruct::set_max, "Axis min/max in data units.")
      .def_prop_rw(
          "tick_min",
          &QpAxisStruct::tick_min,
          &QpAxisStruct::set_tick_min,
          "Min tick location along axis in data units."
      )
      .def_prop_rw(
          "tick_max",
          &QpAxisStruct::tick_max,
          &QpAxisStruct::set_tick_max,
          "Max tick location along axis in data units."
      )
      .def_prop_rw(
          "eval_min",
          &QpAxisStruct::eval_min,
          &QpAxisStruct::set_eval_min,
          "For general use. Not set by quick_plot."
      )
      .def_prop_rw(
          "eval_max",
          &QpAxisStruct::eval_max,
          &QpAxisStruct::set_eval_max,
          "For general use. Not set by quick_plot."
      )
      .def_prop_rw(
          "dtick",
          &QpAxisStruct::dtick,
          &QpAxisStruct::set_dtick,
          "Distance between ticks. In data units. Ticks will be drawn between %min and %max."
      )
      .def_prop_rw(
          "number_offset",
          &QpAxisStruct::number_offset,
          &QpAxisStruct::set_number_offset,
          "Offset from axis line in inches."
      )
      .def_prop_rw(
          "label_offset",
          &QpAxisStruct::label_offset,
          &QpAxisStruct::set_label_offset,
          "Offset from numbers in inches."
      )
      .def_prop_rw(
          "major_tick_len",
          &QpAxisStruct::major_tick_len,
          &QpAxisStruct::set_major_tick_len,
          "In inches."
      )
      .def_prop_rw(
          "minor_tick_len",
          &QpAxisStruct::minor_tick_len,
          &QpAxisStruct::set_minor_tick_len,
          "In inches."
      )
      .def_prop_rw(
          "label_color",
          &QpAxisStruct::label_color,
          &QpAxisStruct::set_label_color,
          "Color of the label."
      )
      .def_prop_rw(
          "major_div",
          &QpAxisStruct::major_div,
          &QpAxisStruct::set_major_div,
          "Actual numbrer of major divisions"
      )
      .def_prop_rw(
          "major_div_nominal",
          &QpAxisStruct::major_div_nominal,
          &QpAxisStruct::set_major_div_nominal,
          "Nominal value."
      )
      .def_prop_rw(
          "minor_div",
          &QpAxisStruct::minor_div,
          &QpAxisStruct::set_minor_div,
          "0 = auto choose."
      )
      .def_prop_rw(
          "minor_div_max",
          &QpAxisStruct::minor_div_max,
          &QpAxisStruct::set_minor_div_max,
          "Max number for auto choose."
      )
      .def_prop_rw(
          "places",
          &QpAxisStruct::places,
          &QpAxisStruct::set_places,
          "Number of places after the decimal point to print."
      )
      .def_prop_rw("type", &QpAxisStruct::type, &QpAxisStruct::set_type, "Or 'LOG', or 'CUSTOM'")
      .def_prop_rw(
          "bounds",
          &QpAxisStruct::bounds,
          &QpAxisStruct::set_bounds,
          "Or 'ZERO_AT_END' or 'ZERO_SYMMETRIC'"
      )
      .def_prop_rw(
          "tick_side",
          &QpAxisStruct::tick_side,
          &QpAxisStruct::set_tick_side,
          "+1 = Draw on the side inside the graph, 0 = both (longer tick), -1 = outside."
      )
      .def_prop_rw(
          "number_side",
          &QpAxisStruct::number_side,
          &QpAxisStruct::set_number_side,
          "+1 = Draw to the side inside the graph, -1 = outside."
      )
      .def_prop_rw("draw_label", &QpAxisStruct::draw_label, &QpAxisStruct::set_draw_label)
      .def_prop_rw("draw_numbers", &QpAxisStruct::draw_numbers, &QpAxisStruct::set_draw_numbers)

      .def("__repr__", [](const QpAxisStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpAxisStruct &self) {
            return QpAxisStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const QpAxisStruct &self, nb::dict &memo) { return QpAxisStruct(self); }
      )
      .def(
          "__eq__",
          [](const QpAxisStruct &self, const QpAxisStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const QpAxisStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D QpAxisStruct arrays are not used in structs/routines
  // 2D QpAxisStruct arrays are not used in structs/routines
  // 3D QpAxisStruct arrays are not used in structs/routines
}

// =============================================================================
// qp_legend_struct
void init_qp_legend_struct(nb::module_ &m, nb::class_<QpLegendStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>>(),
         nb::arg("row_spacing") = nb::none(),
         nb::arg("line_length") = nb::none(),
         nb::arg("text_offset") = nb::none(),
         nb::arg("draw_line") = nb::none(),
         nb::arg("draw_symbol") = nb::none(),
         nb::arg("draw_text") = nb::none()
  )
      .def_prop_rw(
          "row_spacing",
          &QpLegendStruct::row_spacing,
          &QpLegendStruct::set_row_spacing,
          "Spacing between rows."
      )
      .def_prop_rw(
          "line_length",
          &QpLegendStruct::line_length,
          &QpLegendStruct::set_line_length,
          "Length of the line in points."
      )
      .def_prop_rw(
          "text_offset",
          &QpLegendStruct::text_offset,
          &QpLegendStruct::set_text_offset,
          "Horizontal offset in points between the line and the text."
      )
      .def_prop_rw(
          "draw_line",
          &QpLegendStruct::draw_line,
          &QpLegendStruct::set_draw_line,
          "Draw lines?"
      )
      .def_prop_rw(
          "draw_symbol",
          &QpLegendStruct::draw_symbol,
          &QpLegendStruct::set_draw_symbol,
          "Draw symbols?"
      )
      .def_prop_rw(
          "draw_text",
          &QpLegendStruct::draw_text,
          &QpLegendStruct::set_draw_text,
          "Draw text?"
      )

      .def("__repr__", [](const QpLegendStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpLegendStruct &self) {
            return QpLegendStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const QpLegendStruct &self, nb::dict &memo) { return QpLegendStruct(self); }
      )
      .def(
          "__eq__",
          [](const QpLegendStruct &self, const QpLegendStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const QpLegendStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D QpLegendStruct arrays are not used in structs/routines
  // 2D QpLegendStruct arrays are not used in structs/routines
  // 3D QpLegendStruct arrays are not used in structs/routines
}

// =============================================================================
// qp_line_struct
void init_qp_line_struct(nb::module_ &m, nb::class_<QpLineStruct> &cls) {
  cls.def(
         nb::init<std::optional<int>, std::optional<std::string>, std::optional<std::string>>(),
         nb::arg("width") = nb::none(),
         nb::arg("color") = nb::none(),
         nb::arg("pattern") = nb::none()
  )
      .def_prop_rw("width", &QpLineStruct::width, &QpLineStruct::set_width)
      .def_prop_rw("color", &QpLineStruct::color, &QpLineStruct::set_color)
      .def_prop_rw("pattern", &QpLineStruct::pattern, &QpLineStruct::set_pattern)

      .def("__repr__", [](const QpLineStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpLineStruct &self) {
            return QpLineStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const QpLineStruct &self, nb::dict &memo) { return QpLineStruct(self); }
      )
      .def(
          "__eq__",
          [](const QpLineStruct &self, const QpLineStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const QpLineStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D QpLineStruct arrays are not used in structs/routines
  // 2D QpLineStruct arrays are not used in structs/routines
  // 3D QpLineStruct arrays are not used in structs/routines
}

// =============================================================================
// qp_point_struct
void init_qp_point_struct(nb::module_ &m, nb::class_<QpPointStruct> &cls) {
  cls.def(
         nb::init<std::optional<double>, std::optional<double>, std::optional<std::string>>(),
         nb::arg("x") = nb::none(),
         nb::arg("y") = nb::none(),
         nb::arg("units") = nb::none()
  )
      .def_prop_rw("x", &QpPointStruct::x, &QpPointStruct::set_x)
      .def_prop_rw("y", &QpPointStruct::y, &QpPointStruct::set_y)
      .def_prop_rw("units", &QpPointStruct::units, &QpPointStruct::set_units)

      .def("__repr__", [](const QpPointStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpPointStruct &self) {
            return QpPointStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const QpPointStruct &self, nb::dict &memo) { return QpPointStruct(self); }
      )
      .def(
          "__eq__",
          [](const QpPointStruct &self, const QpPointStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const QpPointStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D QpPointStruct arrays are not used in structs/routines
  // 2D QpPointStruct arrays are not used in structs/routines
  // 3D QpPointStruct arrays are not used in structs/routines
}

// =============================================================================
// qp_rect_struct
void init_qp_rect_struct(nb::module_ &m, nb::class_<QpRectStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::string>>(),
         nb::arg("x1") = nb::none(),
         nb::arg("x2") = nb::none(),
         nb::arg("y1") = nb::none(),
         nb::arg("y2") = nb::none(),
         nb::arg("units") = nb::none()
  )
      .def_prop_rw("x1", &QpRectStruct::x1, &QpRectStruct::set_x1)
      .def_prop_rw("x2", &QpRectStruct::x2, &QpRectStruct::set_x2)
      .def_prop_rw("y1", &QpRectStruct::y1, &QpRectStruct::set_y1)
      .def_prop_rw("y2", &QpRectStruct::y2, &QpRectStruct::set_y2)
      .def_prop_rw("units", &QpRectStruct::units, &QpRectStruct::set_units)

      .def("__repr__", [](const QpRectStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpRectStruct &self) {
            return QpRectStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const QpRectStruct &self, nb::dict &memo) { return QpRectStruct(self); }
      )
      .def(
          "__eq__",
          [](const QpRectStruct &self, const QpRectStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const QpRectStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D QpRectStruct arrays are not used in structs/routines
  // 2D QpRectStruct arrays are not used in structs/routines
  // 3D QpRectStruct arrays are not used in structs/routines
}

// =============================================================================
// qp_symbol_struct
void init_qp_symbol_struct(nb::module_ &m, nb::class_<QpSymbolStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::string>,
             std::optional<double>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<int>>(),
         nb::arg("type") = nb::none(),
         nb::arg("height") = nb::none(),
         nb::arg("color") = nb::none(),
         nb::arg("fill_pattern") = nb::none(),
         nb::arg("line_width") = nb::none()
  )
      .def_prop_rw("type", &QpSymbolStruct::type, &QpSymbolStruct::set_type)
      .def_prop_rw(
          "height",
          &QpSymbolStruct::height,
          &QpSymbolStruct::set_height,
          "in points (same as text height)"
      )
      .def_prop_rw("color", &QpSymbolStruct::color, &QpSymbolStruct::set_color)
      .def_prop_rw("fill_pattern", &QpSymbolStruct::fill_pattern, &QpSymbolStruct::set_fill_pattern)
      .def_prop_rw("line_width", &QpSymbolStruct::line_width, &QpSymbolStruct::set_line_width)

      .def("__repr__", [](const QpSymbolStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpSymbolStruct &self) {
            return QpSymbolStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const QpSymbolStruct &self, nb::dict &memo) { return QpSymbolStruct(self); }
      )
      .def(
          "__eq__",
          [](const QpSymbolStruct &self, const QpSymbolStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const QpSymbolStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D QpSymbolStruct arrays are not used in structs/routines
  // 2D QpSymbolStruct arrays are not used in structs/routines
  // 3D QpSymbolStruct arrays are not used in structs/routines
}