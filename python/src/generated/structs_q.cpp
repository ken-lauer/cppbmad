#include "pybmad/generated/structs_q.hpp"

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// qp_axis_struct
void init_qp_axis_struct(py::module &m, py::class_<QpAxisStruct> &cls) {
  cls.def(
         py::init<
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
         py::arg("label") = py::none(),
         py::arg("min") = py::none(),
         py::arg("max") = py::none(),
         py::arg("tick_min") = py::none(),
         py::arg("tick_max") = py::none(),
         py::arg("eval_min") = py::none(),
         py::arg("eval_max") = py::none(),
         py::arg("dtick") = py::none(),
         py::arg("number_offset") = py::none(),
         py::arg("label_offset") = py::none(),
         py::arg("major_tick_len") = py::none(),
         py::arg("minor_tick_len") = py::none(),
         py::arg("label_color") = py::none(),
         py::arg("major_div") = py::none(),
         py::arg("major_div_nominal") = py::none(),
         py::arg("minor_div") = py::none(),
         py::arg("minor_div_max") = py::none(),
         py::arg("places") = py::none(),
         py::arg("type") = py::none(),
         py::arg("bounds") = py::none(),
         py::arg("tick_side") = py::none(),
         py::arg("number_side") = py::none(),
         py::arg("draw_label") = py::none(),
         py::arg("draw_numbers") = py::none()
  )
      .def_property("label", &QpAxisStruct::label, &QpAxisStruct::set_label)
      .def_property(
          "min",
          &QpAxisStruct::min,
          &QpAxisStruct::set_min,
          "Axis min/max in data units."
      )
      .def_property(
          "max",
          &QpAxisStruct::max,
          &QpAxisStruct::set_max,
          "Axis min/max in data units."
      )
      .def_property(
          "tick_min",
          &QpAxisStruct::tick_min,
          &QpAxisStruct::set_tick_min,
          "Min tick location along axis in data units."
      )
      .def_property(
          "tick_max",
          &QpAxisStruct::tick_max,
          &QpAxisStruct::set_tick_max,
          "Max tick location along axis in data units."
      )
      .def_property(
          "eval_min",
          &QpAxisStruct::eval_min,
          &QpAxisStruct::set_eval_min,
          "For general use. Not set by quick_plot."
      )
      .def_property(
          "eval_max",
          &QpAxisStruct::eval_max,
          &QpAxisStruct::set_eval_max,
          "For general use. Not set by quick_plot."
      )
      .def_property(
          "dtick",
          &QpAxisStruct::dtick,
          &QpAxisStruct::set_dtick,
          "Distance between ticks. In data units. Ticks will be drawn between %min and %max."
      )
      .def_property(
          "number_offset",
          &QpAxisStruct::number_offset,
          &QpAxisStruct::set_number_offset,
          "Offset from axis line in inches."
      )
      .def_property(
          "label_offset",
          &QpAxisStruct::label_offset,
          &QpAxisStruct::set_label_offset,
          "Offset from numbers in inches."
      )
      .def_property(
          "major_tick_len",
          &QpAxisStruct::major_tick_len,
          &QpAxisStruct::set_major_tick_len,
          "In inches."
      )
      .def_property(
          "minor_tick_len",
          &QpAxisStruct::minor_tick_len,
          &QpAxisStruct::set_minor_tick_len,
          "In inches."
      )
      .def_property(
          "label_color",
          &QpAxisStruct::label_color,
          &QpAxisStruct::set_label_color,
          "Color of the label."
      )
      .def_property(
          "major_div",
          &QpAxisStruct::major_div,
          &QpAxisStruct::set_major_div,
          "Actual numbrer of major divisions"
      )
      .def_property(
          "major_div_nominal",
          &QpAxisStruct::major_div_nominal,
          &QpAxisStruct::set_major_div_nominal,
          "Nominal value."
      )
      .def_property(
          "minor_div",
          &QpAxisStruct::minor_div,
          &QpAxisStruct::set_minor_div,
          "0 = auto choose."
      )
      .def_property(
          "minor_div_max",
          &QpAxisStruct::minor_div_max,
          &QpAxisStruct::set_minor_div_max,
          "Max number for auto choose."
      )
      .def_property(
          "places",
          &QpAxisStruct::places,
          &QpAxisStruct::set_places,
          "Number of places after the decimal point to print."
      )
      .def_property("type", &QpAxisStruct::type, &QpAxisStruct::set_type, "Or 'LOG', or 'CUSTOM'")
      .def_property(
          "bounds",
          &QpAxisStruct::bounds,
          &QpAxisStruct::set_bounds,
          "Or 'ZERO_AT_END' or 'ZERO_SYMMETRIC'"
      )
      .def_property(
          "tick_side",
          &QpAxisStruct::tick_side,
          &QpAxisStruct::set_tick_side,
          "+1 = Draw on the side inside the graph, 0 = both (longer tick), -1 = outside."
      )
      .def_property(
          "number_side",
          &QpAxisStruct::number_side,
          &QpAxisStruct::set_number_side,
          "+1 = Draw to the side inside the graph, -1 = outside."
      )
      .def_property("draw_label", &QpAxisStruct::draw_label, &QpAxisStruct::set_draw_label)
      .def_property("draw_numbers", &QpAxisStruct::draw_numbers, &QpAxisStruct::set_draw_numbers)

      .def("__repr__", [](const QpAxisStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpAxisStruct &self) {
            return QpAxisStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const QpAxisStruct &self, py::dict &memo) { return QpAxisStruct(self); }
      )

      ;

  // 1D QpAxisStruct arrays are not used in structs/routines
  // 2D QpAxisStruct arrays are not used in structs/routines
  // 3D QpAxisStruct arrays are not used in structs/routines
}

// =============================================================================
// qp_legend_struct
void init_qp_legend_struct(py::module &m, py::class_<QpLegendStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>>(),
         py::arg("row_spacing") = py::none(),
         py::arg("line_length") = py::none(),
         py::arg("text_offset") = py::none(),
         py::arg("draw_line") = py::none(),
         py::arg("draw_symbol") = py::none(),
         py::arg("draw_text") = py::none()
  )
      .def_property(
          "row_spacing",
          &QpLegendStruct::row_spacing,
          &QpLegendStruct::set_row_spacing,
          "Spacing between rows."
      )
      .def_property(
          "line_length",
          &QpLegendStruct::line_length,
          &QpLegendStruct::set_line_length,
          "Length of the line in points."
      )
      .def_property(
          "text_offset",
          &QpLegendStruct::text_offset,
          &QpLegendStruct::set_text_offset,
          "Horizontal offset in points between the line and the text."
      )
      .def_property(
          "draw_line",
          &QpLegendStruct::draw_line,
          &QpLegendStruct::set_draw_line,
          "Draw lines?"
      )
      .def_property(
          "draw_symbol",
          &QpLegendStruct::draw_symbol,
          &QpLegendStruct::set_draw_symbol,
          "Draw symbols?"
      )
      .def_property(
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
          [](const QpLegendStruct &self, py::dict &memo) { return QpLegendStruct(self); }
      )

      ;

  // 1D QpLegendStruct arrays are not used in structs/routines
  // 2D QpLegendStruct arrays are not used in structs/routines
  // 3D QpLegendStruct arrays are not used in structs/routines
}

// =============================================================================
// qp_line_struct
void init_qp_line_struct(py::module &m, py::class_<QpLineStruct> &cls) {
  cls.def(
         py::init<std::optional<int>, std::optional<std::string>, std::optional<std::string>>(),
         py::arg("width") = py::none(),
         py::arg("color") = py::none(),
         py::arg("pattern") = py::none()
  )
      .def_property("width", &QpLineStruct::width, &QpLineStruct::set_width)
      .def_property("color", &QpLineStruct::color, &QpLineStruct::set_color)
      .def_property("pattern", &QpLineStruct::pattern, &QpLineStruct::set_pattern)

      .def("__repr__", [](const QpLineStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpLineStruct &self) {
            return QpLineStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const QpLineStruct &self, py::dict &memo) { return QpLineStruct(self); }
      )

      ;

  // 1D QpLineStruct arrays are not used in structs/routines
  // 2D QpLineStruct arrays are not used in structs/routines
  // 3D QpLineStruct arrays are not used in structs/routines
}

// =============================================================================
// qp_point_struct
void init_qp_point_struct(py::module &m, py::class_<QpPointStruct> &cls) {
  cls.def(
         py::init<std::optional<double>, std::optional<double>, std::optional<std::string>>(),
         py::arg("x") = py::none(),
         py::arg("y") = py::none(),
         py::arg("units") = py::none()
  )
      .def_property("x", &QpPointStruct::x, &QpPointStruct::set_x)
      .def_property("y", &QpPointStruct::y, &QpPointStruct::set_y)
      .def_property("units", &QpPointStruct::units, &QpPointStruct::set_units)

      .def("__repr__", [](const QpPointStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpPointStruct &self) {
            return QpPointStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const QpPointStruct &self, py::dict &memo) { return QpPointStruct(self); }
      )

      ;

  // 1D QpPointStruct arrays are not used in structs/routines
  // 2D QpPointStruct arrays are not used in structs/routines
  // 3D QpPointStruct arrays are not used in structs/routines
}

// =============================================================================
// qp_rect_struct
void init_qp_rect_struct(py::module &m, py::class_<QpRectStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::string>>(),
         py::arg("x1") = py::none(),
         py::arg("x2") = py::none(),
         py::arg("y1") = py::none(),
         py::arg("y2") = py::none(),
         py::arg("units") = py::none()
  )
      .def_property("x1", &QpRectStruct::x1, &QpRectStruct::set_x1)
      .def_property("x2", &QpRectStruct::x2, &QpRectStruct::set_x2)
      .def_property("y1", &QpRectStruct::y1, &QpRectStruct::set_y1)
      .def_property("y2", &QpRectStruct::y2, &QpRectStruct::set_y2)
      .def_property("units", &QpRectStruct::units, &QpRectStruct::set_units)

      .def("__repr__", [](const QpRectStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpRectStruct &self) {
            return QpRectStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const QpRectStruct &self, py::dict &memo) { return QpRectStruct(self); }
      )

      ;

  // 1D QpRectStruct arrays are not used in structs/routines
  // 2D QpRectStruct arrays are not used in structs/routines
  // 3D QpRectStruct arrays are not used in structs/routines
}

// =============================================================================
// qp_symbol_struct
void init_qp_symbol_struct(py::module &m, py::class_<QpSymbolStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<double>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<int>>(),
         py::arg("type") = py::none(),
         py::arg("height") = py::none(),
         py::arg("color") = py::none(),
         py::arg("fill_pattern") = py::none(),
         py::arg("line_width") = py::none()
  )
      .def_property("type", &QpSymbolStruct::type, &QpSymbolStruct::set_type)
      .def_property(
          "height",
          &QpSymbolStruct::height,
          &QpSymbolStruct::set_height,
          "in points (same as text height)"
      )
      .def_property("color", &QpSymbolStruct::color, &QpSymbolStruct::set_color)
      .def_property(
          "fill_pattern",
          &QpSymbolStruct::fill_pattern,
          &QpSymbolStruct::set_fill_pattern
      )
      .def_property("line_width", &QpSymbolStruct::line_width, &QpSymbolStruct::set_line_width)

      .def("__repr__", [](const QpSymbolStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpSymbolStruct &self) {
            return QpSymbolStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const QpSymbolStruct &self, py::dict &memo) { return QpSymbolStruct(self); }
      )

      ;

  // 1D QpSymbolStruct arrays are not used in structs/routines
  // 2D QpSymbolStruct arrays are not used in structs/routines
  // 3D QpSymbolStruct arrays are not used in structs/routines
}