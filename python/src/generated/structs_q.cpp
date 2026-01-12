#include "pybmad/generated/structs_q.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// qp_axis_struct
void init_qp_axis_struct(py::module& m, py::class_<QpAxisProxy>& cls) {
  cls.def(py::init<>())
      // QpAxisProxy.label (0D_NOT_character -
      .def_property("label", &QpAxisProxy::label, &QpAxisProxy::set_label)
      // QpAxisProxy.min (0D_NOT_real - Axis min/max in data units.
      .def_property("min", &QpAxisProxy::min, &QpAxisProxy::set_min)
      // QpAxisProxy.max (0D_NOT_real - Axis min/max in data units.
      .def_property("max", &QpAxisProxy::max, &QpAxisProxy::set_max)
      // QpAxisProxy.tick_min (0D_NOT_real - Min tick location along axis in data units.
      .def_property(
          "tick_min", &QpAxisProxy::tick_min, &QpAxisProxy::set_tick_min)
      // QpAxisProxy.tick_max (0D_NOT_real - Max tick location along axis in data units.
      .def_property(
          "tick_max", &QpAxisProxy::tick_max, &QpAxisProxy::set_tick_max)
      // QpAxisProxy.eval_min (0D_NOT_real - For general use. Not set by quick_plot.
      .def_property(
          "eval_min", &QpAxisProxy::eval_min, &QpAxisProxy::set_eval_min)
      // QpAxisProxy.eval_max (0D_NOT_real - For general use. Not set by quick_plot.
      .def_property(
          "eval_max", &QpAxisProxy::eval_max, &QpAxisProxy::set_eval_max)
      // QpAxisProxy.dtick (0D_NOT_real - Distance between ticks. In data units. Ticks will be drawn between %min and %max.
      .def_property("dtick", &QpAxisProxy::dtick, &QpAxisProxy::set_dtick)
      // QpAxisProxy.number_offset (0D_NOT_real - Offset from axis line in inches.
      .def_property(
          "number_offset",
          &QpAxisProxy::number_offset,
          &QpAxisProxy::set_number_offset)
      // QpAxisProxy.label_offset (0D_NOT_real - Offset from numbers in inches.
      .def_property(
          "label_offset",
          &QpAxisProxy::label_offset,
          &QpAxisProxy::set_label_offset)
      // QpAxisProxy.major_tick_len (0D_NOT_real - In inches.
      .def_property(
          "major_tick_len",
          &QpAxisProxy::major_tick_len,
          &QpAxisProxy::set_major_tick_len)
      // QpAxisProxy.minor_tick_len (0D_NOT_real - In inches.
      .def_property(
          "minor_tick_len",
          &QpAxisProxy::minor_tick_len,
          &QpAxisProxy::set_minor_tick_len)
      // QpAxisProxy.label_color (0D_NOT_character - Color of the label.
      .def_property(
          "label_color",
          &QpAxisProxy::label_color,
          &QpAxisProxy::set_label_color)
      // QpAxisProxy.major_div (0D_NOT_integer - Actual numbrer of major divisions
      .def_property(
          "major_div", &QpAxisProxy::major_div, &QpAxisProxy::set_major_div)
      // QpAxisProxy.major_div_nominal (0D_NOT_integer - Nominal value.
      .def_property(
          "major_div_nominal",
          &QpAxisProxy::major_div_nominal,
          &QpAxisProxy::set_major_div_nominal)
      // QpAxisProxy.minor_div (0D_NOT_integer - 0 = auto choose.
      .def_property(
          "minor_div", &QpAxisProxy::minor_div, &QpAxisProxy::set_minor_div)
      // QpAxisProxy.minor_div_max (0D_NOT_integer - Max number for auto choose.
      .def_property(
          "minor_div_max",
          &QpAxisProxy::minor_div_max,
          &QpAxisProxy::set_minor_div_max)
      // QpAxisProxy.places (0D_NOT_integer - Number of places after the decimal point to print.
      .def_property("places", &QpAxisProxy::places, &QpAxisProxy::set_places)
      // QpAxisProxy.type (0D_NOT_character - Or 'LOG', or 'CUSTOM'
      .def_property("type", &QpAxisProxy::type, &QpAxisProxy::set_type)
      // QpAxisProxy.bounds (0D_NOT_character - Or 'ZERO_AT_END' or 'ZERO_SYMMETRIC'
      .def_property("bounds", &QpAxisProxy::bounds, &QpAxisProxy::set_bounds)
      // QpAxisProxy.tick_side (0D_NOT_integer - +1 = Draw on the side inside the graph, 0 = both (longer tick), -1 = outside.
      .def_property(
          "tick_side", &QpAxisProxy::tick_side, &QpAxisProxy::set_tick_side)
      // QpAxisProxy.number_side (0D_NOT_integer - +1 = Draw to the side inside the graph, -1 = outside.
      .def_property(
          "number_side",
          &QpAxisProxy::number_side,
          &QpAxisProxy::set_number_side)
      // QpAxisProxy.draw_label (0D_NOT_logical -
      .def_property(
          "draw_label", &QpAxisProxy::draw_label, &QpAxisProxy::set_draw_label)
      // QpAxisProxy.draw_numbers (0D_NOT_logical -
      .def_property(
          "draw_numbers",
          &QpAxisProxy::draw_numbers,
          &QpAxisProxy::set_draw_numbers)

      .def("__repr__", [](const QpAxisProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpAxisProxy& self) {
            return QpAxisProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const QpAxisProxy& self, py::dict& memo) {
            return QpAxisProxy(self);
          })

      ;

  // 1D QpAxisProxy arrays are not used in structs/routines
  // 2D QpAxisProxy arrays are not used in structs/routines
  // 3D QpAxisProxy arrays are not used in structs/routines
}

// =============================================================================
// qp_legend_struct
void init_qp_legend_struct(py::module& m, py::class_<QpLegendProxy>& cls) {
  cls.def(py::init<>())
      // QpLegendProxy.row_spacing (0D_NOT_real - Spacing between rows.
      .def_property(
          "row_spacing",
          &QpLegendProxy::row_spacing,
          &QpLegendProxy::set_row_spacing)
      // QpLegendProxy.line_length (0D_NOT_real - Length of the line in points.
      .def_property(
          "line_length",
          &QpLegendProxy::line_length,
          &QpLegendProxy::set_line_length)
      // QpLegendProxy.text_offset (0D_NOT_real - Horizontal offset in points between the line and the text.
      .def_property(
          "text_offset",
          &QpLegendProxy::text_offset,
          &QpLegendProxy::set_text_offset)
      // QpLegendProxy.draw_line (0D_NOT_logical - Draw lines?
      .def_property(
          "draw_line", &QpLegendProxy::draw_line, &QpLegendProxy::set_draw_line)
      // QpLegendProxy.draw_symbol (0D_NOT_logical - Draw symbols?
      .def_property(
          "draw_symbol",
          &QpLegendProxy::draw_symbol,
          &QpLegendProxy::set_draw_symbol)
      // QpLegendProxy.draw_text (0D_NOT_logical - Draw text?
      .def_property(
          "draw_text", &QpLegendProxy::draw_text, &QpLegendProxy::set_draw_text)

      .def(
          "__repr__", [](const QpLegendProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpLegendProxy& self) {
            return QpLegendProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const QpLegendProxy& self, py::dict& memo) {
            return QpLegendProxy(self);
          })

      ;

  // 1D QpLegendProxy arrays are not used in structs/routines
  // 2D QpLegendProxy arrays are not used in structs/routines
  // 3D QpLegendProxy arrays are not used in structs/routines
}

// =============================================================================
// qp_line_struct
void init_qp_line_struct(py::module& m, py::class_<QpLineProxy>& cls) {
  cls.def(py::init<>())
      // QpLineProxy.width (0D_NOT_integer -
      .def_property("width", &QpLineProxy::width, &QpLineProxy::set_width)
      // QpLineProxy.color (0D_NOT_character -
      .def_property("color", &QpLineProxy::color, &QpLineProxy::set_color)
      // QpLineProxy.pattern (0D_NOT_character -
      .def_property("pattern", &QpLineProxy::pattern, &QpLineProxy::set_pattern)

      .def("__repr__", [](const QpLineProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpLineProxy& self) {
            return QpLineProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const QpLineProxy& self, py::dict& memo) {
            return QpLineProxy(self);
          })

      ;

  // 1D QpLineProxy arrays are not used in structs/routines
  // 2D QpLineProxy arrays are not used in structs/routines
  // 3D QpLineProxy arrays are not used in structs/routines
}

// =============================================================================
// qp_point_struct
void init_qp_point_struct(py::module& m, py::class_<QpPointProxy>& cls) {
  cls.def(py::init<>())
      // QpPointProxy.x (0D_NOT_real -
      .def_property("x", &QpPointProxy::x, &QpPointProxy::set_x)
      // QpPointProxy.y (0D_NOT_real -
      .def_property("y", &QpPointProxy::y, &QpPointProxy::set_y)
      // QpPointProxy.units (0D_NOT_character -
      .def_property("units", &QpPointProxy::units, &QpPointProxy::set_units)

      .def("__repr__", [](const QpPointProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpPointProxy& self) {
            return QpPointProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const QpPointProxy& self, py::dict& memo) {
            return QpPointProxy(self);
          })

      ;

  // 1D QpPointProxy arrays are not used in structs/routines
  // 2D QpPointProxy arrays are not used in structs/routines
  // 3D QpPointProxy arrays are not used in structs/routines
}

// =============================================================================
// qp_rect_struct
void init_qp_rect_struct(py::module& m, py::class_<QpRectProxy>& cls) {
  cls.def(py::init<>())
      // QpRectProxy.x1 (0D_NOT_real -
      .def_property("x1", &QpRectProxy::x1, &QpRectProxy::set_x1)
      // QpRectProxy.x2 (0D_NOT_real -
      .def_property("x2", &QpRectProxy::x2, &QpRectProxy::set_x2)
      // QpRectProxy.y1 (0D_NOT_real -
      .def_property("y1", &QpRectProxy::y1, &QpRectProxy::set_y1)
      // QpRectProxy.y2 (0D_NOT_real -
      .def_property("y2", &QpRectProxy::y2, &QpRectProxy::set_y2)
      // QpRectProxy.units (0D_NOT_character -
      .def_property("units", &QpRectProxy::units, &QpRectProxy::set_units)

      .def("__repr__", [](const QpRectProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpRectProxy& self) {
            return QpRectProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const QpRectProxy& self, py::dict& memo) {
            return QpRectProxy(self);
          })

      ;

  // 1D QpRectProxy arrays are not used in structs/routines
  // 2D QpRectProxy arrays are not used in structs/routines
  // 3D QpRectProxy arrays are not used in structs/routines
}

// =============================================================================
// qp_symbol_struct
void init_qp_symbol_struct(py::module& m, py::class_<QpSymbolProxy>& cls) {
  cls.def(py::init<>())
      // QpSymbolProxy.type (0D_NOT_character -
      .def_property("type", &QpSymbolProxy::type, &QpSymbolProxy::set_type)
      // QpSymbolProxy.height (0D_NOT_real - in points (same as text height)
      .def_property(
          "height", &QpSymbolProxy::height, &QpSymbolProxy::set_height)
      // QpSymbolProxy.color (0D_NOT_character -
      .def_property("color", &QpSymbolProxy::color, &QpSymbolProxy::set_color)
      // QpSymbolProxy.fill_pattern (0D_NOT_character -
      .def_property(
          "fill_pattern",
          &QpSymbolProxy::fill_pattern,
          &QpSymbolProxy::set_fill_pattern)
      // QpSymbolProxy.line_width (0D_NOT_integer -
      .def_property(
          "line_width",
          &QpSymbolProxy::line_width,
          &QpSymbolProxy::set_line_width)

      .def(
          "__repr__", [](const QpSymbolProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const QpSymbolProxy& self) {
            return QpSymbolProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const QpSymbolProxy& self, py::dict& memo) {
            return QpSymbolProxy(self);
          })

      ;

  // 1D QpSymbolProxy arrays are not used in structs/routines
  // 2D QpSymbolProxy arrays are not used in structs/routines
  // 3D QpSymbolProxy arrays are not used in structs/routines
}