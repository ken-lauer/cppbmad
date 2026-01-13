#include "pybmad/generated/Bmad_routines_p.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyParseCartesianMap python_parse_cartesian_map(
    CartesianMapProxy& ct_map,
    EleProxy& ele,
    LatProxy& lat,
    std::string delim,
    bool delim_found,
    bool err_flag) {
  Bmad::parse_cartesian_map(ct_map, ele, lat, delim, delim_found, err_flag);
  auto py_result{PyParseCartesianMap{delim, delim_found, err_flag}};
  return py_result;
}
PyParseCylindricalMap python_parse_cylindrical_map(
    CylindricalMapProxy& cl_map,
    EleProxy& ele,
    LatProxy& lat,
    std::string delim,
    bool delim_found,
    bool err_flag) {
  Bmad::parse_cylindrical_map(cl_map, ele, lat, delim, delim_found, err_flag);
  auto py_result{PyParseCylindricalMap{delim, delim_found, err_flag}};
  return py_result;
}
PyParseGenGradMap python_parse_gen_grad_map(
    GenGradMapProxy& gg_map,
    EleProxy& ele,
    LatProxy& lat,
    std::string delim,
    bool delim_found,
    bool err_flag) {
  Bmad::parse_gen_grad_map(gg_map, ele, lat, delim, delim_found, err_flag);
  auto py_result{PyParseGenGradMap{delim, delim_found, err_flag}};
  return py_result;
}
PyParseGridField python_parse_grid_field(
    GridFieldProxy& g_field,
    EleProxy& ele,
    LatProxy& lat,
    std::string delim,
    bool delim_found,
    bool err_flag) {
  Bmad::parse_grid_field(g_field, ele, lat, delim, delim_found, err_flag);
  auto py_result{PyParseGridField{delim, delim_found, err_flag}};
  return py_result;
}
PyParseIntegerList python_parse_integer_list(
    std::string err_str,
    LatProxy& lat,
    IntAlloc1D& int_array,
    bool exact_size,
    std::string delim,
    bool delim_found,
    std::optional<std::string> open_delim,
    std::optional<std::string> separator,
    std::optional<std::string> close_delim,
    std::optional<int> default_value,
    bool is_ok) {
  Bmad::parse_integer_list(
      err_str,
      lat,
      int_array,
      exact_size,
      delim,
      delim_found,
      make_opt_ref(open_delim),
      make_opt_ref(separator),
      make_opt_ref(close_delim),
      make_opt_ref(default_value),
      is_ok);
  auto py_result{PyParseIntegerList{
      err_str,
      exact_size,
      delim,
      delim_found,
      open_delim,
      separator,
      close_delim,
      default_value,
      is_ok}};
  return py_result;
}
PyParseIntegerList2 python_parse_integer_list2(
    std::string err_str,
    LatProxy& lat,
    IntAlloc1D& int_array,
    std::optional<int> num_expected = std::nullopt,
    std::optional<std::string> open_delim = std::nullopt,
    std::optional<std::string> separator = std::nullopt,
    std::optional<std::string> close_delim = std::nullopt,
    std::optional<int> default_value = std::nullopt) {
  auto _result = Bmad::parse_integer_list2(
      err_str,
      lat,
      int_array,
      make_opt_ref(num_expected),
      make_opt_ref(open_delim),
      make_opt_ref(separator),
      make_opt_ref(close_delim),
      make_opt_ref(default_value));
  auto py_result{PyParseIntegerList2{
      _result,
      num_expected,
      open_delim,
      separator,
      close_delim,
      default_value}};
  return py_result;
}
PyParseRealList2 python_parse_real_list2(
    LatProxy& lat,
    std::string err_str,
    RealAlloc1D& real_array,
    std::optional<int> num_expected = std::nullopt,
    std::optional<std::string> open_brace = std::nullopt,
    std::optional<std::string> separator = std::nullopt,
    std::optional<std::string> close_brace = std::nullopt,
    std::optional<double> default_value = std::nullopt,
    std::optional<bool> single_value = std::nullopt) {
  auto _result = Bmad::parse_real_list2(
      lat,
      err_str,
      real_array,
      make_opt_ref(num_expected),
      make_opt_ref(open_brace),
      make_opt_ref(separator),
      make_opt_ref(close_brace),
      make_opt_ref(default_value),
      make_opt_ref(single_value));
  auto py_result{PyParseRealList2{
      _result,
      num_expected,
      open_brace,
      separator,
      close_brace,
      default_value,
      single_value}};
  return py_result;
}
PyParserAddConstant python_parser_add_constant(
    std::string word,
    LatProxy& lat,
    bool redef_is_error) {
  Bmad::parser_add_constant(word, lat, redef_is_error);
  auto py_result{PyParserAddConstant{word, redef_is_error}};
  return py_result;
}
PyParserCallCheck python_parser_call_check(
    std::string word,
    int ix_word,
    std::string delim,
    bool delim_found,
    bool call_found,
    std::optional<bool> err_flag = std::nullopt) {
  Bmad::parser_call_check(
      word, ix_word, delim, delim_found, call_found, make_opt_ref(err_flag));
  auto py_result{PyParserCallCheck{
      word, ix_word, delim, delim_found, call_found, err_flag}};
  return py_result;
}
PyParserFastIntegerRead python_parser_fast_integer_read(
    IntAlloc1D& int_vec,
    EleProxy& ele,
    std::string delim_wanted,
    std::string err_str,
    bool is_ok) {
  Bmad::parser_fast_integer_read(int_vec, ele, delim_wanted, err_str, is_ok);
  auto py_result{PyParserFastIntegerRead{delim_wanted, err_str, is_ok}};
  return py_result;
}
PyParserFileStack python_parser_file_stack(
    std::string how,
    std::optional<std::string> file_name_in = std::nullopt,
    std::optional<bool> finished = std::nullopt,
    std::optional<bool> err = std::nullopt,
    std::optional<bool> open_file = std::nullopt,
    std::optional<bool> abort_on_open_error = std::nullopt) {
  Bmad::parser_file_stack(
      how,
      make_opt_ref(file_name_in),
      make_opt_ref(finished),
      make_opt_ref(err),
      make_opt_ref(open_file),
      make_opt_ref(abort_on_open_error));
  auto py_result{PyParserFileStack{
      how, file_name_in, finished, err, open_file, abort_on_open_error}};
  return py_result;
}
PyParserGetInteger python_parser_get_integer(
    int int_val,
    std::string word,
    int ix_word,
    std::string delim,
    bool delim_found,
    bool err,
    std::optional<std::string> str1 = std::nullopt,
    std::optional<std::string> str2 = std::nullopt) {
  Bmad::parser_get_integer(
      int_val,
      word,
      ix_word,
      delim,
      delim_found,
      err,
      make_opt_ref(str1),
      make_opt_ref(str2));
  auto py_result{PyParserGetInteger{
      int_val, word, ix_word, delim, delim_found, err, str1, str2}};
  return py_result;
}
PyParserGetLogical python_parser_get_logical(
    std::string attrib_name,
    bool this_logic,
    std::string ele_name,
    std::string delim,
    bool delim_found,
    bool err) {
  Bmad::parser_get_logical(
      attrib_name, this_logic, ele_name, delim, delim_found, err);
  auto py_result{PyParserGetLogical{
      attrib_name, this_logic, ele_name, delim, delim_found, err}};
  return py_result;
}
PyParserPrintLine python_parser_print_line(LatProxy& lat, bool end_of_file) {
  Bmad::parser_print_line(lat, end_of_file);
  auto py_result{PyParserPrintLine{end_of_file}};
  return py_result;
}
PyParserReadLrWake python_parser_read_lr_wake(
    EleProxy& ele,
    std::string delim,
    bool delim_found,
    bool err_flag) {
  Bmad::parser_read_lr_wake(ele, delim, delim_found, err_flag);
  auto py_result{PyParserReadLrWake{delim, delim_found, err_flag}};
  return py_result;
}
PyParserReadSrWake python_parser_read_sr_wake(
    EleProxy& ele,
    std::string delim,
    bool delim_found,
    bool err_flag) {
  Bmad::parser_read_sr_wake(ele, delim, delim_found, err_flag);
  auto py_result{PyParserReadSrWake{delim, delim_found, err_flag}};
  return py_result;
}
PyParticleIsMovingBackwards python_particle_is_moving_backwards(
    CoordProxy& orbit,
    bool is_moving_backwards) {
  Bmad::particle_is_moving_backwards(orbit, is_moving_backwards);
  auto py_result{PyParticleIsMovingBackwards{is_moving_backwards}};
  return py_result;
}
PyParticleIsMovingForward python_particle_is_moving_forward(
    CoordProxy& orbit,
    std::optional<int> dir,
    bool is_moving_forward) {
  Bmad::particle_is_moving_forward(orbit, dir, is_moving_forward);
  auto py_result{PyParticleIsMovingForward{is_moving_forward}};
  return py_result;
}
PyParticleRfTime python_particle_rf_time(
    CoordProxy& orbit,
    EleProxy& ele,
    std::optional<bool> reference_active_edge,
    std::optional<double> s_rel,
    std::optional<bool> time_coords,
    std::optional<double> rf_freq,
    std::optional<bool> abs_time,
    long double time) {
  Bmad::particle_rf_time(
      orbit,
      ele,
      reference_active_edge,
      s_rel,
      time_coords,
      rf_freq,
      abs_time,
      time);
  auto py_result{PyParticleRfTime{time}};
  return py_result;
}
PyPatchFlipsPropagationDirection python_patch_flips_propagation_direction(
    double x_pitch,
    double y_pitch,
    bool is_flip) {
  Bmad::patch_flips_propagation_direction(x_pitch, y_pitch, is_flip);
  auto py_result{PyPatchFlipsPropagationDirection{is_flip}};
  return py_result;
}
PyPatchLength python_patch_length(
    EleProxy& patch,
    std::optional<int> ref_coords,
    double length) {
  Bmad::patch_length(patch, ref_coords, length);
  auto py_result{PyPatchLength{length}};
  return py_result;
}
PyPhotonAddToDetectorStatistics python_photon_add_to_detector_statistics(
    CoordProxy& orbit0,
    CoordProxy& orbit,
    EleProxy& ele,
    std::optional<int> ix_pt = std::nullopt,
    std::optional<int> iy_pt = std::nullopt,
    optional_ref<PixelPtProxy> pixel_pt = std::nullopt) {
  Bmad::photon_add_to_detector_statistics(
      orbit0, orbit, ele, make_opt_ref(ix_pt), make_opt_ref(iy_pt), pixel_pt);
  auto py_result{PyPhotonAddToDetectorStatistics{ix_pt, iy_pt}};
  return py_result;
}
PyPhotonTargetCornerCalc python_photon_target_corner_calc(
    EleProxy& aperture_ele,
    double x_lim,
    double y_lim,
    double z_lim,
    EleProxy& source_ele) {
  auto _result = Bmad::photon_target_corner_calc(
      aperture_ele, x_lim, y_lim, z_lim, source_ele);
  auto py_result{PyPhotonTargetCornerCalc{_result, x_lim, y_lim, z_lim}};
  return py_result;
}
PyPhysicalEleEnd python_physical_ele_end(
    int track_end,
    CoordProxy& orbit,
    int ele_orientation,
    std::optional<bool> return_stream_end,
    int physical_end) {
  Bmad::physical_ele_end(
      track_end, orbit, ele_orientation, return_stream_end, physical_end);
  auto py_result{PyPhysicalEleEnd{physical_end}};
  return py_result;
}
PyPointerToSurfaceDisplacementPt python_pointer_to_surface_displacement_pt(
    EleProxy& ele,
    bool nearest,
    double x,
    double y,
    std::optional<int> ix = std::nullopt,
    std::optional<int> iy = std::nullopt,
    std::optional<bool> extend_grid = std::nullopt,
    std::optional<double> xx = std::nullopt,
    std::optional<double> yy = std::nullopt) {
  auto _result = Bmad::pointer_to_surface_displacement_pt(
      ele,
      nearest,
      x,
      y,
      make_opt_ref(ix),
      make_opt_ref(iy),
      extend_grid,
      make_opt_ref(xx),
      make_opt_ref(yy));
  auto py_result{
      PyPointerToSurfaceDisplacementPt{_result, x, y, ix, iy, xx, yy}};
  return py_result;
}
PyPointerToSurfaceSegmentedPt python_pointer_to_surface_segmented_pt(
    EleProxy& ele,
    bool nearest,
    double x,
    double y,
    std::optional<int> ix = std::nullopt,
    std::optional<int> iy = std::nullopt,
    std::optional<bool> extend_grid = std::nullopt,
    std::optional<double> xx = std::nullopt,
    std::optional<double> yy = std::nullopt) {
  auto _result = Bmad::pointer_to_surface_segmented_pt(
      ele,
      nearest,
      x,
      y,
      make_opt_ref(ix),
      make_opt_ref(iy),
      extend_grid,
      make_opt_ref(xx),
      make_opt_ref(yy));
  auto py_result{PyPointerToSurfaceSegmentedPt{_result, x, y, ix, iy, xx, yy}};
  return py_result;
}

void init_Bmad_routines_p(py::module& m) {
  py::class_<PyParseCartesianMap, std::unique_ptr<PyParseCartesianMap>>(
      m, "ParseCartesianMap", "parse_cartesian_map return type")
      .def_readonly("delim", &PyParseCartesianMap::delim)
      .def_readonly("delim_found", &PyParseCartesianMap::delim_found)
      .def_readonly("err_flag", &PyParseCartesianMap::err_flag)
      .def("__len__", [](const PyParseCartesianMap&) { return 3; })
      .def(
          "__getitem__", [](const PyParseCartesianMap& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.delim);
            if (i == 1)
              return py::cast(s.delim_found);
            if (i == 2)
              return py::cast(s.err_flag);
            throw py::index_error();
          });
  m.def(
      "parse_cartesian_map",
      &python_parse_cartesian_map,
      py::arg("ct_map"),
      py::arg("ele"),
      py::arg("lat"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("err_flag"),
      R"""(Subroutine parse_cartesian_map (ct_map, ele, lat, delim, delim_found, err_flag)

  Subroutine to parse a "cartesian_map = {}" construct

  This subroutine is used by bmad_parser and bmad_parser2.
  This subroutine is private to bmad_parser_mod.
  This must read in:
  {type = ,
     dr = ,
     r0 = ,
     pt(i,j,k) = ( (ex_re, ex_im), .... (bz_re, bz_im) )
     .
     .
     . ) },

  )""");
  py::class_<PyParseCylindricalMap, std::unique_ptr<PyParseCylindricalMap>>(
      m, "ParseCylindricalMap", "parse_cylindrical_map return type")
      .def_readonly("delim", &PyParseCylindricalMap::delim)
      .def_readonly("delim_found", &PyParseCylindricalMap::delim_found)
      .def_readonly("err_flag", &PyParseCylindricalMap::err_flag)
      .def("__len__", [](const PyParseCylindricalMap&) { return 3; })
      .def(
          "__getitem__",
          [](const PyParseCylindricalMap& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.delim);
            if (i == 1)
              return py::cast(s.delim_found);
            if (i == 2)
              return py::cast(s.err_flag);
            throw py::index_error();
          });
  m.def(
      "parse_cylindrical_map",
      &python_parse_cylindrical_map,
      py::arg("cl_map"),
      py::arg("ele"),
      py::arg("lat"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("err_flag"),
      R"""(Parameters
  ----------
  cl_map : 
  ele : 
  lat : 
  delim : 
  delim_found : 
  err_flag : 
  )""");
  py::class_<PyParseGenGradMap, std::unique_ptr<PyParseGenGradMap>>(
      m, "ParseGenGradMap", "parse_gen_grad_map return type")
      .def_readonly("delim", &PyParseGenGradMap::delim)
      .def_readonly("delim_found", &PyParseGenGradMap::delim_found)
      .def_readonly("err_flag", &PyParseGenGradMap::err_flag)
      .def("__len__", [](const PyParseGenGradMap&) { return 3; })
      .def("__getitem__", [](const PyParseGenGradMap& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.delim);
        if (i == 1)
          return py::cast(s.delim_found);
        if (i == 2)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "parse_gen_grad_map",
      &python_parse_gen_grad_map,
      py::arg("gg_map"),
      py::arg("ele"),
      py::arg("lat"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("err_flag"),
      R"""(Subroutine parse_gen_grad_map (gg_map, ele, lat, delim, delim_found, err_flag)

  Subroutine to parse a "gen_grad_map = {}" construct

  )""");
  py::class_<PyParseGridField, std::unique_ptr<PyParseGridField>>(
      m, "ParseGridField", "parse_grid_field return type")
      .def_readonly("delim", &PyParseGridField::delim)
      .def_readonly("delim_found", &PyParseGridField::delim_found)
      .def_readonly("err_flag", &PyParseGridField::err_flag)
      .def("__len__", [](const PyParseGridField&) { return 3; })
      .def("__getitem__", [](const PyParseGridField& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.delim);
        if (i == 1)
          return py::cast(s.delim_found);
        if (i == 2)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "parse_grid_field",
      &python_parse_grid_field,
      py::arg("g_field"),
      py::arg("ele"),
      py::arg("lat"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("err_flag"),
      R"""(Parameters
  ----------
  g_field : 
  ele : 
  lat : 
  delim : 
  delim_found : 
  err_flag : 
  )""");
  py::class_<PyParseIntegerList, std::unique_ptr<PyParseIntegerList>>(
      m, "ParseIntegerList", "parse_integer_list return type")
      .def_readonly("err_str", &PyParseIntegerList::err_str)
      .def_readonly("exact_size", &PyParseIntegerList::exact_size)
      .def_readonly("delim", &PyParseIntegerList::delim)
      .def_readonly("delim_found", &PyParseIntegerList::delim_found)
      .def_readonly("open_delim", &PyParseIntegerList::open_delim)
      .def_readonly("separator", &PyParseIntegerList::separator)
      .def_readonly("close_delim", &PyParseIntegerList::close_delim)
      .def_readonly("default_value", &PyParseIntegerList::default_value)
      .def_readonly("is_ok", &PyParseIntegerList::is_ok)
      .def("__len__", [](const PyParseIntegerList&) { return 9; })
      .def("__getitem__", [](const PyParseIntegerList& s, int i) -> py::object {
        if (i < 0)
          i += 9;
        if (i == 0)
          return py::cast(s.err_str);
        if (i == 1)
          return py::cast(s.exact_size);
        if (i == 2)
          return py::cast(s.delim);
        if (i == 3)
          return py::cast(s.delim_found);
        if (i == 4)
          return py::cast(s.open_delim);
        if (i == 5)
          return py::cast(s.separator);
        if (i == 6)
          return py::cast(s.close_delim);
        if (i == 7)
          return py::cast(s.default_value);
        if (i == 8)
          return py::cast(s.is_ok);
        throw py::index_error();
      });
  m.def(
      "parse_integer_list",
      &python_parse_integer_list,
      py::arg("err_str"),
      py::arg("lat"),
      py::arg("int_array"),
      py::arg("exact_size"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("open_delim") = py::none(),
      py::arg("separator") = py::none(),
      py::arg("close_delim") = py::none(),
      py::arg("default_value") = py::none(),
      py::arg("is_ok"),
      R"""(Function parse_integer_list (err_str, lat, int_array, exact_size, delim, delim_found, open_delim,
                                        separator, close_delim, default_value) result (is_ok)

  Routine to parse a list of integers of the form:
     open_delim integer_1 separator integer_2 . . . close_delim
  Example:   "(1.2, 2.3, 4.4, 8.5)"

  Similar to parse_integer_list2 except does not use allocatable array.
  See parse_integer_list2 for more details

  )""");
  py::class_<PyParseIntegerList2, std::unique_ptr<PyParseIntegerList2>>(
      m, "ParseIntegerList2", "parse_integer_list2 return type")
      .def_readonly("num_found", &PyParseIntegerList2::num_found)
      .def_readonly("delim", &PyParseIntegerList2::delim)
      .def_readonly("delim_found", &PyParseIntegerList2::delim_found)
      .def_readonly("is_ok", &PyParseIntegerList2::is_ok)
      .def_readonly("num_expected", &PyParseIntegerList2::num_expected)
      .def_readonly("open_delim", &PyParseIntegerList2::open_delim)
      .def_readonly("separator", &PyParseIntegerList2::separator)
      .def_readonly("close_delim", &PyParseIntegerList2::close_delim)
      .def_readonly("default_value", &PyParseIntegerList2::default_value)
      .def("__len__", [](const PyParseIntegerList2&) { return 9; })
      .def(
          "__getitem__", [](const PyParseIntegerList2& s, int i) -> py::object {
            if (i < 0)
              i += 9;
            if (i == 0)
              return py::cast(s.num_found);
            if (i == 1)
              return py::cast(s.delim);
            if (i == 2)
              return py::cast(s.delim_found);
            if (i == 3)
              return py::cast(s.is_ok);
            if (i == 4)
              return py::cast(s.num_expected);
            if (i == 5)
              return py::cast(s.open_delim);
            if (i == 6)
              return py::cast(s.separator);
            if (i == 7)
              return py::cast(s.close_delim);
            if (i == 8)
              return py::cast(s.default_value);
            throw py::index_error();
          });
  m.def(
      "parse_integer_list2",
      &python_parse_integer_list2,
      py::arg("err_str"),
      py::arg("lat"),
      py::arg("int_array"),
      py::arg("num_expected") = py::none(),
      py::arg("open_delim") = py::none(),
      py::arg("separator") = py::none(),
      py::arg("close_delim") = py::none(),
      py::arg("default_value") = py::none(),
      R"""(Function parse_integer_list2 (err_str, lat, int_array, num_found, delim, delim_found, num_expected,
                                         open_delim, separator, close_delim, default_value) result (is_ok)

  Routine to parse a list of integers of the form
     open_delim integer_1 separator integer_2 . . . close_delim
  Example:   (1, 2, 4, 8)

  Parameters
  ----------
  err_str : unknown
      Error string to print if there is an error.
  lat : LatStruct
      lattice
  int_array : int
      the array to be read in Optional: num_expected = 1     -- integer: number of expected arguments. Used to
      initialize int_array. open_delim   = '('   -- character(1): opening delimeter. separator    = ','   --
      character(1): separating character. close_delim  = ')'   -- character(1): closing delimeter. default_value
      = 0    -- real(rp): inital assignment of int_array elements.
      This parameter is an input/output and is modified in-place. As an output: Array of values.

  Returns
  -------
  is_ok : bool
      Set True if everything is ok.
  num_found : int
      number of elements.
  delim : unknown
      Delimiter found where the parsing of the input line stops.
  delim_found : bool
      Delimiter found? False if end of input command.
  )""");
  py::class_<Bmad::ParseRealList, std::unique_ptr<Bmad::ParseRealList>>(
      m, "ParseRealList", "parse_real_list return type")
      .def_readonly("real_array", &Bmad::ParseRealList::real_array)
      .def_readonly("delim", &Bmad::ParseRealList::delim)
      .def_readonly("delim_found", &Bmad::ParseRealList::delim_found)
      .def_readonly("num_found", &Bmad::ParseRealList::num_found)
      .def_readonly("is_ok", &Bmad::ParseRealList::is_ok)
      .def("__len__", [](const Bmad::ParseRealList&) { return 5; })
      .def(
          "__getitem__", [](const Bmad::ParseRealList& s, int i) -> py::object {
            if (i < 0)
              i += 5;
            if (i == 0)
              return py::cast(s.real_array);
            if (i == 1)
              return py::cast(s.delim);
            if (i == 2)
              return py::cast(s.delim_found);
            if (i == 3)
              return py::cast(s.num_found);
            if (i == 4)
              return py::cast(s.is_ok);
            throw py::index_error();
          });
  m.def(
      "parse_real_list",
      &Bmad::parse_real_list,
      py::arg("lat"),
      py::arg("err_str"),
      py::arg("exact_size"),
      py::arg("open_delim") = py::none(),
      py::arg("separator") = py::none(),
      py::arg("close_delim") = py::none(),
      py::arg("default_value") = py::none(),
      R"""(Function parse_real_list (lat, err_str, real_array, exact_size, delim, delim_found, open_delim,
                                 separator, close_delim, default_value, num_found) result (is_ok)

  Routine to parse a list of reals of the form:
     open_delim real_1 separator real_2 . . . close_delim
  Example:   "(1.2, 2.3, 4.4, 8.5)"

  Similar to parse_real_list2 except does not use allocatable array.
  Also see: parse_real_matrix.

  Parameters
  ----------
  lat : LatStruct
      Lattice
  err_str : unknown
      Error string to print if there is an error.
  exact_size : 
  open_delim : 
  separator : 
  close_delim : 
  default_value : 

  Returns
  -------
  real_array : 
  delim : 
  delim_found : 
  num_found : 

  Notes
  -----
  Related routines:
  parse_real_matrix.
  )""");
  py::class_<PyParseRealList2, std::unique_ptr<PyParseRealList2>>(
      m, "ParseRealList2", "parse_real_list2 return type")
      .def_readonly("num_found", &PyParseRealList2::num_found)
      .def_readonly("delim", &PyParseRealList2::delim)
      .def_readonly("delim_found", &PyParseRealList2::delim_found)
      .def_readonly("is_ok", &PyParseRealList2::is_ok)
      .def_readonly("num_expected", &PyParseRealList2::num_expected)
      .def_readonly("open_brace", &PyParseRealList2::open_brace)
      .def_readonly("separator", &PyParseRealList2::separator)
      .def_readonly("close_brace", &PyParseRealList2::close_brace)
      .def_readonly("default_value", &PyParseRealList2::default_value)
      .def_readonly("single_value", &PyParseRealList2::single_value)
      .def("__len__", [](const PyParseRealList2&) { return 10; })
      .def("__getitem__", [](const PyParseRealList2& s, int i) -> py::object {
        if (i < 0)
          i += 10;
        if (i == 0)
          return py::cast(s.num_found);
        if (i == 1)
          return py::cast(s.delim);
        if (i == 2)
          return py::cast(s.delim_found);
        if (i == 3)
          return py::cast(s.is_ok);
        if (i == 4)
          return py::cast(s.num_expected);
        if (i == 5)
          return py::cast(s.open_brace);
        if (i == 6)
          return py::cast(s.separator);
        if (i == 7)
          return py::cast(s.close_brace);
        if (i == 8)
          return py::cast(s.default_value);
        if (i == 9)
          return py::cast(s.single_value);
        throw py::index_error();
      });
  m.def(
      "parse_real_list2",
      &python_parse_real_list2,
      py::arg("lat"),
      py::arg("err_str"),
      py::arg("real_array"),
      py::arg("num_expected") = py::none(),
      py::arg("open_brace") = py::none(),
      py::arg("separator") = py::none(),
      py::arg("close_brace") = py::none(),
      py::arg("default_value") = py::none(),
      py::arg("single_value") = py::none(),
      R"""(Function parse_real_list2 (lat, err_str, real_array, num_found, delim, delim_found, num_expected,
                             open_delim, separator, close_delim, default_value, single_value) result (is_ok)

  Routine to parse a list of reals of the form:
     open_brace real_1 separator real_2 . . . close_brace
  Example:   "(1.2, 2.3, 4.4, 8.5)"

  Parameters
  ----------
  lat : LatStruct
      lattice
  err_str : unknown
      Error string to print if there is an error.
  real_array : float
      the array to be read in num_expected = 10       -- integer, optional: number of expected arguments Used to
      initialize real_array open_brace   = '('      -- character(1), optional: opening delimeter. separator    =
      ','      -- character(1), optional: separating character close_brace  = ')'      -- character(1),
      optional: closing delimeter default_value = 0.0_rp  -- real(rp), optional: inital assignment of real_array
      elements. single_value = False    -- logical, optional: If true then an array with a single value and no
      braces is accepted.
      This parameter is an input/output and is modified in-place. As an output: Array of values

  Returns
  -------
  is_ok : bool
      Set True if everything is ok
  num_found : int
      number of elements
  delim : unknown
      Delimiter found where the parsing of the input line stops.
  delim_found : bool
      Stopping delimiter found? False if end of input command.

  Notes
  -----
  Related routines:
  pase_real_list parse_real_matrix.
  )""");
  py::class_<PyParserAddConstant, std::unique_ptr<PyParserAddConstant>>(
      m, "ParserAddConstant", "parser_add_constant return type")
      .def_readonly("word", &PyParserAddConstant::word)
      .def_readonly("redef_is_error", &PyParserAddConstant::redef_is_error)
      .def("__len__", [](const PyParserAddConstant&) { return 2; })
      .def(
          "__getitem__", [](const PyParserAddConstant& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.word);
            if (i == 1)
              return py::cast(s.redef_is_error);
            throw py::index_error();
          });
  m.def(
      "parser_add_constant",
      &python_parser_add_constant,
      py::arg("word"),
      py::arg("lat"),
      py::arg("redef_is_error"),
      R"""(Parameters
  ----------
  word : 
  lat : 
  redef_is_error : 
  )""");
  py::class_<PyParserCallCheck, std::unique_ptr<PyParserCallCheck>>(
      m, "ParserCallCheck", "parser_call_check return type")
      .def_readonly("word", &PyParserCallCheck::word)
      .def_readonly("ix_word", &PyParserCallCheck::ix_word)
      .def_readonly("delim", &PyParserCallCheck::delim)
      .def_readonly("delim_found", &PyParserCallCheck::delim_found)
      .def_readonly("call_found", &PyParserCallCheck::call_found)
      .def_readonly("err_flag", &PyParserCallCheck::err_flag)
      .def("__len__", [](const PyParserCallCheck&) { return 6; })
      .def("__getitem__", [](const PyParserCallCheck& s, int i) -> py::object {
        if (i < 0)
          i += 6;
        if (i == 0)
          return py::cast(s.word);
        if (i == 1)
          return py::cast(s.ix_word);
        if (i == 2)
          return py::cast(s.delim);
        if (i == 3)
          return py::cast(s.delim_found);
        if (i == 4)
          return py::cast(s.call_found);
        if (i == 5)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "parser_call_check",
      &python_parser_call_check,
      py::arg("word"),
      py::arg("ix_word"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("call_found"),
      py::arg("err_flag") = py::none(),
      R"""(Subroutine parser_call_check(word, ix_word, delim, delim_found, call_found, err_flag))

  Routine to check if there is a "call::XXX" construct in the input stream.

  )""");
  py::class_<
      Bmad::ParserFastComplexRead,
      std::unique_ptr<Bmad::ParserFastComplexRead>>(
      m, "ParserFastComplexRead", "parser_fast_complex_read return type")
      .def_readonly("cmplx_vec", &Bmad::ParserFastComplexRead::cmplx_vec)
      .def_readonly("delim", &Bmad::ParserFastComplexRead::delim)
      .def_readonly("is_ok", &Bmad::ParserFastComplexRead::is_ok)
      .def("__len__", [](const Bmad::ParserFastComplexRead&) { return 3; })
      .def(
          "__getitem__",
          [](const Bmad::ParserFastComplexRead& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.cmplx_vec);
            if (i == 1)
              return py::cast(s.delim);
            if (i == 2)
              return py::cast(s.is_ok);
            throw py::index_error();
          });
  m.def(
      "parser_fast_complex_read",
      &Bmad::parser_fast_complex_read,
      py::arg("ele"),
      py::arg("err_str"),
      R"""(Function parser_fast_complex_read (cmplx_vec, ele, delim, err_str)  result (is_ok)

  Routine to read an array of complex numbers.

  This routine assumes that the array values are pure numbers in the form "<re>" or "(<re> <im>)"
  where <re> and <im> are real numbers (not expressions) and there are no commas except possibly
  at the end of the array.

  Parameters
  ----------
  ele : EleStruct
      Lattice element associated with the array. Used for error messages.
  err_str : unknown
      String used when printing error messages identifying where in the lattice file the error is occuring.

  Returns
  -------
  cmplx_vec : complex
      Complex vector.
  delim : unknown
      Delimitor at end of array. Must be "," or "}"
  is_ok : bool
      True if everything OK. False otherwise.
  )""");
  py::class_<PyParserFastIntegerRead, std::unique_ptr<PyParserFastIntegerRead>>(
      m, "ParserFastIntegerRead", "parser_fast_integer_read return type")
      .def_readonly("delim_wanted", &PyParserFastIntegerRead::delim_wanted)
      .def_readonly("err_str", &PyParserFastIntegerRead::err_str)
      .def_readonly("is_ok", &PyParserFastIntegerRead::is_ok)
      .def("__len__", [](const PyParserFastIntegerRead&) { return 3; })
      .def(
          "__getitem__",
          [](const PyParserFastIntegerRead& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.delim_wanted);
            if (i == 1)
              return py::cast(s.err_str);
            if (i == 2)
              return py::cast(s.is_ok);
            throw py::index_error();
          });
  m.def(
      "parser_fast_integer_read",
      &python_parser_fast_integer_read,
      py::arg("int_vec"),
      py::arg("ele"),
      py::arg("delim_wanted"),
      py::arg("err_str"),
      py::arg("is_ok"),
      R"""(Function parser_fast_integer_read (int_vec, ele, delim_wanted, err_str)  result (is_ok)


  Returns
  -------
  is_ok
  )""");
  py::class_<
      Bmad::ParserFastRealRead,
      std::unique_ptr<Bmad::ParserFastRealRead>>(
      m, "ParserFastRealRead", "parser_fast_real_read return type")
      .def_readonly("real_vec", &Bmad::ParserFastRealRead::real_vec)
      .def_readonly("delim", &Bmad::ParserFastRealRead::delim)
      .def_readonly("n_real", &Bmad::ParserFastRealRead::n_real)
      .def_readonly("is_ok", &Bmad::ParserFastRealRead::is_ok)
      .def("__len__", [](const Bmad::ParserFastRealRead&) { return 4; })
      .def(
          "__getitem__",
          [](const Bmad::ParserFastRealRead& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.real_vec);
            if (i == 1)
              return py::cast(s.delim);
            if (i == 2)
              return py::cast(s.n_real);
            if (i == 3)
              return py::cast(s.is_ok);
            throw py::index_error();
          });
  m.def(
      "parser_fast_real_read",
      &Bmad::parser_fast_real_read,
      py::arg("ele"),
      py::arg("end_delims"),
      py::arg("err_str"),
      py::arg("exact_size") = py::none(),
      R"""(Function parser_fast_real_read (real_vec, ele, end_delims, delim, err_str, exact_size, n_real)  result (is_ok)

  Routine to read an array of real numbers.

  This routine assumes that the array values are pure numbers in the form "<re1> <re2> ...,"
  where <re1>, <re2>, etc. are real numbers (not expressions) and there are no commas except possibly,
  at the end of the array.

  Note: if end_delim is "," and next character is a delim but not ",", the next character is taken as the delim.

  Parameters
  ----------
  ele : EleStruct
      Lattice element associated with the array. Used for error messages.
  end_delims : unknown
      List of possible ending delimitors.
  err_str : unknown
      String used when printing error messages identifying where in the lattice file the error is occuring.
  exact_size : bool, optional
      If True (default), number of values must match real_vec size.

  Returns
  -------
  real_vec : complex
      Real vector.
  delim : unknown
      Delimitor at end of array.
  is_ok : bool
      True if everything OK. False otherwise.
  n_real : int
      Number of elements found.
  )""");
  py::class_<PyParserFileStack, std::unique_ptr<PyParserFileStack>>(
      m, "ParserFileStack", "parser_file_stack return type")
      .def_readonly("how", &PyParserFileStack::how)
      .def_readonly("file_name_in", &PyParserFileStack::file_name_in)
      .def_readonly("finished", &PyParserFileStack::finished)
      .def_readonly("err", &PyParserFileStack::err)
      .def_readonly("open_file", &PyParserFileStack::open_file)
      .def_readonly(
          "abort_on_open_error", &PyParserFileStack::abort_on_open_error)
      .def("__len__", [](const PyParserFileStack&) { return 6; })
      .def("__getitem__", [](const PyParserFileStack& s, int i) -> py::object {
        if (i < 0)
          i += 6;
        if (i == 0)
          return py::cast(s.how);
        if (i == 1)
          return py::cast(s.file_name_in);
        if (i == 2)
          return py::cast(s.finished);
        if (i == 3)
          return py::cast(s.err);
        if (i == 4)
          return py::cast(s.open_file);
        if (i == 5)
          return py::cast(s.abort_on_open_error);
        throw py::index_error();
      });
  m.def(
      "parser_file_stack",
      &python_parser_file_stack,
      py::arg("how"),
      py::arg("file_name_in") = py::none(),
      py::arg("finished") = py::none(),
      py::arg("err") = py::none(),
      py::arg("open_file") = py::none(),
      py::arg("abort_on_open_error") = py::none(),
      R"""(Subroutine parser_file_stack (how, file_name_in, finished, err, open_file, abort_on_open_error)

  Subroutine to keep track of the files that are opened for reading.
  This subroutine is used by bmad_parser and bmad_parser2.
  This subroutine is not intended for general use.

  )""");
  py::class_<PyParserGetInteger, std::unique_ptr<PyParserGetInteger>>(
      m, "ParserGetInteger", "parser_get_integer return type")
      .def_readonly("int_val", &PyParserGetInteger::int_val)
      .def_readonly("word", &PyParserGetInteger::word)
      .def_readonly("ix_word", &PyParserGetInteger::ix_word)
      .def_readonly("delim", &PyParserGetInteger::delim)
      .def_readonly("delim_found", &PyParserGetInteger::delim_found)
      .def_readonly("err", &PyParserGetInteger::err)
      .def_readonly("str1", &PyParserGetInteger::str1)
      .def_readonly("str2", &PyParserGetInteger::str2)
      .def("__len__", [](const PyParserGetInteger&) { return 8; })
      .def("__getitem__", [](const PyParserGetInteger& s, int i) -> py::object {
        if (i < 0)
          i += 8;
        if (i == 0)
          return py::cast(s.int_val);
        if (i == 1)
          return py::cast(s.word);
        if (i == 2)
          return py::cast(s.ix_word);
        if (i == 3)
          return py::cast(s.delim);
        if (i == 4)
          return py::cast(s.delim_found);
        if (i == 5)
          return py::cast(s.err);
        if (i == 6)
          return py::cast(s.str1);
        if (i == 7)
          return py::cast(s.str2);
        throw py::index_error();
      });
  m.def(
      "parser_get_integer",
      &python_parser_get_integer,
      py::arg("int_val"),
      py::arg("word"),
      py::arg("ix_word"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("err"),
      py::arg("str1") = py::none(),
      py::arg("str2") = py::none(),
      R"""(Parameters
  ----------
  int_val : 
  word : 
  ix_word : 
  delim : 
  delim_found : 
  err : 
  str1 : 
  str2 : 
  )""");
  py::class_<PyParserGetLogical, std::unique_ptr<PyParserGetLogical>>(
      m, "ParserGetLogical", "parser_get_logical return type")
      .def_readonly("attrib_name", &PyParserGetLogical::attrib_name)
      .def_readonly("this_logic", &PyParserGetLogical::this_logic)
      .def_readonly("ele_name", &PyParserGetLogical::ele_name)
      .def_readonly("delim", &PyParserGetLogical::delim)
      .def_readonly("delim_found", &PyParserGetLogical::delim_found)
      .def_readonly("err", &PyParserGetLogical::err)
      .def("__len__", [](const PyParserGetLogical&) { return 6; })
      .def("__getitem__", [](const PyParserGetLogical& s, int i) -> py::object {
        if (i < 0)
          i += 6;
        if (i == 0)
          return py::cast(s.attrib_name);
        if (i == 1)
          return py::cast(s.this_logic);
        if (i == 2)
          return py::cast(s.ele_name);
        if (i == 3)
          return py::cast(s.delim);
        if (i == 4)
          return py::cast(s.delim_found);
        if (i == 5)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "parser_get_logical",
      &python_parser_get_logical,
      py::arg("attrib_name"),
      py::arg("this_logic"),
      py::arg("ele_name"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("err"),
      R"""(Parameters
  ----------
  attrib_name : 
  this_logic : 
  ele_name : 
  delim : 
  delim_found : 
  err : 
  )""");
  m.def(
      "parser_identify_fork_to_element",
      &Bmad::parser_identify_fork_to_element,
      py::arg("lat"),
      R"""(Subroutine parser_identify_fork_to_element (lat)

  Routine to identify the elements the forks in a lattice are branching to.

  This subroutine is used by bmad_parser and bmad_parser2.
  This subroutine is not intended for general use.

  )""");
  m.def(
      "parser_init_custom_elements",
      &Bmad::parser_init_custom_elements,
      py::arg("lat"),
      R"""(Subroutine parser_init_custom_elements (lat)

  )""");
  py::class_<PyParserPrintLine, std::unique_ptr<PyParserPrintLine>>(
      m, "ParserPrintLine", "parser_print_line return type")
      .def_readonly("end_of_file", &PyParserPrintLine::end_of_file)
      .def("__len__", [](const PyParserPrintLine&) { return 1; })
      .def("__getitem__", [](const PyParserPrintLine& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.end_of_file);
        throw py::index_error();
      });
  m.def(
      "parser_print_line",
      &python_parser_print_line,
      py::arg("lat"),
      py::arg("end_of_file"),
      R"""(Subroutine parser_print_line(end_of_file)

  This routine is called when a print statement is found in the lattice file.

  )""");
  py::class_<PyParserReadLrWake, std::unique_ptr<PyParserReadLrWake>>(
      m, "ParserReadLrWake", "parser_read_lr_wake return type")
      .def_readonly("delim", &PyParserReadLrWake::delim)
      .def_readonly("delim_found", &PyParserReadLrWake::delim_found)
      .def_readonly("err_flag", &PyParserReadLrWake::err_flag)
      .def("__len__", [](const PyParserReadLrWake&) { return 3; })
      .def("__getitem__", [](const PyParserReadLrWake& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.delim);
        if (i == 1)
          return py::cast(s.delim_found);
        if (i == 2)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "parser_read_lr_wake",
      &python_parser_read_lr_wake,
      py::arg("ele"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("err_flag"),
      R"""(Subroutine parser_read_lr_wake (ele, delim, delim_found, err_flag)

  Subroutine to read in a long-range wake field from an external file.
  This subroutine is used by bmad_parser and bmad_parser2.

  Parameters
  ----------
  ele : EleStruct
      Element containing wake structure.
      This parameter is an input/output and is modified in-place. As an output: Element with wake information.
  )""");
  m.def(
      "parser_read_old_format_lr_wake",
      &Bmad::parser_read_old_format_lr_wake,
      py::arg("ele"),
      py::arg("lr_file_name"),
      R"""(Subroutine parser_read_old_format_lr_wake (ele, lr_file_name)

  Subroutine to read in a long-range wake field from an external file.
  This subroutine is used by bmad_parser and bmad_parser2.

  Parameters
  ----------
  ele : EleStruct
      Element containing wake structure.
      This parameter is an input/output and is modified in-place. As an output: Element with wake information.
  lr_file_name : unknown
      Name of long-range wake field file.
  )""");
  m.def(
      "parser_read_old_format_sr_wake",
      &Bmad::parser_read_old_format_sr_wake,
      py::arg("ele"),
      py::arg("sr_file_name"),
      R"""(Subroutine parser_read_old_format_sr_wake (ele, sr_file_name)

  Subroutine to read in a short-range wake field from an external file.
  This subroutine is used by bmad_parser and bmad_parser2.

  Parameters
  ----------
  ele : EleStruct
      Element containing wake structure.
      This parameter is an input/output and is modified in-place. As an output: Element with wake information.
  sr_file_name : unknown
      Name of short-range wake field file.
  )""");
  py::class_<PyParserReadSrWake, std::unique_ptr<PyParserReadSrWake>>(
      m, "ParserReadSrWake", "parser_read_sr_wake return type")
      .def_readonly("delim", &PyParserReadSrWake::delim)
      .def_readonly("delim_found", &PyParserReadSrWake::delim_found)
      .def_readonly("err_flag", &PyParserReadSrWake::err_flag)
      .def("__len__", [](const PyParserReadSrWake&) { return 3; })
      .def("__getitem__", [](const PyParserReadSrWake& s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.delim);
        if (i == 1)
          return py::cast(s.delim_found);
        if (i == 2)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "parser_read_sr_wake",
      &python_parser_read_sr_wake,
      py::arg("ele"),
      py::arg("delim"),
      py::arg("delim_found"),
      py::arg("err_flag"),
      R"""(Subroutine parser_read_sr_wake (ele, delim, delim_found, err_flag)

  Subroutine to read in a short-range wake field.
  This subroutine is used by bmad_parser and bmad_parser2.

  Parameters
  ----------
  ele : EleStruct
      Element containing wake structure.
      This parameter is an input/output and is modified in-place. As an output: Element with wake information.
  )""");
  m.def(
      "parser_transfer_control_struct",
      &Bmad::parser_transfer_control_struct,
      py::arg("con_in"),
      py::arg("lord"),
      py::arg("ix_var"),
      R"""(Subroutine parser_transfer_control_struct (con_in, con_out, lord, ix_var)

  Routine to transfer the information from an input control_struct (which stores
  the user input parameters) to a control_struct that will be stored in the lat%control
  or lord%control%ramp for a ramper.

  Parameters
  ----------
  con_in : ControlStruct
      Input control structure.
  lord : EleStruct
      Lord element associated with the control_struct.
  ix_var : int
      If an expression stack evaluates to a constant, this routine will modify the expression stack to evaluate
      to the value of: lord.control.var(ix_var) * constant

  Returns
  -------
  con_out : ControlStruct
      Output control structure.
  )""");
  m.def(
      "particle_in_global_frame",
      &Bmad::particle_in_global_frame,
      py::arg("orb"),
      py::arg("branch"),
      py::arg("in_time_coordinates") = py::none(),
      py::arg("in_body_frame") = py::none(),
      py::arg("w_mat_out") = py::none(),
      py::arg("particle"),
      R"""(Function particle_in_global_frame (orb, in_time_coordinates, in_body_frame, w_mat_out) result (particle)

  Returns the particle in global time coordinates given is coordinates orb in lattice lat.

  Parameters
  ----------
  orb : CoordStruct
      particle in s-coordinates
  branch : BranchStruct
      branch that contains branch.ele(orb.ix_ele)
  in_time_coordinates : bool
      Default is false. If true, orb will taken as in time coordinates.
  in_body_frame : bool
      Default is true. If false, ele offsets will be ignored. Result:
  particle : CoordStruct
      particle in global time coordinates

  Returns
  -------
  particle
  )""");
  py::class_<
      PyParticleIsMovingBackwards,
      std::unique_ptr<PyParticleIsMovingBackwards>>(
      m,
      "ParticleIsMovingBackwards",
      "particle_is_moving_backwards return type")
      .def_readonly(
          "is_moving_backwards",
          &PyParticleIsMovingBackwards::is_moving_backwards)
      .def("__len__", [](const PyParticleIsMovingBackwards&) { return 1; })
      .def(
          "__getitem__",
          [](const PyParticleIsMovingBackwards& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_moving_backwards);
            throw py::index_error();
          });
  m.def(
      "particle_is_moving_backwards",
      &python_particle_is_moving_backwards,
      py::arg("orbit"),
      py::arg("is_moving_backwards"),
      R"""(Parameters
  ----------
  orbit : CoordStruct
      Particle coordinates
  is_moving_backwards : 
  )""");
  py::class_<
      PyParticleIsMovingForward,
      std::unique_ptr<PyParticleIsMovingForward>>(
      m, "ParticleIsMovingForward", "particle_is_moving_forward return type")
      .def_readonly(
          "is_moving_forward", &PyParticleIsMovingForward::is_moving_forward)
      .def("__len__", [](const PyParticleIsMovingForward&) { return 1; })
      .def(
          "__getitem__",
          [](const PyParticleIsMovingForward& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_moving_forward);
            throw py::index_error();
          });
  m.def(
      "particle_is_moving_forward",
      &python_particle_is_moving_forward,
      py::arg("orbit"),
      py::arg("dir") = py::none(),
      py::arg("is_moving_forward"),
      R"""(Parameters
  ----------
  orbit : CoordStruct
      Particle coordinates
  dir : int, optional
      +1 if tracking forward(default) or -1 to return True if tracking backwards.
  is_moving_forward : 
  )""");
  py::class_<PyParticleRfTime, std::unique_ptr<PyParticleRfTime>>(
      m, "ParticleRfTime", "particle_rf_time return type")
      .def_readonly("time", &PyParticleRfTime::time)
      .def("__len__", [](const PyParticleRfTime&) { return 1; })
      .def("__getitem__", [](const PyParticleRfTime& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.time);
        throw py::index_error();
      });
  m.def(
      "particle_rf_time",
      &python_particle_rf_time,
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("reference_active_edge") = py::none(),
      py::arg("s_rel") = py::none(),
      py::arg("time_coords") = py::none(),
      py::arg("rf_freq") = py::none(),
      py::arg("abs_time") = py::none(),
      py::arg("time"),
      R"""(Parameters
  ----------
  orbit : CoordStruct
      Particle coordinates
  ele : EleStruct
      Element being tracked through.
  reference_active_edge : bool
      If True, and ele is a rfcavity or lcavity, use the active edge (edge of the -- logical: If True, and ele
      is a rfcavity or lcavity, use the active edge (edge of the region with non-zero field) as the reference
      point.
  s_rel : float, optional
      Longitudinal position relative to the upstream edge of the element. Needed for relative time tracking when
      the particle is inside the element. Default is 0.
  time_coords : bool, optional
      Default False. If True then orbit is using time based phase space coordinates.
  rf_freq : float, optional
      If present, the returned time shifted by an integer multiple of 1/rf_freq to be in the range
      [-1/2*rf_freq, 1/2*rf_freq]. This is useful to avoid round-off errors.
  abs_time : float, optional
      If False (default) use setting of bmad_com.absolute_time_tracking. If True, use absolute time instead of
      relative time. Ouput:
  time : 
  )""");
  py::class_<
      PyPatchFlipsPropagationDirection,
      std::unique_ptr<PyPatchFlipsPropagationDirection>>(
      m,
      "PatchFlipsPropagationDirection",
      "patch_flips_propagation_direction return type")
      .def_readonly("is_flip", &PyPatchFlipsPropagationDirection::is_flip)
      .def("__len__", [](const PyPatchFlipsPropagationDirection&) { return 1; })
      .def(
          "__getitem__",
          [](const PyPatchFlipsPropagationDirection& s, int i) -> py::object {
            if (i < 0)
              i += 1;
            if (i == 0)
              return py::cast(s.is_flip);
            throw py::index_error();
          });
  m.def(
      "patch_flips_propagation_direction",
      &python_patch_flips_propagation_direction,
      py::arg("x_pitch"),
      py::arg("y_pitch"),
      py::arg("is_flip"),
      R"""(Parameters
  ----------
  x_pitch : float
      Rotaion around y-axis
  y_pitch : float
      Rotation around x-axis.
  is_flip : 
  )""");
  py::class_<PyPatchLength, std::unique_ptr<PyPatchLength>>(
      m, "PatchLength", "patch_length return type")
      .def_readonly("length", &PyPatchLength::length)
      .def("__len__", [](const PyPatchLength&) { return 1; })
      .def("__getitem__", [](const PyPatchLength& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.length);
        throw py::index_error();
      });
  m.def(
      "patch_length",
      &python_patch_length,
      py::arg("patch"),
      py::arg("ref_coords") = py::none(),
      py::arg("length"),
      R"""(Parameters
  ----------
  patch : EleStruct
      Patch element.
  ref_coords : int, optional
      Reference coords to use. entrance_end$, exit_end$ Default is nint(patch.value(ref_coords$)).
  length : 
  )""");
  py::class_<
      Bmad::PhotonAbsorptionAndPhaseShift,
      std::unique_ptr<Bmad::PhotonAbsorptionAndPhaseShift>>(
      m,
      "PhotonAbsorptionAndPhaseShift",
      "photon_absorption_and_phase_shift return type")
      .def_readonly(
          "absorption", &Bmad::PhotonAbsorptionAndPhaseShift::absorption)
      .def_readonly(
          "phase_shift", &Bmad::PhotonAbsorptionAndPhaseShift::phase_shift)
      .def_readonly("err_flag", &Bmad::PhotonAbsorptionAndPhaseShift::err_flag)
      .def(
          "__len__",
          [](const Bmad::PhotonAbsorptionAndPhaseShift&) { return 3; })
      .def(
          "__getitem__",
          [](const Bmad::PhotonAbsorptionAndPhaseShift& s,
             int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.absorption);
            if (i == 1)
              return py::cast(s.phase_shift);
            if (i == 2)
              return py::cast(s.err_flag);
            throw py::index_error();
          });
  m.def(
      "photon_absorption_and_phase_shift",
      &Bmad::photon_absorption_and_phase_shift,
      py::arg("material"),
      py::arg("Energy"),
      R"""(Subroutine photon_absorption_and_phase_shift (material, Energy, absorption, phase_shift, err_flag)

  Routine to calcualte the absorption and phase shift values for a photon with a given
  energy going through a particular material.

  Parameters
  ----------
  material : unknown
      Material name.
  Energy : float
      Photon energy (eV).

  Returns
  -------
  absorption : float
      E_field ~ Exp(-absorption * length)
  phase_shift : float
      E_field Phase shift (radians) per unit length relative to vacuum.
  err_flag : bool
      Set true if material not recognized.
  )""");
  py::class_<
      PyPhotonAddToDetectorStatistics,
      std::unique_ptr<PyPhotonAddToDetectorStatistics>>(
      m,
      "PhotonAddToDetectorStatistics",
      "photon_add_to_detector_statistics return type")
      .def_readonly("ix_pt", &PyPhotonAddToDetectorStatistics::ix_pt)
      .def_readonly("iy_pt", &PyPhotonAddToDetectorStatistics::iy_pt)
      .def("__len__", [](const PyPhotonAddToDetectorStatistics&) { return 2; })
      .def(
          "__getitem__",
          [](const PyPhotonAddToDetectorStatistics& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.ix_pt);
            if (i == 1)
              return py::cast(s.iy_pt);
            throw py::index_error();
          });
  m.def(
      "photon_add_to_detector_statistics",
      &python_photon_add_to_detector_statistics,
      py::arg("orbit0"),
      py::arg("orbit"),
      py::arg("ele"),
      py::arg("ix_pt") = py::none(),
      py::arg("iy_pt") = py::none(),
      py::arg("pixel_pt") = py::none(),
      R"""(Subroutine photon_add_to_detector_statistics (orbit0, orbit, ele, ix_pt, iy_pt, pixel_pt)

  Routine to add photon statistics to the appropriate pixel of a "detector" grid.

  It is assumed that track_to_surface has been called so that the photon is at the
  detector surface and that orbit%vec(1) and %vec(3) are in surface coords (needed for curved detectors).

  Parameters
  ----------
  orbit0 : CoordStruct
      Photon coords at beginning of lattice
  orbit : CoordStruct
      Photon coords at the detector.
  ele : EleStruct
      Element with grid.
      This parameter is an input/output and is modified in-place. As an output: Element with updatted grid.
  pixel_pt : PixelPtStruct, optional
      If present then use this grid point instead of the grid point determined by the (x, y) coords of the
      photon
  )""");
  py::class_<Bmad::PhotonReflection, std::unique_ptr<Bmad::PhotonReflection>>(
      m, "PhotonReflection", "photon_reflection return type")
      .def_readonly("graze_angle_out", &Bmad::PhotonReflection::graze_angle_out)
      .def_readonly("phi_out", &Bmad::PhotonReflection::phi_out)
      .def("__len__", [](const Bmad::PhotonReflection&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::PhotonReflection& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.graze_angle_out);
            if (i == 1)
              return py::cast(s.phi_out);
            throw py::index_error();
          });
  m.def(
      "photon_reflection",
      &Bmad::photon_reflection,
      py::arg("graze_angle_in"),
      py::arg("energy"),
      py::arg("surface"),
      R"""(Subroutine photon_reflection (graze_angle_in, energy, surface, graze_angle_out, phi_out)

  Routine to reflect a photon from a surface including both diffuse and specular reflections.

  Parameters
  ----------
  graze_angle_in : float
      Incident grazing (not polar) angle in radians.
  energy : float
      Photon energy in eV.
  surface : PhotonReflectSurfaceStruct
      surface info

  Returns
  -------
  graze_angle_out : float
      graze_angle in radians.
  phi_out : float
      Azimuthal angle in radians.
  )""");
  m.def(
      "photon_reflection_std_surface_init",
      &Bmad::photon_reflection_std_surface_init,
      R"""(Subroutine photon_reflection_std_surface_init (surface)

  Routine to initialize the standard proton reflection probability tables.
  The standard tables are for 10 nm C film on Al substrate.
  The surface roughness for diffuse scattering is 200 nm and the
  the surface roughness correlation length is 5.5 um.


  Returns
  -------
  surface : 
      photon_reflect_surface_struct
  )""");
  py::class_<
      Bmad::PhotonReflectivity,
      std::unique_ptr<Bmad::PhotonReflectivity>>(
      m, "PhotonReflectivity", "photon_reflectivity return type")
      .def_readonly("p_reflect", &Bmad::PhotonReflectivity::p_reflect)
      .def_readonly("rel_p_specular", &Bmad::PhotonReflectivity::rel_p_specular)
      .def("__len__", [](const Bmad::PhotonReflectivity&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::PhotonReflectivity& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.p_reflect);
            if (i == 1)
              return py::cast(s.rel_p_specular);
            throw py::index_error();
          });
  m.def(
      "photon_reflectivity",
      &Bmad::photon_reflectivity,
      py::arg("angle"),
      py::arg("energy"),
      py::arg("surface"),
      R"""(Subroutine photon_reflectivity (angle, energy, surface, p_reflect, rel_p_specular)

  Routine to evaluate the photon reflectivity.
    probability of absorption          = 1 - p_reflect
    probability of reflection          = p_reflect
    probability of specular reflection = p_reflect * rel_p_specular
    probability of diffuse reflection  = p_reflect * (1 - rel_p_specular)

  Use photon_reflection_std_surface_init or read_surface_reflection_file to get surface info.

  Parameters
  ----------
  angle : float
      Incident grazing angle in radians.
  energy : float
      Photon energy in eV.
  surface : PhotonReflectSurfaceStruct
      surface info

  Returns
  -------
  p_reflect : float
      Reflection probability.
  rel_p_specular : float
      Relative specular reflection probability.
  )""");
  py::class_<
      PyPhotonTargetCornerCalc,
      std::unique_ptr<PyPhotonTargetCornerCalc>>(
      m, "PhotonTargetCornerCalc", "photon_target_corner_calc return type")
      .def_readonly("corner", &PyPhotonTargetCornerCalc::corner)
      .def_readonly("x_lim", &PyPhotonTargetCornerCalc::x_lim)
      .def_readonly("y_lim", &PyPhotonTargetCornerCalc::y_lim)
      .def_readonly("z_lim", &PyPhotonTargetCornerCalc::z_lim)
      .def("__len__", [](const PyPhotonTargetCornerCalc&) { return 4; })
      .def(
          "__getitem__",
          [](const PyPhotonTargetCornerCalc& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.corner);
            if (i == 1)
              return py::cast(s.x_lim);
            if (i == 2)
              return py::cast(s.y_lim);
            if (i == 3)
              return py::cast(s.z_lim);
            throw py::index_error();
          });
  m.def(
      "photon_target_corner_calc",
      &python_photon_target_corner_calc,
      py::arg("aperture_ele"),
      py::arg("x_lim"),
      py::arg("y_lim"),
      py::arg("z_lim"),
      py::arg("source_ele"),
      R"""(Subroutine photon_target_corner_calc (aperture_ele, x_lim, y_lim, z_lim, source_ele, corner)

  Routine to calculate the corner coords in the source_ele ref frame.

  Parameters
  ----------
  aperture_ele : EleStruct
      Element containing the aperture x_lim, y_lim  -- real(rp): Transverse corner points in aperture_ele coord
      frame.
  source_ele : EleStruct
      Photon source element.

  Returns
  -------
  corner : TargetPointStruct
      Corner coords in source_ele ref frame.
  )""");
  m.def(
      "photon_target_setup",
      &Bmad::photon_target_setup,
      py::arg("ele"),
      R"""(Subroutine photon_target_setup (ele)

  Routine to calculate and store the parmeters needed for photon targeting.
  This routine is called by Bmad parsing routines and is not meant for general use.

  Photon initialization with targeting is done by the routine init_photon_from_a_photon_init_ele
  Which is called by init_coord.

  Parameters
  ----------
  ele : EleStruct
      Source element to setup. Element will be of type: sample, diffraction_plate or photon_init.
      This parameter is an input/output and is modified in-place. As an output: Source element with target
      parameters setup in ele.photon.target.
  )""");
  m.def(
      "photon_type",
      &Bmad::photon_type,
      py::arg("ele"),
      R"""(Function photon_type (ele) result (e_type)

  Routine to return the type of photon to be tracked: coherent$ or incoherent$.

  Parameters
  ----------
  ele : EleStruct
      Element being tracked through.

  Returns
  -------
  e_type : int
      coherent$ or incoherent$
  )""");
  py::class_<PyPhysicalEleEnd, std::unique_ptr<PyPhysicalEleEnd>>(
      m, "PhysicalEleEnd", "physical_ele_end return type")
      .def_readonly("physical_end", &PyPhysicalEleEnd::physical_end)
      .def("__len__", [](const PyPhysicalEleEnd&) { return 1; })
      .def("__getitem__", [](const PyPhysicalEleEnd& s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.physical_end);
        throw py::index_error();
      });
  m.def(
      "physical_ele_end",
      &python_physical_ele_end,
      py::arg("track_end"),
      py::arg("orbit"),
      py::arg("ele_orientation"),
      py::arg("return_stream_end") = py::none(),
      py::arg("physical_end"),
      R"""(Parameters
  ----------
  track_end : int
      first_track_edge$, second_track_edge$, surface$, or in_between$
  orbit : CoordStruct
      Particle position.
  ele_orientation : int
      Either 1 = Normal or -1 = element reversed.
  return_stream_end : bool, optional
      If True return the stream end instead of the physical end. Default is False.
  physical_end : 
  )""");
  m.def(
      "point_photon_emission",
      &Bmad::point_photon_emission,
      py::arg("ele"),
      py::arg("param"),
      py::arg("orbit"),
      py::arg("direction"),
      py::arg("max_target_area"),
      py::arg("w_to_surface") = py::none(),
      R"""(Subroutine point_photon_emission (ele, param, orbit, direction, max_target_area, w_to_surface)

  Routine to emit a photon from a point that may be on a surface.
  If there is a downstream target, the emission calc will take this into account.

  Parameters
  ----------
  ele : EleStruct
      Emitting element.
  param : LatParamStruct
      lattice parameters.
  orbit : CoordStruct
      phase-space coords of photon. --   Will be in curved surface coords if there is a curved surface.
      This parameter is an input/output and is modified in-place. As an output: Final phase-space coords
  direction : int
      +1 -> Emit in forward +z direction, -1 -> emit backwards.
  max_target_area : float
      Area of the solid angle photons may be emitted over. max_target_area is used for normalizing the photon
      field. generally will be equal to twopi or fourpi.
  w_to_surface : float, optional
      Rotation matrix for curved surface.
  )""");
  m.def(
      "pointer_to_branch",
      py::overload_cast<EleProxy&>(&Bmad::pointer_to_branch),
      py::arg("ele"),
      R"""(Function pointer_to_branch

  Routine to return a pointer to the lattice branch associated with a given name
  or a given element.

  This routine is an overloaded name for:
    pointer_to_branch_given_ele (ele) result (branch_ptr))
    pointer_to_branch_given_name (branch_name, lat, parameter_is_branch0, blank_branch) result (branch_ptr)

  The lattice branch *associated* with a given element is not necessarily the
  branch where the element is *located*. For example, all lords live in branch #0.
  But the branch associated with a super_lord element is the branch of its slaves.

  To get the branch where the element is located, simply use ele%ix_branch.

  Note: Result is ambiguous if ele argument is associated with multiple branches
  which can happen, for example, with overlay elements.

  Parameters
  ----------
  ele : EleStruct
      Element contained in the branch.
  branch_name : unknown
      May be a branch name or a branch index.
  lat : LatStruct
      Lattice to search.
  parameter_is_branch0 : bool, optional
      If True, 'PARAMETER' is taken to be an alternative name for branch(0). Default is False.
  blank_branch : int, optional
      Branch index if branch_name = ''. Default is blank is an error.

  Returns
  -------
  branch_ptr : BranchStruct
      Pointer to the branch. Nullified if there is no associated branch.

  Notes
  -----
  Overloaded versions:
  )""");
  m.def(
      "pointer_to_branch",
      py::overload_cast<
          std::string,
          LatProxy&,
          std::optional<bool>,
          std::optional<int>>(&Bmad::pointer_to_branch),
      py::arg("branch_name"),
      py::arg("lat"),
      py::arg("parameter_is_branch0") = py::none(),
      py::arg("blank_branch") = py::none(),
      R"""(Function pointer_to_branch

  Routine to return a pointer to the lattice branch associated with a given name
  or a given element.

  This routine is an overloaded name for:
    pointer_to_branch_given_ele (ele) result (branch_ptr))
    pointer_to_branch_given_name (branch_name, lat, parameter_is_branch0, blank_branch) result (branch_ptr)

  The lattice branch *associated* with a given element is not necessarily the
  branch where the element is *located*. For example, all lords live in branch #0.
  But the branch associated with a super_lord element is the branch of its slaves.

  To get the branch where the element is located, simply use ele%ix_branch.

  Note: Result is ambiguous if ele argument is associated with multiple branches
  which can happen, for example, with overlay elements.

  Parameters
  ----------
  ele : EleStruct
      Element contained in the branch.
  branch_name : unknown
      May be a branch name or a branch index.
  lat : LatStruct
      Lattice to search.
  parameter_is_branch0 : bool, optional
      If True, 'PARAMETER' is taken to be an alternative name for branch(0). Default is False.
  blank_branch : int, optional
      Branch index if branch_name = ''. Default is blank is an error.

  Returns
  -------
  branch_ptr : BranchStruct
      Pointer to the branch. Nullified if there is no associated branch.

  Notes
  -----
  Overloaded versions:
  )""");
  m.def(
      "pointer_to_ele",
      py::overload_cast<LatProxy&, int, std::optional<int>>(
          &Bmad::pointer_to_ele),
      py::arg("lat"),
      py::arg("ix_ele"),
      py::arg("ix_branch") = py::none(),
      R"""(Function pointer_to_ele (...)

  Routine to return a pointer to an element.
  pointer_to_ele is an overloaded name for:
      Function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr)
      Function pointer_to_ele2 (lat, ele_loc) result (ele_ptr)
      Function pointer_to_ele3 (lat, ele_name) result (ele_ptr)
      Function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)

  pointer_to_ele4(lat, foreign_ele) is useful when foreign_ele is associated with a separate
  lattice that has an identical layout. pointer_to_ele4 will then return the corresponding
  element in lat.

  Note that using ele_name to locate an element is potentially dangerous if there
  are multiple elements that have the same name. Better in this case is to use:
    lat_ele_locator

  Parameters
  ----------
  lat : LatStruct
      Lattice.
  ix_ele : int
      Index of element in lat.branch(ix_branch).
  ix_branch : int
      Index of the lat.branch(:) containing the element.
  ix_nametable : int
      Nametable index. See above
  ele_loc : LatEleLocStruct
      Location identification.
  ele_name : unknown
      Name or index of element.
  foreign_ele : EleStruct
      Lattice element in another lattice.

  Returns
  -------
  ele_ptr : EleStruct
      Pointer to the element. Nullified if no match or error.

  Notes
  -----
  Related routines:
  pointer_to_slave pointer_to_lord
  Overloaded versions: Function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr), Function
  pointer_to_ele2 (lat, ele_loc) result (ele_ptr), Function pointer_to_ele3 (lat, ele_name) result (ele_ptr),
  Function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)
  )""");
  m.def(
      "pointer_to_ele",
      py::overload_cast<LatProxy&, LatEleLocProxy&>(&Bmad::pointer_to_ele),
      py::arg("lat"),
      py::arg("ele_loc"),
      R"""(Function pointer_to_ele (...)

  Routine to return a pointer to an element.
  pointer_to_ele is an overloaded name for:
      Function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr)
      Function pointer_to_ele2 (lat, ele_loc) result (ele_ptr)
      Function pointer_to_ele3 (lat, ele_name) result (ele_ptr)
      Function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)

  pointer_to_ele4(lat, foreign_ele) is useful when foreign_ele is associated with a separate
  lattice that has an identical layout. pointer_to_ele4 will then return the corresponding
  element in lat.

  Note that using ele_name to locate an element is potentially dangerous if there
  are multiple elements that have the same name. Better in this case is to use:
    lat_ele_locator

  Parameters
  ----------
  lat : LatStruct
      Lattice.
  ix_ele : int
      Index of element in lat.branch(ix_branch).
  ix_branch : int
      Index of the lat.branch(:) containing the element.
  ix_nametable : int
      Nametable index. See above
  ele_loc : LatEleLocStruct
      Location identification.
  ele_name : unknown
      Name or index of element.
  foreign_ele : EleStruct
      Lattice element in another lattice.

  Returns
  -------
  ele_ptr : EleStruct
      Pointer to the element. Nullified if no match or error.

  Notes
  -----
  Related routines:
  pointer_to_slave pointer_to_lord
  Overloaded versions: Function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr), Function
  pointer_to_ele2 (lat, ele_loc) result (ele_ptr), Function pointer_to_ele3 (lat, ele_name) result (ele_ptr),
  Function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)
  )""");
  m.def(
      "pointer_to_ele",
      py::overload_cast<LatProxy&, std::string>(&Bmad::pointer_to_ele),
      py::arg("lat"),
      py::arg("ele_name"),
      R"""(Function pointer_to_ele (...)

  Routine to return a pointer to an element.
  pointer_to_ele is an overloaded name for:
      Function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr)
      Function pointer_to_ele2 (lat, ele_loc) result (ele_ptr)
      Function pointer_to_ele3 (lat, ele_name) result (ele_ptr)
      Function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)

  pointer_to_ele4(lat, foreign_ele) is useful when foreign_ele is associated with a separate
  lattice that has an identical layout. pointer_to_ele4 will then return the corresponding
  element in lat.

  Note that using ele_name to locate an element is potentially dangerous if there
  are multiple elements that have the same name. Better in this case is to use:
    lat_ele_locator

  Parameters
  ----------
  lat : LatStruct
      Lattice.
  ix_ele : int
      Index of element in lat.branch(ix_branch).
  ix_branch : int
      Index of the lat.branch(:) containing the element.
  ix_nametable : int
      Nametable index. See above
  ele_loc : LatEleLocStruct
      Location identification.
  ele_name : unknown
      Name or index of element.
  foreign_ele : EleStruct
      Lattice element in another lattice.

  Returns
  -------
  ele_ptr : EleStruct
      Pointer to the element. Nullified if no match or error.

  Notes
  -----
  Related routines:
  pointer_to_slave pointer_to_lord
  Overloaded versions: Function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr), Function
  pointer_to_ele2 (lat, ele_loc) result (ele_ptr), Function pointer_to_ele3 (lat, ele_name) result (ele_ptr),
  Function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)
  )""");
  m.def(
      "pointer_to_ele",
      py::overload_cast<LatProxy&, EleProxy&>(&Bmad::pointer_to_ele),
      py::arg("lat"),
      py::arg("foreign_ele"),
      R"""(Function pointer_to_ele (...)

  Routine to return a pointer to an element.
  pointer_to_ele is an overloaded name for:
      Function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr)
      Function pointer_to_ele2 (lat, ele_loc) result (ele_ptr)
      Function pointer_to_ele3 (lat, ele_name) result (ele_ptr)
      Function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)

  pointer_to_ele4(lat, foreign_ele) is useful when foreign_ele is associated with a separate
  lattice that has an identical layout. pointer_to_ele4 will then return the corresponding
  element in lat.

  Note that using ele_name to locate an element is potentially dangerous if there
  are multiple elements that have the same name. Better in this case is to use:
    lat_ele_locator

  Parameters
  ----------
  lat : LatStruct
      Lattice.
  ix_ele : int
      Index of element in lat.branch(ix_branch).
  ix_branch : int
      Index of the lat.branch(:) containing the element.
  ix_nametable : int
      Nametable index. See above
  ele_loc : LatEleLocStruct
      Location identification.
  ele_name : unknown
      Name or index of element.
  foreign_ele : EleStruct
      Lattice element in another lattice.

  Returns
  -------
  ele_ptr : EleStruct
      Pointer to the element. Nullified if no match or error.

  Notes
  -----
  Related routines:
  pointer_to_slave pointer_to_lord
  Overloaded versions: Function pointer_to_ele1 (lat, ix_ele, ix_branch) result (ele_ptr), Function
  pointer_to_ele2 (lat, ele_loc) result (ele_ptr), Function pointer_to_ele3 (lat, ele_name) result (ele_ptr),
  Function pointer_to_ele4 (lat, foreign_ele) result (ele_ptr)
  )""");
  py::class_<
      Bmad::PointerToElementAtS,
      std::unique_ptr<Bmad::PointerToElementAtS>>(
      m, "PointerToElementAtS", "pointer_to_element_at_s return type")
      .def_readonly("err_flag", &Bmad::PointerToElementAtS::err_flag)
      .def_readonly("s_eff", &Bmad::PointerToElementAtS::s_eff)
      .def_readonly("position", &Bmad::PointerToElementAtS::position)
      .def_readonly("ele", &Bmad::PointerToElementAtS::ele)
      .def("__len__", [](const Bmad::PointerToElementAtS&) { return 4; })
      .def(
          "__getitem__",
          [](const Bmad::PointerToElementAtS& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.err_flag);
            if (i == 1)
              return py::cast(s.s_eff);
            if (i == 2)
              return py::cast(s.position);
            if (i == 3)
              return py::cast(s.ele);
            throw py::index_error();
          });
  m.def(
      "pointer_to_element_at_s",
      &Bmad::pointer_to_element_at_s,
      py::arg("branch"),
      py::arg("s"),
      py::arg("choose_max"),
      py::arg("print_err") = py::none(),
      R"""(Function pointer_to_element_at_s (branch, s, choose_max, err_flag, s_eff, position) result (ele)

  Function to return a pointer to the element at position s.
  That is, return ele => branch%ele(ix_ele) such that:
  If choose_max = True:
      If s = branch%ele(ix_end_of_branch): ix_ele = ix_end_of_branch
      Else: branch%ele(ix_ele)%s_strat <= s < branch%ele(ix_ele)%s
  If choose_max = False:
      If s = branch%ele(0): ix_ele = 0
      Else: branch%ele(ix_ele)%s_start < s <= branch%ele(ix_ele)%s
  That is, if s corresponds to an element boundary between elements with indexes ix1 and ix2 = ix1 + 1:
      choose_max = True  => ix_ele = ix2
      choose_max = False => ix_ele = ix1

  Also see: element_at_s

  The setting of choose_max only makes a difference when s corresponds to an element boundary.

  Note: For a circular lattice, s is evaluated at the effective s which
  is modulo the branch length:
      s_eff = s - branch_length * floor(s/branch_length)

  Note: If there are multiple elements that are at the given s position due to the presence of
  an element with a negative length, which of the possible elements is actually chosen is ill-defined.

  Parameters
  ----------
  branch : BranchStruct
      Branch to use
  s : float
      Longitudinal position.
  choose_max : bool
      See above.
  print_err : bool, optional
      Print error message if there is an error? Default is True.

  Returns
  -------
  ele : EleStruct
      Pointer to element at s.
  err_flag : bool
      Set True if s is out of bounds. False otherwise.
  s_eff : float
      Effective s. Equal to s with a open lattice. See above.
  position : CoordStruct
      Positional information. .s         -- Same as input s. .ix_ele    -- Same as output ix_ele .location  --
      Location relative to element. Upstream_end$, downstream_end$, or inside$

  Notes
  -----
  Related routines:
  element_at_s The setting of choose_max only makes a difference when s corresponds to an element boundary. For
  a circular lattice s is evaluated at the effective s which s_eff = s - branch_length * floor(s/branch_length)
  If there are multiple elements that are at the given s position due to the presence of an element with a
  negative length which of the possible elements is actually chosen is ill-defined.
  )""");
  m.def(
      "pointer_to_field_ele",
      &Bmad::pointer_to_field_ele,
      py::arg("ele"),
      py::arg("ix_field_ele"),
      py::arg("field_ele"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Element with sum number of associated field elements.
  ix_field_ele : int
      Index of the field element to point to. This index runs from 1 to num_field_eles(ele).
  dz_offset : float
      Longitudinal offset of ele upstream edge from the field ele pointed to.
  field_ele : 
  )""");
  m.def(
      "pointer_to_girder",
      &Bmad::pointer_to_girder,
      py::arg("ele"),
      py::arg("girder"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Element to check.
  ix_slave_back : int
      Index back to ele. That is, pointer_to_slave(girder, ix_slave_back) will point back to ele. Set to -1 if
      no girder present
  girder : 
  )""");
  py::class_<Bmad::PointerToLord, std::unique_ptr<Bmad::PointerToLord>>(
      m, "PointerToLord", "pointer_to_lord return type")
      .def_readonly("control", &Bmad::PointerToLord::control)
      .def_readonly("ix_slave_back", &Bmad::PointerToLord::ix_slave_back)
      .def_readonly("ix_control", &Bmad::PointerToLord::ix_control)
      .def_readonly("ix_ic", &Bmad::PointerToLord::ix_ic)
      .def("__len__", [](const Bmad::PointerToLord&) { return 4; })
      .def(
          "__getitem__", [](const Bmad::PointerToLord& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.control);
            if (i == 1)
              return py::cast(s.ix_slave_back);
            if (i == 2)
              return py::cast(s.ix_control);
            if (i == 3)
              return py::cast(s.ix_ic);
            throw py::index_error();
          });
  m.def(
      "pointer_to_lord",
      &Bmad::pointer_to_lord,
      py::arg("slave"),
      py::arg("ix_lord"),
      py::arg("lord_type") = py::none(),
      py::arg("lord_ptr"),
      R"""(Parameters
  ----------
  slave : EleStruct
      Slave element.
  ix_lord : int
      Index of the lord.
  control : ControlStruct
      Pointer to control info for this lord/slave relationship. Nullified if there is an error.
  ix_slave_back : int
      Index back to the slave. That is, pointer_to_slave(lord_ptr, ix_slave_back) will point back to slave. Set
      to -1 if there is an error or the slave is a slice_slave.
  lord_type : int, optional
      See above.
  ix_control : int
      Index in lat.control(:) array the control argument is at. For ramper lord elements, ix_control is index
      for the lord.control.ramper(:) array.
  ix_ic : int
      Index of the lat.ic(:) element associated with the control argument.
  lord_ptr : 
  )""");
  py::class_<
      Bmad::PointerToMultipassLord,
      std::unique_ptr<Bmad::PointerToMultipassLord>>(
      m, "PointerToMultipassLord", "pointer_to_multipass_lord return type")
      .def_readonly("ix_pass", &Bmad::PointerToMultipassLord::ix_pass)
      .def_readonly("super_lord", &Bmad::PointerToMultipassLord::super_lord)
      .def("__len__", [](const Bmad::PointerToMultipassLord&) { return 2; })
      .def(
          "__getitem__",
          [](const Bmad::PointerToMultipassLord& s, int i) -> py::object {
            if (i < 0)
              i += 2;
            if (i == 0)
              return py::cast(s.ix_pass);
            if (i == 1)
              return py::cast(s.super_lord);
            throw py::index_error();
          });
  m.def(
      "pointer_to_multipass_lord",
      &Bmad::pointer_to_multipass_lord,
      py::arg("ele"),
      py::arg("multi_lord"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Lattice element.
  ix_pass : int
      Multipass turn number. Set to 0 if element is a multipass_lord. Set to -1 if element is not a
      multipass_slave.
  super_lord : EleStruct
      super_lord of the element. Set to NULL if ele is not a super_slave or super_lord. Note: if ele is a
      multipass_lord there are multiple possible super_lord slaves.
  multi_lord : 
  )""");
  m.def(
      "pointer_to_next_ele",
      &Bmad::pointer_to_next_ele,
      py::arg("this_ele"),
      py::arg("offset") = py::none(),
      py::arg("skip_beginning") = py::none(),
      py::arg("follow_fork") = py::none(),
      py::arg("next_ele"),
      R"""(Parameters
  ----------
  this_ele : EleStruct
      Starting element.
  offset : int, optional
      +1 -> return next element, +2 -> element after that, etc. Can be negative. Default = +1.
  skip_beginning : bool, optional
      If True then skip beginning element #0 when wrapping around. Default is False.
  follow_fork : bool, optional
      If True then fork at any fork element. Default is False.
  next_ele : 
  )""");
  py::class_<Bmad::PointerToSlave, std::unique_ptr<Bmad::PointerToSlave>>(
      m, "PointerToSlave", "pointer_to_slave return type")
      .def_readonly("control", &Bmad::PointerToSlave::control)
      .def_readonly("ix_lord_back", &Bmad::PointerToSlave::ix_lord_back)
      .def_readonly("ix_control", &Bmad::PointerToSlave::ix_control)
      .def_readonly("ix_ic", &Bmad::PointerToSlave::ix_ic)
      .def_readonly("slave_ptr", &Bmad::PointerToSlave::slave_ptr)
      .def("__len__", [](const Bmad::PointerToSlave&) { return 5; })
      .def(
          "__getitem__",
          [](const Bmad::PointerToSlave& s, int i) -> py::object {
            if (i < 0)
              i += 5;
            if (i == 0)
              return py::cast(s.control);
            if (i == 1)
              return py::cast(s.ix_lord_back);
            if (i == 2)
              return py::cast(s.ix_control);
            if (i == 3)
              return py::cast(s.ix_ic);
            if (i == 4)
              return py::cast(s.slave_ptr);
            throw py::index_error();
          });
  m.def(
      "pointer_to_slave",
      &Bmad::pointer_to_slave,
      py::arg("lord"),
      py::arg("ix_slave"),
      py::arg("lord_type") = py::none(),
      R"""(Function pointer_to_slave (lord, ix_slave, control, lord_type, ix_lord_back, ix_control, ix_ic) result (slave_ptr)

  Function to point to a slave of a lord.
  Note: Ramper lords do not have any associated slaves (slaves are assigned dynamically at run time).

  If lord_type = all$ (the default) the range for ix_slave is:
    1 to lord%n_slave                                 for "regular" slaves.
    lord%n_slave+1 to lord%n_slave+lord%n_slave_field for field overlap slaves.

  If lord_type = field_lord$, only the field overlap slaves may be accessed and the range for ix_slave is:
    1 to lord%n_slave_field

  Parameters
  ----------
  lord : EleStruct
      Lord element
  ix_slave : int
      Index of the slave in the list of slaves controled by the lord..
  lord_type : int, optional
      See above.

  Returns
  -------
  slave_ptr : EleStruct
      Pointer to the slave. Nullified if there is an error.
  control : ControlStruct
      Pointer to control info for this lord/slave relationship. Nullified if there is an error.
  ix_lord_back : int
      Index back to the lord. That is, pointer_to_lord(slave_ptr, ix_lord_back) will point back to the lord. Set
      to -1 if there is an error.
  ix_control : int
      Index in lat.control(:) array the control argument is at.
  ix_ic : int
      Index of the lat.ic(:) element associated with the control argument.

  Notes
  -----
  Related routines:
  pointer_to_lord pointer_to_super_lord pointer_to_ele num_lords
  )""");
  py::class_<
      Bmad::PointerToSuperLord,
      std::unique_ptr<Bmad::PointerToSuperLord>>(
      m, "PointerToSuperLord", "pointer_to_super_lord return type")
      .def_readonly("control", &Bmad::PointerToSuperLord::control)
      .def_readonly("ix_slave_back", &Bmad::PointerToSuperLord::ix_slave_back)
      .def_readonly("ix_control", &Bmad::PointerToSuperLord::ix_control)
      .def_readonly("ix_ic", &Bmad::PointerToSuperLord::ix_ic)
      .def("__len__", [](const Bmad::PointerToSuperLord&) { return 4; })
      .def(
          "__getitem__",
          [](const Bmad::PointerToSuperLord& s, int i) -> py::object {
            if (i < 0)
              i += 4;
            if (i == 0)
              return py::cast(s.control);
            if (i == 1)
              return py::cast(s.ix_slave_back);
            if (i == 2)
              return py::cast(s.ix_control);
            if (i == 3)
              return py::cast(s.ix_ic);
            throw py::index_error();
          });
  m.def(
      "pointer_to_super_lord",
      &Bmad::pointer_to_super_lord,
      py::arg("slave"),
      py::arg("lord_type") = py::none(),
      py::arg("lord_ptr"),
      R"""(Parameters
  ----------
  slave : EleStruct
      Slave element.
  control : ControlStruct
      Pointer to control info for this lord/slave relationship. Nullified if there is an error.
  ix_slave_back : int
      Index back to the slave. That is, pointer_to_slave(lord_ptr, ix_slave_back) will point back to slave. Set
      to -1 if there is an error or the slave is a slice_slave.
  ix_control : int
      Index in lat.control(:) array the control argument is at. For ramper lord elements, ix_control is index
      for the lord.control.ramper(:) array.
  ix_ic : int
      Index of the lat.ic(:) element associated with the control argument.
  lord_type : int, optional
      If present, only return a super_lord of this type.
  lord_ptr : 
  )""");
  py::class_<
      PyPointerToSurfaceDisplacementPt,
      std::unique_ptr<PyPointerToSurfaceDisplacementPt>>(
      m,
      "PointerToSurfaceDisplacementPt",
      "pointer_to_surface_displacement_pt return type")
      .def_readonly("pt", &PyPointerToSurfaceDisplacementPt::pt)
      .def_readonly("x", &PyPointerToSurfaceDisplacementPt::x)
      .def_readonly("y", &PyPointerToSurfaceDisplacementPt::y)
      .def_readonly("ix", &PyPointerToSurfaceDisplacementPt::ix)
      .def_readonly("iy", &PyPointerToSurfaceDisplacementPt::iy)
      .def_readonly("xx", &PyPointerToSurfaceDisplacementPt::xx)
      .def_readonly("yy", &PyPointerToSurfaceDisplacementPt::yy)
      .def("__len__", [](const PyPointerToSurfaceDisplacementPt&) { return 7; })
      .def(
          "__getitem__",
          [](const PyPointerToSurfaceDisplacementPt& s, int i) -> py::object {
            if (i < 0)
              i += 7;
            if (i == 0)
              return py::cast(s.pt);
            if (i == 1)
              return py::cast(s.x);
            if (i == 2)
              return py::cast(s.y);
            if (i == 3)
              return py::cast(s.ix);
            if (i == 4)
              return py::cast(s.iy);
            if (i == 5)
              return py::cast(s.xx);
            if (i == 6)
              return py::cast(s.yy);
            throw py::index_error();
          });
  m.def(
      "pointer_to_surface_displacement_pt",
      &python_pointer_to_surface_displacement_pt,
      py::arg("ele"),
      py::arg("nearest"),
      py::arg("x"),
      py::arg("y"),
      py::arg("ix") = py::none(),
      py::arg("iy") = py::none(),
      py::arg("extend_grid") = py::none(),
      py::arg("xx") = py::none(),
      py::arg("yy") = py::none(),
      R"""(Function pointer_to_surface_displacement_pt (ele, nearest, x, y, ix, iy, extend_grid, xx, yy) result (pt)

  Routine to point to the grid point struct associated with point (x,y).

  Note: If nearest = True, the grid boundary is a length dr/2 from the boundary grid points.

  Parameters
  ----------
  ele : EleStruct
      Element containing the grid
  nearest : bool
      If True, return pointer to nearest grid point. If False, return pointer to the grid point lower and left
      of (x,y). x, y        -- real(rp): Photon position.
  extend_grid : bool, optional
      If (x,y) past grid pretend (x,y) is at grid boundary. Default is False. ix, iy      -- integer, optional:
      Grid point index.

  Returns
  -------
  pt : GridPointStruct
      Pointer to grid point. Will not be associated if (x,y) outside the grid. xx, yy      -- real(rp),
      optional: Set equal to (x, y) except if (x,y) is outside of the grid. In this case, (xx, yy) will be set
      to be on the nearest grid boundary point.
  )""");
  py::class_<
      PyPointerToSurfaceSegmentedPt,
      std::unique_ptr<PyPointerToSurfaceSegmentedPt>>(
      m,
      "PointerToSurfaceSegmentedPt",
      "pointer_to_surface_segmented_pt return type")
      .def_readonly("pt", &PyPointerToSurfaceSegmentedPt::pt)
      .def_readonly("x", &PyPointerToSurfaceSegmentedPt::x)
      .def_readonly("y", &PyPointerToSurfaceSegmentedPt::y)
      .def_readonly("ix", &PyPointerToSurfaceSegmentedPt::ix)
      .def_readonly("iy", &PyPointerToSurfaceSegmentedPt::iy)
      .def_readonly("xx", &PyPointerToSurfaceSegmentedPt::xx)
      .def_readonly("yy", &PyPointerToSurfaceSegmentedPt::yy)
      .def("__len__", [](const PyPointerToSurfaceSegmentedPt&) { return 7; })
      .def(
          "__getitem__",
          [](const PyPointerToSurfaceSegmentedPt& s, int i) -> py::object {
            if (i < 0)
              i += 7;
            if (i == 0)
              return py::cast(s.pt);
            if (i == 1)
              return py::cast(s.x);
            if (i == 2)
              return py::cast(s.y);
            if (i == 3)
              return py::cast(s.ix);
            if (i == 4)
              return py::cast(s.iy);
            if (i == 5)
              return py::cast(s.xx);
            if (i == 6)
              return py::cast(s.yy);
            throw py::index_error();
          });
  m.def(
      "pointer_to_surface_segmented_pt",
      &python_pointer_to_surface_segmented_pt,
      py::arg("ele"),
      py::arg("nearest"),
      py::arg("x"),
      py::arg("y"),
      py::arg("ix") = py::none(),
      py::arg("iy") = py::none(),
      py::arg("extend_grid") = py::none(),
      py::arg("xx") = py::none(),
      py::arg("yy") = py::none(),
      R"""(Function pointer_to_surface_segmented_pt (ele, nearest, x, y, ix, iy, extend_grid, xx, yy) result (pt)

  Routine to point to the grid point struct associated with point (x,y).

  Note: If nearest = True, the grid boundary is a length dr/2 from the boundary grid points.

  Parameters
  ----------
  ele : EleStruct
      Element containing the grid
  nearest : bool
      If True, return pointer to nearest grid point. If False, return pointer to the grid point lower and left
      of (x,y). x, y        -- real(rp): Photon position.
  extend_grid : bool, optional
      If (x,y) past grid pretend (x,y) is at grid boundary. Default is False. ix, iy      -- integer, optional:
      Grid point index.

  Returns
  -------
  pt : GridPointStruct
      Pointer to grid point. Will not be associated if (x,y) outside the grid. xx, yy      -- real(rp),
      optional: Set equal to (x, y) except if (x,y) is outside of the grid. In this case, (xx, yy) will be set
      to be on the nearest grid boundary point.
  )""");
  m.def(
      "pointer_to_wake_ele",
      &Bmad::pointer_to_wake_ele,
      py::arg("ele"),
      py::arg("wake_ele"),
      R"""(Parameters
  ----------
  ele : EleStruct
      Lattice element.
  delta_s : float
      distance of wake locaiton from beginning of ele.
  wake_ele : 
  )""");
  py::class_<Bmad::PointerToWall3d, std::unique_ptr<Bmad::PointerToWall3d>>(
      m, "PointerToWall3d", "pointer_to_wall3d return type")
      .def_readonly("ds_offset", &Bmad::PointerToWall3d::ds_offset)
      .def_readonly("is_branch_wall", &Bmad::PointerToWall3d::is_branch_wall)
      .def_readonly("wall3d", &Bmad::PointerToWall3d::wall3d)
      .def("__len__", [](const Bmad::PointerToWall3d&) { return 3; })
      .def(
          "__getitem__",
          [](const Bmad::PointerToWall3d& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.ds_offset);
            if (i == 1)
              return py::cast(s.is_branch_wall);
            if (i == 2)
              return py::cast(s.wall3d);
            throw py::index_error();
          });
  m.def(
      "pointer_to_wall3d",
      &Bmad::pointer_to_wall3d,
      py::arg("ele"),
      py::arg("ix_wall") = py::none(),
      R"""(Function pointer_to_wall3d (ele, ix_wall, ds_offset, is_branch_wall) result (wall3d)

  Function to return a pointer to a wall3d structure associated
  with a given lattice element.

  Note: The wall associated with a the vacuum chamber is the branch%wall3d.

  Parameters
  ----------
  ele : EleStruct
      lattice element.
  ix_wall : int, optional
      index in wall3d(:) array. Default is 1.

  Returns
  -------
  wall3d : Wall3DStruct
      Pointer to the associated wall structure. Will be nullified if there is no associated wall.
  ds_offset : float
      Element offset: s(beginning of ele) - s(beginning of wall3d)
  is_branch_wall : bool
      Set True if wall3d points to branch.wall3d.
  )""");
  m.def(
      "polar_to_spinor",
      &Bmad::polar_to_spinor,
      py::arg("polar"),
      py::arg("spinor"),
      R"""(Parameters
  ----------
  polar : SpinPolarStruct
      includes polar phase
  spinor : 
  )""");
  m.def(
      "polar_to_vec",
      &Bmad::polar_to_vec,
      py::arg("polar"),
      py::arg("vec"),
      R"""(Parameters
  ----------
  polar : 
      Spin_polar_struct
  vec : 
  )""");
  py::class_<Bmad::ProjectEmitToXyz, std::unique_ptr<Bmad::ProjectEmitToXyz>>(
      m, "ProjectEmitToXyz", "project_emit_to_xyz return type")
      .def_readonly("sigma_x", &Bmad::ProjectEmitToXyz::sigma_x)
      .def_readonly("sigma_y", &Bmad::ProjectEmitToXyz::sigma_y)
      .def_readonly("sigma_z", &Bmad::ProjectEmitToXyz::sigma_z)
      .def("__len__", [](const Bmad::ProjectEmitToXyz&) { return 3; })
      .def(
          "__getitem__",
          [](const Bmad::ProjectEmitToXyz& s, int i) -> py::object {
            if (i < 0)
              i += 3;
            if (i == 0)
              return py::cast(s.sigma_x);
            if (i == 1)
              return py::cast(s.sigma_y);
            if (i == 2)
              return py::cast(s.sigma_z);
            throw py::index_error();
          });
  m.def(
      "project_emit_to_xyz",
      &Bmad::project_emit_to_xyz,
      py::arg("ring"),
      py::arg("ix"),
      py::arg("mode"),
      R"""(Subroutine project_emit_to_xyz(ring, ix, mode, sigma_x, sigma_y, sigma_z)

  Obtains the projected x, y, and z beamsizes by building the sigma matrix
  from the normal mode emittances and 1-turn transfer matrix.
  These projectes beamsize are what would be seen by instrumentation.

  This method of projecting takes into account transverse and longitudinal coupling.

  This method of obtaining the projected beam sizes is from "Alternitive approach to general
  coupled linear optics" by Andrzej Wolski.

  The normal mode emittances used to generate a beam envelop sigma matrix from the
  1-turn transfer matrix.  The projected sizes are from the 1, 1 3, 3 and 5, 5 elements of
  the sigma matrix.

  Parameters
  ----------
  ring : LatStruct
      the storage ring
  ix : int
      element at which to make the projection
  mode : NormalModesStruct
      normal mode emittances .a.emittance -- real(rp): a-mode emittance .b.emittance -- real(rp): b-mode
      emittance .z.emittance -- real(rp): z-mode emittance .a.tune      -- real(rp): a-mode tune.  Used to
      associate emittances with the proper mode. .b.tune      -- real(rp): b-mode tune.  Used to associate
      emittances with the proper mode.

  Returns
  -------
  sigma_x : float
      projected horizontal beamsize
  sigma_y : float
      projected vertical beamsize
  sigma_z : float
      projected longitudinal beamsize
  )""");
  m.def(
      "psi_prime_sca",
      &Bmad::psi_prime_sca,
      py::arg("t"),
      py::arg("p"),
      py::arg("args"),
      R"""(Subroutine psi_prime_sca(t, p, dpdt, args)

  This wraps the array-valued psi_prime function as a scalar.

  See psi_prime comments for details.

  Parameters
  ----------
  t : float
      time relative to RF bucket
  p : float
      psi(t)
  args : float
      parameters and constants of DEQ

  Returns
  -------
  dpdt : float
      dpsi_dt
  )""");
  m.def(
      "ptc_bookkeeper",
      &Bmad::ptc_bookkeeper,
      py::arg("lat"),
      R"""(Parameters
  ----------
  lat : LatStruct
      Bmad lattice.
  )""");
  m.def(
      "ptc_closed_orbit_calc",
      &Bmad::ptc_closed_orbit_calc,
      py::arg("branch"),
      py::arg("radiation_damping_on") = py::none(),
      R"""(Subroutine ptc_closed_orbit_calc (branch, closed_orbit, radiation_damping_on)

  Routine to calculate the closed orbit of a lattice branch using PTC.
  This routine assumes the associated PTC layout has been crated
  with lat_to_ptc_layout.

  Parameters
  ----------
  branch : BranchStruct
      Branch of a lattice.
  radiation_damping_on : bool, optional
      If True, radiation dampling is included in the calculation. -- logical, optional: If True, radiation
      dampling is included in the calculation. Default is the setting of bmad_com..radiation_damping_on.

  Returns
  -------
  closed_orbit : CoordStruct
      closed_orbit
  )""");
  py::class_<Bmad::PtcEmitCalc, std::unique_ptr<Bmad::PtcEmitCalc>>(
      m, "PtcEmitCalc", "ptc_emit_calc return type")
      .def_readonly("norm_mode", &Bmad::PtcEmitCalc::norm_mode)
      .def_readonly("closed_orb", &Bmad::PtcEmitCalc::closed_orb)
      .def("__len__", [](const Bmad::PtcEmitCalc&) { return 2; })
      .def("__getitem__", [](const Bmad::PtcEmitCalc& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.norm_mode);
        if (i == 1)
          return py::cast(s.closed_orb);
        throw py::index_error();
      });
  m.def(
      "ptc_emit_calc",
      &Bmad::ptc_emit_calc,
      py::arg("ele"),
      py::arg("sigma_mat"),
      R"""(Subroutine ptc_emit_calc (ele, norm_mode, sigma_mat, closed_orb)

  Routine to calculate emittances, etc.

  Note: This routine calls the PTC init_all routine.

  Parameters
  ----------
  ele : EleStruct
      Element at which to evaluate the parameters.

  Returns
  -------
  norm_mode : 
      Normal_modes_struct %a%tune, %b%tune, %z%tune %a%alpha_damp, etc. %a%emittance, etc.
  sigma_map : float
      Sigma matrix (Bmad coordinates).
  closed_orb : CoordStruct
      Closed orbit at ele (Bmad coordinates). Notice: This closed orbit is the closed orbit with radiation on.
  )""");
  m.def(
      "ptc_layouts_resplit",
      &Bmad::ptc_layouts_resplit,
      py::arg("dKL_max"),
      py::arg("l_max"),
      py::arg("l_max_drift_only"),
      py::arg("bend_dorb"),
      py::arg("sex_dx"),
      py::arg("even") = py::none(),
      py::arg("crossover") = py::none(),
      py::arg("crossover_wiggler") = py::none(),
      R"""(Subroutine ptc_layouts_resplit (dKL_max, l_max, l_max_drift_only, bend_dorb, sex_dx,
                                                            even, crossover, crossover_wiggler)

  Routine to resplit (that is, recalculate the number of integration steps for an element)
  For the fibres in all layouts. After doing a resplit, the tune (and any other relavent
  "adjustable" parameters) should be adjusted to the correct values.

  Parameters
  ----------
  dKL_max : float
      Maximum K1 * L quadrupole strength allowed for an integration step. Reasonable value would be something
      like 0.04.
  l_max : float
      Maximum step length. Ignored if set to 0.
  l_max_drift_only : bool
      If True then l_max is only used for splitting drifts. -- logical: If True then l_max is only used for
      splitting drifts.
  bend_dorb : float
      Residual bend orbit error. With some integration methods a zero orbit at the start of the bend will not be
      zero at the end. In this case, bend_dorb sets a maximum allowable orbit deviation. If set to zero, this
      argument will be ignored. A resonable value is 10d-7. Note that the actual orbit deviation is not simply
      related to bend_dorb and can be larger. In any case, lowering bend_dorb (without making it zero) will
      lower the
  sex_dx : float
      To split sextupoles, sex_dx is used as the reference position about which the quadrupole strength is
      calculated. This quadrupole strength is then used with dKL_max to calculate the number of integration
      steps. Set to zero to ignore.
  even : bool, optional
      If True then each fibre  will have an even number of steps. If False then the number of steps will be odd.
      If not present then number of steps is not constrained to be even or odd.
  crossover : int, optional
      crossover(1) sets the maximum number of 2nd order integration steps to use. If the number of steps would
      exceed crossover(1) then integration is switched to 4th order. crossover(2) sets the maximum number of 4th
      order integration steps. If this number is exceeded, 6th order integration is used. Currently the default
      in PTC is [4, 18].
  crossover_wiggler(2) : int, optional
      crossover for wiggler elements. -- integer, optional: crossover for wiggler elements.
  )""");
  m.def(
      "ptc_one_turn_mat_and_closed_orbit_calc",
      &Bmad::ptc_one_turn_mat_and_closed_orbit_calc,
      py::arg("branch"),
      py::arg("pz") = py::none(),
      R"""(Subroutine ptc_one_turn_mat_and_closed_orbit_calc (branch, pz)

  Routine to compute the transfer matrices for the individual elements and closed orbit
  for a lattice branch with closed geometry.

  Note: PTC itself does not compute Twiss parameters. Use twiss_from_mat6 to compute this.

  Parameters
  ----------
  branch : BranchStruct
      Lattice branch.
      This parameter is an input/output and is modified in-place. As an output: Lattice branch containing the
      matrices.
  pz : float, optional
      energy offset around which to calculate the matrices if there is no RF.
  )""");
  m.def(
      "ptc_ran_seed_put",
      &Bmad::ptc_ran_seed_put,
      py::arg("iseed"),
      R"""(Parameters
  ----------
  iseed : int
      0 -> Use system clock.
  )""");
  m.def(
      "ptc_set_rf_state_for_c_normal",
      &Bmad::ptc_set_rf_state_for_c_normal,
      py::arg("nocavity"),
      R"""(Parameters
  ----------
  nocavity : bool
      True -> RF is off and vice versa.
  )""");
  m.def(
      "ptc_set_taylor_order_if_needed",
      &Bmad::ptc_set_taylor_order_if_needed,
      R"""(Subroutine ptc_set_taylor_order_if_needed()

  Routine to see if the taylor_order for PTC needs to be set/changed.
  For example, for a change in bmad_com%taylor_order.

  )""");
  py::class_<Bmad::PtcSpinCalc, std::unique_ptr<Bmad::PtcSpinCalc>>(
      m, "PtcSpinCalc", "ptc_spin_calc return type")
      .def_readonly("norm_mode", &Bmad::PtcSpinCalc::norm_mode)
      .def_readonly("closed_orb", &Bmad::PtcSpinCalc::closed_orb)
      .def("__len__", [](const Bmad::PtcSpinCalc&) { return 2; })
      .def("__getitem__", [](const Bmad::PtcSpinCalc& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.norm_mode);
        if (i == 1)
          return py::cast(s.closed_orb);
        throw py::index_error();
      });
  m.def(
      "ptc_spin_calc",
      &Bmad::ptc_spin_calc,
      py::arg("ele"),
      py::arg("sigma_mat"),
      R"""(Subroutine ptc_spin_calc (ele, norm_mode, sigma_mat, closed_orb)

  Routine to equilibrium polarizations, etc.

  Parameters
  ----------
  ele : EleStruct
      Element at which to evaluate the parameters.

  Returns
  -------
  norm_mode : 
      Normal_modes_struct %a%tune, %b%tune, %z%tune %a%alpha_damp, etc. %a%emittance, etc.
  sigma_map : float
      Sigma matrix (Bmad coordinates).
  closed_orb : CoordStruct
      Closed orbit at ele (Bmad coordinates). Notice: This closed orbit is the closed orbit with radiation on.
  )""");
  py::class_<Bmad::PtcTrackAll, std::unique_ptr<Bmad::PtcTrackAll>>(
      m, "PtcTrackAll", "ptc_track_all return type")
      .def_readonly("track_state", &Bmad::PtcTrackAll::track_state)
      .def_readonly("err_flag", &Bmad::PtcTrackAll::err_flag)
      .def("__len__", [](const Bmad::PtcTrackAll&) { return 2; })
      .def("__getitem__", [](const Bmad::PtcTrackAll& s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.track_state);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "ptc_track_all",
      &Bmad::ptc_track_all,
      py::arg("branch"),
      py::arg("orbit"),
      R"""(Subroutine ptc_track_all (branch, orbit, track_state, err_flag)

  Routine to track from the start to the end of a lattice branch.

  Parameters
  ----------
  branch : LatStruct
      Lat to track through.
  orbit : CoordStruct
      Coordinates at beginning of branch.
      This parameter is an input/output and is modified in-place. As an output: Orbit array.

  Returns
  -------
  track_state : int
      Set to moving_forward$ if everything is OK. Otherwise: set to index of element where particle was lost.
  err_flag : bool
      Set true if particle lost or error. False otherwise
  )""");
  m.def(
      "ptc_transfer_map_with_spin",
      &Bmad::ptc_transfer_map_with_spin,
      py::arg("branch"),
      py::arg("t_map"),
      py::arg("s_map"),
      py::arg("orb0"),
      py::arg("ix1") = py::none(),
      py::arg("ix2") = py::none(),
      py::arg("one_turn") = py::none(),
      py::arg("unit_start") = py::none(),
      R"""(Parameters
  ----------
  branch : BranchStruct
      Lattice branch used in the calculation.
  t_map : TaylorStruct
      Initial orbital map (used when unit_start = False)
      This parameter is an input/output and is modified in-place. As an output: Orbital transfer map.
  s_map : TaylorStruct
      Initial spin map (used when unit_start = False)
      This parameter is an input/output and is modified in-place. As an output: Quaternion spin transfer map.
  orb0 : CoordStruct
      Initial orbit around which the map is made.
  err_flag : bool
      Set True if problem like number overflow, etc.
  ix1 : int, optional
      Element start index for the calculation. Default is 0.
  ix2 : int, optional
      Element end index for the calculation. Default is branch.n_ele_track.
  one_turn : bool, optional
      If present and True, and if ix1 = ix2, and the lattice is circular, then construct the one-turn map from
      ix1 back to ix1. Default = False.
  unit_start : bool, optional
      If present and False then t_map will be used as the starting map instead of the unit map. Default = True
  )""");
  m.def(
      "pwd_mat",
      &Bmad::pwd_mat,
      py::arg("lat"),
      py::arg("t6"),
      py::arg("inductance"),
      py::arg("sig_z"),
      R"""(Function pwd_mat(t6, inductance, sig_z) result (t6_pwd)

  Calculates potential well distortion as RF defocusing.  Calculates t6_pwd=t6.Mpwd,
  where Mpwd is identity with 65 element proportional to the inductance.

  Vpwd = -inductance * lat%param%n_part * e_charge * c_light**3 / SQRT(twopi) / sig_z**3 / omega_RF  !effective RF voltage from PWD
  Mpwd(6,5) = omega_RF * Vpwd / c_light / lat%ele(0)%value(E_TOT$) * branch%ele(i)%value(l$) / lat%param%total_length

  Parameters
  ----------
  lat : LatStruct
      Bunch current in # per bunch .ele(0).value(E_TOT$) -- real(rp): Beam energy
  t6 : float
      1-turn transfer matrix
  inductance : float
      Longitudinal inductance in Henrys.  Something on the order of nH.
  sig_z : float
      Bunch length.

  Returns
  -------
  t6_pwd : float
      1-turn transfer matrix with PWD defocusing applied
  )""");
}
