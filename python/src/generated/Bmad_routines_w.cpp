#include "pybmad/generated/Bmad_routines_w.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

PyWriteLatLine python_write_lat_line(
    std::string line,
    int iu,
    bool end_is_neigh,
    std::optional<bool> do_split = std::nullopt,
    std::optional<bool> ampersand_at_ends = std::nullopt
) {
  Bmad::write_lat_line(line, iu, end_is_neigh, do_split, ampersand_at_ends);
  auto py_result{PyWriteLatLine{line}};
  return py_result;
}

void init_Bmad_routines_w(nb::module_ &m) {
  m.def(
      "w_mat_for_bend_angle",
      &Bmad::w_mat_for_bend_angle,
      nb::arg("angle"),
      nb::arg("ref_tilt"),
      nb::arg("r_vec") = nb::none(),
      R"""(Wrapper for Fortran routine w_mat_for_bend_angle

Parameters
----------
angle : float
    Bending angle.

ref_tilt : float
    Reference tilt.

r_vec : 1D array of float (shape: 3), optional
    Starting position.
    This parameter is an input/output and is modified in-place.
    As an output, r_vec: position with ref_tilt transformation

Returns
-------
w_mat : 2D array of float (shape: 3,3)
    W matrix
)"""
  );
  m.def(
      "w_mat_for_tilt",
      &Bmad::w_mat_for_tilt,
      nb::arg("tilt"),
      nb::arg("return_inverse") = nb::none(),
      R"""(Wrapper for Fortran routine w_mat_for_tilt

Parameters
----------
tilt : float
    pitch angle

return_inverse : bool, optional
    If True, return the inverse matrix. Default is False.

Returns
-------
w_mat : 2D array of float (shape: 3,3)
    Transformation matrix.
)"""
  );
  m.def(
      "w_mat_for_x_pitch",
      &Bmad::w_mat_for_x_pitch,
      nb::arg("x_pitch"),
      nb::arg("return_inverse") = nb::none(),
      R"""(Wrapper for Fortran routine w_mat_for_x_pitch

Parameters
----------
x_pitch : float
    pitch angle

return_inverse : bool, optional
    If True, return the inverse matrix. Default is False.

Returns
-------
w_mat : 2D array of float (shape: 3,3)
    Transformation matrix.
)"""
  );
  m.def(
      "w_mat_for_y_pitch",
      &Bmad::w_mat_for_y_pitch,
      nb::arg("y_pitch"),
      nb::arg("return_inverse") = nb::none(),
      R"""(Wrapper for Fortran routine w_mat_for_y_pitch

Parameters
----------
y_pitch : float
    pitch angle

return_inverse : bool, optional
    If True, return the inverse matrix. Default is False.

Returns
-------
w_mat : 2D array of float (shape: 3,3)
    Transformation matrix.
)"""
  );
  nb::class_<Bmad::Wall3dDRadius>(m, "Wall3dDRadius", "wall3d_d_radius return type")
      .def_ro("perp", &Bmad::Wall3dDRadius::perp)
      .def_ro("ix_section", &Bmad::Wall3dDRadius::ix_section)
      .def_ro("no_wall_here", &Bmad::Wall3dDRadius::no_wall_here)
      .def_ro("origin", &Bmad::Wall3dDRadius::origin)
      .def_ro("radius_wall", &Bmad::Wall3dDRadius::radius_wall)
      .def_ro("err_flag", &Bmad::Wall3dDRadius::err_flag)
      .def_ro("d_radius", &Bmad::Wall3dDRadius::d_radius)
      .def("__len__", [](const Bmad::Wall3dDRadius &) { return 7; })
      .def("__getitem__", [](const Bmad::Wall3dDRadius &s, int i) -> nb::object {
        if (i < 0)
          i += 7;
        if (i == 0)
          return nb::cast(s.perp);
        if (i == 1)
          return nb::cast(s.ix_section);
        if (i == 2)
          return nb::cast(s.no_wall_here);
        if (i == 3)
          return nb::cast(s.origin);
        if (i == 4)
          return nb::cast(s.radius_wall);
        if (i == 5)
          return nb::cast(s.err_flag);
        if (i == 6)
          return nb::cast(s.d_radius);
        throw nb::index_error();
      });
  m.def(
      "wall3d_d_radius",
      &Bmad::wall3d_d_radius,
      nb::arg("position"),
      nb::arg("ele"),
      nb::arg("ix_wall") = nb::none(),
      R"""(                                     no_wall_here, origin, radius_wall, err_flag) result (d_radius)

Routine to calculate the difference radius = particle_radius - wall_radius.
Radiuses are measured along a line from the wall origin with the line passing through
the particle position.
The wall origin itself lies on a line connecting the centers of the bounding sections.

Module needed:
  use wall3d_mod

Parameters
----------
position : 1D array of float
    Particle position in element coordinates. In a patch, with respect to entrance coords. [position(1),
    position(3)] = [x, y] transverse coords. position(5)                = Longitudinal position relative to
    beginning of element. position(6)                = Longitudinal velocity (only +/- sign matters).

ele : EleStruct
    Element with wall

ix_wall : int, optional
    Index of wall in .wall3d(:) array. Default is 1.

Returns
-------
d_radius : float
    r_particle - r_wall

perp : 1D array of float (shape: 3), optional
    Perpendicular normal to the wall.

ix_section : int, optional
    Set to wall slice section particle is in. That is between ix_section and ix_section+1.

no_wall_here : bool, optional
    True if the sub-chamber under consideration does not exist at the longitudinal location of the particle.

origin : 1D array of float (shape: 3), optional
    (x, y, s) origin with respect to the radius is measured. Uses the same coords as position.

radius_wall : float, optional
    Radius of the wall.

err_flag : bool, optional
    Set True if error. (EG noassociated .wall3d), false otherwise.
)"""
  );
  m.def(
      "wall3d_initializer",
      &Bmad::wall3d_initializer,
      nb::arg("wall3d"),
      R"""(Routine to initialize a wall3d_struct
  1) Add vertex points if there is symmetry.
  2) Compute circular and elliptical centers.
  3) Compute spline coefficients, etc.

Parameters
----------
wall3d : Wall3dStruct
    Wall.
    This parameter is an input/output and is modified in-place.
    As an output, wall3d: Initialized wall.

Returns
-------
err : bool
    Set true if there is a problem.
)"""
  );
  m.def(
      "wall3d_section_initializer",
      &Bmad::wall3d_section_initializer,
      nb::arg("section"),
      R"""(Routine to initialize a wall3d_section_struct:
  1) Add vertex points if there is symmetry.
  2) Compute circular and elliptical centers.

Parameters
----------
section : Wall3dSectionStruct
    Wall3d section.
    This parameter is an input/output and is modified in-place.
    As an output, section: Initialized section-section.

Returns
-------
err : bool
    Set true if there is a problem.
)"""
  );
  m.def(
      "wall3d_to_position",
      &Bmad::wall3d_to_position,
      nb::arg("orbit"),
      nb::arg("ele"),
      R"""(Routine to return the suitable postion to be used in calling wall3d_d_radius

This routine assumes that if in a patch the coordinates of orbit are with respect
to the downstream end if orbit%direction*orbit%time_dir = 1 and vice versa.

Parameters
----------
orbit : CoordStruct
    Particle position.

ele : EleStruct
    Element particle is in.

Returns
-------
position : 1D array of float (shape: 6)
    Position used in wall3d_d_radius call.
)"""
  );
  m.def(
      "word_to_value",
      [](std::string word, LatStruct &lat, double value, bool err_flag, EleStruct *ele) {
        auto fn =
            static_cast<void (*)(std::string, LatStruct &, double, bool, optional_ref<EleStruct>)>(
                &Bmad::word_to_value
            );
        return fn(word, lat, value, err_flag, ptr_to_opt_ref(ele));
      },
      nb::arg("word"),
      nb::arg("lat"),
      nb::arg("value"),
      nb::arg("err_flag"),
      nb::arg("ele") = nb::none(),
      R"""(Wrapper for Fortran routine word_to_value

Parameters
----------
word : str

lat : LatStruct

value : float

err_flag : bool

ele : EleStruct, optional
)"""
  );
  m.def(
      "write_ascii_beam_file",
      &Bmad::write_ascii_beam_file,
      nb::arg("file_name"),
      nb::arg("beam"),
      nb::arg("new_file") = nb::none(),
      nb::arg("alive_only") = nb::none(),
      R"""(Routine to write a beam file in ASCII format (version 4).

Parameters
----------
file_name : str
    Name of file.

beam : BeamStruct
    Beam to write

new_file : bool, optional
    New file or append? Default = True.

alive_only : bool, optional
    Only write live (includes pre_born) particles to the file? Default is False.
)"""
  );
  m.def(
      "write_astra_bend",
      &Bmad::write_astra_bend,
      nb::arg("iu"),
      nb::arg("strength"),
      nb::arg("id"),
      nb::arg("d1"),
      nb::arg("d2"),
      nb::arg("d3"),
      nb::arg("d4"),
      R"""(Wrapper for Fortran routine write_astra_bend

Parameters
----------
iu : int

strength : float

id : int

d1 : 1D array of float (shape: 2)

d2 : 1D array of float (shape: 2)

d3 : 1D array of float (shape: 2)

d4 : 1D array of float (shape: 2)
)"""
  );
  m.def(
      "write_astra_ele",
      [](int iu,
         EleStruct &ele,
         int id,
         StrIndexStruct *fieldgrid_names,
         std::optional<int> dimensions) {
        auto fn = static_cast<
            void (*)(int, EleStruct &, int, optional_ref<StrIndexStruct>, std::optional<int>)>(
            &Bmad::write_astra_ele
        );
        return fn(iu, ele, id, ptr_to_opt_ref(fieldgrid_names), dimensions);
      },
      nb::arg("iu"),
      nb::arg("ele"),
      nb::arg("id"),
      nb::arg("fieldgrid_names") = nb::none(),
      nb::arg("dimensions") = nb::none(),
      R"""(Wrapper for Fortran routine write_astra_ele

Parameters
----------
iu : int

ele : EleStruct

id : int

fieldgrid_names : StrIndexStruct, optional

dimensions : int, optional
)"""
  );
  nb::class_<Bmad::WriteAstraFieldGridFile>(
      m,
      "WriteAstraFieldGridFile",
      "write_astra_field_grid_file return type"
  )
      .def_ro("maxfield", &Bmad::WriteAstraFieldGridFile::maxfield)
      .def_ro("err", &Bmad::WriteAstraFieldGridFile::err)
      .def("__len__", [](const Bmad::WriteAstraFieldGridFile &) { return 2; })
      .def("__getitem__", [](const Bmad::WriteAstraFieldGridFile &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.maxfield);
        if (i == 1)
          return nb::cast(s.err);
        throw nb::index_error();
      });
  m.def(
      "write_astra_field_grid_file",
      &Bmad::write_astra_field_grid_file,
      nb::arg("astra_file_unit"),
      nb::arg("ele"),
      nb::arg("dz") = nb::none(),
      R"""(Write 1-D field map files for Astra. The format is:
z field
...

Note: Simplified from write_opal_field_grid_file

Parameters
----------
astra_file_unit : int
    unit number to write to, if > 0 if < 0, nothing is written, and only maxfield is returned

ele : EleStruct
    element to make map

dz : float, optional
    z step size in m. Default: 0.001 m

Returns
-------
maxfield : float
    absolute maximum found for element field scaling

err : bool, optional
    Set True if, say a file could not be opened.
)"""
  );
  nb::class_<Bmad::WriteAstraFieldGridFile3d>(
      m,
      "WriteAstraFieldGridFile3d",
      "write_astra_field_grid_file_3d return type"
  )
      .def_ro("maxfield", &Bmad::WriteAstraFieldGridFile3d::maxfield)
      .def_ro("err", &Bmad::WriteAstraFieldGridFile3d::err)
      .def("__len__", [](const Bmad::WriteAstraFieldGridFile3d &) { return 2; })
      .def("__getitem__", [](const Bmad::WriteAstraFieldGridFile3d &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.maxfield);
        if (i == 1)
          return nb::cast(s.err);
        throw nb::index_error();
      });
  m.def(
      "write_astra_field_grid_file_3d",
      &Bmad::write_astra_field_grid_file_3d,
      nb::arg("base_filename"),
      nb::arg("ele"),
      nb::arg("dz") = nb::none(),
      R"""(Writes 3-D field map files for Astra. The format is:
Nx x[1] x[2] ....... x[Nx-1] x[Nx]
Ny y[1] y[2] ....... y[Ny-1] y[Ny]
Nz z[1] z[2] ....... z[Nz-1] z[Nz]
<field values>
where field values are produced from a loop as in:
do iz = 1, Nz
  do iy = 1, Ny
    write single line: field(:, iy, iz)

Note: similar to write_astra_field_grid_file

Parameters
----------
base_filename : str
    Base filename. Files will be written as: base_filename.ex, .ey, .ez, .bx, .by, .bz If set to '', no files
    will be written

ele : EleStruct
    element to make map

dz : float, optional
    z step size in m. Default: 0.001 m

Returns
-------
maxfield : float
    absolute maximum on-axis field found for element field scaling

err : bool, optional
    Set True if, say a file could not be opened.
)"""
  );
  m.def(
      "write_astra_lattice_file",
      &Bmad::write_astra_lattice_file,
      nb::arg("astra_file_unit"),
      nb::arg("lat"),
      nb::arg("astra_lattice_param"),
      R"""(Subroutine to write an Astra lattice file using the information in a lat_struct.

Parameters
----------
astra_file_unit : int
    unit number to write to

lat : LatStruct
    Holds the lattice information.

Returns
-------
err : bool
    Set True if, say a file could not be opened.
)"""
  );
  m.def(
      "write_beam_file",
      [](std::string file_name,
         BeamStruct &beam,
         std::optional<bool> new_file,
         std::optional<int> file_format,
         LatStruct *lat,
         std::optional<bool> alive_only) {
        auto fn = static_cast<
            void (*)(std::string, BeamStruct &, std::optional<bool>, std::optional<int>, optional_ref<LatStruct>, std::optional<bool>)>(
            &Bmad::write_beam_file
        );
        return fn(file_name, beam, new_file, file_format, ptr_to_opt_ref(lat), alive_only);
      },
      nb::arg("file_name"),
      nb::arg("beam"),
      nb::arg("new_file") = nb::none(),
      nb::arg("file_format") = nb::none(),
      nb::arg("lat") = nb::none(),
      nb::arg("alive_only") = nb::none(),
      R"""(Routine to write a beam file.

A '.h5' suffix will be appended to the created file if hdf5$ format is used and file_name does not
already have a '.h5' or '.hdf5' suffix.

Parameters
----------
file_name : str
    Name of file.

beam : BeamStruct
    Beam to write

new_file : bool, optional
    New file or append? Default = True.

file_format : int, optional
    ascii$, or hdf5$ (default). old_ascii$ (deprecated) is still accepted.

lat : LatStruct, optional
    If present, lattice info will be writen to hdf5 files.

alive_only : bool, optional
    Only write live (includes pre_born) particles to the file? Default is False.
)"""
  );
  m.def(
      "write_beam_floor_positions",
      &Bmad::write_beam_floor_positions,
      nb::arg("file_name"),
      nb::arg("beam"),
      nb::arg("ele"),
      nb::arg("new_file") = nb::none(),
      R"""(Wrapper for Fortran routine write_beam_floor_positions

Parameters
----------
file_name : str
    Name of file.

beam : BeamStruct
    Beam to write

ele : EleStruct
    Element that the beam is at.

new_file : bool, optional
    New file or append? Default = True.
)"""
  );
  m.def(
      "write_binary_cartesian_map",
      &Bmad::write_binary_cartesian_map,
      nb::arg("file_name"),
      nb::arg("ele"),
      nb::arg("cart_map"),
      R"""(Routine to write a binary cartesian_map structure.
Note: The file name should have a ".bin" suffix.

Parameters
----------
file_name : str
    File to create.

ele : EleStruct
    Element associated with the map.

cart_map : CartesianMapStruct
    Cartesian map.

Returns
-------
err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "write_binary_cylindrical_map",
      &Bmad::write_binary_cylindrical_map,
      nb::arg("file_name"),
      nb::arg("ele"),
      nb::arg("cl_map"),
      R"""(Routine to write a binary cylindrical_map structure.
Note: The file name should have a ".bin" suffix.

Parameters
----------
file_name : str
    File to create.

ele : EleStruct
    Element associated with the map.

cl_map : CylindricalMapStruct
    Cylindrical map.

Returns
-------
err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "write_binary_grid_field",
      &Bmad::write_binary_grid_field,
      nb::arg("file_name"),
      nb::arg("ele"),
      nb::arg("g_field"),
      R"""(Routine to write a binary grid_field structure.
Note: The file name should have a ".bin" suffix.

Parameters
----------
file_name : str
    File to create.

ele : EleStruct
    Element associated with the map.

g_field : GridFieldStruct
    Cylindrical map.

Returns
-------
err_flag : bool
    Set True if there is an error. False otherwise.
)"""
  );
  m.def(
      "write_blender_ele",
      &Bmad::write_blender_ele,
      nb::arg("iu"),
      nb::arg("ele"),
      nb::arg("old_format") = nb::none(),
      R"""(Wrapper for Fortran routine write_blender_ele

Parameters
----------
iu : int

ele : EleStruct

old_format : bool, optional
)"""
  );
  m.def(
      "write_blender_lat_layout",
      &Bmad::write_blender_lat_layout,
      nb::arg("file_name"),
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine write_blender_lat_layout

Parameters
----------
file_name : str

lat : LatStruct
)"""
  );
  m.def(
      "write_bmad_lattice_file",
      [](std::string bmad_file, LatStruct &lat, std::optional<int> output_form, CoordStruct *orbit0
      ) {
        auto fn = static_cast<
            bool (*)(std::string, LatStruct &, std::optional<int>, optional_ref<CoordStruct>)>(
            &Bmad::write_bmad_lattice_file
        );
        return fn(bmad_file, lat, output_form, ptr_to_opt_ref(orbit0));
      },
      nb::arg("bmad_file"),
      nb::arg("lat"),
      nb::arg("output_form") = nb::none(),
      nb::arg("orbit0") = nb::none(),
      R"""(Wrapper for Fortran routine write_bmad_lattice_file

Parameters
----------
bmad_file : str
    Name of the output lattice file.

lat : LatStruct
    Holds the lattice information.

output_form : int, optional
    binary$   -> Write grid_field info in binary hdf5 form in separate files. Default. All other fields are
    writen in separate files in ASCII ascii$    -> Fields will be put in separate ASCII files. one_file$ ->
    Everything in one file.

orbit0 : CoordStruct, optional
    Initial orbit. Used to write the inital orbit if the lattice geometry is closed.

Returns
-------
err : bool, optional
    Set True if, say a file could not be opened.
)"""
  );
  m.def(
      "write_digested_bmad_file",
      [](std::string digested_name,
         LatStruct &lat,
         std::optional<int> n_files,
         CharacterAlloc1D *file_names,
         ExtraParsingInfoStruct *extra) {
        auto fn = static_cast<
            bool (*)(std::string, LatStruct &, std::optional<int>, optional_ref<CharacterAlloc1D>, optional_ref<ExtraParsingInfoStruct>)>(
            &Bmad::write_digested_bmad_file
        );
        return fn(digested_name, lat, n_files, ptr_to_opt_ref(file_names), ptr_to_opt_ref(extra));
      },
      nb::arg("digested_name"),
      nb::arg("lat"),
      nb::arg("n_files") = nb::none(),
      nb::arg("file_names") = nb::none(),
      nb::arg("extra") = nb::none(),
      R"""(Wrapper for Fortran routine write_digested_bmad_file

Parameters
----------
digested_name : str
    Name for the digested file.

lat : LatStruct
    Input lat structure.

n_files : int, optional
    Number of original files

file_names : 1D array of str, optional
    Names of the original files used to create the lat structure.

extra : ExtraParsingInfoStruct, optional
    Extra info that can be stored in the digested file.

Returns
-------
err_flag : bool, optional
    Set True if there is a problem. EG: No write permission. Set False if everything is OK.
)"""
  );
  m.def(
      "write_gpt_ele",
      [](int iu,
         EleStruct &ele,
         std::string name,
         int dimensions,
         StrIndexStruct *fieldgrid_names,
         std::optional<bool> only_phasing) {
        auto fn = static_cast<
            void (*)(int, EleStruct &, std::string, int, optional_ref<StrIndexStruct>, std::optional<bool>)>(
            &Bmad::write_gpt_ele
        );
        return fn(iu, ele, name, dimensions, ptr_to_opt_ref(fieldgrid_names), only_phasing);
      },
      nb::arg("iu"),
      nb::arg("ele"),
      nb::arg("name"),
      nb::arg("dimensions"),
      nb::arg("fieldgrid_names") = nb::none(),
      nb::arg("only_phasing") = nb::none(),
      R"""(Wrapper for Fortran routine write_gpt_ele

Parameters
----------
iu : int

ele : EleStruct

name : str

dimensions : int

fieldgrid_names : StrIndexStruct, optional

only_phasing : bool, optional
)"""
  );
  nb::class_<Bmad::WriteGptFieldGridFile1d>(
      m,
      "WriteGptFieldGridFile1d",
      "write_gpt_field_grid_file_1d return type"
  )
      .def_ro("maxfield", &Bmad::WriteGptFieldGridFile1d::maxfield)
      .def_ro("ref_time", &Bmad::WriteGptFieldGridFile1d::ref_time)
      .def_ro("err", &Bmad::WriteGptFieldGridFile1d::err)
      .def("__len__", [](const Bmad::WriteGptFieldGridFile1d &) { return 3; })
      .def("__getitem__", [](const Bmad::WriteGptFieldGridFile1d &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.maxfield);
        if (i == 1)
          return nb::cast(s.ref_time);
        if (i == 2)
          return nb::cast(s.err);
        throw nb::index_error();
      });
  m.def(
      "write_gpt_field_grid_file_1d",
      &Bmad::write_gpt_field_grid_file_1d,
      nb::arg("gpt_file_unit"),
      nb::arg("ele"),
      nb::arg("dz") = nb::none(),
      R"""(Write 1-D field map files for gpt. The format is:
z field
...

Note: Simplified from write_opal_field_grid_file

Parameters
----------
gpt_file_unit : int
    unit number to write to, if > 0 if < 0, nothing is written, and only maxfield is returned

ele : EleStruct
    element to make map

dz : float, optional
    z step size in m. Default: 0.001 m

Returns
-------
maxfield : float
    absolute maximum found for element field scaling

ref_time : float
    time that the field was evaluated at

err : bool, optional
    Set True if, say a file could not be opened.
)"""
  );
  nb::class_<Bmad::WriteGptFieldGridFile2d>(
      m,
      "WriteGptFieldGridFile2d",
      "write_gpt_field_grid_file_2d return type"
  )
      .def_ro("maxfield", &Bmad::WriteGptFieldGridFile2d::maxfield)
      .def_ro("ref_time", &Bmad::WriteGptFieldGridFile2d::ref_time)
      .def_ro("err", &Bmad::WriteGptFieldGridFile2d::err)
      .def("__len__", [](const Bmad::WriteGptFieldGridFile2d &) { return 3; })
      .def("__getitem__", [](const Bmad::WriteGptFieldGridFile2d &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.maxfield);
        if (i == 1)
          return nb::cast(s.ref_time);
        if (i == 2)
          return nb::cast(s.err);
        throw nb::index_error();
      });
  m.def(
      "write_gpt_field_grid_file_2d",
      &Bmad::write_gpt_field_grid_file_2d,
      nb::arg("gpt_file_unit"),
      nb::arg("ele"),
      nb::arg("dr") = nb::none(),
      nb::arg("dz") = nb::none(),
      nb::arg("r_max") = nb::none(),
      R"""(Subroutine to write an GPT lattice file using the information in
a lat_struct. Optionally only part of the lattice can be generated.

Parameters
----------
gpt_file_unit : int
    unit number to write to, if > 0 if < 0, nothing is written, and only maxfield is returned

ele : EleStruct
    element to make map

dr : float, optional
    r step size in m. Default: 0.001 m

dz : float, optional
    z step size in m. Default: 0.001 m

r_max : float, optional
    maximum radius in m. Default: 0.02 m

Returns
-------
maxfield : float
    absolute maximum found for element field scaling

ref_time : float
    time that the field was evaluated at

err : bool, optional
    Set True if, say a file could not be opened.
)"""
  );
  nb::class_<Bmad::WriteGptFieldGridFile3d>(
      m,
      "WriteGptFieldGridFile3d",
      "write_gpt_field_grid_file_3d return type"
  )
      .def_ro("maxfield", &Bmad::WriteGptFieldGridFile3d::maxfield)
      .def_ro("ref_time", &Bmad::WriteGptFieldGridFile3d::ref_time)
      .def_ro("err", &Bmad::WriteGptFieldGridFile3d::err)
      .def("__len__", [](const Bmad::WriteGptFieldGridFile3d &) { return 3; })
      .def("__getitem__", [](const Bmad::WriteGptFieldGridFile3d &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.maxfield);
        if (i == 1)
          return nb::cast(s.ref_time);
        if (i == 2)
          return nb::cast(s.err);
        throw nb::index_error();
      });
  m.def(
      "write_gpt_field_grid_file_3d",
      &Bmad::write_gpt_field_grid_file_3d,
      nb::arg("base_filename"),
      nb::arg("ele"),
      nb::arg("dz") = nb::none(),
      R"""(Writes 3-D field map files for gpt. The format is:

E-fields:
'x', 'y', 'z', 'ExRe', 'EyRe', 'EzRe', 'ExIm ', 'EyIm ', 'EzIm '
H-fields
'x', 'y', 'z', 'HxRe', 'HyRe', 'HzRe', 'HxIm ', 'HyIm ', 'HzIm '

where the fields oscillate as exp(+i \omega t)

Note: similar to write_gpt_field_grid_file

Parameters
----------
base_filename : str
    Base filename. Files will be written as: base_filename_E_ASCII.gpt, _H_ASCII.gpt If set to '', no files
    will be written

ele : EleStruct
    element to make map

dz : float, optional
    z step size in m. Default: 0.001 m

Returns
-------
maxfield : float
    absolute maximum on-axis field found for element field scaling

ref_time : float
    time that the field was evaluated at

err : bool, optional
    Set True if, say a file could not be opened.
)"""
  );
  m.def(
      "write_gpt_lattice_file",
      &Bmad::write_gpt_lattice_file,
      nb::arg("lat"),
      nb::arg("gpt_lat_param"),
      R"""(Subroutine to write a gpt lattice file using the information in a Bmad lattice.

Parameters
----------
lat : LatStruct
    Holds the lattice information.

gpt_lat_param : GptLatParamStruct
    Parameters for constructing the lattice

Returns
-------
err : bool
    Set True if there is an error
)"""
  );
  nb::class_<PyWriteLatLine>(m, "WriteLatLine", "write_lat_line return type")
      .def_ro("line", &PyWriteLatLine::line)
      .def("__len__", [](const PyWriteLatLine &) { return 1; })
      .def("__getitem__", [](const PyWriteLatLine &s, int i) -> nb::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return nb::cast(s.line);
        throw nb::index_error();
      });
  m.def(
      "write_lat_line",
      &python_write_lat_line,
      nb::arg("line"),
      nb::arg("iu"),
      nb::arg("end_is_neigh"),
      nb::arg("do_split") = nb::none(),
      nb::arg("ampersand_at_ends") = nb::none(),
      R"""(Routine to write strings to a lattice file.
This routine will break the string up into multiple lines
if the string is too long and add a continuation character if needed.

If the "line" arg does not represent a full "sentence" (end_is_neigh = False),
then only part of the line may be written and the part not written will be returned.

Parameters
----------
line : str
    String of text.
    This parameter is an input/output and is modified in-place.
    As an output, line: part of the string not written.

iu : int
    Unit number to write to.

end_is_neigh : bool
    If true then write out everything. Otherwise wait for a full line of max_char characters or so.

do_split : bool, optional
    Split line if overlength? Default is True. False is used when line has already been split for expressions
    since the expression splitting routine does a much better job of it.

ampersand_at_ends : bool, optional
    Default True. If False then do not include "&" line continuation

Returns
-------
line : str
    String of text.
    This parameter is an input/output and is modified in-place.
    As an output, line: part of the string not written.
)"""
  );
  m.def(
      "write_lattice_elegant_format",
      &Bmad::write_lattice_elegant_format,
      nb::arg("out_file_name"),
      nb::arg("lat"),
      nb::arg("ref_orbit") = nb::none(),
      nb::arg("use_matrix_model") = nb::none(),
      nb::arg("include_apertures") = nb::none(),
      nb::arg("dr12_drift_max") = nb::none(),
      nb::arg("ix_branch") = nb::none(),
      R"""(Wrapper for Fortran routine write_lattice_elegant_format

Parameters
----------
out_file_name : str
    Name of the mad output lattice file.

lat : LatStruct
    Holds the lattice information.

ref_orbit : 1D array of CoordStruct, optional
    Referece orbit for sad_mult and patch elements. This argument must be present if the lattice has sad_mult
    or patch elements and is being translated to MAD-8 or SAD.

use_matrix_model : bool, optional
    Use a drift-matrix_drift model for wigglers/undulators? [A MAD "matrix" is a 2nd order Taylor map.] This
    switch is ignored for SAD conversion. Default is False -> Use a bend-drift-bend model. Note: sol_quad
    elements always use a drift-matrix-drift model.

include_apertures : bool, optional
    If True (the default), add to the output lattice a zero length collimator element next to any non-
    collimator element that has an aperture. Note: MADX translations for non-drift elements can handle non-
    collimator elements with an aperture so in this case this argument is ignored.

dr12_drift_max : float, optional
    Max deviation for drifts allowed before a correction matrix element is added. Default value is 1d-5. A
    negative number means use default.

ix_branch : int, optional
    Index of lattice branch to use. Default = 0.

Returns
-------
err : bool, optional
    Set True if, say a file could not be opened.
)"""
  );
  m.def(
      "write_lattice_foreign_format",
      &Bmad::write_lattice_foreign_format,
      nb::arg("out_type"),
      nb::arg("out_file_name"),
      nb::arg("lat"),
      nb::arg("ref_orbit") = nb::none(),
      nb::arg("use_matrix_model") = nb::none(),
      nb::arg("include_apertures") = nb::none(),
      nb::arg("dr12_drift_max") = nb::none(),
      nb::arg("ix_branch") = nb::none(),
      R"""(Wrapper for Fortran routine write_lattice_foreign_format

Parameters
----------
out_type : str
    Either 'ELEGANT', 'MAD-8', 'MAD-X', 'SAD', 'OPAL-T', 'PALS', or 'SCIBMAD'.

out_file_name : str
    Name of the mad output lattice file.

lat : LatStruct
    Holds the lattice information.

ref_orbit : 1D array of CoordStruct, optional
    Referece orbit for sad_mult and patch elements. This argument must be present if the lattice has sad_mult
    or patch elements and is being translated to MAD-8 or SAD.

use_matrix_model : bool, optional
    Use a drift-matrix_drift model for wigglers/undulators? [A MAD "matrix" is a 2nd order Taylor map.] This
    switch is ignored for SAD conversion. Default is False -> Use a bend-drift-bend model. Note: sol_quad
    elements always use a drift-matrix-drift model.

include_apertures : bool, optional
    If True (the default), add to the output lattice a zero length collimator element next to any non-
    collimator element that has an aperture. Note: MADX translations for non-drift elements can handle non-
    collimator elements with an aperture so in this case this argument is ignored.

dr12_drift_max : float, optional
    Max deviation for drifts allowed before a correction matrix element is added. Default value is 1d-5. A
    negative number means use default.

ix_branch : int, optional
    Index of lattice branch to use. Default = 0.

Returns
-------
err : bool, optional
    Set True if, say a file could not be opened.
)"""
  );
  m.def(
      "write_lattice_mad_format",
      &Bmad::write_lattice_mad_format,
      nb::arg("out_type"),
      nb::arg("out_file_name"),
      nb::arg("lat"),
      nb::arg("ref_orbit") = nb::none(),
      nb::arg("use_matrix_model") = nb::none(),
      nb::arg("include_apertures") = nb::none(),
      nb::arg("dr12_drift_max") = nb::none(),
      nb::arg("ix_branch") = nb::none(),
      R"""(Wrapper for Fortran routine write_lattice_mad_format

Parameters
----------
out_type : str
    Either 'MAD-8', or 'MAD-X'

out_file_name : str
    Name of the mad output lattice file.

lat : LatStruct
    Holds the lattice information.

ref_orbit : 1D array of CoordStruct, optional
    Referece orbit for sad_mult and patch elements. This argument must be present if the lattice has sad_mult
    or patch elements and is being translated to MAD-8 or SAD.

use_matrix_model : bool, optional
    Use a drift-matrix_drift model for wigglers/undulators? [A MAD "matrix" is a 2nd order Taylor map.] This
    switch is ignored for SAD conversion. Default is False -> Use a bend-drift-bend model. Note: sol_quad
    elements always use a drift-matrix-drift model.

include_apertures : bool, optional
    If True (the default), add to the output lattice a zero length collimator element next to any non-
    collimator element that has an aperture. Note: MADX translations for non-drift elements can handle non-
    collimator elements with an aperture so in this case this argument is ignored.

dr12_drift_max : float, optional
    Max deviation for drifts allowed before a correction matrix element is added. Default value is 1d-5. A
    negative number means use default.

ix_branch : int, optional
    Index of lattice branch to use. Default = 0.

Returns
-------
err : bool, optional
    Set True if, say a file could not be opened.
)"""
  );
  nb::class_<Bmad::WriteLatticePalsFormat>(
      m,
      "WriteLatticePalsFormat",
      "write_lattice_pals_format return type"
  )
      .def_ro("pals_file", &Bmad::WriteLatticePalsFormat::pals_file)
      .def_ro("err_flag", &Bmad::WriteLatticePalsFormat::err_flag)
      .def("__len__", [](const Bmad::WriteLatticePalsFormat &) { return 2; })
      .def("__getitem__", [](const Bmad::WriteLatticePalsFormat &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.pals_file);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "write_lattice_pals_format",
      &Bmad::write_lattice_pals_format,
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine write_lattice_pals_format

Parameters
----------
lat : LatStruct
    Lattice

Returns
-------
pals_file : str
    Pals lattice file name.

err_flag : bool, optional
    Error flag
)"""
  );
  m.def(
      "write_lattice_sad_format",
      &Bmad::write_lattice_sad_format,
      nb::arg("out_file_name"),
      nb::arg("lat"),
      nb::arg("include_apertures") = nb::none(),
      nb::arg("ix_branch") = nb::none(),
      nb::arg("err") = nb::none(),
      R"""(Wrapper for Fortran routine write_lattice_sad_format

Parameters
----------
out_file_name : str

lat : LatStruct

include_apertures : bool, optional

ix_branch : int, optional

err : bool, optional
)"""
  );
  nb::class_<Bmad::WriteLatticeScibmadFormat>(
      m,
      "WriteLatticeScibmadFormat",
      "write_lattice_scibmad_format return type"
  )
      .def_ro("scibmad_file", &Bmad::WriteLatticeScibmadFormat::scibmad_file)
      .def_ro("err_flag", &Bmad::WriteLatticeScibmadFormat::err_flag)
      .def("__len__", [](const Bmad::WriteLatticeScibmadFormat &) { return 2; })
      .def("__getitem__", [](const Bmad::WriteLatticeScibmadFormat &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.scibmad_file);
        if (i == 1)
          return nb::cast(s.err_flag);
        throw nb::index_error();
      });
  m.def(
      "write_lattice_scibmad_format",
      &Bmad::write_lattice_scibmad_format,
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine write_lattice_scibmad_format

Parameters
----------
lat : LatStruct
    Lattice

Returns
-------
scibmad_file : str
    SciBmad lattice file name.

err_flag : bool, optional
    Error flag
)"""
  );
  m.def(
      "write_line_element",
      &Bmad::write_line_element,
      nb::arg("line"),
      nb::arg("iu"),
      nb::arg("ele"),
      nb::arg("lat"),
      R"""(Wrapper for Fortran routine write_line_element

Parameters
----------
line : str

iu : int

ele : EleStruct

lat : LatStruct
)"""
  );
  nb::class_<Bmad::WriteOpalFieldGridFile>(
      m,
      "WriteOpalFieldGridFile",
      "write_opal_field_grid_file return type"
  )
      .def_ro("maxfield", &Bmad::WriteOpalFieldGridFile::maxfield)
      .def_ro("err", &Bmad::WriteOpalFieldGridFile::err)
      .def("__len__", [](const Bmad::WriteOpalFieldGridFile &) { return 2; })
      .def("__getitem__", [](const Bmad::WriteOpalFieldGridFile &s, int i) -> nb::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return nb::cast(s.maxfield);
        if (i == 1)
          return nb::cast(s.err);
        throw nb::index_error();
      });
  m.def(
      "write_opal_field_grid_file",
      &Bmad::write_opal_field_grid_file,
      nb::arg("opal_file_unit"),
      nb::arg("ele"),
      nb::arg("param"),
      R"""(Subroutine to write an OPAL lattice file using the information in
a lat_struct. Optionally only part of the lattice can be generated.

Parameters
----------
opal_file_unit : int
    unit number to write to, if > 0 if < 0, nothing is written, and only maxfield is returned

ele : EleStruct
    element to make map

param : LatParamStruct
    Contains lattice information

Returns
-------
maxfield : float
    absolute maximum found for element field scaling

err : bool, optional
    Set True if, say a file could not be opened.
)"""
  );
  m.def(
      "write_opal_lattice_file",
      &Bmad::write_opal_lattice_file,
      nb::arg("opal_file_unit"),
      nb::arg("lat"),
      R"""(Subroutine to write an OPAL lattice file using the information in
a lat_struct. Optionally only part of the lattice can be generated.

Parameters
----------
opal_file_unit : int
    unit number to write to

lat : LatStruct
    Holds the lattice information.

Returns
-------
err : bool, optional
    Set True if, say a file could not be opened.
)"""
  );
  m.def(
      "write_time_particle_distribution",
      [](int time_file_unit,
         BunchStruct &bunch,
         EleStruct &ele,
         std::optional<std::string> style,
         BranchStruct *branch,
         std::optional<std::string> format) {
        auto fn = static_cast<
            bool (*)(int, BunchStruct &, EleStruct &, std::optional<std::string>, optional_ref<BranchStruct>, std::optional<std::string>)>(
            &Bmad::write_time_particle_distribution
        );
        return fn(time_file_unit, bunch, ele, style, ptr_to_opt_ref(branch), format);
      },
      nb::arg("time_file_unit"),
      nb::arg("bunch"),
      nb::arg("ele"),
      nb::arg("style") = nb::none(),
      nb::arg("branch") = nb::none(),
      nb::arg("format") = nb::none(),
      R"""(Subroutine to write a time-based bunch from a standard Bmad bunch

Note: 'BMAD' style (absolute curvilinear coordinates):
      n_particles_alive
      x/m  m*c^2 \beta_x*\gamma/eV y/m m*c^2\beta_y*\gamma/eV s/m m*c^2\beta_z*\gamma/eV time/s charge/C

      'OPAL' style (absolute curvilinear coordinates):
      n_particles_alive
      x/m  \beta_x*\gamma  y/m \beta_y*\gamma s/m \beta_s*\gamma

      'ASTRA' style (global Cartesian coordinates, first line is the reference particle used for z, pz, and t calculation):
      x/m y/m  z/m  m*c^2 \beta_x*\gamma/eV m*c^2 \beta_y*\gamma/eV m*c^2 \beta_z*\gamma/eV time/ns charge/nC species status

      'GPT' style (global Cartesian coordinates, with header labeling the columns)
      x/m y/m z/m \beta_x*\gamma \beta_y*\gamma \beta_z*\gamma t/s elementary_charge/C charge/elementary_charge

Parameters
----------
time_file_unit : int
    unit number to write to, if > 0

bunch : BunchStruct
    bunch to be written. Particles are drifted to bmad_bunch.t_center for output

ele : EleStruct
    Element being tracked through.

style : str, optional
    Style of output file: 'BMAD' (default), 'OPAL', 'ASTRA', 'GPT'

branch : BranchStruct, optional
    Required for 'ASTRA' style

format : str, optional
    format for numerical output. default: 'es15.7'

Returns
-------
err : bool, optional
    Set True if, say a file could not be opened.
)"""
  );
}
