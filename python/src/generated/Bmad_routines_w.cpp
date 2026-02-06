#include "pybmad/generated/Bmad_routines_w.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

PyWriteLatLine python_write_lat_line(
    std::string line,
    int iu,
    bool end_is_neigh,
    std::optional<bool> do_split = std::nullopt,
    std::optional<bool> scibmad = std::nullopt
) {
  Bmad::write_lat_line(line, iu, end_is_neigh, do_split, scibmad);
  auto py_result{PyWriteLatLine{line}};
  return py_result;
}

void init_Bmad_routines_w(py::module &m) {
  m.def(
      "w_mat_for_bend_angle",
      &Bmad::w_mat_for_bend_angle,
      py::arg("angle"),
      py::arg("ref_tilt"),
      py::arg("r_vec") = py::none(),
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
      py::arg("tilt"),
      py::arg("return_inverse") = py::none(),
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
      py::arg("x_pitch"),
      py::arg("return_inverse") = py::none(),
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
      py::arg("y_pitch"),
      py::arg("return_inverse") = py::none(),
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
  py::class_<Bmad::Wall3dDRadius, std::unique_ptr<Bmad::Wall3dDRadius>>(
      m,
      "Wall3dDRadius",
      "wall3d_d_radius return type"
  )
      .def_readonly("perp", &Bmad::Wall3dDRadius::perp)
      .def_readonly("ix_section", &Bmad::Wall3dDRadius::ix_section)
      .def_readonly("no_wall_here", &Bmad::Wall3dDRadius::no_wall_here)
      .def_readonly("origin", &Bmad::Wall3dDRadius::origin)
      .def_readonly("radius_wall", &Bmad::Wall3dDRadius::radius_wall)
      .def_readonly("err_flag", &Bmad::Wall3dDRadius::err_flag)
      .def_readonly("d_radius", &Bmad::Wall3dDRadius::d_radius)
      .def("__len__", [](const Bmad::Wall3dDRadius &) { return 7; })
      .def("__getitem__", [](const Bmad::Wall3dDRadius &s, int i) -> py::object {
        if (i < 0)
          i += 7;
        if (i == 0)
          return py::cast(s.perp);
        if (i == 1)
          return py::cast(s.ix_section);
        if (i == 2)
          return py::cast(s.no_wall_here);
        if (i == 3)
          return py::cast(s.origin);
        if (i == 4)
          return py::cast(s.radius_wall);
        if (i == 5)
          return py::cast(s.err_flag);
        if (i == 6)
          return py::cast(s.d_radius);
        throw py::index_error();
      });
  m.def(
      "wall3d_d_radius",
      &Bmad::wall3d_d_radius,
      py::arg("position"),
      py::arg("ele"),
      py::arg("ix_wall") = py::none(),
      R"""(Function wall3d_d_radius (position, ele, ix_wall, perp, ix_section,
                                     no_wall_here, origin, radius_wall, err_flag) result (d_radius)

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
      py::arg("wall3d"),
      R"""(Subroutine wall3d_initializer (wall3d, err)

Routine to initialize a wall3d_struct
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
      py::arg("section"),
      R"""(Subroutine wall3d_section_initializer (section, err)

Routine to initialize a wall3d_section_struct:
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
      py::arg("orbit"),
      py::arg("ele"),
      R"""(Function wall3d_to_position (orbit, ele) result (position)

Routine to return the suitable postion to be used in calling wall3d_d_radius

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
      &Bmad::word_to_value,
      py::arg("word"),
      py::arg("lat"),
      py::arg("value"),
      py::arg("err_flag"),
      py::arg("ele") = py::none(),
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
      py::arg("file_name"),
      py::arg("beam"),
      py::arg("new_file") = py::none(),
      py::arg("alive_only") = py::none(),
      R"""(Subroutine write_ascii_beam_file (file_name, beam, new_file, alive_only)

Routine to write a beam file in ASCII format (version 4).

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
      py::arg("iu"),
      py::arg("strength"),
      py::arg("id"),
      py::arg("d1"),
      py::arg("d2"),
      py::arg("d3"),
      py::arg("d4"),
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
  py::class_<Bmad::WriteAstraFieldGridFile, std::unique_ptr<Bmad::WriteAstraFieldGridFile>>(
      m,
      "WriteAstraFieldGridFile",
      "write_astra_field_grid_file return type"
  )
      .def_readonly("maxfield", &Bmad::WriteAstraFieldGridFile::maxfield)
      .def_readonly("err", &Bmad::WriteAstraFieldGridFile::err)
      .def("__len__", [](const Bmad::WriteAstraFieldGridFile &) { return 2; })
      .def("__getitem__", [](const Bmad::WriteAstraFieldGridFile &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.maxfield);
        if (i == 1)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "write_astra_field_grid_file",
      &Bmad::write_astra_field_grid_file,
      py::arg("astra_file_unit"),
      py::arg("ele"),
      py::arg("dz") = py::none(),
      R"""(Subroutine write_astra_field_grid_file (astra_file_unit, ele, maxfield, err)

  Write 1-D field map files for Astra. The format is:
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
  py::class_<Bmad::WriteAstraFieldGridFile3d, std::unique_ptr<Bmad::WriteAstraFieldGridFile3d>>(
      m,
      "WriteAstraFieldGridFile3d",
      "write_astra_field_grid_file_3d return type"
  )
      .def_readonly("maxfield", &Bmad::WriteAstraFieldGridFile3d::maxfield)
      .def_readonly("err", &Bmad::WriteAstraFieldGridFile3d::err)
      .def("__len__", [](const Bmad::WriteAstraFieldGridFile3d &) { return 2; })
      .def("__getitem__", [](const Bmad::WriteAstraFieldGridFile3d &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.maxfield);
        if (i == 1)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "write_astra_field_grid_file_3d",
      &Bmad::write_astra_field_grid_file_3d,
      py::arg("base_filename"),
      py::arg("ele"),
      py::arg("dz") = py::none(),
      R"""(Subroutine write_astra_field_grid_file_3D (base_filename, ele, maxfield, dz, err)

  Writes 3-D field map files for Astra. The format is:
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
      "write_beam_file",
      &Bmad::write_beam_file,
      py::arg("file_name"),
      py::arg("beam"),
      py::arg("new_file") = py::none(),
      py::arg("file_format") = py::none(),
      py::arg("lat") = py::none(),
      py::arg("alive_only") = py::none(),
      R"""(Subroutine write_beam_file (file_name, beam, new_file, file_format, lat, alive_only)

Routine to write a beam file.

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
      py::arg("file_name"),
      py::arg("beam"),
      py::arg("ele"),
      py::arg("new_file") = py::none(),
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
      py::arg("file_name"),
      py::arg("ele"),
      py::arg("cart_map"),
      R"""(Subroutine write_binary_cartesian_map (file_name, ele, cart_map, err_flag)

Routine to write a binary cartesian_map structure.
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
      py::arg("file_name"),
      py::arg("ele"),
      py::arg("cl_map"),
      R"""(Subroutine write_binary_cylindrical_map (file_name, ele, cl_map, err_flag)

Routine to write a binary cylindrical_map structure.
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
      py::arg("file_name"),
      py::arg("ele"),
      py::arg("g_field"),
      R"""(Subroutine write_binary_grid_field (file_name, ele, g_field, err_flag)

Routine to write a binary grid_field structure.
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
      py::arg("iu"),
      py::arg("ele"),
      py::arg("old_format") = py::none(),
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
      py::arg("file_name"),
      py::arg("lat"),
      R"""(Wrapper for Fortran routine write_blender_lat_layout

Parameters
----------
file_name : str

lat : LatStruct
)"""
  );
  m.def(
      "write_bmad_lattice_file",
      &Bmad::write_bmad_lattice_file,
      py::arg("bmad_file"),
      py::arg("lat"),
      py::arg("output_form") = py::none(),
      py::arg("orbit0") = py::none(),
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
  py::class_<Bmad::WriteGptFieldGridFile1d, std::unique_ptr<Bmad::WriteGptFieldGridFile1d>>(
      m,
      "WriteGptFieldGridFile1d",
      "write_gpt_field_grid_file_1d return type"
  )
      .def_readonly("maxfield", &Bmad::WriteGptFieldGridFile1d::maxfield)
      .def_readonly("ref_time", &Bmad::WriteGptFieldGridFile1d::ref_time)
      .def_readonly("err", &Bmad::WriteGptFieldGridFile1d::err)
      .def("__len__", [](const Bmad::WriteGptFieldGridFile1d &) { return 3; })
      .def("__getitem__", [](const Bmad::WriteGptFieldGridFile1d &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.maxfield);
        if (i == 1)
          return py::cast(s.ref_time);
        if (i == 2)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "write_gpt_field_grid_file_1d",
      &Bmad::write_gpt_field_grid_file_1d,
      py::arg("gpt_file_unit"),
      py::arg("ele"),
      py::arg("dz") = py::none(),
      R"""(Subroutine write_gpt_field_grid_file_1D (gpt_file_unit, ele, maxfield, ref_time, dz, err)

  Write 1-D field map files for gpt. The format is:
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
  py::class_<Bmad::WriteGptFieldGridFile2d, std::unique_ptr<Bmad::WriteGptFieldGridFile2d>>(
      m,
      "WriteGptFieldGridFile2d",
      "write_gpt_field_grid_file_2d return type"
  )
      .def_readonly("maxfield", &Bmad::WriteGptFieldGridFile2d::maxfield)
      .def_readonly("ref_time", &Bmad::WriteGptFieldGridFile2d::ref_time)
      .def_readonly("err", &Bmad::WriteGptFieldGridFile2d::err)
      .def("__len__", [](const Bmad::WriteGptFieldGridFile2d &) { return 3; })
      .def("__getitem__", [](const Bmad::WriteGptFieldGridFile2d &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.maxfield);
        if (i == 1)
          return py::cast(s.ref_time);
        if (i == 2)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "write_gpt_field_grid_file_2d",
      &Bmad::write_gpt_field_grid_file_2d,
      py::arg("gpt_file_unit"),
      py::arg("ele"),
      py::arg("dr") = py::none(),
      py::arg("dz") = py::none(),
      py::arg("r_max") = py::none(),
      R"""(Subroutine write_gpt_field_grid_file_2D (gpt_file_unit, ele, maxfield, ref_time, dr, dz,  err)

Subroutine to write an GPT lattice file using the information in
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
  py::class_<Bmad::WriteGptFieldGridFile3d, std::unique_ptr<Bmad::WriteGptFieldGridFile3d>>(
      m,
      "WriteGptFieldGridFile3d",
      "write_gpt_field_grid_file_3d return type"
  )
      .def_readonly("maxfield", &Bmad::WriteGptFieldGridFile3d::maxfield)
      .def_readonly("ref_time", &Bmad::WriteGptFieldGridFile3d::ref_time)
      .def_readonly("err", &Bmad::WriteGptFieldGridFile3d::err)
      .def("__len__", [](const Bmad::WriteGptFieldGridFile3d &) { return 3; })
      .def("__getitem__", [](const Bmad::WriteGptFieldGridFile3d &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.maxfield);
        if (i == 1)
          return py::cast(s.ref_time);
        if (i == 2)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "write_gpt_field_grid_file_3d",
      &Bmad::write_gpt_field_grid_file_3d,
      py::arg("base_filename"),
      py::arg("ele"),
      py::arg("dz") = py::none(),
      R"""(Subroutine write_gpt_field_grid_file_3D (base_filename, ele, maxfield, ref_time, dz, err)

  Writes 3-D field map files for gpt. The format is:

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
  py::class_<PyWriteLatLine, std::unique_ptr<PyWriteLatLine>>(
      m,
      "WriteLatLine",
      "write_lat_line return type"
  )
      .def_readonly("line", &PyWriteLatLine::line)
      .def("__len__", [](const PyWriteLatLine &) { return 1; })
      .def("__getitem__", [](const PyWriteLatLine &s, int i) -> py::object {
        if (i < 0)
          i += 1;
        if (i == 0)
          return py::cast(s.line);
        throw py::index_error();
      });
  m.def(
      "write_lat_line",
      &python_write_lat_line,
      py::arg("line"),
      py::arg("iu"),
      py::arg("end_is_neigh"),
      py::arg("do_split") = py::none(),
      py::arg("scibmad") = py::none(),
      R"""(Subroutine write_lat_line (line, iu, end_is_neigh, do_split)

Routine to write strings to a lattice file.
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

scibmad : bool, optional
    Default False. If True then do not include "&" line continuation

Returns
-------
line : str
    String of text.
    This parameter is an input/output and is modified in-place.
    As an output, line: part of the string not written.
)"""
  );
  m.def(
      "write_lattice_in_elegant_format",
      &Bmad::write_lattice_in_elegant_format,
      py::arg("out_file_name"),
      py::arg("lat"),
      py::arg("ref_orbit") = py::none(),
      py::arg("use_matrix_model") = py::none(),
      py::arg("include_apertures") = py::none(),
      py::arg("dr12_drift_max") = py::none(),
      py::arg("ix_branch") = py::none(),
      R"""(Wrapper for Fortran routine write_lattice_in_elegant_format

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
      "write_lattice_in_foreign_format",
      &Bmad::write_lattice_in_foreign_format,
      py::arg("out_type"),
      py::arg("out_file_name"),
      py::arg("lat"),
      py::arg("ref_orbit") = py::none(),
      py::arg("use_matrix_model") = py::none(),
      py::arg("include_apertures") = py::none(),
      py::arg("dr12_drift_max") = py::none(),
      py::arg("ix_branch") = py::none(),
      R"""(Wrapper for Fortran routine write_lattice_in_foreign_format

Parameters
----------
out_type : str
    Either 'ELEGANT', 'MAD-8', 'MAD-X', 'SAD', or 'OPAL-T', 'SCIBMAD'.

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
      "write_lattice_in_mad_format",
      &Bmad::write_lattice_in_mad_format,
      py::arg("out_type"),
      py::arg("out_file_name"),
      py::arg("lat"),
      py::arg("ref_orbit") = py::none(),
      py::arg("use_matrix_model") = py::none(),
      py::arg("include_apertures") = py::none(),
      py::arg("dr12_drift_max") = py::none(),
      py::arg("ix_branch") = py::none(),
      R"""(Wrapper for Fortran routine write_lattice_in_mad_format

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
  m.def(
      "write_lattice_in_sad_format",
      &Bmad::write_lattice_in_sad_format,
      py::arg("out_file_name"),
      py::arg("lat"),
      py::arg("include_apertures") = py::none(),
      py::arg("ix_branch") = py::none(),
      py::arg("err") = py::none(),
      R"""(Wrapper for Fortran routine write_lattice_in_sad_format

Parameters
----------
out_file_name : str

lat : LatStruct

include_apertures : bool, optional

ix_branch : int, optional

err : bool, optional
)"""
  );
  py::class_<Bmad::WriteLatticeInScibmad, std::unique_ptr<Bmad::WriteLatticeInScibmad>>(
      m,
      "WriteLatticeInScibmad",
      "write_lattice_in_scibmad return type"
  )
      .def_readonly("scibmad_file", &Bmad::WriteLatticeInScibmad::scibmad_file)
      .def_readonly("err_flag", &Bmad::WriteLatticeInScibmad::err_flag)
      .def("__len__", [](const Bmad::WriteLatticeInScibmad &) { return 2; })
      .def("__getitem__", [](const Bmad::WriteLatticeInScibmad &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.scibmad_file);
        if (i == 1)
          return py::cast(s.err_flag);
        throw py::index_error();
      });
  m.def(
      "write_lattice_in_scibmad",
      &Bmad::write_lattice_in_scibmad,
      py::arg("lat"),
      R"""(Wrapper for Fortran routine write_lattice_in_scibmad

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
      py::arg("line"),
      py::arg("iu"),
      py::arg("ele"),
      py::arg("lat"),
      R"""(Wrapper for Fortran routine write_line_element

Parameters
----------
line : str

iu : int

ele : EleStruct

lat : LatStruct
)"""
  );
  py::class_<Bmad::WriteOpalFieldGridFile, std::unique_ptr<Bmad::WriteOpalFieldGridFile>>(
      m,
      "WriteOpalFieldGridFile",
      "write_opal_field_grid_file return type"
  )
      .def_readonly("maxfield", &Bmad::WriteOpalFieldGridFile::maxfield)
      .def_readonly("err", &Bmad::WriteOpalFieldGridFile::err)
      .def("__len__", [](const Bmad::WriteOpalFieldGridFile &) { return 2; })
      .def("__getitem__", [](const Bmad::WriteOpalFieldGridFile &s, int i) -> py::object {
        if (i < 0)
          i += 2;
        if (i == 0)
          return py::cast(s.maxfield);
        if (i == 1)
          return py::cast(s.err);
        throw py::index_error();
      });
  m.def(
      "write_opal_field_grid_file",
      &Bmad::write_opal_field_grid_file,
      py::arg("opal_file_unit"),
      py::arg("ele"),
      py::arg("param"),
      R"""(Subroutine write_opal_field_grid_file (opal_file_unit, ele, param, maxfield, err)

Subroutine to write an OPAL lattice file using the information in
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
      py::arg("opal_file_unit"),
      py::arg("lat"),
      R"""(Subroutine write_opal_lattice_file (opal_file_unit, lat, err)

Subroutine to write an OPAL lattice file using the information in
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
      &Bmad::write_time_particle_distribution,
      py::arg("time_file_unit"),
      py::arg("bunch"),
      py::arg("ele"),
      py::arg("style") = py::none(),
      py::arg("branch") = py::none(),
      py::arg("format") = py::none(),
      R"""(Subroutine write_time_particle_distribution  (time_file_unit, bunch, ele, style, branch, format, err)

Subroutine to write a time-based bunch from a standard Bmad bunch

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
