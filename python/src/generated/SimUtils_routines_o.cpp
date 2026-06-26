#include "pybmad/generated/SimUtils_routines_o.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_o(nb::module_ &m) {
  m.def(
      "omega_to_quat",
      &SimUtils::omega_to_quat,
      nb::arg("omega"),
      R"""(Routine to convert from omega + angle representation to a quaternion.

Parameters
----------
omega : 1D array of float (shape: 3)
    Axis of rotation + magnitude = rotation angle.

Returns
-------
quat : 1D array of float (shape: 0:3)
    Rotation quaternion.
)"""
  );
  m.def(
      "openpmd_species_name",
      &SimUtils::openpmd_species_name,
      nb::arg("species"),
      R"""(Routine to return the openPMD name of a particle species given the Bmad species ID.
Note: the pmd_name does not include the particle charge. For example, if species
corresponds to He+ then the pmd_name will be "He".

Parameters
----------
species : int
    Bmad species ID number.

Returns
-------
pmd_name : str
    Name of the species. Will return 'INVALID!' (= invalid_name) if index is not valid.
)"""
  );
  m.def(
      "ordinal_str",
      &SimUtils::ordinal_str,
      nb::arg("n"),
      R"""(Wrapper for Fortran routine ordinal_str

Parameters
----------
n : int

Returns
-------
str : str
)"""
  );
  m.def(
      "out_io_buffer_get_line",
      &SimUtils::out_io_buffer_get_line,
      nb::arg("ix_line"),
      R"""(Routine to return the nuber of lines in the internal buffer.
See the output_direct documentation for more details.
)"""
  );
  m.def(
      "out_io_buffer_num_lines",
      &SimUtils::out_io_buffer_num_lines,
      R"""(Routine to return the nuber of lines in the internal buffer.
See the output_direct documentation for more details.
)"""
  );
  m.def(
      "out_io_buffer_reset",
      &SimUtils::out_io_buffer_reset,
      R"""(Routine to initialize the buffer used for capturing output.
)"""
  );
  m.def(
      "out_io",
      nb::overload_cast<int, std::string, std::string, int, std::optional<bool>>(&SimUtils::out_io),
      nb::arg("level"),
      nb::arg("routine_name"),
      nb::arg("line"),
      nb::arg("i_num"),
      nb::arg("insert_tag_line") = nb::none(),
      R"""(Wrapper for Fortran routine out_io_int

Parameters
----------
level : int

routine_name : str

line : str

i_num : int

insert_tag_line : bool, optional
)"""
  );
  m
      .def(
          "out_io",
          [](int level,
             std::string routine_name,
             std::string line1,
             std::optional<std::string> line2,
             std::optional<std::string> line3,
             std::optional<std::string> line4,
             std::optional<std::string> line5,
             std::optional<std::string> line6,
             std::optional<std::string> line7,
             std::optional<std::string> line8,
             std::optional<std::string> line9,
             std::optional<std::string> line10,
             std::optional<std::string> line11,
             std::optional<std::string> line12,
             std::optional<FArray1D<Real>> r_array,
             std::optional<FArray1D<Int>> i_array,
             BoolAlloc1D *l_array,
             std::optional<bool> insert_tag_line) {
            auto fn =
                static_cast<void (*)(int, std::string, std::string, std::optional<std::string>, std::optional<std::string>, std::optional<std::string>, std::optional<std::string>, std::optional<std::string>, std::optional<std::string>, std::optional<std::string>, std::optional<std::string>, std::optional<std::string>, std::optional<std::string>, std::optional<std::string>, std::optional<FArray1D<Real>>, std::optional<FArray1D<Int>>, optional_ref<BoolAlloc1D>, std::optional<bool>)>(
                    &SimUtils::out_io
                );
            return fn(
                level,
                routine_name,
                line1,
                line2,
                line3,
                line4,
                line5,
                line6,
                line7,
                line8,
                line9,
                line10,
                line11,
                line12,
                r_array,
                i_array,
                ptr_to_opt_ref(l_array),
                insert_tag_line
            );
          },
          nb::arg("level"),
          nb::arg("routine_name"),
          nb::arg("line1"),
          nb::arg("line2") = nb::none(),
          nb::arg("line3") = nb::none(),
          nb::arg("line4") = nb::none(),
          nb::arg("line5") = nb::none(),
          nb::arg("line6") = nb::none(),
          nb::arg("line7") = nb::none(),
          nb::arg("line8") = nb::none(),
          nb::arg("line9") = nb::none(),
          nb::arg("line10") = nb::none(),
          nb::arg("line11") = nb::none(),
          nb::arg("line12") = nb::none(),
          nb::arg("r_array") = nb::none(),
          nb::arg("i_array") = nb::none(),
          nb::arg("l_array") = nb::none(),
          nb::arg("insert_tag_line") = nb::none(),
          R"""(Wrapper for Fortran routine out_io_line12

Parameters
----------
level : int

routine_name : str

line1 : str

line2 : str, optional

line3 : str, optional

line4 : str, optional

line5 : str, optional

line6 : str, optional

line7 : str, optional

line8 : str, optional

line9 : str, optional

line10 : str, optional

line11 : str, optional

line12 : str, optional

r_array : 1D array of float, optional

i_array : 1D array of int, optional

l_array : 1D array of bool, optional

insert_tag_line : bool, optional
)"""
      );
  m.def(
      "out_io",
      [](int level,
         std::string routine_name,
         CharacterAlloc1D &lines,
         std::optional<FArray1D<Real>> r_array,
         std::optional<FArray1D<Int>> i_array,
         BoolAlloc1D *l_array,
         std::optional<bool> insert_tag_line) {
        auto fn = static_cast<
            void (*)(int, std::string, CharacterAlloc1D &, std::optional<FArray1D<Real>>, std::optional<FArray1D<Int>>, optional_ref<BoolAlloc1D>, std::optional<bool>)>(
            &SimUtils::out_io
        );
        return fn(
            level,
            routine_name,
            lines,
            r_array,
            i_array,
            ptr_to_opt_ref(l_array),
            insert_tag_line
        );
      },
      nb::arg("level"),
      nb::arg("routine_name"),
      nb::arg("lines"),
      nb::arg("r_array") = nb::none(),
      nb::arg("i_array") = nb::none(),
      nb::arg("l_array") = nb::none(),
      nb::arg("insert_tag_line") = nb::none(),
      R"""(Wrapper for Fortran routine out_io_lines

Parameters
----------
level : int

routine_name : str

lines : 1D array of str

r_array : 1D array of float, optional

i_array : 1D array of int, optional

l_array : 1D array of bool, optional

insert_tag_line : bool, optional
)"""
  );
  m.def(
      "out_io",
      nb::overload_cast<int, std::string, std::string, bool, std::optional<bool>>(&SimUtils::out_io
      ),
      nb::arg("level"),
      nb::arg("routine_name"),
      nb::arg("line"),
      nb::arg("l_num"),
      nb::arg("insert_tag_line") = nb::none(),
      R"""(Wrapper for Fortran routine out_io_logical

Parameters
----------
level : int

routine_name : str

line : str

l_num : bool

insert_tag_line : bool, optional
)"""
  );
  m.def(
      "out_io_print_and_capture_setup",
      &SimUtils::out_io_print_and_capture_setup,
      nb::arg("print_on") = nb::none(),
      nb::arg("capture_state") = nb::none(),
      nb::arg("capture_add_null") = nb::none(),
      R"""(Set whether a message from a call to out_io is sent to the terminal for printing and/or captured for program use.

Capture may be desired, for example, to display the output in a separate window or captured output could be passed
to a python process for processing.

The procedure for how a message is handled is as follows:
  First: When out_io is called, the message level is used to determine if anything is to be printed or captured at all.
    When a program is started, everything will pass this test for printing and/or capturing.
    This behavior can be modified by calls to the output_direct routine.
  Second: If a message is to be printed and/or captured (passes the first step), then the internal print_on flag is used
    to determine if printing to the terminal and the internal capture_state flag is used to determine if capture is
    to be done. The initial setting of these flags is print_on = True and capture_state = 'OFF'.
    These internal flags can be set using the print_on and capture_state arguments of this routine.

Notice that whether a message is also written to a file is independent of print and capture logic (see output_direct for more details).

There are two capture modes. buffered (blocked) and unbuffered (unblocked) output.
If a message is to be captured as outlined above, one and only one capture mode is used

Unbuffered output is used when running multithreaded so that the program does not have to wait for output. For example, with a GUI.
With unbuffered output, out_io calls three routines:
  out_io_called(level, routine_name)  ! Called at the start of a message.
  out_io_line(line)                   ! Called for each line of a message.
  out_io_end()                        ! Called at end of a message.
The versions of these routines in the sim_utils library are just dummies.
The idea is that modified versions of these routines can be used to capture the output.

Buffered output uses an internal buffer to store the output.
Output that has been buffered is retrieved by using the routines:
  out_io_buffer_reset
  out_io_buffer_num_lines and
  out_io_buffer_get_line

Parameters
----------
print_on : bool, optional
    If present, set the internal print_on flag to the value of this argument.

capture_state : str, optional
    If present, set the internal capture_state to the value of this argument. Possible values:

capture_add_null : bool, optional
    Is captured output null terminated (for interfacing with C/C++)?
)"""
  );
  m.def(
      "out_io",
      nb::overload_cast<int, std::string, std::string, double, std::optional<bool>>(
          &SimUtils::out_io
      ),
      nb::arg("level"),
      nb::arg("routine_name"),
      nb::arg("line"),
      nb::arg("r_num"),
      nb::arg("insert_tag_line") = nb::none(),
      R"""(Wrapper for Fortran routine out_io_real

Parameters
----------
level : int

routine_name : str

line : str

r_num : float

insert_tag_line : bool, optional
)"""
  );
  m.def(
      "output_direct",
      [](std::optional<int> file_unit,
         std::optional<bool> print_and_capture,
         std::optional<int> min_level,
         std::optional<int> max_level,
         OutIoOutputDirectStruct *set) {
        auto fn = static_cast<
            OutIoOutputDirectStruct (*)(std::optional<int>, std::optional<bool>, std::optional<int>, std::optional<int>, optional_ref<OutIoOutputDirectStruct>)>(
            &SimUtils::output_direct
        );
        return fn(file_unit, print_and_capture, min_level, max_level, ptr_to_opt_ref(set));
      },
      nb::arg("file_unit") = nb::none(),
      nb::arg("print_and_capture") = nb::none(),
      nb::arg("min_level") = nb::none(),
      nb::arg("max_level") = nb::none(),
      nb::arg("set") = nb::none(),
      R"""(Subroutine to set where the output goes when out_io is called.
Output may be sent to the terminal screen, written to a file, and/or captured for program use.

Settings can be made on a message status level by level basis.
See the top of this file for the list of the message status levels.

Once set for a given status level, the settings remain until the next call to
output_direct that cover the same status level.

Parameters
----------
file_unit : int, optional
    Unit number for writing to a file. -1 => No writing (initial default setting).

print_and_capture : bool, optional
    If present then this sets whether output is printed to the terminal and/or captured for program use. Note:
    How output capture works is also set by the out_io_print_and_caputure_setup routine. See the
    out_io_print_and_caputure_setup routine documentation for more details.

min_level : int, optional
    Minimum message status level to apply to. Default is s_blank$

max_level : int, optional
    Maximum message status level to apply to. Default is s_important$

set : OutIoOutputDirectStruct, optional
    If present, use this structure to set where output goes. This structure can be used in place of specifying
    file_unit, etc. One way to use "set" is to first call this routine with the "get" argument to get the
    output direction state.

Returns
-------
get : OutIoOutputDirectStruct, optional
    If present, capture the output direction state before any setting is done.
)"""
  );
}
