#include "pybmad/generated/SimUtils_routines_o.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_o(py::module &m) {
  m.def(
      "omega_to_quat",
      &SimUtils::omega_to_quat,
      py::arg("omega"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function omega_to_quat (omega) result (quat)

Routine to convert from omega + angle representation to a quaternion.

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
      py::arg("species"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function openpmd_species_name (species) result(pmd_name)

Routine to return the openPMD name of a particle species given the Bmad species ID.
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
      py::arg("n"),
      py::arg("str"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Wrapper for Fortran routine ordinal_str

Parameters
----------
n : int

str : str
)"""
  );
  m.def(
      "out_io_buffer_get_line",
      &SimUtils::out_io_buffer_get_line,
      py::arg("ix_line"),
      py::arg("line"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function out_io_buffer_get_line(ix_line) result (line)

Routine to return the nuber of lines in the internal buffer.
See the output_direct documentation for more details.
)"""
  );
  m.def(
      "out_io_buffer_num_lines",
      &SimUtils::out_io_buffer_num_lines,
      py::arg("n_lines"),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Function out_io_buffer_num_lines() result (n_lines)

Routine to return the nuber of lines in the internal buffer.
See the output_direct documentation for more details.
)"""
  );
  m.def(
      "out_io_buffer_reset",
      &SimUtils::out_io_buffer_reset,
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine out_io_buffer_reset ()

Routine to initialize the buffer used for capturing output.
)"""
  );
  m.def(
      "out_io",
      py::overload_cast<int, std::string, std::string, int, std::optional<bool>>(&SimUtils::out_io),
      py::arg("level"),
      py::arg("routine_name"),
      py::arg("line"),
      py::arg("i_num"),
      py::arg("insert_tag_line") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
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
  m.def(
      "out_io",
      py::overload_cast<
          int,
          std::string,
          std::string,
          std::optional<std::string>,
          std::optional<std::string>,
          std::optional<std::string>,
          std::optional<std::string>,
          std::optional<std::string>,
          std::optional<std::string>,
          std::optional<std::string>,
          std::optional<std::string>,
          std::optional<std::string>,
          std::optional<std::string>,
          std::optional<std::string>,
          std::optional<FArray1D<Real>>,
          std::optional<FArray1D<Int>>,
          optional_ref<BoolAlloc1D>,
          std::optional<bool>>(&SimUtils::out_io),
      py::arg("level"),
      py::arg("routine_name"),
      py::arg("line1"),
      py::arg("line2") = py::none(),
      py::arg("line3") = py::none(),
      py::arg("line4") = py::none(),
      py::arg("line5") = py::none(),
      py::arg("line6") = py::none(),
      py::arg("line7") = py::none(),
      py::arg("line8") = py::none(),
      py::arg("line9") = py::none(),
      py::arg("line10") = py::none(),
      py::arg("line11") = py::none(),
      py::arg("line12") = py::none(),
      py::arg("r_array") = py::none(),
      py::arg("i_array") = py::none(),
      py::arg("l_array") = py::none(),
      py::arg("insert_tag_line") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
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
      py::overload_cast<
          int,
          std::string,
          CharacterAlloc1D &,
          std::optional<FArray1D<Real>>,
          std::optional<FArray1D<Int>>,
          optional_ref<BoolAlloc1D>,
          std::optional<bool>>(&SimUtils::out_io),
      py::arg("level"),
      py::arg("routine_name"),
      py::arg("lines"),
      py::arg("r_array") = py::none(),
      py::arg("i_array") = py::none(),
      py::arg("l_array") = py::none(),
      py::arg("insert_tag_line") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
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
      py::overload_cast<int, std::string, std::string, bool, std::optional<bool>>(&SimUtils::out_io
      ),
      py::arg("level"),
      py::arg("routine_name"),
      py::arg("line"),
      py::arg("l_num"),
      py::arg("insert_tag_line") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
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
      py::arg("print_on") = py::none(),
      py::arg("capture_state") = py::none(),
      py::arg("capture_add_null") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
      R"""(Subroutine out_io_print_and_capture_setup (print_on, capture_state, capture_add_null)

Set whether a message from a call to out_io is sent to the terminal for printing and/or captured for program use.

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
      py::overload_cast<int, std::string, std::string, double, std::optional<bool>>(
          &SimUtils::out_io
      ),
      py::arg("level"),
      py::arg("routine_name"),
      py::arg("line"),
      py::arg("r_num"),
      py::arg("insert_tag_line") = py::none(),
      py::call_guard<py::gil_scoped_release>(),
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
}
