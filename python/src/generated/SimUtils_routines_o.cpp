#include "pybmad/generated/SimUtils_routines_o.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_o(py::module &m) {
  m.def(
      "omega_to_quat",
      &SimUtils::omega_to_quat,
      py::arg("omega"),
      R"""(Function omega_to_quat (omega) result (quat)

Routine to convert from omega + angle representation to a quaternion.

Parameters
----------
omega : float
    Axis of rotation + magnitude = rotation angle.

Returns
-------
quat : float
    Rotation quaternion.
)"""
  );
  m.def(
      "openpmd_species_name",
      &SimUtils::openpmd_species_name,
      py::arg("species"),
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
pmd_name : unknown
    Name of the species. Will return 'INVALID!' (= invalid_name) if index is not valid.
)"""
  );
  m.def(
      "ordinal_str",
      &SimUtils::ordinal_str,
      py::arg("n"),
      py::arg("str"),
      R"""(Wrapper for Fortran routine ordinal_str

Parameters
----------
n : 
str : 
)"""
  );
}
