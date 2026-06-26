#include "pybmad/generated/SimUtils_routines_a.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_SimUtils_routines_a(nb::module_ &m) {
  m.def(
      "all_pointer_to_string",
      &SimUtils::all_pointer_to_string,
      nb::arg("a_ptr"),
      nb::arg("err") = nb::none(),
      R"""(Wrapper for Fortran routine all_pointer_to_string

Parameters
----------
a_ptr : AllPointerStruct

err : bool, optional

Returns
-------
str : str
)"""
  );
  m.def(
      "allocate_thread_states",
      &SimUtils::allocate_thread_states,
      R"""(Routine to allocate random number state structures when openMP is used.
)"""
  );
  m.def(
      "anomalous_moment_of",
      &SimUtils::anomalous_moment_of,
      nb::arg("species"),
      R"""(Routine to return the anomolous moment for subatomic species type. Otherwise returns 0.

Parameters
----------
species : int
    Species ID.

Returns
-------
moment : float
    Anomalous moment.
)"""
  );
  m.def(
      "antiparticle",
      &SimUtils::antiparticle,
      nb::arg("species"),
      R"""(Routine to return the antiparticle ID given the particle ID.
For a molecule the anti-species is just the molecude with the charge reversed.

Parameters
----------
species : int
    Particle ID.

Returns
-------
anti_species : int
    Antiparticle ID.
)"""
  );
  m.def(
      "apfft",
      &SimUtils::apfft,
      nb::arg("rdata_in"),
      nb::arg("bounds"),
      nb::arg("window"),
      nb::arg("phase"),
      nb::arg("diag") = nb::none(),
      R"""(Implements the All Phase FFT method for obtaining accurate phase from signal data.

The signal data is truncated to an odd length, and the phase is relative to the central point.
)"""
  );
  nb::class_<SimUtils::ApfftCorr>(m, "ApfftCorr", "apfft_corr return type")
      .def_ro("phase", &SimUtils::ApfftCorr::phase)
      .def_ro("amp", &SimUtils::ApfftCorr::amp)
      .def_ro("freq", &SimUtils::ApfftCorr::freq)
      .def("__len__", [](const SimUtils::ApfftCorr &) { return 3; })
      .def("__getitem__", [](const SimUtils::ApfftCorr &s, int i) -> nb::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return nb::cast(s.phase);
        if (i == 1)
          return nb::cast(s.amp);
        if (i == 2)
          return nb::cast(s.freq);
        throw nb::index_error();
      });
  m.def(
      "apfft_corr",
      &SimUtils::apfft_corr,
      nb::arg("rdata_in"),
      nb::arg("window"),
      nb::arg("bounds") = nb::none(),
      nb::arg("diag") = nb::none(),
      R"""(For real signal rdata_in, computes phase, frequency, and amplitude
of peak found within bounds.  Algorithm is corrected all-phase FFT and should.

This routine finds only one peak:  the largest amplitude within the bound.  Signals with multiple
components can be investigated by varying bounds appropriately.

Parameters
----------
rdata_in : 1D array of float
    signal data.

window : str
    'rec' or 'han' for rectangular or Hann window.

bounds : 1D array of float (shape: 2), optional
    range within which to search for peak.

diag : int, optional
    causes low-level routine apfft_ext to produce a fort.X file where X=9000+fid containing diag information.

Returns
-------
phase : float
    phase of peak found in signal.

amp : float
    amplitude of peak

freq : float
    frequency of peak
)"""
  );
  m.def(
      "apfft_ext",
      &SimUtils::apfft_ext,
      nb::arg("rdata"),
      nb::arg("bounds"),
      nb::arg("window"),
      nb::arg("phase"),
      nb::arg("amp"),
      nb::arg("freq"),
      nb::arg("diag") = nb::none(),
      R"""(Implements the All Phase FFT method for obtaining accurate phase from signal data.

This "extended" apfft subroutine returns the amplitudes and frequency as well, for use
by the corrected apfft subroutine in this module.
)"""
  );
  m.def(
      "asinc",
      &SimUtils::asinc,
      nb::arg("x"),
      nb::arg("nd") = nb::none(),
      R"""(Wrapper for Fortran routine asinc

Parameters
----------
x : float

nd : int, optional
    Derivative order. nd = 0 (default) -> compute arcsin(x) / x NOTE: Currently only nd = 0 and nd = 1 are
    implemented.

Returns
-------
y : float
    nd^th derivative. of arcsin(x)/x
)"""
  );
  m.def(
      "assert_equal",
      &SimUtils::assert_equal,
      nb::arg("int_arr"),
      nb::arg("err_str"),
      R"""(Wrapper for Fortran routine assert_equal

Parameters
----------
int_arr : 1D array of int

err_str : str

Returns
-------
ival : int
)"""
  );
  m.def(
      "atomic_number",
      &SimUtils::atomic_number,
      nb::arg("species"),
      R"""(Routine to return the atomic number Z if species argument corresponds to an atomic particle  or is a proton.
Set to the charge for atomic particles.
Set to zero for molecules.

Parameters
----------
species : int
    Spicies ID.

Returns
-------
atomic_num : int
    Atomic index. Set to zero if a molecule
)"""
  );
  m.def(
      "atomic_species_id",
      &SimUtils::atomic_species_id,
      nb::arg("charge"),
      nb::arg("is_anti"),
      nb::arg("atomic_num"),
      nb::arg("n_nuc"),
      R"""(Routine to return the species ID for an atom

Parameters
----------
charge : int
    Charge of the atom.

is_anti : bool
    Is an anti-atom.

atomic_num : int
    Atomic number.

n_nuc : int
    Number of nucleons.

Returns
-------
species_id : int
    Species ID number.
)"""
  );
  m.def(
      "axis_angle_to_quat",
      &SimUtils::axis_angle_to_quat,
      nb::arg("axis"),
      nb::arg("angle"),
      R"""(Routine to convert from axis + angle representation to a quaternion.

Parameters
----------
axis : 1D array of float (shape: 3)
    Axis of rotation.

angle : float
    angle of rotation.

Returns
-------
quat : 1D array of float (shape: 0:3)
    Rotation quaternion.
)"""
  );
  m.def(
      "axis_angle_to_w_mat",
      &SimUtils::axis_angle_to_w_mat,
      nb::arg("axis"),
      nb::arg("angle"),
      R"""(Routine to construct the 3D rotation matrix w_mat given an axis of rotation
and a rotation angle.

Parameters
----------
axis : 1D array of float (shape: 3)
    Rotation axis. Does not have to be normalized.

angle : float
    Rotation angle in the range [-pi, pi].

Returns
-------
w_mat : 2D array of float (shape: 3,3)
    Rotation matrix
)"""
  );
}
