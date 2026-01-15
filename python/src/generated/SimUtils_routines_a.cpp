#include "pybmad/generated/SimUtils_routines_a.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Pybmad;

void init_SimUtils_routines_a(py::module &m) {
  m.def(
      "allocate_thread_states",
      &SimUtils::allocate_thread_states,
      R"""(Subroutine allocate_thread_states()

Routine to allocate random number state structures when openMP is used.

)"""
  );
  m.def(
      "anomalous_moment_of",
      &SimUtils::anomalous_moment_of,
      py::arg("species"),
      R"""(Function anomalous_moment_of (species) result (moment)

Routine to return the anomolous moment for subatomic species type. Otherwise returns 0.

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
      py::arg("species"),
      R"""(Function antiparticle (species) result (anti_species)

Routine to return the antiparticle ID given the particle ID.
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
      py::arg("rdata_in"),
      py::arg("bounds"),
      py::arg("window"),
      py::arg("phase"),
      py::arg("diag") = py::none(),
      R"""(subroutine apfft(rdata_in, bounds, window, phase, diag)

Implements the All Phase FFT method for obtaining accurate phase from signal data.

The signal data is truncated to an odd length, and the phase is relative to the central point.

)"""
  );
  py::class_<SimUtils::ApfftCorr, std::unique_ptr<SimUtils::ApfftCorr>>(
      m,
      "ApfftCorr",
      "apfft_corr return type"
  )
      .def_readonly("phase", &SimUtils::ApfftCorr::phase)
      .def_readonly("amp", &SimUtils::ApfftCorr::amp)
      .def_readonly("freq", &SimUtils::ApfftCorr::freq)
      .def("__len__", [](const SimUtils::ApfftCorr &) { return 3; })
      .def("__getitem__", [](const SimUtils::ApfftCorr &s, int i) -> py::object {
        if (i < 0)
          i += 3;
        if (i == 0)
          return py::cast(s.phase);
        if (i == 1)
          return py::cast(s.amp);
        if (i == 2)
          return py::cast(s.freq);
        throw py::index_error();
      });
  m.def(
      "apfft_corr",
      &SimUtils::apfft_corr,
      py::arg("rdata_in"),
      py::arg("bounds") = py::none(),
      py::arg("window"),
      py::arg("diag") = py::none(),
      R"""(subroutine apfft_corr(rdata_in, bounds, window, phase, amp, freq, diag)

For real signal rdata_in, computes phase, frequency, and amplitude
of peak found within bounds.  Algorithm is corrected all-phase FFT and should.

This routine finds only one peak:  the largest amplitude within the bound.  Signals with multiple
components can be investigated by varying bounds appropriately.

Parameters
----------
rdata_in : float
    signal data.
bounds : float
    range within which to search for peak.
window : unknown
    'rec' or 'han' for rectangular or Hann window.
diag : int, optional
    causes low-level routine apfft_ext to produce a fort.X file where X=9000+fid containing diag information.

Returns
-------
phase : float
    phase of peak found in signal.
freq : float
    frequency of peak
amp : float
    amplitude of peak
)"""
  );
  m.def(
      "apfft_ext",
      &SimUtils::apfft_ext,
      py::arg("rdata"),
      py::arg("bounds"),
      py::arg("window"),
      py::arg("phase"),
      py::arg("amp"),
      py::arg("freq"),
      py::arg("diag") = py::none(),
      R"""(subroutine apfft_ext(rdata,bounds, window, phase, amp, freq, diag)

Implements the All Phase FFT method for obtaining accurate phase from signal data.

This "extended" apfft subroutine returns the amplitudes and frequency as well, for use
by the corrected apfft subroutine in this module.

)"""
  );
  m.def(
      "asinc",
      &SimUtils::asinc,
      py::arg("x"),
      py::arg("nd") = py::none(),
      py::arg("y"),
      R"""(Parameters
----------
x : 
nd : 
y : 
)"""
  );
  m.def(
      "assert_equal",
      &SimUtils::assert_equal,
      py::arg("int_arr"),
      py::arg("err_str"),
      py::arg("ival"),
      R"""(Parameters
----------
int_arr : 
err_str : 
ival : 
)"""
  );
  m.def(
      "atomic_number",
      &SimUtils::atomic_number,
      py::arg("species"),
      R"""(Function atomic_number(species) result (atomic_num)

Routine to return the atomic number Z if species argument corresponds to an atomic particle  or is a proton.
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
      py::arg("charge"),
      py::arg("is_anti"),
      py::arg("atomic_num"),
      py::arg("n_nuc"),
      R"""(Function atomic_species_id(charge, is_anti, atomic_num, n_nuc) result (species_id)

Routine to return the species ID for an atom

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
      py::arg("axis"),
      py::arg("angle"),
      R"""(Function axis_angle_to_quat (axis, angle) result (quat)

Routine to convert from axis + angle representation to a quaternion.

Parameters
----------
axis : float
    Axis of rotation.
angle : float
    angle of rotation.

Returns
-------
quat : float
    Rotation quaternion.
)"""
  );
  m.def(
      "axis_angle_to_w_mat",
      &SimUtils::axis_angle_to_w_mat,
      py::arg("axis"),
      py::arg("angle"),
      R"""(Subroutine axis_angle_to_w_mat (axis, angle, w_mat)

Routine to construct the 3D rotation matrix w_mat given an axis of rotation
and a rotation angle.

Parameters
----------
axis : float
    Rotation axis. Does not have to be normalized.
angle : float
    Rotation angle in the range [-pi, pi].

Returns
-------
w_mat : float
    Rotation matrix
)"""
  );
}
