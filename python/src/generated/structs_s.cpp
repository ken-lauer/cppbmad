#include "pybmad/generated/structs_s.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// space_charge_common_struct
void init_space_charge_common_struct(
    py::module& m,
    py::class_<SpaceChargeCommonProxy>& cls) {
  cls.def(py::init<>())
      // SpaceChargeCommonProxy.ds_track_step (0D_NOT_real - CSR tracking step size
      .def_property(
          "ds_track_step",
          &SpaceChargeCommonProxy::ds_track_step,
          &SpaceChargeCommonProxy::set_ds_track_step)
      // SpaceChargeCommonProxy.dt_track_step (0D_NOT_real - Time Runge kutta initial step.
      .def_property(
          "dt_track_step",
          &SpaceChargeCommonProxy::dt_track_step,
          &SpaceChargeCommonProxy::set_dt_track_step)
      // SpaceChargeCommonProxy.cathode_strength_cutoff (0D_NOT_real - Cutoff for the cathode field calc.
      .def_property(
          "cathode_strength_cutoff",
          &SpaceChargeCommonProxy::cathode_strength_cutoff,
          &SpaceChargeCommonProxy::set_cathode_strength_cutoff)
      // SpaceChargeCommonProxy.rel_tol_tracking (0D_NOT_real - Relative tolerance for tracking.
      .def_property(
          "rel_tol_tracking",
          &SpaceChargeCommonProxy::rel_tol_tracking,
          &SpaceChargeCommonProxy::set_rel_tol_tracking)
      // SpaceChargeCommonProxy.abs_tol_tracking (0D_NOT_real - Absolute tolerance for tracking.
      .def_property(
          "abs_tol_tracking",
          &SpaceChargeCommonProxy::abs_tol_tracking,
          &SpaceChargeCommonProxy::set_abs_tol_tracking)
      // SpaceChargeCommonProxy.beam_chamber_height (0D_NOT_real - Used in shielding calculation.
      .def_property(
          "beam_chamber_height",
          &SpaceChargeCommonProxy::beam_chamber_height,
          &SpaceChargeCommonProxy::set_beam_chamber_height)
      // SpaceChargeCommonProxy.lsc_sigma_cutoff (0D_NOT_real - Cutoff for the 1-dim longitudinal SC calc. If a bin sigma is < cutoff * sigma_ave then ignore.
      .def_property(
          "lsc_sigma_cutoff",
          &SpaceChargeCommonProxy::lsc_sigma_cutoff,
          &SpaceChargeCommonProxy::set_lsc_sigma_cutoff)
      // SpaceChargeCommonProxy.particle_sigma_cutoff (0D_NOT_real - 3D SC calc cutoff for particles with (x,y,z) position far from the center. Negative or zero means ignore.
      .def_property(
          "particle_sigma_cutoff",
          &SpaceChargeCommonProxy::particle_sigma_cutoff,
          &SpaceChargeCommonProxy::set_particle_sigma_cutoff)
      // SpaceChargeCommonProxy.space_charge_mesh_size (1D_NOT_integer - Gird size for fft_3d space charge calc.
      .def_property_readonly(
          "space_charge_mesh_size",
          &SpaceChargeCommonProxy::space_charge_mesh_size)
      // SpaceChargeCommonProxy.csr3d_mesh_size (1D_NOT_integer - Gird size for CSR.
      .def_property_readonly(
          "csr3d_mesh_size", &SpaceChargeCommonProxy::csr3d_mesh_size)
      // SpaceChargeCommonProxy.n_bin (0D_NOT_integer - Number of bins used
      .def_property(
          "n_bin",
          &SpaceChargeCommonProxy::n_bin,
          &SpaceChargeCommonProxy::set_n_bin)
      // SpaceChargeCommonProxy.particle_bin_span (0D_NOT_integer - Longitudinal particle length / dz_bin
      .def_property(
          "particle_bin_span",
          &SpaceChargeCommonProxy::particle_bin_span,
          &SpaceChargeCommonProxy::set_particle_bin_span)
      // SpaceChargeCommonProxy.n_shield_images (0D_NOT_integer - Chamber wall shielding. 0 = no shielding.
      .def_property(
          "n_shield_images",
          &SpaceChargeCommonProxy::n_shield_images,
          &SpaceChargeCommonProxy::set_n_shield_images)
      // SpaceChargeCommonProxy.sc_min_in_bin (0D_NOT_integer - Minimum number of particles in a bin for sigmas to be valid.
      .def_property(
          "sc_min_in_bin",
          &SpaceChargeCommonProxy::sc_min_in_bin,
          &SpaceChargeCommonProxy::set_sc_min_in_bin)
      // SpaceChargeCommonProxy.lsc_kick_transverse_dependence (0D_NOT_logical -
      .def_property(
          "lsc_kick_transverse_dependence",
          &SpaceChargeCommonProxy::lsc_kick_transverse_dependence,
          &SpaceChargeCommonProxy::set_lsc_kick_transverse_dependence)
      // SpaceChargeCommonProxy.debug (0D_NOT_logical -
      .def_property(
          "debug",
          &SpaceChargeCommonProxy::debug,
          &SpaceChargeCommonProxy::set_debug)
      // SpaceChargeCommonProxy.diagnostic_output_file (0D_NOT_character - If non-blank write a diagnostic (EG wake) file
      .def_property(
          "diagnostic_output_file",
          &SpaceChargeCommonProxy::diagnostic_output_file,
          &SpaceChargeCommonProxy::set_diagnostic_output_file)

      .def(
          "__repr__",
          [](const SpaceChargeCommonProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SpaceChargeCommonProxy& self) {
            return SpaceChargeCommonProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SpaceChargeCommonProxy& self, py::dict& memo) {
            return SpaceChargeCommonProxy(self);
          })

      ;

  // 1D SpaceChargeCommonProxy arrays are not used in structs/routines
  // 2D SpaceChargeCommonProxy arrays are not used in structs/routines
  // 3D SpaceChargeCommonProxy arrays are not used in structs/routines
}

// =============================================================================
// spin_axis_struct
void init_spin_axis_struct(py::module& m, py::class_<SpinAxisProxy>& cls) {
  cls.def(py::init<>())
      // SpinAxisProxy.l (1D_NOT_real - Transverse axis.
      .def_property_readonly("l", &SpinAxisProxy::l)
      // SpinAxisProxy.n0 (1D_NOT_real - Invariant spin axis on closed orbit.
      .def_property_readonly("n0", &SpinAxisProxy::n0)
      // SpinAxisProxy.m (1D_NOT_real - Transverse axis.
      .def_property_readonly("m", &SpinAxisProxy::m)

      .def(
          "__repr__", [](const SpinAxisProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SpinAxisProxy& self) {
            return SpinAxisProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SpinAxisProxy& self, py::dict& memo) {
            return SpinAxisProxy(self);
          })

      ;

  // 1D SpinAxisProxy arrays are not used in structs/routines
  // 2D SpinAxisProxy arrays are not used in structs/routines
  // 3D SpinAxisProxy arrays are not used in structs/routines
}

// =============================================================================
// spin_orbit_map1_struct
void init_spin_orbit_map1_struct(
    py::module& m,
    py::class_<SpinOrbitMap1Proxy>& cls) {
  cls.def(py::init<>())
      // SpinOrbitMap1Proxy.orb_mat (2D_NOT_real - Orbital matrix
      .def_property_readonly("orb_mat", &SpinOrbitMap1Proxy::orb_mat)
      // SpinOrbitMap1Proxy.vec0 (1D_NOT_real - Orbital 0th order map: r_out = mat6 * r_in + vec0
      .def_property_readonly("vec0", &SpinOrbitMap1Proxy::vec0)
      // SpinOrbitMap1Proxy.spin_q (2D_NOT_real - 0th and 1st order quaternion spin map
      .def_property_readonly("spin_q", &SpinOrbitMap1Proxy::spin_q)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return SpinOrbitMap1ProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const SpinOrbitMap1Proxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SpinOrbitMap1Proxy& self) {
            return SpinOrbitMap1Proxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SpinOrbitMap1Proxy& self, py::dict& memo) {
            return SpinOrbitMap1Proxy(self);
          })

      ;

  bind_FTypeArrayND<SpinOrbitMap1ProxyArray1D>(m, "SpinOrbitMap1StructArray1D");
  bind_FTypeAlloc1D<SpinOrbitMap1ProxyAlloc1D>(m, "SpinOrbitMap1StructAlloc1D");
  // 2D SpinOrbitMap1Proxy arrays are not used in structs/routines
  // 3D SpinOrbitMap1Proxy arrays are not used in structs/routines
}

// =============================================================================
// spin_polar_struct
void init_spin_polar_struct(py::module& m, py::class_<SpinPolarProxy>& cls) {
  cls.def(py::init<>())
      // SpinPolarProxy.polarization (0D_NOT_real -
      .def_property(
          "polarization",
          &SpinPolarProxy::polarization,
          &SpinPolarProxy::set_polarization)
      // SpinPolarProxy.theta (0D_NOT_real - Spherical coords: Angle from z-axis.
      .def_property("theta", &SpinPolarProxy::theta, &SpinPolarProxy::set_theta)
      // SpinPolarProxy.phi (0D_NOT_real - Spherical coords: Angle in (x,y) plane.
      .def_property("phi", &SpinPolarProxy::phi, &SpinPolarProxy::set_phi)
      // SpinPolarProxy.xi (0D_NOT_real - Spinor phase angle (See Bmad manual).
      .def_property("xi", &SpinPolarProxy::xi, &SpinPolarProxy::set_xi)

      .def(
          "__repr__",
          [](const SpinPolarProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SpinPolarProxy& self) {
            return SpinPolarProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SpinPolarProxy& self, py::dict& memo) {
            return SpinPolarProxy(self);
          })

      ;

  // 1D SpinPolarProxy arrays are not used in structs/routines
  // 2D SpinPolarProxy arrays are not used in structs/routines
  // 3D SpinPolarProxy arrays are not used in structs/routines
}

// =============================================================================
// strong_beam_struct
void init_strong_beam_struct(py::module& m, py::class_<StrongBeamProxy>& cls) {
  cls.def(py::init<>())
      // StrongBeamProxy.ix_slice (0D_NOT_integer - 0 -> at element center and not at slice.
      .def_property(
          "ix_slice",
          &StrongBeamProxy::ix_slice,
          &StrongBeamProxy::set_ix_slice)
      // StrongBeamProxy.x_center (0D_NOT_real - Strong beam slice center.
      .def_property(
          "x_center",
          &StrongBeamProxy::x_center,
          &StrongBeamProxy::set_x_center)
      // StrongBeamProxy.y_center (0D_NOT_real - Strong beam slice center.
      .def_property(
          "y_center",
          &StrongBeamProxy::y_center,
          &StrongBeamProxy::set_y_center)
      // StrongBeamProxy.x_sigma (0D_NOT_real - Strong beam slice sigma.
      .def_property(
          "x_sigma", &StrongBeamProxy::x_sigma, &StrongBeamProxy::set_x_sigma)
      // StrongBeamProxy.y_sigma (0D_NOT_real - Strong beam slice sigma.
      .def_property(
          "y_sigma", &StrongBeamProxy::y_sigma, &StrongBeamProxy::set_y_sigma)
      // StrongBeamProxy.dx (0D_NOT_real - Particle - beam slice distance.
      .def_property("dx", &StrongBeamProxy::dx, &StrongBeamProxy::set_dx)
      // StrongBeamProxy.dy (0D_NOT_real - Particle - beam slice distance.
      .def_property("dy", &StrongBeamProxy::dy, &StrongBeamProxy::set_dy)

      .def(
          "__repr__",
          [](const StrongBeamProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const StrongBeamProxy& self) {
            return StrongBeamProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const StrongBeamProxy& self, py::dict& memo) {
            return StrongBeamProxy(self);
          })

      ;

  // 1D StrongBeamProxy arrays are not used in structs/routines
  // 2D StrongBeamProxy arrays are not used in structs/routines
  // 3D StrongBeamProxy arrays are not used in structs/routines
}

// =============================================================================
// surface_curvature_struct
void init_surface_curvature_struct(
    py::module& m,
    py::class_<SurfaceCurvatureProxy>& cls) {
  cls.def(py::init<>())
      // SurfaceCurvatureProxy.xy (2D_NOT_real -
      .def_property_readonly("xy", &SurfaceCurvatureProxy::xy)
      // SurfaceCurvatureProxy.spherical (0D_NOT_real -
      .def_property(
          "spherical",
          &SurfaceCurvatureProxy::spherical,
          &SurfaceCurvatureProxy::set_spherical)
      // SurfaceCurvatureProxy.elliptical (1D_NOT_real - Total curvature = elliptical + spherical
      .def_property_readonly("elliptical", &SurfaceCurvatureProxy::elliptical)
      // SurfaceCurvatureProxy.has_curvature (0D_NOT_logical - Dependent var. Will be set by Bmad
      .def_property(
          "has_curvature",
          &SurfaceCurvatureProxy::has_curvature,
          &SurfaceCurvatureProxy::set_has_curvature)

      .def(
          "__repr__",
          [](const SurfaceCurvatureProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceCurvatureProxy& self) {
            return SurfaceCurvatureProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SurfaceCurvatureProxy& self, py::dict& memo) {
            return SurfaceCurvatureProxy(self);
          })

      ;

  // 1D SurfaceCurvatureProxy arrays are not used in structs/routines
  // 2D SurfaceCurvatureProxy arrays are not used in structs/routines
  // 3D SurfaceCurvatureProxy arrays are not used in structs/routines
}

// =============================================================================
// surface_displacement_pt_struct
void init_surface_displacement_pt_struct(
    py::module& m,
    py::class_<SurfaceDisplacementPtProxy>& cls) {
  cls.def(py::init<>())
      // SurfaceDisplacementPtProxy.x0 (0D_NOT_real - Position at center
      .def_property(
          "x0",
          &SurfaceDisplacementPtProxy::x0,
          &SurfaceDisplacementPtProxy::set_x0)
      // SurfaceDisplacementPtProxy.y0 (0D_NOT_real - Position at center
      .def_property(
          "y0",
          &SurfaceDisplacementPtProxy::y0,
          &SurfaceDisplacementPtProxy::set_y0)
      // SurfaceDisplacementPtProxy.z0 (0D_NOT_real -
      .def_property(
          "z0",
          &SurfaceDisplacementPtProxy::z0,
          &SurfaceDisplacementPtProxy::set_z0)
      // SurfaceDisplacementPtProxy.dz_dx (0D_NOT_real -
      .def_property(
          "dz_dx",
          &SurfaceDisplacementPtProxy::dz_dx,
          &SurfaceDisplacementPtProxy::set_dz_dx)
      // SurfaceDisplacementPtProxy.dz_dy (0D_NOT_real -
      .def_property(
          "dz_dy",
          &SurfaceDisplacementPtProxy::dz_dy,
          &SurfaceDisplacementPtProxy::set_dz_dy)
      // SurfaceDisplacementPtProxy.d2z_dxdy (0D_NOT_real -
      .def_property(
          "d2z_dxdy",
          &SurfaceDisplacementPtProxy::d2z_dxdy,
          &SurfaceDisplacementPtProxy::set_d2z_dxdy)

      .def(
          "__repr__",
          [](const SurfaceDisplacementPtProxy& self) {
            return to_string(self);
          })

      .def(
          "__copy__",
          [](const SurfaceDisplacementPtProxy& self) {
            return SurfaceDisplacementPtProxy(
                self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SurfaceDisplacementPtProxy& self, py::dict& memo) {
            return SurfaceDisplacementPtProxy(self);
          })

      ;

  // 1D SurfaceDisplacementPtProxy arrays are not used in structs/routines
  bind_FTypeArrayND<SurfaceDisplacementPtProxyArray2D>(
      m, "SurfaceDisplacementPtStructArray2D");
  // 3D SurfaceDisplacementPtProxy arrays are not used in structs/routines
}

// =============================================================================
// surface_displacement_struct
void init_surface_displacement_struct(
    py::module& m,
    py::class_<SurfaceDisplacementProxy>& cls) {
  cls.def(py::init<>())
      // SurfaceDisplacementProxy.active (0D_NOT_logical -
      .def_property(
          "active",
          &SurfaceDisplacementProxy::active,
          &SurfaceDisplacementProxy::set_active)
      // SurfaceDisplacementProxy.dr (1D_NOT_real -
      .def_property_readonly("dr", &SurfaceDisplacementProxy::dr)
      // SurfaceDisplacementProxy.r0 (1D_NOT_real -
      .def_property_readonly("r0", &SurfaceDisplacementProxy::r0)
      // SurfaceDisplacementProxy.pt (2D_ALLOC_type -
      .def_property_readonly("pt", &SurfaceDisplacementProxy::pt)

      .def(
          "__repr__",
          [](const SurfaceDisplacementProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceDisplacementProxy& self) {
            return SurfaceDisplacementProxy(
                self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SurfaceDisplacementProxy& self, py::dict& memo) {
            return SurfaceDisplacementProxy(self);
          })

      ;

  // 1D SurfaceDisplacementProxy arrays are not used in structs/routines
  // 2D SurfaceDisplacementProxy arrays are not used in structs/routines
  // 3D SurfaceDisplacementProxy arrays are not used in structs/routines
}

// =============================================================================
// surface_h_misalign_pt_struct
void init_surface_h_misalign_pt_struct(
    py::module& m,
    py::class_<SurfaceHMisalignPtProxy>& cls) {
  cls.def(py::init<>())
      // SurfaceHMisalignPtProxy.x0 (0D_NOT_real - Position at center
      .def_property(
          "x0", &SurfaceHMisalignPtProxy::x0, &SurfaceHMisalignPtProxy::set_x0)
      // SurfaceHMisalignPtProxy.y0 (0D_NOT_real - Position at center
      .def_property(
          "y0", &SurfaceHMisalignPtProxy::y0, &SurfaceHMisalignPtProxy::set_y0)
      // SurfaceHMisalignPtProxy.rot_y (0D_NOT_real - rot_t = x-rotation for Bragg and z-rotation for Laue.
      .def_property(
          "rot_y",
          &SurfaceHMisalignPtProxy::rot_y,
          &SurfaceHMisalignPtProxy::set_rot_y)
      // SurfaceHMisalignPtProxy.rot_t (0D_NOT_real - rot_t = x-rotation for Bragg and z-rotation for Laue.
      .def_property(
          "rot_t",
          &SurfaceHMisalignPtProxy::rot_t,
          &SurfaceHMisalignPtProxy::set_rot_t)
      // SurfaceHMisalignPtProxy.rot_y_rms (0D_NOT_real - rot_t = x-rotation for Bragg and z-rotation for Laue.
      .def_property(
          "rot_y_rms",
          &SurfaceHMisalignPtProxy::rot_y_rms,
          &SurfaceHMisalignPtProxy::set_rot_y_rms)
      // SurfaceHMisalignPtProxy.rot_t_rms (0D_NOT_real - rot_t = x-rotation for Bragg and z-rotation for Laue.
      .def_property(
          "rot_t_rms",
          &SurfaceHMisalignPtProxy::rot_t_rms,
          &SurfaceHMisalignPtProxy::set_rot_t_rms)

      .def(
          "__repr__",
          [](const SurfaceHMisalignPtProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceHMisalignPtProxy& self) {
            return SurfaceHMisalignPtProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SurfaceHMisalignPtProxy& self, py::dict& memo) {
            return SurfaceHMisalignPtProxy(self);
          })

      ;

  // 1D SurfaceHMisalignPtProxy arrays are not used in structs/routines
  bind_FTypeArrayND<SurfaceHMisalignPtProxyArray2D>(
      m, "SurfaceHMisalignPtStructArray2D");
  // 3D SurfaceHMisalignPtProxy arrays are not used in structs/routines
}

// =============================================================================
// surface_h_misalign_struct
void init_surface_h_misalign_struct(
    py::module& m,
    py::class_<SurfaceHMisalignProxy>& cls) {
  cls.def(py::init<>())
      // SurfaceHMisalignProxy.active (0D_NOT_logical -
      .def_property(
          "active",
          &SurfaceHMisalignProxy::active,
          &SurfaceHMisalignProxy::set_active)
      // SurfaceHMisalignProxy.dr (1D_NOT_real -
      .def_property_readonly("dr", &SurfaceHMisalignProxy::dr)
      // SurfaceHMisalignProxy.r0 (1D_NOT_real -
      .def_property_readonly("r0", &SurfaceHMisalignProxy::r0)
      // SurfaceHMisalignProxy.pt (2D_ALLOC_type -
      .def_property_readonly("pt", &SurfaceHMisalignProxy::pt)

      .def(
          "__repr__",
          [](const SurfaceHMisalignProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceHMisalignProxy& self) {
            return SurfaceHMisalignProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SurfaceHMisalignProxy& self, py::dict& memo) {
            return SurfaceHMisalignProxy(self);
          })

      ;

  // 1D SurfaceHMisalignProxy arrays are not used in structs/routines
  // 2D SurfaceHMisalignProxy arrays are not used in structs/routines
  // 3D SurfaceHMisalignProxy arrays are not used in structs/routines
}

// =============================================================================
// surface_segmented_pt_struct
void init_surface_segmented_pt_struct(
    py::module& m,
    py::class_<SurfaceSegmentedPtProxy>& cls) {
  cls.def(py::init<>())
      // SurfaceSegmentedPtProxy.x0 (0D_NOT_real - Position at center
      .def_property(
          "x0", &SurfaceSegmentedPtProxy::x0, &SurfaceSegmentedPtProxy::set_x0)
      // SurfaceSegmentedPtProxy.y0 (0D_NOT_real - Position at center
      .def_property(
          "y0", &SurfaceSegmentedPtProxy::y0, &SurfaceSegmentedPtProxy::set_y0)
      // SurfaceSegmentedPtProxy.z0 (0D_NOT_real - Position at center
      .def_property(
          "z0", &SurfaceSegmentedPtProxy::z0, &SurfaceSegmentedPtProxy::set_z0)
      // SurfaceSegmentedPtProxy.dz_dx (0D_NOT_real - Slope at center
      .def_property(
          "dz_dx",
          &SurfaceSegmentedPtProxy::dz_dx,
          &SurfaceSegmentedPtProxy::set_dz_dx)
      // SurfaceSegmentedPtProxy.dz_dy (0D_NOT_real - Slope at center
      .def_property(
          "dz_dy",
          &SurfaceSegmentedPtProxy::dz_dy,
          &SurfaceSegmentedPtProxy::set_dz_dy)

      .def(
          "__repr__",
          [](const SurfaceSegmentedPtProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceSegmentedPtProxy& self) {
            return SurfaceSegmentedPtProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SurfaceSegmentedPtProxy& self, py::dict& memo) {
            return SurfaceSegmentedPtProxy(self);
          })

      ;

  // 1D SurfaceSegmentedPtProxy arrays are not used in structs/routines
  bind_FTypeArrayND<SurfaceSegmentedPtProxyArray2D>(
      m, "SurfaceSegmentedPtStructArray2D");
  // 3D SurfaceSegmentedPtProxy arrays are not used in structs/routines
}

// =============================================================================
// surface_segmented_struct
void init_surface_segmented_struct(
    py::module& m,
    py::class_<SurfaceSegmentedProxy>& cls) {
  cls.def(py::init<>())
      // SurfaceSegmentedProxy.active (0D_NOT_logical -
      .def_property(
          "active",
          &SurfaceSegmentedProxy::active,
          &SurfaceSegmentedProxy::set_active)
      // SurfaceSegmentedProxy.dr (1D_NOT_real -
      .def_property_readonly("dr", &SurfaceSegmentedProxy::dr)
      // SurfaceSegmentedProxy.r0 (1D_NOT_real -
      .def_property_readonly("r0", &SurfaceSegmentedProxy::r0)
      // SurfaceSegmentedProxy.pt (2D_ALLOC_type -
      .def_property_readonly("pt", &SurfaceSegmentedProxy::pt)

      .def(
          "__repr__",
          [](const SurfaceSegmentedProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceSegmentedProxy& self) {
            return SurfaceSegmentedProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SurfaceSegmentedProxy& self, py::dict& memo) {
            return SurfaceSegmentedProxy(self);
          })

      ;

  // 1D SurfaceSegmentedProxy arrays are not used in structs/routines
  // 2D SurfaceSegmentedProxy arrays are not used in structs/routines
  // 3D SurfaceSegmentedProxy arrays are not used in structs/routines
}

// =============================================================================
// spline_struct
void init_spline_struct(py::module& m, py::class_<SplineProxy>& cls) {
  cls.def(py::init<>())
      // SplineProxy.x0 (0D_NOT_real - Point at start of spline
      .def_property("x0", &SplineProxy::x0, &SplineProxy::set_x0)
      // SplineProxy.y0 (0D_NOT_real - Point at start of spline
      .def_property("y0", &SplineProxy::y0, &SplineProxy::set_y0)
      // SplineProxy.x1 (0D_NOT_real - Point at end of spline
      .def_property("x1", &SplineProxy::x1, &SplineProxy::set_x1)
      // SplineProxy.coef (1D_NOT_real - coefficients for cubic spline
      .def_property_readonly("coef", &SplineProxy::coef)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return SplineProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const SplineProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SplineProxy& self) {
            return SplineProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SplineProxy& self, py::dict& memo) {
            return SplineProxy(self);
          })

      ;

  bind_FTypeArrayND<SplineProxyArray1D>(m, "SplineStructArray1D");
  bind_FTypeAlloc1D<SplineProxyAlloc1D>(m, "SplineStructAlloc1D");
  // 2D SplineProxy arrays are not used in structs/routines
  // 3D SplineProxy arrays are not used in structs/routines
}

// =============================================================================
// summation_rdt_struct
void init_summation_rdt_struct(
    py::module& m,
    py::class_<SummationRdtProxy>& cls) {
  cls.def(py::init<>())
      // SummationRdtProxy.h11001 (0D_NOT_complex -
      .def_property(
          "h11001", &SummationRdtProxy::h11001, &SummationRdtProxy::set_h11001)
      // SummationRdtProxy.h00111 (0D_NOT_complex -
      .def_property(
          "h00111", &SummationRdtProxy::h00111, &SummationRdtProxy::set_h00111)
      // SummationRdtProxy.h20001 (0D_NOT_complex -
      .def_property(
          "h20001", &SummationRdtProxy::h20001, &SummationRdtProxy::set_h20001)
      // SummationRdtProxy.h00201 (0D_NOT_complex -
      .def_property(
          "h00201", &SummationRdtProxy::h00201, &SummationRdtProxy::set_h00201)
      // SummationRdtProxy.h10002 (0D_NOT_complex -
      .def_property(
          "h10002", &SummationRdtProxy::h10002, &SummationRdtProxy::set_h10002)
      // SummationRdtProxy.h21000 (0D_NOT_complex -
      .def_property(
          "h21000", &SummationRdtProxy::h21000, &SummationRdtProxy::set_h21000)
      // SummationRdtProxy.h30000 (0D_NOT_complex -
      .def_property(
          "h30000", &SummationRdtProxy::h30000, &SummationRdtProxy::set_h30000)
      // SummationRdtProxy.h10110 (0D_NOT_complex -
      .def_property(
          "h10110", &SummationRdtProxy::h10110, &SummationRdtProxy::set_h10110)
      // SummationRdtProxy.h10020 (0D_NOT_complex -
      .def_property(
          "h10020", &SummationRdtProxy::h10020, &SummationRdtProxy::set_h10020)
      // SummationRdtProxy.h10200 (0D_NOT_complex - 2nd order in K2 moments
      .def_property(
          "h10200", &SummationRdtProxy::h10200, &SummationRdtProxy::set_h10200)
      // SummationRdtProxy.h31000 (0D_NOT_complex -
      .def_property(
          "h31000", &SummationRdtProxy::h31000, &SummationRdtProxy::set_h31000)
      // SummationRdtProxy.h40000 (0D_NOT_complex -
      .def_property(
          "h40000", &SummationRdtProxy::h40000, &SummationRdtProxy::set_h40000)
      // SummationRdtProxy.h20110 (0D_NOT_complex -
      .def_property(
          "h20110", &SummationRdtProxy::h20110, &SummationRdtProxy::set_h20110)
      // SummationRdtProxy.h11200 (0D_NOT_complex -
      .def_property(
          "h11200", &SummationRdtProxy::h11200, &SummationRdtProxy::set_h11200)
      // SummationRdtProxy.h20020 (0D_NOT_complex -
      .def_property(
          "h20020", &SummationRdtProxy::h20020, &SummationRdtProxy::set_h20020)
      // SummationRdtProxy.h20200 (0D_NOT_complex -
      .def_property(
          "h20200", &SummationRdtProxy::h20200, &SummationRdtProxy::set_h20200)
      // SummationRdtProxy.h00310 (0D_NOT_complex -
      .def_property(
          "h00310", &SummationRdtProxy::h00310, &SummationRdtProxy::set_h00310)
      // SummationRdtProxy.h00400 (0D_NOT_complex -
      .def_property(
          "h00400", &SummationRdtProxy::h00400, &SummationRdtProxy::set_h00400)
      // SummationRdtProxy.h22000 (0D_NOT_complex -
      .def_property(
          "h22000", &SummationRdtProxy::h22000, &SummationRdtProxy::set_h22000)
      // SummationRdtProxy.h00220 (0D_NOT_complex -
      .def_property(
          "h00220", &SummationRdtProxy::h00220, &SummationRdtProxy::set_h00220)
      // SummationRdtProxy.h11110 (0D_NOT_complex -
      .def_property(
          "h11110", &SummationRdtProxy::h11110, &SummationRdtProxy::set_h11110)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return SummationRdtProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const SummationRdtProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SummationRdtProxy& self) {
            return SummationRdtProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const SummationRdtProxy& self, py::dict& memo) {
            return SummationRdtProxy(self);
          })

      ;

  bind_FTypeArrayND<SummationRdtProxyArray1D>(m, "SummationRdtStructArray1D");
  bind_FTypeAlloc1D<SummationRdtProxyAlloc1D>(m, "SummationRdtStructAlloc1D");
  // 2D SummationRdtProxy arrays are not used in structs/routines
  // 3D SummationRdtProxy arrays are not used in structs/routines
}