#include "pybmad/generated/structs_s.hpp"

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// space_charge_common_struct
void init_space_charge_common_struct(py::module &m, py::class_<SpaceChargeCommonStruct> &cls) {
  cls.def(py::init<>())
      // SpaceChargeCommonStruct.ds_track_step (0D_NOT_real - CSR tracking step size
      .def_property(
          "ds_track_step",
          &SpaceChargeCommonStruct::ds_track_step,
          &SpaceChargeCommonStruct::set_ds_track_step
      )
      // SpaceChargeCommonStruct.dt_track_step (0D_NOT_real - Time Runge kutta initial step.
      .def_property(
          "dt_track_step",
          &SpaceChargeCommonStruct::dt_track_step,
          &SpaceChargeCommonStruct::set_dt_track_step
      )
      // SpaceChargeCommonStruct.cathode_strength_cutoff (0D_NOT_real - Cutoff for the cathode field
      // calc.
      .def_property(
          "cathode_strength_cutoff",
          &SpaceChargeCommonStruct::cathode_strength_cutoff,
          &SpaceChargeCommonStruct::set_cathode_strength_cutoff
      )
      // SpaceChargeCommonStruct.rel_tol_tracking (0D_NOT_real - Relative tolerance for tracking.
      .def_property(
          "rel_tol_tracking",
          &SpaceChargeCommonStruct::rel_tol_tracking,
          &SpaceChargeCommonStruct::set_rel_tol_tracking
      )
      // SpaceChargeCommonStruct.abs_tol_tracking (0D_NOT_real - Absolute tolerance for tracking.
      .def_property(
          "abs_tol_tracking",
          &SpaceChargeCommonStruct::abs_tol_tracking,
          &SpaceChargeCommonStruct::set_abs_tol_tracking
      )
      // SpaceChargeCommonStruct.beam_chamber_height (0D_NOT_real - Used in shielding calculation.
      .def_property(
          "beam_chamber_height",
          &SpaceChargeCommonStruct::beam_chamber_height,
          &SpaceChargeCommonStruct::set_beam_chamber_height
      )
      // SpaceChargeCommonStruct.lsc_sigma_cutoff (0D_NOT_real - Cutoff for the 1-dim longitudinal
      // SC calc. If a bin sigma is < cutoff * sigma_ave then ignore.
      .def_property(
          "lsc_sigma_cutoff",
          &SpaceChargeCommonStruct::lsc_sigma_cutoff,
          &SpaceChargeCommonStruct::set_lsc_sigma_cutoff
      )
      // SpaceChargeCommonStruct.particle_sigma_cutoff (0D_NOT_real - 3D SC calc cutoff for
      // particles with (x,y,z) position far from the center. Negative or zero means ignore.
      .def_property(
          "particle_sigma_cutoff",
          &SpaceChargeCommonStruct::particle_sigma_cutoff,
          &SpaceChargeCommonStruct::set_particle_sigma_cutoff
      )
      // SpaceChargeCommonStruct.space_charge_mesh_size (1D_NOT_integer - Gird size for fft_3d space
      // charge calc.
      .def_property_readonly(
          "space_charge_mesh_size",
          &SpaceChargeCommonStruct::space_charge_mesh_size
      )
      // SpaceChargeCommonStruct.csr3d_mesh_size (1D_NOT_integer - Gird size for CSR.
      .def_property_readonly("csr3d_mesh_size", &SpaceChargeCommonStruct::csr3d_mesh_size)
      // SpaceChargeCommonStruct.n_bin (0D_NOT_integer - Number of bins used
      .def_property("n_bin", &SpaceChargeCommonStruct::n_bin, &SpaceChargeCommonStruct::set_n_bin)
      // SpaceChargeCommonStruct.particle_bin_span (0D_NOT_integer - Longitudinal particle length /
      // dz_bin
      .def_property(
          "particle_bin_span",
          &SpaceChargeCommonStruct::particle_bin_span,
          &SpaceChargeCommonStruct::set_particle_bin_span
      )
      // SpaceChargeCommonStruct.n_shield_images (0D_NOT_integer - Chamber wall shielding. 0 = no
      // shielding.
      .def_property(
          "n_shield_images",
          &SpaceChargeCommonStruct::n_shield_images,
          &SpaceChargeCommonStruct::set_n_shield_images
      )
      // SpaceChargeCommonStruct.sc_min_in_bin (0D_NOT_integer - Minimum number of particles in a
      // bin for sigmas to be valid.
      .def_property(
          "sc_min_in_bin",
          &SpaceChargeCommonStruct::sc_min_in_bin,
          &SpaceChargeCommonStruct::set_sc_min_in_bin
      )
      // SpaceChargeCommonStruct.lsc_kick_transverse_dependence (0D_NOT_logical -
      .def_property(
          "lsc_kick_transverse_dependence",
          &SpaceChargeCommonStruct::lsc_kick_transverse_dependence,
          &SpaceChargeCommonStruct::set_lsc_kick_transverse_dependence
      )
      // SpaceChargeCommonStruct.debug (0D_NOT_logical -
      .def_property("debug", &SpaceChargeCommonStruct::debug, &SpaceChargeCommonStruct::set_debug)
      // SpaceChargeCommonStruct.diagnostic_output_file (0D_NOT_character - If non-blank write a
      // diagnostic (EG wake) file
      .def_property(
          "diagnostic_output_file",
          &SpaceChargeCommonStruct::diagnostic_output_file,
          &SpaceChargeCommonStruct::set_diagnostic_output_file
      )

      .def("__repr__", [](const SpaceChargeCommonStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SpaceChargeCommonStruct &self) {
            return SpaceChargeCommonStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SpaceChargeCommonStruct &self, py::dict &memo) {
            return SpaceChargeCommonStruct(self);
          }
      )

      ;

  // 1D SpaceChargeCommonStruct arrays are not used in structs/routines
  // 2D SpaceChargeCommonStruct arrays are not used in structs/routines
  // 3D SpaceChargeCommonStruct arrays are not used in structs/routines
}

// =============================================================================
// spin_axis_struct
void init_spin_axis_struct(py::module &m, py::class_<SpinAxisStruct> &cls) {
  cls.def(py::init<>())
      // SpinAxisStruct.l (1D_NOT_real - Transverse axis.
      .def_property_readonly("l", &SpinAxisStruct::l)
      // SpinAxisStruct.n0 (1D_NOT_real - Invariant spin axis on closed orbit.
      .def_property_readonly("n0", &SpinAxisStruct::n0)
      // SpinAxisStruct.m (1D_NOT_real - Transverse axis.
      .def_property_readonly("m", &SpinAxisStruct::m)

      .def("__repr__", [](const SpinAxisStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SpinAxisStruct &self) {
            return SpinAxisStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SpinAxisStruct &self, py::dict &memo) { return SpinAxisStruct(self); }
      )

      ;

  // 1D SpinAxisStruct arrays are not used in structs/routines
  // 2D SpinAxisStruct arrays are not used in structs/routines
  // 3D SpinAxisStruct arrays are not used in structs/routines
}

// =============================================================================
// spin_orbit_map1_struct
void init_spin_orbit_map1_struct(py::module &m, py::class_<SpinOrbitMap1Struct> &cls) {
  cls.def(py::init<>())
      // SpinOrbitMap1Struct.orb_mat (2D_NOT_real - Orbital matrix
      .def_property_readonly("orb_mat", &SpinOrbitMap1Struct::orb_mat)
      // SpinOrbitMap1Struct.vec0 (1D_NOT_real - Orbital 0th order map: r_out = mat6 * r_in + vec0
      .def_property_readonly("vec0", &SpinOrbitMap1Struct::vec0)
      // SpinOrbitMap1Struct.spin_q (2D_NOT_real - 0th and 1st order quaternion spin map
      .def_property_readonly("spin_q", &SpinOrbitMap1Struct::spin_q)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return SpinOrbitMap1StructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const SpinOrbitMap1Struct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SpinOrbitMap1Struct &self) {
            return SpinOrbitMap1Struct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SpinOrbitMap1Struct &self, py::dict &memo) { return SpinOrbitMap1Struct(self); }
      )

      ;

  bind_FTypeArrayND<SpinOrbitMap1StructArray1D>(m, "SpinOrbitMap1StructArray1D");
  bind_FTypeAlloc1D<SpinOrbitMap1StructAlloc1D>(m, "SpinOrbitMap1StructAlloc1D");
  // 2D SpinOrbitMap1Struct arrays are not used in structs/routines
  // 3D SpinOrbitMap1Struct arrays are not used in structs/routines
}

// =============================================================================
// spin_polar_struct
void init_spin_polar_struct(py::module &m, py::class_<SpinPolarStruct> &cls) {
  cls.def(py::init<>())
      // SpinPolarStruct.polarization (0D_NOT_real -
      .def_property(
          "polarization",
          &SpinPolarStruct::polarization,
          &SpinPolarStruct::set_polarization
      )
      // SpinPolarStruct.theta (0D_NOT_real - Spherical coords: Angle from z-axis.
      .def_property("theta", &SpinPolarStruct::theta, &SpinPolarStruct::set_theta)
      // SpinPolarStruct.phi (0D_NOT_real - Spherical coords: Angle in (x,y) plane.
      .def_property("phi", &SpinPolarStruct::phi, &SpinPolarStruct::set_phi)
      // SpinPolarStruct.xi (0D_NOT_real - Spinor phase angle (See Bmad manual).
      .def_property("xi", &SpinPolarStruct::xi, &SpinPolarStruct::set_xi)

      .def("__repr__", [](const SpinPolarStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SpinPolarStruct &self) {
            return SpinPolarStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SpinPolarStruct &self, py::dict &memo) { return SpinPolarStruct(self); }
      )

      ;

  // 1D SpinPolarStruct arrays are not used in structs/routines
  // 2D SpinPolarStruct arrays are not used in structs/routines
  // 3D SpinPolarStruct arrays are not used in structs/routines
}

// =============================================================================
// strong_beam_struct
void init_strong_beam_struct(py::module &m, py::class_<StrongBeamStruct> &cls) {
  cls.def(py::init<>())
      // StrongBeamStruct.ix_slice (0D_NOT_integer - 0 -> at element center and not at slice.
      .def_property("ix_slice", &StrongBeamStruct::ix_slice, &StrongBeamStruct::set_ix_slice)
      // StrongBeamStruct.x_center (0D_NOT_real - Strong beam slice center.
      .def_property("x_center", &StrongBeamStruct::x_center, &StrongBeamStruct::set_x_center)
      // StrongBeamStruct.y_center (0D_NOT_real - Strong beam slice center.
      .def_property("y_center", &StrongBeamStruct::y_center, &StrongBeamStruct::set_y_center)
      // StrongBeamStruct.x_sigma (0D_NOT_real - Strong beam slice sigma.
      .def_property("x_sigma", &StrongBeamStruct::x_sigma, &StrongBeamStruct::set_x_sigma)
      // StrongBeamStruct.y_sigma (0D_NOT_real - Strong beam slice sigma.
      .def_property("y_sigma", &StrongBeamStruct::y_sigma, &StrongBeamStruct::set_y_sigma)
      // StrongBeamStruct.dx (0D_NOT_real - Particle - beam slice distance.
      .def_property("dx", &StrongBeamStruct::dx, &StrongBeamStruct::set_dx)
      // StrongBeamStruct.dy (0D_NOT_real - Particle - beam slice distance.
      .def_property("dy", &StrongBeamStruct::dy, &StrongBeamStruct::set_dy)

      .def("__repr__", [](const StrongBeamStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const StrongBeamStruct &self) {
            return StrongBeamStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const StrongBeamStruct &self, py::dict &memo) { return StrongBeamStruct(self); }
      )

      ;

  // 1D StrongBeamStruct arrays are not used in structs/routines
  // 2D StrongBeamStruct arrays are not used in structs/routines
  // 3D StrongBeamStruct arrays are not used in structs/routines
}

// =============================================================================
// surface_curvature_struct
void init_surface_curvature_struct(py::module &m, py::class_<SurfaceCurvatureStruct> &cls) {
  cls.def(py::init<>())
      // SurfaceCurvatureStruct.xy (2D_NOT_real -
      .def_property_readonly("xy", &SurfaceCurvatureStruct::xy)
      // SurfaceCurvatureStruct.spherical (0D_NOT_real -
      .def_property(
          "spherical",
          &SurfaceCurvatureStruct::spherical,
          &SurfaceCurvatureStruct::set_spherical
      )
      // SurfaceCurvatureStruct.elliptical (1D_NOT_real - Total curvature = elliptical + spherical
      .def_property_readonly("elliptical", &SurfaceCurvatureStruct::elliptical)
      // SurfaceCurvatureStruct.has_curvature (0D_NOT_logical - Dependent var. Will be set by Bmad
      .def_property(
          "has_curvature",
          &SurfaceCurvatureStruct::has_curvature,
          &SurfaceCurvatureStruct::set_has_curvature
      )

      .def("__repr__", [](const SurfaceCurvatureStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceCurvatureStruct &self) {
            return SurfaceCurvatureStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SurfaceCurvatureStruct &self, py::dict &memo) {
            return SurfaceCurvatureStruct(self);
          }
      )

      ;

  // 1D SurfaceCurvatureStruct arrays are not used in structs/routines
  // 2D SurfaceCurvatureStruct arrays are not used in structs/routines
  // 3D SurfaceCurvatureStruct arrays are not used in structs/routines
}

// =============================================================================
// surface_displacement_pt_struct
void init_surface_displacement_pt_struct(
    py::module &m,
    py::class_<SurfaceDisplacementPtStruct> &cls
) {
  cls.def(py::init<>())
      // SurfaceDisplacementPtStruct.x0 (0D_NOT_real - Position at center
      .def_property("x0", &SurfaceDisplacementPtStruct::x0, &SurfaceDisplacementPtStruct::set_x0)
      // SurfaceDisplacementPtStruct.y0 (0D_NOT_real - Position at center
      .def_property("y0", &SurfaceDisplacementPtStruct::y0, &SurfaceDisplacementPtStruct::set_y0)
      // SurfaceDisplacementPtStruct.z0 (0D_NOT_real -
      .def_property("z0", &SurfaceDisplacementPtStruct::z0, &SurfaceDisplacementPtStruct::set_z0)
      // SurfaceDisplacementPtStruct.dz_dx (0D_NOT_real -
      .def_property(
          "dz_dx",
          &SurfaceDisplacementPtStruct::dz_dx,
          &SurfaceDisplacementPtStruct::set_dz_dx
      )
      // SurfaceDisplacementPtStruct.dz_dy (0D_NOT_real -
      .def_property(
          "dz_dy",
          &SurfaceDisplacementPtStruct::dz_dy,
          &SurfaceDisplacementPtStruct::set_dz_dy
      )
      // SurfaceDisplacementPtStruct.d2z_dxdy (0D_NOT_real -
      .def_property(
          "d2z_dxdy",
          &SurfaceDisplacementPtStruct::d2z_dxdy,
          &SurfaceDisplacementPtStruct::set_d2z_dxdy
      )

      .def("__repr__", [](const SurfaceDisplacementPtStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceDisplacementPtStruct &self) {
            return SurfaceDisplacementPtStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SurfaceDisplacementPtStruct &self, py::dict &memo) {
            return SurfaceDisplacementPtStruct(self);
          }
      )

      ;

  // 1D SurfaceDisplacementPtStruct arrays are not used in structs/routines
  bind_FTypeArrayND<SurfaceDisplacementPtStructArray2D>(m, "SurfaceDisplacementPtStructArray2D");
  // 3D SurfaceDisplacementPtStruct arrays are not used in structs/routines
}

// =============================================================================
// surface_displacement_struct
void init_surface_displacement_struct(py::module &m, py::class_<SurfaceDisplacementStruct> &cls) {
  cls.def(py::init<>())
      // SurfaceDisplacementStruct.active (0D_NOT_logical -
      .def_property(
          "active",
          &SurfaceDisplacementStruct::active,
          &SurfaceDisplacementStruct::set_active
      )
      // SurfaceDisplacementStruct.dr (1D_NOT_real -
      .def_property_readonly("dr", &SurfaceDisplacementStruct::dr)
      // SurfaceDisplacementStruct.r0 (1D_NOT_real -
      .def_property_readonly("r0", &SurfaceDisplacementStruct::r0)
      // SurfaceDisplacementStruct.pt (2D_ALLOC_type -
      .def_property_readonly("pt", &SurfaceDisplacementStruct::pt)

      .def("__repr__", [](const SurfaceDisplacementStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceDisplacementStruct &self) {
            return SurfaceDisplacementStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SurfaceDisplacementStruct &self, py::dict &memo) {
            return SurfaceDisplacementStruct(self);
          }
      )

      ;

  // 1D SurfaceDisplacementStruct arrays are not used in structs/routines
  // 2D SurfaceDisplacementStruct arrays are not used in structs/routines
  // 3D SurfaceDisplacementStruct arrays are not used in structs/routines
}

// =============================================================================
// surface_h_misalign_pt_struct
void init_surface_h_misalign_pt_struct(py::module &m, py::class_<SurfaceHMisalignPtStruct> &cls) {
  cls.def(py::init<>())
      // SurfaceHMisalignPtStruct.x0 (0D_NOT_real - Position at center
      .def_property("x0", &SurfaceHMisalignPtStruct::x0, &SurfaceHMisalignPtStruct::set_x0)
      // SurfaceHMisalignPtStruct.y0 (0D_NOT_real - Position at center
      .def_property("y0", &SurfaceHMisalignPtStruct::y0, &SurfaceHMisalignPtStruct::set_y0)
      // SurfaceHMisalignPtStruct.rot_y (0D_NOT_real - rot_t = x-rotation for Bragg and z-rotation
      // for Laue.
      .def_property("rot_y", &SurfaceHMisalignPtStruct::rot_y, &SurfaceHMisalignPtStruct::set_rot_y)
      // SurfaceHMisalignPtStruct.rot_t (0D_NOT_real - rot_t = x-rotation for Bragg and z-rotation
      // for Laue.
      .def_property("rot_t", &SurfaceHMisalignPtStruct::rot_t, &SurfaceHMisalignPtStruct::set_rot_t)
      // SurfaceHMisalignPtStruct.rot_y_rms (0D_NOT_real - rot_t = x-rotation for Bragg and
      // z-rotation for Laue.
      .def_property(
          "rot_y_rms",
          &SurfaceHMisalignPtStruct::rot_y_rms,
          &SurfaceHMisalignPtStruct::set_rot_y_rms
      )
      // SurfaceHMisalignPtStruct.rot_t_rms (0D_NOT_real - rot_t = x-rotation for Bragg and
      // z-rotation for Laue.
      .def_property(
          "rot_t_rms",
          &SurfaceHMisalignPtStruct::rot_t_rms,
          &SurfaceHMisalignPtStruct::set_rot_t_rms
      )

      .def("__repr__", [](const SurfaceHMisalignPtStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceHMisalignPtStruct &self) {
            return SurfaceHMisalignPtStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SurfaceHMisalignPtStruct &self, py::dict &memo) {
            return SurfaceHMisalignPtStruct(self);
          }
      )

      ;

  // 1D SurfaceHMisalignPtStruct arrays are not used in structs/routines
  bind_FTypeArrayND<SurfaceHMisalignPtStructArray2D>(m, "SurfaceHMisalignPtStructArray2D");
  // 3D SurfaceHMisalignPtStruct arrays are not used in structs/routines
}

// =============================================================================
// surface_h_misalign_struct
void init_surface_h_misalign_struct(py::module &m, py::class_<SurfaceHMisalignStruct> &cls) {
  cls.def(py::init<>())
      // SurfaceHMisalignStruct.active (0D_NOT_logical -
      .def_property("active", &SurfaceHMisalignStruct::active, &SurfaceHMisalignStruct::set_active)
      // SurfaceHMisalignStruct.dr (1D_NOT_real -
      .def_property_readonly("dr", &SurfaceHMisalignStruct::dr)
      // SurfaceHMisalignStruct.r0 (1D_NOT_real -
      .def_property_readonly("r0", &SurfaceHMisalignStruct::r0)
      // SurfaceHMisalignStruct.pt (2D_ALLOC_type -
      .def_property_readonly("pt", &SurfaceHMisalignStruct::pt)

      .def("__repr__", [](const SurfaceHMisalignStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceHMisalignStruct &self) {
            return SurfaceHMisalignStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SurfaceHMisalignStruct &self, py::dict &memo) {
            return SurfaceHMisalignStruct(self);
          }
      )

      ;

  // 1D SurfaceHMisalignStruct arrays are not used in structs/routines
  // 2D SurfaceHMisalignStruct arrays are not used in structs/routines
  // 3D SurfaceHMisalignStruct arrays are not used in structs/routines
}

// =============================================================================
// surface_segmented_pt_struct
void init_surface_segmented_pt_struct(py::module &m, py::class_<SurfaceSegmentedPtStruct> &cls) {
  cls.def(py::init<>())
      // SurfaceSegmentedPtStruct.x0 (0D_NOT_real - Position at center
      .def_property("x0", &SurfaceSegmentedPtStruct::x0, &SurfaceSegmentedPtStruct::set_x0)
      // SurfaceSegmentedPtStruct.y0 (0D_NOT_real - Position at center
      .def_property("y0", &SurfaceSegmentedPtStruct::y0, &SurfaceSegmentedPtStruct::set_y0)
      // SurfaceSegmentedPtStruct.z0 (0D_NOT_real - Position at center
      .def_property("z0", &SurfaceSegmentedPtStruct::z0, &SurfaceSegmentedPtStruct::set_z0)
      // SurfaceSegmentedPtStruct.dz_dx (0D_NOT_real - Slope at center
      .def_property("dz_dx", &SurfaceSegmentedPtStruct::dz_dx, &SurfaceSegmentedPtStruct::set_dz_dx)
      // SurfaceSegmentedPtStruct.dz_dy (0D_NOT_real - Slope at center
      .def_property("dz_dy", &SurfaceSegmentedPtStruct::dz_dy, &SurfaceSegmentedPtStruct::set_dz_dy)

      .def("__repr__", [](const SurfaceSegmentedPtStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceSegmentedPtStruct &self) {
            return SurfaceSegmentedPtStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SurfaceSegmentedPtStruct &self, py::dict &memo) {
            return SurfaceSegmentedPtStruct(self);
          }
      )

      ;

  // 1D SurfaceSegmentedPtStruct arrays are not used in structs/routines
  bind_FTypeArrayND<SurfaceSegmentedPtStructArray2D>(m, "SurfaceSegmentedPtStructArray2D");
  // 3D SurfaceSegmentedPtStruct arrays are not used in structs/routines
}

// =============================================================================
// surface_segmented_struct
void init_surface_segmented_struct(py::module &m, py::class_<SurfaceSegmentedStruct> &cls) {
  cls.def(py::init<>())
      // SurfaceSegmentedStruct.active (0D_NOT_logical -
      .def_property("active", &SurfaceSegmentedStruct::active, &SurfaceSegmentedStruct::set_active)
      // SurfaceSegmentedStruct.dr (1D_NOT_real -
      .def_property_readonly("dr", &SurfaceSegmentedStruct::dr)
      // SurfaceSegmentedStruct.r0 (1D_NOT_real -
      .def_property_readonly("r0", &SurfaceSegmentedStruct::r0)
      // SurfaceSegmentedStruct.pt (2D_ALLOC_type -
      .def_property_readonly("pt", &SurfaceSegmentedStruct::pt)

      .def("__repr__", [](const SurfaceSegmentedStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceSegmentedStruct &self) {
            return SurfaceSegmentedStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SurfaceSegmentedStruct &self, py::dict &memo) {
            return SurfaceSegmentedStruct(self);
          }
      )

      ;

  // 1D SurfaceSegmentedStruct arrays are not used in structs/routines
  // 2D SurfaceSegmentedStruct arrays are not used in structs/routines
  // 3D SurfaceSegmentedStruct arrays are not used in structs/routines
}

// =============================================================================
// spline_struct
void init_spline_struct(py::module &m, py::class_<SplineStruct> &cls) {
  cls.def(py::init<>())
      // SplineStruct.x0 (0D_NOT_real - Point at start of spline
      .def_property("x0", &SplineStruct::x0, &SplineStruct::set_x0)
      // SplineStruct.y0 (0D_NOT_real - Point at start of spline
      .def_property("y0", &SplineStruct::y0, &SplineStruct::set_y0)
      // SplineStruct.x1 (0D_NOT_real - Point at end of spline
      .def_property("x1", &SplineStruct::x1, &SplineStruct::set_x1)
      // SplineStruct.coef (1D_NOT_real - coefficients for cubic spline
      .def_property_readonly("coef", &SplineStruct::coef)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return SplineStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const SplineStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SplineStruct &self) {
            return SplineStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SplineStruct &self, py::dict &memo) { return SplineStruct(self); }
      )

      ;

  bind_FTypeArrayND<SplineStructArray1D>(m, "SplineStructArray1D");
  bind_FTypeAlloc1D<SplineStructAlloc1D>(m, "SplineStructAlloc1D");
  // 2D SplineStruct arrays are not used in structs/routines
  // 3D SplineStruct arrays are not used in structs/routines
}

// =============================================================================
// summation_rdt_struct
void init_summation_rdt_struct(py::module &m, py::class_<SummationRdtStruct> &cls) {
  cls.def(py::init<>())
      // SummationRdtStruct.h11001 (0D_NOT_complex -
      .def_property("h11001", &SummationRdtStruct::h11001, &SummationRdtStruct::set_h11001)
      // SummationRdtStruct.h00111 (0D_NOT_complex -
      .def_property("h00111", &SummationRdtStruct::h00111, &SummationRdtStruct::set_h00111)
      // SummationRdtStruct.h20001 (0D_NOT_complex -
      .def_property("h20001", &SummationRdtStruct::h20001, &SummationRdtStruct::set_h20001)
      // SummationRdtStruct.h00201 (0D_NOT_complex -
      .def_property("h00201", &SummationRdtStruct::h00201, &SummationRdtStruct::set_h00201)
      // SummationRdtStruct.h10002 (0D_NOT_complex -
      .def_property("h10002", &SummationRdtStruct::h10002, &SummationRdtStruct::set_h10002)
      // SummationRdtStruct.h21000 (0D_NOT_complex -
      .def_property("h21000", &SummationRdtStruct::h21000, &SummationRdtStruct::set_h21000)
      // SummationRdtStruct.h30000 (0D_NOT_complex -
      .def_property("h30000", &SummationRdtStruct::h30000, &SummationRdtStruct::set_h30000)
      // SummationRdtStruct.h10110 (0D_NOT_complex -
      .def_property("h10110", &SummationRdtStruct::h10110, &SummationRdtStruct::set_h10110)
      // SummationRdtStruct.h10020 (0D_NOT_complex -
      .def_property("h10020", &SummationRdtStruct::h10020, &SummationRdtStruct::set_h10020)
      // SummationRdtStruct.h10200 (0D_NOT_complex - 2nd order in K2 moments
      .def_property("h10200", &SummationRdtStruct::h10200, &SummationRdtStruct::set_h10200)
      // SummationRdtStruct.h31000 (0D_NOT_complex -
      .def_property("h31000", &SummationRdtStruct::h31000, &SummationRdtStruct::set_h31000)
      // SummationRdtStruct.h40000 (0D_NOT_complex -
      .def_property("h40000", &SummationRdtStruct::h40000, &SummationRdtStruct::set_h40000)
      // SummationRdtStruct.h20110 (0D_NOT_complex -
      .def_property("h20110", &SummationRdtStruct::h20110, &SummationRdtStruct::set_h20110)
      // SummationRdtStruct.h11200 (0D_NOT_complex -
      .def_property("h11200", &SummationRdtStruct::h11200, &SummationRdtStruct::set_h11200)
      // SummationRdtStruct.h20020 (0D_NOT_complex -
      .def_property("h20020", &SummationRdtStruct::h20020, &SummationRdtStruct::set_h20020)
      // SummationRdtStruct.h20200 (0D_NOT_complex -
      .def_property("h20200", &SummationRdtStruct::h20200, &SummationRdtStruct::set_h20200)
      // SummationRdtStruct.h00310 (0D_NOT_complex -
      .def_property("h00310", &SummationRdtStruct::h00310, &SummationRdtStruct::set_h00310)
      // SummationRdtStruct.h00400 (0D_NOT_complex -
      .def_property("h00400", &SummationRdtStruct::h00400, &SummationRdtStruct::set_h00400)
      // SummationRdtStruct.h22000 (0D_NOT_complex -
      .def_property("h22000", &SummationRdtStruct::h22000, &SummationRdtStruct::set_h22000)
      // SummationRdtStruct.h00220 (0D_NOT_complex -
      .def_property("h00220", &SummationRdtStruct::h00220, &SummationRdtStruct::set_h00220)
      // SummationRdtStruct.h11110 (0D_NOT_complex -
      .def_property("h11110", &SummationRdtStruct::h11110, &SummationRdtStruct::set_h11110)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return SummationRdtStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const SummationRdtStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SummationRdtStruct &self) {
            return SummationRdtStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SummationRdtStruct &self, py::dict &memo) { return SummationRdtStruct(self); }
      )

      ;

  bind_FTypeArrayND<SummationRdtStructArray1D>(m, "SummationRdtStructArray1D");
  bind_FTypeAlloc1D<SummationRdtStructAlloc1D>(m, "SummationRdtStructAlloc1D");
  // 2D SummationRdtStruct arrays are not used in structs/routines
  // 3D SummationRdtStruct arrays are not used in structs/routines
}