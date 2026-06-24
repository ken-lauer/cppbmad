#include "pybmad/generated/structs_s.hpp"

#include <cstdint>
#include <functional>

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// space_charge_common_struct
void init_space_charge_common_struct(py::module &m, py::class_<SpaceChargeCommonStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<int>>,
             std::optional<std::vector<int>>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<std::string>>(),
         py::arg("ds_track_step") = py::none(),
         py::arg("dt_track_step") = py::none(),
         py::arg("cathode_strength_cutoff") = py::none(),
         py::arg("rel_tol_tracking") = py::none(),
         py::arg("abs_tol_tracking") = py::none(),
         py::arg("beam_chamber_height") = py::none(),
         py::arg("lsc_sigma_cutoff") = py::none(),
         py::arg("particle_sigma_cutoff") = py::none(),
         py::arg("mesh_growth_factor") = py::none(),
         py::arg("mesh_shrink_factor") = py::none(),
         py::arg("space_charge_mesh_size") = py::none(),
         py::arg("csr3d_mesh_size") = py::none(),
         py::arg("n_bin") = py::none(),
         py::arg("particle_bin_span") = py::none(),
         py::arg("n_shield_images") = py::none(),
         py::arg("sc_min_in_bin") = py::none(),
         py::arg("lsc_kick_transverse_dependence") = py::none(),
         py::arg("debug") = py::none(),
         py::arg("diagnostic_output_file") = py::none()
  )
      .def_property(
          "ds_track_step",
          &SpaceChargeCommonStruct::ds_track_step,
          &SpaceChargeCommonStruct::set_ds_track_step,
          "CSR tracking step size"
      )
      .def_property(
          "dt_track_step",
          &SpaceChargeCommonStruct::dt_track_step,
          &SpaceChargeCommonStruct::set_dt_track_step,
          "Time Runge kutta initial step."
      )
      .def_property(
          "cathode_strength_cutoff",
          &SpaceChargeCommonStruct::cathode_strength_cutoff,
          &SpaceChargeCommonStruct::set_cathode_strength_cutoff,
          "Cutoff for the cathode field calc."
      )
      .def_property(
          "rel_tol_tracking",
          &SpaceChargeCommonStruct::rel_tol_tracking,
          &SpaceChargeCommonStruct::set_rel_tol_tracking,
          "Relative tolerance for tracking."
      )
      .def_property(
          "abs_tol_tracking",
          &SpaceChargeCommonStruct::abs_tol_tracking,
          &SpaceChargeCommonStruct::set_abs_tol_tracking,
          "Absolute tolerance for tracking."
      )
      .def_property(
          "beam_chamber_height",
          &SpaceChargeCommonStruct::beam_chamber_height,
          &SpaceChargeCommonStruct::set_beam_chamber_height,
          "Used in shielding calculation."
      )
      .def_property(
          "lsc_sigma_cutoff",
          &SpaceChargeCommonStruct::lsc_sigma_cutoff,
          &SpaceChargeCommonStruct::set_lsc_sigma_cutoff,
          "Cutoff for the 1-dim longitudinal SC calc. If a bin sigma is < cutoff * sigma_ave then "
          "ignore."
      )
      .def_property(
          "particle_sigma_cutoff",
          &SpaceChargeCommonStruct::particle_sigma_cutoff,
          &SpaceChargeCommonStruct::set_particle_sigma_cutoff,
          "3D SC calc cutoff for particles with (x,y,z) position far from the center. Negative or "
          "zero means ignore."
      )
      .def_property(
          "mesh_growth_factor",
          &SpaceChargeCommonStruct::mesh_growth_factor,
          &SpaceChargeCommonStruct::set_mesh_growth_factor,
          "Fractional padding when growing SC mesh (default: 10%). Set to 0 for tight-fit (no "
          "caching speedup)."
      )
      .def_property(
          "mesh_shrink_factor",
          &SpaceChargeCommonStruct::mesh_shrink_factor,
          &SpaceChargeCommonStruct::set_mesh_shrink_factor,
          "Fractional threshold for shrinking SC mesh (default: 10%). Mesh shrinks when bunch "
          "fills < (1-this) of the mesh range."
      )
      .def_property(
          "space_charge_mesh_size",
          py::cpp_function(
              &SpaceChargeCommonStruct::space_charge_mesh_size,
              py::keep_alive<0, 1>()
          ),
          &SpaceChargeCommonStruct::set_space_charge_mesh_size,
          "Gird size for fft_3d space charge calc."
      )
      .def_property(
          "csr3d_mesh_size",
          py::cpp_function(&SpaceChargeCommonStruct::csr3d_mesh_size, py::keep_alive<0, 1>()),
          &SpaceChargeCommonStruct::set_csr3d_mesh_size,
          "Gird size for CSR."
      )
      .def_property(
          "n_bin",
          &SpaceChargeCommonStruct::n_bin,
          &SpaceChargeCommonStruct::set_n_bin,
          "Number of bins used"
      )
      .def_property(
          "particle_bin_span",
          &SpaceChargeCommonStruct::particle_bin_span,
          &SpaceChargeCommonStruct::set_particle_bin_span,
          "Longitudinal particle length / dz_bin"
      )
      .def_property(
          "n_shield_images",
          &SpaceChargeCommonStruct::n_shield_images,
          &SpaceChargeCommonStruct::set_n_shield_images,
          "Chamber wall shielding. 0 = no shielding."
      )
      .def_property(
          "sc_min_in_bin",
          &SpaceChargeCommonStruct::sc_min_in_bin,
          &SpaceChargeCommonStruct::set_sc_min_in_bin,
          "Minimum number of particles in a bin for sigmas to be valid."
      )
      .def_property(
          "lsc_kick_transverse_dependence",
          &SpaceChargeCommonStruct::lsc_kick_transverse_dependence,
          &SpaceChargeCommonStruct::set_lsc_kick_transverse_dependence
      )
      .def_property("debug", &SpaceChargeCommonStruct::debug, &SpaceChargeCommonStruct::set_debug)
      .def_property(
          "diagnostic_output_file",
          &SpaceChargeCommonStruct::diagnostic_output_file,
          &SpaceChargeCommonStruct::set_diagnostic_output_file,
          "If non-blank write a diagnostic (EG wake) file"
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
      .def(
          "__eq__",
          [](const SpaceChargeCommonStruct &self, const SpaceChargeCommonStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const SpaceChargeCommonStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
  cls.def(
         py::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>>(),
         py::arg("l") = py::none(),
         py::arg("n0") = py::none(),
         py::arg("m") = py::none()
  )
      .def_property(
          "l",
          py::cpp_function(&SpinAxisStruct::l, py::keep_alive<0, 1>()),
          &SpinAxisStruct::set_l,
          "Transverse axis."
      )
      .def_property(
          "n0",
          py::cpp_function(&SpinAxisStruct::n0, py::keep_alive<0, 1>()),
          &SpinAxisStruct::set_n0,
          "Invariant spin axis on closed orbit."
      )
      .def_property(
          "m",
          py::cpp_function(&SpinAxisStruct::m, py::keep_alive<0, 1>()),
          &SpinAxisStruct::set_m,
          "Transverse axis."
      )

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
      .def(
          "__eq__",
          [](const SpinAxisStruct &self, const SpinAxisStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const SpinAxisStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D SpinAxisStruct arrays are not used in structs/routines
  // 2D SpinAxisStruct arrays are not used in structs/routines
  // 3D SpinAxisStruct arrays are not used in structs/routines
}

// =============================================================================
// spin_orbit_map1_struct
void init_spin_orbit_map1_struct(py::module &m, py::class_<SpinOrbitMap1Struct> &cls) {
  cls.def(
         py::init<
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>>(),
         py::arg("orb_mat") = py::none(),
         py::arg("vec0") = py::none(),
         py::arg("spin_q") = py::none()
  )
      .def_property(
          "orb_mat",
          py::cpp_function(&SpinOrbitMap1Struct::orb_mat, py::keep_alive<0, 1>()),
          &SpinOrbitMap1Struct::set_orb_mat,
          "Orbital matrix"
      )
      .def_property(
          "vec0",
          py::cpp_function(&SpinOrbitMap1Struct::vec0, py::keep_alive<0, 1>()),
          &SpinOrbitMap1Struct::set_vec0,
          "Orbital 0th order map: r_out = mat6 * r_in + vec0"
      )
      .def_property(
          "spin_q",
          py::cpp_function(&SpinOrbitMap1Struct::spin_q, py::keep_alive<0, 1>()),
          &SpinOrbitMap1Struct::set_spin_q,
          "0th and 1st order quaternion spin map"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return SpinOrbitMap1StructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = SpinOrbitMap1StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const SpinOrbitMap1Struct &self, const SpinOrbitMap1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const SpinOrbitMap1Struct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<SpinOrbitMap1StructArray1D, SpinOrbitMap1StructAlloc1D>(
      m,
      "SpinOrbitMap1StructArray1D",
      "SpinOrbitMap1StructAlloc1D"
  );
  // 2D SpinOrbitMap1Struct arrays are not used in structs/routines
  // 3D SpinOrbitMap1Struct arrays are not used in structs/routines
}

// =============================================================================
// spin_polar_struct
void init_spin_polar_struct(py::module &m, py::class_<SpinPolarStruct> &cls) {
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("polarization") = py::none(),
         py::arg("theta") = py::none(),
         py::arg("phi") = py::none(),
         py::arg("xi") = py::none()
  )
      .def_property(
          "polarization",
          &SpinPolarStruct::polarization,
          &SpinPolarStruct::set_polarization
      )
      .def_property(
          "theta",
          &SpinPolarStruct::theta,
          &SpinPolarStruct::set_theta,
          "Spherical coords: Angle from z-axis."
      )
      .def_property(
          "phi",
          &SpinPolarStruct::phi,
          &SpinPolarStruct::set_phi,
          "Spherical coords: Angle in (x,y) plane."
      )
      .def_property(
          "xi",
          &SpinPolarStruct::xi,
          &SpinPolarStruct::set_xi,
          "Spinor phase angle (See Bmad manual)."
      )

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
      .def(
          "__eq__",
          [](const SpinPolarStruct &self, const SpinPolarStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const SpinPolarStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D SpinPolarStruct arrays are not used in structs/routines
  // 2D SpinPolarStruct arrays are not used in structs/routines
  // 3D SpinPolarStruct arrays are not used in structs/routines
}

// =============================================================================
// strong_beam_struct
void init_strong_beam_struct(py::module &m, py::class_<StrongBeamStruct> &cls) {
  cls.def(
         py::init<
             std::optional<int>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("ix_slice") = py::none(),
         py::arg("x_center") = py::none(),
         py::arg("y_center") = py::none(),
         py::arg("x_sigma") = py::none(),
         py::arg("y_sigma") = py::none(),
         py::arg("dx") = py::none(),
         py::arg("dy") = py::none()
  )
      .def_property(
          "ix_slice",
          &StrongBeamStruct::ix_slice,
          &StrongBeamStruct::set_ix_slice,
          "0 -> at element center and not at slice."
      )
      .def_property(
          "x_center",
          &StrongBeamStruct::x_center,
          &StrongBeamStruct::set_x_center,
          "Strong beam slice center."
      )
      .def_property(
          "y_center",
          &StrongBeamStruct::y_center,
          &StrongBeamStruct::set_y_center,
          "Strong beam slice center."
      )
      .def_property(
          "x_sigma",
          &StrongBeamStruct::x_sigma,
          &StrongBeamStruct::set_x_sigma,
          "Strong beam slice sigma."
      )
      .def_property(
          "y_sigma",
          &StrongBeamStruct::y_sigma,
          &StrongBeamStruct::set_y_sigma,
          "Strong beam slice sigma."
      )
      .def_property(
          "dx",
          &StrongBeamStruct::dx,
          &StrongBeamStruct::set_dx,
          "Particle - beam slice distance."
      )
      .def_property(
          "dy",
          &StrongBeamStruct::dy,
          &StrongBeamStruct::set_dy,
          "Particle - beam slice distance."
      )

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
      .def(
          "__eq__",
          [](const StrongBeamStruct &self, const StrongBeamStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const StrongBeamStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D StrongBeamStruct arrays are not used in structs/routines
  // 2D StrongBeamStruct arrays are not used in structs/routines
  // 3D StrongBeamStruct arrays are not used in structs/routines
}

// =============================================================================
// surface_curvature_struct
void init_surface_curvature_struct(py::module &m, py::class_<SurfaceCurvatureStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::vector<std::vector<double>>>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<bool>>(),
         py::arg("xy") = py::none(),
         py::arg("spherical") = py::none(),
         py::arg("elliptical") = py::none(),
         py::arg("has_curvature") = py::none()
  )
      .def_property(
          "xy",
          py::cpp_function(&SurfaceCurvatureStruct::xy, py::keep_alive<0, 1>()),
          &SurfaceCurvatureStruct::set_xy
      )
      .def_property(
          "spherical",
          &SurfaceCurvatureStruct::spherical,
          &SurfaceCurvatureStruct::set_spherical
      )
      .def_property(
          "elliptical",
          py::cpp_function(&SurfaceCurvatureStruct::elliptical, py::keep_alive<0, 1>()),
          &SurfaceCurvatureStruct::set_elliptical,
          "Total curvature = elliptical + spherical"
      )
      .def_property(
          "has_curvature",
          &SurfaceCurvatureStruct::has_curvature,
          &SurfaceCurvatureStruct::set_has_curvature,
          "Dependent var. Will be set by Bmad"
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
      .def(
          "__eq__",
          [](const SurfaceCurvatureStruct &self, const SurfaceCurvatureStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const SurfaceCurvatureStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("x0") = py::none(),
         py::arg("y0") = py::none(),
         py::arg("z0") = py::none(),
         py::arg("dz_dx") = py::none(),
         py::arg("dz_dy") = py::none(),
         py::arg("d2z_dxdy") = py::none()
  )
      .def_property(
          "x0",
          &SurfaceDisplacementPtStruct::x0,
          &SurfaceDisplacementPtStruct::set_x0,
          "Position at center"
      )
      .def_property(
          "y0",
          &SurfaceDisplacementPtStruct::y0,
          &SurfaceDisplacementPtStruct::set_y0,
          "Position at center"
      )
      .def_property("z0", &SurfaceDisplacementPtStruct::z0, &SurfaceDisplacementPtStruct::set_z0)
      .def_property(
          "dz_dx",
          &SurfaceDisplacementPtStruct::dz_dx,
          &SurfaceDisplacementPtStruct::set_dz_dx
      )
      .def_property(
          "dz_dy",
          &SurfaceDisplacementPtStruct::dz_dy,
          &SurfaceDisplacementPtStruct::set_dz_dy
      )
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
      .def(
          "__eq__",
          [](const SurfaceDisplacementPtStruct &self, const SurfaceDisplacementPtStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const SurfaceDisplacementPtStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
  cls.def(
         py::init<
             std::optional<bool>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>>(),
         py::arg("active") = py::none(),
         py::arg("dr") = py::none(),
         py::arg("r0") = py::none()
  )
      .def_property(
          "active",
          &SurfaceDisplacementStruct::active,
          &SurfaceDisplacementStruct::set_active
      )
      .def_property(
          "dr",
          py::cpp_function(&SurfaceDisplacementStruct::dr, py::keep_alive<0, 1>()),
          &SurfaceDisplacementStruct::set_dr
      )
      .def_property(
          "r0",
          py::cpp_function(&SurfaceDisplacementStruct::r0, py::keep_alive<0, 1>()),
          &SurfaceDisplacementStruct::set_r0
      )
      .def_property_readonly(
          "pt",
          py::cpp_function(&SurfaceDisplacementStruct::pt, py::keep_alive<0, 1>())
      )

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
      .def(
          "__eq__",
          [](const SurfaceDisplacementStruct &self, const SurfaceDisplacementStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const SurfaceDisplacementStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("x0") = py::none(),
         py::arg("y0") = py::none(),
         py::arg("rot_y") = py::none(),
         py::arg("rot_t") = py::none(),
         py::arg("rot_y_rms") = py::none(),
         py::arg("rot_t_rms") = py::none()
  )
      .def_property(
          "x0",
          &SurfaceHMisalignPtStruct::x0,
          &SurfaceHMisalignPtStruct::set_x0,
          "Position at center"
      )
      .def_property(
          "y0",
          &SurfaceHMisalignPtStruct::y0,
          &SurfaceHMisalignPtStruct::set_y0,
          "Position at center"
      )
      .def_property(
          "rot_y",
          &SurfaceHMisalignPtStruct::rot_y,
          &SurfaceHMisalignPtStruct::set_rot_y,
          "rot_t = x-rotation for Bragg and z-rotation for Laue."
      )
      .def_property(
          "rot_t",
          &SurfaceHMisalignPtStruct::rot_t,
          &SurfaceHMisalignPtStruct::set_rot_t,
          "rot_t = x-rotation for Bragg and z-rotation for Laue."
      )
      .def_property(
          "rot_y_rms",
          &SurfaceHMisalignPtStruct::rot_y_rms,
          &SurfaceHMisalignPtStruct::set_rot_y_rms,
          "rot_t = x-rotation for Bragg and z-rotation for Laue."
      )
      .def_property(
          "rot_t_rms",
          &SurfaceHMisalignPtStruct::rot_t_rms,
          &SurfaceHMisalignPtStruct::set_rot_t_rms,
          "rot_t = x-rotation for Bragg and z-rotation for Laue."
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
      .def(
          "__eq__",
          [](const SurfaceHMisalignPtStruct &self, const SurfaceHMisalignPtStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const SurfaceHMisalignPtStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
  cls.def(
         py::init<
             std::optional<bool>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>>(),
         py::arg("active") = py::none(),
         py::arg("dr") = py::none(),
         py::arg("r0") = py::none()
  )
      .def_property("active", &SurfaceHMisalignStruct::active, &SurfaceHMisalignStruct::set_active)
      .def_property(
          "dr",
          py::cpp_function(&SurfaceHMisalignStruct::dr, py::keep_alive<0, 1>()),
          &SurfaceHMisalignStruct::set_dr
      )
      .def_property(
          "r0",
          py::cpp_function(&SurfaceHMisalignStruct::r0, py::keep_alive<0, 1>()),
          &SurfaceHMisalignStruct::set_r0
      )
      .def_property_readonly(
          "pt",
          py::cpp_function(&SurfaceHMisalignStruct::pt, py::keep_alive<0, 1>())
      )

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
      .def(
          "__eq__",
          [](const SurfaceHMisalignStruct &self, const SurfaceHMisalignStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const SurfaceHMisalignStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         py::arg("x0") = py::none(),
         py::arg("y0") = py::none(),
         py::arg("z0") = py::none(),
         py::arg("dz_dx") = py::none(),
         py::arg("dz_dy") = py::none()
  )
      .def_property(
          "x0",
          &SurfaceSegmentedPtStruct::x0,
          &SurfaceSegmentedPtStruct::set_x0,
          "Position at center"
      )
      .def_property(
          "y0",
          &SurfaceSegmentedPtStruct::y0,
          &SurfaceSegmentedPtStruct::set_y0,
          "Position at center"
      )
      .def_property(
          "z0",
          &SurfaceSegmentedPtStruct::z0,
          &SurfaceSegmentedPtStruct::set_z0,
          "Position at center"
      )
      .def_property(
          "dz_dx",
          &SurfaceSegmentedPtStruct::dz_dx,
          &SurfaceSegmentedPtStruct::set_dz_dx,
          "Slope at center"
      )
      .def_property(
          "dz_dy",
          &SurfaceSegmentedPtStruct::dz_dy,
          &SurfaceSegmentedPtStruct::set_dz_dy,
          "Slope at center"
      )

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
      .def(
          "__eq__",
          [](const SurfaceSegmentedPtStruct &self, const SurfaceSegmentedPtStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const SurfaceSegmentedPtStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
  cls.def(
         py::init<
             std::optional<bool>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>>(),
         py::arg("active") = py::none(),
         py::arg("dr") = py::none(),
         py::arg("r0") = py::none()
  )
      .def_property("active", &SurfaceSegmentedStruct::active, &SurfaceSegmentedStruct::set_active)
      .def_property(
          "dr",
          py::cpp_function(&SurfaceSegmentedStruct::dr, py::keep_alive<0, 1>()),
          &SurfaceSegmentedStruct::set_dr
      )
      .def_property(
          "r0",
          py::cpp_function(&SurfaceSegmentedStruct::r0, py::keep_alive<0, 1>()),
          &SurfaceSegmentedStruct::set_r0
      )
      .def_property_readonly(
          "pt",
          py::cpp_function(&SurfaceSegmentedStruct::pt, py::keep_alive<0, 1>())
      )

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
      .def(
          "__eq__",
          [](const SurfaceSegmentedStruct &self, const SurfaceSegmentedStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const SurfaceSegmentedStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
  cls.def(
         py::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<double>>>(),
         py::arg("x0") = py::none(),
         py::arg("y0") = py::none(),
         py::arg("x1") = py::none(),
         py::arg("coef") = py::none()
  )
      .def_property("x0", &SplineStruct::x0, &SplineStruct::set_x0, "Point at start of spline")
      .def_property("y0", &SplineStruct::y0, &SplineStruct::set_y0, "Point at start of spline")
      .def_property("x1", &SplineStruct::x1, &SplineStruct::set_x1, "Point at end of spline")
      .def_property(
          "coef",
          py::cpp_function(&SplineStruct::coef, py::keep_alive<0, 1>()),
          &SplineStruct::set_coef,
          "coefficients for cubic spline"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return SplineStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = SplineStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const SplineStruct &self, const SplineStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const SplineStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<SplineStructArray1D, SplineStructAlloc1D>(
      m,
      "SplineStructArray1D",
      "SplineStructAlloc1D"
  );
  // 2D SplineStruct arrays are not used in structs/routines
  // 3D SplineStruct arrays are not used in structs/routines
}

// =============================================================================
// summation_rdt_struct
void init_summation_rdt_struct(py::module &m, py::class_<SummationRdtStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>>(),
         py::arg("h11001") = py::none(),
         py::arg("h00111") = py::none(),
         py::arg("h20001") = py::none(),
         py::arg("h00201") = py::none(),
         py::arg("h10002") = py::none(),
         py::arg("h21000") = py::none(),
         py::arg("h30000") = py::none(),
         py::arg("h10110") = py::none(),
         py::arg("h10020") = py::none(),
         py::arg("h10200") = py::none(),
         py::arg("h31000") = py::none(),
         py::arg("h40000") = py::none(),
         py::arg("h20110") = py::none(),
         py::arg("h11200") = py::none(),
         py::arg("h20020") = py::none(),
         py::arg("h20200") = py::none(),
         py::arg("h00310") = py::none(),
         py::arg("h00400") = py::none(),
         py::arg("h22000") = py::none(),
         py::arg("h00220") = py::none(),
         py::arg("h11110") = py::none()
  )
      .def_property("h11001", &SummationRdtStruct::h11001, &SummationRdtStruct::set_h11001)
      .def_property("h00111", &SummationRdtStruct::h00111, &SummationRdtStruct::set_h00111)
      .def_property("h20001", &SummationRdtStruct::h20001, &SummationRdtStruct::set_h20001)
      .def_property("h00201", &SummationRdtStruct::h00201, &SummationRdtStruct::set_h00201)
      .def_property("h10002", &SummationRdtStruct::h10002, &SummationRdtStruct::set_h10002)
      .def_property("h21000", &SummationRdtStruct::h21000, &SummationRdtStruct::set_h21000)
      .def_property("h30000", &SummationRdtStruct::h30000, &SummationRdtStruct::set_h30000)
      .def_property("h10110", &SummationRdtStruct::h10110, &SummationRdtStruct::set_h10110)
      .def_property("h10020", &SummationRdtStruct::h10020, &SummationRdtStruct::set_h10020)
      .def_property(
          "h10200",
          &SummationRdtStruct::h10200,
          &SummationRdtStruct::set_h10200,
          "2nd order in K2 moments"
      )
      .def_property("h31000", &SummationRdtStruct::h31000, &SummationRdtStruct::set_h31000)
      .def_property("h40000", &SummationRdtStruct::h40000, &SummationRdtStruct::set_h40000)
      .def_property("h20110", &SummationRdtStruct::h20110, &SummationRdtStruct::set_h20110)
      .def_property("h11200", &SummationRdtStruct::h11200, &SummationRdtStruct::set_h11200)
      .def_property("h20020", &SummationRdtStruct::h20020, &SummationRdtStruct::set_h20020)
      .def_property("h20200", &SummationRdtStruct::h20200, &SummationRdtStruct::set_h20200)
      .def_property("h00310", &SummationRdtStruct::h00310, &SummationRdtStruct::set_h00310)
      .def_property("h00400", &SummationRdtStruct::h00400, &SummationRdtStruct::set_h00400)
      .def_property("h22000", &SummationRdtStruct::h22000, &SummationRdtStruct::set_h22000)
      .def_property("h00220", &SummationRdtStruct::h00220, &SummationRdtStruct::set_h00220)
      .def_property("h11110", &SummationRdtStruct::h11110, &SummationRdtStruct::set_h11110)
      .def_static(
          "new_array1d",
          [](int sz) { return SummationRdtStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = SummationRdtStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const SummationRdtStruct &self, const SummationRdtStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const SummationRdtStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<SummationRdtStructArray1D, SummationRdtStructAlloc1D>(
      m,
      "SummationRdtStructArray1D",
      "SummationRdtStructAlloc1D"
  );
  // 2D SummationRdtStruct arrays are not used in structs/routines
  // 3D SummationRdtStruct arrays are not used in structs/routines
}