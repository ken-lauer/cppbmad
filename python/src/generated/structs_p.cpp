#include "pybmad/generated/structs_p.hpp"

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// photon_element_struct
void init_photon_element_struct(py::module &m, py::class_<PhotonElementStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const SurfaceCurvatureStruct>,
             optional_ref<const PhotonTargetStruct>,
             optional_ref<const PhotonMaterialStruct>,
             optional_ref<const SurfaceSegmentedStruct>,
             optional_ref<const SurfaceHMisalignStruct>,
             optional_ref<const SurfaceDisplacementStruct>,
             optional_ref<const PixelDetecStruct>,
             std::optional<int>,
             optional_ref<const PhotonReflectTableStruct>,
             optional_ref<const PhotonReflectTableStruct>,
             optional_ref<const std::vector<double>>>(),
         py::arg("curvature") = py::none(),
         py::arg("target") = py::none(),
         py::arg("material") = py::none(),
         py::arg("segmented") = py::none(),
         py::arg("h_misalign") = py::none(),
         py::arg("displacement") = py::none(),
         py::arg("pixel") = py::none(),
         py::arg("reflectivity_table_type") = py::none(),
         py::arg("reflectivity_table_sigma") = py::none(),
         py::arg("reflectivity_table_pi") = py::none(),
         py::arg("integrated_init_energy_prob") = py::none()
  )
      // PhotonElementStruct.curvature (0D_NOT_type -
      .def_property(
          "curvature",
          &PhotonElementStruct::curvature,
          &PhotonElementStruct::set_curvature
      )
      // PhotonElementStruct.target (0D_NOT_type -
      .def_property("target", &PhotonElementStruct::target, &PhotonElementStruct::set_target)
      // PhotonElementStruct.material (0D_NOT_type -
      .def_property("material", &PhotonElementStruct::material, &PhotonElementStruct::set_material)
      // PhotonElementStruct.segmented (0D_NOT_type -
      .def_property(
          "segmented",
          &PhotonElementStruct::segmented,
          &PhotonElementStruct::set_segmented
      )
      // PhotonElementStruct.h_misalign (0D_NOT_type -
      .def_property(
          "h_misalign",
          &PhotonElementStruct::h_misalign,
          &PhotonElementStruct::set_h_misalign
      )
      // PhotonElementStruct.displacement (0D_NOT_type -
      .def_property(
          "displacement",
          &PhotonElementStruct::displacement,
          &PhotonElementStruct::set_displacement
      )
      // PhotonElementStruct.pixel (0D_NOT_type -
      .def_property("pixel", &PhotonElementStruct::pixel, &PhotonElementStruct::set_pixel)
      // PhotonElementStruct.reflectivity_table_type (0D_NOT_integer -
      .def_property(
          "reflectivity_table_type",
          &PhotonElementStruct::reflectivity_table_type,
          &PhotonElementStruct::set_reflectivity_table_type
      )
      // PhotonElementStruct.reflectivity_table_sigma (0D_NOT_type - If polarization is ignored use
      // sigma table.
      .def_property(
          "reflectivity_table_sigma",
          &PhotonElementStruct::reflectivity_table_sigma,
          &PhotonElementStruct::set_reflectivity_table_sigma
      )
      // PhotonElementStruct.reflectivity_table_pi (0D_NOT_type -
      .def_property(
          "reflectivity_table_pi",
          &PhotonElementStruct::reflectivity_table_pi,
          &PhotonElementStruct::set_reflectivity_table_pi
      )
      // PhotonElementStruct.init_energy_prob (1D_ALLOC_type - Initial energy probability density
      .def_property_readonly("init_energy_prob", &PhotonElementStruct::init_energy_prob)
      // PhotonElementStruct.integrated_init_energy_prob (1D_ALLOC_real -
      .def_property(
          "integrated_init_energy_prob",
          &PhotonElementStruct::integrated_init_energy_prob,
          &PhotonElementStruct::set_integrated_init_energy_prob
      )

      .def("__repr__", [](const PhotonElementStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PhotonElementStruct &self) {
            return PhotonElementStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const PhotonElementStruct &self, py::dict &memo) { return PhotonElementStruct(self); }
      )

      ;

  // 1D PhotonElementStruct arrays are not used in structs/routines
  // 2D PhotonElementStruct arrays are not used in structs/routines
  // 3D PhotonElementStruct arrays are not used in structs/routines
}

// =============================================================================
// photon_material_struct
void init_photon_material_struct(py::module &m, py::class_<PhotonMaterialStruct> &cls) {
  cls.def(
         py::init<
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<double>>>(),
         py::arg("f0_m1") = py::none(),
         py::arg("f0_m2") = py::none(),
         py::arg("f_0") = py::none(),
         py::arg("f_h") = py::none(),
         py::arg("f_hbar") = py::none(),
         py::arg("f_hkl") = py::none(),
         py::arg("h_norm") = py::none(),
         py::arg("l_ref") = py::none()
  )
      // PhotonMaterialStruct.f0_m1 (0D_NOT_complex - For multilayer_mirror only.
      .def_property("f0_m1", &PhotonMaterialStruct::f0_m1, &PhotonMaterialStruct::set_f0_m1)
      // PhotonMaterialStruct.f0_m2 (0D_NOT_complex - For multilayer_mirror only.
      .def_property("f0_m2", &PhotonMaterialStruct::f0_m2, &PhotonMaterialStruct::set_f0_m2)
      // PhotonMaterialStruct.f_0 (0D_NOT_complex -
      .def_property("f_0", &PhotonMaterialStruct::f_0, &PhotonMaterialStruct::set_f_0)
      // PhotonMaterialStruct.f_h (0D_NOT_complex - Structure factor for H direction.
      .def_property("f_h", &PhotonMaterialStruct::f_h, &PhotonMaterialStruct::set_f_h)
      // PhotonMaterialStruct.f_hbar (0D_NOT_complex - Structure factor for -H direction.
      .def_property("f_hbar", &PhotonMaterialStruct::f_hbar, &PhotonMaterialStruct::set_f_hbar)
      // PhotonMaterialStruct.f_hkl (0D_NOT_complex - = sqrt(f_h * f_hbar)
      .def_property("f_hkl", &PhotonMaterialStruct::f_hkl, &PhotonMaterialStruct::set_f_hkl)
      // PhotonMaterialStruct.h_norm (1D_NOT_real - Normalized H vector for crystals.
      .def_property("h_norm", &PhotonMaterialStruct::h_norm, &PhotonMaterialStruct::set_h_norm)
      // PhotonMaterialStruct.l_ref (1D_NOT_real - Crystal reference orbit displacement vector in
      // element coords.
      .def_property("l_ref", &PhotonMaterialStruct::l_ref, &PhotonMaterialStruct::set_l_ref)

      .def("__repr__", [](const PhotonMaterialStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PhotonMaterialStruct &self) {
            return PhotonMaterialStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const PhotonMaterialStruct &self, py::dict &memo) {
            return PhotonMaterialStruct(self);
          }
      )

      ;

  // 1D PhotonMaterialStruct arrays are not used in structs/routines
  // 2D PhotonMaterialStruct arrays are not used in structs/routines
  // 3D PhotonMaterialStruct arrays are not used in structs/routines
}

// =============================================================================
// photon_reflect_surface_struct
void init_photon_reflect_surface_struct(
    py::module &m,
    py::class_<PhotonReflectSurfaceStruct> &cls
) {
  cls.def(
         py::init<
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>>(),
         py::arg("name") = py::none(),
         py::arg("description") = py::none(),
         py::arg("reflectivity_file") = py::none(),
         py::arg("surface_roughness_rms") = py::none(),
         py::arg("roughness_correlation_len") = py::none(),
         py::arg("ix_surface") = py::none()
  )
      // PhotonReflectSurfaceStruct.name (0D_NOT_character -
      .def_property(
          "name",
          &PhotonReflectSurfaceStruct::name,
          &PhotonReflectSurfaceStruct::set_name
      )
      // PhotonReflectSurfaceStruct.description (0D_NOT_character - Descriptive name
      .def_property(
          "description",
          &PhotonReflectSurfaceStruct::description,
          &PhotonReflectSurfaceStruct::set_description
      )
      // PhotonReflectSurfaceStruct.reflectivity_file (0D_NOT_character -
      .def_property(
          "reflectivity_file",
          &PhotonReflectSurfaceStruct::reflectivity_file,
          &PhotonReflectSurfaceStruct::set_reflectivity_file
      )
      // PhotonReflectSurfaceStruct.table (1D_ALLOC_type -
      .def_property_readonly("table", &PhotonReflectSurfaceStruct::table)
      // PhotonReflectSurfaceStruct.surface_roughness_rms (0D_NOT_real - sigma in Dugan's notation
      .def_property(
          "surface_roughness_rms",
          &PhotonReflectSurfaceStruct::surface_roughness_rms,
          &PhotonReflectSurfaceStruct::set_surface_roughness_rms
      )
      // PhotonReflectSurfaceStruct.roughness_correlation_len (0D_NOT_real - T in Dugan's notation
      .def_property(
          "roughness_correlation_len",
          &PhotonReflectSurfaceStruct::roughness_correlation_len,
          &PhotonReflectSurfaceStruct::set_roughness_correlation_len
      )
      // PhotonReflectSurfaceStruct.ix_surface (0D_NOT_integer -
      .def_property(
          "ix_surface",
          &PhotonReflectSurfaceStruct::ix_surface,
          &PhotonReflectSurfaceStruct::set_ix_surface
      )

      .def("__repr__", [](const PhotonReflectSurfaceStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PhotonReflectSurfaceStruct &self) {
            return PhotonReflectSurfaceStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const PhotonReflectSurfaceStruct &self, py::dict &memo) {
            return PhotonReflectSurfaceStruct(self);
          }
      )

      ;

  // 1D PhotonReflectSurfaceStruct arrays are not used in structs/routines
  // 2D PhotonReflectSurfaceStruct arrays are not used in structs/routines
  // 3D PhotonReflectSurfaceStruct arrays are not used in structs/routines
}

// =============================================================================
// photon_reflect_table_struct
void init_photon_reflect_table_struct(py::module &m, py::class_<PhotonReflectTableStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<std::vector<double>>>,
             std::optional<double>,
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<double>>>(),
         py::arg("angle") = py::none(),
         py::arg("energy") = py::none(),
         py::arg("p_reflect") = py::none(),
         py::arg("max_energy") = py::none(),
         py::arg("p_reflect_scratch") = py::none(),
         py::arg("bragg_angle") = py::none()
  )
      // PhotonReflectTableStruct.angle (1D_ALLOC_real - Vector of angle values for %p_reflect
      .def_property("angle", &PhotonReflectTableStruct::angle, &PhotonReflectTableStruct::set_angle)
      // PhotonReflectTableStruct.energy (1D_ALLOC_real - Vector of energy values for %p_reflect
      .def_property(
          "energy",
          &PhotonReflectTableStruct::energy,
          &PhotonReflectTableStruct::set_energy
      )
      // PhotonReflectTableStruct.int1 (1D_ALLOC_type -
      .def_property_readonly("int1", &PhotonReflectTableStruct::int1)
      // PhotonReflectTableStruct.p_reflect (2D_ALLOC_real - (angle, ev) probability. Log used for
      // smooth surface reflection
      .def_property(
          "p_reflect",
          &PhotonReflectTableStruct::p_reflect,
          &PhotonReflectTableStruct::set_p_reflect
      )
      // PhotonReflectTableStruct.max_energy (0D_NOT_real - maximum energy for this table
      .def_property(
          "max_energy",
          &PhotonReflectTableStruct::max_energy,
          &PhotonReflectTableStruct::set_max_energy
      )
      // PhotonReflectTableStruct.p_reflect_scratch (1D_ALLOC_real - Scratch space
      .def_property(
          "p_reflect_scratch",
          &PhotonReflectTableStruct::p_reflect_scratch,
          &PhotonReflectTableStruct::set_p_reflect_scratch
      )
      // PhotonReflectTableStruct.bragg_angle (1D_ALLOC_real - Bragg angle at energy values.
      .def_property(
          "bragg_angle",
          &PhotonReflectTableStruct::bragg_angle,
          &PhotonReflectTableStruct::set_bragg_angle
      )
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return PhotonReflectTableStructAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1
      )

      .def("__repr__", [](const PhotonReflectTableStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PhotonReflectTableStruct &self) {
            return PhotonReflectTableStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const PhotonReflectTableStruct &self, py::dict &memo) {
            return PhotonReflectTableStruct(self);
          }
      )

      ;

  bind_FTypeArrayND<PhotonReflectTableStructArray1D>(m, "PhotonReflectTableStructArray1D");
  bind_FTypeAlloc1D<PhotonReflectTableStructAlloc1D>(m, "PhotonReflectTableStructAlloc1D");
  // 2D PhotonReflectTableStruct arrays are not used in structs/routines
  // 3D PhotonReflectTableStruct arrays are not used in structs/routines
}

// =============================================================================
// photon_target_struct
void init_photon_target_struct(py::module &m, py::class_<PhotonTargetStruct> &cls) {
  cls.def(
         py::init<
             std::optional<int>,
             std::optional<int>,
             optional_ref<const LatEleLocStruct>,
             optional_ref<const TargetPointStruct>>(),
         py::arg("type") = py::none(),
         py::arg("n_corner") = py::none(),
         py::arg("ele_loc") = py::none(),
         py::arg("center") = py::none()
  )
      // PhotonTargetStruct.type (0D_NOT_integer - or rectangular$
      .def_property("type", &PhotonTargetStruct::type, &PhotonTargetStruct::set_type)
      // PhotonTargetStruct.n_corner (0D_NOT_integer -
      .def_property("n_corner", &PhotonTargetStruct::n_corner, &PhotonTargetStruct::set_n_corner)
      // PhotonTargetStruct.ele_loc (0D_NOT_type -
      .def_property("ele_loc", &PhotonTargetStruct::ele_loc, &PhotonTargetStruct::set_ele_loc)
      // PhotonTargetStruct.corner (1D_NOT_type -
      .def_property_readonly("corner", &PhotonTargetStruct::corner)
      // PhotonTargetStruct.center (0D_NOT_type -
      .def_property("center", &PhotonTargetStruct::center, &PhotonTargetStruct::set_center)

      .def("__repr__", [](const PhotonTargetStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PhotonTargetStruct &self) {
            return PhotonTargetStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const PhotonTargetStruct &self, py::dict &memo) { return PhotonTargetStruct(self); }
      )

      ;

  // 1D PhotonTargetStruct arrays are not used in structs/routines
  // 2D PhotonTargetStruct arrays are not used in structs/routines
  // 3D PhotonTargetStruct arrays are not used in structs/routines
}

// =============================================================================
// pixel_detec_struct
void init_pixel_detec_struct(py::module &m, py::class_<PixelDetecStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<double>>,
             std::optional<int64_t>,
             std::optional<int64_t>,
             std::optional<int64_t>>(),
         py::arg("dr") = py::none(),
         py::arg("r0") = py::none(),
         py::arg("n_track_tot") = py::none(),
         py::arg("n_hit_detec") = py::none(),
         py::arg("n_hit_pixel") = py::none()
  )
      // PixelDetecStruct.dr (1D_NOT_real -
      .def_property("dr", &PixelDetecStruct::dr, &PixelDetecStruct::set_dr)
      // PixelDetecStruct.r0 (1D_NOT_real -
      .def_property("r0", &PixelDetecStruct::r0, &PixelDetecStruct::set_r0)
      // PixelDetecStruct.n_track_tot (0D_NOT_integer8 - How many photons were launched from source
      // element.
      .def_property(
          "n_track_tot",
          &PixelDetecStruct::n_track_tot,
          &PixelDetecStruct::set_n_track_tot
      )
      // PixelDetecStruct.n_hit_detec (0D_NOT_integer8 - How many photons hit the detector.
      .def_property(
          "n_hit_detec",
          &PixelDetecStruct::n_hit_detec,
          &PixelDetecStruct::set_n_hit_detec
      )
      // PixelDetecStruct.n_hit_pixel (0D_NOT_integer8 - How many photons hit the pixel grid of the
      // detector.
      .def_property(
          "n_hit_pixel",
          &PixelDetecStruct::n_hit_pixel,
          &PixelDetecStruct::set_n_hit_pixel
      )
      // PixelDetecStruct.pt (2D_ALLOC_type - Grid of pixels
      .def_property_readonly("pt", &PixelDetecStruct::pt)

      .def("__repr__", [](const PixelDetecStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PixelDetecStruct &self) {
            return PixelDetecStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const PixelDetecStruct &self, py::dict &memo) { return PixelDetecStruct(self); }
      )

      ;

  // 1D PixelDetecStruct arrays are not used in structs/routines
  // 2D PixelDetecStruct arrays are not used in structs/routines
  // 3D PixelDetecStruct arrays are not used in structs/routines
}

// =============================================================================
// pixel_pt_struct
void init_pixel_pt_struct(py::module &m, py::class_<PixelPtStruct> &cls) {
  cls.def(
         py::init<
             std::optional<int64_t>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<double>>,
             optional_ref<const std::vector<double>>>(),
         py::arg("n_photon") = py::none(),
         py::arg("E_x") = py::none(),
         py::arg("E_y") = py::none(),
         py::arg("intensity_x") = py::none(),
         py::arg("intensity_y") = py::none(),
         py::arg("intensity") = py::none(),
         py::arg("orbit") = py::none(),
         py::arg("orbit_rms") = py::none(),
         py::arg("init_orbit") = py::none(),
         py::arg("init_orbit_rms") = py::none()
  )
      // PixelPtStruct.n_photon (0D_NOT_integer8 -
      .def_property("n_photon", &PixelPtStruct::n_photon, &PixelPtStruct::set_n_photon)
      // PixelPtStruct.E_x (0D_NOT_complex -
      .def_property("E_x", &PixelPtStruct::E_x, &PixelPtStruct::set_E_x)
      // PixelPtStruct.E_y (0D_NOT_complex -
      .def_property("E_y", &PixelPtStruct::E_y, &PixelPtStruct::set_E_y)
      // PixelPtStruct.intensity_x (0D_NOT_real -
      .def_property("intensity_x", &PixelPtStruct::intensity_x, &PixelPtStruct::set_intensity_x)
      // PixelPtStruct.intensity_y (0D_NOT_real -
      .def_property("intensity_y", &PixelPtStruct::intensity_y, &PixelPtStruct::set_intensity_y)
      // PixelPtStruct.intensity (0D_NOT_real -
      .def_property("intensity", &PixelPtStruct::intensity, &PixelPtStruct::set_intensity)
      // PixelPtStruct.orbit (1D_NOT_real - x, Vx/c, y, Vy/c, dummy, E - E_ref.
      .def_property("orbit", &PixelPtStruct::orbit, &PixelPtStruct::set_orbit)
      // PixelPtStruct.orbit_rms (1D_NOT_real - RMS statistics.
      .def_property("orbit_rms", &PixelPtStruct::orbit_rms, &PixelPtStruct::set_orbit_rms)
      // PixelPtStruct.init_orbit (1D_NOT_real - Initial orbit at start of lattice statistics.
      .def_property("init_orbit", &PixelPtStruct::init_orbit, &PixelPtStruct::set_init_orbit)
      // PixelPtStruct.init_orbit_rms (1D_NOT_real - Initial orbit at start of lattice RMS
      // statistics.
      .def_property(
          "init_orbit_rms",
          &PixelPtStruct::init_orbit_rms,
          &PixelPtStruct::set_init_orbit_rms
      )

      .def("__repr__", [](const PixelPtStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PixelPtStruct &self) {
            return PixelPtStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const PixelPtStruct &self, py::dict &memo) { return PixelPtStruct(self); }
      )

      ;

  // 1D PixelPtStruct arrays are not used in structs/routines
  bind_FTypeArrayND<PixelPtStructArray2D>(m, "PixelPtStructArray2D");
  // 3D PixelPtStruct arrays are not used in structs/routines
}

// =============================================================================
// pre_tracker_struct
void init_pre_tracker_struct(py::module &m, py::class_<PreTrackerStruct> &cls) {
  cls.def(
         py::init<
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<std::string>>(),
         py::arg("who") = py::none(),
         py::arg("ix_ele_start") = py::none(),
         py::arg("ix_ele_end") = py::none(),
         py::arg("input_file") = py::none()
  )
      // PreTrackerStruct.who (0D_NOT_integer - Can be opal$, or impactt$
      .def_property("who", &PreTrackerStruct::who, &PreTrackerStruct::set_who)
      // PreTrackerStruct.ix_ele_start (0D_NOT_integer -
      .def_property(
          "ix_ele_start",
          &PreTrackerStruct::ix_ele_start,
          &PreTrackerStruct::set_ix_ele_start
      )
      // PreTrackerStruct.ix_ele_end (0D_NOT_integer -
      .def_property("ix_ele_end", &PreTrackerStruct::ix_ele_end, &PreTrackerStruct::set_ix_ele_end)
      // PreTrackerStruct.input_file (0D_NOT_character -
      .def_property("input_file", &PreTrackerStruct::input_file, &PreTrackerStruct::set_input_file)

      .def("__repr__", [](const PreTrackerStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PreTrackerStruct &self) {
            return PreTrackerStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const PreTrackerStruct &self, py::dict &memo) { return PreTrackerStruct(self); }
      )

      ;

  // 1D PreTrackerStruct arrays are not used in structs/routines
  // 2D PreTrackerStruct arrays are not used in structs/routines
  // 3D PreTrackerStruct arrays are not used in structs/routines
}

// =============================================================================
// ptc_normal_form_struct
void init_ptc_normal_form_struct(py::module &m, py::class_<PtcNormalFormStruct> &cls) {
  cls.def(
         py::init<
             optional_ref<const EleStruct>,
             optional_ref<const std::vector<double>>,
             std::optional<bool>>(),
         py::arg("ele_origin") = py::none(),
         py::arg("orb0") = py::none(),
         py::arg("valid_map") = py::none()
  )
      // PtcNormalFormStruct.ele_origin (0D_PTR_type - Element at which the on-turn map was created.
      .def_property(
          "ele_origin",
          &PtcNormalFormStruct::ele_origin,
          &PtcNormalFormStruct::set_ele_origin
      )
      // PtcNormalFormStruct.orb0 (1D_NOT_real - Closed orbit at element.
      .def_property("orb0", &PtcNormalFormStruct::orb0, &PtcNormalFormStruct::set_orb0)
      // PtcNormalFormStruct.valid_map (0D_NOT_logical -
      .def_property(
          "valid_map",
          &PtcNormalFormStruct::valid_map,
          &PtcNormalFormStruct::set_valid_map
      )

      .def("__repr__", [](const PtcNormalFormStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PtcNormalFormStruct &self) {
            return PtcNormalFormStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const PtcNormalFormStruct &self, py::dict &memo) { return PtcNormalFormStruct(self); }
      )

      ;

  // 1D PtcNormalFormStruct arrays are not used in structs/routines
  // 2D PtcNormalFormStruct arrays are not used in structs/routines
  // 3D PtcNormalFormStruct arrays are not used in structs/routines
}