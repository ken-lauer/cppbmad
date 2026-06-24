#include "pybmad/generated/structs_p.hpp"

#include <cstdint>
#include <functional>

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
             std::optional<std::vector<double>>>(),
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
      .def_property(
          "curvature",
          py::cpp_function(&PhotonElementStruct::curvature, py::keep_alive<0, 1>()),
          &PhotonElementStruct::set_curvature
      )
      .def_property(
          "target",
          py::cpp_function(&PhotonElementStruct::target, py::keep_alive<0, 1>()),
          &PhotonElementStruct::set_target
      )
      .def_property(
          "material",
          py::cpp_function(&PhotonElementStruct::material, py::keep_alive<0, 1>()),
          &PhotonElementStruct::set_material
      )
      .def_property(
          "segmented",
          py::cpp_function(&PhotonElementStruct::segmented, py::keep_alive<0, 1>()),
          &PhotonElementStruct::set_segmented
      )
      .def_property(
          "h_misalign",
          py::cpp_function(&PhotonElementStruct::h_misalign, py::keep_alive<0, 1>()),
          &PhotonElementStruct::set_h_misalign
      )
      .def_property(
          "displacement",
          py::cpp_function(&PhotonElementStruct::displacement, py::keep_alive<0, 1>()),
          &PhotonElementStruct::set_displacement
      )
      .def_property(
          "pixel",
          py::cpp_function(&PhotonElementStruct::pixel, py::keep_alive<0, 1>()),
          &PhotonElementStruct::set_pixel
      )
      .def_property(
          "reflectivity_table_type",
          &PhotonElementStruct::reflectivity_table_type,
          &PhotonElementStruct::set_reflectivity_table_type
      )
      .def_property(
          "reflectivity_table_sigma",
          py::cpp_function(&PhotonElementStruct::reflectivity_table_sigma, py::keep_alive<0, 1>()),
          &PhotonElementStruct::set_reflectivity_table_sigma,
          "If polarization is ignored use sigma table."
      )
      .def_property(
          "reflectivity_table_pi",
          py::cpp_function(&PhotonElementStruct::reflectivity_table_pi, py::keep_alive<0, 1>()),
          &PhotonElementStruct::set_reflectivity_table_pi
      )
      .def_property_readonly(
          "init_energy_prob",
          py::cpp_function(&PhotonElementStruct::init_energy_prob, py::keep_alive<0, 1>()),
          "Initial energy probability density"
      )
      .def_property(
          "integrated_init_energy_prob",
          py::cpp_function(
              &PhotonElementStruct::integrated_init_energy_prob,
              py::keep_alive<0, 1>()
          ),
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
      .def(
          "__eq__",
          [](const PhotonElementStruct &self, const PhotonElementStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const PhotonElementStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
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
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>>(),
         py::arg("f0_m1") = py::none(),
         py::arg("f0_m2") = py::none(),
         py::arg("f_0") = py::none(),
         py::arg("f_h") = py::none(),
         py::arg("f_hbar") = py::none(),
         py::arg("f_hkl") = py::none(),
         py::arg("h_norm") = py::none(),
         py::arg("l_ref") = py::none()
  )
      .def_property(
          "f0_m1",
          &PhotonMaterialStruct::f0_m1,
          &PhotonMaterialStruct::set_f0_m1,
          "For multilayer_mirror only."
      )
      .def_property(
          "f0_m2",
          &PhotonMaterialStruct::f0_m2,
          &PhotonMaterialStruct::set_f0_m2,
          "For multilayer_mirror only."
      )
      .def_property("f_0", &PhotonMaterialStruct::f_0, &PhotonMaterialStruct::set_f_0)
      .def_property(
          "f_h",
          &PhotonMaterialStruct::f_h,
          &PhotonMaterialStruct::set_f_h,
          "Structure factor for H direction."
      )
      .def_property(
          "f_hbar",
          &PhotonMaterialStruct::f_hbar,
          &PhotonMaterialStruct::set_f_hbar,
          "Structure factor for -H direction."
      )
      .def_property(
          "f_hkl",
          &PhotonMaterialStruct::f_hkl,
          &PhotonMaterialStruct::set_f_hkl,
          "= sqrt(f_h * f_hbar)"
      )
      .def_property(
          "h_norm",
          py::cpp_function(&PhotonMaterialStruct::h_norm, py::keep_alive<0, 1>()),
          &PhotonMaterialStruct::set_h_norm,
          "Normalized H vector for crystals."
      )
      .def_property(
          "l_ref",
          py::cpp_function(&PhotonMaterialStruct::l_ref, py::keep_alive<0, 1>()),
          &PhotonMaterialStruct::set_l_ref,
          "Crystal reference orbit displacement vector in element coords."
      )

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
      .def(
          "__eq__",
          [](const PhotonMaterialStruct &self, const PhotonMaterialStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const PhotonMaterialStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
      .def_property(
          "name",
          &PhotonReflectSurfaceStruct::name,
          &PhotonReflectSurfaceStruct::set_name
      )
      .def_property(
          "description",
          &PhotonReflectSurfaceStruct::description,
          &PhotonReflectSurfaceStruct::set_description,
          "Descriptive name"
      )
      .def_property(
          "reflectivity_file",
          &PhotonReflectSurfaceStruct::reflectivity_file,
          &PhotonReflectSurfaceStruct::set_reflectivity_file
      )
      .def_property_readonly(
          "table",
          py::cpp_function(&PhotonReflectSurfaceStruct::table, py::keep_alive<0, 1>())
      )
      .def_property(
          "surface_roughness_rms",
          &PhotonReflectSurfaceStruct::surface_roughness_rms,
          &PhotonReflectSurfaceStruct::set_surface_roughness_rms,
          "sigma in Dugan's notation"
      )
      .def_property(
          "roughness_correlation_len",
          &PhotonReflectSurfaceStruct::roughness_correlation_len,
          &PhotonReflectSurfaceStruct::set_roughness_correlation_len,
          "T in Dugan's notation"
      )
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
      .def(
          "__eq__",
          [](const PhotonReflectSurfaceStruct &self, const PhotonReflectSurfaceStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const PhotonReflectSurfaceStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
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
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>>(),
         py::arg("angle") = py::none(),
         py::arg("energy") = py::none(),
         py::arg("p_reflect") = py::none(),
         py::arg("max_energy") = py::none(),
         py::arg("p_reflect_scratch") = py::none(),
         py::arg("bragg_angle") = py::none()
  )
      .def_property(
          "angle",
          py::cpp_function(&PhotonReflectTableStruct::angle, py::keep_alive<0, 1>()),
          &PhotonReflectTableStruct::set_angle,
          "Vector of angle values for %p_reflect"
      )
      .def_property(
          "energy",
          py::cpp_function(&PhotonReflectTableStruct::energy, py::keep_alive<0, 1>()),
          &PhotonReflectTableStruct::set_energy,
          "Vector of energy values for %p_reflect"
      )
      .def_property_readonly(
          "int1",
          py::cpp_function(&PhotonReflectTableStruct::int1, py::keep_alive<0, 1>())
      )
      .def_property(
          "p_reflect",
          py::cpp_function(&PhotonReflectTableStruct::p_reflect, py::keep_alive<0, 1>()),
          &PhotonReflectTableStruct::set_p_reflect,
          "(angle, ev) probability. Log used for smooth surface reflection"
      )
      .def_property(
          "max_energy",
          &PhotonReflectTableStruct::max_energy,
          &PhotonReflectTableStruct::set_max_energy,
          "maximum energy for this table"
      )
      .def_property(
          "p_reflect_scratch",
          py::cpp_function(&PhotonReflectTableStruct::p_reflect_scratch, py::keep_alive<0, 1>()),
          &PhotonReflectTableStruct::set_p_reflect_scratch,
          "Scratch space"
      )
      .def_property(
          "bragg_angle",
          py::cpp_function(&PhotonReflectTableStruct::bragg_angle, py::keep_alive<0, 1>()),
          &PhotonReflectTableStruct::set_bragg_angle,
          "Bragg angle at energy values."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return PhotonReflectTableStructAlloc1D(sz); },
          py::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = PhotonReflectTableStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          py::arg("lbound"),
          py::arg("ubound")
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
      .def(
          "__eq__",
          [](const PhotonReflectTableStruct &self, const PhotonReflectTableStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const PhotonReflectTableStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<PhotonReflectTableStructArray1D, PhotonReflectTableStructAlloc1D>(
      m,
      "PhotonReflectTableStructArray1D",
      "PhotonReflectTableStructAlloc1D"
  );
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
      .def_property(
          "type",
          &PhotonTargetStruct::type,
          &PhotonTargetStruct::set_type,
          "or rectangular$"
      )
      .def_property("n_corner", &PhotonTargetStruct::n_corner, &PhotonTargetStruct::set_n_corner)
      .def_property(
          "ele_loc",
          py::cpp_function(&PhotonTargetStruct::ele_loc, py::keep_alive<0, 1>()),
          &PhotonTargetStruct::set_ele_loc
      )
      .def_property_readonly(
          "corner",
          py::cpp_function(&PhotonTargetStruct::corner, py::keep_alive<0, 1>())
      )
      .def_property(
          "center",
          py::cpp_function(&PhotonTargetStruct::center, py::keep_alive<0, 1>()),
          &PhotonTargetStruct::set_center
      )

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
      .def(
          "__eq__",
          [](const PhotonTargetStruct &self, const PhotonTargetStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const PhotonTargetStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
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
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<int64_t>,
             std::optional<int64_t>,
             std::optional<int64_t>>(),
         py::arg("dr") = py::none(),
         py::arg("r0") = py::none(),
         py::arg("n_track_tot") = py::none(),
         py::arg("n_hit_detec") = py::none(),
         py::arg("n_hit_pixel") = py::none()
  )
      .def_property(
          "dr",
          py::cpp_function(&PixelDetecStruct::dr, py::keep_alive<0, 1>()),
          &PixelDetecStruct::set_dr
      )
      .def_property(
          "r0",
          py::cpp_function(&PixelDetecStruct::r0, py::keep_alive<0, 1>()),
          &PixelDetecStruct::set_r0
      )
      .def_property(
          "n_track_tot",
          &PixelDetecStruct::n_track_tot,
          &PixelDetecStruct::set_n_track_tot,
          "How many photons were launched from source element."
      )
      .def_property(
          "n_hit_detec",
          &PixelDetecStruct::n_hit_detec,
          &PixelDetecStruct::set_n_hit_detec,
          "How many photons hit the detector."
      )
      .def_property(
          "n_hit_pixel",
          &PixelDetecStruct::n_hit_pixel,
          &PixelDetecStruct::set_n_hit_pixel,
          "How many photons hit the pixel grid of the detector."
      )
      .def_property_readonly(
          "pt",
          py::cpp_function(&PixelDetecStruct::pt, py::keep_alive<0, 1>()),
          "Grid of pixels"
      )

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
      .def(
          "__eq__",
          [](const PixelDetecStruct &self, const PixelDetecStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const PixelDetecStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
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
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>>(),
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
      .def_property("n_photon", &PixelPtStruct::n_photon, &PixelPtStruct::set_n_photon)
      .def_property("E_x", &PixelPtStruct::E_x, &PixelPtStruct::set_E_x)
      .def_property("E_y", &PixelPtStruct::E_y, &PixelPtStruct::set_E_y)
      .def_property("intensity_x", &PixelPtStruct::intensity_x, &PixelPtStruct::set_intensity_x)
      .def_property("intensity_y", &PixelPtStruct::intensity_y, &PixelPtStruct::set_intensity_y)
      .def_property("intensity", &PixelPtStruct::intensity, &PixelPtStruct::set_intensity)
      .def_property(
          "orbit",
          py::cpp_function(&PixelPtStruct::orbit, py::keep_alive<0, 1>()),
          &PixelPtStruct::set_orbit,
          "x, Vx/c, y, Vy/c, dummy, E - E_ref."
      )
      .def_property(
          "orbit_rms",
          py::cpp_function(&PixelPtStruct::orbit_rms, py::keep_alive<0, 1>()),
          &PixelPtStruct::set_orbit_rms,
          "RMS statistics."
      )
      .def_property(
          "init_orbit",
          py::cpp_function(&PixelPtStruct::init_orbit, py::keep_alive<0, 1>()),
          &PixelPtStruct::set_init_orbit,
          "Initial orbit at start of lattice statistics."
      )
      .def_property(
          "init_orbit_rms",
          py::cpp_function(&PixelPtStruct::init_orbit_rms, py::keep_alive<0, 1>()),
          &PixelPtStruct::set_init_orbit_rms,
          "Initial orbit at start of lattice RMS statistics."
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
      .def(
          "__eq__",
          [](const PixelPtStruct &self, const PixelPtStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const PixelPtStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
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
      .def_property(
          "who",
          &PreTrackerStruct::who,
          &PreTrackerStruct::set_who,
          "Can be opal$, or impactt$"
      )
      .def_property(
          "ix_ele_start",
          &PreTrackerStruct::ix_ele_start,
          &PreTrackerStruct::set_ix_ele_start
      )
      .def_property("ix_ele_end", &PreTrackerStruct::ix_ele_end, &PreTrackerStruct::set_ix_ele_end)
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
      .def(
          "__eq__",
          [](const PreTrackerStruct &self, const PreTrackerStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const PreTrackerStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
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
             std::optional<std::vector<double>>,
             std::optional<bool>>(),
         py::arg("ele_origin") = py::none(),
         py::arg("orb0") = py::none(),
         py::arg("valid_map") = py::none()
  )
      .def_property(
          "ele_origin",
          py::cpp_function(&PtcNormalFormStruct::ele_origin, py::keep_alive<0, 1>()),
          &PtcNormalFormStruct::set_ele_origin,
          "Element at which the on-turn map was created."
      )
      .def_property(
          "orb0",
          py::cpp_function(&PtcNormalFormStruct::orb0, py::keep_alive<0, 1>()),
          &PtcNormalFormStruct::set_orb0,
          "Closed orbit at element."
      )
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
      .def(
          "__eq__",
          [](const PtcNormalFormStruct &self, const PtcNormalFormStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          py::is_operator()
      )
      .def(
          "__hash__",
          [](const PtcNormalFormStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D PtcNormalFormStruct arrays are not used in structs/routines
  // 2D PtcNormalFormStruct arrays are not used in structs/routines
  // 3D PtcNormalFormStruct arrays are not used in structs/routines
}