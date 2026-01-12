#include "pybmad/generated/structs_p.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// photon_element_struct
void init_photon_element_struct(
    py::module& m,
    py::class_<PhotonElementProxy>& cls) {
  cls.def(py::init<>())
      // PhotonElementProxy.curvature (0D_NOT_type -
      .def_property(
          "curvature",
          &PhotonElementProxy::curvature,
          &PhotonElementProxy::set_curvature)
      // PhotonElementProxy.target (0D_NOT_type -
      .def_property(
          "target",
          &PhotonElementProxy::target,
          &PhotonElementProxy::set_target)
      // PhotonElementProxy.material (0D_NOT_type -
      .def_property(
          "material",
          &PhotonElementProxy::material,
          &PhotonElementProxy::set_material)
      // PhotonElementProxy.segmented (0D_NOT_type -
      .def_property(
          "segmented",
          &PhotonElementProxy::segmented,
          &PhotonElementProxy::set_segmented)
      // PhotonElementProxy.h_misalign (0D_NOT_type -
      .def_property(
          "h_misalign",
          &PhotonElementProxy::h_misalign,
          &PhotonElementProxy::set_h_misalign)
      // PhotonElementProxy.displacement (0D_NOT_type -
      .def_property(
          "displacement",
          &PhotonElementProxy::displacement,
          &PhotonElementProxy::set_displacement)
      // PhotonElementProxy.pixel (0D_NOT_type -
      .def_property(
          "pixel", &PhotonElementProxy::pixel, &PhotonElementProxy::set_pixel)
      // PhotonElementProxy.reflectivity_table_type (0D_NOT_integer -
      .def_property(
          "reflectivity_table_type",
          &PhotonElementProxy::reflectivity_table_type,
          &PhotonElementProxy::set_reflectivity_table_type)
      // PhotonElementProxy.reflectivity_table_sigma (0D_NOT_type - If polarization is ignored use sigma table.
      .def_property(
          "reflectivity_table_sigma",
          &PhotonElementProxy::reflectivity_table_sigma,
          &PhotonElementProxy::set_reflectivity_table_sigma)
      // PhotonElementProxy.reflectivity_table_pi (0D_NOT_type -
      .def_property(
          "reflectivity_table_pi",
          &PhotonElementProxy::reflectivity_table_pi,
          &PhotonElementProxy::set_reflectivity_table_pi)
      // PhotonElementProxy.init_energy_prob (1D_ALLOC_type - Initial energy probability density
      .def_property_readonly(
          "init_energy_prob", &PhotonElementProxy::init_energy_prob)
      // PhotonElementProxy.integrated_init_energy_prob (1D_ALLOC_real -
      .def_property_readonly(
          "integrated_init_energy_prob",
          &PhotonElementProxy::integrated_init_energy_prob)

      .def(
          "__repr__",
          [](const PhotonElementProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PhotonElementProxy& self) {
            return PhotonElementProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const PhotonElementProxy& self, py::dict& memo) {
            return PhotonElementProxy(self);
          })

      ;

  // 1D PhotonElementProxy arrays are not used in structs/routines
  // 2D PhotonElementProxy arrays are not used in structs/routines
  // 3D PhotonElementProxy arrays are not used in structs/routines
}

// =============================================================================
// photon_material_struct
void init_photon_material_struct(
    py::module& m,
    py::class_<PhotonMaterialProxy>& cls) {
  cls.def(py::init<>())
      // PhotonMaterialProxy.f0_m1 (0D_NOT_complex - For multilayer_mirror only.
      .def_property(
          "f0_m1", &PhotonMaterialProxy::f0_m1, &PhotonMaterialProxy::set_f0_m1)
      // PhotonMaterialProxy.f0_m2 (0D_NOT_complex - For multilayer_mirror only.
      .def_property(
          "f0_m2", &PhotonMaterialProxy::f0_m2, &PhotonMaterialProxy::set_f0_m2)
      // PhotonMaterialProxy.f_0 (0D_NOT_complex -
      .def_property(
          "f_0", &PhotonMaterialProxy::f_0, &PhotonMaterialProxy::set_f_0)
      // PhotonMaterialProxy.f_h (0D_NOT_complex - Structure factor for H direction.
      .def_property(
          "f_h", &PhotonMaterialProxy::f_h, &PhotonMaterialProxy::set_f_h)
      // PhotonMaterialProxy.f_hbar (0D_NOT_complex - Structure factor for -H direction.
      .def_property(
          "f_hbar",
          &PhotonMaterialProxy::f_hbar,
          &PhotonMaterialProxy::set_f_hbar)
      // PhotonMaterialProxy.f_hkl (0D_NOT_complex - = sqrt(f_h * f_hbar)
      .def_property(
          "f_hkl", &PhotonMaterialProxy::f_hkl, &PhotonMaterialProxy::set_f_hkl)
      // PhotonMaterialProxy.h_norm (1D_NOT_real - Normalized H vector for crystals.
      .def_property_readonly("h_norm", &PhotonMaterialProxy::h_norm)
      // PhotonMaterialProxy.l_ref (1D_NOT_real - Crystal reference orbit displacement vector in element coords.
      .def_property_readonly("l_ref", &PhotonMaterialProxy::l_ref)

      .def(
          "__repr__",
          [](const PhotonMaterialProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PhotonMaterialProxy& self) {
            return PhotonMaterialProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const PhotonMaterialProxy& self, py::dict& memo) {
            return PhotonMaterialProxy(self);
          })

      ;

  // 1D PhotonMaterialProxy arrays are not used in structs/routines
  // 2D PhotonMaterialProxy arrays are not used in structs/routines
  // 3D PhotonMaterialProxy arrays are not used in structs/routines
}

// =============================================================================
// photon_reflect_surface_struct
void init_photon_reflect_surface_struct(
    py::module& m,
    py::class_<PhotonReflectSurfaceProxy>& cls) {
  cls.def(py::init<>())
      // PhotonReflectSurfaceProxy.name (0D_NOT_character -
      .def_property(
          "name",
          &PhotonReflectSurfaceProxy::name,
          &PhotonReflectSurfaceProxy::set_name)
      // PhotonReflectSurfaceProxy.description (0D_NOT_character - Descriptive name
      .def_property(
          "description",
          &PhotonReflectSurfaceProxy::description,
          &PhotonReflectSurfaceProxy::set_description)
      // PhotonReflectSurfaceProxy.reflectivity_file (0D_NOT_character -
      .def_property(
          "reflectivity_file",
          &PhotonReflectSurfaceProxy::reflectivity_file,
          &PhotonReflectSurfaceProxy::set_reflectivity_file)
      // PhotonReflectSurfaceProxy.table (1D_ALLOC_type -
      .def_property_readonly("table", &PhotonReflectSurfaceProxy::table)
      // PhotonReflectSurfaceProxy.surface_roughness_rms (0D_NOT_real - sigma in Dugan's notation
      .def_property(
          "surface_roughness_rms",
          &PhotonReflectSurfaceProxy::surface_roughness_rms,
          &PhotonReflectSurfaceProxy::set_surface_roughness_rms)
      // PhotonReflectSurfaceProxy.roughness_correlation_len (0D_NOT_real - T in Dugan's notation
      .def_property(
          "roughness_correlation_len",
          &PhotonReflectSurfaceProxy::roughness_correlation_len,
          &PhotonReflectSurfaceProxy::set_roughness_correlation_len)
      // PhotonReflectSurfaceProxy.ix_surface (0D_NOT_integer -
      .def_property(
          "ix_surface",
          &PhotonReflectSurfaceProxy::ix_surface,
          &PhotonReflectSurfaceProxy::set_ix_surface)

      .def(
          "__repr__",
          [](const PhotonReflectSurfaceProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PhotonReflectSurfaceProxy& self) {
            return PhotonReflectSurfaceProxy(
                self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const PhotonReflectSurfaceProxy& self, py::dict& memo) {
            return PhotonReflectSurfaceProxy(self);
          })

      ;

  // 1D PhotonReflectSurfaceProxy arrays are not used in structs/routines
  // 2D PhotonReflectSurfaceProxy arrays are not used in structs/routines
  // 3D PhotonReflectSurfaceProxy arrays are not used in structs/routines
}

// =============================================================================
// photon_reflect_table_struct
void init_photon_reflect_table_struct(
    py::module& m,
    py::class_<PhotonReflectTableProxy>& cls) {
  cls.def(py::init<>())
      // PhotonReflectTableProxy.angle (1D_ALLOC_real - Vector of angle values for %p_reflect
      .def_property_readonly("angle", &PhotonReflectTableProxy::angle)
      // PhotonReflectTableProxy.energy (1D_ALLOC_real - Vector of energy values for %p_reflect
      .def_property_readonly("energy", &PhotonReflectTableProxy::energy)
      // PhotonReflectTableProxy.int1 (1D_ALLOC_type -
      .def_property_readonly("int1", &PhotonReflectTableProxy::int1)
      // PhotonReflectTableProxy.p_reflect (2D_ALLOC_real - (angle, ev) probability. Log used for smooth surface reflection
      .def_property_readonly("p_reflect", &PhotonReflectTableProxy::p_reflect)
      // PhotonReflectTableProxy.max_energy (0D_NOT_real - maximum energy for this table
      .def_property(
          "max_energy",
          &PhotonReflectTableProxy::max_energy,
          &PhotonReflectTableProxy::set_max_energy)
      // PhotonReflectTableProxy.p_reflect_scratch (1D_ALLOC_real - Scratch space
      .def_property_readonly(
          "p_reflect_scratch", &PhotonReflectTableProxy::p_reflect_scratch)
      // PhotonReflectTableProxy.bragg_angle (1D_ALLOC_real - Bragg angle at energy values.
      .def_property_readonly(
          "bragg_angle", &PhotonReflectTableProxy::bragg_angle)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return PhotonReflectTableProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const PhotonReflectTableProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PhotonReflectTableProxy& self) {
            return PhotonReflectTableProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const PhotonReflectTableProxy& self, py::dict& memo) {
            return PhotonReflectTableProxy(self);
          })

      ;

  bind_FTypeArrayND<PhotonReflectTableProxyArray1D>(
      m, "PhotonReflectTableStructArray1D");
  bind_FTypeAlloc1D<PhotonReflectTableProxyAlloc1D>(
      m, "PhotonReflectTableStructAlloc1D");
  // 2D PhotonReflectTableProxy arrays are not used in structs/routines
  // 3D PhotonReflectTableProxy arrays are not used in structs/routines
}

// =============================================================================
// photon_target_struct
void init_photon_target_struct(
    py::module& m,
    py::class_<PhotonTargetProxy>& cls) {
  cls.def(py::init<>())
      // PhotonTargetProxy.type (0D_NOT_integer - or rectangular$
      .def_property(
          "type", &PhotonTargetProxy::type, &PhotonTargetProxy::set_type)
      // PhotonTargetProxy.n_corner (0D_NOT_integer -
      .def_property(
          "n_corner",
          &PhotonTargetProxy::n_corner,
          &PhotonTargetProxy::set_n_corner)
      // PhotonTargetProxy.ele_loc (0D_NOT_type -
      .def_property(
          "ele_loc",
          &PhotonTargetProxy::ele_loc,
          &PhotonTargetProxy::set_ele_loc)
      // PhotonTargetProxy.corner (1D_NOT_type -
      .def_property_readonly("corner", &PhotonTargetProxy::corner)
      // PhotonTargetProxy.center (0D_NOT_type -
      .def_property(
          "center", &PhotonTargetProxy::center, &PhotonTargetProxy::set_center)

      .def(
          "__repr__",
          [](const PhotonTargetProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PhotonTargetProxy& self) {
            return PhotonTargetProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const PhotonTargetProxy& self, py::dict& memo) {
            return PhotonTargetProxy(self);
          })

      ;

  // 1D PhotonTargetProxy arrays are not used in structs/routines
  // 2D PhotonTargetProxy arrays are not used in structs/routines
  // 3D PhotonTargetProxy arrays are not used in structs/routines
}

// =============================================================================
// pixel_detec_struct
void init_pixel_detec_struct(py::module& m, py::class_<PixelDetecProxy>& cls) {
  cls.def(py::init<>())
      // PixelDetecProxy.dr (1D_NOT_real -
      .def_property_readonly("dr", &PixelDetecProxy::dr)
      // PixelDetecProxy.r0 (1D_NOT_real -
      .def_property_readonly("r0", &PixelDetecProxy::r0)
      // PixelDetecProxy.n_track_tot (0D_NOT_integer8 - How many photons were launched from source element.
      .def_property(
          "n_track_tot",
          &PixelDetecProxy::n_track_tot,
          &PixelDetecProxy::set_n_track_tot)
      // PixelDetecProxy.n_hit_detec (0D_NOT_integer8 - How many photons hit the detector.
      .def_property(
          "n_hit_detec",
          &PixelDetecProxy::n_hit_detec,
          &PixelDetecProxy::set_n_hit_detec)
      // PixelDetecProxy.n_hit_pixel (0D_NOT_integer8 - How many photons hit the pixel grid of the detector.
      .def_property(
          "n_hit_pixel",
          &PixelDetecProxy::n_hit_pixel,
          &PixelDetecProxy::set_n_hit_pixel)
      // PixelDetecProxy.pt (2D_ALLOC_type - Grid of pixels
      .def_property_readonly("pt", &PixelDetecProxy::pt)

      .def(
          "__repr__",
          [](const PixelDetecProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PixelDetecProxy& self) {
            return PixelDetecProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const PixelDetecProxy& self, py::dict& memo) {
            return PixelDetecProxy(self);
          })

      ;

  // 1D PixelDetecProxy arrays are not used in structs/routines
  // 2D PixelDetecProxy arrays are not used in structs/routines
  // 3D PixelDetecProxy arrays are not used in structs/routines
}

// =============================================================================
// pixel_pt_struct
void init_pixel_pt_struct(py::module& m, py::class_<PixelPtProxy>& cls) {
  cls.def(py::init<>())
      // PixelPtProxy.n_photon (0D_NOT_integer8 -
      .def_property(
          "n_photon", &PixelPtProxy::n_photon, &PixelPtProxy::set_n_photon)
      // PixelPtProxy.E_x (0D_NOT_complex -
      .def_property("E_x", &PixelPtProxy::E_x, &PixelPtProxy::set_E_x)
      // PixelPtProxy.E_y (0D_NOT_complex -
      .def_property("E_y", &PixelPtProxy::E_y, &PixelPtProxy::set_E_y)
      // PixelPtProxy.intensity_x (0D_NOT_real -
      .def_property(
          "intensity_x",
          &PixelPtProxy::intensity_x,
          &PixelPtProxy::set_intensity_x)
      // PixelPtProxy.intensity_y (0D_NOT_real -
      .def_property(
          "intensity_y",
          &PixelPtProxy::intensity_y,
          &PixelPtProxy::set_intensity_y)
      // PixelPtProxy.intensity (0D_NOT_real -
      .def_property(
          "intensity", &PixelPtProxy::intensity, &PixelPtProxy::set_intensity)
      // PixelPtProxy.orbit (1D_NOT_real - x, Vx/c, y, Vy/c, dummy, E - E_ref.
      .def_property_readonly("orbit", &PixelPtProxy::orbit)
      // PixelPtProxy.orbit_rms (1D_NOT_real - RMS statistics.
      .def_property_readonly("orbit_rms", &PixelPtProxy::orbit_rms)
      // PixelPtProxy.init_orbit (1D_NOT_real - Initial orbit at start of lattice statistics.
      .def_property_readonly("init_orbit", &PixelPtProxy::init_orbit)
      // PixelPtProxy.init_orbit_rms (1D_NOT_real - Initial orbit at start of lattice RMS statistics.
      .def_property_readonly("init_orbit_rms", &PixelPtProxy::init_orbit_rms)

      .def("__repr__", [](const PixelPtProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PixelPtProxy& self) {
            return PixelPtProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const PixelPtProxy& self, py::dict& memo) {
            return PixelPtProxy(self);
          })

      ;

  // 1D PixelPtProxy arrays are not used in structs/routines
  bind_FTypeArrayND<PixelPtProxyArray2D>(m, "PixelPtStructArray2D");
  // 3D PixelPtProxy arrays are not used in structs/routines
}

// =============================================================================
// pre_tracker_struct
void init_pre_tracker_struct(py::module& m, py::class_<PreTrackerProxy>& cls) {
  cls.def(py::init<>())
      // PreTrackerProxy.who (0D_NOT_integer - Can be opal$, or impactt$
      .def_property("who", &PreTrackerProxy::who, &PreTrackerProxy::set_who)
      // PreTrackerProxy.ix_ele_start (0D_NOT_integer -
      .def_property(
          "ix_ele_start",
          &PreTrackerProxy::ix_ele_start,
          &PreTrackerProxy::set_ix_ele_start)
      // PreTrackerProxy.ix_ele_end (0D_NOT_integer -
      .def_property(
          "ix_ele_end",
          &PreTrackerProxy::ix_ele_end,
          &PreTrackerProxy::set_ix_ele_end)
      // PreTrackerProxy.input_file (0D_NOT_character -
      .def_property(
          "input_file",
          &PreTrackerProxy::input_file,
          &PreTrackerProxy::set_input_file)

      .def(
          "__repr__",
          [](const PreTrackerProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PreTrackerProxy& self) {
            return PreTrackerProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const PreTrackerProxy& self, py::dict& memo) {
            return PreTrackerProxy(self);
          })

      ;

  // 1D PreTrackerProxy arrays are not used in structs/routines
  // 2D PreTrackerProxy arrays are not used in structs/routines
  // 3D PreTrackerProxy arrays are not used in structs/routines
}

// =============================================================================
// ptc_normal_form_struct
void init_ptc_normal_form_struct(
    py::module& m,
    py::class_<PtcNormalFormProxy>& cls) {
  cls.def(py::init<>())
      // PtcNormalFormProxy.ele_origin (0D_PTR_type - Element at which the on-turn map was created.
      .def_property(
          "ele_origin",
          &PtcNormalFormProxy::ele_origin,
          &PtcNormalFormProxy::set_ele_origin)
      // PtcNormalFormProxy.orb0 (1D_NOT_real - Closed orbit at element.
      .def_property_readonly("orb0", &PtcNormalFormProxy::orb0)
      // PtcNormalFormProxy.valid_map (0D_NOT_logical -
      .def_property(
          "valid_map",
          &PtcNormalFormProxy::valid_map,
          &PtcNormalFormProxy::set_valid_map)

      .def(
          "__repr__",
          [](const PtcNormalFormProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PtcNormalFormProxy& self) {
            return PtcNormalFormProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const PtcNormalFormProxy& self, py::dict& memo) {
            return PtcNormalFormProxy(self);
          })

      ;

  // 1D PtcNormalFormProxy arrays are not used in structs/routines
  // 2D PtcNormalFormProxy arrays are not used in structs/routines
  // 3D PtcNormalFormProxy arrays are not used in structs/routines
}