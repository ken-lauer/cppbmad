#include "pybmad/generated/structs_p.hpp"

#include <cstdint>
#include <functional>

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

using namespace Pybmad;
namespace nb = nanobind;

// =============================================================================
// parser_controller_struct
void init_parser_controller_struct(nb::module_ &m, nb::class_<ParserControllerStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::vector<double>>,
             std::optional<int>>(),
         nb::arg("name") = nb::none(),
         nb::arg("attrib_name") = nb::none(),
         nb::arg("y_knot") = nb::none(),
         nb::arg("n_stk") = nb::none()
  )
      .def_prop_rw("name", &ParserControllerStruct::name, &ParserControllerStruct::set_name)
      .def_prop_rw(
          "attrib_name",
          &ParserControllerStruct::attrib_name,
          &ParserControllerStruct::set_attrib_name
      )
      .def_prop_ro(
          "stack",
          &ParserControllerStruct::stack,
          nb::keep_alive<0, 1>(),
          "Arithmetic expression stack"
      )
      .def_prop_rw(
          "y_knot",
          &ParserControllerStruct::y_knot,
          &ParserControllerStruct::set_y_knot,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw("n_stk", &ParserControllerStruct::n_stk, &ParserControllerStruct::set_n_stk)
      .def_static(
          "new_array1d",
          [](int sz) { return ParserControllerStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = ParserControllerStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const ParserControllerStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ParserControllerStruct &self) {
            return ParserControllerStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ParserControllerStruct &self, nb::dict &memo) {
            return ParserControllerStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const ParserControllerStruct &self, const ParserControllerStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ParserControllerStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  bind_1d_type_array_pair<ParserControllerStructArray1D, ParserControllerStructAlloc1D>(
      m,
      "ParserControllerStructArray1D",
      "ParserControllerStructAlloc1D"
  );
  // 2D ParserControllerStruct arrays are not used in structs/routines
  // 3D ParserControllerStruct arrays are not used in structs/routines
}

// =============================================================================
// parser_ele_struct
void init_parser_ele_struct(nb::module_ &m, nb::class_<ParserEleStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::string>,
             std::optional<int>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<double>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<std::string>>(),
         nb::arg("ref_name") = nb::none(),
         nb::arg("ix_super_ref_multipass") = nb::none(),
         nb::arg("ele_name") = nb::none(),
         nb::arg("lat_file") = nb::none(),
         nb::arg("offset") = nb::none(),
         nb::arg("ix_line_in_file") = nb::none(),
         nb::arg("ix_count") = nb::none(),
         nb::arg("ele_pt") = nb::none(),
         nb::arg("ref_pt") = nb::none(),
         nb::arg("index") = nb::none(),
         nb::arg("superposition_command_here") = nb::none(),
         nb::arg("superposition_has_been_set") = nb::none(),
         nb::arg("wrap_superimpose") = nb::none(),
         nb::arg("create_jumbo_slave") = nb::none(),
         nb::arg("is_range") = nb::none(),
         nb::arg("default_attrib") = nb::none()
  )
      .def_prop_ro("control", &ParserEleStruct::control, nb::keep_alive<0, 1>())
      .def_prop_ro("field_overlaps", &ParserEleStruct::field_overlaps, nb::keep_alive<0, 1>())
      .def_prop_rw("ref_name", &ParserEleStruct::ref_name, &ParserEleStruct::set_ref_name)
      .def_prop_rw(
          "ix_super_ref_multipass",
          &ParserEleStruct::ix_super_ref_multipass,
          &ParserEleStruct::set_ix_super_ref_multipass,
          "Multipass index for superimpose reference element."
      )
      .def_prop_rw(
          "ele_name",
          &ParserEleStruct::ele_name,
          &ParserEleStruct::set_ele_name,
          "For fork element or superimpose statement."
      )
      .def_prop_ro(
          "names1",
          &ParserEleStruct::names1,
          nb::keep_alive<0, 1>(),
          "Currently just used by feedback element."
      )
      .def_prop_ro(
          "names2",
          &ParserEleStruct::names2,
          nb::keep_alive<0, 1>(),
          "Currently just used by feedback element."
      )
      .def_prop_rw(
          "lat_file",
          &ParserEleStruct::lat_file,
          &ParserEleStruct::set_lat_file,
          "File where element was defined."
      )
      .def_prop_rw("offset", &ParserEleStruct::offset, &ParserEleStruct::set_offset)
      .def_prop_rw(
          "ix_line_in_file",
          &ParserEleStruct::ix_line_in_file,
          &ParserEleStruct::set_ix_line_in_file,
          "Line in file where element was defined."
      )
      .def_prop_rw("ix_count", &ParserEleStruct::ix_count, &ParserEleStruct::set_ix_count)
      .def_prop_rw("ele_pt", &ParserEleStruct::ele_pt, &ParserEleStruct::set_ele_pt)
      .def_prop_rw("ref_pt", &ParserEleStruct::ref_pt, &ParserEleStruct::set_ref_pt)
      .def_prop_rw("index", &ParserEleStruct::index, &ParserEleStruct::set_index)
      .def_prop_rw(
          "superposition_command_here",
          &ParserEleStruct::superposition_command_here,
          &ParserEleStruct::set_superposition_command_here
      )
      .def_prop_rw(
          "superposition_has_been_set",
          &ParserEleStruct::superposition_has_been_set,
          &ParserEleStruct::set_superposition_has_been_set
      )
      .def_prop_rw(
          "wrap_superimpose",
          &ParserEleStruct::wrap_superimpose,
          &ParserEleStruct::set_wrap_superimpose
      )
      .def_prop_rw(
          "create_jumbo_slave",
          &ParserEleStruct::create_jumbo_slave,
          &ParserEleStruct::set_create_jumbo_slave
      )
      .def_prop_rw(
          "is_range",
          &ParserEleStruct::is_range,
          &ParserEleStruct::set_is_range,
          "For girders"
      )
      .def_prop_rw(
          "default_attrib",
          &ParserEleStruct::default_attrib,
          &ParserEleStruct::set_default_attrib,
          "For group/overlay elements: slave attribute"
      )

      .def("__repr__", [](const ParserEleStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ParserEleStruct &self) {
            return ParserEleStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const ParserEleStruct &self, nb::dict &memo) { return ParserEleStruct(self); }
      )
      .def(
          "__eq__",
          [](const ParserEleStruct &self, const ParserEleStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const ParserEleStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D ParserEleStruct arrays are not used in structs/routines
  // 2D ParserEleStruct arrays are not used in structs/routines
  // 3D ParserEleStruct arrays are not used in structs/routines
}

// =============================================================================
// photon_element_struct
void init_photon_element_struct(nb::module_ &m, nb::class_<PhotonElementStruct> &cls) {
  cls.def(
         "__init__",
         [](PhotonElementStruct *self,
            const SurfaceCurvatureStruct *curvature,
            const PhotonTargetStruct *target,
            const PhotonMaterialStruct *material,
            const SurfaceSegmentedStruct *segmented,
            const SurfaceHMisalignStruct *h_misalign,
            const SurfaceDisplacementStruct *displacement,
            const PixelDetecStruct *pixel,
            std::optional<int> reflectivity_table_type,
            const PhotonReflectTableStruct *reflectivity_table_sigma,
            const PhotonReflectTableStruct *reflectivity_table_pi,
            std::optional<std::vector<double>> integrated_init_energy_prob) {
           new (self) PhotonElementStruct(
               ptr_to_opt_ref(curvature),
               ptr_to_opt_ref(target),
               ptr_to_opt_ref(material),
               ptr_to_opt_ref(segmented),
               ptr_to_opt_ref(h_misalign),
               ptr_to_opt_ref(displacement),
               ptr_to_opt_ref(pixel),
               reflectivity_table_type,
               ptr_to_opt_ref(reflectivity_table_sigma),
               ptr_to_opt_ref(reflectivity_table_pi),
               integrated_init_energy_prob
           );
         },
         nb::arg("curvature") = nb::none(),
         nb::arg("target") = nb::none(),
         nb::arg("material") = nb::none(),
         nb::arg("segmented") = nb::none(),
         nb::arg("h_misalign") = nb::none(),
         nb::arg("displacement") = nb::none(),
         nb::arg("pixel") = nb::none(),
         nb::arg("reflectivity_table_type") = nb::none(),
         nb::arg("reflectivity_table_sigma") = nb::none(),
         nb::arg("reflectivity_table_pi") = nb::none(),
         nb::arg("integrated_init_energy_prob") = nb::none()
  )
      .def_prop_rw(
          "curvature",
          &PhotonElementStruct::curvature,
          &PhotonElementStruct::set_curvature,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "target",
          &PhotonElementStruct::target,
          &PhotonElementStruct::set_target,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "material",
          &PhotonElementStruct::material,
          &PhotonElementStruct::set_material,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "segmented",
          &PhotonElementStruct::segmented,
          &PhotonElementStruct::set_segmented,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "h_misalign",
          &PhotonElementStruct::h_misalign,
          &PhotonElementStruct::set_h_misalign,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "displacement",
          &PhotonElementStruct::displacement,
          &PhotonElementStruct::set_displacement,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "pixel",
          &PhotonElementStruct::pixel,
          &PhotonElementStruct::set_pixel,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "reflectivity_table_type",
          &PhotonElementStruct::reflectivity_table_type,
          &PhotonElementStruct::set_reflectivity_table_type
      )
      .def_prop_rw(
          "reflectivity_table_sigma",
          &PhotonElementStruct::reflectivity_table_sigma,
          &PhotonElementStruct::set_reflectivity_table_sigma,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "If polarization is ignored use sigma table."
      )
      .def_prop_rw(
          "reflectivity_table_pi",
          &PhotonElementStruct::reflectivity_table_pi,
          &PhotonElementStruct::set_reflectivity_table_pi,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro(
          "init_energy_prob",
          &PhotonElementStruct::init_energy_prob,
          nb::keep_alive<0, 1>(),
          "Initial energy probability density"
      )
      .def_prop_rw(
          "integrated_init_energy_prob",
          &PhotonElementStruct::integrated_init_energy_prob,
          &PhotonElementStruct::set_integrated_init_energy_prob,
          nb::for_getter(nb::keep_alive<0, 1>())
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
          [](const PhotonElementStruct &self, nb::dict &memo) { return PhotonElementStruct(self); }
      )
      .def(
          "__eq__",
          [](const PhotonElementStruct &self, const PhotonElementStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const PhotonElementStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D PhotonElementStruct arrays are not used in structs/routines
  // 2D PhotonElementStruct arrays are not used in structs/routines
  // 3D PhotonElementStruct arrays are not used in structs/routines
}

// =============================================================================
// photon_material_struct
void init_photon_material_struct(nb::module_ &m, nb::class_<PhotonMaterialStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::complex<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>>(),
         nb::arg("f0_m1") = nb::none(),
         nb::arg("f0_m2") = nb::none(),
         nb::arg("f_0") = nb::none(),
         nb::arg("f_h") = nb::none(),
         nb::arg("f_hbar") = nb::none(),
         nb::arg("f_hkl") = nb::none(),
         nb::arg("h_norm") = nb::none(),
         nb::arg("l_ref") = nb::none()
  )
      .def_prop_rw(
          "f0_m1",
          &PhotonMaterialStruct::f0_m1,
          &PhotonMaterialStruct::set_f0_m1,
          "For multilayer_mirror only."
      )
      .def_prop_rw(
          "f0_m2",
          &PhotonMaterialStruct::f0_m2,
          &PhotonMaterialStruct::set_f0_m2,
          "For multilayer_mirror only."
      )
      .def_prop_rw("f_0", &PhotonMaterialStruct::f_0, &PhotonMaterialStruct::set_f_0)
      .def_prop_rw(
          "f_h",
          &PhotonMaterialStruct::f_h,
          &PhotonMaterialStruct::set_f_h,
          "Structure factor for H direction."
      )
      .def_prop_rw(
          "f_hbar",
          &PhotonMaterialStruct::f_hbar,
          &PhotonMaterialStruct::set_f_hbar,
          "Structure factor for -H direction."
      )
      .def_prop_rw(
          "f_hkl",
          &PhotonMaterialStruct::f_hkl,
          &PhotonMaterialStruct::set_f_hkl,
          "= sqrt(f_h * f_hbar)"
      )
      .def_prop_rw(
          "h_norm",
          &PhotonMaterialStruct::h_norm,
          &PhotonMaterialStruct::set_h_norm,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Normalized H vector for crystals."
      )
      .def_prop_rw(
          "l_ref",
          &PhotonMaterialStruct::l_ref,
          &PhotonMaterialStruct::set_l_ref,
          nb::for_getter(nb::keep_alive<0, 1>()),
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
          [](const PhotonMaterialStruct &self, nb::dict &memo) {
            return PhotonMaterialStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const PhotonMaterialStruct &self, const PhotonMaterialStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const PhotonMaterialStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
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
    nb::module_ &m,
    nb::class_<PhotonReflectSurfaceStruct> &cls
) {
  cls.def(
         nb::init<
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<double>,
             std::optional<double>,
             std::optional<int>>(),
         nb::arg("name") = nb::none(),
         nb::arg("description") = nb::none(),
         nb::arg("reflectivity_file") = nb::none(),
         nb::arg("surface_roughness_rms") = nb::none(),
         nb::arg("roughness_correlation_len") = nb::none(),
         nb::arg("ix_surface") = nb::none()
  )
      .def_prop_rw("name", &PhotonReflectSurfaceStruct::name, &PhotonReflectSurfaceStruct::set_name)
      .def_prop_rw(
          "description",
          &PhotonReflectSurfaceStruct::description,
          &PhotonReflectSurfaceStruct::set_description,
          "Descriptive name"
      )
      .def_prop_rw(
          "reflectivity_file",
          &PhotonReflectSurfaceStruct::reflectivity_file,
          &PhotonReflectSurfaceStruct::set_reflectivity_file
      )
      .def_prop_ro("table", &PhotonReflectSurfaceStruct::table, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "surface_roughness_rms",
          &PhotonReflectSurfaceStruct::surface_roughness_rms,
          &PhotonReflectSurfaceStruct::set_surface_roughness_rms,
          "sigma in Dugan's notation"
      )
      .def_prop_rw(
          "roughness_correlation_len",
          &PhotonReflectSurfaceStruct::roughness_correlation_len,
          &PhotonReflectSurfaceStruct::set_roughness_correlation_len,
          "T in Dugan's notation"
      )
      .def_prop_rw(
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
          [](const PhotonReflectSurfaceStruct &self, nb::dict &memo) {
            return PhotonReflectSurfaceStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const PhotonReflectSurfaceStruct &self, const PhotonReflectSurfaceStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const PhotonReflectSurfaceStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D PhotonReflectSurfaceStruct arrays are not used in structs/routines
  // 2D PhotonReflectSurfaceStruct arrays are not used in structs/routines
  // 3D PhotonReflectSurfaceStruct arrays are not used in structs/routines
}

// =============================================================================
// photon_reflect_table_struct
void init_photon_reflect_table_struct(nb::module_ &m, nb::class_<PhotonReflectTableStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>>(),
         nb::arg("angle") = nb::none(),
         nb::arg("energy") = nb::none(),
         nb::arg("p_reflect") = nb::none(),
         nb::arg("max_energy") = nb::none(),
         nb::arg("p_reflect_scratch") = nb::none(),
         nb::arg("bragg_angle") = nb::none()
  )
      .def_prop_rw(
          "angle",
          &PhotonReflectTableStruct::angle,
          &PhotonReflectTableStruct::set_angle,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Vector of angle values for %p_reflect"
      )
      .def_prop_rw(
          "energy",
          &PhotonReflectTableStruct::energy,
          &PhotonReflectTableStruct::set_energy,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Vector of energy values for %p_reflect"
      )
      .def_prop_ro("int1", &PhotonReflectTableStruct::int1, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "p_reflect",
          &PhotonReflectTableStruct::p_reflect,
          &PhotonReflectTableStruct::set_p_reflect,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "(angle, ev) probability. Log used for smooth surface reflection"
      )
      .def_prop_rw(
          "max_energy",
          &PhotonReflectTableStruct::max_energy,
          &PhotonReflectTableStruct::set_max_energy,
          "maximum energy for this table"
      )
      .def_prop_rw(
          "p_reflect_scratch",
          &PhotonReflectTableStruct::p_reflect_scratch,
          &PhotonReflectTableStruct::set_p_reflect_scratch,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Scratch space"
      )
      .def_prop_rw(
          "bragg_angle",
          &PhotonReflectTableStruct::bragg_angle,
          &PhotonReflectTableStruct::set_bragg_angle,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Bragg angle at energy values."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return PhotonReflectTableStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = PhotonReflectTableStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const PhotonReflectTableStruct &self, nb::dict &memo) {
            return PhotonReflectTableStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const PhotonReflectTableStruct &self, const PhotonReflectTableStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const PhotonReflectTableStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
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
void init_photon_target_struct(nb::module_ &m, nb::class_<PhotonTargetStruct> &cls) {
  cls.def(
         "__init__",
         [](PhotonTargetStruct *self,
            std::optional<int> type,
            std::optional<int> n_corner,
            const LatEleLocStruct *ele_loc,
            const TargetPointStruct *center) {
           new (self)
               PhotonTargetStruct(type, n_corner, ptr_to_opt_ref(ele_loc), ptr_to_opt_ref(center));
         },
         nb::arg("type") = nb::none(),
         nb::arg("n_corner") = nb::none(),
         nb::arg("ele_loc") = nb::none(),
         nb::arg("center") = nb::none()
  )
      .def_prop_rw(
          "type",
          &PhotonTargetStruct::type,
          &PhotonTargetStruct::set_type,
          "or rectangular$"
      )
      .def_prop_rw("n_corner", &PhotonTargetStruct::n_corner, &PhotonTargetStruct::set_n_corner)
      .def_prop_rw(
          "ele_loc",
          &PhotonTargetStruct::ele_loc,
          &PhotonTargetStruct::set_ele_loc,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro("corner", &PhotonTargetStruct::corner, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "center",
          &PhotonTargetStruct::center,
          &PhotonTargetStruct::set_center,
          nb::for_getter(nb::keep_alive<0, 1>())
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
          [](const PhotonTargetStruct &self, nb::dict &memo) { return PhotonTargetStruct(self); }
      )
      .def(
          "__eq__",
          [](const PhotonTargetStruct &self, const PhotonTargetStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const PhotonTargetStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D PhotonTargetStruct arrays are not used in structs/routines
  // 2D PhotonTargetStruct arrays are not used in structs/routines
  // 3D PhotonTargetStruct arrays are not used in structs/routines
}

// =============================================================================
// pixel_detec_struct
void init_pixel_detec_struct(nb::module_ &m, nb::class_<PixelDetecStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<int64_t>,
             std::optional<int64_t>,
             std::optional<int64_t>>(),
         nb::arg("dr") = nb::none(),
         nb::arg("r0") = nb::none(),
         nb::arg("n_track_tot") = nb::none(),
         nb::arg("n_hit_detec") = nb::none(),
         nb::arg("n_hit_pixel") = nb::none()
  )
      .def_prop_rw(
          "dr",
          &PixelDetecStruct::dr,
          &PixelDetecStruct::set_dr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "r0",
          &PixelDetecStruct::r0,
          &PixelDetecStruct::set_r0,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "n_track_tot",
          &PixelDetecStruct::n_track_tot,
          &PixelDetecStruct::set_n_track_tot,
          "How many photons were launched from source element."
      )
      .def_prop_rw(
          "n_hit_detec",
          &PixelDetecStruct::n_hit_detec,
          &PixelDetecStruct::set_n_hit_detec,
          "How many photons hit the detector."
      )
      .def_prop_rw(
          "n_hit_pixel",
          &PixelDetecStruct::n_hit_pixel,
          &PixelDetecStruct::set_n_hit_pixel,
          "How many photons hit the pixel grid of the detector."
      )
      .def_prop_ro("pt", &PixelDetecStruct::pt, nb::keep_alive<0, 1>(), "Grid of pixels")

      .def("__repr__", [](const PixelDetecStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PixelDetecStruct &self) {
            return PixelDetecStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const PixelDetecStruct &self, nb::dict &memo) { return PixelDetecStruct(self); }
      )
      .def(
          "__eq__",
          [](const PixelDetecStruct &self, const PixelDetecStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const PixelDetecStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D PixelDetecStruct arrays are not used in structs/routines
  // 2D PixelDetecStruct arrays are not used in structs/routines
  // 3D PixelDetecStruct arrays are not used in structs/routines
}

// =============================================================================
// pixel_pt_struct
void init_pixel_pt_struct(nb::module_ &m, nb::class_<PixelPtStruct> &cls) {
  cls.def(
         nb::init<
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
         nb::arg("n_photon") = nb::none(),
         nb::arg("E_x") = nb::none(),
         nb::arg("E_y") = nb::none(),
         nb::arg("intensity_x") = nb::none(),
         nb::arg("intensity_y") = nb::none(),
         nb::arg("intensity") = nb::none(),
         nb::arg("orbit") = nb::none(),
         nb::arg("orbit_rms") = nb::none(),
         nb::arg("init_orbit") = nb::none(),
         nb::arg("init_orbit_rms") = nb::none()
  )
      .def_prop_rw("n_photon", &PixelPtStruct::n_photon, &PixelPtStruct::set_n_photon)
      .def_prop_rw("E_x", &PixelPtStruct::E_x, &PixelPtStruct::set_E_x)
      .def_prop_rw("E_y", &PixelPtStruct::E_y, &PixelPtStruct::set_E_y)
      .def_prop_rw("intensity_x", &PixelPtStruct::intensity_x, &PixelPtStruct::set_intensity_x)
      .def_prop_rw("intensity_y", &PixelPtStruct::intensity_y, &PixelPtStruct::set_intensity_y)
      .def_prop_rw("intensity", &PixelPtStruct::intensity, &PixelPtStruct::set_intensity)
      .def_prop_rw(
          "orbit",
          &PixelPtStruct::orbit,
          &PixelPtStruct::set_orbit,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "x, Vx/c, y, Vy/c, dummy, E - E_ref."
      )
      .def_prop_rw(
          "orbit_rms",
          &PixelPtStruct::orbit_rms,
          &PixelPtStruct::set_orbit_rms,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "RMS statistics."
      )
      .def_prop_rw(
          "init_orbit",
          &PixelPtStruct::init_orbit,
          &PixelPtStruct::set_init_orbit,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Initial orbit at start of lattice statistics."
      )
      .def_prop_rw(
          "init_orbit_rms",
          &PixelPtStruct::init_orbit_rms,
          &PixelPtStruct::set_init_orbit_rms,
          nb::for_getter(nb::keep_alive<0, 1>()),
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
          [](const PixelPtStruct &self, nb::dict &memo) { return PixelPtStruct(self); }
      )
      .def(
          "__eq__",
          [](const PixelPtStruct &self, const PixelPtStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const PixelPtStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D PixelPtStruct arrays are not used in structs/routines
  bind_FTypeArrayND<PixelPtStructArray2D>(m, "PixelPtStructArray2D");
  // 3D PixelPtStruct arrays are not used in structs/routines
}

// =============================================================================
// pre_tracker_struct
void init_pre_tracker_struct(nb::module_ &m, nb::class_<PreTrackerStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<std::string>>(),
         nb::arg("who") = nb::none(),
         nb::arg("ix_ele_start") = nb::none(),
         nb::arg("ix_ele_end") = nb::none(),
         nb::arg("input_file") = nb::none()
  )
      .def_prop_rw(
          "who",
          &PreTrackerStruct::who,
          &PreTrackerStruct::set_who,
          "Can be opal$, or impactt$"
      )
      .def_prop_rw(
          "ix_ele_start",
          &PreTrackerStruct::ix_ele_start,
          &PreTrackerStruct::set_ix_ele_start
      )
      .def_prop_rw("ix_ele_end", &PreTrackerStruct::ix_ele_end, &PreTrackerStruct::set_ix_ele_end)
      .def_prop_rw("input_file", &PreTrackerStruct::input_file, &PreTrackerStruct::set_input_file)

      .def("__repr__", [](const PreTrackerStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const PreTrackerStruct &self) {
            return PreTrackerStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const PreTrackerStruct &self, nb::dict &memo) { return PreTrackerStruct(self); }
      )
      .def(
          "__eq__",
          [](const PreTrackerStruct &self, const PreTrackerStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const PreTrackerStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D PreTrackerStruct arrays are not used in structs/routines
  // 2D PreTrackerStruct arrays are not used in structs/routines
  // 3D PreTrackerStruct arrays are not used in structs/routines
}

// =============================================================================
// ptc_normal_form_struct
void init_ptc_normal_form_struct(nb::module_ &m, nb::class_<PtcNormalFormStruct> &cls) {
  cls.def(
         "__init__",
         [](PtcNormalFormStruct *self,
            const EleStruct *ele_origin,
            std::optional<std::vector<double>> orb0,
            std::optional<bool> valid_map) {
           new (self) PtcNormalFormStruct(ptr_to_opt_ref(ele_origin), orb0, valid_map);
         },
         nb::arg("ele_origin") = nb::none(),
         nb::arg("orb0") = nb::none(),
         nb::arg("valid_map") = nb::none()
  )
      .def_prop_rw(
          "ele_origin",
          &PtcNormalFormStruct::ele_origin,
          &PtcNormalFormStruct::set_ele_origin,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Element at which the on-turn map was created."
      )
      .def_prop_rw(
          "orb0",
          &PtcNormalFormStruct::orb0,
          &PtcNormalFormStruct::set_orb0,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Closed orbit at element."
      )
      .def_prop_rw(
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
          [](const PtcNormalFormStruct &self, nb::dict &memo) { return PtcNormalFormStruct(self); }
      )
      .def(
          "__eq__",
          [](const PtcNormalFormStruct &self, const PtcNormalFormStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const PtcNormalFormStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D PtcNormalFormStruct arrays are not used in structs/routines
  // 2D PtcNormalFormStruct arrays are not used in structs/routines
  // 3D PtcNormalFormStruct arrays are not used in structs/routines
}