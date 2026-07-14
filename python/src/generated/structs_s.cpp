#include "pybmad/generated/structs_s.hpp"

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
// seq_ele_struct
void init_seq_ele_struct(nb::module_ &m, nb::class_<SeqEleStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<std::string>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<int>>(),
         nb::arg("name") = nb::none(),
         nb::arg("tag") = nb::none(),
         nb::arg("slice_start") = nb::none(),
         nb::arg("slice_end") = nb::none(),
         nb::arg("type") = nb::none(),
         nb::arg("ix_ele") = nb::none(),
         nb::arg("ix_arg") = nb::none(),
         nb::arg("rep_count") = nb::none(),
         nb::arg("ele_order_reflect") = nb::none(),
         nb::arg("ele_orientation") = nb::none()
  )
      .def_prop_rw(
          "name",
          &SeqEleStruct::name,
          &SeqEleStruct::set_name,
          "name of element, subline, or sublist"
      )
      .def_prop_ro("actual_arg", &SeqEleStruct::actual_arg, nb::keep_alive<0, 1>())
      .def_prop_rw("tag", &SeqEleStruct::tag, &SeqEleStruct::set_tag, "tag name.")
      .def_prop_rw(
          "slice_start",
          &SeqEleStruct::slice_start,
          &SeqEleStruct::set_slice_start,
          "For 'my_line[start:end]' slice constructs."
      )
      .def_prop_rw(
          "slice_end",
          &SeqEleStruct::slice_end,
          &SeqEleStruct::set_slice_end,
          "For 'my_line[start:end]' slice constructs."
      )
      .def_prop_rw(
          "type",
          &SeqEleStruct::type,
          &SeqEleStruct::set_type,
          "LINE$, REPLACEMENT_LINE$, LIST$, ELEMENT$"
      )
      .def_prop_rw(
          "ix_ele",
          &SeqEleStruct::ix_ele,
          &SeqEleStruct::set_ix_ele,
          "if an element: pointer to ELE array if a line or list: pointer to SEQ array"
      )
      .def_prop_rw(
          "ix_arg",
          &SeqEleStruct::ix_arg,
          &SeqEleStruct::set_ix_arg,
          "index in arg list (for replacement lines)"
      )
      .def_prop_rw(
          "rep_count",
          &SeqEleStruct::rep_count,
          &SeqEleStruct::set_rep_count,
          "how many copies of an element"
      )
      .def_prop_rw(
          "ele_order_reflect",
          &SeqEleStruct::ele_order_reflect,
          &SeqEleStruct::set_ele_order_reflect,
          "Travel through ele sequence in reverse order"
      )
      .def_prop_rw(
          "ele_orientation",
          &SeqEleStruct::ele_orientation,
          &SeqEleStruct::set_ele_orientation,
          "element has reverse orientation."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return SeqEleStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = SeqEleStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const SeqEleStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SeqEleStruct &self) {
            return SeqEleStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SeqEleStruct &self, nb::dict &memo) { return SeqEleStruct(self); }
      )
      .def(
          "__eq__",
          [](const SeqEleStruct &self, const SeqEleStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const SeqEleStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<SeqEleStructArray1D, SeqEleStructAlloc1D>(
      m,
      "SeqEleStructArray1D",
      "SeqEleStructAlloc1D"
  );
  // 2D SeqEleStruct arrays are not used in structs/routines
  // 3D SeqEleStruct arrays are not used in structs/routines
}

// =============================================================================
// seq_struct
void init_seq_struct(nb::module_ &m, nb::class_<SeqStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::string>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<int>,
             std::optional<std::string>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<bool>,
             std::optional<bool>>(),
         nb::arg("name") = nb::none(),
         nb::arg("type") = nb::none(),
         nb::arg("ix_list") = nb::none(),
         nb::arg("list_upcount") = nb::none(),
         nb::arg("index") = nb::none(),
         nb::arg("file_name") = nb::none(),
         nb::arg("ix_file_line") = nb::none(),
         nb::arg("multipass") = nb::none(),
         nb::arg("ptc_layout") = nb::none(),
         nb::arg("active") = nb::none()
  )
      .def_prop_rw("name", &SeqStruct::name, &SeqStruct::set_name, "name of sequence")
      .def_prop_ro("ele", &SeqStruct::ele, nb::keep_alive<0, 1>(), "Elements in the sequence")
      .def_prop_ro("dummy_arg", &SeqStruct::dummy_arg, nb::keep_alive<0, 1>())
      .def_prop_ro(
          "corresponding_actual_arg",
          &SeqStruct::corresponding_actual_arg,
          nb::keep_alive<0, 1>()
      )
      .def_prop_rw(
          "type",
          &SeqStruct::type,
          &SeqStruct::set_type,
          "LINE$, REPLACEMENT_LINE$ or LIST$"
      )
      .def_prop_rw(
          "ix_list",
          &SeqStruct::ix_list,
          &SeqStruct::set_ix_list,
          "Current index for lists"
      )
      .def_prop_rw("list_upcount", &SeqStruct::list_upcount, &SeqStruct::set_list_upcount)
      .def_prop_rw(
          "index",
          &SeqStruct::index,
          &SeqStruct::set_index,
          "Alphabetical order sorted index"
      )
      .def_prop_rw(
          "file_name",
          &SeqStruct::file_name,
          &SeqStruct::set_file_name,
          "File where sequence is defined"
      )
      .def_prop_rw(
          "ix_file_line",
          &SeqStruct::ix_file_line,
          &SeqStruct::set_ix_file_line,
          "Line number in file where sequence is defined"
      )
      .def_prop_rw("multipass", &SeqStruct::multipass, &SeqStruct::set_multipass)
      .def_prop_rw(
          "ptc_layout",
          &SeqStruct::ptc_layout,
          &SeqStruct::set_ptc_layout,
          "Put in separate PTC layout"
      )
      .def_prop_rw(
          "active",
          &SeqStruct::active,
          &SeqStruct::set_active,
          "Used to prevent infinite loops."
      )
      .def_static(
          "new_array1d",
          [](int sz) { return SeqStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = SeqStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const SeqStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SeqStruct &self) {
            return SeqStruct(self); // under-the-hood fortran copy
          }
      )
      .def("__deepcopy__", [](const SeqStruct &self, nb::dict &memo) { return SeqStruct(self); })
      .def(
          "__eq__",
          [](const SeqStruct &self, const SeqStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const SeqStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<SeqStructArray1D, SeqStructAlloc1D>(
      m,
      "SeqStructArray1D",
      "SeqStructAlloc1D"
  );
  // 2D SeqStruct arrays are not used in structs/routines
  // 3D SeqStruct arrays are not used in structs/routines
}

// =============================================================================
// space_charge_common_struct
void init_space_charge_common_struct(nb::module_ &m, nb::class_<SpaceChargeCommonStruct> &cls) {
  cls.def(
         nb::init<
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
         nb::arg("ds_track_step") = nb::none(),
         nb::arg("dt_track_step") = nb::none(),
         nb::arg("cathode_strength_cutoff") = nb::none(),
         nb::arg("rel_tol_tracking") = nb::none(),
         nb::arg("abs_tol_tracking") = nb::none(),
         nb::arg("beam_chamber_height") = nb::none(),
         nb::arg("lsc_sigma_cutoff") = nb::none(),
         nb::arg("particle_sigma_cutoff") = nb::none(),
         nb::arg("space_charge_mesh_size") = nb::none(),
         nb::arg("csr3d_mesh_size") = nb::none(),
         nb::arg("n_bin") = nb::none(),
         nb::arg("particle_bin_span") = nb::none(),
         nb::arg("n_shield_images") = nb::none(),
         nb::arg("sc_min_in_bin") = nb::none(),
         nb::arg("lsc_kick_transverse_dependence") = nb::none(),
         nb::arg("debug") = nb::none(),
         nb::arg("diagnostic_output_file") = nb::none()
  )
      .def_prop_rw(
          "ds_track_step",
          &SpaceChargeCommonStruct::ds_track_step,
          &SpaceChargeCommonStruct::set_ds_track_step,
          "CSR tracking step size"
      )
      .def_prop_rw(
          "dt_track_step",
          &SpaceChargeCommonStruct::dt_track_step,
          &SpaceChargeCommonStruct::set_dt_track_step,
          "Time Runge kutta initial step."
      )
      .def_prop_rw(
          "cathode_strength_cutoff",
          &SpaceChargeCommonStruct::cathode_strength_cutoff,
          &SpaceChargeCommonStruct::set_cathode_strength_cutoff,
          "Cutoff for the cathode field calc."
      )
      .def_prop_rw(
          "rel_tol_tracking",
          &SpaceChargeCommonStruct::rel_tol_tracking,
          &SpaceChargeCommonStruct::set_rel_tol_tracking,
          "Relative tolerance for tracking."
      )
      .def_prop_rw(
          "abs_tol_tracking",
          &SpaceChargeCommonStruct::abs_tol_tracking,
          &SpaceChargeCommonStruct::set_abs_tol_tracking,
          "Absolute tolerance for tracking."
      )
      .def_prop_rw(
          "beam_chamber_height",
          &SpaceChargeCommonStruct::beam_chamber_height,
          &SpaceChargeCommonStruct::set_beam_chamber_height,
          "Used in shielding calculation."
      )
      .def_prop_rw(
          "lsc_sigma_cutoff",
          &SpaceChargeCommonStruct::lsc_sigma_cutoff,
          &SpaceChargeCommonStruct::set_lsc_sigma_cutoff,
          "Cutoff for the 1-dim longitudinal SC calc. If a bin sigma is < cutoff * sigma_ave then "
          "ignore."
      )
      .def_prop_rw(
          "particle_sigma_cutoff",
          &SpaceChargeCommonStruct::particle_sigma_cutoff,
          &SpaceChargeCommonStruct::set_particle_sigma_cutoff,
          "3D SC calc cutoff for particles with (x,y,z) position far from the center. Negative or "
          "zero means ignore."
      )
      .def_prop_rw(
          "space_charge_mesh_size",
          &SpaceChargeCommonStruct::space_charge_mesh_size,
          &SpaceChargeCommonStruct::set_space_charge_mesh_size,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Gird size for fft_3d space charge calc."
      )
      .def_prop_rw(
          "csr3d_mesh_size",
          &SpaceChargeCommonStruct::csr3d_mesh_size,
          &SpaceChargeCommonStruct::set_csr3d_mesh_size,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Gird size for CSR."
      )
      .def_prop_rw(
          "n_bin",
          &SpaceChargeCommonStruct::n_bin,
          &SpaceChargeCommonStruct::set_n_bin,
          "Number of bins used"
      )
      .def_prop_rw(
          "particle_bin_span",
          &SpaceChargeCommonStruct::particle_bin_span,
          &SpaceChargeCommonStruct::set_particle_bin_span,
          "Longitudinal particle length / dz_bin"
      )
      .def_prop_rw(
          "n_shield_images",
          &SpaceChargeCommonStruct::n_shield_images,
          &SpaceChargeCommonStruct::set_n_shield_images,
          "Chamber wall shielding. 0 = no shielding."
      )
      .def_prop_rw(
          "sc_min_in_bin",
          &SpaceChargeCommonStruct::sc_min_in_bin,
          &SpaceChargeCommonStruct::set_sc_min_in_bin,
          "Minimum number of particles in a bin for sigmas to be valid."
      )
      .def_prop_rw(
          "lsc_kick_transverse_dependence",
          &SpaceChargeCommonStruct::lsc_kick_transverse_dependence,
          &SpaceChargeCommonStruct::set_lsc_kick_transverse_dependence
      )
      .def_prop_rw("debug", &SpaceChargeCommonStruct::debug, &SpaceChargeCommonStruct::set_debug)
      .def_prop_rw(
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
          [](const SpaceChargeCommonStruct &self, nb::dict &memo) {
            return SpaceChargeCommonStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const SpaceChargeCommonStruct &self, const SpaceChargeCommonStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_spin_axis_struct(nb::module_ &m, nb::class_<SpinAxisStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>>(),
         nb::arg("l") = nb::none(),
         nb::arg("n0") = nb::none(),
         nb::arg("m") = nb::none()
  )
      .def_prop_rw(
          "l",
          &SpinAxisStruct::l,
          &SpinAxisStruct::set_l,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Transverse axis."
      )
      .def_prop_rw(
          "n0",
          &SpinAxisStruct::n0,
          &SpinAxisStruct::set_n0,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Invariant spin axis on closed orbit."
      )
      .def_prop_rw(
          "m",
          &SpinAxisStruct::m,
          &SpinAxisStruct::set_m,
          nb::for_getter(nb::keep_alive<0, 1>()),
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
          [](const SpinAxisStruct &self, nb::dict &memo) { return SpinAxisStruct(self); }
      )
      .def(
          "__eq__",
          [](const SpinAxisStruct &self, const SpinAxisStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
// spin_eigen_struct
void init_spin_eigen_struct(nb::module_ &m, nb::class_<SpinEigenStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<std::complex<double>>>,
             std::optional<std::complex<double>>>(),
         nb::arg("vec") = nb::none(),
         nb::arg("val") = nb::none()
  )
      .def_prop_rw(
          "vec",
          &SpinEigenStruct::vec,
          &SpinEigenStruct::set_vec,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw("val", &SpinEigenStruct::val, &SpinEigenStruct::set_val)
      .def_static(
          "new_array1d",
          [](int sz) { return SpinEigenStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = SpinEigenStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const SpinEigenStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SpinEigenStruct &self) {
            return SpinEigenStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SpinEigenStruct &self, nb::dict &memo) { return SpinEigenStruct(self); }
      )
      .def(
          "__eq__",
          [](const SpinEigenStruct &self, const SpinEigenStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const SpinEigenStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<SpinEigenStructArray1D, SpinEigenStructAlloc1D>(
      m,
      "SpinEigenStructArray1D",
      "SpinEigenStructAlloc1D"
  );
  // 2D SpinEigenStruct arrays are not used in structs/routines
  // 3D SpinEigenStruct arrays are not used in structs/routines
}

// =============================================================================
// spin_matching_struct
void init_spin_matching_struct(nb::module_ &m, nb::class_<SpinMatchingStruct> &cls) {
  cls.def(
         "__init__",
         [](SpinMatchingStruct *self,
            const SpinAxisStruct *axis,
            std::optional<std::vector<double>> dn_dpz,
            std::optional<std::vector<double>> alpha,
            std::optional<std::vector<double>> beta,
            std::optional<std::vector<double>> orb0,
            std::optional<std::vector<std::vector<double>>> M_1turn,
            std::optional<std::vector<std::vector<double>>> M_ele,
            std::optional<std::vector<double>> sq_ele,
            std::optional<std::vector<double>> sq_1turn,
            std::optional<bool> valid) {
           new (self) SpinMatchingStruct(
               ptr_to_opt_ref(axis),
               dn_dpz,
               alpha,
               beta,
               orb0,
               M_1turn,
               M_ele,
               sq_ele,
               sq_1turn,
               valid
           );
         },
         nb::arg("axis") = nb::none(),
         nb::arg("dn_dpz") = nb::none(),
         nb::arg("alpha") = nb::none(),
         nb::arg("beta") = nb::none(),
         nb::arg("orb0") = nb::none(),
         nb::arg("M_1turn") = nb::none(),
         nb::arg("M_ele") = nb::none(),
         nb::arg("sq_ele") = nb::none(),
         nb::arg("sq_1turn") = nb::none(),
         nb::arg("valid") = nb::none()
  )
      .def_prop_rw(
          "axis",
          &SpinMatchingStruct::axis,
          &SpinMatchingStruct::set_axis,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro("eigen", &SpinMatchingStruct::eigen, nb::keep_alive<0, 1>())
      .def_prop_rw(
          "dn_dpz",
          &SpinMatchingStruct::dn_dpz,
          &SpinMatchingStruct::set_dn_dpz,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Invariant spin derivative"
      )
      .def_prop_rw(
          "alpha",
          &SpinMatchingStruct::alpha,
          &SpinMatchingStruct::set_alpha,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Alpha vector"
      )
      .def_prop_rw(
          "beta",
          &SpinMatchingStruct::beta,
          &SpinMatchingStruct::set_beta,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Beta vector"
      )
      .def_prop_rw(
          "orb0",
          &SpinMatchingStruct::orb0,
          &SpinMatchingStruct::set_orb0,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Closed orbit"
      )
      .def_prop_rw(
          "M_1turn",
          &SpinMatchingStruct::M_1turn,
          &SpinMatchingStruct::set_M_1turn,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "1-turn matrix"
      )
      .def_prop_rw(
          "M_ele",
          &SpinMatchingStruct::M_ele,
          &SpinMatchingStruct::set_M_ele,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Transfer matrix through element."
      )
      .def_prop_rw(
          "sq_ele",
          &SpinMatchingStruct::sq_ele,
          &SpinMatchingStruct::set_sq_ele,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "sq_1turn",
          &SpinMatchingStruct::sq_1turn,
          &SpinMatchingStruct::set_sq_1turn,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw("valid", &SpinMatchingStruct::valid, &SpinMatchingStruct::set_valid)
      .def_static(
          "new_array1d",
          [](int sz) { return SpinMatchingStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = SpinMatchingStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const SpinMatchingStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SpinMatchingStruct &self) {
            return SpinMatchingStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SpinMatchingStruct &self, nb::dict &memo) { return SpinMatchingStruct(self); }
      )
      .def(
          "__eq__",
          [](const SpinMatchingStruct &self, const SpinMatchingStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const SpinMatchingStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<SpinMatchingStructArray1D, SpinMatchingStructAlloc1D>(
      m,
      "SpinMatchingStructArray1D",
      "SpinMatchingStructAlloc1D"
  );
  // 2D SpinMatchingStruct arrays are not used in structs/routines
  // 3D SpinMatchingStruct arrays are not used in structs/routines
}

// =============================================================================
// spin_orbit_map1_struct
void init_spin_orbit_map1_struct(nb::module_ &m, nb::class_<SpinOrbitMap1Struct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<std::vector<double>>>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<std::vector<double>>>>(),
         nb::arg("orb_mat") = nb::none(),
         nb::arg("vec0") = nb::none(),
         nb::arg("spin_q") = nb::none()
  )
      .def_prop_rw(
          "orb_mat",
          &SpinOrbitMap1Struct::orb_mat,
          &SpinOrbitMap1Struct::set_orb_mat,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Orbital matrix"
      )
      .def_prop_rw(
          "vec0",
          &SpinOrbitMap1Struct::vec0,
          &SpinOrbitMap1Struct::set_vec0,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Orbital 0th order map: r_out = mat6 * r_in + vec0"
      )
      .def_prop_rw(
          "spin_q",
          &SpinOrbitMap1Struct::spin_q,
          &SpinOrbitMap1Struct::set_spin_q,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "0th and 1st order quaternion spin map"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return SpinOrbitMap1StructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = SpinOrbitMap1StructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const SpinOrbitMap1Struct &self, nb::dict &memo) { return SpinOrbitMap1Struct(self); }
      )
      .def(
          "__eq__",
          [](const SpinOrbitMap1Struct &self, const SpinOrbitMap1Struct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_spin_polar_struct(nb::module_ &m, nb::class_<SpinPolarStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("polarization") = nb::none(),
         nb::arg("theta") = nb::none(),
         nb::arg("phi") = nb::none(),
         nb::arg("xi") = nb::none()
  )
      .def_prop_rw(
          "polarization",
          &SpinPolarStruct::polarization,
          &SpinPolarStruct::set_polarization
      )
      .def_prop_rw(
          "theta",
          &SpinPolarStruct::theta,
          &SpinPolarStruct::set_theta,
          "Spherical coords: Angle from z-axis."
      )
      .def_prop_rw(
          "phi",
          &SpinPolarStruct::phi,
          &SpinPolarStruct::set_phi,
          "Spherical coords: Angle in (x,y) plane."
      )
      .def_prop_rw(
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
          [](const SpinPolarStruct &self, nb::dict &memo) { return SpinPolarStruct(self); }
      )
      .def(
          "__eq__",
          [](const SpinPolarStruct &self, const SpinPolarStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_strong_beam_struct(nb::module_ &m, nb::class_<StrongBeamStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<int>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("ix_slice") = nb::none(),
         nb::arg("x_center") = nb::none(),
         nb::arg("y_center") = nb::none(),
         nb::arg("x_sigma") = nb::none(),
         nb::arg("y_sigma") = nb::none(),
         nb::arg("dx") = nb::none(),
         nb::arg("dy") = nb::none()
  )
      .def_prop_rw(
          "ix_slice",
          &StrongBeamStruct::ix_slice,
          &StrongBeamStruct::set_ix_slice,
          "0 -> at element center and not at slice."
      )
      .def_prop_rw(
          "x_center",
          &StrongBeamStruct::x_center,
          &StrongBeamStruct::set_x_center,
          "Strong beam slice center."
      )
      .def_prop_rw(
          "y_center",
          &StrongBeamStruct::y_center,
          &StrongBeamStruct::set_y_center,
          "Strong beam slice center."
      )
      .def_prop_rw(
          "x_sigma",
          &StrongBeamStruct::x_sigma,
          &StrongBeamStruct::set_x_sigma,
          "Strong beam slice sigma."
      )
      .def_prop_rw(
          "y_sigma",
          &StrongBeamStruct::y_sigma,
          &StrongBeamStruct::set_y_sigma,
          "Strong beam slice sigma."
      )
      .def_prop_rw(
          "dx",
          &StrongBeamStruct::dx,
          &StrongBeamStruct::set_dx,
          "Particle - beam slice distance."
      )
      .def_prop_rw(
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
          [](const StrongBeamStruct &self, nb::dict &memo) { return StrongBeamStruct(self); }
      )
      .def(
          "__eq__",
          [](const StrongBeamStruct &self, const StrongBeamStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_surface_curvature_struct(nb::module_ &m, nb::class_<SurfaceCurvatureStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<std::vector<std::vector<double>>>,
             std::optional<double>,
             std::optional<std::vector<double>>,
             std::optional<bool>>(),
         nb::arg("xy") = nb::none(),
         nb::arg("spherical") = nb::none(),
         nb::arg("elliptical") = nb::none(),
         nb::arg("has_curvature") = nb::none()
  )
      .def_prop_rw(
          "xy",
          &SurfaceCurvatureStruct::xy,
          &SurfaceCurvatureStruct::set_xy,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "spherical",
          &SurfaceCurvatureStruct::spherical,
          &SurfaceCurvatureStruct::set_spherical
      )
      .def_prop_rw(
          "elliptical",
          &SurfaceCurvatureStruct::elliptical,
          &SurfaceCurvatureStruct::set_elliptical,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Total curvature = elliptical + spherical"
      )
      .def_prop_rw(
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
          [](const SurfaceCurvatureStruct &self, nb::dict &memo) {
            return SurfaceCurvatureStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const SurfaceCurvatureStruct &self, const SurfaceCurvatureStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
    nb::module_ &m,
    nb::class_<SurfaceDisplacementPtStruct> &cls
) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("x0") = nb::none(),
         nb::arg("y0") = nb::none(),
         nb::arg("z0") = nb::none(),
         nb::arg("dz_dx") = nb::none(),
         nb::arg("dz_dy") = nb::none(),
         nb::arg("d2z_dxdy") = nb::none()
  )
      .def_prop_rw(
          "x0",
          &SurfaceDisplacementPtStruct::x0,
          &SurfaceDisplacementPtStruct::set_x0,
          "Position at center"
      )
      .def_prop_rw(
          "y0",
          &SurfaceDisplacementPtStruct::y0,
          &SurfaceDisplacementPtStruct::set_y0,
          "Position at center"
      )
      .def_prop_rw("z0", &SurfaceDisplacementPtStruct::z0, &SurfaceDisplacementPtStruct::set_z0)
      .def_prop_rw(
          "dz_dx",
          &SurfaceDisplacementPtStruct::dz_dx,
          &SurfaceDisplacementPtStruct::set_dz_dx
      )
      .def_prop_rw(
          "dz_dy",
          &SurfaceDisplacementPtStruct::dz_dy,
          &SurfaceDisplacementPtStruct::set_dz_dy
      )
      .def_prop_rw(
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
          [](const SurfaceDisplacementPtStruct &self, nb::dict &memo) {
            return SurfaceDisplacementPtStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const SurfaceDisplacementPtStruct &self, const SurfaceDisplacementPtStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_surface_displacement_struct(nb::module_ &m, nb::class_<SurfaceDisplacementStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<bool>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>>(),
         nb::arg("active") = nb::none(),
         nb::arg("dr") = nb::none(),
         nb::arg("r0") = nb::none()
  )
      .def_prop_rw(
          "active",
          &SurfaceDisplacementStruct::active,
          &SurfaceDisplacementStruct::set_active
      )
      .def_prop_rw(
          "dr",
          &SurfaceDisplacementStruct::dr,
          &SurfaceDisplacementStruct::set_dr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "r0",
          &SurfaceDisplacementStruct::r0,
          &SurfaceDisplacementStruct::set_r0,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro("pt", &SurfaceDisplacementStruct::pt, nb::keep_alive<0, 1>())

      .def("__repr__", [](const SurfaceDisplacementStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceDisplacementStruct &self) {
            return SurfaceDisplacementStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SurfaceDisplacementStruct &self, nb::dict &memo) {
            return SurfaceDisplacementStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const SurfaceDisplacementStruct &self, const SurfaceDisplacementStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_surface_h_misalign_pt_struct(nb::module_ &m, nb::class_<SurfaceHMisalignPtStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("x0") = nb::none(),
         nb::arg("y0") = nb::none(),
         nb::arg("rot_y") = nb::none(),
         nb::arg("rot_t") = nb::none(),
         nb::arg("rot_y_rms") = nb::none(),
         nb::arg("rot_t_rms") = nb::none()
  )
      .def_prop_rw(
          "x0",
          &SurfaceHMisalignPtStruct::x0,
          &SurfaceHMisalignPtStruct::set_x0,
          "Position at center"
      )
      .def_prop_rw(
          "y0",
          &SurfaceHMisalignPtStruct::y0,
          &SurfaceHMisalignPtStruct::set_y0,
          "Position at center"
      )
      .def_prop_rw(
          "rot_y",
          &SurfaceHMisalignPtStruct::rot_y,
          &SurfaceHMisalignPtStruct::set_rot_y,
          "rot_t = x-rotation for Bragg and z-rotation for Laue."
      )
      .def_prop_rw(
          "rot_t",
          &SurfaceHMisalignPtStruct::rot_t,
          &SurfaceHMisalignPtStruct::set_rot_t,
          "rot_t = x-rotation for Bragg and z-rotation for Laue."
      )
      .def_prop_rw(
          "rot_y_rms",
          &SurfaceHMisalignPtStruct::rot_y_rms,
          &SurfaceHMisalignPtStruct::set_rot_y_rms,
          "rot_t = x-rotation for Bragg and z-rotation for Laue."
      )
      .def_prop_rw(
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
          [](const SurfaceHMisalignPtStruct &self, nb::dict &memo) {
            return SurfaceHMisalignPtStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const SurfaceHMisalignPtStruct &self, const SurfaceHMisalignPtStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_surface_h_misalign_struct(nb::module_ &m, nb::class_<SurfaceHMisalignStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<bool>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>>(),
         nb::arg("active") = nb::none(),
         nb::arg("dr") = nb::none(),
         nb::arg("r0") = nb::none()
  )
      .def_prop_rw("active", &SurfaceHMisalignStruct::active, &SurfaceHMisalignStruct::set_active)
      .def_prop_rw(
          "dr",
          &SurfaceHMisalignStruct::dr,
          &SurfaceHMisalignStruct::set_dr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "r0",
          &SurfaceHMisalignStruct::r0,
          &SurfaceHMisalignStruct::set_r0,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro("pt", &SurfaceHMisalignStruct::pt, nb::keep_alive<0, 1>())

      .def("__repr__", [](const SurfaceHMisalignStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceHMisalignStruct &self) {
            return SurfaceHMisalignStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SurfaceHMisalignStruct &self, nb::dict &memo) {
            return SurfaceHMisalignStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const SurfaceHMisalignStruct &self, const SurfaceHMisalignStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_surface_segmented_pt_struct(nb::module_ &m, nb::class_<SurfaceSegmentedPtStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<double>>(),
         nb::arg("x0") = nb::none(),
         nb::arg("y0") = nb::none(),
         nb::arg("z0") = nb::none(),
         nb::arg("dz_dx") = nb::none(),
         nb::arg("dz_dy") = nb::none()
  )
      .def_prop_rw(
          "x0",
          &SurfaceSegmentedPtStruct::x0,
          &SurfaceSegmentedPtStruct::set_x0,
          "Position at center"
      )
      .def_prop_rw(
          "y0",
          &SurfaceSegmentedPtStruct::y0,
          &SurfaceSegmentedPtStruct::set_y0,
          "Position at center"
      )
      .def_prop_rw(
          "z0",
          &SurfaceSegmentedPtStruct::z0,
          &SurfaceSegmentedPtStruct::set_z0,
          "Position at center"
      )
      .def_prop_rw(
          "dz_dx",
          &SurfaceSegmentedPtStruct::dz_dx,
          &SurfaceSegmentedPtStruct::set_dz_dx,
          "Slope at center"
      )
      .def_prop_rw(
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
          [](const SurfaceSegmentedPtStruct &self, nb::dict &memo) {
            return SurfaceSegmentedPtStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const SurfaceSegmentedPtStruct &self, const SurfaceSegmentedPtStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_surface_segmented_struct(nb::module_ &m, nb::class_<SurfaceSegmentedStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<bool>,
             std::optional<std::vector<double>>,
             std::optional<std::vector<double>>>(),
         nb::arg("active") = nb::none(),
         nb::arg("dr") = nb::none(),
         nb::arg("r0") = nb::none()
  )
      .def_prop_rw("active", &SurfaceSegmentedStruct::active, &SurfaceSegmentedStruct::set_active)
      .def_prop_rw(
          "dr",
          &SurfaceSegmentedStruct::dr,
          &SurfaceSegmentedStruct::set_dr,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_rw(
          "r0",
          &SurfaceSegmentedStruct::r0,
          &SurfaceSegmentedStruct::set_r0,
          nb::for_getter(nb::keep_alive<0, 1>())
      )
      .def_prop_ro("pt", &SurfaceSegmentedStruct::pt, nb::keep_alive<0, 1>())

      .def("__repr__", [](const SurfaceSegmentedStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const SurfaceSegmentedStruct &self) {
            return SurfaceSegmentedStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const SurfaceSegmentedStruct &self, nb::dict &memo) {
            return SurfaceSegmentedStruct(self);
          }
      )
      .def(
          "__eq__",
          [](const SurfaceSegmentedStruct &self, const SurfaceSegmentedStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
// str_index_struct
void init_str_index_struct(nb::module_ &m, nb::class_<StrIndexStruct> &cls) {
  cls.def(
         nb::init<std::optional<std::vector<int>>, std::optional<int>, std::optional<int>>(),
         nb::arg("index") = nb::none(),
         nb::arg("n_min") = nb::none(),
         nb::arg("n_max") = nb::none()
  )
      .def_prop_ro("name", &StrIndexStruct::name, nb::keep_alive<0, 1>(), "Array of names.")
      .def_prop_rw(
          "index",
          &StrIndexStruct::index,
          &StrIndexStruct::set_index,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "Sorted index for names(:) array. names(an_index(i)) is in alphabetical order."
      )
      .def_prop_rw("n_min", &StrIndexStruct::n_min, &StrIndexStruct::set_n_min)
      .def_prop_rw(
          "n_max",
          &StrIndexStruct::n_max,
          &StrIndexStruct::set_n_max,
          "Use only names(n_min:n_max) part of array."
      )

      .def("__repr__", [](const StrIndexStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const StrIndexStruct &self) {
            return StrIndexStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const StrIndexStruct &self, nb::dict &memo) { return StrIndexStruct(self); }
      )
      .def(
          "__eq__",
          [](const StrIndexStruct &self, const StrIndexStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const StrIndexStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D StrIndexStruct arrays are not used in structs/routines
  // 2D StrIndexStruct arrays are not used in structs/routines
  // 3D StrIndexStruct arrays are not used in structs/routines
}

// =============================================================================
// spline_struct
void init_spline_struct(nb::module_ &m, nb::class_<SplineStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<double>,
             std::optional<double>,
             std::optional<std::vector<double>>>(),
         nb::arg("x0") = nb::none(),
         nb::arg("y0") = nb::none(),
         nb::arg("x1") = nb::none(),
         nb::arg("coef") = nb::none()
  )
      .def_prop_rw("x0", &SplineStruct::x0, &SplineStruct::set_x0, "Point at start of spline")
      .def_prop_rw("y0", &SplineStruct::y0, &SplineStruct::set_y0, "Point at start of spline")
      .def_prop_rw("x1", &SplineStruct::x1, &SplineStruct::set_x1, "Point at end of spline")
      .def_prop_rw(
          "coef",
          &SplineStruct::coef,
          &SplineStruct::set_coef,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "coefficients for cubic spline"
      )
      .def_static(
          "new_array1d",
          [](int sz) { return SplineStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = SplineStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const SplineStruct &self, nb::dict &memo) { return SplineStruct(self); }
      )
      .def(
          "__eq__",
          [](const SplineStruct &self, const SplineStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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
void init_summation_rdt_struct(nb::module_ &m, nb::class_<SummationRdtStruct> &cls) {
  cls.def(
         nb::init<
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
         nb::arg("h11001") = nb::none(),
         nb::arg("h00111") = nb::none(),
         nb::arg("h20001") = nb::none(),
         nb::arg("h00201") = nb::none(),
         nb::arg("h10002") = nb::none(),
         nb::arg("h21000") = nb::none(),
         nb::arg("h30000") = nb::none(),
         nb::arg("h10110") = nb::none(),
         nb::arg("h10020") = nb::none(),
         nb::arg("h10200") = nb::none(),
         nb::arg("h31000") = nb::none(),
         nb::arg("h40000") = nb::none(),
         nb::arg("h20110") = nb::none(),
         nb::arg("h11200") = nb::none(),
         nb::arg("h20020") = nb::none(),
         nb::arg("h20200") = nb::none(),
         nb::arg("h00310") = nb::none(),
         nb::arg("h00400") = nb::none(),
         nb::arg("h22000") = nb::none(),
         nb::arg("h00220") = nb::none(),
         nb::arg("h11110") = nb::none()
  )
      .def_prop_rw("h11001", &SummationRdtStruct::h11001, &SummationRdtStruct::set_h11001)
      .def_prop_rw("h00111", &SummationRdtStruct::h00111, &SummationRdtStruct::set_h00111)
      .def_prop_rw("h20001", &SummationRdtStruct::h20001, &SummationRdtStruct::set_h20001)
      .def_prop_rw("h00201", &SummationRdtStruct::h00201, &SummationRdtStruct::set_h00201)
      .def_prop_rw("h10002", &SummationRdtStruct::h10002, &SummationRdtStruct::set_h10002)
      .def_prop_rw("h21000", &SummationRdtStruct::h21000, &SummationRdtStruct::set_h21000)
      .def_prop_rw("h30000", &SummationRdtStruct::h30000, &SummationRdtStruct::set_h30000)
      .def_prop_rw("h10110", &SummationRdtStruct::h10110, &SummationRdtStruct::set_h10110)
      .def_prop_rw("h10020", &SummationRdtStruct::h10020, &SummationRdtStruct::set_h10020)
      .def_prop_rw(
          "h10200",
          &SummationRdtStruct::h10200,
          &SummationRdtStruct::set_h10200,
          "2nd order in K2 moments"
      )
      .def_prop_rw("h31000", &SummationRdtStruct::h31000, &SummationRdtStruct::set_h31000)
      .def_prop_rw("h40000", &SummationRdtStruct::h40000, &SummationRdtStruct::set_h40000)
      .def_prop_rw("h20110", &SummationRdtStruct::h20110, &SummationRdtStruct::set_h20110)
      .def_prop_rw("h11200", &SummationRdtStruct::h11200, &SummationRdtStruct::set_h11200)
      .def_prop_rw("h20020", &SummationRdtStruct::h20020, &SummationRdtStruct::set_h20020)
      .def_prop_rw("h20200", &SummationRdtStruct::h20200, &SummationRdtStruct::set_h20200)
      .def_prop_rw("h00310", &SummationRdtStruct::h00310, &SummationRdtStruct::set_h00310)
      .def_prop_rw("h00400", &SummationRdtStruct::h00400, &SummationRdtStruct::set_h00400)
      .def_prop_rw("h22000", &SummationRdtStruct::h22000, &SummationRdtStruct::set_h22000)
      .def_prop_rw("h00220", &SummationRdtStruct::h00220, &SummationRdtStruct::set_h00220)
      .def_prop_rw("h11110", &SummationRdtStruct::h11110, &SummationRdtStruct::set_h11110)
      .def_static(
          "new_array1d",
          [](int sz) { return SummationRdtStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = SummationRdtStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
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
          [](const SummationRdtStruct &self, nb::dict &memo) { return SummationRdtStruct(self); }
      )
      .def(
          "__eq__",
          [](const SummationRdtStruct &self, const SummationRdtStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
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