#include "pybmad/generated/structs_i.hpp"

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
// interval1_coef_struct
void init_interval1_coef_struct(nb::module_ &m, nb::class_<Interval1CoefStruct> &cls) {
  cls.def(
         nb::init<std::optional<double>, std::optional<double>, std::optional<double>>(),
         nb::arg("c0") = nb::none(),
         nb::arg("c1") = nb::none(),
         nb::arg("n_exp") = nb::none()
  )
      .def_prop_rw("c0", &Interval1CoefStruct::c0, &Interval1CoefStruct::set_c0)
      .def_prop_rw("c1", &Interval1CoefStruct::c1, &Interval1CoefStruct::set_c1)
      .def_prop_rw("n_exp", &Interval1CoefStruct::n_exp, &Interval1CoefStruct::set_n_exp)
      .def_static(
          "new_array1d",
          [](int sz) { return Interval1CoefStructAlloc1D(sz); },
          nb::arg("sz") = 0
      )
      .def_static(
          "new_array1d_bounds",
          [](int lbound, int ubound) {
            auto cnt = Interval1CoefStructAlloc1D();
            cnt.resize_bounds(lbound, ubound);
            return cnt;
          },
          nb::arg("lbound"),
          nb::arg("ubound")
      )

      .def("__repr__", [](const Interval1CoefStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Interval1CoefStruct &self) {
            return Interval1CoefStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const Interval1CoefStruct &self, nb::dict &memo) { return Interval1CoefStruct(self); }
      )
      .def(
          "__eq__",
          [](const Interval1CoefStruct &self, const Interval1CoefStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const Interval1CoefStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  bind_1d_type_array_pair<Interval1CoefStructArray1D, Interval1CoefStructAlloc1D>(
      m,
      "Interval1CoefStructArray1D",
      "Interval1CoefStructAlloc1D"
  );
  // 2D Interval1CoefStruct arrays are not used in structs/routines
  // 3D Interval1CoefStruct arrays are not used in structs/routines
}

// =============================================================================
// ibs_lifetime_struct
void init_ibs_lifetime_struct(nb::module_ &m, nb::class_<IbsLifetimeStruct> &cls) {
  cls.def(
         nb::init<std::optional<double>, std::optional<double>, std::optional<double>>(),
         nb::arg("Tlx") = nb::none(),
         nb::arg("Tly") = nb::none(),
         nb::arg("Tlp") = nb::none()
  )
      .def_prop_rw("Tlx", &IbsLifetimeStruct::Tlx, &IbsLifetimeStruct::set_Tlx)
      .def_prop_rw("Tly", &IbsLifetimeStruct::Tly, &IbsLifetimeStruct::set_Tly)
      .def_prop_rw("Tlp", &IbsLifetimeStruct::Tlp, &IbsLifetimeStruct::set_Tlp)

      .def("__repr__", [](const IbsLifetimeStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const IbsLifetimeStruct &self) {
            return IbsLifetimeStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const IbsLifetimeStruct &self, nb::dict &memo) { return IbsLifetimeStruct(self); }
      )
      .def(
          "__eq__",
          [](const IbsLifetimeStruct &self, const IbsLifetimeStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const IbsLifetimeStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D IbsLifetimeStruct arrays are not used in structs/routines
  // 2D IbsLifetimeStruct arrays are not used in structs/routines
  // 3D IbsLifetimeStruct arrays are not used in structs/routines
}

// =============================================================================
// ibs_maxratio_struct
void init_ibs_maxratio_struct(nb::module_ &m, nb::class_<IbsMaxratioStruct> &cls) {
  cls.def(
         nb::init<std::optional<double>, std::optional<double>, std::optional<double>>(),
         nb::arg("rx") = nb::none(),
         nb::arg("ry") = nb::none(),
         nb::arg("r_p") = nb::none()
  )
      .def_prop_rw("rx", &IbsMaxratioStruct::rx, &IbsMaxratioStruct::set_rx)
      .def_prop_rw("ry", &IbsMaxratioStruct::ry, &IbsMaxratioStruct::set_ry)
      .def_prop_rw("r_p", &IbsMaxratioStruct::r_p, &IbsMaxratioStruct::set_r_p)

      .def("__repr__", [](const IbsMaxratioStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const IbsMaxratioStruct &self) {
            return IbsMaxratioStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const IbsMaxratioStruct &self, nb::dict &memo) { return IbsMaxratioStruct(self); }
      )
      .def(
          "__eq__",
          [](const IbsMaxratioStruct &self, const IbsMaxratioStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const IbsMaxratioStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D IbsMaxratioStruct arrays are not used in structs/routines
  // 2D IbsMaxratioStruct arrays are not used in structs/routines
  // 3D IbsMaxratioStruct arrays are not used in structs/routines
}

// =============================================================================
// ibs_sim_param_struct
void init_ibs_sim_param_struct(nb::module_ &m, nb::class_<IbsSimParamStruct> &cls) {
  cls.def(
         nb::init<
             std::optional<double>,
             std::optional<int>,
             std::optional<bool>,
             std::optional<double>,
             std::optional<double>,
             std::optional<bool>,
             std::optional<double>,
             std::optional<std::string>>(),
         nb::arg("tau_a") = nb::none(),
         nb::arg("clog_to_use") = nb::none(),
         nb::arg("set_dispersion") = nb::none(),
         nb::arg("eta_set") = nb::none(),
         nb::arg("etap_set") = nb::none(),
         nb::arg("do_pwd") = nb::none(),
         nb::arg("inductance") = nb::none(),
         nb::arg("formula") = nb::none()
  )
      .def_prop_rw(
          "tau_a",
          &IbsSimParamStruct::tau_a,
          &IbsSimParamStruct::set_tau_a,
          "horizontal damping rate (needed for coulomb log tail cut)"
      )
      .def_prop_rw(
          "clog_to_use",
          &IbsSimParamStruct::clog_to_use,
          &IbsSimParamStruct::set_clog_to_use,
          "see multi_coulomb_log subroutine for valid settings.  Set to 1 to disable tail-cut.  "
          "Set to 1 for linacs."
      )
      .def_prop_rw(
          "set_dispersion",
          &IbsSimParamStruct::set_dispersion,
          &IbsSimParamStruct::set_set_dispersion,
          "True: add vertical dispersion to transfer matrix.  Valid for kubo method."
      )
      .def_prop_rw(
          "eta_set",
          &IbsSimParamStruct::eta_set,
          &IbsSimParamStruct::set_eta_set,
          "If set_dispersion, then this value is used to add y-z coupling to the transfer matrix."
      )
      .def_prop_rw(
          "etap_set",
          &IbsSimParamStruct::etap_set,
          &IbsSimParamStruct::set_etap_set,
          "If set_dispersion, then this value is used to add y-z coupling to the transfer matrix."
      )
      .def_prop_rw(
          "do_pwd",
          &IbsSimParamStruct::do_pwd,
          &IbsSimParamStruct::set_do_pwd,
          "If true, then use potential well distortion to calculate bunch lengths.  If false, "
          "bunch length is proportional to energy spread."
      )
      .def_prop_rw(
          "inductance",
          &IbsSimParamStruct::inductance,
          &IbsSimParamStruct::set_inductance,
          "Inductive part of impedance for pwd calc."
      )
      .def_prop_rw(
          "formula",
          &IbsSimParamStruct::formula,
          &IbsSimParamStruct::set_formula,
          "Which IBS formulation to use.  See subroutine ibs1 for a list. real(rp) :: fake_3HC = "
          "-1   ! If greater than zero, divide growth rates by this factor."
      )

      .def("__repr__", [](const IbsSimParamStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const IbsSimParamStruct &self) {
            return IbsSimParamStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const IbsSimParamStruct &self, nb::dict &memo) { return IbsSimParamStruct(self); }
      )
      .def(
          "__eq__",
          [](const IbsSimParamStruct &self, const IbsSimParamStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const IbsSimParamStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D IbsSimParamStruct arrays are not used in structs/routines
  // 2D IbsSimParamStruct arrays are not used in structs/routines
  // 3D IbsSimParamStruct arrays are not used in structs/routines
}

// =============================================================================
// ibs_struct
void init_ibs_struct(nb::module_ &m, nb::class_<IbsStruct> &cls) {
  cls.def(
         nb::init<std::optional<double>, std::optional<double>, std::optional<double>>(),
         nb::arg("inv_Ta") = nb::none(),
         nb::arg("inv_Tb") = nb::none(),
         nb::arg("inv_Tz") = nb::none()
  )
      .def_prop_rw("inv_Ta", &IbsStruct::inv_Ta, &IbsStruct::set_inv_Ta)
      .def_prop_rw("inv_Tb", &IbsStruct::inv_Tb, &IbsStruct::set_inv_Tb)
      .def_prop_rw("inv_Tz", &IbsStruct::inv_Tz, &IbsStruct::set_inv_Tz)

      .def("__repr__", [](const IbsStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const IbsStruct &self) {
            return IbsStruct(self); // under-the-hood fortran copy
          }
      )
      .def("__deepcopy__", [](const IbsStruct &self, nb::dict &memo) { return IbsStruct(self); })
      .def(
          "__eq__",
          [](const IbsStruct &self, const IbsStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const IbsStruct &self) {
            return std::hash<std::uintptr_t>{}(
                reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr())
            );
          }
      )

      ;

  // 1D IbsStruct arrays are not used in structs/routines
  // 2D IbsStruct arrays are not used in structs/routines
  // 3D IbsStruct arrays are not used in structs/routines
}