#include "pybmad/generated/structs_k.hpp"

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
// kv_beam_init_struct
void init_kv_beam_init_struct(nb::module_ &m, nb::class_<KvBeamInitStruct> &cls) {
  cls.def(
         nb::init<std::optional<std::vector<int>>, std::optional<int>, std::optional<double>>(),
         nb::arg("part_per_phi") = nb::none(),
         nb::arg("n_I2") = nb::none(),
         nb::arg("A") = nb::none()
  )
      .def_prop_rw(
          "part_per_phi",
          &KvBeamInitStruct::part_per_phi,
          &KvBeamInitStruct::set_part_per_phi,
          nb::for_getter(nb::keep_alive<0, 1>()),
          "number of particles per angle variable."
      )
      .def_prop_rw("n_I2", &KvBeamInitStruct::n_I2, &KvBeamInitStruct::set_n_I2, "number of I2")
      .def_prop_rw("A", &KvBeamInitStruct::A, &KvBeamInitStruct::set_A, "A = I1/e")

      .def("__repr__", [](const KvBeamInitStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const KvBeamInitStruct &self) {
            return KvBeamInitStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const KvBeamInitStruct &self, nb::dict &memo) { return KvBeamInitStruct(self); }
      )
      .def(
          "__eq__",
          [](const KvBeamInitStruct &self, const KvBeamInitStruct &other) {
            return self.get_fortran_ptr() == other.get_fortran_ptr();
          },
          nb::is_operator()
      )
      .def(
          "__hash__",
          [](const KvBeamInitStruct &self) {
            return std::hash<std::uintptr_t>{
            }(reinterpret_cast<std::uintptr_t>(self.get_fortran_ptr()));
          }
      )

      ;

  // 1D KvBeamInitStruct arrays are not used in structs/routines
  // 2D KvBeamInitStruct arrays are not used in structs/routines
  // 3D KvBeamInitStruct arrays are not used in structs/routines
}