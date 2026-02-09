#include "pybmad/generated/structs_k.hpp"

#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// kv_beam_init_struct
void init_kv_beam_init_struct(py::module &m, py::class_<KvBeamInitStruct> &cls) {
  cls.def(
         py::init<optional_ref<const std::vector<int>>, std::optional<int>, std::optional<double>>(
         ),
         py::arg("part_per_phi") = py::none(),
         py::arg("n_I2") = py::none(),
         py::arg("A") = py::none()
  )
      .def_property(
          "part_per_phi",
          &KvBeamInitStruct::part_per_phi,
          &KvBeamInitStruct::set_part_per_phi,
          "number of particles per angle variable."
      )
      .def_property("n_I2", &KvBeamInitStruct::n_I2, &KvBeamInitStruct::set_n_I2, "number of I2")
      .def_property("A", &KvBeamInitStruct::A, &KvBeamInitStruct::set_A, "A = I1/e")

      .def("__repr__", [](const KvBeamInitStruct &self) { return to_string(self); })

      .def(
          "__copy__",
          [](const KvBeamInitStruct &self) {
            return KvBeamInitStruct(self); // under-the-hood fortran copy
          }
      )
      .def(
          "__deepcopy__",
          [](const KvBeamInitStruct &self, py::dict &memo) { return KvBeamInitStruct(self); }
      )

      ;

  // 1D KvBeamInitStruct arrays are not used in structs/routines
  // 2D KvBeamInitStruct arrays are not used in structs/routines
  // 3D KvBeamInitStruct arrays are not used in structs/routines
}