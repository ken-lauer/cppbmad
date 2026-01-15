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
  cls.def(py::init<>())
      // KvBeamInitStruct.part_per_phi (1D_NOT_integer - number of particles per angle variable.
      .def_property_readonly("part_per_phi", &KvBeamInitStruct::part_per_phi)
      // KvBeamInitStruct.n_I2 (0D_NOT_integer - number of I2
      .def_property("n_I2", &KvBeamInitStruct::n_I2, &KvBeamInitStruct::set_n_I2)
      // KvBeamInitStruct.A (0D_NOT_real - A = I1/e
      .def_property("A", &KvBeamInitStruct::A, &KvBeamInitStruct::set_A)

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