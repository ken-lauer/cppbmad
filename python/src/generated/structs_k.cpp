#include "pybmad/generated/structs_k.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// kv_beam_init_struct
void init_kv_beam_init_struct(py::module& m, py::class_<KvBeamInitProxy>& cls) {
  cls.def(py::init<>())
      // KvBeamInitProxy.part_per_phi (1D_NOT_integer - number of particles per angle variable.
      .def_property_readonly("part_per_phi", &KvBeamInitProxy::part_per_phi)
      // KvBeamInitProxy.n_I2 (0D_NOT_integer - number of I2
      .def_property("n_I2", &KvBeamInitProxy::n_I2, &KvBeamInitProxy::set_n_I2)
      // KvBeamInitProxy.A (0D_NOT_real - A = I1/e
      .def_property("A", &KvBeamInitProxy::A, &KvBeamInitProxy::set_A)

      .def(
          "__repr__",
          [](const KvBeamInitProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const KvBeamInitProxy& self) {
            return KvBeamInitProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const KvBeamInitProxy& self, py::dict& memo) {
            return KvBeamInitProxy(self);
          })

      ;

  // 1D KvBeamInitProxy arrays are not used in structs/routines
  // 2D KvBeamInitProxy arrays are not used in structs/routines
  // 3D KvBeamInitProxy arrays are not used in structs/routines
}