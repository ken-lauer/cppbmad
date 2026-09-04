#include "pybmad/generated/Tao_routines_d.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_Tao_routines_d(nb::module_ &m) {
  m.def(
      "deallocate_node_components",
      &Tao::deallocate_node_components,
      nb::arg("node"),
      R"""(Routine to deallocate the allocatable components of a tree node.

Note: This is needed since gfortran does not deallocate the allocatable components of the
elements of a pointer array when the array itself is deallocated. Without this there is a
memory leak every time an expression is evaluated.

Parameters
----------
node : TaoEvalNodeStruct
    Node to clean up.
    This parameter is an input/output and is modified in-place.
    As an output, node: Node with allocatable components deallocated.
)"""
  );
}
