#include "pybmad/generated/CppBmadExtra_routines_s.hpp"

namespace nb = nanobind;
using namespace nanobind::literals;
using namespace Pybmad;

void init_CppBmadExtra_routines_s(nb::module_ &m) {
  m.def(
      "set_ele_misalignments",
      &CppBmadExtra::set_ele_misalignments,
      nb::arg("ele"),
      nb::arg("x_offset"),
      nb::arg("y_offset"),
      nb::arg("z_offset"),
      nb::arg("x_pitch"),
      nb::arg("y_pitch"),
      nb::arg("tilt"),
      nb::arg("check_free") = nb::none(),
      R"""(Set the six misalignment attributes on ele in one call and mark bookkeeping
flags once at the end. Equivalent to six set_ele_attribute calls minus the
string parsing and minus the per-call attribute_set_bookkeeping work.
Caller must invoke lattice_bookkeeper before tracking.

The trailing underscore in the name marks this as a pybmad-side helper, not
a routine from Bmad proper.

When check_free is true (default), all six attributes are verified to
resolve via attribute_index() to their expected slots AND to satisfy
attribute_free(..., dependent_attribs_free = .true.). The check is
all-or-nothing: if any attribute fails, the function returns .false. and ele
is not modified. With check_free = .false. the writes proceed unconditionally;
use only when the ele key is known to support all six attributes.

Parameters
----------
ele : EleStruct
    Element to modify.
    This parameter is an input/output and is modified in-place.
    As an output, ele: Values stored and bookkeeping flags marked stale

x_offset : float
    Values to write.  Tilt is 'roll' only for bends.

y_offset : float
    Values to write.  Tilt is 'roll' only for bends.

z_offset : float
    Values to write.  Tilt is 'roll' only for bends.

x_pitch : float
    Values to write.  Tilt is 'roll' only for bends.

y_pitch : float
    Values to write.  Tilt is 'roll' only for bends.

tilt : float
    Values to write.  Tilt is 'roll' only for bends.

check_free : bool, optional
    Default .true.. If true, validate slot layout and attribute freeness before any write.

Returns
-------
ok : bool
    .true. on success, .false. if check_free caught a non-standard or non-free attribute.
)"""
  );
}
