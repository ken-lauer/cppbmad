// Hand-written nanobind bindings exposing Tao's hooks to Python as read/write
// properties on a singleton holder: `tao.hooks.<name> = callable`,
// `tao.hooks.<name> = None` to clear, and reading `tao.hooks.<name>` returns the
// current callable or None.
//
// Each property setter wraps the Python callable in a std::function that, on
// every invocation, acquires the GIL, hands the Fortran-owned proxies to Python
// as non-owning references (rv_policy::reference), calls the Python callable, and
// translates its return value into the hook's out-parameter. Python exceptions
// are caught and reported here so they never unwind into the C trampoline (and
// thus never into Fortran). See plans/hook-plan.md.
//
// The proxies passed to the callback are borrowed and valid only for the
// duration of the call; do not stash them.

#include "pybmad/tao_hooks.hpp"

#include <nanobind/stl/optional.h>

#include "pybmad/hook_util.hpp"

using Pybmad::hooks::obj_to_bool;
using Pybmad::hooks::ref;
using Pybmad::hooks::report_error;

namespace {

// Holds the current Python callable for each hook (for readback via the getter).
struct TaoHooks {
  nb::object lattice_calc = nb::none();
  nb::object optimizer = nb::none();
  nb::object merit_var = nb::none();
  nb::object merit_data = nb::none();
};

} // namespace

void Pybmad::init_tao_hooks(nb::module_ &m) {
  nb::class_<TaoHooks> cls(m, "TaoHooks", "Registry of Tao hook callbacks.");

  cls.def_prop_rw(
      "lattice_calc",
      [](TaoHooks &h) { return h.lattice_calc; },
      [](TaoHooks &h, std::optional<nb::object> value) {
        nb::object fn = value ? *value : nb::none();
        h.lattice_calc = fn;
        if (fn.is_none()) {
          Tao::clear_lattice_calc_hook();
          return;
        }
        Tao::set_lattice_calc_hook([fn](bool calc_ok) -> bool {
          nb::gil_scoped_acquire acq;
          try {
            return obj_to_bool(fn(calc_ok));
          } catch (nb::python_error &e) {
            report_error(e);
            return calc_ok;
          }
        });
      }
  );

  cls.def_prop_rw(
      "optimizer",
      [](TaoHooks &h) { return h.optimizer; },
      [](TaoHooks &h, std::optional<nb::object> value) {
        nb::object fn = value ? *value : nb::none();
        h.optimizer = fn;
        if (fn.is_none()) {
          Tao::clear_optimizer_hook();
          return;
        }
        Tao::set_optimizer_hook([fn](bool abort) -> bool {
          nb::gil_scoped_acquire acq;
          try {
            return obj_to_bool(fn(abort));
          } catch (nb::python_error &e) {
            report_error(e);
            return abort;
          }
        });
      }
  );

  cls.def_prop_rw(
      "merit_var",
      [](TaoHooks &h) { return h.merit_var; },
      [](TaoHooks &h, std::optional<nb::object> value) {
        nb::object fn = value ? *value : nb::none();
        h.merit_var = fn;
        if (fn.is_none()) {
          Tao::clear_merit_var_hook();
          return;
        }
        Tao::set_merit_var_hook([fn](int i_uni, int j_var, Bmad::TaoVarStruct &var) {
          nb::gil_scoped_acquire acq;
          try {
            fn(i_uni, j_var, ref(var));
          } catch (nb::python_error &e) {
            report_error(e);
          }
        });
      }
  );

  cls.def_prop_rw(
      "merit_data",
      [](TaoHooks &h) { return h.merit_data; },
      [](TaoHooks &h, std::optional<nb::object> value) {
        nb::object fn = value ? *value : nb::none();
        h.merit_data = fn;
        if (fn.is_none()) {
          Tao::clear_merit_data_hook();
          return;
        }
        Tao::set_merit_data_hook(
            [fn](int i_uni, int j_data, Bmad::TaoDataStruct &datum, bool valid_value_set) -> bool {
              nb::gil_scoped_acquire acq;
              try {
                return obj_to_bool(fn(i_uni, j_data, ref(datum), valid_value_set));
              } catch (nb::python_error &e) {
                report_error(e);
                return valid_value_set;
              }
            }
        );
      }
  );

  m.attr("hooks") = nb::cast(TaoHooks{});
}
