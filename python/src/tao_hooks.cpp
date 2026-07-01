// Hand-written nanobind bindings exposing Tao's hooks to Python as
// tao.set_<x>_hook(callable) / tao.clear_<x>_hook().
//
// Each set_* wraps the Python callable in a std::function that, on every
// invocation, acquires the GIL, hands the Fortran-owned proxies to Python as
// non-owning references (rv_policy::reference), calls the Python callable, and
// translates its return value into the hook's out-parameter. Python exceptions
// are caught and reported here so they never unwind into the C trampoline (and
// thus never into Fortran). See plans/hook-plan.md.
//
// The proxies passed to the callback are borrowed and valid only for the
// duration of the call; do not stash them.

#include "pybmad/tao_hooks.hpp"

#include "pybmad/hook_util.hpp"

using Pybmad::hooks::obj_to_bool;
using Pybmad::hooks::ref;
using Pybmad::hooks::report_error;

void Pybmad::init_tao_hooks(nb::module_ &m) {
  m.def("set_lattice_calc_hook", [](nb::callable fn) {
    Tao::set_lattice_calc_hook([fn](bool calc_ok) -> bool {
      nb::gil_scoped_acquire acq;
      try {
        return obj_to_bool(fn(calc_ok));
      } catch (nb::python_error &e) {
        report_error(e);
        return calc_ok;
      }
    });
  });
  m.def("clear_lattice_calc_hook", []() { Tao::clear_lattice_calc_hook(); });

  m.def("set_optimizer_hook", [](nb::callable fn) {
    Tao::set_optimizer_hook([fn](bool abort) -> bool {
      nb::gil_scoped_acquire acq;
      try {
        return obj_to_bool(fn(abort));
      } catch (nb::python_error &e) {
        report_error(e);
        return abort;
      }
    });
  });
  m.def("clear_optimizer_hook", []() { Tao::clear_optimizer_hook(); });

  m.def("set_merit_var_hook", [](nb::callable fn) {
    Tao::set_merit_var_hook([fn](int i_uni, int j_var, Bmad::TaoVarStruct &var) {
      nb::gil_scoped_acquire acq;
      try {
        fn(i_uni, j_var, ref(var));
      } catch (nb::python_error &e) {
        report_error(e);
      }
    });
  });
  m.def("clear_merit_var_hook", []() { Tao::clear_merit_var_hook(); });

  m.def("set_merit_data_hook", [](nb::callable fn) {
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
  });
  m.def("clear_merit_data_hook", []() { Tao::clear_merit_data_hook(); });
}
