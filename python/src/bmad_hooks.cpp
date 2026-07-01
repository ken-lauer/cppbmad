// Hand-written nanobind bindings exposing Bmad's tracking hooks to Python as
// read/write properties on a singleton holder: `bmad.hooks.<name> = callable`,
// `bmad.hooks.<name> = None` to clear, and reading `bmad.hooks.<name>` returns
// the current callable or None.
//
// Each property setter wraps the Python callable in a std::function that acquires
// the GIL, hands the Fortran-owned proxies (and coord-array views) to Python as
// non-owning references, calls the Python callable, and translates its return
// value into the hook's out-parameters. Out-parameters are passed to Python by
// value and written back from the callable's return value:
//   * a callable that returns None leaves all out-parameters unchanged,
//   * otherwise it returns the out-parameters in signature order (a tuple when
//     there is more than one).
// Optional arguments (absent in Fortran) are passed to Python as None.
// Python exceptions are caught and reported so they never unwind into Fortran.
// The proxies/views are borrowed and valid only for the duration of the call.

#include "pybmad/bmad_hooks.hpp"

#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/tuple.h>

#include "pybmad/hook_util.hpp"

using Pybmad::hooks::ref;
using Pybmad::hooks::report_error;

namespace {

// Holds the current Python callable for each hook (for readback via the getter).
struct BmadHooks {
  nb::object time_runge_kutta_periodic_kick = nb::none();
  nb::object track1_bunch = nb::none();
  nb::object track1_custom = nb::none();
  nb::object track_many = nb::none();
  nb::object track1_postprocess = nb::none();
  nb::object track1_preprocess = nb::none();
  nb::object track1_spin_custom = nb::none();
  nb::object track1_wake = nb::none();
  nb::object wall_hit_handler_custom = nb::none();
};

} // namespace

void Pybmad::init_bmad_hooks(nb::module_ &m) {
  nb::class_<BmadHooks> cls(m, "BmadHooks", "Registry of Bmad tracking-hook callbacks.");

  cls.def_prop_rw(
      "time_runge_kutta_periodic_kick",
      [](BmadHooks &h) { return h.time_runge_kutta_periodic_kick; },
      [](BmadHooks &h, std::optional<nb::object> value) {
        nb::object fn = value ? *value : nb::none();
        h.time_runge_kutta_periodic_kick = fn;
        if (fn.is_none()) {
          Bmad::clear_time_runge_kutta_periodic_kick_hook();
          return;
        }
        Bmad::set_time_runge_kutta_periodic_kick_hook([fn](
                                                          Bmad::CoordStruct &orbit,
                                                          Bmad::EleStruct &ele,
                                                          Bmad::LatParamStruct &param,
                                                          double &stop_time,
                                                          int &init_needed
                                                      ) {
          nb::gil_scoped_acquire acq;
          try {
            nb::object r = fn(ref(orbit), ref(ele), ref(param), stop_time, init_needed);
            if (!r.is_none()) {
              auto t = nb::cast<std::pair<double, int>>(r);
              stop_time = t.first;
              init_needed = t.second;
            }
          } catch (nb::python_error &e) {
            report_error(e);
          }
        });
      }
  );

  cls.def_prop_rw(
      "track1_bunch",
      [](BmadHooks &h) { return h.track1_bunch; },
      [](BmadHooks &h, std::optional<nb::object> value) {
        nb::object fn = value ? *value : nb::none();
        h.track1_bunch = fn;
        if (fn.is_none()) {
          Bmad::clear_track1_bunch_hook();
          return;
        }
        Bmad::set_track1_bunch_hook([fn](
                                        Bmad::BunchStruct &bunch,
                                        Bmad::EleStruct &ele,
                                        bool &err,
                                        Bmad::CoordStructArray1D *centroid,
                                        int *direction,
                                        bool &finished,
                                        Bmad::BunchTrackStruct *bunch_track
                                    ) {
          nb::gil_scoped_acquire acq;
          try {
            nb::object rc = centroid ? ref(*centroid) : nb::none();
            nb::object rd = direction ? nb::cast(*direction) : nb::none();
            nb::object rb = bunch_track ? ref(*bunch_track) : nb::none();
            nb::object r = fn(ref(bunch), ref(ele), err, rc, rd, finished, rb);
            if (!r.is_none()) {
              auto t = nb::cast<std::pair<bool, bool>>(r);
              err = t.first;
              finished = t.second;
            }
          } catch (nb::python_error &e) {
            report_error(e);
          }
        });
      }
  );

  cls.def_prop_rw(
      "track1_custom",
      [](BmadHooks &h) { return h.track1_custom; },
      [](BmadHooks &h, std::optional<nb::object> value) {
        nb::object fn = value ? *value : nb::none();
        h.track1_custom = fn;
        if (fn.is_none()) {
          Bmad::clear_track1_custom_hook();
          return;
        }
        Bmad::set_track1_custom_hook([fn](
                                         Bmad::CoordStruct &start_orb,
                                         Bmad::EleStruct &ele,
                                         Bmad::LatParamStruct &param,
                                         bool &err_flag,
                                         bool &finished,
                                         Bmad::TrackStruct *track
                                     ) {
          nb::gil_scoped_acquire acq;
          try {
            nb::object rt = track ? ref(*track) : nb::none();
            nb::object r = fn(ref(start_orb), ref(ele), ref(param), err_flag, finished, rt);
            if (!r.is_none()) {
              auto t = nb::cast<std::pair<bool, bool>>(r);
              err_flag = t.first;
              finished = t.second;
            }
          } catch (nb::python_error &e) {
            report_error(e);
          }
        });
      }
  );

  cls.def_prop_rw(
      "track_many",
      [](BmadHooks &h) { return h.track_many; },
      [](BmadHooks &h, std::optional<nb::object> value) {
        nb::object fn = value ? *value : nb::none();
        h.track_many = fn;
        if (fn.is_none()) {
          Bmad::clear_track_many_hook();
          return;
        }
        Bmad::set_track_many_hook([fn](
                                      bool &finished,
                                      Bmad::LatStruct &lat,
                                      Bmad::CoordStructArray1D &orbit,
                                      int ix_start,
                                      int ix_end,
                                      int direction,
                                      int *ix_branch,
                                      int *track_state
                                  ) {
          nb::gil_scoped_acquire acq;
          try {
            nb::object rib = ix_branch ? nb::cast(*ix_branch) : nb::none();
            nb::object rts = track_state ? nb::cast(*track_state) : nb::none();
            nb::object r =
                fn(finished, ref(lat), ref(orbit), ix_start, ix_end, direction, rib, rts);
            if (!r.is_none())
              finished = nb::cast<bool>(r);
          } catch (nb::python_error &e) {
            report_error(e);
          }
        });
      }
  );

  cls.def_prop_rw(
      "track1_postprocess",
      [](BmadHooks &h) { return h.track1_postprocess; },
      [](BmadHooks &h, std::optional<nb::object> value) {
        nb::object fn = value ? *value : nb::none();
        h.track1_postprocess = fn;
        if (fn.is_none()) {
          Bmad::clear_track1_postprocess_hook();
          return;
        }
        Bmad::set_track1_postprocess_hook([fn](
                                              Bmad::CoordStruct &start_orb,
                                              Bmad::EleStruct &ele,
                                              Bmad::LatParamStruct &param,
                                              Bmad::CoordStruct &end_orb
                                          ) {
          nb::gil_scoped_acquire acq;
          try {
            fn(ref(start_orb), ref(ele), ref(param), ref(end_orb));
          } catch (nb::python_error &e) {
            report_error(e);
          }
        });
      }
  );

  cls.def_prop_rw(
      "track1_preprocess",
      [](BmadHooks &h) { return h.track1_preprocess; },
      [](BmadHooks &h, std::optional<nb::object> value) {
        nb::object fn = value ? *value : nb::none();
        h.track1_preprocess = fn;
        if (fn.is_none()) {
          Bmad::clear_track1_preprocess_hook();
          return;
        }
        Bmad::set_track1_preprocess_hook([fn](
                                             Bmad::CoordStruct &start_orb,
                                             Bmad::EleStruct &ele,
                                             Bmad::LatParamStruct &param,
                                             bool &err_flag,
                                             bool &finished,
                                             bool &radiation_included,
                                             Bmad::TrackStruct *track
                                         ) {
          nb::gil_scoped_acquire acq;
          try {
            nb::object rt = track ? ref(*track) : nb::none();
            nb::object r =
                fn(ref(start_orb),
                   ref(ele),
                   ref(param),
                   err_flag,
                   finished,
                   radiation_included,
                   rt);
            if (!r.is_none()) {
              auto t = nb::cast<std::tuple<bool, bool, bool>>(r);
              err_flag = std::get<0>(t);
              finished = std::get<1>(t);
              radiation_included = std::get<2>(t);
            }
          } catch (nb::python_error &e) {
            report_error(e);
          }
        });
      }
  );

  cls.def_prop_rw(
      "track1_spin_custom",
      [](BmadHooks &h) { return h.track1_spin_custom; },
      [](BmadHooks &h, std::optional<nb::object> value) {
        nb::object fn = value ? *value : nb::none();
        h.track1_spin_custom = fn;
        if (fn.is_none()) {
          Bmad::clear_track1_spin_custom_hook();
          return;
        }
        Bmad::set_track1_spin_custom_hook([fn](
                                              Bmad::CoordStruct &start_orb,
                                              Bmad::EleStruct &ele,
                                              Bmad::LatParamStruct &param,
                                              Bmad::CoordStruct &end_orb,
                                              bool &err_flag,
                                              bool *make_quaternion
                                          ) {
          nb::gil_scoped_acquire acq;
          try {
            nb::object rmq = make_quaternion ? nb::cast(*make_quaternion) : nb::none();
            nb::object r = fn(ref(start_orb), ref(ele), ref(param), ref(end_orb), err_flag, rmq);
            if (!r.is_none()) {
              // Either err_flag alone, or (err_flag, make_quaternion).
              if (nb::isinstance<nb::tuple>(r)) {
                auto t = nb::cast<std::pair<bool, bool>>(r);
                err_flag = t.first;
                if (make_quaternion)
                  *make_quaternion = t.second;
              } else {
                err_flag = nb::cast<bool>(r);
              }
            }
          } catch (nb::python_error &e) {
            report_error(e);
          }
        });
      }
  );

  cls.def_prop_rw(
      "track1_wake",
      [](BmadHooks &h) { return h.track1_wake; },
      [](BmadHooks &h, std::optional<nb::object> value) {
        nb::object fn = value ? *value : nb::none();
        h.track1_wake = fn;
        if (fn.is_none()) {
          Bmad::clear_track1_wake_hook();
          return;
        }
        Bmad::set_track1_wake_hook(
            [fn](Bmad::BunchStruct &bunch, Bmad::EleStruct &ele, bool &finished) {
              nb::gil_scoped_acquire acq;
              try {
                nb::object r = fn(ref(bunch), ref(ele), finished);
                if (!r.is_none())
                  finished = nb::cast<bool>(r);
              } catch (nb::python_error &e) {
                report_error(e);
              }
            }
        );
      }
  );

  cls.def_prop_rw(
      "wall_hit_handler_custom",
      [](BmadHooks &h) { return h.wall_hit_handler_custom; },
      [](BmadHooks &h, std::optional<nb::object> value) {
        nb::object fn = value ? *value : nb::none();
        h.wall_hit_handler_custom = fn;
        if (fn.is_none()) {
          Bmad::clear_wall_hit_handler_custom_hook();
          return;
        }
        Bmad::set_wall_hit_handler_custom_hook(
            [fn](Bmad::CoordStruct &orb, Bmad::EleStruct &ele, double s) {
              nb::gil_scoped_acquire acq;
              try {
                fn(ref(orb), ref(ele), s);
              } catch (nb::python_error &e) {
                report_error(e);
              }
            }
        );
      }
  );

  m.attr("hooks") = nb::cast(BmadHooks{});
}
