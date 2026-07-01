// Hand-written nanobind bindings exposing Bmad's tracking hooks to Python as
// bmad.set_<x>_hook(callable) / bmad.clear_<x>_hook().
//
// Each set_* wraps the Python callable in a std::function that acquires the GIL,
// hands the Fortran-owned proxies (and coord-array views) to Python as
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

#include <nanobind/stl/pair.h>
#include <nanobind/stl/tuple.h>

#include "pybmad/hook_util.hpp"

using Pybmad::hooks::ref;
using Pybmad::hooks::report_error;

void Pybmad::init_bmad_hooks(nb::module_ &m) {
  m.def("set_time_runge_kutta_periodic_kick_hook", [](nb::callable fn) {
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
  });
  m.def("clear_time_runge_kutta_periodic_kick_hook", []() {
    Bmad::clear_time_runge_kutta_periodic_kick_hook();
  });

  m.def("set_track1_bunch_hook", [](nb::callable fn) {
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
  });
  m.def("clear_track1_bunch_hook", []() { Bmad::clear_track1_bunch_hook(); });

  m.def("set_track1_custom_hook", [](nb::callable fn) {
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
  });
  m.def("clear_track1_custom_hook", []() { Bmad::clear_track1_custom_hook(); });

  m.def("set_track_many_hook", [](nb::callable fn) {
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
        nb::object r = fn(finished, ref(lat), ref(orbit), ix_start, ix_end, direction, rib, rts);
        if (!r.is_none())
          finished = nb::cast<bool>(r);
      } catch (nb::python_error &e) {
        report_error(e);
      }
    });
  });
  m.def("clear_track_many_hook", []() { Bmad::clear_track_many_hook(); });

  m.def("set_track1_postprocess_hook", [](nb::callable fn) {
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
  });
  m.def("clear_track1_postprocess_hook", []() { Bmad::clear_track1_postprocess_hook(); });

  m.def("set_track1_preprocess_hook", [](nb::callable fn) {
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
            fn(ref(start_orb), ref(ele), ref(param), err_flag, finished, radiation_included, rt);
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
  });
  m.def("clear_track1_preprocess_hook", []() { Bmad::clear_track1_preprocess_hook(); });

  m.def("set_track1_spin_custom_hook", [](nb::callable fn) {
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
  });
  m.def("clear_track1_spin_custom_hook", []() { Bmad::clear_track1_spin_custom_hook(); });

  m.def("set_track1_wake_hook", [](nb::callable fn) {
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
  });
  m.def("clear_track1_wake_hook", []() { Bmad::clear_track1_wake_hook(); });

  m.def("set_wall_hit_handler_custom_hook", [](nb::callable fn) {
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
  });
  m.def("clear_wall_hit_handler_custom_hook", []() { Bmad::clear_wall_hit_handler_custom_hook(); });
}
