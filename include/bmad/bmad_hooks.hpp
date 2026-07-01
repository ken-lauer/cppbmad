#pragma once

// Hand-written C++ interface to Bmad's tracking hook / custom procedure
// pointers (bmad_routine_interface).
//
// Each hook `<x>` can be given a C++ callable via `Bmad::set_<x>_hook(...)`;
// `Bmad::clear_<x>_hook()` restores Bmad's default behavior. The callable
// receives non-owning proxy views of the Fortran derived-type arguments, valid
// only for the duration of the call, and out-parameters by reference (write into
// them to affect Bmad). Optional derived-type / integer / logical arguments are
// passed as pointers that are null when the argument is absent; optional and
// required `coord_struct` arrays are passed as live `CoordStructArray1D` views.
//
// Callables must not throw across the boundary; exceptions are caught and
// swallowed at the trampoline (see src/hooks/bmad_hooks.cpp). This layer has no
// nanobind dependency. See plans/hook-plan.md.

#include <functional>

#include "bmad/generated/proxy.hpp"

namespace Bmad {

// (orbit, ele, param, stop_time [inout], init_needed [inout])
using TimeRungeKuttaHook =
    std::function<void(CoordStruct &, EleStruct &, LatParamStruct &, double &, int &)>;
void set_time_runge_kutta_periodic_kick_hook(TimeRungeKuttaHook fn);
void clear_time_runge_kutta_periodic_kick_hook();

// (bunch, ele, err [out], centroid [optional], direction [optional], finished [out], bunch_track
// [optional])
using Track1BunchHook = std::function<void(
    BunchStruct &,
    EleStruct &,
    bool &,
    CoordStructArray1D *,
    int *,
    bool &,
    BunchTrackStruct *
)>;
void set_track1_bunch_hook(Track1BunchHook fn);
void clear_track1_bunch_hook();

// (start_orb, ele, param, err_flag [out], finished [out], track [optional])
using Track1CustomHook = std::function<
    void(CoordStruct &, EleStruct &, LatParamStruct &, bool &, bool &, TrackStruct *)>;
void set_track1_custom_hook(Track1CustomHook fn);
void clear_track1_custom_hook();

// (finished [out], lat, orbit, ix_start, ix_end, direction, ix_branch [optional], track_state
// [optional])
using TrackManyHook =
    std::function<void(bool &, LatStruct &, CoordStructArray1D &, int, int, int, int *, int *)>;
void set_track_many_hook(TrackManyHook fn);
void clear_track_many_hook();

// (start_orb, ele, param, end_orb)
using Track1PostprocessHook =
    std::function<void(CoordStruct &, EleStruct &, LatParamStruct &, CoordStruct &)>;
void set_track1_postprocess_hook(Track1PostprocessHook fn);
void clear_track1_postprocess_hook();

// (start_orb, ele, param, err_flag [out], finished [out], radiation_included [out], track
// [optional])
using Track1PreprocessHook = std::function<
    void(CoordStruct &, EleStruct &, LatParamStruct &, bool &, bool &, bool &, TrackStruct *)>;
void set_track1_preprocess_hook(Track1PreprocessHook fn);
void clear_track1_preprocess_hook();

// (start_orb, ele, param, end_orb, err_flag [out], make_quaternion [optional, inout])
using Track1SpinCustomHook = std::function<
    void(CoordStruct &, EleStruct &, LatParamStruct &, CoordStruct &, bool &, bool *)>;
void set_track1_spin_custom_hook(Track1SpinCustomHook fn);
void clear_track1_spin_custom_hook();

// (bunch, ele, finished [out])
using Track1WakeHook = std::function<void(BunchStruct &, EleStruct &, bool &)>;
void set_track1_wake_hook(Track1WakeHook fn);
void clear_track1_wake_hook();

// (orb, ele, s)
using WallHitHook = std::function<void(CoordStruct &, EleStruct &, double)>;
void set_wall_hit_handler_custom_hook(WallHitHook fn);
void clear_wall_hit_handler_custom_hook();

} // namespace Bmad
