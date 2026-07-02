// Hand-written C++ glue between Bmad's Fortran hook trampolines
// (src/hooks/bmad_hooks.f90) and C++ callables declared in
// include/bmad/bmad_hooks.hpp.
//
// For each hook we keep a file-local std::function slot and an extern "C"
// trampoline matching the Fortran-side C ABI. The trampoline wraps each void*
// argument in a non-owning proxy (and each coord array in a live
// CoordStructArray1D view), invokes the std::function inside a try/catch(...) so
// no exception unwinds into Fortran, and writes back any out-parameters.

#include "bmad/bmad_hooks.hpp"

#include <cstddef>
#include <optional>
#include <utility>

#include "bmad/hook_util.hpp"

using Bmad::make_fortran_proxy;
using Bmad::hooks::log_hook_error;
using Bmad::hooks::make_coord_view;

// --- Fortran registration entry points -------------------------------------
extern "C" {
void bmad_hook_register_time_runge_kutta_periodic_kick(
    void (*)(void *, void *, void *, double *, int *)
);
void bmad_hook_register_track1_bunch(
    void (*)(void *, void *, int *, void *, int, int, std::size_t, void *, int *, void *)
);
void bmad_hook_register_track1_custom(void (*)(void *, void *, void *, int *, int *, void *));
void bmad_hook_register_track_many(
    void (*)(int *, void *, void *, int, int, std::size_t, int, int, int, void *, void *)
);
void bmad_hook_register_track1_postprocess(void (*)(void *, void *, void *, void *));
void bmad_hook_register_track1_preprocess(
    void (*)(void *, void *, void *, int *, int *, int *, void *)
);
void bmad_hook_register_track1_spin_custom(void (*)(void *, void *, void *, void *, int *, void *));
void bmad_hook_register_track1_wake(void (*)(void *, void *, int *));
void bmad_hook_register_wall_hit_handler_custom(void (*)(void *, void *, double));
}

namespace {

Bmad::TimeRungeKuttaHook g_time_rk;
extern "C" void
tramp_time_rk(void *orbit, void *ele, void *param, double *stop_time, int *init_needed) {
  if (!g_time_rk)
    return;
  try {
    auto po = make_fortran_proxy<Bmad::CoordStruct>(orbit);
    auto pe = make_fortran_proxy<Bmad::EleStruct>(ele);
    auto pp = make_fortran_proxy<Bmad::LatParamStruct>(param);
    g_time_rk(po, pe, pp, *stop_time, *init_needed);
  } catch (...) {
    log_hook_error("time_runge_kutta_periodic_kick");
  }
}

Bmad::Track1BunchHook g_track1_bunch;
extern "C" void tramp_track1_bunch(
    void *bunch,
    void *ele,
    int *err,
    void *cdata,
    int clb,
    int cub,
    std::size_t cesize,
    void *dir,
    int *finished,
    void *btrack
) {
  if (!g_track1_bunch)
    return;
  try {
    auto pb = make_fortran_proxy<Bmad::BunchStruct>(bunch);
    auto pe = make_fortran_proxy<Bmad::EleStruct>(ele);
    Bmad::CoordStructArray1D centroid = make_coord_view(cdata, clb, cub, cesize);
    Bmad::CoordStructArray1D *centroid_p = cdata ? &centroid : nullptr;
    std::optional<Bmad::BunchTrackStruct> bt_storage;
    Bmad::BunchTrackStruct *bt_p = nullptr;
    if (btrack) {
      bt_storage.emplace(make_fortran_proxy<Bmad::BunchTrackStruct>(btrack));
      bt_p = &*bt_storage;
    }
    bool e = (*err != 0);
    bool f = (*finished != 0);
    g_track1_bunch(pb, pe, e, centroid_p, static_cast<int *>(dir), f, bt_p);
    *err = e ? 1 : 0;
    *finished = f ? 1 : 0;
  } catch (...) {
    log_hook_error("track1_bunch");
    *err = 1;
  }
}

Bmad::Track1CustomHook g_track1_custom;
extern "C" void tramp_track1_custom(
    void *start_orb,
    void *ele,
    void *param,
    int *err_flag,
    int *finished,
    void *track
) {
  if (!g_track1_custom)
    return;
  try {
    auto ps = make_fortran_proxy<Bmad::CoordStruct>(start_orb);
    auto pe = make_fortran_proxy<Bmad::EleStruct>(ele);
    auto pp = make_fortran_proxy<Bmad::LatParamStruct>(param);
    std::optional<Bmad::TrackStruct> tk_storage;
    Bmad::TrackStruct *tk_p = nullptr;
    if (track) {
      tk_storage.emplace(make_fortran_proxy<Bmad::TrackStruct>(track));
      tk_p = &*tk_storage;
    }
    bool e = (*err_flag != 0);
    bool f = (*finished != 0);
    g_track1_custom(ps, pe, pp, e, f, tk_p);
    *err_flag = e ? 1 : 0;
    *finished = f ? 1 : 0;
  } catch (...) {
    log_hook_error("track1_custom");
    *err_flag = 1;
  }
}

Bmad::TrackManyHook g_track_many;
extern "C" void tramp_track_many(
    int *finished,
    void *lat,
    void *odata,
    int olb,
    int oub,
    std::size_t oesize,
    int ix_start,
    int ix_end,
    int direction,
    void *ix_branch,
    void *track_state
) {
  if (!g_track_many)
    return;
  try {
    auto pl = make_fortran_proxy<Bmad::LatStruct>(lat);
    Bmad::CoordStructArray1D orbit = make_coord_view(odata, olb, oub, oesize);
    bool f = (*finished != 0);
    g_track_many(
        f,
        pl,
        orbit,
        ix_start,
        ix_end,
        direction,
        static_cast<int *>(ix_branch),
        static_cast<int *>(track_state)
    );
    *finished = f ? 1 : 0;
  } catch (...) {
    log_hook_error("track_many");
  }
}

Bmad::Track1PostprocessHook g_track1_postprocess;
extern "C" void tramp_track1_postprocess(void *start_orb, void *ele, void *param, void *end_orb) {
  if (!g_track1_postprocess)
    return;
  try {
    auto ps = make_fortran_proxy<Bmad::CoordStruct>(start_orb);
    auto pe = make_fortran_proxy<Bmad::EleStruct>(ele);
    auto pp = make_fortran_proxy<Bmad::LatParamStruct>(param);
    auto pend = make_fortran_proxy<Bmad::CoordStruct>(end_orb);
    g_track1_postprocess(ps, pe, pp, pend);
  } catch (...) {
    log_hook_error("track1_postprocess");
  }
}

Bmad::Track1PreprocessHook g_track1_preprocess;
extern "C" void tramp_track1_preprocess(
    void *start_orb,
    void *ele,
    void *param,
    int *err_flag,
    int *finished,
    int *radiation_included,
    void *track
) {
  if (!g_track1_preprocess)
    return;
  try {
    auto ps = make_fortran_proxy<Bmad::CoordStruct>(start_orb);
    auto pe = make_fortran_proxy<Bmad::EleStruct>(ele);
    auto pp = make_fortran_proxy<Bmad::LatParamStruct>(param);
    std::optional<Bmad::TrackStruct> tk_storage;
    Bmad::TrackStruct *tk_p = nullptr;
    if (track) {
      tk_storage.emplace(make_fortran_proxy<Bmad::TrackStruct>(track));
      tk_p = &*tk_storage;
    }
    bool e = (*err_flag != 0);
    bool f = (*finished != 0);
    bool r = (*radiation_included != 0);
    g_track1_preprocess(ps, pe, pp, e, f, r, tk_p);
    *err_flag = e ? 1 : 0;
    *finished = f ? 1 : 0;
    *radiation_included = r ? 1 : 0;
  } catch (...) {
    log_hook_error("track1_preprocess");
  }
}

Bmad::Track1SpinCustomHook g_track1_spin_custom;
extern "C" void tramp_track1_spin_custom(
    void *start_orb,
    void *ele,
    void *param,
    void *end_orb,
    int *err_flag,
    void *make_quaternion
) {
  if (!g_track1_spin_custom)
    return;
  try {
    auto ps = make_fortran_proxy<Bmad::CoordStruct>(start_orb);
    auto pe = make_fortran_proxy<Bmad::EleStruct>(ele);
    auto pp = make_fortran_proxy<Bmad::LatParamStruct>(param);
    auto pend = make_fortran_proxy<Bmad::CoordStruct>(end_orb);
    int *mq_i = static_cast<int *>(make_quaternion);
    bool mq_storage = mq_i ? (*mq_i != 0) : false;
    bool *mq_p = mq_i ? &mq_storage : nullptr;
    bool e = (*err_flag != 0);
    g_track1_spin_custom(ps, pe, pp, pend, e, mq_p);
    *err_flag = e ? 1 : 0;
    if (mq_i)
      *mq_i = mq_storage ? 1 : 0;
  } catch (...) {
    log_hook_error("track1_spin_custom");
    *err_flag = 1;
  }
}

Bmad::Track1WakeHook g_track1_wake;
extern "C" void tramp_track1_wake(void *bunch, void *ele, int *finished) {
  if (!g_track1_wake)
    return;
  try {
    auto pb = make_fortran_proxy<Bmad::BunchStruct>(bunch);
    auto pe = make_fortran_proxy<Bmad::EleStruct>(ele);
    bool f = (*finished != 0);
    g_track1_wake(pb, pe, f);
    *finished = f ? 1 : 0;
  } catch (...) {
    log_hook_error("track1_wake");
  }
}

Bmad::WallHitHook g_wall_hit;
extern "C" void tramp_wall_hit(void *orb, void *ele, double s) {
  if (!g_wall_hit)
    return;
  try {
    auto po = make_fortran_proxy<Bmad::CoordStruct>(orb);
    auto pe = make_fortran_proxy<Bmad::EleStruct>(ele);
    g_wall_hit(po, pe, s);
  } catch (...) {
    log_hook_error("wall_hit_handler_custom");
  }
}

} // namespace

// --- set_*/clear_* -----------------------------------------------------------

#define BMAD_HOOK_SETCLEAR(name, HookType, slot, regfn, tramp) \
  void Bmad::set_##name##_hook(HookType fn) {                  \
    slot = std::move(fn);                                      \
    regfn(slot ? &tramp : nullptr);                            \
  }                                                            \
  void Bmad::clear_##name##_hook() {                           \
    slot = nullptr;                                            \
    regfn(nullptr);                                            \
  }

BMAD_HOOK_SETCLEAR(
    time_runge_kutta_periodic_kick,
    TimeRungeKuttaHook,
    g_time_rk,
    bmad_hook_register_time_runge_kutta_periodic_kick,
    tramp_time_rk
)
BMAD_HOOK_SETCLEAR(
    track1_bunch,
    Track1BunchHook,
    g_track1_bunch,
    bmad_hook_register_track1_bunch,
    tramp_track1_bunch
)
BMAD_HOOK_SETCLEAR(
    track1_custom,
    Track1CustomHook,
    g_track1_custom,
    bmad_hook_register_track1_custom,
    tramp_track1_custom
)
BMAD_HOOK_SETCLEAR(
    track_many,
    TrackManyHook,
    g_track_many,
    bmad_hook_register_track_many,
    tramp_track_many
)
BMAD_HOOK_SETCLEAR(
    track1_postprocess,
    Track1PostprocessHook,
    g_track1_postprocess,
    bmad_hook_register_track1_postprocess,
    tramp_track1_postprocess
)
BMAD_HOOK_SETCLEAR(
    track1_preprocess,
    Track1PreprocessHook,
    g_track1_preprocess,
    bmad_hook_register_track1_preprocess,
    tramp_track1_preprocess
)
BMAD_HOOK_SETCLEAR(
    track1_spin_custom,
    Track1SpinCustomHook,
    g_track1_spin_custom,
    bmad_hook_register_track1_spin_custom,
    tramp_track1_spin_custom
)
BMAD_HOOK_SETCLEAR(
    track1_wake,
    Track1WakeHook,
    g_track1_wake,
    bmad_hook_register_track1_wake,
    tramp_track1_wake
)
BMAD_HOOK_SETCLEAR(
    wall_hit_handler_custom,
    WallHitHook,
    g_wall_hit,
    bmad_hook_register_wall_hit_handler_custom,
    tramp_wall_hit
)

#undef BMAD_HOOK_SETCLEAR
