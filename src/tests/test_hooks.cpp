#include <bmad.hpp>
#include <bmad/bmad_hooks.hpp>

#include <string>
#include <vector>

#include "doctest.h"

using namespace Bmad;

// Exercises the track1 pre/post-process hooks. Both must fire when an element is
// tracked, and an observer that leaves the out-parameters at their defaults must
// NOT abort tracking -- so track1_postprocess still runs afterwards. (err_flag is
// a pure output Bmad passes uninitialized; the trampolines default it to false so
// an observer does not accidentally signal an error. See src/hooks/.)
TEST_CASE("track1 preprocess/postprocess hooks") {
  auto lat = bmad_parser("data/csr_example/lat.bmad").lat;
  auto branch = lat.branch()[0];
  auto lat_param = lat.param().value();

  std::vector<std::string> pre, post;

  Bmad::set_track1_preprocess_hook(
      [&](CoordStruct &, EleStruct &ele, LatParamStruct &, bool &, bool &, bool &, TrackStruct *) {
        // Observe only: leaving err_flag/finished at their defaults lets Bmad track
        // the element normally (and lets track1_postprocess fire afterwards).
        pre.push_back(ele.name());
      }
  );
  Bmad::set_track1_postprocess_hook(
      [&](CoordStruct &, EleStruct &ele, LatParamStruct &, CoordStruct &) {
        post.push_back(ele.name());
      }
  );

  // On-axis starting orbit (survives the lattice).
  auto orbit = CoordStructAlloc1D();
  reallocate_coord(orbit, lat, 0);
  auto orb0{orbit[0]};
  auto ele0{branch.ele()[0]};
  init_coord(orb0, {}, ele0, Bmad::DOWNSTREAM_END);
  lattice_bookkeeper(lat);

  SUBCASE("single element via track1") {
    auto ele = branch.ele()[1];
    track1(orb0, ele, lat_param);

    CHECK(pre.size() == 1);
    CHECK(post.size() == 1); // postprocess fires only if preprocess did not abort
    CHECK(pre[0] == ele.name());
    CHECK(post[0] == ele.name());
  }

  SUBCASE("all elements via track_all") {
    auto result = track_all(lat, orbit);

    CHECK(result.track_state == Bmad::MOVING_FORWARD);
    CHECK(pre.size() == static_cast<size_t>(branch.n_ele_track()));
    CHECK(post.size() == pre.size());
  }

  Bmad::clear_track1_preprocess_hook();
  Bmad::clear_track1_postprocess_hook();
}
