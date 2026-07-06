// Same as simple_track_all, but installs track1 pre/post-process hooks that
// observe each element as it is tracked. The hook arguments are live, non-owning
// proxy views valid only for the duration of the call.
//
// Run from the repository root so the default (bundled) lattice path resolves, or
// pass a lattice file as the first argument:
//   ./debug/hook_track_all [lattice.bmad]

#include <bmad.hpp>
#include <bmad/bmad_hooks.hpp>

#include <iostream>
#include <string>
#include <vector>

using namespace Bmad;

int main(int argc, char *argv[]) {
  std::string lat_file = (argc > 1) ? argv[1] : "data/csr_example/lat.bmad";

  std::cout << "Parsing lattice: " << lat_file << "\n";
  auto parsed = bmad_parser(lat_file);
  if (parsed.err_flag) {
    std::cerr << "Failed to parse lattice: " << lat_file << "\n";
    return 1;
  }
  auto lat = parsed.lat;
  auto branch = lat.branch()[0];

  // Observe each element before it is tracked. Leave err_flag/finished at their
  // defaults so Bmad tracks the element normally; setting finished=true would make
  // the hook replace the tracking (and skip track1_postprocess).
  std::vector<std::string> preprocessed;
  Bmad::set_track1_preprocess_hook(
      [&](CoordStruct &, EleStruct &ele, LatParamStruct &, bool &, bool &, bool &, TrackStruct *) {
        preprocessed.push_back(ele.name());
      }
  );

  // Observe the exit orbit after each element (end_orb is a live proxy).
  Bmad::set_track1_postprocess_hook(
      [](CoordStruct &, EleStruct &ele, LatParamStruct &, CoordStruct &end_orb) {
        std::cout << "postprocess " << ele.name() << ": x = " << end_orb.vec()[0] << "\n";
      }
  );

  auto orbit = CoordStructAlloc1D();
  reallocate_coord(orbit, lat, 0);
  auto orb0{orbit[0]};
  auto ele0{branch.ele()[0]};
  init_coord(orb0, FixedArray1D<double, 6>{1e-4, 0, 0, 0, 0, 0}, ele0, Bmad::DOWNSTREAM_END);

  auto result = track_all(lat, orbit);

  // Clear by using the following:
  // Bmad::clear_track1_preprocess_hook();
  // Bmad::clear_track1_postprocess_hook();

  std::cout << "\nPreprocess hook fired for " << preprocessed.size() << " element(s).\n";
  if (result.track_state != Bmad::MOVING_FORWARD)
    std::cout << "Particle lost at element: " << result.track_state << "\n";
  return 0;
}
