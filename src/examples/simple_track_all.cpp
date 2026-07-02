// Parse a lattice, initialize a single particle at the start of the root branch,
// track it through every element, and print the orbit at each element.
//
// Run from the repository root so the default (bundled) lattice path resolves, or
// pass a lattice file as the first argument:
//   ./debug/simple_track_all [lattice.bmad]

#include <bmad.hpp>

#include <iostream>
#include <string>

using namespace Bmad;
using namespace SimUtils; // for operator<< on std::vector

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

  // Build and initialize the starting orbit for the root branch.
  auto orbit = CoordStructAlloc1D();
  reallocate_coord(orbit, lat, 0);
  auto orb0{orbit[0]};
  auto ele0{branch.ele()[0]};
  init_coord(orb0, {1e-4, 0, 0, 0, 0, 0}, ele0, Bmad::DOWNSTREAM_END);

  auto result = track_all(lat, orbit);

  int ix_end = branch.n_ele_track() - 1;
  if (result.track_state != Bmad::MOVING_FORWARD) {
    std::cout << "Particle lost at element: " << result.track_state << "\n";
    ix_end = result.track_state;
  }

  for (int ie = 0; ie <= ix_end; ++ie) {
    std::cout << branch.ele()[ie].name() << " " << orbit[ie].vec().to_vector() << "\n";
  }
  return 0;
}
