// Parse a lattice with coherent-synchrotron-radiation (CSR) and space-charge
// tracking enabled, generate a beam, track its centroid and then the full beam
// relative to that centroid, and print / save the resulting distribution.
//
// Run from the repository root so the default (bundled) lattice path resolves, or
// pass a lattice file as the first argument:
//   ./debug/csr [lattice.bmad]

#include <bmad.hpp>

#include <fstream>
#include <iomanip>
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
  ran_seed_put(123456);

  // Beam parameters.
  auto beam_init{BeamInitStruct()};
  beam_init.set_a_norm_emit(4e-12);
  beam_init.set_b_norm_emit(4e-12);
  beam_init.set_dPz_dz(0);
  beam_init.set_sig_z(0.3e-3);
  beam_init.set_sig_pz(0e-20);
  beam_init.set_bunch_charge(0.01e-10);
  beam_init.set_n_particle(1000);
  beam_init.set_n_bunch(1);

  // Enable CSR + space charge.
  auto bmad_com{get_bmad_com()};
  bmad_com.set_csr_and_space_charge_on(true);

  auto space_charge_com{get_space_charge_com()};
  space_charge_com.set_ds_track_step(0.1);
  space_charge_com.set_n_bin(400);
  space_charge_com.set_beam_chamber_height(0.02);
  space_charge_com.set_n_shield_images(0);
  space_charge_com.set_particle_bin_span(8);

  auto ele0{lat.ele()[0]};
  auto lat_param{lat.param().value()};
  auto [beam, beam_err, beam_set] = init_beam_distribution(ele0, lat_param, beam_init);
  if (beam_err) {
    std::cerr << "init_beam_distribution failed\n";
    return 1;
  }

  // Compute the initial bunch centroid.
  auto particles = beam.bunch()[0].particle();
  size_t n_particles = particles.size();
  auto ave{FixedArray1D<double, 6>()};
  for (int i = 0; i < 6; ++i) {
    double sum = 0.0;
    for (const auto &p : particles)
      sum += p.vec()[i];
    ave[i] = (n_particles > 0) ? (sum / n_particles) : 0.0;
  }

  // Track the centroid through the lattice, then track the full beam relative to it.
  auto centroid{CoordStructAlloc1D()};
  reallocate_coord(centroid, lat, 0);
  auto centroid0{centroid[0]};
  init_coord(centroid0, ave, ele0, Bmad::DOWNSTREAM_END);
  track_all(lat, centroid);

  auto [beam1, beam1_err, beam1_set] = init_beam_distribution(ele0, lat_param, beam_init);
  track_beam(lat, beam1, std::nullopt, std::nullopt, centroid);

  auto first = beam1.bunch()[0].particle()[0].vec().to_vector();
  std::cout << "First particle after tracking: " << first << "\n";

  // Save the full distribution.
  std::ofstream out("csr.dat");
  out.setf(std::ios::fixed, std::ios::floatfield);
  out.precision(8);
  int idx = 1;
  for (const auto &part : beam1.bunch()[0].particle()) {
    out << std::setw(8) << idx << " " << part.vec().to_vector() << "\n";
    idx++;
  }
  std::cout << "Wrote csr.dat (" << (idx - 1) << " particles)\n";
  return 0;
}
