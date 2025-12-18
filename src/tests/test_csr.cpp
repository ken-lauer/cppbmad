#include <bmad.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace Bmad;
using namespace SimUtils;

int main() {
  auto lat = bmad_parser("data/csr_example/lat.bmad").lat;
  ran_seed_put(123456);

  auto beam_init{BeamInitProxy()};
  beam_init.set_a_norm_emit(4e-12);
  beam_init.set_b_norm_emit(4e-12);
  beam_init.set_dPz_dz(0);
  beam_init.set_sig_z(0.3e-3);
  beam_init.set_sig_pz(0e-20);
  beam_init.set_bunch_charge(0.01e-10);
  beam_init.set_n_particle(1000);
  beam_init.set_n_bunch(1);

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
  auto [beam, err_flag, beam_init_set] =
      init_beam_distribution(ele0, lat_param, beam_init);

  // First bunch and its particles
  auto bunch = beam.bunch()[0];
  auto particles = bunch.particle();
  size_t n_particles = particles.size();

  // Calculate the average (centroid)
  auto ave{FixedArray1D<double, 6>()};
  for (int i = 0; i < 6; ++i) {
    double sum = 0.0;
    for (const auto& p : particles) {
      sum += p.vec()[i];
    }
    ave[i] = (n_particles > 0) ? (sum / n_particles) : 0.0;
  }

  auto centroid{CoordProxyAlloc1D()};
  reallocate_coord(centroid, lat, 0);
  auto centroid0{centroid[0]};
  init_coord(centroid0, ave, ele0, Bmad::DOWNSTREAM_END);

  track_all(lat, centroid);

  auto [beam1, err_flag1, beam_init_set1] =
      init_beam_distribution(ele0, lat_param, beam_init);
  std::cout << "init_beam_distribution err=" << err_flag1 << "\n";

  auto track_res = track_beam(lat, beam1, std::nullopt, std::nullopt, centroid);

  std::cout << "track_beam result=" << track_res << "\n";

  std::cout << "First particle coords at end of lattice:\n";
  std::cout << beam1.bunch()[0].particle()[0].vec().to_vector() << "\n";

  std::ofstream out("csr.dat");
  out.setf(std::ios::fixed, std::ios::floatfield);
  out.precision(8);
  int idx = 1;
  for (const auto& part : beam1.bunch()[0].particle()) {
    out << std::setw(8) << idx << " " << part.vec().to_vector() << "\n";
    idx++;
  }
  out.close();
  return 0;
}
