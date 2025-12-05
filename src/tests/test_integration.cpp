#include <fstream>
#include <iostream>
#include "bmad.hpp"
#include "bmad/generated/proxy.hpp"
// #include "bmad/fortran_arrays.hpp"
// #include "bmad/generated/routines.hpp"

using namespace Bmad;
using namespace SimUtils;

extern "C" int tao_c_init_tao(const char*);
extern "C" int tao_is_initialized();
extern "C" char* tao_c_out_io_buffer_get_line(int n);
extern "C" int tao_c_out_io_buffer_num_lines();

void test_allocatable() {
  auto ctr{EleProxyAllocatable1D(0, 10)};

  ctr.resize(0, 5);
  ctr[0].set_name("foo0");
  ctr[1].set_name("foo1");
  ctr[2].set_name("foo2");
  ctr[3].set_name("foo3");
  std::cout << "lbound is: " << ctr.view().lower_bound() << "\n";

  for (auto ele : ctr) {
    std::cout << ele.get_fortran_ptr() << " " << ele.name() << "\n";
  }
}

void test_real_container() {
  auto a{RealAllocatable1D()};
  auto b{RealAllocatable1D()};

  a.resize(0, 20);
  b.resize(0, 20);
  int n{10}, isn{1}, ierr{0};

  auto aview{a.view()}, bview{b.view()};
  aview[0] = 1.;
  aview[1] = 2.;
  aview[2] = 3.;
  bview[0] = 1.;
  bview[1] = 2.;
  bview[2] = 3.;

  fft1(a, b, n, isn, ierr);
  std::cout << "a=" << a.view().to_vector() << "\n";
  std::cout << "b=" << b.view().to_vector() << "\n";
  std::cout << "ierr=" << ierr << "\n";

  // TODO: ierr not coming back when isn!=1
}

void track_test() {
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

  auto centroid{CoordProxyAllocatable1D()};
  reallocate_coord_lat(centroid, lat, 0);
  auto centroid0{centroid[0]};
  init_coord1(centroid0, ave, ele0, Bmad::DOWNSTREAM_END);

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
}

void test_struct_access() {
  auto prox{AllEncompassingProxy()};
  auto cpx = std::complex<double>(1.0, 2.0);

  prox.set_complex_dp_0d(cpx);
  std::cout << "real: " << prox.complex_dp_0d().real() << "vs" << cpx.real()
            << "\n";
  std::cout << "imag: " << prox.complex_dp_0d().imag() << "vs" << cpx.imag()
            << "\n";
}

int main() {
  try {
    std::cout << "tao_is_initialized() = " << tao_is_initialized() << "\n";

    test_allocatable();
    test_struct_access();
    test_real_container();

    // auto initret = tao_c_init_tao(
    //     "-init_file $ACC_ROOT_DIR/bmad-doc/tao_examples/optics_matching/tao.init -noplot");
    // auto initret = tao_c_init_tao("-lat lattices/branch_fork.bmad -noplot");
    auto initret = tao_c_init_tao("-lat data/csr_example/lat.bmad -noplot");
    std::cout << "tao_c_init_tao returned: " << initret << "\n";
    std::cout << "tao_is_initialized() = " << tao_is_initialized() << "\n";

    std::cout << "Lines: " << tao_c_out_io_buffer_num_lines() << "\n";
    for (int i = 1; i < tao_c_out_io_buffer_num_lines(); i++) {
      std::cout << tao_c_out_io_buffer_get_line(i) << "\n";
    }
    std::cout << "Tao universes available: " << tao_get_n_universes() << "\n";

    Tao::TaoLatticeIndexProxy tao_lattice(1, Tao::LatticeType::MODEL);
    auto lattice_proxy = *tao_lattice;

    auto lat = lattice_proxy.lat();

    auto n_branches = lat.branch().size();
    std::cout << "Number of branches: " << n_branches << std::endl;

    for (const auto branch : lat.branch()) { // 1D_ALLOC_TYPE
      std::cout << "Number of elements in branch " << branch.ix_branch() << ": "
                << branch.ele().size() << std::endl;

      for (auto ele : branch.ele()) {
        std::cout << ele.ix_ele() << " @ " << ele.s() << ": " << ele.name()
                  << " (" << ele.type() << ")\n";
      }

      // Direct access still works
      // auto element_direct =
      //     *Bmad::TaoElementProxy(1, Bmad::LatticeType::MODEL, 0, 5);
    }
    track_test();

  } catch (const Bmad::BmadException& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
