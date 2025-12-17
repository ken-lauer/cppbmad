#include <bmad.hpp>
#include <iostream>
// #include "bmad/fortran_arrays.hpp"
// #include "bmad/generated/routines.hpp"

using namespace Bmad;
using namespace SimUtils;

extern "C" int tao_c_init_tao(const char*);
extern "C" int tao_is_initialized();
extern "C" char* tao_c_out_io_buffer_get_line(int n);
extern "C" int tao_c_out_io_buffer_num_lines();

void test_allocatable() {
  auto ctr{EleProxyAlloc1D(0, 10)};

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
  auto a{RealAlloc1D()};
  auto b{RealAlloc1D()};

  a.resize(0, 20);
  b.resize(0, 20);
  int n{10}, isn{1};

  auto aview{a.view()}, bview{b.view()};
  aview[0] = 1.;
  aview[1] = 2.;
  aview[2] = 3.;
  bview[0] = 1.;
  bview[1] = 2.;
  bview[2] = 3.;

  auto ierr = fft1(a, b, n, isn);
  std::cout << "a=" << a.view().to_vector() << "\n";
  std::cout << "b=" << b.view().to_vector() << "\n";
  std::cout << "ierr=" << ierr << "\n";

  // TODO: ierr not coming back when isn!=1
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

  } catch (const Bmad::BmadException& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
