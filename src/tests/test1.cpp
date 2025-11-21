#include <iostream>
#include <string>
#include "bmad_std_typedef.h"
// #include "pybmad.hpp"
#include "tao_proxies.hpp"

using namespace Bmad;

extern "C" int tao_c_init_tao(const char*);
extern "C" int tao_is_initialized();
extern "C" char* tao_c_out_io_buffer_get_line(int n);
extern "C" int tao_c_out_io_buffer_num_lines();

int main() {
  try {
    std::cout << "tao_is_initialized() = " << tao_is_initialized() << "\n";

    // auto initret = tao_c_init_tao("-init_file
    // $ACC_ROOT_DIR/bmad-doc/tao_examples/optics_matching/tao.init");
    auto initret = tao_c_init_tao("-lat lattices/branch_fork.bmad");
    std::cout << "tao_c_init_tao returned: " << initret << "\n";
    std::cout << "tao_is_initialized() = " << tao_is_initialized() << "\n";

    std::cout << "Lines: " << tao_c_out_io_buffer_num_lines() << "\n";
    for (int i = 1; i < tao_c_out_io_buffer_num_lines(); i++) {
      std::cout << tao_c_out_io_buffer_get_line(i) << "\n";
    }
    std::cout << "Tao universes available: " << tao_get_n_universes() << "\n";

    tao::TaoLatticeIndexProxy tao_lattice(1, tao::LatticeType::MODEL);
    auto lattice_proxy = *tao_lattice;

    auto lat = lattice_proxy.lat();

    auto n_branches = lat.branch().size();
    std::cout << "Number of branches: " << n_branches << std::endl;

    for (const auto branch : lat.branch()) { // 1D_ALLOC_TYPE
      // int n_elements = branch.get_n_elements();
      std::cout << "Number of elements in branch " << branch.ix_branch() << ": "
                << branch.ele().size() << std::endl;

      for (auto ele : branch.ele()) {
        std::cout << ele.ix_ele() << " @ " << ele.s() << ": " << ele.name()
                  << " (" << ele.type() << ")\n";
        // std::cout << "absolute_time_tracking(ele) = "
        //           << absolute_time_tracking(ele) << "\n";
        // auto cpp_ele = ele.deepcopy();
        // std::cout << "-- deepcopy result:" << cpp_ele->name << "\n";
        // std::cout << element.spin_taylor_ref_orb_in().to_vector() << "\n";
        // std::cout << cpp_ele->spin_taylor_ref_orb_in << "\n";
        //
        // json js;
        // to_json(js, *cpp_ele);
        // std::cout << js;
        //
        // for (const auto taylor : ele.taylor()) { // 1D_NOT_TYPE
        //   std::cout << "taylor ref " << taylor.ref()
        //             << ", n_terms=" << taylor.term().size() << " \n";
        //   for (const auto term : taylor.term()) { // 1D_PTR_TYPE
        //     std::cout << "    term: coef=" << term.coef() << " = "
        //               << term.expn().to_vector() << "\n";
        //   }
        // }
        // std::cout << " attrs           : " << ele.value().to_vector() <<
        // "\n";
        //  TODO the deepcopy is off-by-1 with an extra 0 at the start, I think
        // std::cout << " attrs (deepcopy):  " << cpp_ele->value << "\n";
      }

      // Direct access still works
      // auto element_direct =
      //     *tao::TaoElementProxy(1, tao::LatticeType::MODEL, 0, 5);
    }

  } catch (const tao::TaoException& e) {
    std::cerr << "Tao error: " << e.what() << std::endl;
    return 1;
  }

  // dspline_len(0.0, double s_chord1, CPP_spline &spline, std::optional<double>
  // dtheta_ref)
  return 0;
}
