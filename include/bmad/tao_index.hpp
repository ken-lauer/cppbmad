#pragma once

#include "bmad/fortran_arrays.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/proxy_base.hpp"

#include <complex>
#include <memory>
#include <string>

extern "C" {

int tao_get_n_universes();
void* tao_c_get_universe_ptr(int ix_uni);
void* tao_c_get_tao_lattice_ptr(int ix_uni, int ix_lat);
void* tao_c_get_lattice_ptr(int ix_uni, int ix_lat);
void* tao_c_get_branch_ptr(int ix_uni, int ix_lat, int ix_branch);
void* tao_c_get_element_ptr(int ix_uni, int ix_lat, int ix_branch, int ix_ele);

int tao_lat_get_n_branches(void* lat_ptr);
void* tao_lat_get_branch_ptr(void* lat_ptr, int ix_branch);
int tao_branch_get_n_elements(void* branch_ptr);
void* tao_branch_get_element_ptr(void* branch_ptr, int ix_ele);
}

namespace Tao {

enum class LatticeType : int { MODEL = 1, DESIGN = 2, BASE = 3 };

class TaoUniverseIndexProxy;
class TaoLatticeIndexProxy;
class TaoBranchIndexProxy;

class TaoElementIndexProxy {
 private:
  int ix_uni_, ix_lat_, ix_branch_, ix_ele_;

  void* get_fortran_ptr_() const {
    void* ele_ptr =
        tao_c_get_element_ptr(ix_uni_, ix_lat_, ix_branch_, ix_ele_);
    if (!ele_ptr) {
      throw Bmad::NullPointerException(
          "TaoElementIndexProxy dereference for ix_uni=" +
          std::to_string(ix_uni_) + " ix_lat=" + std::to_string(ix_lat_) +
          " ix_branch=" + std::to_string(ix_branch_) +
          " ix_ele=" + std::to_string(ix_ele_) + "");
    }
    return ele_ptr;
  }

 public:
  TaoElementIndexProxy(
      int ix_uni,
      LatticeType lattice_type,
      int ix_branch,
      int ix_ele)
      : ix_uni_(ix_uni),
        ix_lat_(static_cast<int>(lattice_type)),
        ix_branch_(ix_branch),
        ix_ele_(ix_ele) {}

  Bmad::EleProxy operator*() const {
    return Bmad::EleProxy(get_fortran_ptr_());
  }

  std::unique_ptr<Bmad::EleProxy> operator->() const {
    return std::make_unique<Bmad::EleProxy>(get_fortran_ptr_());
  }
};

class TaoBranchIndexProxy {
 private:
  int ix_uni_, ix_lat_, ix_branch_;

  void* get_fortran_ptr_() const {
    void* branch_ptr = tao_c_get_branch_ptr(ix_uni_, ix_lat_, ix_branch_);
    if (!branch_ptr) {
      throw Bmad::NullPointerException(
          "TaoBranchIndexProxy dereference for [" + std::to_string(ix_uni_) +
          "," + std::to_string(ix_lat_) + "," + std::to_string(ix_branch_) +
          "]");
    }
    return branch_ptr;
  }

 public:
  TaoBranchIndexProxy(int ix_uni, LatticeType lattice_type, int ix_branch)
      : ix_uni_(ix_uni),
        ix_lat_(static_cast<int>(lattice_type)),
        ix_branch_(ix_branch) {}

  Bmad::BranchProxy operator*() const {
    return Bmad::BranchProxy(get_fortran_ptr_());
  }

  std::unique_ptr<Bmad::BranchProxy> operator->() const {
    return std::make_unique<Bmad::BranchProxy>(get_fortran_ptr_());
  }

  TaoElementIndexProxy get_element(int ix_ele) const {
    return TaoElementIndexProxy(
        ix_uni_, static_cast<LatticeType>(ix_lat_), ix_branch_, ix_ele);
  }
};

class TaoLatticeIndexProxy {
 private:
  int ix_uni_, ix_lat_;

  void* get_fortran_ptr_() const {
    void* lat_ptr = tao_c_get_tao_lattice_ptr(ix_uni_, ix_lat_);
    if (!lat_ptr) {
      throw Bmad::NullPointerException(
          "TaoLatticeIndexProxy dereference for [" + std::to_string(ix_uni_) +
          "," + std::to_string(ix_lat_) + "]");
    }
    return lat_ptr;
  }

 public:
  TaoLatticeIndexProxy(int ix_uni, LatticeType lattice_type)
      : ix_uni_(ix_uni), ix_lat_(static_cast<int>(lattice_type)) {}

  TaoBranchIndexProxy get_branch(int ix_branch) const {
    return TaoBranchIndexProxy(
        ix_uni_, static_cast<LatticeType>(ix_lat_), ix_branch);
  }

  Bmad::TaoLatticeProxy operator*() const {
    return Bmad::TaoLatticeProxy(get_fortran_ptr_());
  }
  std::unique_ptr<Bmad::TaoLatticeProxy> operator->() const {
    return std::make_unique<Bmad::TaoLatticeProxy>(**this);
  }
};

class TaoUniverseIndexProxy {
 private:
  int ix_uni_;

  void* get_fortran_ptr_() const {
    void* uni_ptr = tao_c_get_universe_ptr(ix_uni_);
    if (!uni_ptr) {
      throw Bmad::NullPointerException(
          "TaoUniverseIndexProxy dereference for universe " +
          std::to_string(ix_uni_));
    }
    return uni_ptr;
  }

 public:
  explicit TaoUniverseIndexProxy(int ix_uni) : ix_uni_(ix_uni) {}

  Bmad::TaoUniverseProxy operator*() const;
  std::unique_ptr<Bmad::TaoUniverseProxy> operator->() const {
    return std::make_unique<Bmad::TaoUniverseProxy>(get_fortran_ptr_());
  }

  TaoLatticeIndexProxy get_lattice(LatticeType lattice_type) const {
    return TaoLatticeIndexProxy(ix_uni_, lattice_type);
  }
};

} // namespace Tao
