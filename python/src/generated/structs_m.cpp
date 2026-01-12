#include "pybmad/generated/structs_m.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// mode3_struct
void init_mode3_struct(py::module& m, py::class_<Mode3Proxy>& cls) {
  cls.def(py::init<>())
      // Mode3Proxy.v (2D_NOT_real -
      .def_property_readonly("v", &Mode3Proxy::v)
      // Mode3Proxy.a (0D_NOT_type -
      .def_property("a", &Mode3Proxy::a, &Mode3Proxy::set_a)
      // Mode3Proxy.b (0D_NOT_type -
      .def_property("b", &Mode3Proxy::b, &Mode3Proxy::set_b)
      // Mode3Proxy.c (0D_NOT_type -
      .def_property("c", &Mode3Proxy::c, &Mode3Proxy::set_c)
      // Mode3Proxy.x (0D_NOT_type -
      .def_property("x", &Mode3Proxy::x, &Mode3Proxy::set_x)
      // Mode3Proxy.y (0D_NOT_type -
      .def_property("y", &Mode3Proxy::y, &Mode3Proxy::set_y)

      .def("__repr__", [](const Mode3Proxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const Mode3Proxy& self) {
            return Mode3Proxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const Mode3Proxy& self, py::dict& memo) {
            return Mode3Proxy(self);
          })

      ;

  // 1D Mode3Proxy arrays are not used in structs/routines
  // 2D Mode3Proxy arrays are not used in structs/routines
  // 3D Mode3Proxy arrays are not used in structs/routines
}

// =============================================================================
// mode_info_struct
void init_mode_info_struct(py::module& m, py::class_<ModeInfoProxy>& cls) {
  cls.def(py::init<>())
      // ModeInfoProxy.stable (0D_NOT_logical - Is the mode stable?
      .def_property(
          "stable", &ModeInfoProxy::stable, &ModeInfoProxy::set_stable)
      // ModeInfoProxy.tune (0D_NOT_real - 'fractional' tune in radians
      .def_property("tune", &ModeInfoProxy::tune, &ModeInfoProxy::set_tune)
      // ModeInfoProxy.emit (0D_NOT_real - Emittance (unnormalized).
      .def_property("emit", &ModeInfoProxy::emit, &ModeInfoProxy::set_emit)
      // ModeInfoProxy.chrom (0D_NOT_real - Chromaticity.
      .def_property("chrom", &ModeInfoProxy::chrom, &ModeInfoProxy::set_chrom)
      // ModeInfoProxy.sigma (0D_NOT_real - Beam size.
      .def_property("sigma", &ModeInfoProxy::sigma, &ModeInfoProxy::set_sigma)
      // ModeInfoProxy.sigmap (0D_NOT_real - Beam divergence.
      .def_property(
          "sigmap", &ModeInfoProxy::sigmap, &ModeInfoProxy::set_sigmap)

      .def(
          "__repr__", [](const ModeInfoProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const ModeInfoProxy& self) {
            return ModeInfoProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const ModeInfoProxy& self, py::dict& memo) {
            return ModeInfoProxy(self);
          })

      ;

  // 1D ModeInfoProxy arrays are not used in structs/routines
  // 2D ModeInfoProxy arrays are not used in structs/routines
  // 3D ModeInfoProxy arrays are not used in structs/routines
}

// =============================================================================
// mad_energy_struct
void init_mad_energy_struct(py::module& m, py::class_<MadEnergyProxy>& cls) {
  cls.def(py::init<>())
      // MadEnergyProxy.total (0D_NOT_real -
      .def_property("total", &MadEnergyProxy::total, &MadEnergyProxy::set_total)
      // MadEnergyProxy.beta (0D_NOT_real - normalized velocity: v/c
      .def_property("beta", &MadEnergyProxy::beta, &MadEnergyProxy::set_beta)
      // MadEnergyProxy.gamma (0D_NOT_real - relativistic factor: 1/sqrt(1-beta^2)
      .def_property("gamma", &MadEnergyProxy::gamma, &MadEnergyProxy::set_gamma)
      // MadEnergyProxy.kinetic (0D_NOT_real - kinetic energy
      .def_property(
          "kinetic", &MadEnergyProxy::kinetic, &MadEnergyProxy::set_kinetic)
      // MadEnergyProxy.p0c (0D_NOT_real - particle momentum
      .def_property("p0c", &MadEnergyProxy::p0c, &MadEnergyProxy::set_p0c)
      // MadEnergyProxy.particle (0D_NOT_integer - particle species
      .def_property(
          "particle", &MadEnergyProxy::particle, &MadEnergyProxy::set_particle)

      .def(
          "__repr__",
          [](const MadEnergyProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MadEnergyProxy& self) {
            return MadEnergyProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const MadEnergyProxy& self, py::dict& memo) {
            return MadEnergyProxy(self);
          })

      ;

  // 1D MadEnergyProxy arrays are not used in structs/routines
  // 2D MadEnergyProxy arrays are not used in structs/routines
  // 3D MadEnergyProxy arrays are not used in structs/routines
}

// =============================================================================
// mad_map_struct
void init_mad_map_struct(py::module& m, py::class_<MadMapProxy>& cls) {
  cls.def(py::init<>())
      // MadMapProxy.k (1D_NOT_real - 0th order map.
      .def_property_readonly("k", &MadMapProxy::k)
      // MadMapProxy.r (2D_NOT_real - 1st order map.
      .def_property_readonly("r", &MadMapProxy::r)
      // MadMapProxy.t (3D_NOT_real - 2nd order map.
      .def_property_readonly("t", &MadMapProxy::t)

      .def("__repr__", [](const MadMapProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const MadMapProxy& self) {
            return MadMapProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const MadMapProxy& self, py::dict& memo) {
            return MadMapProxy(self);
          })

      ;

  // 1D MadMapProxy arrays are not used in structs/routines
  // 2D MadMapProxy arrays are not used in structs/routines
  // 3D MadMapProxy arrays are not used in structs/routines
}