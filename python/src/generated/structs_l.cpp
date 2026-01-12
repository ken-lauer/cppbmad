#include "pybmad/generated/structs_l.hpp"
#include "bmad/generated/proxy.hpp"
#include "bmad/generated/to_string.hpp"
#include "bmad/to_string.hpp"
#include "pybmad/arrays.hpp"

using namespace Pybmad;
namespace py = pybind11;

// =============================================================================
// lat_ele_loc_struct
void init_lat_ele_loc_struct(py::module& m, py::class_<LatEleLocProxy>& cls) {
  cls.def(py::init<>())
      // LatEleLocProxy.ix_ele (0D_NOT_integer -
      .def_property(
          "ix_ele", &LatEleLocProxy::ix_ele, &LatEleLocProxy::set_ix_ele)
      // LatEleLocProxy.ix_branch (0D_NOT_integer -
      .def_property(
          "ix_branch",
          &LatEleLocProxy::ix_branch,
          &LatEleLocProxy::set_ix_branch)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return LatEleLocProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const LatEleLocProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatEleLocProxy& self) {
            return LatEleLocProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const LatEleLocProxy& self, py::dict& memo) {
            return LatEleLocProxy(self);
          })

      ;

  bind_FTypeArrayND<LatEleLocProxyArray1D>(m, "LatEleLocStructArray1D");
  bind_FTypeAlloc1D<LatEleLocProxyAlloc1D>(m, "LatEleLocStructAlloc1D");
  // 2D LatEleLocProxy arrays are not used in structs/routines
  // 3D LatEleLocProxy arrays are not used in structs/routines
}

// =============================================================================
// lat_ele_order1_struct
void init_lat_ele_order1_struct(
    py::module& m,
    py::class_<LatEleOrder1Proxy>& cls) {
  cls.def(py::init<>())
      // LatEleOrder1Proxy.ix_branch (0D_NOT_integer - Branch index
      .def_property(
          "ix_branch",
          &LatEleOrder1Proxy::ix_branch,
          &LatEleOrder1Proxy::set_ix_branch)
      // LatEleOrder1Proxy.ix_order (0D_NOT_integer - Order index. -1 -> Unique in lattice, 0 -> unique in branch.
      .def_property(
          "ix_order",
          &LatEleOrder1Proxy::ix_order,
          &LatEleOrder1Proxy::set_ix_order)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return LatEleOrder1ProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const LatEleOrder1Proxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatEleOrder1Proxy& self) {
            return LatEleOrder1Proxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const LatEleOrder1Proxy& self, py::dict& memo) {
            return LatEleOrder1Proxy(self);
          })

      ;

  bind_FTypeArrayND<LatEleOrder1ProxyArray1D>(m, "LatEleOrder1StructArray1D");
  bind_FTypeAlloc1D<LatEleOrder1ProxyAlloc1D>(m, "LatEleOrder1StructAlloc1D");
  // 2D LatEleOrder1Proxy arrays are not used in structs/routines
  // 3D LatEleOrder1Proxy arrays are not used in structs/routines
}

// =============================================================================
// lat_ele_order_array_struct
void init_lat_ele_order_array_struct(
    py::module& m,
    py::class_<LatEleOrderArrayProxy>& cls) {
  cls.def(py::init<>())
      // LatEleOrderArrayProxy.ele (1D_ALLOC_type -
      .def_property_readonly("ele", &LatEleOrderArrayProxy::ele)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) {
            return LatEleOrderArrayProxyAlloc1D(lbound, sz);
          },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def(
          "__repr__",
          [](const LatEleOrderArrayProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatEleOrderArrayProxy& self) {
            return LatEleOrderArrayProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const LatEleOrderArrayProxy& self, py::dict& memo) {
            return LatEleOrderArrayProxy(self);
          })

      ;

  bind_FTypeArrayND<LatEleOrderArrayProxyArray1D>(
      m, "LatEleOrderArrayStructArray1D");
  bind_FTypeAlloc1D<LatEleOrderArrayProxyAlloc1D>(
      m, "LatEleOrderArrayStructAlloc1D");
  // 2D LatEleOrderArrayProxy arrays are not used in structs/routines
  // 3D LatEleOrderArrayProxy arrays are not used in structs/routines
}

// =============================================================================
// lat_ele_order_struct
void init_lat_ele_order_struct(
    py::module& m,
    py::class_<LatEleOrderProxy>& cls) {
  cls.def(py::init<>())
      // LatEleOrderProxy.branch (1D_ALLOC_type -
      .def_property_readonly("branch", &LatEleOrderProxy::branch)

      .def(
          "__repr__",
          [](const LatEleOrderProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatEleOrderProxy& self) {
            return LatEleOrderProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const LatEleOrderProxy& self, py::dict& memo) {
            return LatEleOrderProxy(self);
          })

      ;

  // 1D LatEleOrderProxy arrays are not used in structs/routines
  // 2D LatEleOrderProxy arrays are not used in structs/routines
  // 3D LatEleOrderProxy arrays are not used in structs/routines
}

// =============================================================================
// lat_param_struct
void init_lat_param_struct(py::module& m, py::class_<LatParamProxy>& cls) {
  cls.def(py::init<>())
      // LatParamProxy.n_part (0D_NOT_real - Particles/bunch (for BeamBeam elements).
      .def_property(
          "n_part", &LatParamProxy::n_part, &LatParamProxy::set_n_part)
      // LatParamProxy.total_length (0D_NOT_real - total_length of branch. Warning: branch may not start at s = 0.
      .def_property(
          "total_length",
          &LatParamProxy::total_length,
          &LatParamProxy::set_total_length)
      // LatParamProxy.unstable_factor (0D_NOT_real - If positive: Growth rate/turn if unstable in closed branches or |orbit-aperture|/aperture if particle hits wall. Zero otherwise.
      .def_property(
          "unstable_factor",
          &LatParamProxy::unstable_factor,
          &LatParamProxy::set_unstable_factor)
      // LatParamProxy.t1_with_RF (2D_NOT_real - Full 1-turn matrix with RF on.
      .def_property_readonly("t1_with_RF", &LatParamProxy::t1_with_RF)
      // LatParamProxy.t1_no_RF (2D_NOT_real - Full 1-turn matrix with RF off.
      .def_property_readonly("t1_no_RF", &LatParamProxy::t1_no_RF)
      // LatParamProxy.spin_tune (0D_NOT_real - Closed orbit spin tune.
      .def_property(
          "spin_tune", &LatParamProxy::spin_tune, &LatParamProxy::set_spin_tune)
      // LatParamProxy.particle (0D_NOT_integer - Reference particle: positron$, electron$, etc. Call lattice_bookkeeper if this is changed.
      .def_property(
          "particle", &LatParamProxy::particle, &LatParamProxy::set_particle)
      // LatParamProxy.default_tracking_species (0D_NOT_integer - Default particle type to use in tracking.
      .def_property(
          "default_tracking_species",
          &LatParamProxy::default_tracking_species,
          &LatParamProxy::set_default_tracking_species)
      // LatParamProxy.geometry (0D_NOT_integer - open$ or closed$
      .def_property(
          "geometry", &LatParamProxy::geometry, &LatParamProxy::set_geometry)
      // LatParamProxy.ixx (0D_NOT_integer - Integer for general use
      .def_property("ixx", &LatParamProxy::ixx, &LatParamProxy::set_ixx)
      // LatParamProxy.stable (0D_NOT_logical - is closed lat stable?
      .def_property(
          "stable", &LatParamProxy::stable, &LatParamProxy::set_stable)
      // LatParamProxy.live_branch (0D_NOT_logical - Should tracking be done on the branch?
      .def_property(
          "live_branch",
          &LatParamProxy::live_branch,
          &LatParamProxy::set_live_branch)
      // LatParamProxy.g1_integral (0D_NOT_real - Approximate |g| (bending strength) integral of branch.
      .def_property(
          "g1_integral",
          &LatParamProxy::g1_integral,
          &LatParamProxy::set_g1_integral)
      // LatParamProxy.g2_integral (0D_NOT_real - Approximate g^2 integral of branch.
      .def_property(
          "g2_integral",
          &LatParamProxy::g2_integral,
          &LatParamProxy::set_g2_integral)
      // LatParamProxy.g3_integral (0D_NOT_real - Approximate g^2 integral of branch.
      .def_property(
          "g3_integral",
          &LatParamProxy::g3_integral,
          &LatParamProxy::set_g3_integral)
      // LatParamProxy.bookkeeping_state (0D_NOT_type - Overall status for the branch.
      .def_property(
          "bookkeeping_state",
          &LatParamProxy::bookkeeping_state,
          &LatParamProxy::set_bookkeeping_state)
      // LatParamProxy.beam_init (0D_NOT_type - For beam initialization.
      .def_property(
          "beam_init", &LatParamProxy::beam_init, &LatParamProxy::set_beam_init)

      .def(
          "__repr__", [](const LatParamProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatParamProxy& self) {
            return LatParamProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const LatParamProxy& self, py::dict& memo) {
            return LatParamProxy(self);
          })

      ;

  // 1D LatParamProxy arrays are not used in structs/routines
  // 2D LatParamProxy arrays are not used in structs/routines
  // 3D LatParamProxy arrays are not used in structs/routines
}

// =============================================================================
// lat_struct
void init_lat_struct(py::module& m, py::class_<LatProxy>& cls) {
  cls.def(py::init<>())
      // LatProxy.use_name (0D_NOT_character - Name of lat given by USE statement
      .def_property("use_name", &LatProxy::use_name, &LatProxy::set_use_name)
      // LatProxy.lattice (0D_NOT_character - Lattice
      .def_property("lattice", &LatProxy::lattice, &LatProxy::set_lattice)
      // LatProxy.machine (0D_NOT_character - Name of the machine the lattice is for ('LHC', etc).
      .def_property("machine", &LatProxy::machine, &LatProxy::set_machine)
      // LatProxy.input_file_name (0D_NOT_character - Name of the lattice input file
      .def_property(
          "input_file_name",
          &LatProxy::input_file_name,
          &LatProxy::set_input_file_name)
      // LatProxy.title (0D_NOT_character - General title
      .def_property("title", &LatProxy::title, &LatProxy::set_title)
      // LatProxy.print_str (1D_ALLOC_character - Saved print statements.
      .def_property_readonly("print_str", &LatProxy::print_str)
      // LatProxy.constant (1D_ALLOC_type - Constants defined in the lattice
      .def_property_readonly("constant", &LatProxy::constant)
      // LatProxy.a (0D_PTR_type - Tunes (fractional part), etc.
      .def_property("a", &LatProxy::a, &LatProxy::set_a)
      // LatProxy.b (0D_PTR_type - Tunes (fractional part), etc.
      .def_property("b", &LatProxy::b, &LatProxy::set_b)
      // LatProxy.z (0D_PTR_type - Tunes (fractional part), etc.
      .def_property("z", &LatProxy::z, &LatProxy::set_z)
      // LatProxy.param (0D_PTR_type - Parameters
      .def_property("param", &LatProxy::param, &LatProxy::set_param)
      // LatProxy.lord_state (0D_NOT_type - lord bookkeeping status.
      .def_property(
          "lord_state", &LatProxy::lord_state, &LatProxy::set_lord_state)
      // LatProxy.ele_init (0D_NOT_type - For use by any program
      .def_property("ele_init", &LatProxy::ele_init, &LatProxy::set_ele_init)
      // LatProxy.ele (1D_PTR_type - Array of elements [=> branch(0)].
      .def_property_readonly("ele", &LatProxy::ele)
      // LatProxy.branch (1D_ALLOC_type - Branch(0:) array
      .def_property_readonly("branch", &LatProxy::branch)
      // LatProxy.control (1D_ALLOC_type - Control list
      .def_property_readonly("control", &LatProxy::control)
      // LatProxy.particle_start (0D_PTR_type - Starting particle_coords.
      .def_property(
          "particle_start",
          &LatProxy::particle_start,
          &LatProxy::set_particle_start)
      // LatProxy.beam_init (0D_NOT_type - Beam initialization.
      .def_property("beam_init", &LatProxy::beam_init, &LatProxy::set_beam_init)
      // LatProxy.pre_tracker (0D_NOT_type - For OPAL/IMPACT-T
      .def_property(
          "pre_tracker", &LatProxy::pre_tracker, &LatProxy::set_pre_tracker)
      // LatProxy.custom (1D_ALLOC_real - Custom attributes.
      .def_property_readonly("custom", &LatProxy::custom)
      // LatProxy.version (0D_NOT_integer - Version number
      .def_property("version", &LatProxy::version, &LatProxy::set_version)
      // LatProxy.n_ele_track (0D_PTR_integer - Number of lat elements to track through.
      .def_property(
          "n_ele_track", &LatProxy::n_ele_track, &LatProxy::set_n_ele_track)
      // LatProxy.n_ele_max (0D_PTR_integer - Index of last valid element in %ele(:) array
      .def_property("n_ele_max", &LatProxy::n_ele_max, &LatProxy::set_n_ele_max)
      // LatProxy.n_control_max (0D_NOT_integer - Last index used in control_array
      .def_property(
          "n_control_max",
          &LatProxy::n_control_max,
          &LatProxy::set_n_control_max)
      // LatProxy.n_ic_max (0D_NOT_integer - Last index used in ic_array
      .def_property("n_ic_max", &LatProxy::n_ic_max, &LatProxy::set_n_ic_max)
      // LatProxy.input_taylor_order (0D_NOT_integer - As set in the input file
      .def_property(
          "input_taylor_order",
          &LatProxy::input_taylor_order,
          &LatProxy::set_input_taylor_order)
      // LatProxy.ic (1D_ALLOC_integer - Index to %control(:) from slaves.
      .def_property_readonly("ic", &LatProxy::ic)
      // LatProxy.photon_type (0D_NOT_integer - Or coherent$. For X-ray simulations.
      .def_property(
          "photon_type", &LatProxy::photon_type, &LatProxy::set_photon_type)
      // LatProxy.creation_hash (0D_NOT_integer - Set by bmad_parser. creation_hash will vary if any of the lattice files are modified.
      .def_property(
          "creation_hash",
          &LatProxy::creation_hash,
          &LatProxy::set_creation_hash)
      // LatProxy.ramper_slave_bookkeeping (0D_NOT_integer -
      .def_property(
          "ramper_slave_bookkeeping",
          &LatProxy::ramper_slave_bookkeeping,
          &LatProxy::set_ramper_slave_bookkeeping)
      .def_static(
          "new_array1d",
          [](int sz, int lbound) { return LatProxyAlloc1D(lbound, sz); },
          py::arg("sz"),
          py::arg("lbound") = 1)

      .def("__repr__", [](const LatProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LatProxy& self) {
            return LatProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const LatProxy& self, py::dict& memo) { return LatProxy(self); })

      ;

  bind_FTypeArrayND<LatProxyArray1D>(m, "LatStructArray1D");
  bind_FTypeAlloc1D<LatProxyAlloc1D>(m, "LatStructAlloc1D");
  // 2D LatProxy arrays are not used in structs/routines
  // 3D LatProxy arrays are not used in structs/routines
}

// =============================================================================
// linac_normal_mode_struct
void init_linac_normal_mode_struct(
    py::module& m,
    py::class_<LinacNormalModeProxy>& cls) {
  cls.def(py::init<>())
      // LinacNormalModeProxy.i2_E4 (0D_NOT_real - Integral: g^2 * gamma^4
      .def_property(
          "i2_E4",
          &LinacNormalModeProxy::i2_E4,
          &LinacNormalModeProxy::set_i2_E4)
      // LinacNormalModeProxy.i3_E7 (0D_NOT_real - Integral: g^3 * gamma^7
      .def_property(
          "i3_E7",
          &LinacNormalModeProxy::i3_E7,
          &LinacNormalModeProxy::set_i3_E7)
      // LinacNormalModeProxy.i5a_E6 (0D_NOT_real - Integral: (g^3 * H_a) * gamma^6
      .def_property(
          "i5a_E6",
          &LinacNormalModeProxy::i5a_E6,
          &LinacNormalModeProxy::set_i5a_E6)
      // LinacNormalModeProxy.i5b_E6 (0D_NOT_real - Integral: (g^3 * H_b) * gamma^6
      .def_property(
          "i5b_E6",
          &LinacNormalModeProxy::i5b_E6,
          &LinacNormalModeProxy::set_i5b_E6)
      // LinacNormalModeProxy.sig_E1 (0D_NOT_real - Energy spread after 1 pass (eV)
      .def_property(
          "sig_E1",
          &LinacNormalModeProxy::sig_E1,
          &LinacNormalModeProxy::set_sig_E1)
      // LinacNormalModeProxy.a_emittance_end (0D_NOT_real - a mode emittance at end of linac
      .def_property(
          "a_emittance_end",
          &LinacNormalModeProxy::a_emittance_end,
          &LinacNormalModeProxy::set_a_emittance_end)
      // LinacNormalModeProxy.b_emittance_end (0D_NOT_real - b mode emittance at end of linac
      .def_property(
          "b_emittance_end",
          &LinacNormalModeProxy::b_emittance_end,
          &LinacNormalModeProxy::set_b_emittance_end)

      .def(
          "__repr__",
          [](const LinacNormalModeProxy& self) { return to_string(self); })

      .def(
          "__copy__",
          [](const LinacNormalModeProxy& self) {
            return LinacNormalModeProxy(self); // under-the-hood fortran copy
          })
      .def(
          "__deepcopy__",
          [](const LinacNormalModeProxy& self, py::dict& memo) {
            return LinacNormalModeProxy(self);
          })

      ;

  // 1D LinacNormalModeProxy arrays are not used in structs/routines
  // 2D LinacNormalModeProxy arrays are not used in structs/routines
  // 3D LinacNormalModeProxy arrays are not used in structs/routines
}