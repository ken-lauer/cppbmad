#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Bmad_routines_e(py::module& m);

struct PyEAccelField {
  double field;
};
struct PyEleFullName {
  std::string str;
};
struct PyEleHasConstantDsDtRef {
  bool is_const;
};

struct PyEleHasNonzeroKick {
  EleProxy ele;
  bool has_kick;
};
struct PyEleHasNonzeroOffset {
  bool has_offset;
};
struct PyEleLocName {
  std::string str;
};
struct PyEleNametableIndex {
  int ix_nt;
};
struct PyEleRfStepIndex {
  int ix_step;
};
struct PyEleUniqueName {
  std::string unique_name;
};
struct PyEleValueHasChanged {
  bool has_changed;
};

struct PyEmFieldDerivatives {
  EmFieldProxy dfield;
  double s_pos;
  bool local_ref_frame;
  std::optional<bool> grid_allow_s_out_of_bounds;
  std::optional<double> rf_time;
};
struct PyEnteringElement {
  bool is_entering;
};
struct PyEqAcKicker {
  bool is_eq;
};
struct PyEqAcKickerFreq {
  bool is_eq;
};
struct PyEqAcKickerTime {
  bool is_eq;
};
struct PyEqAnormalMode {
  bool is_eq;
};
struct PyEqApertureParam {
  bool is_eq;
};
struct PyEqAperturePoint {
  bool is_eq;
};
struct PyEqApertureScan {
  bool is_eq;
};
struct PyEqBeam {
  bool is_eq;
};
struct PyEqBeamInit {
  bool is_eq;
};
struct PyEqBmadCommon {
  bool is_eq;
};
struct PyEqBookkeepingState {
  bool is_eq;
};
struct PyEqBpmPhaseCoupling {
  bool is_eq;
};
struct PyEqBranch {
  bool is_eq;
};
struct PyEqBunch {
  bool is_eq;
};
struct PyEqBunchParams {
  bool is_eq;
};
struct PyEqCartesianMap {
  bool is_eq;
};
struct PyEqCartesianMapTerm {
  bool is_eq;
};
struct PyEqCartesianMapTerm1 {
  bool is_eq;
};
struct PyEqComplexTaylor {
  bool is_eq;
};
struct PyEqComplexTaylorTerm {
  bool is_eq;
};
struct PyEqControl {
  bool is_eq;
};
struct PyEqControlRamp1 {
  bool is_eq;
};
struct PyEqControlVar1 {
  bool is_eq;
};
struct PyEqController {
  bool is_eq;
};
struct PyEqCoord {
  bool is_eq;
};
struct PyEqCoordArray {
  bool is_eq;
};
struct PyEqCylindricalMap {
  bool is_eq;
};
struct PyEqCylindricalMapTerm {
  bool is_eq;
};
struct PyEqCylindricalMapTerm1 {
  bool is_eq;
};
struct PyEqEle {
  bool is_eq;
};
struct PyEqEllipseBeamInit {
  bool is_eq;
};
struct PyEqEmField {
  bool is_eq;
};
struct PyEqEmTaylor {
  bool is_eq;
};
struct PyEqEmTaylorTerm {
  bool is_eq;
};
struct PyEqExpressionAtom {
  bool is_eq;
};
struct PyEqFloorPosition {
  bool is_eq;
};
struct PyEqGenGrad1 {
  bool is_eq;
};
struct PyEqGenGradMap {
  bool is_eq;
};
struct PyEqGridBeamInit {
  bool is_eq;
};
struct PyEqGridField {
  bool is_eq;
};
struct PyEqGridFieldPt {
  bool is_eq;
};
struct PyEqGridFieldPt1 {
  bool is_eq;
};
struct PyEqHighEnergySpaceCharge {
  bool is_eq;
};
struct PyEqInterval1Coef {
  bool is_eq;
};
struct PyEqKvBeamInit {
  bool is_eq;
};
struct PyEqLat {
  bool is_eq;
};
struct PyEqLatEleLoc {
  bool is_eq;
};
struct PyEqLatParam {
  bool is_eq;
};
struct PyEqLinacNormalMode {
  bool is_eq;
};
struct PyEqMode3 {
  bool is_eq;
};
struct PyEqModeInfo {
  bool is_eq;
};
struct PyEqNormalModes {
  bool is_eq;
};
struct PyEqPhotonElement {
  bool is_eq;
};
struct PyEqPhotonMaterial {
  bool is_eq;
};
struct PyEqPhotonReflectSurface {
  bool is_eq;
};
struct PyEqPhotonReflectTable {
  bool is_eq;
};
struct PyEqPhotonTarget {
  bool is_eq;
};
struct PyEqPixelDetec {
  bool is_eq;
};
struct PyEqPixelPt {
  bool is_eq;
};
struct PyEqPreTracker {
  bool is_eq;
};
struct PyEqRadInt1 {
  bool is_eq;
};
struct PyEqRadIntAllEle {
  bool is_eq;
};
struct PyEqRadIntBranch {
  bool is_eq;
};
struct PyEqRadMap {
  bool is_eq;
};
struct PyEqRadMapEle {
  bool is_eq;
};
struct PyEqRamperLord {
  bool is_eq;
};
struct PyEqSpaceChargeCommon {
  bool is_eq;
};
struct PyEqSpinPolar {
  bool is_eq;
};
struct PyEqSpline {
  bool is_eq;
};
struct PyEqStrongBeam {
  bool is_eq;
};
struct PyEqSurfaceCurvature {
  bool is_eq;
};
struct PyEqSurfaceDisplacement {
  bool is_eq;
};
struct PyEqSurfaceDisplacementPt {
  bool is_eq;
};
struct PyEqSurfaceHMisalign {
  bool is_eq;
};
struct PyEqSurfaceHMisalignPt {
  bool is_eq;
};
struct PyEqSurfaceSegmented {
  bool is_eq;
};
struct PyEqSurfaceSegmentedPt {
  bool is_eq;
};
struct PyEqTargetPoint {
  bool is_eq;
};
struct PyEqTaylor {
  bool is_eq;
};
struct PyEqTaylorTerm {
  bool is_eq;
};
struct PyEqTrack {
  bool is_eq;
};
struct PyEqTrackPoint {
  bool is_eq;
};
struct PyEqTwiss {
  bool is_eq;
};
struct PyEqWake {
  bool is_eq;
};
struct PyEqWakeLr {
  bool is_eq;
};
struct PyEqWakeLrMode {
  bool is_eq;
};
struct PyEqWakeSr {
  bool is_eq;
};
struct PyEqWakeSrMode {
  bool is_eq;
};
struct PyEqWakeSrZLong {
  bool is_eq;
};
struct PyEqWall3d {
  bool is_eq;
};
struct PyEqWall3dSection {
  bool is_eq;
};
struct PyEqWall3dVertex {
  bool is_eq;
};
struct PyEqXyDisp {
  bool is_eq;
};
struct PyEqualSignHere {
  std::string delim;
  bool is_here;
};
struct PyEquivalentTaylorAttributes {
  bool equiv;
};
struct PyEtdiv {
  double A;
  double B;
  double C;
  double D;
  double E;
  double F;
};
struct PyExpectOneOf {
  std::string delim;
  bool delim_found;
  bool is_ok;
};
