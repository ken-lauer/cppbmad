#pragma once
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pybmad/arrays.hpp"
#include "pybmad/util.hpp"

namespace py = pybind11;

void init_Tao_routines_t(py::module& m);

struct PyTaoAllocateDataArray {
  int n_data;
  std::optional<bool> exact;
};
struct PyTaoAllocateV1Var {
  int n_v1;
  bool save_old;
};
struct PyTaoAllocateVarArray {
  bool default_good_user;
};
struct PyTaoBeamEmitCalc {
  double emit;
};
struct PyTaoBranchIndex {
  int ix_this;
};
struct PyTaoCalcDataAtSPts {
  double comp_sign;
};

struct PyTaoChangeEle {
  bool err_flag;
  bool update;
};
struct PyTaoChromCalcNeeded {
  std::string data_type;
  std::string data_source;
  bool do_chrom;
};
struct PyTaoClipCmd {
  double value1;
  double value2;
};
struct PyTaoCmdHistoryRecord {
  std::string cmd;
};

struct PyTaoCommand {
  bool err_is_fatal;
  bool err;
};
struct PyTaoConstraintTypeName {
  std::string datum_name;
};
struct PyTaoCurveEleRef {
  bool point_to_ele_ref;
};
struct PyTaoCurveIxUni {
  int ix_uni;
};
struct PyTaoCurveName {
  std::string curve_name;
};
struct PyTaoD2D1Name {
  std::string d2_d1_name;
};
struct PyTaoD2DataStuffit {
  std::string d2_name;
  int n_d1_data;
};
struct PyTaoDataCheck {
  bool err;
};
struct PyTaoDataSanityCheck {
  bool is_valid;
};
struct PyTaoDatumHasAssociatedEle {
  int has_associated_ele;
};
struct PyTaoDatumName {
  std::string datum_name;
};
struct PyTaoDrawCurveData {
  bool have_data;
};
struct PyTaoDrawEleForFloorPlan {
  double offset1;
  double offset2;
};
struct PyTaoDrawHistogramData {
  bool have_data;
};
struct PyTaoEleShapeInfo : public Tao::TaoEleShapeInfo {
  double y1;
  double y2;
  std::optional<int> ix_shape_min;
  PyTaoEleShapeInfo(
      Tao::TaoEleShapeInfo _base,
      double y1,
      double y2,
      std::optional<int> ix_shape_min)
      : Tao::TaoEleShapeInfo(std::move(_base)),
        y1(y1),
        y2(y2),
        ix_shape_min(ix_shape_min) {}
};
struct PyTaoEvaluateLatOrBeamData : public Tao::TaoEvaluateLatOrBeamData {
  std::string default_source;
  PyTaoEvaluateLatOrBeamData(
      Tao::TaoEvaluateLatOrBeamData _base,
      std::string default_source)
      : Tao::TaoEvaluateLatOrBeamData(std::move(_base)),
        default_source(default_source) {}
};
struct PyTaoEvaluateTune {
  double q_val;
};
struct PyTaoGraphName {
  std::string graph_name;
};
struct PyTaoInitBeamInUniverse {
  std::string track_start;
  std::string track_end;
  double comb_ds_save;
};
struct PyTaoInitDataInUniverse {
  int n_d2_add;
  std::optional<bool> keep_existing_data;
};
struct PyTaoInitLattice {
  std::string lat_file;
  bool err_flag;
};
struct PyTaoInitPlotting {
  std::string plot_file;
};
struct PyTaoInjectParticle {
  int ix_branch;
};

struct PyTaoIsValidName {
  std::string why_invalid;
  bool is_valid;
};
struct PyTaoKeyInfoToStr {
  int ix_key;
  int ix_min_key;
  int ix_max_key;
  std::string key_str;
  std::string header_str;
};
struct PyTaoLatEmitCalc {
  double emit;
};
struct PyTaoLatSigmaCalcNeeded {
  std::string data_type;
  std::string data_source;
  bool do_lat_sigma;
};
struct PyTaoLoadThisDatum {
  double datum_value;
  bool valid_value;
  std::optional<std::string> why_invalid;
};
struct PyTaoLocateElements : public Tao::TaoLocateElements {
  std::optional<bool> above_ubound_is_err;
  PyTaoLocateElements(
      Tao::TaoLocateElements _base,
      std::optional<bool> above_ubound_is_err)
      : Tao::TaoLocateElements(std::move(_base)),
        above_ubound_is_err(above_ubound_is_err) {}
};

struct PyTaoMerit {
  bool calc_ok;
  double this_merit;
};

struct PyTaoNextWord {
  std::string word;
  std::string line;
};
struct PyTaoOneTurnMapCalcNeeded {
  std::string data_type;
  std::string data_source;
  bool do_one_turn_map;
};

struct PyTaoOpenFile {
  int iunit;
  std::string file;
};

struct PyTaoOpenScratchFile {
  bool err;
  int iu;
};
struct PyTaoOptimizationStatus {
  std::string why_str;
};
struct PyTaoParamValueAtS : public Tao::TaoParamValueAtS {
  std::string dat_name;
  double value;
  PyTaoParamValueAtS(
      Tao::TaoParamValueAtS _base,
      std::string dat_name,
      double value)
      : Tao::TaoParamValueAtS(std::move(_base)),
        dat_name(dat_name),
        value(value) {}
};

struct PyTaoParseCommandArgs {
  bool error;
  std::optional<std::string> cmd_line;
};
struct PyTaoPointerToDatumEle : public Tao::TaoPointerToDatumEle {
  std::string ele_name;
  PyTaoPointerToDatumEle(Tao::TaoPointerToDatumEle _base, std::string ele_name)
      : Tao::TaoPointerToDatumEle(std::move(_base)), ele_name(ele_name) {}
};
struct PyTaoPointerToEleShape : public Tao::TaoPointerToEleShape {
  std::optional<int> ix_shape_min;
  PyTaoPointerToEleShape(
      Tao::TaoPointerToEleShape _base,
      std::optional<int> ix_shape_min)
      : Tao::TaoPointerToEleShape(std::move(_base)),
        ix_shape_min(ix_shape_min) {}
};

struct PyTaoPointerToUniverseStr {
  TaoUniverseProxy u;
  std::string string;
};
struct PyTaoRadIntCalcNeeded {
  std::string data_type;
  std::string data_source;
  bool do_rad_int;
};
struct PyTaoReExecute {
  std::string string;
  bool err;
};
struct PyTaoReadCmd {
  std::string which;
  std::string file;
};
struct PyTaoReadPhaseSpaceIndex {
  int ix_ps;
};
struct PyTaoRemoveBlankCharacters {
  std::string str;
};
struct PyTaoSetCalculateCmd {
  std::optional<std::string> switch_;
};
struct PyTaoSetDataCmd {
  std::optional<bool> silent;
};
struct PyTaoSetElementsCmd {
  bool update;
};
struct PyTaoSetFloorPlanAxisLabel {
  std::string which;
};
struct PyTaoSpinMatricesCalcNeeded {
  std::string data_type;
  std::string data_source;
  bool do_calc;
};
struct PyTaoSrdtCalcNeeded {
  std::string data_type;
  std::string data_source;
  int do_srdt;
};

struct PyTaoSubinUniNumber {
  std::string name_out;
  bool ok;
};
struct PyTaoToChangeNumber {
  std::string num_str;
  int n_size;
  std::string abs_or_rel;
  bool err;
};
struct PyTaoToInt {
  std::string str;
  int i_int;
  bool err;
};
struct PyTaoToPhaseAndCouplingReading
    : public Tao::TaoToPhaseAndCouplingReading {
  std::string why_invalid;
  PyTaoToPhaseAndCouplingReading(
      Tao::TaoToPhaseAndCouplingReading _base,
      std::string why_invalid)
      : Tao::TaoToPhaseAndCouplingReading(std::move(_base)),
        why_invalid(why_invalid) {}
};
struct PyTaoTooManyParticlesLost {
  bool no_beam;
};
struct PyTaoUniverseIndex {
  int i_this_uni;
};
struct PyTaoVar1Name {
  std::string var1_name;
};
struct PyTaoVarAttribName {
  std::string var_attrib_name;
};
struct PyTaoWaveCmd {
  bool err_flag;
};
struct PyTaoXScaleGraph {
  double x_min;
  double x_max;
  std::optional<bool> include_wall;
  std::optional<bool> have_scaled;
};
