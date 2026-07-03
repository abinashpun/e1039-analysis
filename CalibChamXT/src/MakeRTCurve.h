#ifndef _MAKE_RT_CURVE__H_
#define _MAKE_RT_CURVE__H_
#include <map>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <fun4all/SubsysReco.h>
class CalibParam;
class CalibData;

class MakeRTCurve {
  typedef enum {
    ST_ID_D3P = 4,
    ST_ID_D3M = 5,
    ST_ID_D23 = 6
  } StationID_t;
  static const int N_RT_PT = 9; // N of points of R-T curve.
  int m_iter;
  std::string m_dir_name_out;
  
  CalibParam cal_par;
  CalibData  cal_dat;

  unsigned int m_n_evt_all;
  unsigned int m_n_evt_ana;
  unsigned int m_n_trk_all;
  unsigned int m_n_trk_ana;

  int m_verb;

 public:
  MakeRTCurve(const int iter);
  virtual ~MakeRTCurve();

  void Init();
  void AnalyzeListOfFiles(const char* fn_list);
  void AnalyzeFile(const char* fname);

  void DrawHistEvent();
  void DrawHistHit();
  void ExtractRT();
  void WriteRT();
  void DrawCalibResult();

  void Verbosity(const int verb) { m_verb = verb; }
  int  Verbosity() const  { return m_verb; }

 private:
  void EvalD2XY(const Tracklet* trk, double& x, double& y);
  void EvalD3XY(const Tracklet* trk, double& x, double& y);
};

#endif // _MAKE_RT_CURVE__H_
