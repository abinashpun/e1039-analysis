#ifndef __CALIB_DATA_H__
#define __CALIB_DATA_H__
#include <string>
#include <map>
class TH1;
class TH2;
class Tracklet;
class CalibParam;
class FitRTDist;

/// Class to hold calibration data.
class CalibData {
  static const int T_BIN =  30; ///< N of T bins for drawing.
  static const int R_BIN = 100; ///< N of R bins for drawing.

  static const int ST_ID_MIN = 1; // 1=D0, 2=D1, 3=D2, 4=D3p, 5=D3m, 6=D23, 7=Global
  static const int ST_ID_MAX = 7;

  CalibParam* cal_par;

  /// Historgrams about events and tracks
  TH1* h1_stat;
  TH1* h1_st_id;
  TH2* h2_ntrk;
  TH2* h2_nhit;
  TH2* h2_mom;
  TH2* h2_rchi2;
  TH2* h2_x0;
  TH2* h2_y0;
  TH2* h2_tx;
  TH2* h2_ty;
  TH2* h2_x_st;
  TH2* h2_y_st;

  /// Histograms for R-T
  TH1* h1_time_wide[99];
  TH1* h1_time [99];
  TH2* h2_tp_wp[99]; // track position, wire position
  TH2* h2_RT   [99];
  TH2* h2_TR   [99];
  TH2* h2_TR_ex[99];
  TH2* h2_Resi [99];
  TH2* h2_td_dd[99]; // x-axis: track distance, y-axis: drift distance

 public:
  CalibData();
  virtual ~CalibData();

  TH2* GetHistRT(const int ipl) const { return h2_RT[ipl]; }

  void Init(CalibParam* ptr);
  void FillEventInfo(const int rec_stat, const std::map<int, int> list_n_trk);
  void FillTracklet(const Tracklet* trk);
  void FillTrackletHits(const Tracklet* trk);
  void FillHit(const int det_id, const double drift_dist, const double tdc_time, const double track_pos, const double wire_pos);
  void DrawHistEvent(const std::string dir_out);
  void DrawHistHit  (const std::string dir_out);

 private:
  void DrawIn1D(TH2* h2, const std::string label, const std::string dir_out, int iy_lo=0, int iy_hi=0);
};

#endif // __CALIB_DATA_H__
