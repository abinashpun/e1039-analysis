#ifndef __CALIB_PARAM_H__
#define __CALIB_PARAM_H__
#include <string>
class TGraph;
class RTCurve;

/// Class to hold the calibration parameters.
class CalibParam {
  static const int N_PL = 30;
  static const double DT_RT; // time interval between R-T points.

  bool m_ana_pl[N_PL];

  bool fix_time_window;
  TGraph* m_gr_t2r_in[N_PL]; // Input R-T curve
  RTCurve* m_rtc[N_PL];

  double T_MIN[N_PL];
  double T_MAX[N_PL];
  double R_MAX[N_PL];
  double T0[N_PL]; ///< T0 of the fixed time window.
  double T1[N_PL]; ///< T1 of the fixed time window.

 public:
  CalibParam();
  virtual ~CalibParam();

  int GetNumPlanes() const { return N_PL; }

  void SetAnaPlanes(const bool d0, const bool d1, const bool d2, const bool d3p, const bool d3m);
  bool GetAnaPlane(const int i_pl) { return m_ana_pl[i_pl]; }

  bool TimeWindowIsFixed() const { return fix_time_window; }
  RTCurve* GetRTCurve(const int i_pl) const { return m_rtc[i_pl]; }

  double GetTMin(const int i_pl) const { return T_MIN[i_pl]; }
  double GetTMax(const int i_pl) const { return T_MAX[i_pl]; }
  double GetRMax(const int i_pl) const { return R_MAX[i_pl]; }
  double GetT0  (const int i_pl) const { return T0   [i_pl]; }
  double GetT1  (const int i_pl) const { return T1   [i_pl]; }

  TGraph*   GetGraphT2R (const int i_pl) const { return m_gr_t2r_in[i_pl]; }

  void Init(const int n_rt_pt);
  void  ReadRTParam(const std::string fname);
  void WriteRTParam(const std::string dir_name, const std::string fname);
  void WriteRTGraph(const std::string dir_name, const std::string fname);
  void ReadTimeWindow(const std::string fname);

  static double ZOfStationID(const int st_id);

 private:
};

#endif // __CALIB_PARAM_H__
