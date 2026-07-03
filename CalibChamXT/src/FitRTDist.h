#ifndef __FIT_RT_DIST_H__
#define __FIT_RT_DIST_H__
#include <string>
#include <vector>
class TH2;
class RTCurve;

/// Class for fitting R-T histogram
class FitRTDist {
 protected:
  bool m_fix_t1t0;
  double m_t1; ///< Fixed T1 value.  Valid only when m_fix_t1t0 = true.
  double m_t0; ///< Fixed T0 value.  Valid only when m_fix_t1t0 = true.

  static int    m_n_pt; ///< N of R-T points.
  static int    m_n_par; ///< N of free parameters
  static TH2*   m_h2_RT;
  static double m_t_min;
  static double m_t_max;
  static double m_r_max;

  int m_verb;

 public:
  FitRTDist();
  virtual ~FitRTDist();

  static double FCN(const double* xx);
  static void CalcChi2(const double* pars, double& chi2, int& ndf);

  void FixT1T0(const double t1, const double t0);
  int  DoFit(const int n_pt, TH2* h2, double r_max, TGraph* gr_init, RTCurve* rtc);

  void Verbosity(const int verb) { m_verb = verb; }
  int  Verbosity() const  { return m_verb; }

 private:
  static void SetRTCurve(const double* pars, RTCurve* rtc);
};

#endif // __FIT_RT_DIST_H__
