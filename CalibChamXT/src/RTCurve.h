#ifndef __RT_CURVE_H__
#define __RT_CURVE_H__
class TSpline3;
class TGraph;

/// Class to represent R-T curve.
class RTCurve {
  int m_n_pt; ///< N of R-T points
  double* m_array_r;
  double* m_array_t;
  double  m_dr; ///< Width in R

  bool m_t2r_changed;
  bool m_t2dr_changed;
  TSpline3* m_spl_t2r; // Must access via GetT2RSpline().
  TGraph*   m_gr_t2dr; // Must access via GetT2DRGraph().

 public:
  RTCurve(const int n_pt);
  virtual ~RTCurve();

  int GetNumPoints() const { return m_n_pt; }
  void SetPoint(const int i_pt, const double r, const double t);
  void GetPoint(const int i_pt, double& r, double& t) const;

  void   SetRWidth(const double dr);
  double GetRWidth() const;

  double GetTMin() const { return m_array_t[0]; }
  double GetTMax() const { return m_array_t[m_n_pt-1]; }
  double GetRMax() const { return m_array_r[0]; }

  double EvalR(const double t);

  TSpline3* GetT2RSpline();
  TGraph*   GetT2DRGraph();

 private:
  TSpline3* CreateT2RSpline() const;
  TGraph*   CreateT2DRGraph() const;
};

#endif // __RT_CURVE_H__
