#include <cassert>
#include <TSpline.h>
#include <TGraph.h>
#include "RTCurve.h"
using namespace std;

RTCurve::RTCurve(const int n_pt)
  : m_n_pt(n_pt)
  , m_array_r(new double[n_pt])
  , m_array_t(new double[n_pt])
  , m_t2r_changed(true)
  , m_t2dr_changed(true)
  , m_spl_t2r(0)
  , m_gr_t2dr(0)
{ 
  ;
}

RTCurve::~RTCurve()
{
  delete m_array_t;
  delete m_array_r;
  if (m_spl_t2r) delete m_spl_t2r;
  if (m_gr_t2dr) delete m_gr_t2dr;
}

void RTCurve::SetPoint(const int i_pt, const double r, const double t)
{
  assert(0 <= i_pt && i_pt < m_n_pt);
  m_array_r[i_pt] = r;
  m_array_t[i_pt] = t;
  m_t2r_changed = true;
}

void RTCurve::GetPoint(const int i_pt, double& r, double& t) const
{
  assert(0 <= i_pt && i_pt < m_n_pt);
  r = m_array_r[i_pt];
  t = m_array_t[i_pt];
}

void RTCurve::SetRWidth(const double dr)
{
  m_dr = dr; 
  m_t2dr_changed = true;
}

double RTCurve::GetRWidth() const
{
  return m_dr;
}

double RTCurve::EvalR(const double t)
{
  if (t < GetTMin()) return GetRMax();
  if (t > GetTMax()) return 0;
  double r = GetT2RSpline()->Eval(t);
  if (r > GetRMax()) return GetRMax();
  if (r < 0        ) return 0;
  return r;
}

TSpline3* RTCurve::GetT2RSpline()
{
  if (m_t2r_changed) {
    if (m_spl_t2r) delete m_spl_t2r;
    m_spl_t2r = CreateT2RSpline();
    m_t2r_changed = false;
  }
  return m_spl_t2r;
}

TGraph* RTCurve::GetT2DRGraph()
{
  if (m_t2dr_changed) {
    if (m_gr_t2dr) delete m_gr_t2dr;
    m_gr_t2dr = CreateT2DRGraph();
    m_t2dr_changed = false;
  }
  return m_gr_t2dr;
}

TSpline3* RTCurve::CreateT2RSpline() const
{
  return new TSpline3("spl", m_array_t, m_array_r, m_n_pt, "b1e1", 0, 0);
}

TGraph* RTCurve::CreateT2DRGraph() const
{
  TGraph* gr = new TGraph();
  gr->SetPoint(0, GetTMin(), m_dr);
  return gr;
}
