#include <iostream>
#include <sstream>
#include <TMath.h>
#include <TH2D.h>
#include <TSpline.h>
#include <TGraph.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include "RTCurve.h"
#include "FitRTDist.h"
using namespace std;

int    FitRTDist::m_n_pt  = 0;
int    FitRTDist::m_n_par = 0;
TH2*   FitRTDist::m_h2_RT = 0;
double FitRTDist::m_t_min = 0;
double FitRTDist::m_t_max = 0;
double FitRTDist::m_r_max = 0;

FitRTDist::FitRTDist()
  : m_fix_t1t0(false)
  , m_t1(0)
  , m_t0(0)
  , m_verb(0)
{ 
  ;
}

FitRTDist::~FitRTDist()
{
  //if (m_minuit) delete m_minuit;
}

double FitRTDist::FCN(const double* xx)
{
  double chi2;
  int ndf;
  CalcChi2(xx, chi2, ndf);
  return chi2;
}

/**
 * N of parameters = m_n_pt + 1
 * pars[0] = T1
 *     [1] = T0
 *     [2] = DR (constant over T)
 *     [3]...[m_n_pt] = R at each non-edge T points
 *
 * Note that the T position of each point is not fixed but depends on T1 & T0.
 */
void FitRTDist::CalcChi2(const double* pars, double& chi2, int& ndf)
{
  RTCurve rtc(m_n_pt);
  SetRTCurve(pars, &rtc);
  double dr = rtc.GetRWidth();
  
  chi2 = 0;
  ndf  = 0;
  double r_width = m_h2_RT->GetYaxis()->GetBinWidth(1);
  int n_r = m_h2_RT->GetNbinsY();
  for (int it = m_h2_RT->GetNbinsX(); it > 0; it--) {
    double cont_t  = m_h2_RT->Integral(it, it, 1, n_r);
    if (cont_t < 10) continue;
    double t_cent = m_h2_RT->GetXaxis()->GetBinCenter(it);
    double r_func = rtc.EvalR(t_cent);
    for (int ir = n_r; ir > 0; ir--) {
      double cont = m_h2_RT->GetBinContent(it, ir);
      if (cont <= 0) continue;
      double r_cent  = m_h2_RT->GetYaxis()->GetBinCenter(ir);
      double cont_func = cont_t * r_width * exp( -pow( (r_cent-r_func)/dr, 2 )/2 ) / (sqrt(2*TMath::Pi()) * dr);
      double chi2p = pow(cont - cont_func, 2) / cont;
      //cout << "XXX " << it << " " << ir << " : " << dr << " " << cont_t << " " << cont_func << " " << cont << " " << chi2p << endl;
      chi2 += chi2p;
      ndf++;
    }
  }
  ndf -= m_n_par;
}

void FitRTDist::FixT1T0(const double t1, const double t0)
{
  m_fix_t1t0 = true;
  m_t1 = t1;
  m_t0 = t0;
}

int FitRTDist::DoFit(const int n_pt, TH2* h2, double r_max, TGraph* gr_init, RTCurve* rtc)
{
  const double DR = 0.04;
  m_n_pt  = n_pt;
  m_h2_RT = h2;
  m_r_max = r_max;
  m_t_min = h2->GetXaxis()->GetXmin();
  m_t_max = h2->GetXaxis()->GetXmax();

  ROOT::Math::Minimizer* m_mini = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  m_mini->SetMaxFunctionCalls(1000000);
  m_mini->SetTolerance(0.001);
  m_mini->SetPrintLevel(m_verb > 0  ?  m_verb - 1  :  0);

  double t1, t0;
  if (m_fix_t1t0) {
    t1 = m_t1;
    t0 = m_t0;
    m_mini->SetFixedVariable(0, "T1", t1);
    m_mini->SetFixedVariable(1, "T1", t0);
  } else {
    t1 = gr_init->GetX()[0];
    t0 = gr_init->GetX()[gr_init->GetN()-1];
    m_mini->SetLimitedVariable(0, "T1", t1, 0.1, m_t_min, m_t_max);
    m_mini->SetLimitedVariable(1, "T0", t0, 0.1, m_t_min, m_t_max);
  }
  m_mini->SetLimitedVariable(2, "DR", DR, 0.01, 0, 1);
  if (m_verb > 0) cout << "    T: " << t1 << "   " << t0 << "   " << DR << endl;

  ostringstream oss;

  for (int i_pt = 1; i_pt < m_n_pt - 1; i_pt++) {
    double t = t1 + (t0 - t1) * (i_pt + 1) / m_n_pt;
    double r = gr_init->Eval(t);
    oss.str("");   oss << "R" << i_pt+1;
    double r_lo = r - 3 * DR; // 3 sigma at max
    double r_hi = r + 3 * DR; // 3 sigma at max
    if (r_lo < 0    ) r_lo = 0;
    if (r_hi > r_max) r_hi = r_max;
    m_mini->SetLimitedVariable(i_pt+2, oss.str().c_str(), r, 0.001, r_lo, r_hi);
    if (m_verb > 0) cout << "    " << i_pt+1 << ": " << t << " --- " << r_lo << " < " << r << " < " << r_hi << "\n";
  }

  m_n_par = m_mini->NFree();
  if (m_verb > 0) cout << "    N_par " << m_n_par << "\n";

  ROOT::Math::Functor functor(&FitRTDist::FCN, m_n_par);
  m_mini->SetFunction(functor);

  bool mini_result = m_mini->Minimize();
  if (! mini_result) {
    cerr << "    WARNING:  Bad fit result (ret = " << mini_result << ").  Please inspect it." << endl;
    //<< "Please fix the problem.  Abort" << endl;
    //exit(1);
  }

  if (m_verb > 1) {
    double min_value = m_mini->MinValue();
    cout << "    MinValue = " << min_value << "\n";
  }
  const double* list_val = m_mini->X();
  const double* list_err = m_mini->Errors();
  //double c_ij = m_mini->CovMatrix(i, j);

  if (m_verb > 1) {
    cout << "    T: " << list_val[0] << "   " << list_val[1] << "   " << list_val[2] << "\n";
    for (int i_pt = 1; i_pt < m_n_pt - 1; i_pt++) {
      cout << "    " << i_pt+1 << ": " << list_val[i_pt + 2] << "\n";
    }
  }

  double chi2;
  int    ndf ;
  CalcChi2(list_val, chi2, ndf);
  cout << "    chi2 / ndf = " << chi2 << " / " << ndf << " = " << chi2/ndf << endl;

  SetRTCurve(list_val, rtc);

  delete m_mini;
  return 0;
}

void FitRTDist::SetRTCurve(const double* pars, RTCurve* rtc)
{
  double t1 = pars[0];
  double t0 = pars[1];
  rtc->SetRWidth(pars[2]);
  rtc->SetPoint(0, m_r_max, t1);
  for (int i_pt = 1; i_pt < m_n_pt - 1; i_pt++) {
    double t = t1 + (t0 - t1) * (i_pt + 1) / m_n_pt;
    rtc->SetPoint(i_pt, pars[i_pt+2], t);
  }
  rtc->SetPoint(m_n_pt-1, 0, t0);
}
