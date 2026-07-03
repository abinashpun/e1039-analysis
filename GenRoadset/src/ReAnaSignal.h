#ifndef _RE_ANA_SIGNAL__H_
#define _RE_ANA_SIGNAL__H_
#include "AnaSignal.h"
class TFile;
class TH1;
class TH2;

/// Class to analyze the GMC data with a roadset applied.
class ReAnaSignal : public AnaSignal {
 protected:
  TFile* m_file_out;
  TH1* m_h1_mass_tb;
  TH1* m_h1_mass_trig;
  TH1* m_h1_xF_tb;
  TH1* m_h1_xF_trig;
  TH1* m_h1_x1_tb;
  TH1* m_h1_x1_trig;
  TH1* m_h1_x2_tb;
  TH1* m_h1_x2_trig;
  TH1* m_h1_mom_tb;
  TH1* m_h1_mom_trig;
  TH1* m_h1_phi_tb;
  TH1* m_h1_phi_trig;
  TH1* m_h1_theta_tb;
  TH1* m_h1_theta_trig;

  TH2* m_h2_xF_tb;
  TH2* m_h2_xF_trig;
  TH2* m_h2_x1_tb;
  TH2* m_h2_x1_trig;
  TH2* m_h2_x2_tb;
  TH2* m_h2_x2_trig;
  TH2* m_h2_mom_tb;
  TH2* m_h2_mom_trig;
  TH2* m_h2_phi_tb;
  TH2* m_h2_phi_trig;
  TH2* m_h2_theta_tb;
  TH2* m_h2_theta_trig;

 public:
  ReAnaSignal(const std::string label="re_ana_signal");
  virtual ~ReAnaSignal();

  virtual void Init();
  virtual void End();
  virtual void Analyze();

 protected:
  virtual void ProcessOneEvent();

  void DrawOneVar(const char* name, TH1* h1_tb, TH1* h1_trig);
  //void DrawOneVar(TTree* tree, const char* name, const char* var, const int N, const double X0, const double X1);
  void DrawOneVar2D(const char* name, TH2* h2_tb, TH2* h2_trig);
};

#endif // _RE_ANA_SIGNAL__H_
