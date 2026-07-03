#include <iostream>
#include <iomanip>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include "RoadMap.h"
#include "RoadList.h"
#include "UtilRoad.h"
#include "ReAnaSignal.h"
using namespace std;

ReAnaSignal::ReAnaSignal(const std::string label)
  : AnaSignal(label)
{
  ;
}

ReAnaSignal::~ReAnaSignal()
{
  ;
}

void ReAnaSignal::Init()
{
  ostringstream oss;
  oss << m_dir_out << "/output.root";
  m_file_out = new TFile(oss.str().c_str(), "RECREATE");
  m_h1_mass_tb   = new TH1D("h1_mass_tb"  , "", 18, 1, 10);
  m_h1_mass_trig = new TH1D("h1_mass_trig", "", 18, 1, 10);
  m_h1_mass_tb  ->Sumw2();
  m_h1_mass_trig->Sumw2();
  m_h1_xF_tb   = new TH1D("h1_xF_tb"  , "", 20, -1, 1);
  m_h1_xF_trig = new TH1D("h1_xF_trig", "", 20, -1, 1);
  m_h1_xF_tb  ->Sumw2();
  m_h1_xF_trig->Sumw2();
  m_h1_x1_tb   = new TH1D("h1_x1_tb"  , "", 20, 0, 1);
  m_h1_x1_trig = new TH1D("h1_x1_trig", "", 20, 0, 1);
  m_h1_x1_tb  ->Sumw2();
  m_h1_x1_trig->Sumw2();
  m_h1_x2_tb   = new TH1D("h1_x2_tb"  , "", 20, 0, 1);
  m_h1_x2_trig = new TH1D("h1_x2_trig", "", 20, 0, 1);
  m_h1_x2_tb  ->Sumw2();
  m_h1_x2_trig->Sumw2();

  m_h1_mom_tb     = new TH1D("h1_mom_tb"  , "", 20, 20, 120);
  m_h1_mom_trig   = new TH1D("h1_mom_trig", "", 20, 20, 120);
  m_h1_mom_tb    ->Sumw2();
  m_h1_mom_trig  ->Sumw2();
  m_h1_phi_tb     = new TH1D("h1_phi_tb"  , "", 20, -TMath::Pi(), TMath::Pi());
  m_h1_phi_trig   = new TH1D("h1_phi_trig", "", 20, -TMath::Pi(), TMath::Pi());
  m_h1_phi_tb    ->Sumw2();
  m_h1_phi_trig  ->Sumw2();
  m_h1_theta_tb   = new TH1D("h1_theta_tb"  , "", 20, 0, 0.1);
  m_h1_theta_trig = new TH1D("h1_theta_trig", "", 20, 0, 0.1);
  m_h1_theta_tb  ->Sumw2();
  m_h1_theta_trig->Sumw2();

  m_h2_xF_tb   = new TH2D("h2_xF_tb"  , "", 20, -1, 1, 6, 2.5, 8.5);
  m_h2_xF_trig = new TH2D("h2_xF_trig", "", 20, -1, 1, 6, 2.5, 8.5);
  m_h2_xF_tb  ->Sumw2();
  m_h2_xF_trig->Sumw2();
  m_h2_x1_tb   = new TH2D("h2_x1_tb"  , "", 20, 0, 1, 6, 2.5, 8.5);
  m_h2_x1_trig = new TH2D("h2_x1_trig", "", 20, 0, 1, 6, 2.5, 8.5);
  m_h2_x1_tb  ->Sumw2();
  m_h2_x1_trig->Sumw2();
  m_h2_x2_tb   = new TH2D("h2_x2_tb"  , "", 20, 0, 1, 6, 2.5, 8.5);
  m_h2_x2_trig = new TH2D("h2_x2_trig", "", 20, 0, 1, 6, 2.5, 8.5);
  m_h2_x2_tb  ->Sumw2();
  m_h2_x2_trig->Sumw2();

  m_h2_mom_tb     = new TH2D("h2_mom_tb"  , "", 20, 20, 120, 6, 2.5, 8.5);
  m_h2_mom_trig   = new TH2D("h2_mom_trig", "", 20, 20, 120, 6, 2.5, 8.5);
  m_h2_mom_tb    ->Sumw2();
  m_h2_mom_trig  ->Sumw2();
  m_h2_phi_tb     = new TH2D("h2_phi_tb"  , "", 20, -TMath::Pi(), TMath::Pi(), 6, 2.5, 8.5);
  m_h2_phi_trig   = new TH2D("h2_phi_trig", "", 20, -TMath::Pi(), TMath::Pi(), 6, 2.5, 8.5);
  m_h2_phi_tb    ->Sumw2();
  m_h2_phi_trig  ->Sumw2();
  m_h2_theta_tb   = new TH2D("h2_theta_tb"  , "", 20, 0, 0.1, 6, 2.5, 8.5);
  m_h2_theta_trig = new TH2D("h2_theta_trig", "", 20, 0, 0.1, 6, 2.5, 8.5);
  m_h2_theta_tb  ->Sumw2();
  m_h2_theta_trig->Sumw2();
}

void ReAnaSignal::End()
{
  m_file_out->cd();
  m_file_out->Write();
  m_file_out->Close();
}

void ReAnaSignal::ProcessOneEvent()
{
  double weight = m_sig_data->weight;
  double mass   = m_sig_data->mass;
  //double pT     = m_sig_data->pT;
  double xF     = m_sig_data->xF;
  double x1     = m_sig_data->x1;
  double x2     = m_sig_data->x2;
  double mom    = m_sig_data->mom;
  double phi    = m_sig_data->phi;
  double theta  = m_sig_data->theta;
  int road_pos  = m_sig_data->road_pos;
  int road_neg  = m_sig_data->road_neg;

  if (road_pos * road_neg >= 0) return;

  m_h1_mass_tb->Fill(mass, weight);
  if (m_mass_lo <= mass && mass <= m_mass_hi) {
    m_h1_xF_tb   ->Fill(xF   , weight);
    m_h1_x1_tb   ->Fill(x1   , weight);
    m_h1_x2_tb   ->Fill(x2   , weight);
    m_h1_mom_tb  ->Fill(mom  , weight);
    m_h1_phi_tb  ->Fill(phi  , weight);
    m_h1_theta_tb->Fill(theta, weight);
  }
  m_h2_xF_tb   ->Fill(xF   , mass, weight);
  m_h2_x1_tb   ->Fill(x1   , mass, weight);
  m_h2_x2_tb   ->Fill(x2   , mass, weight);
  m_h2_mom_tb  ->Fill(mom  , mass, weight);
  m_h2_phi_tb  ->Fill(phi  , mass, weight);
  m_h2_theta_tb->Fill(theta, mass, weight);

  RoadMap* map_pos = road_pos > 0  ?  &m_road_map_pos_top  :  &m_road_map_pos_bot;
  RoadMap* map_neg = road_neg > 0  ?  &m_road_map_neg_top  :  &m_road_map_neg_bot;
  if (map_pos->Find(road_pos) && map_neg->Find(road_neg)) {
    m_h1_mass_trig->Fill(mass, weight);
    if (m_mass_lo <= mass && mass <= m_mass_hi) {
      m_h1_xF_trig   ->Fill(xF   , weight);
      m_h1_x1_trig   ->Fill(x1   , weight);
      m_h1_x2_trig   ->Fill(x2   , weight);
      m_h1_mom_trig  ->Fill(mom  , weight);
      m_h1_phi_trig  ->Fill(phi  , weight);
      m_h1_theta_trig->Fill(theta, weight);
    }
    m_h2_xF_trig   ->Fill(xF   , mass, weight);
    m_h2_x1_trig   ->Fill(x1   , mass, weight);
    m_h2_x2_trig   ->Fill(x2   , mass, weight);
    m_h2_mom_trig  ->Fill(mom  , mass, weight);
    m_h2_phi_trig  ->Fill(phi  , mass, weight);
    m_h2_theta_trig->Fill(theta, mass, weight);
  }
}

void ReAnaSignal::Analyze()
{
  DrawOneVar("mass" , m_h1_mass_tb , m_h1_mass_trig );
  DrawOneVar("xF"   , m_h1_xF_tb   , m_h1_xF_trig   );
  DrawOneVar("x1"   , m_h1_x1_tb   , m_h1_x1_trig   );
  DrawOneVar("x2"   , m_h1_x2_tb   , m_h1_x2_trig   );
  DrawOneVar("mom"  , m_h1_mom_tb  , m_h1_mom_trig  );
  DrawOneVar("phi"  , m_h1_phi_tb  , m_h1_phi_trig  );
  DrawOneVar("theta", m_h1_theta_tb, m_h1_theta_trig);

  DrawOneVar2D("xF"   , m_h2_xF_tb   , m_h2_xF_trig   );
  DrawOneVar2D("x1"   , m_h2_x1_tb   , m_h2_x1_trig   );
  DrawOneVar2D("x2"   , m_h2_x2_tb   , m_h2_x2_trig   );
  DrawOneVar2D("mom"  , m_h2_mom_tb  , m_h2_mom_trig  );
  DrawOneVar2D("phi"  , m_h2_phi_tb  , m_h2_phi_trig  );
  DrawOneVar2D("theta", m_h2_theta_tb, m_h2_theta_trig);
}

void ReAnaSignal::DrawOneVar(const char* name, TH1* h1_tb, TH1* h1_trig)
{
  TCanvas* c1 = new TCanvas("c1", "");
  c1->SetGrid();
  gStyle->SetOptStat(0000);

  h1_tb  ->SetLineColor(kBlack);
  h1_trig->SetLineColor(kRed  );

  ostringstream oss;
  oss << "GMC;" << name << ";Weighted yield";
  THStack hs("hs", oss.str().c_str());
  hs.Add(h1_tb  , "E1");
  hs.Add(h1_trig, "E1");
  hs.Draw("nostack");

  TLegend leg (0.75, 0.80, 0.99, 0.99);
  leg.AddEntry(h1_tb  , "All T+B or B+T", "l");
  leg.AddEntry(h1_trig, "Triggered"     , "l");
  leg.SetTextFont(22);
  leg.SetBorderSize(1);
  leg.SetFillColor(0);
  leg.Draw();

  oss.str("");
  oss << m_dir_out << "/h1_" << name << ".png";
  c1->SaveAs(oss.str().c_str());

  c1->SetLogy(true);
  
  oss.str("");
  oss << m_dir_out << "/h1_" << name << "_log.png";
  c1->SaveAs(oss.str().c_str());
  
  c1->SetLogy(false);

  TH1* h1_ratio = (TH1*)h1_trig->Clone("h1_ratio");
  oss.str("");
  oss << "GMC;" << name << ";Trigger acceptance";
  h1_ratio->SetTitle(oss.str().c_str());
  for (int ib = 1; ib <= h1_ratio->GetNbinsX(); ib++) {
    double n_tb   = h1_tb  ->GetBinContent(ib);
    double e_tb   = h1_tb  ->GetBinError  (ib);
    double n_trig = h1_trig->GetBinContent(ib);
    double e_trig = h1_trig->GetBinError  (ib);
    double frac   = 0;
    double err    = 0;
    if (n_tb > 0) {
      double n_nont = n_tb - n_trig; // non-triggered
      double e_nont = sqrt( pow(e_tb,2) - pow(e_trig,2) );
      frac = n_trig / n_tb;
      err  = sqrt( pow(n_trig*e_nont, 2) + pow(e_trig*n_nont, 2) ) / pow(n_tb, 2);
    }
    h1_ratio->SetBinContent(ib, frac);
    h1_ratio->SetBinError  (ib, err );
  }
  h1_ratio->Draw("E1");
  h1_ratio->GetYaxis()->SetRangeUser(0, 1);
  oss.str("");
  oss << m_dir_out << "/h1_" << name << "_ratio.png";
  c1->SaveAs(oss.str().c_str());

  delete h1_ratio;
  delete c1;
}

void ReAnaSignal::DrawOneVar2D(const char* name, TH2* h2_tb, TH2* h2_trig)
{
  TCanvas* c1 = new TCanvas("c1", "");
  c1->SetGrid();
  gStyle->SetOptStat(0000);

  ostringstream oss;
  oss << "GMC;" << name << ";Trigger acceptance";
  THStack hs_ratio("hs", oss.str().c_str());
  TLegend leg_ratio(0.89, 0.70, 0.99, 0.99);
  TH1* h1_ratio[99];

  for (int iy = 1; iy <= h2_tb->GetNbinsY(); iy++) {
    double m_lo = h2_tb->GetYaxis()->GetBinLowEdge(iy  );
    double m_hi = h2_tb->GetYaxis()->GetBinLowEdge(iy+1);
    TH1* h1_tb   = h2_tb  ->ProjectionX("h1_tb"  , iy, iy);
    TH1* h1_trig = h2_trig->ProjectionX("h1_trig", iy, iy);
    h1_tb  ->SetLineColor(kBlack);
    h1_trig->SetLineColor(kRed  );

    oss.str("");
    oss << "GMC: M=" << m_lo << "-" << m_hi << " GeV;" << name << ";Weighted yield";
    THStack hs("hs", oss.str().c_str());
    hs.Add(h1_tb  , "E1");
    hs.Add(h1_trig, "E1");
    hs.Draw("nostack");

    TLegend leg(0.75, 0.80, 0.99, 0.99);
    leg.AddEntry(h1_tb  , "All T+B or B+T", "l");
    leg.AddEntry(h1_trig, "Triggered"     , "l");
    leg.SetTextFont(22);
    leg.SetBorderSize(1);
    leg.SetFillColor(0);
    leg.Draw();

    oss.str("");
    oss << m_dir_out << "/h1_" << name << "_m" << iy << ".png";
    c1->SaveAs(oss.str().c_str());
    
    c1->SetLogy(true);
    
    oss.str("");
    oss << m_dir_out << "/h1_" << name << "_m" << iy << "_log.png";
    c1->SaveAs(oss.str().c_str());
    
    c1->SetLogy(false);

    oss.str("");
    oss << "h1_ratio_m" << iy;
    h1_ratio[iy] = (TH1*)h1_trig->Clone(oss.str().c_str());
    oss.str("");
    oss << "GMC: M=" << m_lo << "-" << m_hi << " GeV;" << name << ";Trigger acceptance";
    h1_ratio[iy]->SetTitle(oss.str().c_str());
    for (int ib = 1; ib <= h1_ratio[iy]->GetNbinsX(); ib++) {
      double n_tb   = h1_tb  ->GetBinContent(ib);
      double e_tb   = h1_tb  ->GetBinError  (ib);
      double n_trig = h1_trig->GetBinContent(ib);
      double e_trig = h1_trig->GetBinError  (ib);
      double frac   = 0;
      double err    = 0;
      if (n_tb > 0) {
        double n_nont = n_tb - n_trig; // non-triggered
        double e_nont = sqrt( pow(e_tb,2) - pow(e_trig,2) );
        frac = n_trig / n_tb;
        err  = sqrt( pow(n_trig*e_nont, 2) + pow(e_trig*n_nont, 2) ) / pow(n_tb, 2);
      }
      h1_ratio[iy]->SetBinContent(ib, frac);
      h1_ratio[iy]->SetBinError  (ib, err );
    }
    //h1_ratio[iy]->Draw("E1");
    //h1_ratio[iy]->GetYaxis()->SetRangeUser(0, 1);
    //oss.str("");
    //oss << m_dir_out << "/h1_" << name << "_ratio_m" << iy << ".png";
    //c1->SaveAs(oss.str().c_str());

    h1_ratio[iy]->SetLineColor  (iy+1);
    h1_ratio[iy]->SetMarkerColor(iy+1);
    h1_ratio[iy]->SetMarkerStyle(21);

    oss.str("");
    oss << m_lo << "-" << m_hi;
    leg_ratio.AddEntry(h1_ratio[iy], oss.str().c_str(), "l");
    hs_ratio .Add     (h1_ratio[iy], "E1");

    delete h1_tb;
    delete h1_trig;
  }

  hs_ratio.SetMinimum(0.0);
  hs_ratio.SetMaximum(1.0);
  hs_ratio .Draw("nostack");
  leg_ratio.SetHeader("Mass");
  leg_ratio.SetTextFont(22);
  leg_ratio.SetBorderSize(1);
  leg_ratio.SetFillColor(0);
  leg_ratio.Draw();
  oss.str("");
  oss << m_dir_out << "/h1_" << name << "_ratio.png";
  c1->SaveAs(oss.str().c_str());

  for (int iy = 1; iy <= h2_tb->GetNbinsY(); iy++) delete h1_ratio[iy];
  delete c1;
}

//void ReAnaSignal::DrawOneVar(TTree* tree, const char* name, const char* var, const int N, const double X0, const double X1)
//{
//  static TCanvas* c1 = 0;
//  if (! c1) {
//    c1 = new TCanvas("c1", "");
//    c1->SetGrid();
//    c1->SetLogy(true);
//    gStyle->SetOptStat(0000);
//  }
//
//  TH1* h1_all = new TH1D("h1_all", "", N, X0, X1);
//  TH1* h1_acc = new TH1D("h1_acc", "", N, X0, X1);
//  TH1* h1_tb  = new TH1D("h1_tb" , "", N, X0, X1);
//  h1_all->Sumw2();
//  h1_acc->Sumw2();
//  h1_tb ->Sumw2();
//  tree->Project("h1_all", var, "weight");
//  tree->Project("h1_acc", var, "weight * (road_pos*road_neg != 0)");
//  tree->Project("h1_tb" , var, "weight * (road_pos*road_neg <  0)");
//
//  h1_all->SetLineColor(kBlack);
//  h1_acc->SetLineColor(kBlue);
//  h1_tb ->SetLineColor(kRed);
//
//  ostringstream oss;
//  oss << "GMC;" << name << ";Weighted yield";
//  THStack hs("hs", oss.str().c_str());
//  hs.Add(h1_all, "E1");
//  hs.Add(h1_acc, "E1");
//  hs.Add(h1_tb , "E1");
//  hs.Draw("nostack");
//
//  TLegend leg (0.75, 0.75, 0.99, 0.99);
//  leg.AddEntry(h1_all, "All generated", "l");
//  leg.AddEntry(h1_acc, "In hodo acceptance", "l");
//  leg.AddEntry(h1_tb , "T+B or B+T", "l");
//  leg.SetTextFont(22);
//  leg.SetBorderSize(1);
//  leg.SetFillColor(0);
//  leg.Draw();
//
//  oss.str("");
//  oss << m_dir_out << "/h1_" << name << ".png";
//  c1->SaveAs(oss.str().c_str());
//
//  delete h1_all;
//  delete h1_acc;
//  delete h1_tb ;
//}

