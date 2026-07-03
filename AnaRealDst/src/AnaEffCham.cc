#include <iomanip>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <THStack.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <interface_main/SQRun.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQHitVector.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <geom_svc/GeomSvc.h>
#include <UtilAna/UtilSQHit.h>
#include "AnaEffCham.h"
using namespace std;

AnaEffCham::AnaEffCham(const ChamType_t type) : m_type(type)
{
  switch (m_type) {
  case D0 :  m_pl0 =  1;  break;
  case D1 :  m_pl0 =  7;  break;
  case D2 :  m_pl0 = 13;  break;
  case D3p:  m_pl0 = 19;  break;
  case D3m:  m_pl0 = 25;  break;
  }
}

int AnaEffCham::Init(PHCompositeNode* topNode)
{
  n_evt_all = n_evt_trig = n_evt_nhit = 0;

  ostringstream oss;
  GeomSvc* geom = GeomSvc::instance();
  ofs.open("output.txt");
  f_out = new TFile("output.root", "RECREATE");

  h1_eff_all = new TH1D("h1_eff_all", "", N_PL, 0, N_PL);
  h1_eff_ok  = new TH1D("h1_eff_ok" , "", N_PL, 0, N_PL);

  for (int i_pl = 0; i_pl < N_PL; i_pl++) {
    int det_id = m_pl0 + i_pl;
    string name = geom->getDetectorName(det_id);
    int   n_ele = geom->getPlaneNElements(det_id); 
    list_pl_id   .push_back(det_id);
    list_pl_name .push_back(name);
    list_pl_n_ele.push_back(n_ele);
    cout << "  " << setw(6) << name << " " << det_id << " " << n_ele << endl;

    h1_eff_all->GetXaxis()->SetBinLabel(i_pl+1, name.c_str());
    h1_eff_ok ->GetXaxis()->SetBinLabel(i_pl+1, name.c_str());

    oss.str("");
    oss << "h1_nhit_" << det_id << "_" << name;
    h1_nhit[i_pl] = new TH1D(oss.str().c_str(), "", 10, -0.5, 9.5);
    oss.str("");
    oss << name << ";N of hits/plane/event;Hit count";
    h1_nhit[i_pl]->SetTitle(oss.str().c_str());

    oss.str("");
    oss << "h1_ele_" << det_id << "_" << name;
    h1_ele[i_pl] = new TH1D(oss.str().c_str(), "", n_ele, 0.5, n_ele+0.5);
    oss.str("");
    oss << name << ";Element ID;Hit count";
    h1_ele[i_pl]->SetTitle(oss.str().c_str());

    const double DT = 20/9.0; // 4/9 ns per single count of Taiwan TDC
    const int NT = 1000;
    const double T0 = 0.5*DT;
    const double T1 = (NT+0.5)*DT;
    oss.str("");
    oss << "h1_time_" << det_id << "_" << name;
    h1_time[i_pl] = new TH1D(oss.str().c_str(), "", NT, T0, T1);
    oss.str("");
    oss << name << ";tdcTime;Hit count";
    h1_time[i_pl]->SetTitle(oss.str().c_str());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaEffCham::InitRun(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaEffCham::process_event(PHCompositeNode* topNode)
{
  SQEvent* event       = findNode::getClass<SQEvent    >(topNode, "SQEvent");
  SQHitVector* hit_vec = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if (!event || !hit_vec) return Fun4AllReturnCodes::ABORTEVENT;
  //int run_id   = event->get_run_id();
  //int spill_id = event->get_spill_id();
  //int event_id = event->get_event_id();
  n_evt_all++;

  ///
  /// Event selection
  ///
  if (! event->get_trigger(SQEvent::NIM1)) { // depends on the setting of run(s) you analyze
    return Fun4AllReturnCodes::EVENT_OK;
  }
  n_evt_trig++;

  //shared_ptr<SQHitVector> hv_h1(UtilSQHit::FindFirstHits(hit_vec, "H1T"));
  //shared_ptr<SQHitVector> hv_h2(UtilSQHit::FindFirstHits(hit_vec, "H2T"));
  shared_ptr<SQHitVector> hv_h3(UtilSQHit::FindFirstHits(hit_vec, "H3T"));
  shared_ptr<SQHitVector> hv_h4(UtilSQHit::FindFirstHits(hit_vec, "H4T"));
  if (//hv_h1->size() != 1 ||
      //hv_h2->size() != 1 ||
      hv_h3->size() != 1 ||
      hv_h4->size() != 1   ) return Fun4AllReturnCodes::EVENT_OK;
  n_evt_nhit++;

  ///
  /// Analysis of detector hits
  ///
  for (int i_pl = 0; i_pl < N_PL; i_pl++) {
    shared_ptr<SQHitVector> hv(UtilSQHit::FindHits(hit_vec, list_pl_id[i_pl]));

    bool is_eff = hv->size() > 0; // very simple judgment for now
    if (is_eff) h1_eff_ok->Fill(i_pl);
    h1_eff_all->Fill(i_pl);

    h1_nhit[i_pl]->Fill(hv->size());
    for (SQHitVector::ConstIter it = hv->begin(); it != hv->end(); it++) {
      int    ele  = (*it)->get_element_id();
      double time = (*it)->get_tdc_time  ();
      h1_ele [i_pl]->Fill(ele );
      h1_time[i_pl]->Fill(time);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaEffCham::End(PHCompositeNode* topNode)
{
  ostringstream oss;

  ofs << "N of events:\n"
      << "  All        = " << n_evt_all << "\n"
      << "  Trigger OK = " << n_evt_trig << "\n"
      << "  n_hit OK   = " << n_evt_nhit << "\n"
      << endl;

  TCanvas* c1 = new TCanvas("c1", "");
  c1->SetGrid();

  h1_eff_all->SetLineColor(kRed);
  h1_eff_all->SetLineWidth(2);
  THStack hs("hs", ";;N of hits");
  hs.Add(h1_eff_all);
  hs.Add(h1_eff_ok );
  hs.Draw("nostack");
  c1->SaveAs("h1_eff.png");

  TEfficiency* teff_pl = new TEfficiency(*h1_eff_ok, *h1_eff_all);
  teff_pl->SetMarkerStyle(21);
  teff_pl->SetMarkerColor(kRed);
  teff_pl->SetLineColor  (kRed);
  teff_pl->Draw("APE1");
  teff_pl->SetTitle(";Plane;Plane efficiency");
  c1->SaveAs("teff_pl.png");

  for (int i_pl = 0; i_pl < N_PL; i_pl++) {
    h1_nhit[i_pl]->Draw();
    oss.str("");
    oss << h1_nhit[i_pl]->GetName() << ".png";
    c1->SaveAs(oss.str().c_str());

    h1_ele[i_pl]->Draw();
    oss.str("");
    oss << h1_ele[i_pl]->GetName() << ".png";
    c1->SaveAs(oss.str().c_str());

    h1_time[i_pl]->Draw();
    oss.str("");
    oss << h1_time[i_pl]->GetName() << ".png";
    c1->SaveAs(oss.str().c_str());
  }
  delete c1;

  f_out->cd();
  f_out->Write();
  f_out->Close();
  ofs.close();
  
  return Fun4AllReturnCodes::EVENT_OK;
}
