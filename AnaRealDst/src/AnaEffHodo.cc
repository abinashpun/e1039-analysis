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
#include "AnaEffHodo.h"
using namespace std;

AnaEffHodo::AnaEffHodo()
{
  ;
}

int AnaEffHodo::Init(PHCompositeNode* topNode)
{
  n_evt_all = n_evt_trig = n_evt_nhit = 0;

  ostringstream oss;
  ofs.open("output.txt");
  f_out = new TFile("output.root", "RECREATE");

  h1_eff_all = new TH1D("h1_eff_all", "", 1, 0, 1);
  h1_eff_ok  = new TH1D("h1_eff_ok" , "", 1, 0, 1);

  h1_nhit = new TH1D("h1_nhit", ";N of hits/plane/event;Hit count", 10, -0.5, 9.5);
  h1_ele  = new TH1D("h1_ele", ";Element ID;Hit count", 23, 0.5, 23.5);
  
  const double DT = 20/9.0; // 4/9 ns per single count of Taiwan TDC
  const int NT = 1000;
  const double T0 = 0.5*DT;
  const double T1 = (NT+0.5)*DT;
  h1_time = new TH1D("h1_time", ";tdcTime;Hit count", NT, T0, T1);

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaEffHodo::InitRun(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaEffHodo::process_event(PHCompositeNode* topNode)
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

  shared_ptr<SQHitVector> hv_h1t(UtilSQHit::FindFirstHits(hit_vec, "H1T"));
  shared_ptr<SQHitVector> hv_h1b(UtilSQHit::FindFirstHits(hit_vec, "H1B"));
  shared_ptr<SQHitVector> hv_h2t(UtilSQHit::FindFirstHits(hit_vec, "H2T"));
  shared_ptr<SQHitVector> hv_h2b(UtilSQHit::FindFirstHits(hit_vec, "H2B"));
  shared_ptr<SQHitVector> hv_h3t(UtilSQHit::FindFirstHits(hit_vec, "H3T"));
  shared_ptr<SQHitVector> hv_h3b(UtilSQHit::FindFirstHits(hit_vec, "H3B"));
  shared_ptr<SQHitVector> hv_h4t(UtilSQHit::FindFirstHits(hit_vec, "H4T"));
  shared_ptr<SQHitVector> hv_h4b(UtilSQHit::FindFirstHits(hit_vec, "H4B"));
  if (hv_h1t->size() + hv_h1b->size() != 1 ||
      hv_h2t->size() + hv_h2b->size() != 1 ||
      hv_h3t->size() + hv_h3b->size() != 1 ||
      hv_h4t->size() + hv_h4b->size() != 1 
    ) return Fun4AllReturnCodes::EVENT_OK;

  n_evt_nhit++;

  ///
  /// Analysis of detector hits
  ///
  shared_ptr<SQHitVector> hv_h2l(UtilSQHit::FindFirstHits(hit_vec, "H2L"));
  shared_ptr<SQHitVector> hv_h2r(UtilSQHit::FindFirstHits(hit_vec, "H2R"));

  int nhit = hv_h2l->size() + hv_h2r->size();
  h1_nhit->Fill(nhit);
  for (SQHitVector::ConstIter it = hv_h2l->begin(); it != hv_h2l->end(); it++) {
    h1_ele ->Fill((*it)->get_element_id());
    h1_time->Fill((*it)->get_tdc_time  ());
  }
  for (SQHitVector::ConstIter it = hv_h2r->begin(); it != hv_h2r->end(); it++) {
    h1_ele ->Fill((*it)->get_element_id());
    h1_time->Fill((*it)->get_tdc_time  ());
  }

  bool is_eff = nhit > 0; // very simple judgment for now
  if (is_eff) h1_eff_ok->Fill(0);
  h1_eff_all->Fill(0);

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaEffHodo::End(PHCompositeNode* topNode)
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

  h1_nhit->Draw();
  c1->SaveAs("h1_nhit.png");
  
  h1_ele->Draw();
  c1->SaveAs("h1_ele.png");
  
  h1_time->Draw();
  c1->SaveAs("h1_time.png");

  delete c1;

  f_out->cd();
  f_out->Write();
  f_out->Close();
  ofs.close();
  
  return Fun4AllReturnCodes::EVENT_OK;
}
