#include <iomanip>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
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
#include <UtilAna/UtilHist.h>
#include "AnaChamHit.h"
using namespace std;

AnaChamHit::AnaChamHit(const std::string& name)
  : SubsysReco  (name)
  , m_sq_evt    (0)
  , m_sq_hit_vec(0)
  , m_file_name ("output.root")
  , m_file      (0)
  , m_tree      (0)
{
  ;
}

AnaChamHit::~AnaChamHit()
{
  ;
}

int AnaChamHit::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaChamHit::InitRun(PHCompositeNode* topNode)
{
  m_sq_evt     = findNode::getClass<SQEvent    >(topNode, "SQEvent");
  m_sq_hit_vec = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if (!m_sq_evt || !m_sq_hit_vec) return Fun4AllReturnCodes::ABORTEVENT;

  m_file = new TFile(m_file_name.c_str(), "RECREATE");
  m_tree = new TTree("tree", "Created by AnaChamHit");
  m_tree->Branch("event"   , &m_evt);
  m_tree->Branch("hit_list", &m_hit_list);

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaChamHit::process_event(PHCompositeNode* topNode)
{
  if (! m_sq_evt->get_trigger(SQEvent::MATRIX1)) {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_evt.run_id    = m_sq_evt->get_run_id();
  m_evt.spill_id  = m_sq_evt->get_spill_id();
  m_evt.event_id  = m_sq_evt->get_event_id();
  m_evt.fpga_bits = (m_sq_evt->get_trigger() >> SQEvent::MATRIX1) & 0x1f;
  m_evt.nim_bits  = (m_sq_evt->get_trigger() >> SQEvent::NIM1   ) & 0x1f;
  m_evt.D1 = m_evt.D2 = m_evt.D3p = m_evt.D3m = 0;
  
  GeomSvc* geom = GeomSvc::instance();
  m_hit_list.clear();
  for (SQHitVector::Iter it = m_sq_hit_vec->begin(); it != m_sq_hit_vec->end(); it++) {
    SQHit* hit = *it;
    int det_id = hit->get_detector_id();
    if (! geom->isChamber(det_id)) continue;

    if      ( 0 < det_id && det_id <=  6) m_evt.D1++;
    else if (12 < det_id && det_id <= 18) m_evt.D2++;
    else if (18 < det_id && det_id <= 24) m_evt.D3p++;
    else if (24 < det_id && det_id <= 30) m_evt.D3m++;

    HitData hd;
    hd.det_id     = det_id;
    hd.ele_id     = hit->get_element_id();
    hd.tdc_time   = hit->get_tdc_time();
    hd.drift_dist = hit->get_drift_distance();
    m_hit_list.push_back(hd);
  }

  m_tree->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaChamHit::End(PHCompositeNode* topNode)
{
  m_file->cd();
  m_file->Write();
  m_file->Close();  
  return Fun4AllReturnCodes::EVENT_OK;
}

void AnaChamHit::AnalyzeTree(TChain* tree)
{
  cout << "N of trees = " << tree->GetNtrees() << endl;

  TH1* h1_ele [31];
  TH1* h1_time[31];
  TH1* h1_dist[31];
  GeomSvc* geom = GeomSvc::instance();
  ostringstream oss;
  for (int pl = 1; pl <= 30; pl++) {
    if (6 < pl && pl <= 12) continue;
    string name  = geom->getDetectorName  (pl);
    int   n_ele  = geom->getPlaneNElements(pl);
    double space = geom->getPlaneSpacing  (pl);
    oss.str("");
    oss << "h1_ele_" << pl;
    h1_ele[pl] = new TH1D(oss.str().c_str(), "", n_ele, 0.5, n_ele+0.5);
    oss.str("");
    oss << setw(2) << pl << ":" << name << ";Element ID;Hit count";
    h1_ele[pl]->SetTitle(oss.str().c_str());

    const double DT = 12/9.0; // 4/9 ns per single count of Taiwan TDC
    oss.str("");
    oss << "h1_time_" << pl;
    h1_time[pl] = new TH1D(oss.str().c_str(), "", 1000, 400.5*DT, 1400.5*DT);
    oss.str("");
    oss << setw(2) << pl << ":" << name << ";tdcTime (ns);Hit count";
    h1_time[pl]->SetTitle(oss.str().c_str());

    oss.str("");
    oss << "h1_dist_" << pl;
    h1_dist[pl] = new TH1D(oss.str().c_str(), "", 110, 0.0, 0.55*space);
    oss.str("");
    oss << setw(2) << pl << ":" << name << ";Drift distance (cm);Hit count";
    h1_dist[pl]->SetTitle(oss.str().c_str());
  }
  
  EventData* evt = 0;
  HitList* hit_list = 0;
  tree->SetBranchAddress("event"   , &evt);
  tree->SetBranchAddress("hit_list", &hit_list);

  int n_ent = tree->GetEntries();
  cout << "N of entries = " << n_ent << endl;
  for (int i_ent = 0; i_ent < n_ent; i_ent++) {
    if ((i_ent+1) % (n_ent/10) == 0) cout << "  " << 10*(i_ent+1)/(n_ent/10) << "%" << flush;
    tree->GetEntry(i_ent);
    for (auto it = hit_list->begin(); it != hit_list->end(); it++) {
      HitData* hd = &(*it);
      short det_id = hd->det_id;
      if (6 < det_id && det_id <= 12) continue;
      short ele_id = hd->ele_id;
      double time = hd->tdc_time;
      double dist = hd->drift_dist;
      h1_ele [det_id]->Fill(ele_id);
      h1_time[det_id]->Fill(time);
      h1_dist[det_id]->Fill(dist);
    }
  }
  
  gSystem->mkdir("result", true);
  TCanvas* c1 = new TCanvas("c1", "");
  c1->SetGrid();
  //c1->SetLogy(true);

  oss << setfill('0');
  for (int pl = 1; pl <= 30; pl++) {
    if (6 < pl && pl <= 12) continue;
    h1_ele[pl]->Draw();
    oss.str("");
    oss << "result/h1_ele_" << setw(2) << pl << ".png";
    c1->SaveAs(oss.str().c_str());

    h1_time[pl]->Draw();
    UtilHist::AutoSetRange(h1_time[pl]);
    oss.str("");
    oss << "result/h1_time_" << setw(2) << pl << ".png";
    c1->SaveAs(oss.str().c_str());

    h1_dist[pl]->Draw();
    oss.str("");
    oss << "result/h1_dist_" << setw(2) << pl << ".png";
    c1->SaveAs(oss.str().c_str());
  }

  delete c1;
}

