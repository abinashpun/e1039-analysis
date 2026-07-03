#include <iomanip>
#include <TSystem.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLegend.h>
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
#include "AnaTrigSignal.h"
using namespace std;

/// List of detectors that you want to analyze
const vector<string> AnaTrigSignal::m_list_det_name = { "BeforeInhMatrix", "AfterInhMatrix", "AfterInhNIM" };

AnaTrigSignal::AnaTrigSignal(const std::string& name)
  : SubsysReco(name)
  , m_file(0)
  , m_tree(0)
{
  ;
}

int AnaTrigSignal::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTrigSignal::InitRun(PHCompositeNode* topNode)
{
  ///
  /// Input
  ///
  m_evt     = findNode::getClass<SQEvent    >(topNode, "SQEvent");
  m_hit_vec = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if (!m_evt || !m_hit_vec) return Fun4AllReturnCodes::ABORTEVENT;

  ///
  /// Output
  ///
  gSystem->mkdir("result", true);

  m_file = new TFile("result/output.root", "RECREATE");
  m_tree = new TTree("tree", "Created by AnaTrigSignal");
  m_tree->Branch("det_name", &b_det_name, "det_name/C");
  m_tree->Branch("det"     , &b_det     ,      "det/I");
  m_tree->Branch("ele"     , &b_ele     ,      "ele/I");
  m_tree->Branch("time"    , &b_time    ,     "time/D");

  ostringstream oss;
  GeomSvc* geom = GeomSvc::instance();
  for (unsigned int i_det = 0; i_det < m_list_det_name.size(); i_det++) {
    string name = m_list_det_name[i_det];
    int id = geom->getDetectorID(name);
    if (id <= 0) {
      cerr << "!ERROR!  AnaTrigSignal::InitRun():  Invalid ID (" << id << ").  Probably the detector name that you specified in 'list_det_name' (" << name << ") is not valid.  Abort." << endl;
      exit(1);
    }
    m_list_det_id.push_back(id);

    const double DT = 12/9.0; // 4/9 ns per single count of Taiwan TDC
    int    NT = 1500;
    double T0 =    0.5*DT;
    double T1 = 1500.5*DT;
        
    oss.str("");
    oss << "h2_time_ele_" << name;
    m_h2_time_ele[i_det] = new TH2D(oss.str().c_str(), "", NT, T0, T1,  N_ELE, 0.5, N_ELE+0.5);
    oss.str("");
    oss << name << ";TDC time (ns);Element ID";
    m_h2_time_ele[i_det]->SetTitle(oss.str().c_str());

    oss.str("");
    oss << "h1_nhit_" << name;
    m_h1_nhit[i_det] = new TH1D(oss.str().c_str(), "", 10, -0.5, 9.5);
    oss.str("");
    oss << name << ";N of hits/plane/event;Hit count";
    m_h1_nhit[i_det]->SetTitle(oss.str().c_str());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTrigSignal::process_event(PHCompositeNode* topNode)
{
  //int spill_id = m_evt->get_spill_id();
  //int event_id = m_evt->get_event_id();

  ///
  /// Event selection
  ///
  //if (! m_evt->get_trigger(SQEvent::NIM2)) {
  //  return Fun4AllReturnCodes::EVENT_OK;
  //}

  ///
  /// Get & fill the hit info
  ///
  for (unsigned int i_det = 0; i_det < m_list_det_name.size(); i_det++) {
    strncpy(b_det_name, m_list_det_name[i_det].c_str(), sizeof(b_det_name));
    b_det = m_list_det_id[i_det];
    shared_ptr<SQHitVector> hv(UtilSQHit::FindHits(m_hit_vec, b_det));
    for (SQHitVector::ConstIter it = hv->begin(); it != hv->end(); it++) {
      b_ele  = (*it)->get_element_id();
      b_time = (*it)->get_tdc_time  ();
      m_tree->Fill();

      m_h2_time_ele[i_det]->Fill(b_time, b_ele);
    }
    m_h1_nhit[i_det]->Fill(hv->size());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTrigSignal::End(PHCompositeNode* topNode)
{
  ostringstream oss;
  TCanvas* c1 = new TCanvas("c1", "");
  c1->SetGrid();
  for (unsigned int i_det = 0; i_det < m_list_det_id.size(); i_det++) {
    string name = m_list_det_name[i_det];

    gStyle->SetOptStat(0);
    UtilHist::AutoSetRangeX(m_h2_time_ele[i_det]);
    m_h2_time_ele[i_det]->Draw("colz");
    oss.str("");
    oss << "result/" << m_h2_time_ele[i_det]->GetName() << ".png";
    c1->SaveAs(oss.str().c_str());

    TH1* h1_ele = m_h2_time_ele[i_det]->ProjectionY("h1_ele");
    h1_ele->Draw();
    oss.str("");
    oss << "result/h1_ele_" << name << ".png";
    c1->SaveAs(oss.str().c_str());
    delete h1_ele;

    oss.str("");
    oss << name << ";TDC time (ns);Hit count";
    THStack hs("hs", oss.str().c_str());
    TLegend leg(0.9, 0.5, 0.99, 0.9);
    for (int i_ele = 1; i_ele <= N_ELE; i_ele++) {
      oss.str("");
      oss << "h1_time_" << i_ele;
      TH1* h1_time = m_h2_time_ele[i_det]->ProjectionX(oss.str().c_str(), i_ele, i_ele);
      h1_time->SetLineColor(i_ele);
      hs.Add(h1_time);
      oss.str("");
      oss << "#" << i_ele;
      leg.AddEntry(h1_time, oss.str().c_str(), "l");
    }
    hs.Draw("nostack");
    leg.Draw();
    oss.str("");
    oss << "result/h1_time_" << name << ".png";
    c1->SaveAs(oss.str().c_str());

    m_h1_nhit[i_det]->Draw();
    oss.str("");
    oss << "result/" << m_h1_nhit[i_det]->GetName() << ".png";
    c1->SaveAs(oss.str().c_str());
  }
  delete c1;

  m_file->cd();
  m_file->Write();
  m_file->Close();
  
  return Fun4AllReturnCodes::EVENT_OK;
}
