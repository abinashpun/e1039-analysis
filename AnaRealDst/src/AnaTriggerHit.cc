#include <iomanip>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
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
#include "AnaTriggerHit.h"
using namespace std;

/// List of detectors that you want to analyze
const vector<string> AnaTriggerHit::m_list_det_name = {
  "H1T", "H1B", "H2T", "H2B", "H3T", "H3B", "H4T", "H4B"
};

AnaTriggerHit::AnaTriggerHit(const std::string& name)
  : SubsysReco(name)
  , m_intime(false)
  , m_level(1)
  , m_file(0)
  , m_tree(0)
{
  ;
}

int AnaTriggerHit::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTriggerHit::InitRun(PHCompositeNode* topNode)
{
  ///
  /// Input
  ///
  m_evt     = findNode::getClass<SQEvent    >(topNode, "SQEvent");
  m_hit_vec = findNode::getClass<SQHitVector>(topNode, "SQTriggerHitVector");
  if (!m_evt || !m_hit_vec) return Fun4AllReturnCodes::ABORTEVENT;

  ///
  /// Output
  ///
  gSystem->mkdir("result", true);

  m_file = new TFile("result/output.root", "RECREATE");
  m_tree  = new TTree("tree", "Created by AnaTriggerHit");
  m_tree->Branch("det_name", &b_det_name, "det_name/C");
  m_tree->Branch("det_id"  , &b_det_id  ,   "det_id/I");
  m_tree->Branch("ele_id"  , &b_ele_id  ,   "ele_id/I");
  m_tree->Branch("time"    , &b_time    ,     "time/D");

  ostringstream oss;
  GeomSvc* geom = GeomSvc::instance();
  m_list_det_id.clear();
  for (unsigned int i_det = 0; i_det < m_list_det_name.size(); i_det++) {
    string name = m_list_det_name[i_det];
    int id = geom->getDetectorID(name);
    if (id <= 0) {
      cerr << "!ERROR!  AnaTriggerHit::InitRun():  Invalid ID (" << id << ").  Probably the detector name that you specified in 'm_list_det_name' (" << name << ") is not valid.  Abort." << endl;
      exit(1);
    }
    m_list_det_id.push_back(id);
    int n_ele = geom->getPlaneNElements(id);
    cout << "  " << setw(6) << name << " = " << id << endl;

    oss.str("");
    oss << "h1_ele_id_" << name;
    m_h1_ele[i_det] = new TH1D(oss.str().c_str(), "", n_ele, 0.5, n_ele+0.5);
    oss.str("");
    oss << "V1495 Level-" << m_level << " " << name << ";Element ID;Hit count";
    m_h1_ele[i_det]->SetTitle(oss.str().c_str());

    oss.str("");
    oss << "h1_nhit_" << name;
    m_h1_nhit[i_det] = new TH1D(oss.str().c_str(), "", 20, -0.5, 19.5);
    oss.str("");
    oss << "V1495 Level-" << m_level << " " << name << ";N of hits/plane/event;Hit count";
    m_h1_nhit[i_det]->SetTitle(oss.str().c_str());

    oss.str("");
    oss << "h1_time_" << name;
    m_h1_time[i_det] = new TH1D(oss.str().c_str(), "", 2000, -0.5, 1999.5);
    oss.str("");
    oss << "V1495 Level-" << m_level << " " << name << ";TDC time (ns);Hit count";
    m_h1_time[i_det]->SetTitle(oss.str().c_str());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTriggerHit::process_event(PHCompositeNode* topNode)
{
  //int spill_id = m_evt->get_spill_id();
  //int event_id = m_evt->get_event_id();

  ///
  /// Event selection
  ///
//  if (! m_evt->get_trigger(SQEvent::NIM2)) {
//    return Fun4AllReturnCodes::EVENT_OK;
//  }

  ///
  /// Get & fill the hit info
  ///
  for (unsigned int i_det = 0; i_det < m_list_det_name.size(); i_det++) {
    strncpy(b_det_name, m_list_det_name[i_det].c_str(), sizeof(b_det_name));
    b_det_id = m_list_det_id[i_det];
    shared_ptr<SQHitVector> hv(UtilSQHit::FindHits(m_hit_vec, b_det_id));
    for (SQHitVector::ConstIter it = hv->begin(); it != hv->end(); it++) {
      if ((*it)->get_level() != m_level) continue;
      if (m_intime && ! (*it)->is_in_time()) continue;

      b_ele_id = (*it)->get_element_id();
      b_time   = (*it)->get_tdc_time  ();
      m_tree->Fill();

      m_h1_ele [i_det]->Fill(b_ele_id);
      m_h1_time[i_det]->Fill(b_time);
    }
    m_h1_nhit[i_det]->Fill(hv->size());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTriggerHit::End(PHCompositeNode* topNode)
{
  ostringstream oss;
  TCanvas* c1 = new TCanvas("c1", "");
  c1->SetGrid();
  for (unsigned int i_det = 0; i_det < m_list_det_id.size(); i_det++) {
    m_h1_ele[i_det]->Draw();
    oss.str("");
    oss << "result/" << m_h1_ele[i_det]->GetName() << ".png";
    c1->SaveAs(oss.str().c_str());

    m_h1_nhit[i_det]->Draw();
    oss.str("");
    oss << "result/" << m_h1_nhit[i_det]->GetName() << ".png";
    c1->SaveAs(oss.str().c_str());

    m_h1_time[i_det]->Draw();
    oss.str("");
    oss << "result/" << m_h1_time[i_det]->GetName() << ".png";
    c1->SaveAs(oss.str().c_str());
  }
  delete c1;

  m_file->cd();
  m_file->Write();
  m_file->Close();
  
  return Fun4AllReturnCodes::EVENT_OK;
}
