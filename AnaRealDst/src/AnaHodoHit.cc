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
#include "AnaHodoHit.h"
using namespace std;

/// List of detectors that you want to analyze
const vector<string> AnaHodoHit::m_list_det_name = { "H1T", "H1B", "H2T", "H2B" };

AnaHodoHit::AnaHodoHit(const std::string& name)
  : SubsysReco(name)
  , m_file(0)
  , m_tree(0)
{
  ;
}

int AnaHodoHit::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaHodoHit::InitRun(PHCompositeNode* topNode)
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
  m_tree = new TTree("tree", "Created by AnaHodoHit");
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
      cerr << "!ERROR!  AnaHodoHit::InitRun():  Invalid ID (" << id << ").  Probably the detector name that you specified in 'list_det_name' (" << name << ") is not valid.  Abort." << endl;
      exit(1);
    }
    m_list_det_id.push_back(id);
    int n_ele = geom->getPlaneNElements(id);
    cout << "  " << setw(6) << name << " = " << id << endl;

    oss.str("");
    oss << "h1_ele_" << name;
    m_h1_ele[i_det] = new TH1D(oss.str().c_str(), "", n_ele, 0.5, n_ele+0.5);
    oss.str("");
    oss << name << ";Element ID;Hit count";
    m_h1_ele[i_det]->SetTitle(oss.str().c_str());

    oss.str("");
    oss << "h1_nhit_" << name;
    m_h1_nhit[i_det] = new TH1D(oss.str().c_str(), "", 10, -0.5, 9.5);
    oss.str("");
    oss << name << ";N of hits/plane/event;Hit count";
    m_h1_nhit[i_det]->SetTitle(oss.str().c_str());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaHodoHit::process_event(PHCompositeNode* topNode)
{
  //int spill_id = m_evt->get_spill_id();
  //int event_id = m_evt->get_event_id();

  ///
  /// Event selection
  ///
  if (! m_evt->get_trigger(SQEvent::NIM2)) {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  /// Example of event selection using H3T and H3B hodo hits
  static int did_h3t = 0;
  static int did_h3b = 0;
  if (did_h3t == 0) { // Resolve IDs only once.
    GeomSvc* geom = GeomSvc::instance();
    did_h3t = geom->getDetectorID("H3T");
    did_h3b = geom->getDetectorID("H3B");
    cout << "H3T = " << did_h3t << ", H3B = " << did_h3b << endl;
  }
  shared_ptr<SQHitVector> hv_h3t(UtilSQHit::FindHits(m_hit_vec, did_h3t));
  shared_ptr<SQHitVector> hv_h3b(UtilSQHit::FindHits(m_hit_vec, did_h3b));
  if (hv_h3t->size() + hv_h3b->size() != 1) return Fun4AllReturnCodes::EVENT_OK;

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

      m_h1_ele[i_det]->Fill(b_ele);
    }
    m_h1_nhit[i_det]->Fill(hv->size());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaHodoHit::End(PHCompositeNode* topNode)
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
  }
  delete c1;

  m_file->cd();
  m_file->Write();
  m_file->Close();
  
  return Fun4AllReturnCodes::EVENT_OK;
}
