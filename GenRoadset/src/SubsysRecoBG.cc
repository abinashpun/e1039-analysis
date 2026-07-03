#include <iomanip>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <interface_main/SQRun.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQMCEvent.h>
#include <interface_main/SQHitVector.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <UtilAna/UtilSQHit.h>
#include "SubsysRecoBG.h"
using namespace std;

SubsysRecoBG::SubsysRecoBG(const std::string &name)
  : SubsysReco(name)
  , m_mode(DEFAULT)
  , mi_evt(0)
  , mi_mc_evt(0)
  , mi_vec_hit(0)
  , mo_file(0)
  , mo_tree(0)
{
  ;
}

int SubsysRecoBG::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int SubsysRecoBG::InitRun(PHCompositeNode* topNode)
{
  mi_evt     = findNode::getClass<SQEvent    >(topNode, "SQEvent");
  mi_vec_hit = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if (!mi_evt || !mi_vec_hit) return Fun4AllReturnCodes::ABORTEVENT;

  if (m_mode == FULL_BG) {
    mi_mc_evt = findNode::getClass<SQMCEvent>(topNode, "SQMCEvent");
    if (!mi_mc_evt) return Fun4AllReturnCodes::ABORTEVENT;
  }

  mo_file = new TFile("bg_data.root", "RECREATE");
  mo_tree = new TTree("bg_tree", "Created by SubsysRecoBG");
  mo_tree->Branch("bg_data", &mo_bg);

  return Fun4AllReturnCodes::EVENT_OK;
}

int SubsysRecoBG::process_event(PHCompositeNode* topNode)
{
  mo_bg.run   = mi_evt->get_run_id();
  mo_bg.evt   = mi_evt->get_event_id();
  mo_bg.fpga1 = mi_evt->get_trigger(SQEvent::MATRIX1);

  if (m_mode == FULL_BG) {
    mo_bg.pot_rfp00 = mi_mc_evt->get_cross_section();
  }
  mo_bg.inte_rfp00 = mi_evt->get_qie_rf_intensity(0);
  mo_bg.inte_max   = 0;
  for (int ii = -8; ii <= 8; ii++) {
    int inte = mi_evt->get_qie_rf_intensity(ii);
    if (inte > mo_bg.inte_max) mo_bg.inte_max = inte;
  }

  ExtractHits(mi_vec_hit, "H1T", mo_bg.h1t);
  ExtractHits(mi_vec_hit, "H2T", mo_bg.h2t);
  ExtractHits(mi_vec_hit, "H3T", mo_bg.h3t);
  ExtractHits(mi_vec_hit, "H4T", mo_bg.h4t);
  ExtractHits(mi_vec_hit, "H1B", mo_bg.h1b);
  ExtractHits(mi_vec_hit, "H2B", mo_bg.h2b);
  ExtractHits(mi_vec_hit, "H3B", mo_bg.h3b);
  ExtractHits(mi_vec_hit, "H4B", mo_bg.h4b);

  mo_tree->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

int SubsysRecoBG::End(PHCompositeNode* topNode)
{
  mo_file->cd();
  mo_file->Write();
  mo_file->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}

void SubsysRecoBG::ExtractHits(const SQHitVector* hit_vec, const std::string det_name, std::vector<int>& list_ele)
{
  list_ele.clear();
  shared_ptr<SQHitVector> hv(UtilSQHit::FindFirstHits(hit_vec, det_name, true));
  for (SQHitVector::Iter it = hv->begin(); it != hv->end(); it++) {
    list_ele.push_back( (*it)->get_element_id() );
  }
}
