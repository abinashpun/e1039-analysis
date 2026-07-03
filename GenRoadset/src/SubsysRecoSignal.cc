#include <iomanip>
#include <TFile.h>
#include <TTree.h>
#include <interface_main/SQMCEvent.h>
#include <interface_main/SQTrackVector.h>
#include <interface_main/SQDimuonVector.h>
#include <interface_main/SQHitVector.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <UtilAna/UtilTrack.h>
#include <UtilAna/UtilDimuon.h>
#include <UtilAna/UtilTrigger.h>
#include <geom_svc/GeomSvc.h>
#include "SubsysRecoSignal.h"
using namespace std;

SubsysRecoSignal::SubsysRecoSignal(const std::string &name)
  : SubsysReco(name)
{
  ;
}

int SubsysRecoSignal::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int SubsysRecoSignal::InitRun(PHCompositeNode* topNode)
{
  mi_evt_true = findNode::getClass<SQMCEvent     >(topNode, "SQMCEvent");
  mi_vec_trk  = findNode::getClass<SQTrackVector >(topNode, "SQTruthTrackVector");
  mi_vec_dim  = findNode::getClass<SQDimuonVector>(topNode, "SQTruthDimuonVector");
  mi_vec_hit  = findNode::getClass<SQHitVector   >(topNode, "SQHitVector");
  if (!mi_evt_true || !mi_vec_trk || !mi_vec_dim || !mi_vec_hit) return Fun4AllReturnCodes::ABORTEVENT;

  mo_file = new TFile("signal_tree.root", "RECREATE");
  mo_tree = new TTree("signal_tree", "Created by SubsysRecoSignal");
  mo_tree->Branch("signal_data", &mo_sig);

  return Fun4AllReturnCodes::EVENT_OK;
}

int SubsysRecoSignal::process_event(PHCompositeNode* topNode)
{
  static unsigned int n_evt = 0;
  if    (++n_evt % 100000 == 0) cout << n_evt << endl;
  else if (n_evt %  10000 == 0) cout << " . " << flush;

  mo_sig.weight = mi_evt_true->get_weight();

  for (unsigned int ii = 0; ii < mi_vec_dim->size(); ii++) {
    SQDimuon* dim = mi_vec_dim->at(ii);
    int id_trk_pos = dim->get_track_id_pos();
    int id_trk_neg = dim->get_track_id_neg();
    SQTrack* trk_pos = UtilTrack::FindTrackByID(mi_vec_trk, id_trk_pos, true);
    SQTrack* trk_neg = UtilTrack::FindTrackByID(mi_vec_trk, id_trk_neg, true);
    if (! trk_pos || ! trk_neg) {
      cerr << "Invalid true tracks!  Abort." << endl;
      exit(1);
    }
    
    shared_ptr<SQHitVector> hv_trk_pos(UtilTrack::FindHodoHitsOfTrack(mi_vec_hit, id_trk_pos));
    shared_ptr<SQHitVector> hv_trk_neg(UtilTrack::FindHodoHitsOfTrack(mi_vec_hit, id_trk_neg));
    int road_pos = UtilTrigger::ExtractRoadID(hv_trk_pos.get());
    int road_neg = UtilTrigger::ExtractRoadID(hv_trk_neg.get());
    //cout << "  pos: " << id_trk_pos << " " << hv_trk_pos->size() << " " << road_pos << "\n"
    //     << "  neg: " << id_trk_neg << " " << hv_trk_neg->size() << " " << road_neg << endl;
    //if (road_pos * road_neg >= 0) continue // Keep only T+B/B+T pairs

    UtilDimuon::CalcVar(dim, mo_sig.mass, mo_sig.pT, mo_sig.x1, mo_sig.x2, mo_sig.xF);
    mo_sig.dz       = dim->get_pos().Z();
    mo_sig.mom      = dim->get_mom().P();
    mo_sig.phi      = dim->get_mom().Phi();
    mo_sig.theta    = dim->get_mom().Theta();
    mo_sig.road_pos = road_pos;
    mo_sig.road_neg = road_neg;
    mo_tree->Fill();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int SubsysRecoSignal::End(PHCompositeNode* topNode)
{
  mo_file->cd();
  mo_file->Write();
  mo_file->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}
