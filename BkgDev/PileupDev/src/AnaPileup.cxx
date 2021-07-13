/**
 * \class AnaPileup
 * \ module for track quality assurance
 * \author Abinash Pun
 *
 * 
 */


#include "AnaPileup.h"

#include <interface_main/SQHit.h>
#include <interface_main/SQHit_v1.h>
#include <interface_main/SQMCHit_v1.h>
#include <interface_main/SQHitMap_v1.h>
#include <interface_main/SQHitVector_v1.h>
#include <interface_main/SQEvent_v1.h>
#include <interface_main/SQRun_v1.h>
#include <interface_main/SQSpill_v1.h>
#include <interface_main/SQSpillMap_v1.h>

#include <ktracker/SRecEvent.h>
#include <geom_svc/GeomSvc.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4HitDefs.h>
#include <g4main/PHG4VtxPoint.h>

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include <cstring>
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <limits>
#include <tuple>

#include <boost/lexical_cast.hpp>

#define NDET 62
//#define LogDebug(exp)		std::cout<<"DEBUG: "  <<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
#define LogError(exp)		std::cout<<"ERROR: "  <<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
#define LogWarning(exp)	    std::cout<<"WARNING: "<<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl

using namespace std;

AnaPileup::AnaPileup(const std::string& name) :
  SubsysReco(name),
  _hit_container_type("Vector"),
  _event(0),
  _run_header(nullptr),
  _spill_map(nullptr),
  _event_header(nullptr),
  _hit_map(nullptr),
  _hit_vector(nullptr),
  _out_name("eval.root")
{
 
}

int AnaPileup::Init(PHCompositeNode* topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaPileup::InitRun(PHCompositeNode* topNode) {

  ResetEvalVars();
  InitEvalTree();

  p_geomSvc = GeomSvc::instance();

  int ret = GetNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaPileup::process_event(PHCompositeNode* topNode) {
  
  //int ret = Fun4AllReturnCodes::ABORTRUN;

//  if(_recEvent) {    
    //ret = RecoEval(topNode);
    Eval(topNode);
  //  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  //}

  ++_event;

  return Fun4AllReturnCodes::EVENT_OK;
}


//play ground for Abi==============================================
int AnaPileup::Eval(PHCompositeNode* topNode)
{
  ResetEvalVars();
	
  if(_truth) {
    for(auto iter=_truth->GetPrimaryParticleRange().first;
	iter!=_truth->GetPrimaryParticleRange().second;
	++iter) {
      PHG4Particle * par = iter->second; 
     
      int trk_id = par->get_track_id();

      pid[n_tracks] = par->get_pid();
      int vtx_id =  par->get_vtx_id();
		
      PHG4VtxPoint* vtx = _truth->GetVtx(vtx_id);
      gvx[n_tracks] = vtx->get_x();
      gvy[n_tracks] = vtx->get_y();
      gvz[n_tracks] = vtx->get_z();

      TVector3 mom_truth(par->get_px(), par->get_py(), par->get_pz());
      gpx[n_tracks] = par->get_px();
      gpy[n_tracks] = par->get_py();
      gpz[n_tracks] = par->get_pz();
      gpt[n_tracks] = mom_truth.Pt();
      geta[n_tracks] = mom_truth.Eta();
      gphi[n_tracks] = mom_truth.Phi();

      ++n_tracks;
      //if(n_tracks>=10000) break;
             
    }//truth loop
	

  }//truth condition
 
    int h1Bhit = 0;
    int h2Bhit = 0;
    int h3Bhit = 0;
    int h4Bhit = 0;

    int h1Thit = 0;
    int h2Thit = 0;
    int h3Thit = 0;
    int h4Thit = 0;

    int d1xhit = 0;
    int d2xhit = 0;
    int d3xhit = 0;
    
    for(int ihit=0; ihit<_hit_vector->size(); ++ihit) {
    SQHit *sqhit = _hit_vector->at(ihit);
    int sq_detid = sqhit->get_detector_id();

     /// ST.1 HODOS 
     if (sq_detid == 31){
	 new((*pos_H1B)[h1Bhit]) TVector3(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
         new((*mom_H1B)[h1Bhit]) TVector3(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
         ++h1Bhit;
     }
     if (sq_detid == 32){
         new((*pos_H1T)[h1Thit]) TVector3(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
         new((*mom_H1T)[h1Thit]) TVector3(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
         ++h1Thit;
     }
     
     /// ST.2 HODOS
     if (sq_detid == 37){
         new((*pos_H2B)[h2Bhit]) TVector3(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
         new((*mom_H2B)[h2Bhit]) TVector3(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
         ++h2Bhit;
     }
     if (sq_detid == 38){
         new((*pos_H2T)[h2Thit]) TVector3(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
         new((*mom_H2T)[h2Thit]) TVector3(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
         ++h2Thit;
     }


     /// ST.3 HODOS
     if (sq_detid == 39){
         new((*pos_H3B)[h3Bhit]) TVector3(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
         new((*mom_H3B)[h3Bhit]) TVector3(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());      
         ++h3Bhit;
     }
     if (sq_detid == 40){
         new((*pos_H3T)[h3Thit]) TVector3(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
         new((*mom_H3T)[h3Thit]) TVector3(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
         ++h3Thit;
     }

     /// ST.4 HODOS 
     if (sq_detid == 45){
         new((*pos_H4B)[h4Bhit]) TVector3(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
         new((*mom_H4B)[h4Bhit]) TVector3(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
         ++h4Bhit;
     }
     if (sq_detid == 46){
         new((*pos_H4T)[h4Thit]) TVector3(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
         new((*mom_H4T)[h4Thit]) TVector3(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
         ++h4Thit;
     }

     /// DOX CHAMBER
     if (sq_detid ==3 ) {
         new((*pos_D0X)[d1xhit]) TVector3(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
         new((*mom_D0X)[d1xhit]) TVector3(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
         ++d1xhit;
     }

     /// D2XP CHAMBER
     if (sq_detid ==15 ) {
         new((*pos_D2Xp)[d2xhit]) TVector3(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
         new((*mom_D2Xp)[d2xhit]) TVector3(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
         ++d2xhit;
     }

     /// D3MXP OR DMPXP CHAMBER
     if (sq_detid == 27 ||sq_detid == 21 ) {
         new((*pos_D3mXp_OR_D3pXp)[d3xhit]) TVector3(sqhit->get_truth_x(), sqhit->get_truth_y(), sqhit->get_truth_z());
         new((*mom_D3mXp_OR_D3pXp)[d3xhit]) TVector3(sqhit->get_truth_px(), sqhit->get_truth_py(), sqhit->get_truth_pz());
         ++d3xhit;
     }
 
  }

///Reco information
      SRecTrack* Best_recTrack = NULL;
      n_recTracks = _recEvent->getNTracks();
	
  
       if(n_recTracks>0){
	 for(int i = 0; i < n_recTracks; ++i) {
      	  SRecTrack* recTrack = &_recEvent->getTrack(i);                      

	TVector3 rec_mom = recTrack->getTargetMom();
	nhits[i] =recTrack->getNHits();	
	charge[i] = recTrack->getCharge();

	TVector3 rec_vtx = recTrack->getTargetPos();
 
	rec_vx[i]  = rec_vtx.X();
	rec_vy[i]  = rec_vtx.Y();
	rec_vz[i]  = rec_vtx.Z();
    		
	rec_px[i]  = rec_mom.Px();
	rec_py[i]  = rec_mom.Py();
	rec_pz[i]  = rec_mom.Pz();
	rec_pt[i]  = rec_mom.Pt();
	rec_eta[i] = rec_mom.Eta();
	rec_phi[i] = rec_mom.Phi();

	nhits_st1[i] = recTrack->getNHitsInStation(1);
	nhits_st2[i] = recTrack->getNHitsInStation(2);
        nhits_st3[i] = recTrack->getNHitsInStation(3);
	    
	chisq_st1[i] = recTrack->getChisq();
	prob_st1[i] = recTrack->getProb();                
	quality[i] = recTrack->getQuality(); 

	}

      }//if best reco track


      
  _qa_tree->Fill();

  (*pos_H1B).Delete();
  (*pos_H2B).Delete();
  (*pos_H3B).Delete();
  (*pos_H4B).Delete();

  (*pos_H1T).Delete();
  (*pos_H2T).Delete();
  (*pos_H3T).Delete();
  (*pos_H4T).Delete();

  (*mom_H1B).Delete();
  (*mom_H2B).Delete();
  (*mom_H3B).Delete();
  (*mom_H4B).Delete();
  
  (*mom_H1T).Delete();
  (*mom_H2T).Delete();
  (*mom_H3T).Delete();
  (*mom_H4T).Delete();

  (*pos_D0X).Delete();
  (*pos_D2Xp).Delete();
  (*pos_D3mXp_OR_D3pXp).Delete();

  (*mom_D0X).Delete();
  (*mom_D2Xp).Delete();
  (*mom_D3mXp_OR_D3pXp).Delete();

  return Fun4AllReturnCodes::EVENT_OK;
}

///===========================
int AnaPileup::End(PHCompositeNode* topNode) {
  if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
    std::cout << "AnaPileup::End" << std::endl;

  PHTFileServer::get().cd(_out_name.c_str());
  _qa_tree->Write();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaPileup::InitEvalTree() {
 
  PHTFileServer::get().open(_out_name.c_str(), "RECREATE");

  mom_H1T = new TClonesArray("TVector3");
  pos_H1T = new TClonesArray("TVector3");
  mom_H1B = new TClonesArray("TVector3");
  pos_H1B = new TClonesArray("TVector3");

  mom_H2T = new TClonesArray("TVector3");
  pos_H2T = new TClonesArray("TVector3");
  mom_H2B = new TClonesArray("TVector3");
  pos_H2B = new TClonesArray("TVector3");

  mom_H3T = new TClonesArray("TVector3");
  pos_H3T = new TClonesArray("TVector3");
  mom_H3B = new TClonesArray("TVector3");
  pos_H3B = new TClonesArray("TVector3");

  mom_H4T = new TClonesArray("TVector3");
  pos_H4T = new TClonesArray("TVector3");
  mom_H4B = new TClonesArray("TVector3");
  pos_H4B = new TClonesArray("TVector3");


  mom_D0X = new TClonesArray("TVector3");
  pos_D0X = new TClonesArray("TVector3");

  mom_D2Xp = new TClonesArray("TVector3");
  pos_D2Xp = new TClonesArray("TVector3");

  mom_D3mXp_OR_D3pXp = new TClonesArray("TVector3");
  pos_D3mXp_OR_D3pXp = new TClonesArray("TVector3");

/*  mom_H4B = new TClonesArray("TVector3");
  pos_H4B = new TClonesArray("TVector3");
*/  

  ///For the pileup QA tree======================
  _qa_tree = new TTree("pileup_ana", "Analsysis of reconstruction and simulation from piling up bkg interactions");

  _qa_tree->Branch("n_tracks",      &n_tracks,           "n_tracks/I");
  _qa_tree->Branch("n_recTracks",   &n_recTracks,        "n_recTracks/I");

  ///Generated Truth info
  _qa_tree->Branch("pid",           pid,                 "pid[n_tracks]/I");
  _qa_tree->Branch("gvx",           gvx,                 "gvx[n_tracks]/F");
  _qa_tree->Branch("gvy",           gvy,                 "gvy[n_tracks]/F");
  _qa_tree->Branch("gvz",           gvz,                 "gvz[n_tracks]/F");
  _qa_tree->Branch("gpx",           gpx,                 "gpx[n_tracks]/F");
  _qa_tree->Branch("gpy",           gpy,                 "gpy[n_tracks]/F");
  _qa_tree->Branch("gpz",           gpz,                 "gpz[n_tracks]/F");
  _qa_tree->Branch("gpt",           gpt,                 "gpt[n_tracks]/F");
  _qa_tree->Branch("geta",          geta,                "geta[n_tracks]/F");
  _qa_tree->Branch("gphi",          gphi,                "gphi[n_tracks]/F");

  ///HODO HITS
  _qa_tree->Branch("mom_H1T", &mom_H1T, 256000, 99);
  _qa_tree->Branch("pos_H1T", &pos_H1T, 256000, 99);
  _qa_tree->Branch("mom_H1B", &mom_H1B, 256000, 99);
  _qa_tree->Branch("pos_H1B", &pos_H1B, 256000, 99);

  _qa_tree->Branch("mom_H2T", &mom_H2T, 256000, 99);
  _qa_tree->Branch("pos_H2T", &pos_H2T, 256000, 99);
  _qa_tree->Branch("mom_H2B", &mom_H2B, 256000, 99);
  _qa_tree->Branch("pos_H2B", &pos_H2B, 256000, 99);

  _qa_tree->Branch("mom_H3T", &mom_H3T, 256000, 99);
  _qa_tree->Branch("pos_H3T", &pos_H3T, 256000, 99);
  _qa_tree->Branch("mom_H3B", &mom_H3B, 256000, 99);
  _qa_tree->Branch("pos_H3B", &pos_H3B, 256000, 99);

  _qa_tree->Branch("mom_H4T", &mom_H4T, 256000, 99);
  _qa_tree->Branch("pos_H4T", &pos_H4T, 256000, 99);
  _qa_tree->Branch("mom_H4B", &mom_H4B, 256000, 99);
  _qa_tree->Branch("pos_H4B", &pos_H4B, 256000, 99);

  ///CHAMBER HITS
  _qa_tree->Branch("mom_D0X", &mom_D0X, 256000, 99);
  _qa_tree->Branch("pos_D0X", &pos_D0X, 256000, 99);

  _qa_tree->Branch("mom_D2Xp", &mom_D2Xp, 256000, 99);
  _qa_tree->Branch("pos_D2Xp", &pos_D2Xp, 256000, 99);

  _qa_tree->Branch("mom_D3mXp_OR_D3pXp", &mom_D3mXp_OR_D3pXp, 256000, 99);
  _qa_tree->Branch("pos_D3mXp_OR_D3pXp", &pos_D3mXp_OR_D3pXp, 256000, 99);

  ///Reco info in vertex
  _qa_tree->Branch("rec_vx",            rec_vx,                  "rec_vx[n_recTracks]/F");
  _qa_tree->Branch("rec_vy",            rec_vy,                  "rec_vy[n_recTracks]/F");
  _qa_tree->Branch("rec_vz",            rec_vz,                  "rec_vz[n_recTracks]/F");
  _qa_tree->Branch("rec_px",            rec_px,                  "rec_px[n_recTracks]/F");
  _qa_tree->Branch("rec_py",            rec_py,                  "rec_py[n_recTracks]/F");
  _qa_tree->Branch("rec_pz",            rec_pz,                  "rec_pz[n_recTracks]/F");
  _qa_tree->Branch("rec_pt",            rec_pt,                  "rec_pt[n_recTracks]/F");
  _qa_tree->Branch("rec_eta",           rec_eta,                 "rec_eta[n_recTracks]/F");
  _qa_tree->Branch("rec_phi",           rec_phi,                 "rec_phi[n_recTracks]/F");

  _qa_tree->Branch("chisq_st1",      chisq_st1,            "chisq_st1[n_recTracks]/F");
  _qa_tree->Branch("prob_st1",       prob_st1,             "prob_st1[n_recTracks]/F");
  _qa_tree->Branch("quality",        quality,             "quality[n_recTracks]/F");

  _qa_tree->Branch("nhits",         nhits,               "nhits[n_recTracks]/I");
  _qa_tree->Branch("nhits_st1",         nhits_st1,               "nhits_st1[n_recTracks]/I");
  _qa_tree->Branch("nhits_st2",         nhits_st2,               "nhits_st2[n_recTracks]/I");
  _qa_tree->Branch("nhits_st3",         nhits_st3,               "nhits_st3[n_recTracks]/I");
  _qa_tree->Branch("charge",        charge,              "charge[n_recTracks]/I");

  mom_H1T->BypassStreamer();
  pos_H1T->BypassStreamer(); 
  mom_H1B->BypassStreamer();
  pos_H1B->BypassStreamer();

  mom_H2T->BypassStreamer();
  pos_H2T->BypassStreamer();
  mom_H2B->BypassStreamer();
  pos_H2B->BypassStreamer();

  mom_H2T->BypassStreamer();
  pos_H2T->BypassStreamer();
  mom_H2B->BypassStreamer();
  pos_H2B->BypassStreamer();

  mom_H2T->BypassStreamer();
  pos_H2T->BypassStreamer();
  mom_H2B->BypassStreamer();
  pos_H2B->BypassStreamer();
 
  mom_D0X->BypassStreamer();
  pos_D0X->BypassStreamer();

  mom_D2Xp->BypassStreamer();
  pos_D2Xp->BypassStreamer();

  mom_D3mXp_OR_D3pXp->BypassStreamer();
  pos_D3mXp_OR_D3pXp->BypassStreamer();

  return 0;
}

int AnaPileup::ResetEvalVars() {

  run_id = std::numeric_limits<int>::max();
  spill_id = std::numeric_limits<int>::max();
  event_id = std::numeric_limits<int>::max();
  emu_trigger = 0;


    n_tracks = 0;
   
  for(int i=0; i<100; ++i) {
    pid[i]       = std::numeric_limits<int>::max();
    gvx[i]        = std::numeric_limits<float>::max();
    gvy[i]        = std::numeric_limits<float>::max();
    gvz[i]        = std::numeric_limits<float>::max();
    gpx[i]        = std::numeric_limits<float>::max();
    gpy[i]        = std::numeric_limits<float>::max();
    gpz[i]        = std::numeric_limits<float>::max();
    gpt[i]        = std::numeric_limits<float>::max();
    geta[i]       = std::numeric_limits<float>::max();
    gphi[i]       = std::numeric_limits<float>::max();
    nhits[i]      = std::numeric_limits<int>::max();
    charge[i]     = std::numeric_limits<int>::max();

    rec_vx[i]         = std::numeric_limits<float>::max();
    rec_vy[i]         = std::numeric_limits<float>::max();
    rec_vz[i]         = std::numeric_limits<float>::max();
    rec_px[i]         = std::numeric_limits<float>::max();
    rec_py[i]         = std::numeric_limits<float>::max();
    rec_pz[i]         = std::numeric_limits<float>::max();
    rec_pt[i]         = std::numeric_limits<float>::max();
    rec_eta[i]        = std::numeric_limits<float>::max();
    rec_phi[i]        = std::numeric_limits<float>::max();


  }

  return 0;
}

int AnaPileup::GetNodes(PHCompositeNode* topNode) {

  _run_header = findNode::getClass<SQRun>(topNode, "SQRun");
  if (!_run_header) {
    LogError("!_run_header");
    //return Fun4AllReturnCodes::ABORTEVENT;
  }

  _spill_map = findNode::getClass<SQSpillMap>(topNode, "SQSpillMap");
  if (!_spill_map) {
    LogError("!_spill_map");
    //return Fun4AllReturnCodes::ABORTEVENT;
  }

  _event_header = findNode::getClass<SQEvent>(topNode, "SQEvent");
  if (!_event_header) {
    LogError("!_event_header");
    //return Fun4AllReturnCodes::ABORTEVENT;
  }

  if(_hit_container_type.find("Map") != std::string::npos) {
    _hit_map = findNode::getClass<SQHitMap>(topNode, "SQHitMap");
    if (!_hit_map) {
      LogError("!_hit_map");
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if(_hit_container_type.find("Vector") != std::string::npos) {
    _hit_vector = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
    if (!_hit_vector) {
      LogError("!_hit_vector");
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  _truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truth) {
    LogError("!_truth");
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _recEvent = findNode::getClass<SRecEvent>(topNode, "SRecEvent");
  if (!_recEvent) {
    LogError("!_recEvent");
    //return Fun4AllReturnCodes::ABORTEVENT;
  }
 
  g4hc_d1x  = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D1X");
  g4hc_d2xp  = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D2Xp");
  g4hc_d3px = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D3pXp");
  g4hc_d3mx = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D3mXp");
  if (! g4hc_d1x) g4hc_d1x = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D0X");

  if ( !g4hc_d1x || !g4hc_d3px || !g4hc_d3mx) {
    cout << "Failed at getting nodes: "<< g4hc_d1x << " " << g4hc_d3px << " " << g4hc_d3mx << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

//hodoscope for the acceptance study 
  g4hc_h1t  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H1T");
  g4hc_h1b  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H1B");
  g4hc_h2t  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H2T");
  g4hc_h2b  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H2B");
  g4hc_h3t  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H3T");
  g4hc_h3b  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H3B");
  g4hc_h4t  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H4T");
  g4hc_h4b  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H4B");

  if (!g4hc_h1t || !g4hc_h1b || !g4hc_h2t || !g4hc_h2b ||
      !g4hc_h3t || !g4hc_h3b || !g4hc_h4t || !g4hc_h4b   ) {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

///Prop tubes hits
  g4hc_p1y1  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1Y1");
  g4hc_p1y2  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1Y2");
  g4hc_p1x1  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1X1");
  g4hc_p1x2  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P1X2");
  g4hc_p2x1  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2X1");
  g4hc_p2x2  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2X2");
  g4hc_p2y1  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2Y1");
  g4hc_p2y2  = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_P2Y2");

  if (!g4hc_p1y1 || !g4hc_p1y2 || !g4hc_p1x1 || !g4hc_p1x2 ||
      !g4hc_p2x1 || !g4hc_p2x2 || !g4hc_p2y1 || !g4hc_p2y2   ) {
    return Fun4AllReturnCodes::ABORTEVENT;
  }



  return Fun4AllReturnCodes::EVENT_OK;
}

/*
//For finding g4hit information in stations (following Kenichi's truth node maker)
bool AnaPileup::FindG4HitAtStation(const int trk_id, const PHG4HitContainer* g4hc, TVector3* pos, TLorentzVector* mom)
{
  //const double M_MU = 0.1056583745;
  PHG4HitContainer::ConstRange range = g4hc->getHits();
  for (PHG4HitContainer::ConstIterator it = range.first; it != range.second; it++) {
    PHG4Hit* hit = it->second;
    if (hit->get_trkid() == trk_id) {
      pos->SetXYZ (hit->get_x(0)     , hit->get_y(0)     , hit->get_z(0)           );
      mom->SetXYZM(hit->get_px(0),     hit->get_py(0),     hit->get_pz(0), M_MU);
      return true;
    }
  }
  return false;
}


//Function for finding best reco track
SRecTrack* AnaPileup::FindBestMomRecTrack(SRecEvent *recEvent,  const float true_TargetP)
{
  double dP = 100.;
  double hold_dP = 99999.;
  
  SRecTrack* Best_recTrack =  NULL;
  for(int itrack=0; itrack<recEvent->getNTracks(); ++itrack){   
    if (hold_dP>dP) hold_dP = dP;
    SRecTrack *recTrack = &recEvent->getTrack(itrack);
    dP = fabs(true_TargetP - recTrack->getTargetMom().Mag());
   
    //Finding out best match track in terms of energy
    if(dP-hold_dP<0.) Best_recTrack = recTrack;  
  }
  return Best_recTrack;
  
}


//Function to find common hit ids for reco and truth tracks
int AnaPileup::FindCommonHitIDs(vector<int>& hitidvec1, vector<int>& hitidvec2)
{
  //This function assumes the input vectors have been sorted
  auto iter = hitidvec1.begin();
  auto jter = hitidvec2.begin();

  int nCommon = 0;
  while(iter != hitidvec1.end() && jter != hitidvec2.end()) {
    if(*iter < *jter) {
      ++iter;
    } else {
      if(!(*jter < *iter)) {
        ++nCommon;
        ++iter;
      }
      ++jter;
    }
  }

  return nCommon;
}


//functions for the acceptance
bool AnaPileup::FindG4HitAtHodo(const int trk_id, const PHG4HitContainer* g4hc, TVector3* pos, TLorentzVector* mom)
{
  //const double M_MU = 0.1056583745; 
  PHG4HitContainer::ConstRange range = g4hc->getHits();
  for (PHG4HitContainer::ConstIterator it = range.first; it != range.second; it++) {
    PHG4Hit* hit = it->second;
    if (hit->get_trkid() == trk_id) {
      pos->SetXYZ (hit->get_x(0)     , hit->get_y(0)     , hit->get_z(0) );
      mom->SetXYZM(hit->get_px(0),     hit->get_py(0),     hit->get_pz(0), M_MU);
     return true;
    }
  }
  return false;
}


bool AnaPileup::FindG4HitAtProp(const int trk_id, const PHG4HitContainer* g4hc)
{
  PHG4HitContainer::ConstRange range = g4hc->getHits();
  for (PHG4HitContainer::ConstIterator it = range.first; it != range.second; it++) {
    PHG4Hit* hit = it->second;
    if (hit->get_trkid() == trk_id) {
      return true;
    }
  }
  return false;
}
*/
