/**
 * \class AnaTrkQA
 * \ module for track quality assurance
 * \author Haiwang, modified by Abinash Pun
 *
 * 
 */


#include "AnaTrkQA.h"

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

#include <cstring>
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <limits>
#include <tuple>

#include <boost/lexical_cast.hpp>

#define NDET 62
#define LogDebug(exp)		std::cout<<"DEBUG: "  <<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
#define LogError(exp)		std::cout<<"ERROR: "  <<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
#define LogWarning(exp)	    std::cout<<"WARNING: "<<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl

using namespace std;

AnaTrkQA::AnaTrkQA(const std::string& name) :
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

int AnaTrkQA::Init(PHCompositeNode* topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTrkQA::InitRun(PHCompositeNode* topNode) {

  ResetEvalVars();
  InitEvalTree();

  p_geomSvc = GeomSvc::instance();

  int ret = GetNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTrkQA::process_event(PHCompositeNode* topNode) {
  
  int ret = Fun4AllReturnCodes::ABORTRUN;

  if(_recEvent) {    
    ret = RecoEval(topNode);
    ret = RecoEvalv2(topNode);
    if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  }

  ++_event;

  return ret;
}

int AnaTrkQA::RecoEval(PHCompositeNode* topNode)
{
  if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
    std::cout << "Entering AnaTrkQA::RecoEval: " << _event << std::endl;
 
  ResetEvalVars();

  if(_spill_map) {
    auto spill_info = _spill_map->get(spill_id);
    if(spill_info) {
      target_pos = spill_info->get_target_pos();
    } else {
      LogWarning("");
    }
  }

  if(_event_header) {
    event_id    = _event_header->get_event_id();
    emu_trigger = _event_header->get_trigger();
    spill_id    = _event_header->get_spill_id();
    run_id      = _event_header->get_run_id();
  }

  PHG4HitContainer *D1X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D0X");
  if (!D1X_hits)
    D1X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D1X");

  if (!D1X_hits)
    {
      if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
        cout << Name() << " Could not locate g4 hit node " << "G4HIT_D0X or G4HIT_D1X" << endl;
    }

  std::map<int, int> parID_nhits_dc;
  std::map<int, int> parID_nhits_hodo;
  std::map<int, int> parID_nhits_prop;
  std::map<int, int> parID_nhits_dp;

  std::map<int, std::map<int, int> > parID_detid_elmid;

  typedef std::tuple<int, int> ParDetPair;
  std::map<ParDetPair, int> parID_detID_ihit;

  std::map<int, int> hitID_ihit;

  // If using vector, index it first
  if(_hit_vector) {
    for(int ihit=0; ihit<_hit_vector->size(); ++ihit) {
      SQHit *hit = _hit_vector->at(ihit);
      int hitID = hit->get_hit_id();
      hitID_ihit[hitID]      = ihit;

      if(_truth) {
      	int track_id = hit->get_track_id();
      	int det_id = hit->get_detector_id();
      	parID_detID_ihit[std::make_tuple(track_id, det_id)] = ihit;

      	if(hit->get_detector_id() >= 1 and hit->get_detector_id() <=30) {
	  if(parID_nhits_dc.find(track_id)!=parID_nhits_dc.end())
	    parID_nhits_dc[track_id] = parID_nhits_dc[track_id]+1;
	  else
	    parID_nhits_dc[track_id] = 1;
      	}
      	if(hit->get_detector_id() >= 31 and hit->get_detector_id() <=46) {
	  if(parID_nhits_hodo.find(track_id)!=parID_nhits_hodo.end())
	    parID_nhits_hodo[track_id] = parID_nhits_hodo[track_id]+1;
	  else
	    parID_nhits_hodo[track_id] = 1;
      	}
      	if(hit->get_detector_id() >= 47 and hit->get_detector_id() <=54) {
	  if(parID_nhits_prop.find(track_id)!=parID_nhits_prop.end())
	    parID_nhits_prop[track_id] = parID_nhits_prop[track_id]+1;
	  else
	    parID_nhits_prop[track_id] = 1;
      	}
      	if(hit->get_detector_id() >= 55 and hit->get_detector_id() <=62) {
	  if(parID_nhits_dp.find(track_id)!=parID_nhits_dp.end())
	    parID_nhits_dp[track_id] = parID_nhits_dp[track_id]+1;
	  else
	    parID_nhits_dp[track_id] = 1;
      	}
      }
    }
  }

  typedef std::tuple<int, int> ParRecoPair;
  std::map<ParRecoPair, int> parID_recID_nHit;

  typedef std::tuple<int, int> TrkIDNHit;
  std::map<int, TrkIDNHit> parID_bestRecID;
  std::map<int, TrkIDNHit> recID_bestParID;

  if(!_recEvent) {
    LogInfo("!_recEvent");
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if(!_truth) {
    LogInfo("!_truth");
    return Fun4AllReturnCodes::ABORTRUN;
  }

  krecstat = _recEvent->getRecStatus();

  for(int itrack=0; itrack<_recEvent->getNTracks(); ++itrack){
    SRecTrack recTrack = _recEvent->getTrack(itrack);

    if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT) {
      cout << "--------- itrack: " << itrack << " ---------"<< endl;
    }

    // Fill map parID_recID => nHit
    for(int ihit=0; ihit<recTrack.getNHits();++ihit) {
      int hitID = recTrack.getHitIndex(ihit);

      if(Verbosity() >= Fun4AllBase::VERBOSITY_MORE) {
	LogDebug("hitID: " << hitID);
      }

      // signed hitID to hitID
      hitID = abs(hitID);

      //! TODO change back to map?
      if(hitID_ihit.find(hitID)==hitID_ihit.end()) continue;

      SQHit *hit = _hit_vector->at(hitID_ihit[hitID]);

      if(!hit) {
	if(Verbosity() >= Fun4AllBase::VERBOSITY_MORE) {
	  LogWarning("!hit");
	}
      }

                         
      // TODO better way to exclude hodo and prop hits?
      // this is try to exclude hodo and prop hits
      if(hit->get_detector_id() > 30) continue;  	
                       
      int parID = hit->get_track_id();
      //LogInfo(parID);
      if(parID > 9999) continue;
      //LogInfo(parID);

      ParRecoPair key = std::make_tuple(parID, itrack);

      if(parID_recID_nHit.find(key)!=parID_recID_nHit.end())
	parID_recID_nHit[key] = parID_recID_nHit[key]+1;
      else
	parID_recID_nHit[key] = 1;
    }
  }

  for(auto iter=parID_recID_nHit.begin();
      iter!=parID_recID_nHit.end(); ++iter) {
    int parID = std::get<0>(iter->first);
    int recID = std::get<1>(iter->first);
    int nHit  = iter->second;

    if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT) {
      LogInfo("");
      std::cout
	<< parID << " : " << recID
	<< " => " << nHit
	<< std::endl;
    }

    if(parID_bestRecID.find(parID)!=parID_bestRecID.end()) {
      int nHit_current_best  = std::get<1>(parID_bestRecID[parID]);
      if (nHit > nHit_current_best)
	parID_bestRecID[parID] = std::make_tuple(recID, nHit);
    }
    else
      parID_bestRecID[parID] = std::make_tuple(recID, nHit);

    if(recID_bestParID.find(recID)!=recID_bestParID.end()) {
      int nHit_current_best  = std::get<1>(recID_bestParID[recID]);
      if (nHit > nHit_current_best)
	recID_bestParID[recID] = std::make_tuple(parID, nHit);
    }
    else
      recID_bestParID[recID] = std::make_tuple(parID, nHit);
  }

  // Fill Reco Tree
  for(int itrack=0; itrack<_recEvent->getNTracks(); ++itrack){
    SRecTrack recTrack = _recEvent->getTrack(itrack);

    if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT) {
      cout << "--------- itrack: " << itrack << " ---------"<< endl;
    }

    nhits[n_tracks] = recTrack.getNHits();
    charge[n_tracks] = recTrack.getCharge();
    TVector3 rec_vtx = recTrack.getTargetPos();
    vx[n_tracks]  = rec_vtx.X();
    vy[n_tracks]  = rec_vtx.Y();
    vz[n_tracks]  = rec_vtx.Z();
    TVector3 rec_mom = recTrack.getTargetMom();
    px[n_tracks]  = rec_mom.Px();
    py[n_tracks]  = rec_mom.Py();
    pz[n_tracks]  = rec_mom.Pz();
    pt[n_tracks]  = rec_mom.Pt();
    eta[n_tracks] = rec_mom.Eta();
    phi[n_tracks] = rec_mom.Phi();
    chisq_st1[n_tracks] = recTrack.getChisq();
    prob_st1[n_tracks] = recTrack.getProb();                
    quality[n_tracks] = recTrack.getQuality(); 
                     
    {   
      double tx, ty, tz;
      recTrack.getMomentumSt1(tx, ty, tz);
      px_st1[n_tracks] = tx;
      py_st1[n_tracks] = ty;
      pz_st1[n_tracks] = tz;
 
      double x, y;
      recTrack.getPositionSt1(x, y);		
      x_st1[n_tracks] = x;		
      y_st1[n_tracks] = y;

    }

    int trackID = itrack;

    rec_id[n_tracks] = trackID;

    if(recID_bestParID.find(trackID)!=recID_bestParID.end()) {
      ntruhits[n_tracks] = std::get<1>(recID_bestParID[trackID]);
      int parID = std::get<0>(recID_bestParID[trackID]);

      PHG4Particle * par = _truth->GetParticle(parID);

      par_id[n_tracks] = parID;

      pid[n_tracks] = par->get_pid();

      int vtx_id =  par->get_vtx_id();
      PHG4VtxPoint* vtx = _truth->GetVtx(vtx_id);
      gvx[n_tracks] = vtx->get_x();
      gvy[n_tracks] = vtx->get_y();
      gvz[n_tracks] = vtx->get_z();

      TVector3 mom(par->get_px(), par->get_py(), par->get_pz());
      gpx[n_tracks] = par->get_px();
      gpy[n_tracks] = par->get_py();
      gpz[n_tracks] = par->get_pz();
      gpt[n_tracks] = mom.Pt();
      geta[n_tracks] = mom.Eta();
      gphi[n_tracks] = mom.Phi();

      // trackID + detID -> SQHit -> PHG4Hit -> momentum
      for(int det_id=1; det_id<=12; ++det_id) {
	auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	if(iter != parID_detID_ihit.end()) {
	  if(verbosity >= Fun4AllBase::VERBOSITY_A_LOT) {
	    LogDebug("ihit: " << iter->second);
	  }
	  SQHit *hit = _hit_vector->at(iter->second);
	  if(verbosity >= Fun4AllBase::VERBOSITY_A_LOT) {
	    hit->identify();
	  }

	  if(hit and D1X_hits) {
	    PHG4Hit* g4hit =  D1X_hits->findHit(hit->get_g4hit_id());
	    if (g4hit) {
	      if(verbosity >= 2) {
		g4hit->identify();
	      }
	      gx_st1[n_tracks]  = g4hit->get_x(0);
	      gy_st1[n_tracks]  = g4hit->get_y(0);
	      gz_st1[n_tracks]  = g4hit->get_z(0);
	      if(det_id==1){
		sqx_st1[n_tracks] = hit->get_truth_x();
		sqy_st1[n_tracks] = hit->get_truth_y();
		sqz_st1[n_tracks] = hit->get_truth_z();
	      }                             
                      
	      double x0 = hit->get_truth_x()-g4hit->get_px(0)/g4hit->get_pz(0) *hit->get_truth_z();
	      double y0 = hit->get_truth_y()-g4hit->get_py(0)/g4hit->get_pz(0) *hit->get_truth_z();

	      double x00 = x_st1[n_tracks]-px_st1[n_tracks]/pz_st1[n_tracks] *hit->get_truth_z();
	      double y00 = y_st1[n_tracks]-py_st1[n_tracks]/pz_st1[n_tracks] *hit->get_truth_z();

	      if(det_id==1){
		sq_decID[n_tracks]=hit->get_detector_id();
		sq_pos_st1[n_tracks]=hit->get_pos();
		
		sq_drift_st1[n_tracks] = hit->get_drift_distance();
               
		rec_drift_st1[n_tracks] = p_geomSvc->getDCA(hit->get_detector_id(), hit->get_element_id(),  px_st1[n_tracks]/pz_st1[n_tracks], py_st1[n_tracks]/pz_st1[n_tracks], x00, y00);
		
	      }

	      gpx_st1[n_tracks] = g4hit->get_px(0)/1.;
	      gpy_st1[n_tracks] = g4hit->get_py(0)/1.;
	      gpz_st1[n_tracks] = g4hit->get_pz(0)/1.;
                       
                                        

	      break;
	    }
	  }
	}
      }
      gnhits[n_tracks] =
	parID_nhits_dc[parID] +
	parID_nhits_hodo[parID] +
	parID_nhits_prop[parID] +
	parID_nhits_dp[parID];

      gndc[n_tracks] = parID_nhits_dc[parID];
      gnhodo[n_tracks] = parID_nhits_hodo[parID];
      gnprop[n_tracks] = parID_nhits_prop[parID];
      gndp[n_tracks] = parID_nhits_dp[parID];
    }
    ++n_tracks;
    if(n_tracks>=1000) break;
  }

  if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT) LogInfo("parID_recID_nHit mapping finished");

  _tout_reco->Fill();

  if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
    std::cout << "Leaving AnaTrkQA::RecoEval: " << _event << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//play ground for Abi========================
int AnaTrkQA::RecoEvalv2(PHCompositeNode* topNode)
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


      //st1. true hit info
      TVector3 g_pos_st1;
      TLorentzVector g_mom_st1;
      int st1_hit = FindHitAtStation(trk_id, g4hc_d1x, &g_pos_st1, &g_mom_st1);
		
      gx_st1[n_tracks]  = g_pos_st1.X();
      gy_st1[n_tracks]  = g_pos_st1.Y();
      gz_st1[n_tracks]  = g_pos_st1.Z();
     		
      gpx_st1[n_tracks] = g_mom_st1.Px();
      gpy_st1[n_tracks] = g_mom_st1.Py();
      gpz_st1[n_tracks] = g_mom_st1.Pz();

      //st3. trute hit info (work in progress: get it from either d3pXp or d3mXp)
      TVector3 g_pos_st3;
      TLorentzVector g_mom_st3;
      int d3phit = FindHitAtStation(trk_id, g4hc_d3px, &g_pos_st3, &g_mom_st3);
      int d3mhit = FindHitAtStation(trk_id, g4hc_d3mx, &g_pos_st3, &g_mom_st3);
      //create branches to put st3. g4hit info
                 



      //=================SQ Hit information for st1. and st3.
      SQHit *sqhit_st1 = NULL; //just incase need for later on
      SQHit *sqhit_st3 = NULL;
      // Get the sqhit position at st1. , st3 and drift distances
      if(_hit_vector) {
	for(int ihit=0; ihit<_hit_vector->size(); ++ihit) {
	  SQHit *hit = _hit_vector->at(ihit);                    

	  //st. 1
	  if(hit->get_track_id() == trk_id && hit->get_detector_id() == 1) {
	    sqhit_st1 = hit;
     
	    sqx_st1[n_tracks] = hit->get_truth_x();
	    sqy_st1[n_tracks] = hit->get_truth_y();
	    sqz_st1[n_tracks] = hit->get_truth_z();
	    sq_pos_st1[n_tracks]=hit->get_pos();
	    sq_drift_st1[n_tracks] = hit->get_drift_distance();
	  }
		
	  //st. 3
	  if(hit->get_track_id() == trk_id && (hit->get_detector_id() == 27/*D3mXp*/||hit->get_detector_id() == 21/*D3pXp*/)){
	    sqhit_st3 = hit;
	    //create branch to store st.3 sqhit information
	  }
			
	}
      }


           
      //====Find best reco track in terms of momentum======
 
      float deltaP = 100.;
      float hold_deltaP = 9999.;
            
      SRecTrack* Best_recTrack = NULL;
      n_recTracks = _recEvent->getNTracks();
      //Loop over the SRecTracks
      for(int itrack=0; itrack<_recEvent->getNTracks(); ++itrack){
			
	if (hold_deltaP>deltaP) hold_deltaP = deltaP;

	SRecTrack* recTrack = &_recEvent->getTrack(itrack);

	deltaP = fabs(mom_truth.Mag()- recTrack->getMomentumVertex().Mag());
	//Finding out best match track in terms of energy
	if(deltaP-hold_deltaP<0.) Best_recTrack = recTrack;
			
      }//reco track loop

      //===========================Now fill all info from Best Reco track 

      if (Best_recTrack){
	TVector3 rec_mom = Best_recTrack->getTargetMom();
	cout<<" reco momentum Px,Py,Pz: "<<rec_mom.Px()<<", "<<rec_mom.Py()<<", "<<rec_mom.Pz()<<endl;
	nhits[n_tracks] = Best_recTrack->getNHits();
	charge[n_tracks] = Best_recTrack->getCharge();
	TVector3 rec_vtx = Best_recTrack->getTargetPos();
 
	vx[n_tracks]  = rec_vtx.X();
	vy[n_tracks]  = rec_vtx.Y();
	vz[n_tracks]  = rec_vtx.Z();
    		
	px[n_tracks]  = rec_mom.Px();
	py[n_tracks]  = rec_mom.Py();
	pz[n_tracks]  = rec_mom.Pz();
	pt[n_tracks]  = rec_mom.Pt();
	eta[n_tracks] = rec_mom.Eta();
	phi[n_tracks] = rec_mom.Phi();
    
	chisq_st1[n_tracks] = Best_recTrack->getChisq();
	prob_st1[n_tracks] = Best_recTrack->getProb();                
	quality[n_tracks] = Best_recTrack->getQuality(); 


		
	//st. 1 reco values   
	double tx, ty, tz;
	Best_recTrack->getMomentumSt1(tx, ty, tz);
	px_st1[n_tracks] = tx;
	py_st1[n_tracks] = ty;
	pz_st1[n_tracks] = tz;
	cout<<" In station 1, pz: "<<pz_st1[n_tracks]<<endl;

	double x, y;
	Best_recTrack->getPositionSt1(x, y);		
	x_st1[n_tracks] = x;		
	y_st1[n_tracks] = y;

	// Pull distribution work 
	double state_vec00_st1  = Best_recTrack->getStateVector(0)[0][0];
	double cov00_st1 = Best_recTrack->getCovariance(0)[0][0];
	double true_mom_st1 = sqrt(gpx_st1[n_tracks]*gpx_st1[n_tracks]+gpy_st1[n_tracks]*gpy_st1[n_tracks]+gpz_st1[n_tracks]*gpz_st1[n_tracks]);
	pull_state00[n_tracks] = (fabs(state_vec00_st1) - 1./true_mom_st1)/sqrt(cov00_st1); 
	//work for other state vectors too


	//calculating drift distance for reco track
	if(sqhit_st1 ) {
	
	  double x00 = x_st1[n_tracks]-px_st1[n_tracks]/pz_st1[n_tracks] *sqhit_st1->get_truth_z();
	  double y00 = y_st1[n_tracks]-py_st1[n_tracks]/pz_st1[n_tracks] *sqhit_st1->get_truth_z();

	  rec_drift_st1[n_tracks] = p_geomSvc->getDCA(sqhit_st1->get_detector_id(), sqhit_st1->get_element_id(),  px_st1[n_tracks]/pz_st1[n_tracks], py_st1[n_tracks]/pz_st1[n_tracks], x00, y00);
	  cout<<" x00: "<<x00<<" reco drift: "<<rec_drift_st1[n_tracks]<<endl;

	}    
      }//if best track
	 
      ++n_tracks;
      if(n_tracks>=1000) break;
             
    }//truth loop
	

  }//truth condition
 
  _qa_tree->Fill();


  return Fun4AllReturnCodes::EVENT_OK;
}

//===========================
int AnaTrkQA::End(PHCompositeNode* topNode) {
  if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
    std::cout << "AnaTrkQA::End" << std::endl;

  PHTFileServer::get().cd(_out_name.c_str());

  _tout_reco->Write();
  _qa_tree->Write();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTrkQA::InitEvalTree() {
 
  PHTFileServer::get().open(_out_name.c_str(), "RECREATE");


  _tout_reco = new TTree("Reco", "Reco Eval");

  _tout_reco->Branch("krecstat",      &krecstat,           "krecstat/I");
  _tout_reco->Branch("n_tracks",      &n_tracks,           "n_tracks/I");
  _tout_reco->Branch("par_id",        par_id,              "par_id[n_tracks]/I");
  _tout_reco->Branch("rec_id",        rec_id,              "rec_id[n_tracks]/I");
  _tout_reco->Branch("pid",           pid,                 "pid[n_tracks]/I");
  _tout_reco->Branch("gvx",           gvx,                 "gvx[n_tracks]/F");
  _tout_reco->Branch("gvy",           gvy,                 "gvy[n_tracks]/F");
  _tout_reco->Branch("gvz",           gvz,                 "gvz[n_tracks]/F");
  _tout_reco->Branch("gpx",           gpx,                 "gpx[n_tracks]/F");
  _tout_reco->Branch("gpy",           gpy,                 "gpy[n_tracks]/F");
  _tout_reco->Branch("gpz",           gpz,                 "gpz[n_tracks]/F");
  _tout_reco->Branch("gpt",           gpt,                 "gpt[n_tracks]/F");
  _tout_reco->Branch("geta",          geta,                "geta[n_tracks]/F");
  _tout_reco->Branch("gphi",          gphi,                "gphi[n_tracks]/F");
  _tout_reco->Branch("gx_st1",        gx_st1,              "gx_st1[n_tracks]/F");
  _tout_reco->Branch("gy_st1",        gy_st1,              "gy_st1[n_tracks]/F");
  _tout_reco->Branch("gz_st1",        gz_st1,              "gz_st1[n_tracks]/F");
  _tout_reco->Branch("gpx_st1",       gpx_st1,             "gpx_st1[n_tracks]/F");
  _tout_reco->Branch("gpy_st1",       gpy_st1,             "gpy_st1[n_tracks]/F");
  _tout_reco->Branch("gpz_st1",       gpz_st1,             "gpz_st1[n_tracks]/F");
  _tout_reco->Branch("gnhits",        gnhits,              "gnhits[n_tracks]/I");
  _tout_reco->Branch("gndc",          gndc,                "gndc[n_tracks]/I");
  _tout_reco->Branch("gnhodo",        gnhodo,              "gnhodo[n_tracks]/I");
  _tout_reco->Branch("gnprop",        gnprop,              "gnprop[n_tracks]/I");
  _tout_reco->Branch("gndp",          gndp,                "gndp[n_tracks]/I");
 
  _tout_reco->Branch("sq_pos_st1",      sq_pos_st1,         "sq_pos_st1[n_tracks]/F");
  _tout_reco->Branch("sq_drift_st1",    sq_drift_st1,       "sq_drift_st1[n_tracks]/F");
  _tout_reco->Branch("sq_decID",        sq_decID,           "sq_decID[n_tracks]/F");
  _tout_reco->Branch("rec_drift_st1",   rec_drift_st1,      "rec_drift_st1[n_tracks]/F"); 
  _tout_reco->Branch("sqx_st1",        sqx_st1,              "sqx_st1[n_tracks]/F");
  _tout_reco->Branch("sqy_st1",        sqy_st1,              "sqy_st1[n_tracks]/F");
  _tout_reco->Branch("sqz_st1",        sqz_st1,              "sqz_st1[n_tracks]/F");
  _tout_reco->Branch("chisq_st1",      chisq_st1,            "chisq_st1[n_tracks]/F");
  _tout_reco->Branch("prob_st1",       prob_st1,             "prob_st1[n_tracks]/F");
  _tout_reco->Branch("quality",        quality,             "quality[n_tracks]/F");


  //_tout_reco->Branch("gelmid",        gelmid,              "gelmid[n_tracks][128]/I");

  _tout_reco->Branch("ntruhits",      ntruhits,            "ntruhits[n_tracks]/I");
  _tout_reco->Branch("nhits",         nhits,               "nhits[n_tracks]/I");
  _tout_reco->Branch("charge",        charge,              "charge[n_tracks]/I");
  _tout_reco->Branch("vx",            vx,                  "vx[n_tracks]/F");
  _tout_reco->Branch("vy",            vy,                  "vy[n_tracks]/F");
  _tout_reco->Branch("vz",            vz,                  "vz[n_tracks]/F");
  _tout_reco->Branch("px",            px,                  "px[n_tracks]/F");
  _tout_reco->Branch("py",            py,                  "py[n_tracks]/F");
  _tout_reco->Branch("pz",            pz,                  "pz[n_tracks]/F");
  _tout_reco->Branch("pt",            pt,                  "pt[n_tracks]/F");
  _tout_reco->Branch("eta",           eta,                 "eta[n_tracks]/F");
  _tout_reco->Branch("phi",           phi,                 "phi[n_tracks]/F");
  _tout_reco->Branch("x_st1",         x_st1,               "x_st1[n_tracks]/F");
  _tout_reco->Branch("y_st1",         y_st1,               "y_st1[n_tracks]/F");
  _tout_reco->Branch("px_st1",        px_st1,              "px_st1[n_tracks]/F");
  _tout_reco->Branch("py_st1",        py_st1,              "py_st1[n_tracks]/F");
  _tout_reco->Branch("pz_st1",        pz_st1,              "pz_st1[n_tracks]/F");



  // For the QA tree======================
  _qa_tree = new TTree("QA_ana", "QA analsysis of reconstruction and simulation");

  _qa_tree->Branch("krecstat",      &krecstat,           "krecstat/I");
  _qa_tree->Branch("n_tracks",      &n_tracks,           "n_tracks/I");
  _qa_tree->Branch("n_recTracks",   &n_recTracks,        "n_recTracks/I");
  _qa_tree->Branch("par_id",        par_id,              "par_id[n_tracks]/I");
  _qa_tree->Branch("rec_id",        rec_id,              "rec_id[n_tracks]/I");
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
  _qa_tree->Branch("gx_st1",        gx_st1,              "gx_st1[n_tracks]/F");
  _qa_tree->Branch("gy_st1",        gy_st1,              "gy_st1[n_tracks]/F");
  _qa_tree->Branch("gz_st1",        gz_st1,              "gz_st1[n_tracks]/F");
  _qa_tree->Branch("gpx_st1",       gpx_st1,             "gpx_st1[n_tracks]/F");
  _qa_tree->Branch("gpy_st1",       gpy_st1,             "gpy_st1[n_tracks]/F");
  _qa_tree->Branch("gpz_st1",       gpz_st1,             "gpz_st1[n_tracks]/F");
  _qa_tree->Branch("gnhits",        gnhits,              "gnhits[n_tracks]/I");
  _qa_tree->Branch("gndc",          gndc,                "gndc[n_tracks]/I");
  _qa_tree->Branch("gnhodo",        gnhodo,              "gnhodo[n_tracks]/I");
  _qa_tree->Branch("gnprop",        gnprop,              "gnprop[n_tracks]/I");
  _qa_tree->Branch("gndp",          gndp,                "gndp[n_tracks]/I");
 
  _qa_tree->Branch("sq_pos_st1",      sq_pos_st1,         "sq_pos_st1[n_tracks]/F");
  _qa_tree->Branch("sq_drift_st1",    sq_drift_st1,       "sq_drift_st1[n_tracks]/F");
  _qa_tree->Branch("sq_decID",        sq_decID,           "sq_decID[n_tracks]/F");
  _qa_tree->Branch("rec_drift_st1",   rec_drift_st1,      "rec_drift_st1[n_tracks]/F"); 
  _qa_tree->Branch("sqx_st1",        sqx_st1,              "sqx_st1[n_tracks]/F");
  _qa_tree->Branch("sqy_st1",        sqy_st1,              "sqy_st1[n_tracks]/F");
  _qa_tree->Branch("sqz_st1",        sqz_st1,              "sqz_st1[n_tracks]/F");
  _qa_tree->Branch("chisq_st1",      chisq_st1,            "chisq_st1[n_tracks]/F");
  _qa_tree->Branch("prob_st1",       prob_st1,             "prob_st1[n_tracks]/F");
  _qa_tree->Branch("quality",        quality,             "quality[n_tracks]/F");


  //_qa_tree->Branch("gelmid",        gelmid,              "gelmid[n_tracks][128]/I");

  _qa_tree->Branch("ntruhits",      ntruhits,            "ntruhits[n_tracks]/I");
  _qa_tree->Branch("nhits",         nhits,               "nhits[n_tracks]/I");
  _qa_tree->Branch("charge",        charge,              "charge[n_tracks]/I");
  _qa_tree->Branch("vx",            vx,                  "vx[n_tracks]/F");
  _qa_tree->Branch("vy",            vy,                  "vy[n_tracks]/F");
  _qa_tree->Branch("vz",            vz,                  "vz[n_tracks]/F");
  _qa_tree->Branch("px",            px,                  "px[n_tracks]/F");
  _qa_tree->Branch("py",            py,                  "py[n_tracks]/F");
  _qa_tree->Branch("pz",            pz,                  "pz[n_tracks]/F");
  _qa_tree->Branch("pt",            pt,                  "pt[n_tracks]/F");
  _qa_tree->Branch("eta",           eta,                 "eta[n_tracks]/F");
  _qa_tree->Branch("phi",           phi,                 "phi[n_tracks]/F");
  _qa_tree->Branch("x_st1",         x_st1,               "x_st1[n_tracks]/F");
  _qa_tree->Branch("y_st1",         y_st1,               "y_st1[n_tracks]/F");
  _qa_tree->Branch("px_st1",        px_st1,              "px_st1[n_tracks]/F");
  _qa_tree->Branch("py_st1",        py_st1,              "py_st1[n_tracks]/F");
  _qa_tree->Branch("pz_st1",        pz_st1,              "pz_st1[n_tracks]/F");

  _qa_tree->Branch("pull_state00",   pull_state00,       "pull_state00[n_tracks]/F");

  return 0;
}

int AnaTrkQA::ResetEvalVars() {
  run_id = std::numeric_limits<int>::max();
  spill_id = std::numeric_limits<int>::max();
  target_pos = std::numeric_limits<float>::max();
  event_id = std::numeric_limits<int>::max();
  emu_trigger = 0;
  krecstat = std::numeric_limits<int>::max();

  n_hits = 0;
  for(int i=0; i<10000; ++i) {
    detector_id[i]    = std::numeric_limits<short>::max();
    element_id[i]     = std::numeric_limits<short>::max();
    hodo_mask[i]      = std::numeric_limits<short>::max();
    drift_distance[i] = std::numeric_limits<float>::max();
    pos[i]            = std::numeric_limits<float>::max();
    detector_z[i]     = std::numeric_limits<float>::max();

    truth_x[i]       = std::numeric_limits<float>::max();
    truth_y[i]       = std::numeric_limits<float>::max();
    truth_z[i]       = std::numeric_limits<float>::max();
    truth_pos[i]     = std::numeric_limits<float>::max();
  }

  n_tracks = 0;
  for(int i=0; i<1000; ++i) {
    rec_id[i]     = std::numeric_limits<int>::max();
    par_id[i]     = std::numeric_limits<int>::max();
    pid[i]        = std::numeric_limits<int>::max();
    gvx[i]        = std::numeric_limits<float>::max();
    gvy[i]        = std::numeric_limits<float>::max();
    gvz[i]        = std::numeric_limits<float>::max();
    gpx[i]        = std::numeric_limits<float>::max();
    gpy[i]        = std::numeric_limits<float>::max();
    gpz[i]        = std::numeric_limits<float>::max();
    gpt[i]        = std::numeric_limits<float>::max();
    geta[i]       = std::numeric_limits<float>::max();
    gphi[i]       = std::numeric_limits<float>::max();
    gnhits[i]     = std::numeric_limits<int>::max();
    gx_st1[i]     = std::numeric_limits<float>::max();
    gy_st1[i]     = std::numeric_limits<float>::max();
    gz_st1[i]     = std::numeric_limits<float>::max();
    gpx_st1[i]    = std::numeric_limits<float>::max();
    gpy_st1[i]    = std::numeric_limits<float>::max();
    gpz_st1[i]    = std::numeric_limits<float>::max();
    gndc[i]       = std::numeric_limits<int>::max();
    gnhodo[i]     = std::numeric_limits<int>::max();
    gnprop[i]     = std::numeric_limits<int>::max();
    gndp[i]       = std::numeric_limits<int>::max();

    for(int j=0; j<NDET+1; ++j) {
      gelmid[i][j] = std::numeric_limits<int>::max();
    }

    ntruhits[i]   = std::numeric_limits<int>::max();
    nhits[i]      = std::numeric_limits<int>::max();
    charge[i]     = std::numeric_limits<int>::max();
    vx[i]         = std::numeric_limits<float>::max();
    vy[i]         = std::numeric_limits<float>::max();
    vz[i]         = std::numeric_limits<float>::max();
    px[i]         = std::numeric_limits<float>::max();
    py[i]         = std::numeric_limits<float>::max();
    pz[i]         = std::numeric_limits<float>::max();
    pt[i]         = std::numeric_limits<float>::max();
    eta[i]        = std::numeric_limits<float>::max();
    phi[i]        = std::numeric_limits<float>::max();
    x_st1[i]     = std::numeric_limits<float>::max();
    y_st1[i]     = std::numeric_limits<float>::max();
    px_st1[i]     = std::numeric_limits<float>::max();
    py_st1[i]     = std::numeric_limits<float>::max();
    pz_st1[i]     = std::numeric_limits<float>::max();
  }

  return 0;
}

int AnaTrkQA::GetNodes(PHCompositeNode* topNode) {

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


  //Abi Add
 
  g4hc_d1x  = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D1X");
  g4hc_d3px = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D3pXp");
  g4hc_d3mx = findNode::getClass<PHG4HitContainer      >(topNode, "G4HIT_D3mXp");
  if (! g4hc_d1x) g4hc_d1x = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D0X");

  if ( !g4hc_d1x || !g4hc_d3px || !g4hc_d3mx) {
    cout << "Failed at getting nodes: "<< g4hc_d1x << " " << g4hc_d3px << " " << g4hc_d3mx << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }


  return Fun4AllReturnCodes::EVENT_OK;
}


//Adapted from Kenichi's truth node maker's function
bool AnaTrkQA::FindHitAtStation(const int trk_id, const PHG4HitContainer* g4hc, TVector3* pos, TLorentzVector* mom)
{
  const double M_MU = 0.1056583745; // GeV
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


