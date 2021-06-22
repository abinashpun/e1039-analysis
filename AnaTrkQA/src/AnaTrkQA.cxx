/**
 * \class AnaTrkQA
 * \ module for track quality assurance
 * \author Abinash Pun
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
#include <interface_main/SQDimuonVector_v1.h>

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
#include <TVector3.h>
#include <TLorentzVector.h>

#include <cstring>
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <limits>
#include <tuple>

#include <boost/lexical_cast.hpp>

//#define NDET 62
//#define LogDebug(exp)		std::cout<<"DEBUG: "  <<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
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
  legacyContainer(true),
  _out_name("eval.root")
{
 
}

int AnaTrkQA::Init(PHCompositeNode* topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTrkQA::InitRun(PHCompositeNode* topNode) {

  event_id = 0;

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
    //ret = RecoEval(topNode);
    ret = TruthRecoEval(topNode);
    if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  }

  DimuonInfo(topNode);

  ++event_id;

  return ret;
}


//play ground for Abi==============================================
int AnaTrkQA::TruthRecoEval(PHCompositeNode* topNode)
{
    ResetEvalVars();

    //std::cout << "event ID = " << event_id << std::endl;

  if(_truth) {
    for(auto iter=_truth->GetPrimaryParticleRange().first;
	iter!=_truth->GetPrimaryParticleRange().second;
	++iter) {
      PHG4Particle * par = iter->second; 
     
      int trk_id = par->get_track_id();

      pid = par->get_pid();

// << "pid = " << pid << std::endl;

      int vtx_id =  par->get_vtx_id();

      // generated momentum and postion at target
      PHG4VtxPoint* vtx = _truth->GetVtx(vtx_id);

      gpos->SetXYZ(vtx->get_x(), vtx->get_y(), vtx->get_z());
      gmom->SetXYZ(par->get_px(), par->get_py(), par->get_pz());
       
 /// G4Hits at different stations----  (G4Hit information are not stored in tree for now. Can be extended to use later on)   
      TVector3 g_pos_st1;
      TLorentzVector g_mom_st1;
      bool st1hit = FindG4HitAtStation(trk_id, g4hc_d1x, &g_pos_st1, &g_mom_st1);
      if(st1hit){
	
	gx_st1  = g_pos_st1.X();
	gy_st1  = g_pos_st1.Y();
	gz_st1  = g_pos_st1.Z();
     		
	gpx_st1 = g_mom_st1.Px();
	gpy_st1 = g_mom_st1.Py();
	gpz_st1 = g_mom_st1.Pz();
      }

     TVector3 g_pos_st2;
     TLorentzVector g_mom_st2;
     bool st2hit =  FindG4HitAtStation(trk_id, g4hc_d2xp, &g_pos_st2, &g_mom_st2);
    
     TVector3 g_pos_st3;
     TLorentzVector g_mom_st3;

      bool st3hit = FindG4HitAtStation(trk_id, g4hc_d3px, &g_pos_st3, &g_mom_st3)|| FindG4HitAtStation(trk_id, g4hc_d3mx, &g_pos_st3, &g_mom_st3);
		        
     bool prophit = (FindG4HitAtProp(trk_id,g4hc_p1x1)||FindG4HitAtProp(trk_id,g4hc_p1x2)||FindG4HitAtProp(trk_id,g4hc_p1y1)||FindG4HitAtProp(trk_id,g4hc_p1y2)) && (FindG4HitAtProp(trk_id,g4hc_p2x1)||FindG4HitAtProp(trk_id,g4hc_p2x2)||FindG4HitAtProp(trk_id,g4hc_p2y1)||FindG4HitAtProp(trk_id,g4hc_p2y2));//truth having prop hits at station 4

/// Detector acceptance: Truth particle passing through all the 4 stations------------------
    if(st1hit && st2hit && st3hit && prophit){
        ac_mom = gmom;
    }


///==========Implementing functions to catch best track===
      SRecTrack* Best_recTrack = NULL;
      n_recTracks = _recEvent->getNTracks();
	
     // if(_recEvent->getNTracks()>0) Best_recTrack = FindBestMomRecTrack(_recEvent, mom_truth.Mag());      
  
///@Hit matching condition for choosing reconstructed track ---------------   
       vector<int> sqhit_idvec;
       map<int, vector<int> > rtrkid_hitidvec;

	///fill sqhit hit_id vector    
	 for(int ihit=0; ihit<_hit_vector->size(); ++ihit) {
	 	SQHit *sqhit = _hit_vector->at(ihit);
	   	int sq_detid = sqhit->get_detector_id();
           	if(sq_detid > nChamberPlanes || (sq_detid >= 7 && sq_detid <= 12)) continue;
	   	if(trk_id != sqhit->get_track_id()) continue;
		sqhit_idvec.push_back(sqhit->get_hit_id());
				
	  }
	
         sort(sqhit_idvec.begin(), sqhit_idvec.end());

	///fill map of reco track id and hitindex
	 if(n_recTracks>0){
	 	for(int i = 0; i < n_recTracks; ++i) {
      		 	SRecTrack* recTrack = &_recEvent->getTrack(i);                      
      		 	rtrkid_hitidvec[i] = vector<int>();
      
       			int n_rhits = recTrack->getNHits();
      			for(int j = 0; j < n_rhits; ++j) { 
        		rtrkid_hitidvec[i].push_back(fabs(recTrack->getHitIndex(j)));
			
      			}

      		sort(rtrkid_hitidvec[i].begin(), rtrkid_hitidvec[i].end());
    		}
         }
	
///Now try to find the matching reco id
	double m_matching_threshold = 0.75;
	int rtrkid = -1;
     	unsigned int n_match = 0;
	
	for(auto it = rtrkid_hitidvec.begin(); it != rtrkid_hitidvec.end(); ++it) {
      		int n_match_new = FindCommonHitIDs(sqhit_idvec, it->second);
      		if(n_match < n_match_new) {
        		n_match = n_match_new;
        		rtrkid = it->first;
      		}
    	}

    	if(rtrkid >= 0 && double(n_match)/double(sqhit_idvec.size()) > m_matching_threshold) {
      	Best_recTrack = &_recEvent->getTrack(rtrkid);
    	}
///@Hitmatching...ends...........................................


      if(Best_recTrack){

    // reconstructed momentum and postion @ target
	*rec_mom = Best_recTrack->getTargetMom();
	*rec_pos = Best_recTrack->getTargetPos();

	nhits = Best_recTrack->getNHits();
	charge = Best_recTrack->getCharge();

	nhits_st1 = Best_recTrack->getNHitsInStation(1);
	nhits_st2 = Best_recTrack->getNHitsInStation(2);
    nhits_st3 = Best_recTrack->getNHitsInStation(3);
	    
	chisq_st1 = Best_recTrack->getChisq();
	prob_st1 = Best_recTrack->getProb();
	quality = Best_recTrack->getQuality();
		
	//Old way of getting st. 1 reco values   
	/*double tx, ty, tz;
	Best_recTrack->getMomentumSt1(tx, ty, tz);
	px_st1[n_tracks] = tx;
	py_st1[n_tracks] = ty;
	pz_st1[n_tracks] = tz;
		
	double x, y;
	Best_recTrack->getPositionSt1(x, y);		
        x_st1[n_tracks] = x;		
	y_st1[n_tracks] = y;
        */
      }//if best reco track



      ///=================Digitized hit information at different stations and corresponding reco values (if any)

      if(_hit_vector) {
	for(int ihit=0; ihit<_hit_vector->size(); ++ihit) {
	  SQHit *hit = _hit_vector->at(ihit);
           
           int detid = hit->get_detector_id();      
           if(detid > nChamberPlanes || (detid >= 7 && detid <= 12)) continue; //only for the chamber hits
	  
	  /// D0X is considered as st1
	  if(hit->get_track_id() == trk_id && hit->get_detector_id() ==3 ) {

	    // truth info station 1
	    sq_mom_st1->SetXYZ(hit->get_truth_px(), hit->get_truth_py(), hit->get_truth_pz());
	    sq_pos_st1->SetXYZ(hit->get_truth_x(), hit->get_truth_y(), hit->get_truth_z());
	    sq_drift_st1 = hit->get_drift_distance();
 
	    ///if the best reco track available
	    if (Best_recTrack){						
	      double sq_z_st1 = hit->get_truth_z();
	      int rec_index_st1 = Best_recTrack->getNearestNode(sq_z_st1);
			
	      double z_rec_st1 = Best_recTrack->getZ(rec_index_st1);
						
	      if(fabs(sq_z_st1- z_rec_st1>1.)) continue;//to avoid mismatch of nodes
			
	      double p_rec_st1 = fabs(1./Best_recTrack->getStateVector(rec_index_st1)[0][0]);
	      double tx_rec_st1 = Best_recTrack->getStateVector(rec_index_st1)[1][0];
	      double ty_rec_st1 = Best_recTrack->getStateVector(rec_index_st1)[2][0];
	      double x_rec_st1 = Best_recTrack->getStateVector(rec_index_st1)[3][0];
	      double y_rec_st1 = Best_recTrack->getStateVector(rec_index_st1)[4][0];

	      double x0_st1 = x_rec_st1 - tx_rec_st1 *z_rec_st1;
	      double y0_st1 = y_rec_st1 - ty_rec_st1 *z_rec_st1;

	      double rec_p_st1 =  p_rec_st1;
	      double rec_pz_st1 = p_rec_st1/sqrt(1.+tx_rec_st1*tx_rec_st1+ty_rec_st1*ty_rec_st1);
	      double rec_px_st1 = rec_pz_st1* tx_rec_st1;
	      double rec_py_st1 = rec_pz_st1* ty_rec_st1;

	      // reco. info station 1
	      rec_pos_st1->SetXYZ(x_rec_st1, y_rec_st1, z_rec_st1);
	      rec_mom_st1->SetXYZ(rec_px_st1, rec_py_st1, rec_pz_st1);
	      rec_drift_st1 = p_geomSvc->getDCA(hit->get_detector_id(), hit->get_element_id(),tx_rec_st1, ty_rec_st1, x0_st1,y0_st1);

	      double px = hit->get_truth_px();
	      double py = hit->get_truth_py();
	      double pz = hit->get_truth_pz();

	      // Pull distribution work
	      double cov00_st1 = Best_recTrack->getCovariance(rec_index_st1)[0][0];
	      double mom_st1 = sqrt(px* px + py* py + pz* pz);
	      pull_q2p_st1 = (fabs(Best_recTrack->getStateVector(rec_index_st1)[0][0]) - 1./mom_st1)/sqrt(cov00_st1);
		
	  
	    }///best reco condition

	  }///st.1 work ends

	  ///==========

	  ///st. 2 is now D2Xp
	  if(hit->get_track_id() == trk_id && hit->get_detector_id() ==15 ) {

	    // truth info  station 2
	    sq_pos_st2->SetXYZ(hit->get_truth_x(), hit->get_truth_y(), hit->get_truth_z());
	    sq_mom_st2->SetXYZ(hit->get_truth_px(), hit->get_truth_py(), hit->get_truth_pz());
	    sq_drift_st2 = hit->get_drift_distance();
 
	    ///if the best reco track available
	    if (Best_recTrack){						
	      double sq_z = hit->get_truth_z();
	      int rec_index = Best_recTrack->getNearestNode(sq_z);
			
	      double rec_z = Best_recTrack->getZ(rec_index);
						
	      if(fabs(sq_z- rec_z>1.)) continue;///to avid mismatch of node
			
	      double p_rec = fabs(1./Best_recTrack->getStateVector(rec_index)[0][0]);
	      double tx_rec = Best_recTrack->getStateVector(rec_index)[1][0];
	      double ty_rec = Best_recTrack->getStateVector(rec_index)[2][0];
	      double x_rec = Best_recTrack->getStateVector(rec_index)[3][0];
	      double y_rec = Best_recTrack->getStateVector(rec_index)[4][0];

	      double x0 = x_rec - tx_rec *rec_z;
	      double y0 = y_rec - ty_rec *rec_z;

	      double rec_p_st2 =  p_rec;
	      double rec_pz_st2 = p_rec/sqrt(1.+tx_rec*tx_rec+ty_rec*ty_rec);
	      double rec_px_st2 = rec_pz_st2* tx_rec;
	      double rec_py_st2 = rec_pz_st2* ty_rec;

	      // reco. info station 2
	      rec_pos_st2->SetXYZ(x_rec, y_rec, rec_z);
	      rec_mom_st2->SetXYZ(rec_px_st2, rec_py_st2, rec_pz_st2);
  	      rec_drift_st2 = p_geomSvc->getDCA(hit->get_detector_id(), hit->get_element_id(),tx_rec, ty_rec, x0,y0);

  	      double px = hit->get_truth_px();
	      double py = hit->get_truth_py();
	      double pz = hit->get_truth_pz();

	      ///Pull distribution work 
	      double cov00_st2 = Best_recTrack->getCovariance(rec_index)[0][0];
	      double mom_st2 = sqrt(px* px + py* py + pz* pz);
	      pull_q2p_st2 = (fabs(Best_recTrack->getStateVector(rec_index)[0][0]) - 1./mom_st2)/sqrt(cov00_st2);
	    }///if best reco track

	  }//st2. work done

	  //=======================
	  			 	

	  //st. 3 is now D3mXp(id = 27) or D3pXp (id=21)
	  if(hit->get_track_id() == trk_id && (hit->get_detector_id() == 27 ||hit->get_detector_id() == 21)) {          		

	  // truth info station 3
	  sq_pos_st3->SetXYZ(hit->get_truth_x(), hit->get_truth_y(), hit->get_truth_z());
	  sq_mom_st3->SetXYZ(hit->get_truth_px(), hit->get_truth_py(), hit->get_truth_pz());
	  sq_drift_st3 = hit->get_drift_distance();
 
	  ///if the best reco track available
	  if (Best_recTrack){						
	  double sq_z = hit->get_truth_z();
	  int rec_index = Best_recTrack->getNearestNode(sq_z);
			
	  double rec_z = Best_recTrack->getZ(rec_index);
						
	  if(fabs(sq_z- rec_z>1.)) continue;///to avoid mismatch of nodes
			
	  double p_rec = fabs(1./Best_recTrack->getStateVector(rec_index)[0][0]);
	  double tx_rec = Best_recTrack->getStateVector(rec_index)[1][0];
	  double ty_rec = Best_recTrack->getStateVector(rec_index)[2][0];
	  double x_rec = Best_recTrack->getStateVector(rec_index)[3][0];
	  double y_rec = Best_recTrack->getStateVector(rec_index)[4][0];

	  double x0 = x_rec - tx_rec *rec_z;
	  double y0 = y_rec - ty_rec *rec_z;

	  double rec_p_st3 =  p_rec;
	  double rec_pz_st3 = p_rec/sqrt(1.+tx_rec*tx_rec+ty_rec*ty_rec);
	  double rec_px_st3 = rec_pz_st3* tx_rec;
	  double rec_py_st3 = rec_pz_st3* ty_rec;

      // reco. info station 3
      rec_mom_st3->SetXYZ(rec_px_st3, rec_py_st3, rec_pz_st3);
      rec_pos_st3->SetXYZ(x_rec, y_rec, rec_z);
	  rec_drift_st3 = p_geomSvc->getDCA(hit->get_detector_id(), hit->get_element_id(),tx_rec, ty_rec, x0,y0);

	  double px = hit->get_truth_px();
	  double py = hit->get_truth_py();
	  double pz = hit->get_truth_pz();

	  // Pull distribution work 
	  double cov00_st3 = Best_recTrack->getCovariance(rec_index)[0][0];
	  double mom_st3 = sqrt(px* px + py* py + pz* pz);
	  pull_q2p_st3 = (fabs(Best_recTrack->getStateVector(rec_index)[0][0]) - 1./mom_st3)/sqrt(cov00_st3);

	  //std::cout << "pull_q2p_st3 = " << pull_q2p_st3 << std::endl;

	  }//if best reco track
	  }//st3. work done
	  
	  //=======================		
			
	}//sqhit vector loop

      }//if hit vector

      tree1->Fill();

      ++n_tracks;
      if(n_tracks>=100) break;
             
    }//truth loop
  }//truth condition

  return Fun4AllReturnCodes::EVENT_OK;
}

//dimuon info
int AnaTrkQA:: DimuonInfo(PHCompositeNode* topNode){

  int nDimuons = dimuonVector->size();
  int nRecDimuons =  legacyContainer ? _recEvent->getNDimuons() : (recDimuonVector ? recDimuonVector->size() : -1);

  for(int i = 0; i < nDimuons; ++i)
  {
    SQDimuon* dimuon = dimuonVector->at(i);
    mass = dimuon->get_mom().M();
    *vtx = dimuon->get_pos();
    *pmom = dimuon->get_mom_pos().Vect();
    *nmom = dimuon->get_mom_neg().Vect();
    //xsec = mcEvent->get_cross_section();
    if(AllChamberPlaneHits(dimuon->get_track_id_pos(), _hit_vector) && AllChamberPlaneHits(dimuon->get_track_id_neg(), _hit_vector)) mass_acc =  mass;
    //if(!AllChamberPlaneHits(dimuon->get_track_id_pos(), _hit_vector)) Fun4AllReturnCodes::DISCARDEVENT;
    //if(!AllChamberPlaneHits(dimuon->get_track_id_neg(), _hit_vector)) Fun4AllReturnCodes::DISCARDEVENT; 
    //std::cout << "mass : " << mass << std::endl;

    int recid = dimuon->get_rec_dimuon_id();

    if(recid >= 0 && recid < nRecDimuons)
    {
      SRecDimuon* recDimuon = legacyContainer ? &(_recEvent->getDimuon(recid)) : dynamic_cast<SRecDimuon*>(recDimuonVector->at(recid));
      rec_mass = recDimuon->mass;
      *rec_pmom = recDimuon->p_pos.Vect();
      *rec_nmom = recDimuon->p_neg.Vect();
      *rec_ppos = recDimuon->vtx_pos;
      *rec_npos = recDimuon->vtx_neg;
      *rec_vtx  = recDimuon->vtx;
    }
///@temporary modification by Abi to work with lecayContainer
    if(recid <0 && nDimuons==nRecDimuons)
    {
      SRecDimuon* recDimuon = &(_recEvent->getDimuon(0));
      rec_mass = recDimuon->mass;
      *rec_pmom = recDimuon->p_pos.Vect();
      *rec_nmom = recDimuon->p_neg.Vect();
      *rec_ppos = recDimuon->vtx_pos;
      *rec_npos = recDimuon->vtx_neg;
      *rec_vtx  = recDimuon->vtx;
    }
    tree2->Fill();
    }

return Fun4AllReturnCodes::EVENT_OK;
}

///===========================
int AnaTrkQA::End(PHCompositeNode* topNode) {
  if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
    std::cout << "AnaTrkQA::End" << std::endl;

  PHTFileServer::get().cd(_out_name.c_str());
  tree1->Write();
  tree2->Write();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTrkQA::InitEvalTree() {
 
  PHTFileServer::get().open(_out_name.c_str(), "RECREATE");

  ///For track tree ======================
  tree1 = new TTree("track", "track info");

  tree1->Branch("event_id",      &event_id,           "event_id/I");
  tree1 -> Branch("pid",           &pid,                 "pid/I");

  ///Generated Truth info
  tree1->Branch("gpos",          "TVector3",       &gpos);
  tree1->Branch("gmom",          "TVector3",       &gmom);

  ///Accepted truth info
  tree1->Branch("ac_mom",       "TVector3",         &ac_mom);

  ///Reco info in vertex
  tree1->Branch("rec_pos",           "TVector3",           &rec_pos);
  tree1->Branch("rec_mom",           "TVector3",           &rec_mom);

 
  ///station 1 truth and reco info
  tree1->Branch("sq_mom_st1",           "TVector3",           &sq_mom_st1);
  tree1->Branch("sq_pos_st1",           "TVector3",           &sq_pos_st1);
  tree1->Branch("sq_drift_st1",         &sq_drift_st1,        "sq_drift_st1/F");

  tree1->Branch("rec_mom_st1",           "TVector3",           &rec_mom_st1);
  tree1->Branch("rec_pos_st1",           "TVector3",           &rec_pos_st1);
  tree1->Branch("rec_drift_st1",         &rec_drift_st1,       "rec_drift_st1/F");

 
 /// Station 2 truth and reco info
  tree1->Branch("sq_mom_st2",           "TVector3",           &sq_mom_st2);
  tree1->Branch("sq_pos_st2",           "TVector3",           &sq_pos_st2);
  tree1->Branch("sq_drift_st2",         &sq_drift_st2,        "sq_drift_st2/F");

  tree1->Branch("rec_mom_st2",           "TVector3",           &rec_mom_st2);
  tree1->Branch("rec_pos_st2",           "TVector3",           &rec_pos_st2);
  tree1->Branch("rec_drift_st2",         &rec_drift_st2,       "rec_drift_st2/F");


/// Station 3 truth and reco info
  tree1->Branch("sq_mom_st3",           "TVector3",           &sq_mom_st3);
  tree1->Branch("sq_pos_st3",           "TVector3",           &sq_pos_st3);
  tree1->Branch("sq_drift_st3",         &sq_drift_st3,        "sq_drift_st3/F");

  tree1->Branch("rec_mom_st3",           "TVector3",           &rec_mom_st3);
  tree1->Branch("rec_pos_st3",           "TVector3",           &rec_pos_st3);
  tree1->Branch("rec_drift_st3",         &rec_drift_st3,       "rec_drift_st3/F");


///quality info of reconstructed tracks
  tree1->Branch("pull_q2p_st1",      &pull_q2p_st1,      "pull_q2p_st1/F");
  tree1->Branch("pull_q2p_st2",      &pull_q2p_st2,      "pull_q2p_st2/F");
  tree1->Branch("pull_q2p_st3",      &pull_q2p_st3,      "pull_q2p_st3/F");

  tree1->Branch("chisq_st1",      &chisq_st1,            "chisq_st1/F");
  tree1->Branch("prob_st1",       &prob_st1,             "prob_st1/F");
  tree1->Branch("quality",        &quality,             "quality/F");

  tree1->Branch("nhits",         &nhits,               "nhits/I");
  tree1->Branch("nhits_st1",     &nhits_st1,           "nhits_st1/I");
  tree1->Branch("nhits_st2",     &nhits_st2,           "nhits_st2/I");
  tree1->Branch("nhits_st3",     &nhits_st3,           "nhits_st3/I");
  tree1->Branch("charge",        &charge,              "charge/I");

  // dimuon tree
  tree2 = new TTree("dimuon", "dimuon info");
  tree2 -> Branch("event_id", &event_id, "event_id/I");
  tree2 -> Branch("n_tracks",      &n_tracks,           "n_tracks/I");
  tree2 -> Branch("n_recTracks",   &n_recTracks,        "n_recTracks/I");
  tree2->Branch("mass", &mass, "mass/D");
  tree2->Branch("mass_acc", &mass_acc, "mass_acc/D");
  tree2->Branch("rec_mass", &rec_mass, "rec_mass/D");
  tree2->Branch("vtx", &vtx, 256000, 99);
  tree2->Branch("pmom", &pmom, 256000, 99);
  tree2->Branch("nmom", &nmom, 256000, 99);
  tree2->Branch("rec_pmom", &rec_pmom, 256000, 99);
  tree2->Branch("rec_nmom", &rec_nmom, 256000, 99);
  tree2->Branch("rec_ppos", &rec_ppos, 256000, 99);
  tree2->Branch("rec_npos", &rec_npos, 256000, 99);
  tree2->Branch("rec_vtx", &rec_vtx, 256000, 99);

  return 0;
}

int AnaTrkQA::ResetEvalVars() {
  run_id = -9999;
  spill_id = -9999;
  target_pos = -9999.;
  //event_id = std::numeric_limits<int>::max();;
  emu_trigger = 0;
  krecstat = -9999;

    n_hits = 0;
    detector_id    = 9999;
    element_id     = 9999;
    hodo_mask      = 9999;
    drift_distance = 9999;
    pos            = 9999;
    detector_z     = 9999;

    truth_x       = 9999.;
    truth_y       = 9999.;
    truth_z       = 9999.;
    truth_pos     = 9999.;

    n_tracks = 0;

    rec_id    = 9999;
    par_id     = 9999;
    pid       = 9999;
    gnhits     = 9999;
    gx_st1     = 9999.;
    gy_st1     = 9999.;
    gz_st1     = 9999.;
    gpx_st1    = 9999.;
    gpy_st1    = 9999.;
    gpz_st1    = 9999.;
    gndc       = 9999;
    gnhodo     = 9999;
    gnprop     = 9999;
    gndp       = 9999;



    /*for(int j=0; j<NDET+1; ++j) {
      gelmid[i][j] = std::numeric_limits<int>::max();;
    }
*/

    nhits      = 9999;
    charge     = 9999;
  /*  x_st1[i]     = std::numeric_limits<float>::max();
    y_st1[i]     = std::numeric_limits<float>::max();
    px_st1[i]     = std::numeric_limits<float>::max();
    py_st1[i]     = std::numeric_limits<float>::max();
    pz_st1[i]     = std::numeric_limits<float>::max();
*/

    sq_drift_st1     = 9999.;
    rec_drift_st1     = 9999.;
    sq_drift_st2     = 9999.;
    rec_drift_st2     = 9999.;
    sq_drift_st3     = 9999.;
    rec_drift_st3     = 9999.;

    nhits_st1 = 9999.;
    nhits_st2 = 9999.;
    nhits_st3 = 9999.;


    rec_mom_st3 = new TVector3(9999., 9999., 9999.);
    rec_pos_st3 = new TVector3(9999., 9999., 9999.);
    sq_mom_st3 = new TVector3(9999., 9999., 9999.);
    sq_pos_st3 = new TVector3(9999., 9999., 9999.);
    rec_mom_st2 = new TVector3(9999., 9999., 9999.);
    rec_pos_st2 = new TVector3(9999., 9999., 9999.);
    sq_mom_st2 = new TVector3(9999., 9999., 9999.);
    sq_pos_st2 = new TVector3(9999., 9999., 9999.);
    rec_pos_st1 = new TVector3(9999., 9999., 9999.);
    rec_mom_st1 = new TVector3(9999., 9999., 9999.);
    sq_mom_st1 = new TVector3(9999., 9999., 9999.);
    sq_pos_st1 = new TVector3(9999., 9999., 9999.);
    rec_pos = new TVector3(9999., 9999., 9999.);
    rec_mom = new TVector3(9999., 9999., 9999.);
    ac_mom = new TVector3(9999., 9999., 9999.);
    gpos = new TVector3(9999., 9999., 9999.);
    gmom = new TVector3(9999., 9999., 9999.);

    // dimuon info
    mass = -9999.;
    mass_acc = -9999.;
    vtx = new TVector3(-9999., -9999., -9999.);
    pmom = new TVector3(-9999., -9999., -9999.);
    nmom = new TVector3(-9999., -9999., -9999.);

    rec_mass = -9999.;
    rec_pmom = new TVector3(-9999., -9999., -9999.);
    rec_nmom = new TVector3(-9999., -9999., -9999.);
    rec_ppos = new TVector3(-9999., -9999., -9999.);
    rec_npos = new TVector3(-9999., -9999., -9999.);
    rec_vtx = new TVector3(-9999., -9999., -9999.);

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

 // dimuon info
 dimuonVector = findNode::getClass<SQDimuonVector>(topNode, "SQTruthDimuonVector");
 recDimuonVector = findNode::getClass<SQDimuonVector>(topNode, "SQRecDimuonVector");
 if(!dimuonVector)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
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


//For finding g4hit information in stations (following Kenichi's truth node maker)
bool AnaTrkQA::FindG4HitAtStation(const int trk_id, const PHG4HitContainer* g4hc, TVector3* pos, TLorentzVector* mom)
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
SRecTrack* AnaTrkQA::FindBestMomRecTrack(SRecEvent *recEvent,  const float true_TargetP)
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
int AnaTrkQA::FindCommonHitIDs(vector<int>& hitidvec1, vector<int>& hitidvec2)
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
bool AnaTrkQA::FindG4HitAtHodo(const int trk_id, const PHG4HitContainer* g4hc)
{
  //const double M_MU = 0.1056583745; 
  PHG4HitContainer::ConstRange range = g4hc->getHits();
  for (PHG4HitContainer::ConstIterator it = range.first; it != range.second; it++) {
    PHG4Hit* hit = it->second;
    if (hit->get_trkid() == trk_id) {
      return true;
    }
  }
  return false;
}


bool AnaTrkQA::FindG4HitAtProp(const int trk_id, const PHG4HitContainer* g4hc)
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

bool AnaTrkQA::AllChamberPlaneHits(const int trk_id, SQHitVector *_hit_vector )
{
  bool D0Xhit, D0Xphit, D0Vhit, D0Vphit, D0Uhit, D0Uphit, D2Xhit, D2Xphit, D2Uhit, D2Uphit, D2Vhit, D2Vphit, D3pXhit, D3pXphit, D3pUhit, D3pUphit, D3pVhit, D3pVphit, D3mXhit, D3mXphit, D3mUhit, D3mUphit, D3mVhit, D3mVphit;

  D0Xhit = D0Xphit = D0Vhit =D0Vphit= D0Uhit= D0Uphit= D2Xhit= D2Xphit= D2Uhit= D2Uphit= D2Vhit= D2Vphit= D3pXhit= D3pXphit= D3pUhit= D3pUphit= D3pVhit= D3pVphit= D3mXhit= D3mXphit= D3mUhit= D3mUphit= D3mVhit= D3mVphit = false;
 

        for(int ihit=0; ihit<_hit_vector->size(); ++ihit) {
          SQHit *sqhit = _hit_vector->at(ihit);
            if(sqhit->get_track_id() == trk_id){
              if(sqhit->get_detector_id() == 1 ) D0Uhit=true;
              if(sqhit->get_detector_id() == 2 ) D0Uphit=true;
              if(sqhit->get_detector_id() == 3 ) D0Xphit=true;
              if(sqhit->get_detector_id() == 4 ) D0Xhit=true;
              if(sqhit->get_detector_id() == 5 ) D0Vhit=true;
              if(sqhit->get_detector_id() == 6 ) D0Vphit=true;

              if(sqhit->get_detector_id() == 13) D2Vhit=true;
              if(sqhit->get_detector_id() == 14) D2Vphit=true;
              if(sqhit->get_detector_id() == 15) D2Xphit=true;
              if(sqhit->get_detector_id() == 16) D2Xhit=true;
              if(sqhit->get_detector_id() == 17) D2Uhit=true;
              if(sqhit->get_detector_id() == 18) D2Uphit=true;

              if(sqhit->get_detector_id() == 19) D3pVphit=true;
              if(sqhit->get_detector_id() == 20) D3pVhit=true;
              if(sqhit->get_detector_id() == 21) D3pXphit=true;
              if(sqhit->get_detector_id() == 22) D3pXhit=true;
              if(sqhit->get_detector_id() == 23) D3pUphit=true;
              if(sqhit->get_detector_id() == 24) D3pUhit=true;

              if(sqhit->get_detector_id() == 25) D3mVphit=true;
              if(sqhit->get_detector_id() == 26) D3mVhit=true;
              if(sqhit->get_detector_id() == 27) D3mXphit=true;
              if(sqhit->get_detector_id() == 28) D3mXhit=true;
              if(sqhit->get_detector_id() == 29) D3mUphit=true;
              if(sqhit->get_detector_id() == 30) D3mUhit=true;
        }    
      }

 st1hit = false;
 st2hit = false;
 st3phit = false;
 st3mhit = false;
 st3hit = false;

 st1hit = D0Uhit && D0Uphit && D0Xphit && D0Xhit && D0Vhit && D0Vphit;
 st2hit = D2Uhit && D2Uphit && D2Xphit && D2Xhit && D2Vhit && D2Vphit;
 st3phit = D3pUhit && D3pUphit && D3pXphit && D3pXhit && D3pVhit && D3pVphit;
 st3mhit = D3mUhit && D3mUphit && D3mXphit && D3mXhit && D3mVhit && D3mVphit;
 st3hit = st3phit || st3mhit;

  if (st3phit && st3mhit) return false;
  else if(st1hit && st2hit && st3hit) return true;
  else return false;
}

