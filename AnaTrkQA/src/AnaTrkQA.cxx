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
    //ret = RecoEval(topNode);
    ret = RecoEvalv2(topNode);
    if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  }

  ++_event;

  return ret;
}


//play ground for Abi==============================================
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


      //st1. true g4 hit info (with the D0X) [just filled the tree but sqhit info is used. Can be removed later]
      TVector3 g_pos_st1;
      TLorentzVector g_mom_st1;
      if(FindHitAtStation(trk_id, g4hc_d1x, &g_pos_st1, &g_mom_st1)){;
	
	gx_st1[n_tracks]  = g_pos_st1.X();
	gy_st1[n_tracks]  = g_pos_st1.Y();
	gz_st1[n_tracks]  = g_pos_st1.Z();
     		
	gpx_st1[n_tracks] = g_mom_st1.Px();
	gpy_st1[n_tracks] = g_mom_st1.Py();
	gpz_st1[n_tracks] = g_mom_st1.Pz();
      }

      //st3. trute hit info [can be filled to tree but for now sqhit infor is used]
      TVector3 g_pos_st3;
      TLorentzVector g_mom_st3;
      if (FindHitAtStation(trk_id, g4hc_d3px, &g_pos_st3, &g_mom_st3)|| FindHitAtStation(trk_id, g4hc_d3mx, &g_pos_st3, &g_mom_st3)){
		  //work in progress...(for now sqhit is   )        
      } 

		
      //==========Implementing function to catch best track in terms of momemtum===
      SRecTrack* Best_recTrack = NULL;
      n_recTracks[n_tracks]= _recEvent->getNTracks();
	
      if(_recEvent->getNTracks()>0) Best_recTrack = FindBestRecTrack(_recEvent, mom_truth.Mag());

      if(Best_recTrack){

	TVector3 rec_mom = Best_recTrack->getTargetMom();
		
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


		
	//st. 1 reco values [Just filled in tree but not used..other "rec_*" branches are filled inside sqhit vector loop   
	double tx, ty, tz;
	Best_recTrack->getMomentumSt1(tx, ty, tz);
	px_st1[n_tracks] = tx;
	py_st1[n_tracks] = ty;
	pz_st1[n_tracks] = tz;
		
	double x, y;
	Best_recTrack->getPositionSt1(x, y);		
	x_st1[n_tracks] = x;		
	y_st1[n_tracks] = y;

      }//if best track


      //=================SQ Hit information
      SQHit *sqhit_st1 = NULL;
      SQHit *sqhit_st3 = NULL;
      // Get the sqhit position at st1. , st3 and drift distances
      if(_hit_vector) {
	for(int ihit=0; ihit<_hit_vector->size(); ++ihit) {
	  SQHit *hit = _hit_vector->at(ihit);                    


	  //st. 1 for now D0X
	  if(hit->get_track_id() == trk_id && hit->get_detector_id() ==3 ) {
	    sqhit_st1 = hit;
            		
	    sq_px_st1[n_tracks] = hit->get_truth_px();
	    sq_py_st1[n_tracks] = hit->get_truth_py();
	    sq_pz_st1[n_tracks] = hit->get_truth_pz();
		
	    sq_x_st1[n_tracks] = hit->get_truth_x();
	    sq_y_st1[n_tracks] = hit->get_truth_y();
	    sq_z_st1[n_tracks] = hit->get_truth_z();
	    sq_pos_st1[n_tracks]=hit->get_pos();
	    sq_drift_st1[n_tracks] = hit->get_drift_distance();
 
	    //iF the best reco track available
	    if (Best_recTrack){						
	      double sq_z = hit->get_truth_z();
	      int rec_index = Best_recTrack->getNearestNode(sq_z);
			
	      double rec_z = Best_recTrack->getZ(rec_index);
						
	      if(fabs(sq_z- rec_z>1.)) continue;//to avid mismatch of node
			
	      double p_rec = fabs(1./Best_recTrack->getStateVector(rec_index)[0][0]);
	      double tx_rec = Best_recTrack->getStateVector(rec_index)[1][0];
	      double ty_rec = Best_recTrack->getStateVector(rec_index)[2][0];
	      double x_rec = Best_recTrack->getStateVector(rec_index)[3][0];
	      double y_rec = Best_recTrack->getStateVector(rec_index)[4][0];

	      double x0 = x_rec - tx_rec *rec_z;
	      double y0 = y_rec - ty_rec *rec_z;

	      rec_p_st1[n_tracks] =  p_rec;
	      rec_pz_st1[n_tracks] = p_rec/sqrt(1.+tx_rec*tx_rec+ty_rec*ty_rec);
	      rec_px_st1[n_tracks] = rec_pz_st1[n_tracks]* tx_rec;
	      rec_py_st1[n_tracks] = rec_pz_st1[n_tracks]* ty_rec;

	      rec_x_st1[n_tracks] = x_rec;
	      rec_y_st1[n_tracks] = y_rec;
							
	      rec_drift_st1[n_tracks] = p_geomSvc->getDCA(hit->get_detector_id(), hit->get_element_id(),tx_rec, ty_rec, x0,y0);			
	      // Pull distribution work 
	      double cov00_st1 = Best_recTrack->getCovariance(rec_index)[0][0];
	      double sq_mom_st1 = sqrt(sq_px_st1[n_tracks]*sq_px_st1[n_tracks]+sq_py_st1[n_tracks]*sq_py_st1[n_tracks]+sq_pz_st1[n_tracks]*sq_pz_st1[n_tracks]);
	      pull_q2p_st1[n_tracks] = (fabs(Best_recTrack->getStateVector(rec_index)[0][0]) - 1./sq_mom_st1)/sqrt(cov00_st1);		        
		
	  
	    }//best reco condition
	  }//st.1 work ends

	  //==========

	  //st. 2 for now D2Xp
	  if(hit->get_track_id() == trk_id && hit->get_detector_id() ==15 ) {
			
            		
	    sq_px_st2[n_tracks] = hit->get_truth_px();
	    sq_py_st2[n_tracks] = hit->get_truth_py();
	    sq_pz_st2[n_tracks] = hit->get_truth_pz();
		
	    sq_x_st2[n_tracks] = hit->get_truth_x();
	    sq_y_st2[n_tracks] = hit->get_truth_y();
	    sq_z_st2[n_tracks] = hit->get_truth_z();
	    sq_pos_st2[n_tracks]=hit->get_pos();
	    sq_drift_st2[n_tracks] = hit->get_drift_distance();
 
	    //iF the best reco track available
	    if (Best_recTrack){						
	      double sq_z = hit->get_truth_z();
	      int rec_index = Best_recTrack->getNearestNode(sq_z);
			
	      double rec_z = Best_recTrack->getZ(rec_index);
						
	      if(fabs(sq_z- rec_z>1.)) continue;//to avid mismatch of node
			
	      double p_rec = fabs(1./Best_recTrack->getStateVector(rec_index)[0][0]);
	      double tx_rec = Best_recTrack->getStateVector(rec_index)[1][0];
	      double ty_rec = Best_recTrack->getStateVector(rec_index)[2][0];
	      double x_rec = Best_recTrack->getStateVector(rec_index)[3][0];
	      double y_rec = Best_recTrack->getStateVector(rec_index)[4][0];

	      double x0 = x_rec - tx_rec *rec_z;
	      double y0 = y_rec - ty_rec *rec_z;

	      rec_p_st2[n_tracks] =  p_rec;
	      rec_pz_st2[n_tracks] = p_rec/sqrt(1.+tx_rec*tx_rec+ty_rec*ty_rec);
	      rec_px_st2[n_tracks] = rec_pz_st2[n_tracks]* tx_rec;
	      rec_py_st2[n_tracks] = rec_pz_st2[n_tracks]* ty_rec;


	      rec_x_st2[n_tracks] = x_rec;
	      rec_y_st2[n_tracks] = y_rec;
							
	      rec_drift_st2[n_tracks] = p_geomSvc->getDCA(hit->get_detector_id(), hit->get_element_id(),tx_rec, ty_rec, x0,y0);			

	      // Pull distribution work 
	      double cov00_st2 = Best_recTrack->getCovariance(rec_index)[0][0];
	      double sq_mom_st2 = sqrt(sq_px_st2[n_tracks]*sq_px_st2[n_tracks]+sq_py_st2[n_tracks]*sq_py_st2[n_tracks]+sq_pz_st2[n_tracks]*sq_pz_st2[n_tracks]);
	      pull_q2p_st2[n_tracks] = (fabs(Best_recTrack->getStateVector(rec_index)[0][0]) - 1./sq_mom_st2)/sqrt(cov00_st2);		        		
	  
	    }//if best reco track	
	
	  }//st2. work done

	  //=======================
	  			 	

	  //st. 3 for now D3mXp(id = 27) or D3pXp (id=21)
	  if(hit->get_track_id() == trk_id && (hit->get_detector_id() == 27 ||hit->get_detector_id() == 27)) {          		
	  sq_px_st3[n_tracks] = hit->get_truth_px();
	  sq_py_st3[n_tracks] = hit->get_truth_py();
	  sq_pz_st3[n_tracks] = hit->get_truth_pz();
		
	  sq_x_st3[n_tracks] = hit->get_truth_x();
	  sq_y_st3[n_tracks] = hit->get_truth_y();
	  sq_z_st3[n_tracks] = hit->get_truth_z();
	  sq_pos_st3[n_tracks]=hit->get_pos();
	  sq_drift_st3[n_tracks] = hit->get_drift_distance();
 
	  //iF the best reco track available
	  if (Best_recTrack){						
	  double sq_z = hit->get_truth_z();
	  int rec_index = Best_recTrack->getNearestNode(sq_z);
			
	  double rec_z = Best_recTrack->getZ(rec_index);
						
	  if(fabs(sq_z- rec_z>1.)) continue;//to avid mismatch of node
			
	  double p_rec = fabs(1./Best_recTrack->getStateVector(rec_index)[0][0]);
	  double tx_rec = Best_recTrack->getStateVector(rec_index)[1][0];
	  double ty_rec = Best_recTrack->getStateVector(rec_index)[2][0];
	  double x_rec = Best_recTrack->getStateVector(rec_index)[3][0];
	  double y_rec = Best_recTrack->getStateVector(rec_index)[4][0];

	  double x0 = x_rec - tx_rec *rec_z;
	  double y0 = y_rec - ty_rec *rec_z;

	  rec_p_st3[n_tracks] =  p_rec;
	  rec_pz_st3[n_tracks] = p_rec/sqrt(1.+tx_rec*tx_rec+ty_rec*ty_rec);
	  rec_px_st3[n_tracks] = rec_pz_st3[n_tracks]* tx_rec;
	  rec_py_st3[n_tracks] = rec_pz_st3[n_tracks]* ty_rec;

	  rec_x_st3[n_tracks] = x_rec;
	  rec_y_st3[n_tracks] = y_rec;
							
	  rec_drift_st3[n_tracks] = p_geomSvc->getDCA(hit->get_detector_id(), hit->get_element_id(),tx_rec, ty_rec, x0,y0);			

	  // Pull distribution work 
	  double cov00_st3 = Best_recTrack->getCovariance(rec_index)[0][0];
	  double sq_mom_st3 = sqrt(sq_px_st3[n_tracks]*sq_px_st3[n_tracks]+sq_py_st3[n_tracks]*sq_py_st3[n_tracks]+sq_pz_st3[n_tracks]*sq_pz_st3[n_tracks]);
	  pull_q2p_st3[n_tracks] = (fabs(Best_recTrack->getStateVector(rec_index)[0][0]) - 1./sq_mom_st3)/sqrt(cov00_st3);		        		
	  
	  }//if best reco track	
	
	  }//st3. work done
	  
	  //=======================		
			
	}//sqhit vector loop
      }//if hit vector

	 
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
  _qa_tree->Write();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTrkQA::InitEvalTree() {
 
  PHTFileServer::get().open(_out_name.c_str(), "RECREATE");

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
  _qa_tree->Branch("sq_x_st1",        sq_x_st1,              "sq_x_st1[n_tracks]/F");
  _qa_tree->Branch("sq_y_st1",        sq_y_st1,              "sq_y_st1[n_tracks]/F");
  _qa_tree->Branch("sq_z_st1",        sq_z_st1,              "sq_z_st1[n_tracks]/F");
  _qa_tree->Branch("sq_px_st1",        sq_px_st1,              "sq_px_st1[n_tracks]/F");
  _qa_tree->Branch("sq_py_st1",        sq_py_st1,              "sq_py_st1[n_tracks]/F");
  _qa_tree->Branch("sq_pz_st1",        sq_pz_st1,              "sq_pz_st1[n_tracks]/F");


  _qa_tree->Branch("rec_drift_st1",   rec_drift_st1,      "rec_drift_st1[n_tracks]/F");
  _qa_tree->Branch("rec_px_st1",      rec_px_st1,      "rec_px_st1[n_tracks]/F");
  _qa_tree->Branch("rec_py_st1",      rec_py_st1,      "rec_py_st1[n_tracks]/F");
  _qa_tree->Branch("rec_pz_st1",      rec_pz_st1,      "rec_pz_st1[n_tracks]/F");
  _qa_tree->Branch("rec_x_st1",      rec_x_st1,      "rec_x_st1[n_tracks]/F");
  _qa_tree->Branch("rec_y_st1",      rec_y_st1,      "rec_y_st1[n_tracks]/F");


  _qa_tree->Branch("sq_pos_st2",      sq_pos_st2,         "sq_pos_st2[n_tracks]/F");
  _qa_tree->Branch("sq_drift_st2",    sq_drift_st2,       "sq_drift_st2[n_tracks]/F");
  _qa_tree->Branch("sq_x_st2",        sq_x_st2,              "sq_x_st2[n_tracks]/F");
  _qa_tree->Branch("sq_y_st2",        sq_y_st2,              "sq_y_st2[n_tracks]/F");
  _qa_tree->Branch("sq_z_st2",        sq_z_st2,              "sq_z_st2[n_tracks]/F");
  _qa_tree->Branch("sq_px_st2",        sq_px_st2,              "sq_px_st2[n_tracks]/F");
  _qa_tree->Branch("sq_py_st2",        sq_py_st2,              "sq_py_st2[n_tracks]/F");
  _qa_tree->Branch("sq_pz_st2",        sq_pz_st2,              "sq_pz_st2[n_tracks]/F");


  _qa_tree->Branch("rec_drift_st2",   rec_drift_st2,      "rec_drift_st2[n_tracks]/F");
  _qa_tree->Branch("rec_px_st2",      rec_px_st2,      "rec_px_st2[n_tracks]/F");
  _qa_tree->Branch("rec_py_st2",      rec_py_st2,      "rec_py_st2[n_tracks]/F");
  _qa_tree->Branch("rec_pz_st2",      rec_pz_st2,      "rec_pz_st2[n_tracks]/F");
  _qa_tree->Branch("rec_x_st2",      rec_x_st2,      "rec_x_st2[n_tracks]/F");
  _qa_tree->Branch("rec_y_st2",      rec_y_st2,      "rec_y_st2[n_tracks]/F");

  _qa_tree->Branch("sq_pos_st3",      sq_pos_st3,         "sq_pos_st3[n_tracks]/F");
  _qa_tree->Branch("sq_drift_st3",    sq_drift_st3,       "sq_drift_st3[n_tracks]/F");
  _qa_tree->Branch("sq_x_st3",        sq_x_st3,              "sq_x_st3[n_tracks]/F");
  _qa_tree->Branch("sq_y_st3",        sq_y_st3,              "sq_y_st3[n_tracks]/F");
  _qa_tree->Branch("sq_z_st3",        sq_z_st3,              "sq_z_st3[n_tracks]/F");
  _qa_tree->Branch("sq_px_st3",        sq_px_st3,              "sq_px_st3[n_tracks]/F");
  _qa_tree->Branch("sq_py_st3",        sq_py_st3,              "sq_py_st3[n_tracks]/F");
  _qa_tree->Branch("sq_pz_st3",        sq_pz_st3,              "sq_pz_st3[n_tracks]/F");


  _qa_tree->Branch("rec_drift_st3",   rec_drift_st3,      "rec_drift_st3[n_tracks]/F");
  _qa_tree->Branch("rec_px_st3",      rec_px_st3,      "rec_px_st3[n_tracks]/F");
  _qa_tree->Branch("rec_py_st3",      rec_py_st3,      "rec_py_st3[n_tracks]/F");
  _qa_tree->Branch("rec_pz_st3",      rec_pz_st3,      "rec_pz_st3[n_tracks]/F");
  _qa_tree->Branch("rec_x_st3",      rec_x_st3,      "rec_x_st3[n_tracks]/F");
  _qa_tree->Branch("rec_y_st3",      rec_y_st3,      "rec_y_st3[n_tracks]/F");


  _qa_tree->Branch("pull_q2p_st1",      pull_q2p_st1,      "pull_q2p_st1[n_tracks]/F");
  _qa_tree->Branch("pull_q2p_st2",      pull_q2p_st2,      "pull_q2p_st2[n_tracks]/F");
  _qa_tree->Branch("pull_q2p_st3",      pull_q2p_st3,      "pull_q2p_st3[n_tracks]/F");


  _qa_tree->Branch("chisq_st1",      chisq_st1,            "chisq_st1[n_tracks]/F");
  _qa_tree->Branch("prob_st1",       prob_st1,             "prob_st1[n_tracks]/F");
  _qa_tree->Branch("quality",        quality,             "quality[n_tracks]/F");

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

    sq_x_st1[i]     = std::numeric_limits<float>::max();
    sq_y_st1[i]     = std::numeric_limits<float>::max();
    sq_z_st1[i]     = std::numeric_limits<float>::max();
    sq_px_st1[i]     = std::numeric_limits<float>::max();
    sq_py_st1[i]     = std::numeric_limits<float>::max();
    sq_pz_st1[i]     = std::numeric_limits<float>::max();
    sq_pos_st1[i]     = std::numeric_limits<float>::max();
    sq_drift_st1[i]     = std::numeric_limits<float>::max(); 
 
    rec_x_st1[i]     = std::numeric_limits<float>::max();
    rec_y_st1[i]     = std::numeric_limits<float>::max();
    rec_px_st1[i]     = std::numeric_limits<float>::max();
    rec_py_st1[i]     = std::numeric_limits<float>::max();
    rec_pz_st1[i]     = std::numeric_limits<float>::max();
    rec_drift_st1[i]     = std::numeric_limits<float>::max();
 

    sq_x_st2[i]     = std::numeric_limits<float>::max();
    sq_y_st2[i]     = std::numeric_limits<float>::max();
    sq_z_st2[i]     = std::numeric_limits<float>::max();
    sq_px_st2[i]     = std::numeric_limits<float>::max();
    sq_py_st2[i]     = std::numeric_limits<float>::max();
    sq_pz_st2[i]     = std::numeric_limits<float>::max();
    sq_pos_st2[i]     = std::numeric_limits<float>::max();
    sq_drift_st2[i]     = std::numeric_limits<float>::max();

    rec_x_st2[i]     = std::numeric_limits<float>::max();
    rec_y_st2[i]     = std::numeric_limits<float>::max();
    rec_px_st2[i]     = std::numeric_limits<float>::max();
    rec_py_st2[i]     = std::numeric_limits<float>::max();
    rec_pz_st2[i]     = std::numeric_limits<float>::max();
    rec_drift_st2[i]     = std::numeric_limits<float>::max();
 

    sq_x_st3[i]     = std::numeric_limits<float>::max();
    sq_y_st3[i]     = std::numeric_limits<float>::max();
    sq_z_st3[i]     = std::numeric_limits<float>::max();
    sq_px_st3[i]     = std::numeric_limits<float>::max();
    sq_py_st3[i]     = std::numeric_limits<float>::max();
    sq_pz_st3[i]     = std::numeric_limits<float>::max();
    sq_pos_st3[i]     = std::numeric_limits<float>::max();
    sq_drift_st3[i]     = std::numeric_limits<float>::max();

    rec_x_st3[i]     = std::numeric_limits<float>::max();
    rec_y_st3[i]     = std::numeric_limits<float>::max();
    rec_px_st3[i]     = std::numeric_limits<float>::max();
    rec_py_st3[i]     = std::numeric_limits<float>::max();
    rec_pz_st3[i]     = std::numeric_limits<float>::max();
    rec_drift_st3[i]     = std::numeric_limits<float>::max();

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


//For finding g4hit information in stations (following Kenichi's truth node maker)
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


//Function for finding best reco track
SRecTrack* AnaTrkQA::FindBestRecTrack(SRecEvent *recEvent,  const float true_TargetP)
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
  cout<<" best track inside function: "<<Best_recTrack->getTargetMom().Mag()<<"True P"<<true_TargetP<<endl;
  return Best_recTrack;
  
}

