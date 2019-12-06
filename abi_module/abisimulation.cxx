/*============================================================================= 
  Author: Abinash Pun
  Goal: Trying to make the vertex distribution from PYTHIA more realistic
  (may be add some realistic distribution function)
  Strategy: write an class similar to photons (but with cmake instead of autogen.sh, follow the example module from Haiwang)
  ============================================================================*/
#include "abisimulation.h"
#include <TLorentzVector.h>

//#include <geom_svc/GeomSvc.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Particlev2.h>
#include <g4main/PHG4HitDefs.h>
#include <g4main/PHG4VtxPoint.h>
#include <E906LegacyGen/SQDimuonTruthInfoContainer.h>
#include  <E906LegacyGen/SQMCDimuon.h>
#include <TFile.h>
#include <TTree.h>

#include <cstring>
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <limits>
#include <tuple>
#include <cassert>
#include <iostream>


#include <boost/lexical_cast.hpp>
#define LogError(exp)		std::cout<<"ERROR: "  <<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl

using namespace std;


abisimulation::abisimulation(const std::string &name)
  : SubsysReco("ABISIMULATION")
{
  outfilename = name;
  //initialize global variables to -999 so that they have a spot in memory
  //initialize_to_zero(); 
  ResetVars();

  //add other initializers here 
}

int abisimulation::Init(PHCompositeNode *topnode)
{
  file = new TFile(outfilename.c_str(), "RECREATE");
  InitTree(); 
  return Fun4AllReturnCodes::EVENT_OK;

}


int abisimulation::InitRun(PHCompositeNode *topNode) {

 
 

  //p_geomSvc = GeomSvc::instance();

  int ret = GetNodes(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return Fun4AllReturnCodes::EVENT_OK;
}


int abisimulation::process_event(PHCompositeNode *topnode)
{
  
  
   _dimuoninfo = findNode::getClass<SQDimuonTruthInfoContainer>(topnode, "DimuonInfo");
    if (!_dimuoninfo) {
      LogError("! _dimuon");
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    if(_dimuoninfo->get_Dimuon_xs()<=0)  return Fun4AllReturnCodes::ABORTEVENT;
    
    dimuon_xs = _dimuoninfo->get_Dimuon_xs();    
    dimuon_CosThetaCS =_dimuoninfo->get_Dimuon_cosThetaCS() ;
    dimuon_phiCS =_dimuoninfo->get_Dimuon_phiCS() ;	
    dimuon_m =_dimuoninfo->get_Dimuon_m() ;
    dimuon_pt =_dimuoninfo->get_Dimuon_pt();
    dimuon_xf =  _dimuoninfo->get_Dimuon_xF();
    // cout<<"Inside analysis module: dimuon cross section: "<<_dimuoninfo->get_Dimuon_xs()<<endl;
   // }

  int suni = 0;
  if(_truth) {
    for(auto iter=_truth->GetPrimaryParticleRange().first;
    	iter!=_truth->GetPrimaryParticleRange().second;
    	++iter) {

      // for(auto iter=_truth->GetSecondaryParticleRange().first;
      // 	iter!=_truth->GetSecondaryParticleRange().second;
      //  	++iter) {
  
      PHG4Particle * par= iter->second;
      
      if(suni>1) continue;
      int vtx_id =  par->get_vtx_id();
      PHG4VtxPoint* vtx = _truth->GetVtx(vtx_id);
    // const PHG4Particle * last_primary =
    //   _truth->GetMap().rbegin()->second;
    //   const PHG4VtxPoint* vtx = _truth->GetPrimaryVtx(last_primary->get_vtx_id());
      // assert(primary_vtx);
  
    truthpid[suni] = par->get_pid();
    truth_vtxx[suni] = vtx->get_x();
    truth_vtxy[suni] = vtx->get_y();
    truth_vtxz[suni] = vtx->get_z();
    truthpx[suni] = par->get_px();
    truthpy[suni] = par->get_py();
    truthpz[suni] = par->get_pz();
    truthe[suni] = par->get_e();
    
 

    // cout<<"Inside Analysis module; track_id: "<<par->get_track_id()<<" PID: "<<par->get_pid()<<" name: "<<par->get_name()<<" Energy: "<<par->get_e()<<" suni: "<<suni<<endl;
    
    suni++;
   
  
    //Calculation  from scratch
    //@@
    
    if (suni>1){ 
       if(truthpid[0]==truthpid[1]) continue;
      	TLorentzVector mup, mum;
	mup.SetPxPyPzE(truthpx[0],truthpy[0],truthpz[0], truthe[0]);
	mum.SetPxPyPzE(truthpx[1],truthpy[1],truthpz[1], truthe[1]);
	// TLorentzVector dimuon = mup+mum;
	// dimuon_theta = TMath::Cos(dimuon.Theta());


	Double_t mp = 0.938;
	Double_t ebeam = 120.;

	TLorentzVector p_beam(0., 0., sqrt(ebeam*ebeam - mp*mp), ebeam);
	TLorentzVector p_target(0., 0., 0., mp);

	TLorentzVector p_cms = p_beam + p_target;
	TLorentzVector p_sum = mup+mum;


	double fMass = p_sum.M();
	double fpT = p_sum.Perp();

	// fx1 = (p_target*p_sum)/(p_target*p_cms);
	// fx2 = (p_beam*p_sum)/(p_beam*p_cms);

	Double_t s = p_cms.M2();
	Double_t sqrts = p_cms.M();
	TVector3 bv_cms = p_cms.BoostVector();
	p_sum.Boost(-bv_cms);

	double fxF = 2.*p_sum.Pz()/sqrts/(1. - fMass*fMass/s);
	double fCosTh = 2.*(mum.E()*mup.Pz() - mup.E()*mum.Pz())/fMass/TMath::Sqrt(fMass*fMass + fpT*fpT);
	double fPhi = TMath::ATan(2.*TMath::Sqrt(fMass*fMass + fpT*fpT)/fMass*(mum.Px()*mup.Py() - mup.Px()*mum.Py())/(mup.Px()*mup.Px() - mum.Px()*mum.Px() + mup.Py()*mup.Py() - mum.Py()*mum.Py()));


	//dimuon_CosThetaCS = fCosTh;
	//dimuon_phiCS = fPhi;
	
	//dimuon_m = fMass;
	dimuon_pt = fpT;
	dimuon_xf = fxF;
	//	cout<<"Inside analysis module: dimuon, mass, phi, costheta: "<<fMass<<" "<<fPhi<<" "<<fCosTh<<endl;

  }
    
    //@@
    
    }
    
   
	if(truthpid[0]==-13||truthpid[0]==13) {
      truth_tree->Fill();
	}//To filter the empty events not passing the genereateDimuons conditions
 


  }
 


  return 0;
}

int abisimulation::End(PHCompositeNode *topnode)
{
  std::cout << " DONE PROCESSING " << endl;

  file->Write();
  file->Close();
  return 0;
}

void abisimulation::ResetVars()
{
  truth_vtxx[2]= {0};
  truth_vtxy[2]={0};
  truth_vtxz[2]={0};
  truthpid[2] = {0};
  truthpx[2]={0};
  truthpy[2]={0};
  truthpz[2]={0};
  truthe[2]={0};
  dimuon_CosThetaCS = 0.;
  dimuon_phiCS = 0.;
  dimuon_xs = 1.;
  dimuon_m= 0.;
  dimuon_pt = 0.;
  dimuon_xf = 0.;
}

void abisimulation::InitTree()
{
  truth_tree = new TTree("truthtree","a tree with all truth information from generator");
  truth_tree->Branch("truth_vtxx", truth_vtxx, "truth_vtxx[2]/F");
  truth_tree->Branch("truth_vtxy", truth_vtxy, "truth_vtxy[2]/F");
  truth_tree->Branch("truth_vtxz", truth_vtxz, "truth_vtxz[2]/F");
  truth_tree->Branch("truthpx", truthpx, "truthpx[2]/F");
  truth_tree->Branch("truthpy", truthpy, "truthpy[2]/F");
  truth_tree->Branch("truthpz", truthpz, "truthpz[2]/F");
  truth_tree->Branch("truthe", truthe, "truthe[2]/F");
  truth_tree->Branch("truthpid", truthpid, "truthpid[2]/I");
  truth_tree->Branch("dimuon_CosThetaCS", &dimuon_CosThetaCS, "dimuon_CosThetaCS/F");
  truth_tree->Branch("dimuon_phiCS", &dimuon_phiCS, "dimuon_phiCS/F");
  truth_tree->Branch("dimuon_xs", &dimuon_xs, "dimuon_xs/F");
  truth_tree->Branch("dimuon_m", &dimuon_m, "dimuon_m/F");
  truth_tree->Branch("dimuon_xf", &dimuon_xf, "dimuon_xf/F");
  truth_tree->Branch("dimuon_pt", &dimuon_pt, "dimuon_pt/F");

}

 int abisimulation::GetNodes(PHCompositeNode* topNode)
{

  _truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truth) {
    LogError("!_truth");
    return Fun4AllReturnCodes::ABORTEVENT;
  }
 
  return Fun4AllReturnCodes::EVENT_OK;
}
