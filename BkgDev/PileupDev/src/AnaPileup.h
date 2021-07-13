/**
 * \class AnaPileup
 * \Analysis module for bkg pile up study
 * \author  Abinash Pun
 *
 * Created: 07-05-2020
 */

#ifndef _H_AnaPileup_H_
#define _H_AnaPileup_H_

// ROOT
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TMatrixD.h>
#include <TVector3.h>

// Fun4All includes
#include <fun4all/SubsysReco.h>

// STL includes
#include <vector>
#include <string>
#include <iostream>
#include <list>
#include <map>

class TVector3;
class TLorentzVector;
class TClonesArray;

class SQRun;
class SQSpillMap;
class SQEvent;
class SQHitMap;
class SQHitVector;
class SQHit;

class PHG4TruthInfoContainer;
class PHG4HitContainer;
class SRecEvent;
class SRecTrack;
class GeomSvc;

class TFile;
class TTree;

class AnaPileup: public SubsysReco {

 public:

  AnaPileup(const std::string &name = "AnaPileup.root");
  virtual ~AnaPileup() {
  }

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  int InitEvalTree();
  int ResetEvalVars();

  const std::string& get_hit_container_choice() const {
    return _hit_container_type;
  }

  void set_hit_container_choice(const std::string& hitContainerChoice) {
    _hit_container_type = hitContainerChoice;
  }

  const std::string& get_out_name() const {
    return _out_name;
  }

  void set_out_name(const std::string& outName) {
    _out_name = outName;
  }

 private:

  int GetNodes(PHCompositeNode *topNode);

  int Eval(PHCompositeNode *topNode); 
 
  std::string _hit_container_type;

  size_t _event;
  SQRun* _run_header;
  SQSpillMap * _spill_map;

  SQEvent * _event_header;
  SQHitMap *_hit_map;
  SQHitVector *_hit_vector;

  PHG4TruthInfoContainer* _truth;
  SRecEvent* _recEvent;

  PHG4HitContainer *g4hc_d1x;
  PHG4HitContainer *g4hc_d2xp;
  PHG4HitContainer *g4hc_d3px;
  PHG4HitContainer *g4hc_d3mx;

  PHG4HitContainer *g4hc_h1t;
  PHG4HitContainer *g4hc_h1b;
  PHG4HitContainer *g4hc_h2t;
  PHG4HitContainer *g4hc_h2b;
  PHG4HitContainer *g4hc_h3t;
  PHG4HitContainer *g4hc_h3b;
  PHG4HitContainer *g4hc_h4t;
  PHG4HitContainer *g4hc_h4b;

  PHG4HitContainer *g4hc_p1y1;
  PHG4HitContainer *g4hc_p1y2;
  PHG4HitContainer *g4hc_p1x1;
  PHG4HitContainer *g4hc_p1x2;
  PHG4HitContainer *g4hc_p2x1;
  PHG4HitContainer *g4hc_p2x2;
  PHG4HitContainer *g4hc_p2y1;
  PHG4HitContainer *g4hc_p2y2;

  std::string _out_name;
	
  TTree* _tout_reco;
  TTree* _qa_tree;
  TFile *file;

  int run_id;
  int spill_id;
  int event_id;
  unsigned short emu_trigger;

/// truth info
  int n_tracks;
  int pid[9999];

//generated info at vtx
  float gvx[9999];
  float gvy[9999];
  float gvz[9999];
  float gpx[9999];
  float gpy[9999];
  float gpz[9999];

  float gpt[9999];
  float geta[9999];
  float gphi[9999];

//hodoscope sqhit info
  TClonesArray* pos_H1T;
  TClonesArray* pos_H1B;
  TClonesArray* mom_H1T;
  TClonesArray* mom_H1B;

  TClonesArray* pos_H2T;
  TClonesArray* pos_H2B;
  TClonesArray* mom_H2T;
  TClonesArray* mom_H2B;

  TClonesArray* pos_H3T;
  TClonesArray* pos_H3B;
  TClonesArray* mom_H3T;
  TClonesArray* mom_H3B;

  TClonesArray* pos_H4T;
  TClonesArray* pos_H4B;
  TClonesArray* mom_H4T;
  TClonesArray* mom_H4B;

  TClonesArray* pos_H1X;
  TClonesArray* mom_H1X;

  TClonesArray* pos_H2X;
  TClonesArray* mom_H2X;

  TClonesArray* pos_H3X;
  TClonesArray* mom_H3X;
  
  TClonesArray* pos_H4X;
  TClonesArray* mom_H4X;

//DCC sqhit info
  TClonesArray* pos_D0U;
  TClonesArray* pos_D0Up;
  TClonesArray* pos_D0X;
  TClonesArray* pos_D0Xp;
  TClonesArray* pos_D0V;
  TClonesArray* pos_D0Vp;

  TClonesArray* pos_D2U;
  TClonesArray* pos_D2Up;
  TClonesArray* pos_D2X;
  TClonesArray* pos_D2Xp;
  TClonesArray* pos_D2V;
  TClonesArray* pos_D2Vp;

  TClonesArray* pos_D3pU;
  TClonesArray* pos_D3pUp;
  TClonesArray* pos_D3pX;
  TClonesArray* pos_D3pXp;
  TClonesArray* pos_D3pV;
  TClonesArray* pos_D3pVp;

  TClonesArray* pos_D3mU;
  TClonesArray* pos_D3mUp;
  TClonesArray* pos_D3mX;
  TClonesArray* pos_D3mXp;
  TClonesArray* pos_D3mV;
  TClonesArray* pos_D3mVp;
  TClonesArray*pos_D3mXp_OR_D3pXp;


  TClonesArray*mom_D0X;
  TClonesArray*mom_D2Xp;
  TClonesArray*mom_D3mXp_OR_D3pXp;

//=========

///reco info

  int n_recTracks;
  int n_hits;

  int nhits[10000];
  int charge[10000];
  int nhits_st1[10000];
  int nhits_st2[10000];
  int nhits_st3[10000];

  float rec_vx[10000];
  float rec_vy[10000];
  float rec_vz[10000];
  float rec_px[10000];
  float rec_py[10000];
  float rec_pz[10000];
  float rec_pt[10000];
  float rec_eta[10000];
  float rec_phi[10000];

  float chisq_st1[10000];
  float prob_st1[10000];
  float quality[10000];

  GeomSvc *p_geomSvc;
};


#endif /* _H_AnaPileup_H_ */
