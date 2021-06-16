/**
 * \class AnaTrkQA
 * \Analysis module for track QA
 * \author  Abinash Pun
 *
 * Created: 07-05-2020
 */

#ifndef _H_AnaTrkQA_H_
#define _H_AnaTrkQA_H_

// ROOT
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <ktracker/SRecEvent.h>

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

class SQRun;
class SQSpillMap;
class SQEvent;
class SQHitMap;
class SQHitVector;
class SQHit;
class SQDimuonVector;

class PHG4TruthInfoContainer;
class PHG4HitContainer;
class SRecEvent;
class SRecTrack;
class GeomSvc;

class TFile;
class TTree;

class AnaTrkQA: public SubsysReco {

 public:

  AnaTrkQA(const std::string &name = "AnaTrkQA.root");
  virtual ~AnaTrkQA() {
  }

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void set_legacy_rec_container(bool b) { legacyContainer = b; }

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

  int TruthRecoEval(PHCompositeNode *topNode);
  int DimuonInfo(PHCompositeNode* topNode);

  bool legacyContainer;
  
  bool FindG4HitAtStation(const int trk_id, const PHG4HitContainer* g4hc, TVector3* pos, TLorentzVector* mom);
  int  FindCommonHitIDs(std::vector<int>& hitidvec1, std::vector<int>& hitidvec2);
  SRecTrack* FindBestMomRecTrack(SRecEvent *recEvent, const float true_P); 
  bool FindG4HitAtHodo(const int trk_id, const PHG4HitContainer* g4hc);
  bool FindG4HitAtProp(const int trk_id, const PHG4HitContainer* g4hc);
 
  std::string _hit_container_type;

  size_t _event;
  SQRun* _run_header;
  SQSpillMap * _spill_map;

  SQEvent * _event_header;
  SQHitMap *_hit_map;
  SQHitVector *_hit_vector;

  PHG4TruthInfoContainer* _truth;
  SRecEvent* _recEvent;

  // dimuon vectors
  SQDimuonVector* dimuonVector;
  SQDimuonVector* recDimuonVector;

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
  TTree* tree1;
  TTree* tree2;
  TFile* file;

  int run_id;
  int spill_id;
  float target_pos;
  int event_id;
  int krecstat;
  int kalman_stat;
  unsigned short emu_trigger;

  int n_hits;
  int nhits_st1;
  int nhits_st2;
  int nhits_st3;
  int hit_id;
  int detector_id;
  int element_id;
  int hodo_mask;
  float drift_distance;
  float pos;
  float detector_z;

  float truth_x;
  float truth_y;
  float truth_z;
  float truth_pos;

  int n_tracks;
  int n_recTracks;
  int rec_id;
  int par_id;
  int pid;

  float gx_st1;
  float gy_st1;
  float gz_st1;

  float sq_drift_st1;
  float sq_decID;

  float sq_drift_st2;
  float sq_drift_st3;

  float chisq_st1;
  float prob_st1;
  float quality;
  float fit_time;
  
  float pull_q2p_st1;
  float pull_q2p_st2;
  float pull_q2p_st3;

  float rec_drift_st1;
  float rec_drift_st2;
  float rec_drift_st3;

  TVector3* gpos;
  TVector3* gmom;
  TVector3* rec_mom;
  TVector3* rec_pos;

  TVector3* ac_mom;

  TVector3* sq_pos_st1;
  TVector3* sq_mom_st1;
  TVector3* sq_pos_st2;
  TVector3* sq_mom_st2;
  TVector3* sq_pos_st3;
  TVector3* sq_mom_st3;
  TVector3* rec_mom_st1;
  TVector3* rec_pos_st1;
  TVector3* rec_pos_st2;
  TVector3* rec_mom_st2;
  TVector3* rec_pos_st3;
  TVector3* rec_mom_st3;

// dimuon info

double mass;
TVector3* vtx;
TVector3* pmom;
TVector3* nmom;

double rec_mass;
TVector3* rec_pmom;
TVector3* rec_nmom;
TVector3* rec_ppos;
TVector3* rec_npos;
TVector3* rec_vtx;

  float gpx_st1;
  float gpy_st1;
  float gpz_st1;
  int gnhits;
  int gndc;
  int gnhodo;
  int gnprop;
  int gndp;
  int ntruhits;
  int nhits;
  int charge;
  float pull_state00[100];

  int gelmid[1000][128];



  GeomSvc *p_geomSvc;
};


#endif /* _H_AnaTrkQA_H_ */
