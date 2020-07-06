/**
 * \class AnaTrkQA
 * \Analysis module for track QA
 * \author Haiwang, Updated by Abinash Pun
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

// Fun4All includes
#include <fun4all/SubsysReco.h>

// STL includes
#include <vector>
#include <string>
#include <iostream>
#include <list>
#include <map>


class SQRun;
class SQSpillMap;

class SQEvent;
class SQHitMap;
class SQHitVector;

class PHG4TruthInfoContainer;

class SRecEvent;

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

  int TruthEval(PHCompositeNode *topNode);
  int RecoEval(PHCompositeNode *topNode);

  std::string _hit_container_type;

  size_t _event;

  SQRun* _run_header;
  SQSpillMap * _spill_map;

  SQEvent * _event_header;
  SQHitMap *_hit_map;
  SQHitVector *_hit_vector;

  PHG4TruthInfoContainer* _truth;

  SRecEvent* _recEvent;

  std::string _out_name;
	
  TTree* _tout_reco;
  TFile *file;

  int run_id;
  int spill_id;
  float target_pos;
  int event_id;
  int krecstat;
  unsigned short emu_trigger;

  int n_hits;
  int hit_id[10000];
  int detector_id[10000];
  int element_id[10000];
  int hodo_mask[10000];
  float drift_distance[10000];
  float pos[10000];
  float detector_z[10000];

  float truth_x[10000];
  float truth_y[10000];
  float truth_z[10000];
  float truth_pos[10000];

  int n_tracks;
  int rec_id[1000];
  int par_id[1000];
  int pid[1000];
  float gvx[1000];
  float gvy[1000];
  float gvz[1000];
  float gpx[1000];
  float gpy[1000];
  float gpz[1000];
  float gx_st1[1000];
  float gy_st1[1000];
  float gz_st1[1000];
        
       
  float sq_pos_st1[1000];
  float sq_drift_st1[1000];
  float sq_decID[1000];
  float rec_drift_st1[100];        
  float sqx_st1[1000];
  float sqy_st1[1000];
  float sqz_st1[1000];
  float chisq_st1[1000];
  float prob_st1[1000];
  float quality[1000];
       
  float gpx_st1[1000];
  float gpy_st1[1000];
  float gpz_st1[1000];
  float gpt[1000];
  float geta[1000];
  float gphi[1000];
  int gnhits[1000];
  int gndc[1000];
  int gnhodo[1000];
  int gnprop[1000];
  int gndp[1000];
  int ntruhits[1000];
  int nhits[1000];
  int charge[1000];
  float vx[1000];
  float vy[1000];
  float vz[1000];
  float px[1000];
  float py[1000];
  float pz[1000];
  float pt[1000];
  float eta[1000];
  float phi[1000];
  float x_st1[1000];
  float y_st1[1000];
  float px_st1[1000];
  float py_st1[1000];
  float pz_st1[1000];

  int gelmid[1000][128];



  GeomSvc *p_geomSvc;
};


#endif /* _H_AnaTrkQA_H_ */
