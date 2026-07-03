#ifndef _SUBSYS_RECO_BG__H_
#define _SUBSYS_RECO_BG__H_
#include <vector>
#include <fun4all/SubsysReco.h>
#include "TreeData.h"
class TFile;
class TTree;
class SQEvent;
class SQMCEvent;
class SQHitVector;

class SubsysRecoBG: public SubsysReco {
  typedef enum { DEFAULT, FULL_BG } Mode_t;
  Mode_t m_mode;

  /// Input
  SQEvent*     mi_evt;
  SQMCEvent*   mi_mc_evt;
  SQHitVector* mi_vec_hit;

  /// Output
  TFile* mo_file;
  TTree* mo_tree;
  BgData mo_bg;

 public:
  SubsysRecoBG(const std::string &name="SubsysRecoBG");
  virtual ~SubsysRecoBG() {;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void FullBgMode() { m_mode = FULL_BG; }

 protected:
  void ExtractHits(const SQHitVector* hit_vec, const std::string det_name, std::vector<int>& list_ele);
};

#endif // _SUBSYS_RECO_BG__H_
