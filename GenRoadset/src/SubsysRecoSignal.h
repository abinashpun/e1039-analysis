#ifndef _SUBSYS_RECO_SIGNAL__H_
#define _SUBSYS_RECO_SIGNAL__H_
#include <map>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <fun4all/SubsysReco.h>
#include "TreeData.h"
class TFile;
class TTree;
class SQMCEvent;
class SQTrackVector;
class SQDimuonVector;
class SQHitVector;

class SubsysRecoSignal: public SubsysReco {
  /// Input
  SQMCEvent     * mi_evt_true;
  SQTrackVector * mi_vec_trk;
  SQDimuonVector* mi_vec_dim;
  SQHitVector   * mi_vec_hit;

  /// Output
  TFile*     mo_file;
  TTree*     mo_tree;
  SignalData mo_sig;

 public:
  SubsysRecoSignal(const std::string &name="SubsysRecoSignal");
  virtual ~SubsysRecoSignal() {;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
};

#endif // _SUBSYS_RECO_SIGNAL__H_
