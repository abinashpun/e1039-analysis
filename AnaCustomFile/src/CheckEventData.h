#ifndef _CHECK_EVENT_DATA__H_
#define _CHECK_EVENT_DATA__H_
#include <iostream>
#include <map>
#include <fun4all/SubsysReco.h>
class SQEvent;
class SQHitVector;
class TFile;
class TTree;

/// A SubsysReco module to look into simulated hits.
class CheckEventData: public SubsysReco {
  /// Input
  SQEvent* m_evt;
  SQHitVector* m_hit_vec;

  /// Output
  TFile* m_file;
  TTree* m_tree;

 public:
  CheckEventData();
  virtual ~CheckEventData() {;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
};

#endif // _CHECK_EVENT_DATA__H_
