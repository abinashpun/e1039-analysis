#ifndef _ANA_CHAM_HIT__H_
#define _ANA_CHAM_HIT__H_
#include <fun4all/SubsysReco.h>
#include "TreeData.h"
class TFile;
class TTree;
class TChain;
class SQEvent;
class SQHitVector;

/// An example class to analyze hodoscope hits in E1039 DST file.
class AnaChamHit: public SubsysReco {
  SQEvent* m_sq_evt;
  SQHitVector* m_sq_hit_vec;

  std::string m_file_name;
  TFile*      m_file;
  TTree*      m_tree;
  EventData   m_evt;
  HitList     m_hit_list;

 public:
  AnaChamHit(const std::string& name="AnaChamHit");
  virtual ~AnaChamHit();
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  void SetOutputFileName(const std::string name) { m_file_name = name; }

  static void AnalyzeTree(TChain* tree);
};

#endif // _ANA_CHAM_HIT__H_
