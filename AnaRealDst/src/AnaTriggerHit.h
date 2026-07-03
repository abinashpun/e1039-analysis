#ifndef _ANA_TRIGGER_HIT__H_
#define _ANA_TRIGGER_HIT__H_
#include <fun4all/SubsysReco.h>
class TFile;
class TTree;
class TH1;
class TH2;
class SQEvent;
class SQHitVector;

/// A simple SubsysReco module to analyze the trigger (i.e. V1495 TDC) hit.
/**
 * It outputs small numbers of trees and histograms.
 * You can modify either/both of them as you need.
 */
class AnaTriggerHit: public SubsysReco {
  static const std::vector<std::string> m_list_det_name;
  std::vector<int> m_list_det_id;

  bool m_intime; ///< Require in-time or not.
  int m_level; ///< Level 0 (A), 1 (B) or 2 (C)
  SQEvent* m_evt;
  SQHitVector* m_hit_vec;

  TFile* m_file;
  TTree* m_tree;
  char   b_det_name[16];
  int    b_det_id;
  int    b_ele_id;
  double b_time;
  TH1*   m_h1_ele[99];
  TH1*   m_h1_nhit[99];
  TH1*   m_h1_time[99];

 public:
  AnaTriggerHit(const std::string& name="AnaTriggerHit");
  virtual ~AnaTriggerHit() {;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void UseInTime(const bool val) { m_intime = val; }
  int  UseInTime() const  { return m_intime; }
  void SetLevel(const int level) { m_level = level; }
  int  GetLevel() const   { return m_level; }
};

#endif // _ANA_TRIGGER_HIT__H_
