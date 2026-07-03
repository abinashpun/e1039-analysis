#ifndef _ANA_HODO_HIT__H_
#define _ANA_HODO_HIT__H_
#include <fun4all/SubsysReco.h>
class TFile;
class TTree;
class TH1;
class TH2;
class SQEvent;
class SQHitVector;

/// An example class to analyze hodoscope hits in E1039 DST file.
class AnaHodoHit: public SubsysReco {
  static const std::vector<std::string> m_list_det_name;
  std::vector<int> m_list_det_id;

  SQEvent* m_evt;
  SQHitVector* m_hit_vec;

  TFile* m_file;
  TTree* m_tree;
  char   b_det_name[16];
  int    b_det;
  int    b_ele;
  double b_time;
  TH1*   m_h1_ele[99];
  TH1*   m_h1_nhit[99];

 public:
  AnaHodoHit(const std::string& name="AnaHodoHit");
  virtual ~AnaHodoHit() {;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
};

#endif // _ANA_HODO_HIT__H_
