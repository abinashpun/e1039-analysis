#ifndef _ANA_TRIG_SIGNAL__H_
#define _ANA_TRIG_SIGNAL__H_
#include <fun4all/SubsysReco.h>
class TFile;
class TTree;
class TH1;
class TH2;
class SQEvent;
class SQHitVector;

/// An example class to analyze hodoscope hits in E1039 DST file.
class AnaTrigSignal: public SubsysReco {
  static const int N_ELE = 5;
  static const std::vector<std::string> m_list_det_name;
  std::vector<int> m_list_det_id;

  SQEvent* m_evt;
  SQHitVector* m_hit_vec;

  TFile* m_file;
  TTree* m_tree;
  char   b_det_name[32];
  int    b_det;
  int    b_ele;
  double b_time;
  TH2*   m_h2_time_ele[99];
  TH1*   m_h1_nhit[99];

 public:
  AnaTrigSignal(const std::string& name="AnaTrigSignal");
  virtual ~AnaTrigSignal() {;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
};

#endif // _ANA_TRIG_SIGNAL__H_
