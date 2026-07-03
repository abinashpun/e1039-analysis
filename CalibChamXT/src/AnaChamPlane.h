#ifndef _ANA_CHAM_PLANE_H_
#define _ANA_CHAM_PLANE_H_
#include <fstream>
#include <fun4all/SubsysReco.h>
class TFile;
class TH1;
class SQEvent;
class SQHitVector;

class AnaChamPlane: public SubsysReco {
  std::string m_plane_name;
  int m_plane_id;
  int m_n_ele;
  int m_n_evt_all;
  int m_n_hit_all;

  SQEvent*     m_evt;
  SQHitVector* m_hit_vec;

  TFile* m_file_out;
  TH1* m_h1_nhit;
  TH1* m_h1_ele ;
  TH1* m_h1_pos ;
  TH1* m_h1_time;
  TH1* m_h1_dist;

 public:
  AnaChamPlane(const std::string& plane_name);
  virtual ~AnaChamPlane() {;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
};

#endif // _ANA_CHAM_PLANE_H_
