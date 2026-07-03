#ifndef _EXTRACT_TDC_DIST__H_
#define _EXTRACT_TDC_DIST__H_
#include <fstream>
#include <fun4all/SubsysReco.h>
#include <interface_main/SQEvent.h>
class TFile;
class TH1;

class ExtractTdcDist: public SubsysReco {
 public:
  static const int N_PL = 6;

 private:
  typedef enum { UNDEF, D0, D1, D2, D3p, D3m } ChamType_t;

  std::string m_name;
  ChamType_t m_type;
  int m_pl0;
  std::string m_dir_out;

  SQEvent::TriggerMask m_trig_mask;

  std::vector<int>         m_list_pl_id;
  std::vector<std::string> m_list_pl_name;
  std::vector<int>         m_list_pl_n_ele;

  unsigned int m_n_evt_all;
  unsigned int m_n_evt_trig;
  unsigned int m_n_evt_ana;

  std::ofstream m_ofs;
  TFile* m_file;
  TH1* m_h1_nhit[N_PL];
  TH1* m_h1_ele [N_PL];
  TH1* m_h1_time[N_PL];
  TH1* m_h1_time_wide[N_PL];

 public:
  ExtractTdcDist(const std::string& name="d0");
  virtual ~ExtractTdcDist() {;}

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
};

#endif // _EXTRACT_TDC_DIST__H_
