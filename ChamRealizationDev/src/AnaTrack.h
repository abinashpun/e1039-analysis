#ifndef _ANA_TRACK_H_
#define _ANA_TRACK_H_
#include <fstream>
#include <map>
#include <fun4all/SubsysReco.h>
class TFile;
class TH1;
class SQEvent;
class SQHitVector;
class SRecEvent;
class SQTrackVector;
class TrackletVector;

class AnaTrack: public SubsysReco {
  static const int N_PL = 30; // N of all chamber planes

  SQEvent*        m_evt;
  SRecEvent*      m_srec;
  SQHitVector*    m_vec_hit;
  SQTrackVector*  m_vec_trk_true;
  TrackletVector* m_vec_trklet;

  TFile* m_file_out;

  TH1* m_h1_ntrk;
  TH1* m_h1_nhit; // N of hits/tracklet
  TH1* m_h1_mom;
  TH1* m_h1_rchi2;
  TH1* m_h1_x0;
  TH1* m_h1_y0;
  TH1* m_h1_tx;
  TH1* m_h1_ty;

  TH1* m_h1_ele [N_PL+1];
  TH1* m_h1_pos [N_PL+1];
  TH1* m_h1_time[N_PL+1]; // TDC (not drift) time
  TH1* m_h1_dd  [N_PL+1]; // drift distance
  TH1* m_h1_td  [N_PL+1]; // track distance = track position - wire position

 public:
  AnaTrack(const std::string& name="AnaTrack");
  virtual ~AnaTrack() {;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 private:
  typedef std::map<int, int> IdMap_t; // For now the key is not ID but index.
  void FindTrackRelation (IdMap_t& id_map);

  void DrawTrackPar();
  void DrawOneGroup(const std::string name, TH1* h1[], const int i_pl0, const int i_pl1);
};

#endif // _ANA_TRACK_H_
