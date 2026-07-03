#ifndef _ANA_TRACK__H_
#define _ANA_TRACK__H_
#include <fun4all/SubsysReco.h>
#include "TreeData.h"
class TFile;
class TH2;
class TTree;
class TChain;
class SQEvent;
class SQHitVector;
class SQTrackVector;

/// An example class to analyze hodoscope hits in E1039 DST file.
class AnaTrack: public SubsysReco {
  SQEvent*     m_sq_evt;
  SQHitVector* m_sq_hit_vec;
  SQTrackVector*  m_sq_trk_vec;

  std::string m_file_name;
  TFile*      m_file;
  TTree*      m_tree;
  EventData   m_evt;
  TrackList   m_trk_list;
  
 public:
  AnaTrack(const std::string& name="AnaTrack");
  virtual ~AnaTrack();
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  void SetOutputFileName(const std::string name) { m_file_name = name; }

  static void AnalyzeTree(TChain* tree);
  static void DrawHistIn1D(TH2* h2, const std::string label, const std::string dir_out);
};

#endif // _ANA_TRACK__H_
