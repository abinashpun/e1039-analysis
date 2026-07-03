#ifndef _ANA_SIM_DIMUON__H_
#define _ANA_SIM_DIMUON__H_
#include <map>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <fun4all/SubsysReco.h>
#include "TreeData.h"
class TFile;
class TTree;
class SQEvent;
class SRecEvent;
class SQMCEvent;
class SQTrackVector;
class SQDimuonVector;

/// An example class to analyze the simulated uDST file.
class AnaSimDimuon: public SubsysReco {
  /// Input
  SQEvent       * mi_evt;
  SRecEvent     * mi_srec;
  SQMCEvent     * mi_evt_true;
  SQTrackVector * mi_vec_trk;
  SQDimuonVector* mi_vec_dim;

  /// Output
  TFile* mo_file;
  TTree* mo_tree;
  EventData  mo_evt;
  TrackList  mo_trk_true;
  TrackList  mo_trk_reco;
  DimuonList mo_dim_true;
  DimuonList mo_dim_reco;

  std::map<int, int> m_proc_id_cnt;
  std::map<int, int> m_part_id_1_cnt;
  std::map<int, int> m_part_id_2_cnt;
  std::map<int, int> m_part_id_3_cnt;
  std::map<int, int> m_part_id_4_cnt;
  std::map<std::string, int> m_part_id_12_cnt;

 public:
  AnaSimDimuon(const std::string& name="AnaSimDimuon");
  virtual ~AnaSimDimuon() {;}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 private:
  typedef std::map<int, int> IdMap_t; // For now the key is not ID but index.
  void FindDimuonRelation(IdMap_t& id_map);
};

#endif // _ANA_SIM_DIMUON__H_
