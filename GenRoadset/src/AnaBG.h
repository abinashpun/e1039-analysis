#ifndef _ANA_BG__H_
#define _ANA_BG__H_
#include <map>
#include "RoadMap.h"
#include "RoadList.h"
#include "TreeData.h"
#include "AnaBase.h"
class TFile;
class TH1;
class AnaSignal;

class AnaBG : public AnaBase {
 protected:
  /// N of RFs per spill
  /** Expected counts per spill = [N of "fired" events] / [N of all events] * [N of filled RFs]
   *  N of filled RFs = 186e6 = 588 * 369000 * 6 / 7
   */
  static constexpr double N_RF_PER_SPILL = 186e6; 
  typedef std::vector<int> EleList;

  int  m_n_evt_used;
  int  m_n_evt_tb;
  int  m_n_evt_fired;

  RoadMap* m_road_map_pos_top;
  RoadMap* m_road_map_pos_bot;
  RoadMap* m_road_map_neg_top;
  RoadMap* m_road_map_neg_bot;

  ///
  /// Inputs
  ///
  BgData* m_bg_data;

  ///
  /// Outputs
  ///
  TFile* m_file_out;
  TH1* m_h1_inte_max;
  TH1* m_h1_inte;

  ///
  /// Parameters
  ///
  int m_inte_cut;

 public:
  AnaBG(const std::string label="ana_bg");
  virtual ~AnaBG();
  virtual void Init();
  virtual void End();
  virtual void ReadEvents(const char* fname);
  virtual void Analyze();

  void SetInteCut(const int inte_cut) { m_inte_cut = inte_cut; }
  int  GetInteCut()            { return m_inte_cut; }

  void SetRoads(AnaSignal* ana_signal);
  void SetRoads(RoadMap* pos_top, RoadMap* pos_bot, RoadMap* neg_top, RoadMap* neg_bot);

 protected:
  virtual void ProcessOneEvent();
  void FindAllRoads(const EleList* h1, const EleList* h2, const EleList* h3, const EleList* h4, const int tb, RoadMap* road_map);
  void DrawInteMax();
};

#endif // _ANA_BG__H_
