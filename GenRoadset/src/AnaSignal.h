#ifndef _ANA_SIGNAL__H_
#define _ANA_SIGNAL__H_
#include "RoadMap.h"
#include "RoadList.h"
#include "TreeData.h"
#include "AnaBase.h"
class TTree;

class AnaSignal : public AnaBase {
 protected:
  double  m_mass_lo;
  double  m_mass_hi;
  RoadMap m_road_map_pos_top;
  RoadMap m_road_map_pos_bot;
  RoadMap m_road_map_neg_top;
  RoadMap m_road_map_neg_bot;

  SignalData* m_sig_data;

 public:
  AnaSignal(const std::string label="ana_signal");
  virtual ~AnaSignal();

  void SetMassRange(const double mass_lo, const double mass_hi);

  virtual void Init();
  virtual void End();
  virtual void ReadEvents(const char* fname="signal_tree.root");
  virtual void Analyze();

  RoadMap * GetRoadMapPosTop () { return &m_road_map_pos_top ; }
  RoadMap * GetRoadMapPosBot () { return &m_road_map_pos_bot ; }
  RoadMap * GetRoadMapNegTop () { return &m_road_map_neg_top ; }
  RoadMap * GetRoadMapNegBot () { return &m_road_map_neg_bot ; }

 protected:
  virtual void ProcessOneEvent();
};

#endif // _ANA_SIGNAL__H_
