#ifndef _RE_ANA_FULL_BG__H_
#define _RE_ANA_FULL_BG__H_
#include "AnaBG.h"

class ReAnaBG : public AnaBG {
 protected:
  RoadMap  m_road_map_top;
  RoadMap  m_road_map_bot;

 public:
  ReAnaBG(const std::string label="re_ana_bg");
  virtual ~ReAnaBG();
  virtual void Init();
  virtual void End();
  virtual void Analyze();

  RoadMap * GetRoadMapTop () { return &m_road_map_top ; }
  RoadMap * GetRoadMapBot () { return &m_road_map_bot ; }

 protected:
  virtual void ProcessOneEvent();
};

#endif // _RE_ANA_FULL_BG__H_
