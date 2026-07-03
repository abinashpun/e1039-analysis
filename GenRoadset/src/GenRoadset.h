#ifndef _GEN_ROADSET__H_
#define _GEN_ROADSET__H_
#include "RoadMap.h"
#include "RoadList.h"
#include "AnaBase.h"
class AnaSignal;

class GenRoadset : public AnaBase {
 protected:
  RoadMap* m_map_pos_top;
  RoadMap* m_map_pos_bot;
  RoadMap* m_map_neg_top;
  RoadMap* m_map_neg_bot;

  RoadList m_list_pos_top;
  RoadList m_list_pos_bot;
  RoadList m_list_neg_top;
  RoadList m_list_neg_bot;

 public:
  GenRoadset(const std::string label="gen_rs");
  virtual ~GenRoadset();

  void SetMap(AnaSignal* ana_signal);

  void SymmetrizeRoads();

  void GenerateRoadSetBySignal(const double frac_cut, const std::string rs_id);
  void GenerateRoadSetByFoM(const double frac_cut, const std::string rs_id);
  void GenerateSymmetricRoadSetByFoM(const double frac_cut, const std::string rs_id);

  //RoadMap * GetRoadMapPosTop () { return &m_map_pos_top ; }
  //RoadMap * GetRoadMapPosBot () { return &m_map_pos_bot ; }
  //RoadMap * GetRoadMapNegTop () { return &m_map_neg_top ; }
  //RoadMap * GetRoadMapNegBot () { return &m_map_neg_bot ; }

  //RoadList* GetRoadListPosTop() { return &m_list_pos_top; }
  //RoadList* GetRoadListPosBot() { return &m_list_pos_bot; }
  //RoadList* GetRoadListNegTop() { return &m_list_neg_top; }
  //RoadList* GetRoadListNegBot() { return &m_list_neg_bot; }

 protected:
  void DisableRoadsByWeightFraction(RoadList* road_list, const double frac);
  void DisableRoadsByFoM(RoadList* road_list, const double frac);

  int FlipRoadLeftRight(const int road);
  int FlipRoadTopBottom(const int road);
};

#endif /* _GEN_ROADSET__H_ */
