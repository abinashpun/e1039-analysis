#ifndef _DRAW_ROADSET__H_
#define _DRAW_ROADSET__H_
#include "RoadMap.h"
#include "AnaBase.h"
class TTree;

class DrawRoadset : public AnaBase {
 protected:
  typedef std::pair<int, int> PlaneConn_t;
  typedef std::map<PlaneConn_t, int> PlaneConnCount_t;

  std::string m_rs_id;
  double  m_cnt_max;
  RoadMap m_road_map_pos_top;
  RoadMap m_road_map_pos_bot;
  RoadMap m_road_map_neg_top;
  RoadMap m_road_map_neg_bot;

 public:
  DrawRoadset(const std::string rs_id);
  virtual ~DrawRoadset();

  void SetMaxCount(const double cnt_max) { m_cnt_max = cnt_max; }

  RoadMap* GetRoadMapPosTop () { return &m_road_map_pos_top ; }
  RoadMap* GetRoadMapPosBot () { return &m_road_map_pos_bot ; }
  RoadMap* GetRoadMapNegTop () { return &m_road_map_neg_top ; }
  RoadMap* GetRoadMapNegBot () { return &m_road_map_neg_bot ; }

  void Draw();

 protected:
  void DrawOneSet(RoadMap* map);
  void DrawConnection(const PlaneConnCount_t* list_cnt, const int plane, const int count_max);
  void DrawElement  (const int st, const int n_ele, const double y_label);
  void DrawElementV1(const int st, const int ele1, const int ele2);
};

#endif // _DRAW_ROADSET__H_
