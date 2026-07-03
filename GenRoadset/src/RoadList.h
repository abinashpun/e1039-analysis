#ifndef _ROAD_LIST__H_
#define _ROAD_LIST__H_
#include <vector>
#include "RoadListBase.h"
class RoadInfo;
class RoadMap;

/// Class to hold an ordered set (i.e. vector) of roads.
class RoadList : public RoadListBase {
  typedef std::vector<RoadInfo*> InfoVec;
  InfoVec m_vec;

 public:
  RoadList();
  virtual ~RoadList();

  int  Size() const { return m_vec.size(); }
  RoadInfo*       Get(int i)       { return m_vec.at(i); }
  const RoadInfo* Get(int i) const { return m_vec.at(i); }
  void Add(RoadInfo* info);
  void Add(const RoadMap* road_map);

  void SortBySignal();
  void SortByBG();
  void SortByFoM();
};

#endif // _ROAD_LIST__H_
