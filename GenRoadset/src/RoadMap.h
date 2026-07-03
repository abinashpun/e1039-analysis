#ifndef _ROAD_MAP__H_
#define _ROAD_MAP__H_
#include <map>
#include "RoadListBase.h"
class RoadInfo;

/// Class to hold a non-ordered set (i.e. map) of roads.
class RoadMap : public RoadListBase {
  typedef std::map<int, RoadInfo*> InfoMap;
  InfoMap m_map;

 public:
  typedef InfoMap::iterator Iter;
  typedef InfoMap::const_iterator ConstIter;

  RoadMap();
  virtual ~RoadMap();

  int       Size () const { return m_map.size (); }
  ConstIter Begin() const { return m_map.begin(); }
  ConstIter End  () const { return m_map.end  (); }
  Iter      Begin()       { return m_map.begin(); }
  Iter      End  ()       { return m_map.end  (); }

  RoadInfo* Find(const int road) const;
  bool      Contain(const RoadMap* map) const;
  bool      ContainEnabled(const RoadMap* map) const;

  void Add(const int road, const double weight, const int count=1);
  void Add(const RoadInfo* info);
  void Add(const RoadMap* map);
  void Scale(const double val);

  void AddBG(const int road, const double weight, const int count=1);
  void AddBG(const RoadInfo* info);
  void AddBG(const RoadMap* map);
  void ScaleBG(const double val);
};

#endif // _ROAD_MAP__H_
