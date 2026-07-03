#ifndef _UTIL_ROAD__H_
#define _UTIL_ROAD__H_
#include <iostream>
class RoadMap;
class RoadList;

namespace UtilRoad {
  //void MakeOrderedList(const RoadMap* road_map, RoadList* road_list);
  void PrintList(RoadList* road_list, int i_begin=0, int i_end=0, std::ostream& os=std::cout);

  void WriteToFile(const RoadMap* pos_top, const RoadMap* pos_bot, const RoadMap* neg_top, const RoadMap* neg_bot, const std::string id);
  void ReadFromFile(RoadMap* pos_top, RoadMap* pos_bot, RoadMap* neg_top, RoadMap* neg_bot, const std::string id);
}

#endif // _UTIL_ROAD__H_
