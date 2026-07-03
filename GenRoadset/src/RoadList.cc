#include <iostream>
#include <iomanip>
#include <algorithm>
#include "RoadInfo.h"
#include "RoadMap.h"
#include "RoadList.h"
using namespace std;

RoadList::RoadList()
{
  ;
}

RoadList::~RoadList()
{
  m_vec.clear();
}

void RoadList::Add(RoadInfo* info)
{
  m_vec.push_back(info);
}

void RoadList::Add(const RoadMap* road_map)
{
  for (RoadMap::ConstIter it = road_map->Begin(); it != road_map->End(); it++) Add(it->second);
}

void RoadList::SortBySignal()
{
  sort(m_vec.begin(), m_vec.end(),
       [](RoadInfo* a, RoadInfo* b) { return a->GetWeight() > b->GetWeight(); }
    );
}

void RoadList::SortByBG()
{
  sort(m_vec.begin(), m_vec.end(),
       [](RoadInfo* a, RoadInfo* b) { return a->GetWeightBG() > b->GetWeightBG(); }
    );
}

/// Sort the list by Figure of Merit, i.e. "signal / sqrt(bg)".
void RoadList::SortByFoM()
{
  sort(m_vec.begin(), m_vec.end(),
       [](RoadInfo* a, RoadInfo* b) { return a->GetFoM() > b->GetFoM(); }
    );
}
