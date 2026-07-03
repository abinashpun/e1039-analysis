//#include <iostream>
//#include <iomanip>
//#include <algorithm>
#include "RoadInfo.h"
#include "RoadMap.h"
using namespace std;

RoadMap::RoadMap()
{
  ;
}

RoadMap::~RoadMap()
{
  for (Iter it = m_map.begin(); it != m_map.end(); it++) {
    if (it->second) delete it->second;
  }
  m_map.clear();
}

RoadInfo* RoadMap::Find(const int road) const
{
  return m_map.find(road) != m_map.end()  ?  m_map.at(road)  :  0;
}

bool RoadMap::Contain(const RoadMap* map) const
{
  for (RoadMap::ConstIter it = map->Begin(); it != map->End(); it++) {
    if (m_map.find(it->first) != m_map.end()) return true;
  }
  return false;
}

bool RoadMap::ContainEnabled(const RoadMap* map) const
{
  for (RoadMap::ConstIter it = map->Begin(); it != map->End(); it++) {
    if (! it->second->GetOnOff()) continue;
    if (m_map.find(it->first) != m_map.end()) return true;
  }
  return false;
}

void RoadMap::Add(const int road, const double weight, const int count)
{
  RoadInfo* info = Find(road);
  if (! info) {
    if (m_frozen) return;
    info = new RoadInfo(road);
    m_map[road] = info;
  }
  info->Add(weight, count);
}

void RoadMap::Add(const RoadInfo* info)
{
  Add(info->GetRoadID(), info->GetWeight(), info->GetCount());
}

void RoadMap::Add(const RoadMap* map)
{
  for (RoadMap::ConstIter it = map->Begin(); it != map->End(); it++) Add(it->second);
}

void RoadMap::Scale(const double val)
{
  for (RoadMap::Iter it = m_map.begin(); it != m_map.end(); it++) it->second->Scale(val);
}

void RoadMap::AddBG(const int road, const double weight, const int count)
{
  RoadInfo* info = Find(road);
  if (! info) {
    if (m_frozen) return;
    info = new RoadInfo(road);
    m_map[road] = info;
  }
  info->AddBG(weight, count);
}

void RoadMap::AddBG(const RoadInfo* info)
{
  AddBG(info->GetRoadID(), info->GetWeightBG(), info->GetCountBG());
}

void RoadMap::AddBG(const RoadMap* map)
{
  for (RoadMap::ConstIter it = map->Begin(); it != map->End(); it++) AddBG(it->second);
}

void RoadMap::ScaleBG(const double val)
{
  for (RoadMap::Iter it = m_map.begin(); it != m_map.end(); it++) it->second->ScaleBG(val);
}

