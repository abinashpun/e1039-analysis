#include <iostream>
#include <iomanip>
#include <algorithm>
#include "RoadInfo.h"
using namespace std;

RoadInfo::RoadInfo(const int road_id)
  : m_road_id  (road_id)
  , m_count    (0)
  , m_weight   (0)
  , m_count_bg (0)
  , m_weight_bg(0)
  , m_on_off   (true)
{
  ;
}

RoadInfo::RoadInfo(const RoadInfo& obj)
  : m_road_id  (obj.GetRoadID  ())
  , m_count    (obj.GetCount   ())
  , m_weight   (obj.GetWeight  ())
  , m_count_bg (obj.GetCountBG ())
  , m_weight_bg(obj.GetWeightBG())
  , m_on_off   (obj.GetOnOff   ())
{
  ;
}

//RoadInfo::RoadInfo(const int road_id, const int count, const double weight, const bool on_off)
//  : m_road_id(road_id)
//  , m_count  (count  )
//  , m_weight (weight )
//  , m_on_off (on_off )
//{
//  ;
//}

RoadInfo::~RoadInfo()
{
  ;
}

void RoadInfo::Add(const double weight, const int count)
{
  m_count  += count ;
  m_weight += weight;
}

void RoadInfo::Scale(const double val)
{
  m_weight *= val;
}

void RoadInfo::AddBG(const double weight, const int count)
{
  m_count_bg  += count ;
  m_weight_bg += weight;
}

void RoadInfo::ScaleBG(const double val)
{
  m_weight_bg *= val;
}

//bool RoadInfo::operator<(const RoadInfo& obj) const
//{
//  return m_weight < obj.GetWeight() || m_road_id < obj.GetRoadID();
//  //return m_count < obj.GetCount() || m_road_id < obj.GetRoadID();
//}
