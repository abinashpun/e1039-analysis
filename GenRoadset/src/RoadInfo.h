#ifndef _ROAD_INFO__H_
#define _ROAD_INFO__H_
#include <vector>
#include <cmath>
//#include <map>

/// Class to hold one road.
class RoadInfo {
  int    m_road_id;
  int    m_count;
  double m_weight; //< signal
  int    m_count_bg;
  double m_weight_bg;
  bool   m_on_off;

 public:
  RoadInfo(const int road_id=0);
  RoadInfo(const RoadInfo& obj);
  //RoadInfo(const int road_id, const int count=0, const double weight=0, const bool on_off=true);
  virtual ~RoadInfo();

  void SetRoadID(const int road_id) { m_road_id = road_id; }
  int  GetRoadID() const { return m_road_id; }

  void   Add      (const double weight, const int count=1);
  void   SetCount (const int count    ) { m_count   = count  ; }
  void   SetWeight(const double weight) { m_weight  = weight ; }
  int    GetCount () const { return m_count  ; }
  double GetWeight() const { return m_weight ; }
  void   Scale(const double val);

  void   AddBG      (const double weight, const int count=1);
  void   SetCountBG (const int count    ) { m_count_bg  = count  ; }
  void   SetWeightBG(const double weight) { m_weight_bg = weight ; }
  int    GetCountBG () const { return m_count_bg  ; }
  double GetWeightBG() const { return m_weight_bg ; }
  void   ScaleBG(const double val);

  void   SetOnOff (const bool on_off) { m_on_off  = on_off ; }
  double GetOnOff () const { return m_on_off ; }

  double GetFoM() const { return m_weight / sqrt(m_weight_bg); }
  
  //bool operator<(const RoadInfo& obj) const;
};

#endif // _ROAD_INFO__H_
