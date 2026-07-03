#ifndef _ROAD_LIST_BASE__H_
#define _ROAD_LIST_BASE__H_
#include <string>

/// Base class to hold a set (list or map) of roads.
class RoadListBase {
 protected:
  int m_pn; //< +1 = pos, -1 = neg, 0 = unset
  int m_tb; //< +1 = top, -1 = bot, 0 = unset
  bool m_frozen;

  RoadListBase();

 public:
  virtual ~RoadListBase();

  void SetPNTB  (const int pn, const int tb) { m_pn = pn; m_tb = tb; }
  void SetPosNeg(const int pn) { m_pn = pn; }
  void SetTopBot(const int tb) { m_tb = tb; }
  int  GetPosNeg() const { return m_pn; }
  int  GetTopBot() const { return m_tb; }
  std::string GetPosNegStr() const { return m_pn>0 ? "Pos" : (m_pn<0 ? "Neg" : "!PN"); }
  std::string GetTopBotStr() const { return m_tb>0 ? "Top" : (m_tb<0 ? "Bot" : "!TB"); }

  void Frozen() { m_frozen = true; }
  bool IsFrozen() { return m_frozen; }
};

#endif // _ROAD_LIST_BASE__H_
