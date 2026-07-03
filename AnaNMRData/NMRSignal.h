#ifndef _NMR_SIGNAL__H_
#define _NMR_SIGNAL__H_

class NMRSignal {
  int m_evt_num;
  std::vector<double> m_list_pt;
public:
  NMRSignal() : m_evt_num(0) {;}
  virtual ~NMRSignal() {;}

  void SetEventNum(const int evt_num) { m_evt_num = evt_num; }
  int  GetEventNum() const { return m_evt_num; }

  int GetNumPoint() { return m_list_pt.size(); }
  double GetPoint(const int idx) { return m_list_pt[idx]; }

  void AddPoint(const double val) { m_list_pt.push_back(val); }
  void ClearPoints() { m_list_pt.clear(); }
};

#endif // _NMR_SIGNAL__H_
