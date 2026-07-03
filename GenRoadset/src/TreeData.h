#ifndef _TREE_DATA__H_
#define _TREE_DATA__H_
#include <vector>
#include <phool/PHObject.h>
#include <TLorentzVector.h>

struct SignalData {
  double weight;
  double mass;
  double pT;
  double xF;
  double x1;
  double x2;
  double dz;
  double mom;
  double phi;
  double theta;
  int    road_pos;
  int    road_neg;

  SignalData();
  virtual ~SignalData() {;}
  ClassDef(SignalData, 1);
};

struct BgData {
  int    run;
  int    evt;
  bool   fpga1;
  double pot_rfp00; ///< In unit of N of protons.
  int    inte_rfp00; ///< In unit of QIE count.
  int    inte_max;
  std::vector<int> h1t;
  std::vector<int> h2t;
  std::vector<int> h3t;
  std::vector<int> h4t;
  std::vector<int> h1b;
  std::vector<int> h2b;
  std::vector<int> h3b;
  std::vector<int> h4b;

  BgData();
  virtual ~BgData() {;}
  ClassDef(BgData, 1);
};

#endif /* _TREE_DATA__H_ */
