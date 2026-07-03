R__LOAD_LIBRARY(libGenRoadset)
#include "Common.h"
using namespace std;

void DrawRS()
{
  string rs_id;
  double mass_lo;
  double mass_hi;
  int    inte_cut;
  double frac_cut;
  string list_sig;
  string list_bg;
  GetParams(rs_id, mass_lo, mass_hi, inte_cut, frac_cut, list_sig, list_bg);
  
  DrawRoadset* draw = new DrawRoadset(rs_id);
  //draw->SetMaxCount(20); // default = 20
  UtilRoad::ReadFromFile(
    draw->GetRoadMapPosTop(), 
    draw->GetRoadMapPosBot(), 
    draw->GetRoadMapNegTop(), 
    draw->GetRoadMapNegBot(), 
    rs_id
    );
  draw->Draw();

  delete draw;
  cout << "E N D" << endl;
  exit(0);
}
