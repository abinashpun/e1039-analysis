R__LOAD_LIBRARY(libGenRoadset)
#include "Common.h"
using namespace std;

void ApplyRoadset2Signal()
{
  string rs_id;
  double mass_lo;
  double mass_hi;
  int    inte_cut;
  double frac_cut;
  string list_sig;
  string list_bg;
  GetParams(rs_id, mass_lo, mass_hi, inte_cut, frac_cut, list_sig, list_bg);

  ReAnaSignal* ana_sig = new ReAnaSignal("app_sig_rs" + rs_id);
  ana_sig->SetInputTreeName  ("signal_tree");
  ana_sig->SetInputBranchName("signal_data");
  ana_sig->SetMassRange(mass_lo, mass_hi);
  ana_sig->Init();
  UtilRoad::ReadFromFile(
    ana_sig->GetRoadMapPosTop(), 
    ana_sig->GetRoadMapPosBot(), 
    ana_sig->GetRoadMapNegTop(), 
    ana_sig->GetRoadMapNegBot(), 
    rs_id
    );
  ana_sig->ReadEventsFromFileList(list_sig.c_str()); 
  ana_sig->Analyze();
  ana_sig->End();

  delete ana_sig;
  cout << "E N D" << endl;
  exit(0);
}
