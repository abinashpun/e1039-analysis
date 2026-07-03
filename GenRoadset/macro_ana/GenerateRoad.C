R__LOAD_LIBRARY(libGenRoadset)
#include "Common.h"
using namespace std;

/// Macro to generate a roadset based on the S/N FoM.
void GenerateRoad()
{
  const bool make_symm = true;
  string rs_id;
  double mass_lo;
  double mass_hi;
  int    inte_cut;
  double frac_cut;
  string list_sig;
  string list_bg;
  GetParams(rs_id, mass_lo, mass_hi, inte_cut, frac_cut, list_sig, list_bg);

  cout << "\nRead and analyze the signal files." << endl;
  AnaSignal* ana_sig = new AnaSignal();
  ana_sig->SetInputTreeName  ("signal_tree");
  ana_sig->SetInputBranchName("signal_data");
  ana_sig->SetMassRange(mass_lo, mass_hi);
  ana_sig->Init();
  ana_sig->ReadEventsFromFileList(list_sig.c_str());
  ana_sig->Analyze();
  ana_sig->End();

  cout << "\nRead and analyze the BG files." << endl;
  AnaBG* ana_bg = new AnaBG();
  ana_bg->SetInputTreeName  ("bg_tree");
  ana_bg->SetInputBranchName("bg_data");
  ana_bg->SetInteCut(inte_cut);
  ana_bg->SetRoads(ana_sig);
  ana_bg->Init();
  ana_bg->ReadEventsFromFileList(list_bg.c_str());
  ana_bg->Analyze();
  ana_bg->End();

  cout << "\nGenerate a roadset." << endl;
  GenRoadset* gen_rs = new GenRoadset("gen_rs" + rs_id);
  gen_rs->SetMap(ana_sig);
  if (make_symm) {
    gen_rs->SymmetrizeRoads();
    gen_rs->GenerateSymmetricRoadSetByFoM(frac_cut, rs_id);
  } else {
    gen_rs->GenerateRoadSetByFoM(frac_cut, rs_id);
  }

  delete gen_rs;
  delete ana_bg;
  delete ana_sig;
  cout << "E N D" << endl;
  exit(0);
}

/// Macro to generate a roadset based on the GMC signal yield.
void GenerateRoadBySignal(const string rs_id="901")
{
  const double  mass_lo = 4.0;
  const double  mass_hi = 5.0;
  const double frac_cut = 0.05;

  AnaSignal* ana_sig = new AnaSignal();
  ana_sig->SetInputTreeName  ("signal_tree");
  ana_sig->SetInputBranchName("signal_data");
  ana_sig->SetMassRange(mass_lo, mass_hi);
  ana_sig->Init();
  ana_sig->ReadEventsFromFileList("list_signal.txt");
  ana_sig->Analyze();

  GenRoadset* gen_rs = new GenRoadset("gen_rs" + rs_id);
  gen_rs->SetMap(ana_sig);
  gen_rs->GenerateRoadSetBySignal(frac_cut, rs_id);

  ana_sig->End();

  delete gen_rs;
  delete ana_sig;
  cout << "E N D" << endl;
  exit(0);
}
