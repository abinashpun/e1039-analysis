R__LOAD_LIBRARY(libGenRoadset)
#include "Common.h"
using namespace std;

void ApplyRoadset2BG(const int inte_cut=0)
{
  string rs_id;
  double mass_lo;
  double mass_hi;
  int    inte_cut_gen;
  double frac_cut;
  string list_sig;
  string list_bg;
  GetParams(rs_id, mass_lo, mass_hi, inte_cut_gen, frac_cut, list_sig, list_bg);

  RoadMap map_pos_top;
  RoadMap map_pos_bot;
  RoadMap map_neg_top;
  RoadMap map_neg_bot;
  UtilRoad::ReadFromFile(&map_pos_top, &map_pos_bot, &map_neg_top, &map_neg_bot, rs_id);

  ostringstream oss;
  oss << "app_bg_rs" << rs_id << "_inte" << setfill('0') << setw(4) << inte_cut;
  ReAnaBG* ana_bg = new ReAnaBG(oss.str().c_str());
  ana_bg->SetInputTreeName  ("bg_tree");
  ana_bg->SetInputBranchName("bg_data");
  ana_bg->SetInteCut(inte_cut);
  ana_bg->SetRoads(&map_pos_top, &map_pos_bot, &map_neg_top, &map_neg_bot);
  ana_bg->Init();
  ana_bg->ReadEventsFromFileList(list_bg.c_str()); 
  //ana_bg->ReadEventsFromFileList("list_bg.txt");
//  if (rs_id == "103") ana_bg->ReadEventsFromFileList("list_bg_NIM3e906.txt");
  //else if (rs_id == "113") ana_bg->ReadEventsFromFileList("list_bg_fullsimRun06.txt");
  //ana_bg->ReadEventsFromFileList("list_bg_fullsimRun06.txt");
  //ana_bg->ReadEventsFromFileList("list_bg_NIM3e906.txt");  
//  if(rs_id == "123" || rs_id == "129") ana_bg->ReadEventsFromFileList("list_bg_fullsimRun06_reverseKMAG.txt");
//  if(rs_id == "113" || rs_id == "119") ana_bg->ReadEventsFromFileList("list_bg_fullsimRun06.txt"); 

  ana_bg->Analyze();
  ana_bg->End();

  delete ana_bg;
  cout << "E N D" << endl;
  exit(0);
}
