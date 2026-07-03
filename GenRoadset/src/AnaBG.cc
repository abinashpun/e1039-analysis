#include <iostream>
#include <fstream>
#include <iomanip>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TLegend.h>
#include <TH1D.h>
#include <UtilAna/UtilTrigger.h>
#include <UtilAna/UtilBeam.h>
#include "RoadInfo.h"
#include "RoadMap.h"
#include "RoadList.h"
#include "UtilRoad.h"
#include "AnaSignal.h"
#include "AnaBG.h"
using namespace std;

AnaBG::AnaBG(const std::string label)
  : AnaBase(label)
  , m_n_evt_used (0)
  , m_n_evt_tb   (0)
  , m_n_evt_fired(0)
  , m_road_map_pos_top(0)
  , m_road_map_pos_bot(0)
  , m_road_map_neg_top(0)
  , m_road_map_neg_bot(0)
  , m_bg_data(0)
  , m_inte_cut(0)
{
  ;
}

AnaBG::~AnaBG() 
{
}

void AnaBG::Init()
{
  if (!m_road_map_pos_top) {
    cout << "AnaBG::Init():  ERROR  No road given.  Abort." << endl;
    exit(1);
  }

  int n_val;
  double* list_val;
  UtilBeam::ListOfRfValues(n_val, list_val);

  ostringstream oss;
  oss << m_dir_out << "/output.root";
  m_file_out = new TFile(oss.str().c_str(), "RECREATE");
 // m_h1_inte_max = new TH1D("h1_inte_max"  , "", 200, 0, 20e3);
  m_h1_inte_max = new TH1D("h1_inte_max"  , "", 1200, 0, 130e3);
 // m_h1_inte = new TH1D("h1_inte"  , "Run 06 NIM3 ", n_val-1, list_val);
  m_h1_inte = new TH1D("h1_inte"  , "", 1200, 0, 1300e3);

  m_ofs << "Intensity cut = " << m_inte_cut << endl;
}

void AnaBG::End()
{
  m_file_out->cd();
  m_file_out->Write();
  m_file_out->Close();
}

void AnaBG::ReadEvents(const char* fname)
{
  gROOT->cd();
  TFile* file = new TFile(fname);
  if (! file->IsOpen()) {
    cout << "Cannot open the file, '" << fname << "'.  Abort." << endl;
    exit(1);
  }
  TTree* tree = (TTree*)file->Get(m_tree_name.c_str());
  if (! tree) {
    cout << "Cannot get the tree, '" << m_tree_name.c_str() << "'.  Abort." << endl;
    exit(1);
  }


  tree->SetBranchAddress(m_branch_name.c_str(), &m_bg_data);

  for (int i_ent = 0; i_ent < tree->GetEntries(); i_ent++) {
    tree->GetEntry(i_ent);
    ProcessOneEvent();
  }
  file->Close();
  delete file;
}

void AnaBG::ProcessOneEvent()
{
  //int  run        =  m_bg_data->run;
  //int  evt        =  m_bg_data->evt;
  //bool fpga1      =  m_bg_data->fpga1;
  int  inte_rfp00 =  m_bg_data->inte_rfp00;
  int  inte_max   =  m_bg_data->inte_max;
  EleList* h1t    = &m_bg_data->h1t;
  EleList* h2t    = &m_bg_data->h2t;
  EleList* h3t    = &m_bg_data->h3t;
  EleList* h4t    = &m_bg_data->h4t;
  EleList* h1b    = &m_bg_data->h1b;
  EleList* h2b    = &m_bg_data->h2b;
  EleList* h3b    = &m_bg_data->h3b;
  EleList* h4b    = &m_bg_data->h4b;

  m_h1_inte->Fill(inte_rfp00);

  if (inte_max == 0) return; // Some readout error
  m_n_evt_used++;

  m_h1_inte_max->Fill(inte_max);

//  if (m_inte_cut > 0 && inte_max >= m_inte_cut) return;

  if(m_inte_cut > 0 && inte_rfp00 >= m_inte_cut) return; //modification by Abi

  if (h1t->size() == 0 || h1b->size() == 0 || 
      h2t->size() == 0 || h2b->size() == 0 || 
      h3t->size() == 0 || h3b->size() == 0 || 
      h4t->size() == 0 || h4b->size() == 0   ) return;
  m_n_evt_tb++;
  
  RoadMap map_top;
  RoadMap map_bot;
  FindAllRoads(h1t, h2t, h3t, h4t, +1, &map_top);
  FindAllRoads(h1b, h2b, h3b, h4b, -1, &map_bot);
  
  bool fired = false;
  if (map_top.Contain(m_road_map_pos_top) &&
      map_bot.Contain(m_road_map_neg_bot)   ) { // T+B
    fired = true;
    for (RoadMap::Iter it = map_top.Begin(); it != map_top.End(); it++) {
      if (m_road_map_pos_top->Find(it->first)) m_road_map_pos_top->AddBG(it->second);
    }
    for (RoadMap::Iter it = map_bot.Begin(); it != map_bot.End(); it++) {
      if (m_road_map_neg_bot->Find(it->first)) m_road_map_neg_bot->AddBG(it->second);
    }
  }
  if (map_top.Contain(m_road_map_neg_top) &&
      map_bot.Contain(m_road_map_pos_bot)   ) { // B+T
    fired = true;
    for (RoadMap::Iter it = map_top.Begin(); it != map_top.End(); it++) {
      if (m_road_map_neg_top->Find(it->first)) m_road_map_neg_top->AddBG(it->second);
    }
    for (RoadMap::Iter it = map_bot.Begin(); it != map_bot.End(); it++) {
      if (m_road_map_pos_bot->Find(it->first)) m_road_map_pos_bot->AddBG(it->second);
    }
  }

  if (fired) m_n_evt_fired++;
}

void AnaBG::Analyze()
{
  cout << "AnaBG::Analyze():" << endl;

  m_ofs << "Event Counts:\n"
        << "N of analyzed events:\n"
        << "  Used  " << m_n_evt_used << "\n"
        << "  T+B   " << m_n_evt_tb << "\n"
        << "  Fired " << m_n_evt_fired << "\n \n"
        << "Expected counts per spill:\n"
        << "  T+B   " << N_RF_PER_SPILL * m_n_evt_tb    / m_n_evt_used << "\n"
        << "  Fired " << N_RF_PER_SPILL * m_n_evt_fired / m_n_evt_used << "\n"
        << endl;

  DrawInteMax();

  m_road_map_pos_top->ScaleBG(N_RF_PER_SPILL / m_n_evt_used);
  m_road_map_pos_bot->ScaleBG(N_RF_PER_SPILL / m_n_evt_used);
  m_road_map_neg_top->ScaleBG(N_RF_PER_SPILL / m_n_evt_used);
  m_road_map_neg_bot->ScaleBG(N_RF_PER_SPILL / m_n_evt_used);

  ///
  /// Regularize BG, namely set its weight to 1 (per spill), so that FoM is always finite.
  ///
  for (RoadMap::Iter it = m_road_map_pos_top->Begin(); it != m_road_map_pos_top->End(); it++) {
    if (it->second->GetWeightBG() == 0) it->second->AddBG(1);
  }
  for (RoadMap::Iter it = m_road_map_pos_bot->Begin(); it != m_road_map_pos_bot->End(); it++) {
    if (it->second->GetWeightBG() == 0) it->second->AddBG(1);
  }
  for (RoadMap::Iter it = m_road_map_neg_top->Begin(); it != m_road_map_neg_top->End(); it++) {
    if (it->second->GetWeightBG() == 0) it->second->AddBG(1);
  }
  for (RoadMap::Iter it = m_road_map_neg_bot->Begin(); it != m_road_map_neg_bot->End(); it++) {
    if (it->second->GetWeightBG() == 0) it->second->AddBG(1);
  }
}

void AnaBG::SetRoads(AnaSignal* ana_signal)
{
  m_road_map_pos_top = ana_signal->GetRoadMapPosTop();
  m_road_map_pos_bot = ana_signal->GetRoadMapPosBot();
  m_road_map_neg_top = ana_signal->GetRoadMapNegTop();
  m_road_map_neg_bot = ana_signal->GetRoadMapNegBot();
  m_road_map_pos_top->Frozen();
  m_road_map_pos_bot->Frozen();
  m_road_map_neg_top->Frozen();
  m_road_map_neg_bot->Frozen();
}

void AnaBG::SetRoads(RoadMap* pos_top, RoadMap* pos_bot, RoadMap* neg_top, RoadMap* neg_bot)
{
  m_road_map_pos_top = pos_top;
  m_road_map_pos_bot = pos_bot;
  m_road_map_neg_top = neg_top;
  m_road_map_neg_bot = neg_bot;
}

void AnaBG::FindAllRoads(const EleList* h1, const EleList* h2, const EleList* h3, const EleList* h4, const int tb, RoadMap* road_map)
{
  if (abs(tb) != 1) return;
  for (EleList::const_iterator it1 = h1->begin(); it1 != h1->end(); it1++) {
  for (EleList::const_iterator it2 = h2->begin(); it2 != h2->end(); it2++) {
  for (EleList::const_iterator it3 = h3->begin(); it3 != h3->end(); it3++) {
  for (EleList::const_iterator it4 = h4->begin(); it4 != h4->end(); it4++) {
    int road = UtilTrigger::Hodo2Road(*it1, *it2, *it3, *it4, tb);
    road_map->AddBG(road, 1.0);
  }
  }
  }
  }
}


void AnaBG::DrawInteMax()
{
  //int bin_cut = m_h1_inte_max->FindBin(m_inte_cut) - 1;
  //double frac_acc = m_h1_inte_max->Integral(0, bin_cut)*1.0 / m_h1_inte_max->Integral()*1.;
  int bin_cut = m_h1_inte->FindBin(m_inte_cut) - 1;
  double frac_acc = m_h1_inte->Integral(0, bin_cut)*1.0 / m_h1_inte->Integral()*1.;

  m_ofs << "NIM3 event fraction = " << frac_acc << " @ inte_rfp00 < " << m_inte_cut<<" m_h1_inte_max->Integral(0, bin_cut): "<< m_h1_inte_max->Integral(0, bin_cut)<<" m_h1_inte_max->Integral(): "<<m_h1_inte_max->Integral() << endl;

  TCanvas* c1 = new TCanvas("c1", "");
  c1->SetGrid();
  c1->SetLogy(true);
  gStyle->SetOptStat(0000);

  m_h1_inte_max->SetTitle("NIM3;Max intensity;Yield");

  m_h1_inte_max->GetXaxis()->SetRangeUser(0, 5000);
  m_h1_inte_max->Draw();
  ostringstream oss;
  oss << m_dir_out << "/h1_inte_max.png";
  c1->SaveAs(oss.str().c_str());

  m_h1_inte_max->GetXaxis()->SetRange();
  m_h1_inte_max->Rebin(10);
  m_h1_inte_max->Draw();
  oss.str("");
  oss << m_dir_out << "/h1_inte_max_all.png";
  c1->SaveAs(oss.str().c_str());

  delete c1;
}
