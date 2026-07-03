#include <iostream>
#include <iomanip>
#include <TFile.h>
#include <TTree.h>
#include "RoadInfo.h"
#include "RoadMap.h"
#include "RoadList.h"
#include "UtilRoad.h"
#include "AnaSignal.h"
using namespace std;

AnaSignal::AnaSignal(const std::string label)
  : AnaBase(label)
  , m_mass_lo(0.0)
  , m_mass_hi(100.0)
  , m_sig_data(0)
{
  m_road_map_pos_top.SetPNTB(+1, +1);
  m_road_map_pos_bot.SetPNTB(+1, -1);
  m_road_map_neg_top.SetPNTB(-1, +1);
  m_road_map_neg_bot.SetPNTB(-1, -1);
}

AnaSignal::~AnaSignal()
{
  ;
}

void AnaSignal::SetMassRange(const double mass_lo, const double mass_hi)
{
  m_mass_lo = mass_lo;
  m_mass_hi = mass_hi;
}

void AnaSignal::Init()
{
  ;
}

void AnaSignal::End()
{
  ;
}

void AnaSignal::ReadEvents(const char* fname)
{
  //gROOT->cd();
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

  tree->SetBranchAddress(m_branch_name.c_str(), &m_sig_data);

  for (int i_ent = 0; i_ent < tree->GetEntries(); i_ent++) {
    tree->GetEntry(i_ent);
    ProcessOneEvent();
  }

  file->Close();
  delete file;
}

void AnaSignal::ProcessOneEvent()
{
  double weight = m_sig_data->weight;
  double mass   = m_sig_data->mass;
  int road_pos  = m_sig_data->road_pos;
  int road_neg  = m_sig_data->road_neg;

  if (mass < m_mass_lo || mass > m_mass_hi) return;
  if (road_pos * road_neg >= 0) return; // Require T+B or B+T

  RoadMap* map = &(road_pos > 0 ? m_road_map_pos_top : m_road_map_pos_bot);
  map->Add(road_pos, weight);
  
  map = &(road_neg > 0 ? m_road_map_neg_top : m_road_map_neg_bot);
  map->Add(road_neg, weight);
}

void AnaSignal::Analyze()
{
  cout << "AnaSignal::Analyze():" << endl;
  cout << "  Do nothing for now.  Probably it will scale the signal weight by the simulated integrated luminosity." << endl;

  //m_road_list_pos_top.reserve(m_road_map_pos_top.size());
  //m_road_list_pos_bot.reserve(m_road_map_pos_bot.size());
  //m_road_list_neg_top.reserve(m_road_map_neg_top.size());
  //m_road_list_neg_bot.reserve(m_road_map_neg_bot.size());

  //UtilRoad::MakeOrderedList(&m_road_map_pos_top, &m_road_list_pos_top);
  //UtilRoad::MakeOrderedList(&m_road_map_pos_bot, &m_road_list_pos_bot);
  //UtilRoad::MakeOrderedList(&m_road_map_neg_top, &m_road_list_neg_top);
  //UtilRoad::MakeOrderedList(&m_road_map_neg_bot, &m_road_list_neg_bot);

//  m_road_list_pos_top.Add(&m_road_map_pos_top);
//  m_road_list_pos_bot.Add(&m_road_map_pos_bot);
//  m_road_list_neg_top.Add(&m_road_map_neg_top);
//  m_road_list_neg_bot.Add(&m_road_map_neg_bot);
  //m_road_list_pos_top.SortBySignal();
  //m_road_list_pos_bot.SortBySignal();
  //m_road_list_neg_top.SortBySignal();
  //m_road_list_neg_bot.SortBySignal();

  //cout << "  pos top: " << m_road_list_pos_top.Size() << endl
  //     << "  pos bot: " << m_road_list_pos_bot.Size() << endl
  //     << "  neg top: " << m_road_list_neg_top.Size() << endl
  //     << "  neg bot: " << m_road_list_neg_bot.Size() << endl
  //     << endl;

  //m_ofs << "Pos Top: First Five\n";
  //UtilRoad::PrintList(&m_road_list_pos_top, +1, +5, m_ofs);
  //m_ofs << "          Last Five\n";
  //UtilRoad::PrintList(&m_road_list_pos_top, -5, -1, m_ofs);
  //m_ofs << "Pos Bot: First Five\n";
  //UtilRoad::PrintList(&m_road_list_pos_bot, +1, +5, m_ofs);
  //m_ofs << "          Last Five\n";
  //UtilRoad::PrintList(&m_road_list_pos_bot, -5, -1, m_ofs);
  //m_ofs << "Neg Top: First Five\n";
  //UtilRoad::PrintList(&m_road_list_neg_top, +1, +5, m_ofs);
  //m_ofs << "          Last Five\n";
  //UtilRoad::PrintList(&m_road_list_neg_top, -5, -1, m_ofs);
  //m_ofs << "Neg Bot: First Five\n";
  //UtilRoad::PrintList(&m_road_list_neg_bot, +1, +5, m_ofs);
  //m_ofs << "          Last Five\n";
  //UtilRoad::PrintList(&m_road_list_neg_bot, -5, -1, m_ofs);
}
