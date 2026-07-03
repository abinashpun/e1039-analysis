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
#include "RoadMap.h"
#include "RoadList.h"
#include "UtilRoad.h"
#include "ReAnaBG.h"
using namespace std;

ReAnaBG::ReAnaBG(const std::string label)
  : AnaBG(label)
{
  ;
}

ReAnaBG::~ReAnaBG() 
{
  ;
}

void ReAnaBG::Init()
{
  if (!m_road_map_pos_top) {
    cout << "AnaBG::Init():  ERROR  No road given.  Abort." << endl;
    exit(1);
  }
}

void ReAnaBG::End()
{
  ;
}

void ReAnaBG::ProcessOneEvent()
{
  //int  run        =  m_bg_data->run;
  //int  evt        =  m_bg_data->evt;
  //bool fpga1      =  m_bg_data->fpga1;
  //int  inte_rfp00 =  m_bg_data->inte_rfp00;
  int  inte_max   =  m_bg_data->inte_max;
  EleList* h1t    = &m_bg_data->h1t;
  EleList* h2t    = &m_bg_data->h2t;
  EleList* h3t    = &m_bg_data->h3t;
  EleList* h4t    = &m_bg_data->h4t;
  EleList* h1b    = &m_bg_data->h1b;
  EleList* h2b    = &m_bg_data->h2b;
  EleList* h3b    = &m_bg_data->h3b;
  EleList* h4b    = &m_bg_data->h4b;

  if (inte_max == 0) return; // Some readout error
  m_n_evt_used++;

  if (m_inte_cut > 0 && inte_max >= m_inte_cut) return;
  if (h1t->size() == 0 || h1b->size() == 0 || 
      h2t->size() == 0 || h2b->size() == 0 || 
      h3t->size() == 0 || h3b->size() == 0 || 
      h4t->size() == 0 || h4b->size() == 0   ) return;
  m_n_evt_tb++;
  
  RoadMap map_top;
  RoadMap map_bot;
  FindAllRoads(h1t, h2t, h3t, h4t, +1, &map_top);
  FindAllRoads(h1b, h2b, h3b, h4b, -1, &map_bot);
  
  if ((map_top.ContainEnabled(m_road_map_pos_top) && 
       map_bot.ContainEnabled(m_road_map_neg_bot)   ) || // T+B ||
      (map_top.ContainEnabled(m_road_map_neg_top) &&
       map_bot.ContainEnabled(m_road_map_pos_bot)   )) { // B+T
    m_n_evt_fired++;
        
    m_road_map_top.AddBG(&map_top);
    m_road_map_bot.AddBG(&map_bot);
  }
}

void ReAnaBG::Analyze()
{
  cout << "ReAnaBG::Analyze():" << endl;

  m_ofs << "Event Counts:\n"
        << "N of analyzed events:\n"
        << "  Used  " << m_n_evt_used << "\n"
        << "  T+B   " << m_n_evt_tb << "\n"
        << "  Fired " << m_n_evt_fired << "\n \n"
        << "Expected counts per spill:\n"
        << "  T+B   " << N_RF_PER_SPILL * m_n_evt_tb    / m_n_evt_used << "\n"
        << "  Fired " << N_RF_PER_SPILL * m_n_evt_fired / m_n_evt_used << "\n"
        << endl;

  RoadList road_list_top;
  RoadList road_list_bot;
  road_list_top.Add(&m_road_map_top);
  road_list_bot.Add(&m_road_map_bot);
  road_list_top.SortByBG();
  road_list_bot.SortByBG();
  m_ofs << "  N of fired roads @ Top: " << road_list_top.Size() << endl
        << "  N of fired roads @ Bot: " << road_list_bot.Size() << endl
        << endl;

  m_ofs << "Top: First Five\n";
  UtilRoad::PrintList(&road_list_top, +1, +5, m_ofs);
  m_ofs << "      Last Five\n";
  UtilRoad::PrintList(&road_list_top, -5, -1, m_ofs);
  m_ofs << "Bot: First Five\n";
  UtilRoad::PrintList(&road_list_bot, +1, +5, m_ofs);
  m_ofs << "      Last Five\n";
  UtilRoad::PrintList(&road_list_bot, -5, -1, m_ofs);
}
