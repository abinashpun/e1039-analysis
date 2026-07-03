#include <iostream>
#include <sstream>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TPaletteAxis.h>
#include <TLine.h>
#include <TMarker.h>
#include <TText.h>
#include <UtilAna/UtilTrigger.h>
#include "RoadInfo.h"
#include "RoadMap.h"
#include "UtilRoad.h"
#include "DrawRoadset.h"
using namespace std;

DrawRoadset::DrawRoadset(const std::string rs_id)
  : AnaBase("draw_rs" + rs_id)
  , m_rs_id(rs_id)
  , m_cnt_max(20)
{
  m_road_map_pos_top.SetPNTB(+1, +1);
  m_road_map_pos_bot.SetPNTB(+1, -1);
  m_road_map_neg_top.SetPNTB(-1, +1);
  m_road_map_neg_bot.SetPNTB(-1, -1);
}

DrawRoadset::~DrawRoadset()
{
  ;
}

void DrawRoadset::Draw()
{
  DrawOneSet(&m_road_map_pos_top);
  DrawOneSet(&m_road_map_pos_bot);
  DrawOneSet(&m_road_map_neg_top);
  DrawOneSet(&m_road_map_neg_bot);
}

void DrawRoadset::DrawOneSet(RoadMap* map)
{
  PlaneConnCount_t list_cnt_12;
  PlaneConnCount_t list_cnt_23;
  PlaneConnCount_t list_cnt_34;
  for (RoadMap::ConstIter it = map->Begin(); it != map->End(); it++) {
    int id = it->first;
    int h1, h2, h3, h4, tb;
    UtilTrigger::Road2Hodo(id, h1, h2, h3, h4, tb);
    list_cnt_12[PlaneConn_t(h1, h2)]++;
    list_cnt_23[PlaneConn_t(h2, h3)]++;
    list_cnt_34[PlaneConn_t(h3, h4)]++;
  }

  TCanvas* c1 = new TCanvas("c1", "", 600, 960);
  c1->SetMargin(0.05, 0.15, 0.05, 0.05); // (l, r, b, t)
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);

  TH2* h2_pa = new TH2D("h2_pa", "", 1,0,1, 1,0,1);
  h2_pa->Fill(0.0, 0.0); // "pa" will be null if "h2_pa" is empty!!
  h2_pa->SetMinimum(0);
  h2_pa->SetMaximum(m_cnt_max);
  h2_pa->Draw("colz");
  c1->Update();
  TPaletteAxis* pa = (TPaletteAxis*)h2_pa->GetListOfFunctions()->FindObject("palette");

  const double DY = 16 + 2; // N of H234T paddles + margin

  ostringstream oss;
  oss << "Roadset " << m_rs_id << " : " << map->GetPosNegStr() << " " << map->GetTopBotStr()
      << ";Station;Paddle";
  TH2* h2_fr = new TH2D("h2_fr", oss.str().c_str(), 1, 0.7, 4.3, 1, -DY/2, DY/2);
  h2_fr->Draw("AH");
  pa->Draw();

  int cnt_max = 0;
  for (PlaneConnCount_t::const_iterator it = list_cnt_12.begin(); it != list_cnt_12.end(); it++) {
    if (it->second > cnt_max) cnt_max = it->second;
  }
  for (PlaneConnCount_t::const_iterator it = list_cnt_23.begin(); it != list_cnt_23.end(); it++) {
    if (it->second > cnt_max) cnt_max = it->second;
  }
  for (PlaneConnCount_t::const_iterator it = list_cnt_34.begin(); it != list_cnt_34.end(); it++) {
    if (it->second > cnt_max) cnt_max = it->second;
  }
  cout << "Max count of inter-station connections = " << cnt_max << endl;
  if (cnt_max > m_cnt_max) {
    cout << "  cnt_max > m_cnt_max.  Consider changing m_cnt_max." << endl;
  }

  DrawConnection(&list_cnt_12, 1, m_cnt_max);
  DrawConnection(&list_cnt_23, 2, m_cnt_max);
  DrawConnection(&list_cnt_34, 3, m_cnt_max);

  DrawElement(1, 23, -DY/2-0.5);
  DrawElement(2, 16, -DY/2-0.5);
  DrawElement(3, 16, -DY/2-0.5);
  DrawElement(4, 16, -DY/2-0.5);

  oss.str("");
  oss << m_dir_out << "/road_count_" << map->GetPosNegStr() << "_" << map->GetTopBotStr() << ".png";
  c1->SaveAs(oss.str().c_str());
  delete h2_fr;
  delete h2_pa;
  delete c1;
}

void DrawRoadset::DrawConnection(const PlaneConnCount_t* list_cnt, const int plane, const int count_max)
{
  double cent1 = (plane==1 ? (23+1)/2.0 : (16+1)/2.0);
  double cent2 = (16+1)/2.0; // "plane+1" never be "1".
  double step1 = (plane==1 ? 0.7 : 1.0);
  double step2 = 1.0;
  
  for (PlaneConnCount_t::const_iterator it = list_cnt->begin(); it != list_cnt->end(); it++) {
    int h1  = it->first.first;
    int h2  = it->first.second;
    int cnt = it->second;
    int color = TColor::GetPalette()[254 * cnt / count_max];
    TLine line;
    line.SetLineColor(color);
    line.DrawLine(plane, step1*(h1-cent1), plane+1, step2*(h2-cent2));
  }
}

void DrawRoadset::DrawElement(const int st, const int n_ele, const double y_label)
{
  double cent = (n_ele + 1) / 2.0;
  double step = st==1 ? 0.7 : 1.0;
  TMarker marker;
  marker.SetMarkerStyle(21);
  for (int ele = 1; ele <= n_ele; ele++) {
    marker.DrawMarker(st, step * (ele - cent));
  }
  TText text;
  text.SetTextAlign(22); // center center
  text.DrawText(st, step * (1     - cent) - 0.5, "1");

  ostringstream oss;
  oss << n_ele;
  text.DrawText(st, step * (n_ele - cent) + 0.5, oss.str().c_str());

  oss.str("");
  oss << "St." << st;
  text.DrawText(st, y_label, oss.str().c_str());
}

void DrawRoadset::DrawElementV1(const int st, const int ele1, const int ele2)
{
  double ele_cnt = (ele1 + ele2)/2.0;
  TMarker marker;
  marker.SetMarkerStyle(21);
  for (int ele = ele1; ele <= ele2; ele++) {
    marker.DrawMarker(st, ele - ele_cnt);
  }
  ostringstream oss;
  oss << ele1;
  TText text;
  text.SetTextAlign(22); // center center
  text.DrawText(st, ele1 - ele_cnt - 0.5, oss.str().c_str());
  oss.str("");
  oss << ele2;
  text.DrawText(st, ele2 - ele_cnt + 0.5, oss.str().c_str());

  oss.str("");
  oss << "St." << st;
  text.DrawText(st, -10.3, oss.str().c_str());
}
