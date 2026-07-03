#include <iomanip>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <interface_main/SQHitVector.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <geom_svc/GeomSvc.h>
#include <UtilAna/UtilSQHit.h>
#include "ExtractTdcDist.h"
using namespace std;

ExtractTdcDist::ExtractTdcDist(const std::string& name)
  : SubsysReco("ExtractTdcDist"+name)
  , m_name    (name)
  , m_n_evt_all (0)
  , m_n_evt_trig(0)
  , m_n_evt_ana (0)
{
  if      (name == "d0" ) m_type = D0 ;
  else if (name == "d1" ) m_type = D1 ;
  else if (name == "d2" ) m_type = D2 ;
  else if (name == "d3p") m_type = D3p;
  else if (name == "d3m") m_type = D3m;
  else {
    cout << "ExtractTdcDist:  Bad name (" << name << ").  Abort." << endl;
    exit(1);
  }
}

///
/// Event selection
///   NIM1: (H1X||H1Y)&&(H2X||H2Y)&&(H3X||H4Y1)&&(H4X||H4Y2)
///   NIM2: (H1X||H1Y)&&(H2X||H2Y)
///   NIM3: Random
///   NIM4: (H2X||H2Y)&&(H4X||H4Y2)
int ExtractTdcDist::Init(PHCompositeNode* topNode)
{
  double t0_lo;
  double t0_hi;
  switch (m_type) {
  case D0 :  m_pl0 =  1; m_trig_mask = SQEvent::NIM2; t0_lo = 1000; t0_hi = 1400; break;
  case D1 :  m_pl0 =  7; m_trig_mask = SQEvent::NIM2; t0_lo =  850; t0_hi = 1250; break;
  case D2 :  m_pl0 = 13; m_trig_mask = SQEvent::NIM4; t0_lo =  850; t0_hi = 1250; break;
  case D3p:  m_pl0 = 19; m_trig_mask = SQEvent::NIM4; t0_lo =  740; t0_hi = 1140; break;
  case D3m:  m_pl0 = 25; m_trig_mask = SQEvent::NIM4; t0_lo =  850; t0_hi = 1250; break;
  }

  m_dir_out = "tdc_dist_" + m_name;
  gSystem->mkdir(m_dir_out.c_str(), true);

  GeomSvc* geom = GeomSvc::instance();
  m_ofs.open( (m_dir_out+"/output.txt").c_str() );
  m_file = new TFile( (m_dir_out+"/output.root").c_str(), "RECREATE");

  ostringstream oss;
  oss << setfill('0');

  for (int i_pl = 0; i_pl < N_PL; i_pl++) {
    int det_id = m_pl0 + i_pl;
    string name = geom->getDetectorName(det_id);
    int   n_ele = geom->getPlaneNElements(det_id); 
    m_list_pl_id   .push_back(det_id);
    m_list_pl_name .push_back(name);
    m_list_pl_n_ele.push_back(n_ele);
    cout << "  " << setw(6) << name << " " << det_id << " " << n_ele << endl;

    oss.str("");
    oss << "h1_nhit_" << setw(2) << det_id;
    m_h1_nhit[i_pl] = new TH1D(oss.str().c_str(), "", 50, -0.5, 49.5);
    oss.str("");
    oss << name << ";N of hits/plane/event;Hit count";
    m_h1_nhit[i_pl]->SetTitle(oss.str().c_str());

    oss.str("");
    oss << "h1_ele_" << setw(2) << det_id;
    m_h1_ele[i_pl] = new TH1D(oss.str().c_str(), "", n_ele, 0.5, n_ele+0.5);
    oss.str("");
    oss << name << ";Element ID;Hit count";
    m_h1_ele[i_pl]->SetTitle(oss.str().c_str());

    const double DT = 20/9.0; // 4/9 ns per single count of Taiwan TDC
    int i_t_lo = (int)(t0_lo / DT);
    int i_t_hi = (int)(t0_hi / DT) + 1;
    double t_lo = i_t_lo * DT;
    double t_hi = i_t_hi * DT;
    int    t_n  = i_t_hi - i_t_lo;
    oss.str("");
    oss << "h1_time_" << setw(2) << det_id;
    m_h1_time[i_pl] = new TH1D(oss.str().c_str(), "", t_n, t_lo, t_hi);
    oss.str("");
    oss << name << ";tdcTime;Hit count";
    m_h1_time[i_pl]->SetTitle(oss.str().c_str());

    t_lo =    0 * DT;
    t_hi = 1000 * DT;
    t_n  = 1000;
    oss.str("");
    oss << "h1_time_wide_" << setw(2) << det_id;
    m_h1_time_wide[i_pl] = new TH1D(oss.str().c_str(), "", t_n, t_lo, t_hi);
    oss.str("");
    oss << name << ";tdcTime;Hit count";
    m_h1_time_wide[i_pl]->SetTitle(oss.str().c_str());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int ExtractTdcDist::InitRun(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int ExtractTdcDist::process_event(PHCompositeNode* topNode)
{
  SQEvent* event       = findNode::getClass<SQEvent    >(topNode, "SQEvent");
  SQHitVector* hit_vec = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if (!event || !hit_vec) return Fun4AllReturnCodes::ABORTEVENT;
  m_n_evt_all++;

  if (! event->get_trigger(m_trig_mask)) return Fun4AllReturnCodes::EVENT_OK;
  m_n_evt_trig++;

  bool good_event = true;
  int n_hit_tot = 0;
  shared_ptr<SQHitVector> list_hv[N_PL];
  for (int i_pl = 0; i_pl < N_PL; i_pl++) {
    list_hv[i_pl] = shared_ptr<SQHitVector>(UtilSQHit::FindHits(hit_vec, m_list_pl_id[i_pl]));

    int n_hit = list_hv[i_pl]->size();
    m_h1_nhit[i_pl]->Fill(n_hit);
    if (n_hit > 1) {
      good_event = false;
      break;
    }
    n_hit_tot += n_hit;
  }
  if (!good_event || n_hit_tot < 5) return Fun4AllReturnCodes::EVENT_OK;

  m_n_evt_ana++;

  ///
  /// Analysis of detector hits
  ///
  for (int i_pl = 0; i_pl < N_PL; i_pl++) {
    for (auto it = list_hv[i_pl]->begin(); it != list_hv[i_pl]->end(); it++) {
      int    ele  = (*it)->get_element_id();
      double time = (*it)->get_tdc_time  ();
      m_h1_ele      [i_pl]->Fill(ele );
      m_h1_time     [i_pl]->Fill(time);
      m_h1_time_wide[i_pl]->Fill(time);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int ExtractTdcDist::End(PHCompositeNode* topNode)
{
  ostringstream oss;

  m_ofs << "N of events:\n"
        << "  All        = " << m_n_evt_all << "\n"
        << "  Trigger OK = " << m_n_evt_trig << "\n"
        << "  Analyzed   = " << m_n_evt_ana << "\n"
        << endl;

  TCanvas* c1 = new TCanvas("c1", "");
  c1->SetGrid();
  for (int i_pl = 0; i_pl < N_PL; i_pl++) {
    m_h1_nhit[i_pl]->Draw();
    oss.str("");
    oss << m_dir_out << "/" << m_h1_nhit[i_pl]->GetName() << ".png";
    c1->SaveAs(oss.str().c_str());
  
    m_h1_ele[i_pl]->Draw();
    oss.str("");
    oss << m_dir_out << "/" << m_h1_ele[i_pl]->GetName() << ".png";
    c1->SaveAs(oss.str().c_str());
  
    m_h1_time[i_pl]->Draw();
    oss.str("");
    oss << m_dir_out << "/" << m_h1_time[i_pl]->GetName() << ".png";
    c1->SaveAs(oss.str().c_str());

    m_h1_time_wide[i_pl]->Draw();
    oss.str("");
    oss << m_dir_out << "/" << m_h1_time_wide[i_pl]->GetName() << ".png";
    c1->SaveAs(oss.str().c_str());
  }
  delete c1;

  m_file->cd();
  m_file->Write();
  m_file->Close();
  m_ofs.close();
  
  return Fun4AllReturnCodes::EVENT_OK;
}
