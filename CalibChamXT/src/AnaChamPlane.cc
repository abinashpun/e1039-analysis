#include <iomanip>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <THStack.h>
#include <TEfficiency.h>
#include <interface_main/SQRun.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQHitVector.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <geom_svc/GeomSvc.h>
#include <UtilAna/UtilSQHit.h>
#include "AnaChamPlane.h"
using namespace std;

AnaChamPlane::AnaChamPlane(const std::string& plane_name)
  : SubsysReco("AnaChamPlane"+plane_name)
  , m_plane_name(plane_name)
  , m_plane_id  (0)
  , m_n_ele     (0)
  , m_n_evt_all (0)
  , m_n_hit_all (0)
  , m_evt       (0)
  , m_hit_vec   (0)
  , m_file_out  (0)
  , m_h1_nhit   (0)
  , m_h1_ele    (0)
  , m_h1_pos    (0)
  , m_h1_time   (0)
  , m_h1_dist   (0)
{
  ;
}

int AnaChamPlane::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaChamPlane::InitRun(PHCompositeNode* topNode)
{
  m_evt     = findNode::getClass<SQEvent    >(topNode, "SQEvent");
  m_hit_vec = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if (!m_evt || !m_hit_vec) return Fun4AllReturnCodes::ABORTEVENT;

  GeomSvc* geom = GeomSvc::instance();
  m_plane_id = geom->getDetectorID(m_plane_name);
  if (m_plane_id <= 0) {
    cout << "!!ERROR!!  " << Name() << " is not given a proper plane name.\n";
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  m_n_ele = geom->getPlaneNElements(m_plane_id);

  m_file_out = new TFile(TString::Format("output_%s.root", Name().c_str()).Data(), "RECREATE");

  m_h1_nhit = new TH1D(
    TString::Format("h1_nhit_%d_%s", m_plane_id, m_plane_name.c_str()).Data(), 
    TString::Format("#%d: %s;N of hits/plane/event;Hit count", m_plane_id, m_plane_name.c_str()).Data(), 
    m_n_ele+1, -0.5, m_n_ele+0.5);
  
  m_h1_ele = new TH1D(
    TString::Format("h1_ele_%d_%s", m_plane_id, m_plane_name.c_str()).Data(), 
    TString::Format("#%d: %s;Element ID;Hit count", m_plane_id, m_plane_name.c_str()).Data(), 
    m_n_ele, 0.5, m_n_ele+0.5);
  
  double plane_width = geom->getPlaneScaleX(m_plane_id);
  m_h1_pos = new TH1D(
    TString::Format("h1_pos_%d_%s", m_plane_id, m_plane_name.c_str()).Data(), 
    TString::Format("#%d: %s;Position (cm);Hit count", m_plane_id, m_plane_name.c_str()).Data(), 
    60, -0.6*plane_width, 0.6*plane_width);
  
  const double DT = 20/9.0; // 4/9 ns per single count of Taiwan TDC
  const int NT = 450;
  const double T0 = 225.5*DT;
  const double T1 = 675.5*DT;
  m_h1_time = new TH1D(
    TString::Format("h1_time_%d_%s", m_plane_id, m_plane_name.c_str()).Data(), 
    TString::Format("#%d: %s;tdcTime;Hit count", m_plane_id, m_plane_name.c_str()).Data(), 
    NT, T0, T1);

  m_h1_dist = new TH1D(
    TString::Format("h1_dist_%d_%s", m_plane_id, m_plane_name.c_str()).Data(), 
    TString::Format("#%d: %s;Drift dist;Hit count", m_plane_id, m_plane_name.c_str()).Data(), 
    120, 0.0, 0.6 * geom->getCellWidth(m_plane_id));
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaChamPlane::process_event(PHCompositeNode* topNode)
{
  //int run_id   = m_evt->get_run_id();
  //int spill_id = m_evt->get_spill_id();
  //int event_id = m_evt->get_event_id();
  m_n_evt_all++;

  shared_ptr<SQHitVector> hv(UtilSQHit::FindHits(m_hit_vec, m_plane_id));
  m_n_hit_all += hv->size();

  m_h1_nhit->Fill(hv->size());
  for (SQHitVector::ConstIter it = hv->begin(); it != hv->end(); it++) {
    int    ele  = (*it)->get_element_id();
    double pos  = (*it)->get_pos();
    double time = (*it)->get_tdc_time();
    double dist = (*it)->get_drift_distance();
    m_h1_ele ->Fill(ele );
    m_h1_pos ->Fill(pos );
    m_h1_time->Fill(time);
    m_h1_dist->Fill(dist);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaChamPlane::End(PHCompositeNode* topNode)
{
  std::ofstream ofs(TString::Format("output_%s.txt", Name().c_str()).Data());
  ofs << "N of events = " << m_n_evt_all << "\n"
      << "N of hits   = " << m_n_hit_all << "\n";
  ofs.close();

  m_file_out->cd();
  m_file_out->Write();
  m_file_out->Close();
  
  return Fun4AllReturnCodes::EVENT_OK;
}
