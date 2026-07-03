#include <sstream>
#include <iomanip>
#include <TSystem.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <interface_main/SQRun.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQHardEvent.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <geom_svc/GeomSvc.h>
//#include <UtilAna/UtilSQHit.h>
#include "AnaHardEvent.h"
using namespace std;

AnaHardEvent::AnaHardEvent(const std::string& name)
  : SubsysReco(name)
  , m_file(0)
  , m_h1_n_ttdc(0)
{
  ;
}

int AnaHardEvent::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaHardEvent::InitRun(PHCompositeNode* topNode)
{
  ///
  /// Input
  ///
  m_evt      = findNode::getClass<SQEvent    >(topNode, "SQEvent");
  m_hard_evt = findNode::getClass<SQHardEvent>(topNode, "SQHardEvent");
  if (!m_evt || !m_hard_evt) return Fun4AllReturnCodes::ABORTEVENT;

  ///
  /// Output
  ///
  gSystem->mkdir("result", true);

  m_ofs.open("result/log.txt");

  m_file = new TFile("result/output.root", "RECREATE");
  m_h1_n_ttdc = new TH1D("h1_n_ttdc", ";N of Taiwan TDCs", 100, -0.5, 99.5);

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaHardEvent::process_event(PHCompositeNode* topNode)
{
  int sp_id  = m_evt->get_spill_id();
  int evt_id = m_evt->get_event_id();
  //if (! m_evt->get_trigger(SQEvent::NIM2)) {
  //  return Fun4AllReturnCodes::EVENT_OK;
  //}

  int coda_evt_id = m_hard_evt->get_coda_event_id();
  short n_ttdc    = m_hard_evt->get_n_board_taiwan();
  //short n_v1495   = m_hard_evt->get_n_board_v1495();
  if (evt_id == 27294) cout << "E " << evt_id << " C " << coda_evt_id << " S " << sp_id << " N " << n_ttdc << endl;

  m_h1_n_ttdc->Fill(n_ttdc);
  if (n_ttdc != 82) m_ofs << evt_id << "\t" << coda_evt_id << "\t" << n_ttdc << "\n";

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaHardEvent::End(PHCompositeNode* topNode)
{
  ostringstream oss;
  TCanvas* c1 = new TCanvas("c1", "");
  c1->SetGrid();

  m_h1_n_ttdc->Draw();
  c1->SetLogy(true);
  c1->SaveAs("result/h1_n_ttdc.png");

  delete c1;
  m_file->cd();
  m_file->Write();
  m_file->Close();
  m_ofs.close();
  
  return Fun4AllReturnCodes::EVENT_OK;
}
