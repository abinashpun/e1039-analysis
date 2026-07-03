#include <iomanip>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <interface_main/SQRun.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQHitVector.h>
#include <interface_main/SQTrackVector.h>
#include <ktracker/SRecEvent.h>
#include <ktracker/FastTracklet.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <geom_svc/GeomSvc.h>
#include <UtilAna/UtilSQHit.h>
#include "AnaTrack.h"
using namespace std;

AnaTrack::AnaTrack(const std::string& name)
  : SubsysReco(name)
  , m_evt         (0)
  , m_srec        (0)
  , m_vec_hit     (0)
  , m_vec_trk_true(0)
  , m_vec_trklet  (0)
  , m_file_out    (0)
//  , m_h1_ele      (0)
//  , m_h1_pos      (0)
//  , m_h1_time     (0)
//  , m_h1_dd       (0)
//  , m_h1_td       (0)
{
  ;
}

int AnaTrack::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTrack::InitRun(PHCompositeNode* topNode)
{
  m_evt          = findNode::getClass<SQEvent       >(topNode, "SQEvent");
  m_vec_hit      = findNode::getClass<SQHitVector   >(topNode, "SQHitVector");
  m_vec_trk_true = findNode::getClass<SQTrackVector >(topNode, "SQTruthTrackVector");
  m_srec         = findNode::getClass<SRecEvent     >(topNode, "SRecEvent");
  m_vec_trklet   = findNode::getClass<TrackletVector>(topNode, "TrackletVector");
  if (!m_evt || !m_vec_hit || !m_vec_trk_true || !m_srec || !m_vec_trklet) return Fun4AllReturnCodes::ABORTEVENT;

  GeomSvc* geom = GeomSvc::instance();

  gSystem->mkdir(Name().c_str(), true);
  m_file_out = new TFile(TString::Format("%s/output.root", Name().c_str()).Data(), "RECREATE");

  m_h1_ntrk  = new TH1D("h1_ntrk" , ";N of tracks/event;N of tracks", 100, -0.5, 199.5);
  m_h1_nhit  = new TH1D("h1_nhit" ,   ";N of hits/track;N of tracks",  18,  0.5,  18.5);
  m_h1_mom   = new TH1D("h1_mom"  ,    ";Momentum (GeV);N of tracks", 120,  0.0, 120.0);
  m_h1_rchi2 = new TH1D("h1_rchi2",      ";#chi^{2}/NDF;N of tracks", 100,  0.0,   3.0);
  m_h1_x0    = new TH1D("h1_x0"   ,        ";x_{0} (cm);N of tracks", 100, -100,   100);
  m_h1_y0    = new TH1D("h1_y0"   ,        ";y_{0} (cm);N of tracks", 100, -100,   100);
  m_h1_tx    = new TH1D("h1_tx"   ,             ";t_{x};N of tracks", 100, -0.1,   0.1);
  m_h1_ty    = new TH1D("h1_ty"   ,             ";t_{y};N of tracks", 100, -0.1,   0.1);

  for (int i_pl = 1; i_pl <= N_PL; i_pl++) {
    string name_pl     = geom->getDetectorName(i_pl);
    int    n_ele       = geom->getPlaneNElements(i_pl);
    double plane_width = geom->getPlaneScaleX(i_pl);
    double cell_width  = geom->getCellWidth(i_pl);

    m_h1_ele[i_pl] = new TH1D(
      TString::Format("h1_ele_%d_%s", i_pl, name_pl.c_str()).Data(), 
      TString::Format("#%d: %s;Element ID;Yield", i_pl, name_pl.c_str()).Data(), 
      n_ele, 0.5, n_ele+0.5);
    
    m_h1_pos[i_pl] = new TH1D(
      TString::Format("h1_pos_%d_%s", i_pl, name_pl.c_str()).Data(), 
      TString::Format("#%d: %s;Position (cm);Yield", i_pl, name_pl.c_str()).Data(), 
      60, -0.6*plane_width, 0.6*plane_width);
    
    m_h1_time[i_pl] = new TH1D(
      TString::Format("h1_time_%d_%s", i_pl, name_pl.c_str()).Data(), 
      TString::Format("#%d: %s;tdcTime;Yield", i_pl, name_pl.c_str()).Data(), 
      150, 500, 2000);

    m_h1_dd[i_pl] = new TH1D(
      TString::Format("h1_dd_%d_%s", i_pl, name_pl.c_str()).Data(), 
      TString::Format("#%d: %s;Drift dist;Yield", i_pl, name_pl.c_str()).Data(), 
      120, 0.0, 0.6 * cell_width);

    m_h1_td[i_pl] = new TH1D(
      TString::Format("h1_td_%d_%s", i_pl, name_pl.c_str()).Data(), 
      TString::Format("#%d: %s;Track dist;Yield", i_pl, name_pl.c_str()).Data(), 
      120, 0.0, 0.6 * cell_width);
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTrack::process_event(PHCompositeNode* topNode)
{
  int n_trk = m_vec_trklet->size();
  m_h1_ntrk->Fill(n_trk);
  if (Verbosity() >= 2) {
    cout << "AnaTrack:\n"
         << "  N of true tracks = " << m_vec_trk_true->size() << "\n"
         << "  N of SRec tracks = " << m_srec->getNTracks() << "\n"
         << "  N of tracklets   = " << n_trk << endl;
  }

  for (int i_trk = 0; i_trk < n_trk; i_trk++) {
    Tracklet* trk = m_vec_trklet->at(i_trk);
    if (Verbosity() >= 2) {
      cout << "  Track " << i_trk << ": stationID = " << trk->stationID << endl;
    }

    int n_hit = trk->getNHits();
    int ndf = n_hit - 4; // Correct only when KMag is off
    double rchi2 = trk->chisq / ndf;
    m_h1_nhit ->Fill(n_hit      );
    m_h1_mom  ->Fill(1/trk->invP);
    m_h1_rchi2->Fill(rchi2      );
    m_h1_x0   ->Fill(trk->x0    );
    m_h1_y0   ->Fill(trk->y0    );
    m_h1_tx   ->Fill(trk->tx    );
    m_h1_ty   ->Fill(trk->ty    );

    for (std::list<SignedHit>::iterator it = trk->hits.begin(); it != trk->hits.end(); ++it) {
      if(it->hit.index < 0) continue;
      int    det_id     = it->hit.detectorID;
      int    ele_id     = it->hit.elementID;
      double pos        = it->hit.pos;
      double tdc_time   = it->hit.tdcTime;
      double drift_dist = it->hit.driftDistance;
      double track_dist = trk->getExpPositionW(det_id) - pos; // track position - wire position
      //if (det_id == 1) cout << " Z " << trk->getExpPositionW(det_id) << " " << it->hit.pos << " " << r << " " << t << " " << TMAX[det_id-1]-t << endl;
      m_h1_ele [det_id]->Fill(ele_id);
      m_h1_pos [det_id]->Fill(pos);
      m_h1_time[det_id]->Fill(tdc_time);
      m_h1_dd  [det_id]->Fill(drift_dist);
      m_h1_td  [det_id]->Fill(track_dist);
    }
  }

  //for (SQHitVector::ConstIter it = hv->begin(); it != hv->end(); it++) {
  //  int    ele  = (*it)->get_element_id();
  //  double pos  = (*it)->get_pos();
  //  double time = (*it)->get_tdc_time();
  //  double dist = (*it)->get_drift_distance();
  //  m_h1_ele ->Fill(ele );
  //  m_h1_pos ->Fill(pos );
  //  m_h1_time->Fill(time);
  //  m_h1_dist->Fill(dist);
  //}

//  IdMap_t id_trk_t2r; // [idx_true] -> idx_reco
//  FindTrackRelation(id_trk_t2r);
  for (unsigned int i_true = 0; i_true < m_vec_trk_true->size(); i_true++) {
    //if (id_trk_t2r[i_true] == 0) continue;
    
    SQTrack* trk_true = m_vec_trk_true->at(i_true);
    int charge = trk_true->get_charge();
    TVector3 pos_vtx = trk_true->get_pos_vtx();
    TVector3 pos_st1 = trk_true->get_pos_st1();
    TVector3 pos_st3 = trk_true->get_pos_st3();
    TLorentzVector mom_vtx = trk_true->get_mom_vtx();
    TLorentzVector mom_st1 = trk_true->get_mom_st1();
    TLorentzVector mom_st3 = trk_true->get_mom_st3();
    if (Verbosity() >= 3) {
      cout << "  True track #" << i_true << ": " << charge << "\n"
           << "    pos (" << pos_vtx.X() << ", " << pos_vtx.Y() << " " << pos_vtx.Z() << ")"
           <<        " (" << pos_st1.X() << ", " << pos_st1.Y() << " " << pos_st1.Z() << ")"
           <<        " (" << pos_st3.X() << ", " << pos_st3.Y() << " " << pos_st3.Z() << ")\n"
           << "    mom (" << mom_vtx.X() << ", " << mom_vtx.Y() << ", " << mom_vtx.Z() << ")"
           <<        " (" << mom_st1.X() << ", " << mom_st1.Y() << ", " << mom_st1.Z() << ")"
           <<        " (" << mom_st3.X() << ", " << mom_st3.Y() << ", " << mom_st3.Z() << ")"
           << endl;
    }
//  
//      TrackData tdr;
//      SRecTrack* trk_reco = &m_srec->getTrack(id_trk_t2r[i_true]);
//      tdr.charge  = trk_reco->getCharge();
//      tdr.pos_vtx = trk_reco->getVertex();
//      tdr.mom_vtx = trk_reco->getMomentumVertex();
//
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTrack::End(PHCompositeNode* topNode)
{
  DrawTrackPar();

  DrawOneGroup("ele_d0" , m_h1_ele,  1,  6);
  DrawOneGroup("ele_d2" , m_h1_ele, 13, 18);
  DrawOneGroup("ele_d3p", m_h1_ele, 19, 24);
  DrawOneGroup("ele_d3m", m_h1_ele, 25, 30);

  DrawOneGroup("pos_d0" , m_h1_pos,  1,  6);
  DrawOneGroup("pos_d2" , m_h1_pos, 13, 18);
  DrawOneGroup("pos_d3p", m_h1_pos, 19, 24);
  DrawOneGroup("pos_d3m", m_h1_pos, 25, 30);

  DrawOneGroup("time_d0" , m_h1_time,  1,  6);
  DrawOneGroup("time_d2" , m_h1_time, 13, 18);
  DrawOneGroup("time_d3p", m_h1_time, 19, 24);
  DrawOneGroup("time_d3m", m_h1_time, 25, 30);

  DrawOneGroup("dd_d0" , m_h1_dd,  1,  6);
  DrawOneGroup("dd_d2" , m_h1_dd, 13, 18);
  DrawOneGroup("dd_d3p", m_h1_dd, 19, 24);
  DrawOneGroup("dd_d3m", m_h1_dd, 25, 30);

  DrawOneGroup("td_d0" , m_h1_td,  1,  6);
  DrawOneGroup("td_d2" , m_h1_td, 13, 18);
  DrawOneGroup("td_d3p", m_h1_td, 19, 24);
  DrawOneGroup("td_d3m", m_h1_td, 25, 30);

  //std::ofstream ofs(TString::Format("%s/output.txt", Name().c_str()).Data());
  //ofs << "N of events = " << m_n_evt_all << "\n";
  //ofs.close();

  m_file_out->cd();
  m_file_out->Write();
  m_file_out->Close();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void AnaTrack::FindTrackRelation(IdMap_t& id_map)
{
  id_map.clear();
  for (unsigned int i_true = 0; i_true < m_vec_trk_true->size(); i_true++) {
    SQTrack* trk_true = m_vec_trk_true->at(i_true);
    int     ch_true = trk_true->get_charge();
    double mom_true = trk_true->get_mom_vtx().Mag();

    int i_reco_best = -1;
    double mom_diff_best;
    for (int i_reco = 0; i_reco < m_srec->getNTracks(); i_reco++) {
      SRecTrack* trk_reco = &m_srec->getTrack(i_reco);
      if (trk_reco->getCharge() != ch_true) continue;
      double mom_diff = fabs(trk_reco->getMomentumVertex().Mag() - mom_true);
      if (i_reco_best < 0 || mom_diff < mom_diff_best) {
        i_reco_best   = i_reco;
        mom_diff_best = mom_diff;
      }
    }
    id_map[i_true] = i_reco_best;
  }
}

void AnaTrack::DrawTrackPar()
{
  TCanvas* c1 = new TCanvas("c1", "");
  c1->SetGrid();
  
  m_h1_nhit ->Draw();  c1->SaveAs((Name()+"/h1_track_nhit.png" ).c_str());
  m_h1_mom  ->Draw();  c1->SaveAs((Name()+"/h1_track_mom.png"  ).c_str());
  m_h1_rchi2->Draw();  c1->SaveAs((Name()+"/h1_track_rchi2.png").c_str());
  m_h1_x0   ->Draw();  c1->SaveAs((Name()+"/h1_track_x0.png"   ).c_str());
  m_h1_y0   ->Draw();  c1->SaveAs((Name()+"/h1_track_y0.png"   ).c_str());
  m_h1_tx   ->Draw();  c1->SaveAs((Name()+"/h1_track_tx.png"   ).c_str());
  m_h1_ty   ->Draw();  c1->SaveAs((Name()+"/h1_track_ty.png"   ).c_str());

  delete c1;
}

void AnaTrack::DrawOneGroup(const std::string name, TH1* h1[], const int i_pl0, const int i_pl1)
{
  ostringstream oss;
  oss << ";" << h1[i_pl0]->GetXaxis()->GetTitle() << ";" << h1[i_pl0]->GetYaxis()->GetTitle();
  THStack* hs  = new THStack("hs" , oss.str().c_str());
  TLegend* leg = new TLegend(0.85, 0.70, 0.99, 0.99);

  for (int i_pl = i_pl0; i_pl <= i_pl1; i_pl++) {
    h1[i_pl]->SetLineColor(i_pl - i_pl0 + 1);
    hs->Add(h1[i_pl]);
    leg->AddEntry(h1[i_pl], h1[i_pl]->GetTitle(), "l");
  }

  TCanvas* c1 = new TCanvas("c1", "");
  c1->SetGrid();
  hs->Draw("nostack");
  leg->Draw();

  oss.str("");
  oss << Name() << "/h1_" << name << ".png";
  c1->SaveAs(oss.str().c_str());

  delete c1;
  delete leg;
  delete hs;
}
