#include <iomanip>
#include <sstream>
#include <TH2D.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>
#include <THStack.h>
#include <TGraph.h>
#include <TLegend.h>
#include <geom_svc/GeomSvc.h>
#include <UtilAna/UtilHist.h>
#include <ktracker/FastTracklet.h>
#include "CalibParam.h"
#include "CalibData.h"
using namespace std;

CalibData::CalibData()
  : cal_par(0)
{
  ;
}

CalibData::~CalibData()
{
  ;
}

void CalibData::Init(CalibParam* ptr)
{
  cal_par = ptr;

  h1_stat  = new TH1D("h1_stat" ,      ";Reco. status;N of tracks",  25, -24.5, 0.5);

  int    id_n = ST_ID_MAX - ST_ID_MIN + 1;
  double id_l = ST_ID_MIN - 0.5;
  double id_h = ST_ID_MAX + 0.5;
  h1_st_id = new TH1D("h1_st_id",                   ";Station ID;N of tracks", id_n, id_l, id_h);
  h2_ntrk  = new TH2D("h2_ntrk" , ";N of tracks/event;Station ID;N of tracks", 200, -0.5, 199.5,  id_n, id_l, id_h);
  h2_nhit  = new TH2D("h2_nhit" ,   ";N of hits/track;Station ID;N of tracks",  18,  0.5,  18.5,  id_n, id_l, id_h);
  h2_mom   = new TH2D("h2_mom"  ,    ";Momentum (GeV);Station ID;N of tracks",  22,  0.0, 120.0,  id_n, id_l, id_h);
  h2_rchi2 = new TH2D("h2_rchi2",      ";#chi^{2}/NDF;Station ID;N of tracks",  50,  0.0,  10.0,  id_n, id_l, id_h);
  h2_x0    = new TH2D("h2_x0"   ,        ";x_{0} (cm);Station ID;N of tracks",  100,  -500, 500,  id_n, id_l, id_h);
  h2_y0    = new TH2D("h2_y0"   ,        ";y_{0} (cm);Station ID;N of tracks",  100,  -500, 500,  id_n, id_l, id_h);
  h2_tx    = new TH2D("h2_tx"   ,             ";t_{x};Station ID;N of tracks",  100,  -1.0, 1.0,  id_n, id_l, id_h);
  h2_ty    = new TH2D("h2_ty"   ,             ";t_{y};Station ID;N of tracks",  100,  -1.0, 1.0,  id_n, id_l, id_h);
  h2_x_st  = new TH2D("h2_x_st" ,       ";x_{St} (cm);Station ID;N of tracks",  100,  -250, 250,  id_n, id_l, id_h);
  h2_y_st  = new TH2D("h2_y_st" ,       ";y_{St} (cm);Station ID;N of tracks",  100,  -250, 250,  id_n, id_l, id_h);

  GeomSvc* geom = GeomSvc::instance();
  ostringstream oss;
  ostringstream oss2;
  oss2 << setfill('0');
  for (int ip = 0; ip < cal_par->GetNumPlanes(); ip++) {
    string det_name = geom->getDetectorName(ip+1);
    double det_dy   = geom->getPlaneScaleX (ip+1);
    double t_min = cal_par->GetTMin(ip);
    double t_max = cal_par->GetTMax(ip);
    double r_max = cal_par->GetRMax(ip);

    const double t_width = 4.0/9.0 * 3; // 4/9 ns * integer
    oss.str("");
    oss << "h1_time_wide_" << det_name;
    oss2.str("");
    oss2 << "#" << setw(2) << (ip+1) << ": " << det_name << ";TDC time (ns);Yields";
    h1_time_wide[ip] = new TH1D(oss.str().c_str(), oss2.str().c_str(),  2000, 0, 2000*t_width);

    int bin_t_lo = (int)(t_min / t_width);
    int bin_t_hi = (int)(t_max / t_width) + 1;
    oss.str("");
    oss << "h1_time_" << det_name;
    oss2.str("");
    oss2 << "#" << setw(2) << (ip+1) << ": " << det_name << ";TDC time (ns);Yields";
    h1_time[ip] = new TH1D(oss.str().c_str(), oss2.str().c_str(),  bin_t_hi-bin_t_lo, t_width*(bin_t_lo)-0.01, t_width*(bin_t_hi)-0.01);

    oss.str("");
    oss << "h2_tp_wp_" << det_name;
    oss2.str("");
    oss2 << "#" << setw(2) << (ip+1) << ": " << det_name << ";Track position (cm);Wire position (cm);";
    h2_tp_wp[ip] = new TH2D(oss.str().c_str(), oss2.str().c_str(),  120, -0.6*det_dy, 0.6*det_dy,  120, -0.6*det_dy, 0.6*det_dy);

    oss.str("");
    oss << "h2_rt_" << det_name;
    oss2.str("");
    oss2 << "#" << setw(2) << (ip+1) << ": " << det_name << ";TDC time (ns);R (cm);";
    h2_RT[ip] = new TH2D(oss.str().c_str(), oss2.str().c_str(),  T_BIN, t_min, t_max, R_BIN, 0, 1.2*r_max );

    oss.str("");
    oss << "h2_resi_" << det_name;
    oss2.str("");
    oss2 << "#" << setw(2) << (ip+1) << ": " << det_name << ";TDC time (ns);Residual (cm);";
    h2_Resi[ip] = new TH2D(oss.str().c_str(), oss2.str().c_str(),  T_BIN, t_min, t_max, R_BIN, -0.5, 0.5);

    oss.str("");
    oss << "h2_tr_" << det_name;
    oss2.str("");
    oss2 << "#" << setw(2) << (ip+1) << ": " << det_name << ";TDC time (ns);R (cm);";
    h2_TR[ip] = new TH2D(oss.str().c_str(), oss2.str().c_str(),  R_BIN, 0, 1.2*r_max, T_BIN, t_min, t_max );

    oss.str("");
    oss << "h2_tr_ex_" << det_name;
    oss2.str("");
    oss2 << "#" << setw(2) << (ip+1) << ": " << det_name << ";Track distance (cm);TDC time (ns);";
    h2_TR_ex[ip] = new TH2D(oss.str().c_str(), oss2.str().c_str(),  R_BIN*2, -1.2*r_max, 1.2*r_max, T_BIN, t_min, t_max );

    oss.str("");
    oss << "h2_td_dd_" << det_name;
    oss2.str("");
    oss2 << "#" << setw(2) << (ip+1) << ": " << det_name << ";Track distance (cm);Drift distance (cm);Yields";
    h2_td_dd[ip] = new TH2D(oss.str().c_str(), oss2.str().c_str(),  R_BIN*2, -1.2*r_max, 1.2*r_max,    R_BIN, 0, 1.2*r_max);
  }
}

void CalibData::FillEventInfo(const int rec_stat, const std::map<int, int> list_n_trk)
{
  h1_stat->Fill(rec_stat);

  for (int st_id = ST_ID_MIN; st_id <= ST_ID_MAX; st_id++) {
    int n_trk = list_n_trk.find(st_id) != list_n_trk.end()  ?  list_n_trk.at(st_id)  :  0;
    h2_ntrk->Fill(n_trk, st_id);
  }
}

void CalibData::FillTracklet(const Tracklet* trk)
{
  int st_id = trk->stationID;
  int n_hit = trk->getNHits();
  int ndf = n_hit - (st_id == 7  ?  5  :  4); // Incorrect, when KMag is off
  double rchi2 = ndf > 0  ?  trk->chisq / ndf  :  0;

  double z_st = CalibParam::ZOfStationID(st_id);
  double x_st = trk->x0 + z_st * trk->tx;
  double y_st = trk->y0 + z_st * trk->ty;

  h1_st_id->Fill(st_id);
  h2_nhit ->Fill(n_hit      , st_id);
  h2_mom  ->Fill(1/trk->invP, st_id);
  h2_rchi2->Fill(rchi2      , st_id);
  h2_x0   ->Fill(trk->x0    , st_id);
  h2_y0   ->Fill(trk->y0    , st_id);
  h2_tx   ->Fill(trk->tx    , st_id);
  h2_ty   ->Fill(trk->ty    , st_id);
  h2_x_st ->Fill(x_st       , st_id);
  h2_y_st ->Fill(y_st       , st_id);
}

void CalibData::FillTrackletHits(const Tracklet* trk)
{
  for(auto it = trk->hits.begin(); it != trk->hits.end(); ++it) {
    if(it->hit.index < 0) continue; // Probably not necessary.
    //int sign = it->hit.sign;
    int    det_id     = it->hit.detectorID;
    int    ele_id     = it->hit.elementID;
    double drift_dist = it->hit.driftDistance;
    double tdc_time   = it->hit.tdcTime;
    double track_pos  = trk->getExpPositionW(det_id);
    double wire_pos   = it->hit.pos; // GeomSvc::instance()->getMeasurement(det_id, ele_id);
    FillHit(det_id, drift_dist, tdc_time, track_pos, wire_pos);
    //cout << "D " << det_id << " " << it->hit.pos << " " << tdc_time << " " << drift_dist << " " << track_dist << endl;
  }
}

/// Fill the hit information.
/**
 * Here "track_dist" is the signed distance.
 */
void CalibData::FillHit(const int det_id, const double drift_dist, const double tdc_time, const double track_pos, const double wire_pos)
{
  double track_dist = track_pos - wire_pos;
  h1_time_wide[det_id-1]->Fill(tdc_time);
  h1_time     [det_id-1]->Fill(tdc_time);
  h2_tp_wp    [det_id-1]->Fill(track_pos, wire_pos);
  h2_RT       [det_id-1]->Fill(tdc_time, fabs(track_dist));
  h2_Resi     [det_id-1]->Fill(tdc_time, fabs(drift_dist) - fabs(track_dist));
  h2_TR       [det_id-1]->Fill(fabs(track_dist), tdc_time);
  h2_TR_ex    [det_id-1]->Fill(track_dist, tdc_time);
  h2_td_dd    [det_id-1]->Fill(track_dist, drift_dist);
}

void CalibData::DrawHistEvent(const string dir_out)
{
  cout << "DrawHistEvent()" << endl;
  ostringstream oss;
  TCanvas* c1 = new TCanvas("c1", "", 800, 600);
  c1->SetGrid();
  gStyle->SetOptStat(0);
  
  //h1_stat->Draw();
  //oss.str("");
  //oss << dir_out << "/h1_rec_stat.png";
  //c1->SaveAs(oss.str().c_str());

  h1_st_id->Draw();
  oss.str("");
  oss << dir_out << "/h1_st_id.png";
  c1->SaveAs(oss.str().c_str());
  
  const int iy_lo = 1; // 4 = D3p, 5 = D3m, 6 = D23
  const int iy_hi = 5;
  DrawIn1D(h2_ntrk , "ntrk" , dir_out, iy_lo, iy_hi);
  DrawIn1D(h2_nhit , "nhit" , dir_out, iy_lo, iy_hi);
  //DrawIn1D(h2_mom  , "mom"  , dir_out, iy_lo, iy_hi);
  DrawIn1D(h2_rchi2, "rchi2", dir_out, iy_lo, iy_hi);
  //DrawIn1D(h2_x0   , "x0"   , dir_out, iy_lo, iy_hi);
  //DrawIn1D(h2_y0   , "y0"   , dir_out, iy_lo, iy_hi);
  DrawIn1D(h2_tx   , "tx"   , dir_out, iy_lo, iy_hi);
  DrawIn1D(h2_ty   , "ty"   , dir_out, iy_lo, iy_hi);
  DrawIn1D(h2_x_st , "x_st" , dir_out, iy_lo, iy_hi);
  DrawIn1D(h2_y_st , "y_st" , dir_out, iy_lo, iy_hi);

  delete c1;
}

void CalibData::DrawHistHit(const string dir_out)
{
  cout << "DrawHistHit()" << endl;
  ostringstream oss;
  TCanvas* c1 = new TCanvas("c1", "", 800, 600);
  c1->SetGrid();
  gStyle->SetOptStat(0);
  
  GeomSvc* geom = GeomSvc::instance();
  
  oss << setfill('0');
  for (int ip = 0; ip < cal_par->GetNumPlanes(); ip++) {
    if (! cal_par->GetAnaPlane(ip)) continue;
    string det_name = geom->getDetectorName(ip+1);
    //cout << "  Plane " << ip+1 << endl;

    oss.str("");
    oss << setw(2) << (ip+1) << "_" << det_name << ".png";
    string suffix = oss.str();

    UtilHist::AutoSetRange(h1_time_wide[ip]);
    h1_time_wide[ip]->Draw();
    c1->SaveAs( (dir_out + "/h1_time_wide_" + suffix).c_str() );

    h1_time[ip]->Draw();
    c1->SaveAs( (dir_out + "/h1_time_" + suffix).c_str() );

    //h2_tp_wp[ip]->Draw("colz"); // To check if any outlier exists.
    //c1->SaveAs( (dir_out + "/h2_tp_wp_" + suffix).c_str() );

    TH1* h1_tp = h2_tp_wp[ip]->ProjectionX("h1_tp");
    h1_tp->Draw();
    c1->SaveAs( (dir_out + "/h1_tp_" + suffix).c_str() );
    delete h1_tp;

    //TH1* h1_wp = h2_tp_wp[ip]->ProjectionY("h1_wp");
    //h1_wp->Draw();
    //c1->SaveAs( (dir_out + "/h1_wp_" + suffix).c_str() );
    //delete h1_wp;
    
    h2_TR_ex[ip]->Draw("colz");
    c1->SaveAs( (dir_out + "/h2_tr_ex_" + suffix).c_str() );
    
    h2_Resi[ip]->Draw("colz");
    c1->SaveAs( (dir_out + "/h2_resi_" + suffix).c_str() );
    
    h2_td_dd[ip]->Draw("colz");
    c1->SaveAs( (dir_out + "/h2_td_dd_" + suffix).c_str() );
    
    TH1* h1_td = h2_td_dd[ip]->ProjectionX("h1_td");
    h1_td->Draw();
    c1->SaveAs( (dir_out + "/h1_td_" + suffix).c_str() );
    delete h1_td;
    
    TH1* h1_dd = h2_td_dd[ip]->ProjectionY("h1_dd");
    h1_dd->Draw();
    c1->SaveAs( (dir_out + "/h1_dd_" + suffix).c_str() );
    delete h1_dd;
  }
  delete c1;
}

void CalibData::DrawIn1D(TH2* h2, const string label, const string dir_out, int iy_lo, int iy_hi)
{
  if (iy_lo <= 0) iy_lo = 1;
  if (iy_hi <= 0) iy_hi = h2->GetNbinsY();

  ostringstream oss;
  oss << h2->GetTitle() << ";" << h2->GetXaxis()->GetTitle() << ";" << h2->GetZaxis()->GetTitle();
  THStack hs("hs", oss.str().c_str());
  TLegend leg(0.90, 0.90-0.03*(iy_hi-iy_lo+2), 0.99, 0.99);
  leg.SetHeader("St. ID");
  for (int iy = iy_lo; iy <= iy_hi; iy++) {
    oss.str("");
    oss << h2->GetName() << "_" << iy;
    TH1* h1 = h2->ProjectionX(oss.str().c_str(), iy, iy);
    h1->SetLineColor(iy - iy_lo + 1);
    hs.Add(h1);
    oss.str("");   
    //oss << h2->GetYaxis()->GetTitle() << h2->GetYaxis()->GetBinCenter(iy);
    oss << h2->GetYaxis()->GetBinCenter(iy);
    leg.AddEntry(h1, oss.str().c_str(), "l");
  }
  
  TCanvas* c1 = new TCanvas("c1_din1d", "", 800, 600);
  c1->SetGrid();
  c1->SetLogy();
  gStyle->SetOptStat(0);
  hs.Draw("nostack");
  leg.Draw();
  
  oss.str("");
  oss << dir_out << "/h1_" << label << ".png";
  c1->SaveAs(oss.str().c_str());
  delete c1;
}
