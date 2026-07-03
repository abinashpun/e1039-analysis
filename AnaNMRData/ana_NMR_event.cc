#include "NMRDataManager.h"

void ana_NMR_event()
{
  NMRDataManager man;
  //man.Verb(1); // 0 = silent, 1 = some, 2 = every event
  gSystem->mkdir("result", true);

  ///
  /// Read event file(s).
  ///
  man.ReadEventFile("2024-02-08_19h02m57s");
  man.ReadEventFile("2024-02-08_22h13m47s");
  man.PrintKeyList();

  ///
  /// Set the reference time.  Set 'utime0' to 0 to pick up the time of the 1st data.
  ///
  TTimeStamp ts0(2024, 2, 8, 18, 0, 0);
  int utime0 = ts0.GetSec() + TTimeStamp::GetZoneOffset();

  ///
  /// Create graphs/histograms and fill them in the event loop.
  ///
  TGraph* gr_pol  = new TGraph();
  TGraph* gr_pres = new TGraph();

  int n_evt = man.GetNumEvent();
  for (int i_evt = 0; i_evt < n_evt; i_evt++) {
    NMREvent* evt = man.GetEvent(i_evt);
    int   evt_num = evt->GetInt   ("EventNum");
    double    pol = evt->GetDouble("Polarization");
    double   area = evt->GetDouble("Area");
    string qmeter = evt->GetString("QMeterName");
    double   pres = evt->GetDouble("He4Press");
    //string YMDhms = NMRUtil::ConvEventNum(evt_num);

    if (utime0 == 0) utime0 = evt_num;
    if (evt_num < utime0) continue;

    int n_pt = gr_pol->GetN();
    double time = (evt_num - utime0) / 3600.0;
    gr_pol ->SetPoint(n_pt, time, pol );
    gr_pres->SetPoint(n_pt, time, pres);
  }

  ///
  /// Draw graphs/histograms.
  ///
  TCanvas* c1 = new TCanvas("c1", "");
  c1->SetGrid(true);

  gr_pol->SetTitle( (";Hours since "+NMRUtil::ConvEventNum(utime0)+";Polarization").c_str() );
  gr_pol->SetMarkerStyle(7);
  gr_pol->Draw("AP");
  c1->SaveAs("result/gr_pol.png");

  gr_pres->SetTitle( (";Hours since "+NMRUtil::ConvEventNum(utime0)+";Pressure (Torr)").c_str() );
  gr_pres->SetMarkerStyle(7);
  gr_pres->Draw("AP");
  gr_pres->GetYaxis()->SetRangeUser(21.0, 26.0);
  c1->SaveAs("result/gr_pres.png");

  delete c1;
  exit(0);
}
