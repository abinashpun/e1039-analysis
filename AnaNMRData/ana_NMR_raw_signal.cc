#include "NMRDataManager.h"

void ana_NMR_raw_signal()
{
  NMRDataManager man;
  man.Verb(0); // 0 = silent, 1 = some, 2 = every event
  gSystem->mkdir("result", true);

  ///
  /// Read event file(s).
  ///
  //man.ReadEventFile("2024-02-08_19h02m57s");
  //man.ReadEventFile("2024-02-08_22h13m47s");
  man.ReadEventFile("2024-02-12_16h05m57s");
  //man.ReadEventFile("2024-02-13_19h58m19s");

  int ScanSweeps; // N of sweeps per measurement
  int ScanSteps;  // N of frequency steps per scan
  double RFFreq;  // Center frequency in MHz
  double RFMod;   // Half width of frequency range in kHz
  man.GetScanParam(ScanSweeps, ScanSteps, RFFreq, RFMod);

  ///
  /// Set the reference time.  Set 'utime0' to 0 to pick up the time of the 1st data.
  ///
  //TTimeStamp ts0(2024, 2, 8, 18, 0, 0);
  TTimeStamp ts0(2024, 2, 11,  18, 0, 0); // For "2024-02-12_16h05m57s"
  //TTimeStamp ts0(2024, 2, 12,  18, 0, 0); // For "2024-02-12_16h05m57s"
  int utime0 = ts0.GetSec() + TTimeStamp::GetZoneOffset();

  ///
  /// Create graphs/histograms and fill them in the event loop.
  ///
  TGraph* gr_first = new TGraph();
  TGraph* gr_last  = new TGraph();
  TGraph* gr_diff  = new TGraph();

  int n_evt = man.GetNumEvent();
  int i_evt_ok = 0;
  for (int i_evt = 0; i_evt < n_evt; i_evt++) {
    NMREvent* evt = man.GetEvent(i_evt);
    int   evt_num = evt->GetInt("EventNum");
    int      n_pt = evt->GetInt("ScanSteps");
    if (utime0 == 0) utime0 = evt_num;
    if (evt_num < utime0) continue;

    NMRSignal* sig = man.GetRawSignal(evt_num);
    double val0 = sig->GetPoint(0);
    double val1 = sig->GetPoint(n_pt-1);
    double hour = (evt_num - utime0) / 3600.0;
    gr_first->SetPoint(gr_first->GetN(), hour, val0);
    gr_last ->SetPoint(gr_last ->GetN(), hour, val1);
    gr_diff ->SetPoint(gr_diff ->GetN(), hour, val0 - val1);
  }

  ///
  /// Draw graphs/histograms.
  ///
  TCanvas* c1 = new TCanvas("c1", "", 1600, 600);
  c1->SetGrid(true);
  c1->SetMargin(0.05, 0.05, 0.1, 0.1); // (l, r, b, t)
  TGaxis::SetMaxDigits(3);

  gr_last ->SetMarkerStyle(7);
  gr_last ->SetMarkerColor(kBlue);
  gr_first->SetMarkerStyle(7);
  gr_first->SetMarkerColor(kRed);

  TMultiGraph* mg = new TMultiGraph();
  string title_x = "Hours since " + NMRUtil::ConvEventNum(utime0);
  mg->SetTitle((";" + title_x + ";RawSignal").c_str());
  mg->Add(gr_first, "P");
  mg->Add(gr_last , "P");
  mg->Draw("A");
  mg->GetYaxis()->SetTitleOffset(0.6);

  TLegend* leg = new TLegend(0.85, 0.90, 0.99, 0.99);
  leg->AddEntry(gr_first, "First point", "p");
  leg->AddEntry(gr_last ,  "Last point", "p");
  leg->Draw();

  c1->SaveAs("result/gr_voltage_vs_time.png");

  gr_diff->SetTitle((";" + title_x + ";V_{First} - V_{Last}").c_str());
  gr_diff->GetYaxis()->SetTitleOffset(0.6);
  gr_diff->SetMarkerStyle(7);
  gr_diff->SetMarkerColor(kRed);
  gr_diff->Draw("AP");
  c1->SaveAs("result/gr_diff_vs_time.png");

  delete c1;
  exit(0);
}
