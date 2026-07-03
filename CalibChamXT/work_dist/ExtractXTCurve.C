/// ExtractXTCurve.C:  ROOT macro to extract the X-T curve.
R__LOAD_LIBRARY(libCalibChamXT)

int ExtractXTCurve(const string cham = "d0", const int run=3824)
{
  recoConsts* rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", run);

  int det_id_1st;
  double T1, T0;
  if      (cham == "d0" ) { det_id_1st =  1; T1 = 1130; T0 = 1330; }
  else if (cham == "d1" ) { det_id_1st =  7; T1 = 1700; T0 = 1800; }
  else if (cham == "d2" ) { det_id_1st = 13; T1 =  850; T0 = 1300; }
  else if (cham == "d3p") { det_id_1st = 19; T1 =  740; T0 = 1300; }
  else if (cham == "d3m") { det_id_1st = 25; T1 =  850; T0 = 1300; }
  else {
    cout << "Bad chamber name (" << cham << ").  Abort." << endl;
    exit(1);
  }

  string fname = "tdc_dist_"+cham+"/output.root";
  TFile* file = new TFile(fname.c_str());

  TCanvas* c1 = new TCanvas("c1", "");
  c1->SetGrid();

  string dir_out = "result_" + cham;
  gSystem->mkdir(dir_out.c_str(), true);

  ofstream ofs( (dir_out + "/xt_curve.tsv").c_str() );

  ostringstream oss;
  oss << setfill('0');
  for (int det_id = det_id_1st; det_id < det_id_1st+6; det_id++) {
    oss.str("");
    oss << "h1_time_" << setw(2) << det_id;
    TH1* h1 = (TH1*)file->Get(oss.str().c_str());
    if (! h1) {
      cout << "Cannot find '" << oss.str() << "'.  Abort." << endl;
      exit(1);
    }

    string det_name = h1->GetTitle();
    double x_max = GeomSvc::instance()->getCellWidth(det_id) / 2;
    cout << det_id << " " << det_name << " " << x_max << endl;

    TGraph* gr = new TGraph();

    int bin1 = h1->GetXaxis()->FindBin(T1);
    int bin0 = h1->GetXaxis()->FindBin(T0);
    double cont_tot = h1->Integral(bin1, bin0);
    for (int bin = bin1; bin <= bin0; bin++) {
      double cont = h1->Integral(bin1, bin);
      double x = x_max * (1 - cont / cont_tot);
      double t = h1->GetXaxis()->GetBinCenter(bin);
      //cout << t << "\t" << x << "\n";
      gr->SetPoint(bin - bin1, t, x);
    }

    oss.str("");
    oss << det_name << ";T (ns);X (cm)";
    gr->SetTitle(oss.str().c_str());
    gr->SetMarkerStyle(7);
    gr->SetMarkerColor(kRed);
    gr->SetLineColor  (kRed);
    gr->Draw("ACP");

    oss.str("");
    oss << dir_out << "/gr_xt_" << setw(2) << det_id << ".png";
    c1->SaveAs(oss.str().c_str());

    const int T_STEP = 2;
    for (int t = T1; t <= T0; t += T_STEP) {
      double x  = gr->Eval(t);
      double dx = 0.03;
      double dt = dx * (T0 - T1) / x_max / 2; // Rough estimate
      ofs << det_name << "\t" << t << "\t" << x << "\t" << dt << "\t" << dx << "\n";
    }
    
    delete gr;
  }

  delete c1;
  ofs.close();
  exit(0);
}
