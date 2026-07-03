/// Macro to draw the anticipated precision of TSSA of J/psi and psi'.
/**
 * One week of data taking is assumed.
 * delta_{AN}^{1w} = 1 / (df * pol * sqrt(R * L^1w))
 */
const double lumi_1w  = 7e4; // integrated luminosity in 1 week [/pb]
const double df       = 0.176; // dilution_factor
const double pol      = 0.8;   // polarization
TFile* file;
TCanvas* c1;
void DrawOneParticle(const char* part);
void DrawOneVar(const char* part, const char* var, const double x_min, const double x_max);

void DrawStatExp()
{
  cout << "Luminosity   = " << lumi_1w << " /pb for one week\n"
       << "Duty factor  = " << df << "\n"
       << "Polarization = " << pol << endl;
  c1 = new TCanvas("c1", "");
  c1->SetGrid();
  gSystem->mkdir("result_asym", true);
  DrawOneParticle("jpsi");
  DrawOneParticle("psip");
  exit(0);
}

void DrawOneParticle(const char* part)
{
  ostringstream oss;
  oss << part << "/result/hist.root";
  file = new TFile(oss.str().c_str());
  DrawOneVar(part, "x2", 0.0, 0.2);
  DrawOneVar(part, "x1", 0.2, 1.0);
  DrawOneVar(part, "xF", 0.2, 1.0);
  DrawOneVar(part, "pT", 0.0, 4.0);
  delete file;
}

void DrawOneVar(const char* part, const char* var, const double x_min, const double x_max)
{
  ostringstream oss;
  oss << "h1_" << var;
  TH1* h1 = (TH1*)file->Get(oss.str().c_str());
  h1->Draw("E1");
  oss.str("");
  oss << "result_asym/h1_" << part << "_" << var << "_input.png";
  c1->SaveAs(oss.str().c_str());

  h1->Scale(lumi_1w);
  h1->GetYaxis()->SetTitle("Yield per week");
  oss.str("");
  oss << "result_asym/h1_" << part << "_" << var << "_exp.png";
  c1->SaveAs(oss.str().c_str());

  int    n_pt = 0;
  double val_x[99];
  double val_y[99];
  double err_x[99];
  double err_y[99];
  for (int ib = 1; ib <= h1->GetNbinsX(); ib++) {
    double cont = h1->GetBinContent(ib);
    if (cont == 0) continue;
    double err = 1 / sqrt(cont) / df / pol;
    if (err > 0.3) continue;
    val_x[n_pt] = h1->GetBinCenter(ib);
    val_y[n_pt] = 0;
    err_x[n_pt] = h1->GetBinWidth(ib) / 2;
    err_y[n_pt] = err;
    n_pt++;
  }
  TGraphErrors* gr = new TGraphErrors(n_pt, val_x, val_y, err_x, err_y);
  gr->SetMarkerStyle(21);
  gr->SetMarkerColor(kRed);
  gr->SetLineColor(kRed);

  oss.str("");
  oss << ";" << var << ";Anticipated statistical presicion of A_{N}: #delta_{AN}";
  c1->DrawFrame(x_min, -0.3, x_max, 0.3, oss.str().c_str());
  gr->Draw("PE1same");

  oss.str("");
  oss << "result_asym/gr_" << part << "_tssa_err_" << var << ".png";
  c1->SaveAs(oss.str().c_str());

  delete gr;
}
