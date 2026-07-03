void DrawMassDist()
{
  gSystem->mkdir("result_asym", true);
  TFile* file_jpsi = new TFile("jpsi/result/hist.root");
  TFile* file_psip = new TFile("psip/result/hist.root");
  TFile* file_dy   = new TFile("dy/result/hist.root");
  TH1* h1_jpsi = (TH1*)file_jpsi->Get("h1_mass_rec");
  TH1* h1_psip = (TH1*)file_psip->Get("h1_mass_rec");
  TH1* h1_dy   = (TH1*)file_dy  ->Get("h1_mass_rec");

  TH1* h1_dy_jpsi = (TH1*)h1_dy->Clone("h1_dy_jpsi");
  h1_dy_jpsi->Add(h1_jpsi);

  TH1* h1_dy_psip = (TH1*)h1_dy->Clone("h1_dy_psip");
  h1_dy_psip->Add(h1_psip);

  TH1* h1_all = (TH1*)h1_dy->Clone("h1_all");
  h1_all->Add(h1_jpsi);
  h1_all->Add(h1_psip);

  h1_jpsi   ->SetLineColor(kRed);
  h1_psip   ->SetLineColor(kBlue);
  h1_dy     ->SetLineColor(kBlack);
  h1_dy_jpsi->SetLineColor(kRed);
  h1_dy_psip->SetLineColor(kBlue);
  h1_all    ->SetLineColor(kGreen);
  h1_all    ->SetLineWidth(2);
  h1_all    ->SetTitle(";Reconstructed mass;Antcipated rate (pb)");

  TCanvas* c1 = new TCanvas("c1", "");
  c1->SetGrid();

  h1_all    ->Draw("HIST");
  h1_dy     ->Draw("HISTsame");
  h1_dy_psip->Draw("HISTsame");
  h1_dy_jpsi->Draw("HISTsame");

  c1->SaveAs("result_asym/h1_mass_tot.png");

  c1->SetLogy(true);
  c1->SaveAs("result_asym/h1_mass_tot_log.png");

  exit(0);
}
