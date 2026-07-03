using namespace std;
const double PI = acos(-1);
const int N_PHI = 20;
const int N_X2  = 10;
const double X2_MIN = 0.00;
const double X2_MAX = 0.20;

void asymmetry(const char* fname="sim_tree.root", const char* tname="tree")
{
  TFile* file = new TFile(fname);
  TTree* tree = (TTree*)file->Get(tname);

  const double PI = TMath::Pi();
  TH1D *h1_x2  = new TH1D("h1_x2", "GMC J/psi;x2;Yield", N_X2, X2_MIN, X2_MAX); 
  //TH1D *h1_m   = new TH1D("h1_m","",10, 0,10);
  TH2D *h2_phi_x2 = new TH2D("h2_phi_x2", "GMC J/psi;#phi;x2;Yield", N_PHI, -PI, PI,  N_X2, X2_MIN, X2_MAX);

  tree->Project("h1_x2"    , "dim_true.x2"               , "dim_true.pdg_id==443");
  //tree->Project("h1_m"     , "dim_true.mom.M()"          , "dim_true.pdg_id==443");
  tree->Project("h2_phi_x2", "dim_true.x2:dim_true.phi_s", "dim_true.pdg_id==443");
  //tree->Project("h2_phi_x2", "dim_true.x2:dim_true.phi_s_up", "dim_true.pdg_id==443");

  gSystem->mkdir("result_asym", true);
  ofstream ofs("result_asym/result.txt");

  TCanvas* c1 = new TCanvas("c1","");
  c1->SetGrid();
  h1_x2->Draw();
  c1->SaveAs("result_asym/h1_x2.png");
  //h1_m->Draw();
  //c1->SaveAs("result_asym/h1_m.png");
  h2_phi_x2->Draw("colz");
  c1->SaveAs("result_asym/h2_phi_x2.png");

  ostringstream oss;
  oss << setfill('0');
  TH1* h1_asym = new TH1D("h1_asym", "GMC J/psi;x2;A_{N}", N_X2, X2_MIN, X2_MAX);//Sivers asymetry An
  for (int i_x2 = 1; i_x2 <= N_X2; i_x2++) {
    TH1* h1_phi = h2_phi_x2->ProjectionX("h1_phi", i_x2, i_x2);
    double cnt_tot = h1_phi->Integral();
    h1_phi->Draw();
    h1_phi->SetMinimum(0);
    oss.str("");
    oss << h1_phi->GetTitle() << " : x2 = " << h2_phi_x2->GetYaxis()->GetBinCenter(i_x2);
    h1_phi->SetTitle(oss.str().c_str());
    oss.str("");
    oss << "result_asym/h1_phi_" << setw(2) << i_x2 << ".png";
    c1->SaveAs(oss.str().c_str());
    delete h1_phi;

    double x2 = X2_MIN + (X2_MAX - X2_MIN) * (i_x2 - 0.5) / N_X2;
    double x2_width = (X2_MAX - X2_MIN) / N_X2;
    float sum   = 0;
    float sum_N = 0;
    //float new_err=0;
    //for (int i_phi = 1; i_phi <= N_PHI; i_phi++) {
    //  float phi_angle = -PI + (i_phi - 0.5)*2.0*PI/N_PHI;
    //  float a = h2_phi_x2->GetBinContent(i_phi, i_x2);
    //  float b = h2_phi_x2->GetBinError  (i_phi, i_x2);
    //  sum     += a*sin(phi_angle);
    //  sum_N   += a;
    //  new_err += pow(b * (sin(phi_angle) * sum_N - sum) / pow(sum_N, 2), 2);
    //}
    //double asym_val = 0;
    //double asym_err = 0;
    //if (sum_N > 0) {
    //  asym_val = 2 * sum / sum_N;
    //  asym_err = TMath::Sqrt(new_err);
    //  h1_asym->SetBinContent(i_x2, asym_val);
    //  h1_asym->SetBinError  (i_x2, asym_err);
    //}

    bool stat_ok = true;
    for (int i_phi = 1; i_phi <= N_PHI; i_phi++) {
      float phi = -PI + (i_phi - 0.5)*2.0*PI/N_PHI;
      float cont = h2_phi_x2->GetBinContent(i_phi, i_x2);
      sum     += cont*sin(phi);
      sum_N   += cont;
//      if (cont <= 0) stat_ok = false;
    }
    if (sum == 0) stat_ok = false;
    double asym_val = 0;
    double asym_err = 0;
    if (stat_ok > 0) {
      asym_val = 2 * sum / sum_N;
      for (int i_phi = 1; i_phi <= N_PHI; i_phi++) {
        float phi = -PI + (i_phi - 0.5)*2.0*PI/N_PHI;
        float err = h2_phi_x2->GetBinError(i_phi, i_x2);
        asym_err += pow(err * (sin(phi) * sum_N - sum) / pow(sum_N, 2), 2);
      }
      asym_err = 2 * sqrt(asym_err);
    }

    h1_asym->SetBinContent(i_x2, asym_val);
    h1_asym->SetBinError  (i_x2, asym_err);
    ofs << i_x2 << "\t" << x2 << "\t" << x2_width/2 << "\t" << cnt_tot << "\t" << asym_val << "\t" << asym_err << endl;
  }
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(true);
  h1_asym->GetYaxis()->SetRangeUser(-0.1, 0.1);
  h1_asym->Draw("E1");
  h1_asym->Fit("pol0", "EM");
  c1->SaveAs("result_asym/h1_asym.png");

  ofs.close();
}
