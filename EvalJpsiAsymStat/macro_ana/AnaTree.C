/// macro_ana/AnaTree.C:  Macro to analyze the tree created by Fun4All.C
R__LOAD_LIBRARY(EvalJpsiAsymStat)

using namespace std;
TFile* file_out;
string mode; // jpsi, psip, dy
string generator; // legacy, pythia
TFile* file;
TTree* tree;
TCanvas* c1;
double inte_lumi;
void MakeHist1D(const string name, const string var_x, const string title_x, const int n_x, const double x_lo, const double x_hi);
double GetInteLumi(const char* fn_lumi="lumi_tot.txt");

///
/// Main
///
void AnaTree(const char* mode_in="jpsi", const char* generator_in="legacy", const char* fname="sim_tree.root", const char* tname="tree")
{
  mode = mode_in;
  generator = generator_in;
  file = new TFile(fname);
  tree = (TTree*)file->Get(tname);
  gSystem->mkdir("result", true);
  c1 = new TCanvas("c1", "");
  c1->SetGrid();
  gStyle->SetOptStat(0);

  file_out = new TFile("result/hist.root", "RECREATE");

  inte_lumi = GetInteLumi();
  cout << "Integrated luminosity = " << inte_lumi << endl;

  c1->SetLogy(true);
  tree->Draw("n_dim_true");
  c1->SaveAs("result/h1_true_n_dim.png");
  tree->Draw("rec_stat", "n_dim_true==1");
  c1->SaveAs("result/h1_rec_stat.png");
  tree->Draw("weight", "n_dim_true==1");
  c1->SaveAs("result/h1_weight.png");

  c1->SetLogy(false);
  const double PI = TMath::Pi();
  //MakeHist1D("vtx_x", "dim_true.pos.X()"   , "X-vertex", 100, -2, 2);
  //MakeHist1D("vtx_y", "dim_true.pos.Y()"   , "Y-vertex", 100, -2, 2);
  MakeHist1D("vtx_z", "dim_true.pos.Z()"   , "Z-vertex", 100, -310, -290);
  //MakeHist1D("px"      , "dim_true.mom.X()"   , "px"  , 50, -5, 5);
  //MakeHist1D("py"      , "dim_true.mom.Y()"   , "py"  , 50, -5, 5);
  MakeHist1D("pz"      , "dim_true.mom.Z()"   , "pz"  , 30, 30, 120);
  MakeHist1D("pT"      , "dim_true.mom.Perp()", "pT"  , 10, 0,   5);
  MakeHist1D("mass"    , "dim_true.mom.M()"   , "Mass"        , 50, 1, 6);
  MakeHist1D("mass_rec", "dim_reco.mom.M()"   , "Reco. Mass"  , 50, 1, 6);
  //MakeHist1D("phi"     , "dim_true.mom.Phi()" , "#phi", 30, -PI, PI);
  MakeHist1D("x1"      , "dim_true.x1"        , "x1"  , 10, 0.0, 1.0);
  MakeHist1D("x2"      , "dim_true.x2"        , "x2"  , 10, 0.0, 0.2);
  MakeHist1D("xF"      , "dim_true.xF"        , "xF"  , 10, 0.0, 1.0);
  //MakeHist1D("cs_costh", "dim_true.costh"     , "cos#theta in C-S frame", 50, -1, 1);
  //MakeHist1D("cs_phi"  , "dim_true.phi"       , "#phi in C-S frame"     , 50, 0, 2*PI);

  file_out->Write();
  file_out->Close();
  exit(0);
}

///
/// Sub-Functions
///
void MakeHist1D(const string name, const string var_x, const string title_x, const int n_x, const double x_lo, const double x_hi)
{
  static double mass_min = -1;
  static double mass_max = -1;
  static double ww_br    = -1; // weight for branching ratio
  if (mass_min < 0) {
    if (mode == "jpsi") {
      mass_min = 3.097 - 0.010;
      mass_max = 3.097 + 0.010;
      ww_br    = (generator == "pythia" ? 0.0593 : 1.0);
    } else if (mode == "psip") {
      mass_min = 3.686 - 0.010;
      mass_max = 3.686 + 0.010;
      ww_br    = (generator == "pythia" ? 7.3e-3 : 1.0);
    } else if (mode == "dy") {
      mass_min =  0.0;
      mass_max = 15.0;
      ww_br    =  1.0;
    } else {
      cout << "Mode '" << mode << "' is not supported.  Abort." << endl;
      exit(1);
    }
    cout << "mass_min = " << mass_min << "\n"
         << "mass_max = " << mass_max << "\n"
         << "ww_br    = " << ww_br << "\n";
  }

  ostringstream oss;
  oss << "h1_" << name;
  string name_h1 = oss.str();
  TH1* h1 = new TH1D(name_h1.c_str(), "", n_x, x_lo, x_hi);

  oss.str("");
  oss << ";" << title_x << ";Rate (pb)";
  h1->SetTitle(oss.str().c_str());

  oss.str("");
  oss << "weight * " << ww_br << " * (n_dim_true == 1 && " << mass_min << " < dim_true.mom.M() && dim_true.mom.M() < " << mass_max << ")";
  string par_cut = oss.str();

  tree->Project(name_h1.c_str(), var_x.c_str(), par_cut.c_str());
  h1->Scale(1/inte_lumi);
  h1->SetLineColor(kRed);
  h1->Draw("E1");

  oss.str("");
  oss << "result/" << name_h1 << ".png";
  c1->SaveAs(oss.str().c_str());
}

double GetInteLumi(const char* fn_lumi)
{
  ifstream ifs(fn_lumi);
  if (!ifs.is_open()) {
    cout << "GetInteLumi():  Cannot open '" << fn_lumi << "'.  Just return 1.0." << endl;
    return 1.0;
  }
  double val;
  ifs >> val;
  ifs.close();
  return val;
}
