/// Fun4SimTree.C:  Macro to analyze the simulated tree created by Fun4SimMicroDst.C.
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
R__LOAD_LIBRARY(libana_sim_dst)
#endif

using namespace std;
TCanvas* c1;
void DrawDimTrueKin(TTree* tr);
void DrawDimRecoKin(TTree* tr);
void DrawTrkTrueKin(TTree* tr);
void DrawTrueVar(TTree* tr, const string varname, const string title_x, const int n_x, const double x_lo, const double x_hi);
void AnaEvents(TTree* tr);

void Fun4SimTree(const char* fname="sim_tree.root", const char* tname="tree")
{
  //gSystem->Load("libana_sim_dst.so");

  TFile* file = new TFile(fname);
  TTree* tr = (TTree*)file->Get(tname);

  gSystem->mkdir("result", true);
  c1 = new TCanvas("c1", "");
  c1->SetGrid();
  c1->SetLogy(true);
  
  DrawDimTrueKin(tr);
  DrawDimRecoKin(tr);
  DrawTrkTrueKin(tr);

  AnaEvents(tr);

  exit(0);
}

///
/// Functions below
///
void DrawDimTrueKin(TTree* tr)
{
  tr->Draw("n_dim_true");
  c1->SaveAs("result/h1_true_n_dim.png");
  tr->Draw("n_dim_reco");
  c1->SaveAs("result/h1_reco_n_dim.png");

  const double PI = TMath::Pi();
  DrawTrueVar(tr, "dim_true.pdg_id"    , "True dimuon PDG ID", 1000, 0, 0);
  DrawTrueVar(tr, "dim_true.mom.X()"   , "True dimuon px (GeV)", 100, -5,   5);
  DrawTrueVar(tr, "dim_true.mom.Y()"   , "True dimuon py (GeV)", 100, -5,   5);
  DrawTrueVar(tr, "dim_true.mom.Z()"   , "True dimuon pz (GeV)", 100,  0, 100);
  DrawTrueVar(tr, "dim_true.mom.M()"   , "True dimuon mass (GeV)", 100, 0, 5);
  DrawTrueVar(tr, "dim_true.mom.Eta()" , "True dimuon #eta", 110, 0, 11);
  DrawTrueVar(tr, "dim_true.mom.Phi()" , "True dimuon #phi", 100, -PI, PI);
  DrawTrueVar(tr, "dim_true.x1"        , "True x1", 50, 0, 1);
  DrawTrueVar(tr, "dim_true.x2"        , "True x2", 50, 0, 1);
}

void DrawDimRecoKin(TTree* tr)
{
  tr->Draw("rec_stat"); // cf. GlobalConsts.h.
  c1->SaveAs("result/h1_rec_stat.png");
  
  tr->Draw("trig_bits", "rec_stat==0");
  c1->SaveAs("result/h1_trig_bits.png");

  tr->Draw("dim_reco.mom.M()", "rec_stat==0");
  c1->SaveAs("result/h1_dim_reco_mass.png");

  tr->Draw("dim_reco.x1", "rec_stat==0");
  c1->SaveAs("result/h1_dim_reco_x1.png");

  tr->Draw("dim_reco.x2", "rec_stat==0");
  c1->SaveAs("result/h1_dim_reco_x2.png");
}


void DrawTrkTrueKin(TTree* tr)
{
  DrawTrueVar(tr, "dim_true.mom_pos.X()", "True px (GeV) of mu+", 100, -5, 5);
  DrawTrueVar(tr, "dim_true.mom_pos.Y()", "True py (GeV) of mu+", 100, -5, 5);
  DrawTrueVar(tr, "dim_true.mom_pos.Z()", "True pz (GeV) of mu+", 100,  0, 100);
  DrawTrueVar(tr, "dim_true.mom_neg.X()", "True px (GeV) of mu-", 100, -5, 5);
  DrawTrueVar(tr, "dim_true.mom_neg.Y()", "True py (GeV) of mu-", 100, -5, 5);
  DrawTrueVar(tr, "dim_true.mom_neg.Z()", "True pz (GeV) of mu-", 100,  0, 100);

  THStack* hs;
  TH1* h1_all = new TH1D("h1_all", "", 100, -1, 1);
  TH1* h1_rec = new TH1D("h1_rec", "", 100, -1, 1);
  tr->Project("h1_all", "(dim_true.mom_pos.Z() - dim_true.mom_neg.Z())/(dim_true.mom_pos.Z() + dim_true.mom_neg.Z())");
  tr->Project("h1_rec", "(dim_true.mom_pos.Z() - dim_true.mom_neg.Z())/(dim_true.mom_pos.Z() + dim_true.mom_neg.Z())", "rec_stat==0");
  hs = new THStack("hs", "J/#psi GMC;gpz+gpz (GeV) of tracks;N of tracks");
  hs->Add(h1_all);
  hs->Add(h1_rec);
  h1_rec->SetLineColor(kRed);
  hs->Draw("nostack");
  c1->SaveAs("result/h1_trk_true_pz_asym.png");
}

void DrawTrueVar(TTree* tr, const string varname, const string title_x, const int n_x, const double x_lo, const double x_hi)
{
  TH1* h1_all = new TH1D("h1_all", "", n_x, x_lo, x_hi);
  TH1* h1_rec = new TH1D("h1_rec", "", n_x, x_lo, x_hi);
  tr->Project("h1_all", varname.c_str());
  tr->Project("h1_rec", varname.c_str(), "rec_stat==0");

  ostringstream oss;
  oss << "J/#psi GMC;" << title_x << ";Yield";
  THStack hs("hs", oss.str().c_str());
  hs.Add(h1_all);
  hs.Add(h1_rec);
  h1_rec->SetLineColor(kRed);
  hs.Draw("nostack");

  oss.str("");
  oss << "result/h1_";
  for (string::const_iterator it = varname.begin(); it != varname.end(); it++) {
    switch (*it) { // modify bad chars for file name
    case '.': case '*': case '/': oss << '_'; break;
    case '(': case ')': case ' ': /* omit */ break;
    default: oss << *it;
    }
  }
  oss << ".png";
  c1->SaveAs(oss.str().c_str());

  delete h1_all;
  delete h1_rec;
}

void AnaEvents(TTree* tr)
{
  typedef map<int, int> IntCount_t;
  IntCount_t id_cnt;
  DimuonList* list_dim = new DimuonList();
  tr->SetBranchAddress("dim_true", &list_dim);

  int n_ent = tr->GetEntries();
  cout << "AnaEvents(): n = " << n_ent << endl;
  for (int i_ent = 0; i_ent < n_ent; i_ent++) {
    if ((i_ent+1) % (n_ent/10) == 0) cout << "  " << 100*(i_ent+1)/n_ent << "%" << flush;
    tr->GetEntry(i_ent);
    for (DimuonList::iterator it = list_dim->begin(); it != list_dim->end(); it++) {
      DimuonData* dd = &(*it);
      int pdg_id = dd->pdg_id;
      if (id_cnt.find(pdg_id) == id_cnt.end()) id_cnt[pdg_id] = 1;
      else                                     id_cnt[pdg_id]++;
    }
  }
  cout << endl;
  for (IntCount_t::iterator it = id_cnt.begin(); it != id_cnt.end(); it++) {
    cout << setw(10) << it->first << "  " << setw(10) << it->second << endl;
  }
}
