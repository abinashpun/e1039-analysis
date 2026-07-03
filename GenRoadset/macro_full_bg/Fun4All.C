#include <top/G4_Beamline.C>
#include <top/G4_Target.C>
#include <top/G4_InsensitiveVolumes.C>
#include <top/G4_SensitiveDetectors.C>
R__LOAD_LIBRARY(libg4dst)
R__LOAD_LIBRARY(libSQPrimaryGen)
R__LOAD_LIBRARY(libGenRoadset)
using namespace std;

int Fun4All(const string fname, const int n_evt=0, const int KMag_polarity=+1, const double KMag_scale=1.0)
{
  recoConsts *rc = recoConsts::instance();
  Fun4AllServer *se = Fun4AllServer::instance();

  ///
  /// Global parameters
  ///
  const int run_id = 5433; // To select the plane geometry.
  double FMAGSTR = -1.044; // -1.054; F&KMAGSTR were changed on 2024-06-28
  double KMAGSTR = -1.025 * KMag_scale * KMag_polarity; // -0.951
  rc->set_IntFlag   ("RUNNUMBER", run_id);
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);

  ///
  /// Event generator
  ///
  const int gen_switch = 2;
  SQPileupGen* sq_gen = new SQPileupGen();
  switch (gen_switch) {
  case 1: { // Function
    TF1* intensity_profile = new TF1("intensity_profile", "[0]*pow(x,[1])*exp(-[2]*x+exp(-[3]*x))+[4]", 0, 5000);
    intensity_profile->SetParameters(6.35, 1.38, 4.9e-3, 4.7e-3, 178.8);
    sq_gen->set_beam_intensity_profile(intensity_profile);
    break; }
  case 2: { // Histogram
/*    TFile* file = new TFile("e906_rf00.root");
    TH1D* hprof = (TH1D*)file->Get("h1_rf00");
    for (int ibin = 1; ibin <= hprof->GetNbinsX(); ++ibin) {
      hprof->SetBinContent(ibin, hprof->GetBinContent(ibin) * hprof->GetBinWidth(ibin));
    }

    TFile* file = new TFile("run6nim3.root");
    TH1D* hprof = (TH1D*)file->Get("h1_inte");
    sq_gen->set_inhibit_threshold(200000);//3000-5000
*/

    TFile* file = new TFile("nim3pot_run6.root");
    TH1F* hprof = (TH1F*)file->Get("h_pot_p00");
    sq_gen->set_inhibit_threshold(1300e3);//3000-5000

    //for (int ibin = 1; ibin <= hprof->GetNbinsX(); ++ibin) {
      //hprof->SetBinContent(ibin, hprof->GetBinContent(ibin) * hprof->GetBinWidth(ibin));
    //}

    sq_gen->set_proton_coeff(1.0);//For Run 15361-15370, thredhold = 3000 and coeff = 0.033
    sq_gen->set_beam_intensity_profile_histo(hprof);
    break; }
  case 3: {// Constant
    sq_gen->setBucketSize(5000);
    break; }
  }
  sq_gen->setExtInputFile(fname);
  se->registerSubsystem(sq_gen);

  ///
  /// Detector setting
  ///
  PHG4Reco *g4Reco = new PHG4Reco();
  g4Reco->set_field_map(
      rc->get_CharFlag("fMagFile")+" "+
      rc->get_CharFlag("kMagFile")+" "+
      Form("%f",FMAGSTR) + " " +
      Form("%f",KMAGSTR) + " " +
      "5.0",
      PHFieldConfig::RegionalConst);
  g4Reco->SetWorldSizeX(1000);
  g4Reco->SetWorldSizeY(1000);
  g4Reco->SetWorldSizeZ(5000);
  g4Reco->SetWorldShape("G4BOX");
  g4Reco->SetWorldMaterial("G4_AIR"); //G4_Galactic, G4_AIR
  g4Reco->SetPhysicsList("FTFP_BERT");

  SetupInsensitiveVolumes(g4Reco);
  SetupBeamline(g4Reco);
  SetupTarget(g4Reco);
  SetupSensitiveDetectors(g4Reco);

  se->registerSubsystem(g4Reco);

  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  SQDigitizer *digitizer = new SQDigitizer("DPDigitizer", 0);
  se->registerSubsystem(digitizer);

  TruthNodeMaker* tnm = new TruthNodeMaker();
  se->registerSubsystem(tnm);

  SubsysRecoBG* sr_bg = new SubsysRecoBG();
  sr_bg->FullBgMode();
  se->registerSubsystem(sr_bg);

  ///
  /// Input, output and execution
  ///
  Fun4AllInputManager* man_in = new Fun4AllDummyInputManager("DUMMY");
  se->registerInputManager(man_in);

  const bool count_only_good_events = true;
  se->run(n_evt, count_only_good_events);
  
  se->End();
  se->PrintTimer();
  rc->WriteToFile("recoConsts.tsv");
  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
  return 0;
}
