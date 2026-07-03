#include <top/G4_Beamline.C>
#include <top/G4_Target.C>
#include <top/G4_InsensitiveVolumes.C>
#include <top/G4_SensitiveDetectors.C>
#include <MacroCommon.h>
R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libg4testbench)
R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libg4dst)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libChamRealizationDev)
using namespace std;

/// Fun4All macro to generate events and carry out the track reconstruction.
/**
 * The chamber realization is *not* carried out in this step.
 * It takes 10 minutes to generate 10k events by default.
 */
int Fun4AllStep1(const int n_evt=10)
{
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", RUNNUMBER);
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  ///
  /// Event generator
  ///
  PHG4SimpleEventGenerator *genp = new PHG4SimpleEventGenerator("MUP");
  //genp->Verbosity(2);
  //genp->set_seed(123);
  const int n_mu = 1;
  genp->add_particles("mu+", n_mu);  // mu+,e+,proton,pi+,Upsilon
  genp->set_vertex_distribution_mean (0, 0, 300);
  genp->set_vertex_distribution_width(0, 0,   0);
  genp->set_vertex_size_parameters(0, 0);
  genp->set_pxpypz_range(-1,-1, +1,+1, +40,+40); // for FMAGSTR < 0
  se->registerSubsystem(genp);

  ///
  /// Detector setting
  ///
  PHG4Reco *g4Reco = new PHG4Reco();
  g4Reco->set_field_map();
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

  // save truth info to the Node Tree
  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  SQDigitizer *digitizer = new SQDigitizer("DPDigitizer", 0);
  se->registerSubsystem(digitizer);

  /// Save only events that are in the geometric acceptance.
  SQGeomAcc* geom_acc = new SQGeomAcc();
  geom_acc->SetMuonMode(SQGeomAcc::SINGLE); // PAIR, PAIR_TBBT, SINGLE, SINGLE_T, etc.
  geom_acc->SetPlaneMode(SQGeomAcc::HODO_CHAM); // HODO, CHAM or HODO_CHAM
  geom_acc->SetNumOfH1EdgeElementsExcluded(4); // Exclude 4 elements at H1 edges
  se->registerSubsystem(geom_acc);

  // Make SQ nodes for truth info
  se->registerSubsystem(new TruthNodeMaker());

  ///
  /// Reconstruction and analysis
  ///
  SQReco* reco = new SQReco();
  //reco->Verbosity(2);
  reco->set_evt_reducer_opt("none");
  reco->set_enable_eval_dst(true); // To save Tracklets
  se->registerSubsystem(reco);

  AnaTrack* ana = new AnaTrack("AnaTrack_raw");
  //ana->Verbosity(2);
  se->registerSubsystem(ana);

  ///
  /// Input, output and execution
  ///
  Fun4AllInputManager* man_in = new Fun4AllDummyInputManager("DUMMY");
  se->registerInputManager(man_in);

  Fun4AllDstOutputManager *man_out = new Fun4AllDstOutputManager("DSTOUT", "DST.root");
  se->registerOutputManager(man_out);
  man_out->AddNode("SQEvent");
  man_out->AddNode("SQHitVector");
  man_out->AddNode("SRecEvent");
  man_out->AddNode("SQMCEvent");
  man_out->AddNode("SQTruthTrackVector");
  man_out->AddNode("TrackletVector");

  const bool count_only_good_events = true;
  se->run(n_evt, count_only_good_events);

  PHGeomUtility::ExportGeomtry(se->topNode(),"geom.root");

  se->End();
  se->PrintTimer();
  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
  return 0;
}
