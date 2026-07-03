#include <MacroCommon.h>
R__LOAD_LIBRARY(libfun4all)
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libChamRealizationDev)
using namespace std;

/// Fun4All macro to carry out the chamber realization and then the track reconstruction.
int Fun4AllStep2(const int n_evt=0, const char* fn_dst="DST.root")
{
  ///
  /// Global parameters
  ///
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", RUNNUMBER);
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);

  Fun4AllServer *se = Fun4AllServer::instance();

  ///
  /// Reconstruction and analysis
  ///
  SQChamberRealization* cal_cr = new SQChamberRealization();
  //cal_cr->Verbosity(2);
  cal_cr->SetChamEff(0.94, 0.94, 0.94, 0.94, 0.94); // (D0, D1, D2, D3p, D3m)
  //cal_cr->ScaleChamReso(2, 2, 2, 2, 2); // (D0, D1, D2, D3p, D3m)
  //cal_cr->FixChamReso(0.04, 0.04, 0.04, 0.04, 0.04); // (D0, D1, D2, D3p, D3m) in cm
  se->registerSubsystem(cal_cr);

  SQReco* reco = new SQReco();
  //reco->Verbosity(2);
  reco->use_geom_io_node(true); // Use geom.root.
  reco->set_evt_reducer_opt("none");
  reco->set_enable_eval_dst(true); // To save Tracklets
  se->registerSubsystem(reco);

  AnaTrack* ana = new AnaTrack("AnaTrack_real");
  //ana->Verbosity(2);
  se->registerSubsystem(ana);

  ///
  /// Input, output and execution
  ///
  Fun4AllInputManager* man_in = new Fun4AllDstInputManager("DSTIN");
  se->registerInputManager(man_in);

  man_in->fileopen(fn_dst);
  se->run(n_evt);
  
  se->End();
  se->PrintTimer();
  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
  return 0;
}
