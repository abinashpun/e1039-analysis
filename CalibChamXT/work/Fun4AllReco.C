/// Fun4AllReco.C:  Fun4all macro to run the reconstruction.
// /seaquest/users/apun/abi_project/data_manage/e1039-data-mgt_test/RecoE1039Data.C
#include <top/G4_InsensitiveVolumes.C>
#include <top/G4_SensitiveDetectors.C>
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libphg4hit) // Need be included in libevt_filter
R__LOAD_LIBRARY(libevt_filter)
R__LOAD_LIBRARY(libcalibrator)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libCalibChamXT)

int Fun4AllReco(const int run, const int iter, const char* fn_dst, const int n_evt)
{
  recoConsts* rc = recoConsts::instance();
  rc->init("cosmic");
  rc->set_IntFlag("RUNNUMBER", run);
  rc->set_BoolFlag("COARSE_MODE", false);
  rc->set_DoubleFlag("KMAGSTR", 0.);
  rc->set_DoubleFlag("FMAGSTR", 0.);
  rc->set_IntFlag("MaxHitsDC0" , 12); // default = 100
  rc->set_IntFlag("MaxHitsDC1" , 12); // default = 100
  rc->set_IntFlag("MaxHitsDC2" , 12); // default = 100
  rc->set_IntFlag("MaxHitsDC3p", 12); // default = 100
  rc->set_IntFlag("MaxHitsDC3m", 12); // default = 100
  //rc->Print();

  Fun4AllServer* se = Fun4AllServer::instance();
  //se->Verbosity(1);

  PHG4Reco* g4reco = new PHG4Reco();
  g4reco->set_field_map();
  g4reco->SetWorldSizeX(1000);
  g4reco->SetWorldSizeY(1000);
  g4reco->SetWorldSizeZ(5000);
  g4reco->SetWorldShape("G4BOX");
  g4reco->SetWorldMaterial("G4_AIR"); //G4_Galactic, G4_AIR
  g4reco->SetPhysicsList("FTFP_BERT");
  SetupInsensitiveVolumes(g4reco);
  SetupSensitiveDetectors(g4reco);
  se->registerSubsystem(g4reco);

  FilterByTrigger* fbt = new FilterByTrigger();
  fbt->SetNimBits(1,1,0,1,0);
  fbt->EnableOutput();
  se->registerSubsystem(fbt);

  CalibDriftDist* cal_dd = new CalibDriftDist();
  cal_dd->Verbosity(1);
  if (iter == 1) {
    cal_dd->SetResolution(0.03, 0.03, 0.03, 0.03, 0.03); // (D0, D1, D2, D3p, D3m) in cm.
  //} else {
  //cal_dd->ReadParamFromFile(fn_in_time, fn_xt_curve);
  }
  se->registerSubsystem(cal_dd);

  if (iter == 1) {
    se->registerSubsystem(new AnaChamPlane("D2V" ));
    se->registerSubsystem(new AnaChamPlane("D2Vp"));
    se->registerSubsystem(new AnaChamPlane("D2Xp"));
    se->registerSubsystem(new AnaChamPlane("D2X" ));
    se->registerSubsystem(new AnaChamPlane("D2U" ));
    se->registerSubsystem(new AnaChamPlane("D2Up"));

    se->registerSubsystem(new AnaChamPlane("D3pV" ));
    se->registerSubsystem(new AnaChamPlane("D3pVp"));
    se->registerSubsystem(new AnaChamPlane("D3pX" ));
    se->registerSubsystem(new AnaChamPlane("D3pXp"));
    se->registerSubsystem(new AnaChamPlane("D3pU" ));
    se->registerSubsystem(new AnaChamPlane("D3pUp"));

    se->registerSubsystem(new AnaChamPlane("D3mV" ));
    se->registerSubsystem(new AnaChamPlane("D3mVp"));
    se->registerSubsystem(new AnaChamPlane("D3mX" ));
    se->registerSubsystem(new AnaChamPlane("D3mXp"));
    se->registerSubsystem(new AnaChamPlane("D3mU" ));
    se->registerSubsystem(new AnaChamPlane("D3mUp"));
  }

  auto reco = new SQTrackletReco();
  //reco->Verbosity(9);
  reco->set_legacy_rec_container(false);
  reco->set_enable_KF(false); // default = true
  reco->set_evt_reducer_opt("none"); // was "e", default = aoc, "" = aoc
  reco->set_enable_eval(true);
  reco->set_enable_eval_dst(true);
  reco->add_eval_list(0); // D0, D1
  reco->add_eval_list(1); // D2
  reco->add_eval_list(2); // D3p, D3m
  reco->add_eval_list(3); // D2+3, back-partial tracklets
  reco->add_eval_list(4); // D1+2+3
  reco->drop_empty_event(true);
  se->registerSubsystem(reco);

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTIN");
  in->Verbosity(0);
  in->fileopen(fn_dst);
  se->registerInputManager(in);

  Fun4AllDstOutputManager* out = new Fun4AllDstOutputManager("DSTOUT", "DSTreco.root");
  se->registerOutputManager(out);
  out->AddNode("SQEvent");
  out->AddNode("TrackletVector");

  se->run(n_evt);
  se->End();
  //se->PrintTimer();
  rc->WriteToFile("recoConsts.tsv");
  delete se;
  return 0;
}
