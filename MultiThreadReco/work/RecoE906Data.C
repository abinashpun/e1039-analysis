#include <TSystem.h>
#include <top/G4_InsensitiveVolumes.C>
#include <top/G4_SensitiveDetectors.C>
//R__LOAD_LIBRARY(libinterface_main)
//R__LOAD_LIBRARY(libfun4all)
//R__LOAD_LIBRARY(libg4detectors)
//R__LOAD_LIBRARY(libg4eval)
R__LOAD_LIBRARY(libcalibrator)
//R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(MultiThreadReco)

int RecoE906Data(const int run_id, const int spill_id, const string fname_in, const string fname_out, const int n_evt=100)
{
  ROOT::EnableThreadSafety();

  recoConsts* rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", 6); // To select the plane geometry for E906 Run 6.
  rc->set_DoubleFlag("FMAGSTR", -1.044);
  rc->set_DoubleFlag("KMAGSTR", -1.025);
  rc->set_CharFlag("TRIGGER_L1", "78");
  rc->set_CharFlag("EventReduceOpts", "aoc");
  rc->set_DoubleFlag("RejectWinDC0" , 0.12);
  rc->set_DoubleFlag("RejectWinDC1" , 0.25);
  rc->set_DoubleFlag("RejectWinDC2" , 0.15);
  rc->set_DoubleFlag("RejectWinDC3p", 0.16);
  rc->set_DoubleFlag("RejectWinDC3m", 0.14);
  rc->set_IntFlag("MaxHitsDC0" , 70); // 350);
  rc->set_IntFlag("MaxHitsDC1" , 70); // 350);
  rc->set_IntFlag("MaxHitsDC2" , 34); // 170);
  rc->set_IntFlag("MaxHitsDC3p", 28); // 140);
  rc->set_IntFlag("MaxHitsDC3m", 28); // 140);
  //rc->Print();

  Fun4AllServer* se = Fun4AllServer::instance();
  //se->Verbosity(100);

  ///
  /// Geometry.  You might create and use "geom.root" instead.
  ///
  PHG4Reco* g4reco = new PHG4Reco();
  g4reco->set_field_map();
  g4reco->SetWorldSizeX(1000);
  g4reco->SetWorldSizeY(1000);
  g4reco->SetWorldSizeZ(5000);
  g4reco->SetWorldShape("G4BOX");
  g4reco->SetWorldMaterial("G4_AIR");
  g4reco->SetPhysicsList("FTFP_BERT");
  SetupInsensitiveVolumes(g4reco);
  SetupSensitiveDetectors(g4reco);
  se->registerSubsystem(g4reco);

  ///
  /// Calibrator
  ///
  CalibHitElementPos* cal_ele_pos = new CalibHitElementPos();
  se->registerSubsystem(cal_ele_pos);

  CalibDriftDist* cal_drift_dist = new CalibDriftDist();
  cal_drift_dist->Verbosity(1);
  se->registerSubsystem(cal_drift_dist);

  ///
  /// Reconstruction
  ///
  auto reco = new MultiThreadReco();
  //reco->Verbosity(100);
  reco->SetOutputFileName(fname_out);
  //KScheduler::SetNumThreads(16); // default = 16
  //KScheduler::SetInputPipeDepth(32); // default = 32
  //KScheduler::SetOutputPipeDepth(32); // default = 32
  KScheduler::Verbose(1); // default = 0
  KScheduler::PrintFreq(1000); // default = 100
  //KScheduler::SaveNum(5000); // default = 5000
  //reco->setInputTy(SQReco::E1039); // default = E1039
  reco->set_enable_KF(true); // default = true
  //reco->setFitterTy(SQReco::KFREF); // default = KFREF
  se->registerSubsystem(reco);

  //VertexFit* vtx_fit = new VertexFit();
  //se->registerSubsystem(vtx_fit);

  ///
  /// Input, output and run.
  ///
  Fun4AllSRawEventInputManager* in = new Fun4AllSRawEventInputManager("SRawEventIM");
  in->Verbosity(0);
  in->enable_E1039_translation();
  in->set_tree_name("save");
  in->set_branch_name("rawEvent");
  in->fileopen(fname_in);
  se->registerInputManager(in);

  se->run(n_evt);
  //PHGeomUtility::ExportGeomtry(se->topNode(), "geom.root");
  //rc->WriteToFile("recoConsts.tsv");
  se->End();
  delete se;
  gSystem->Exit(0);
  return 0;
}
