R__LOAD_LIBRARY(libAnaCustomFile)

/**
 * Fun4All macro to analyze a Track-QA v1 file.
 */
int AnaTrackQAFile(
  const int n_evt=10,
  const string fn_in="/project/ptgroup/spinquest/MC_merge_files/singMup_x2y2z300_750K/merged_trackQA.root")
{
  recoConsts* rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", 1); // To select the right plane geometry.
  rc->set_DoubleFlag("FMAGSTR", -1.054);
  rc->set_DoubleFlag("KMAGSTR", -0.951);
  rc->set_BoolFlag("COARSE_MODE", true);
  rc->Print();

  Fun4AllServer* se = Fun4AllServer::instance();
  se->Verbosity(0);

  auto calib_hit_pos = new CalibHitPos();
  //calib_hit_pos->Verbosity(1); // 0 = none, 1 = old+new pos
  se->registerSubsystem(calib_hit_pos);
  
  auto check_evt = new CheckEventData();
  check_evt->Verbosity(0); // 0 = none, 1 = det:ele, 2 = det:ele[dist]
  se->registerSubsystem(check_evt);

  SQReco* reco = new SQReco();
  reco->Verbosity(2); // 0 = none, 2 = some, 4 = lot
  reco->set_legacy_rec_container(false); // false = SQTrackVector, true = SRecEvent
  reco->set_geom_file_name("geom.root");
  reco->set_enable_KF(true);
  reco->setInputTy(SQReco::E1039);
  reco->setFitterTy(SQReco::KFREF);
  reco->set_evt_reducer_opt("none");
  reco->set_enable_eval(true);
  reco->set_eval_file_name("eval.root");
  se->registerSubsystem(reco);

  auto in = new Fun4AllTrackQAv1InputManager("IN");
  in->Verbosity(0);
  in->fileopen(fn_in);
  se->registerInputManager(in);

  se->run(n_evt);

  se->End();
  se->PrintTimer();
  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
  return 0;
}
