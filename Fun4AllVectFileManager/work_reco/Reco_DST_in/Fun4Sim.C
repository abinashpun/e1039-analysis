R__LOAD_LIBRARY(libcalibrator)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libfun4all_vect_file_manager)

int Fun4Sim(const int n_evt = 10){
  const double FMAGSTR = -1.044;
  const double KMAGSTR = -1.025;
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", 5433); /// The geometry is selected based on run number.
  rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
  rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
  rc->set_DoubleFlag("SIGX_BEAM", 0.3);
  rc->set_DoubleFlag("SIGY_BEAM", 0.3);
  rc->set_DoubleFlag("Z_UPSTREAM", -700.);
  

  Fun4AllServer* se = Fun4AllServer::instance();
  //se->setRun(run_id);
  se->Verbosity(100);

  SQReco* reco = new SQReco();
  reco->Verbosity(1);
  reco->set_legacy_rec_container(true); // default = true
  reco->set_geom_file_name((string)gSystem->Getenv("E1039_RESOURCE") + "/geometry/geom_run005433.root");
  reco->set_enable_KF(true);
  reco->setInputTy(SQReco::E1039);
  reco->setFitterTy(SQReco::KFREF);
  reco->set_evt_reducer_opt("none");
  reco->set_enable_eval_dst(true);
  for (int ii = 0; ii <= 3; ii++) reco->add_eval_list(ii);
  reco->set_enable_eval(true);
  reco->set_eval_file_name("eval.root");
  se->registerSubsystem(reco);

  SQVertexing* vtx = new SQVertexing();
  vtx->Verbosity(100);
  vtx->set_legacy_rec_container(true); // default = false
  //vtx->set_single_retracking(true);
  se->registerSubsystem(vtx);

 Fun4AllInputManager *in = new Fun4AllDstInputManager("DST");
 se->registerInputManager(in);
 in->fileopen("/seaquest/users/mhossain/HitMatrix/e1039-analysis/SimChainDev/data/JPsi_vec_test/1/out/DST.root"); //This DST file contains the hit information from simulation or real data.
  // DST output manager
  Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", "DST.root");
  se->registerOutputManager(out);
  se->run();
  se->End();
  //rc->WriteToFile("recoConsts.tsv");
  delete se;
  exit(0); //return 0;
}
