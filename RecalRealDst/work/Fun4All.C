R__LOAD_LIBRARY(geom_svc)
R__LOAD_LIBRARY(calibrator)
R__LOAD_LIBRARY(RecalRealDst)

int Fun4All(const int run_id, const int spill_id, const std::string fn_in, const std::string fn_out)
{
  recoConsts* rc = recoConsts::instance();
  rc->set_CharFlag("AlignmentMille", "config/align_mille_v10_a.txt");
  rc->set_CharFlag("AlignmentHodo", "");
  rc->set_CharFlag("AlignmentProp", "");
  rc->set_CharFlag("Calibration", "");

  Fun4AllServer* se = Fun4AllServer::instance();
  se->Verbosity(0);
  se->setRun(run_id);

  CalibHitElementPos* cal_ele_pos = new CalibHitElementPos();
  cal_ele_pos->CalibTriggerHit(false);
  se->registerSubsystem(cal_ele_pos);

  CalibDriftDist* cal_dd = new CalibDriftDist();
  se->registerSubsystem(cal_dd);

  /// In case you use your own SubsysReco module under `src/`.
  //CalibDriftDistMod* cal_dd = new CalibDriftDistMod();
  //se->registerSubsystem(cal_dd);

  Fun4AllInputManager* im = new Fun4AllDstInputManager("DSTIN");
  im->Verbosity(0);
  im->fileopen(fn_in.c_str());
  se->registerInputManager(im);

  Fun4AllDstOutputManager* om_dst = new Fun4AllDstOutputManager("DSTOUT", fn_out);
  se->registerOutputManager(om_dst);

  /// In case you output the SRawEvent file.
  //Fun4AllSRawEventOutputManager* om_sraw = new Fun4AllSRawEventOutputManager("SRAWOUT");
  //om_sraw->SetFileName(fn_out);
  //se->registerOutputManager(om_sraw);

  se->run();
  se->End();
  //rc->Print();
  delete se;
  //gSystem->Exit(0);
  return 0;
}
