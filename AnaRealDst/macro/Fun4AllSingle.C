R__LOAD_LIBRARY(AnaRealDst)

/// Main function.
int Fun4AllSingle(const int run_id=5518, const int spill_id=1898254, const int nevent=0)
{
  Fun4AllServer* se = Fun4AllServer::instance();
  se->setRun(run_id);
  //se->Verbosity(1);

  //se->registerSubsystem(new AnaHodoHit());
  //se->registerSubsystem(new AnaTriggerHit());
  //se->registerSubsystem(new AnaRealDst());
  //se->registerSubsystem(new AnaEffHodo());
  //se->registerSubsystem(new AnaEffCham(AnaEffCham::D3p));
  se->registerSubsystem(new AnaHardEvent());

  string fn_in = UtilMine::GetDstFilePath(run_id, spill_id);
  cout << "DST file = " << fn_in << endl;
  Fun4AllInputManager *in = new Fun4AllDstInputManager("DST");
  in->fileopen(fn_in);
  se->registerInputManager(in);

  se->run(nevent);
  se->End();
  delete se;
  return 0;
}
