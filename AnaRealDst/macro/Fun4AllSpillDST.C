R__LOAD_LIBRARY(AnaRealDst)

/// Fun4All macro to analyze spill-by-spill DST files of the given run.
int Fun4AllSpillDST(const int run_id=5598, const int nevent=0)
{
  Fun4AllServer* se = Fun4AllServer::instance();
  se->setRun(run_id);
  //se->Verbosity(1);

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DST");
  se->registerInputManager(in);

  //se->registerSubsystem(new AnaHodoHit());
  //se->registerSubsystem(new AnaTriggerHit());
  se->registerSubsystem(new AnaTrigSignal());
  //se->registerSubsystem(new AnaRealDst());
  //se->registerSubsystem(new AnaEffHodo());
  //se->registerSubsystem(new AnaEffCham(AnaEffCham::D3p));
  //se->registerSubsystem(new AnaHardEvent());

  auto list_dst = UtilMine::GetListOfSpillDSTs(run_id);
  for (auto it = list_dst.begin(); it != list_dst.end(); it++) {
    string fn_in = *it;
    cout << "DST file = " << fn_in << endl;
    in->fileopen(fn_in);
    se->run(nevent);
  }

  se->End();
  delete se;
  exit(0); //return 0;
}
