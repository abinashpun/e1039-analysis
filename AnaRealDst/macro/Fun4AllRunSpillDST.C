R__LOAD_LIBRARY(AnaRealDst)

/// Fun4All macro to analyze spill-by-spill DST files of multiple runs.
int Fun4AllRunSpillDST(const char* fn_list="list_run.txt", const int nevent=0)
{
  vector<int> list_run_id;
  ifstream ifs(fn_list);
  int run_id;
  while (ifs >> run_id) list_run_id.push_back(run_id);
  ifs.close();
  if (list_run_id.size() == 0) return 0;

  Fun4AllServer* se = Fun4AllServer::instance();
  se->setRun(list_run_id[0]);
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

  for (auto it1 = list_run_id.begin(); it1 != list_run_id.end(); it1++) {
    run_id = *it1;
    cout << "Run " << run_id << endl;
    auto list_dst = UtilMine::GetListOfSpillDSTs(run_id);
    for (auto it2 = list_dst.begin(); it2 != list_dst.end(); it2++) {
      string fn_in = *it2;
      cout << "  DST file = " << fn_in << endl;
      in->fileopen(fn_in);
      se->run(nevent);
    }
  }
  se->End();
  delete se;
  exit(0); //return 0;
}
