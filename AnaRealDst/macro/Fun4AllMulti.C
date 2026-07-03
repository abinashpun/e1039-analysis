R__LOAD_LIBRARY(AnaRealDst)

/// Function to analyze multiple DST files listed in `list_run_spill.txt`.
int Fun4AllMulti(const char* fn_list="list_run_spill.txt")
{
  Fun4AllServer* se = Fun4AllServer::instance();
  //se->Verbosity(1);

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DST");
  se->registerInputManager(in);

  //se->registerSubsystem(new AnaEffHodo());
  //se->registerSubsystem(new AnaEffCham(AnaEffCham::D3p));
  se->registerSubsystem(new AnaHardEvent());

  ifstream ifs(fn_list);
  int run_id, spill_id;
  while (ifs >> run_id >> spill_id) {
    string fn_in = UtilMine::GetDstFilePath(run_id, spill_id);
    cout << "Run " << run_id << " Spill " << spill_id << ": " << fn_in << endl;
    in->fileopen(fn_in);
    se->setRun(run_id);
    se->run();
  }
  ifs.close();

  se->End();
  delete se;
  return 0;
}
