//R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(RecoData2024)

/// Fun4All macro to analyze spill-by-spill DST files of multiple runs.
int Fun4All(const int run_id, const string dir_in, const string fn_list, const int n_evt=0)
{
  recoConsts* rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", run_id);
  rc->set_DoubleFlag("FMAGSTR", -1.044);
  rc->set_DoubleFlag("KMAGSTR", -1.025);
  //rc->set_CharFlag("AlignmentMille", "config/align_mille_v10.txt");
  //rc->set_CharFlag("AlignmentHodo", "");
  //rc->set_CharFlag("AlignmentProp", "");
  //rc->set_CharFlag("Calibration", "");
  
  Fun4AllServer* se = Fun4AllServer::instance();
  se->setRun(run_id);
  //se->Verbosity(1);

  //SQVertexing* vtx = new SQVertexing();
  //vtx->Verbosity(99);
  //se->registerSubsystem(vtx);
  
  //se->registerSubsystem(new AnaDimuon());
  se->registerSubsystem(new AnaTrack());

  vector<string> list_in;
  ifstream ifs(fn_list);
  int spill;
  string bname;
  while (ifs >> spill >> bname) list_in.push_back(dir_in+"/"+bname+".root");
  ifs.close();

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DST");
  se->registerInputManager(in);
  for (auto it = list_in.begin(); it != list_in.end(); it++) {
    cout << "Input = " << *it << endl;
    in->fileopen(*it);
    se->run(n_evt);
  }
  se->End();
  delete se;
  exit(0); //return 0;
}
