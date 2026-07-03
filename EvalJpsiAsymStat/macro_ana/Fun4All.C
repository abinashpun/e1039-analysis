/// macro_ana/Fun4All.C:  Fun4all macro to analyze the simDST files for J/psi TSSA.
R__LOAD_LIBRARY(EvalJpsiAsymStat)

int Fun4All(const int n_dst_ana=0, const char* fn_list_dst="list_dst.txt")
{
  Fun4AllServer* se = Fun4AllServer::instance();
  //se->Verbosity(1);
  Fun4AllInputManager *in = new Fun4AllDstInputManager("SimDst");
  se->registerInputManager(in);

  se->registerSubsystem(new AnaSimRunInfo());
  se->registerSubsystem(new AnaSimDimuon());

  vector<string> list_dst;
  string fn_dst;
  ifstream ifs(fn_list_dst);
  while (ifs >> fn_dst) list_dst.push_back(fn_dst);
  ifs.close();

  int n_dst = list_dst.size();
  cout << "N of DSTs: all = " << n_dst;
  if (n_dst_ana > 0 && n_dst > n_dst_ana) n_dst = n_dst_ana;
  cout << ", to be analyzed = " << n_dst << endl;
  for (int i_dst = 0; i_dst < n_dst; i_dst++) {
    string fn_dst = list_dst[i_dst];
    cout << "DST: " << i_dst << "/" << n_dst << ": " << fn_dst << endl;
    in->fileopen(fn_dst);
    se->run();
  }

  se->End();
  delete se;
  return 0;
}
