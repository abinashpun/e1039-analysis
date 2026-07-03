/// Fun4AllReco.C:  Fun4all macro to run the reconstruction.
// /seaquest/users/apun/abi_project/data_manage/e1039-data-mgt_test/RecoE1039Data.C
R__LOAD_LIBRARY(libg4detectors)
R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libCalibChamXT)

int Fun4AllAna(const int run, const int iter, const char* fn_list)
{
  recoConsts* rc = recoConsts::instance();
  rc->init("cosmic");
  rc->set_IntFlag("RUNNUMBER", run);
  rc->set_BoolFlag("COARSE_MODE", false);
  rc->set_DoubleFlag("KMAGSTR", 0.);
  rc->set_DoubleFlag("FMAGSTR", 0.);

  Fun4AllServer* se = Fun4AllServer::instance();

  SRMakeRTCurve* ana = new SRMakeRTCurve(iter);
  se->registerSubsystem(ana);

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTIN");
  in->Verbosity(0);
  se->registerInputManager(in);

  ifstream ifs(fn_list);
  string fname;
  vector<string> list_fname;
  while (ifs >> fname) list_fname.push_back(fname);
  ifs.close();

  unsigned int n_file = list_fname.size();
  for (unsigned int i_file = 0; i_file < n_file; i_file++) {
    fname = list_fname[i_file];
    cout << "Input " << i_file << " / " << n_file << ": " << fname << endl;
    in->fileopen(fname.c_str());
    se->run();
  }

  gErrorIgnoreLevel = 1111;
  se->End();
  delete se;
  return 0;
}
