/// Fun4All.C:  Fun4all macro to execute "ExtractTdcDist".
R__LOAD_LIBRARY(libCalibChamXT)

int Fun4All(const int run, const char* fn_dst, const int n_evt=0)
{
  recoConsts* rc = recoConsts::instance();
  rc->init("cosmic");
  rc->set_IntFlag("RUNNUMBER", run);

  Fun4AllServer* se = Fun4AllServer::instance();
  //se->Verbosity(1);

  se->registerSubsystem(new ExtractTdcDist("d0"));
  se->registerSubsystem(new ExtractTdcDist("d2"));
  se->registerSubsystem(new ExtractTdcDist("d3p"));
  se->registerSubsystem(new ExtractTdcDist("d3m"));

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTIN");
  in->Verbosity(0);
  in->fileopen(fn_dst);
  se->registerInputManager(in);

  se->run(n_evt);
  se->End();
  //rc->WriteToFile("recoConsts.tsv");
  delete se;
  exit(0);
}
