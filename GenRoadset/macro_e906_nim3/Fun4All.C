R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(libGenRoadset)
using namespace std;

int Fun4All(const int run, const string fname, const int nevent=0)
{
  Fun4AllServer* se = Fun4AllServer::instance();
  //se->Verbosity(1);
  //const bool use_trigger_hit = false;

  Fun4AllSRawEventInputManager* in = new Fun4AllSRawEventInputManager("SRawEventIM");
  in->Verbosity(0);
  in->enable_E1039_translation();
  se->registerInputManager(in);

  E906SpillSelector* e906ss = new E906SpillSelector();
  e906ss->EnableOutput();
  se->registerSubsystem(e906ss);

  FilterByTrigger* fbt = new FilterByTrigger();
  fbt->SetNimBits(0,0,1,0,0);
  fbt->EnableOutput();
  se->registerSubsystem(fbt);

  SubsysRecoBG* sr_bg = new SubsysRecoBG();
  se->registerSubsystem(sr_bg);

  in->fileopen(fname);
  se->run(nevent);
  se->End();
  delete se;
  return 0;
}
