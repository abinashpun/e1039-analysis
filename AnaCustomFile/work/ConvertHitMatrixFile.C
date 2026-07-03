R__LOAD_LIBRARY(libAnaCustomFile)

/// Fun4All macro to convert a hit-matrix file to DST and SRawEvent files.
/**
 * The input parameters like file name are given via macro argument.
 * You can change them in command line or in this file itself;
 * @code
 * root -b 'ConvertHitMatrixFile.C(0)' // or
 * root -b 'ConvertHitMatrixFile.C(0, "/path/to/input/file")'
 * @endcode
 *
 * If you don't need one of the two outputs, you can comment out `out1` or `out2` to save time and disk space.
 */
int ConvertHitMatrixFile(
  const int n_evt=10,
  const string fn_in="/project/ptgroup/Kenichi/tmp/2022-03-20/output_merged.root",
  const string tree_name="tree_nim3",
  const string branch_name="hit_matrix",
  const string fn_out1="DST.root",
  const string fn_out2="srawevent.root")
{
  recoConsts* rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", 1); // To select the right plane geometry.
  //rc->Print();

  Fun4AllServer* se = Fun4AllServer::instance();
  se->Verbosity(0);

  auto calib_hit_pos = new CalibHitPos();
  //calib_hit_pos->Verbosity(1); // 0 = none, 1 = old+new pos
  se->registerSubsystem(calib_hit_pos);

  /// Use this for debug.
  //auto check_evt = new CheckEventData();
  //check_evt->Verbosity(0); // 0 = none, 1 = det:ele, 2 = det:ele[dist]
  //se->registerSubsystem(check_evt);

  auto in = new Fun4AllHitMatrixInputManager("IN");
  in->Verbosity(0);
  in->SetTreeName(tree_name);
  in->SetBranchName(branch_name);
  in->fileopen(fn_in);
  se->registerInputManager(in);

  auto out1 = new Fun4AllDstOutputManager("DSTOUT", fn_out1);
  se->registerOutputManager(out1);

  auto out2 = new Fun4AllSRawEventOutputManagerV2("SRAWEVENTOUT", fn_out2);
  se->registerOutputManager(out2);

  se->run(n_evt);

  se->End();
  se->PrintTimer();
  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
  return 0;
}
