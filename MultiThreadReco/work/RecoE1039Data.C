/// ROOT macro to run the online reconstruction.
/** 
 * Based on `RunScheduler.C`.
 */
//R__LOAD_LIBRARY(libinterface_main)
R__LOAD_LIBRARY(libcalibrator)
//R__LOAD_LIBRARY(libktracker)
R__LOAD_LIBRARY(MultiThreadReco)

/**
 * Not tested at all.
 */
int RecoE1039Data(const int run_id, const int spill_id, TString infile, TString outfile, const int n_evt=100)
{
    ROOT::EnableThreadSafety();
    std::cout << "Start of kTThreadTest" << std::endl;

    const bool cosmic = true;
    recoConsts* rc = recoConsts::instance();
    rc->set_IntFlag("RUNNUMBER", run_id);
    rc->set_CharFlag("EventReduceOpts", "a"); // Better set "none" if possible
    if (cosmic) {
        rc->init("cosmic");
        rc->set_BoolFlag("COARSE_MODE", true);
        rc->set_DoubleFlag("KMAGSTR", 0.);
        rc->set_DoubleFlag("FMAGSTR", 0.);
    } else {
      const double FMAGSTR = -1.054;
      const double KMAGSTR = -0.951;
      rc->set_DoubleFlag("FMAGSTR", FMAGSTR);
      rc->set_DoubleFlag("KMAGSTR", KMAGSTR);
    }
    rc->Print();

    // Start the scheduler call
    std::cout << "Initializing Scheduler" << std::endl;
    //KJob      ::Verbose(1);
    //KScheduler::Verbose(1);
    KScheduler* ksched = new KScheduler(infile, outfile);
    ksched->UseTrackletReco(true);
    ksched->Init();

    std::cout << "Master call to run the threads and join" << std::endl;
    ksched->runThreads();

    std::cout << "All threads dead! all done... bye bye!\n";

    delete ksched;
    return 0;
}
