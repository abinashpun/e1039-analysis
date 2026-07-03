R__LOAD_LIBRARY(interface_main)

int WriteDstTree()
{
  string fn_dst = "DSTx.root";
  TFile* file = new TFile(fn_dst.c_str(), "RECREATE");

  ///
  /// Create the RUN node.
  ///
  TTree* tree_run = new TTree("T1", "Created by WriteDstTree.C.");
  SQRun* sq_run = new SQRun_v1();
  tree_run->Branch("RUN.SQRun", &sq_run);
  sq_run->set_run_id(1234);
  sq_run->set_unix_time_begin(10);
  sq_run->set_unix_time_end  (20);
  tree_run->Fill();

  ///
  /// Create the DST node.
  ///
  TTree* tree = new TTree("T", "Created by WriteDstTree.C.");
  SQEvent*     sq_evt          = new SQEvent_v2();
  SQHitVector* sq_hit_vec      = new SQHitVector_v1();
  SQHitVector* sq_trig_hit_vec = new SQHitVector_v1();
  tree->Branch("DST.SQEvent"           , &sq_evt);
  tree->Branch("DST.SQHitVector"       , &sq_hit_vec);
  tree->Branch("DST.SQTriggerHitVector", &sq_trig_hit_vec);

  const int n_evt = 100;
  for (int i_evt = 0; i_evt < n_evt; i_evt++) {
    sq_evt->Reset();
    sq_evt->set_event_id(i_evt + 1);
    sq_evt->set_trigger(i_evt % 8);

    sq_hit_vec->Reset();
    const size_t n_hit = 10;
    for (size_t i_hit = 0; i_hit < n_hit; i_hit++) {
      SQHit* hit = new SQHit_v1();
      hit->set_hit_id     (i_hit+1);
      hit->set_detector_id(i_hit+2);
      hit->set_element_id (i_hit+10);
      hit->set_tdc_time   (i_hit+900.);
      hit->set_drift_distance(0);
      hit->set_pos           (0);
      sq_hit_vec->push_back(hit);
      delete hit;
    }

    sq_trig_hit_vec->Reset();
    const size_t n_trig_hit = 10;
    for (size_t i_hit = 0; i_hit < n_trig_hit; i_hit++) {
      SQHit* hit = new SQHit_v1();
      hit->set_hit_id     (i_hit+1);
      hit->set_detector_id(i_hit+31);
      hit->set_element_id (i_hit+1);
      hit->set_tdc_time   (i_hit+500.);
      hit->set_drift_distance(0);
      hit->set_pos           (0);
      sq_trig_hit_vec->push_back(hit);
      delete hit;
    }

    tree->Fill();
  }

  file->Write();
  file->Close();

  exit(0);
}
