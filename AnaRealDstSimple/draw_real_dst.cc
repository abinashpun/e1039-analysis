/// ROOT macro to look into the E1039 DST file in a simple way.
/** 
 * This macro uses TTree::Draw() to plot the contents of the E1039 DST file,
 * rather than the Fun4All event loop.  Anyone can try it for quick checks,
 * but not recommended to use it for serious analyses because it cannot 
 * make use of the mapping and calibration info of the e1039-core framework.
 *
 * It still needs the library file of the e1039-core framework to interpret 
 * the SQ data structure.  On "seaquestdaq01" and "e1039-monitor", you can 
 * just execute the following command before running this macro;
 *   source /data2/e1039/this-e1039.sh
 *
 * The variable name that you can use in TTree::Draw() is given in the 
 * e1039-core header files, e.g. 
 *   /data2/e1039/core/include/interface_main/SQHit.h
 */
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
R__LOAD_LIBRARY(libinterface_main)
#endif

TTree* GetTree(const int run, const int spill=0);
void PrintTreeElements(TTree* tree);

/// Main function.
void draw_real_dst(const int run=318, const int spill=0)
{
  gSystem->Load("libinterface_main.so");

  TTree* tree = GetTree(run, spill);
  PrintTreeElements(tree); // Call this function to find all variable names.
  tree->SetAlias("event", "DST.SQEvent");
  tree->SetAlias("hit"  , "DST.SQHitVector._vector");
  //tree->SetAlias("trk"  , "DST.SRecEvent.fAllTracks");

  TCanvas* c1 = new TCanvas("c1", "");
  c1->SetGrid();

  //// Example 1:  Number of hits per event
  //tree->Draw("DST.SQHitVector.@_vector.size()", "");
  //c1->SaveAs("n_hits.png");

  //// Example 2:  TDC time with trigger-type and detector selections
  tree->Draw("hit.get_tdc_time()",
             "event._trigger == 0x4 && hit.get_detector_id() == 13");
  c1->SaveAs("tdc_time.png");

  //// Example 3:  Element ID with trigger-type and detector selections
  tree->Draw("hit.get_element_id()",
             "event._trigger == 0x4 && hit.get_detector_id() == 13");
  c1->SaveAs("ele_id.png");

  //// Example 3:  Tracks.  Only works on DST afte reconstruction.
  //tree->Draw("trk.fChisq");
  //c1->SaveAs("trk_chi2.png");

  exit(0);
}

/// Sub-function to open the ROOT file and get the tree.
TTree* GetTree(const int run, const int spill)
{
  ostringstream oss;
  oss << setfill('0') << setw(6) << run;
  string run6 = oss.str();
  oss.str("");
  oss << setw(9) << spill;
  string spill9 = oss.str();

  oss.str("");
  oss << "/data2/e1039/dst/run_" << run6 << "/run_" << run6 << "_spill_" << spill9 << "_spin.root";
  string fn_in = oss.str();
  cout << "DST = " << fn_in << endl;
  TFile* file = new TFile(fn_in.c_str());
  if (! file->IsOpen()) {
    cout << "  Cannot open the DST file.  Abort." << endl;
    exit(1);
  }
  TTree* tree = (TTree*)file->Get("T");
  if (! tree) {
    cout << "  Cannot get the tree object.  Abort." << endl;
    exit(1);
  }
  return tree;
}

/// Sub-function to just print out the list of tree variables.
void PrintTreeElements(TTree* tree)
{
  TObjArray* arr = tree->GetListOfLeaves();
  int n_ent = arr->GetEntries();
  for (int i_ent = 0; i_ent < n_ent; i_ent++) {
    TObject* obj = arr->At(i_ent);
    cout << setw(2) << i_ent << " " << obj->ClassName() << " " << obj->GetName() << "\n";
  }
  cout << endl;
  
  //TLeaf* leaf = tree->GetLeaf("DST.SQHitVector._vector");
  //leaf->Print();
  //TObject* cl = leaf->GetValuePointer();
  //cout << leaf->GetTypeName() << " " << cl->GetName() << endl;
  //leaf->PrintValue(10);
  //exit(0);
}
