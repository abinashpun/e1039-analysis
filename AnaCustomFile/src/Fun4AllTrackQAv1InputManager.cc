#include "Fun4AllTrackQAv1InputManager.h"
#include <interface_main/SQHit_v1.h>
#include <interface_main/SQHitVector_v1.h>
#include <interface_main/SQEvent_v1.h>
#include <interface_main/SQRun_v1.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllSyncManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllUtils.h>
#include <ffaobjects/RunHeader.h>
#include <ffaobjects/SyncObjectv2.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/recoConsts.h>
#include <cstdlib>
#include <memory>
#include <TFile.h>
#include <TTree.h>
using namespace std;

Fun4AllTrackQAv1InputManager::Fun4AllTrackQAv1InputManager(const string& name, const string& topnodename)
  : Fun4AllInputManager(name, "")
  , segment(-999)
  , isopen(0)
  , events_total(0)
  , events_thisfile(0)
  , topNodeName(topnodename)
  , _tree_name("QA_ana")
  , run_header(nullptr)
  , event_header(nullptr)
  , hit_vec(nullptr)
{
  Fun4AllServer* se = Fun4AllServer::instance();
  topNode = se->topNode(topNodeName.c_str());
  PHNodeIterator iter(topNode);

  PHCompositeNode* runNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode) {
    runNode = new PHCompositeNode("RUN");
    topNode->addNode(runNode);
  }

  PHCompositeNode* eventNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!eventNode) {
    eventNode = new PHCompositeNode("DST");
    topNode->addNode(eventNode);
  }

  run_header = new SQRun_v1();
  runNode->addNode(new PHIODataNode<PHObject>(run_header, "SQRun", "PHObject"));
  
  event_header = new SQEvent_v1();
  eventNode->addNode(new PHIODataNode<PHObject>(event_header, "SQEvent", "PHObject"));
  
  hit_vec = new SQHitVector_v1();
  eventNode->addNode(new PHIODataNode<PHObject>(hit_vec, "SQHitVector", "PHObject"));
  
  syncobject = new SyncObjectv2();
}

Fun4AllTrackQAv1InputManager::~Fun4AllTrackQAv1InputManager()
{
  if(isopen) fileclose();
  delete syncobject;
}

int Fun4AllTrackQAv1InputManager::fileopen(const string &filenam)
{
  if (isopen) {
    cout << "Closing currently open file "
         << filename
         << " and opening " << filenam << endl;
    fileclose();
  }
  filename = filenam;
  if (verbosity > 0) {
    cout << ThisName << ": opening file " << filename.c_str() << endl;
  }
  
  events_thisfile = 0;

  _fin = TFile::Open(filenam.c_str(), "READ");
  if (! _fin->IsOpen()) {
    cerr << "!!ERROR!! Failed at opening the input file (" << filenam << ").  Exit.\n";
    return -1;
  }
  _tin = (TTree*) _fin->Get(_tree_name.c_str());
  if (!_tin) {
    cerr << "!!ERROR!! Failed at getting the input tree (" << _tree_name << ").  Exit.\n";
    return -1;
  }

  _tin->SetBranchAddress("nHits"        , &nHits);
  _tin->SetBranchAddress("detectorID"   , &detectorID);
  _tin->SetBranchAddress("elementID"    , &elementID);
  _tin->SetBranchAddress("driftDistance", &driftDistance);

  segment = 0;
  isopen = 1;
  AddToFileOpened(filename); // add file to the list of files which were opened
  return 0;
}

int Fun4AllTrackQAv1InputManager::run(const int nevents)
{
readagain:
  if (!isopen) {
    if (filelist.empty()) {
      if (verbosity > 0) {
        cout << Name() << ": No Input file open" << endl;
      }
      return -1;
    } else {
      if (OpenNextFile()) {
        cout << Name() << ": No Input file from filelist opened" << endl;
        return -1;
      }
    }
  }
  if (verbosity > 3) {
    cout << "Getting Event from " << Name() << endl;
  }
  
  if (events_thisfile >= _tin->GetEntries()) {
    fileclose();
    goto readagain;
  }
  
  _tin->GetEntry(events_thisfile);
  events_thisfile++;
  events_total++;
  
  int run_id   = 0;
  int spill_id = 0;
  int event_id = events_total;
  
  SetRunNumber                (run_id);
  mySyncManager->PrdfEvents   (events_thisfile);
  mySyncManager->SegmentNumber(spill_id);
  mySyncManager->CurrentEvent (0);
  
  syncobject->RunNumber       (run_id);
  syncobject->EventCounter    (events_thisfile);
  syncobject->SegmentNumber   (spill_id);
  syncobject->EventNumber     (event_id);

  run_header->set_run_id(run_id);

  event_header->set_run_id  (run_id);
  event_header->set_spill_id(spill_id);
  event_header->set_event_id(event_id);
  event_header->set_data_quality(0);
  event_header->set_trigger(0);
  event_header->set_qie_turn_id(0);
  event_header->set_qie_rf_id  (0);
  for(int i = -16; i < 16; ++i) event_header->set_qie_rf_intensity(i, 0);

  for (auto i_hit = 0; i_hit < nHits; ++i_hit) {
    SQHit* hit = new SQHit_v1();
    hit->set_hit_id(i_hit + 1);
    hit->set_detector_id(detectorID[i_hit]);
    hit->set_element_id (elementID [i_hit]);
    hit->set_tdc_time   (0);
    hit->set_drift_distance(driftDistance[i_hit]);
    hit->set_in_time(true);
    hit->set_pos(0);
    hit_vec->push_back(hit);
    delete hit;
  }

  // check if the local SubsysReco discards this event
  if (RejectEvent() != Fun4AllReturnCodes::EVENT_OK) {
    ResetEvent();
    goto readagain;
  }

  return 0;
}

int Fun4AllTrackQAv1InputManager::fileclose()
{
  if (!isopen) {
    cout << Name() << ": fileclose: No Input file open" << endl;
    return -1;
  }
  
  _fin->Close();
  isopen = 0;
  
  // if we have a file list, move next entry to top of the list
  // or repeat the same entry again
  if (!filelist.empty()) {
    if (repeat) {
      filelist.push_back(*(filelist.begin()));
      if (repeat > 0) {
        repeat--;
      }
    }
    filelist.pop_front();
  }
  
  return 0;
}

void
Fun4AllTrackQAv1InputManager::Print(const string &what) const
{
  Fun4AllInputManager::Print(what);
  return ;
}

int
Fun4AllTrackQAv1InputManager::OpenNextFile()
{
  while (!filelist.empty()) {
    list<string>::const_iterator iter = filelist.begin();
    if (verbosity) {
      cout << PHWHERE << " opening next file: " << *iter << endl;
    }
    if (fileopen(*iter)) {
      cout << PHWHERE << " could not open file: " << *iter << endl;
      filelist.pop_front();
    } else {
      return 0;
    }
  }
  return -1;
}

int
Fun4AllTrackQAv1InputManager::ResetEvent()
{
  syncobject->Reset();
  return 0;
}

int
Fun4AllTrackQAv1InputManager::PushBackEvents(const int i)
{
  cerr << "!!ERROR!!  PushBackEvents():  Not implemented yet." << endl;
  return -1;
}

int
Fun4AllTrackQAv1InputManager::GetSyncObject(SyncObject **mastersync)
{
  // here we copy the sync object from the current file to the
  // location pointed to by mastersync. If mastersync is a 0 pointer
  // the syncobject is cloned. If mastersync allready exists the content
  // of syncobject is copied
  if (!(*mastersync)) {
    if (syncobject) *mastersync = syncobject->clone();
  } else {
    *(*mastersync) = *syncobject; // copy syncobject content
  }
  return Fun4AllReturnCodes::SYNC_OK;
}

int
Fun4AllTrackQAv1InputManager::SyncIt(const SyncObject *mastersync)
{
  if (!mastersync) {
    cout << PHWHERE << Name() << " No MasterSync object, cannot perform synchronization" << endl;
    cout << "Most likely your first file does not contain a SyncObject and the file" << endl;
    cout << "opened by the Fun4AllDstInputManager with Name " << Name() << " has one" << endl;
    cout << "Change your macro and use the file opened by this input manager as first input" << endl;
    cout << "and you will be okay. Fun4All will not process the current configuration" << endl << endl;
    return Fun4AllReturnCodes::SYNC_FAIL;
  }
  int iret = syncobject->Different(mastersync);
  if (iret) {
    cout << "big problem" << endl;
    exit(1);
  }
  return Fun4AllReturnCodes::SYNC_OK;
}
