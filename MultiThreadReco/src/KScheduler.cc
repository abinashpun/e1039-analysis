/*
 * KScheduler.cxx
 * Scheduler for KEventJobs in the KJobQueue
 * Author: Noah Wuerfel, nwuerfel@umich.edu
 * created: 10/14/2020
 */
#include <iostream>
#include <algorithm>
#include <cmath>
#include <interface_main/SQHitVector.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQRun.h>
#include <interface_main/SQSpillMap.h>
#include <interface_main/SQTrackVector.h>
#include <phool/recoConsts.h>
#include <phfield/PHFieldConfig_v3.h>
#include <phfield/PHFieldUtility.h>
#include <phgeom/PHGeomUtility.h>
#include <ktracker/KalmanFastTrackletting.h>
#include "KJob.h"
#include "KScheduler.h"
using namespace std;

//sure there's a better way for passing compound args...
// thread_local may not be right?.... root is annoying
struct workerArg{
    KScheduler* kschdPtr;
    unsigned threadId;
};

// FUCK this was tricky need to initialize blecause static
// need static to acess via the member threads because pthreads are typically
// static code
int KScheduler::verb = 0;
int KScheduler::n_threads = 16;
int KScheduler::input_pipe_depth  = 32; // play with iobuffer to load memory vs do computation...
int KScheduler::output_pipe_depth = 32;
int KScheduler::print_freq = 100; // num events per timeupdate TODO post fraction need update from fReaderThread
int KScheduler::save_num = 5000;
bool KScheduler::save_raw_evt = false;

//bool KScheduler::use_e906_data = true;
TString KScheduler::inputFilename="";
TString KScheduler::outputFilename="";
int KScheduler::completedEvents = 0;

KScheduler::KScheduler(TString inFile, TString outFile) 
  : use_tracklet_reco(false)
{
    this->setInputFilename(inFile);
    this->setOutputFilename(outFile);
}

void KScheduler::Init(PHField* phfield, TGeoManager* t_geo_manager, KalmanFitter* kfitter, SQGenFit::GFFitter* gfitter, SQGenFit::GFField* gfield)
{
    std::cout << "Initialization of KScheduler ..." << std::endl;
    std::cout << "================================" << std::endl;

    // init services
    std::cout << "starting KScheduler services\n";
    //GeomSvc* p_geomSvc = GeomSvc::instance();
    recoConsts* rc = recoConsts::instance();

    //Initialize event reducer
    string evt_red_opt = rc->get_CharFlag("EventReduceOpts");
    for (int i = 0; i < n_threads; i++) {
      EventReducer* eventReducer = 0;
      if (evt_red_opt != "none") {
        eventReducer = new EventReducer(evt_red_opt);
        assert(eventReducer);
      }
      eventReducerQueue.push(eventReducer);
      vec_eventReducer.push_back(eventReducer);
    }

    //Initialize the kfasttrackers
    KalmanFastTracking* kFastTracker = 0;
    for (int i = 0; i < n_threads; i++) {
        if (! use_tracklet_reco) kFastTracker = new KalmanFastTracking    (phfield, t_geo_manager, false);
        else                     kFastTracker = new KalmanFastTrackletting(phfield, t_geo_manager, false);
        assert(kFastTracker);
        kFastTrkQueue.push(kFastTracker);
        vec_kFastTrk.push_back(kFastTracker);
    }

    m_kfitter = kfitter;
    m_gfitter = gfitter;

    // build TClonesArrays for tracklet outputs
    TClonesArray* trackletArray = 0;
    for (int i = 0; i < n_threads; i++) {
        trackletArray= new TClonesArray("Tracklet",1000);
        trackletArray->BypassStreamer();
        assert(trackletArray);
        kTrackArrQueue.push(trackletArray);
    }

    // worker threads
    workThreadArr.resize(n_threads, 0);

    // net time for alljobs
    avgTimer = new TStopwatch();
    totalTimer = new TStopwatch();

    // lock for wArg structure when initializing structs
    wArgMutex = new TMutex();
    fReaderMutex = new TMutex();

    // pipeline inside
    newJobQueuePutMutex = new TMutex();
    newJobQueueTakeMutex = new TMutex();

    // pipeline evReducer
    evRedQueuePutMutex = new TMutex();
    evRedQueueTakeMutex = new TMutex();

    // pipeline finder
    kFTrkQueuePutMutex = new TMutex();
    kFTrkQueueTakeMutex = new TMutex();

    // pipeline out
    cmpJobQueuePutMutex = new TMutex();
    cmpJobQueueTakeMutex = new TMutex();

    //tracklet arr
    ktrkQueueMutex = new TMutex();
   
//    jobSem = new TSemaphore(MAXJOBS);
    njqFSem = new TSemaphore(0);
    njqESem = new TSemaphore(input_pipe_depth);

    // evReducer
    erqFSem = new TSemaphore(n_threads);
    erqESem = new TSemaphore(0);

    // finder
    kftqFSem = new TSemaphore(n_threads);
    kftqESem = new TSemaphore(0);

    // trackletArrays
    ktrkqFSem = new TSemaphore(n_threads);
    ktrkqESem = new TSemaphore(0);

    // output-stage pipeline
    cjqFSem = new TSemaphore(0);
    cjqESem = new TSemaphore(output_pipe_depth);

    // inits TThreads to 0
    fRDPtr = 0;
    fRPPtr = 0;

    return;
};

KScheduler::~KScheduler(){

    int i = 0;

    //TODO kill remaning threads and clean memory for them
    /*
    for(int i = 0; i<n_threads; i++){
        delete jobMutexArr[i];
    }
    */
    
    delete fRDPtr;
    fRDPtr = 0;
    delete fRPPtr;
    fRPPtr = 0;
    for(i=0;i<n_threads;i++){
        assert(workThreadArr[i]);
        delete workThreadArr[i];
        workThreadArr[i]=0;
    }

    delete wArgMutex;
    delete fReaderMutex;

    // first stage of pipeline
    delete newJobQueuePutMutex;
    delete newJobQueueTakeMutex;

    // eventReducer
    delete evRedQueuePutMutex;
    delete evRedQueueTakeMutex;
    // TODO CLEAN UP EVENTREDUCERS
    EventReducer* er = 0;
    for(i=0;i<n_threads;i++){
        er = eventReducerQueue.front(); 
        eventReducerQueue.pop();
        delete er;
    }
    //delete eventReducerQueue;
    
    // finder
    delete kFTrkQueuePutMutex;
    delete kFTrkQueueTakeMutex;

    KalmanFastTracking* ft = 0;
    for(i=0; i<n_threads; i++){
        ft = kFastTrkQueue.front();
        kFastTrkQueue.pop();
        delete ft;
    }

    TClonesArray* trkArr= 0;
    for(i=0; i<n_threads; i++){
        trkArr = kTrackArrQueue.front();
        assert(trkArr);
        kTrackArrQueue.pop();
        delete trkArr;
    }

    // second stage of pipeline
    delete cmpJobQueuePutMutex;
    delete cmpJobQueueTakeMutex;

    //tracklet arr mutex
    delete ktrkQueueMutex;
    
    // first stage of pipeline sems
    delete njqFSem;
    delete njqESem;

    // ev reducer
    delete erqFSem;
    delete erqESem;

    // finder
    delete kftqFSem;
    delete kftqESem;

    // tracklet arrays sems
    delete ktrkqFSem;
    delete ktrkqESem;

    //output stage of pipeline sems
    delete cjqFSem;
    delete cjqESem;

    // timer
    delete avgTimer;
    TThread::Printf("Total Time and average per event below:");
    // something with printing or formatting is wrong here...
    // yeah the formatting is stupid.. I pulled this from TStopwatch::Print()...
//    TThread::Printf("Total: %d, AvgPerEntry: %d",ttimeelapsed,ttimeelapsed/(double)completedEvents);
    Double_t  realt = totalTimer->RealTime();
    Double_t  cput  = totalTimer->CpuTime();
    Double_t avgt = totalTimer->RealTime()/(double)completedEvents;
 
    Int_t  hours = Int_t(realt / 3600);
    realt -= hours * 3600;
    Int_t  min   = Int_t(realt / 60);
    realt -= min * 60;
    //Int_t  sec   = Int_t(realt);

    Int_t  avhours = Int_t(avgt / 3600);
    avgt -= avhours * 3600;
    Int_t  avmin   = Int_t(avgt / 60);
    avgt -= avmin * 60;
    //Int_t  avsec   = Int_t(avgt);
 
    if (realt < 0) realt = 0;
    if (cput  < 0) cput  = 0;
    if (totalTimer->Counter() > 1) {
          TThread::Printf("Real time %d:%02d:%06.3f, CP time %.3f, %d slices", 
            hours, min, realt, cput, totalTimer->Counter());
    } 
    else {
          TThread::Printf("Real time %d:%02d:%06.3f, CP time %.3f", 
            hours, min, realt, cput);
          TThread::Printf("Average time per event: %d:%02d:%06.3f", 
            avhours, avmin, avgt);
    }


//    totalTimer->Print("m");

//   TThread::Printf("total events: %i",completedEvents);
    delete totalTimer;
    

    return;
}

Int_t KScheduler::runThreads(){
    Int_t ret; 
    std::cout << "KScheduler spawning threads..." << std::endl;
//    std::cout << "KScheduler spawning readerThread..." << std::endl;
//    ret = this->startReaderThread();
//    assert(ret == 0);
    std::cout << "KScheduler spawning reaperThread..." << std::endl;
    ret = this->startReaperThread();
    assert(ret == 0);
    std::cout << "KScheduler spawning worker threads..." << std::endl;
    ret = this->startWorkerThreads();
    assert(ret == 0);
    return 0;
}

Int_t KScheduler::endThreads()
{
  //TThread::Join(fRDPtr->GetId(),NULL);
  //delete fRDPtr;
  //fRDPtr = 0;
  TThread::Join(fRPPtr->GetId(),NULL);
  delete fRPPtr;
  fRPPtr = 0;
  TThread* wPtr = 0;
  for (int i=0; i < n_threads; i++) {
    wPtr = workThreadArr[i];
    TThread::Join(wPtr->GetId(),NULL);
  }
  return 0;
}

// static getters and setters
TString KScheduler::getInputFilename(){
    return inputFilename;
}
void KScheduler::setInputFilename(TString name){
    inputFilename = name;
}
// static getters and setters
TString KScheduler::getOutputFilename(){
    return outputFilename;
}

void KScheduler::setOutputFilename(TString name){
    outputFilename = name;
}

void KScheduler::postCompletedEvent(){
  if(Verbose() > 0 && completedEvents % print_freq == 0){
    avgTimer->Stop();
    TThread::Printf("completed: %i events, last %i in time shown below:",
                    completedEvents, print_freq);
    // uhoh gap here?
    avgTimer->Print("u");
    avgTimer->Start();
    //TThread::Printf("\n");
    
  }
  completedEvents++;
}

// allocates memory for KJOBs - need to delete them in reaperThread
/**
 * This function is no longer used.
 * Use MultiThreadReco::process_event() instead.
 */
//void* KScheduler::fReaderThread(void* readerArgPtr){
//
//    int i;
//
//    // wack... outdated instructions for ROOT?
//    KScheduler* kschd = (KScheduler*) readerArgPtr;
//    TString filename = kschd->getInputFilename();
//    if(Verbose() > 0){
//        TThread::Printf("Start of fReaderThread - ifile is: %s\n", filename.Data());
//    }
//
//    // to be filled from Tree
//    // or from converter in e1039
//    SRawEvent* localEventPtr=0;
//
//    // in case we have E1039 data
//    SQEvent* SQEvPtr = 0;
//    SQHitVector* SQHitVecPtr = 0;
//    SQHitVector* SQTrigHitVecPtr = 0;
//
//    // rules of synch for ROOT obj are still not clear, so i just sync
//    // everything.. Not sure if i need to use global mutex for aLL root obj or
//    // local ones are good enough with the EnableThreadSafety feature of R6...
//    // if its global i'm boned unless i set up a sleeping system for queue
//    // filling... 
//    //
//    // GetEntries does make a new TObject... BUT IN ROOT6
//    // ENABLETHREADSAFETYCLAIMSTHAT THIS IS ACCEPTABLEEEEEEEEE.... yikes
//
//
//    // just lock everything for now... performace issues are better than race
//    // cond.
//    kschd->fReaderMutex->Lock();
//    TFile* inputFile = new TFile(filename, "READ"); 
//    TTree* dataTree = 0;
//    if(use_e906_data){
//        dataTree = (TTree*) inputFile->Get("save");
//        dataTree->SetBranchAddress("rawEvent",&localEventPtr);
//    }
//    // e1039 data formatting...
//    else{
//        dataTree = (TTree*) inputFile->Get("T");
//        dataTree->SetBranchAddress("DST.SQEvent",&SQEvPtr);
//        dataTree->SetBranchAddress("DST.SQHitVector", &SQHitVecPtr);
//        dataTree->SetBranchAddress("DST.SQTriggerHitVector", &SQTrigHitVecPtr);
//    }
//    int nEvents = dataTree->GetEntries();
//    TThread::Printf("fReaderThread: starting read of %i events...\n",nEvents);
//    kschd->fReaderMutex->UnLock();
//
//    // KJob pointer to add to queue 
//    KJob* newKJobPtr = NULL;
//
//    bool copy = true;
//
//    // TThread sched_yield? for IO block?
//    for(i=0; i<nEvents; i++){
//        dataTree->GetEntry(i);
//
//        // convert SQHitVector to SRawEvent
//        // TODO memory leak from new alloc in BuildSRawEvent?
//        if(! use_e906_data){
//            localEventPtr = KScheduler::BuildSRawEvent(SQEvPtr, SQHitVecPtr, SQTrigHitVecPtr);    
//            copy = false;
//        }
//
//        // makes a copy of the SRawEvent
//        // takes a mutex when allocating memory for the event
//        newKJobPtr = new KJob(i,localEventPtr,kschd,copy);
//
//        // need shared queue class
//        // for now theres just one mutex....
//        kschd->njqESem->Wait();
//        kschd->newJobQueuePutMutex->Lock();
//        kschd->newJobQueue.push(newKJobPtr);
//        kschd->newJobQueuePutMutex->UnLock();
//        kschd->njqFSem->Post();
//    }
//
//    TThread::Printf("fReaderThread: finished reading all events... dumping poison...");
//
//    for(i=0;i<n_threads;i++){
//        if(Verbose() > 0){
//            TThread::Printf("poisoning...");
//        }
//        newKJobPtr = new KJob(true);
//        kschd->njqESem->Wait();
//        kschd->newJobQueuePutMutex->Lock();
//        kschd->newJobQueue.push(newKJobPtr);
//        kschd->newJobQueuePutMutex->UnLock();
//        kschd->njqFSem->Post();
//    }
//    return 0;
//}

void KScheduler::PushEvent(SRawEvent* sraw, bool copy)
{
  static int job_id = 0;
  KJob* newKJobPtr = new KJob(job_id++, sraw, this, copy);

  // need shared queue class
  // for now theres just one mutex....
  njqESem->Wait();
  newJobQueuePutMutex->Lock();
  newJobQueue.push(newKJobPtr);
  newJobQueuePutMutex->UnLock();
  njqFSem->Post();
}

void KScheduler::PushPoison()
{
  TThread::Printf("fReaderThread: finished reading all events... dumping poison...");
  for (int i = 0; i < n_threads; i++) {
    if(Verbose() > 2) TThread::Printf("  poisoning...");
    KJob* newKJobPtr = new KJob(true);
    njqESem->Wait();
    newJobQueuePutMutex->Lock();
    newJobQueue.push(newKJobPtr);
    newJobQueuePutMutex->UnLock();
    njqFSem->Post();
  }
}


// reaper thread
void* KScheduler::fReaperThread(void* reaperArg){

    KScheduler* kschd = (KScheduler*) reaperArg;
    //recoConsts* rc = recoConsts::instance();

//    TThread::Printf("worker #%u got ptr:%p",dwArgPtr->threadId,wArgPtr);
    TThread::Printf("Starting fReaper thread\n");

    // output fields for saveTree->Fill
    TString outputFilename = kschd->getOutputFilename();
    SRawEvent* outputRawEventPtr = 0;
    SRecEvent* outputRecEventPtr = 0;
    TClonesArray* outputTracklets = new TClonesArray("Tracklet",1000);
    int nTracklets = 0;
    int sizeTrackArr = 0;

    // setup the output files...
    TThread::Printf("opening output file\n");
    TFile* saveFile = new TFile(outputFilename.Data(), "recreate"); 
    TTree* saveTree = new TTree("save","save");
    if (save_raw_evt) saveTree->Branch("rawEvent", &outputRawEventPtr, 256000, 99);
    saveTree->Branch("recEvent", &outputRecEventPtr, 256000, 99);
    //saveTree->Branch("time", &time, "time/D");
    saveTree->Branch("outputTracklets", &outputTracklets, 256000, 99);
    saveTree->Branch("nTracklets", &nTracklets, "nTracklets/I");
    saveTree->Branch("sizeTrackArr", &sizeTrackArr, "nTracklets/I");

    // try to reap jobs forever
    int poisonPills = 0;
    bool running = true;
    while(running){
        // try to acquire a job from the queue...
        kschd->cjqFSem->Wait();
        kschd->cmpJobQueuePutMutex->Lock();
        KJob* tCompleteJobPtr = kschd->cmpJobQueue.front();
        kschd->cmpJobQueue.pop();
        kschd->cmpJobQueuePutMutex->UnLock();
        kschd->cjqESem->Post();

        // check for poison
        if(tCompleteJobPtr->isPoison){
          // cleanup job and die... 
          if(Verbose() > 1) TThread::Printf("ReaperThread got a poison pill...");
          delete tCompleteJobPtr;
          poisonPills++;
          if (poisonPills == n_threads) {
            TThread::Printf("ReaperThread caught all pills...");
            break;
          } else {
            continue;
          }
        }
        
        // check for halted thread
        if(tCompleteJobPtr->p_JobStatus==HALTED){
           // job failed for some reason... 
           delete tCompleteJobPtr;
           continue;
        }
        
        if(Verbose() > 1){
            tCompleteJobPtr->jobMutex->Lock();
            TThread::Printf("fReaper gets jobId: %i, eventID: %i", 
                tCompleteJobPtr->jobId, tCompleteJobPtr->evData->getEventID());
            tCompleteJobPtr->jobMutex->UnLock();
        }

        // TODO CHECK JOB STATUS 
        // otherwise no atomic way to safely delete... because i'd be deleting
        // mutex too...
        // TODO PRAY .... is there a race conidition here? I'm 99% sure the
        // pipeline enforces enough serialism here tobe safe
        //
        // TODO UNLESS there's a 3rd party daemon that checsk events or jobs...
        // maybe when doing time scheduling? ... careful here...

        /// Perform the track fitting here, because KalmanFitter/GFFitter are not thread-safe.
        kschd->DoTrackFitting(tCompleteJobPtr->tracklets, tCompleteJobPtr->recEvData);

        outputRawEventPtr = tCompleteJobPtr->evData;
        outputRecEventPtr = tCompleteJobPtr->recEvData;
        //outputTracklets = tCompleteJobPtr->tracklets;

        //outputTracklets = new TClonesArray(*(tCompleteJobPtr->tracklets));
        *outputTracklets = *(tCompleteJobPtr->tracklets);

        TClonesArray& ref_output = *outputTracklets;
        nTracklets = tCompleteJobPtr->nTracklets;
        sizeTrackArr = outputTracklets->GetEntries();
        assert(outputRawEventPtr);
        assert(outputRecEventPtr);
        assert(outputTracklets);
        if(Verbose() > 1){
            TThread::Printf("got tracklets: %i for disk for event: %i\n",outputTracklets->GetEntries(), tCompleteJobPtr->evData->getEventID());
            TThread::Printf("outputTracklets pointer is: %p\n", outputTracklets);
        
            if(outputTracklets->GetEntries()>0){
                Tracklet* printer = (Tracklet*) ref_output[0];
                TThread::Printf("first tracklet for disk is%p\n",printer);
                if(printer) printer->print();
            }
        }

        saveTree->Fill();
        if(saveTree->GetEntries() % save_num == 0){
            TThread::Printf("fReaper saving another %i jobs", save_num);
            saveTree->AutoSave("SaveSelf");
        }
        outputTracklets->Clear();
        // need to reuse these job objects because tracklets are large...
        delete tCompleteJobPtr;
        kschd->postCompletedEvent();
    }

    TThread::Printf("ReaperThread attempting to save tree");
    // save outputs
    saveTree->AutoSave("SaveSelf");
    saveFile->cd();
    saveTree->Write();
    saveFile->Close();
    
    // cleanup
    // oh no how do i clean these all up...
    // outputTracklets? do i need to delete shallow copy?
    return 0;
}

// need to delete the memory for the worker thread for now
void* KScheduler::fWorkerThread(void* wArgPtr)
{
    //workerArg lwArg;
    workerArg* dwArgPtr = (workerArg*) wArgPtr;
    unsigned threadId = dwArgPtr->threadId;
    KScheduler* kschd = dwArgPtr->kschdPtr;
    //recoConsts* rc = recoConsts::instance();
    //KJob* tCompleteJobPtr = NULL;

    // worker tool pointers
    //int nEvents_dimuon=0;
//    std::list<Tracklet>& rec_tracklets = NULL;
//    TClonesArray& arr_tracklets = NULL;
//    TThread::Printf("worker #%u got ptr:%p",dwArgPtr->threadId,wArgPtr);
    TThread::Printf("Starting worker thread%u\n",threadId);

    // try to get jobs forever
    bool running = true;
    while(running){

        // try to acquire a job from the queue...
        kschd->njqFSem->Wait();
        kschd->newJobQueuePutMutex->Lock();
        KJob* tCompleteJobPtr = kschd->newJobQueue.front();
        assert(tCompleteJobPtr);
        kschd->newJobQueue.pop();
        kschd->newJobQueuePutMutex->UnLock();
        kschd->njqESem->Post();

        // TODO check for poison pill
        if(tCompleteJobPtr->isPoison){
          if(Verbose() > 2) TThread::Printf("WorkerThread got a poison pill...");
          running = false;
          //put job in complete queue to kill next part of pipeline
          kschd->cjqESem->Wait();
          kschd->cmpJobQueuePutMutex->Lock();
          kschd->cmpJobQueue.push(tCompleteJobPtr);
          kschd->cmpJobQueuePutMutex->UnLock();
          kschd->cjqFSem->Post();
          break;
        }
       
        // check poison stream so no crash
        assert(tCompleteJobPtr->evData);

        // try to acquire a job from the queue...
        /// TODO stuff on the job..
        // semaphore enforces synchronization... (I hope...)
        if(Verbose() > 1){
            TThread::Printf("Worker %u gets jobId: %i, eventID: %i\n", 
                threadId, tCompleteJobPtr->jobId, tCompleteJobPtr->evData->getEventID());
        }

        EventReducer* evReducer = kschd->vec_eventReducer.at(threadId);

        //// acquire an eventReducer...
        //kschd->erqFSem->Wait();
        //kschd->evRedQueuePutMutex->Lock();
        //EventReducer* evReducer = kschd->eventReducerQueue.front();
        //assert(evReducer);
        //kschd->eventReducerQueue.pop();
        //kschd->evRedQueuePutMutex->UnLock();
        //kschd->erqESem->Post();
        //
        //if(Verbose() > 1){
        //    TThread::Printf("Worker %u gets evReducer: %p\n",threadId,evReducer);
        //}

        // reduce the event...
        int n_red = evReducer->reduceEvent(tCompleteJobPtr->evData);
        if(Verbose() > 2){
          int n_hits = tCompleteJobPtr->evData->getNChamberHitsAll();
          TThread::Printf("Worker %u reduced %u hits (%u -> %u).\n", threadId, n_red, n_hits+n_red, n_hits);
        }

        // TODO update hit hinfo for the SQhitvector? needs a lot of
        // bookkeeping...

        //// put eventReducer back in the queue...
        //kschd->erqESem->Wait();
        //kschd->evRedQueuePutMutex->Lock();
        //kschd->eventReducerQueue.push(evReducer);
        //kschd->evRedQueuePutMutex->UnLock();
        //kschd->erqFSem->Post();

        KalmanFastTracking* kFastTracker = kschd->vec_kFastTrk.at(threadId);

        //// now the same for the fast tracker
        //kschd->kftqFSem->Wait();
        //kschd->kFTrkQueuePutMutex->Lock();
        //KalmanFastTracking* kFastTracker = kschd->kFastTrkQueue.front();
        ////kFastTracker->Verbosity(99);
        //assert(kFastTracker);
        //kschd->kFastTrkQueue.pop();
        //kschd->kFTrkQueuePutMutex->UnLock();
        //kschd->kftqESem->Post();

        // do something with the tracker
        if(Verbose() > 1) TThread::Printf("Worker %u gets kFastTracker: %p\n", threadId, kFastTracker);

        // set the event
        int recStatus = kFastTracker->setRawEvent(tCompleteJobPtr->evData); 
        if((recStatus != 0 && Verbose() > 1) || Verbose() > 2) TThread::Printf("kFastTrackRecStatus: %i", recStatus);
        tCompleteJobPtr->recEvData->setRecStatus(recStatus);
        tCompleteJobPtr->recEvData->setRawEvent(tCompleteJobPtr->evData);
        std::list<Tracklet>& rec_tracklets = kFastTracker->getFinalTracklets(); 
        if(Verbose() > 1) {
          TThread::Printf("Worker %u completed setRawEvent: %p, trackletsize = %i\n", threadId, kFastTracker, (int)rec_tracklets.size());
          // tCompleteJobPtr->tracklets->GetEntries()
          //TThread::Printf("job pointer for tracklets is%p\n",tCompleteJobPtr->tracklets);
        }

        TClonesArray* _tracklets = tCompleteJobPtr->tracklets; // Output
        int nTracklets = 0;
        int nFittedTracks = 0;        
        for(std::list<Tracklet>::iterator iter = rec_tracklets.begin(); iter != rec_tracklets.end(); ++iter){
            iter->calcChisq();

            //KalmanFitter*       kfitter = kschd->GetKFitter();
            //SQGenFit::GFFitter* gfitter = kschd->GetGFitter();
            //SRecTrack recTrack;
            //bool fitOK = false;
            //if      (kfitter) fitOK = kschd->fitTrackCand(*iter, kfitter, recTrack);
            //else if (gfitter) fitOK = kschd->fitTrackCand(*iter, gfitter, recTrack);
            //
            //if (!fitOK) {
            //  recTrack = iter->getSRecTrack(kfitter != 0);
            //  recTrack.setKalmanStatus(-1);
            //  //fillRecTrack(recTrack);
            //} else {
            //  ++nFittedTracks;
            //}

            new((*_tracklets)[nTracklets]) Tracklet(*iter);
            ++nTracklets;
            tCompleteJobPtr->nTracklets++;
            //tCompleteJobPtr->recEvData->insertTrack(recTrack);
        }
        if(Verbose() > 1) {
          TThread::Printf("Worker %u: Fitter:  nTracklets = %i, nFittedTracks = %i\n", threadId, nTracklets, nFittedTracks);
        }

        /// Keep all trackets.  idx = 0 (D0/D1), 1 (D2), 2 (D3p/m), 3 (D2+3), 4 (D1+2+3)
//        for (int idx = 0; idx <= 4; idx++) {
//            std::list<Tracklet>& tracklets_temp = kFastTracker->getTrackletList(idx);
//            for(std::list<Tracklet>::iterator iter = tracklets_temp.begin(); iter != tracklets_temp.end(); ++iter){
//                iter->calcChisq();
//                //TODO tracklets
//                new(arr_tracklets[nTracklets++]) Tracklet(*iter);
//                tCompleteJobPtr->nTracklets++;
//            }
//        }
//        if(Verbose() > 0){
//            TThread::Printf("arr_tracklet has %i entries for eventID: %i\n",arr_tracklets.GetEntries(), 
//                tCompleteJobPtr->evData->getEventID());
//            if(arr_tracklets.GetEntries() > 0){
//              Tracklet* printer = (Tracklet*) arr_tracklets[0]; 
//                TThread::Printf("first tracklet is:%p for eventID: %i\n",printer, tCompleteJobPtr->evData->getEventID());
//                if(printer)
//                    printer->print();
//            }
//       }

        //// put tracker back in queue...
        //kschd->kftqESem->Wait();
        //kschd->kFTrkQueuePutMutex->Lock();
        //kschd->kFastTrkQueue.push(kFastTracker);
        //kschd->kFTrkQueuePutMutex->UnLock();
        //kschd->kftqFSem->Post();
      
        //put job in complete queue
        kschd->cjqESem->Wait();
        kschd->cmpJobQueuePutMutex->Lock();
        kschd->cmpJobQueue.push(tCompleteJobPtr);
        kschd->cmpJobQueuePutMutex->UnLock();
        kschd->cjqFSem->Post();
    }

    // figure out right place to do this... malloced in startWorkerThread
    delete dwArgPtr;
    return 0;
}

// takes mem for thred
//Int_t KScheduler::startReaderThread(){
//    // check threadstatus
//    std::cout << "Booting fReaderThread" << std::endl;
//    if(!fRDPtr){
////        fRDPtr = new TThread("fReaderThread", (void(*)(void*)) &fReaderThread, (void*) this );
//        fRDPtr = new TThread("fReaderThread", 
//            (TThread::VoidRtnFunc_t) &fReaderThread, (void*) this );
//
//        fRDPtr->Run();
//        return 0;
//    }
//    return 1;
//}
// frees mem for thred
//Int_t KScheduler::stopReaderThread(){
//    if(fRDPtr){
//        TThread::Delete(fRDPtr);
//        delete fRDPtr;
//        fRDPtr = 0;
//        return 0;
//    }
//    return 1;
//}

// takes mem for thred
Int_t KScheduler::startReaperThread(){
    // check threadstatus
    std::cout << "Booting fReaperThread" << std::endl;
    if(!fRPPtr){
//        fRPPtr = new TThread("fReaperThread", (void(*)(void*)) &fReaperThread, (void*) this );
        fRPPtr = new TThread("fReaperThread", 
            (TThread::VoidRtnFunc_t) &fReaperThread, (void*) this );

        fRPPtr->Run();
        return 0;
    }
    return 1;
}
Int_t KScheduler::stopReaperThread(){
    if(fRPPtr){
        TThread::Delete(fRPPtr);
        delete fRPPtr;
        fRPPtr = 0;
        return 0;
    }
    return 1;
}

// takes mem for argptr, needs to be freed by the Kthread
// TODO FIX INPUT ARG TO THIS NO SLOT ID NEEDED
Int_t KScheduler::startWorkerThread(unsigned threadId){
    // check threadstatus
    std::cout << "Booting fWorkerThread:" << threadId << std::endl;
    std::string threadnm = "workerThread" + std::to_string(threadId);
    const char* formatTnm = threadnm.c_str();
    TThread* thisThread;
    workerArg* wArgPtr = new workerArg;
    wArgPtr->kschdPtr = this;
    wArgPtr->threadId = threadId;

    Printf("thread %u Gets Ptr:%p", threadId, wArgPtr);
    if(!workThreadArr[threadId]){
//        workThreadArr[threadId] = new TThread(formatTnm, (void(*)(void*)) &fWorkerThread, 
//            (void*) wArgPtr );
        workThreadArr[threadId] = new TThread(formatTnm, 
            (TThread::VoidRtnFunc_t) &fWorkerThread, (void*) wArgPtr);
        thisThread = workThreadArr[threadId];
        assert(thisThread != 0);
        thisThread->Run();
        return 0;
    }
    return 1;
}

// wrapper
Int_t KScheduler::startWorkerThreads(){
    Int_t ret;
    //TThread* thisThread;
    for(int i = 0; i < n_threads; i++){
        ret = startWorkerThread(i);
        // delay in print here is enough for the wArg to work out from run...
        // TODO NEED REAL SYNC FOR THAT (sig after set warg?)
        // THIS IS REALLY BAD TODO
        std::cout << "started thread:" << i << std::endl;
        assert(ret == 0);
    }
    return 0;
}

// TODO this is a little more confusing with worker ids
/*
Int_t KScheduler::stopWorkerThread(){
    if(workThreadArr[]){
        TThread::Delete(fRPPtr);
        delete fRPPtr;
        fRPPtr = 0;
        return 0;
    }
    return 1;

}
*/

void KScheduler::DoTrackFitting(TClonesArray* tracklets, SRecEvent* srec)
{
  int n_trk = tracklets->GetEntries();
  int n_trk_fit = 0;
  for (int i_trk = 0; i_trk < n_trk; i_trk++) {
    Tracklet* tracklet = (Tracklet*)tracklets->At(i_trk);

    KalmanFitter*       kfitter = GetKFitter();
    SQGenFit::GFFitter* gfitter = GetGFitter();
    SRecTrack recTrack;
    bool fitOK = false;
    if      (kfitter) fitOK = fitTrackCand(*tracklet, kfitter, recTrack);
    else if (gfitter) fitOK = fitTrackCand(*tracklet, gfitter, recTrack);
    
    if (!fitOK) {
      recTrack = tracklet->getSRecTrack(kfitter != 0);
      recTrack.setKalmanStatus(-1);
    } else {
      n_trk_fit++;
    }
    srec->insertTrack(recTrack);
  }
  if ((Verbose() > 0 && n_trk > 0) || Verbose() > 1) {
    int run_id = srec->getRunID();
    int spill_id = srec->getSpillID();
    int event_id = srec->getEventID();
    cout << "DoTrackFitting: run " << run_id << " spill " << spill_id << " event " << event_id << " | " << n_trk_fit << " / " << n_trk << endl;
  }
}

/**
 * Copied from SQReco...
 */
bool KScheduler::fitTrackCand(Tracklet& tracklet, KalmanFitter* fitter, SRecTrack& strack)
{
  KalmanTrack kmtrk;
  kmtrk.setTracklet(tracklet);

  if(kmtrk.getNodeList().empty()) 
  {
    LogDebug("kmtrk nodelist empty");
    return false;
  }

  if(fitter->processOneTrack(kmtrk) == 0)
  {
    LogDebug("kFitter failed to converge");
    return false;
  }

  fitter->updateTrack(kmtrk);//update after fitting

  if(!kmtrk.isValid()) 
  {
    LogDebug("kmtrk quality cut failed");
    return false;
  }

  //SRecTrack strack = kmtrk.getSRecTrack();
  strack = kmtrk.getSRecTrack();

  //Set trigger road ID
  TriggerRoad road(tracklet);
  strack.setTriggerRoad(road.getRoadID());

  //Set prop tube slopes
  strack.setNHitsInPT(tracklet.seg_x.getNHits(), tracklet.seg_y.getNHits());
  strack.setPTSlope(tracklet.seg_x.a, tracklet.seg_y.a);
  strack.setKalmanStatus(1);
  //fillRecTrack(strack);
  return true;
}

/**
 * Copied from SQReco...
 */
bool KScheduler::fitTrackCand(Tracklet& tracklet, SQGenFit::GFFitter* fitter, SRecTrack& strack)
{
  SQGenFit::GFTrack gftrk;
  gftrk.setTracklet(tracklet);

  int fitOK = fitter->processTrack(gftrk);
  if(fitOK != 0)
  {
    LogDebug("gFitter failed to converge.");
    return false;
  }

  //if(Verbosity() > Fun4AllBase::VERBOSITY_A_LOT)
  //{
  //  gftrk.postFitUpdate();
  //  gftrk.print(2);
  //}

  //TODO: A gtrack quality cut?

  //SRecTrack strack = gftrk.getSRecTrack();
  strack = gftrk.getSRecTrack();

  //Set trigger road ID
  TriggerRoad road(tracklet);
  strack.setTriggerRoad(road.getRoadID());

  //Set prop tube slopes
  strack.setNHitsInPT(tracklet.seg_x.getNHits(), tracklet.seg_y.getNHits());
  strack.setPTSlope(tracklet.seg_x.a, tracklet.seg_y.a);

  //fillRecTrack(strack);
  return true;
}
