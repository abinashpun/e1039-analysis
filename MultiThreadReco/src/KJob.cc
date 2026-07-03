#include <iostream>
#include <algorithm>
#include <cmath>
#include <ktracker/SRawEvent.h>
#include <ktracker/SRecEvent.h>
#include "KScheduler.h"
#include "KJob.h"

int KJob::verb = 0;

// primary instantiator...
// evPtr is pointer to copy of event that doesn't get clobbered on get entry
KJob::KJob(int jobId, SRawEvent* evPtr, KScheduler* universe, bool copy)
  : jobId(jobId)
  , universe(universe)
{
    if(Verbose() > 0)
        TThread::Printf("starting KJOB: %i\n", jobId);
    
    jobMutex = new TMutex();
    // make some big memory allocations...
    // check them... somehow....
    // need DEEP copy because event is clobbered by GetEntries
    // src - dest
    if(copy){
        evData = new SRawEvent(*evPtr);
    }
    // 1039 i produce new upstream
    else{
        evData = evPtr;
    }
    recEvData = new SRecEvent();
    assert(evData);
    if(!evData){
        std::cout << "failed to allocate memory for SRawEvent... smaller INPUT_PIPE_DEPTH?\n";
    }
    p_JobStatus=STARTED;
    jobTimer = new TStopwatch();
    isPoison = false;

    // get an available, clean TClonesArr
    universe->ktrkqFSem->Wait();
    universe->ktrkQueueMutex->Lock();
    tracklets = universe->kTrackArrQueue.front();
    universe->kTrackArrQueue.pop();
    universe->ktrkQueueMutex->UnLock();
    universe->ktrkqESem->Post();
    nTracklets = 0;

    return;
};


// poisoned jobs to kill downstream threads
KJob::KJob(bool poisoned)
{
    assert(poisoned &&"shouldn't call this if not poisoned...");
    isPoison = poisoned;
    return;
}

// failing to initialize everything in poison pill was causing crashouts
KJob::~KJob()
{
    if(isPoison)
        return;
    
    delete jobMutex;
    delete evData;
    delete recEvData;
    //p_jobStatus
    if(Verbose() > 0){
        TThread::Printf("Kjob: %i completed in time:\n",jobId);
        jobTimer->Print("u");
    }
    delete jobTimer;
    // clear and yield tracklet array
    // TODO check if "C" needed in arg for tracklets
    tracklets->Clear();
    universe->ktrkqESem->Wait();
    universe->ktrkQueueMutex->Lock();
    universe->kTrackArrQueue.push(tracklets);
    universe->ktrkQueueMutex->UnLock();
    universe->ktrkqFSem->Post();
    return;
};
