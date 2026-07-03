#ifndef _KJOB_H_
#define _KJOB_H_
#include <list>
#include <map>
#include <queue>
#include <new>
#include <TClonesArray.h>
#include <TStopwatch.h>
#include <TThread.h>
#include <TSemaphore.h>
#include <TMutex.h>
#include <TCondition.h>
#include <GlobalConsts.h>
class SRawEvent;
class SRecEvent;
class KScheduler;

enum KJobStatus {STARTED, RUNNING, HALTED, COMPLETE};
enum KJobQueueStatus {UNINITIALIZED, INITIALIZED, OPEN, FULL};

class KJob {
    static int verb;
    
public:
    KJob(int jobId, SRawEvent* evPtr, KScheduler* universe, bool copy);
    KJob(bool poisoned);
    ~KJob();

    static void Verbose(const int a) { verb = a; }
    static int  Verbose() { return verb; }

    //TODO needs a mutex     // TMUTEX 
    //KJobStatus getJobStatus();
    int jobId;
    KScheduler* universe;
    int nTracklets = 0;
    bool isPoison;
    // TODO probably dont need a jobmutex
    TMutex* jobMutex;
    TStopwatch* jobTimer;
    SRawEvent* evData;

    // tracklets
    TClonesArray* tracklets; 

    SRecEvent* recEvData;
    KJobStatus p_JobStatus;
};

#endif 
