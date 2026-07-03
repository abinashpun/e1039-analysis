#ifndef Fun4AllVectEventInputManager_H_
#define Fun4AllVectEventInputManager_H_
#include <fun4all/Fun4AllInputManager.h>
#include <string>
#include <vector>
#include <map>

class PHCompositeNode;
class SyncObject;

class TFile;
class TTree;
class SQRun;
class SQSpillMap;
class SQEvent;
class SQHitVector;

class Fun4AllVectEventInputManager : public Fun4AllInputManager {
public:
    Fun4AllVectEventInputManager(const std::string &name = "DUMMY", const std::string &topnodename = "TOP");
    virtual ~Fun4AllVectEventInputManager();

    int fileopen(const std::string &filename);
    int fileclose();
    int run(const int nevents = 0);
    int isOpen() { return isopen; }

    void Print(const std::string &what = "ALL") const;
    int ResetEvent();
    int PushBackEvents(const int i);
    int GetSyncObject(SyncObject **mastersync);
    int SyncIt(const SyncObject *mastersync);

    const std::string& get_branch_name() const { return _branch_name; }
    void set_branch_name(const std::string& branchName) { _branch_name = branchName; }

    const std::string& get_tree_name() const { return _tree_name; }
    void set_tree_name(const std::string& treeName) { _tree_name = treeName; }

protected:
    int OpenNextFile();
    void VectToE1039();
    int segment;
    int isopen;
    int events_total;
    int events_thisfile;
    std::string topNodeName;
    PHCompositeNode *topNode;
    SyncObject* syncobject;

    std::string _tree_name;
    std::string _branch_name;

    TFile* _fin;
    TTree* _tin;
int eventID;
int spillID;
int runID;

std::vector<int>* detectorIDs = nullptr;
std::vector<int>* elementIDs = nullptr;
std::vector<double>* driftDistances = nullptr;
std::vector<double>* tdcTimes = nullptr;
std::vector<bool>* hitsInTime = nullptr;

std::vector<int>* triggerDetectorIDs = nullptr;
std::vector<int>* triggerElementIDs = nullptr;
std::vector<double>* triggerTdcTimes = nullptr;
std::vector<int>* triggerDriftDistances = nullptr;
std::vector<bool>* triggerHitsInTime = nullptr;

int fpgaTriggers[5] = {0};
int nimTriggers[5] = {0};
int rfIntensities[33] = {0};


   SQRun*       run_header;
   SQSpillMap*  spill_map;
   SQEvent*     event_header;
   SQHitVector* hit_vec;
   SQHitVector* trig_hit_vec;
};

#endif /* __Fun4AllVectEventInputManager_H_ */
