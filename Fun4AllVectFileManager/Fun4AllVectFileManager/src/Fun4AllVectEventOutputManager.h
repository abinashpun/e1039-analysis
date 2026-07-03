#ifndef _FUN4ALL_UNIVERSAL_OUTPUT_MANAGER__H_
#define _FUN4ALL_UNIVERSAL_OUTPUT_MANAGER__H_

#include <fun4all/Fun4AllOutputManager.h>
#include <string>
#include <vector>
#include <map>
#include <TStopwatch.h>

class TFile;
class TTree;
class PHCompositeNode;
class SQEvent;
class SQSpillMap;
class SQHitVector;


class Fun4AllVectEventOutputManager : public Fun4AllOutputManager {
public:
    Fun4AllVectEventOutputManager(const std::string &myname = "UNIVERSALOUT");
    virtual ~Fun4AllVectEventOutputManager();

    void SetTreeName(const std::string& name) { m_tree_name = name; }
    void SetFileName(const std::string& name) { m_file_name = name; }
    virtual int Write(PHCompositeNode* startNode);
    void ResetBranches();
protected:
    int OpenFile(PHCompositeNode* startNode);
    void CloseFile();

private:
    std::string m_tree_name;
    std::string m_file_name;
    std::string m_dir_base;
    bool m_dimuon_mode;
   
    TFile* m_file;
    TTree* m_tree;

    SQEvent* m_evt;
    SQSpillMap* m_sp_map;
    SQHitVector* m_hit_vec;
    SQHitVector* m_trig_hit_vec;

int runID;
int spillID;
int eventID;
int rfID;
int turnID;
int rfIntensities[33];
int fpgaTriggers[5] = {0};
int nimTriggers[5] = {0};

std::vector<int> detectorIDs;
std::vector<int> elementIDs;
std::vector<double> tdcTimes;
std::vector<double> driftDistances;
std::vector<bool> hitsInTime;

std::vector<int> triggerDetectorIDs;
std::vector<int> triggerElementIDs;
std::vector<double> triggerTdcTimes;
std::vector<double> triggerDriftDistances;
std::vector<bool> triggerHitsInTime;

};

#endif /* _FUN4ALL_VECT_EVENT_OUTPUT_MANAGER__H_ */

