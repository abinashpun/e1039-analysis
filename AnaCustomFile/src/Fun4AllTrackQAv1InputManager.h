#ifndef _FUN4ALL_TRACK_QAV1_INPUT_MANAGER__H_
#define _FUN4ALL_TRACK_QAV1_INPUT_MANAGER__H_
#include <string>
#include <map>
#include <fun4all/Fun4AllInputManager.h>
class PHCompositeNode;
class SyncObject;
class TFile;
class TTree;
class SQRun;
class SQEvent;
class SQHitVector;

class Fun4AllTrackQAv1InputManager : public Fun4AllInputManager
{
public:
  Fun4AllTrackQAv1InputManager(const std::string &name="DUMMY", const std::string &topnodename="TOP");
  virtual ~Fun4AllTrackQAv1InputManager();
  int fileopen(const std::string &filenam);
  int fileclose();
  int run(const int nevents=0);
  int isOpen() {return isopen;}

  void Print(const std::string &what = "ALL") const;
  int ResetEvent();
  int PushBackEvents(const int i);
  int GetSyncObject(SyncObject **mastersync);
  int SyncIt(const SyncObject *mastersync);

  const std::string& GetTreeName() const { return _tree_name; }
  void SetTreeName(const std::string& treeName) { _tree_name = treeName; }
  
protected:
  int OpenNextFile();

  int segment;
  int isopen;
  int events_total;
  int events_thisfile;
  std::string topNodeName;
  PHCompositeNode* topNode;
  SyncObject* syncobject;

  /// Input
  std::string _tree_name;
  TFile* _fin;
  TTree* _tin;
  int nHits;
  int detectorID[500];
  int elementID [500];
  float driftDistance[500]; // [nHits]

  /// Output
  SQRun*       run_header;
  SQEvent*     event_header;
  SQHitVector* hit_vec;
};

#endif // _FUN4ALL_TRACK_QAV1_INPUT_MANAGER__H_
