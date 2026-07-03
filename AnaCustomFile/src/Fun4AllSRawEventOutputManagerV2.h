#ifndef __FUN4ALL_SRAWEVENT_OUTPUT_MANAGER_V2_H__
#define __FUN4ALL_SRAWEVENT_OUTPUT_MANAGER_V2_H__
#include <fun4all/Fun4AllOutputManager.h>
#include <string>
#include <vector>
class TFile;
class TTree;
class PHCompositeNode;
class SRawEvent;
class SQEvent;
class SQSpillMap;
class SQHitVector;

class Fun4AllSRawEventOutputManagerV2: public Fun4AllOutputManager
{
  std::string m_file_name;
  std::string m_tree_name;
  std::string m_branch_name;
  TFile* m_file;
  TTree* m_tree;
  SRawEvent* m_sraw;
  
  SQEvent* m_evt;
  SQSpillMap* m_sp_map;
  SQHitVector* m_hit_vec;
  SQHitVector* m_trig_hit_vec;

 public:
  Fun4AllSRawEventOutputManagerV2(const std::string &myname="SRAWEVENTOUT", const std::string &file_name="srawevent.root");
  virtual ~Fun4AllSRawEventOutputManagerV2();

  void SetTreeName  (const std::string name) { m_tree_name   = name; }
  void SetBranchName(const std::string name) { m_branch_name = name; }

  virtual int Write(PHCompositeNode *startNode);

 protected:
  void OpenFile();
  void CloseFile();
};

#endif /* __FUN4ALL_SRAWEVENT_OUTPUT_MANAGER_V2_H__ */
