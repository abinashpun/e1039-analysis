#ifndef _ANA_BASE__H_
#define _ANA_BASE__H_
#include <string>
#include <fstream>

class AnaBase {
 protected:
  std::string m_label;
  std::string m_tree_name;
  std::string m_branch_name;
  std::string m_dir_out;
  std::ofstream m_ofs;

  AnaBase(const std::string label);

 public:
  virtual ~AnaBase();

  virtual void SetInputTreeName  (const std::string name) { m_tree_name   = name; }
  virtual void SetInputBranchName(const std::string name) { m_branch_name = name; }

  virtual void Init();
  virtual void Analyze();
  virtual void End();
  virtual void ReadEventsFromFileList(const char* fn_list, const int n_max=0);
  virtual void ReadEvents(const char* fname);

 protected:
  virtual void ProcessOneEvent();
};

#endif // _ANA_BASE__H_ 
