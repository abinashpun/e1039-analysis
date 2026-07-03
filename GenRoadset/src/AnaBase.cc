#include <iostream>
#include <sstream>
#include <TSystem.h>
#include "AnaBase.h"
using namespace std;

AnaBase::AnaBase(const std::string label)
  : m_label      (label)
  , m_tree_name  ("default_tree_name")
  , m_branch_name("default_branch_name")
{
  ostringstream oss;
  oss << "result/" << m_label;
  m_dir_out = oss.str();
  gSystem->mkdir(oss.str().c_str(), true);
  oss << "/log.txt";
  m_ofs.open(oss.str().c_str());
}

AnaBase::~AnaBase()
{
  m_ofs.close();

  ostringstream oss;
  oss << m_dir_out << "/log.txt";
  string fname = oss.str();
  ifstream ifs(fname.c_str());
  cout << "+----------------------------------------------------------------+\n"
       << "| Messages in '" << fname << "':\n"
       << "+----------------------------------------------------------------+\n"
       << string(istreambuf_iterator<char>(ifs), istreambuf_iterator<char>())
       << "+----------------------------------------------------------------+\n"
       << endl;
  ifs.close();
}

void AnaBase::Init()
{
  ;
}

void AnaBase::Analyze()
{
  ;
}

void AnaBase::End()
{
  ;
}

void AnaBase::ReadEventsFromFileList(const char* fn_list, const int n_max)
{
  vector<string> list_fname;
  ifstream ifs(fn_list);
  string fname;
  while (ifs >> fname) list_fname.push_back(fname);
  ifs.close();
  int n_fname = list_fname.size();

  m_ofs << "N of input files = " << n_fname << endl;
  cout  << "N of input files = " << n_fname << endl;
  for (int ii = 0; ii < n_fname; ii++) {
    if (n_max > 0 && ii == n_max) break;
    cout << "." << flush;
    if ((ii+1) % 50 == 0) cout << " " << (ii+1) << endl;
    //cout << "File " << setw(3) << ii << " / " << n_fname << endl;
    ReadEvents(list_fname[ii].c_str());
  }
  cout << endl;
}

void AnaBase::ReadEvents(const char* fname)
{
  ;
}

void AnaBase::ProcessOneEvent()
{
  ;
}
