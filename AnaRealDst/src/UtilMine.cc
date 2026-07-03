#include <iostream>
#include <sstream>
#include <iomanip>
#include <TSystem.h>
#include <UtilAna/UtilOnline.h>
#include "UtilMine.h"
using namespace std;

std::string UtilMine::GetDstFilePath(const int run_id, const int spill_id, const std::string dir_dst)
{
  ostringstream oss;
  if (dir_dst != "") oss << dir_dst;
  else               oss << UtilOnline::GetDstFileDir();
  oss << setfill('0') << "/run_" << setw(6) << run_id;
  oss << "/run_" << setw(6) << run_id << "_spill_" << setw(9) << spill_id << "_spin.root";
  return oss.str();
}

vector<string> UtilMine::GetListOfSpillDSTs(const int run, const string dir_dst)
{
  ostringstream oss;
  if (dir_dst != "") oss << dir_dst;
  else               oss << UtilOnline::GetDstFileDir();
  oss << setfill('0') << "/run_" << setw(6) << run;
  string dir_run = oss.str();
  cout << "DST directory = " << dir_run << endl;

  vector<string> list_dst;

  void* dirp = gSystem->OpenDirectory(dir_run.c_str());
  if (dirp == 0) return list_dst; // The directory does not exist.

  const char* name_char;
  while ((name_char = gSystem->GetDirEntry(dirp))) {
    string name = name_char;
    int length = name.length();
    if (length < 10 ||
        name.substr(0, 4) != "run_" ||
        name.substr(length-10, 10) != "_spin.root") continue;
    cout << "  DST file = " << name << endl;
    list_dst.push_back(dir_run+"/"+name);
  }
  cout << "  N of DST files = " << list_dst.size() << endl;
  gSystem->FreeDirectory(dirp);
  sort(list_dst.begin(), list_dst.end());
  return list_dst;
}
