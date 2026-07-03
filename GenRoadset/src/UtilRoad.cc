#include <fstream>
#include <sstream>
#include <iomanip>
#include <TSystem.h>
#include "RoadMap.h"
#include "RoadList.h"
#include "RoadInfo.h"
#include "UtilRoad.h"
using namespace std;

//void UtilRoad::MakeOrderedList(const RoadMap* road_map, RoadList* road_list)
//{
//  //road_list.clear();
//  //road_list.reserve(road_map.size());
//  //cout << "N " << road_map->Size() << endl;
//  for (RoadMap::ConstIter it = road_map->Begin(); it != road_map->End(); it++) {
//    //cout << "  I " << it->first << " " << it->second << endl;
//    //int road = it->first;
//    //it->second.SetRoadID(road);
//    road_list->Add(it->second);
//  }
//  road_list->Sort();
//}

void UtilRoad::PrintList(RoadList* road_list, int i_begin, int i_end, std::ostream& os)
{
  int nn = road_list->Size();
  if (nn == 0) return;

  if      (i_begin == 0) i_begin = 1;
  else if (i_begin <  0) i_begin = nn + i_begin + 1;

  if      (i_end == 0) i_end = nn;
  else if (i_end <  0) i_end = nn + i_end + 1;

  if (i_begin <  1) i_begin =  1;
  if (i_end   <  1) i_end   =  1;
  if (i_begin > nn) i_begin = nn;
  if (i_end   > nn) i_end   = nn;

  for (int ii = i_begin; ii <= i_end; ii++) {
    RoadInfo* info = road_list->Get(ii - 1);
    os << "  " << setw(5) << ii << ": "
       << setw(6) << info->GetRoadID() << " "
       << setw(5) << info->GetCount () << "  "
       << setw(6) << info->GetWeight() << "  "
       << setw(5) << info->GetCountBG () << "  "
       << setw(6) << info->GetWeightBG() << "  "
       << setw(6) << info->GetFoM() << endl;
  }
}

void UtilRoad::WriteToFile(const RoadMap* pos_top, const RoadMap* pos_bot, const RoadMap* neg_top, const RoadMap* neg_bot, const string id)
{
  const string dir_base = "roadset";
  string dir_id = dir_base + "/" + id;
  gSystem->mkdir(dir_id.c_str(), true);

  string        name[4] = { "pos_top", "pos_bot", "neg_top", "neg_bot" };
  const RoadMap* map[4] = {  pos_top ,  pos_bot ,  neg_top ,  neg_bot  };
  for (int ii = 0; ii < 4; ii++) {
    ostringstream oss;
    oss << dir_id << "/" << name[ii] << ".txt";
    ofstream ofs(oss.str());
    for (RoadMap::ConstIter it = map[ii]->Begin(); it != map[ii]->End(); it++) {
      const RoadInfo* info = it->second;
      if (! info->GetOnOff()) continue;
      ofs << info->GetRoadID() << "\t" << info->GetCount() << "\t" << info->GetWeight() << "\n";
    }
    ofs.close();
  }
}

void UtilRoad::ReadFromFile(RoadMap* pos_top, RoadMap* pos_bot, RoadMap* neg_top, RoadMap* neg_bot, const std::string id)
{
  cout << "UtilRoad::ReadFromFile():\n";
  const string dir_base = "roadset";

  string  name[4] = { "pos_top", "pos_bot", "neg_top", "neg_bot" };
  RoadMap* map[4] = {  pos_top ,  pos_bot ,  neg_top ,  neg_bot  };
  for (int ii = 0; ii < 4; ii++) {
    ostringstream oss;
    oss << dir_base << "/" << id << "/" << name[ii] << ".txt";
    cout << "  " << oss.str() << ": ";
    ifstream ifs(oss.str());
    if (! ifs.is_open()) {
      cout << "Cannot open.  Abort." << endl;
      exit(1);
    }
    int road_id, count;
    double weight;
    while (ifs >> road_id >> count >> weight) map[ii]->Add(road_id, count, weight);
    ifs.close();
    cout << map[ii]->Size() << endl;
  }
}
