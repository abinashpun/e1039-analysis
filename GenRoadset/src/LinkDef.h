#ifdef __CINT__
#pragma link off all class;
#pragma link off all function;
#pragma link off all global;

// TreeData.(h|cc)
#pragma link C++ class SignalData+;
#pragma link C++ class BgData+;

#pragma link C++ class SubsysRecoSignal-!;
#pragma link C++ class SubsysRecoBG-!;
#pragma link C++ class E906SpillSelector-!;
#pragma link C++ class FilterByTrigger-!;

#pragma link C++ class RoadInfo-!;
#pragma link C++ class RoadList-!;
#pragma link C++ class RoadMap-!;

#pragma link C++ class AnaSignal-!;
#pragma link C++ class ReAnaSignal-!;

#pragma link C++ class AnaBG-!;
#pragma link C++ class ReAnaBG-!;

#pragma link C++ class GenRoadset-!;
#pragma link C++ class DrawRoadset-!;

#endif
