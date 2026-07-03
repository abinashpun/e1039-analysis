#ifndef _NMR_UTIL__H_
#define _NMR_UTIL__H_

namespace NMRUtil {

void ConvEventNum(const int event_num, int& YYYY, int& MM, int& DD, int& hh, int& mm, int& ss)
{
  TTimeStamp ts(event_num - TTimeStamp::GetZoneOffset());
  int YYYYMMDD = ts.GetDate();
  int   hhmmss = ts.GetTime();
  YYYY =  YYYYMMDD / 10000;
  MM   = (YYYYMMDD % 10000) / 100;
  DD   =  YYYYMMDD %   100;
  hh   =  hhmmss / 10000;
  mm   = (hhmmss % 10000) / 100;
  ss   =  hhmmss %   100;
}

std::string ConvEventNum(const int event_num)
{
  int YYYY, MM, DD, hh, mm, ss;
  ConvEventNum(event_num, YYYY, MM, DD, hh, mm, ss);
  ostringstream oss;
  oss << setfill('0') << YYYY << "-" << setw(2) << MM << "-" << setw(2) << DD << " "
      << setw(2) << hh << ":" << setw(2) << mm << ":" << setw(2) << ss;
  return oss.str();
}

};

#endif // _NMR_UTIL__H_
