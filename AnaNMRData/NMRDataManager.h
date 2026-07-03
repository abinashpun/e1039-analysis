#ifndef _NMR_DATA_MANAGER__H_
#define _NMR_DATA_MANAGER__H_
#include "NMRUtil.h"
#include "NMREvent.h"
#include "NMRSignal.h"

class NMRDataManager {
  int m_verb;
  std::string m_dir_base;
  std::vector<std::string> m_list_label;
  std::vector<NMREvent> m_list_evt;
  std::map<int, NMRSignal> m_list_raw_sig;
  std::map<int, NMRSignal> m_list_poly_sig;
  double m_raw_sig_min;
  double m_raw_sig_max;
  double m_poly_sig_min;
  double m_poly_sig_max;

public:
  NMRDataManager();
  virtual ~NMRDataManager() {;}

  void Verb(const int verb) { m_verb = verb; }
  int  Verb() const { return m_verb; }
  void SetDirBase(const std::string dir_base) { m_dir_base = dir_base; }
  void PrintKeyList(std::ostream& os=std::cout);

  int GetNumEvent() { return m_list_evt.size(); }
  NMREvent* GetEvent(const int index) { return &m_list_evt[index]; }
  NMREvent* FindEvent(const int evt_num);
  void ClearEvent() { m_list_evt.clear(); }
  void ReadEventFile(const std::string label);
  void ReadRawSignalFile (const std::string label);
  void ReadPolySignalFile(const std::string label);
  void ReadRawSignalFileListed ();
  void ReadPolySignalFileListed();

  NMRSignal* GetRawSignal (const int evt_num);
  NMRSignal* GetPolySignal(const int evt_num);
  void GetRawSignalRange (double& min, double& max);
  void GetPolySignalRange(double& min, double& max);

  void GetScanParam(int& ScanSweeps, int& ScanSteps, double& RFFreq, double& RFMod);
};

NMRDataManager::NMRDataManager()
  : m_verb(0)
  , m_dir_base("/data2/e1039_data/target_data/NMR")
  , m_raw_sig_min (+1e10)
  , m_raw_sig_max (-1e10)
  , m_poly_sig_min(+1e10)
  , m_poly_sig_max(-1e10)
{
  ;
}

void NMRDataManager::PrintKeyList(std::ostream& os=std::cout)
{
  int n_key = NMREvent::GetNumKey();
  for (int i_key = 0; i_key < n_key; i_key++) {
    os << "  " << i_key << "=" << NMREvent::GetKey(i_key);
  }
  os << endl;
}

NMREvent* NMRDataManager::FindEvent(const int evt_num)
{
  int idx = NMREvent::FindKey("EventNum");
  for (int i_evt = 0; i_evt < m_list_evt.size(); i_evt++) {
    NMREvent* ne = GetEvent(i_evt);
    if (ne->GetInt(idx) == evt_num) return ne;
  }
  return 0;
}

void NMRDataManager::ReadEventFile(const std::string label)
{
  std::string fname = m_dir_base + "/" + label + ".csv";
  ifstream ifs(fname.c_str());
  if (! ifs.is_open()) {
    cerr << "!!ERROR!!  Cannot find the event file for '" << label << "': " << fname << ".  Abort." << endl;
    exit(1);
  }
  m_list_label.push_back(label);
  
  bool has_key = (NMREvent::GetNumKey() > 0);
  int i_key = 0; // Used only when has_key is true.

  if (Verb() > 0 && ! has_key) cout << "Keys:";
  std::string line;
  std::getline(ifs, line);
  istringstream iss(line);
  std::string key;
  while (std::getline(iss, key, ',')) {
    if (has_key) {
      if (NMREvent::GetKey(i_key) != key) {
        cerr << "!!ERROR!!  Key mismatch: " << NMREvent::GetKey(i_key) << " != " << key << ".  Abort." << endl;
      }
      i_key++;
    } else {
      NMREvent::AddKey(key);
    }
    if (Verb() > 0 && ! has_key) cout << "  " << key;
  }
  if (Verb() > 0 && ! has_key) cout << endl;
  
  while (std::getline(ifs, line)) {
    if (Verb() > 1) cout << "Event:";
    NMREvent ne;
    istringstream iss(line);
    ne.SetLine(line);
    for (int ii = 0; ii < NMREvent::GetNumKey(); ii++) {
      std::string value;
      std::getline(iss, value, ',');
      ne.AddValue(value);
      if (Verb() > 1) cout << "  " << value;
    }
    m_list_evt.push_back(ne);
    if (Verb() > 1) cout << endl;
  }
  ifs.close();
}

void NMRDataManager::ReadRawSignalFile(const std::string label)
{
  std::string fname = m_dir_base + "/" + label + "-RawSignal.csv";
  ifstream ifs(fname.c_str());
  if (! ifs.is_open()) {
    cerr << "!!ERROR!!  Cannot find the RawSignal file for '" << label << "': " << fname << ".  Abort." << endl;
    exit(1);
  }

  std::string line;
  while (std::getline(ifs, line)) {
    if (Verb() > 1) cout << "RawSignal:";
    istringstream iss(line);
    std::string value;
    std::getline(iss, value, ',');
    int evt_num = atoi(value.c_str());
    NMRSignal* ne = &m_list_raw_sig[evt_num];
    ne->SetEventNum(evt_num);

    while (std::getline(iss, value, ',')) {
      double val = atof(value.c_str());
      ne->AddPoint(val);
      if (val < m_raw_sig_min) m_raw_sig_min = val;
      if (val > m_raw_sig_max) m_raw_sig_max = val;
      if (Verb() > 1) cout << "  " << value;
    }
    if (Verb() > 1) cout << endl;
  }
  ifs.close();
}

void NMRDataManager::ReadPolySignalFile(const std::string label)
{
  std::string fname = m_dir_base + "/" + label + "-PolySignal.csv";
  ifstream ifs(fname.c_str());
  if (! ifs.is_open()) {
    cerr << "!!ERROR!!  Cannot find the PolySignal file for '" << label << "': " << fname << ".  Abort." << endl;
    exit(1);
  }

  std::string line;
  while (std::getline(ifs, line)) {
    if (Verb() > 1) cout << "PolySignal:";
    istringstream iss(line);
    std::string value;
    std::getline(iss, value, ',');
    int evt_num = atoi(value.c_str());
    NMRSignal* ne = &m_list_poly_sig[evt_num];
    ne->SetEventNum(evt_num);

    while (std::getline(iss, value, ',')) {
      double val = atof(value.c_str());
      ne->AddPoint(val);
      if (val < m_poly_sig_min) m_poly_sig_min = val;
      if (val > m_poly_sig_max) m_poly_sig_max = val;
      if (Verb() > 1) cout << "  " << value;
    }
    if (Verb() > 1) cout << endl;
  }
  ifs.close();
}

void NMRDataManager::ReadRawSignalFileListed()
{
  if (Verb() > 0) cout << "Reading RawSignal file(s)..." << endl;
  for (auto it = m_list_label.begin(); it != m_list_label.end(); it++) {
    ReadRawSignalFile(*it);
  }
  if (m_list_raw_sig.size() == 0) {
    cerr << "!!ERROR!!  No RawSignal record was found.  Abort." << endl;
    exit(1);
  }
}

void NMRDataManager::ReadPolySignalFileListed()
{
  if (Verb() > 0) cout << "Reading PolySignal file(s)..." << endl;
  for (auto it = m_list_label.begin(); it != m_list_label.end(); it++) {
    ReadPolySignalFile(*it);
  }
  if (m_list_poly_sig.size() == 0) {
    cerr << "!!ERROR!!  No PolySignal record was found.  Abort." << endl;
    exit(1);
  }
}

NMRSignal* NMRDataManager::GetRawSignal(const int evt_num)
{
  if (m_list_raw_sig.size() == 0) ReadRawSignalFileListed();
  if (m_list_raw_sig.find(evt_num) == m_list_raw_sig.end()) return 0;
  return &m_list_raw_sig[evt_num];
}

NMRSignal* NMRDataManager::GetPolySignal(const int evt_num)
{
  if (m_list_poly_sig.size() == 0) ReadPolySignalFileListed();
  if (m_list_poly_sig.find(evt_num) == m_list_poly_sig.end()) return 0;
  return &m_list_poly_sig[evt_num];
}

void NMRDataManager::GetRawSignalRange(double& min, double& max)
{
  if (m_list_raw_sig.size() == 0) ReadRawSignalFileListed();
  min = m_raw_sig_min;
  max = m_raw_sig_max;
}

void NMRDataManager::GetPolySignalRange(double& min, double& max)
{
  if (m_list_poly_sig.size() == 0) ReadPolySignalFileListed();
  min = m_poly_sig_min;
  max = m_poly_sig_max;
}

void NMRDataManager::GetScanParam(int& ScanSweeps, int& ScanSteps, double& RFFreq, double& RFMod)
{
  if (GetNumEvent() == 0) { 
    cerr << "!!ERROR!!  GetScanParam():  No event is available.  Abort." << endl;
    exit(1);
  }
  NMREvent* ne = GetEvent(0);
  ScanSweeps = ne->GetInt   ("ScanSweeps");
  ScanSteps  = ne->GetInt   ("ScanSteps");
  RFFreq     = ne->GetDouble("RFFreq");
  RFMod      = ne->GetDouble("RFMod");
  if (Verb() > 0) cout << "ScanSweeps = " << ScanSweeps << ", ScanSteps = " << ScanSteps << ", RFFreq = " << RFFreq << ", RFMod = " << RFMod << endl;
}

#endif // _NMR_DATA_MANAGER__H_
