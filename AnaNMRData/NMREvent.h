#ifndef _NMR_EVENT__H_
#define _NMR_EVENT__H_

class NMREvent {
  static std::vector<std::string> m_list_key;
  static std::map<std::string, int> m_map_key;

  std::string m_line;
  std::vector<std::string> m_list_val;
public:
  NMREvent() : m_line("") {;}
  virtual ~NMREvent() {;}

  static int GetNumKey() { return m_list_key.size(); }
  static std::string GetKey(const int idx) { return m_list_key[idx]; }
  static void AddKey(const std::string key);
  static int FindKey(const std::string key);

  std::string GetLine() { return m_line; }
  void        SetLine(const std::string line) { m_line = line; }

  std::string GetString(const int idx) { return      m_list_val[idx]; }
  int         GetInt   (const int idx) { return atoi(m_list_val[idx].c_str()); }
  double      GetDouble(const int idx) { return atof(m_list_val[idx].c_str()); }

  std::string GetString(const std::string key) { return GetString(FindKey(key)); }
  int         GetInt   (const std::string key) { return GetInt   (FindKey(key)); }
  double      GetDouble(const std::string key) { return GetDouble(FindKey(key)); }

  void AddValue(const std::string value) { m_list_val.push_back(value); }
  void ClearValue() { m_list_val.clear(); }
};

std::vector<std::string> NMREvent::m_list_key;
std::map<std::string, int> NMREvent::m_map_key;

void NMREvent::AddKey(const std::string key)
{
  m_list_key.push_back(key);
  m_map_key[key] = m_list_key.size() - 1;
}

int NMREvent::FindKey(const std::string key)
{
  if (m_map_key.find(key) == m_map_key.end()) {
    cerr << "!!ERROR!!  Cannot find '" << key << "' in the key list.  Abort." << endl;
    exit(1);
  }
  return m_map_key[key]; 
}

#endif // _NMR_EVENT__H_
