#include <cstdlib>
#include <string>
#include <iostream>
//#include <iomanip>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <phool/phool.h>
#include <phool/getClass.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIOManager.h>
#include <phool/PHNodeIterator.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQSpillMap.h>
#include <interface_main/SQHitVector.h>
#include <ktracker/SRawEvent.h>
#include <ktracker/UtilSRawEvent.h>
#include "Fun4AllSRawEventOutputManagerV2.h"
using namespace std;

Fun4AllSRawEventOutputManagerV2::Fun4AllSRawEventOutputManagerV2(const string &myname, const std::string &file_name)
  : Fun4AllOutputManager(myname)
  , m_file_name(file_name)
  , m_tree_name("save")
  , m_branch_name("rawEvent")
  , m_file(0)
  , m_tree(0)
  , m_sraw(0)
  , m_evt(0)
  , m_sp_map(0)
  , m_hit_vec(0)
  , m_trig_hit_vec(0)
{
  ;
}

Fun4AllSRawEventOutputManagerV2::~Fun4AllSRawEventOutputManagerV2()
{
  CloseFile();
  if (m_sraw) delete m_sraw;
}

int Fun4AllSRawEventOutputManagerV2::Write(PHCompositeNode *startNode)
{
  if (! m_evt) {
    m_evt          = findNode::getClass<SQEvent    >(startNode, "SQEvent");
    m_sp_map       = findNode::getClass<SQSpillMap >(startNode, "SQSpillMap");
    m_hit_vec      = findNode::getClass<SQHitVector>(startNode, "SQHitVector");
    m_trig_hit_vec = findNode::getClass<SQHitVector>(startNode, "SQTriggerHitVector");
    if (!m_evt) {
      cout << PHWHERE << "Cannot find the SQEvent data nodes.  Abort." << endl;
      exit(1);
    }
  }
  if (! m_file) OpenFile();

  UtilSRawEvent::SetEvent(m_sraw, m_evt);
  if (m_sp_map      ) UtilSRawEvent::SetSpill     (m_sraw, m_sp_map->get(m_evt->get_spill_id()));
  if (m_hit_vec     ) UtilSRawEvent::SetHit       (m_sraw, m_hit_vec);
  if (m_trig_hit_vec) UtilSRawEvent::SetTriggerHit(m_sraw, m_trig_hit_vec);
  m_tree->Fill();
  return 0;
}

void Fun4AllSRawEventOutputManagerV2::OpenFile()
{
  if (Verbosity() > 0) cout << "Fun4AllSRawEventOutputManagerV2::OpenFile(): " << m_file_name << endl;
  m_file = new TFile(m_file_name.c_str(), "RECREATE");
  if (!m_file->IsOpen()) {
    cout << PHWHERE << "Could not open " << m_file_name << ".  Abort." << endl;
    exit(1);
  }
  if (!m_sraw) m_sraw = new SRawEvent();
  m_tree = new TTree(m_tree_name.c_str(), "");
  m_tree->Branch(m_branch_name.c_str(), &m_sraw);
}

void Fun4AllSRawEventOutputManagerV2::CloseFile()
{
  if (Verbosity() > 0) cout << "Fun4AllSRawEventOutputManagerV2::CloseFile(): " << m_file_name << endl;
  if (! m_file) return;
  m_file->Write();
  m_file->Close();
  delete m_file;
  m_file = 0;
}
