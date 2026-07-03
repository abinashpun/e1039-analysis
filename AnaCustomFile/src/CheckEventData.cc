#include <iomanip>
#include <TFile.h>
#include <TTree.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQHitVector.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <UtilAna/UtilSQHit.h>
#include "CheckEventData.h"
using namespace std;

CheckEventData::CheckEventData()
  : SubsysReco("CheckEventData")
  , m_evt     (0)
  , m_hit_vec (0)
  , m_file    (0)
  , m_tree    (0)
{
  ;
}

int CheckEventData::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int CheckEventData::InitRun(PHCompositeNode* topNode)
{
  m_evt     = findNode::getClass<SQEvent    >(topNode, "SQEvent");
  m_hit_vec = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if (!m_evt || !m_hit_vec) return Fun4AllReturnCodes::ABORTEVENT;
  return Fun4AllReturnCodes::EVENT_OK;
}

int CheckEventData::process_event(PHCompositeNode* topNode)
{
  if (Verbosity() > 0) {
    cout << m_hit_vec->size() << " | ";
    for (SQHitVector::ConstIter it = m_hit_vec->begin(); it != m_hit_vec->end(); it++) {
      SQHit* hit = *it;
      cout << " " << hit->get_detector_id() << ":" << hit->get_element_id();
      if (Verbosity() > 1) cout << "[" << hit->get_drift_distance() << "]";
    }
    cout << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CheckEventData::End(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
