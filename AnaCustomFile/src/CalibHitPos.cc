#include <iomanip>
#include <interface_main/SQHitVector.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/recoConsts.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <geom_svc/GeomSvc.h>
#include "CalibHitPos.h"
using namespace std;

CalibHitPos::CalibHitPos(const std::string& name)
  : SubsysReco(name)
  , m_vec_hit(0)
{
  ;
}

CalibHitPos::~CalibHitPos()
{
  ;
}

int CalibHitPos::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int CalibHitPos::InitRun(PHCompositeNode* topNode)
{
  m_vec_hit = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
  if (!m_vec_hit) return Fun4AllReturnCodes::ABORTEVENT;
  return Fun4AllReturnCodes::EVENT_OK;
}

int CalibHitPos::process_event(PHCompositeNode* topNode)
{
  GeomSvc* geom = GeomSvc::instance();
  for (SQHitVector::Iter it = m_vec_hit->begin(); it != m_vec_hit->end(); it++) {
    SQHit* hit = *it;
    int det = hit->get_detector_id();
    if (!geom->isChamber(det) && !geom->isPropTube(det)) continue;

    int ele = hit->get_element_id();
    double pos = geom->getMeasurement(det, ele);
    if (Verbosity() > 0) cout << "CalibHitPos: " << det << " " << ele << ": " << hit->get_pos() << " -> " << pos << endl;
    hit->set_pos(pos);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int CalibHitPos::End(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
