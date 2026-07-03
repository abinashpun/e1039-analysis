#ifndef __CALIB_HIT_POS_H__
#define __CALIB_HIT_POS_H__
#include <fun4all/SubsysReco.h>
class SQHitVector;

/// SubsysReco module to update the hit position, i.e. SQHit::set_pos().
class CalibHitPos: public SubsysReco {
  SQHitVector* m_vec_hit;

 public:
  CalibHitPos(const std::string &name="CalibHitPos");
  virtual ~CalibHitPos();
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
};

#endif // __CALIB_HIT_POS_H__
