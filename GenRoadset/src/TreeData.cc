#include "TreeData.h"
using namespace std;

SignalData::SignalData()
  : weight(1.0)
  , mass  (0.0)
  , pT    (0.0)
  , xF    (0.0)
  , x1    (0.0)
  , x2    (0.0)
  , dz    (0.0)
  , mom   (0.0)
  , phi   (0.0)
  , theta (0.0)
  , road_pos(0)
  , road_neg(0)
{
  ;
}

BgData::BgData()
  : run       (0)
  , evt       (0)
  , fpga1     (false)
  , pot_rfp00 (0)
  , inte_rfp00(0)
  , inte_max  (0)
{
  ;
}
