#include <iomanip>
#include <TSystem.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include "KJob.h"
#include "KScheduler.h"
#include "MultiThreadReco.h"
using namespace std;

MultiThreadReco::MultiThreadReco(const std::string& name)
  : SQReco(name)
  , m_out_file_name("multi_thread_reco.root")
  , m_save_raw_evt(false)
  , m_ksc(0)
{
  ;
}

MultiThreadReco::~MultiThreadReco()
{
  if (m_ksc) delete m_ksc;
}

int MultiThreadReco::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int MultiThreadReco::InitRun(PHCompositeNode* topNode)
{
  //int ret = MakeNodes(topNode);
  //if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = InitField(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = InitGeom(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  if (_enable_KF) {
    if (_fitter_type == SQReco::LEGACY) {
      _kfitter = new KalmanFitter(_phfield, _t_geo_manager);
      _kfitter->setControlParameter(50, 0.001);
    } else {
      _gfitter = new SQGenFit::GFFitter();
      if      (_fitter_type == SQReco::KF    ) _gfitter->init(_gfield, "KalmanFitter");
      else if (_fitter_type == SQReco::KFREF ) _gfitter->init(_gfield, "KalmanFitterRefTrack");
      else if (_fitter_type == SQReco::DAF   ) _gfitter->init(_gfield, "DafSimple");
      else if (_fitter_type == SQReco::DAFREF) _gfitter->init(_gfield, "DafRef");
    }
  }

  m_ksc = new KScheduler("", m_out_file_name);
  //m_ksc->UseTrackletReco(true);
  m_ksc->SaveRawEvent(m_save_raw_evt);
  m_ksc->Init(_phfield, _t_geo_manager, _kfitter, _gfitter, _gfield);
  m_ksc->runThreads();

  return Fun4AllReturnCodes::EVENT_OK;
}

int MultiThreadReco::process_event(PHCompositeNode* topNode)
{
  bool copy = true;
  if (_input_type == SQReco::E1039) {
    _rawEvent = BuildSRawEvent();
    copy = false;
  }

  m_ksc->PushEvent(_rawEvent, copy);

  //if(_input_type == SQReco::E1039) {
  //  delete _rawEvent;
  //  _rawEvent = 0;
  //}

  return Fun4AllReturnCodes::EVENT_OK;
}

int MultiThreadReco::End(PHCompositeNode* topNode)
{
  m_ksc->PushPoison();
  m_ksc->endThreads();
  
  return Fun4AllReturnCodes::EVENT_OK;
}
