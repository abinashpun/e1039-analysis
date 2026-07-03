#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <phool/phool.h>
#include <phool/getClass.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIOManager.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQSpillMap.h>
#include <interface_main/SQHitVector.h>
#include <interface_main/SQDimuonVector.h>
#include <ktracker/SRecEvent.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include "Fun4AllVectEventOutputManager.h"
using namespace std;

Fun4AllVectEventOutputManager::Fun4AllVectEventOutputManager(const std::string &myname)
  : Fun4AllOutputManager(myname),
    m_file(0),
    m_tree(0),
    m_tree_name("tree"),
    m_file_name("output.root"),
    m_evt(0),
    m_sp_map(0),
    m_hit_vec(0)
{
 ;
}

Fun4AllVectEventOutputManager::~Fun4AllVectEventOutputManager() {
    CloseFile();
}

int Fun4AllVectEventOutputManager::OpenFile(PHCompositeNode* startNode) {
std::cout << "Fun4AllVectEventOutputManager::OpenFile(): Attempting to open file: " << m_file_name << " with tree: " << m_tree_name << std::endl;
m_file = new TFile(m_file_name.c_str(), "RECREATE");

if (!m_file || m_file->IsZombie()) {
    std::cerr << "Error: Could not create file " << m_file_name << std::endl;
    exit(1);
} else {
    std::cout << "File " << m_file->GetName() << " opened successfully." << std::endl;
}

m_tree = new TTree(m_tree_name.c_str(), "Tree for storing events");
if (!m_tree) {
    std::cerr << "Error: Could not create tree " << m_tree_name << std::endl;
    exit(1);
} else {
    std::cout << "Tree " << m_tree->GetName() << " created successfully." << std::endl;
}
    
m_tree->Branch("runID", &runID, "runID/I");
m_tree->Branch("spillID", &spillID, "spillID/I");
m_tree->Branch("eventID", &eventID, "eventID/I");
m_tree->Branch("rfID", &rfID, "rfID/I");
m_tree->Branch("turnID", &turnID, "turnID/I");
m_tree->Branch("rfIntensities", rfIntensities, "rfIntensities[33]/I");
m_tree->Branch("fpgaTriggers", fpgaTriggers, "fpgaTriggers[5]/I");
m_tree->Branch("nimTriggers", nimTriggers, "nimTriggers[5]/I");

m_tree->Branch("detectorIDs", &detectorIDs);
m_tree->Branch("elementIDs", &elementIDs);
m_tree->Branch("tdcTimes", &tdcTimes);
m_tree->Branch("driftDistances", &driftDistances);
m_tree->Branch("hitsInTime", &hitsInTime);

m_tree->Branch("triggerDetectorIDs", &triggerDetectorIDs);
m_tree->Branch("triggerElementIDs", &triggerElementIDs);
m_tree->Branch("triggerTdcTimes", &triggerTdcTimes);
m_tree->Branch("triggerDriftDistances", &triggerDriftDistances);
m_tree->Branch("triggerHitsInTime", &triggerHitsInTime);


    m_evt = findNode::getClass<SQEvent>(startNode, "SQEvent");
    m_hit_vec = findNode::getClass<SQHitVector>(startNode, "SQHitVector");
    m_trig_hit_vec = findNode::getClass<SQHitVector>(startNode, "SQTriggerHitVector");

    if (!m_evt || !m_hit_vec || !m_trig_hit_vec) return Fun4AllReturnCodes::ABORTEVENT;
 return Fun4AllReturnCodes::EVENT_OK;
}
int Fun4AllVectEventOutputManager::Write(PHCompositeNode* startNode) {
	if (!m_file || !m_tree) {
		OpenFile(startNode);
	}
	ResetBranches();
	runID = m_evt->get_run_id();
	spillID = m_evt->get_spill_id();
	rfID = m_evt->get_qie_rf_id();
	eventID = m_evt->get_event_id();
	turnID = m_evt->get_qie_turn_id();

	fpgaTriggers[0] = m_evt->get_trigger(SQEvent::MATRIX1);
	fpgaTriggers[1] = m_evt->get_trigger(SQEvent::MATRIX2);
	fpgaTriggers[2] = m_evt->get_trigger(SQEvent::MATRIX3);
	fpgaTriggers[3] = m_evt->get_trigger(SQEvent::MATRIX4);
	fpgaTriggers[4] = m_evt->get_trigger(SQEvent::MATRIX5);

	nimTriggers[0] = m_evt->get_trigger(SQEvent::NIM1);
	nimTriggers[1] = m_evt->get_trigger(SQEvent::NIM2);
	nimTriggers[2] = m_evt->get_trigger(SQEvent::NIM3);
	nimTriggers[3] = m_evt->get_trigger(SQEvent::NIM4);
	nimTriggers[4] = m_evt->get_trigger(SQEvent::NIM5);


	for (int i = -16; i < 16; ++i) {
    // cout << "intensity index: i" << i+16 << endl;
    	rfIntensities[i + 16] = m_evt->get_qie_rf_intensity(i);
}

if (m_hit_vec) {
    for (int ihit = 0; ihit < m_hit_vec->size(); ++ihit) {
        SQHit* hit = m_hit_vec->at(ihit);
        detectorIDs.push_back(hit->get_detector_id());
        elementIDs.push_back(hit->get_element_id());
        tdcTimes.push_back(hit->get_tdc_time());
        driftDistances.push_back(hit->get_drift_distance());
        // cout << "get drift distance: " << hit->get_drift_distance() << endl;
        hitsInTime.push_back(hit->is_in_time());
    }
}

if (m_trig_hit_vec) {
    for (int ihit = 0; ihit < m_trig_hit_vec->size(); ++ihit) {
        SQHit* hit = m_trig_hit_vec->at(ihit);
        triggerDetectorIDs.push_back(hit->get_detector_id());
        triggerElementIDs.push_back(hit->get_element_id());
        triggerTdcTimes.push_back(hit->get_tdc_time());
        triggerDriftDistances.push_back(hit->get_drift_distance());
        triggerHitsInTime.push_back(hit->is_in_time());
    }
}

	m_tree->Fill();
	return 0;
}

void Fun4AllVectEventOutputManager::CloseFile() {
    if (!m_file) return;
    std::cout << "Fun4AllVectEventOutputManager::CloseFile(): Closing file: " << m_file_name << std::endl;
    m_file->Write();
    m_file->Close();
    delete m_file;
    m_file = nullptr;
}

void Fun4AllVectEventOutputManager::ResetBranches() {
    detectorIDs.clear();
    elementIDs.clear();
    tdcTimes.clear();
    driftDistances.clear();
    hitsInTime.clear();

    triggerDetectorIDs.clear();
    triggerElementIDs.clear();
    triggerTdcTimes.clear();
    triggerDriftDistances.clear();
    triggerHitsInTime.clear();
}

