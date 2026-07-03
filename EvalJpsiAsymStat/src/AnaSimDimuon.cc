#include <iomanip>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzRotation.h>
#include <interface_main/SQEvent.h>
#include <interface_main/SQMCEvent.h>
#include <interface_main/SQTrackVector.h>
#include <interface_main/SQDimuonVector.h>
#include <ktracker/SRecEvent.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <UtilAna/UtilDimuon.h>
#include "UtilAsym.h"
#include "AnaSimDimuon.h"
using namespace std;

AnaSimDimuon::AnaSimDimuon(const std::string& name)
  : SubsysReco(name)
{
  ;
}

int AnaSimDimuon::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaSimDimuon::InitRun(PHCompositeNode* topNode)
{
  mi_evt      = findNode::getClass<SQEvent       >(topNode, "SQEvent");
  mi_srec     = findNode::getClass<SRecEvent     >(topNode, "SRecEvent");
  mi_evt_true = findNode::getClass<SQMCEvent     >(topNode, "SQMCEvent");
  mi_vec_trk  = findNode::getClass<SQTrackVector >(topNode, "SQTruthTrackVector");
  mi_vec_dim  = findNode::getClass<SQDimuonVector>(topNode, "SQTruthDimuonVector");
  if (!mi_evt || !mi_srec || !mi_evt_true || !mi_vec_trk || !mi_vec_dim) {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  mo_file = new TFile("sim_tree.root", "RECREATE");
  mo_tree = new TTree("tree", "Created by AnaSimDimuon");
  //mo_tree->Branch("sqevt"   ,  mi_evt);
  mo_tree->Branch("evt"     , &mo_evt);
  //mo_tree->Branch("trk_true", &mo_trk_true);
  //mo_tree->Branch("trk_reco", &mo_trk_reco);
  mo_tree->Branch("dim_true", &mo_dim_true);
  mo_tree->Branch("dim_reco", &mo_dim_reco);

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaSimDimuon::process_event(PHCompositeNode* topNode)
{
  ///
  /// Event info
  ///
  mo_evt.proc_id = mi_evt_true->get_process_id();
  for (int ii = 0; ii < 4; ii++) {
    mo_evt.par_id [ii] = mi_evt_true->get_particle_id      (ii);
    mo_evt.par_mom[ii] = mi_evt_true->get_particle_momentum(ii);
  }
  mo_evt.weight     = mi_evt_true->get_weight();
  mo_evt.nim1       = mi_evt->get_trigger(SQEvent::NIM1);
  mo_evt.fpga1      = mi_evt->get_trigger(SQEvent::MATRIX1);
  mo_evt.trig_bits  = mi_evt->get_trigger();
  mo_evt.rec_stat   = mi_srec->getRecStatus();
  mo_evt.n_dim_true = mi_vec_dim->size();
  mo_evt.n_dim_reco = mi_srec->getNDimuons();

  int pol_sign = (mi_evt->get_event_id() % 2 == 0  ?  +1  : -1);

  m_proc_id_cnt   [mi_evt_true->get_process_id()  ]++;
  m_part_id_1_cnt [mi_evt_true->get_particle_id(0)]++;
  m_part_id_2_cnt [mi_evt_true->get_particle_id(1)]++;
  m_part_id_3_cnt [mi_evt_true->get_particle_id(2)]++;
  m_part_id_4_cnt [mi_evt_true->get_particle_id(3)]++;

  ostringstream oss;
  oss << mi_evt_true->get_particle_id(0) << " " << mi_evt_true->get_particle_id(1);
  m_part_id_12_cnt[oss.str()]++;

  ///
  /// Dimuon info
  ///
  IdMap_t id_dim_t2r;
  FindDimuonRelation(id_dim_t2r);
  mo_dim_true.clear();
  mo_dim_reco.clear();
  for (unsigned int ii = 0; ii < mi_vec_dim->size(); ii++) {
    SQDimuon* dim = mi_vec_dim->at(ii);
    DimuonData dd;
    //dd.pdg_id  = dim->get_pdg_id();
    dd.pos     = dim->get_pos();
    dd.mom     = dim->get_mom();
    dd.mom_pos = dim->get_mom_pos();
    dd.mom_neg = dim->get_mom_neg();
    //if (fabs(dd.mom.M() - 3.097) < 0.001) dd.pdg_id = 443;
    UtilDimuon::CalcVar(dim, dd.mass, dd.pT, dd.x1, dd.x2, dd.xF);
    UtilDimuon::Lab2CollinsSoper(dim, dd.costh, dd.phi_s);

    //double x1mod, x2mod;
    //double phi_s, phi_s_up;
    //UtilAsym::CalcAngle(dd.mom_pos, dd.mom_neg,       +1, x1mod, x2mod, phi_s_up);
    //UtilAsym::CalcAngle(dd.mom_pos, dd.mom_neg, pol_sign, x1mod, x2mod, phi_s);
    //if (dd.pT > 2) cout << "X " << dd.x1 << " " << x1mod << " "  << dd.x2 << " " << x2mod << " "  << dd.phi_s << " " << phi_s << endl;
    //dd.x1       = x1mod;
    //dd.x2       = x2mod;
    //dd.phi_s    = phi_s;
    //dd.phi_s_up = phi_s_up;
    mo_dim_true.push_back(dd);

    //double costh_dim, phi_dim;
    //UtilDimuon::Lab2CollinsSoper(dim, costh_dim, phi_dim);
    //cout << "X "<<dd.x1<<" "<<dd.x2<<" "<<x1mod<<" "<<x2mod<<endl;
    //cout << "P "<<phi_s<<" "<<phi_dim<<" "<<costh_dim<<endl;

    DimuonData ddr;
    if (id_dim_t2r[ii] >= 0) {
      SRecDimuon dim_reco = mi_srec->getDimuon(id_dim_t2r[ii]);
      ddr.pos     = dim_reco.vtx;
      ddr.mom     = dim_reco.p_pos + dim_reco.p_neg;
      ddr.mom_pos = dim_reco.p_pos;
      ddr.mom_neg = dim_reco.p_neg;
      ddr.x1      = dim_reco.x1;
      ddr.x2      = dim_reco.x2;
    }
    mo_dim_reco.push_back(ddr);
  }
  mo_tree->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaSimDimuon::End(PHCompositeNode* topNode)
{
  ostringstream oss;
  oss << "Process ID:\n";
  for (map<int, int>::iterator it = m_proc_id_cnt.begin(); it != m_proc_id_cnt.end(); it++) {
    oss << setw(10) << it->first << "\t" << setw(10) << it->second << "\n";
  }
  oss << "Particle-1&2 IDs:\n";
  for (map<string, int>::iterator it = m_part_id_12_cnt.begin(); it != m_part_id_12_cnt.end(); it++) {
    oss << setw(10) << it->first << "\t" << setw(10) << it->second << "\n";
  }
  oss << "Particle-1 ID:\n";
  for (map<int, int>::iterator it = m_part_id_1_cnt.begin(); it != m_part_id_1_cnt.end(); it++) {
    oss << setw(10) << it->first << "\t" << setw(10) << it->second << "\n";
  }
  oss << "Particle-2 ID:\n";
  for (map<int, int>::iterator it = m_part_id_2_cnt.begin(); it != m_part_id_2_cnt.end(); it++) {
    oss << setw(10) << it->first << "\t" << setw(10) << it->second << "\n";
  }
  oss << "Particle-3 ID:\n";
  for (map<int, int>::iterator it = m_part_id_3_cnt.begin(); it != m_part_id_3_cnt.end(); it++) {
    oss << setw(10) << it->first << "\t" << setw(10) << it->second << "\n";
  }
  oss << "Particle-4 ID:\n";
  for (map<int, int>::iterator it = m_part_id_4_cnt.begin(); it != m_part_id_4_cnt.end(); it++) {
    oss << setw(10) << it->first << "\t" << setw(10) << it->second << "\n";
  }
  cout << "\n" << oss.str() << endl;
 
  ofstream ofs("pid_count.txt");
  ofs << oss.str() << endl;
  ofs.close();

  mo_file->cd();
  mo_file->Write();
  mo_file->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}

void AnaSimDimuon::FindDimuonRelation(IdMap_t& id_map)
{
  id_map.clear();
  for (unsigned int i_true = 0; i_true < mi_vec_dim->size(); i_true++) {
    SQDimuon* dim_true = mi_vec_dim->at(i_true);
    double mass_true = dim_true->get_mom().M();

    int i_reco_best = -1;
    double mass_diff_best;
    for (int i_reco = 0; i_reco < mi_srec->getNDimuons(); i_reco++) {
      SRecDimuon dim_reco = mi_srec->getDimuon(i_reco);
      double mass_diff = fabs(dim_reco.mass - mass_true);
      if (i_reco_best < 0 || mass_diff < mass_diff_best) {
        i_reco_best   = i_reco;
        mass_diff_best = mass_diff;
      }
    }
    id_map[i_true] = i_reco_best;
  }
}
