#ifndef __ABISIMULATION_H__
#define __ABISIMULATION_H__

#include <fun4all/SubsysReco.h>
#include <vector>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TTree.h>

class PHCompositeNode;
class RawClusterContainer;
class RawCluster;
class SvtxTrackMap;
class JetMap;
class JetEvalStack;
class JetRecoEval;
class SvtxTrackEval;
class PHG4TruthInfoContainer;
class PHHepMCGenEvent;
//class RawTowerGeom;
//class RawTowerGeomContainer;
class GlobalVertex;
class PHG4HitContainer;
class SQDimuonTruthInfoContainer;
class MCDimuon;

class abisimulation : public SubsysReco
{
 public:
  abisimulation(const std::string &name = "abivertex.root");

  int Init(PHCompositeNode *);
  int InitRun(PHCompositeNode *);
  int process_event(PHCompositeNode *);
  int End(PHCompositeNode *);

 
 private:
  void ResetVars();
  void InitTree();
  int GetNodes(PHCompositeNode *topNode);
  PHG4TruthInfoContainer* _truth;
  SQDimuonTruthInfoContainer* _dimuoninfo;
  //MCDimuon* _dimuoninfo;
  TFile *file;
  TTree *truth_tree;
  std::string outfilename;

  //all truth particles
  /* float truthpx, truthpy, truthpz; */
  /* float truthp; */
  /* float truthphi; */
  /* float trutheta; */
  /* float truthpt; */
  /* float truthenergy; */
  /* int truthpid; */
  /* int numparticlesinevent; */
  /* int process_id; */
 
  //Vertex variable
  float truth_vtxx[2],truth_vtxy[2],truth_vtxz[2], truthpx[2], truthpy[2], truthpz[2] ,dimuon_CosThetaCS,dimuon_phiCS,dimuon_xs, truthe[2],dimuon_m;
  int truthpid[2];


  const float pi = 3.1415926;
};

#endif  // __ABIVERTEX_H__
