/*Author Abinash Pun (NMSU)
* For track QA plots from output of AnaTrkQA module
*/
#include <TStyle.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <iostream>
#include <TFile.h>
#include <TF1.h>
#include "TEventList.h"
#include <TH2F.h>
#include <math.h>
#include <TAxis.h>
#include "TEventList.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <climits>
#include <cfloat>
#include <cmath>

using namespace std;

void do_dp_kinmetics(const TH2D*dpvsp){

 TCanvas *c7= new TCanvas("c7","c7");
 Double_t sigma2var[16];
 Double_t sigma2var_err[16];
 Double_t variable[16];
 Double_t mean_diff_err[16];
 Double_t mean_diff[16];
 TH1D *projection[16];
 char name[16];
 Double_t maxbin,center, suni, maxx, minx;
 TF1* gausfit=new TF1("gausfit","gaus",0,1000);
 gStyle->SetOptFit();  
  TF1* res_fit =new TF1("res_fit", "[0]/sqrt(x)+[1]+[2]/x",5.,100.);  
  c7->Divide(4,4);

  for(int i=0;i<16; i++){

    variable[i] = 5*i+10;
    Int_t lowerbin = dpvsp->GetXaxis()->FindBin(5*i+12.5);
    Int_t upperbin = dpvsp->GetXaxis()->FindBin(5*i+17.5);
    sprintf(name,"%i<P<%i",5*i+12,5*i+17);
    projection[i] = (TH1D*)dpvsp->ProjectionY(name,lowerbin,upperbin,"");
    c7->cd(i+1);
    maxbin = projection[i]->GetMaximumBin();
    center = projection[i]->GetBinCenter(maxbin);
    
    //choose the range of the gaussfit according to need
    if(i<4) suni=.5;
    else if (i<7) suni = 0.75;
    else if( i<8) suni = 1.;
    else if(i<12) suni = 1.5;
    else suni =2.0;
    
    maxx = center + suni;
    minx = center - suni;
    projection[i]->Draw();
    projection[i]->Fit("gausfit","Q",",",minx,maxx);
    sigma2var[i]=gausfit->GetParameter(2)/variable[i];
    sigma2var_err[i]=gausfit->GetParError(2)/variable[i];
    mean_diff[i] = gausfit->GetParameter(1);
    mean_diff_err[i] = gausfit->GetParError(1);    
}
 
  TGraphErrors* sigma2p_p = new TGraphErrors(16,variable,sigma2var,0,sigma2var_err);
  TCanvas *c8 = new TCanvas("c8","c8");
  sigma2p_p->Draw("AP");
  sigma2p_p->Fit("res_fit","Q","",5,90);
  c8->SaveAs("plots/mom_reso.png");

 //Mean difference vs truth plot
  TGraphErrors* dp_p = new TGraphErrors(16,variable,mean_diff,0,mean_diff_err);
  TCanvas *c9 = new TCanvas("c9","c9");
  dp_p->Draw("AP");
  dp_p->SetMarkerStyle(8);
  TLine * line = new TLine(4.,0,92.,0);
  line->Draw("same");
  line->SetLineStyle(9);
  line->SetLineWidth(2);
  dp_p->SetName("plots/deltap_p.png"); 
}

//=========

void plot_trackQA(const char* infile = "trackQA.root"){
 float mu_mass = .105658;

 TFile *fup = new TFile(infile,"read");

 TTree *Reco_eval = (TTree*) fup->Get("QA_ana");

 //tree info
 float sq_x_st1, sq_y_st1, sq_px_st1, sq_py_st1, sq_pz_st1,gnhits, gndc, ntruthhits;
 float rec_x_st1, rec_y_st1, rec_px_st1, rec_py_st1, rec_pz_st1;
 float rec_drift_st1, sq_drift_st1, chisq_st1, prob_st1, pull_q2p_st1; 

 float sq_x_st2, sq_y_st2, sq_px_st2, sq_py_st2, sq_pz_st2;
 float rec_x_st2, rec_y_st2, rec_px_st2, rec_py_st2, rec_pz_st2;
 float rec_drift_st2, sq_drift_st2, chisq_st2, prob_st2, pull_q2p_st2;  

 float sq_x_st3, sq_y_st3, sq_px_st3, sq_py_st3, sq_pz_st3;
 float rec_x_st3, rec_y_st3, rec_px_st3, rec_py_st3, rec_pz_st3;
 float rec_drift_st3, sq_drift_st3, chisq_st3, prob_st3, pull_q2p_st3; 

 int pid, n_recTracks, nhits;


//Setting the branch address 
 Reco_eval->SetBranchAddress("pid", &pid);
 Reco_eval->SetBranchAddress("n_recTracks", &n_recTracks);
 Reco_eval->SetBranchAddress("chisq_st1", &chisq_st1); 
 Reco_eval->SetBranchAddress("prob_st1", &prob_st1);

 Reco_eval->SetBranchAddress("sq_x_st1", &sq_x_st1);
 Reco_eval->SetBranchAddress("sq_y_st1", &sq_y_st1);
 
 Reco_eval->SetBranchAddress("sq_px_st1", &sq_px_st1);
 Reco_eval->SetBranchAddress("sq_py_st1", &sq_py_st1);
 Reco_eval->SetBranchAddress("sq_pz_st1", &sq_pz_st1);

 Reco_eval->SetBranchAddress("rec_x_st1", &rec_x_st1);
 Reco_eval->SetBranchAddress("rec_y_st1", &rec_y_st1);

 Reco_eval->SetBranchAddress("rec_px_st1", &rec_px_st1);
 Reco_eval->SetBranchAddress("rec_py_st1", &rec_py_st1);
 Reco_eval->SetBranchAddress("rec_pz_st1", &rec_pz_st1);

 Reco_eval->SetBranchAddress("rec_drift_st1", &rec_drift_st1);
 Reco_eval->SetBranchAddress("sq_drift_st1", &sq_drift_st1);

 Reco_eval->SetBranchAddress("pull_q2p_st1", &pull_q2p_st1);

//for st2
 Reco_eval->SetBranchAddress("sq_x_st2", &sq_x_st2);
 Reco_eval->SetBranchAddress("sq_y_st2", &sq_y_st2);
 
 Reco_eval->SetBranchAddress("sq_px_st2", &sq_px_st2);
 Reco_eval->SetBranchAddress("sq_py_st2", &sq_py_st2);
 Reco_eval->SetBranchAddress("sq_pz_st2", &sq_pz_st2);

 Reco_eval->SetBranchAddress("rec_x_st2", &rec_x_st2);
 Reco_eval->SetBranchAddress("rec_y_st2", &rec_y_st2);

 Reco_eval->SetBranchAddress("rec_px_st2", &rec_px_st2);
 Reco_eval->SetBranchAddress("rec_py_st2", &rec_py_st2);
 Reco_eval->SetBranchAddress("rec_pz_st2", &rec_pz_st2);

 Reco_eval->SetBranchAddress("rec_drift_st2", &rec_drift_st2);
 Reco_eval->SetBranchAddress("sq_drift_st2", &sq_drift_st2);

 Reco_eval->SetBranchAddress("pull_q2p_st2", &pull_q2p_st2);


//for st3
 Reco_eval->SetBranchAddress("sq_x_st3", &sq_x_st3);
 Reco_eval->SetBranchAddress("sq_y_st3", &sq_y_st3);
 
 Reco_eval->SetBranchAddress("sq_px_st3", &sq_px_st3);
 Reco_eval->SetBranchAddress("sq_py_st3", &sq_py_st3);
 Reco_eval->SetBranchAddress("sq_pz_st3", &sq_pz_st3);

 Reco_eval->SetBranchAddress("rec_x_st3", &rec_x_st3);
 Reco_eval->SetBranchAddress("rec_y_st3", &rec_y_st3);

 Reco_eval->SetBranchAddress("rec_px_st3", &rec_px_st3);
 Reco_eval->SetBranchAddress("rec_py_st3", &rec_py_st3);
 Reco_eval->SetBranchAddress("rec_pz_st3", &rec_pz_st3);

 Reco_eval->SetBranchAddress("rec_drift_st3", &rec_drift_st3);
 Reco_eval->SetBranchAddress("sq_drift_st3", &sq_drift_st3);

 Reco_eval->SetBranchAddress("pull_q2p_st3", &pull_q2p_st3);

 //some local variable
 float reco_pT_st1, truth_pT_st1, reco_ptot_st1, truth_ptot_st1; 
 float reco_e_st1, truth_e_st1;

 float reco_pT_st2, truth_pT_st2, reco_ptot_st2, truth_ptot_st2; 
 float reco_e_st2, truth_e_st2;

 float reco_pT_st3, truth_pT_st3, reco_ptot_st3, truth_ptot_st3; 
 float reco_e_st3, truth_e_st3;
 
//Histograms initializatins
 
//for st1.
 TH1D *ddrift_st1 = new TH1D("ddrift_st1","DCA_residual_st1", 200, -.01, .01); 
 TH1D *dE2E_st1 = new TH1D("dE2E_st1","(RecoE-truthE)/truthE at st1", 240, -.3, .3);//deltaE/E
 TH1D *dP2P_st1 = new TH1D("dP2P_st1","(RecoP-truthP)/truthP at st1", 240, -.3, .3);//deltaP/P

 TH1D *dpT_st1 = new TH1D("dpT_st1","(RecopT-truthpT) at st1", 240, -.3, .3);
 TH1D *dpx_st1 = new TH1D("dpx_st1","(Recopx-truthpx) at st1", 240, -.3, .3);
 TH1D *dpy_st1 = new TH1D("dpy_st1","(Recopy-truthpy) at st1", 240, -.3, .3);

 TH1D* pull_00_st1 = new TH1D("pull_00_st1","q/p pull at st1",240,-4.0,4.0);

 TH1D *nTracks = new TH1D("nTracks","nTracks",5,0,5);
 TH1D *Chisq_st1 = new TH1D("chisq","Chisq_st1",100,0,100.);
 TH1D *Prob_st1 = new TH1D("Prob_st1","Prob_st1",100,0,1.);

 TH2D *ddrift_p_st1 = new TH2D("ddrift_p","ddrift vs truthP at st1",140,5.,95.,200,-0.01,0.01);
 TH2D *dE_E_st1 = new TH2D("dE_E","deltaE vs truthE at st1",180,0.,100.,100,-5,5.);
 TH2D *dP_P_st1 = new TH2D("dP_P","deltaP vs truthP at st1",180,0.,100.,100,-5,5.);


//for st2.
 TH1D *ddrift_st2 = new TH1D("ddrift_st2","DCA_residual_st2", 200, -.01, .01); 
 TH1D *dE2E_st2 = new TH1D("dE2E_st2","(RecoE-truthE)/truthE at st2", 240, -.3, .3);//deltaE/E
 TH1D *dP2P_st2 = new TH1D("dP2P_st2","(RecoP-truthP)/truthP at st2", 240, -.3, .3);//deltaP/P

 TH1D *dpT_st2 = new TH1D("dpT_st2","(RecopT-truthpT) at st2", 240, -.3, .3);
 TH1D *dpx_st2 = new TH1D("dpx_st2","(Recopx-truthpx) at st2", 240, -.3, .3);
 TH1D *dpy_st2 = new TH1D("dpy_st2","(Recopy-truthpy) at st2", 240, -.3, .3);

 TH1D* pull_00_st2 = new TH1D("pull_00_st2","q/p pull at st2",240,-4.0,4.0);

//for st3
 TH1D *ddrift_st3 = new TH1D("ddrift_st3","DCA_residual_st3", 200, -.01, .01); 
 TH1D *dE2E_st3 = new TH1D("dE2E_st3","(RecoE-truthE)/truthE at st3", 240, -.3, .3);//deltaE/E
 TH1D *dP2P_st3 = new TH1D("dP2P_st3","(RecoP-truthP)/truthP at st3", 240, -.3, .3);//deltaP/P

 TH1D *dpT_st3 = new TH1D("dpT_st3","(RecopT-truthpT) at st3", 240, -.3, .3);
 TH1D *dpx_st3 = new TH1D("dpx_st3","(Recopx-truthpx) at st3", 240, -.3, .3);
 TH1D *dpy_st3 = new TH1D("dpy_st3","(Recopy-truthpy) at st3", 240, -.3, .3);

 TH1D* pull_00_st3 = new TH1D("pull_00_st3","q/p pull at st3",240,-4.0,4.0);

//=====loop over the entries

 for( Int_t i=0; i<Reco_eval->GetEntries(); i++) {
   Reco_eval->GetEntry( i );   

//calculation of variables st1
   truth_pT_st1 = sqrt(sq_px_st1*sq_px_st1+sq_py_st1*sq_py_st1);
   reco_pT_st1 = sqrt(rec_px_st1*rec_px_st1+rec_py_st1*rec_py_st1);
   truth_e_st1 =  sqrt(sq_px_st1*sq_px_st1 + sq_py_st1*sq_py_st1 + sq_pz_st1*sq_pz_st1 + mu_mass*mu_mass);
   reco_e_st1 =  sqrt(rec_px_st1*rec_px_st1 + rec_py_st1*rec_py_st1 + rec_pz_st1*rec_pz_st1 + mu_mass*mu_mass);
   reco_ptot_st1 = sqrt(rec_px_st1*rec_px_st1 + rec_py_st1*rec_py_st1 + rec_pz_st1*rec_pz_st1);
   truth_ptot_st1 = sqrt(sq_px_st1*sq_px_st1 + sq_py_st1*sq_py_st1 + sq_pz_st1*sq_pz_st1);

//calculation of variables st2
   truth_pT_st2 = sqrt(sq_px_st2*sq_px_st2+sq_py_st2*sq_py_st2);
   reco_pT_st2 = sqrt(rec_px_st2*rec_px_st2+rec_py_st2*rec_py_st2);
   truth_e_st2 =  sqrt(sq_px_st2*sq_px_st2 + sq_py_st2*sq_py_st2 + sq_pz_st2*sq_pz_st2 + mu_mass*mu_mass);
   reco_e_st2 =  sqrt(rec_px_st2*rec_px_st2 + rec_py_st2*rec_py_st2 + rec_pz_st2*rec_pz_st2 + mu_mass*mu_mass);
   reco_ptot_st2 = sqrt(rec_px_st2*rec_px_st2 + rec_py_st2*rec_py_st2 + rec_pz_st2*rec_pz_st2);
   truth_ptot_st2 = sqrt(sq_px_st2*sq_px_st2 + sq_py_st2*sq_py_st2 + sq_pz_st2*sq_pz_st2);


//calculation of variables st3
   truth_pT_st3 = sqrt(sq_px_st3*sq_px_st3+sq_py_st3*sq_py_st3);
   reco_pT_st3 = sqrt(rec_px_st3*rec_px_st3+rec_py_st3*rec_py_st3);
   truth_e_st3 =  sqrt(sq_px_st3*sq_px_st3 + sq_py_st3*sq_py_st3 + sq_pz_st3*sq_pz_st3 + mu_mass*mu_mass);
   reco_e_st3 =  sqrt(rec_px_st3*rec_px_st3 + rec_py_st3*rec_py_st3 + rec_pz_st3*rec_pz_st3 + mu_mass*mu_mass);
   reco_ptot_st3 = sqrt(rec_px_st3*rec_px_st3 + rec_py_st3*rec_py_st3 + rec_pz_st3*rec_pz_st3);
   truth_ptot_st3 = sqrt(sq_px_st3*sq_px_st3 + sq_py_st3*sq_py_st3 + sq_pz_st3*sq_pz_st3);

  //Fill the histograms
   nTracks->Fill(n_recTracks);
   if(n_recTracks>0){
    Chisq_st1->Fill(chisq_st1);
    Prob_st1->Fill(prob_st1);
    ddrift_st1->Fill(rec_drift_st1-sq_drift_st1); 
    dE2E_st1->Fill((reco_e_st1-truth_e_st1)/truth_e_st1);
    dP2P_st1->Fill((reco_ptot_st1-truth_ptot_st1)/truth_ptot_st1); 

    pull_00_st1->Fill(pull_q2p_st1);
   
    dpT_st1->Fill((reco_pT_st1-truth_pT_st1));
    dpx_st1->Fill((rec_px_st1-sq_px_st1));
    dpy_st1->Fill((rec_py_st1-sq_py_st1));   

    dE_E_st1->Fill(truth_e_st1,(reco_e_st1-truth_e_st1));  
    ddrift_p_st1->Fill(truth_ptot_st1,rec_drift_st1-sq_drift_st1);  
    dP_P_st1->Fill(truth_ptot_st1,reco_ptot_st1-truth_ptot_st1); 
    
    //st2. hist fill
    ddrift_st2->Fill(rec_drift_st2-sq_drift_st2); 
    dE2E_st2->Fill((reco_e_st2-truth_e_st2)/truth_e_st2);
    dP2P_st2->Fill((reco_ptot_st2-truth_ptot_st2)/truth_ptot_st2); 

    pull_00_st2->Fill(pull_q2p_st2);
   
    dpT_st2->Fill((reco_pT_st2-truth_pT_st2));
    dpx_st2->Fill((rec_px_st2-sq_px_st2));
    dpy_st2->Fill((rec_py_st2-sq_py_st2));   


   //st2. hist fill
    ddrift_st3->Fill(rec_drift_st3-sq_drift_st3); 
    dE2E_st3->Fill((reco_e_st3-truth_e_st3)/truth_e_st3);
    dP2P_st3->Fill((reco_ptot_st3-truth_ptot_st3)/truth_ptot_st3); 

    pull_00_st3->Fill(pull_q2p_st3);
   
    dpT_st3->Fill((reco_pT_st3-truth_pT_st3));
    dpx_st3->Fill((rec_px_st3-sq_px_st3));
    dpy_st3->Fill((rec_py_st3-sq_py_st3)); 

   }
  }

//fitting functions
 TF1* gausfit=new TF1("gausfit","gaus",0,1000);
 gStyle->SetOptFit();
 double maxbin, center, maxx, minx;
 

 gSystem->mkdir("plots", true);

 //Draw the histograms

 TCanvas *cdpT = new TCanvas("cdpT","cdpT");
 dpT_st1->Draw();
 cdpT->SaveAs("plots/dpT_st1.png");

 TCanvas *cdpx = new TCanvas("cdpx","cdpx");
 dpx_st1->Draw();
 cdpx->SaveAs("plots/dpx_st1.png");


 TCanvas *cdpy = new TCanvas("cdpy","cdpy");
 dpy_st1->Draw();
 cdpy->SaveAs("plots/dpy_st1.png");

 TCanvas *cpull0 = new TCanvas("cpull0","cpull0");
 pull_00_st1->Draw();
 cpull0->SaveAs("plots/pull_00_st1.png");


 TCanvas *c1 = new TCanvas("c1","c1");
 dP2P_st1->Draw();
 maxbin = dP2P_st1->GetMaximumBin();
 center = dP2P_st1->GetBinCenter(maxbin);
//choose the fit range carefully 
 maxx = center + 0.02;
 minx = center - 0.02;
 dP2P_st1->Fit("gausfit","Q","",minx,maxx);
 //c1->SetLogy();
 c1->SaveAs("plots/dP2P_st1.png");


 TCanvas *c2 = new TCanvas("c2","c2");
 maxbin = dE2E_st1->GetMaximumBin();
 center = dE2E_st1->GetBinCenter(maxbin);
//choose the fit range carefully 
 maxx = center + 0.02;
 minx = center - 0.02;
 dE2E_st1->Draw();
 //c2->SetLogy();
 dE2E_st1->Fit("gausfit","Q","",minx,maxx);
 c2->SaveAs("plots/dE2E_st1.png"); 

 TCanvas *c3 = new TCanvas("c3","c3");
 ddrift_st1->Draw();
 c3->SetLogy();
 c3->SaveAs("plots/ddrift_st1.png");
 

 TCanvas *chi =  new TCanvas("chi","chi");
 Chisq_st1->Draw();
 chi->SetLogy();
 chi->SaveAs("plots/chisq_st1.png");

 TCanvas *pb =  new TCanvas("pb","pb");
 Prob_st1->Draw();
 pb->SetLogy();
 pb->SaveAs("plots/prob_st1.png");


 TCanvas *ntrk = new TCanvas("ntrk","ntrk");
 nTracks->Draw();
 ntrk->SaveAs("plots/nTracks.png");

//2D plots
 TCanvas *c4 = new TCanvas("c4","c4");
 dP_P_st1->Draw("colz");
 c4->SaveAs("plots/dP_P_st1.png");

 TCanvas *c5 = new TCanvas("c5","c5");
 dE_E_st1->Draw("colz");
 c5->SaveAs("plots/dE_E_st1.png");

 TCanvas *c6 = new TCanvas("c6","c6");
 ddrift_p_st1->Draw("colz"); 
 c6->SaveAs("plots/ddrift_p_st1.png"); 


//====st2 plots

 TCanvas *cdpT_st2 = new TCanvas("cdpT_st2","cdpT_st2");
 dpT_st2->Draw();
 cdpT_st2->SaveAs("plots/dpT_st2.png");

 TCanvas *cdpx_st2 = new TCanvas("cdpx_st2","cdpx_st2");
 dpx_st2->Draw();
 cdpx_st2->SaveAs("plots/dpx_st2.png");


 TCanvas *cdpy_st2 = new TCanvas("cdpy_st2","cdpy_st2");
 dpy_st2->Draw();
 cdpy_st2->SaveAs("plots/dpy_st2.png");

 TCanvas *cpull0_st2 = new TCanvas("cpull0_st2","cpull0_st2");
 pull_00_st2->Draw();
 cpull0_st2->SaveAs("plots/pull_00_st2.png");


 TCanvas *c1_st2 = new TCanvas("c1_st2","c1_st2");
 dP2P_st2->Draw();
 maxbin = dP2P_st2->GetMaximumBin();
 center = dP2P_st2->GetBinCenter(maxbin);
//choose the fit range carefully 
 maxx = center + 0.02;
 minx = center - 0.02;
 dP2P_st2->Fit("gausfit","Q","",minx,maxx);
 //c1->SetLogy();
 c1_st2->SaveAs("plots/dP2P_st2.png");


 TCanvas *c2_st2 = new TCanvas("c2_st2","c2_st2");
 maxbin = dE2E_st2->GetMaximumBin();
 center = dE2E_st2->GetBinCenter(maxbin);
//choose the fit range carefully 
 maxx = center + 0.02;
 minx = center - 0.02;
 dE2E_st2->Draw();
 //c2->SetLogy();
 dE2E_st2->Fit("gausfit","Q","",minx,maxx);
 c2_st2->SaveAs("plots/dE2E_st2.png"); 

 TCanvas *c3_st2 = new TCanvas("c3_st2","c3_st2");
 ddrift_st2->Draw();
 c3_st2->SetLogy();
 c3_st2->SaveAs("plots/ddrift_st2.png");



//====st3 plots

 TCanvas *cdpT_st3 = new TCanvas("cdpT_st3","cdpT_st3");
 dpT_st3->Draw();
 cdpT_st3->SaveAs("plots/dpT_st3.png");

 TCanvas *cdpx_st3 = new TCanvas("cdpx_st3","cdpx_st3");
 dpx_st3->Draw();
 cdpx_st3->SaveAs("plots/dpx_st3.png");


 TCanvas *cdpy_st3 = new TCanvas("cdpy_st3","cdpy_st3");
 dpy_st3->Draw();
 cdpy_st3->SaveAs("plots/dpy_st3.png");

 TCanvas *cpull0_st3 = new TCanvas("cpull0_st3","cpull0_st3");
 pull_00_st3->Draw();
 cpull0_st3->SaveAs("plots/pull_00_st3.png");


 TCanvas *c1_st3 = new TCanvas("c1_st3","c1_st3");
 dP2P_st3->Draw();
 maxbin = dP2P_st3->GetMaximumBin();
 center = dP2P_st3->GetBinCenter(maxbin);
//choose the fit range carefully 
 maxx = center + 0.02;
 minx = center - 0.02;
 dP2P_st2->Fit("gausfit","Q","",minx,maxx);
 //c1->SetLogy();
 c1_st3->SaveAs("plots/dP2P_st3.png");


 TCanvas *c2_st3 = new TCanvas("c2_st3","c2_st3");
 maxbin = dE2E_st3->GetMaximumBin();
 center = dE2E_st3->GetBinCenter(maxbin);
//choose the fit range carefully 
 maxx = center + 0.02;
 minx = center - 0.02;
 dE2E_st3->Draw();
 //c2->SetLogy();
 dE2E_st3->Fit("gausfit","Q","",minx,maxx);
 c2_st3->SaveAs("plots/dE2E_st3.png"); 

 TCanvas *c3_st3 = new TCanvas("c3_st3","c3_st3");
 ddrift_st3->Draw();
 c3_st3->SetLogy();
 c3_st3->SaveAs("plots/ddrift_st3.png");

//=======Kinematic dependence work
 if(Reco_eval->GetEntries()>=50000)do_dp_kinmetics(dP_P_st1); //only for higer statistics

}


