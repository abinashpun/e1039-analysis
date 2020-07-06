//Author: Abinash Pun (NMSU)
//Plotting macro for track QA

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

 TTree *Reco_eval = (TTree*) fup->Get("Reco");

 //tree info
 float gx_st1, gy_st1, gpx_st1, gpy_st1, gpz_st1,gnhits, gndc, ntruthhits;
 float x_st1, y_st1, z_st1, px_st1, py_st1, pz_st1;
 float rec_drift_st1, sq_drift_st1, chisq_st1, prob_st1; 

 int krecstat,pid, n_tracks, nhits;


//Setting the branch address 
 Reco_eval->SetBranchAddress("pid", &pid);
 Reco_eval->SetBranchAddress("n_tracks", &n_tracks);
 Reco_eval->SetBranchAddress("chisq_st1", &chisq_st1); 
 Reco_eval->SetBranchAddress("prob_st1", &prob_st1);

 Reco_eval->SetBranchAddress("gx_st1", &gx_st1);
 Reco_eval->SetBranchAddress("gy_st1", &gy_st1);
 
 Reco_eval->SetBranchAddress("gpx_st1", &gpx_st1);
 Reco_eval->SetBranchAddress("gpy_st1", &gpy_st1);
 Reco_eval->SetBranchAddress("gpz_st1", &gpz_st1);

 Reco_eval->SetBranchAddress("x_st1", &x_st1);
 Reco_eval->SetBranchAddress("y_st1", &y_st1);

 Reco_eval->SetBranchAddress("px_st1", &px_st1);
 Reco_eval->SetBranchAddress("py_st1", &py_st1);
 Reco_eval->SetBranchAddress("pz_st1", &pz_st1);

 Reco_eval->SetBranchAddress("rec_drift_st1", &rec_drift_st1);
 Reco_eval->SetBranchAddress("sq_drift_st1", &sq_drift_st1);

 //some local variable
 float reco_pt_st1, truth_pt_st1, reco_ptot_st1, truth_ptot_st1; 
 float reco_e_st1, truth_e_st1;
 
//Histograms initializatins
 
 TH1D *ddrift = new TH1D("ddrift","DCA_residual_st1", 200, -.2, .2); 
 TH1D *dE = new TH1D("dE","(RecoE-truthE)/truthE", 240, -.3, .3);//deltaE/E
 TH1D *dP = new TH1D("dP","(RecoP-truthP)/truthP", 240, -.3, .3);//deltaP/P
 TH1D *nTracks = new TH1D("nTracks","nTracks",5,0,5);
 TH1D *chisq = new TH1D("chisq","chisq_st1",100,0,100.);
 TH1D *prob = new TH1D("prob","prob_st1",100,0,1.);

 TH2D *ddrift_p = new TH2D("ddrift_p","ddrift vs truthP",140,5.,95.,200,-0.2,0.2);
 TH2D *dE_E = new TH2D("dE_E","deltaE vs truthE",180,0.,100.,100,-5,5.);
 TH2D *dP_P = new TH2D("dP_P","deltaP vs truthP",180,0.,100.,100,-5,5.);

 for( Int_t i=1; i<Reco_eval->GetEntries(); i++) {
   Reco_eval->GetEntry( i );   

//calculation of variables
   truth_pt_st1 = sqrt(gpx_st1*gpx_st1+gpy_st1*gpy_st1);
   reco_pt_st1 = sqrt(px_st1*px_st1+py_st1*py_st1);
   truth_e_st1 =  sqrt(gpx_st1*gpx_st1 + gpy_st1*gpy_st1 + gpz_st1*gpz_st1 + mu_mass*mu_mass);
   reco_e_st1 =  sqrt(px_st1*px_st1 + py_st1*py_st1 + pz_st1*pz_st1 + mu_mass*mu_mass);
   reco_ptot_st1 = sqrt(px_st1*px_st1 + py_st1*py_st1 + pz_st1*pz_st1);
   truth_ptot_st1 = sqrt(gpx_st1*gpx_st1 + gpy_st1*gpy_st1 + gpz_st1*gpz_st1);

  //Fill the histograms
   nTracks->Fill(n_tracks);
   if(n_tracks>0){
    chisq->Fill(chisq_st1);
    prob->Fill(prob_st1);
    ddrift->Fill(rec_drift_st1-sq_drift_st1); 
    dE->Fill((reco_e_st1-truth_e_st1)/truth_e_st1);
    dP->Fill((reco_ptot_st1-truth_ptot_st1)/truth_ptot_st1); 

    dE_E->Fill(truth_e_st1,(reco_e_st1-truth_e_st1));  
    ddrift_p->Fill(truth_ptot_st1,rec_drift_st1-sq_drift_st1);  
    dP_P->Fill(truth_ptot_st1,reco_ptot_st1-truth_ptot_st1); 
   }
  }

//fitting functions
 TF1* gausfit=new TF1("gausfit","gaus",0,1000);
 gStyle->SetOptFit();
 double maxbin, center, maxx, minx;
 

 gSystem->mkdir("plots", true);

 //Draw the histograms

 TCanvas *c1 = new TCanvas("c1","c1");
 dP->Draw();
 maxbin = dP->GetMaximumBin();
 center = dP->GetBinCenter(maxbin);
//choose the fit range carefully 
 maxx = center + 0.02;
 minx = center - 0.02;
 dP->Fit("gausfit","Q","",minx,maxx);
 //c1->SetLogy();
 c1->SaveAs("plots/dP.png");


 TCanvas *c2 = new TCanvas("c2","c2");
 maxbin = dE->GetMaximumBin();
 center = dE->GetBinCenter(maxbin);
//choose the fit range carefully 
 maxx = center + 0.02;
 minx = center - 0.02;
 dE->Draw();
 //c2->SetLogy();
 dE->Fit("gausfit","Q","",minx,maxx);
 c2->SaveAs("plots/dE.png"); 

 TCanvas *c3 = new TCanvas("c3","c3");
 ddrift->Draw();
 c3->SetLogy();
 c3->SaveAs("plots/ddrift.png");
 

 TCanvas *chi =  new TCanvas("chi","chi");
 chisq->Draw();
 chi->SetLogy();
 chi->SaveAs("plots/chisq_st1.png");

 TCanvas *pb =  new TCanvas("pb","pb");
 prob->Draw();
 pb->SetLogy();
 pb->SaveAs("plots/prob_st1.png");


 TCanvas *ntrk = new TCanvas("ntrk","ntrk");
 nTracks->Draw();
 ntrk->SaveAs("plots/nTracks.png");

//2D plots
 TCanvas *c4 = new TCanvas("c4","c4");
 dP_P->Draw("colz");
 c4->SaveAs("plots/dP_P.png");

 TCanvas *c5 = new TCanvas("c5","c5");
 dE_E->Draw("colz");
 c5->SaveAs("plots/dE_E.png");

 TCanvas *c6 = new TCanvas("c6","c6");
 ddrift_p->Draw("colz"); 
 c6->SaveAs("plots/ddrift_p.png"); 


//=======Kinematic dependence work
 do_dp_kinmetics(dP_P);

}



