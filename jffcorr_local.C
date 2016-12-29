#include <iostream>
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "assert.h"
#include <fstream>
#include "TMath.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include "mixing_tree.h"
#include <TCut.h>
#include <vector>
#include "TCanvas.h"

using namespace std;

#define pthat_samples 10

#define nCBins 4
#define nptBins 10
#define npfbins 21

Bool_t isPbPb = true;

int phi_nbins = 20;
int eta_nbins = 18;

int jt_nbins = 42;
int cs_nbins = 14;

int mypbin, mycbin, myptbin, myrefptbin;

Double_t jt_bin_bounds[43] = {80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280.,290.,300.,310.,320.,330.,340.,350.,360.,370.,380.,390.,400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,500.};

Double_t jtpt_bins[11] = {80.,90.,100.,110.,120.,140.,170.,220.,280.,370.,500.};

Double_t cs_bin_bounds[14] = {0., 4., 6., 8., 10., 12., 14., 16., 18., 20., 22., 24., 26., 40.};

float CBins[nCBins+1] = {0, 20, 60, 100, 200};

char saythis[500];

int pthats[11] = {0,15,30,50,80,120,170,220,280,370,500};

TString cent[4] = {"0","1","2","3"};

double calo_jtpt, calo_refpt, calo_jteta, calo_jtphi, calo_refeta, calo_refphi, closure, pt_size, genpt_size, eta, phi, dr, genpt, geneta, genphi, trk_pt, chg, ratio, corr_factor, corr_factor_it1, corr_factor_it2, pt_hat, vz, hiBin, adjust, corr_calo_jtpt, pthat_weight;

//Auxillary functions defined below
void ReadFileList(std::vector<TString> &my_file_names, TString file_of_names, bool debug=false);

void FitGauss(TH1F* hist_p, Bool_t isPbPb, Float_t& mean, Float_t& meanErr, Float_t& res, Float_t& resErr);

///// CUTS ///////
const double etacut = 1.6;
const double pTmincut = 80.;
const double trketacut = 2.4;
const double trk_ptcut = 2.;

////main loop starts here////

void jffcorr_local(){	

  ///////// histograms ////////////

  TH2F *h_cs_reco[nCBins][nptBins];
  TH2F *h_cs_ref[nCBins][nptBins];

  TH2F *h_jt_closure[nCBins];
  TH2F *h_jt_closure_q[nCBins];
  TH2F *h_jt_closure_g[nCBins];
  
  TH1D* h_gen[nCBins];
  TH2D* h_reco_ref[nCBins];
 
  for (int ibin=0;ibin<nCBins;ibin++){

  	sprintf(saythis,"h_jt_closure_cent%d",ibin);
    h_jt_closure[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,100,0.5,1.5);
    h_jt_closure[ibin]->Sumw2();

    sprintf(saythis,"h_jt_closure_q_cent%d",ibin);
    h_jt_closure_q[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,100,0.5,1.5);  
    h_jt_closure_q[ibin]->Sumw2();

    sprintf(saythis,"h_jt_closure_g_cent%d",ibin);
    h_jt_closure_g[ibin] = new TH2F(saythis, "", jt_nbins-1,jt_bin_bounds,100,0.5,1.5);
  	h_jt_closure_g[ibin]->Sumw2(); 

    sprintf(saythis,"h_gen_cent%d",ibin);
    h_gen[ibin] = new TH1D(saythis,"",420,80.,500.);
    h_gen[ibin]->Sumw2();

    sprintf(saythis,"h_reco_ref_cent%d",ibin);
    h_reco_ref[ibin] = new TH2D(saythis,"",420,80.,500.,420,80.,500.);
    h_reco_ref[ibin]->Sumw2();

  	for (int ibin2=0;ibin2<nptBins;ibin2++){

      sprintf(saythis,"h_cs_reco_cent%d_pt%d",ibin,ibin2);
      h_cs_reco[ibin][ibin2] = new TH2F(saythis,"",40,0.,40.,100,0.5,1.5);	
      h_cs_reco[ibin][ibin2]->Sumw2();

      sprintf(saythis,"h_cs_ref_cent%d_pt%d",ibin,ibin2);
      h_cs_ref[ibin][ibin2] = new TH2F(saythis,"",40,0.,40.,100,0.5,1.5); 
      h_cs_ref[ibin][ibin2]->Sumw2(); 

    }   
  }

  ///////////////// centrality reweighting ///////////////////////

  TFile *cen = TFile::Open("VertexHiNcollFits.root");

  TF1 *f_cen = (TF1*)cen->Get("xfit_hi")->Clone("f_cen");

  ///////////////// correction parameters //////////////////

  TF1 *f1_par0[4];
  TF1 *f1_par1[4];

  TF1 *f1_pol0[4];

  TF1 *f1_itpol2[4]; 

  for (ibin = 0; ibin < 4; ibin++){

    TFile *f_param_a0 = TFile::Open("fitparam_a0.root");

    f1_par0[ibin] = (TF1*)f_param_a0->Get((TString)("f1_par0_cent"+cent[ibin]))->Clone((TString)("f1_par0_cent"+cent[ibin]));
  
    TFile *f_param_a1 = TFile::Open("fitparam_a1.root");

    f1_par1[ibin] = (TF1*)f_param_a1->Get((TString)("f1_par1_cent"+cent[ibin]))->Clone((TString)("f1_par1_cent"+cent[ibin]));

    TFile *f_param_pol0 = TFile::Open("fitparam_pol0.root");

    f1_pol0[ibin] = (TF1*)f_param_pol0->Get((TString)("f1_pt_cent"+cent[ibin]))->Clone((TString)("f1_pt_cent"+cent[ibin]));
/*
    TFile *f_param_it = TFile::Open("fitparam_it.root");

    f1_itpol2[ibin] = (TF1*)f_param_it->Get((TString)("f1_itpol2_cent"+cent[ibin]))->Clone((TString)("f1_itpol2_cent"+cent[ibin]));
*/  
  }

  /////open files //////

  TFile *my_file = TFile::Open("unzippedSkim_full.root");
  
  TTree *inp_tree = (TTree*)my_file->Get("unzipMixTree");
  mixing_tree *my_primary = new mixing_tree(inp_tree);
  std::cout << "Successfully retrieved tree from input file!" << std::endl;
  Long64_t n_jets = my_primary->fChain->GetEntriesFast();
  
    
  //// Loop over all reco jets ////

  for (int jet = 0; jet < n_jets; jet++){

    if (jet%500000==0) cout<<jet<<endl;

    my_primary->fChain->GetEntry(jet);

    calo_jteta = my_primary->jteta;  
    if(fabs(calo_jteta) >= etacut) continue ;  

    calo_jtpt = my_primary->jtpt;
    if (calo_jtpt <= pTmincut) continue;
            
    calo_refpt = my_primary->refpt;
            
    int refparton_flavor = my_primary->refparton_flavor;
            
    int nCSGT2 = my_primary->nCScand;
    //cout<<nCSGT2<<endl;

    //// centrality weight

    hiBin = my_primary->hiBin;
        
    if (hiBin == 0 ){continue; }
                 
    double weight_cen = f_cen->Eval(hiBin);

    for (int cbin = 0; cbin < nCBins; cbin++){ 

      if (hiBin > CBins[cbin] && hiBin <= CBins[cbin+1]){

        mycbin = cbin; 

      }

    }

    ///// pthat weight

    pthat_weight = my_primary->weight;

    ///// pt bin

    for (int ptbin = 0; ptbin < 10; ptbin++){

      if (calo_jtpt > jtpt_bins[ptbin] && calo_jtpt <= jtpt_bins[ptbin+1]){
        
        myptbin = ptbin; 
      }
    }

    for (int refptbin = 0; refptbin < 10; refptbin++){

      if (calo_refpt > jtpt_bins[refptbin] && calo_refpt <= jtpt_bins[refptbin+1]){
        
        myrefptbin = refptbin; 
      }
    }

    //////apply residual correction

    double p0_cs = f1_par0[mycbin]->Eval(calo_jtpt);

    double p1_cs = f1_par1[mycbin]->Eval(calo_jtpt);

    corr_factor = ((p1_cs[mycbin])*nCSGT2)+p0_cs[mycbin];

    calo_jtpt = calo_jtpt*(1./corr_factor);  

    ///////pol0 correction
            
    double pol0_corr = f1_pol0[mycbin]->GetParameter(0);

    calo_jtpt = calo_jtpt*(1./pol0_corr);

    /////// closure ///////// 

    closure = calo_jtpt/calo_refpt;

    ///////  filling histos ////

    h_reco_ref[mycbin]->Fill(calo_jtpt,calo_refpt,pthat_weight*weight_cen);
              
    h_cs_reco[mycbin][myptbin]->Fill(nCSGT2,closure,pthat_weight*weight_cen);

    h_cs_ref[mycbin][myrefptbin]->Fill(nCSGT2,closure,pthat_weight*weight_cen); 
              
    h_jt_closure[mycbin]->Fill(calo_refpt,closure,pthat_weight*weight_cen);
              
              
    if(refparton_flavor == -999) continue;
               
    if (fabs(refparton_flavor) == 21)h_jt_closure_g[mycbin]->Fill(calo_refpt,closure,pthat_weight*weight_cen);
    else h_jt_closure_q[mycbin]->Fill(calo_refpt,closure,pthat_weight*weight_cen);

  }

  cout<<"Ready to write histograms"<<endl;
  TFile f1("/home/dhanush/Documents/JEC/local/histos/closure_histos_corr_ncs_pol0.root", "RECREATE");   

  for (int ibin=0;ibin<nCBins;ibin++){

    h_jt_closure[ibin]->Write();
    h_jt_closure_q[ibin]->Write();
    h_jt_closure_g[ibin]->Write();

    h_reco_ref[ibin]->Write();
    h_gen[ibin]->Write();

    for (int ibin3=0;ibin3< 10;ibin3++){

      h_cs_reco[ibin][ibin3]->Write();
      h_cs_ref[ibin][ibin3]->Write();
 
    }

  }  

  f1.Close();

  cout<<"done"<<endl;

}

