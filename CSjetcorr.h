#ifndef CSJETCORR_H
#define CSJETCORR_H

#ifndef ROOT_TH2F
#include "TH2F.h"
#endif

#ifndef ROOT_TH1F
#include "TH1F.h"
#endif

#include <iostream>

class CSjetcorr{

	public:
		TFile* f_param_a0;
		TFile* f_param_a1;
		TFile* f_param_pol0;
		TFile* f_param_ratio;
		CSjetcorr();
		double applyCSjetcorr(int nCSGT2, double calo_jtpt, int hibin);

	private: 
        double p0_cs;
        double p1_cs;
        double corr_factor;
        double pol0_corr;
        double pt_bin;
        Double_t jt_bin_bounds[43];

	    TString cent[4];
     	TF1 *f1_par0[4];
        TF1 *f1_par1[4];
        TF1 *f1_pol0[4];
        TProfile *h_ratio[4];

};

CSjetcorr::CSjetcorr(){
    cent[4] = {"0","1","2","3"};

	TFile* f_param_a0 = TFile::Open("fitparam_a0.root","read");
	TFile* f_param_a1 = TFile::Open("fitparam_a1.root","read");
	TFile* f_param_pol0 = TFile::Open("fitparam_pol0.root","read");
	TFile* f_param_ratio = TFile::Open("fitparam_ratio.root","read");

    for (int ibin = 0; ibin < 4; ibin++){

      f1_par0[ibin] = (TF1*)f_param_a0->Get((TString)("f1_par0_cent"+cent[ibin]))->Clone((TString)("f1_par0_cent"+cent[ibin]));

      f1_par1[ibin] = (TF1*)f_param_a1->Get((TString)("f1_par1_cent"+cent[ibin]))->Clone((TString)("f1_par1_cent"+cent[ibin]));

      f1_pol0[ibin] = (TF1*)f_param_pol0->Get((TString)("f1_pt_cent"+cent[ibin]))->Clone((TString)("f1_pt_cent"+cent[ibin]));

      h_ratio[ibin] = (TProfile*)f_param_ratio->Get((TString)("h_ratio_cent"+cent[ibin]))->Clone((TString)("h_ratio_cent"+cent[ibin]));
      
    }


}

double CSjetcorr::applyCSjetcorr(int nCSGT2, double calo_jtpt, int mycbin) {

	jt_bin_bounds[43]; = {80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 260., 270., 280.,290.,300.,310.,320.,330.,340.,350.,360.,370.,380.,390.,400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,500.};

    //////apply residual correction

    p0_cs = f1_par0[mycbin]->Eval(calo_jtpt);

    p1_cs = f1_par1[mycbin]->Eval(calo_jtpt);

    corr_factor = ((p1_cs[mycbin])*nCSGT2)+p0_cs[mycbin];

    calo_jtpt = calo_jtpt*(1./corr_factor);  

    ///////pol0 correction
            
    pol0_corr = f1_pol0[mycbin]->GetParameter(0);

    calo_jtpt = calo_jtpt*(1./pol0_corr); 

    ///// corrections

    for (int pt_bin = 0; pt_bin < 42; pt_bin++){

      if (calo_jtpt > jt_bin_bounds[pt_bin] && calo_jtpt <= jt_bin_bounds[pt_bin+1]){
        
        mypt_bin = pt_bin; 
      }
    }    

    corr_calo_jtpt = calo_jtpt*((h_ratio[mycbin]->GetBinContent(mypt_bin+1))); 

    return corr_calo_jtpt; 
}
#endif