#include <iostream>
#include "math.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "hlt.hh"
#include <fstream>
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include <vector>

#include "DM_1DRatio.hh"

int main(){

  gROOT->Reset();
  const int r2B[4] = {11, 6, 6, 4};
  float c1B[] = {0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.9, 0.95, 1.0, 1.2};
  float c2B[] = {0.50, 0.575, 0.65, 0.75, 0.85, .950, 1.2};
  float c3B[] = {0.50, 0.575, 0.65, 0.75, 0.85, .950, 1.2};
  float c4B[] = {0.50, 0.60, 0.70, .950, 1.2};
  std::vector<float*> v;
  v.push_back(c1B);
  v.push_back(c2B);
  v.push_back(c3B);
  v.push_back(c4B);
  TH1F* h_in[6];
  TFile* f = new TFile("/Users/cmorgoth/Software/git/DM_Signal/combine/DMm100Vu_combine_rsq_cat3.root");
  h_in[0] = (TH1F*)f->Get("z_rsq");
  h_in[1] = (TH1F*)f->Get("w_rsq");
  h_in[2] = (TH1F*)f->Get("dy_rsq");
  h_in[3] = (TH1F*)f->Get("tt_rsq");
  h_in[4] = (TH1F*)f->Get("signal_rsq");
  h_in[5] = (TH1F*)f->Get("data_obs");
  
  int nbins = h_in[0]->GetNbinsX();
  int bin_index;
  
  switch(nbins){
  case 11:
    bin_index = 0;
    break;
  case 6:
    bin_index = 1;
    break;
  case 4:
    bin_index = 3;
    break;
  default:
    std::cout << "Invalid Binning!!!" << std::endl;
    break;
  }


  TH1F* h_stack[6];//Histograms no error for stack plot(Same order as from f->Get)
  for(int i = 0; i < 6; i++){
    TString hn(Form("h_stack_%d",i));
    h_stack[i] = new TH1F(hn, hn, nbins, v.at(bin_index));
    for(int j = 1; j <= nbins; j++){
      h_stack[i]->SetBinContent(j,h_in[i]->GetBinContent(j));
    }
  }
  
  
  
  return 0;	
}
