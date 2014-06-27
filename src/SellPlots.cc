#include <iostream>
#include "math.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TAxis.h" 
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
  TFile* f = new TFile("/Users/cmorgoth/Software/git/DM_Signal/combine/DMm100Vu_combine_rsq_cat4.root");
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
  
  THStack* stack1 = new THStack("stack1", "");
  h_stack[3]->SetFillColor(kPink+9);//ttbar
  h_stack[2]->SetFillColor(kViolet+9);//dy
  h_stack[0]->SetFillColor(kYellow-4);//z
  h_stack[1]->SetFillColor(kSpring+4);//w
  h_stack[4]->SetLineColor(kRed-4);//Signal
  h_stack[4]->SetLineStyle(2);//Signal
  h_stack[4]->Scale(10);
  
  TLegend* leg = new TLegend(0.72,0.69,0.895,0.895);
  
  
  leg->AddEntry(h_stack[0],"Z(#nu#bar{#nu}) + jets","f");
  leg->AddEntry(h_stack[1],"W + jets","f");
  leg->AddEntry(h_stack[2],"Z/#gamma^{*}(ll) + jets","f");
  leg->AddEntry(h_stack[3],"t #bar{t} + jets","f");
  leg->AddEntry(h_in[5],"Data","lep");
  leg->AddEntry(h_stack[4],"pp #rightarrow #chi#chi + jets","l");
  leg->AddEntry((TObject*)0,"m_{#chi} = 100 GeV","");
  leg->AddEntry((TObject*)0,"#sigma = 10 pb","");
  
  
  
  stack1->Add(h_stack[3]);
  stack1->Add(h_stack[2]);
  stack1->Add(h_stack[1]);
  stack1->Add(h_stack[0]);
  
  StackSignal(stack1, h_stack[4], h_in[5], "blah1", "blah2", "DMm100Vu_cat4", "blah4", leg);
  
  return 0;	
}
