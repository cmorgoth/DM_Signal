/*
Generates the input for combine-tool limit setting code
including data card and the root file with the histograms
for the systematics, signal is scale to the correct LUMI
and includes the ISR, PDF, JES systematic, as well as the
PU re-weighting, and the HLT weight for the turn-on curve
*/
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
  
  TCanvas* ca = new TCanvas("c","c", 640, 640);
  TFile* f2 = new TFile("trigger/hlt_eff_HTMHT_0mu_BOX_0bTag_Final.root");
  TEfficiency* hlt = (TEfficiency*)f2->Get("Eff2d");

  TH1F* dy[4];
  TH1F* z[4];
  TH1F* w[4];
  TH1F* tt[4];
  TH1F* bkg[4];
  TH1F* data[4];
  //Getting MR categories plots from Prediction
  TString dys, zs, ws, tts, bkgs, datas;
  double tt_N[4];//Total contribution #
  double dy_N[4];//Total contribution #
  double z_N[4];//Total contribution #
  double w_N[4];//Total contribution #
  double data_N[4];//Total contribution #
  //TFile* in = new TFile("/Users/cmorgoth/Software/git/BkgPredictionDM/Pred_Files/MR_Cat_PredV2.root");
  //TFile* in = new TFile("/Users/cmorgoth/Software/git/BkgPredictionDM/Pred_Files/MR_Cat_PredV2_Nominal.root");
  //TFile* in = new TFile("/Users/cmorgoth/Software/git/BkgPredictionDM/Pred_Files/MR_Cat_PredV2_ISR_Off.root");
  
  TFile* in = new TFile("/Users/cmorgoth/Software/git/BkgPredictionDM/PredFilesFinal/MR_Cat_PredV2_NEW_kF.root");
  
  for(int i = 0; i < 4; i++){
    dys = TString(Form("cat%d_dy_Pred",i+1));
    zs = TString(Form("cat%d_z_Pred",i+1));
    ws = TString(Form("cat%d_w_Pred",i+1));
    tts = TString(Form("cat%d_tt_Pred",i+1));
    bkgs = TString(Form("cat%d_1D_0mu_Box_Pred_sys",i+1));
    datas = TString(Form("Data_cat%d_1D_0mu_Box",i+1));
    dy[i] = (TH1F*)in->Get(dys);
    z[i] = (TH1F*)in->Get(zs);
    w[i] = (TH1F*)in->Get(ws);
    tt[i] = (TH1F*)in->Get(tts);
    bkg[i] = (TH1F*)in->Get(bkgs);
    data[i] = (TH1F*)in->Get(datas);
    dy_N[i] = dy[i]->Integral();
    z_N[i] = z[i]->Integral();
    w_N[i] = w[i]->Integral();
    tt_N[i] = tt[i]->Integral();
    data_N[i] = data[i]->Integral();
  }

  //Creating Histos for Shape Analysis
  
  TH1F* dy_up[4];
  TH1F* z_up[4];
  TH1F* w_up[4];
  TH1F* tt_up[4];
  
  TH1F* dy_down[4];
  TH1F* z_down[4];
  TH1F* w_down[4];
  TH1F* tt_down[4];
  
  TString bkgn;
  for(int i = 0; i < 4; i++){
    bkgn = TString(Form("dy_cat%d",i));
    dy_up[i] = new TH1F(bkgn, bkgn, r2B[i], v.at(i));
    bkgn = TString(Form("dy_down_cat%d",i));
    dy_down[i] = new TH1F(bkgn, bkgn, r2B[i], v.at(i));
    
    bkgn = TString(Form("z_up_cat%d",i));
    z_up[i] = new TH1F(bkgn, bkgn, r2B[i], v.at(i));
    bkgn = TString(Form("z_down_cat%d",i));
    z_down[i] = new TH1F(bkgn, bkgn, r2B[i], v.at(i));

    bkgn = TString(Form("w_up_cat%d",i));
    w_up[i] = new TH1F(bkgn, bkgn, r2B[i], v.at(i));
    bkgn = TString(Form("w_down_cat%d",i));
    w_down[i] = new TH1F(bkgn, bkgn, r2B[i], v.at(i));

    bkgn = TString(Form("tt_up_cat%d",i));
    tt_up[i] = new TH1F(bkgn, bkgn, r2B[i], v.at(i));
    bkgn = TString(Form("tt_down_cat%d",i));
    tt_down[i] = new TH1F(bkgn, bkgn, r2B[i], v.at(i));

    double min_norm = 0.0001;
    for(int j = 1; j <= dy_up[i]->GetNbinsX(); j++){
      dy_up[i]->SetBinContent(j,dy[i]->GetBinContent(j)+dy[i]->GetBinError(j));
      if(dy[i]->GetBinContent(j)-dy[i]->GetBinError(j) > 0.0){
	dy_down[i]->SetBinContent(j,dy[i]->GetBinContent(j)-dy[i]->GetBinError(j));
      }else{
	dy_down[i]->SetBinContent(j,min_norm);
      }
      
      z_up[i]->SetBinContent(j,z[i]->GetBinContent(j)+z[i]->GetBinError(j));
      if(z[i]->GetBinContent(j)-z[i]->GetBinError(j) > 0.0){
	z_down[i]->SetBinContent(j,z[i]->GetBinContent(j)-z[i]->GetBinError(j));
      }else{
	z_down[i]->SetBinContent(j,min_norm);
      }

      w_up[i]->SetBinContent(j,w[i]->GetBinContent(j)+w[i]->GetBinError(j));
      if(w[i]->GetBinContent(j)-w[i]->GetBinError(j) > 0.0){
	w_down[i]->SetBinContent(j,w[i]->GetBinContent(j)-w[i]->GetBinError(j));
      }else{
	w_down[i]->SetBinContent(j,min_norm);
      }
      
      tt_up[i]->SetBinContent(j,tt[i]->GetBinContent(j)+tt[i]->GetBinError(j));
      if(tt[i]->GetBinContent(j)-tt[i]->GetBinError(j) > 0.0){
	tt_down[i]->SetBinContent(j,tt[i]->GetBinContent(j)-tt[i]->GetBinError(j));
      }else{
	tt_down[i]->SetBinContent(j,min_norm);
      }
      
    }
    
  }

  TH1F* h_rsq[24][4];
  
  TH1F* h_rsq_ISR_up[24][4];
  TH1F* h_rsq_ISR_down[24][4];
  
  TH1F* h_rsq_JES_up[24][4];
  TH1F* h_rsq_JES_down[24][4];
  
  TH1F* h_rsq_PDF_up[24][4];
  TH1F* h_rsq_PDF_down[24][4];
  
  TH1F* s_up[24][4];
  TH1F* s_down[24][4];

  TH1F* pdf_acc[24][4];
  
  TString sn;
  for(int j = 0; j < 24; j++){
    for(int i = 0; i < 4; i++){
      sn = TString(Form("signal%d_cat%d",j,i));
      h_rsq[j][i] = new TH1F(sn, sn, r2B[i], v.at(i));
      
      h_rsq_ISR_up[j][i] = new TH1F(sn+"ISR_up", sn+"ISR_up", r2B[i], v.at(i));
      h_rsq_ISR_down[j][i] = new TH1F(sn+"ISR_down", sn+"ISR_down", r2B[i], v.at(i));
      
      h_rsq_JES_up[j][i] = new TH1F(sn+"JES_up", sn+"JES_up", r2B[i], v.at(i));
      h_rsq_JES_down[j][i] = new TH1F(sn+"JES_down", sn+"JES_down", r2B[i], v.at(i));
      
      //h_rsq_PDF_up[j][i] = new TH1F(sn+"PDF_up", sn+"PDF_up", r2B[i], v.at(i));
      //h_rsq_PDF_down[j][i] = new TH1F(sn+"PDF_down", sn+"PDF_down", r2B[i], v.at(i));
      
      s_up[j][i] = new TH1F(sn+"_up", sn+"_up", r2B[i], v.at(i));
      s_down[j][i] = new TH1F(sn+"_down", sn+"_down", r2B[i], v.at(i));
    }
  }

  //Here the program starts
  
  //std::ifstream mfile0("list_of_files_v2.list");
  std::ifstream mfile0("list_DM_BIS.list");
  std::ofstream outfile("eff_table_normal_R2_0p5_MR_200_Dphi_B_2p5_New.tex");
  
  outfile << "\\begin{table}[htdp]\n\\caption{default}\n\\begin{center}\n\\begin{tabular}{|c|c|}\n\\hline\n";
  
  std::string fname0;
  std::cout.precision(16);
  int xs_counter = 0;
  
  TFile* f_acc;
  
  if (mfile0.is_open()){
    while ( mfile0.good() ){
      mfile0 >> fname0;
      if(mfile0.eof())break;
      std::cout << fname0 << std::endl;
      int low_ = fname0.find("DMm");
      int high_ = fname0.find("_testMC_0.root") - low_;
      
      std::string dm_sample = fname0.substr(low_,high_);
      std::cout << "============ " << dm_sample << " ==================" << std::endl;
      TString pdf_acc_name = dm_sample.c_str() ;
      pdf_acc_name = "AccBIS/"+ pdf_acc_name + "_Acc.root";
      std::cout << "============ " << pdf_acc_name << " ==================" << std::endl;
      
      f_acc = new TFile(pdf_acc_name);
      for(int cat = 0; cat < 4; cat++){
	TString sname = Form("PDF_SYS_cat%d",cat+1);
	pdf_acc[xs_counter][cat] = (TH1F*)f_acc->Get(sname);
      }
      
      TFile* f = new TFile(fname0.c_str());
      TTree* eff = (TTree*)f->Get("effTree");
      TTree* out = (TTree*)f->Get("outTree");
      
      double mr[4], rsq[4], Jet_PT[20], Jet_Eta[20], Jet_Phi[20], metCorrX[4], metCorrY[4];
      double mr_up[4], rsq_up[4], Jet_PT_up[20], Jet_Eta_up[20], Jet_Phi_up[20], metCorrX_up[4], metCorrY_up[4];
      double mr_down[4], rsq_down[4], Jet_PT_down[20], Jet_Eta_down[20], Jet_Phi_down[20], metCorrX_down[4], metCorrY_down[4];
      
      double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
      double pTHem1_up, pTHem2_up, etaHem1_up, etaHem2_up, phiHem1_up, phiHem2_up;
      double pTHem1_down, pTHem2_down, etaHem1_down, etaHem2_down, phiHem1_down, phiHem2_down;
      
      int btag, box, N_Jets;
      double pu_w, ISR, ISR_up, ISR_down;
      
      double Npassed_In, Npassed_ISR, Npassed_ISR_up, Npassed_ISR_down;
      
      eff->SetBranchStatus("*", 0);
      eff->SetBranchStatus("Npassed_In", 1);
      eff->SetBranchStatus("Npassed_ISR", 1);
      eff->SetBranchStatus("Npassed_ISR_up", 1);
      eff->SetBranchStatus("Npassed_ISR_down", 1);
      
      eff->SetBranchAddress("Npassed_In", &Npassed_In);
      eff->SetBranchAddress("Npassed_ISR", &Npassed_ISR);
      eff->SetBranchAddress("Npassed_ISR_up", &Npassed_ISR_up);
      eff->SetBranchAddress("Npassed_ISR_down", &Npassed_ISR_down);
      
      
      int N_eff = eff->GetEntries();
      int Gen_Evts = 0;
      int Gen_Evts_isr = 0;
      int Gen_Evts_isrUp = 0;
      int Gen_Evts_isrDown = 0;
      
      for(int i = 0; i < N_eff; i++){
	eff->GetEntry(i);
	Gen_Evts += Npassed_In;
	Gen_Evts_isr += Npassed_ISR;
	Gen_Evts_isrUp += Npassed_ISR_up;
	Gen_Evts_isrDown += Npassed_ISR_down;
      }
      
      std::cout << "Gen_Events: " << Gen_Evts << std::endl;
      std::cout << "Gen_ISR: " << Gen_Evts_isr << std::endl;
      std::cout << "Gen_ISR_Up: " << Gen_Evts_isrUp << std::endl;
      std::cout << "Gen_ISR_Down: " << Gen_Evts_isrDown << std::endl;
      
      out->SetBranchStatus("*", 0);
      
      out->SetBranchStatus("pu_w", 1);
      out->SetBranchStatus("ISR", 1);
      out->SetBranchStatus("ISR_up", 1);
      out->SetBranchStatus("ISR_down", 1);
      
      out->SetBranchStatus("MR", 1);
      out->SetBranchStatus("MR_up", 1);
      out->SetBranchStatus("MR_down", 1);

      out->SetBranchStatus("RSQ",1);
      out->SetBranchStatus("RSQ_up",1);
      out->SetBranchStatus("RSQ_down",1);
      
      out->SetBranchStatus("nBtag", 1);
      out->SetBranchStatus("BOX_NUM",1);
      out->SetBranchStatus("N_Jets",1);
      
      out->SetBranchStatus("Jet_PT",1);
      out->SetBranchStatus("Jet_PT_up",1);
      out->SetBranchStatus("Jet_PT_down",1);
      
      out->SetBranchStatus("Jet_Phi",1);
      out->SetBranchStatus("Jet_Phi_up",1);
      out->SetBranchStatus("Jet_Phi_down",1);
      
      out->SetBranchStatus("Jet_Eta",1);
      out->SetBranchStatus("Jet_Eta_up",1);
      out->SetBranchStatus("Jet_Eta_down",1);
      
      out->SetBranchStatus("pTHem1",1);
      out->SetBranchStatus("pTHem1_up",1);
      out->SetBranchStatus("pTHem1_down",1);
      
      out->SetBranchStatus("pTHem2",1);
      out->SetBranchStatus("pTHem2_up",1);
      out->SetBranchStatus("pTHem2_down",1);
       
      out->SetBranchStatus("etaHem1",1);
      out->SetBranchStatus("etaHem1_up",1);
      out->SetBranchStatus("etaHem1_down",1);
      
      out->SetBranchStatus("etaHem2",1);
      out->SetBranchStatus("etaHem2_up",1);
      out->SetBranchStatus("etaHem2_down",1);
      
      out->SetBranchStatus("phiHem1",1);
      out->SetBranchStatus("phiHem1_up",1);
      out->SetBranchStatus("phiHem1_down",1);
      
      out->SetBranchStatus("phiHem2",1);
      out->SetBranchStatus("phiHem2_up",1);
      out->SetBranchStatus("phiHem2_down",1);
      
      out->SetBranchStatus("metCorrX",1);
      out->SetBranchStatus("metCorrX_up",1);
      out->SetBranchStatus("metCorrX_down",1);

      out->SetBranchStatus("metCorrY",1);
      out->SetBranchStatus("metCorrY_up",1);
      out->SetBranchStatus("metCorrY_down",1);

      ///////////////////////////////
      ///////////Addresses///////////
      ///////////////////////////////
      
      out->SetBranchAddress("pu_w", &pu_w);
      out->SetBranchAddress("ISR", &ISR);
      out->SetBranchAddress("ISR_up", &ISR_up);
      out->SetBranchAddress("ISR_down", &ISR_down);

      out->SetBranchAddress("MR", mr);
      out->SetBranchAddress("MR_up", mr_up);
      out->SetBranchAddress("MR_down", mr_down);
      
      out->SetBranchAddress("RSQ", rsq);
      out->SetBranchAddress("RSQ_up", rsq_up);
      out->SetBranchAddress("RSQ_down", rsq_down);
      
      out->SetBranchAddress("nBtag", &btag);
      out->SetBranchAddress("BOX_NUM", &box);
      out->SetBranchAddress("N_Jets", &N_Jets);
      
      out->SetBranchAddress("Jet_PT", Jet_PT);
      out->SetBranchAddress("Jet_PT_up", Jet_PT_up);
      out->SetBranchAddress("Jet_PT_down", Jet_PT_down);
      
      out->SetBranchAddress("Jet_Phi", Jet_Phi);
      out->SetBranchAddress("Jet_Phi_up", Jet_Phi_up);
      out->SetBranchAddress("Jet_Phi_down", Jet_Phi_down);
      
      out->SetBranchAddress("Jet_Eta", Jet_Eta);
      out->SetBranchAddress("Jet_Eta_up", Jet_Eta_up);
      out->SetBranchAddress("Jet_Eta_down", Jet_Eta_down);
      
      out->SetBranchAddress("pTHem1", &pTHem1);
      out->SetBranchAddress("pTHem1_up", &pTHem1_up);
      out->SetBranchAddress("pTHem1_down", &pTHem1_down);
      
      out->SetBranchAddress("pTHem2", &pTHem2);
      out->SetBranchAddress("pTHem2_up", &pTHem2_up);
      out->SetBranchAddress("pTHem2_down", &pTHem2_down);
      
      out->SetBranchAddress("etaHem1", &etaHem1);
      out->SetBranchAddress("etaHem1_up", &etaHem1_up);
      out->SetBranchAddress("etaHem1_down", &etaHem1_down);
      
      out->SetBranchAddress("etaHem2", &etaHem2);
      out->SetBranchAddress("etaHem2_up", &etaHem2_up);
      out->SetBranchAddress("etaHem2_down", &etaHem2_down);
      
      out->SetBranchAddress("phiHem1", &phiHem1);
      out->SetBranchAddress("phiHem1_up", &phiHem1_up);
      out->SetBranchAddress("phiHem1_down", &phiHem1_down);
      
      out->SetBranchAddress("phiHem2", &phiHem2);
      out->SetBranchAddress("phiHem2_up", &phiHem2_up);
      out->SetBranchAddress("phiHem2_down", &phiHem2_down);
      
      out->SetBranchAddress("metCorrX", metCorrX);
      out->SetBranchAddress("metCorrX_up", metCorrX_up);
      out->SetBranchAddress("metCorrX_down", metCorrX_down);
      
      out->SetBranchAddress("metCorrY", metCorrY);
      out->SetBranchAddress("metCorrY_up", metCorrY_up);
      out->SetBranchAddress("metCorrY_down", metCorrY_down);
      
      int N_out = out->GetEntries();
      //double Lumi = 18.51;//PromptReco
      double Lumi = 18.836;//Jan22Rereco
      //double scaleF = Lumi*1000./Gen_Evts;//Scale to 1 pb
      double scaleF = Lumi*1000./Gen_Evts_isr;
      double scaleF_up = Lumi*1000./Gen_Evts_isrUp;
      double scaleF_down = Lumi*1000./Gen_Evts_isrDown;
      

      double N_passed = 0.0;
      double N_passed_ISR = 0.0;
      double N_passed_ISR_up = 0.0;
      double N_passed_ISR_down = 0.0;
      double N_passed_JES_up = 0.0;
      double N_passed_JES_down = 0.0;
      
      for(int j = 0; j < N_out; j++){
	out->GetEntry(j);
	double hlt_w = HLTscale(mr[2], rsq[2], hlt);
	
	//Nominal
	TLorentzVector j1;
	TLorentzVector j2;
	j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
	j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
	double Dphi = j1.DeltaPhi(j2);
	
	if(mr[2] >= 200.0 && rsq[2] >= 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5){
	  if(mr[2] > 200.0 && mr[2] <= 300.0 ){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_rsq[xs_counter][0]->Fill(rsq[2], hlt_w*scaleF*pu_w*ISR);
	  }else if(mr[2] > 300.0 && mr[2] <= 400.0){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_rsq[xs_counter][1]->Fill(rsq[2], hlt_w*scaleF*pu_w*ISR);
	  }else if(mr[2] > 400.0 && mr[2] <= 600.0){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_rsq[xs_counter][2]->Fill(rsq[2], hlt_w*scaleF*pu_w*ISR);
	  }else if(mr[2] > 600.0 && mr[2] <= 3500.0){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_rsq[xs_counter][3]->Fill(rsq[2], hlt_w*scaleF*pu_w*ISR);
	  }
	  N_passed += hlt_w*pu_w;
	  N_passed_ISR += hlt_w*pu_w*ISR;
	}
	
	//ISR UP
	if(mr[2] >= 200.0 && rsq[2] >= 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5){
	  if(mr[2] > 200.0 && mr[2] <= 300.0 ){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_rsq_ISR_up[xs_counter][0]->Fill(rsq[2], hlt_w*scaleF_up*pu_w*ISR_up);
	  }else if(mr[2] > 300.0 && mr[2] <= 400.0){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_rsq_ISR_up[xs_counter][1]->Fill(rsq[2], hlt_w*scaleF_up*pu_w*ISR_up);
	  }else if(mr[2] > 400.0 && mr[2] <= 600.0){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_rsq_ISR_up[xs_counter][2]->Fill(rsq[2], hlt_w*scaleF_up*pu_w*ISR_up);
	  }else if(mr[2] > 600.0 && mr[2] <= 3500.0){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_rsq_ISR_up[xs_counter][3]->Fill(rsq[2], hlt_w*scaleF_up*pu_w*ISR_up);
	  }
	  N_passed_ISR_up += hlt_w*pu_w*ISR_up;
	}
	
	//ISR DOWN
	if(mr[2] >= 200.0 && rsq[2] >= 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5){
	  if(mr[2] > 200.0 && mr[2] <= 300.0 ){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_rsq_ISR_down[xs_counter][0]->Fill(rsq[2], hlt_w*scaleF_down*pu_w*ISR_down);
	  }else if(mr[2] > 300.0 && mr[2] <= 400.0){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_rsq_ISR_down[xs_counter][1]->Fill(rsq[2], hlt_w*scaleF_down*pu_w*ISR_down);
	  }else if(mr[2] > 400.0 && mr[2] <= 600.0){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_rsq_ISR_down[xs_counter][2]->Fill(rsq[2], hlt_w*scaleF_down*pu_w*ISR_down);
	  }else if(mr[2] > 600.0 && mr[2] <= 3500.0){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_rsq_ISR_down[xs_counter][3]->Fill(rsq[2], hlt_w*scaleF_down*pu_w*ISR_down);
	  }
	  N_passed_ISR_down += hlt_w*pu_w*ISR_down;
	}
	
	//JES UP
	j1.SetPtEtaPhiE(pTHem1_up, etaHem1_up, phiHem1_up, pTHem1_up*cosh(etaHem1_up));//Hemisphere
	j2.SetPtEtaPhiE(pTHem2_up, etaHem2_up, phiHem2_up, pTHem2_up*cosh(etaHem2_up));//Hemisphere
	Dphi = j1.DeltaPhi(j2);
	if(mr_up[2] >= 200.0 && rsq_up[2] >= 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5){
	  if(mr_up[2] > 200.0 && mr_up[2] <= 300.0 ){
	    if(rsq_up[2] > 1.2)rsq_up[2] = 1.1;
	    h_rsq_JES_up[xs_counter][0]->Fill(rsq_up[2], hlt_w*scaleF*pu_w*ISR);
	  }else if(mr_up[2] > 300.0 && mr_up[2] <= 400.0){
	    if(rsq_up[2] > 1.2)rsq_up[2] = 1.1;
	    h_rsq_JES_up[xs_counter][1]->Fill(rsq_up[2], hlt_w*scaleF*pu_w*ISR);
	  }else if(mr_up[2] > 400.0 && mr_up[2] <= 600.0){
	    if(rsq_up[2] > 1.2)rsq_up[2] = 1.1;
	    h_rsq_JES_up[xs_counter][2]->Fill(rsq_up[2], hlt_w*scaleF*pu_w*ISR);
	  }else if(mr_up[2] > 600.0 && mr_up[2] <= 3500.0){
	    if(rsq_up[2] > 1.2)rsq_up[2] = 1.1;
	    h_rsq_JES_up[xs_counter][3]->Fill(rsq_up[2], hlt_w*scaleF*pu_w*ISR);
	  }
	}
	
	//JES DOWN
	j1.SetPtEtaPhiE(pTHem1_down, etaHem1_down, phiHem1_down, pTHem1_down*cosh(etaHem1_down));//Hemisphere
	j2.SetPtEtaPhiE(pTHem2_down, etaHem2_down, phiHem2_down, pTHem2_down*cosh(etaHem2_down));//Hemisphere
	Dphi = j1.DeltaPhi(j2);
	if(mr_down[2] >= 200.0 && rsq_down[2] >= 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5){
	  if(mr_down[2] > 200.0 && mr_down[2] <= 300.0 ){
	    if(rsq_down[2] > 1.2)rsq_down[2] = 1.1;
	    h_rsq_JES_down[xs_counter][0]->Fill(rsq_down[2], hlt_w*scaleF*pu_w*ISR);
	  }else if(mr_down[2] > 300.0 && mr_down[2] <= 400.0){
	    if(rsq_down[2] > 1.2)rsq_down[2] = 1.1;
	    h_rsq_JES_down[xs_counter][1]->Fill(rsq_down[2], hlt_w*scaleF*pu_w*ISR);
	  }else if(mr_down[2] > 400.0 && mr_down[2] <= 600.0){
	    if(rsq_down[2] > 1.2)rsq_down[2] = 1.1;
	    h_rsq_JES_down[xs_counter][2]->Fill(rsq_down[2], hlt_w*scaleF*pu_w*ISR);
	  }else if(mr_down[2] > 600.0 && mr_down[2] <= 3500.0){
	    if(rsq_down[2] > 1.2)rsq_down[2] = 1.1;
	    h_rsq_JES_down[xs_counter][3]->Fill(rsq_down[2], hlt_w*scaleF*pu_w*ISR);
	  }
	}
	
      }
      
      double sample_eff_old = N_passed/Gen_Evts;
      double sample_eff = N_passed_ISR/Gen_Evts_isr;
      double sample_eff_up = N_passed_ISR_up/Gen_Evts_isrUp;
      double sample_eff_down = N_passed_ISR_down/Gen_Evts_isrDown;
      
      std::cout << "Sample Eff OLD: " << sample_eff_old*100 << "%" << std::endl;
      std::cout << "Sample Eff: " << sample_eff*100 << "%" << std::endl;
      std::cout << "Sample Eff Up: " << sample_eff_up*100 << "%" << std::endl;
      std::cout << "Sample Eff Down: " << sample_eff_down*100 << "%" << std::endl;
      
      outfile << dm_sample << " & " << sample_eff*100 << "\\%" << "\\" << "\\" << "\n";
      outfile << "\\hline" << std::endl;
      
      for(int i = 0; i < 4; i++){
	for(int j = 1; j <= s_up[xs_counter][i]->GetNbinsX(); j++){
	  s_up[xs_counter][i]->SetBinContent(j, h_rsq[xs_counter][i]->GetBinContent(j)+h_rsq[xs_counter][i]->GetBinError(j));
	  s_down[xs_counter][i]->SetBinContent(j, h_rsq[xs_counter][i]->GetBinContent(j)-h_rsq[xs_counter][i]->GetBinError(j));
	}
      }
      
      /////1D files
      TFile* fo;
      for(int i = 0; i < 4; i++){
	TString data_card_name1 = dm_sample.c_str();
	TString s1 = data_card_name1;
	s1 = s1+Form("_combine_rsq_cat%d.root",i+1);
	
	data_card_name1 = "CombineBIS/" + data_card_name1 + Form("_rsq_cat%d.txt",i+1);
	std::ofstream data_card_f1(data_card_name1);
	data_card_f1 << "imax 1\njmax 4\nkmax 6\n";
	data_card_f1 << "------------------------------------------------------------------------------------------\n";
	data_card_f1 << "shapes * *\t" << s1 << "\t\t$PROCESS\t$PROCESS_$SYSTEMATIC\n";
	data_card_f1 << "------------------------------------------------------------------------------------------\n";
	data_card_f1 << "Observation\t" << data_N[i] << "\n";
	data_card_f1 << "------------------------------------------------------------------------------------------\n";
	data_card_f1 << "bin\t\tb1\t\tb1\t\tb1\t\tb1\t\tb1\n";
	data_card_f1 << "process\t\tsignal_rsq\ttt_rsq\t\tdy_rsq\t\tz_rsq\t\tw_rsq\n";
	data_card_f1 << "process\t\t0\t\t1\t\t2\t\t3\t\t4\n";
	data_card_f1 << "rate\t\t"<< h_rsq[xs_counter][i]->Integral() <<"\t\t"<< tt_N[i] <<"\t\t" << dy_N[i] << "\t\t" << z_N[i] << "\t\t" << w_N[i] << "\n";
	data_card_f1 << "------------------------------------------------------------------------------------------\n";
	data_card_f1 << "lumi\tlnN\t1.026\t\t1.0\t\t1.0\t\t1.0\t\t1.0\n";
	//data_card_f1 << "alpha\tshape\t1\t\t-\t\t-\t\t-\t\t-\n";
	data_card_f1 << "Isr\tshape\t1\t\t-\t\t-\t\t-\t\t-\n";
	data_card_f1 << "Jes\tshape\t1\t\t-\t\t-\t\t-\t\t-\n";
	data_card_f1 << "Pdf\tshape\t1\t\t-\t\t-\t\t-\t\t-\n";
	data_card_f1 << "beta\tshape\t-\t\t1\t\t-\t\t-\t\t-\n";
	data_card_f1 << "gamma\tshape\t-\t\t-\t\t1\t\t1\t\t1\n";
	data_card_f1.close();

	fo = new TFile("CombineBIS/"+s1, "RECREATE");
	
	data[i]->Write("data_obs");
	
	h_rsq_PDF_up[xs_counter][i] = new TH1F( *h_rsq[xs_counter][i] );
	h_rsq_PDF_down[xs_counter][i] = new TH1F( *h_rsq[xs_counter][i] );
	for(int bin = 1; bin <= h_rsq_PDF_up[xs_counter][i]->GetNbinsX(); bin++){
	  double pdf_up = (h_rsq_PDF_up[xs_counter][i]->GetBinContent(bin))*(1.0 + pdf_acc[xs_counter][i]->GetBinContent(bin));
	  double pdf_down = (h_rsq_PDF_up[xs_counter][i]->GetBinContent(bin))*(1.0 - pdf_acc[xs_counter][i]->GetBinContent(bin));
	  h_rsq_PDF_up[xs_counter][i]->SetBinContent(bin, pdf_up);
	  h_rsq_PDF_down[xs_counter][i]->SetBinContent(bin, pdf_down);
	}
	
	
	h_rsq[xs_counter][i]->Write("signal_rsq");
	//s_up[xs_counter][i]->Write("signal_rsq_alphaUp");
	//s_down[xs_counter][i]->Write("signal_rsq_alphaDown");
	
	h_rsq_ISR_up[xs_counter][i]->Write("signal_rsq_IsrUp");
	h_rsq_ISR_down[xs_counter][i]->Write("signal_rsq_IsrDown");

	h_rsq_JES_up[xs_counter][i]->Write("signal_rsq_JesUp");
	h_rsq_JES_down[xs_counter][i]->Write("signal_rsq_JesDown");
	
	
	//h_rsq_PDF_up[xs_counter][i]->Scale(1.02);
	//h_rsq_PDF_down[xs_counter][i]->Scale(0.98);
	h_rsq_PDF_up[xs_counter][i]->Write("signal_rsq_PdfUp");
	h_rsq_PDF_down[xs_counter][i]->Write("signal_rsq_PdfDown");
	
	dy[i]->Write("dy_rsq");
	dy_up[i]->Write("dy_rsq_gammaUp");
	dy_down[i]->Write("dy_rsq_gammaDown");
	
	z[i]->Write("z_rsq");
	z_up[i]->Write("z_rsq_gammaUp");
	z_down[i]->Write("z_rsq_gammaDown");

	w[i]->Write("w_rsq");
	w_up[i]->Write("w_rsq_gammaUp");
	w_down[i]->Write("w_rsq_gammaDown");
	
	tt[i]->Write("tt_rsq");
	tt_up[i]->Write("tt_rsq_betaUp");
	tt_down[i]->Write("tt_rsq_betaDown");

	bkg[i]->Write("Pred");
	fo->Close();
	
      }
      
      xs_counter++;
      delete out;
      delete eff;
    }
    
  }else{
    std::cout << "Unable to open the file" << std::endl;
  }
  mfile0.close();
  outfile << "\\end{tabular}\n\\end{center}\n\\label{default}\n\\end{table}\n";
  outfile.close();
  
  return 0;
}
