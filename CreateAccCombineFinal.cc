/*
This code obatains the acceptance 
including the PDF systematic,
Calculates the systematic error band
for the pdf acceptance
it outpus a ROOT file with a histogram containing
the error band
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
#include <cmath>

int main(){

  //gROOT->Reset();


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
  
  //TFile* in = new TFile("/Users/cmorgoth/Software/git/BkgPredictionDM/PredFilesFinal/MR_Cat_PredV2_NEW_kF.root");
  TFile* in = new TFile("/Users/cmorgoth/Software/git/BkgPredictionDM/PredFilesAN/MR_Cat_PredV2_NEW_kF.root");
  
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
  
  TH1F* h_acc[24][4];
  
  TH1F* h_acc_pdf_CTEQ[24][45][4];
  TH1F* h_acc_pdf_MRST[24][31][4];
  TH1F* h_acc_pdf_NNPDF[24][101][4];
  
  TH1F* h_acc_pdf_CTEQ_up[24][4];
  TH1F* h_acc_pdf_MRST_up[24][4];
  TH1F* h_acc_pdf_NNPDF_up[24][4];
  
  TH1F* h_acc_pdf_CTEQ_down[24][4];
  TH1F* h_acc_pdf_MRST_down[24][4];
  TH1F* h_acc_pdf_NNPDF_down[24][4];

  TH1F* h_acc_CENTRAL[24][4];
  
  TH1F* h_acc_PDF_up[24][4];
  TH1F* h_acc_PDF_down[24][4];
  
  TH1F* h_acc_ISR_up[24][4];
  TH1F* h_acc_ISR_down[24][4];
  
  TH1F* h_acc_JES_up[24][4];
  TH1F* h_acc_JES_down[24][4];
  //Here the program starts
  
  //std::ifstream mfile0("list_of_files_v2.list");
  std::ifstream mfile0("list_DM_BIS_PDF.list");
  //std::ifstream mfile0("list_DM_BugFixed.list");
  std::ofstream outfile("eff_table_normal_R2_0p5_MR_200_Dphi_B_2p5_New.tex");

  std::string fname0;
  std::cout.precision(16);
  int xs_counter = 0;
   
  if (mfile0.is_open()){
    while ( mfile0.good() ){
      mfile0 >> fname0;
      if(mfile0.eof())break;
      std::cout << fname0 << std::endl;
      int low_ = fname0.find("DMm");
      int high_ = fname0.find("_testMC_0.root") - low_;
      
      std::string dm_sample = fname0.substr(low_,high_);
   
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

      int nCTEQ66, nMRST2006NNLO, nNNPDF10100;
      double wCTEQ66[100], wMRST2006NNLO[100], wNNPDF10100[150];
      
      double N_CTEQ66[100], N_MRST2006NNLO[100], N_NNPDF10100[150];
      double N_CTEQ66_isr[100], N_MRST2006NNLO_isr[100], N_NNPDF10100_isr[150];
      double N_pdf_CTEQ66_isr_up, N_pdf_MRST2006NNLO_isr_up, N_pdf_NNPDF10100_isr_up;
      double N_pdf_CTEQ66_isr_down, N_pdf_MRST2006NNLO_isr_down, N_pdf_NNPDF10100_isr_down;
      
      eff->SetBranchStatus("*", 0);
      eff->SetBranchStatus("Npassed_In", 1);
      eff->SetBranchStatus("Npassed_ISR", 1);
      eff->SetBranchStatus("Npassed_ISR_up", 1);
      eff->SetBranchStatus("Npassed_ISR_down", 1);
      
      eff->SetBranchStatus("nCTEQ66", 1);
      eff->SetBranchStatus("N_pdf_CTEQ66", 1);
      eff->SetBranchStatus("N_pdf_CTEQ66_isr", 1);
      eff->SetBranchStatus("N_pdf_CTEQ66_isr_up", 1);
      eff->SetBranchStatus("N_pdf_CTEQ66_isr_down", 1);
      
      eff->SetBranchStatus("nMRST2006NNLO", 1);
      eff->SetBranchStatus("N_pdf_MRST2006NNLO", 1);
      eff->SetBranchStatus("N_pdf_MRST2006NNLO_isr", 1);
      eff->SetBranchStatus("N_pdf_MRST2006NNLO_isr_up", 1);
      eff->SetBranchStatus("N_pdf_MRST2006NNLO_isr_down", 1);
      
      eff->SetBranchStatus("nNNPDF10100", 1);
      eff->SetBranchStatus("N_pdf_NNPDF10100", 1);
      eff->SetBranchStatus("N_pdf_NNPDF10100_isr", 1);
      eff->SetBranchStatus("N_pdf_NNPDF10100_isr_up", 1);
      eff->SetBranchStatus("N_pdf_NNPDF10100_isr_down", 1);

      //Addresses 
      eff->SetBranchAddress("Npassed_In", &Npassed_In);
      eff->SetBranchAddress("Npassed_ISR", &Npassed_ISR);
      eff->SetBranchAddress("Npassed_ISR_up", &Npassed_ISR_up);
      eff->SetBranchAddress("Npassed_ISR_down", &Npassed_ISR_down);
      
      eff->SetBranchAddress("nCTEQ66", &nCTEQ66);
      eff->SetBranchAddress("N_pdf_CTEQ66",  N_CTEQ66);
      eff->SetBranchAddress("N_pdf_CTEQ66_isr", N_CTEQ66_isr);
      eff->SetBranchAddress("N_pdf_CTEQ66_isr_up", &N_pdf_CTEQ66_isr_up);
      eff->SetBranchAddress("N_pdf_CTEQ66_isr_down", &N_pdf_CTEQ66_isr_down);
      
      eff->SetBranchAddress("nMRST2006NNLO", &nMRST2006NNLO);
      eff->SetBranchAddress("N_pdf_MRST2006NNLO", N_MRST2006NNLO);
      eff->SetBranchAddress("N_pdf_MRST2006NNLO_isr", N_MRST2006NNLO_isr); 
      eff->SetBranchAddress("N_pdf_MRST2006NNLO_isr_up", &N_pdf_MRST2006NNLO_isr_up);
      eff->SetBranchAddress("N_pdf_MRST2006NNLO_isr_down", &N_pdf_MRST2006NNLO_isr_down);
      
      eff->SetBranchAddress("nNNPDF10100", &nNNPDF10100);
      eff->SetBranchAddress("N_pdf_NNPDF10100", N_NNPDF10100);
      eff->SetBranchAddress("N_pdf_NNPDF10100_isr", N_NNPDF10100_isr);
      eff->SetBranchAddress("N_pdf_NNPDF10100_isr_up", &N_pdf_NNPDF10100_isr_up);
      eff->SetBranchAddress("N_pdf_NNPDF10100_isr_down", &N_pdf_NNPDF10100_isr_down);
      
      int N_eff = eff->GetEntries();
      int Gen_Evts = 0;
      int Gen_Evts_isr = 0;
      int Gen_Evts_isrUp = 0;
      int Gen_Evts_isrDown = 0;
      
      double N_pdf_CTEQ66_Tot[100];
      double N_pdf_CTEQ66_isr_Tot[100];
      double N_pdf_MRST2006NNLO_Tot[100];
      double N_pdf_MRST2006NNLO_isr_Tot[100];
      double N_pdf_NNPDF10100_Tot[150];
      double N_pdf_NNPDF10100_isr_Tot[150];
      
      for(int k = 0; k < 100; k++){
	N_pdf_CTEQ66_Tot[k] = 0.0;
	N_pdf_CTEQ66_isr_Tot[k] = 0.0;
	N_pdf_MRST2006NNLO_Tot[k] = 0.0;
	N_pdf_MRST2006NNLO_isr_Tot[k] = 0.0;
      }
      for(int k = 0; k < 150; k++){
	N_pdf_NNPDF10100_Tot[k] = 0.0;
	N_pdf_NNPDF10100_isr_Tot[k] = 0.0;
      }
      
      for(int i = 0; i < N_eff; i++){
	eff->GetEntry(i);
	Gen_Evts += Npassed_In;
	Gen_Evts_isr += Npassed_ISR;
	Gen_Evts_isrUp += Npassed_ISR_up;
	Gen_Evts_isrDown += Npassed_ISR_down;
	std::cout << "nCTEQ66: " << nCTEQ66 << std::endl; 
	for(int l = 0; l < 45; l++){
	  N_pdf_CTEQ66_Tot[l] += N_CTEQ66[l];
	  N_pdf_CTEQ66_isr_Tot[l] += N_CTEQ66_isr[l];
	}
	for(int l = 0; l < 31; l++){
	  N_pdf_MRST2006NNLO_Tot[l] += N_MRST2006NNLO[l];
	  N_pdf_MRST2006NNLO_isr_Tot[l] += N_MRST2006NNLO_isr[l];
	}
	for(int l = 0; l < 101; l++){
	  N_pdf_NNPDF10100_Tot[l] += N_NNPDF10100[l];
	  N_pdf_NNPDF10100_isr_Tot[l] += N_NNPDF10100_isr[l];
	}
      }
      
      TString DM_Name = dm_sample.c_str();
      TString sn;      
      for(int i = 0; i < 4; i++){
	sn = TString(Form("_acc_cat%d",i+1));
	h_acc[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	for(int l = 0; l < 45; l++){
	  sn = TString(Form("_acc_pdfCTEQ_cat%d_index%d",i+1,l));
	  h_acc_pdf_CTEQ[xs_counter][l][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	}
	sn = TString(Form("_acc_pdfCTEQ_up_cat%d",i+1));
	h_acc_pdf_CTEQ_up[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	sn = TString(Form("_acc_pdfCTEQ_down_cat%d",i+1));
	h_acc_pdf_CTEQ_down[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	
	for(int l = 0; l < 31; l++){
	  sn = TString(Form("_acc_pdfMRST_cat%d_index%d",i+1, l));
	  h_acc_pdf_MRST[xs_counter][l][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	}
	sn = TString(Form("_acc_pdfMRST_up_cat%d",i+1));
	h_acc_pdf_MRST_up[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	sn = TString(Form("_acc_pdfMRST_down_cat%d",i+1));
	h_acc_pdf_MRST_down[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	
	for(int l = 0; l < 101; l++){
	  sn = TString(Form("_acc_pdfNNPDF_cat%d_index%d",i+1, l));
	  h_acc_pdf_NNPDF[xs_counter][l][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	}
	sn = TString(Form("_acc_pdfNNPDF_up_cat%d",i+1));
	h_acc_pdf_NNPDF_up[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	sn = TString(Form("_acc_pdfNNPDF_down_cat%d",i+1));
	h_acc_pdf_NNPDF_down[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));

	sn = TString(Form("_acc_CENTRAL_cat%d",i+1));
	h_acc_CENTRAL[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	sn = TString(Form("_acc_TOTAL_PDF_up_cat%d",i+1));
	h_acc_PDF_up[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	sn = TString(Form("_acc_TOTAL_PDF_down_cat%d",i+1));
	h_acc_PDF_down[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	
	//JES
	sn = TString(Form("_acc_JES_up_cat%d",i+1));
	h_acc_JES_up[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	sn = TString(Form("_acc_JES_down_cat%d",i+1));
	h_acc_JES_down[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));

	//ISR
	sn = TString(Form("_acc_ISR_up_cat%d",i+1));
	h_acc_ISR_up[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
	sn = TString(Form("_acc_ISR_down_cat%d",i+1));
	h_acc_ISR_down[xs_counter][i] = new TH1F(DM_Name+sn, DM_Name+sn, r2B[i], v.at(i));
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

      out->SetBranchStatus("nCTEQ66", 1);
      out->SetBranchStatus("wCTEQ66", 1);
      out->SetBranchStatus("nMRST2006NNLO", 1);
      out->SetBranchStatus("wMRST2006NNLO", 1);
      out->SetBranchStatus("nNNPDF10100", 1);
      out->SetBranchStatus("wNNPDF10100", 1);
      

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
      
      out->SetBranchAddress("nCTEQ66", &nCTEQ66);
      out->SetBranchAddress("wCTEQ66", wCTEQ66);
      out->SetBranchAddress("nMRST2006NNLO", &nMRST2006NNLO);
      out->SetBranchAddress("wMRST2006NNLO", wMRST2006NNLO);
      out->SetBranchAddress("nNNPDF10100", &nNNPDF10100);
      out->SetBranchAddress("wNNPDF10100", wNNPDF10100);
      
      int N_out = out->GetEntries();
      //double Lumi = 18.51;//PromptReco
      double Lumi = 18.836;//Jan22Rereco
      //double scaleF = Lumi*1000./Gen_Evts;//Scale to 1 pb
      double scaleF = Lumi*1000./Gen_Evts_isr;
      double scaleF_up = Lumi*1000./Gen_Evts_isrUp;
      double scaleF_down = Lumi*1000./Gen_Evts_isrDown;
      double lumi_SF = Lumi*1000.;
      

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
	  //double isr_up =  hlt_w*pu_w*lumi_SF*ISR_up*(wCTEQ66[0]/N_pdf_CTEQ66_isr_up + wMRST2006NNLO[0]/N_pdf_MRST2006NNLO_isr_up)/2.0;
	  //double isr_down = hlt_w*pu_w*lumi_SF*ISR_down*(wCTEQ66[0]/N_pdf_CTEQ66_isr_down + wMRST2006NNLO[0]/N_pdf_MRST2006NNLO_isr_down)/2.0;
	  //double ff_w = hlt_w*pu_w*ISR*lumi_SF;
	  double isr_up =  hlt_w*pu_w*ISR_up*(wCTEQ66[0]/N_pdf_CTEQ66_isr_up + wMRST2006NNLO[0]/N_pdf_MRST2006NNLO_isr_up)/2.0;
	  double isr_down = hlt_w*pu_w*ISR_down*(wCTEQ66[0]/N_pdf_CTEQ66_isr_down + wMRST2006NNLO[0]/N_pdf_MRST2006NNLO_isr_down)/2.0;
	  double ff_w = hlt_w*pu_w*ISR;
	  if(mr[2] > 200.0 && mr[2] <= 300.0 ){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_acc[xs_counter][0]->Fill(rsq[2], ff_w/Gen_Evts_isr);
	    h_acc_ISR_up[xs_counter][0]->Fill(rsq[2], isr_up);
	    h_acc_ISR_down[xs_counter][0]->Fill(rsq[2], isr_down);
	    for(int l = 0; l < nCTEQ66; l++)h_acc_pdf_CTEQ[xs_counter][l][0]->Fill(rsq[2], ff_w*wCTEQ66[l]/N_pdf_CTEQ66_isr_Tot[l]);
	    for(int l = 0; l < nMRST2006NNLO; l++)h_acc_pdf_MRST[xs_counter][l][0]->Fill(rsq[2], ff_w*wMRST2006NNLO[l]/N_pdf_MRST2006NNLO_isr_Tot[l]);
	    for(int l = 0; l < nNNPDF10100; l++)h_acc_pdf_NNPDF[xs_counter][l][0]->Fill(rsq[2], ff_w*wNNPDF10100[l]/N_pdf_NNPDF10100_isr_Tot[l]);
	  }else if(mr[2] > 300.0 && mr[2] <= 400.0){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_acc[xs_counter][1]->Fill(rsq[2], ff_w/Gen_Evts_isr);
	    h_acc_ISR_up[xs_counter][1]->Fill(rsq[2], isr_up);
	    h_acc_ISR_down[xs_counter][1]->Fill(rsq[2], isr_down);
	    for(int l = 0; l < nCTEQ66; l++)h_acc_pdf_CTEQ[xs_counter][l][1]->Fill(rsq[2], ff_w*wCTEQ66[l]/N_pdf_CTEQ66_isr_Tot[l]);
	    for(int l = 0; l < nMRST2006NNLO; l++)h_acc_pdf_MRST[xs_counter][l][1]->Fill(rsq[2], ff_w*wMRST2006NNLO[l]/N_pdf_MRST2006NNLO_isr_Tot[l]);
	    for(int l = 0; l < nNNPDF10100; l++)h_acc_pdf_NNPDF[xs_counter][l][1]->Fill(rsq[2], ff_w*wNNPDF10100[l]/N_pdf_NNPDF10100_isr_Tot[l]);
	  }else if(mr[2] > 400.0 && mr[2] <= 600.0){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_acc[xs_counter][2]->Fill(rsq[2], ff_w/Gen_Evts_isr);
	    h_acc_ISR_up[xs_counter][2]->Fill(rsq[2], isr_up);
	    h_acc_ISR_down[xs_counter][2]->Fill(rsq[2], isr_down);
	    for(int l = 0; l < nCTEQ66; l++)h_acc_pdf_CTEQ[xs_counter][l][2]->Fill(rsq[2], ff_w*wCTEQ66[l]/N_pdf_CTEQ66_isr_Tot[l]);
	    for(int l = 0; l < nMRST2006NNLO; l++)h_acc_pdf_MRST[xs_counter][l][2]->Fill(rsq[2], ff_w*wMRST2006NNLO[l]/N_pdf_MRST2006NNLO_isr_Tot[l]);
	    for(int l = 0; l < nNNPDF10100; l++)h_acc_pdf_NNPDF[xs_counter][l][2]->Fill(rsq[2], ff_w*wNNPDF10100[l]/N_pdf_NNPDF10100_isr_Tot[l]);
	  }else if(mr[2] > 600.0 && mr[2] <= 3500.0){
	    if(rsq[2] > 1.2)rsq[2] = 1.1;
	    h_acc[xs_counter][3]->Fill(rsq[2], ff_w/Gen_Evts_isr);
	    h_acc_ISR_up[xs_counter][3]->Fill(rsq[2], isr_up);
	    h_acc_ISR_down[xs_counter][3]->Fill(rsq[2], isr_down);
	    for(int l = 0; l < nCTEQ66; l++)h_acc_pdf_CTEQ[xs_counter][l][3]->Fill(rsq[2], ff_w*wCTEQ66[l]/N_pdf_CTEQ66_isr_Tot[l]);
	    for(int l = 0; l < nMRST2006NNLO; l++)h_acc_pdf_MRST[xs_counter][l][3]->Fill(rsq[2], ff_w*wMRST2006NNLO[l]/N_pdf_MRST2006NNLO_isr_Tot[l]);
	    for(int l = 0; l < nNNPDF10100; l++)h_acc_pdf_NNPDF[xs_counter][l][3]->Fill(rsq[2], ff_w*wNNPDF10100[l]/N_pdf_NNPDF10100_isr_Tot[l]);
	  }
	  N_passed += hlt_w*pu_w;
	  N_passed_ISR += hlt_w*pu_w*ISR;
	}

	//JES UP
	j1.SetPtEtaPhiE(pTHem1_up, etaHem1_up, phiHem1_up, pTHem1_up*cosh(etaHem1_up));//Hemisphere
	j2.SetPtEtaPhiE(pTHem2_up, etaHem2_up, phiHem2_up, pTHem2_up*cosh(etaHem2_up));//Hemisphere
	Dphi = j1.DeltaPhi(j2);
	if(mr_up[2] >= 200.0 && rsq_up[2] >= 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5){
	  //double weight = hlt_w*pu_w*lumi_SF*ISR*(wCTEQ66[0]/N_pdf_CTEQ66_isr_Tot[0] + wMRST2006NNLO[0]/N_pdf_MRST2006NNLO_isr_Tot[0])/2.0;
	  double weight = hlt_w*pu_w*ISR*(wCTEQ66[0]/N_pdf_CTEQ66_isr_Tot[0] + wMRST2006NNLO[0]/N_pdf_MRST2006NNLO_isr_Tot[0])/2.0;
	  if(mr_up[2] > 200.0 && mr_up[2] <= 300.0 ){
	    if(rsq_up[2] > 1.2)rsq_up[2] = 1.1;
	    h_acc_JES_up[xs_counter][0]->Fill(rsq_up[2], weight);
	  }else if(mr_up[2] > 300.0 && mr_up[2] <= 400.0){
	    if(rsq_up[2] > 1.2)rsq_up[2] = 1.1;
	    h_acc_JES_up[xs_counter][1]->Fill(rsq_up[2], weight);
	  }else if(mr_up[2] > 400.0 && mr_up[2] <= 600.0){
	    if(rsq_up[2] > 1.2)rsq_up[2] = 1.1;
	    h_acc_JES_up[xs_counter][2]->Fill(rsq_up[2], weight);
	  }else if(mr_up[2] > 600.0 && mr_up[2] <= 3500.0){
	    if(rsq_up[2] > 1.2)rsq_up[2] = 1.1;
	    h_acc_JES_up[xs_counter][3]->Fill(rsq_up[2], weight);
	  }
	}
	
	//JES DOWN
	j1.SetPtEtaPhiE(pTHem1_down, etaHem1_down, phiHem1_down, pTHem1_down*cosh(etaHem1_down));//Hemisphere
	j2.SetPtEtaPhiE(pTHem2_down, etaHem2_down, phiHem2_down, pTHem2_down*cosh(etaHem2_down));//Hemisphere
	Dphi = j1.DeltaPhi(j2);
	if(mr_down[2] >= 200.0 && rsq_down[2] >= 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5){
	  //double weight = hlt_w*pu_w*lumi_SF*ISR*(wCTEQ66[0]/N_pdf_CTEQ66_isr_Tot[0] + wMRST2006NNLO[0]/N_pdf_MRST2006NNLO_isr_Tot[0])/2.0;
	  double weight = hlt_w*pu_w*ISR*(wCTEQ66[0]/N_pdf_CTEQ66_isr_Tot[0] + wMRST2006NNLO[0]/N_pdf_MRST2006NNLO_isr_Tot[0])/2.0;
	  if(mr_down[2] > 200.0 && mr_down[2] <= 300.0 ){
	    if(rsq_down[2] > 1.2)rsq_down[2] = 1.1;
	    h_acc_JES_down[xs_counter][0]->Fill(rsq_down[2], weight);
	  }else if(mr_down[2] > 300.0 && mr_down[2] <= 400.0){
	    if(rsq_down[2] > 1.2)rsq_down[2] = 1.1;
	    h_acc_JES_down[xs_counter][1]->Fill(rsq_down[2], weight);
	  }else if(mr_down[2] > 400.0 && mr_down[2] <= 600.0){
	    if(rsq_down[2] > 1.2)rsq_down[2] = 1.1;
	    h_acc_JES_down[xs_counter][2]->Fill(rsq_down[2], weight);
	  }else if(mr_down[2] > 600.0 && mr_down[2] <= 3500.0){
	    if(rsq_down[2] > 1.2)rsq_down[2] = 1.1;
	    h_acc_JES_down[xs_counter][3]->Fill(rsq_down[2], weight);
	  }
	}
      }

      
      
      double sample_eff_old = N_passed/Gen_Evts;
      double sample_eff = N_passed_ISR/Gen_Evts_isr;
      double sample_eff_up = N_passed_ISR_up/Gen_Evts_isrUp;
      double sample_eff_down = N_passed_ISR_down/Gen_Evts_isrDown;

      TString data_card_name1 = dm_sample.c_str();
      TString s1 = data_card_name1;
      
      //CTEQ66 Error Calculation
      int npair = 22;
      for(int cat = 0; cat < 4; cat++){
	for(int bin = 1; bin <= h_acc_pdf_CTEQ[xs_counter][0][cat]->GetNbinsX(); bin++){
	  double wplus = 0.0;
	  double wminus = 0.0;
	  for(int l = 0; l < npair; l++){
	    double wa = 0.0;
	    if(h_acc_pdf_CTEQ[xs_counter][2*l+1][cat]->GetBinContent(bin)){
	      wa = h_acc_pdf_CTEQ[xs_counter][2*l+1][cat]->GetBinContent(bin)/h_acc_pdf_CTEQ[xs_counter][0][cat]->GetBinContent(bin) - 1.0;
	    }
	    double wb = 0.0;
	    if(h_acc_pdf_CTEQ[xs_counter][2*l+2][cat]->GetBinContent(bin)){
	      wb = h_acc_pdf_CTEQ[xs_counter][2*l+2][cat]->GetBinContent(bin)/h_acc_pdf_CTEQ[xs_counter][0][cat]->GetBinContent(bin) - 1.0;
	    }
	    
	    if(!std::isfinite(wa) || !std::isfinite(wb))continue;
	    
	    if (wa>wb) {
	      if (wa<0.) wa = 0.;
	      if (wb>0.) wb = 0.;
	      wplus += wa*wa;
	      wminus += wb*wb;
	    } else {
	      if (wb<0.) wb = 0.;
	      if (wa>0.) wa = 0.;
	      wplus += wb*wb;
	      wminus += wa*wa;
	    }
	    
	  }
	  wplus = sqrt(wplus);
	  wminus = sqrt(wminus);
	  std::cout << "CAT: " << cat+1 << " Bin: " << bin << " DeltaUp: " << wplus*100.0 << " % " <<
	    " DeltaDown: " << wminus*100.0 << " %"<< std::endl;
	  double acc_up = h_acc_pdf_CTEQ[xs_counter][0][cat]->GetBinContent(bin)+wplus*h_acc_pdf_CTEQ[xs_counter][0][cat]->GetBinContent(bin);
	  double acc_down = h_acc_pdf_CTEQ[xs_counter][0][cat]->GetBinContent(bin)-wminus*h_acc_pdf_CTEQ[xs_counter][0][cat]->GetBinContent(bin);
	  h_acc_pdf_CTEQ_up[xs_counter][cat]->SetBinContent(bin, acc_up);
	  h_acc_pdf_CTEQ_down[xs_counter][cat]->SetBinContent(bin, acc_down);
	}
      }
      
      //MRST Error Calculation
      npair = 15;
      for(int cat = 0; cat < 4; cat++){
	for(int bin = 1; bin <= h_acc_pdf_MRST[xs_counter][0][cat]->GetNbinsX(); bin++){
	  double wplus = 0.0;
	  double wminus = 0.0;
	  for(int l = 0; l < npair; l++){
	    double wa = 0.0;
	    if(h_acc_pdf_MRST[xs_counter][2*l+1][cat]->GetBinContent(bin)){
	      wa = h_acc_pdf_MRST[xs_counter][2*l+1][cat]->GetBinContent(bin)/h_acc_pdf_MRST[xs_counter][0][cat]->GetBinContent(bin) - 1.0;
	    }
	    double wb = 0.0;
	    if(h_acc_pdf_MRST[xs_counter][2*l+2][cat]->GetBinContent(bin)){
	      wb = h_acc_pdf_MRST[xs_counter][2*l+2][cat]->GetBinContent(bin)/h_acc_pdf_MRST[xs_counter][0][cat]->GetBinContent(bin) - 1.0;
	    }
	    
	    if(!std::isfinite(wa) || !std::isfinite(wb))continue;
	    
	    if (wa>wb) {
	      if (wa<0.) wa = 0.;
	      if (wb>0.) wb = 0.;
	      wplus += wa*wa;
	      wminus += wb*wb;
	    } else {
	      if (wb<0.) wb = 0.;
	      if (wa>0.) wa = 0.;
	      wplus += wb*wb;
	      wminus += wa*wa;
	    }
	    
	  }
	  wplus = sqrt(wplus);
	  wminus = sqrt(wminus);
	  std::cout << "CAT: " << cat+1 << " Bin: " << bin << " DeltaUp: " << wplus*100.0 << " % " <<
	    " DeltaDown: " << wminus*100.0 << " %"<< std::endl;
	  double acc_up = h_acc_pdf_MRST[xs_counter][0][cat]->GetBinContent(bin)+wplus*h_acc_pdf_MRST[xs_counter][0][cat]->GetBinContent(bin);
	  double acc_down = h_acc_pdf_MRST[xs_counter][0][cat]->GetBinContent(bin)-wminus*h_acc_pdf_MRST[xs_counter][0][cat]->GetBinContent(bin);
	  h_acc_pdf_MRST_up[xs_counter][cat]->SetBinContent(bin, acc_up);
	  h_acc_pdf_MRST_down[xs_counter][cat]->SetBinContent(bin, acc_down);
	}
      }
      
      for(int cat = 0; cat < 4; cat++){
	for(int bin = 1; bin <= h_acc_pdf_MRST[xs_counter][0][cat]->GetNbinsX(); bin++){
	  double cteqUp = h_acc_pdf_CTEQ_up[xs_counter][cat]->GetBinContent(bin);
	  double cteqDown = h_acc_pdf_CTEQ_down[xs_counter][cat]->GetBinContent(bin);
	  double mrstUp = h_acc_pdf_MRST_up[xs_counter][cat]->GetBinContent(bin);
	  double mrstDown = h_acc_pdf_MRST_down[xs_counter][cat]->GetBinContent(bin);
	  
	  double MAX = -99.0;
	  double MIN = -99.0;
	  if(cteqUp > mrstUp){
	    MAX = cteqUp;
	  }else{
	    MAX = mrstUp;
	  }
	  
	  if(cteqUp < mrstUp){
	    MIN = cteqDown;
	  }else{
	    MIN = mrstDown;
	  }
	  
	  h_acc_CENTRAL[xs_counter][cat]->SetBinContent(bin, (MIN+MAX)/2.0);
	  h_acc_PDF_up[xs_counter][cat]->SetBinContent(bin, (MAX - (MIN+MAX)/2.0)/((MIN+MAX)/2.0));
	  h_acc_PDF_down[xs_counter][cat]->SetBinContent(bin, ((MIN+MAX)/2.0 - MIN)/((MIN+MAX)/2.0));
	}
      }
      

      TFile* fo = new TFile("AccBIS/"+s1+"_Acc.root", "RECREATE");
      TString dummy;
      for(int i = 0; i < 4; i++){
	//h_acc[xs_counter][i]->Write();
	
	//h_acc_pdf_CTEQ[xs_counter][0][i]->Write();
	//h_acc_pdf_CTEQ_up[xs_counter][i]->Write();
	//h_acc_pdf_CTEQ_down[xs_counter][i]->Write();
	
	//h_acc_pdf_MRST[xs_counter][0][i]->Write();
	//h_acc_pdf_MRST_up[xs_counter][i]->Write();
	//h_acc_pdf_MRST_down[xs_counter][i]->Write();

	//h_acc_CENTRAL[xs_counter][i]->Write();
	dummy = Form("PDF_SYS_cat%d",i+1);
	h_acc_PDF_up[xs_counter][i]->Write(dummy);
	//h_acc_PDF_down[xs_counter][i]->Write();
	
	//h_acc_JES_up[xs_counter][i]->Write();
	//h_acc_JES_down[xs_counter][i]->Write();
	
	//h_acc_ISR_up[xs_counter][i]->Write();
	//h_acc_ISR_down[xs_counter][i]->Write();
	}
	
      /*
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
	data_card_f1 << "rate\t\t"<< h_acc_CENTRAL[xs_counter][i]->Integral() <<"\t\t"<< tt_N[i] <<"\t\t" << dy_N[i] << "\t\t" << z_N[i] << "\t\t" << w_N[i] << "\n";
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
	
	h_acc_CENTRAL[xs_counter][i]->Write("signal_rsq");
		
	h_acc_ISR_up[xs_counter][i]->Write("signal_rsq_IsrUp");
	h_acc_ISR_down[xs_counter][i]->Write("signal_rsq_IsrDown");

	h_acc_JES_up[xs_counter][i]->Write("signal_rsq_JesUp");
	h_acc_JES_down[xs_counter][i]->Write("signal_rsq_JesDown");
	
	h_acc_PDF_up[xs_counter][i]->Write("signal_rsq_PdfUp");
	h_acc_PDF_down[xs_counter][i]->Write("signal_rsq_PdfDown");
	
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
	h_acc[xs_counter][i]->Write();
	fo->Close();
	
	}*/
    }
    std::cout << " deb 1" << std::endl;
    xs_counter++;
  }else{
    std::cout << "Unable to open the file" << std::endl;
  }
  std::cout << " deb 2" << std::endl;
  mfile0.close();
  std::cout << " deb 2.1" << std::endl;
  return 0;
  std::cout << " deb 2.2" << std::endl;
}
