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

int main(){

  gROOT->Reset();
  //5x5
  //int MR_nb = 5;
  //int R2_nb = 5;
  //4x4
  int MR_nb = 4;
  int R2_nb = 4;


  TCanvas* ca = new TCanvas("c","c", 640, 640);
  //TFile* f2 = new TFile("trigger/hlt_eff_SignleElePD_Final.root");
  TFile* f2 = new TFile("trigger/hlt_eff_HTMHT_0mu_BOX_0bTag_Final.root");
  TEfficiency* hlt = (TEfficiency*)f2->Get("Eff2d");
  TH2F* mr_rsq;

  
  //Getting 2D plots from Prediction
  TFile* in = new TFile("/Users/cmorgoth/Software/git/BkgPredictionDM/Pred_Files/BkgPred_ttMC_LO_RunAB_2sys.root");
  
  TH1F* tt_1D = (TH1F*)in->Get("tt_1D");
  TH1F* tt_1D_alphaUp = (TH1F*)in->Get("tt_1D_alphaUp");
  TH1F* tt_1D_alphaDown = (TH1F*)in->Get("tt_1D_alphaDown");
  TH1F* z_1D = (TH1F*)in->Get("z_1D");
  TH1F* z_1D_alphaUp = (TH1F*)in->Get("z_1D_alphaUp");
  TH1F* z_1D_alphaDown = (TH1F*)in->Get("z_1D_alphaDown");
  TH1F* dy_1D = (TH1F*)in->Get("dy_1D");
  TH1F* dy_1D_alphaUp = (TH1F*)in->Get("dy_1D_alphaUp");
  TH1F* dy_1D_alphaDown = (TH1F*)in->Get("dy_1D_alphaDown");
  TH1F* w_1D = (TH1F*)in->Get("w_1D");
  TH1F* w_1D_alphaUp = (TH1F*)in->Get("w_1D_alphaUp");
  TH1F* w_1D_alphaDown = (TH1F*)in->Get("w_1D_alphaDown");
  TH1F* data_obs = (TH1F*)in->Get("data_obs");

  double tt_N = tt_1D->Integral();
  double dy_N = dy_1D->Integral();
  double z_N = z_1D->Integral();
  double w_N = w_1D->Integral();

  //Getting 1D plots from Prediction
  TFile* in_rsq = new TFile("/Users/cmorgoth/Software/git/BkgPredictionDM/Pred_Files/Bkg_Pred_1D_Rsq_ttMC_LO_RunAB_2sys.root");
  
  TH1F* tt_rsq = (TH1F*)in_rsq->Get("tt_rsq");
  TH1F* tt_rsq_alphaUp = (TH1F*)in_rsq->Get("tt_rsq_betaUp");
  TH1F* tt_rsq_alphaDown = (TH1F*)in_rsq->Get("tt_rsq_betaDown");
  TH1F* z_rsq = (TH1F*)in_rsq->Get("z_rsq");
  TH1F* z_rsq_alphaUp = (TH1F*)in_rsq->Get("z_rsq_deltaUp");
  TH1F* z_rsq_alphaDown = (TH1F*)in_rsq->Get("z_rsq_deltaDown");
  TH1F* dy_rsq = (TH1F*)in_rsq->Get("dy_rsq");
  TH1F* dy_rsq_alphaUp = (TH1F*)in_rsq->Get("dy_rsq_gammaUp");
  TH1F* dy_rsq_alphaDown = (TH1F*)in_rsq->Get("dy_rsq_gammaDown");
  TH1F* w_rsq = (TH1F*)in_rsq->Get("w_rsq");
  TH1F* w_rsq_alphaUp = (TH1F*)in_rsq->Get("w_rsq_zetaUp");
  TH1F* w_rsq_alphaDown = (TH1F*)in_rsq->Get("w_rsq_zetaDown");
  TH1F* data_rsq = (TH1F*)in_rsq->Get("data_rsq");

  double tt_N1 = tt_rsq->Integral();
  double dy_N1 = dy_rsq->Integral();
  double z_N1 = z_rsq->Integral();
  double w_N1 = w_rsq->Integral();


  
  //Here the program starts
  //double xsec[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
  std::ifstream mfile0("list_of_files_v2.list");
  std::ofstream outfile("eff_table_normal_R2_0p5_MR_200_Dphi_B_2p5.tex");
  
  outfile << "\\begin{table}[htdp]\n\\caption{default}\n\\begin{center}\n\\begin{tabular}{|c|c|}\n\\hline\n";
  
  std::string fname0;
  std::cout.precision(16);
  int xs_counter = 0;
  TH2F* h_2d[24];
  TH1F* h_rsq[24];
  TH1F* h_rsq_up[24];
  TH1F* h_rsq_down[24];
  
  TH2F* h_met_deltaPhi[24];
  TH2F* h_r2_deltaPhi[24];
  TH2F* h_r2_alpha[24];
  TH2F* h_r2_alphaT[24];
  TH2F* h_met_alpha[24];
  TH2F* h_met_alphaT[24];

  
  TH1F* h_M2[24];
  TH1F* h_M2_overPhi[24];
  
  TH2F* h_met_deltaPhi_B[24];
  TH2F* h_r2_deltaPhi_B[24];
  TH2F* h_r2_alpha_B[24];
  TH2F* h_r2_alphaT_B[24];
  TH2F* h_met_alpha_B[24];
  TH2F* h_met_alphaT_B[24];
  TH1F* h_R_DM[24];

  TH1F* h_1d[24];
  TH1F* h_1d_up[24];
  TH1F* h_1d_down[24];
  TH1F* h_alpha[24];
  TH1F* h_alphaT[24];
  TH1F* h_r2[24];
  TH1F* h_met[24];
  TH1F* h_theta[24];
  TH1F* h_Dphi[24];

  TH1F* h_alpha_B[24];
  TH1F* h_alphaT_B[24];
  TH1F* h_r2_B[24];
  TH1F* h_met_B[24];
  TH1F* h_theta_B[24];
  TH1F* h_Dphi_B[24];
  TH1F* h_beta_r[24];
  TH1F* h_R_DM_B[24];
  
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
      
      double mr[4], rsq[4], Jet_PT[20], Jet_Eta[20], Jet_Phi[20], Npassed_In, metCorrX[4], metCorrY[4];
      double pTHem1, pTHem2, etaHem1, etaHem2, phiHem1, phiHem2;
      int btag, box, N_Jets;
      
      eff->SetBranchStatus("*", 0);
      eff->SetBranchStatus("Npassed_In", 1);
      eff->SetBranchAddress("Npassed_In", &Npassed_In);
      
      int N_eff = eff->GetEntries();
      int Gen_Evts = 0;
      
      for(int i = 0; i < N_eff; i++){
	eff->GetEntry(i);
	Gen_Evts += Npassed_In;
      }
      
      std::cout << "Gen_Events: " << Gen_Evts << std::endl;
      
      out->SetBranchStatus("*", 0);
      out->SetBranchStatus("MR", 1);
      out->SetBranchStatus("RSQ",1);
      out->SetBranchStatus("nBtag", 1);
      out->SetBranchStatus("BOX_NUM",1);
      out->SetBranchStatus("N_Jets",1);
      out->SetBranchStatus("Jet_PT",1);
      out->SetBranchStatus("Jet_Phi",1);
      out->SetBranchStatus("Jet_Eta",1);
      out->SetBranchStatus("pTHem1",1);
      out->SetBranchStatus("pTHem2",1);
      out->SetBranchStatus("etaHem1",1);
      out->SetBranchStatus("etaHem2",1);
      out->SetBranchStatus("phiHem1",1);
      out->SetBranchStatus("phiHem2",1);
      out->SetBranchStatus("metCorrX",1);
      out->SetBranchStatus("metCorrY",1);
      out->SetBranchAddress("MR", mr);
      out->SetBranchAddress("RSQ", rsq);
      out->SetBranchAddress("nBtag", &btag);
      out->SetBranchAddress("BOX_NUM", &box);
      out->SetBranchAddress("N_Jets", &N_Jets);
      out->SetBranchAddress("Jet_PT", Jet_PT);
      out->SetBranchAddress("Jet_Phi", Jet_Phi);
      out->SetBranchAddress("Jet_Eta", Jet_Eta);
      out->SetBranchAddress("pTHem1", &pTHem1);
      out->SetBranchAddress("pTHem2", &pTHem2);
      out->SetBranchAddress("etaHem1", &etaHem1);
      out->SetBranchAddress("etaHem2", &etaHem2);
      out->SetBranchAddress("phiHem1", &phiHem1);
      out->SetBranchAddress("phiHem2", &phiHem2);
      out->SetBranchAddress("metCorrX", metCorrX);
      out->SetBranchAddress("metCorrY", metCorrY);
      
      int N_out = out->GetEntries();
      //5x5
      //int n_mr = 5;
      //int n_rsq = 5;
      //double mr_bin[] = {200, 300, 400, 600, 900, 3500};
      //double rsq_bin[] = {0.5, 0.6, 0.725, 0.85, 1.1, 2.5};
      //4x4
      int n_mr = 4;
      int n_rsq = 4;
      double mr_bin[] = {200, 300, 400, 600, 3500};
      double rsq_bin[] = {0.5, 0.6, 0.725, 0.85, 2.5};
      
      mr_rsq = new TH2F("mr_rsq", "mr_rsq", n_mr, mr_bin, n_rsq, rsq_bin);

      TString s_c = dm_sample.c_str();
      TString s_up = dm_sample.c_str();
      TString s_down = dm_sample.c_str();
      TString s_rsq = dm_sample.c_str();
      s_c = s_c + "_c";
      s_up = s_up + "_up";
      s_down = s_down + "s_down";
      
      h_2d[xs_counter] = new TH2F(dm_sample.c_str(), dm_sample.c_str(), n_mr, mr_bin, n_rsq, rsq_bin);
      
      h_rsq[xs_counter] = new TH1F(s_rsq + "_rsq_c", s_rsq + "_rsq_c", n_rsq, rsq_bin);
      h_rsq_up[xs_counter] = new TH1F(s_rsq + "_rsq_up", s_rsq + "_rsq_up", n_rsq, rsq_bin);
      h_rsq_down[xs_counter] = new TH1F(s_rsq + "_rsq_down", s_rsq + "_rsq_down", n_rsq, rsq_bin);
      
      h_1d[xs_counter] = new TH1F(s_c, s_c, n_mr*n_rsq, 0, n_mr*n_rsq);
      h_1d_up[xs_counter] = new TH1F(s_up, s_up, n_mr*n_rsq, 0, n_mr*n_rsq);
      h_1d_down[xs_counter] = new TH1F(s_down, s_down, n_mr*n_rsq, 0, n_mr*n_rsq);
      
      double Lumi = 18.5;
      //double xsec = 1;
      //double scaleF = Lumi*xsec[xs_counter]*1000./Gen_Evts;
      double scaleF = Lumi*1000./Gen_Evts;//Scale to 1 pb
      //double scaleF = 1;
      //std::cout << dm_sample << " " << xsec[xs_counter] << std::endl;
      
      //std::cout << "scaleF : " << scaleF << std::endl;
      
      TString dmsampleTS = dm_sample.c_str();
      TString met_dphi_name = dm_sample.c_str();
      TString r2_dphi_name = dm_sample.c_str();
      TString r2_alpha_name = dm_sample.c_str();
      TString r2_alphaT_name = dm_sample.c_str();
      TString met_alpha_name = dm_sample.c_str();
      TString met_alphaT_name = dm_sample.c_str();
      TString R_DM_name = dm_sample.c_str();
      TString M2_name = dm_sample.c_str();
      TString M2_name_OverPhi = dm_sample.c_str();
      met_dphi_name += "_met_dphi"; 
      r2_dphi_name += "_r2_dphi"; 
      r2_alpha_name += "_r2_alpha"; 
      r2_alphaT_name += "_r2_alphaT"; 
      met_alpha_name += "_met_alpha"; 
      R_DM_name += "_R_DM"; 
      M2_name += "_M2";
      M2_name_OverPhi += "_M2_OverPhi";

      h_met_deltaPhi[xs_counter] = new TH2F(met_dphi_name, met_dphi_name, 100, -TMath::Pi(), TMath::Pi(), 100, 0.0, 1500);
      h_r2_deltaPhi[xs_counter] = new TH2F(r2_dphi_name, r2_dphi_name, 100, -TMath::Pi(), TMath::Pi(), 100, 0.0, 2.5);
      h_r2_alpha[xs_counter] = new TH2F(r2_alpha_name, r2_alpha_name, 100, -0.2, 5, 100, 0.0, 2.5);
      h_r2_alphaT[xs_counter] = new TH2F(r2_alphaT_name, r2_alphaT_name, 100, -0.2, 5, 100, 0.0, 2.5);
      h_met_alpha[xs_counter] = new TH2F(met_alpha_name, met_alpha_name, 100, -0.2, 5, 100, 0.0, 1500.);
      h_met_alphaT[xs_counter] = new TH2F(met_alphaT_name, met_alphaT_name, 100, -0.2, 5, 100, 0.0, 1500.);
      h_alpha[xs_counter] = new TH1F(met_alpha_name+"_1d", met_alpha_name+"_1d", 100, -0.2, 5);
      h_alphaT[xs_counter] = new TH1F(met_alphaT_name+"_1d", met_alphaT_name+"_1d", 100, -0.2, 5);
      h_r2[xs_counter] = new TH1F(dmsampleTS+"_r2", dmsampleTS+"_r2",100, 0.0, 2.5);
      h_met[xs_counter] = new TH1F(dmsampleTS+"_met", dmsampleTS+"_met",100, 0.0, 1500.);
      h_theta[xs_counter] = new TH1F(dmsampleTS+"_theta", dmsampleTS+"_theta", 100, 0, TMath::Pi());
      h_Dphi[xs_counter] = new TH1F(dmsampleTS+"_Dphi", dmsampleTS+"_Dphi", 100, -TMath::Pi(), TMath::Pi());
      h_R_DM[xs_counter] = new TH1F(R_DM_name, R_DM_name, 100, -1000, 10000);
      h_M2[xs_counter] = new TH1F(M2_name, M2_name, 1000, -1000, 1000);
      h_M2_overPhi[xs_counter] = new TH1F(M2_name_OverPhi, M2_name_OverPhi, 100, 0, 2000);
      
      //Boosted Histos
      h_met_deltaPhi_B[xs_counter] = new TH2F(met_dphi_name+"_B", met_dphi_name+"_B", 100, -TMath::Pi(), TMath::Pi(), 100, 0.0, 1500);
      h_r2_deltaPhi_B[xs_counter] = new TH2F(r2_dphi_name+"_B", r2_dphi_name+"_B", 100, -TMath::Pi(), TMath::Pi(), 100, 0.0, 2.5);
      h_r2_alpha_B[xs_counter] = new TH2F(r2_alpha_name+"_B", r2_alpha_name+"_B", 100, -0.2, 5, 100, 0.0, 2.5);
      h_r2_alphaT_B[xs_counter] = new TH2F(r2_alphaT_name+"_B", r2_alphaT_name+"_B", 100, -0.2, 5, 100, 0.0, 2.5);
      h_met_alpha_B[xs_counter] = new TH2F(met_alpha_name+"_B", met_alpha_name+"_B", 100, -0.2, 5, 100, 0.0, 1500.);
      h_met_alphaT_B[xs_counter] = new TH2F(met_alphaT_name+"_B", met_alphaT_name+"_B", 100, -0.2, 5, 100, 0.0, 1500.);
      h_alpha_B[xs_counter] = new TH1F(met_alpha_name+"_B"+"_1d", met_alpha_name+"_B"+"_1d", 100, -0.2, 5);
      h_alphaT_B[xs_counter] = new TH1F(met_alphaT_name+"_B"+"_1d", met_alphaT_name+"_B"+"_1d", 100, -0.2, 5);
      h_r2_B[xs_counter] = new TH1F(dmsampleTS+"_r2_B", dmsampleTS+"_r2_B",100, 0.0, 2.5);
      h_met_B[xs_counter] = new TH1F(dmsampleTS+"_met_B", dmsampleTS+"_met_B",100, 0.0, 1500.);
      h_theta_B[xs_counter] = new TH1F(dmsampleTS+"_theta_B", dmsampleTS+"_theta_B", 100, 0, TMath::Pi());
      h_Dphi_B[xs_counter] = new TH1F(dmsampleTS+"_Dphi_B", dmsampleTS+"_Dphi_B", 100, -TMath::Pi(), TMath::Pi());
      h_beta_r[xs_counter] = new TH1F(dmsampleTS+"_beta_r", dmsampleTS+"_beta_r", 200, -1.5, 1.5);
      h_R_DM_B[xs_counter] = new TH1F(R_DM_name+"_B", R_DM_name+"_B", 100, -1000, 10000);
      
      double N_passed = 0.0;
      for(int j = 0; j < N_out; j++){
	out->GetEntry(j);
	double hlt_w = HLTscale(mr[2], rsq[2], hlt);
	
	TLorentzVector j1;
	TLorentzVector j2;
	TLorentzVector sumJ;
	
	j1.SetPtEtaPhiE(Jet_PT[0], Jet_Eta[0], Jet_Phi[0], Jet_PT[0]*cosh(Jet_Eta[0]));
	j2.SetPtEtaPhiE(Jet_PT[1], Jet_Eta[1], Jet_Phi[1], Jet_PT[1]*cosh(Jet_Eta[1]));
	//j1.SetPtEtaPhiE(pTHem1, etaHem1, phiHem1, pTHem1*cosh(etaHem1));//Hemisphere
        //j2.SetPtEtaPhiE(pTHem2, etaHem2, phiHem2, pTHem2*cosh(etaHem2));//Hemisphere
	
	sumJ = j1+j2;
	//std::cout << "=========beta: " << sumJ.Beta() << "  ===============" << std::endl;
	double cm_beta = sumJ.Beta();
	TVector3 boost_dir = sumJ.Vect().Unit();
	//TVector3 CM_Boost(cm_beta*boost_dir.Px(),cm_beta*boost_dir.Py(),cm_beta*boost_dir.Pz());
	TVector3 CM_Boost(0.0,cm_beta*boost_dir.Py(),cm_beta*boost_dir.Pz());

	double Mt = sumJ.Mt();
	double Minv = sumJ.M();
	double alpha = j2.Pt()/Minv;
	//double alphaT = j2.Pt()/Mt;
	double alphaT = j2.Pt()/sqrt(sumJ.E()*sumJ.E() - (sumJ.Px()*sumJ.Px() + sumJ.Py()*sumJ.Py()) );
	double Dphi = j1.DeltaPhi(j2);
	double Dtheta = j1.Angle(j2.Vect());
	//std::cout << "PhiJ1: " << Jet_Phi[0] << " PhiJ1: " << Jet_Phi[1] <<  " normal: " << deltaPhi << " vector: " << Dphi << std::endl;
	double Met = sqrt(metCorrX[2]*metCorrX[2] + metCorrY[2]*metCorrY[2]);
	double Rdm = Met/(1-cos(Dtheta));
	
	h_met_deltaPhi[xs_counter]->Fill(Dphi,Met,scaleF);
	h_r2_deltaPhi[xs_counter]->Fill(Dphi,rsq[2],scaleF);
	h_r2_alpha[xs_counter]->Fill(alpha,rsq[2],scaleF);
	h_r2_alphaT[xs_counter]->Fill(alphaT,rsq[2],scaleF);
	h_met_alpha[xs_counter]->Fill(alpha,Met,scaleF);
	h_met_alphaT[xs_counter]->Fill(alphaT,Met,scaleF);
	h_alpha[xs_counter]->Fill(alpha,scaleF);
	h_alphaT[xs_counter]->Fill(alphaT,scaleF);
	h_r2[xs_counter]->Fill(rsq[2],scaleF);
	h_met[xs_counter]->Fill(Met,scaleF);
	h_theta[xs_counter]->Fill(Dtheta, scaleF);
	h_Dphi[xs_counter]->Fill(Dphi, scaleF);
	h_R_DM[xs_counter]->Fill(Rdm, scaleF);
	//Boost to the Razor Frame
	double beta_L = (j1.Pz() +  j2.Pz())/(j1.Vect().Mag() +  j2.Vect().Mag());
	double PT_t = (j1+j2).Pt();
	double beta_T_Star_X = (j1+j2).Px()/PT_t;
	double beta_T_Star_Y = (j1+j2).Py()/PT_t;
	double gammaL = 1.0/sqrt(1-beta_L*beta_L);
	TVector3 beta_hat_R(beta_T_Star_X, beta_T_Star_Y, 0);
	TVector3 JSum_T(sumJ.Px(), sumJ.Py(),0);
	double beta_R_Mag = (gammaL*(j1.Vect().Mag() -  j2.Vect().Mag())-gammaL*beta_L*(j1.Pz() -  j2.Pz()))/(beta_hat_R.Dot(JSum_T));
	
	TLorentzVector j1_B(j1);//j1 boosted
	TLorentzVector j2_B(j2);//j2 boosted
	TLorentzVector sumJB;
	
	
	
	j1_B.Boost(0.0, 0.0, -beta_L);
	j2_B.Boost(0.0, 0.0, -beta_L);
	j1_B.Boost(-beta_R_Mag*beta_T_Star_X, -beta_R_Mag*beta_T_Star_Y, 0.0);
	j2_B.Boost(beta_R_Mag*beta_T_Star_X, beta_R_Mag*beta_T_Star_Y, 0.0);
	
	//j1_B.Boost(CM_Boost);
	//j2_B.Boost(CM_Boost);

	//std::cout << "============beta: " << beta_L << " beta_R_Mag: "<< beta_R_Mag << " ==============="<<std::endl;
	//std::cout << "j1.E: " << j1.E() << " j1.Px: " << j1.Px() << " j2.Py: " << j1.Py() << " j1.Pz: " << j1.Pz() << std::endl;
	//std::cout << "j1_B.E: " << j1_B.E() << " j1_B.Px: " << j1_B.Px() << " j2.Py: " << j1_B.Py() << " j1_B.Pz: " << j1_B.Pz() << std::endl;

	//std::cout << "j2.E: " << j2.E() << " j2.Px: " << j2.Px() << " j2.Py: " << j2.Py() << " j2.Pz: " << j2.Pz() << std::endl;
	//std::cout << "j2_B.E: " << j2_B.E() << " j2_B.Px: " << j2_B.Px() << " j2.Py: " << j2_B.Py() << " j2_B.Pz: " << j2_B.Pz() << std::endl;
	

	sumJB = j1_B+j2_B;
	double Mt_B = sumJB.Mt();
	double Minv_B = sumJB.M();
	double alpha_B = j2_B.Pt()/Minv_B;
	double alphaT_B = j2_B.Pt()/sqrt(sumJB.E()*sumJB.E() - (sumJB.Px()*sumJB.Px() + sumJB.Py()*sumJB.Py()) );
	//double alphaT_B = j2_B.Pt()/Mt_B;
	double Dphi_B = j1_B.DeltaPhi(j2_B);
	double Dtheta_B = j1_B.Angle(j2_B.Vect());
	double Rdm_B = Met/(1-cos(Dtheta_B));
	//Filling Boosted Histos
	h_M2[xs_counter]->Fill(Mt_B, scaleF);
	h_M2_overPhi[xs_counter]->Fill(Mt_B/(1-cos(Dphi_B)), scaleF);
	h_met_deltaPhi_B[xs_counter]->Fill(Dphi_B,Met,scaleF);
	h_r2_deltaPhi_B[xs_counter]->Fill(Dphi_B,rsq[2],scaleF);
	h_r2_alpha_B[xs_counter]->Fill(alpha_B,rsq[2],scaleF);
	h_r2_alphaT_B[xs_counter]->Fill(alphaT_B,rsq[2],scaleF);
	h_met_alpha_B[xs_counter]->Fill(alpha_B,Met,scaleF);
	h_met_alphaT_B[xs_counter]->Fill(alphaT_B,Met,scaleF);
	h_alpha_B[xs_counter]->Fill(alpha_B,scaleF);
	h_alphaT_B[xs_counter]->Fill(alphaT_B,scaleF);
	h_r2_B[xs_counter]->Fill(rsq[2],scaleF);
	h_met_B[xs_counter]->Fill(Met,scaleF);
	h_theta_B[xs_counter]->Fill(Dtheta_B, scaleF);
	h_Dphi_B[xs_counter]->Fill(Dphi_B, scaleF);
	h_beta_r[xs_counter]->Fill(beta_L, scaleF);
	h_R_DM_B[xs_counter]->Fill(Rdm_B, scaleF);
	
	//if( mr[2]*rsq[2] > 100 && mr[2] > 200 && btag == 0 && box == 0 /*&& deltaPhi < 2.5*/)mr_rsq->Fill(mr[2], rsq[2], hlt_w);
	if(mr[2] > 200.0 && rsq[2] > 0.5 && btag == 0 && box == 0 && fabs(Dphi) < 2.5){
	  //if(btag == 0 && box == 0 && alphaT_B > 0.5 ){
	  //mr_rsq->Fill(mr[2], rsq[2], hlt_w);
	  h_2d[xs_counter]->Fill(mr[2], rsq[2], hlt_w*scaleF);
	  h_rsq[xs_counter]->Fill(rsq[2], hlt_w*scaleF);
	  N_passed += hlt_w;
	}
      }
      
      //Making Error Plots and Unwrapping
      int b_ctr = 1;
      for(int m = 1; m <= n_mr; m++){
	for(int n = 1; n <= n_rsq; n++){
	  int G_bin = h_2d[xs_counter]->GetBin(m,n,0);
	  double b_c = h_2d[xs_counter]->GetBinContent(m,n);
	  double b_up = h_2d[xs_counter]->GetBinContent(m,n) + h_2d[xs_counter]->GetBinError(G_bin);
	  double b_down = h_2d[xs_counter]->GetBinContent(m,n) - h_2d[xs_counter]->GetBinError(G_bin);
	  h_1d[xs_counter]->SetBinContent(b_ctr, b_c);
	  h_1d_up[xs_counter]->SetBinContent(b_ctr, b_up);
	  h_1d_down[xs_counter]->SetBinContent(b_ctr, b_down);
	  b_ctr++;
	}
      }
      
      //Making erro Plots 1D RSQ
      for(int z = 1; z <= n_rsq; z++){
	double bv = h_rsq[xs_counter]->GetBinContent(z);
	double be = h_rsq[xs_counter]->GetBinError(z);
	h_rsq_up[xs_counter]->SetBinContent(z,bv + be);
	h_rsq_down[xs_counter]->SetBinContent(z,bv - be);
      }


      ca->cd();
      h_2d[xs_counter]->SetStats(0);
      gPad->SetLogy(0);
      gPad->SetLogx(0);
      //gPad->SetLogz();
      h_2d[xs_counter]->Draw("colztext");
      ca->SaveAs("Plots/Signal/"+met_dphi_name+".pdf");
      ca->SaveAs("Plots/Signal/"+met_dphi_name+".png");
      
      h_met_deltaPhi[xs_counter]->SetStats(0);
      gPad->SetLogy(0);
      h_met_deltaPhi[xs_counter]->SetXTitle("#Delta#phi");
      h_met_deltaPhi[xs_counter]->SetYTitle("#slash{E}_{T}");
      h_met_deltaPhi[xs_counter]->SetMaximum(30);
      h_met_deltaPhi[xs_counter]->Draw("colz");
      gPad->SetLogy(0);
      ca->SaveAs("Plots/Dphi/"+met_dphi_name+".pdf");
      ca->SaveAs("Plots/Dphi/"+met_dphi_name+".png");
      
      h_r2_deltaPhi[xs_counter]->SetStats(0);
      h_r2_deltaPhi[xs_counter]->SetXTitle("#Delta#phi");
      h_r2_deltaPhi[xs_counter]->SetYTitle("R^{2}");
      h_r2_deltaPhi[xs_counter]->SetMaximum(30);
      h_r2_deltaPhi[xs_counter]->Draw("colz");
      ca->SaveAs("Plots/Dphi/"+r2_dphi_name+".pdf");
      ca->SaveAs("Plots/Dphi/"+r2_dphi_name+".png");
      
      
      h_Dphi[xs_counter]->SetStats(0);
      h_Dphi[xs_counter]->GetYaxis()->SetRangeUser(1,300);
      h_Dphi[xs_counter]->SetXTitle("#Delta#phi");
      h_Dphi[xs_counter]->Draw("");
      gPad->SetLogy(1);
      ca->SaveAs("Plots/Dphi/"+dmsampleTS+"_Dphi_1d.pdf");
      ca->SaveAs("Plots/Dphi/"+dmsampleTS+"_Dphi_1d.png");

      gPad->SetLogy(0);
      h_r2_alpha[xs_counter]->SetXTitle("#alpha");
      h_r2_alpha[xs_counter]->SetYTitle("R^{2}");
      h_r2_alpha[xs_counter]->SetStats(0);
      h_r2_alpha[xs_counter]->Draw("colz");
      gPad->SetLogy(0);
      ca->SaveAs("Plots/alpha/"+r2_alpha_name+".pdf");
      ca->SaveAs("Plots/alpha/"+r2_alpha_name+".png");
      
      h_r2_alphaT[xs_counter]->SetStats(0);
      h_r2_alphaT[xs_counter]->SetXTitle("#alpha_{T}");
      h_r2_alphaT[xs_counter]->SetYTitle("R^{2}");
      h_r2_alphaT[xs_counter]->Draw("colz");
      ca->SaveAs("Plots/alpha/"+r2_alphaT_name+".pdf");
      ca->SaveAs("Plots/alpha/"+r2_alphaT_name+".png");
      
      
      h_met_alpha[xs_counter]->SetStats(0);
      h_met_alpha[xs_counter]->SetXTitle("#alpha");
      h_met_alpha[xs_counter]->SetYTitle("#slash{E}_{T}");
      h_met_alpha[xs_counter]->Draw("colz");
      ca->SaveAs("Plots/alpha/"+met_alpha_name+".pdf");
      ca->SaveAs("Plots/alpha/"+met_alpha_name+".png");
      
      h_met_alphaT[xs_counter]->SetStats(0);
      h_met_alphaT[xs_counter]->SetXTitle("#alpha_{T}");
      h_met_alphaT[xs_counter]->SetYTitle("#slash{E}_{T}");
      h_met_alphaT[xs_counter]->Draw("colz");
      ca->SaveAs("Plots/alpha/"+met_alphaT_name+".pdf");
      ca->SaveAs("Plots/alpha/"+met_alphaT_name+".png");

      h_alpha[xs_counter]->SetStats(0);
      h_alpha[xs_counter]->SetXTitle("#alpha");
      h_alpha[xs_counter]->Draw();
      gPad->SetLogy(1);
      ca->SaveAs("Plots/alpha/"+dmsampleTS+"_alpha_1d.pdf");
      ca->SaveAs("Plots/alpha/"+dmsampleTS+"_alpha_1d.png");
      
      
      h_alphaT[xs_counter]->SetStats(0);
      h_alphaT[xs_counter]->SetXTitle("#alpha_{T}");
      h_alphaT[xs_counter]->Draw();
      ca->SaveAs("Plots/alpha/"+dmsampleTS+"_alphaT_1d.pdf");
      ca->SaveAs("Plots/alpha/"+dmsampleTS+"_alphaT_1d.png");
      
      h_r2[xs_counter]->SetStats(0);
      h_r2[xs_counter]->SetXTitle("R^{2}");
      h_r2[xs_counter]->Draw();
      ca->SaveAs("Plots/R2/"+dmsampleTS+"_r2_1d.pdf");
      ca->SaveAs("Plots/R2/"+dmsampleTS+"_r2_1d.png");
      
      h_met[xs_counter]->SetStats(0);
      h_met[xs_counter]->SetXTitle("#slash{E}_{T}");
      h_met[xs_counter]->Draw();
      ca->SaveAs("Plots/met/"+dmsampleTS+"_met_1d.pdf");
      ca->SaveAs("Plots/met/"+dmsampleTS+"_met_1d.png");
      
      h_theta[xs_counter]->SetStats(0);
      h_theta[xs_counter]->SetXTitle("#Delta#theta");
      h_theta[xs_counter]->Draw();
      ca->SaveAs("Plots/theta/"+dmsampleTS+"_theta_1d.pdf");
      ca->SaveAs("Plots/theta/"+dmsampleTS+"_theta_1d.png");
      
      h_R_DM[xs_counter]->SetStats(0);
      h_R_DM[xs_counter]->SetXTitle("R_{DM}");
      h_R_DM[xs_counter]->Draw();
      ca->SaveAs("Plots/Rdm/"+dmsampleTS+"_Rdm_1d.pdf");
      ca->SaveAs("Plots/Rdm/"+dmsampleTS+"_Rdm_1d.png");

      h_M2[xs_counter]->SetStats(0);
      h_M2[xs_counter]->SetXTitle("M2");
      h_M2[xs_counter]->Draw();
      ca->SaveAs("Plots/Rdm/"+dmsampleTS+"_M2_1d.pdf");
      ca->SaveAs("Plots/Rdm/"+dmsampleTS+"_M2_1d.png");

      h_M2_overPhi[xs_counter]->SetStats(0);
      h_M2_overPhi[xs_counter]->SetXTitle("M2");
      h_M2_overPhi[xs_counter]->Draw();
      ca->SaveAs("Plots/Rdm/"+dmsampleTS+"_M2_OverPhi_1d.pdf");
      ca->SaveAs("Plots/Rdm/"+dmsampleTS+"_M2_OverPhi_1d.png");

      /////////////////////////
      //Boosted Histos
      ////////////////////////
      
      h_beta_r[xs_counter]->SetStats(0);
      h_beta_r[xs_counter]->Draw();
      ca->SaveAs("Plots/Boost/"+dmsampleTS+"_beta_r_1d.pdf");
      ca->SaveAs("Plots/Boost/"+dmsampleTS+"_beta_r_1d.png");

      h_theta_B[xs_counter]->SetStats(0);
      h_theta_B[xs_counter]->SetXTitle("#Delta#theta_{B}");
      h_theta_B[xs_counter]->Draw();
      ca->SaveAs("Plots/theta/"+dmsampleTS+"_theta_1d_B.pdf");
      ca->SaveAs("Plots/theta/"+dmsampleTS+"_theta_1d_B.png");

      gPad->SetLogy(0);
      h_met_deltaPhi_B[xs_counter]->SetXTitle("#Delta#phi_{R}");
      h_met_deltaPhi_B[xs_counter]->SetYTitle("#slash{E}_{T}");
      h_met_deltaPhi_B[xs_counter]->SetStats(0);
      h_met_deltaPhi_B[xs_counter]->SetMaximum(30);
      h_met_deltaPhi_B[xs_counter]->Draw("colz");
      gPad->SetLogy(0);
      ca->SaveAs("Plots/Dphi/"+met_dphi_name+"_B.pdf");
      ca->SaveAs("Plots/Dphi/"+met_dphi_name+"_B.png");
      
      h_r2_deltaPhi_B[xs_counter]->SetXTitle("#Delta#phi_{R}");
      h_r2_deltaPhi_B[xs_counter]->SetYTitle("R^{2}");
      h_r2_deltaPhi_B[xs_counter]->Draw("colz");
      h_r2_deltaPhi_B[xs_counter]->SetMaximum(30);
      h_r2_deltaPhi_B[xs_counter]->SetStats(0);
      ca->SaveAs("Plots/Dphi/"+r2_dphi_name+"_B.pdf");
      ca->SaveAs("Plots/Dphi/"+r2_dphi_name+"_B.png");

      gPad->SetLogy();
      h_Dphi_B[xs_counter]->SetXTitle("#Delta#phi_{R}");
      h_Dphi_B[xs_counter]->SetStats(0);
      h_Dphi_B[xs_counter]->GetYaxis()->SetRangeUser(1,300);
      h_Dphi_B[xs_counter]->Draw("");
      gPad->SetLogy();
      ca->SaveAs("Plots/Dphi/"+dmsampleTS+"_Dphi_1d_B.pdf");
      ca->SaveAs("Plots/Dphi/"+dmsampleTS+"_Dphi_1d_B.png");
      
      gPad->SetLogy(0);
      h_r2_alpha_B[xs_counter]->SetStats(0);
      h_r2_alpha_B[xs_counter]->SetXTitle("#alpha^{B}");
      h_r2_alpha_B[xs_counter]->SetYTitle("R^{2}");
      h_r2_alpha_B[xs_counter]->Draw("colz");
      gPad->SetLogy(0);
      ca->SaveAs("Plots/alpha/"+r2_alpha_name+"_B.pdf");
      ca->SaveAs("Plots/alpha/"+r2_alpha_name+"_B.png");
      
      h_r2_alphaT_B[xs_counter]->SetStats(0);
      h_r2_alphaT_B[xs_counter]->SetXTitle("#alpha^{B}_{T}");
      h_r2_alphaT_B[xs_counter]->SetYTitle("R^{2}");
      h_r2_alphaT_B[xs_counter]->Draw("colz");
      ca->SaveAs("Plots/alpha/"+r2_alphaT_name+"_B.pdf");
      ca->SaveAs("Plots/alpha/"+r2_alphaT_name+"_B.png");
      
      h_met_alpha_B[xs_counter]->SetStats(0);
      h_met_alpha_B[xs_counter]->SetXTitle("#alpha^{B}");
      h_met_alpha_B[xs_counter]->SetYTitle("#slash{E}_{T}");
      h_met_alpha_B[xs_counter]->Draw("colz");
      ca->SaveAs("Plots/alpha/"+met_alpha_name+"_B.pdf");
      ca->SaveAs("Plots/alpha/"+met_alpha_name+"_B.png");
      
      h_met_alphaT_B[xs_counter]->SetStats(0);
      h_met_alphaT_B[xs_counter]->SetXTitle("#alpha^{B}_{T}");
      h_met_alphaT_B[xs_counter]->SetYTitle("#slash{E}_{T}");
      h_met_alphaT_B[xs_counter]->Draw("colz");
      ca->SaveAs("Plots/alpha/"+met_alphaT_name+"_B.pdf");
      ca->SaveAs("Plots/alpha/"+met_alphaT_name+"_B.png");

      gPad->SetLogy();
      h_alpha_B[xs_counter]->SetStats(0);
      h_alpha_B[xs_counter]->SetXTitle("#alpha^{B}");
      h_alpha_B[xs_counter]->Draw();
      gPad->SetLogy();
      ca->SaveAs("Plots/alpha/"+dmsampleTS+"_alpha_1d_B.pdf");
      ca->SaveAs("Plots/alpha/"+dmsampleTS+"_alpha_1d_B.png");
      
      h_alphaT_B[xs_counter]->SetStats(0);
      h_alphaT_B[xs_counter]->SetXTitle("#alpha^{B}_{T}");
      h_alphaT_B[xs_counter]->Draw();
      ca->SaveAs("Plots/alpha/"+dmsampleTS+"_alphaT_1d_B.pdf");
      ca->SaveAs("Plots/alpha/"+dmsampleTS+"_alphaT_1d_B.png");
      
      h_R_DM_B[xs_counter]->SetStats(0);
      h_R_DM_B[xs_counter]->SetXTitle("R_{DM}");
      h_R_DM_B[xs_counter]->Draw();
      ca->SaveAs("Plots/Rdm/"+dmsampleTS+"_Rdm_1d_B.pdf");
      ca->SaveAs("Plots/Rdm/"+dmsampleTS+"_Rdm_1d_B.png");
      

      double sample_eff = N_passed/Gen_Evts;
      
      std::cout << "Sample Eff: " << sample_eff*100 << "%" <<std::endl;
      outfile << dm_sample << " & " << sample_eff*100 << "\\%" << "\\" << "\\" << "\n";
      outfile << "\\hline" << std::endl;
      //outfile << dm_sample << " " << sample_eff*100 << "\n";
      
      TString s = h_2d[xs_counter]->GetName();
      s = s+"_combine.root";
      TString data_card_name = dm_sample.c_str();
      data_card_name = "combine/" + data_card_name + ".txt";
      std::ofstream data_card_f(data_card_name);
      data_card_f << "imax 1\njmax 4\nkmax 6\n";
      data_card_f << "------------------------------------------------------------------------------------------\n";
      data_card_f << "shapes * *\t" << s << "\t\t$PROCESS\t$PROCESS_$SYSTEMATIC\n";
      data_card_f << "------------------------------------------------------------------------------------------\n";
      data_card_f << "Observation\t" << data_obs->Integral() << "\n";
      data_card_f << "------------------------------------------------------------------------------------------\n";
      data_card_f << "bin\t\tb1\t\tb1\t\tb1\t\tb1\t\tb1\n";
      data_card_f << "process\t\tsignal_1D\ttt_1D\t\tdy_1D\t\tz_1D\t\tw_1D\n";
      data_card_f << "process\t\t0\t\t1\t\t2\t\t3\t\t4\n";
      data_card_f << "rate\t\t"<< h_2d[xs_counter]->Integral() <<"\t\t"<< tt_N <<"\t\t" << dy_N << "\t\t" << z_N << "\t\t" << w_N << "\n";
      data_card_f << "------------------------------------------------------------------------------------------\n";
      data_card_f << "lumi\tlnN\t1.026\t\t1.0\t\t1.0\t\t1.0\t\t1.0\n";
      data_card_f << "alpha\tshape\t1\t\t-\t\t-\t\t-\t\t-\n";
      data_card_f << "beta\tshape\t-\t\t1\t\t-\t\t-\t\t-\n";
      data_card_f << "gamma\tshape\t-\t\t-\t\t1\t\t-\t\t-\n";
      data_card_f << "delta\tshape\t-\t\t-\t\t-\t\t1\t\t-\n";
      data_card_f << "zeta\tshape\t-\t\t-\t\t-\t\t-\t\t1\n";
      data_card_f.close();

      /////1D files
      TString s1 = h_2d[xs_counter]->GetName();
      s1 = s1+"_combine_rsq.root";
      TString data_card_name1 = dm_sample.c_str();
      data_card_name1 = "combine/" + data_card_name1 + "_rsq.txt";
      std::ofstream data_card_f1(data_card_name1);
      data_card_f1 << "imax 1\njmax 4\nkmax 6\n";
      data_card_f1 << "------------------------------------------------------------------------------------------\n";
      data_card_f1 << "shapes * *\t" << s1 << "\t\t$PROCESS\t$PROCESS_$SYSTEMATIC\n";
      data_card_f1 << "------------------------------------------------------------------------------------------\n";
      data_card_f1 << "Observation\t" << data_rsq->Integral() << "\n";
      data_card_f1 << "------------------------------------------------------------------------------------------\n";
      data_card_f1 << "bin\t\tb1\t\tb1\t\tb1\t\tb1\t\tb1\n";
      data_card_f1 << "process\t\tsignal_rsq\ttt_rsq\t\tdy_rsq\t\tz_rsq\t\tw_rsq\n";
      data_card_f1 << "process\t\t0\t\t1\t\t2\t\t3\t\t4\n";
      data_card_f1 << "rate\t\t"<< h_rsq[xs_counter]->Integral() <<"\t\t"<< tt_N1 <<"\t\t" << dy_N1 << "\t\t" << z_N1 << "\t\t" << w_N1 << "\n";
      data_card_f1 << "------------------------------------------------------------------------------------------\n";
      data_card_f1 << "lumi\tlnN\t1.026\t\t1.0\t\t1.0\t\t1.0\t\t1.0\n";
      data_card_f1 << "alpha\tshape\t1\t\t-\t\t-\t\t-\t\t-\n";
      data_card_f1 << "beta\tshape\t-\t\t1\t\t-\t\t-\t\t-\n";
      data_card_f1 << "gamma\tshape\t-\t\t-\t\t1\t\t-\t\t-\n";
      data_card_f1 << "delta\tshape\t-\t\t-\t\t-\t\t1\t\t-\n";
      data_card_f1 << "zeta\tshape\t-\t\t-\t\t-\t\t-\t\t1\n";
      data_card_f1.close();

      
      xs_counter++;
      delete out;
      delete eff;
      delete mr_rsq;
    }
    
  }else{
    std::cout << "Unable to open the file" << std::endl;
  }
  mfile0.close();
  outfile << "\\end{tabular}\n\\end{center}\n\\label{default}\n\\end{table}\n";
  outfile.close();
  
  
  TFile* n_f[24];
  
  for(int k = 0; k < 24; k++){
    TString s = h_2d[k]->GetName();
    s = s+"_combine.root";
    TString s2 = "combine/";
    s = s2+s;
    std::cout << s << std::endl;
    n_f[k] = new TFile(s, "RECREATE");
    h_2d[k]->Write("signal_2D");
    data_obs->Write("data_obs");
    h_1d[k]->Write("signal_1D");
    h_1d_up[k]->Write("signal_1D_alphaUp");
    h_1d_down[k]->Write("signal_1D_alphaDown");
    
    tt_1D->Write();
    tt_1D_alphaUp->Write("tt_1D_betaUp");
    tt_1D_alphaDown->Write("tt_1D_betaDown");
    
    dy_1D->Write();
    dy_1D_alphaUp->Write("dy_1D_gammaUp");
    dy_1D_alphaDown->Write("dy_1D_gammaDown");
      
    z_1D->Write();
    z_1D_alphaUp->Write("z_1D_deltaUp");
    z_1D_alphaDown->Write("z_1D_deltaDown");
    
    w_1D->Write();
    w_1D_alphaUp->Write("w_1D_zetaUp");
    w_1D_alphaDown->Write("w_1D_zetaDown");
    n_f[k]->Close();
  }

  TFile* n_f_rsq[24];
  
  for(int k = 0; k < 24; k++){
    TString s = h_2d[k]->GetName();
    s = s+"_combine_rsq.root";
    TString s2 = "combine/";
    s = s2+s;
    std::cout << s << std::endl;
    n_f_rsq[k] = new TFile(s, "RECREATE");
    data_rsq->Write("data_obs");
    h_rsq[k]->Write("signal_rsq");
    h_rsq_up[k]->Write("signal_rsq_alphaUp");
    h_rsq_down[k]->Write("signal_rsq_alphaDown");
    
    tt_rsq->Write();
    tt_rsq_alphaUp->Write("tt_rsq_betaUp");
    tt_rsq_alphaDown->Write("tt_rsq_betaDown");
    
    dy_rsq->Write();
    dy_rsq_alphaUp->Write("dy_rsq_gammaUp");
    dy_rsq_alphaDown->Write("dy_rsq_gammaDown");
      
    z_rsq->Write();
    z_rsq_alphaUp->Write("z_rsq_deltaUp");
    z_rsq_alphaDown->Write("z_rsq_deltaDown");
    
    w_rsq->Write();
    w_rsq_alphaUp->Write("w_rsq_zetaUp");
    w_rsq_alphaDown->Write("w_rsq_zetaDown");
    n_f_rsq[k]->Close();
  }
  
  return 0;
}
