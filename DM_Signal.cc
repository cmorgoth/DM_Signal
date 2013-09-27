
{

  gROOT->Reset();

  TFile* f2 = new TFile("trigger/hlt_eff_SignleElePD_Final.root");
  TEfficiency* hlt = (TEfficiency*)f2->Get("Eff2d");
    
  
  std::cout.precision(16);
  TFile* f = new TFile("data/DMm100AVu_testMC_0.root");
  TTree* eff = (TTree*)f->Get("effTree");
  TTree* out = (TTree*)f->Get("outTree");

  double mr[4], rsq[4], Npassed_In;
  int btag, box;

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
  out->SetBranchAddress("MR", mr);
  out->SetBranchAddress("RSQ", rsq);
  out->SetBranchAddress("nBtag", &btag);
  out->SetBranchAddress("BOX_NUM", &box);
  
  int N_out = out->GetEntries();
  
  int n_mr = 4;
  int n_rsq = 4;
  double mr_bin[] = {200, 400, 650, 1000, 3500};
  double rsq_bin[] = {0.0, 0.65, 0.8, 0.9, 2.5};
  
  TH2F* mr_rsq = new TH2F("mr_rsq", "mr_rsq", n_mr, mr_bin, n_rsq, rsq_bin);

  double Lumi = 19.6;
  double xsec = 1;
  double scaleF = Lumi*xsec/Gen_Evts;
  
  std::cout << "scaleF : " << scaleF << std::endl;
  
  for(int j = 0; j < N_out; j++){
    out->GetEntry(j);
    if(mr[2] > 200. && rsq[2] > 0.50 && btag == 0 && box == 0 )mr_rsq->Fill(mr[2], rsq[2]);
  }
  
  double sample_eff = mr_rsq->Integral()/Gen_Evts;

  std::cout << "Sample Eff: " << sample_eff*100 << "%" <<std::endl;

  TFile* n_f = new TFile("newFile.root", "RECREATE");
  
  mr_rsq->Write();
  n_f->Write();
  
}
