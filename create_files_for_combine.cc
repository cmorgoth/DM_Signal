{
  gROOT->Reset();
  
  TH1F *bkg_1D = new TH1F("bkg_1D","bkg_1D",16,0,16);
  TH1F *bkg_1D_alphaUp = new TH1F("bkg_1D_alphaUp","bkg_1D_alphaUp",16,0,16);
  TH1F *bkg_1D_alphaDown = new TH1F("bkg_1D_alphaDown","bkg_1D_alphaDown",16,0,16);
  TH1F *data_1D_unwrapt = new TH1F("data_unwrapt","data_Unwrapt_16bins",16,0,16);

  
  TH1F *tt_1D = new TH1F("tt_1D","tt_1D",16,0,16);
  TH1F *tt_1D_alphaUp = new TH1F("tt_1D_alphaUp","tt_1D_alphaUp",16,0,16);
  TH1F *tt_1D_alphaDown = new TH1F("tt_1D_alphaDown","tt_1D_alphaDown",16,0,16);

  TH1F *dy_1D = new TH1F("dy_1D","dy_1D",16,0,16);
  TH1F *dy_1D_alphaUp = new TH1F("dy_1D_alphaUp","dy_1D_alphaUp",16,0,16);
  TH1F *dy_1D_alphaDown = new TH1F("dy_1D_alphaDown","dy_1D_alphaDown",16,0,16);

  TH1F *z_1D = new TH1F("z_1D","z_1D",16,0,16);
  TH1F *z_1D_alphaUp = new TH1F("z_1D_alphaUp","z_1D_alphaUp",16,0,16);
  TH1F *z_1D_alphaDown = new TH1F("z_1D_alphaDown","z_1D_alphaDown",16,0,16);

  TH1F *w_1D = new TH1F("w_1D","w_1D",16,0,16);
  TH1F *w_1D_alphaUp = new TH1F("w_1D_alphaUp","w_1D_alphaUp",16,0,16);
  TH1F *w_1D_alphaDown = new TH1F("w_1D_alphaDown","w_1D_alphaDown",16,0,16);
  

  int counter = 0;
  for ( int i=1; i<5; i++)
  {
  for (int j=1; j<5; j++)
  	{
     counter++;
      double results = bkg->GetBinContent(i,j);
      double results2 = bkg_alphaUp->GetBinContent(i,j);
      double results3 = bkg_alphaDown->GetBinContent(i,j);
      double results_data = data_0b_0mu->GetBinContent(i,j);
      //ttbar
      double tt_c = tt_0b_0mu->GetBinContent(i,j);
      double tt_e = tt_0b_0mu->GetBinError(i,j);
      double tt_u = tt_0b_0mu->GetBinErrorUp(tt_0b_0mu->GetBin(i,j));
      double tt_l = tt_0b_0mu->GetBinErrorLow(tt_0b_0mu->GetBin(i,j));
      //dy
      double dy_c = p_0b_0mu_dy->GetBinContent(i,j);
      double dy_e = p_0b_0mu_dy->GetBinError(i,j);
      double dy_u = p_0b_0mu_dy->GetBinErrorUp(p_0b_0mu_dy->GetBin(i,j));
      double dy_l = p_0b_0mu_dy->GetBinErrorLow(p_0b_0mu_dy->GetBin(i,j));
      //z
      double z_c = p_0b_0mu_Z->GetBinContent(i,j);
      double z_e = p_0b_0mu_Z->GetBinError(i,j);
      double z_u = p_0b_0mu_Z->GetBinErrorUp(p_0b_0mu_Z->GetBin(i,j));
      double z_l = p_0b_0mu_Z->GetBinErrorLow(p_0b_0mu_Z->GetBin(i,j));
      //w
      double w_c = p_0b_0mu_W->GetBinContent(i,j);
      double w_e = p_0b_0mu_W->GetBinError(i,j);
      double w_u = p_0b_0mu_W->GetBinErrorUp(p_0b_0mu_W->GetBin(i,j));
      double w_l = p_0b_0mu_W->GetBinErrorLow(p_0b_0mu_W->GetBin(i,j));
      std::cout << "error: " << w_e << " w_l: " << w_l << " w_u: " << w_l << " w_c-w_l: " << w_c-w_l << std::endl;
      std::cout << "error: " << z_e << " z_l: " << z_l << " z_u: " << z_l << " z_c-z_l: " << z_c-z_l << std::endl;
      std::cout << "error: " << dy_e << " dy_l: " << dy_l << " dy_u: " << dy_l << " dy_c-dy_l: " << dy_c-dy_l << std::endl;
      std::cout << "error: " << tt_e << " tt_l: " << tt_l << " tt_u: " << tt_l << " tt_c-tt_l: " << tt_c-tt_l << std::endl;
      //unwrapping 2D into 1D
      //tt
      tt_1D->SetBinContent(counter,tt_c);
      tt_1D_alphaUp->SetBinContent(counter,tt_c+tt_u);
      tt_1D_alphaDown->SetBinContent(counter,tt_c-tt_l);
      //dy
      dy_1D->SetBinContent(counter,dy_c);
      dy_1D_alphaUp->SetBinContent(counter,dy_c+dy_u);
      dy_1D_alphaDown->SetBinContent(counter,dy_c-dy_l);
      //z
      z_1D->SetBinContent(counter,z_c);
      z_1D_alphaUp->SetBinContent(counter,z_c+z_u);
      z_1D_alphaDown->SetBinContent(counter,z_c-z_l);
      //w
      w_1D->SetBinContent(counter,w_c);
      w_1D_alphaUp->SetBinContent(counter,w_c+w_u);
      w_1D_alphaDown->SetBinContent(counter,w_c-w_l);
      
      data_1D_unwrapt->SetBinContent(counter,results_data);
      bkg_1D->SetBinContent(counter,results);
      bkg_1D_alphaUp->SetBinContent(counter,results2);
      bkg_1D_alphaDown->SetBinContent(counter,results3);
    }
  }
  //bkg_file_1DRsq->Write("BkgPred_Unwrapt");
  //data_1D_unwrapt->Write("data_Unwrapt");
  //bkg_1D->Write("BkgPred_Unwrapt");
  //bkg_1D_alphaUp->Write("BkgPred_Unwrapt_alphaUp");
  //bkg_1D_alphaDown->Write("BkgPred_Unwrapt_alphaDown");
  
  tt_1D->Write();
  tt_1D_alphaUp->Write();
  tt_1D_alphaDown->Write();

  dy_1D->Write();
  dy_1D_alphaUp->Write();
  dy_1D_alphaDown->Write();

  z_1D->Write();
  z_1D_alphaUp->Write();
  z_1D_alphaDown->Write();

  w_1D->Write();
  w_1D_alphaUp->Write();
  w_1D_alphaDown->Write();
}
