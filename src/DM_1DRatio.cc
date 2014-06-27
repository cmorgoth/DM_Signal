#include "DM_1DRatio.hh"

int StackSignal(THStack* stack, TH1F* h1, TH1F* h2, TString h1Name = "h1Name", TString h2Name = "h2Name", TString fname = "default_name", TString type = "defaultType", TLegend* leg = NULL){
  
  TCanvas* C = new TCanvas("C", "C	", 400, 500);
  C->cd();  
  
  TString label;
  
  TPad *pad1 = new TPad("pad1","pad1",0,0,1,1);
  pad1->SetBottomMargin(0.1);
  pad1->SetLeftMargin(0.13);
  pad1->Draw();
  pad1->cd();
  h1->SetStats(0);
  h1->SetTitle("");
  //h2->SetMarkerSize(0.6);
  h2->SetLineWidth(2);
  h1->SetLineWidth(2);
  stack->SetMinimum(0.1*h2->GetBinContent(h2->GetNbinsX()));
  stack->SetMaximum(4*h2->GetBinContent(1));
  
  stack->Draw();
  stack->GetHistogram()->GetXaxis()->CenterTitle();
  stack->GetHistogram()->GetXaxis()->SetTitleOffset(1.2);
  stack->GetHistogram()->GetXaxis()->SetTitle("R^{2}");
  stack->GetHistogram()->GetYaxis()->CenterTitle();
  stack->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
  stack->GetHistogram()->GetYaxis()->SetTitle("Events");
  h2->SetXTitle("Events");
  
  h1->Draw("same");
  h2->Draw("pesame");
  h1->Draw("same");
  pad1->Update();
  C->cd();
  
  pad1->SetLogy();
  C->Update();

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextSize(0.03);
  t->DrawLatex(0.22,0.95,"CMS Preliminary:");
  t->DrawLatex(0.42,0.95,"#sqrt{s} = 8 TeV,");
  t->DrawLatex(0.62,0.95,"#int L dt = 18.5 fb^{-1}");
  leg->SetTextSize(.022);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->Draw();
  C->cd();
  C->SaveAs(fname + ".C");
  C->SaveAs(fname + ".root");
  C->SaveAs(fname + ".pdf");
  C->SaveAs(fname + ".png");
  
  delete leg;
  delete C;
    
  return 0;
  
};
