#include  "../include/root_simu.hpp"
#include  "../include/matrixSingleTop.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH2F.h>


using namespace std;


int main (){

TFile* File = TFile::Open("/home/sane/Stage_M1/LIVMass-master/results/Modulation_Temporelle_13TeV_b1.root");
TFile* tchannel = new TFile ("/home/sane/Stage_M1/LIVMass-master/results/histo/hist_tchannel.root","RECREATE");


  //auto legend = new TLegend(0.70,0.75,0.98,0.95);

  //TCanvas* Signalt = new TCanvas ("", "", 450, 450);
  TGraph* Parti = (TGraph*)File->Get("tX");
  TGraph* AntiParti = (TGraph*)File->Get("tbarX");
  double* Particle = Parti->GetY();
  double* AntiParticle = AntiParti->GetY();



  /*TH1F* h1sigt = new TH1F ( "", "", 24 , 0, 24);
  TH1F* h2bggammat = new TH1F ( "", "", 24 , 0, 24);
  TH1F* h3bgttbart = new TH1F ("","",24 , 0, 24);
  TH1F* h4bgMultit = new TH1F ("","", 24, 0, 24);
  TH1F* h5ftsignt = new TH1F ("","",24 ,0,24);
  TH1F* h6errort = new TH1F ("","",24,0,24);*/

  //h1sigt->SetStats(kFALSE);
  //h2bggammat->SetStats(kFALSE);
  //h3bgttbart->SetStats(kFALSE);
  //h4bgMultit->SetStats(kFALSE);
  //h5ftsignt->SetStats(kFALSE);
  //h6errort->SetStats(kFALSE);
  /*h1sigt->GetXaxis()->SetTitle("t (h)");
  h1sigt->GetYaxis()->SetTitle("Events");
  h1sigt->GetXaxis()->CenterTitle(kTRUE);
  h1sigt->GetYaxis()->CenterTitle(kTRUE);
  h1sigt->GetYaxis()->SetRangeUser(70000,80000);*/



  double sigmatott = sqrt(8446+18298+39392+9048);
  cout<<sigmatott<<endl;
  int k=0;

  /*for (double i = 1 ; i<=24 ; i++)
  {
    h1sigt->SetBinContent(i,8446+18298+39392+9048);
    h3bgttbart->SetBinContent(i,8446+18298+39392);
    h4bgMultit->SetBinContent(i,8446);
    h2bggammat->SetBinContent(i,8446+18298);
    h5ftsignt->SetBinContent(i, 8446+18298+39392+9048+9048*Particle[k]*6.14);
    //h1sigt->SetBinError(i,sigmatott);
    h6errort->SetBinError(i, sigmatott);
    h6errort->SetBinContent(i, 8446+18298+39392+9048);
    k+=2;
  }
  h1sigt->SetFillColor(kRed+1);
  h2bggammat->SetFillColor(kBlue+1);
  h3bgttbart->SetFillColor(kGreen+2);
  h4bgMultit->SetFillColor(kOrange-3);
  h5ftsignt->SetLineColor(kBlack);
  legend->AddEntry(h1sigt,"t channel", "f");
  legend->AddEntry(h2bggammat, "W/Z/#gamma* + jets", "f");
  legend->AddEntry(h3bgttbart,"t#bar{t}/tW", "f");
  legend->AddEntry(h4bgMultit,"Multijet", "f");
  legend->AddEntry(h5ftsignt,"bx = 6.14 GeV","l");
  h5ftsignt->GetYaxis()->SetRangeUser(0,120000);



  h1sigt->Draw("SAME");
  h3bgttbart->Draw("SAME");
  h2bggammat->Draw("SAME");
  h4bgMultit->Draw("SAME");
  h6errort->Draw("SAME");
  h5ftsignt->Draw("SAME");
  legend->Draw("SAME");


  Signalt->SaveAs("Histo_t.pdf");*/


//______________Save t_channel to rootFile____________//

  TH1F* h1sigthist = new TH1F ( "", "", 24 , 0, 24);
  TH1F* h2bggammathist = new TH1F ( "", "", 24 , 0, 24);
  TH1F* h3bgttbarthist = new TH1F ("","",24 , 0, 24);
  TH1F* h4bgMultithist = new TH1F ("","", 24, 0, 24);
  TH1F* h5ftsignthist = new TH1F ("","",24 ,0,24);
  TH1F* asimovt = new TH1F ("","",24,0,24);
  //TH1F* h0background = new TH1F ("","", 24 , 0 , 24);
  /*TH1F* h1sigterror = new TH1F ( "", "", 24 , 0, 24);
  TH1F* h2bggammaterror = new TH1F ( "", "", 24 , 0, 24);
  TH1F* h3bgttbarterror = new TH1F ("","",24 , 0, 24);
  TH1F* h4bgMultiterror = new TH1F ("","", 24, 0, 24);
  TH1F* h0backgrounderror = new TH1F ("","", 24 , 0 , 24);*/

  /*double error1 = sqrt(9048);
  double error2 = sqrt (18298);
  double error3 = sqrt (39392);
  double error4 = sqrt (8446);
  double error0 = sqrt(18298+8446+39392);*/

  tchannel->cd();

  k=0;
  for (double i = 1 ; i<=24 ; i++)
  {
    asimovt->SetBinContent(i,9048+18298+39392+8446);
    h1sigthist->SetBinContent(i,9048);
    h2bggammathist->SetBinContent(i,18298);
    h3bgttbarthist->SetBinContent(i,39392);
    h4bgMultithist->SetBinContent(i,8446);
    h5ftsignthist->SetBinContent(i,Particle[k]);
    /*h0background->SetBinContent(i,18298+8446+39392);
    h1sigthist->SetBinError(i,error1);
    h3bgttbarthist->SetBinError(i,error3);
    h4bgMultithist->SetBinError(i,error4);
    h2bggammathist->SetBinError(i,error2);

    h1sigterror->SetBinError(i,error1);
    h0background->SetBinError(i, error0);
    h3bgttbarterror->SetBinError(i,error3);
    h4bgMultiterror->SetBinError(i,error4);
    h2bggammaterror->SetBinError(i,error2);*/

    k+=2;
  }
  asimovt->Write("Asimov");
  h1sigthist->Write("t_channel");
  h3bgttbarthist->Write("ttbar");
  h4bgMultithist->Write("Multijet");
  h2bggammathist->Write("Electroweak");
  h5ftsignthist->Write("Modulation_t");
  //h0background->Write("background");
  /*h1sigterror->Write("t_channel_error");
  h3bgttbarterror->Write("ttbar_error");
  h4bgMultiterror->Write("Multijet_error");
  h2bggammaterror->Write("Electroweak_error");*/

  tchannel->Close();


//_______________tbar channel sensibility plot_________________//
  TFile* tbarchannel = new TFile ("/home/sane/Stage_M1/LIVMass-master/results/histo/hist_tbarchannel.root","RECREATE");

  /*TCanvas* Signaltbar = new TCanvas("","", 450, 450);
  TH1F* h1sigtbar = new TH1F ( "", "", 24 , 0, 24);
  TH1F* h2bggammatbar = new TH1F ( "", "", 24 , 0, 24);
  TH1F* h3bgttbartbar = new TH1F ("","",24 , 0, 24);
  TH1F* h4bgMultitbar = new TH1F ("","", 24 , 0, 24);
  TH1F* h5ftsigntbar = new TH1F ("","", 24,0 ,24);
  TH1F* h6errortbar = new TH1F ("","",24,0,24);
  auto legendtbar = new TLegend(0.70,0.75,0.98,0.95);


  h1sigtbar->SetStats(kFALSE);
  h2bggammatbar->SetStats(kFALSE);
  h3bgttbartbar->SetStats(kFALSE);
  h4bgMultitbar->SetStats(kFALSE);
  h5ftsigntbar->SetStats(kFALSE);
  h6errortbar->SetStats(kFALSE);

  h1sigtbar->GetXaxis()->SetTitle("t (h)");
  h1sigtbar->GetYaxis()->SetTitle("Events");
  h1sigtbar->GetXaxis()->CenterTitle(kTRUE);
  h1sigtbar->GetYaxis()->CenterTitle(kTRUE);

  double sigmatottbar = sqrt(8404+16233+39618+5739);
  k=0;
  for (double i = 1 ; i<=24 ; i++)
  {
    h1sigtbar->SetBinContent(i,8404+16233+39618+5739);
    h3bgttbartbar->SetBinContent(i,8404+16233+39618);
    h4bgMultitbar->SetBinContent(i,8404);
    h2bggammatbar->SetBinContent(i,8404+16233);
    h5ftsigntbar->SetBinContent(i,8404+16233+39618+5739+5739*AntiParticle[k]*10);
    h6errortbar->SetBinContent( i,8404+16233+39618+5739);
    h6errortbar->SetBinError(i,sigmatottbar);
    k+=2;
  }
  h1sigtbar->SetFillColor(kRed+1);
  h2bggammatbar->SetFillColor(kBlue+1);
  h3bgttbartbar->SetFillColor(kGreen+2);
  h4bgMultitbar->SetFillColor(kOrange-3);
  h5ftsigntbar->SetLineColor(kBlack);

  legendtbar->AddEntry(h1sigt,"#bar{t} channel", "f");
  legendtbar->AddEntry(h2bggammat, "W/Z/#gamma* + jets", "f");
  legendtbar->AddEntry(h3bgttbart,"t#bar{t}/tW", "f");
  legendtbar->AddEntry(h4bgMultit,"Multijet", "f");
  legendtbar->AddEntry(h5ftsignt,"bx = 10 GeV","l");

  h5ftsigntbar->GetYaxis()->SetRangeUser(0,120000);


  h1sigtbar->Draw("SAME");

  h3bgttbartbar->Draw("SAME");

  h2bggammatbar->Draw("SAME");
  h4bgMultitbar->Draw("SAME");
  h6errortbar->Draw("SAME");
  h5ftsigntbar->Draw("SAME");
  legendtbar->Draw("SAME");
  Signaltbar->SaveAs("Histo_tbar.pdf");*/



//______________Save tbar_channel to rootFile_____________//
TH1F* h1sigthistbar = new TH1F ( "", "", 24 , 0, 24);
TH1F* h2bggammathistbar = new TH1F ( "", "", 24 , 0, 24);
TH1F* h3bgttbarthistbar = new TH1F ("","",24 , 0, 24);
TH1F* h4bgMultithistbar = new TH1F ("","", 24, 0, 24);
TH1F* h5ftsignthistbar = new TH1F ("","",24 ,0,24);
TH1F* asimovtbar = new TH1F ("","",24 ,0,24);

  tbarchannel->cd();

  k=0;
  for (double i = 1 ; i<=24 ; i++)
  {
    asimovtbar->SetBinContent(i,5739+39618+8404+16233);
    h1sigthistbar->SetBinContent(i,5739);
    h3bgttbarthistbar->SetBinContent(i,39618);
    h4bgMultithistbar->SetBinContent(i,8404);
    h2bggammathistbar->SetBinContent(i,16233);
    h5ftsignthistbar->SetBinContent(i, AntiParticle[k]);
    k+=2;
  }
  asimovtbar->Write("Asimov");
  h1sigthistbar->Write("tbar_channel");
  h3bgttbarthistbar->Write("ttbar");
  h4bgMultithistbar->Write("Multijet");
  h2bggammathistbar->Write("Electroweak");
  h5ftsignthistbar->Write("Modulation_tbar");
  tbarchannel->Close();

  return 0;

}
