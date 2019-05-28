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

double epjets=140000;
double mupjets=300000;
double eptw=352000;
double mutw=600000;
double epmulti=56000;
double mupmulti=146500;
double ept=74000;
double mut=143000;


TFile* File = TFile::Open("/home/sane/Stage_M1/LIVMass-master/results/Modulation_Temporelle_13TeV_b5.root");
TFile* tchannel = new TFile ("hist_tchannel.root","update");
TFile* tbarchannel = new TFile ("hist_tbarchannel.root","update");

auto legend = new TLegend(0.70,0.75,0.98,0.95);


  TCanvas* Signalt = new TCanvas ("", "", 450, 450);
  TGraph* Parti = (TGraph*)File->Get("tX");
  TGraph* AntiParti = (TGraph*)File->Get("tbarX");
  double* Particle = Parti->GetY();
  double* AntiParticle = AntiParti->GetY();



  TH1F* h1sigt = new TH1F ( "", "", 24 , 0, 24);
  TH1F* h2bggammat = new TH1F ( "", "", 24 , 0, 24);
  TH1F* h3bgttbart = new TH1F ("","",24 , 0, 24);
  TH1F* h4bgMultit = new TH1F ("","", 24, 0, 24);
  TH1F* h5ftsignt = new TH1F ("","",24 ,0,24);
  TH1F* h6errort = new TH1F ("","",24,0,24);

  h1sigt->SetStats(kFALSE);
  h2bggammat->SetStats(kFALSE);
  h3bgttbart->SetStats(kFALSE);
  h4bgMultit->SetStats(kFALSE);
  h5ftsignt->SetStats(kFALSE);
  h6errort->SetStats(kFALSE);
  h1sigt->GetXaxis()->SetTitle("t (h)");
  h1sigt->GetYaxis()->SetTitle("Events");
  h1sigt->GetXaxis()->CenterTitle(kTRUE);
  h1sigt->GetYaxis()->CenterTitle(kTRUE);


  double sigmatott = sqrt((sqrt(epjets)*sqrt(epjets)+sqrt(mupjets)*sqrt(mupjets)+sqrt(eptw)*sqrt(eptw)+sqrt(mutw)*sqrt(mutw)+sqrt(epmulti)*sqrt(epmulti)+sqrt(mupmulti)*sqrt(mupmulti)+sqrt(ept)*sqrt(ept)+sqrt(mut)*sqrt(mut))/24);
  int k=0;
  for (double i = 1 ; i<=24 ; i++)
  {
    h1sigt->SetBinContent(i,8448+18298+39666+9040);
    h3bgttbart->SetBinContent(i,8448+18298+39666);
    h4bgMultit->SetBinContent(i,8448);
    h2bggammat->SetBinContent(i,8448+18298);
    h5ftsignt->SetBinContent(i, 8448+18298+39666+9040+9040*Particle[k]);
    h6errort->SetBinContent(i, 8448+18298+39666+9040);
    h6errort->SetBinError(i, sigmatott);
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
  legend->AddEntry(h5ftsignt,"bx = 5 GeV","l");
  h5ftsignt->GetYaxis()->SetRangeUser(0,120000);



  h1sigt->Draw("SAME");
  h3bgttbart->Draw("SAME");
  h2bggammat->Draw("SAME");
  h4bgMultit->Draw("SAME");
  h6errort->Draw("SAME");
  h5ftsignt->Draw("SAME");
  legend->Draw("SAME");


  Signalt->SaveAs("Histo_t.pdf");


//______________Save t_channel to rootFile____________//
  tchannel->cd();
  h1sigt->Write("Sum_background");
  k=0;
  for (double i = 1 ; i<=24 ; i++)
  {
    h1sigt->SetBinContent(i,9040);
    h3bgttbart->SetBinContent(i,39666);
    h4bgMultit->SetBinContent(i,8448);
    h2bggammat->SetBinContent(i,18298);
    h5ftsignt->SetBinContent(i, 9040*Particle[k]);
    k+=2;
  }
  h1sigt->Write("t_channel");
  h3bgttbart->Write("ttbar/tw");
  h4bgMultit->Write("Multijet");
  h2bggammat->Write("Z/W/gamma_jets");
  h5ftsignt->Write("Modulation_t");



//_______________tbar channel sensibility plot_________________//


  TCanvas* Signaltbar = new TCanvas("","", 450, 450);
  TH1F* h1sigtbar = new TH1F ( "", "", 24 , 0, 24);
  TH1F* h2bggammatbar = new TH1F ( "", "", 24 , 0, 24);
  TH1F* h3bgttbartbar = new TH1F ("","",24 , 0, 24);
  TH1F* h4bgMultitbar = new TH1F ("","", 24 , 0, 24);
  TH1F* h5ftsigntbar = new TH1F ("","", 24,0 ,24);
  TH1F* h6errortbar = new TH1F ("","",24,0,24);
  TH2F* testsignal = new TH2F ( "","", 24,0,24, 200,-50,+50);
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

  double sigmatottbar = sqrt(8416+16250+39708+5750);
  k=0;
  for (double i = 1 ; i<=24 ; i++)
  {
    h1sigtbar->SetBinContent(i,8416+16250+39708+5750);
    h3bgttbartbar->SetBinContent(i,8416+16250+39708);
    h4bgMultitbar->SetBinContent(i,8416);
    h2bggammatbar->SetBinContent(i,8416+16250);
    h5ftsigntbar->SetBinContent(i,8416+16250+39708+5750+5750*AntiParticle[k]);
    h6errortbar->SetBinContent( i,8416+16250+39708+5750);
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
  legendtbar->AddEntry(h5ftsignt,"bx = 5 GeV","l");

  h5ftsigntbar->GetYaxis()->SetRangeUser(0,120000);


  h1sigtbar->Draw("SAME");

  h3bgttbartbar->Draw("SAME");

  h2bggammatbar->Draw("SAME");
  h4bgMultitbar->Draw("SAME");
  h6errortbar->Draw("SAME");
  h5ftsigntbar->Draw("SAME");
  legendtbar->Draw("SAME");
  Signaltbar->SaveAs("Histo_tbar.pdf");



//______________Save tbar_channel to rootFile_____________//
  tbarchannel->cd();
  h1sigtbar->Write("Sum_background");
  k=0;
  for (double i = 1 ; i<=24 ; i++)
  {
    h1sigtbar->SetBinContent(i,5750);
    h3bgttbartbar->SetBinContent(i,39708);
    h4bgMultitbar->SetBinContent(i,8416);
    h2bggammatbar->SetBinContent(i,16250);
    testsignal->SetBinContent(i, 5750*AntiParticle[k]);
    cout<<testsignal->GetBinContent(i)<<endl;
    k+=2;
  }
  h1sigtbar->Write("tbar_channel");
  h3bgttbartbar->Write("ttbar/tw");
  h4bgMultitbar->Write("Multijet");
  h2bggammatbar->Write("Z/W/gamma_jets");
  h5ftsignt->Write("Modulation_tbar");

  return 0;

}
