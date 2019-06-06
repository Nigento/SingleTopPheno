#include  "../include/root_simu.hpp"
#include  "../include/matrixSingleTop.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1F.h>
#include <TLatex.h>

using namespace std;
int main ()
{

  TFile* Filet = TFile::Open("/home/sane/Stage_M1/LIVMass-master/results/histo/hist_tchannel_Xcarre.root");
  TCanvas* feuille = new TCanvas ("", "", 450, 450);
  TFile* Filetbar = TFile::Open("/home/sane/Stage_M1/LIVMass-master/results/histo/hist_tbarchannel_Xcarre.root");
  //TFile* ft = TFile::Open("/home/sane/Stage_M1/LIVMass-master/results/Modulation_Temporelle_13TeV.root");

  feuille->SetLeftMargin(0.14);
  feuille->SetRightMargin(0.03);
  auto legend = new TLegend(0.65,0.75,0.90,0.85);

  double range = 7;
  double systttbar = 1.0;
  double systmulti = 0.974;
  double systsingletop = 1.0;
  double systelectro = 1.0;
  double lumi = 1.0 ;

  TH1F* datt = (TH1F*)Filet->Get("Asimov");
  TH1F* Sigt = (TH1F*)Filet->Get("Modulation_t");
  TH1F* ttbart = (TH1F*)Filet->Get("ttbar");
  TH1F* electroweakt = (TH1F*)Filet->Get("Electroweak");
  TH1F* multijett = (TH1F*)Filet->Get("Multijet");
  TH1F* SingleTopSMt = (TH1F*)Filet->Get("t_channel");
  TH1F* chiCarre = new TH1F ("Chi Carre", "" , range*2/0.0001, -range, range);
  TF1* fithist = new TF1 ("pol2","pol2", -range,range);


  TH1F* dattbar = (TH1F*)Filetbar->Get("Asimov");
  TH1F* Sigtbar = (TH1F*)Filetbar->Get("Modulation_tbar");
  TH1F* ttbartbar = (TH1F*)Filetbar->Get("ttbar");
  TH1F* electroweaktbar = (TH1F*)Filetbar->Get("Electroweak");
  TH1F* SingleTopSMtbar = (TH1F*)Filetbar->Get("tbar_channel");
  TH1F* multijettbar = (TH1F*)Filetbar->Get("Multijet");

  double asimovt = datt->GetBinContent(1);
  double ttbartbkg = ttbart->GetBinContent(1)*systttbar;
  double multijettbkg = multijett->GetBinContent(1)*systmulti;
  double SingleTopSMtbkg = SingleTopSMt->GetBinContent(1)*systsingletop;
  double electroweaktbkg = electroweakt->GetBinContent(1)*systelectro;

  double asimovtbar = dattbar->GetBinContent(1);
  double ttbartbarbkg = ttbartbar->GetBinContent(1)*systttbar;
  double multijettbarbkg = multijettbar->GetBinContent(1)*systmulti;
  double electroweaktbarbkg = electroweaktbar->GetBinContent(1)*systelectro;
  double SingleTopSMtbarbkg = SingleTopSMtbar->GetBinContent(1)*systsingletop;

  double sigmatott = sqrt(lumi*(ttbartbkg+multijettbkg+electroweaktbkg+SingleTopSMtbkg));
  double sigmatottbar = sqrt(lumi*(ttbartbarbkg+multijettbarbkg+electroweaktbarbkg+SingleTopSMtbarbkg));

  chiCarre->GetXaxis()->SetTitle("b_{#mu} (GeV)");
  chiCarre->GetYaxis()->SetTitle("#Delta#chi^{2}");
  chiCarre->GetXaxis()->CenterTitle(kTRUE);
  chiCarre->GetYaxis()->CenterTitle(kTRUE);
  chiCarre->SetStats(kFALSE);

  double dividet = 1.0/(sigmatott*sigmatott*24);
  double dividetbar = 1.0/(sigmatottbar*sigmatottbar*24);
  int bin = 1;


  for (double i = -range ; i<range-0.0001; i+=0.0001)
  {
    double tchannel = 0;
    double tbarchannel = 0;

    for (int p=1; p<26; p++)
    {

      //bsign += (datt->GetBinContent(p)-lumi*((i*Sig->GetBinContent(p)+1)*SingleTopSMt->GetBinContent(p)+ttbart->GetBinContent(p)+electroweakt->GetBinContent(p)+multijett->GetBinContent(p)))*(datt->GetBinContent(p)-lumi*((i*Sig->GetBinContent(p)+1)*SingleTopSMt->GetBinContent(p)+ttbart->GetBinContent(p)+electroweakt->GetBinContent(p)+multijett->GetBinContent(p)));
      tchannel += (asimovt-lumi*((i*Sigt->GetBinContent(p)+1)*SingleTopSMtbkg+ttbartbkg+electroweaktbkg+multijettbkg))*(asimovt-lumi*((i*Sigt->GetBinContent(p)+1)*SingleTopSMtbkg+ttbartbkg+electroweaktbkg+multijettbkg));
      tbarchannel += (asimovtbar-lumi*((i*Sigtbar->GetBinContent(p)+1)*SingleTopSMtbarbkg+ttbartbarbkg+electroweaktbarbkg+multijettbarbkg))*(asimovtbar-lumi*((i*Sigtbar->GetBinContent(p)+1)*SingleTopSMtbarbkg+ttbartbarbkg+electroweaktbarbkg+multijettbarbkg));
    }
    chiCarre->SetBinContent(bin,tchannel*dividet+tchannel*dividetbar);
    //chiCarre->SetBinContent(bin,tchannel*dividet);
    bin++;
  }
  TFile* Xroot = new TFile("XCarre.root","RECREATE");
  Xroot->cd();
  chiCarre->Fit(fithist,"R");
  cout<<"Delta Chi at 1 GeV = "<<fithist->GetX(chiCarre->GetBinContent(chiCarre->GetMinimumBin())+1)<<endl;

  //chiCarre->SetLineColor(kBLACK);
  chiCarre->Draw();
  legend->AddEntry(chiCarre, "#chi_{tot}^{2} = #chi_{t}^{2} + #chi_{#bar{t}}^{2}", "l");


/*  chiCarre = new TH1F ("Chi Carre", "" , range*2/0.0001, -range, range);
  lumi=1.02;
  sigmatot = sqrt(lumi*(ttbartbkg+multijettbkg+electroweaktbkg+SingleTopSMtbkg));

  double buff=0;
  bin = 1;

  for (double i = -range ; i<range-0.0001; i+=0.0001)
  {
    double bsign = 0;

    for (int p=1; p<26; p++)
    {
      bsign += (asimov-lumi*((i*Sig->GetBinContent(p)+1)*SingleTopSMtbkg+ttbartbkg+electroweaktbkg+multijettbkg))*(asimov-lumi*((i*Sig->GetBinContent(p)+1)*SingleTopSMtbkg+ttbartbkg+electroweaktbkg+multijettbkg));
    }
    chiCarre->SetBinContent(bin,bsign*divide);
    //cout<<chiCarre->GetBinContent(bin)<<endl;
    bin++;
  }
  //cout<<chiCarre->GetBinContent(chiCarre->GetMinimumBin())<<endl;
  //cout<<"Delta Chi at 1 GeV = "<<fithist->GetX(chiCarre->GetBinContent(chiCarre->GetMinimumBin())+1)<<endl;

  double min = chiCarre->GetBinContent(chiCarre->GetMinimumBin());
  cout<<"min="<<min<<endl;
  for (int bin = 1 ; bin<=140001; bin++)
  {
    buff = chiCarre->GetBinContent(bin);
    chiCarre->SetBinContent(bin,buff-min);

  }

  chiCarre->SetLineColor(kBlue);
  chiCarre->Draw("SAME");
  chiCarre->Write("XCarre+");
  legend->AddEntry(chiCarre, "Luminosity = 1.02", "l");
  chiCarre = new TH1F ("Chi Carre", "" , range*2/0.0001, -range, range);
  buff=0;

  lumi = 0.98;
  sigmatot = sqrt(lumi*(ttbartbkg+multijettbkg+electroweaktbkg+SingleTopSMtbkg));
  bin=1;
  for (double i = -range ; i<range-0.0001; i+=0.0001)
  {
    double bsign = 0;

    for (int p=1; p<26; p++)
    {
      bsign += (asimov-lumi*((i*Sig->GetBinContent(p)+1)*SingleTopSMtbkg+ttbartbkg+electroweaktbkg+multijettbkg))*(asimov-lumi*((i*Sig->GetBinContent(p)+1)*SingleTopSMtbkg+ttbartbkg+electroweaktbkg+multijettbkg));
    }
    chiCarre->SetBinContent(bin,bsign*divide);
    bin++;
  }
  cout<<"Delta Chi at 1 GeV = "<<fithist->GetX(chiCarre->GetBinContent(chiCarre->GetMinimumBin())+1)<<endl;
  min = chiCarre->GetBinContent(chiCarre->GetMinimumBin());

  for (int bin = 1 ; bin<=140001; bin++)
  {
    buff = chiCarre->GetBinContent(bin);
    chiCarre->SetBinContent(bin,buff-min);
  }

  chiCarre->SetLineColor(kRed);
  chiCarre->Draw("SAME");
  chiCarre->Write("XCarre-");
  legend->AddEntry(chiCarre, "Luminosity = 0.98", "l");*/
  legend->Draw();

  feuille->SaveAs("XCarre.pdf");
  Xroot->Close();
  return 0;
}
