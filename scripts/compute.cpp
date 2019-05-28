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
  double asym;

  TFile* File = TFile::Open("/home/sane/Stage_M1/LIVMass-master/results/Modulation_Temporelle_13TeV_b60.root");
  TFile* Save = new TFile ("Asymmetry.root","update");
  TCanvas * modul = new TCanvas("Asymmetry", "Simple graph", 450 , 450);
  TH2F* axe = new TH2F (" ", " ", 49 ,0,24,100,0,0.4 );
  auto legend = new TLegend(0.75,0.15,0.95,0.35);



  axe->SetStats(kFALSE);
  axe->Draw();
  axe->GetXaxis()->SetTitle("t (h)");
  axe->GetYaxis()->SetTitle("Asymmetry (t)");
  axe->GetYaxis()->CenterTitle(kTRUE);
  axe->GetXaxis()->CenterTitle(kTRUE);
  axe->GetYaxis()->SetTitleOffset(1.35);
  Save->cd();

  TGraph* Parti = (TGraph*)File->Get("tX");
  TGraph* AntiParti = (TGraph*)File->Get("tbarX");
  TGraph* Asym = new TGraph (48);
  double* Particle = Parti->GetY();
  double* AntiParticle = AntiParti->GetY();

  double Xsect = 12.26;
  double Xsectbar = 7.004;
  double l = 0;
for (int i = 0; i<=48; i++ )
  {
    asym = ((Xsect - Xsectbar)+Particle[i]*Xsect - Xsectbar*AntiParticle[i])/((Xsect + Xsectbar)+Particle[i]*Xsect + AntiParticle[i]*Xsectbar);
    Asym->SetPoint(i,l,asym);
    cout<<"Asymx="<<asym<<endl;
    l+=0.5;
  }
  Asym->SetLineColor(kOrange);
  Asym->Write("tX");
  legend->AddEntry(Asym,"bx = 10 GeV", "l");
  Asym->Draw("Same");

  Parti = (TGraph*)File->Get("tY");
  AntiParti = (TGraph*)File->Get("tbarY");
  Particle = Parti->GetY();
  AntiParticle = AntiParti->GetY();
  Asym = new TGraph (48);

  l=0;
  for (int i = 0; i <=48; i++)
  {
    asym = ((Xsect - Xsectbar)+Particle[i]*Xsect - Xsectbar*AntiParticle[i])/((Xsect + Xsectbar)+Particle[i]*Xsect + AntiParticle[i]*Xsectbar);
    Asym->SetPoint(i,l,asym);
    l+=0.5;
  }
  Asym->SetLineColor(kRed);
  Asym->Write("tY");
  legend->AddEntry(Asym,"by = 10 GeV", "l");
  Asym->Draw("Same");


  Parti = (TGraph*)File->Get("tZ");
  AntiParti = (TGraph*)File->Get("tbarZ");
  Particle = Parti->GetY();
  AntiParticle = AntiParti->GetY();
  Asym = new TGraph (48);

l=0;
  for (int i = 0; i <=48; i++ )
  {

    asym = ((Xsect - Xsectbar)+Particle[i]*Xsect - Xsectbar*AntiParticle[i])/((Xsect + Xsectbar)+Particle[i]*Xsect + AntiParticle[i]*Xsectbar);
    Asym->SetPoint(i,l,asym);
    l+=0.5;
  }
  Asym->SetLineColor(kBlue);
  Asym->Write("tZ");
  legend->AddEntry(Asym,"bz = 10 GeV", "l");
  Asym->Draw("Same");


  Parti = (TGraph*)File->Get("tT");
  AntiParti = (TGraph*)File->Get("tbarT");
  Particle = Parti->GetY();
  AntiParticle = AntiParti->GetY();
  Asym = new TGraph (48);

l=0;
  for (int i = 0; i <= 48; i++)
  {
    asym = ((Xsect - Xsectbar)+Particle[i]*Xsect - Xsectbar*AntiParticle[i])/((Xsect + Xsectbar)+Particle[i]*Xsect + AntiParticle[i]*Xsectbar);
    Asym->SetPoint(i,l,asym);
    l+=0.5;
  }
  Asym->SetLineColor(kGreen);
  Asym->Write("tT");
  legend->AddEntry(Asym,"bt = 10 GeV", "l");
  legend->Draw("Same");
  Asym->Draw("Same");

  Parti = (TGraph*)File->Get("tT");
  Asym = new TGraph (48);
  l=0;
  for (int i = 0; i <= 48; i++)
  {
    asym = (Xsect - Xsectbar)/(Xsect + Xsectbar);
    Asym->SetPoint(i,l,asym);
    l+=0.5;
  }
  Asym->SetLineColor(kBlack);
  Asym->Write("SM");
  legend->AddEntry(Asym,"SM", "l");
  legend->Draw("Same");
  Asym->Draw("Same");

modul->SaveAs("Asymmetry_t.pdf");

  return 0;
}
