#include  "../include/root_simu.hpp"
#include  "../include/matrixSingleTop.hpp"
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <TGraph.h>
#include <TH2F.h>
#include <TLegend.h>

using namespace std;

int main (){

  // Some needed variables
  int choice;
  TString diagram;
  int nbtqprecut = 0;
  int nbtqbarprecut = 0;
  int nbtq = 0;
  int nbtqbar =0;
  int nbtw = 0;
  double SM = 0;
  double ratio =0;
  TLorentzVector SME (0. ,0. ,0. ,0.);
  TLorentzVector Wtq (0., 0., 0., 0.);
  TLorentzVector Wtqbar (0., 0., 0., 0.);
  TLorentzVector Wtw (0., 0., 0., 0.);
  TLorentzVector Wt (0., 0., 0., 0.);
  TLorentzVector WtAv(0., 0., 0., 0.);

  //Choose the channel and diagram we want to study.
  cout<<"We consider an energy in the center-of-mass of 13 TeV for both t-channel and tw-channel"<<endl;
  do {
    cout<<"Choose which channel you wish to study within the 4 diagrams available :"<<endl
    <<"t - t-channel through b and q diffusion"<<endl
    <<"tw - tW-channel"<<endl
    <<"tbar - tbar-channel through antib and q diffusion"<<endl
    <<"tbarw - tbarW-channel"<<endl;
    cin>>diagram;}
  while (diagram !="t" && diagram !="tw" && diagram !="tbar" && diagram !="tbarw" );

  root_simu* go;

  //Extract the corresponding rootfile
  if (diagram == "t")
    go = new root_simu("SingleToptchan.root");
  if (diagram == "tw")
    go = new root_simu("SingleToptwchan.root");
  if (diagram == "tbar")
    go = new root_simu("SingleToptbarchan.root");
  if (diagram == "tbarw")
    go = new root_simu("SingleToptbarwchan.root");
  go->Loop();
  int Events = go->fChain->GetEntriesFast();

  //Declare the variable to sum the values of the matrix elements for each iteration
  for (int i=0; i<Events; i++){

    if (go->nature!=0) {
      MatrixSingleTop *M = new MatrixSingleTop(go->pMother1[i], go->b[i], go->pt[i], go->p3[i], go->pMother1_PID[i], go->p3_PID[i], go->nature[i], go->pTmu[i], go->etamu[i], go->pTelec[i], go->etaElec[i]);
      choice = go->nature[i];
        //Treat t-channel in the case of tq and tbarq
        if (choice == 1) {
          nbtqprecut++;
          if (go->etaElec[i]<1.479 && go->pTelec[i]>35){

            if (diagram == "t")
            {
              SME = M->getMbq_mu();
            }
            if (diagram == "tbar")
            {
              SME = -1.0*M->getMbq_mu();
            }

            SM  = M->getMbq();
            ratio = 1/SM;
            Wtq+=SME*ratio;//Sum each calculated change ratio at each iteration
            nbtq++;         //Count the number of time we calculate a charge ratio to calculate the mean value later
          }
           else {if (go->etamu[i]<2.4 && go->pTmu[i]>26){
             if (diagram == "t")
             {
               SME = M->getMbq_mu();
             }
             if (diagram == "tbar")
             {
               SME = -1.0*M->getMbq_mu();
             }
            SM  = M->getMbq();
            ratio = 1/SM;
            Wtq+=SME*ratio;
            nbtq++;
          }}
        }

        //Treat t-channel in the case of tqbar and tbarqbar
        if (choice == 2 ){
          nbtqbarprecut++;
          if (go->etaElec[i]<1.479 && go->pTelec[i]>35){
            if (diagram == "t")
            {
              SME = M->getMbqbar_mu();
            }
            if (diagram == "tbar")
            {
              SME = -1.0*M->getMbqbar_mu();
            }
            SM  = M->getMbqbar();
            /*cout<<"SMtbarq ="<<SM<<endl;
            for (int i=0; i<4;i++)
            {
              cout<<"SM"<<i<<"="<<SME[i]<<endl;
            }*/
            ratio = 1/SM;
            Wtqbar+=SME*ratio;
            nbtqbar++;
            }
          else {if (go->etamu[i]<2.4 && go->pTmu[i]>26){
            if (diagram == "t")
            {
              SME = M->getMbqbar_mu();
            }
            if (diagram == "tbar")
            {
              SME = -1.0*M->getMbqbar_mu();
            }

            SM  = M->getMbqbar();
            ratio = 1/SM;
            Wtqbar+=SME*ratio;
            nbtqbar++;
          }}
        }

        //Ttreat tW-channel in the case of t and tbar
        if(choice == 3){
          if (go->etaElec[i]<2.4 && go->pTelec[i]>20){
            if (diagram == "tbarw")
            {
              SME = -1.0*M->getMbg_mu();
            }
            if (diagram == "tw")
            {
              SME = M->getMbg_mu();
            }
            SM  = M->getMbg();
            ratio = 1/SM;
            Wtw+=SME*ratio;
            nbtw++;}

          else {if (go->etamu[i]<2.4 && go->pTmu[i]>20){
            if (diagram == "tbarw")
            {
              SME = -1.0*M->getMbg_mu();
            }
            if (diagram == "tw")
            {
              SME = M->getMbg_mu();
            }
            SM  = M->getMbg();
            ratio = 1/SM;
            Wtw+=SME*ratio;
            nbtw++;
          }}
        }

    delete M;

  }
}

//Calcuate the charge ratio considering we have different diagram for a single generation
//When we simule tq and so tqbar, we have different cross-sect for tq and tqbar. The following expression weights each calculated matrix element.
double ratiotq = 1.0/nbtq;
double ratiotqbar=1.0/nbtqbar;
double ratiotw = 1.0/nbtw;
double pourctq = 1.0*nbtqprecut/1000000;
double pourctqbar = 1.0*nbtqbarprecut/1000000;
//If Wtw != 0 then Wtq and Wtqbar should be 0

/*cout<<"Somme des W(t)"<<endl;
for (int i=0;i<4;i++){cout<<Wtq[i]<<endl;}
cout<<"nb="<<nb<<endl;*/

if (diagram == "t"){
  WtAv =pourctq*Wtq*ratiotq + pourctqbar*Wtqbar*ratiotqbar;
  //WtAv = Wtq*ratiotq;
  //WtAv = Wtqbar*ratiotqbar;
  ofstream tchan("Charge_Ratio_t-channel.txt");
  for (int i=0; i<4; i++) {tchan<<WtAv[i]<<endl;}
  tchan.close();
}

if (diagram == "tw" ){
  WtAv =Wtw*ratiotw;
  ofstream twchan("Charge_Ratio_tw-channel.txt");
  for (int i=0; i<4; i++) {twchan<<WtAv[i]<<endl;}
  twchan.close();
}
if (diagram == "tbar" ){
  WtAv =pourctq*Wtq*ratiotq +pourctqbar*Wtqbar*ratiotqbar;
  //WtAv = Wtqbar*ratiotqbar;
  //WtAv = Wtq*ratiotq;

  ofstream twchan("Charge_Ratio_tbar-channel.txt");
  for (int i=0; i<4; i++) {twchan<<WtAv[i]<<endl;}
  twchan.close();
}
if (diagram == "tbarw" ){
  WtAv =Wtw*ratiotw;

  ofstream tbarwchan("Charge_Ratio_tbarw-channel.txt");
  for (int i=0; i<4; i++) {tbarwchan<<WtAv[i]<<endl;}
  tbarwchan.close();
}

//________________Macro creating the f(t) graph using root_______________//

double ft = 0;
double xi_1 = sin(46.309)*cos(101.290);
double xi_2 = sin(101.2790);
double xi_3 = cos(46.309)*cos(101.2790);



TLorentzVector b_mu0 (5, 0., 0., 0.);
TLorentzVector b_mu1 (0., 5, 0., 0.);
TLorentzVector b_mu2 (0., 0., 5, 0.);
TLorentzVector b_mu3 (0., 0., 0., 5);

double min,max;
if (diagram == "t" || diagram == "tbar")
{min = -0.13; max = 0.13;}

if (diagram == "tw" || diagram == "tbarw")
{min = -0.5; max = 0.5  ;}

TFile* File = new TFile ("Modulation_Temporelle_13TeV_b5.root", "update");
TCanvas* modul = new TCanvas("Modulation", "Simple graph", 600 , 600);
TH2F* axe = new TH2F (" ", " ", 49 ,0,24,49,min, max);
auto legend = new TLegend(0.75,0.75,0.95,0.95);
TGraph* graph = new TGraph (48);
TString name;
modul->SetLeftMargin(0.14);
modul->SetRightMargin(0.03);
axe->SetStats(kFALSE);
axe->Draw();
axe->GetXaxis()->SetTitle("t (h)");
axe->GetYaxis()->SetTitle("f(t) = SME/SM -1");
axe->GetXaxis()->CenterTitle(kTRUE);
axe->GetYaxis()->CenterTitle(kTRUE);
//axe->GetYaxis()->SetTitleOffset(0.8);
File->cd();

//Save points to Modulation_Temporelle.root and plot the graph at the same time//

//X = 10 Benchmark
int i = 0;
for (double t=0; t<=24 ; t+=0.5)
{
  ft = b_mu0[3]*WtAv[3] - b_mu0[0]*(xi_1*cos(2*M_PI/24*t)+xi_2*sin(2*M_PI/24*t))*WtAv[2]+b_mu0[1]*(xi_2*cos(2*M_PI/24*t)-xi_1*sin(2*M_PI/24*t))*WtAv[2]+b_mu0[2]*xi_3*WtAv[2];
  graph->SetPoint(i, t, ft);
  i++;
}
graph->SetLineColor(kOrange);
graph->Write(diagram+"X");
legend->AddEntry(graph,"bx = 10 GeV", "l");
graph->Draw("Same");

//Y = 10 Benchmark
i = 0;
graph = new TGraph (48);
for (double t=0; t<=24 ; t+=0.5)
{
  ft = b_mu1[3]*WtAv[3] - b_mu1[0]*(xi_1*cos(2*M_PI/24*t)+xi_2*sin(2*M_PI/24*t))*WtAv[2]+b_mu1[1]*(xi_2*cos(2*M_PI/24*t)-xi_1*sin(2*M_PI/24*t))*WtAv[2]+b_mu1[2]*xi_3*WtAv[2];
  graph->SetPoint(i, t, ft);
  i++;
}
graph->SetLineColor(kRed);
graph->Write(diagram+"Y");
legend->AddEntry(graph,"by = 10 GeV", "l");
graph->Draw("Same");

//Z = 10 Benchmark
i = 0;
graph = new TGraph (48);
for (double t=0; t<=24 ; t+=0.5)
{
  ft = b_mu2[3]*WtAv[3] - b_mu2[0]*(xi_1*cos(2*M_PI/24*t)+xi_2*sin(2*M_PI/24*t))*WtAv[2]+b_mu2[1]*(xi_2*cos(2*M_PI/24*t)-xi_1*sin(2*M_PI/24*t))*WtAv[2]+b_mu2[2]*xi_3*WtAv[2];
  graph->SetPoint(i, t, ft);
  i++;
}
graph->SetLineColor(kBlue);
legend->AddEntry(graph,"bz = 10 GeV", "l");
graph->Write(diagram+"Z");
graph->Draw("Same");

//t = 10 Benchmark
i = 0;
graph = new TGraph (48);
for (double t=0; t<=24 ; t+=0.5)
{
  ft = b_mu3[3]*WtAv[3] - b_mu3[0]*(xi_1*cos(2*M_PI/24*t)+xi_2*sin(2*M_PI/24*t))*WtAv[2]+b_mu3[1]*(xi_2*cos(2*M_PI/24*t)-xi_1*sin(2*M_PI/24*t))*WtAv[2]+b_mu3[2]*xi_3*WtAv[2];
  graph->SetPoint(i, t, ft);
  i++;
}
graph->SetLineColor(kGreen);
graph->Write(diagram+"T");
legend->AddEntry(graph,"bt = 10 GeV", "l");
graph->Draw("Same");


/*graph = new TGraph (48);
for (double t=0; t<=24 ; t+=0.5)
{
  ft = 0;
  graph->SetPoint(i, t, ft);
  i++;
}*/
/*graph->SetLineColor(kBlack);
graph->Write(diagram+"SM");
legend->AddEntry(graph,"SM", "l");
graph->Draw("Same");*/


legend->Draw("Same");
modul->SaveAs("Modulation_"+diagram+".pdf");
return 0;
}
