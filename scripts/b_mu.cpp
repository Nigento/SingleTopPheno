#include  "../include/root_simu.hpp"
#include  "../include/matrixSingleTop.hpp"

#include <iostream>
#include <vector>

using namespace std;

int main (){



  // Some needed variables
  int choice;
  int diagram;
  int nbtq = 0;
  int nbtqbar =0;
  int nbtw = 0;
  double SM = 0;
  double tqXsec = 0;
  double tqbarXsec = 0;
  TLorentzVector SME (0. ,0. ,0. ,0.);
  TLorentzVector Wtq (0., 0., 0., 0.);
  TLorentzVector Wtqbar (0., 0., 0., 0.);
  TLorentzVector Wtw (0., 0., 0., 0.);
  TLorentzVector Wt (0., 0., 0., 0.);
  TLorentzVector WtAv(0., 0., 0., 0.);

  //Choose the channel and diagram we want to study.
  cout<<"We consider an energy in the center-of-mass of 13 TeV for both t-channel and tw-channel"<<endl;
  do {
    cout<<"Choose which channel you wish to study within the 6 diagram available :"<<endl
    <<"1 - t-channel through b and q diffusion"<<endl
    <<"2 - t-channel through b and antiq diffusion"<<endl
    <<"3 - tW-channel"<<endl
    <<"4 - tbar-channel through antib and q diffusion"<<endl
    <<"5 - tbar-channel through antib and antiq diffusion"<<endl
    <<"6 - tbarW-channel"<<endl;
    cin>>diagram;}
  while (diagram !=1 && diagram !=2 && diagram !=3 && diagram !=4 && diagram !=5 && diagram !=6);

  root_simu* go;

  //Extract the corresponding rootfile
  if (diagram == 1)
    go = new root_simu("SingleTopTQ.root");
  if (diagram == 2)
    go = new root_simu("SingleTopTQbar.root");
  if (diagram == 3)
    go = new root_simu("");
  if (diagram == 4)
    go = new root_simu("");
  if (diagram == 5)
    go = new root_simu("");
  if (diagram == 6)
    go = new root_simu("");

  go->Loop();
  int Events = go->fChain->GetEntriesFast();

  //Declare the variable to sum the values of the matrix elements for each iteration
  for (int i=0; i<Events; i++){
    if (go->nature!=0) {
      MatrixSingleTop *M = new MatrixSingleTop(go->pMother1[i], go->b[i], go->pt[i], go->p3[i], go->pMother1_PID[i], go->p3_PID[i], go->nature[i], go->pTmu[i], go->etamu[i], go->pTelec[i], go->etaElec[i]);
      choice = go->nature[i];

        //Treat t-channel in the case of tq and tbarq
        if (choice == 1) {
          if (go->etaElec[i]<1.479 && go->pTelec[i]>35){
            SM  = M->getMbq();
            SME = M->getMbq_mu();
            Wtq+=SME*(1/SM);//Sum each calculated change ratio at each iteration
            nbtq++;         //Count the number of time we calculate a charge ratio to calculate the mean value later
          }
           else {if (go->etamu[i]<2.4 && go->pTmu[i]>26){
            SM  = M->getMbq();
            SME = M->getMbq_mu();
            Wtq+=SME*(1/SM);
            nbtq++;
          }}
        }

        //Treat t-channel in the case of tqbar and tbarqbar
        if (choice == 2 ){
          if (go->etaElec[i]<1.479 && go->pTelec[i]>35){
            SM  = M->getMbqbar();
            SME = M->getMbqbar_mu();

            Wtqbar+=SME*(1/SM);

            ntqbar++;
            }
          else {if (go->etamu[i]<2.4 && go->pTmu[i]>26){
            SM  = M->getMbqbar();
            SME = M->getMbqbar_mu();
            Wtqbar+=SME*(1/SM);
            nbtqbar++;
          }}
        }

        //Ttreat tW-channel in the case of t and tbar
        if(choice == 3){
          if (go->etaElec[i]<2.4 && go->pTelec[i]>20){
            SM  = M->getMbg();
            SME = M->getMbg_mu();
            Wtw+=SME*(1/SM);
            nbtw++;}
          else {if (go->etamu[i]<2.4 && go->pTmu[i]>20){
            SM  = M->getMbg();
            SME = M->getMbg_mu();
            Wtw+=SME*(1/SM);
            nbtw++;
          }}
        }

    delete M;

  }
}

//Calcuate the charge ratio considering we have different diagram for a single generation
//When we simule tq and so tqbar, we have different cross-sect for tq and tqbar. The following expression weights each calculated matrix element.

//If Wtw != 0 then Wtq and Wtqbar should be 0
WtAv = (tqXsec*Wtq/nbtq + tqbarXsec*wtqbar/nbtqbar + Wtw/nbtw)

/*
cout<<"Somme des W(t)"<<endl;
for (int i=0;i<4;i++){cout<<Wt[i]<<endl;}
cout<<"nb="<<nb<<endl;*/

/*TLorentzVector b_mu (0., 0., 0., 0.);
cout<<"Choose the condition on b_mu you want to apply:"<<endl
<<"b_0=";cin>>b_mu[0];cout<<endl;
cout<<"b_1";cin>>b_mu[1];cout<<endl;
cout<<"b_2";cin>>b_mu[2];cout<<endl;
cout<<"b_3";cin>>b_mu[3];cout<<endl;*/

cout<<"W(t)/nb="<<endl;
cout<<"Moyenne de W(t)"<<endl;
for (int i=0;i<4;i++){cout<<WtAv[i]<<endl;}


ofstream file("Charge_Ratio_t-channel_bq.txt");

for (int i=0; i<4; i++) {file<<WtAv[i]<<endl;}
file.close();
return 0;
}
