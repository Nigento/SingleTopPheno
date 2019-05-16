
#include "../include/matrixSingleTop.hpp"
#include <iostream>
#include <vector>

using namespace std;

double ckm[3][3] = {{ 0.97446, 0.22452, 0.00365},{ 0.22438,0.97359,0.04214 },{0.00896,0.04133,0.999105}};
/*
//______________ ckm Elements_____________//
TMatrixD ckm(3,3);
//______u_____________________c_________________t___________/
ckm(0,0) = 0.97427; ckm[0][1] = 0.22534; ckm[0][2] = 0.00351;  // d
ckm[1][0] = 0.2252;  ckm[1][1] = 0.97344; ckm[1][2] = 0.0412;   // s
ckm[2][0] = 0.00867; ckm[2][1] = 0.0404;  ckm[2][2] = 0.999146; // b
*/
// Product of ckm elements for t and b particule
const double ckm1  = ckm[2][2]*ckm[2][2];

//____________Matrix elements SM single top_______________//

//*************** Contructor ******************//

MatrixSingleTop::MatrixSingleTop (TLorentzVector pMother1_user, TLorentzVector b_user, TLorentzVector pt_user, TLorentzVector p3_user, int pMother1_PID_user, int p3_PID_user, int nature_user, double pTmu_user, double etamu_user, double pTelec_user, double etaElec_user){

  //Particule identity
  pMother1_PID = pMother1_PID_user;
  p3_PID     = p3_PID_user;
  nature     = nature_user;
  pTelec     = pTelec_user;
  etaElec    = etaElec_user;
  pTmu       = pTmu_user;
  etamu      = etamu_user;
  //Creation
  pMother1 = pMother1_user;
  b      = b_user;
  pt     = pt_user;
  p3     = p3_user;

  if (pMother1_PID==2||pMother1_PID==1) pMother1_PID=0;
  if (pMother1_PID==4||pMother1_PID==3) pMother1_PID=1;
  if (pMother1_PID==5||pMother1_PID==6) pMother1_PID=2;
  if (p3_PID==2||p3_PID==1) p3_PID=0;
  if (p3_PID==4||p3_PID==3) p3_PID=1;
  if (p3_PID==5||p3_PID==6) p3_PID=2;

    if(nature == 1) {Mbq    = calculateMbq(); Mbq_mu = calculateMbq_mu();}

    if(nature == 2) {Mbqbar = calculateMbqbar(); Mbqbar_mu = calculateMbqbar_mu();}

    if(nature == 3){Mbg = calculateMbg(); Mbg_mu = calculateMbg_mu();}

}
MatrixSingleTop::~MatrixSingleTop(){};
//_____________Function Mbg___________//

double MatrixSingleTop::calculateMbg(){

  //Preliminar

  double s = normMinkowski2(pt + p3);
  double t = normMinkowski2(b - pt);
  double frac1 = mt2/mW2;
  double frac2 = s/(t-mt2);
  double frac3 = (mt2-2*mW2+t)/s;
  double frac4 = (mt2-mW2)/(t-mt2);
  double frac5 = mt2/(t-mt2);
  double frac6 = (mt2-mW2)/s;
  double constant = gS2*gW2*ckm1/24;

  //return Mbg
  return  constant*(-2*frac1-(frac1+2)*(frac2+frac3+2*frac4*(frac5+frac6+1)));
}

//______________Function Mbq______________//
/*void MatrixSingleTop::checkckm()
{double ckm2 = ckm [pMother1_PID][p3_PID]*ckm [pMother1_PID][p3_PID];
cout<<"pMother1_PID="<<pMother1_PID<<" and b_PID="<<b_PID<<" and p3_PID="<<p3_PID<<" and pt_PID="<<pt_PID<<endl;}
*/
double MatrixSingleTop::calculateMbq(){

  //Preliminar

  double s = normMinkowski2(pt + p3);
  double t = normMinkowski2(b - pt);
  double ckm2 = ckm [pMother1_PID][p3_PID]*ckm [pMother1_PID][p3_PID];

  //Mbq

  double Mbq = 0.25*gW4*ckm1*ckm2*s*(s - mt2)/((t - mW2)*(t - mW2));

  return Mbq;
}
//______________Function Mbqbar______________//
double MatrixSingleTop::calculateMbqbar(){

  //Preliminar

  double u = normMinkowski2(pMother1 - pt);
  double t = normMinkowski2(b - pt);
  double ckm2 = ckm [pMother1_PID][p3_PID]*ckm [pMother1_PID][p3_PID];
  //Mbqbar

  return 0.25*gW4*ckm1*ckm2*u*(u-mt2)/((t-mW2)*(t-mW2));
}

/*****************************************************/

//___________________Matrix elements SME for single top_______________________//

    //___________Function Mbg_mu__________//

    TLorentzVector MatrixSingleTop::calculateMbg_mu(){

      //Preliminar

      double s = normMinkowski2(pt + p3);
      double t = normMinkowski2(pMother1 - pt);
      double constant = ckm1*gW2*gS2/12;
      double frac1    = 1/(mW2*s);
      double frac2    = mt2*mt2/((t-mt2)*(t-mt2)*(t-mt2));
      double frac3    = mt2/mW2;
      double frac4    = 1/(mW2*s*(t-mt2));
      double frac5    = mt2/((t-mt2)*(t-mt2));
      double frac6    = mt2/s;
      double frac7    = s/mW2;
      double frac8    = mW2/mt2;
      double frac9    = (8*mW2-12*mt2)/s;
      double frac10   = (4*mW2-4*mt2)/s;
      double frac11   = (mt2+3*s)/mW2;

      //return Mbg_mu
      return -constant*(frac1*(pMother1*(mt2-2*mW2)+p3*mt2+pt*t)+8*frac2*(frac3*p3-p3-b)+frac4*(b*(2*mt2*mt2-4*mW2*mW2-mt2*s)
      +pMother1*(3*mt2*mt2-5*mt2*mW2)+(s+mt2-mW2)*(4*mt2*p3-s*pt)+ mW2*s*(2*p3+pt))+frac5*(b*(3*frac3-3)
      +2*pMother1*(frac3*mt2/s-frac6+frac7-frac8)+p3*(4*mt2*frac3/s + frac9 + 9*frac3 - 13) + pt*(frac10 - frac11 + 2*s/mt2 - 5)));
    }

    //____________Function Mbq_mu____________//
    TLorentzVector MatrixSingleTop::calculateMbq_mu(){

      //Preliminar

      double t = normMinkowski2(b - pt);
      double ckm2= ckm[pMother1_PID][p3_PID]*ckm[pMother1_PID][p3_PID];
      double scal = 0.25*gW4*ckm1*ckm2*2/((t-mW2)*(t-mW2));



      //return Mbq_mu
      return scal*p3;
    }

      //_________Function Mbqbar_mu_______//
      TLorentzVector MatrixSingleTop::calculateMbqbar_mu(){

        //Preliminar

        double t = normMinkowski2(b - pt);
        double ckm2= ckm[pMother1_PID][p3_PID]*ckm[pMother1_PID][p3_PID];
        double scal = -0.25*gW4*ckm1*ckm2*2/((t - mW2)*(t - mW2));
        //TLorentzVector test (1.,1.,1.,1.);

        //return Mbqbar_mu
        return scal*pMother1;
      }


/*************************************************/
//_________________Show Outstream_______________//

//Show SM

double MatrixSingleTop::getMbq(){return Mbq;}
double MatrixSingleTop::getMbg(){return Mbg;}
double MatrixSingleTop::getMbqbar(){return Mbqbar;}

//Show SME

TLorentzVector MatrixSingleTop::getMbg_mu(){return Mbg_mu;}
TLorentzVector MatrixSingleTop::getMbq_mu(){return Mbq_mu;}
TLorentzVector MatrixSingleTop::getMbqbar_mu(){return Mbqbar_mu;}

//__________________Some Tools__________________//

double MatrixSingleTop::normMinkowski2(TLorentzVector p){
    return p.Mag2();
}
