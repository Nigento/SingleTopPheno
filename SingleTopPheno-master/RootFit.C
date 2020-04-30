#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <TFile.h>
#include <TH1F.h>
#include <RooHistPdf.h>

using namespace std;
using namespace RooFit;

void RootFit()
{


TFile* File = TFile::Open("/home/sane/Stage_M1/LIVMass-master/results/histo/hist_tchannel.root");
TH1F* dat = (TH1F*)File->Get("Asimov");
TH1F* Sig = (TH1F*)File->Get("Modulation_t");
TH1F* Bkg = (TH1F*)File->Get("background");
TH1F* SingleTopSM = (TH1F*)File->Get("t_channel");

RooRealVar t("t","t",dat->GetXaxis()->GetXmin(),dat->GetXaxis()->GetXmax());

RooDataHist* data = new RooDataHist("data","data",t, dat);

  RooDataHist * sig = new RooDataHist("SingleTopSME","SingleTopSME",t, Sig);
  RooDataHist * bkg = new RooDataHist("BKG","BKG",t, Bkg);
  RooDataHist * singletop = new RooDataHist("SingleTopSM","SingleTopSM",t, SingleTopSM);

  RooHistPdf sigpdf("sigpdf","sigpdf",t,*sig) ;
  RooHistPdf bkgpdf("pdfBkg","pdfBkg",t,*bkg) ;
  RooHistPdf SingleToppdf("pdfSingleTop","pdfSingleTop",t,*singletop) ;

  double yieldsig = Sig->Integral();
  double yieldBkg = Bkg->Integral();
  double yieldsingletop = SingleTopSM->Integral();



  RooRealVar nsig("nsig", "nsig",yieldsig,0,dat->Integral());
  RooRealVar nBkg("nBkg", "nBkg",yieldBkg,0,dat->Integral());
  RooRealVar nSingleTop("nSingleTop", "nSingleTop",yieldsingletop,0,dat->Integral());

//Varying widths
  /*RooRealVar sigsigma("sig-sigma","sig-sigma",0.15, "" );
  RooRealVar bkgsigma("bkg-sigma","bkg-sigma",0.15, "" );
  RooRealVar TTsigma("TT-sigma","TT-sigma",0.4, "" );*/
  RooRealVar b_mu ("b_mu","b_mu", -20, 20, "GeV");
  RooRealVar musign("mu_signal", "mu_signal",-5,6 ,"");
  RooRealVar meanVal("Mean", "Mean",1,"");

  //constant mean=1 and width =0.15

  //RooGaussian f1constraint ("f1constraint","f1constraint",nsig,RooConst(1),); //Enlever la contrainte sur le signal
  RooGaussian f2constraint ("f2constraint","f2constraint",nBkg,RooConst(1),RooConst(0.15));
  RooGaussian f3constraint ("f3constraint","f3constraint",nSingleTop,RooConst(1),RooConst(0.40));


  RooProdPdf m1("m1","m1",RooArgSet(sigpdf,f1constraint)) ;
  RooProdPdf m2("m2","m2",RooArgSet(bkgpdf,f2constraint)) ;
  RooProdPdf m3("m3","m3",RooArgSet(SingleToppdf,f3constraint)) ;



  RooArgList st;
  st.add(m1);
  st.add(m2);
  st.add(m3);

  RooArgList yields;
  yields.add(nsig);
  yields.add(nBkg);
  yields.add(nSingleTop);

  //cout<<nsig;


RooAddPdf *pdf = new RooAddPdf("pdf"," pdf", st, yields);

RooFitResult *fr= pdf->fitTo(*data,Save());

}
