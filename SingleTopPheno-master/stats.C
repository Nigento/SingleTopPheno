#include "/root-CERN/root/include/RooStats/HistFactory/Measurement.h"
#include "/root-CERN/root/include/RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include <sstream>
#include <string>

using namespace std;
using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;

void stats(){

RooStats::HistFactory::Measurement meast("bX_tchan", "stats");

            string inputFileName = "/home/sane/Stage_M1/LIVMass-master/results/histo/hist_tchannel_histfactory.root";
            meast.SetOutputFilePrefix("/home/sane/Stage_M1/LIVMass-master/histFactory/example_");
            meast.SetExportOnly(false);

            //meas.SetPOI("SigXsecOverSM");
            meast.SetPOI("muSME");
            meast.AddConstantParam("alpha_syst1");
            meast.AddConstantParam("Lumi");

            meast.SetLumi(1);
            //meast.SetLumiRelErr(0.0001);

            meast.SetExportOnly(false);

            //____________Set the channel with data_________//
            HistFactory::Channel tchannel("tchannel");
            tchannel.SetData("Asimov", inputFileName);
            //tchannel.SetStatErrorConfig(0.05, "Poisson");

            //______________Signal Histogram______________///
            HistFactory::Sample SMEtchan("Modulation_t", "Modulation_t", inputFileName);
            SMEtchan.AddNormFactor("muSME", 1, 0, 3);
            //SMEtchan.AddOverallSys("syst1",0.95,1.05) ;
            tchannel.AddSample(SMEtchan);


            //______________Background Histogram___________//


            HistFactory::Sample tchan("t_channel", "t_channel", inputFileName);
            tchan.ActivateStatError("t_channel_error",inputFileName);
            //tchan.AddOverallSys("syst1",0.90, 1.10);
            tchan.AddNormFactor("t-chan", 1, 0, 3);
            tchannel.AddSample(tchan);

            HistFactory::Sample bkgdttbartchan("ttbar", "ttbar", inputFileName);
            bkgdttbartchan.ActivateStatError("ttbar_error",inputFileName) ;
            //bkgdttbartchan.AddOverallSys("syst2",0.99, 1.01);
            tchannel.AddSample(bkgdttbartchan);

            HistFactory::Sample bkgdMultijettchan("Multijet", "Multijet", inputFileName);
            bkgdMultijettchan.ActivateStatError("Multijet_error",inputFileName);
            //bkgdMultijettchan.AddOverallSys("syst3",0.6, 1.4);
            tchannel.AddSample(bkgdMultijettchan);

            HistFactory::Sample bkgdZWtchan("Electroweak", "Electroweak", inputFileName);
            bkgdZWtchan.ActivateStatError("Electroweak_error",inputFileName);
            //bkgdZWtchan.AddOverallSys("syst4",0.80, 1.20);
            tchannel.AddSample(bkgdZWtchan);

            /*HistFactory::Sample background("background", "background", inputFileName);
            background.ActivateStatError("ttbar_error",inputFileName);
            //background.AddOverallSys("syst2",0.8,1.2);
            tchannel.AddSample(background);*/

            meast.AddChannel(tchannel);
            meast.CollectHistograms();
            meast.PrintTree();
            MakeModelAndMeasurementFast(meast);

/*RooStats::HistFactory::Measurement meastbar("bX_tbarchan", "stats");

                        ostringstream soutbarchan;
                        soutbarchan<<"histFactory/tbarchan";
                        meastbar.SetOutputFilePrefix(soutbarchan.str());
                        meastbar.SetExportOnly(false);

                        //meas.SetPOI("SigXsecOverSM");
                        meastbar.AddPOI("muSME");

                        meastbar.SetLumi(1);
                        meastbar.SetLumiRelErr(0.02);

                        RooStats::HistFactory::Channel tbarchannel("tbarchannel");

                        tbarchannel.SetData("Asimov", "/home/sane/Stage_M1/LIVMass-master/results/histo/hist_tbarchannel.root");

                        RooStats::HistFactory::Sample tbarchan("tbar_channel", "tbar_channel", "/home/sane/Stage_M1/LIVMass-master/results/histo/hist_tbarchannel.root");
                        tbarchan.AddNormFactor("SigXsecOverSM", 1, -1, 1);
                        tbarchan.AddOverallSys("syst1",  0.95, 1.05);
                        tbarchannel.AddSample(tbarchan);

                        RooStats::HistFactory::Sample bkgdttbartbarchan("ttbar", "ttbar", "/home/sane/Stage_M1/LIVMass-master/results/histo/hist_tbarchannel.root");
                        bkgdttbartbarchan.AddOverallSys("syst2", 0.96, 1.4 );
                        tbarchannel.AddSample(bkgdttbartbarchan);

                        RooStats::HistFactory::Sample bkgdMultijettbarchan("Multijet", "Multijet", "/home/sane/Stage_M1/LIVMass-master/results/histo/hist_tbarchannel.root");
                        bkgdMultijettbarchan.AddOverallSys("syst3", 0.6, 1.4 );
                        tbarchannel.AddSample(bkgdMultijettbarchan);

                        RooStats::HistFactory::Sample bkgdZWtbarchan("Electroweak", "Electroweak", "/home/sane/Stage_M1/LIVMass-master/results/histo/hist_tbarchannel.root");
                        bkgdZWtbarchan.AddOverallSys("syst4", 0.7, 1.3 );
                        tbarchannel.AddSample(bkgdZWtbarchan);

                        RooStats::HistFactory::Sample SMEtbar("Modulation_tbar", "Modulation_tbar", "/home/sane/Stage_M1/LIVMass-master/results/histo/hist_tbarchannel.root");
                        SMEtbar.AddNormFactor("muSME", 1, -1, 1);
                        tbarchannel.AddSample(SMEtbar);


                        meastbar.AddChannel(tbarchannel);
                        meastbar.CollectHistograms();
                        meastbar.PrintTree();
                        MakeModelAndMeasurementFast(meastbar);*/
                        cout<<"Test programme"<<endl;
}
