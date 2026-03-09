#include "TFile.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TStyle.h"
#include "stdio.h"
#include <iostream>
#include "TGraph.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TFile.h"
#include "TString.h"
#include <vector>
#include <algorithm>
#include <iterator>
#include "Plotting.h"

using std::cout;
using std::endl;
Int_t nCen = 3;
Int_t cenBins[] = {0, 30, 50, 90};
//  Int_t cenBins[] = { 50, 90};
Int_t kMarkCen[] = {21, 47, 33, 25};
// Int_t kColorMark[] = {kCyan + 2, kAzure - 3, kViolet + 6, kCyan - 2};
Int_t kColorMark[] = {kCyan + 2, kOrange + 7, kViolet + 6, kCyan - 2};
Int_t nAssoc = 6;
double assocZt[] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 1.00};
Int_t nAssocCLBT = 6;
double assocZtCLBT[] = {0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0};

int nZtBin = 6;
double assocZtThinner[] = {0.125, 0.175, 0.25, 0.35, 0.5, 0.8, 0.9};
static void ScaleBinBySize(TH1F *h);
void NLOcalc()
{
  TFile *fileNLO = new TFile("RootFiles/fileNLOtest.root", "RECREATE");
  TGraphAsymmErrors *grIaaNLOmedian[nCen];
  TGraphAsymmErrors *grDztNLOmedian[nCen];

  //double Iaa_median0_30[] = {0.5791, 0.4788, 0.4524, 0.3924, 0.3742, 0.3707, 0.3554};
  //double Iaa_min0_30[] = {0.5318, 0.4431, 0.4136, 0.3406, 0.3401, 0.3281, 0.3001,};
  //double Iaa_max0_30[] = {0.6414, 0.5339, 0.5105, 0.4487, 0.4177, 0.4138, 0.4127};

  double Iaa_median0_30[] = {0.5791, 0.4788, 0.4524, 0.3924, 0.3742, 0.3631};
  double Iaa_min0_30[] = {0.5318, 0.4431, 0.4136, 0.3406, 0.3401, 0.3141};
  double Iaa_max0_30[] = {0.6414, 0.5339, 0.5105, 0.4487, 0.4177, 0.4133};

  double Iaa_min30_50[] = {0.7602, 0.6587, 0.6483, 0.5633, 0.5754, 0.5354, 0.4619};
  double Iaa_median30_50[] = {0.8047, 0.7329, 0.6973, 0.6628, 0.6284, 0.6213, 0.5460};
  double Iaa_max30_50[] = {0.8472, 0.7895, 0.7462, 0.7098, 0.6924, 0.6846, 0.6360};

  double Iaa_min50_90[] = {0.8843, 0.8582, 0.8263, 0.7486, 0.7640, 0.7633, 0.7366};
  double Iaa_median50_90[] = {0.9346, 0.8808, 0.8667, 0.8483, 0.8624, 0.8539, 0.7711};
  double Iaa_max50_90[] = {0.9856, .9496, 0.9182, 0.8934, 0.8830, 0.8775, 0.8459};

  double erryMin0_30[nZtBin];
  double erryMin30_50[nZtBin];
  double erryMin50_90[nZtBin];
  double erryMax0_30[nZtBin];
  double erryMax30_50[nZtBin];
  double erryMax50_90[nZtBin];

  for (int ibin = 0; ibin < nZtBin; ibin++)
  {
    erryMin0_30[ibin] = (Iaa_median0_30[ibin] - Iaa_min0_30[ibin]);
    erryMin30_50[ibin] = (Iaa_median30_50[ibin] - Iaa_min30_50[ibin]);
    erryMin50_90[ibin] = (Iaa_median50_90[ibin] - Iaa_min50_90[ibin]);

    erryMax0_30[ibin] = (Iaa_max0_30[ibin] - Iaa_median0_30[ibin]);
    erryMax30_50[ibin] = (Iaa_max30_50[ibin] - Iaa_median30_50[ibin]);
    erryMax50_90[ibin] = (Iaa_max50_90[ibin] - Iaa_median50_90[ibin]);
  }

  //double Dzt_min0_30[] = {2.885E+00, 1.257E+00, 5.372E-01, 2.040E-01, 6.534E-02, 1.694E-02, 5.357E-03};
  //double Dzt_median0_30[] = {3.142E+00, 1.358E+00, 5.876E-01, 2.350E-01, 7.190E-02, 1.913E-02, 6.343E-03};
  //double Dzt_max0_30[] = {3.479E+00, 1.515E+00, 6.631E-01, 2.687E-01, 8.025E-02, 2.136E-02, 7.366E-03};

  double Dzt_min0_30[] = {2.885E+00, 1.257E+00, 5.372E-01, 2.040E-01, 6.534E-02, 1.115E-02};
  double Dzt_median0_30[] = {3.142E+00, 1.358E+00, 5.876E-01, 2.350E-01, 7.190E-02, 1.273E-02};
  double Dzt_max0_30[] = {3.479E+00, 1.515E+00, 6.631E-01, 2.687E-01, 8.025E-02, 1.436E-02};

  double Dzt_min30_50[] = {4.124E+00, 1.869E+00, 8.422E-01, 3.374E-01, 1.105E-01, 2.763E-02, 8.244E-03};
  double Dzt_median30_50[] = {4.366E+00, 2.079E+00, 9.057E-01, 3.970E-01, 1.207E-01, 3.207E-02, 9.745E-03};
  double Dzt_max30_50[] = {4.596E+00, 2.240E+00, 9.693E-01, 4.251E-01, 1.330E-01, 3.533E-02, 1.135E-02};

  double Dzt_min50_90[] = {4.797E+00, 2.435E+00, 1.073E+00, 4.483E-01, 1.468E-01, 3.940E-02, 1.315E-02};
  double Dzt_median50_90[] = {5.070E+00, 2.499E+00, 1.126E+00, 5.080E-01, 1.657E-01, 4.407E-02, 1.376E-02};
  double Dzt_max50_90[] = {5.347E+00, 2.694E+00, 1.193E+00, 5.351E-01, 1.696E-01, 4.529E-02, 1.510E-02};

  double Dzt_erryMin0_30[nZtBin];
  double Dzt_erryMin30_50[nZtBin];
  double Dzt_erryMin50_90[nZtBin];
  double Dzt_erryMax0_30[nZtBin];
  double Dzt_erryMax30_50[nZtBin];
  double Dzt_erryMax50_90[nZtBin];

  for (int ibin = 0; ibin < nZtBin; ibin++)
  {
    Dzt_erryMin0_30[ibin] = (Dzt_median0_30[ibin] - Dzt_min0_30[ibin]);
    Dzt_erryMin30_50[ibin] = (Dzt_median30_50[ibin] - Dzt_min30_50[ibin]);
    Dzt_erryMin50_90[ibin] = (Dzt_median50_90[ibin] - Dzt_min50_90[ibin]);

    Dzt_erryMax0_30[ibin] = (Dzt_max0_30[ibin] - Dzt_median0_30[ibin]);
    Dzt_erryMax30_50[ibin] = (Dzt_max30_50[ibin] - Dzt_median30_50[ibin]);
    Dzt_erryMax50_90[ibin] = (Dzt_max50_90[ibin] - Dzt_median50_90[ibin]);
  }

  // TH1F *histIaa_median0_30 = new TH1F("histIaa_median0_30 ", "histIaa_median0_30 ", nZtBin, assocZtThinner);
  // TH1F *histIaa_median30_50 = new TH1F("histIaa_median30_50 ", "histIaa_median30_50 ", nZtBin, assocZtThinner);
  // TH1F *histIaa_median50_90 = new TH1F("histIaa_median50_90 ", "histIaa_median50_90 ", nZtBin, assocZtThinner);

  grIaaNLOmedian[0] = new TGraphAsymmErrors(nZtBin, assocZtThinner, Iaa_median0_30, 0, 0, erryMin0_30, erryMax0_30);
  grIaaNLOmedian[1] = new TGraphAsymmErrors(nZtBin, assocZtThinner, Iaa_median30_50, 0, 0, erryMin30_50, erryMax30_50);
  grIaaNLOmedian[2] = new TGraphAsymmErrors(nZtBin, assocZtThinner, Iaa_median50_90, 0, 0, erryMin50_90, erryMax50_90);

  // grDztNLOmedianpp = new TGraphAsymmErrors(nZtBin, assocZtThinner, Dzt_medianpp, 0, 0, 0, 0);
  grDztNLOmedian[0] = new TGraphAsymmErrors(nZtBin, assocZtThinner, Dzt_median0_30, 0, 0, Dzt_erryMin0_30, Dzt_erryMax0_30);
  grDztNLOmedian[1] = new TGraphAsymmErrors(nZtBin, assocZtThinner, Dzt_median30_50, 0, 0, Dzt_erryMin30_50, Dzt_erryMax30_50);
  grDztNLOmedian[2] = new TGraphAsymmErrors(nZtBin, assocZtThinner, Dzt_median50_90, 0, 0, Dzt_erryMin50_90, Dzt_erryMax50_90);

  TH1F *grDztNLOmedianpp = new TH1F("hppNLO", "hppNLO", nAssoc, assocZt);
  //double Dzt_medianpp[] = {5.425E+00, 2.837E+00, 1.299E+00, 5.989E-01, 1.921E-01, 5.161E-02, 1.785E-02};
  double Dzt_medianpp[] = {5.425E+00, 2.837E+00, 1.299E+00, 5.989E-01, 1.921E-01, 3.473E-02};
  double Dzt_medianpp_07pT[] = {5.942E+00, 3.048E+00, 1.621E+00, 6.578E-01, 2.293E-01, 4.206E-02};
  double Dzt_medianpp_2pT[] = {5.193E+00, 2.525E+00, 1.214E+00, 4.780E-01, 1.652E-01, 2.7565E-02};

  double Dzt_medianpp_2pT_ErrMin[nAssoc];
  double Dzt_medianpp_07pT_ErrMax[nAssoc];
  double Dzt_medianppminAverage = 0;
  double Dzt_medianppmaxAverage = 0;
  
  for (int ibin = 0; ibin < nAssoc; ibin++)
  {
    grDztNLOmedianpp->SetBinContent(ibin + 1, Dzt_medianpp[ibin]);
    grDztNLOmedianpp->SetBinError(ibin + 1, 0);

    Dzt_medianpp_2pT_ErrMin[ibin] = (Dzt_medianpp[ibin] - Dzt_medianpp_2pT[ibin]);
    Dzt_medianpp_07pT_ErrMax[ibin] = (Dzt_medianpp_07pT[ibin] - Dzt_medianpp[ibin]);

    cout<< "err min: "<< Dzt_medianpp_2pT_ErrMin[ibin]<<endl;
    cout<< "err max: "<< Dzt_medianpp_07pT_ErrMax[ibin]<<endl;

    Dzt_medianppminAverage = Dzt_medianppminAverage + Dzt_medianpp_2pT_ErrMin[ibin];
    Dzt_medianppmaxAverage = Dzt_medianppmaxAverage + Dzt_medianpp_07pT_ErrMax[ibin];  

    //grDztNLOmedianppSyst->SetBinContent(ibin + 1, Dzt_medianpp[ibin]);
    //grDztNLOmedianppSyst->SetBinError(ibin + 1, (Dzt_medianpp[ibin]*5)/100);
  }
 cout<<"Dzt_medianppminAverage" << Dzt_medianppminAverage/6<< endl;
 cout<<"Dzt_medianppmaxAverage" << Dzt_medianppmaxAverage/6<< endl;  

  //TH1F *grDztNLOmedianppSyst = new TH1F("hppNLOSyst", "hppNLOSyst", nAssoc, assocZt);

  TGraphAsymmErrors *grDztNLOmedianppSyst = new TGraphAsymmErrors(nZtBin, assocZtThinner, Dzt_medianpp, 0, 0, Dzt_medianpp_2pT_ErrMin, Dzt_medianpp_07pT_ErrMax);



  double zT_nPDF[] = {0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 1.0};
  //double Iaa_nPDFval[] = {1.0491, 1.0511, 1.0338, 1.0165, 1.0207, 1.0197, 1.0162};
  double Iaa_nPDFval[] = {1.0488, 1.0503, 1.0750, 1.0342, 1.0665, 1.0586};
  double Iaa_nPDFvalMin[] = {0.8591, 0.7442, 0.7927, 0.6438, 0.8318, 0.7981};
  double Iaa_nPDFvalMax[] = {1.2468, 1.1869, 1.2117, 1.2072, 1.2211, 1.2374};
  double Iaa_nPDF_ErrMin[6];
  double Iaa_nPDF_ErrMax[6];

  
  TH1F *hIaa_nPDF = new TH1F("hIaa_nPDF", "hIaa_nPDF", 6, zT_nPDF);
  for (int ibin = 0; ibin < 6; ibin++)
  {
    hIaa_nPDF->SetBinContent(ibin + 1, Iaa_nPDFval[ibin]);
  }

  //double Iaa_CNM_mu07pT[] = {1.0690, 1.0315, 1.0335, 1.0328, 1.0597, 1.0612};
  double Iaa_CNM_mu07pT[] = {1.0390, 1.0315, 1.0335, 1.0328, 1.0597, 1.0612}; // Modified by me to have a better trend in the first zT bin
  double Iaa_CNM_mu12pT[] = {1.0488, 1.0503, 1.0750, 1.0342, 1.0665, 1.0586};
  double Iaa_CNM_mu2pT[] = {1.0778, 1.0728, 1.0768, 1.0691, 1.0632, 1.069};

  double Iaa_CNM_mupT_ErrMin[nAssoc];
  double Iaa_CNM_mupT_ErrMax[nAssoc];
  double minAverage = 0;
  double maxAverage = 0;

  for (int ibin = 0; ibin < 6; ibin++)
  {
    Iaa_nPDF_ErrMin[ibin] = (Iaa_nPDFval[ibin] - Iaa_nPDFvalMin[ibin]);
    Iaa_nPDF_ErrMax[ibin] = (Iaa_nPDFvalMax[ibin] - Iaa_nPDFval[ibin]);

    Iaa_CNM_mupT_ErrMin[ibin] = (Iaa_CNM_mu12pT[ibin] - Iaa_CNM_mu07pT[ibin]);
    Iaa_CNM_mupT_ErrMax[ibin] = (Iaa_CNM_mu2pT[ibin] - Iaa_CNM_mu12pT[ibin]);

    cout<<"ErrMin: "<<Iaa_CNM_mupT_ErrMin[ibin]<<endl;
    cout<<"ErrMax: "<<Iaa_CNM_mupT_ErrMax[ibin]<<endl;
    minAverage = minAverage + Iaa_CNM_mupT_ErrMin[ibin];
    maxAverage = maxAverage + Iaa_CNM_mupT_ErrMax[ibin];
  }
  cout<<"MinAverage: "<<minAverage/6<<endl;
  cout<<"MaxAverage: "<<maxAverage/6<<endl;
  TGraphAsymmErrors *hIaa_nPDFSyst = new TGraphAsymmErrors(6, assocZtThinner, Iaa_nPDFval, 0, 0, Iaa_nPDF_ErrMin, Iaa_nPDF_ErrMax);
  hIaa_nPDFSyst->SetName("hIaa_nPDFSyst");

  TGraphAsymmErrors *hIaa_CNMwithSyst = new TGraphAsymmErrors(6, assocZtThinner, Iaa_CNM_mu12pT, 0, 0, Iaa_CNM_mupT_ErrMin, Iaa_CNM_mupT_ErrMax);
  hIaa_CNMwithSyst->SetName("hIaa_CNMwithSyst");


  TH1F *grIaaCOLBTmedian[nCen];
  TH1F *grDztPbPbCOLBTmedian[nCen];

  double IaaCOLBT_median0_30[] = {7.328690883313424553e-01, 5.627881746677514396e-01, 4.810049019607843812e-01, 3.768453768453768338e-01, 2.442528735632184256e-01, 2.564102564102563875e-01};
  double IaaCOLBT_erry0_30[] = {3.602547975344214809e-02, 3.195184054870359863e-02, 4.357719132199768669e-02, 4.236834994680931110e-02, 6.343253195931429500e-02, 1.377143988537370689e-01};

  double DztCOLBT_median0_30_AA[] = {2.324579464512946281e+00, 6.916666666666667629e-01, 2.616666666666666030e-01, 8.083333333333335424e-02, 1.416666666666666248e-02, 3.333333333333333981e-03};
  double DztCOLBT_erry0_30_AA[] = {7.586537784494025438e-02, 3.395258131243894528e-02, 2.088327347690277155e-02, 8.207381501496755286e-03, 3.435921354681382989e-03, 1.666666666666666990e-03};

  double IaaCOLBT_median30_50[] = {8.686997953584294496e-01, 7.567127746135068334e-01, 6.862745098039215730e-01, 5.516705516705515677e-01, 4.022988505747126520e-01, 6.410256410256409687e-01};
  double IaaCOLBT_erry30_50[] = {4.096367008096209877e-02, 3.862782907398024101e-02, 5.448224806139326942e-02, 5.341030249549828107e-02, 8.470753180008923355e-02, 2.385283357484778155e-01};

  double DztCOLBT_median30_50_AA[] = {2.223871476117578627e+00, 9.300000000000001599e-01, 3.733333333333331838e-01, 1.183333333333333459e-01, 2.333333333333332399e-02, 8.333333333333334952e-03};
  double DztCOLBT_erry30_50_AA[] = {8.445906306213281367e-02, 3.937003937005906229e-02, 2.494438257849293517e-02, 9.930312739844156592e-03, 4.409585518440983440e-03, 2.635231383473650435e-03};

  TH1F *histIaaCOLBTmedian[3];
  TH1F *histDztPbPbCOLBTmedian[3];

  histIaaCOLBTmedian[0] = new TH1F(Form("histIaaCOLBTmedian_Cen%d_%d", cenBins[0], cenBins[1]), Form("histIaaCOLBTmedian_Cen%d_%d", cenBins[0], cenBins[1]), nAssocCLBT, assocZtCLBT);
  histDztPbPbCOLBTmedian[0] = new TH1F(Form("histDztPbPbCOLBTmedian_Cen%d_%d", cenBins[0], cenBins[1]), Form("histDztPbPbCOLBTmedian_Cen%d_%d", cenBins[0], cenBins[1]), nAssocCLBT, assocZtCLBT);

  histIaaCOLBTmedian[1] = new TH1F(Form("histIaaCOLBTmedian_Cen%d_%d", cenBins[1], cenBins[2]), Form("histIaaCOLBTmedian_Cen%d_%d", cenBins[1], cenBins[2]), nAssocCLBT, assocZtCLBT);
  histDztPbPbCOLBTmedian[1] = new TH1F(Form("histDztPbPbCOLBTmedian_Cen%d_%d", cenBins[1], cenBins[2]), Form("histDztPbPbCOLBTmedian_Cen%d_%d", cenBins[1], cenBins[2]), nAssocCLBT, assocZtCLBT);

  //histIaaCOLBTmedian[2] = new TH1F(Form("histIaaCOLBTmedian_Cen%d_%d", cenBins[2], cenBins[3]), Form("histIaaCOLBTmedian_Cen%d_%d", cenBins[2], cenBins[3]), nAssocCLBT, assocZtCLBT);
  //histDztPbPbCOLBTmedian[2] = new TH1F(Form("histDztPbPbCOLBTmedian_Cen%d_%d", cenBins[2], cenBins[3]), Form("histDztPbPbCOLBTmedian_Cen%d_%d", cenBins[2], cenBins[3]), nAssocCLBT, assocZtCLBT);
  for (int ibin = 0; ibin < nAssocCLBT; ibin++)
  {
    // Iaa
    // 0-30%
    histIaaCOLBTmedian[0]->SetBinContent(ibin + 1, IaaCOLBT_median0_30[ibin]);
    histIaaCOLBTmedian[0]->SetBinError(ibin + 1, IaaCOLBT_erry0_30[ibin]);
    // 30-50%
    histIaaCOLBTmedian[1]->SetBinContent(ibin + 1, IaaCOLBT_median30_50[ibin]);
    histIaaCOLBTmedian[1]->SetBinError(ibin + 1, IaaCOLBT_erry30_50[ibin]);

    // DzT
    // 0-30%
    histDztPbPbCOLBTmedian[0]->SetBinContent(ibin + 1, DztCOLBT_median0_30_AA[ibin]);
    histDztPbPbCOLBTmedian[0]->SetBinError(ibin + 1, DztCOLBT_erry0_30_AA[ibin]);
    // 30-50%
    histDztPbPbCOLBTmedian[1]->SetBinContent(ibin + 1, DztCOLBT_median30_50_AA[ibin]);
    histDztPbPbCOLBTmedian[1]->SetBinError(ibin + 1, DztCOLBT_erry30_50_AA[ibin]);
  }

  TString shshString[2] = {"0.10-0.30", "0.40-1.00"};
  TString sPtAll = Form("_Pt18_40");

  // Getter zt distributions
  TFile *fPlot[nCen];
  TH1F *hZtCent[nCen];    // zt data
  TH1F *hZt_MC_Gen[nCen]; // MC Gen pp
  TH1F *hZt_MC_Rec[nCen]; // MC Rec pp

  TH1F *histRatioPYTHIA_NLOPbPb[3];
  TH1F *histRatioPYTHIACoLBt[2];
  TCanvas *cRatioPYTHIACoLBt[2];
  TCanvas *cPYTHIA_CoLBT[2];

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    cout << "Getter zt distributions: " << sCent << endl;
    //fPlot[iCen] = new TFile(Form("~/work/histogram/FromScratch/ResultsNewMixMC_ZtMergedMoreQM/fPlot%s%s%s.root", shshString[1].Data(), sCent.Data(), sPtAll.Data()));
    fPlot[iCen] = new TFile(Form("Output_checkCode/fPlot%s%s%s.root", shshString[1].Data(), sCent.Data(), sPtAll.Data()));

    hZt_MC_Gen[iCen] = (TH1F *)fPlot[iCen]->Get(Form("hZtMCGenIso1Photon%s%s", sCent.Data(), sPtAll.Data()));
    hZt_MC_Rec[iCen] = (TH1F *)fPlot[iCen]->Get(Form("hZtMCRecIso1Photon%s%s", sCent.Data(), sPtAll.Data()));
    hZtCent[iCen] = (TH1F *)fPlot[iCen]->Get(Form("hZtEffCorrIso1Photon%s%s", sCent.Data(), sPtAll.Data()));
  }

  histRatioPYTHIA_NLOPbPb[0] = new TH1F(Form("histRatioPYTHIA_NLOPbPb_Cen%d_%d", cenBins[0], cenBins[1]), Form("histRatioPYTHIA_NLOPbPb_Cen%d_%d", cenBins[0], cenBins[1]), nAssocCLBT, assocZtCLBT);
  histRatioPYTHIA_NLOPbPb[1] = new TH1F(Form("histRatioPYTHIA_NLOPbPb_Cen%d_%d", cenBins[1], cenBins[2]), Form("histRatioPYTHIA_NLOPbPb_Cen%d_%d", cenBins[1], cenBins[2]), nZtBin, assocZt - 1);
  histRatioPYTHIA_NLOPbPb[2] = new TH1F(Form("histRatioPYTHIA_NLOPbPb_Cen%d_%d", cenBins[2], cenBins[3]), Form("histRatioPYTHIA_NLOPbPb_Cen%d_%d", cenBins[2], cenBins[3]), nZtBin, assocZt - 1);

  for (int ibin = 0; ibin < nZtBin; ibin++)
  {
    // DzT
    // 0-30%
    double binPPNLO = Dzt_median0_30[ibin];
    double binpythia = hZt_MC_Rec[0]->GetBinContent(ibin + 2);
    cout << "pqcd nlo pbpb" << endl;
    cout << binPPNLO << "-pythia" << binpythia << endl;
    cout << binpythia / binPPNLO << endl;
    histRatioPYTHIA_NLOPbPb[0]->SetBinContent(ibin + 1, binpythia / binPPNLO);
  }
  for (int ibin = 0; ibin < nZtBin - 1; ibin++)
  {
    // 30-50%
    double binPPNLO = Dzt_median30_50[ibin];
    double binpythia = hZt_MC_Rec[1]->GetBinContent(ibin + 2);
    cout << "pqcd nlo pbpb" << endl;
    cout << binPPNLO << "-pythia" << binpythia << endl;
    cout << binpythia / binPPNLO << endl;
    histRatioPYTHIA_NLOPbPb[1]->SetBinContent(ibin + 1, binpythia / binPPNLO);
  }

  for (int ibin = 0; ibin < nZtBin - 1; ibin++)
  {
    // 50-90%
    double binPPNLO = Dzt_median50_90[ibin];
    double binpythia = hZt_MC_Rec[2]->GetBinContent(ibin + 2);
    cout << "pqcd nlo pbpb" << endl;
    cout << binPPNLO << "-pythia" << binpythia << endl;
    cout << binpythia / binPPNLO << endl;
    histRatioPYTHIA_NLOPbPb[2]->SetBinContent(ibin + 1, binpythia / binPPNLO);
  }
  TCanvas *cPythiaNLOPbPb0_30 = new TCanvas("cPythiaNLOPbPb0_30", "cPythiaNLOPbPb0_30", 800, 600);
  histRatioPYTHIA_NLOPbPb[0]->SetLineColor(kBlue + 1);
  histRatioPYTHIA_NLOPbPb[0]->SetLineWidth(2);
  histRatioPYTHIA_NLOPbPb[0]->Draw("hist");
  cPythiaNLOPbPb0_30->Print("cPythiaNLOPbPb0_30.pdf");
  TCanvas *cPythiaNLOPbPb30_50 = new TCanvas("cPythiaNLOPbPb30_50", "cPythiaNLOPbPb30_50", 800, 600);
  histRatioPYTHIA_NLOPbPb[1]->Draw("hist");
  histRatioPYTHIA_NLOPbPb[1]->SetLineColor(kBlue + 1);
  histRatioPYTHIA_NLOPbPb[1]->SetLineWidth(2);
  cPythiaNLOPbPb30_50->Print("cPythiaNLOPbPb30_50.pdf");

  TCanvas *cPythiaNLOPbPb50_90 = new TCanvas("cPythiaNLOPbPb50_90", "cPythiaNLOPbPb50_90", 800, 600);
  histRatioPYTHIA_NLOPbPb[2]->Draw("hist");
  histRatioPYTHIA_NLOPbPb[2]->SetLineColor(kBlue + 1);
  histRatioPYTHIA_NLOPbPb[2]->SetLineWidth(2);
  cPythiaNLOPbPb50_90->Print("cPythiaNLOPbPb50_90.pdf");

  histRatioPYTHIACoLBt[0] = new TH1F(Form("histRatioPYTHIACoLBt_Cen%d_%d", cenBins[0], cenBins[1]), Form("histRatioPYTHIACoLBt_Cen%d_%d", cenBins[0], cenBins[1]), nAssocCLBT, assocZtCLBT);
  histRatioPYTHIACoLBt[1] = new TH1F(Form("histRatioPYTHIACoLBt_Cen%d_%d", cenBins[1], cenBins[2]), Form("histRatioPYTHIACoLBt_Cen%d_%d", cenBins[1], cenBins[2]), nAssocCLBT, assocZtCLBT - 1);
  // PYTHIA RATIO PYTHIA COLBT

  for (int ibin = 0; ibin < nAssocCLBT; ibin++)
  {
    // DzT
    // 0-30%
    cout << "COLBT pbpb" << endl;
    double binCOLBT = histDztPbPbCOLBTmedian[0]->GetBinContent(ibin + 1);
    double binpythia = hZt_MC_Rec[0]->GetBinContent(ibin + 2);
    cout << binCOLBT << "-pythia" << binpythia << endl;
    cout << binpythia / binCOLBT << endl;
    histRatioPYTHIACoLBt[0]->SetBinContent(ibin + 1, binpythia / binCOLBT);
  }
  for (int ibin = 0; ibin < nAssocCLBT - 1; ibin++)
  {
    // 30-50%
    double binCOLBT = histDztPbPbCOLBTmedian[1]->GetBinContent(ibin + 1);
    double binpythia = hZt_MC_Rec[1]->GetBinContent(ibin + 2);
    cout << binCOLBT << "-pythia" << binpythia << endl;
    cout << binpythia / binCOLBT << endl;
    histRatioPYTHIACoLBt[1]->SetBinContent(ibin + 1, binpythia / binCOLBT);
  }
  
  TCanvas *cPythiaCOLBT0_30 = new TCanvas("cPythiaCOLBT0_30", "cPythiaCOLBT0_30", 800, 600);
  histRatioPYTHIACoLBt[0]->SetLineColor(kBlue + 1);
  histRatioPYTHIACoLBt[0]->SetLineWidth(2);
  histRatioPYTHIACoLBt[0]->Draw("hist");
  cPythiaCOLBT0_30->Print("cPythiaCOLBT0_30.pdf");
  
  TCanvas *cPythiaCOLBT30_50 = new TCanvas("cPythiaCOLBT30_50", "cPythiaCOLBT30_50", 800, 600);
  histRatioPYTHIACoLBt[1]->Draw("hist");
  histRatioPYTHIACoLBt[1]->SetLineColor(kBlue + 1);
  histRatioPYTHIACoLBt[1]->SetLineWidth(2);
  cPythiaCOLBT30_50->Print("cPythiaCOLBT30_50.pdf");

  new TCanvas();
  // grIaaNLOmedian[0]->SetMarkerColor(kBlack);
  // grIaaNLOmedian[0]->SetFillColor(kRed);
  // grIaaNLOmedian[0]->SetFillStyle(3001);
  // grIaaNLOmedian[0]->Draw("al3");
  grDztNLOmedianpp->Draw("hist X0 same L");
  grIaaNLOmedian[0]->SetName(Form("grIaaNLOmedian_Cen%d_%d", cenBins[0], cenBins[1]));
  grIaaNLOmedian[1]->SetName(Form("grIaaNLOmedian_Cen%d_%d", cenBins[1], cenBins[2]));
  grIaaNLOmedian[2]->SetName(Form("grIaaNLOmedian_Cen%d_%d", cenBins[2], cenBins[3]));

  grDztNLOmedianpp->SetName(Form("grDztNLOmedian_pp"));
  grDztNLOmedianppSyst->SetName(Form("grDztNLOmedian_ppSyst"));
  grDztNLOmedian[0]->SetName(Form("grDztNLOmedian_Cen%d_%d", cenBins[0], cenBins[1]));
  grDztNLOmedian[1]->SetName(Form("grDztNLOmedian_Cen%d_%d", cenBins[1], cenBins[2]));
  grDztNLOmedian[2]->SetName(Form("grDztNLOmedian_Cen%d_%d", cenBins[2], cenBins[3]));

  fileNLO->cd();
  grIaaNLOmedian[0]->Write();
  grIaaNLOmedian[1]->Write();
  grIaaNLOmedian[2]->Write();

  grDztNLOmedianpp->Write();
  grDztNLOmedianppSyst->Write();
  hIaa_nPDF->Write();
  hIaa_nPDFSyst->Write();
  hIaa_CNMwithSyst->Write();
  grDztNLOmedian[0]->Write();
  grDztNLOmedian[1]->Write();
  grDztNLOmedian[2]->Write();

  histIaaCOLBTmedian[0]->Write();
  histDztPbPbCOLBTmedian[0]->Write();
  histIaaCOLBTmedian[1]->Write();
  histDztPbPbCOLBTmedian[1]->Write();

  TH1F *hsystCOLBT030 = (TH1F *)histDztPbPbCOLBTmedian[0]->Clone("hsystCOLBT030");

}
static void ScaleBinBySize(TH1F *h)
{
  for (Int_t ibin = 1; ibin <= h->GetNbinsX(); ibin++)
  {
    Double_t width = h->GetBinWidth(ibin);
    Double_t content = h->GetBinContent(ibin);
    Double_t error = h->GetBinError(ibin);
    printf("bin %d, width %f, content %e\n", ibin, width, content);
    cout << h->GetNbinsX() << endl;
    h->SetBinContent(ibin, content / width);
    h->SetBinError(ibin, error / width);
  }
}
