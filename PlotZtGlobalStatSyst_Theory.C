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

Int_t kMarkCen[] = {21, 20, 71, 25};
// Int_t kColorMark[] = {kAzure + 2, kOrange + 8, kViolet + 7, kCyan - 2};
// Int_t kColorMarkFill[] = {kAzure + 6, kOrange + 7, kViolet + 6, kCyan - 2};
// Int_t kColorMark[] = {kRed-4, kCyan +4, kViolet + 7,kAzure + 2, };
// Int_t kColorMarkFill[] = {kOrange+8,  kCyan +3, kViolet + 6,kAzure + 6,kCyan - 2};
Int_t kColorMark[] = {kRed - 4, kViolet + 7, kAzure - 3, kAzure + 2,};
Int_t kColorMarkFill[] = {kOrange + 8, kViolet + 6, kAzure + 2, kAzure + 6, kCyan - 2};
Int_t nAssoc = 7;

int nZtBinThin = 28;
double assocZtThinner[] = {0., 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.18, 0.20, 0.21, 0.22, 0.23, 0.24, 0.28, 0.30, 0.40, 0.60, 0.80, 1.00, 1.05};
TH1F *shiftBinXinHisto(TH1F *h1, TH1F *h1Shifted, int iCen);

void PlotZtGlobalStatSyst_Theory(float ptMin = 18, float ptMax = 40, bool Mirror = true, TString sMixed = "Mixed", TString shshBkg = "0.40-1.00", TString dirPlot = "zTresults_PbPb_TheorycheckCode", TString dirInputRef = "Output_checkCode", bool b0_30 = true)
{

  TString sMirror;
  if (Mirror)
  {
    sMirror = "Mirror";
  }
  else
  {
    sMirror = "NoMirror";
  }
  // Define centrality bins and systematics directory
  Int_t nCen;
  std::vector<Int_t> cenBins;
  TString dirSyst;
  if (b0_30)
  {
    nCen = 3;
    cenBins.push_back(0);
    cenBins.push_back(30);
    cenBins.push_back(50);
    cenBins.push_back(90);
    dirSyst = "Systematics_checkCode0_30";
  }
  else if (!b0_30)
  {
    nCen = 4;
    cenBins.push_back(0);
    cenBins.push_back(10);
    cenBins.push_back(30);
    cenBins.push_back(50);
    cenBins.push_back(90);
    dirSyst = "Systematics_checkCode";
  }

  TString shshString[2] = {"0.10-0.30", shshBkg};
  TString sPtAll = Form("_Pt%2.0f_%2.0f", ptMin, ptMax);

  // Define where to save the results
  TString processline = Form(".! mkdir -pv %s", dirPlot.Data());
  gROOT->ProcessLine(processline.Data());

  // Getter zT distributions
  TFile *fPlot[nCen]; // root file with data and MC

  TH1F *hZtCent[nCen];    // zT data
  TH1F *hZt_MC_Gen[nCen]; // MC Gen pp
  TH1F *hZt_MC_Rec[nCen]; // MC Rec pp
  TH1F *h3[nCen];

  // Getter pp-pPb data
  TH1F *hDztpp, *hDztpPb, *hIpPbpp;
  TH1F *hDztpp_stat, *hDztpPb_stat, *hIpPbpp_stat;
  TH1F *hDztpp_systErrGet, *hDztpPb_systErrGet, *hIpPbpp_systErrGet;
  TH1F *hDztpp_syst, *hDztpPb_syst, *hIpPbpp_syst;

  TFile *fileNLO = new TFile("RootFiles/fileNLOtest.root ");

  // Getter NLO pQCD calculations
  TGraphAsymmErrors *grIaaNLOmedian[nCen];
  TH1F *grDztNLOmedianpp;
  //TH1F *grDztNLOmedianppSyst;
  TGraphAsymmErrors *grDztNLOmedianppSyst;
  TGraphAsymmErrors *grDztNLOmedian[nCen];
  TGraphAsymmErrors *grIcpNLOmedian[nCen];
  // Getter COLBT calculations
  TH1F *histIaaCOLBTmedian[nCen];
  TH1F *histDztPbPbCOLBTmedian[nCen];
  TH1F *histIaaCOLBTmedianSyst[2];
  TH1F *histDztPbPbCOLBTmedianSyst[2];
  // Test PYTHIA COLBT ratio
  TH1F *histRatioPYTHIACoLBt[2];
  TCanvas *cRatioPYTHIACoLBt[2];
  TCanvas *cPYTHIA_CoLBT[2];

  cout << "Getter zt distributions: data and models" << endl;

  ///////////////////////////////////////////
  /////////// data: Dzt pPb and pp /////////
  /////////////////////////////////////////

  TFile *finputpp_pPbPaper_DztpPb = new TFile("TheoryCalculations/HEPData-ins1798523-v1-root.root");
  TDirectory *dirPaper_DztpPb = (TDirectory *)finputpp_pPbPaper_DztpPb->Get("Figure 5 Top Panel");

  hDztpPb = (TH1F *)dirPaper_DztpPb->Get("Hist1D_y1");
  hDztpPb_stat = (TH1F *)dirPaper_DztpPb->Get("Hist1D_y1_e1");
  hDztpPb_systErrGet = (TH1F *)dirPaper_DztpPb->Get("Hist1D_y1_e2");
  hDztpPb_syst = (TH1F *)hDztpPb->Clone("hDztpPb_syst");
  for (int ibin = 0; ibin < hDztpPb->GetNbinsX(); ibin++)
  {
    hDztpPb->SetBinError(ibin + 1, hDztpPb_stat->GetBinContent(ibin + 1));
    hDztpPb_syst->SetBinError(ibin + 1, hDztpPb_systErrGet->GetBinContent(ibin + 1));
  }
  cout << "Address DztpPb " << endl;
  cout << hDztpPb_stat << endl;
  cout << hDztpPb_syst << endl;
  PlotStyle(hDztpPb, 24, 1, kBlack, kBlack, "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
  PlotStyle(hDztpPb_syst, 24, 1, kBlack, kGray + 1, "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", true);

  hDztpp = (TH1F *)dirPaper_DztpPb->Get("Hist1D_y2");
  hDztpp_stat = (TH1F *)dirPaper_DztpPb->Get("Hist1D_y2_e1");
  hDztpp_systErrGet = (TH1F *)dirPaper_DztpPb->Get("Hist1D_y2_e2");
  hDztpp_syst = (TH1F *)hDztpp->Clone("hDztpp_syst");
  for (int ibin = 0; ibin < hDztpp->GetNbinsX(); ibin++)
  {
    hDztpp->SetBinError(ibin + 1, hDztpp_stat->GetBinContent(ibin + 1));
    hDztpp_syst->SetBinError(ibin + 1, hDztpp_systErrGet->GetBinContent(ibin + 1));
  }
  cout << "Address Dztpp " << endl;
  cout << hDztpp_stat << endl;
  cout << hDztpp_syst << endl;
  PlotStyle(hDztpp, 20, 1, kBlack, kBlack, "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
  PlotStyle(hDztpp_syst, 20, 1, kBlack, kBlack, "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);

  ////////////////////////////////////////////
  /////////// data: IpA, pPb and pp /////////
  //////////////////////////////////////////
  TFile *finputpp_pPbPaper_IpPbpp = new TFile("TheoryCalculations/HEPData-ins1798523-v1-Figure_5_Bottom_Panel.root");
  TDirectory *dirPaper_IpPbpp = (TDirectory *)finputpp_pPbPaper_IpPbpp->Get("Figure 5 Bottom Panel");
  hIpPbpp = (TH1F *)dirPaper_IpPbpp->Get("Hist1D_y1");
  hIpPbpp_stat = (TH1F *)dirPaper_IpPbpp->Get("Hist1D_y1_e1");
  hIpPbpp_systErrGet = (TH1F *)dirPaper_IpPbpp->Get("Hist1D_y1_e2");
  hIpPbpp_syst = (TH1F *)hIpPbpp->Clone("hIpPbpp_syst");
  for (int ibin = 0; ibin < hIpPbpp->GetNbinsX(); ibin++)
  {
    hIpPbpp->SetBinError(ibin + 1, hIpPbpp_stat->GetBinContent(ibin + 1));
    hIpPbpp_syst->SetBinError(ibin + 1, hIpPbpp_systErrGet->GetBinContent(ibin + 1));
  }
  cout << "Address Iaa = pPb/pp" << endl;
  cout << hIpPbpp_stat << endl;
  cout << hIpPbpp_syst << endl;
  PlotStyle(hIpPbpp, 72, 1, kBlack, kBlack, "#it{z}_{T}", "p-Pb/pp", false);
  PlotStyle(hIpPbpp_syst, 72, 1, kBlack, kBlack, "#it{z}_{T}", "p-Pb/pp", false);

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);

    /////////////////////////////////
    /////////// data: PbPb /////////
    ///////////////////////////////
    cout << "Get Data and MC and set Plot style" << sCent << endl;
    fPlot[iCen] = new TFile(Form("%s/fPlot%s%s%s.root", dirInputRef.Data(), shshString[1].Data(), sCent.Data(), sPtAll.Data()));
    hZt_MC_Gen[iCen] = (TH1F *)fPlot[iCen]->Get(Form("hZtMCGenIso1Photon%s%s", sCent.Data(), sPtAll.Data()));
    hZt_MC_Rec[iCen] = (TH1F *)fPlot[iCen]->Get(Form("hZtMCRecIso1Photon%s%s", sCent.Data(), sPtAll.Data()));
    hZtCent[iCen] = (TH1F *)fPlot[iCen]->Get(Form("hZtEffCorrIso1Photon%s%s", sCent.Data(), sPtAll.Data()));

    PlotStyle(hZt_MC_Gen[iCen], 72, 1, kBlack, kBlack, "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
    PlotStyle(hZt_MC_Rec[iCen], 21, 1, kOrange + 7, kOrange + 7, "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
    PlotStyle(hZtCent[iCen], kMarkCen[iCen], 1, kColorMark[iCen], kColorMarkFill[iCen], "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
    h3[iCen] = (TH1F *)hZtCent[iCen]->Clone(Form("h3%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    h3[iCen]->Divide(hZt_MC_Gen[iCen]);

    /////////////////////////////////////////////////
    /////////// theory: NLO + qhat (eLoss) /////////
    ///////////////////////////////////////////////

    cout << "Get NLO calclulations and set Plot style: " << sCent << endl;
    grIaaNLOmedian[iCen] = (TGraphAsymmErrors *)fileNLO->Get(Form("grIaaNLOmedian%s", sCent.Data()));
    cout << grIaaNLOmedian[iCen]->GetName() << endl;
    grIaaNLOmedian[iCen]->SetLineWidth(8);
    grIaaNLOmedian[iCen]->SetLineColor(kPink + 4);
    grIaaNLOmedian[iCen]->SetFillColorAlpha(kPink - 4, 0.25);
    grIaaNLOmedian[iCen]->SetFillStyle(1001);

    grDztNLOmedian[iCen] = (TGraphAsymmErrors *)fileNLO->Get(Form("grDztNLOmedian%s", sCent.Data()));
    cout << grDztNLOmedian[iCen]->GetName() << endl;
    grDztNLOmedian[iCen]->SetLineWidth(8);
    grDztNLOmedian[iCen]->SetLineColor(kPink + 4);
    grDztNLOmedian[iCen]->SetFillColorAlpha(kPink - 4, 0.25);
    grDztNLOmedian[iCen]->SetFillStyle(1001);

    grIcpNLOmedian[iCen] = (TGraphAsymmErrors *)fileNLO->Get(Form("grDztNLOmedian%s", sCent.Data()));
    cout << grDztNLOmedian[iCen]->GetName() << endl;
    grIcpNLOmedian[iCen]->SetLineWidth(8);
    grIcpNLOmedian[iCen]->SetFillStyle(1001);
    grIcpNLOmedian[iCen]->SetLineColor(kPink + 4);
    grIcpNLOmedian[iCen]->SetFillColorAlpha(kMagenta - 9, 0.40);
  }
  
  // Set x error, in case it helps in plotting, it does not seem
  Float_t widthPqcd[] = {0.025, 0.025, 0.05, 0.05, 0.10, 0.10};
  for (Int_t iCen = 0; iCen < nCen; iCen++)
  {
    cout << iCen << endl;
    Int_t nPoints = grIaaNLOmedian[iCen]->GetN();
    for (Int_t ibin = 0; ibin < nPoints; ibin++)
    {
      grIaaNLOmedian[iCen]->SetPointEXhigh(ibin, widthPqcd[ibin]);
      grIaaNLOmedian[iCen]->SetPointEXlow(ibin, widthPqcd[ibin]);
      grDztNLOmedian[iCen]->SetPointEXhigh(ibin, widthPqcd[ibin]);
      grDztNLOmedian[iCen]->SetPointEXlow(ibin, widthPqcd[ibin]);
      grIcpNLOmedian[iCen]->SetPointEXhigh(ibin, widthPqcd[ibin]);
      grIcpNLOmedian[iCen]->SetPointEXlow(ibin, widthPqcd[ibin]);
    }
  }
  ////////////////////////////////////
  /////////// theory: COLBt /////////
  //////////////////////////////////
  for (int iCen = 0; iCen < 2; iCen++)
  {
    histIaaCOLBTmedian[iCen] = (TH1F *)fileNLO->Get(Form("histIaaCOLBTmedian_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    histIaaCOLBTmedian[iCen]->SetLineWidth(4);
    histIaaCOLBTmedian[iCen]->SetFillStyle(1001);
    histIaaCOLBTmedian[iCen]->SetLineColor(kCyan + 3);
    histIaaCOLBTmedian[iCen]->SetFillColorAlpha(kCyan + 2, 0.20);

    histIaaCOLBTmedianSyst[iCen] = (TH1F *)histIaaCOLBTmedian[iCen]->Clone(Form("histIaaCOLBTmedianSyst_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]));

    histDztPbPbCOLBTmedian[iCen] = (TH1F *)fileNLO->Get(Form("histDztPbPbCOLBTmedian_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    histDztPbPbCOLBTmedian[iCen]->SetLineWidth(3);
    histDztPbPbCOLBTmedian[iCen]->SetFillStyle(1001);
    histDztPbPbCOLBTmedian[iCen]->SetLineColor(kCyan + 3);
    histDztPbPbCOLBTmedian[iCen]->SetFillColorAlpha(kCyan + 2, 0.20);

    histDztPbPbCOLBTmedianSyst[iCen] = (TH1F *)histDztPbPbCOLBTmedian[iCen]->Clone(Form("histDztPbPbCOLBTmedianSyst_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    //////////////////////////////////////////
    //////// Test ratio PYTHIA/COLBT ////////
    ////////////////////////////////////////
    histRatioPYTHIACoLBt[iCen] = (TH1F *)hZt_MC_Rec[iCen]->Clone(Form("histRatioPYTHIACoLBt_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    histRatioPYTHIACoLBt[iCen]->Divide(histDztPbPbCOLBTmedian[iCen]);
    cRatioPYTHIACoLBt[iCen] = new TCanvas(Form("cRatioPYTHIACoLBtCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("cRatioPYTHIACoLBtCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), 800, 600);

    histRatioPYTHIACoLBt[iCen]->Draw("hist");
    cRatioPYTHIACoLBt[iCen]->Print(dirPlot + Form("/cRatioPYTHIACoLBtCen%d_%d.pdf", cenBins[iCen], cenBins[iCen + 1]));

    cPYTHIA_CoLBT[iCen] = new TCanvas(Form("cPYTHIA_CoLBTCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("cPYTHIA_CoLBTCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), 800, 600);
    histDztPbPbCOLBTmedian[iCen]->Draw("hist same ");
    hZt_MC_Rec[iCen]->Draw("hist same ");
    cPYTHIA_CoLBT[iCen]->Print(dirPlot + Form("/cPYTHIA_CoLBTCen%d_%d.pdf", cenBins[iCen], cenBins[iCen + 1]));
  }

  /////////////////////////////////////////////
  /////////// theory: nPDF only //////////////
  ///////////////////////////////////////////

  TH1F *hIaa_nPDF  = (TH1F *)fileNLO->Get(Form("hIaa_nPDF"));
  TGraphAsymmErrors *hIaa_nPDFSyst  = (TGraphAsymmErrors *)fileNLO->Get(Form("hIaa_nPDFSyst"));
  hIaa_nPDF->SetLineColor(kOrange-3);
  hIaa_nPDF->SetLineWidth(4);
  hIaa_nPDF->SetLineStyle(1);
  hIaa_nPDFSyst->SetLineWidth(0);
  hIaa_nPDFSyst->SetLineStyle(0);
  hIaa_nPDFSyst->SetLineColor(kOrange - 3);
  hIaa_nPDFSyst->SetFillColorAlpha(kOrange-2, 0.25);
  
  TGraphAsymmErrors *hIaa_CNMwithSyst  = (TGraphAsymmErrors *)fileNLO->Get(Form("hIaa_CNMwithSyst"));
  hIaa_CNMwithSyst->SetLineWidth(0);
  hIaa_CNMwithSyst->SetLineStyle(7);
  hIaa_CNMwithSyst->SetLineColor(kOrange-3);
  hIaa_CNMwithSyst->SetFillStyle(1001);
  hIaa_CNMwithSyst->SetFillColorAlpha(kGray+3, 0.45);
  

  /////////////////////////////////////
  /////////// theory: NLO pp /////////
  ///////////////////////////////////

  grDztNLOmedianpp = (TH1F *)fileNLO->Get(Form("grDztNLOmedian_pp"));
  grDztNLOmedianpp->SetLineWidth(8);
  grDztNLOmedianpp->SetLineColor(kSpring-6);
  grDztNLOmedianpp->SetMarkerColor(kSpring-6);
  grDztNLOmedianpp->SetMarkerStyle(1);
  grDztNLOmedianpp->SetLineStyle(1);

  grDztNLOmedianppSyst = (TGraphAsymmErrors *)fileNLO->Get(Form("grDztNLOmedian_ppSyst"));
  grDztNLOmedianppSyst->SetLineWidth(4);
  grDztNLOmedianppSyst->SetLineColor(kSpring-7);
  grDztNLOmedianppSyst->SetMarkerColor(kSpring-7);
  grDztNLOmedianppSyst->SetMarkerStyle(1);
  grDztNLOmedianppSyst->SetLineStyle(1);
  grDztNLOmedianppSyst->SetFillStyle(1001);
  grDztNLOmedianppSyst->SetFillColorAlpha(kSpring-7, 0.3);

  TCanvas *cPP_NLO = new TCanvas("cPP_NLO", "cPP_NLO", 800, 600);
  gPad->SetLogy();
  // grDztNLOmedianpp->SetFillStyle(5000);
  // grDztNLOmedianpp->SetFillColorAlpha(16, 0.4);
  grDztNLOmedianpp->Draw("histp");
  grDztNLOmedianppSyst->Draw("histp");
  hZt_MC_Gen[0]->Draw("hist same");
  cPP_NLO->Print(dirPlot + Form("/Checkpp_NLO.pdf"));

  //////////////////////////////////////////
  //////// Get total systematics //////////
  ////////////////////////////////////////

  cout << "Getter zt distribution systematics" << endl;
  TFile *fSystFile = new TFile(Form("%s/fAllSystFile%s%s%s%s.root", dirSyst.Data(), sMixed.Data(), shshBkg.Data(), sMirror.Data(), sPtAll.Data()));
  cout << dirSyst << endl;
  TH1F *hSystZt[nCen];
  TH1F *hsyst[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    hSystZt[iCen] = (TH1F *)fSystFile->Get(Form("hsystCen%d%d", cenBins[iCen], cenBins[iCen + 1]));
    PlotStyle(hSystZt[iCen], kMarkCen[iCen], 1, kColorMark[iCen], kColorMarkFill[iCen], "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", true);
    hsyst[iCen] = (TH1F *)hSystZt[iCen]->Clone(Form("hsyst%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    hsyst[iCen]->Divide(hZt_MC_Gen[iCen]);
    for (int ibin = 0; ibin < nAssoc; ibin++)
    {
      cout << cenBins[iCen] << "-" << cenBins[iCen + 1] << endl;
      cout << ibin << "-" << ibin + 1 << " systematics: " << hSystZt[iCen]->GetBinContent(ibin + 1) << " pm " << hSystZt[iCen]->GetBinError(ibin + 1) << endl;
    }
  }

  //////////////////////////////////////////////////////////
  //////// Compute Iaa with pp from theory ////////////////
  /////// Different binning between data and theory //////
  ///////////////////////////////////////////////////////

  Int_t nAssocNLO = 6;
  double assocZtNLO[] = {0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 1.0};
  TH1F *hPbPb_NLO[nCen];
  TH1F *hPbPbPYTHIA_NLO[nCen];
  TH1F *hPbPb_NLOSyst[nCen];

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    hPbPb_NLO[iCen] = new TH1F(Form("hPbPb%d_%d_NLO", cenBins[iCen], cenBins[iCen + 1]), Form("hPbPb%d_%d_NLO", cenBins[iCen], cenBins[iCen + 1]), nAssocNLO, assocZtNLO);
    hPbPb_NLOSyst[iCen] = new TH1F(Form("hPbPb%d_%d_NLOSyst", cenBins[iCen], cenBins[iCen + 1]), Form("hPbPb%d_%d_NLOSyst", cenBins[iCen], cenBins[iCen + 1]), nAssocNLO, assocZtNLO);

    hPbPbPYTHIA_NLO[iCen] = new TH1F(Form("hPbPbPYTHIA%d_%d_NLO", cenBins[iCen], cenBins[iCen + 1]), Form("hPbPbPYTHIA%d_%d_NLO", cenBins[iCen], cenBins[iCen + 1]), nAssocNLO, assocZtNLO);
    // hPbPbPYTHIA_NLOSyst[iCen] = new TH1F(Form("hPbPbPYTHIA%d_%d_NLOSyst", cenBins[iCen], cenBins[iCen + 1]), Form("hPbPb%d_%d_NLOSyst", cenBins[iCen], cenBins[iCen + 1]), nAssocNLO, assocZtNLO);
    for (int ibin = 0; ibin < nAssocNLO; ibin++)
    {
      cout << "Cent: " << iCen << endl;
      // cout << grDztNLOmedianpp->GetBinCenter(ibin)<<endl;
      // cout << hPbPb_NLOSyst[iCen]->GetBinCenter(ibin)<<endl;
      // cout << hPbPb_NLO[iCen]->GetBinCenter(ibin ) << "____" << hZtCent[iCen]->GetBinCenter(ibin + 1) << endl;
      cout << hPbPb_NLO[iCen]->GetBinCenter(ibin + 1) << ")" << hPbPb_NLO[iCen]->GetBinContent(ibin + 1) << "____" << grDztNLOmedianpp->GetBinContent(ibin + 1) << " Ratio: " << hPbPb_NLO[iCen]->GetBinContent(ibin + 1) / grDztNLOmedianpp->GetBinContent(ibin + 1) << endl;

      hPbPb_NLO[iCen]->SetBinContent(ibin + 1, hZtCent[iCen]->GetBinContent(ibin + 1));
      hPbPb_NLO[iCen]->SetBinError(ibin + 1, hZtCent[iCen]->GetBinError(ibin + 1));

      cout << hPbPb_NLO[iCen]->GetBinContent(ibin + 1) << "____" << grDztNLOmedianpp->GetBinContent(ibin + 1) << " Ratio: " << hPbPb_NLO[iCen]->GetBinContent(ibin + 1) / grDztNLOmedianpp->GetBinContent(ibin + 1) << endl;
      hPbPbPYTHIA_NLO[iCen]->SetBinContent(ibin + 1, hZt_MC_Rec[iCen]->GetBinContent(ibin + 1));
      hPbPbPYTHIA_NLO[iCen]->SetBinError(ibin + 1, hZt_MC_Rec[iCen]->GetBinError(ibin + 1));
      hPbPb_NLOSyst[iCen]->SetBinContent(ibin + 1, hSystZt[iCen]->GetBinContent(ibin + 1));
      hPbPb_NLOSyst[iCen]->SetBinError(ibin + 1, hSystZt[iCen]->GetBinError(ibin + 1));
      cout << hPbPb_NLOSyst[iCen]->GetBinError(ibin + 1) << endl;
      cout << "ERRORE pp: " << grDztNLOmedianpp->GetBinError(ibin + 1) << endl;
      cout << "Cent: " << iCen << endl;
    }
    cout << "Centrality" << iCen << endl;
    hPbPb_NLO[iCen]->Divide(grDztNLOmedianpp);
    hPbPbPYTHIA_NLO[0]->Divide(histDztPbPbCOLBTmedian[0]); // different number of bins
    hPbPb_NLOSyst[iCen]->Divide(grDztNLOmedianpp);

    PlotStyle(hPbPb_NLO[iCen], kMarkCen[iCen], 1, kColorMark[iCen], kColorMarkFill[iCen], "#it{z}_{T}", "Ratio", false);
    // PlotStyle(hPbPb_NLO[iCen], kMarkCen[iCen], 1, kColorMark[iCen], kColorMarkFill[iCen], "#it{z}_{T}", "#it{I}_{NLO} = Pb#font[122]{-}Pb/NLO pQCD", false);
    PlotStyle(hPbPb_NLOSyst[iCen], kMarkCen[iCen], 1, kColorMark[iCen], kColorMarkFill[iCen], "#it{z}_{T}", "Ratio", true);
    // PlotStyle(hPbPb_NLOSyst[iCen], kMarkCen[iCen], 1, kColorMark[iCen], kColorMarkFill[iCen], "#it{z}_{T}", "#it{I}_{NLO} = Pb#font[122]{-}Pb/NLO pQCD", true);
  }

  hPbPb_NLO[0]->GetXaxis()->SetRangeUser(0., 1.0);
  hPbPb_NLO[1]->GetXaxis()->SetRangeUser(0., 0.6);
  hPbPb_NLO[2]->GetXaxis()->SetRangeUser(0., 0.6);
  hPbPb_NLOSyst[0]->GetXaxis()->SetRangeUser(0., 1.0);
  hPbPb_NLOSyst[1]->GetXaxis()->SetRangeUser(0., 0.6);
  hPbPb_NLOSyst[2]->GetXaxis()->SetRangeUser(0., 0.6);

  hPbPbPYTHIA_NLO[0]->GetXaxis()->SetRangeUser(0., 1.0);
  hPbPbPYTHIA_NLO[1]->GetXaxis()->SetRangeUser(0., 0.6);
  hPbPbPYTHIA_NLO[2]->GetXaxis()->SetRangeUser(0., 0.6);

  ////////////////////////////////////////////////////////
  //////// Plotting final results and theory /////////////
  ///////////////////////////////////////////////////////

  TCanvas *cPbPbPYTHIA_NLORatio[nCen];
  TCanvas *cPbPb_NLORatio[nCen];
  TLatex *latPbPb_NLO[nCen];
  TLegend *legPbPb_NLOratio[nCen];
  TLegend *legRatioPbPb[nCen];
  // TH1F *hGeneralRatio = new TH1F("hGeneralRatio", "hGeneralRatio", nZtBinThin, assocZtThinner);
  TH1F *hGeneralRatio = new TH1F("hGeneralRatio", "hGeneralRatio", 105, 0., 1.05);
  PlotStyle(hGeneralRatio, 20, 1, kWhite, kWhite, " #it{z}_{T} ", " #it{I}_{pQCD} = Pb#font[122]{-}Pb (data) / pp (NLO pQCD)", false);
  TGraph *lineNLO = DrawLine(lineNLO, 0, 0.5, 1.2, 0.5);
  TGraph *lineNLO1 = DrawLine(lineNLO1, 0, 1, 1.2, 1);
  TLegend *legRatioNLO = LegStd(legRatioNLO, 0.12, 0.80, 0.32, 0.88);
  TLegend *legTitleRatio[nCen];
  TLegend *legTitleRatioSelection[nCen];

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    cPbPb_NLORatio[iCen] = canvasStdIaa(Form("cPbPb_NLORatioCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), 1, 1);
    // new TCanvas(Form("cPbPb_NLORatioCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("cPbPb_NLORatioCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), 800, 600);
    cPbPb_NLORatio[iCen]->cd(iCen + 1);
    hGeneralRatio->GetYaxis()->SetRangeUser(0.0, 2.5);
    hGeneralRatio->SetTitle(" ");
    hGeneralRatio->GetYaxis()->SetTitleSize(0.046);
    hGeneralRatio->GetYaxis()->SetTitleOffset(1.);
    hGeneralRatio->GetYaxis()->SetLabelSize(0.045);
    hGeneralRatio->GetXaxis()->SetRangeUser(0.05, 1.05);
    hGeneralRatio->GetXaxis()->SetLabelSize(0.045);
    hGeneralRatio->GetXaxis()->SetNdivisions(505);
    hGeneralRatio->GetXaxis()->SetTitleSize(0.055);
    hGeneralRatio->Draw("histsame");
    hPbPb_NLOSyst[iCen]->Draw("E2Psame ");
    hPbPb_NLO[iCen]->Draw("EPX0same");
    grIaaNLOmedian[iCen]->Draw("pl3 same");
    //hIaa_nPDFSyst->Draw("hist same c e3l");
    hIaa_nPDFSyst->Draw("pl3 same");
    hIaa_CNMwithSyst->Draw("pl3 same");
    hIaa_nPDF->Draw("l same ");
    legTitleRatio[iCen] = LegStd(legTitleRatio[iCen], 0.12, 0.91, 0.12, 0.94);
    legTitleRatio[iCen]->SetHeader("ALICE, Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
    legTitleRatioSelection[iCen] = LegStd(legTitleRatioSelection[iCen], 0.56, 0.84, 0.56, 0.96);
    legTitleRatioSelection[iCen]->SetTextSize(0.036);
    legTitleRatioSelection[iCen]->AddEntry("", "|#Delta#it{#varphi}|#color[0]{..}>#color[0]{..}#frac{3}{5}#color[0]{.}#it{#pi}, |#it{#eta}^{#it{#gamma}}| < 0.67", "");
    legTitleRatioSelection[iCen]->AddEntry("", "18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#it{#gamma}}#color[0]{.}<#color[0]{.}40#color[0]{.}GeV/#it{c}#color[0]{.};#color[0]{..}#it{p}_{T}^{h} > 0.5 GeV/#it{c}", "");
    if (iCen == 0)
    {
      legRatioPbPb[iCen] = LegStd(legRatioPbPb[iCen], 0.550, 0.66, 0.840, 0.82);
      legRatioPbPb[iCen]->SetTextSize(0.036);
      legRatioPbPb[iCen]->AddEntry(hIaa_nPDF, "#it{I}_{AA}, NLO pQCD#color[0]{.}+#color[0]{.}CNM", "l");
      legRatioPbPb[iCen]->AddEntry(hIaa_CNMwithSyst, "0.7#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}#it{#mu}#color[0]{.}<#color[0]{..}2#color[0]{.}#it{p}_{T}^{#gamma}", "f");
      legRatioPbPb[iCen]->AddEntry(hIaa_nPDFSyst, "nPDF uncert. EPPS21", "f");
      // legRatioPbPb[iCen]->AddEntry((TObject *)0, "M. Xie et al.", "");
      legPbPb_NLOratio[iCen] = LegStd(legPbPb_NLOratio[iCen], 0.20, 0.66, 0.50, 0.82);
      legPbPb_NLOratio[iCen]->SetTextSize(0.036);
      legPbPb_NLOratio[iCen]->AddEntry(hPbPb_NLOSyst[iCen], Form("#it{I}_{pQCD} %d#font[122]{-}%d%%", cenBins[iCen], cenBins[iCen + 1]), "efp");
      legPbPb_NLOratio[iCen]->AddEntry(grDztNLOmedian[0], "#it{I}_{AA}, NLO pQCD#color[0]{.}+#color[0]{.}#Delta#it{E}_{loss}", "lf");
      // legPbPb_NLOratio[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      //  legPbPb_NLOratio[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
      histIaaCOLBTmedianSyst[0]->SetLineStyle(4);
      histIaaCOLBTmedianSyst[0]->Draw("X0sameE3 ");
      histIaaCOLBTmedian[0]->SetFillStyle(0);
      histIaaCOLBTmedian[0]->SetLineStyle(4);
      histIaaCOLBTmedian[0]->Draw("histsameL ");
      legPbPb_NLOratio[iCen]->AddEntry(histIaaCOLBTmedianSyst[0], "#it{I}_{AA}, CoLBT-hydro", "lf");
      // legPbPb_NLOratio[iCen]->AddEntry((TObject *)0, "X. N. Wang et al.", "");
    }
    else if (iCen == 1)
    {
      legRatioPbPb[iCen] = LegStd(legRatioPbPb[iCen], 0.550, 0.66, 0.840, 0.82);
      legRatioPbPb[iCen]->SetTextSize(0.036);
      legRatioPbPb[iCen]->AddEntry(hIaa_nPDF, "#it{I}_{AA}, NLO pQCD#color[0]{.}+#color[0]{.}CNM", "l");
      legRatioPbPb[iCen]->AddEntry(hIaa_CNMwithSyst, "0.7#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}#it{#mu}#color[0]{.}<#color[0]{..}2#color[0]{.}#it{p}_{T}^{#gamma}", "f");
      legRatioPbPb[iCen]->AddEntry(hIaa_nPDFSyst, "nPDF uncert. EPPS21", "f");
      // legRatioPbPb[iCen]->AddEntry((TObject *)0, "M. Xie et al.", "");
      legPbPb_NLOratio[iCen] = LegStd(legPbPb_NLOratio[iCen], 0.20, 0.66, 0.50, 0.82);
      legPbPb_NLOratio[iCen]->SetTextSize(0.036);
      legPbPb_NLOratio[iCen]->AddEntry(hPbPb_NLOSyst[iCen], Form("#it{I}_{pQCD} %d#font[122]{-}%d%%", cenBins[iCen], cenBins[iCen + 1]), "efp");
      legPbPb_NLOratio[iCen]->AddEntry(grDztNLOmedian[0], "#it{I}_{AA}, NLO pQCD#color[0]{.}+#color[0]{.}#Delta#it{E}_{loss}", "lf");
      // legPbPb_NLOratio[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      //  legPbPb_NLOratio[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
      histIaaCOLBTmedianSyst[1]->SetLineStyle(4);
      histIaaCOLBTmedianSyst[1]->Draw("X0sameE3 ");
      histIaaCOLBTmedian[1]->SetFillStyle(0);
      histIaaCOLBTmedian[1]->SetLineStyle(4);
      histIaaCOLBTmedian[1]->Draw("histsameL ");
      legPbPb_NLOratio[iCen]->AddEntry(histIaaCOLBTmedianSyst[1], "#it{I}_{AA}, CoLBT-hydro", "lf");
      // legPbPb_NLOratio[iCen]->AddEntry((TObject *)0, "X. N. Wang et al.", "");
    }
    else if (iCen == 2)
    {
      legRatioPbPb[iCen] = LegStd(legRatioPbPb[iCen], 0.550, 0.66, 0.840, 0.82);
      legRatioPbPb[iCen]->SetTextSize(0.036);
      legRatioPbPb[iCen]->AddEntry(hIaa_nPDF, "#it{I}_{AA}, NLO pQCD#color[0]{.}+#color[0]{.}CNM", "l");
      legRatioPbPb[iCen]->AddEntry(hIaa_CNMwithSyst, "0.7#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}#it{#mu}#color[0]{.}<#color[0]{..}2#color[0]{.}#it{p}_{T}^{#gamma}", "f");
      legRatioPbPb[iCen]->AddEntry(hIaa_nPDFSyst, "nPDF uncert. EPPS21", "f");
      // legRatioPbPb[iCen]->AddEntry((TObject *)0, "M. Xie et al.", "");
      legPbPb_NLOratio[iCen] = LegStd(legPbPb_NLOratio[iCen], 0.20, 0.715, 0.50, 0.82);
      legPbPb_NLOratio[iCen]->SetTextSize(0.036);
      legPbPb_NLOratio[iCen]->AddEntry(hPbPb_NLOSyst[iCen], Form("#it{I}_{pQCD} %d#font[122]{-}%d%%", cenBins[iCen], cenBins[iCen + 1]), "efp");
      legPbPb_NLOratio[iCen]->AddEntry(grDztNLOmedian[0], "#it{I}_{AA}, NLO pQCD#color[0]{.}+#color[0]{.}#Delta#it{E}_{loss}", "lf");
      // legPbPb_NLOratio[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      //  legPbPb_NLOratio[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
    }
    // latPbPb_NLO[iCen] = LatexStdISO(latPbPb_NLO[iCen], 0.440, 0.92, 0.04, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true);
    legRatioPbPb[iCen]->Draw("same");
    legPbPb_NLOratio[iCen]->Draw("same");
    lineNLO->Draw("l");
    lineNLO1->Draw("l");
    legRatioNLO->Draw("same");
    legTitleRatio[iCen]->Draw("same");
    legTitleRatioSelection[iCen]->Draw("same");
    cPbPb_NLORatio[iCen]->Print(dirPlot + Form("/Ipqcd_Cen%d_%d.pdf", cenBins[iCen], cenBins[iCen + 1]));
    cPbPbPYTHIA_NLORatio[iCen] = canvasStd(Form("cPbPbPYTHIA_NLORatioCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), 1, 1);
    hPbPbPYTHIA_NLO[iCen]->Draw("HIST");
    cPbPbPYTHIA_NLORatio[iCen]->Print(dirPlot + Form("/RatioPYTHIANLOCen%d_%d.pdf", cenBins[iCen], cenBins[iCen + 1]));
  }

  /////////////////////////////////////////////
  // -------- Ratio 1 canvas 3 pad --------- //
  /////////////////////////////////////////////

  TCanvas *cPbPb_NLORatioTogetherTest = new TCanvas("cPbPb_NLORatioTogetherTest", "", 2000, 800);
  TLegend *legPbPb_NLOratioIaaCNMonly[nCen];
  // Number of PADS
  const Int_t Nx = 3;
  const Int_t Ny = 1;

  // Margins
  Float_t lMargin = 0.08;
  Float_t rMargin = 0.03;
  Float_t bMargin = 0.15;
  Float_t tMargin = 0.03;

  // Canvas setup
  CanvasPartition(cPbPb_NLORatioTogetherTest, Nx, Ny, lMargin, rMargin, bMargin, tMargin);

  TH1F *hIaaXplot[nCen];
  TH1F *hIaaXplotSyst[nCen];

  TPad *pad[Nx][Ny];

  for (Int_t iCen = 0; iCen < Nx; iCen++)
  {
    for (Int_t j = 0; j < Ny; j++)
    {
      cPbPb_NLORatioTogetherTest->cd(0);

      // Get the pads previously created.
      pad[iCen][j] = (TPad *)cPbPb_NLORatioTogetherTest->FindObject(TString::Format("pad_%d_%d", iCen, j).Data());
      pad[iCen][j]->Draw();
      pad[iCen][j]->SetFillStyle(4000);
      pad[iCen][j]->SetFrameFillStyle(4000);
      pad[iCen][j]->cd();

      // Size factors
      Float_t xFactor = pad[0][0]->GetAbsWNDC() / pad[iCen][j]->GetAbsWNDC();
      Float_t yFactor = pad[0][0]->GetAbsHNDC() / pad[iCen][j]->GetAbsHNDC();
      TH1F *hFrame = (TH1F *)hGeneralRatio->Clone(TString::Format("h_%d_%d", iCen, j).Data());

      hFrame->GetYaxis()->SetRangeUser(-0.50, 2.05);

      // Format for y axis
      // hFrame->GetYaxis()->SetLabelFont(43);
      hFrame->GetYaxis()->SetLabelSize(0.045);
      hFrame->GetYaxis()->SetLabelOffset(0.02);
      // hFrame->GetYaxis()->SetTitleFont(43);
      hFrame->GetYaxis()->SetTitleSize(0.05);
      hFrame->GetYaxis()->SetTitleOffset(1.8);

      hFrame->GetYaxis()->CenterTitle();

      // Format for x axis
      hFrame->GetXaxis()->SetLabelFont(43);
      hFrame->GetXaxis()->SetLabelSize(30);
      hFrame->GetXaxis()->SetLabelOffset(0.02);
      hFrame->GetXaxis()->SetTitleFont(43);
      hFrame->GetXaxis()->SetTitleSize(34);
      hFrame->GetXaxis()->SetTitleOffset(1.5);
      hFrame->GetXaxis()->SetNdivisions(505);

      // Draw cloned histogram with individual settings
      hFrame->Draw();

      hIaaXplot[iCen] = (TH1F *)hPbPb_NLO[iCen]->Clone(Form("hIaaXplot%d_%d", iCen, iCen + 1));
      hIaaXplotSyst[iCen] = (TH1F *)hPbPb_NLOSyst[iCen]->Clone(Form("hIaaXplotSyst%d_%d", iCen, iCen + 1));

      PlotStyle(hIaaXplot[iCen], 20, 1, kGray + 3, kColorMarkFill[iCen], "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
      PlotStyle(hIaaXplotSyst[iCen], 20, 1, kGray + 3, 16, "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", true);

      hIaaXplotSyst[iCen]->SetFillColorAlpha(16, 0.4);
      hIaaXplotSyst[iCen]->Draw("E2Psame ");
      hIaaXplot[iCen]->Draw("EPX0same");
      grIaaNLOmedian[iCen]->SetLineStyle(1);
      grIaaNLOmedian[iCen]->SetLineWidth(3);
      grIaaNLOmedian[iCen]->Draw("l3 same");
      hIaa_nPDF->SetLineStyle(1);
      hIaa_nPDF->SetLineColor(kOrange - 3);
      hIaa_nPDF->SetLineWidth(4);
      hIaa_nPDFSyst->Draw("E3L same");
      hIaa_CNMwithSyst->SetLineWidth(0);
      hIaa_CNMwithSyst->Draw("pl3 same");
      hIaa_nPDF->Draw("HIST L same ");
      lineNLO->Draw("l");
      lineNLO1->Draw("l");
      legRatioNLO->Draw("same");

      if (iCen == 0)
      {
        TLatex *latTitle = new TLatex();
        latTitle->SetTextFont(42);
        latTitle->SetTextSize(0.045);
        latTitle->SetNDC();
        latTitle->DrawLatex(0.24, 0.94, Form("#font[42]{ALICE} Pb#font[122]{-}Pb,#color[0]{.}#sqrt{#it{s}_{NN}} = 5.02 TeV"));
        latTitle->DrawLatex(0.840, 0.94, Form("#bf{0#font[122]{-}30%%}"));
        histIaaCOLBTmedianSyst[0]->SetLineStyle(4);
        histIaaCOLBTmedianSyst[0]->Draw("X0sameE3 ");
        histIaaCOLBTmedian[0]->SetLineStyle(4);
        histIaaCOLBTmedian[0]->SetLineWidth(4);
        histIaaCOLBTmedian[0]->Draw("hist sameL ");

        legPbPb_NLOratioIaaCNMonly[iCen] = LegStd(legPbPb_NLOratioIaaCNMonly[iCen], 0.215, 0.17, 0.77, 0.34);
        legPbPb_NLOratioIaaCNMonly[iCen]->SetTextSize(0.04);
        legPbPb_NLOratioIaaCNMonly[iCen]->AddEntry(hIaa_nPDF, "#it{I}_{AA}, NLO pQCD#color[0]{.}+#color[0]{.}CNM", "l");
        legPbPb_NLOratioIaaCNMonly[iCen]->AddEntry(hIaa_CNMwithSyst, "0.7#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}#it{#mu}#color[0]{.}<#color[0]{..}2#color[0]{.}#it{p}_{T}^{#gamma}", "f");
        legPbPb_NLOratioIaaCNMonly[iCen]->AddEntry(hIaa_nPDFSyst, "nPDF uncert. EPPS21", "f");
        legPbPb_NLOratioIaaCNMonly[iCen]->Draw("same");
      }
      else if (iCen == 1)
      {
        TLatex *latTitle = new TLatex();
        latTitle->SetTextFont(42);
        latTitle->SetTextSize(0.055);
        latTitle->SetNDC();
        latTitle->DrawLatex(0.78, 0.94, Form("#bf{30#font[122]{-}50%%}"));
      }
      else if (iCen == 2)
      {
        TLatex *latTitle = new TLatex();
        latTitle->SetTextFont(42);
        latTitle->SetTextSize(0.05);
        latTitle->SetNDC();
        latTitle->DrawLatex(0.72, 0.94, Form("#bf{50#font[122]{-}90%%}"));
        // latTitle->DrawLatex(0.14, 0.26, Form("#font[42]{ALICE} Pb#font[122]{-}Pb,#color[0]{..}#sqrt{#it{s}_{NN}} = 5.02 TeV"));
        latTitle->DrawLatex(0.04, 0.26 - 0 * 0.055, Form("|#Delta#it{#varphi}| > #frac{3}{5}#color[0]{.}#it{#pi}, |#it{#eta}^{#it{#gamma}}| < 0.67 "));
        latTitle->DrawLatex(0.04, 0.26 - 1 * 0.055, Form("%2.0f <#color[0]{..}#it{p}_{T}^{#it{#gamma}} < %2.0f GeV/#it{c}#color[0]{.};#color[0]{..}#it{p}_{T}^{h} > 0.5 GeV/#it{c} ", ptMin, ptMax));
      }
      if (iCen == 1)
      {
        //legRatioPbPb[iCen] = LegStd(legRatioPbPb[iCen], 0.01, 0.225, 0.50, 0.27);
        //legRatioPbPb[iCen]->SetTextSize(0.05);
        //legRatioPbPb[iCen]->AddEntry(hIaaXplotSyst[iCen], " Data ", "lfpe");
        // legRatioPbPb[iCen]->AddEntry(hIaa_nPDF, "#it{I}_{AA}, NLO pQCD + CNM", "lf");
        // legRatioPbPb[iCen]->AddEntry((TObject *)0, "M. Xie et al.", "");
        legPbPb_NLOratio[iCen] = LegStd(legPbPb_NLOratio[iCen], 0.01, 0.17, 0.56, 0.34);
        legPbPb_NLOratio[iCen]->SetTextSize(0.05);
        legPbPb_NLOratio[iCen]->AddEntry(hIaaXplotSyst[iCen], "Data ", "lfpe");
        //legPbPb_NLOratio[iCen]->AddEntry(hIaa_CNMwithSyst, "NLO pQCD#color[0]{.}+#color[0]{.}CNM, 0.7#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}#it{#mu}#color[0]{.}<#color[0]{..}2#color[0]{.}#it{p}_{T}^{#gamma}", "lf");
        //legPbPb_NLOratio[iCen]->AddEntry(hIaa_nPDFSyst, "#it{I}_{AA}, NLO pQCD#color[0]{.}+#color[0]{.}CNM", "lf");
        legPbPb_NLOratio[iCen]->AddEntry(grIaaNLOmedian[iCen], "#it{I}_{AA}, NLO pQCD +#color[0]{.}#Delta#it{E}_{loss}", "lf");

        // histIaaCOLBTmedianSyst[1]->SetFillColorAlpha(kAzure-3, 0.20);
        histIaaCOLBTmedianSyst[1]->SetLineStyle(4);
        // histIaaCOLBTmedianSyst[1]->SetLineColor(kAzure-1);
        histIaaCOLBTmedianSyst[1]->Draw("X0sameE3 ");
        // histIaaCOLBTmedian[1]->SetLineColor(kAzure-1);
        histIaaCOLBTmedian[1]->SetLineStyle(4);
        histIaaCOLBTmedian[1]->SetLineWidth(3);
        histIaaCOLBTmedian[1]->Draw("histsameL ");
        legPbPb_NLOratio[iCen]->AddEntry(histIaaCOLBTmedianSyst[iCen], "#it{I}_{AA}, CoLBT-hydro", "lf");
        // legPbPb_NLOratio[iCen]->AddEntry((TObject *)0, "X. N. Wang et al.", "");
        // legRatioPbPb[iCen]->Draw("same");
        legPbPb_NLOratio[iCen]->Draw("same");
      }
    }
  }
  cPbPb_NLORatioTogetherTest->cd();
  cPbPb_NLORatioTogetherTest->Print(dirPlot + Form("/Ipqcd_AllCen.pdf"));

  TCanvas *cPbPb_NLORatioAll = canvasStdIaa(Form("cPbPb_NLORatioAll"), 1, 1);
  TLatex *latPbPb_NLOAll;
  TLegend *legPbPb_NLOratioAll = LegStd(legPbPb_NLOratioAll, 0.55, 0.66, 0.82, 0.82);
  TH1F *hPbPb_NLO_Shifted[nCen];
  TH1F *hPbPb_NLOSyst_Shifted[nCen];
  // legPbPb_NLOratioAll->SetNColumns(3);
  legPbPb_NLOratioAll->SetTextSize(0.036);
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    // hGeneralRatio->GetYaxis()->SetNdivisions(505);
    // hGeneralRatio->GetYaxis()->SetNdivisions(505);
    //  hGeneralRatio->GetXaxis()->SetRangeUser(0.15, 1.0);
    hPbPb_NLO_Shifted[iCen] = new TH1F(Form("hPbPb_NLO_Shifted_Cen%d_%d", iCen, iCen + 1), Form("hPbPb_NLO_Shifted_Cen%d_%d", iCen, iCen + 1), 105, 0, 1.05);
    hPbPb_NLOSyst_Shifted[iCen] = new TH1F(Form("hPbPb_NLOSyst_Shifted_Cen%d_%d", iCen, iCen + 1), Form("hPbPb_NLOSyst_Shifted_Cen%d_%d", iCen, iCen + 1), 105, 0, 1.05);
    PlotStyle(hPbPb_NLO_Shifted[iCen], kMarkCen[iCen], 1, kColorMark[iCen], kColorMarkFill[iCen], "#it{z}_{T}", "Ratio", false);
    PlotStyle(hPbPb_NLOSyst_Shifted[iCen], kMarkCen[iCen], 1, kColorMark[iCen], kColorMarkFill[iCen], "#it{z}_{T}", "Ratio", true);
    hGeneralRatio->SetTitle(" ");
    hGeneralRatio->GetYaxis()->CenterTitle(false);
    hGeneralRatio->GetYaxis()->SetRangeUser(0, 2.5);
    hGeneralRatio->GetYaxis()->SetTitleSize(0.046);
    hGeneralRatio->GetYaxis()->SetTitleOffset(1.);
    hGeneralRatio->GetYaxis()->SetLabelSize(0.040);
    hGeneralRatio->GetXaxis()->SetRangeUser(0.05, 1.05);
    hGeneralRatio->GetXaxis()->SetNdivisions(505);
    hGeneralRatio->GetXaxis()->SetLabelSize(0.040);
    hGeneralRatio->GetXaxis()->SetTitleSize(0.055);
    hGeneralRatio->SetLineWidth(0);
    hGeneralRatio->Draw("histsame");
    ///////// Shift a bit the distributions //////
    // hPbPb_NLO_Shifted[iCen] = shiftBinXinHisto(hPbPb_NLO[iCen], hPbPb_NLO_Shifted[iCen], iCen);
    // hPbPb_NLOSyst_Shifted[iCen] = shiftBinXinHisto(hPbPb_NLOSyst[iCen], hPbPb_NLOSyst_Shifted[iCen], iCen);
    ////
    // hPbPb_NLOSyst_Shifted[iCen]->Draw("E2Psame");
    // hPbPb_NLO_Shifted[iCen]->Draw("EPX0same");
    // hPbPb_NLOSyst[iCen]->GetYaxis()->SetTitle("#it{I}_{pQCD}");
    hPbPb_NLOSyst[iCen]->Draw("E2Psame");
    hPbPb_NLO[iCen]->Draw("EPX0same");
    //  grIaaNLOmedian[iCen]->Draw("pl3 same");
    // hPbPb_NLOSyst[iCen]->SetLineColor(kWhite);
    legPbPb_NLOratioAll->AddEntry(hPbPb_NLOSyst[iCen], Form("%d#font[122]{-}%d%%", cenBins[iCen], cenBins[iCen + 1]), "epf");
    // legPbPb_NLOratioAll->AddEntry(hPbPb_NLOSyst[iCen], "syst. unc.", "f");
  }
  // latPbPb_NLOAll = LatexStdISORatio(latPbPb_NLOAll, 0.40, 0.92, 0.04, cenBins[0], cenBins[1], ptMin, ptMax, false);
  legPbPb_NLOratioAll->Draw("same");
  legTitleRatio[0]->Draw("same");
  legTitleRatioSelection[0]->Draw("same");
  // legRatioNLO->Draw("same");
  lineNLO->Draw("l");
  lineNLO1->Draw("l");
  cPbPb_NLORatioAll->Print(dirPlot + Form("/I_NLOCenAll.pdf"));

  //////////////////////////////////////////////////
  ////////// Comparison with I_pPbpp  //////////////
  //////////////////////////////////////////////////

  TCanvas *cPbPb_NLORatioAll_IpPbpp = canvasStdIaa(Form("cPbPb_NLORatioAll_IpPbpp"), 1, 1);
  TLatex *latPbPb_NLOAll_IpPbpp;
  TLegend *legPbPb_NLOratioAll_IpPbpp = LegStd(legPbPb_NLOratioAll_IpPbpp, 0.52, 0.48, 0.96, 0.78);
  double assocZtThinnerPbpp[] = {0., 0.06, 0.08, 0.10, 0.107, 0.142, 0.15, 0.19, 0.20, 0.253, 0.30, 0.337, 0.40, 0.45, 0.60, 0.80, 1.00, 1.05};
  TH1F *hGeneralRatiopPbpp = new TH1F("hGeneralRatiopPbpp", "hGeneralRatiopPbpp", 17, assocZtThinnerPbpp);
  PlotStyle(hGeneralRatiopPbpp, 20, 1, kWhite, kWhite, " #it{z}_{T} ", " #it{I}_{pQCD}, #it{I}_{pA}", false);
  legPbPb_NLOratioAll_IpPbpp->SetNColumns(2);
  legPbPb_NLOratioAll_IpPbpp->SetTextSize(0.034);
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    // hGeneralRatiopPbpp->GetYaxis()->SetNdivisions(505);
    // hGeneralRatio->GetYaxis()->SetNdivisions(505);
    //  hGeneralRatio->GetXaxis()->SetRangeUser(0.15, 1.0);
    hGeneralRatiopPbpp->SetTitle(" ");
    hGeneralRatiopPbpp->GetYaxis()->CenterTitle(false);
    hGeneralRatiopPbpp->GetYaxis()->SetRangeUser(-0.01, 2.9);
    hGeneralRatiopPbpp->GetYaxis()->SetTitleSize(0.04);
    hGeneralRatiopPbpp->GetYaxis()->SetTitleOffset(1.);
    hGeneralRatiopPbpp->GetYaxis()->SetLabelSize(0.040);
    hGeneralRatiopPbpp->GetXaxis()->SetRangeUser(0.05, 1.05);
    hGeneralRatiopPbpp->GetXaxis()->SetLabelSize(0.040);
    hGeneralRatiopPbpp->GetXaxis()->SetTitleSize(0.055);
    hGeneralRatiopPbpp->SetLineWidth(0);
    hGeneralRatiopPbpp->Draw("histsame");
    hPbPb_NLOSyst[iCen]->Draw("E2Psame ");
    hPbPb_NLO[iCen]->Draw("EPX0same");
    // hPbPb_NLOSyst_Shifted[iCen]->Draw("E2Psame ");
    // hPbPb_NLO_Shifted[iCen]->Draw("EPX0same");
    //  grIaaNLOmedian[iCen]->Draw("pl3 same");
    //  hPbPb_NLOSyst[iCen]->SetLineColor(kWhite);
    legPbPb_NLOratioAll_IpPbpp->AddEntry(hPbPb_NLOSyst[iCen], Form("#it{I}_{pQCD}, %d#font[122]{-}%d%%", cenBins[iCen], cenBins[iCen + 1]), "epf");
    // legPbPb_NLOratioAll_IpPbpp->AddEntry(hPbPb_NLOSyst[iCen], "syst. unc.", "f");
  }

  hIpPbpp_syst->SetFillStyle(0);
  hIpPbpp_syst->SetLineColor(kBlack);
  hIpPbpp_syst->SetLineWidth(1);
  hIpPbpp_syst->Draw("E2same");
  hIpPbpp->Draw("EPX0same");
  legPbPb_NLOratioAll_IpPbpp->AddEntry(hIpPbpp_syst, "#it{I}_{pA}, p#font[122]{-}Pb/pp", "epf");
  // legPbPb_NLOratioAll_IpPbpp->AddEntry(hIpPbpp_syst, "syst. unc.", "epf");
  legPbPb_NLOratioAll_IpPbpp->AddEntry((TObject *)0, "12 <#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}< 40 GeV/#it{c}", "");
  latPbPb_NLOAll_IpPbpp = LatexStdISORatioNoPbPb(latPbPb_NLOAll_IpPbpp, 0.440, 0.92, 0.04, cenBins[0], cenBins[1], ptMin, ptMax, false);
  // latPbPb_NLOAll_IpPbpp = new TLatex();
  // latPbPb_NLOAll_IpPbpp->SetTextFont(42);
  // latPbPb_NLOAll_IpPbpp->SetTextSize(texSize);
  // latPbPb_NLOAll_IpPbpp->SetNDC();
  // latPbPb_NLOAll_IpPbpp->DrawLatex(xpos, ypos, Form("#font[42]{ALICE}"));
  //   latPbPb_NLOAll_IpPbpp->DrawLatex(xpos, ypos - 0.06, Form("#bf{%d#font[122]{-}%d%%} Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV ", cenMin, cenMax));
  // latPbPb_NLOAll_IpPbpp->DrawLatex(xpos, ypos - 2 * 0.06, Form("|#Delta#it{#varphi}| > #frac{3}{5} #it{#pi}, |#it{#eta}^{#it{#gamma}}| < 0.67 "));
  // latPbPb_NLOAll_IpPbpp->DrawLatex(xpos, ypos - 3 * 0.06, Form("%2.0f < #it{p}_{T}^{#it{#gamma}} < %2.0f GeV/#it{c} ; #it{p}_{T}^{h} > 0.5 GeV/#it{c} ", ptMin, ptMax));

  legPbPb_NLOratioAll_IpPbpp->Draw("same");
  legRatioNLO->Draw("same");
  lineNLO->Draw("l");
  lineNLO1->Draw("l");
  cPbPb_NLORatioAll_IpPbpp->Print(dirPlot + Form("/I_NLOCenAll_IpPbpp.pdf"));

  /////////////////////////////////////////////////////////
  ////////// Comparison with I_aa nPDF only  //////////////
  ////////////////////////////////////////////////////////

  TCanvas *cPbPb_NLORatioAll_IaanPDF = canvasStdIaa(Form("cPbPb_NLORatioAll_IaanPDF"), 1, 1);
  TLatex *latPbPb_NLOAll_IaanPDF;
  TLegend *legPbPb_NLOratioAll_IaanPDF = LegStd(legPbPb_NLOratioAll_IaanPDF, 0.54, 0.48, 0.96, 0.740);
  legPbPb_NLOratioAll_IaanPDF->SetNColumns(2);
  legPbPb_NLOratioAll_IaanPDF->SetTextSize(0.034);
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    hGeneralRatio->SetTitle(" ");
    hGeneralRatio->GetYaxis()->CenterTitle(false);
    hGeneralRatio->GetYaxis()->SetRangeUser(-0.01, 2.8);
    hGeneralRatio->GetYaxis()->SetTitleSize(0.046);
    hGeneralRatio->GetYaxis()->SetTitleOffset(1.);
    hGeneralRatio->GetYaxis()->SetLabelSize(0.040);
    hGeneralRatio->GetXaxis()->SetRangeUser(0.05, 1.05);
    hGeneralRatio->GetXaxis()->SetLabelSize(0.040);
    hGeneralRatio->GetXaxis()->SetTitleSize(0.055);
    hGeneralRatio->SetLineWidth(0);
    hGeneralRatio->Draw("histsame");
    hPbPb_NLOSyst[iCen]->Draw("E2Psame ");
    hPbPb_NLO[iCen]->Draw("EPX0same");
    // grIaaNLOmedian[iCen]->Draw("pl3 same");
    // hPbPb_NLOSyst[iCen]->SetLineColor(kWhite);
    legPbPb_NLOratioAll_IaanPDF->AddEntry(hPbPb_NLO[iCen], Form("#it{I}_{pQCD}, %d#font[122]{-}%d%% stat.", cenBins[iCen], cenBins[iCen + 1]), "ep");
    legPbPb_NLOratioAll_IaanPDF->AddEntry(hPbPb_NLOSyst[iCen], "syst. unc.", "f");
  }
  hIaa_nPDF->Draw("hist same c");
  legPbPb_NLOratioAll_IaanPDF->AddEntry(hIaa_nPDF, "#it{I}_{AA}, pQCD + CNM", "l");
  latPbPb_NLOAll_IaanPDF = LatexStdISORatio(latPbPb_NLOAll_IaanPDF, 0.440, 0.92, 0.04, cenBins[0], cenBins[1], ptMin, ptMax, false);
  legPbPb_NLOratioAll_IaanPDF->Draw("same");
  legRatioNLO->Draw("same");
  lineNLO->Draw("l");
  lineNLO1->Draw("l");
  cPbPb_NLORatioAll_IaanPDF->Print(dirPlot + Form("/I_NLOCenAll_IaanPDF.pdf"));

  ///////////////////////////////////////////////////////////////////////////
  ////////// Comparison with I_aa central, I_pPbpp and I_CNM ///////////////
  /////////////////////////////////////////////////////////////////////////

  TCanvas *cPbPb_NLORatio0_10_Iaa_CNM_pA = canvasStdIaa(Form("cPbPb_NLORatio0_10_Iaa_CNM_pA"), 1, 1);
  TLegend *legPbPb_NLORatio0_10_Iaa_CNM_pA = LegStd(legPbPb_NLOratioAll_IaanPDF, 0.52, 0.52, 0.92, 0.740);
  legPbPb_NLORatio0_10_Iaa_CNM_pA->SetNColumns(2);
  legPbPb_NLORatio0_10_Iaa_CNM_pA->SetTextSize(0.034);

  hGeneralRatio->SetTitle(" ");
  hGeneralRatio->GetYaxis()->CenterTitle(false);
  hGeneralRatio->GetYaxis()->SetRangeUser(-0.001, 3.5);
  hGeneralRatio->GetYaxis()->SetTitleSize(0.046);
  hGeneralRatio->GetYaxis()->SetTitleOffset(1.);
  hGeneralRatio->GetYaxis()->SetLabelSize(0.040);
  hGeneralRatio->GetXaxis()->SetRangeUser(-0.01, 1.05);
  hGeneralRatio->GetXaxis()->SetLabelSize(0.040);
  hGeneralRatio->GetXaxis()->SetTitleSize(0.055);
  hGeneralRatio->SetLineWidth(0);
  hGeneralRatio->GetYaxis()->SetTitle("#it{I}_{pQCD}, #it{I}_{pA}");
  hGeneralRatio->Draw("histsame");
  hPbPb_NLOSyst[0]->Draw("E2Psame ");
  hPbPb_NLO[0]->Draw("EPX0same");
  hIpPbpp_syst->SetFillStyle(0);
  hIpPbpp_syst->SetLineColor(kBlack);
  hIpPbpp_syst->SetLineWidth(2);
  hIpPbpp_syst->Draw("E2same");
  hIpPbpp->Draw("EPX0same");
  hIaa_nPDFSyst->Draw("pl3 same");
  hIaa_CNMwithSyst->Draw("pl3 same");
  hIaa_nPDF->SetLineWidth(4);
  hIaa_nPDF->Draw("hist l same");
  //hIaa_nPDF->Draw("hist same c ");
  // grIaaNLOmedian[iCen]->Draw("pl3 same");
  // hPbPb_NLOSyst[iCen]->SetLineColor(kWhite);
  // legPbPb_NLORatio0_10_Iaa_CNM_pA->SetTextSize(0.038);
  // legPbPb_NLORatio0_10_Iaa_CNM_pA->AddEntry(hPbPb_NLO[0], Form("#it{I}_{pQCD}, %d#font[122]{-}%d%% stat.", cenBins[0], cenBins[1]), "ep");
  // legPbPb_NLORatio0_10_Iaa_CNM_pA->AddEntry(hPbPb_NLOSyst[0], "syst. unc.", "f");
  legPbPb_NLORatio0_10_Iaa_CNM_pA->AddEntry(hIpPbpp, "#it{I}_{pA}, p#font[122]{-}Pb/pp stat.", "ep");
  legPbPb_NLORatio0_10_Iaa_CNM_pA->AddEntry(hIpPbpp_syst, "syst. unc.", "f");

  // legPbPb_NLORatio0_10_Iaa_CNM_pA->Draw("same");
  legRatioNLO->Draw("same");
  lineNLO->Draw("l");
  lineNLO1->Draw("l");

  TLegend *legIaaALICE030_gamhad = LegStd(legIaaALICE030_gamhad, 0.12, 0.92, 0.34, 0.96);
  legIaaALICE030_gamhad->SetTextSize(0.036);
  legIaaALICE030_gamhad->SetHeader("ALICE#color[0]{..}#sqrt{#it{s}_{NN}} = 5.02 TeV");
  // legIaaALICE030_gamhad->AddEntry((TObject *)0, "18 <#color[0]{..}#it{p}_{T}^{#it{#gamma}} < 40 GeV/#it{c}#color[0]{.};#color[0]{..}#it{p}_{T}^{h} > 0.5 GeV/#it{c}", "");
  legIaaALICE030_gamhad->Draw("same");

  TLegend *legIaaALICE030pp = LegStd(legIaaALICE030pp, 0.45, 0.86, 0.55, 0.97);
  legIaaALICE030pp->SetTextSize(0.036);
  legIaaALICE030pp->SetHeader(" #color[0]{.}#it{p}_{T}^{iso, ch} < 1.5 GeV/#it{c}, #bf{#it{R} = 0.2}");
  legIaaALICE030pp->AddEntry((TObject *)0, "#bf{18 <#color[0]{..}#it{p}_{T}^{#it{#gamma}} < 40 GeV/#it{c}}#color[0]{.};#color[0]{..}#it{p}_{T}^{h} > 0.5 GeV/#it{c}", "");
  legIaaALICE030pp->Draw("same");

  TLegend *legIaaALICE030Marker = LegStd(legIaaALICE030Marker, 0.45, 0.75, 0.70, 0.805);
  legIaaALICE030Marker->SetTextSize(0.032);
  legIaaALICE030Marker->AddEntry(hPbPb_NLOSyst[0], "0#font[122]{-}30% Pb#font[122]{-}Pb", "epf");
  legIaaALICE030Marker->Draw("same");
  TLegend *legIaaALICE030MarkerIaa_nPDF = LegStd(legIaaALICE030MarkerIaa_nPDF, 0.66, 0.70, 0.92, 0.86);
  legIaaALICE030MarkerIaa_nPDF->SetTextSize(0.032);
  //legIaaALICE030Marker->AddEntry(hIaa_nPDFSyst, "NLO pQCD + CNM", "lf");
  legIaaALICE030MarkerIaa_nPDF->AddEntry(hIaa_nPDF, "#it{I}_{AA}, NLO pQCD + CNM", "l");
  legIaaALICE030MarkerIaa_nPDF->AddEntry(hIaa_CNMwithSyst, "0.7#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}#it{#mu}#color[0]{.}<#color[0]{..}2#it{p}_{T}^{#gamma} ", "f");
  legIaaALICE030MarkerIaa_nPDF->AddEntry(hIaa_nPDFSyst, "nPDF uncert. EPPS21", "f");
  legIaaALICE030MarkerIaa_nPDF->Draw("same");

  TLegend *legIaaALICEppPb = LegStd(legIaaALICEppPb, 0.45, 0.59, 0.55, 0.70);
  legIaaALICEppPb->SetTextSize(0.036);
  // legIaaALICEppPb->SetHeader("Phys Rev C 102 (2020) 044908");
  //  legIaaALICEppPb->AddEntry((TObject *)0, "pp & p#font[122]{-}Pb", "");
  legIaaALICEppPb->AddEntry((TObject *)0, "#color[0]{.}#it{p}_{T}^{iso, ch} < 1.5 GeV/#it{c}, #bf{#it{R} = 0.4}", "");
  legIaaALICEppPb->AddEntry((TObject *)0, "#bf{12 <#color[0]{..}#it{p}_{T}^{#it{#gamma}} < 40 GeV/#it{c}}#color[0]{.}; 0.5 <#color[0]{..}#it{p}_{T}^{h} < 10 GeV/#it{c}", "");
  legIaaALICEppPb->Draw("same");
  TLegend *legIaaALICEppPbMarker = LegStd(legIaaALICEppPbMarker, 0.45, 0.51, 0.70, 0.565);
  // legIaaALICEppPbMarker->SetNColumns(2);
  legIaaALICEppPbMarker->SetTextSize(0.032);
  // legIaaALICEppPbMarker->AddEntry(hIpPbpp, "#it{I}_{pA}, p#font[122]{-}Pb/pp stat.", "ep");
  legIaaALICEppPbMarker->AddEntry(hIpPbpp_syst, "p#font[122]{-}Pb/pp", "epf");
  legIaaALICEppPbMarker->Draw("same");
  TLegend *legIaaALICEppPbPaper = LegStd(legIaaALICEppPbPaper, 0.45, 0.45, 0.55, 0.51);
  legIaaALICEppPbPaper->SetTextSize(0.032);
  legIaaALICEppPbPaper->AddEntry("", "Phys. Rev. C 102 (2020) 044908", "");
  legIaaALICEppPbPaper->Draw("same");

  cPbPb_NLORatio0_10_Iaa_CNM_pA->Print(dirPlot + Form("/Ipqcd_Cen0_30_Iaa_CNM_pA.pdf"));

  TH1F *hGeneral = new TH1F("hGeneral", "hGeneral", nZtBinThin, assocZtThinner);
  PlotStyle(hGeneral, 20, 1, kWhite, kWhite, "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
  hSystZt[0]->GetXaxis()->SetRangeUser(0.10, 1.0);
  hZtCent[0]->GetXaxis()->SetRangeUser(0.10, 1.0);
  hSystZt[1]->GetXaxis()->SetRangeUser(0.10, 0.6);
  hZtCent[1]->GetXaxis()->SetRangeUser(0.10, 0.6);
  hSystZt[2]->GetXaxis()->SetRangeUser(0.10, 0.6);
  hZtCent[2]->GetXaxis()->SetRangeUser(0.10, 0.6);
  hZt_MC_Gen[0]->GetXaxis()->SetRangeUser(0.10, 1.0);
  hZt_MC_Gen[1]->GetXaxis()->SetRangeUser(0.10, 0.6);
  hZt_MC_Gen[2]->GetXaxis()->SetRangeUser(0.10, 0.6);
  h3[0]->GetXaxis()->SetRangeUser(0.10, 1.0);
  h3[1]->GetXaxis()->SetRangeUser(0.10, 0.6);
  h3[2]->GetXaxis()->SetRangeUser(0.10, 0.6);
  hsyst[0]->GetXaxis()->SetRangeUser(0.10, 1.0);
  hsyst[1]->GetXaxis()->SetRangeUser(0.10, 0.6);
  hsyst[2]->GetXaxis()->SetRangeUser(0.10, 0.6);

  TLegend *legZTData = LegStd(legZTData, 0.7, 0.50, 0.85, 0.70);
  TCanvas *cDiffCent = new TCanvas("cDiffCent", "cDiffCent", 800, 600);
  TH1F *hZtDiffCent[nCen];
  TH1F *hZtDiffCentSys[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    hZtDiffCent[iCen] = (TH1F *)hZtCent[iCen]->Clone(Form("hZt%d_%d %%", cenBins[iCen], cenBins[iCen + 1]));
    hZtDiffCentSys[iCen] = (TH1F *)hSystZt[iCen]->Clone(Form("hZtSystm%d_%d %%", cenBins[iCen], cenBins[iCen + 1]));
    gPad->SetLogy();
    PlotStyle(hZtDiffCent[iCen], kMarkCen[iCen], 1, kColorMark[iCen], kColorMarkFill[iCen], "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
    hZtDiffCent[iCen]->SetTitle(" ");
    hZtDiffCentSys[iCen]->SetTitle(" ");
    hZtDiffCent[iCen]->GetYaxis()->SetRangeUser(1e-3, 90);
    hZtDiffCentSys[iCen]->SetMarkerStyle(kMarkCen[iCen]);
    hZtDiffCentSys[iCen]->SetMarkerColor(kColorMark[iCen]);
    hZtDiffCentSys[iCen]->SetFillStyle(0);
    hGeneral->GetYaxis()->SetRangeUser(2 * 1e-3, 99);
    hGeneral->Draw("hist");
    hZtDiffCentSys[iCen]->Draw("E2Psame");
    hZtDiffCent[iCen]->Draw("EPX0same");
    legZTData->AddEntry(hZtDiffCent[iCen], Form("%d-%d %%", cenBins[iCen], cenBins[iCen + 1]), "ep");
  }

  TLatex *lat = LatexStdISORatio(lat, 0.450, 0.84, 0.045, cenBins[0], cenBins[0], ptMin, ptMax, false);
  legZTData->Draw();
  cDiffCent->Print(dirPlot + Form("/ztDiffMethMixed_Cen%d_%d%s.pdf", cenBins[2], cenBins[3], sPtAll.Data()));

  TLegend *legZTPbPbpp[nCen];
  TCanvas *cPbPb[nCen];
  TLatex *lat0 = new TLatex();
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    legZTPbPbpp[iCen] = LegStd(legZTPbPbpp[iCen], 0.180, 0.18, 0.400, 0.30);
    cPbPb[iCen] = new TCanvas(Form("cPbPbCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("cPbPbCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), 800, 600);
    gPad->SetLogy();
    legZTPbPbpp[iCen]->AddEntry(hZt_MC_Gen[iCen], Form("PYTHIA 8 pp "), "l");
    legZTPbPbpp[iCen]->AddEntry(hZtCent[iCen], "Pb#font[122]{-}Pb stat. unc. ", "ep");
    hGeneral->GetYaxis()->SetRangeUser(2 * 1e-3, 99);
    hGeneral->Draw("histSAME");
    hSystZt[iCen]->Draw("E2Psame ");
    hZtCent[iCen]->Draw("EPX0same");
    hZt_MC_Gen[iCen]->SetFillStyle(0);
    hZt_MC_Gen[iCen]->SetLineWidth(10);
    hZt_MC_Gen[iCen]->SetLineColorAlpha(kBlack, 1);
    hZt_MC_Gen[iCen]->Draw("sameE3");
    legZTPbPbpp[iCen]->Draw("same");
    lat0 = LatexStdISO(lat0, 0.45, 0.84, 0.05, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true);
    cPbPb[iCen]->Print(dirPlot + Form("/ztPbPbOverlapppCen%d_%d%s.pdf", cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
  }

  TCanvas *cDiffCent_Pythia[nCen];
  TLegend *legdiffcen[nCen];
  TLegend *legdiffcenpp[nCen];
  TLatex *latCent[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    cDiffCent_Pythia[iCen] = new TCanvas(Form("cDiffCent_Pythia%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("cDiffCent_Pythia%d_%d", cenBins[iCen], cenBins[iCen + 1]), 800, 600);
    legdiffcen[iCen] = LegStd(legdiffcen[iCen], 0.60, 0.615, 0.80, 0.705);
    legdiffcenpp[iCen] = LegStd(legdiffcenpp[iCen], 0.170, 0.12, 0.350, 0.37);
    cDiffCent_Pythia[iCen]->SetTopMargin(0.015);
    cDiffCent_Pythia[iCen]->SetRightMargin(0.02);
    cDiffCent_Pythia[iCen]->SetLeftMargin(0.15);
    cDiffCent_Pythia[iCen]->SetBottomMargin(0.11);
    // hSystZt[iCen]->GetXaxis()->SetLabelSize(0.028);
    // hSystZt[iCen]->GetXaxis()->SetTitleSize(0.032);
    // hSystZt[iCen]->GetXaxis()->SetTitle("#font[12]{{z}_{T}}");
    gPad->SetLogy();
    hGeneral->SetTitle(" ");
    hGeneral->GetXaxis()->SetLabelSize(0.04);
    hGeneral->GetXaxis()->SetLabelOffset(0.01);
    hGeneral->GetXaxis()->SetTitleOffset(1);
    hGeneral->GetYaxis()->SetLabelSize(0.04);
    hGeneral->GetXaxis()->SetTitleSize(0.045);
    hGeneral->GetYaxis()->SetTitleSize(0.045);
    hGeneral->GetYaxis()->SetTitleOffset(1.4);
    // hSystZt[iCen]->SetFillStyle(0);
    hGeneral->Draw("histSAME");
    hSystZt[iCen]->Draw("E2P same ");
    hZtCent[iCen]->Draw("EPX0same");
    hSystZt[iCen]->SetLineColor(kWhite);
    legdiffcen[iCen]->AddEntry(hZtCent[iCen], "Pb#font[122]{-}Pb stat. unc. ", "ep");
    legdiffcen[iCen]->AddEntry(hSystZt[iCen], "Pb#font[122]{-}Pb syst. unc. ", "f");
    legdiffcen[iCen]->Draw("same");
    legdiffcenpp[iCen]->Draw("same");
    TLatex *latdiffCent1 = LatexStdISO(latdiffCent1, 0.30, 0.94, 0.04, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true);
    // cDiffCent_Pythia[iCen]->Print(dirPlot + Form("/withModelNoppztPbPb_pp_Cen%d_%d%s.pdf", cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
    //cDiffCent_Pythia[iCen]->Print(dirPlot + Form("/ztPbPb_Cen%d_%d%s.pdf", cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
  }

  TCanvas *cDiffCent_Theory[nCen];
  TLegend *legdiffcen_Theory[nCen];
  TLatex *latCent_Theory[nCen];
  TLegend *legZTPbPbppSingleMod_Theory[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    // = canvasStdIaa("cAllZtNLOpQCD", 1, 1);
    cDiffCent_Theory[iCen] = canvasStdIaa(Form("cDiffCent_Theory%d_%d", cenBins[iCen], cenBins[iCen + 1]), 1, 1);
    legdiffcen_Theory[iCen] = LegStd(legdiffcen_Theory[iCen], 0.550, 0.66, 0.900, 0.70);
    cDiffCent_Theory[iCen]->SetTopMargin(0.015);
    cDiffCent_Theory[iCen]->SetRightMargin(0.02);
    cDiffCent_Theory[iCen]->SetLeftMargin(0.15);
    cDiffCent_Theory[iCen]->SetBottomMargin(0.11);
    // hSystZt[iCen]->GetXaxis()->SetLabelSize(0.028);
    // hSystZt[iCen]->GetXaxis()->SetTitleSize(0.032);
    // hSystZt[iCen]->GetXaxis()->SetTitle("#font[12]{{z}_{T}}");
    gPad->SetLogy();
    hGeneral->SetTitle(" ");
    // hGeneral->GetYaxis()->SetRangeUser(2*1e-4, 99);
    // hGeneral->GetXaxis()->SetLabelSize(0.04);
    // hGeneral->GetXaxis()->SetLabelOffset(0.01);
    // hGeneral->GetXaxis()->SetTitleOffset(1);
    // hGeneral->GetYaxis()->SetLabelSize(0.04);
    // hGeneral->GetXaxis()->SetTitleSize(0.045);
    // hGeneral->GetYaxis()->SetTitleSize(0.045);
    // hGeneral->GetYaxis()->SetTitleOffset(1.4);
    hGeneral->GetYaxis()->SetTitleSize(0.055);
    // hGeneral->GetXaxis()->SetRangeUser(-0.05, 1.1);
    // hGeneral->GetXaxis()->SetNdivisions(505);
    hGeneral->GetXaxis()->SetLabelSize(0.040);
    hGeneral->GetXaxis()->SetTitleSize(0.055);
    hGeneral->GetYaxis()->SetRangeUser(1e-3, 40);
    hGeneral->GetYaxis()->SetLabelSize(0.04);
    hGeneral->GetYaxis()->SetTitleSize(0.045);
    hGeneral->GetYaxis()->SetTitleOffset(1.0);
    // hSystZt[iCen]->SetFillStyle(0);
    hGeneral->Draw("histSAME");
    hSystZt[iCen]->Draw("E2P same ");
    hZtCent[iCen]->Draw("EPX0same");
    hSystZt[iCen]->SetLineColor(kWhite);
    // legdiffcen_Theory[iCen]->AddEntry(hZtCent[iCen], "Pb#font[122]{-}Pb stat. unc. ", "ep");
    legdiffcen_Theory[iCen]->AddEntry(hSystZt[iCen], "Data", "epf");
    // legdiffcen_Theory[iCen]->Draw("same");
    grDztNLOmedian[iCen]->Draw("pl3 same");
    grDztNLOmedianpp->SetLineStyle(1);
    grDztNLOmedianpp->SetLineWidth(3);
    grDztNLOmedianppSyst->SetLineStyle(1);
    grDztNLOmedianppSyst->SetLineWidth(3);
    grDztNLOmedianppSyst->Draw("pl3 same");
    grDztNLOmedianpp->Draw("pl same ");

    if (iCen == 0)
    {
      histDztPbPbCOLBTmedianSyst[0]->SetLineStyle(4);
      histDztPbPbCOLBTmedian[0]->SetLineStyle(4);
      histDztPbPbCOLBTmedianSyst[0]->Draw("X0sameE3 ");
      histDztPbPbCOLBTmedian[0]->SetFillStyle(0);
      histDztPbPbCOLBTmedian[0]->Draw(" HISTSAMEL ");
    }
    else if (iCen == 1)
    {
      histDztPbPbCOLBTmedianSyst[1]->SetLineStyle(4);
      histDztPbPbCOLBTmedian[1]->SetLineStyle(4);
      histDztPbPbCOLBTmedianSyst[1]->Draw("X0sameE3 ");
      histDztPbPbCOLBTmedian[1]->SetFillStyle(0);
      histDztPbPbCOLBTmedian[1]->Draw(" HISTSAMEL ");
    }

    if (iCen == 0)
    {
      legZTPbPbppSingleMod_Theory[iCen] = LegStd(legZTPbPbppSingleMod_Theory[iCen], 0.12, 0.16, 0.350, 0.41);

      legZTPbPbppSingleMod_Theory[iCen]->AddEntry(hSystZt[iCen], "Data", "epf");
      legZTPbPbppSingleMod_Theory[iCen]->AddEntry(grDztNLOmedianppSyst, "NLO pQCD, pp", "lf");
      legZTPbPbppSingleMod_Theory[iCen]->AddEntry(grDztNLOmedian[0], "NLO pQCD#color[0]{.}+#color[0]{.}#Delta#it{E}_{loss}", "lf");
      // legZTPbPbppSingleMod_Theory[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      //  legZTPbPbppSingleMod_Theory[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
      legZTPbPbppSingleMod_Theory[iCen]->AddEntry(histDztPbPbCOLBTmedianSyst[0], "CoLBT-hydro", "lf");
      // legZTPbPbppSingleMod_Theory[iCen]->AddEntry((TObject *)0, "X. N. Wang et al.", "");
    }
    else if (iCen == 1)
    {
      legZTPbPbppSingleMod_Theory[iCen] = LegStd(legZTPbPbppSingleMod_Theory[iCen], 0.12, 0.16, 0.350, 0.41);
      legZTPbPbppSingleMod_Theory[iCen]->AddEntry(hSystZt[iCen], "Data", "epf");
      legZTPbPbppSingleMod_Theory[iCen]->AddEntry(grDztNLOmedianppSyst, "NLO pQCD, pp", "lf");
      legZTPbPbppSingleMod_Theory[iCen]->AddEntry(grDztNLOmedian[1], "NLO pQCD#color[0]{.}+#color[0]{.}#Delta#it{E}_{loss}", "lf");
      // legZTPbPbppSingleMod_Theory[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      //  legZTPbPbppSingleMod_Theory[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
      legZTPbPbppSingleMod_Theory[iCen]->AddEntry(histDztPbPbCOLBTmedianSyst[1], "CoLBT-hydro", "lf");
      // legZTPbPbppSingleMod_Theory[iCen]->AddEntry((TObject *)0, "X. N. Wang et al.", "");
    }
    else if (iCen == 2)
    {
      legZTPbPbppSingleMod_Theory[iCen] = LegStd(legZTPbPbppSingleMod_Theory[iCen], 0.12, 0.16, 0.350, 0.34);
      legZTPbPbppSingleMod_Theory[iCen]->AddEntry(hSystZt[iCen], "Data", "epf");
      legZTPbPbppSingleMod_Theory[iCen]->AddEntry(grDztNLOmedianppSyst, "NLO pQCD, pp", "lf");
      legZTPbPbppSingleMod_Theory[iCen]->AddEntry(grDztNLOmedian[2], "NLO pQCD#color[0]{.}+#color[0]{.}#Delta#it{E}_{loss}", "lf");
      // legZTPbPbppSingleMod_Theory[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      //  legZTPbPbppSingleMod_Theory[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
    }
    legZTPbPbppSingleMod_Theory[iCen]->Draw("same");

    TLatex *latdiffCent_Theory = LatexStdISO(latdiffCent_Theory, 0.500, 0.92, 0.04, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true);
    // cDiffCent_Pythia[iCen]->Print(dirPlot + Form("/withModelNoppztPbPb_pp_Cen%d_%d%s.pdf", cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
    cDiffCent_Theory[iCen]->Print(dirPlot + Form("/Dzt_Cen%d_%d%s.pdf", cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
  }

  TCanvas *cDiffCent_TheoryAllTest = new TCanvas("cDiffCent_TheoryAllTest", "", 1920, 900);
  cDiffCent_TheoryAllTest->SetFillStyle(4000);

  // Margins
  Float_t lMargin2 = 0.08;
  Float_t rMargin2 = 0.02;
  Float_t bMargin2 = 0.12;
  Float_t tMargin2 = 0.05;

  // Canvas setup
  CanvasPartition(cDiffCent_TheoryAllTest, Nx, Ny, lMargin2, rMargin2, bMargin2, tMargin2);

  TPad *padZt[Nx][Ny];
  TH1F *hDztXplot[nCen];
  TH1F *hDztXplotSyst[nCen];

  for (Int_t iCen = 0; iCen < Nx; iCen++)
  {
    for (Int_t j = 0; j < Ny; j++)
    {
      cDiffCent_TheoryAllTest->cd(0);
      // Get the pads previously created.
      padZt[iCen][j] = (TPad *)cDiffCent_TheoryAllTest->FindObject(TString::Format("pad_%d_%d", iCen, j).Data());
      padZt[iCen][j]->Draw();
      padZt[iCen][j]->SetFillStyle(4000);
      padZt[iCen][j]->SetFrameFillStyle(4000);
      padZt[iCen][j]->cd();
      gPad->SetLogy();

      // Size factors
      Float_t xFactor2 = padZt[0][0]->GetAbsWNDC() / padZt[iCen][j]->GetAbsWNDC();
      Float_t yFactor2 = padZt[0][0]->GetAbsHNDC() / padZt[iCen][j]->GetAbsHNDC();

      TH1F *hFrame2 = (TH1F *)hGeneral->Clone(TString::Format("h_%d_%d", iCen, j).Data());

      hFrame2->GetYaxis()->SetRangeUser(1e-3, 30);
      hFrame2->GetYaxis()->SetLabelSize(0.045);
      hFrame2->GetYaxis()->SetLabelOffset(0.02);
      // hFrame2->GetYaxis()->SetTitleFont(43);
      hFrame2->GetYaxis()->SetTitleSize(0.05);
      hFrame2->GetYaxis()->SetTitleOffset(1.8);

      // Format for x axis
      hFrame2->GetXaxis()->SetRangeUser(-0.01, 1.0005);
      hFrame2->GetXaxis()->SetLabelFont(43);
      hFrame2->GetXaxis()->SetLabelSize(30);
      hFrame2->GetXaxis()->SetLabelOffset(0.02);
      hFrame2->GetXaxis()->SetTitleFont(43);
      hFrame2->GetXaxis()->SetTitleSize(34);
      hFrame2->GetXaxis()->SetTitleOffset(1.25);

      hFrame2->Draw();

      // Plot overlap

      hDztXplot[iCen] = (TH1F *)hZtCent[iCen]->Clone(Form("hDztXplot%d_%d", iCen, iCen + 1));
      hDztXplotSyst[iCen] = (TH1F *)hSystZt[iCen]->Clone(Form("hDztXplotSyst%d_%d", iCen, iCen + 1));

      PlotStyle(hDztXplot[iCen], 20, 1, kGray + 3, kColorMarkFill[iCen], "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
      PlotStyle(hDztXplotSyst[iCen], 20, 1, kGray + 3, 16, "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", true);
      hDztXplotSyst[iCen]->Draw("E2P same ");
      hDztXplot[iCen]->Draw("EPX0same");

      grDztNLOmedian[iCen]->SetLineWidth(3);
      grDztNLOmedian[iCen]->Draw("pl3 same");
      grDztNLOmedianpp->SetLineStyle(1);
      grDztNLOmedianpp->SetLineWidth(3);
      grDztNLOmedianppSyst->SetLineStyle(1);

      //grDztNLOmedianppSyst->Draw("X0sameE3");
      grDztNLOmedianppSyst->Draw("pl3 same");
      grDztNLOmedianpp->Draw("hist l same ");
      if (iCen == 0 || iCen == 1)
      {
        // histDztPbPbCOLBTmedianSyst[iCen]->SetFillColorAlpha(kAzure-3, 0.20);
        histDztPbPbCOLBTmedianSyst[iCen]->SetLineStyle(4);
        histDztPbPbCOLBTmedianSyst[iCen]->SetLineWidth(3);
        // histDztPbPbCOLBTmedianSyst[iCen]->SetLineColor(kAzure-1);
        histDztPbPbCOLBTmedianSyst[iCen]->Draw("X0sameE3 ");
        histDztPbPbCOLBTmedian[iCen]->SetFillStyle(0);
        histDztPbPbCOLBTmedian[iCen]->SetLineStyle(4);
        // histDztPbPbCOLBTmedian[iCen]->SetLineColor(kAzure-1);
        histDztPbPbCOLBTmedian[iCen]->SetLineWidth(3);
        histDztPbPbCOLBTmedian[iCen]->Draw(" HISTSAMEL ");
      }
      if (iCen == 1)
      {
        legZTPbPbppSingleMod_Theory[iCen] = LegStd(legZTPbPbppSingleMod_Theory[iCen], 0.01, 0.16, 0.66, 0.40);
        legZTPbPbppSingleMod_Theory[iCen]->SetTextSize(0.048);
        legZTPbPbppSingleMod_Theory[iCen]->AddEntry(hDztXplotSyst[iCen], Form(" Data "), "epf");
        legZTPbPbppSingleMod_Theory[iCen]->AddEntry(grDztNLOmedianppSyst, "NLO pQCD, pp", "lf");
        legZTPbPbppSingleMod_Theory[iCen]->AddEntry(grDztNLOmedian[1], "NLO pQCD#color[0]{.}+#color[0]{.}#Delta#it{E}_{loss}", "lf");
        legZTPbPbppSingleMod_Theory[iCen]->AddEntry(histDztPbPbCOLBTmedianSyst[1], "CoLBT#font[122]{-}hydro", "lf");
        legZTPbPbppSingleMod_Theory[iCen]->Draw("same");
      }
      if (iCen == 0)
      {
        TLatex *latTitle = new TLatex();
        latTitle->SetTextFont(42);
        latTitle->SetTextSize(0.042);
        latTitle->SetNDC();
        latTitle->DrawLatex(0.24, 0.94, Form("#font[42]{ALICE} Pb#font[122]{-}Pb,#color[0]{.}#sqrt{#it{s}_{NN}} = 5.02 TeV"));
        latTitle->DrawLatex(0.860, 0.94, Form("#bf{0#font[122]{-}30%%}"));
        // latTitle->DrawLatex(0.20, 0.94 - 2*0.055, Form("%2.0f <#color[0]{..}#it{p}_{T}^{#it{#gamma}} < %2.0f GeV/#it{c}#color[0]{..};#color[0]{..}#it{p}_{T}^{h} > 0.5 GeV/#it{c} ", ptMin, ptMax));
      }
      if (iCen == 1)
      {

        TLatex *latTitle = new TLatex();
        latTitle->SetTextFont(42);
        latTitle->SetTextSize(0.052);
        latTitle->SetNDC();
        latTitle->DrawLatex(0.80, 0.94, Form("#bf{30#font[122]{-}50%%}"));
      }

      if (iCen == 2)
      {
        TLatex *latTitle = new TLatex();
        latTitle->SetTextFont(42);
        latTitle->SetTextSize(0.05);
        latTitle->SetNDC();
        latTitle->DrawLatex(0.75, 0.94, Form("#bf{50#font[122]{-}90%%}"));
        // latTitle->DrawLatex(0.14, 0.26, Form("#font[42]{ALICE} Pb#font[122]{-}Pb,#color[0]{..}#sqrt{#it{s}_{NN}} = 5.02 TeV"));
        latTitle->DrawLatex(0.04, 0.29 - 1 * 0.055, Form("|#Delta#it{#varphi}| >#color[0]{..}#frac{3}{5}#color[0]{.}#it{#pi}, |#it{#eta}^{#it{#gamma}}| < 0.67 "));
        latTitle->DrawLatex(0.04, 0.28 - 2 * 0.054, Form("%2.0f <#color[0]{..}#it{p}_{T}^{#it{#gamma}}#color[0]{.}< %2.0f GeV/#it{c}#color[0]{.};#color[0]{..}#it{p}_{T}^{h} > 0.5 GeV/#it{c} ", ptMin, ptMax));
      }
    }
  }

  cDiffCent_TheoryAllTest->Print(dirPlot + Form("/Dzt_AllCen_Theory.pdf"));

  TCanvas *cPbPbppRatio = new TCanvas(Form("cPbPbpp_Ipythia"), Form("cPbPbpp_Ipythia"), 3 * 800, 1 * 800);
  cPbPbppRatio->Divide(3, 1);
  TLatex *lat1[nCen];
  TPad *pad1[nCen];
  TPad *pad2[nCen];

  PlotStyle(hGeneralRatio, 20, 1, kWhite, kWhite, " #it{z}_{T} ", "Ratio", false);

  TCanvas *cPbPbppRatioSingle[nCen];
  TPad *pad1Single[nCen];
  TPad *pad2Single[nCen];
  TLatex *lat1Single[nCen];
  TLegend *legZTPbPbppSingle[nCen];
  TLegend *legZTPbPbppSingleMod[nCen];
  TLegend *legRatioSingleMod[nCen];
  TLegend *legNLO[nCen];
  //  PlotStyle(hGeneralRatio, 20, 1, kWhite, kWhite,"#it{z}_{T}", "1 / #it{N}^{ #it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    legNLO[iCen] = LegStd(legNLO[iCen], 0.16, 0.80, 0.35, 0.995);
    legNLO[iCen]->SetTextSize(0.06);
    cPbPbppRatioSingle[iCen] = new TCanvas(Form("cPbPbppSingle%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("cPbPbppSingle%d_%d", cenBins[iCen], cenBins[iCen + 1]), 600, 800);
    pad1Single[iCen] = new TPad(Form("pad1Single%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("pad1Single%d_%d", cenBins[iCen], cenBins[iCen + 1]), 0, 0.35, 1, 1.0);
    pad1Single[iCen]->SetTopMargin(0.02);
    pad1Single[iCen]->SetBottomMargin(0); // Upper and lower plot are joined
    pad1Single[iCen]->SetLeftMargin(0.15);
    pad1Single[iCen]->SetRightMargin(0.02);
    // pad1[iCen]->SetGridx();         // Vertical grid
    // pad1[iCen]->SetTicks(0,0);
    pad1Single[iCen]->Draw(); // Draw the upper pad: pad1
    pad1Single[iCen]->cd();
    gPad->SetLogy();
    gPad->SetTickx();

    hGeneral->GetYaxis()->SetRangeUser(2 * 1e-4, 99);
    hGeneral->GetYaxis()->SetLabelSize(0.05);
    hGeneral->GetYaxis()->SetTitleSize(0.055);
    hGeneral->GetYaxis()->SetTitleOffset(1.1);
    hGeneral->GetYaxis()->SetLabelOffset(0.01);
    hGeneral->Draw("same");
    // histDztPbPbCOLBTmedian[0]->SetMarkerStyle(0);
    // histDztPbPbCOLBTmedianSyst[0]->SetMarkerStyle(0);
    hSystZt[iCen]->Draw("E2P same");
    hZtCent[iCen]->Draw("EPX0same");
    hZt_MC_Gen[0]->SetLineWidth(5);
    hZt_MC_Gen[0]->Draw("HISTSAMEC ");
    grDztNLOmedian[iCen]->Draw("pl3 same");
    if (iCen == 0)
    {
      histDztPbPbCOLBTmedianSyst[0]->Draw("X0sameE3 ");
      histDztPbPbCOLBTmedian[0]->SetFillStyle(0);
      histDztPbPbCOLBTmedian[0]->Draw(" HISTSAMEL ");
    }
    else if (iCen == 1)
    {
      histDztPbPbCOLBTmedianSyst[1]->Draw("X0sameE3 ");
      histDztPbPbCOLBTmedian[1]->SetFillStyle(0);
      histDztPbPbCOLBTmedian[1]->Draw(" HISTSAMEL ");
    }

    // hZt_MC_Rec[iCen]->Draw("same");
    hSystZt[iCen]->SetTitle("");
    hSystZt[iCen]->GetXaxis()->SetLabelSize(0);
    hSystZt[iCen]->GetXaxis()->SetTitleSize(0);
    // hSystZt[iCen]->GetXaxis()->SetTickSize(0.015);
    lat1Single[iCen] = LatexStdISORatio(lat1Single[iCen], 0.32, 0.92, 0.045, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true);
    //  lat1->SetNDC();
    //  lat->SetTextSize(32);
    legZTPbPbppSingle[iCen] = LegStd(legZTPbPbppSingle[iCen], 0.650, 0.5, 0.830, 0.665);
    legZTPbPbppSingle[iCen]->AddEntry(hZt_MC_Gen[iCen], Form("PYTHIA 8 pp"), "l");
    legZTPbPbppSingle[iCen]->AddEntry(hZtCent[iCen], "Pb#font[122]{-}Pb stat. unc. ", "ep");
    legZTPbPbppSingle[iCen]->AddEntry(hSystZt[iCen], "Pb#font[122]{-}Pb syst. unc. ", "f");

    if (iCen == 0)
    {
      legZTPbPbppSingleMod[iCen] = LegStd(legZTPbPbppSingleMod[iCen], 0.170, 0.02, 0.350, 0.305);
      legZTPbPbppSingleMod[iCen]->AddEntry(grDztNLOmedian[0], "NLO pQCD,", "lf");
      legZTPbPbppSingleMod[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      // legZTPbPbppSingleMod[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
      legZTPbPbppSingleMod[iCen]->AddEntry(histDztPbPbCOLBTmedianSyst[0], "CoLBT-hydro,", "lf");
      // legZTPbPbppSingleMod[iCen]->AddEntry((TObject *)0, "X. N. Wang et al.", "");
    }
    else if (iCen == 1)
    {
      legZTPbPbppSingleMod[iCen] = LegStd(legZTPbPbppSingleMod[iCen], 0.170, 0.02, 0.350, 0.305);
      legZTPbPbppSingleMod[iCen]->AddEntry(grDztNLOmedian[1], "NLO pQCD,", "lf");
      legZTPbPbppSingleMod[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      // legZTPbPbppSingleMod[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
      legZTPbPbppSingleMod[iCen]->AddEntry(histDztPbPbCOLBTmedianSyst[1], "CoLBT-hydro,", "lf");
      // legZTPbPbppSingleMod[iCen]->AddEntry((TObject *)0, "X. N. Wang et al.", "");
    }
    else if (iCen == 2)
    {
      legZTPbPbppSingleMod[iCen] = LegStd(legZTPbPbppSingleMod[iCen], 0.170, 0.02, 0.350, 0.195);
      legZTPbPbppSingleMod[iCen]->AddEntry(grDztNLOmedian[2], "NLO pQCD,", "lf");
      legZTPbPbppSingleMod[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      // legZTPbPbppSingleMod[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
    }
    legZTPbPbppSingle[iCen]->Draw("same");
    legZTPbPbppSingleMod[iCen]->Draw("same");
    // legRatioSingleMod[iCen] = LegStd(legRatioSingleMod[iCen], 0.40, 0.80, 0.80, 0.995);
    legRatioSingleMod[iCen] = LegStd(legRatioSingleMod[iCen], 0.45, 0.88, 0.65, 0.995);
    legRatioSingleMod[iCen]->SetTextSize(0.06);
    cPbPbppRatioSingle[iCen]->cd(); // Go back to the main canvas before defining pad2
    pad2Single[iCen] = new TPad(Form("pad2Single%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("pad2Single%d_%d", cenBins[iCen], cenBins[iCen + 1]), 0, 0.0, 1, 0.35);
    pad2Single[iCen]->SetTopMargin(0);
    pad2Single[iCen]->SetLeftMargin(0.15);
    pad2Single[iCen]->SetBottomMargin(0.22);
    pad2Single[iCen]->SetRightMargin(0.02);
    // pad2[iCen] ->SetGridx(); // vertical grid
    pad2Single[iCen]->Draw();
    pad2Single[iCen]->cd();
    pad2Single[iCen]->SetTicks();

    h3[iCen]->Sumw2();
    h3[iCen]->SetStats(0);
    hsyst[iCen]->SetLineColor(kColorMark[iCen]);
    hsyst[iCen]->SetMarkerColor(kColorMark[iCen]);
    hGeneralRatio->GetYaxis()->SetRangeUser(0, 1.64);
    if (iCen == 2)
      hGeneralRatio->GetYaxis()->SetRangeUser(0, 1.64);
    hGeneralRatio->GetYaxis()->SetTitleSize(0.1);
    hGeneralRatio->GetYaxis()->SetTitleOffset(0.7);
    hGeneralRatio->GetYaxis()->SetLabelSize(0.082);
    hGeneralRatio->GetYaxis()->SetLabelOffset(0.01);
    hGeneralRatio->GetYaxis()->CenterTitle(true);

    // X axis ratio plot settings
    hGeneralRatio->SetTitle("");
    hGeneralRatio->GetYaxis()->SetNdivisions(505);
    hGeneralRatio->GetXaxis()->SetTitleFont(42);
    hGeneralRatio->GetXaxis()->SetTitleSize(0.09);
    hGeneralRatio->GetXaxis()->SetTitleOffset(1);
    hGeneralRatio->GetXaxis()->SetLabelFont(42); // Absolute font size in pixel (precision 3)
    hGeneralRatio->GetXaxis()->SetLabelSize(0.08);
    hGeneralRatio->GetXaxis()->SetLabelOffset(0.01);
    hsyst[iCen]->SetMinimum(-0.1); // Define Y ..
    hsyst[iCen]->SetMaximum(1.2);
    hGeneralRatio->Draw("hist same");
    hsyst[iCen]->Draw("samee2");
    h3[iCen]->Draw("X0SAMEpe");
    grIaaNLOmedian[iCen]->Draw("pl3 same");
    hIaa_nPDF->Draw("hist same c");
    // legPbPb_NLOratioAll_IaanPDF->AddEntry(hIaa_nPDF, "#it{I}_{AA} + CNM", "l");
    legNLO[iCen]->AddEntry(h3[iCen], "#frac{Pb#font[122]{-}Pb}{PYTHIA 8 pp}", "ep");
    // legRatioSingleMod[iCen]->AddEntry(grIaaNLOmedian[iCen], "#it{I}_{AA}: NLO pQCD (X. N. Wang and M. Xie)", "lf");
    if (iCen == 0)
    {

      histIaaCOLBTmedianSyst[0]->Draw("X0sameE3 ");
      histIaaCOLBTmedian[0]->SetFillStyle(0);
      histIaaCOLBTmedian[0]->Draw("histsameL ");
      legRatioSingleMod[iCen]->AddEntry(histIaaCOLBTmedianSyst[0], "#it{I}_{AA}: CoLBT-hydro (X. N. Wang, et al.)", "lf");
    }
    else if (iCen == 1)
    {
      histIaaCOLBTmedianSyst[1]->Draw("X0sameE3 ");
      histIaaCOLBTmedian[1]->SetFillStyle(0);
      histIaaCOLBTmedian[1]->Draw("histsameL ");
      legRatioSingleMod[iCen]->AddEntry(histIaaCOLBTmedianSyst[1], "#it{I}_{AA}: CoLBT-hydro (X. N. Wang, et al.)", "lf");
    }
    // else if (iCen == 2)
    //   legRatioSingleMod[iCen]->SetY1(0.8975);

    legNLO[iCen]->Draw("same");
    legRatioSingleMod[iCen]->Draw("same");
    //  cPbPbpp->Update();
    //  pad2->Update();
    //  TLine line(0.0, 1., 1.2, 1.);

    // Y axis ratio plot settings

    hsyst[iCen]->SetTitle("");
    hsyst[iCen]->GetYaxis()->SetNdivisions(510);
    hsyst[iCen]->GetXaxis()->SetLabelSize(16);
    // hsyst[iCen]->GetYaxis()->SetLabelSize(0.1);
    hsyst[iCen]->GetYaxis()->SetTitleSize(0.08);
    hsyst[iCen]->GetYaxis()->SetTitle("#frac{Data}{PYTHIA}");
    hsyst[iCen]->GetYaxis()->CenterTitle(true);
    hsyst[iCen]->GetYaxis()->SetTitleOffset(0.52);

    // X axis ratio plot settings

    hsyst[iCen]->GetXaxis()->SetTitleSize(20);
    hsyst[iCen]->GetXaxis()->SetTitleFont(43);
    hsyst[iCen]->GetXaxis()->SetTitleOffset(1.2);
    hsyst[iCen]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hsyst[iCen]->GetXaxis()->SetLabelSize(20);
    hsyst[iCen]->GetXaxis()->SetTitle("#it{z}_{T}/#Delta#it{z}_{T}");
    // h3[iCen]->GetXaxis()->SetTitle("z_{#rm{T}} = p_{#rm{T}}^{had}}/p_{#rm{T}}^{#gamma}}");

    TGraph *linea = DrawLine(linea, 0, 0.5, 1.2, 0.5);
    linea->Draw("l");

    TGraph *lineb = DrawLine(lineb, 0, 1, 1.2, 1);
    lineb->Draw("l");
    cPbPbppRatioSingle[iCen]->Print(dirPlot + Form("/SingleztPbPb_pp_Cen%d_%d%s.pdf", cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
  }

  TCanvas *cPbPbppRatioSingleNLOpQCD[nCen];
  TPad *pad1SingleNLOpQCD[nCen];
  TPad *pad2SingleNLOpQCD[nCen];
  TLatex *lat1SingleNLOpQCD[nCen];

  TLegend *legZTPbPbppSingleModNLO[nCen];
  TLegend *legNLOpQCD[nCen];
  TLegend *legZTPbPbppSingleNLOpQCD[nCen];
  TLegend *legRatioSinglePadIpQCDIaaCNM[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    // legNLOpQCD[iCen] = LegStd(legNLOpQCD[iCen], 0.16, 0.80, 0.35, 0.995);
    // legNLOpQCD[iCen]->SetTextSize(0.06);
    legRatioSinglePadIpQCDIaaCNM[iCen] = LegStd(legRatioSinglePadIpQCDIaaCNM[iCen], 0.320, 0.70, 0.92, 0.98);
    legRatioSinglePadIpQCDIaaCNM[iCen]->SetTextSize(0.056);
    legRatioSinglePadIpQCDIaaCNM[iCen]->SetNColumns(2);
    cPbPbppRatioSingleNLOpQCD[iCen] = new TCanvas(Form("cPbPbppSingleNLOpQCD%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("cPbPbppSingleNLOpQCD%d_%d", cenBins[iCen], cenBins[iCen + 1]), 600, 800);
    pad1SingleNLOpQCD[iCen] = new TPad(Form("pad1SingleNLOpQCD%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("pad1SingleNLOpQCD%d_%d", cenBins[iCen], cenBins[iCen + 1]), 0, 0.35, 1, 1.0);
    pad1SingleNLOpQCD[iCen]->SetTopMargin(0.02);
    pad1SingleNLOpQCD[iCen]->SetBottomMargin(0); // Upper and lower plot are joined
    pad1SingleNLOpQCD[iCen]->SetLeftMargin(0.15);
    pad1SingleNLOpQCD[iCen]->SetRightMargin(0.02);
    // pad1[iCen]->SetGridx();         // Vertical grid
    // pad1[iCen]->SetTicks(0,0);
    pad1SingleNLOpQCD[iCen]->Draw(); // Draw the upper pad: pad1
    pad1SingleNLOpQCD[iCen]->cd();
    gPad->SetLogy();
    gPad->SetTickx();

    hGeneral->GetYaxis()->SetRangeUser(2 * 1e-4, 99);
    hGeneral->GetYaxis()->SetLabelSize(0.05);
    hGeneral->GetYaxis()->SetTitleSize(0.055);
    hGeneral->GetYaxis()->SetTitleOffset(1.1);
    hGeneral->GetYaxis()->SetLabelOffset(0.01);
    hGeneral->Draw("same");
    hSystZt[iCen]->Draw("E2P same");
    hZtCent[iCen]->Draw("EPX0same");
    // grDztNLOmedianpp->SetLineWidth(5);
    grDztNLOmedianpp->Draw("histsame");
    grDztNLOmedian[iCen]->Draw("pl3 same");
    if (iCen == 0)
    {
      histDztPbPbCOLBTmedianSyst[0]->Draw("X0sameE3 ");
      histDztPbPbCOLBTmedian[0]->SetFillStyle(0);
      histDztPbPbCOLBTmedian[0]->Draw(" HISTSAMEL ");
    }
    else if (iCen == 1)
    {
      histDztPbPbCOLBTmedianSyst[1]->Draw("X0sameE3 ");
      histDztPbPbCOLBTmedian[1]->SetFillStyle(0);
      histDztPbPbCOLBTmedian[1]->Draw(" HISTSAMEL ");
    }

    // hZt_MC_Rec[iCen]->Draw("same");
    hSystZt[iCen]->SetTitle("");
    hSystZt[iCen]->GetXaxis()->SetLabelSize(0);
    hSystZt[iCen]->GetXaxis()->SetTitleSize(0);
    // hSystZt[iCen]->GetXaxis()->SetTickSize(0.015);
    lat1SingleNLOpQCD[iCen] = LatexStdISORatio(lat1SingleNLOpQCD[iCen], 0.32, 0.92, 0.045, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true);
    //  lat1->SetNDC();
    //  lat->SetTextSize(32);
    legZTPbPbppSingleNLOpQCD[iCen] = LegStd(legZTPbPbppSingleNLOpQCD[iCen], 0.650, 0.5, 0.830, 0.665);
    // legZTPbPbppSingleNLOpQCD[iCen]->AddEntry(grDztNLOmedianpp, Form("NLO pQCD pp"), "l");
    legZTPbPbppSingleNLOpQCD[iCen]->AddEntry(hZtCent[iCen], "Pb#font[122]{-}Pb stat. unc.", "ep");
    legZTPbPbppSingleNLOpQCD[iCen]->AddEntry(hSystZt[iCen], "Pb#font[122]{-}Pb syst. unc.", "f");

    if (iCen == 0)
    {
      legZTPbPbppSingleModNLO[iCen] = LegStd(legZTPbPbppSingleModNLO[iCen], 0.170, 0.02, 0.350, 0.315);
      legZTPbPbppSingleModNLO[iCen]->AddEntry(grDztNLOmedianppSyst, Form("NLO pQCD, pp"), "lf");
      legZTPbPbppSingleModNLO[iCen]->AddEntry(grDztNLOmedian[0], "NLO pQCD#color[0]{.}+#color[0]{.}#Delta#it{E}_{loss}, Pb#font[122]{-}Pb", "lf");
      legZTPbPbppSingleModNLO[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      // legZTPbPbppSingleModNLO[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
      legZTPbPbppSingleModNLO[iCen]->AddEntry(histDztPbPbCOLBTmedianSyst[0], "CoLBT-hydro, Pb#font[122]{-}Pb,", "lf");
      // legZTPbPbppSingleModNLO[iCen]->AddEntry((TObject *)0, "X. N. Wang et al.", "");
    }
    else if (iCen == 1)
    {
      legZTPbPbppSingleModNLO[iCen] = LegStd(legZTPbPbppSingleModNLO[iCen], 0.170, 0.02, 0.350, 0.315);
      legZTPbPbppSingleModNLO[iCen]->AddEntry(grDztNLOmedianppSyst, Form("NLO pQCD, pp"), "lf");
      legZTPbPbppSingleModNLO[iCen]->AddEntry(grDztNLOmedian[1], "NLO pQCD#color[0]{.}+#color[0]{.}#Delta#it{E}_{loss}, Pb#font[122]{-}Pb", "lf");
      legZTPbPbppSingleModNLO[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      // legZTPbPbppSingleModNLO[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
      legZTPbPbppSingleModNLO[iCen]->AddEntry(histDztPbPbCOLBTmedianSyst[1], "CoLBT-hydro, Pb#font[122]{-}Pb,", "lf");
      // legZTPbPbppSingleModNLO[iCen]->AddEntry((TObject *)0, "X. N. Wang et al.", "");
    }
    else if (iCen == 2)
    {
      legZTPbPbppSingleModNLO[iCen] = LegStd(legZTPbPbppSingleModNLO[iCen], 0.170, 0.134, 0.350, 0.315);
      legZTPbPbppSingleModNLO[iCen]->AddEntry(grDztNLOmedianppSyst, Form("NLO pQCD, pp"), "lf");
      legZTPbPbppSingleModNLO[iCen]->AddEntry(grDztNLOmedian[2], "NLO pQCD#color[0]{.}+#color[0]{.}#Delta#it{E}_{loss}, Pb#font[122]{-}Pb", "lf");
      legZTPbPbppSingleModNLO[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      // legZTPbPbppSingleModNLO[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
    }

    legZTPbPbppSingleNLOpQCD[iCen]->Draw("same");
    legZTPbPbppSingleModNLO[iCen]->Draw("same");
    cPbPbppRatioSingleNLOpQCD[iCen]->cd(); // Go back to the main canvas before defining pad2
    pad2SingleNLOpQCD[iCen] = new TPad(Form("pad2SingleNLOpQCD%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("pad2SingleNLOpQCD%d_%d", cenBins[iCen], cenBins[iCen + 1]), 0, 0.0, 1, 0.35);
    pad2SingleNLOpQCD[iCen]->SetTopMargin(0);
    pad2SingleNLOpQCD[iCen]->SetLeftMargin(0.15);
    pad2SingleNLOpQCD[iCen]->SetBottomMargin(0.22);
    pad2SingleNLOpQCD[iCen]->SetRightMargin(0.02);
    // pad2[iCen] ->SetGridx(); // vertical grid
    pad2SingleNLOpQCD[iCen]->Draw();
    pad2SingleNLOpQCD[iCen]->cd();
    pad2SingleNLOpQCD[iCen]->SetTicks();

    h3[iCen]->Sumw2();
    h3[iCen]->SetStats(0);
    hsyst[iCen]->SetLineColor(kColorMark[iCen]);
    hsyst[iCen]->SetMarkerColor(kColorMark[iCen]);
    hGeneralRatio->GetYaxis()->SetRangeUser(0.05, 2.1);
    // hGeneralRatio->GetYaxis()->SetRangeUser(0, 1.64);
    // if (iCen == 2)
    //   hGeneralRatio->GetYaxis()->SetRangeUser(0, 1.64);
    hGeneralRatio->GetYaxis()->SetTitleSize(0.1);
    hGeneralRatio->GetYaxis()->SetTitleOffset(0.7);
    hGeneralRatio->GetYaxis()->SetLabelSize(0.082);
    hGeneralRatio->GetYaxis()->SetLabelOffset(0.01);
    hGeneralRatio->GetYaxis()->CenterTitle(true);

    // X axis ratio plot settings
    hGeneralRatio->SetTitle("");
    hGeneralRatio->GetYaxis()->SetNdivisions(505);
    // hGeneralRatio->GetYaxis()->SetTitle("Pb#font[122]{-}Pb/NLO pQCD pp");
    hGeneralRatio->GetYaxis()->SetTitle("#it{I}_{pQCD}, #it{I}_{AA}");
    hGeneralRatio->GetXaxis()->SetTitleFont(42);
    hGeneralRatio->GetXaxis()->SetTitleSize(0.09);
    hGeneralRatio->GetXaxis()->SetTitleOffset(1);
    hGeneralRatio->GetXaxis()->SetLabelFont(42); // Absolute font size in pixel (precision 3)
    hGeneralRatio->GetXaxis()->SetLabelSize(0.08);
    hGeneralRatio->GetXaxis()->SetLabelOffset(0.01);
    hsyst[iCen]->SetMinimum(-0.1); // Define Y ..
    hsyst[iCen]->SetMaximum(1.2);
    hGeneralRatio->Draw("hist same");
    hPbPb_NLOSyst[iCen]->Draw("E2Psame ");
    hPbPb_NLO[iCen]->Draw("EPX0same");
    grIaaNLOmedian[iCen]->Draw("pl3 same");
    hIaa_nPDF->Draw("hist same c");
    legRatioSinglePadIpQCDIaaCNM[iCen]->AddEntry(h3[iCen], "#it{I}_{pQCD}", "ep");
    legRatioSinglePadIpQCDIaaCNM[iCen]->AddEntry(hIaa_nPDFSyst, "#it{I}_{AA}, NLO pQCD + CNM", "lf");
    legRatioSinglePadIpQCDIaaCNM[iCen]->AddEntry(grIaaNLOmedian[iCen], "#it{I}_{AA}, NLO pQCD#color[0]{.}+#color[0]{.}#Delta#it{E}_{loss}", "lf");
    if (iCen == 0)
    {
      histIaaCOLBTmedianSyst[0]->Draw("X0sameE3 ");
      histIaaCOLBTmedian[0]->SetFillStyle(0);
      histIaaCOLBTmedian[0]->Draw("histsameL ");
      legRatioSinglePadIpQCDIaaCNM[iCen]->AddEntry(histIaaCOLBTmedianSyst[0], "#it{I}_{AA}, CoLBT-hydro", "lf");
    }
    else if (iCen == 1)
    {
      histIaaCOLBTmedianSyst[1]->Draw("X0sameE3 ");
      histIaaCOLBTmedian[1]->SetFillStyle(0);
      histIaaCOLBTmedian[1]->Draw("histsameL ");
      legRatioSinglePadIpQCDIaaCNM[iCen]->AddEntry(histIaaCOLBTmedianSyst[1], "#it{I}_{AA}, CoLBT-hydro", "lf");
    }

    // legNLOpQCD[iCen]->Draw("same");
    //   cPbPbpp->Update();
    legRatioSinglePadIpQCDIaaCNM[iCen]->Draw("same");

    TGraph *linea = DrawLine(linea, 0, 0.5, 1.2, 0.5);
    linea->Draw("l");
    TGraph *lineb = DrawLine(lineb, 0, 1, 1.2, 1);
    lineb->Draw("l");
    cPbPbppRatioSingleNLOpQCD[iCen]->Print(dirPlot + Form("/SingleztPbPb_pp_Cen%d_%d%s_ppNLOpQCD.pdf", cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
  }
  /////////////////////////////////////////////////////////
  ////////////// Dzt all centralities /////////////////////
  /////////////////////////////////////////////////////////

  TCanvas *cAllZt = canvasStd("cAllZt", 1, 1);
  TLegend *legdiffcenZtPYTHIA = LegStd(legdiffcenZtPYTHIA, 0.12, 0.30, 0.44, 0.35);
  TLegend *legdiffcenZt = LegStd(legdiffcenZt, 0.12, 0.16, 0.35, 0.35); //0.12, 0.16, 0.350, 0.41
  gPad->SetLogy();
  hGeneral->GetYaxis()->SetTitleSize(0.055);
  hGeneral->GetXaxis()->SetRangeUser(0.05, 1.05);
  hGeneral->GetXaxis()->SetLabelSize(0.040);
  hGeneral->GetXaxis()->SetTitleSize(0.055);
  hGeneral->GetYaxis()->SetRangeUser(2 * 1e-4, 99);
  hGeneral->GetYaxis()->SetLabelSize(0.04);
  hGeneral->GetYaxis()->SetTitleSize(0.045);
  hGeneral->GetYaxis()->SetTitleOffset(1.0);
  hGeneral->Draw("same");
  hZt_MC_Gen[0]->SetFillStyle(0);
  hZt_MC_Gen[0]->Draw("HISTSAMEC ");
  legdiffcenZtPYTHIA->AddEntry(hZt_MC_Gen[0], "PYTHIA 8 pp ", "l");
  // legdiffcenZt->SetNColumns(2);
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    legdiffcenZt->AddEntry(hSystZt[iCen], Form("%d#font[122]{-}%d%%", cenBins[iCen], cenBins[iCen + 1]), "epf");
    // legdiffcenZt->AddEntry(hSystZt[iCen], Form("syst. unc"), "f");
  }
  legdiffcenZtPYTHIA->Draw("SAME");
  legdiffcenZt->Draw("SAME");
  hSystZt[0]->Draw("samee2");
  hZtCent[0]->Draw("EP X0 same");
  hSystZt[2]->Draw("samee2");
  hZtCent[2]->Draw("EP X0 same");
  hSystZt[1]->Draw("samee2");
  hZtCent[1]->Draw("EP X0 same");
  TLatex *latAliceZtOnly = LatexStdISORatio(latAliceZtOnly, 0.420, 0.92, 0.040, cenBins[0], cenBins[1], ptMin, ptMax, false);
  cAllZt->Print(dirPlot + Form("/ZtAllCent030%s.pdf", sPtAll.Data()));

  /////////////////////////////////////////////////////////////////
  ////////////// Dzt all centralities, pp and NLO pQCD ////////////
  /////////////////////////////////////////////////////////////////

  TCanvas *cAllZtNLOpQCD = canvasStdIaa("cAllZtNLOpQCD", 1, 1);
  cAllZtNLOpQCD->SetTopMargin(0.015);
  cAllZtNLOpQCD->SetRightMargin(0.02);
  cAllZtNLOpQCD->SetLeftMargin(0.15);
  cAllZtNLOpQCD->SetBottomMargin(0.11);
  TLegend *legdiffcenZtNLOpQCD = LegStd(legdiffcenZtNLOpQCD, 0.12, 0.35, 0.35, 0.41);
  gPad->SetLogy();
  hGeneral->GetYaxis()->SetTitleSize(0.055);
  // hGeneral->GetXaxis()->SetRangeUser(-0.05, 1.1);
  hGeneral->GetXaxis()->SetNdivisions(505);
  hGeneral->GetXaxis()->SetLabelSize(0.040);
  hGeneral->GetXaxis()->SetTitleSize(0.055);
  hGeneral->GetYaxis()->SetRangeUser(1e-3, 40);
  hGeneral->GetYaxis()->SetLabelSize(0.04);
  hGeneral->GetYaxis()->SetTitleSize(0.045);
  hGeneral->GetYaxis()->SetTitleOffset(1.0);
  hGeneral->Draw("same");
  //grDztNLOmedianppSyst->Draw("X0sameE3");
  grDztNLOmedianppSyst->Draw("pl3 same");
  grDztNLOmedianpp->Draw("pl same ");
  legdiffcenZtNLOpQCD->AddEntry(grDztNLOmedianppSyst, "NLO pQCD, pp", "lf");
  legdiffcenZtNLOpQCD->Draw("SAME");
  hSystZt[0]->Draw("samee2");
  hZtCent[0]->Draw("EP X0 same");
  hSystZt[2]->Draw("samee2");
  hZtCent[2]->Draw("EP X0 same");
  hSystZt[1]->Draw("samee2");
  hZtCent[1]->Draw("EP X0 same");
  legdiffcenZt->Draw("SAME");
  TLatex *latAliceZtOnlyNLOpQCD = LatexStdISORatio(latAliceZtOnlyNLOpQCD, 0.420, 0.92, 0.040, cenBins[0], cenBins[1], ptMin, ptMax, false);
  cAllZtNLOpQCD->Print(dirPlot + Form("/DztAllCen.pdf"));

  TCanvas *cAllZtNLOpQCD_pPbpp = canvasStdIaa("cAllZtNLOpQCD_pPbpp", 1, 1);
  TLegend *legdiffcenZtpPbpp = LegStd(legdiffcenZtpPbpp, 0.14, 0.35, 0.5025, 0.4);
  legdiffcenZtpPbpp->SetNColumns(2);
  TLegend *legdiffcenZtpPbppNLO = LegStd(legdiffcenZtpPbppNLO, 0.55, 0.20, 0.85, 0.30);
  TH1F *hGeneralDztpPbpp = new TH1F("hGeneralDztpPbpp", "hGeneralDztpPbpp", 17, assocZtThinnerPbpp);
  PlotStyle(hGeneralDztpPbpp, 20, 1, kWhite, kWhite, " #it{z}_{T} ", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
  gPad->SetLogy();
  hGeneralDztpPbpp->GetYaxis()->SetTitleSize(0.055);
  hGeneralDztpPbpp->GetXaxis()->SetRangeUser(0.05, 1.05);
  hGeneralDztpPbpp->GetXaxis()->SetLabelSize(0.040);
  hGeneralDztpPbpp->GetXaxis()->SetTitleSize(0.055);
  hGeneralDztpPbpp->GetYaxis()->SetRangeUser(2 * 1e-4, 99);
  hGeneralDztpPbpp->GetYaxis()->SetLabelSize(0.04);
  hGeneralDztpPbpp->GetYaxis()->SetTitleSize(0.045);
  hGeneralDztpPbpp->GetYaxis()->SetTitleOffset(1.0);
  hGeneralDztpPbpp->SetTitle("");
  hGeneralDztpPbpp->Draw("same");
  grDztNLOmedianpp->Draw("HISTSAMEC ");
  legdiffcenZtpPbpp->AddEntry(hDztpp, "pp stat.    ", "ep");
  legdiffcenZtpPbpp->AddEntry(hDztpp_syst, "syst. unc.", "f");
  legdiffcenZtpPbpp->Draw("SAME");
  legdiffcenZtpPbppNLO->AddEntry(grDztNLOmedianpp, "NLO pQCD, pp", "l");
  legdiffcenZtpPbppNLO->Draw("SAME");
  hDztpp_syst->SetFillStyle(0);
  hDztpp_syst->SetLineColor(kBlack);
  hDztpp_syst->SetLineWidth(1);
  hDztpp_syst->Draw("samee2");
  hDztpp->Draw("EP X0 same");
  hSystZt[0]->Draw("samee2");
  hZtCent[0]->Draw("EP X0 same");
  hSystZt[2]->Draw("samee2");
  hZtCent[2]->Draw("EP X0 same");
  hSystZt[1]->Draw("samee2");
  hZtCent[1]->Draw("EP X0 same");
  legdiffcenZt->Draw("SAME");
  TLatex *latAliceZtOnlyNLOpQCD_pPbpp = LatexStdISORatioNoPbPb(latAliceZtOnlyNLOpQCD_pPbpp, 0.420, 0.92, 0.040, cenBins[0], cenBins[1], ptMin, ptMax, false);
  cAllZtNLOpQCD_pPbpp->Print(dirPlot + Form("/ZtAllCent_pp%s.pdf", sPtAll.Data()));

  /////////////////////////////////////////////////////////////////
  ////////////// Dzt central, pp and NLO pQCD /////////////////////
  /////////////////////////////////////////////////////////////////

  TCanvas *c010ZtNLOpQCD_pPbpp = canvasStdIaa("c010ZtNLOpQCD_pPbpp", 1, 1);
  // TLegend *leg010ZtNLOpQCD_pPbpp = LegStd(leg010ZtNLOpQCD_pPbpp, 0.14, 0.20, 0.56, 0.33);
  // leg010ZtNLOpQCD_pPbpp->SetNColumns(2);
  // TLegend *leg010ZtNLOpQCD_pPbppNLO = LegStd(leg010ZtNLOpQCD_pPbppNLO, 0.55, 0.20, 0.85, 0.30);
  //  cAllZtNLOpQCD->cd();
  gPad->SetLogy();
  hGeneralDztpPbpp->GetYaxis()->SetTitleSize(0.055);
  hGeneralDztpPbpp->GetXaxis()->SetRangeUser(0.05, 1.05);
  hGeneralDztpPbpp->GetXaxis()->SetLabelSize(0.040);
  hGeneralDztpPbpp->GetXaxis()->SetTitleSize(0.055);
  hGeneralDztpPbpp->GetYaxis()->SetRangeUser(1e-3, 99);
  hGeneralDztpPbpp->GetYaxis()->SetLabelSize(0.04);
  hGeneralDztpPbpp->GetYaxis()->SetTitleSize(0.045);
  hGeneralDztpPbpp->GetYaxis()->SetTitleOffset(1.0);
  hGeneralDztpPbpp->SetTitle("");
  hGeneralDztpPbpp->Draw("same");
  grDztNLOmedianpp->SetMarkerStyle(20);
  grDztNLOmedianpp->SetMarkerSize(1.1);
  grDztNLOmedianpp->SetLineWidth(2);
  hSystZt[0]->Draw("samee2");
  hZtCent[0]->Draw("EP X0 same");
  TLegend *legDztALICE030_gamhad = LegStd(legDztALICE030_gamhad, 0.12, 0.90, 0.30, 0.95);
  legDztALICE030_gamhad->SetTextSize(0.038);
  legDztALICE030_gamhad->SetHeader("ALICE,#color[0]{..}#sqrt{#it{s}_{NN}} = 5.02 TeV");
  // legDztALICE030_gamhad->AddEntry((TObject *)0, "18#color[0]{.}<#color[0]{..}#it{p}_{T}^{#it{#gamma}}#color[0]{.}<#color[0]{..}40 GeV/#it{c}#color[0]{.};#color[0]{..}#it{p}_{T}^{h} > 0.5 GeV/#it{c}", "");
  legDztALICE030_gamhad->Draw("same");
  TLegend *legDztALICE030_pp = LegStd(legDztALICE030_pp, 0.52, 0.84, 0.52, 0.95);
  // legDztALICE030Marker->SetNColumns(2);
  legDztALICE030_pp->SetTextSize(0.036);
  legDztALICE030_pp->SetHeader("#color[0]{..}#it{p}_{T}^{iso, ch} < 1.5 GeV/#it{c}, #bf{#it{R} = 0.2}");
  legDztALICE030_pp->AddEntry((TObject *)0, "#bf{18#color[0]{.}<#color[0]{..}#it{p}_{T}^{#it{#gamma}}#color[0]{.}<#color[0]{..}40 GeV/#it{c}}#color[0]{.};#color[0]{..}#it{p}_{T}^{h} > 0.5 GeV/#it{c}", "");
  legDztALICE030_pp->Draw("same");
  TLegend *legDztALICE030Marker = LegStd(legDztALICE030Marker, 0.52, 0.70, 0.780, 0.82);
  legDztALICE030Marker->SetTextSize(0.036);
  legDztALICE030Marker->AddEntry(hPbPb_NLOSyst[0], "0#font[122]{-}30% Pb#font[122]{-}Pb", "epf");
  legDztALICE030Marker->AddEntry(grDztNLOmedianppSyst, "NLO pQCD, pp", "epf");
  legDztALICE030Marker->Draw("same");

  TLegend *legDztALICEppPb = LegStd(legDztALICEppPb, 0.12, 0.28, 0.12, 0.43);
  legDztALICEppPb->SetTextSize(0.036);
  legDztALICEppPb->SetHeader("Phys Rev C 102 (2020) 044908");
  legDztALICEppPb->AddEntry((TObject *)0, "#color[0]{..}#it{p}_{T}^{iso, ch} < 1.5 GeV/#it{c}, #bf{#it{R} = 0.4}", "");
  legDztALICEppPb->AddEntry((TObject *)0, "#bf{12#color[0]{.}<#color[0]{..}#it{p}_{T}^{#it{#gamma}}#color[0]{.}<#color[0]{..}40 GeV/#it{c}}#color[0]{.};#color[0]{..}0.5 <#color[0]{..}#it{p}_{T}^{h} < 10 GeV/#it{c}", "");
  legDztALICEppPb->Draw("same");
  TLegend *legDztALICEppPbMarker = LegStd(legDztALICEppPbMarker, 0.135, 0.16, 0.38, 0.27);
  // legDztALICEppPbMarker->SetNColumns(2);
  legDztALICEppPbMarker->SetTextSize(0.036);
  // legDztALICEppPbMarker->AddEntry(hDztpp, "pp stat.    ", "ep");
  legDztALICEppPbMarker->AddEntry(hDztpp_syst, "pp", "epf");
  // legDztALICEppPbMarker->AddEntry(hDztpPb, "p#font[122]{-}Pb stat.", "ep");
  legDztALICEppPbMarker->AddEntry(hDztpPb_syst, "p#font[122]{-}Pb", "epf");
  legDztALICEppPbMarker->Draw("same");

  hDztpp_syst->SetFillStyle(0);
  hDztpp_syst->SetLineColor(kBlack);
  hDztpp_syst->SetLineWidth(1);
  hDztpp_syst->Draw("samee2");
  hDztpp->Draw("EP X0 same");
  hDztpPb_syst->SetLineColor(kGray + 1);
  hDztpPb_syst->SetLineWidth(1);
  hDztpPb_syst->Draw("samee2");
  hDztpPb->Draw("EP X0 same");
  grDztNLOmedianppSyst->SetMarkerStyle(20);
  grDztNLOmedianppSyst->SetMarkerColor(kSpring-6);
  grDztNLOmedianppSyst->SetLineWidth(0);
  grDztNLOmedianppSyst->Draw("samepe3");
  grDztNLOmedianpp->Draw("EP X0 same");
  // TLatex *lat010ZtNLOpQCD_pPbpp = LatexStdISORatioNoPbPb(lat010ZtNLOpQCD_pPbpp, 0.420, 0.92, 0.040, cenBins[0], cenBins[1], ptMin, ptMax, false);
  c010ZtNLOpQCD_pPbpp->Print(dirPlot + Form("/Dzt_Cent0_30_pQCD_pp%s.pdf", sPtAll.Data()));

  /////////////////////////////////////////////////////////////////
  ////////////// Dzt periph, pp and NLO pQCD /////////////////////
  /////////////////////////////////////////////////////////////////

  TCanvas *c5090ZtNLOpQCD_pPbpp = canvasStdIaa("c5090ZtNLOpQCD_pPbpp", 1, 1);
  TLegend *leg5090ZtNLOpQCD_pPbpp = LegStd(leg5090ZtNLOpQCD_pPbpp, 0.14, 0.20, 0.5025, 0.33);
  leg5090ZtNLOpQCD_pPbpp->SetNColumns(2);
  TLegend *leg5090ZtNLOpQCD_pPbppNLO = LegStd(leg5090ZtNLOpQCD_pPbppNLO, 0.55, 0.20, 0.85, 0.30);
  // cAllZtNLOpQCD->cd();
  gPad->SetLogy();
  hGeneralDztpPbpp->GetYaxis()->SetTitleSize(0.055);
  hGeneralDztpPbpp->GetXaxis()->SetRangeUser(0.05, 1.05);
  hGeneralDztpPbpp->GetXaxis()->SetLabelSize(0.040);
  hGeneralDztpPbpp->GetXaxis()->SetTitleSize(0.055);
  hGeneralDztpPbpp->GetYaxis()->SetRangeUser(2 * 1e-4, 99);
  hGeneralDztpPbpp->GetYaxis()->SetLabelSize(0.04);
  hGeneralDztpPbpp->GetYaxis()->SetTitleSize(0.045);
  hGeneralDztpPbpp->GetYaxis()->SetTitleOffset(1.0);
  hGeneralDztpPbpp->SetTitle("");
  hGeneralDztpPbpp->Draw("same");
  grDztNLOmedianpp->Draw("EPSAME ");
  hSystZt[2]->Draw("samee2");
  hZtCent[2]->Draw("EP X0 same");
  leg5090ZtNLOpQCD_pPbpp->AddEntry(hDztpp, "pp stat.    ", "ep");
  leg5090ZtNLOpQCD_pPbpp->AddEntry(hDztpp_syst, "syst. unc.", "f");
  leg5090ZtNLOpQCD_pPbpp->AddEntry(hZtCent[2], Form("%d#font[122]{-}%d%% stat. ", cenBins[2], cenBins[3]), "ep");
  leg5090ZtNLOpQCD_pPbpp->AddEntry(hSystZt[2], Form("syst. unc"), "f");
  leg5090ZtNLOpQCD_pPbpp->Draw("SAME");
  leg5090ZtNLOpQCD_pPbppNLO->AddEntry(grDztNLOmedianpp, "NLO pQCD, pp", "l");
  leg5090ZtNLOpQCD_pPbppNLO->Draw("SAME");
  hDztpp_syst->SetFillStyle(0);
  hDztpp_syst->SetLineColor(kBlack);
  hDztpp_syst->SetLineWidth(1);
  hDztpp_syst->Draw("samee2");
  hDztpp->Draw("EP X0 same");
  TLatex *lat5090ZtNLOpQCD_pPbpp = LatexStdISORatioNoPbPb(lat5090ZtNLOpQCD_pPbpp, 0.420, 0.92, 0.040, cenBins[0], cenBins[1], ptMin, ptMax, false);
  c5090ZtNLOpQCD_pPbpp->Print(dirPlot + Form("/Zt50_90CentpQCD_pp%s.pdf", sPtAll.Data()));

  // TCanvas *cRatioSuppres = new TCanvas("cRatioSuppres", "cRatioSuppres", 800, 600);
  TCanvas *cRatioSuppres = canvasStd("cRatioSuppres", 1, 1);
  TLegend *legdiffcenRatio = LegStd(legdiffcenRatio, 0.61, 0.52, 0.96, 0.72);
  legdiffcenRatio->SetTextSize(0.038);
  legdiffcenRatio->SetNColumns(2);
  TLegend *legRatioPYTHIA = LegStd(legRatioPYTHIA, 0.14, 0.8, 0.35, 0.92);
  // legRatioPYTHIA->SetHeader("#it{I}_{PYTHIA} = #frac{Pb#font[122]{-}Pb}{PYTHIA 8 pp} ");
  legRatioPYTHIA->SetHeader("#it{I}_{PYTHIA} = #frac{Pb#font[122]{-}Pb}{PYTHIA 8 pp} ");
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    PlotStyle(hsyst[iCen], kMarkCen[iCen], 1, kColorMark[iCen], kColorMarkFill[iCen], "#it{z}_{T}", "#it{I}_{PYTHIA} = Pb#font[122]{-}Pb/PYTHIA", true);
    PlotStyle(h3[iCen], kMarkCen[iCen], 1, kColorMark[iCen], kColorMarkFill[iCen], "#it{z}_{T}", "#it{I}_{PYTHIA}= Pb#font[122]{-}Pb/PYTHIA", false);
    hsyst[iCen]->GetYaxis()->SetRangeUser(-0.01, 1.8);
    hsyst[iCen]->GetYaxis()->SetTitleSize(0.05);
    hsyst[iCen]->GetYaxis()->SetTitleOffset(0.85);
    hsyst[iCen]->GetYaxis()->SetLabelSize(0.035);
    // hsyst[iCen]->SetFillStyle(0);
    hGeneralRatio->SetLineWidth(0);
    hsyst[iCen]->SetLineColor(kWhite);
    hGeneralRatio->GetYaxis()->SetTitle("#it{I}_{PYTHIA}");
    // hGeneralRatio->GetYaxis()->SetTitle("Pb#font[122]{-}Pb/pp PHYTIA");
    hGeneralRatio->GetYaxis()->CenterTitle(false);
    hGeneralRatio->GetYaxis()->SetRangeUser(-0.01, 2.3);
    hGeneralRatio->GetYaxis()->SetTitleSize(0.055);
    hGeneralRatio->GetYaxis()->SetTitleOffset(1.);
    hGeneralRatio->GetYaxis()->SetLabelSize(0.040);
    hGeneralRatio->GetXaxis()->SetRangeUser(0.05, 1.05);
    hGeneralRatio->GetXaxis()->SetLabelSize(0.040);
    hGeneralRatio->GetXaxis()->SetTitleSize(0.055);
    hGeneralRatio->Draw("same");
    hsyst[iCen]->Draw("sameE2");
    h3[iCen]->Draw(" PEX0same");
    legdiffcenRatio->AddEntry(h3[iCen], Form("%d#font[122]{-}%d%% stat.", cenBins[iCen], cenBins[iCen + 1]), "ep");
    legdiffcenRatio->AddEntry(hsyst[iCen], "syst. unc.", "f");
  }

  TGraph *line = DrawLine(line, 0, 0.5, 1.2, 0.5);
  line->Draw("l");

  TGraph *line1 = DrawLine(line1, 0, 1, 1.2, 1);
  line1->Draw("l");
  // cRatioSuppres->cd();
  legRatioPYTHIA->Draw("SAME");
  legdiffcenRatio->Draw("SAME");
  TLatex *latdiffCentRatio = LatexStdISORatio(latdiffCentRatio, 0.50, 0.94, 0.040, cenBins[0], cenBins[1], ptMin, ptMax, false);
  cRatioSuppres->Print(dirPlot + Form("/RatioAllCent030%s.pdf", sPtAll.Data()));

  //////////////////////////////////////////////////////////
  //////// Comparison with other experiments //////////////
  ///////////// CMS, STAR and PHENIX /////////////////////
  ///////////////////////////////////////////////////////

  TFile *fSTAR = new TFile("OtherExpResults/HEPData-ins1442357-v1-Table_3,_au.root");
  TDirectory *dir = (TDirectory *)fSTAR->Get("Table 3, au");
  TGraphErrors *grpp = (TGraphErrors *)dir->Get("Graph1D_y1");

  double zTSTAR[7] = {0.126, 0.23, 0.334, 0.447, 0.54, 0.65, 0.75};
  double DzTSTAR[7] = {7.34329, 1.36828, 0.600217, 0.184328, 0.0643105, 0.0321283, 0.0278139};
  double DzTSTARBox[7] = {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02};
  /*for (int ibin = 0; ibin < 7; ibin++)
   {
   DzTSTARBox[ibin] = (zTSTAR[ibin + 1] - zTSTAR[ibin]) / 2;
   if (ibin == 6)
   DzTSTARBox[ibin] = DzTSTARBox[ibin - 1];
   }*/
  double DzTSTAR_Stat[7] = {0.59479, 0.557656, 0.109973, 0.0291929, 0.0345208, 0.0297178, 0.0467018};
  double DzTSTAR_SysUP[7] = {0.713706, 0.388336, 0.0975644, 0.0193355, 0.00836037, 0.00417668, 0.00361581};
  double DzTSTAR_SysDOWN[7] = {2.00138, 0.159477, 0.194714, 0.0547921, 0.0192932, 0.00963849, 0.00834417};

  TGraphAsymmErrors *grSTARSys = new TGraphAsymmErrors(7, zTSTAR, DzTSTAR, DzTSTARBox, DzTSTARBox, DzTSTAR_SysDOWN, DzTSTAR_SysUP);
  TGraphErrors *grSTAR = new TGraphErrors(7, zTSTAR, DzTSTAR, 0, DzTSTAR_Stat);
  grSTARSys->SetMarkerStyle(20);
  grSTARSys->SetMarkerSize(1.2);
  grSTARSys->SetMarkerColor(kTeal + 4);
  grSTARSys->SetLineColor(kTeal + 4);
  grSTAR->SetLineColor(kTeal + 4);
  grSTAR->SetLineWidth(2);
  grSTAR->SetMarkerStyle(20);
  grSTAR->SetMarkerSize(1.2);

  grSTAR->SetMarkerColor(kTeal + 4);

  TFile *fPHENIX = new TFile("OtherExpResults/HEPData-ins1442357-v1-Table_3,_au.root");

  double xiTPHENIX[7] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
  double zTPHENIX[7] = {0};
  double DzTPHENIX[6] = {0.00282, 0.0239, 0.149, 0.287, 0.593, 0.757};
  double DzTPHENIXBox[6] = {0.02, 0.02, 0.02, 0.02, 0.02, 0.02};
  double DzTPHENIX_Stat[6] = {0.00411, 0.00735, 0.0209, 0.0452, 0.0736, 0.18};
  double DzTPHENIX_Sys[6] = {0.00437, 0.00994, 0.0334, 0.0781, 0.0799, 0.194};

  for (int ibin = 0; ibin < 7; ibin++)
  {
    zTPHENIX[6 - ibin] = 1. / (TMath::Exp(xiTPHENIX[ibin]));
  }

  TH1F *hzTPHENIX = new TH1F("hzTPHENIX", "hzTPHENIX", 6, zTPHENIX);
  TH1F *hzTPHENIXSys = new TH1F("hzTPHENIXSys", "hzTPHENIXSys", 6, zTPHENIX);

  for (int ibin = 0; ibin < 6; ibin++)
  {
    hzTPHENIX->SetBinContent(ibin + 1, DzTPHENIX[5 - ibin]);
    hzTPHENIX->SetBinError(ibin + 1, DzTPHENIX_Stat[5 - ibin]);
    hzTPHENIXSys->SetBinContent(ibin + 1, DzTPHENIX[5 - ibin]);
    hzTPHENIXSys->SetBinError(ibin + 1, DzTPHENIX_Sys[5 - ibin]);
  }
  PlotStyle(hzTPHENIX, 20, 1.2, kAzure + 2, kWhite, "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
  TH1F *hGeneralCMS = new TH1F("hGeneral", "hGeneral", 10, -0.05, 1.1);
  PlotStyle(hGeneralCMS, 71, 1, kWhite, kWhite, "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}|d#it{z}_{T}", false);

  // Iaa STAR and PHENIX

  TH1F *hGeneralIaa = new TH1F("hGeneralIaa", "hGeneralIaa", 10, -0.05, 1.05);
  // PlotStyle(hGeneralIaa, 20, 0, kWhite, kWhite, "#it{z}_{T}", "#it{I}_{NLO pQCD}, #it{I}_{AA}", false);
  //  PlotStyle(hGeneralIaa, 20, 0, kWhite, kWhite, "#it{z}_{T}", "#it{I}_{PYTHIA}, #it{I}_{AA}", false);
  hGeneralIaa->SetDirectory(0);
  double IaaSTAR[7] = {0.734072, 0.435526, 0.415985, 0.273841, 0.429596, 0.13371, 0.281805};
  double IaaSTARBox[7] = {};
  for (int ibin = 0; ibin < 7; ibin++)
  {
    IaaSTARBox[ibin] = (zTSTAR[ibin + 1] - zTSTAR[ibin]) / 2;
    if (ibin == 6)
      IaaSTARBox[ibin] = IaaSTARBox[ibin - 1];
  }
  double IaaSTAR_Stat[7] = {0.068754, 0.179768, 0.0828287, 0.0555868, 0.346568, 0.128112, 0.483996};
  double IaaSTAR_SysUP[7] = {0.114323, 0.132715, 0.0687309, 0.03971, 0.0621083, 0.0173306, 0.105935};
  double IaaSTAR_SysDOWN[7] = {0.208097, 0.056488, 0.13582, 0.0868022, 0.130846, 0.0429338, 0.120614};

  TGraphAsymmErrors *grIaaSTARSys = new TGraphAsymmErrors(7, zTSTAR, IaaSTAR, DzTSTARBox, DzTSTARBox, IaaSTAR_SysDOWN, IaaSTAR_SysUP);
  TGraphErrors *grIaaSTAR = new TGraphErrors(7, zTSTAR, IaaSTAR, 0, IaaSTAR_Stat);
  grIaaSTAR->SetMarkerStyle(24);
  grIaaSTAR->SetMarkerColor(kPink + 8);
  grIaaSTAR->SetLineColor(kPink + 8);
  grIaaSTAR->SetLineWidth(2);
  grIaaSTARSys->SetMarkerStyle(24);
  grIaaSTARSys->SetMarkerColor(kPink + 8);
  grIaaSTARSys->SetLineColor(0);
  grIaaSTARSys->SetFillColor(kPink - 4);
  grIaaSTARSys->SetFillColorAlpha(kPink - 4, 0.30);
  grIaaSTARSys->SetLineWidth(0);
  grIaaSTARSys->SetMarkerSize(1.2);
  // grIaaSTARSys->SetFillStyle(1001);
  // grIaaSTARSys->SetFillStyle(0);
  grIaaSTARSys->SetTitle(" ");

  TH1F *hIaaPHENIX = new TH1F("hIaaPHENIX", "hIaaPHENIX", 6, zTPHENIX);
  TH1F *hIaaPHENIXSys = new TH1F("hIaaPHENIXSys", "hIaaPHENIXSys", 6, zTPHENIX);

  double IaaPHENIX[7] = {0.173, 0.336, 1.04, 1.61, 1.79, 1.17};
  double IaaPHENIXBox[7] = {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02};
  double IaaPHENIX_Stat[7] = {0.253, 0.135, 0.185, 0.236, 0.209, 0.288};
  double IaaPHENIX_Sys[7] = {0.267, 0.141, 0.215, 0.381, 0.167, 0.275};

  for (int ibin = 0; ibin < 6; ibin++)
  {
    hIaaPHENIX->SetBinContent(ibin + 1, IaaPHENIX[5 - ibin]);
    hIaaPHENIX->SetBinError(ibin + 1, IaaPHENIX_Stat[5 - ibin]);
    hIaaPHENIXSys->SetBinContent(ibin + 1, IaaPHENIX[5 - ibin]);
    hIaaPHENIXSys->SetBinError(ibin + 1, IaaPHENIX_Sys[5 - ibin]);
  }
  PlotStyle(hIaaPHENIX, 20, 1, kOrange + 2, kOrange + 2, "#it{z}_{T}", "#it{I}_{PYTHIA}, #it{I}_{AA}", false);
  PlotStyle(hIaaPHENIXSys, 20, 1, kOrange - 8, kOrange - 8, "#it{z}_{T}", "#it{I}_{PYTHIA}, #it{I}_{AA}", true);
  // PlotStyle(hPbPb_NLOSyst[0], 21, 1, kRed+1, kRed+1, "#it{z}_{T}", "Ratio", true);

  TCanvas *cPHENIX = new TCanvas("cPHENIX", "cPHENIX", 1 * 800, 1 * 600);
  cPHENIX->Divide(1, 1);
  cPHENIX->cd(1);
  cPHENIX->cd(1)->SetTopMargin(0.015);
  cPHENIX->cd(1)->SetRightMargin(0.02);
  cPHENIX->cd(1)->SetLeftMargin(0.12);
  cPHENIX->cd(1)->SetBottomMargin(0.11);
  hGeneralIaa->SetTitle(" ");
  hGeneralIaa->GetXaxis()->SetTitle("#it{z}_{T}");
  hGeneralIaa->GetYaxis()->SetTitle("#it{I}_{pQCD}, #it{I}_{AA}");
  hGeneralIaa->GetYaxis()->SetTitleSize(0.052);
  hGeneralIaa->GetXaxis()->SetTitleSize(0.05);
  hGeneralIaa->GetYaxis()->SetLabelSize(0.04);
  hGeneralIaa->GetXaxis()->SetLabelSize(0.045);
  hGeneralIaa->GetYaxis()->SetTitleOffset(1.25);
  hGeneralIaa->GetYaxis()->SetLabelOffset(0.02);
  hGeneralIaa->GetXaxis()->SetLabelOffset(0.02);
  hGeneralIaa->GetYaxis()->SetRangeUser(-0.3, 3.8);
  hGeneralIaa->GetYaxis()->SetTickLength(0.02);
  hGeneralIaa->GetXaxis()->SetTickLength(0.02);

  hGeneralIaa->SetLineWidth(0);
  hGeneralIaa->Draw("histsame");
  for (int ibin = 0; ibin < hPbPb_NLOSyst[0]->GetNbinsX() + 1; ++ibin)
  {
    cout << "bin cent:" << hPbPb_NLOSyst[0]->GetBinCenter(ibin) << endl;
    cout << "bin content:" << hPbPb_NLOSyst[0]->GetBinContent(ibin) << endl;
    cout << "bin error:" << hPbPb_NLOSyst[0]->GetBinError(ibin) << endl;
  }

  hIaaPHENIXSys->Draw("samee2");
  hIaaPHENIX->Draw("EP X0 same");
  grIaaSTARSys->Draw("ep2 same");
  grIaaSTAR->Draw("ep same");
  hPbPb_NLOSyst[0]->Draw("E2Psame");
  hPbPb_NLO[0]->Draw("EP X0 same");

  TGraph *linePHENIX = DrawLine(linePHENIX, -0.05, 0.5, 1.05, 0.5);
  linePHENIX->Draw("l");
  TGraph *linePHENIX1 = DrawLine(linePHENIX1, -0.05, 1, 1.05, 1);
  linePHENIX1->Draw("l");

  TLegend *legIaaALICE1_Head = LegStd(legIaaALICE1_Head, 0.280, 0.90, 0.160, 0.94);
  legIaaALICE1_Head->SetTextSize(0.038);
  legIaaALICE1_Head->SetHeader("Pb#font[122]{-}Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  legIaaALICE1_Head->Draw("same");
  TLegend *legIaaALICE1_gamhad = LegStd(legIaaALICE1_gamhad, 0.13, 0.82, 0.46, 0.87);
  legIaaALICE1_gamhad->SetTextSize(0.034);
  legIaaALICE1_gamhad->AddEntry(hPbPb_NLOSyst[0], "ALICE,#color[0]{.}#it{#gamma}^{iso}#font[122]{-}had. 0#font[122]{-}30%", "epf");
  legIaaALICE1_gamhad->Draw("same");
  TLegend *legIaaALICE1_gamhadSelect = LegStd(legIaaALICE1_gamhadSelect, 0.14, 0.77, 0.16, 0.82);
  legIaaALICE1_gamhadSelect->SetTextSize(0.034);
  legIaaALICE1_gamhadSelect->AddEntry((TObject *)0, "18#color[0]{.}<#color[0]{..}#it{p}_{T}^{#it{#gamma}}#color[0]{.}<#color[0]{..}40 GeV/#it{c}#color[0]{.};#color[0]{..}#it{p}_{T}^{h}#color[0]{.}>#color[0]{..}0.5 GeV/#it{c}", "");
  legIaaALICE1_gamhadSelect->Draw("same");
  // TLegend *legIaaALICE1_gamhadMarker = LegStd(legIaaALICE1_gamhadMarker, 0.18, 0.71, 0.50, 0.76);
  // legIaaALICE1_gamhadMarker->SetNColumns(2);
  // legIaaALICE1_gamhadMarker->SetTextSize(0.032);
  // legIaaALICE1_gamhadMarker->AddEntry(hPbPb_NLO[0], "stat. unc.", "epf");
  // legIaaALICE1_gamhadMarker->AddEntry(hPbPb_NLOSyst[0], "syst. unc.", "f");
  // legIaaALICE1_gamhadMarker->Draw("same");

  TLegend *legIaaSTARPHENIX_Head = LegStd(legIaaSTARPHENIX_Head, 0.58, 0.90, 0.78, 0.94);
  legIaaSTARPHENIX_Head->SetTextSize(0.0380);
  legIaaSTARPHENIX_Head->SetHeader("Au#font[122]{-}Au #sqrt{#it{s}_{NN}} = 200 GeV");
  legIaaSTARPHENIX_Head->Draw("same");

  TLegend *legIaaALICE1_STAR = LegStd(legIaaALICE1_STAR, 0.55, 0.82, 0.84, 0.87);
  legIaaALICE1_STAR->SetTextSize(0.034);
  legIaaALICE1_STAR->AddEntry(grIaaSTARSys, "STAR,#color[0]{.}#it{#gamma}^{dir}#font[122]{-}had. 0#font[122]{-}12% ", "epf");
  legIaaALICE1_STAR->Draw("same");
  TLegend *legIaaALICE1_STARSelect = LegStd(legIaaALICE1_STARSelect, 0.560, 0.72, 0.560, 0.82);
  legIaaALICE1_STARSelect->SetTextSize(0.034);
  legIaaALICE1_STARSelect->AddEntry((TObject *)0, "12#color[0]{.}<#color[0]{..}#it{p}_{T}^{#it{#gamma}}#color[0]{.}<#color[0]{..}20 GeV/#it{c}#color[0]{.};#color[0]{..}#it{p}_{T}^{h}#color[0]{.}>#color[0]{..}1.2 GeV/#it{c}", "");
  legIaaALICE1_STARSelect->AddEntry((TObject *)0, "PL B 760 (2016) 689-696", "");
  legIaaALICE1_STARSelect->Draw("same");
  // TLegend *legIaaALICE1_STARMarker = LegStd(legIaaALICE1_STARMarker, 0.560, 0.66, 0.950, 0.71);
  // legIaaALICE1_STARMarker->SetNColumns(2);
  // legIaaALICE1_STARMarker->SetTextSize(0.032);
  // legIaaALICE1_STARMarker->AddEntry(grIaaSTAR, "stat. unc.", "ep");
  // legIaaALICE1_STARMarker->AddEntry(grIaaSTARSys, "syst. unc.", "f");
  // legIaaALICE1_STARMarker->Draw("same");

  TLegend *legIaaALICE1_PHENIX = LegStd(legIaaALICE1_PHENIX, 0.55, 0.63, 0.85, 0.68);
  legIaaALICE1_PHENIX->SetTextSize(0.034);
  legIaaALICE1_PHENIX->AddEntry(hIaaPHENIXSys, "PHENIX,#color[0]{.}#it{#gamma}^{dir}#font[122]{-}had. 0#font[122]{-}40%", "epf");
  legIaaALICE1_PHENIX->Draw("same");
  TLegend *legIaaALICE1_PHENIXSelect = LegStd(legIaaALICE1_PHENIXSelect, 0.56, 0.53, 0.560, 0.63);
  legIaaALICE1_PHENIXSelect->SetTextSize(0.034);
  legIaaALICE1_PHENIXSelect->AddEntry((TObject *)0, "5#color[0]{.}<#color[0]{..}#it{p}_{T}^{#it{#gamma}}#color[0]{.}<#color[0]{..}9 GeV/#it{c}#color[0]{.};#color[0]{..}0.5#color[0]{.}<#color[0]{..}#it{p}_{T}^{h}#color[0]{.}<#color[0]{..}7 GeV/#it{c}", "");
  legIaaALICE1_PHENIXSelect->AddEntry((TObject *)0, "PRL 111, 032301 (2013)", "");
  legIaaALICE1_PHENIXSelect->Draw("same");
  // TLegend *legIaaALICE1_PHENIXMarker = LegStd(legIaaALICE1_PHENIXMarker, 0.56, 0.44, 0.950, 0.49);
  // legIaaALICE1_PHENIXMarker->SetNColumns(2);
  // legIaaALICE1_PHENIXMarker->SetTextSize(0.032);
  // legIaaALICE1_PHENIXMarker->AddEntry(hIaaPHENIX, "stat. unc.", "ep");
  // legIaaALICE1_PHENIXMarker->AddEntry(hIaaPHENIXSys, "syst. unc.", "f");
  // legIaaALICE1_PHENIXMarker->Draw("same");

  // cPHENIX->cd(2);
  //  TLegend *legGH = LegStd(legGH, 0, 0.95, 0.00, 1);
  //  legGH->SetTextSize(0.045);
  //  legGH->SetHeader("#bf{#it{#gamma}#font[122]{-}hadron}");
  //  legGH->Draw("same");

  // TLegend *legALICE = LegStd(legALICE, 0, 0.74, 0.20, 0.99);
  // legALICE->SetTextSize(0.045);
  // legALICE->SetHeader("ALICE");
  // legALICE->AddEntry((TObject *)0, "#bf{#it{#gamma}#font[122]{-}hadron} Pb#font[122]{-}Pb #sqrt{#it{s}_{NN}} = 5.02 TeV", "");
  // legALICE->AddEntry((TObject *)0, "|#Delta #it{#varphi}| < #frac{3}{5}#pi, |#it{#eta}^{#it{#gamma}}| < 0.67", "");
  // legALICE->AddEntry((TObject *)0, "18 < #it{p}_{T}^{ #it{#gamma}} < 40 GeV/#it{c} ; #it{p}_{T}^{h} > 1 GeV/#it{c}", "");
  // legALICE->Draw("same");

  //  //TLatex *latALICEIaa = LatexStdISORatio(latALICEIaa, 0.400, 0.920, 0.045, cenBins[0], cenBins[1], ptMin, ptMax, true);
  //  TLegend *legALICE3 = LegStd(legALICE3, 0.60, 0.50, 0.84, 0.70);
  //  legALICE3->SetTextSize(0.042);
  //  legALICE3->AddEntry(hPbPb_NLO[0], "ALICE, stat. unc.", "ep");
  //  // legALICE3->AddEntry(h3[0], "ALICE, stat. unc.", "ep");
  //  legALICE3->AddEntry(grIaaSTAR, "STAR, stat. unc.", "ep");
  //  legALICE3->AddEntry(hIaaPHENIX, "PHENIX, stat. unc.", "ep");
  //  legALICE3->AddEntry(hPbPb_NLOSyst[0], " syst. unc.", "f");
  //  // legALICE3->AddEntry(hsyst[0], " syst. unc.", "f");
  //  legALICE3->Draw("same");

  //
  //
  // TLegend *legSTAR1 = LegStd(legSTAR1, 0, 0.44, 0.20, 0.69);
  // legSTAR1->SetTextSize(0.045);
  //
  ////legSTAR1->SetHeader("STAR, Phys.Lett.B 760 (2016) 689-696");
  // legSTAR1->AddEntry((TObject *)0, "#bf{0#font[122]{-}12%} Au#font[122]{-}Au, #sqrt{#it{s}_{NN}} = 200 GeV ", "");
  // legSTAR1->AddEntry((TObject *)0, "|#Delta#it{#varphi} #font[122]{-} #it{#pi}| #leq 1.4 ", "");
  // legSTAR1->AddEntry((TObject *)0, "12 < #it{p}_{T}^{ #it{#gamma}} < 20 GeV/#it{c} ; #it{p}_{T}^{h} > 1.2 GeV/#it{c} ", "");
  // legSTAR1->Draw("same");
  //
  // TLegend *legPHENIX = LegStd(legPHENIX, 0, 0.14, 0.20, 0.39);
  // legPHENIX->SetTextSize(0.045);
  ////legPHENIX->SetHeader("PHENIX, PRL 111, 032301 (2013)");
  // legPHENIX->AddEntry((TObject *)0, "#bf{0#font[122]{-}40%} Au#font[122]{-}Au, #sqrt{#it{s}_{NN}} = 200 GeV ", "");
  // legPHENIX->AddEntry((TObject *)0, "|#Delta#it{#varphi} #font[122]{-} #it{#pi}| < #it{#pi}/2, |#it{y}| < 0.35 ", "");
  // legPHENIX->AddEntry((TObject *)0, "5 < #it{p}_{T}^{#it{#gamma}} < 9 GeV/#it{c} ; 0.5 < #it{p}_{T}^{h} < 7 GeV/#it{c} ", "");
  // legPHENIX->Draw("same");
  //
  // cPHENIX->Print(dirPlot + Form("/ppNLOpQCD_DztIaa_ALICESTARPHENIX%s.pdf", sPtAll.Data()));
  cPHENIX->Print(dirPlot + Form("/ppNLOpQCD_Iaa_ALICESTARPHENIX%s.pdf", sPtAll.Data()));

  double xiTCMS[9] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5};
  double zTCMS[9] = {0};

  double DzTCMS0_10[8] = {0.125, 0.379, 0.717, 1.3, 1.68, 2.5, 3.08, 2.21};
  double DzTCMSBox0_10[8] = {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02};
  double DzTCMS_Stat0_10[8] = {0.026, 0.041, 0.057, 0.073, 0.097, 0.171, 0.295, 0.312};
  double DzTCMS_Sys0_10[8] = {0.02, 0.054, 0.11, 0.141, 0.206, 0.306, 0.395, 0.32};

  double IaaCMS0_10[8] = {0.581, 0.637, 0.797, 0.786, 0.813, 1.15, 1.66, 2.0};
  double IaaCMSBox0_10[8] = {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02};
  double IaaCMS_Stat0_10[8] = {0.12, 0.07, 0.052, 0.045, 0.048, 0.08, 0.16, 0.284};
  double IaaCMS_Sys0_10[8] = {0.093, 0.09, 0.099, 0.087, 0.1, 0.136, 0.196, 0.716};

  for (int ibin = 0; ibin < 9; ibin++)
  {
    zTCMS[8 - ibin] = 1. / (TMath::Exp(xiTCMS[ibin]));
  }
  TH1F *hzTCMS = new TH1F("hzTCMS", "hzTCMS", 8, zTCMS);
  TH1F *hzTCMSSys = new TH1F("hzTCMSSys", "hzTCMSSys", 8, zTCMS);

  for (int ibin = 0; ibin < 8; ibin++)
  {
    hzTCMS->SetBinContent(ibin + 1, DzTCMS0_10[7 - ibin]);
    hzTCMS->SetBinError(ibin + 1, DzTCMS_Stat0_10[7 - ibin]);
    hzTCMSSys->SetBinContent(ibin + 1, DzTCMS0_10[7 - ibin]);
    hzTCMSSys->SetBinError(ibin + 1, DzTCMS_Sys0_10[7 - ibin]);
  }
  PlotStyle(hzTCMS, 21, 1.2, kPink - 5, kPink - 5, "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
  PlotStyle(hzTCMSSys, 21, 1.2, kPink + 2, kPink + 2, "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", true);

  TH1F *hIaaCMS = new TH1F("hIaaCMS", "hIaaCMS", 8, zTCMS);
  TH1F *hIaaCMSSys = new TH1F("hIaaCMSSys", "hIaaCMSSys", 8, zTCMS);
  for (int ibin = 0; ibin < 8; ibin++)
  {
    hIaaCMS->SetBinContent(ibin + 1, IaaCMS0_10[7 - ibin]);
    hIaaCMS->SetBinError(ibin + 1, IaaCMS_Stat0_10[7 - ibin]);
    cout << IaaCMS0_10[7 - ibin] << endl;
    hIaaCMSSys->SetBinContent(ibin + 1, IaaCMS0_10[7 - ibin]);
    hIaaCMSSys->SetBinError(ibin + 1, IaaCMS_Sys0_10[7 - ibin]);
  }
  PlotStyle(hIaaCMS, 24, 1.1, kOrange + 4, kOrange + 2, "#it{z}_{T}", "#it{I}_{pQCD}, #it{I}_{AA}", false);
  PlotStyle(hIaaCMSSys, 24, 1.1, kOrange + 4, kOrange - 3, "#it{z}_{T}", "#it{I}_{pQCD}, #it{I}_{AA}", true);

  double xiTCMS_Zhad[11] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
  double zTCMS_Zhad[11] = {0};

  double DzTCMS0_10_Zhad[10] = {0.0196, 0.0598, 0.157, 0.4, 0.875, 1.673, 3.81, 2.78, 1.13, 0.484};
  double DzTCMSBox0_10_Zhad[10] = {};
  double DzTCMS_Stat0_10_Zhad[10] = {0.0043, 0.0083, 0.015, 0.029, 0.08, 0.711, 0.47, 0.42, 0.25, 0.136};
  double DzTCMS_Sys0_10_Zhad[10] = {0.0028, 0.0057, 0.012, 0.026, 0.053, 0.097, 0.25, 0.71, 0.12, 0.064};

  double IaaCMS0_10_Zhad[10] = {0.612, 0.397, 0.399, 0.507, 0.682, 0.864, 1.441, 1.419, 1.519, 2.32};
  double IaaCMSBox0_10_Zhad[10] = {};
  double IaaCMS_Stat0_10_Zhad[10] = {0.139, 0.056, 0.04, 0.037, 0.063, 0.125, 0.179, 0.214, 0.347, 0.66};
  double IaaCMS_Sys0_10_Zhad[10] = {0.054, 0.026, 0.071, 0.029, 0.04, 0.05, 0.082, 0.088, 0.099, 0.2};

  for (int ibin = 0; ibin < 11; ibin++)
  {
    zTCMS_Zhad[10 - ibin] = 1. / (TMath::Exp(xiTCMS_Zhad[ibin]));
    cout << zTCMS_Zhad[10 - ibin] << endl;
  }
  TH1F *hzTCMS_Zhad = new TH1F("hzTCMS_Zhad", "hzTCMS_Zhad", 10, zTCMS_Zhad);
  TH1F *hzTCMSSys_Zhad = new TH1F("hzTCMSSys_Zhad", "hzTCMSSys_Zhad", 10, zTCMS_Zhad);

  for (int ibin = 0; ibin < 10; ibin++)
  {
    hzTCMS_Zhad->SetBinContent(ibin + 1, DzTCMS0_10_Zhad[9 - ibin]);
    hzTCMS_Zhad->SetBinError(ibin + 1, DzTCMS_Stat0_10_Zhad[9 - ibin]);
    hzTCMSSys_Zhad->SetBinContent(ibin + 1, DzTCMS0_10_Zhad[9 - ibin]);
    hzTCMSSys_Zhad->SetBinError(ibin + 1, DzTCMS_Sys0_10_Zhad[9 - ibin]);
  }
  PlotStyle(hzTCMS_Zhad, 28, 1.5, kTeal + 3, kTeal + 3, "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
  PlotStyle(hzTCMSSys_Zhad, 28, 1.5, kTeal + 3, kTeal + 3, "#it{z}_{T}", "1 / #it{N}^{#it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", true);

  TH1F *hIaaCMS_Zhad = new TH1F("hIaaCMS_Zhad", "hIaaCMS_Zhad", 10, zTCMS_Zhad);
  TH1F *hIaaCMSSys_Zhad = new TH1F("hIaaCMSSys_Zhad", "hIaaCMSSys_Zhad", 10, zTCMS_Zhad);
  for (int ibin = 0; ibin < 10; ibin++)
  {
    hIaaCMS_Zhad->SetBinContent(ibin + 1, IaaCMS0_10_Zhad[9 - ibin]);
    hIaaCMS_Zhad->SetBinError(ibin + 1, IaaCMS_Stat0_10_Zhad[9 - ibin]);
    hIaaCMSSys_Zhad->SetBinContent(ibin + 1, IaaCMS0_10_Zhad[9 - ibin]);
    hIaaCMSSys_Zhad->SetBinError(ibin + 1, IaaCMS_Sys0_10_Zhad[9 - ibin]);
  }
  PlotStyle(hIaaCMS_Zhad, kFullCross, 1.5, kTeal - 6, kTeal - 6, "#it{z}_{T}", "#it{I}_{pQCD}, #it{I}_{AA}", false);
  // PlotStyle(hGeneralCMSIaa, 71, 1, kBlack, "#it{z}_{T}", "#it{I}_{PYTHIA}, #it{I}_{AA}", false);
  PlotStyle(hIaaCMSSys_Zhad, kFullCross, 1.5, kTeal - 5, kTeal - 5, "#it{z}_{T}", "#it{I}_{pQCD}, #it{I}_{AA}", true);

  TCanvas *cCMS_Zhad = new TCanvas("cCMS_Zhad", "cCMS_Zhad", 1 * 800, 1 * 600);
  cCMS_Zhad->Divide(1, 1);
  //  cCMS_Zhad->cd(1);
  //  cCMS_Zhad->cd(1)->SetTopMargin(0.015);
  //  cCMS_Zhad->cd(1)->SetRightMargin(0.02);
  //  cCMS_Zhad->cd(1)->SetLeftMargin(0.12);
  //  cCMS_Zhad->cd(1)->SetBottomMargin(0.11);
  //  gPad->SetLogy();
  //  hGeneralCMS->SetTitle(" ");
  //  hGeneralCMS->GetYaxis()->SetTitleSize(0.045);
  //  hGeneralCMS->GetXaxis()->SetTitleSize(0.05);
  //  hGeneralCMS->GetYaxis()->SetLabelSize(0.04);
  //  hGeneralCMS->GetXaxis()->SetLabelSize(0.045);
  //  hGeneralCMS->GetYaxis()->SetLabelOffset(0.02);
  //  hGeneralCMS->GetXaxis()->SetLabelOffset(0.02);
  //
  //  hGeneralCMS->GetYaxis()->SetRangeUser(3e-3, 10);
  //  hGeneralCMS->GetXaxis()->SetRangeUser(-0.02, 1.08);
  //  // hGeneral->Draw("hist");
  //  hGeneralCMS->SetLineWidth(0);
  //  hGeneralCMS->Draw("hist");
  //  hzTCMSSys_Zhad->SetFillStyle(0);
  //  hzTCMSSys_Zhad->SetLineWidth(2);
  //  hzTCMSSys_Zhad->Draw("samee2");
  //  hzTCMS_Zhad->Draw("EP X0 same");
  //  hzTCMSSys->SetFillStyle(0);
  //  hzTCMSSys->SetLineWidth(2);
  //  hzTCMSSys->Draw("samee2");
  //  hzTCMS->Draw("EP X0 same");
  //  hSystZt[0]->SetMarkerColor(kRed + 1);
  //  hSystZt[0]->SetFillColor(kRed + 1);
  //  hSystZt[0]->SetLineColor(kRed + 1);
  //  hSystZt[0]->SetLineWidth(2);
  //  hZtCent[0]->SetMarkerColor(kRed + 1);
  //  hZtCent[0]->SetLineColor(kRed + 1);
  //  hZtCent[0]->SetMarkerSize(1.2);
  //  hSystZt[0]->SetFillStyle(0);
  //  hSystZt[0]->Draw("samee2");
  //  hZtCent[0]->Draw("EP X0 same");
  //  //TLatex *latALICEcms_Zhad = LatexStdISORatio(latALICEcms_Zhad, 0.380, 0.92, 0.045, cenBins[0], cenBins[1], ptMin, ptMax, true);
  //  //TLegend *legALICE2_Zhad = LegStd(legALICE2_Zhad, 0.140, 0.16, 0.40, 0.40);
  //  TLegend *legALICE2_Zhad = LegStd(legALICE2_Zhad, 0.45, 0.7, 0.9, 0.95);
  //  legALICE2_Zhad->SetTextSize(0.045);
  //  legALICE2_Zhad->AddEntry(hZtCent[0], "ALICE, #bf{#it{#gamma}#font[122]{-}hadron} stat. unc. ", "ep");
  //  legALICE2_Zhad->AddEntry(hzTCMS, "CMS, #bf{#it{#gamma}#font[122]{-}jet}, stat. unc. ", "ep");
  //  legALICE2_Zhad->AddEntry(hzTCMS_Zhad, "CMS, #bf{#it{Z}#font[122]{-}hadron}, stat. unc. ", "ep");
  //  legALICE2_Zhad->AddEntry(hSystZt[0], " syst. unc. ", "f");
  //  legALICE2_Zhad->Draw("same");

  cCMS_Zhad->cd(1);
  cCMS_Zhad->cd(1)->SetTopMargin(0.015);
  cCMS_Zhad->cd(1)->SetRightMargin(0.02);
  cCMS_Zhad->cd(1)->SetLeftMargin(0.12);
  cCMS_Zhad->cd(1)->SetBottomMargin(0.11);
  // hIaaCMSSys_Zhad->SetFillStyle(0);
  // hIaaCMSSys_Zhad->SetLineWidth(2);
  // hGeneralIaa->SetTitle(" ");
  hGeneralIaa->GetYaxis()->SetRangeUser(0, 3.1);
  ////hGeneralIaa->GetXaxis()->SetRangeUser(0.0, 1.1);
  // hGeneralIaa->GetYaxis()->SetTitleSize(0.052);
  // hGeneralIaa->GetYaxis()->SetTitleOffset(1.25);
  // hGeneralIaa->GetXaxis()->SetTitleSize(0.045);
  // hGeneralIaa->GetYaxis()->SetLabelSize(0.04);
  // hGeneralIaa->GetXaxis()->SetLabelSize(0.04);
  // hGeneralIaa->SetLineWidth(0);
  hGeneralIaa->Draw("hist");
  // hIaaCMSSys->SetLineWidth(2);
  // hIaaCMSSys->SetFillStyle(0);
  hIaaCMSSys->Draw("samee2");
  hIaaCMS->Draw("EP X0 same");
  hIaaCMSSys_Zhad->Draw("samee2");
  hIaaCMS_Zhad->Draw("EP X0 same");
  // hPbPb_NLO[0]->SetLineColor(kRed + 1);
  // hPbPb_NLO[0]->SetMarkerColor(kRed + 1);
  // hPbPb_NLO[0]->SetMarkerSize(1.3);
  // hPbPb_NLO[0]->SetLineWidth(2);
  // hPbPb_NLOSyst[0]->SetLineWidth(2);
  // hPbPb_NLOSyst[0]->SetFillColor(kRed + 1);
  // hPbPb_NLOSyst[0]->SetLineColor(kRed + 1);
  // hPbPb_NLOSyst[0]->SetMarkerSize(1.3);
  hPbPb_NLOSyst[0]->Draw("E2Psame ");
  hPbPb_NLO[0]->Draw("EPX0same");
  // h3[0]->SetLineColor(kRed + 1);
  // h3[0]->SetMarkerColor(kRed + 1);
  // h3[0]->SetMarkerSize(1.3);
  // hsyst[0]->SetLineWidth(2);
  // hsyst[0]->SetFillColor(kRed + 1);
  // hsyst[0]->SetLineColor(kRed + 1);
  // hsyst[0]->SetFillStyle(0);
  // hsyst[0]->SetMarkerSize(1.3);
  // hsyst[0]->Draw("samee2");
  // h3[0]->Draw("EP X0 same");
  TLegend *legALICE2_Head = LegStd(legALICE2_Head, 0.20, 0.90, 0.30, 0.94);
  legALICE2_Head->SetTextSize(0.040);
  legALICE2_Head->SetHeader("Pb#font[122]{-}Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  legALICE2_Head->Draw("same");
  TLegend *legALICE2_Zhad = LegStd(legALICE2_Zhad, 0.55, 0.90, 0.88, 0.94);
  legALICE2_Zhad->SetTextSize(0.034);
  legALICE2_Zhad->AddEntry(hPbPb_NLOSyst[0], "ALICE,#color[0]{.}#it{#gamma}^{iso}#font[122]{-}had. 0#font[122]{-}30%", "epf");
  legALICE2_Zhad->Draw("same");
  TLegend *legALICE2_ZhadSelect = LegStd(legALICE2_ZhadSelect, 0.56, 0.85, 0.580, 0.90);
  legALICE2_ZhadSelect->SetTextSize(0.034);
  legALICE2_ZhadSelect->AddEntry((TObject *)0, "18 <#color[0]{..}#it{p}_{T}^{#it{#gamma}} < 40 GeV/#it{c}#color[0]{..};#color[0]{..}#it{p}_{T}^{h} > 0.5 GeV/#it{c}", "");
  legALICE2_ZhadSelect->Draw("same");
  // TLegend *legALICE2_ZhadMarker = LegStd(legALICE2_ZhadMarker, 0.20, 0.73, 0.50, 0.78);
  // legALICE2_ZhadMarker->SetNColumns(2);
  // legALICE2_ZhadMarker->SetTextSize(0.032);
  // legALICE2_ZhadMarker->AddEntry(hPbPb_NLO[0], "stat. unc.", "ep");
  // legALICE2_ZhadMarker->AddEntry(hPbPb_NLOSyst[0], "syst. unc.", "f");
  // legALICE2_ZhadMarker->Draw("same");

  // TLegend *legALICE2_CMSjet = LegStd(legALICE2_CMSjet, 0.60, 0.74, 0.62, 0.98);
  TLegend *legALICE2_CMSjet = LegStd(legALICE2_CMSjet, 0.55, 0.775, 0.88, 0.815);
  legALICE2_CMSjet->SetTextSize(0.034);
  legALICE2_CMSjet->AddEntry(hIaaCMSSys, "CMS,#color[0]{.}#it{#gamma}^{iso}#font[122]{-}jet 0#font[122]{-}10%", "epf");
  legALICE2_CMSjet->Draw("same");
  TLegend *legALICE2_CMSjetSelect = LegStd(legALICE2_CMSjetSelect, 0.56, 0.675, 0.580, 0.775);
  legALICE2_CMSjetSelect->SetTextSize(0.034);
  legALICE2_CMSjetSelect->AddEntry((TObject *)0, "#it{p}_{T}^{#it{#gamma}} > 60 GeV/#it{c}#color[0]{..};#color[0]{..}#it{p}_{T}^{h} > 1 GeV/#it{c}", "");
  legALICE2_CMSjetSelect->AddEntry((TObject *)0, "PRL 121 (2018) 24, 242301", "");
  legALICE2_CMSjetSelect->Draw("same");
  // TLegend *legALICE2_CMSjetMarker = LegStd(legALICE2_CMSjetMarker, 0.62, 0.67, 0.92, 0.72);
  // legALICE2_CMSjetMarker->SetNColumns(2);
  // legALICE2_CMSjetMarker->SetTextSize(0.032);
  // legALICE2_CMSjetMarker->AddEntry(hIaaCMS, "stat. unc.", "ep");
  // legALICE2_CMSjetMarker->AddEntry(hIaaCMSSys, "syst. unc.", "f");
  // legALICE2_CMSjetMarker->Draw("same");

  TLegend *legALICE2_CMSZhad = LegStd(legALICE2_CMSZhad, 0.55, 0.585, 0.88, 0.625);
  legALICE2_CMSZhad->SetTextSize(0.034);
  legALICE2_CMSZhad->AddEntry(hIaaCMSSys_Zhad, "CMS,#color[0]{.}Z^{0}#font[122]{-}had. 0#font[122]{-}30%", "epf");
  legALICE2_CMSZhad->Draw("same");
  TLegend *legALICE2_CMSZhadSelect = LegStd(legALICE2_CMSZhadSelect, 0.56, 0.48, 0.580, 0.58);
  legALICE2_CMSZhadSelect->SetTextSize(0.034);
  legALICE2_CMSZhadSelect->AddEntry((TObject *)0, "#it{p}_{T}^{Z} > 30 GeV/#it{c}#color[0]{..};#color[0]{..}#it{p}_{T}^{h} > 1 GeV/#it{c}", "");
  legALICE2_CMSZhadSelect->AddEntry((TObject *)0, "PRL 128 (2022) 12, 122301", "");
  legALICE2_CMSZhadSelect->Draw("same");
  // legALICE2_CMSZhad->AddEntry((TObject *)0, "PRL 128 (2022) 12, 122301", "");
  // TLegend *legALICE2_CMSZhadMarker = LegStd(legALICE2_CMSZhadMarker, 0.62, 0.39, 0.92, 0.44);
  // legALICE2_CMSZhadMarker->SetNColumns(2);
  // legALICE2_CMSZhadMarker->SetTextSize(0.032);
  // legALICE2_CMSZhadMarker->AddEntry(hIaaCMS_Zhad, "stat. unc.", "ep");
  // legALICE2_CMSZhadMarker->AddEntry(hIaaCMSSys_Zhad, "syst. unc.", "f");
  // legALICE2_CMSZhadMarker->Draw("same");

  TGraph *lineCMS = DrawLine(lineCMS, -0.05, 0.5, 1.1, 0.5);
  lineCMS->Draw("l");
  TGraph *lineCMS1 = DrawLine(lineCMS1, -0.05, 1, 1.1, 1);
  lineCMS1->Draw("l");
  cCMS_Zhad->Print(dirPlot + Form("/ppNLOpQCD_Iaa_ZhadALICECMS%s.pdf", sPtAll.Data()));

  ////////////////
  /// ICP///////
  ///////////////
  ///

  // Recover data Icp calculated in Systematic calculation stage
  //
  TFile *fIcp = TFile::Open("Systematics_checkCode0_30/fIcpWithSystMixed0.40-1.00Mirror_Pt18_40.root");
  TH1F *hIcpSta[nCen - 1];
  TH1F *hIcpSys[nCen - 1];

  for (Int_t iCen = 0; iCen < nCen - 1; iCen++)
  {
    hIcpSta[iCen] = (TH1F *)fIcp->Get(Form("hIcp_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    hIcpSys[iCen] = (TH1F *)fIcp->Get(Form("hIcpSyst_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    // printf("Icp: icen %d, sta %p sys %p\n", iCen, hIcpSta[iCen], hIcpSys[iCen]);
  }

  // Calculate Icp in theory for NLO pQCD
  //
  for (Int_t iCen = 0; iCen < nCen - 1; iCen++)
  {
    Int_t nPoints = grIcpNLOmedian[iCen]->GetN();
    for (Int_t ibin = 0; ibin < nPoints; ibin++)
    {
      Float_t zT = grIcpNLOmedian[iCen]->GetPointX(ibin);
      Float_t dzTCen = grIcpNLOmedian[iCen]->GetPointY(ibin);
      Float_t dzTPer = grIcpNLOmedian[nCen - 1]->GetPointY(ibin);
      Float_t iCP = 0;
      if (dzTPer > 0)
        iCP = dzTCen / dzTPer;

      Float_t errCenHigh = grIcpNLOmedian[iCen]->GetErrorYhigh(ibin);
      Float_t errCenLow = grIcpNLOmedian[iCen]->GetErrorYlow(ibin);
      Float_t errPerHigh = grIcpNLOmedian[nCen - 1]->GetErrorYhigh(ibin);
      Float_t errPerLow = grIcpNLOmedian[nCen - 1]->GetErrorYlow(ibin);
      Float_t errIcpHigh = errCenHigh * errCenHigh / (dzTPer * dzTPer) + errPerHigh * errPerHigh * dzTCen * dzTCen / (dzTPer * dzTPer * dzTPer * dzTPer);
      Float_t errIcpLow = errCenHigh * errCenLow / (dzTPer * dzTPer) + errPerLow * errPerLow * dzTCen * dzTCen / (dzTPer * dzTPer * dzTPer * dzTPer);
      errIcpHigh = TMath::Sqrt(errIcpHigh);
      errIcpLow = TMath::Sqrt(errIcpLow);

      printf("iCen %d Point %d, zT %1.3f, Dzt Cen %1.2e, Per %1.2e, Icp %1.2f --- Err Cen high %1.2e low %1.2e Per high %1.2e low %1.2e Icp high %1.2e low %1.2e\n",
             iCen, ibin,
             zT, dzTCen, dzTPer, iCP,
             errCenHigh, errCenLow,
             errPerHigh, errPerLow,
             errIcpHigh, errIcpLow);
      printf("\t X Err Cen high %1.2e low %1.2e \n",
             grIcpNLOmedian[iCen]->GetErrorXhigh(ibin), grIcpNLOmedian[iCen]->GetErrorXlow(ibin));

      grIcpNLOmedian[iCen]->SetPoint(ibin, zT, iCP);
      grIcpNLOmedian[iCen]->SetPointError(ibin, 0, 0, errIcpHigh, errIcpLow);
    }
    grIcpNLOmedian[iCen]->SetLineStyle(0);
    grIcpNLOmedian[iCen]->SetLineColorAlpha(0, 0);
  }

  gStyle->SetPadTopMargin(0.015);
  gStyle->SetPadRightMargin(0.015);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);

  TCanvas *cICP = new TCanvas(Form("cICP"), Form("cICP"), 1 * 800, 1 * 600);
  cICP->Divide(1, 1);

  gPad->SetTickx();
  gPad->SetTicky();

  double rangeIcpMin[2] = {0.1, 0.1};
  double rangeIcpMax[2] = {0.6, 0.6};

  TH1F *hGeneralIcp = new TH1F("hGeneralIcp", "hGeneral Icp", 100, 0, 2);
  PlotStyle(hGeneralIcp, 20, 1, kWhite, kWhite, "#it{z}_{T}", "#it{I}_{CP}", false);
  hGeneralIcp->SetTitle(" ");
  hGeneralIcp->GetXaxis()->SetRangeUser(0.05, 0.65);
  hGeneralIcp->GetYaxis()->SetRangeUser(-0.1, 2);
  hGeneralIcp->GetYaxis()->SetTitleSize(0.05);
  hGeneralIcp->GetXaxis()->SetTitleSize(0.05);
  hGeneralIcp->GetYaxis()->SetLabelSize(0.05);
  hGeneralIcp->GetXaxis()->SetLabelSize(0.05);
  hGeneralIcp->GetXaxis()->SetTitleOffset(1.0);
  hGeneralIcp->GetYaxis()->SetTitleOffset(1.);

  hGeneralIcp->Draw("");

  TGraph *lineIcp1 = DrawLine(line, 0, 1, 1.2, 1);
  lineIcp1->Draw("same");
  TGraph *lineIcp05 = DrawLine(line1, 0, 0.5, 1.2, 0.5);
  lineIcp05->Draw("same");

  hIcpSys[0]->SetFillColorAlpha(kOrange + 8, 0.30);
  hIcpSys[0]->SetMarkerColor(kRed - 4);
  hIcpSta[0]->SetMarkerColor(kRed - 4);
  hIcpSta[0]->SetLineColor(kRed - 4);
  hIcpSta[1]->SetMarkerColor(kViolet + 7);
  hIcpSta[1]->SetLineColor(kViolet + 7);
  hIcpSys[1]->SetMarkerColor(kViolet + 7);
  hIcpSys[1]->SetFillColorAlpha(kViolet + 6, 0.30);

  TLegend *legIcp[2];
  legIcp[0] = LegStd(legIcp[0], 0.12, 0.80, 0.30, 0.98);
  legIcp[1] = LegStd(legIcp[1], 0.40, 0.80, 0.58, 0.98);

  grIcpNLOmedian[0]->SetFillColorAlpha(kGray, 0.60);

  for (int iCen = 0; iCen < nCen - 1; iCen++)
  {
    hIcpSys[iCen]->GetXaxis()->SetRangeUser(rangeIcpMin[iCen], rangeIcpMax[iCen]);
    hIcpSta[iCen]->GetXaxis()->SetRangeUser(rangeIcpMin[iCen], rangeIcpMax[iCen]);

    grIcpNLOmedian[iCen]->Draw("same E3");

    hIcpSys[iCen]->Draw("same E2");
    hIcpSta[iCen]->Draw("same EX0");

    legIcp[iCen]->SetHeader(Form("%d#font[122]{-}%d%%  / %d#font[122]{-}%d%%:",
                                 cenBins[iCen], cenBins[iCen + 1], cenBins[nCen - 1], cenBins[nCen]));
    //    legIcp[iCen]->AddEntry(hIcpSta[iCen], "Data stat. unc.", "ep");
    legIcp[iCen]->AddEntry(hIcpSys[iCen], Form("Data"), "epf");
    legIcp[iCen]->AddEntry(grIcpNLOmedian[iCen], Form("NLO pQCD#color[0]{.}+#color[0]{.}#Delta#it{E}_{loss}"), "f");
    legIcp[iCen]->Draw("same");
  }

  TLatex *ALICEtexIcp2 = new TLatex();
  ALICEtexIcp2->SetTextFont(42);
  ALICEtexIcp2->SetTextSize(0.04);
  ALICEtexIcp2->SetNDC();
  // ALICEtexIcp2->DrawLatex(0.56, 0.86, Form("#it{This Thesis}"));
  ALICEtexIcp2->DrawLatex(0.68, 0.94, Form("ALICE"));
  ALICEtexIcp2->DrawLatex(0.68, 0.94 - 0.055, Form("Pb#font[122]{-}Pb,#color[0]{..}#sqrt{#it{s}_{NN}} = 5.02 TeV "));
  ALICEtexIcp2->DrawLatex(0.68, 0.94 - 2 * 0.06, Form("|#Delta#it{#varphi}| > #frac{3}{5}#it{#pi}, |#it{#eta}^{#it{#gamma}}| < 0.67 "));
  ALICEtexIcp2->DrawLatex(0.68, 0.94 - 3 * 0.06, Form("%2.0f < #it{p}_{T}^{#it{#gamma}} < %2.0f GeV/#it{c} ", ptMin, ptMax));
  ALICEtexIcp2->DrawLatex(0.68, 0.94 - 4 * 0.06, Form("#it{p}_{T}^{h} > 0.5 GeV/#it{c} "));

  cICP->Print(dirPlot + Form("/Icp_Data_NLOpQCD%s.pdf", sPtAll.Data()));
}

TH1F *shiftBinXinHisto(TH1F *h1, TH1F *h1Shifted, int iCen)
{
  // Assuming you already have:
  // TH1F* hRed;   // e.g. central collisions
  // TH1F* hBlue;  // e.g. peripheral collisions

  // Create a clone of hRed to shift its bin content
  // h1Shifted = (TH1F *)h1->Clone("h1Shifted");
  // h1Shifted->Reset(); // Clear contents before refilling
  int a = 0;
  const int nBins = h1->GetNbinsX();
  for (int i = 1; i <= nBins; ++i)
  {
    double x = h1->GetBinCenter(i);
    double y = h1->GetBinContent(i);
    double err = h1->GetBinError(i);

    // Shift bin center by +0.01 units and set the new value
    if (iCen == 0)
      a = -1;
    else if (iCen == 1)
      a = 1;
    else if (iCen == 2)
      a = 0;
    cout << "a value_______" << a << endl;
    int newBin = h1Shifted->FindBin(x + (a * 0.01));
    cout << (a * 0.01) << endl;
    h1Shifted->SetBinContent(newBin, y);
    h1Shifted->SetBinError(newBin, err);
  }
  return h1Shifted;
}
