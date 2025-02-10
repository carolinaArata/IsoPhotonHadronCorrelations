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
// Int_t nCen = 4;
// Int_t cenBins[] = {0, 10, 30, 50, 90};
//  Int_t cenBins[] = { 50, 90};
Int_t kMarkCen[] = {21, 20, 71, 25};
// Int_t kColorMark[] = {kCyan + 2, kAzure - 3, kViolet + 6, kCyan - 2};
Int_t kColorMark[] = {kAzure + 2, kOrange + 8, kViolet + 7, kCyan - 2};
Int_t kColorMarkFill[] = {kAzure + 5, kOrange + 7, kViolet + 6, kCyan - 2};
Int_t nAssoc = 7;

int nZtBinThin = 9;
double assocZtThinner[] = {0, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80, 1.00, 1.05};


void PlotZtCentCopyNew(float ptMin = 18, float ptMax = 40, bool Mirror = true, TString sMixed = "Mixed", TString shshBkg = "0.40-1.00", TString dirPlot = "~/work/histogram/zTfunct_PbPb_ppZtCheckCode", bool b0_30 = true)
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

  Int_t nCen;
  std::vector<Int_t> cenBins;
  TString dirSyst;
  if (b0_30)
  {
    nCen = 3;
    // Int_t cenBins[] = {0, 10, 30, 50, 90};
    cenBins.push_back(0);
    cenBins.push_back(30);
    cenBins.push_back(50);
    cenBins.push_back(90);
    dirSyst = "~/work/histogram/Systematics_checkCode0_30";
  }
  else if (!b0_30)
  {
    nCen = 4;
    cenBins.push_back(0);
    cenBins.push_back(10);
    cenBins.push_back(30);
    cenBins.push_back(50);
    cenBins.push_back(90);
    dirSyst = "~/work/histogram/Systematics_checkCode";
  }

  TString shshString[2] = {"0.10-0.30", shshBkg};
  TString sPtAll = Form("_Pt%2.0f_%2.0f", ptMin, ptMax);

  // Getter zt distributions
  TFile *fPlot[nCen];
  TH1F *hZtCent[nCen];    // zt data
  TH1F *hZt_MC_Gen[nCen]; // MC Gen pp
  TH1F *hZt_MC_Rec[nCen]; // MC Rec pp
  TH1F *h3[nCen];
  // pQCD NLO
  TGraphAsymmErrors *grIaaNLOmedian[nCen];
  TH1F *grDztNLOmedianpp;
  TGraphAsymmErrors *grDztNLOmedian[nCen];
  TGraphAsymmErrors *grIcpNLOmedian[nCen];
  // COLBT
  TH1F *histIaaCOLBTmedian[nCen];
  TH1F *histDztPbPbCOLBTmedian[nCen];
  TH1F *histIaaCOLBTmedianSyst[2];
  TH1F *histDztPbPbCOLBTmedianSyst[2];
  TH1F *histRatioPYTHIACoLBt[2];
  TCanvas *cRatioPYTHIACoLBt[2];
  TCanvas *cPYTHIA_CoLBT[2];

  TFile *fileNLO = new TFile(" ~/work/histogram/IsoPhotonHadronCorrelations/fileNLO.root ");

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    cout << "Getter zt distributions: " << sCent << endl;
    fPlot[iCen] = new TFile(Form("~/work/histogram/FromScratch/checkCode/fPlot%s%s%s.root", shshString[1].Data(), sCent.Data(), sPtAll.Data()));

    hZt_MC_Gen[iCen] = (TH1F *)fPlot[iCen]->Get(Form("hZtMCGenIso1Photon%s%s", sCent.Data(), sPtAll.Data()));
    hZt_MC_Rec[iCen] = (TH1F *)fPlot[iCen]->Get(Form("hZtMCRecIso1Photon%s%s", sCent.Data(), sPtAll.Data()));
    hZtCent[iCen] = (TH1F *)fPlot[iCen]->Get(Form("hZtEffCorrIso1Photon%s%s", sCent.Data(), sPtAll.Data()));
    grIaaNLOmedian[iCen] = (TGraphAsymmErrors *)fileNLO->Get(Form("grIaaNLOmedian%s", sCent.Data()));
    grIaaNLOmedian[iCen]->SetLineWidth(8);
    grIaaNLOmedian[iCen]->SetLineColor(kPink + 4);
    grIaaNLOmedian[iCen]->SetFillColorAlpha(kMagenta - 7, 0.75);
    grIaaNLOmedian[iCen]->SetFillStyle(3008);

    grDztNLOmedian[iCen] = (TGraphAsymmErrors *)fileNLO->Get(Form("grDztNLOmedian%s", sCent.Data()));
    grDztNLOmedian[iCen]->SetLineWidth(8);
    grDztNLOmedian[iCen]->SetLineColor(kPink + 4);
    grDztNLOmedian[iCen]->SetFillColorAlpha(kMagenta - 7, 0.75);
    grDztNLOmedian[iCen]->SetFillStyle(3008);

    grIcpNLOmedian[iCen] = (TGraphAsymmErrors *)fileNLO->Get(Form("grDztNLOmedian%s", sCent.Data()));
    grIcpNLOmedian[iCen]->SetLineWidth(8);
    grIcpNLOmedian[iCen]->SetFillStyle(3008);

    PlotStyle(hZt_MC_Gen[iCen], 72, 1, kBlack, kBlack, "#it{z}_{T}", "1 / #it{N}^{ #it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d #it{z}_{T}", false);
    PlotStyle(hZt_MC_Rec[iCen], 21, 1, kOrange + 7, kOrange + 7, "#it{z}_{T}", "1 / #it{N}^{ #it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d #it{z}_{T}", false);
    PlotStyle(hZtCent[iCen], kMarkCen[iCen], 1, kColorMark[iCen], kColorMarkFill[iCen], "#it{z}_{T}", "1 / #it{N}^{ #it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d #it{z}_{T}", false);
    h3[iCen] = (TH1F *)hZtCent[iCen]->Clone(Form("h3%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    h3[iCen]->Divide(hZt_MC_Gen[iCen]);
  }
  gSystem->Exec(Form("mkdir %s", dirPlot.Data()));
  cout << "Getter zt systematics : " << endl;

  for (int iCen = 0; iCen < 2; iCen++)
  {
    histIaaCOLBTmedian[iCen] = (TH1F *)fileNLO->Get(Form("histIaaCOLBTmedian_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    histIaaCOLBTmedian[iCen]->SetLineWidth(8);
    histIaaCOLBTmedian[iCen]->SetLineColor(kTeal + 3);
    histIaaCOLBTmedian[iCen]->SetFillColorAlpha(kGreen + 1, 0.40);

    histIaaCOLBTmedianSyst[iCen] = (TH1F *)histIaaCOLBTmedian[iCen]->Clone(Form("histIaaCOLBTmedianSyst_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]));

    histDztPbPbCOLBTmedian[iCen] = (TH1F *)fileNLO->Get(Form("histDztPbPbCOLBTmedian_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    histDztPbPbCOLBTmedian[iCen]->SetLineWidth(8);
    histDztPbPbCOLBTmedian[iCen]->SetLineColor(kTeal + 3);
    histDztPbPbCOLBTmedian[iCen]->SetFillColorAlpha(kGreen + 1, 0.40);

    histDztPbPbCOLBTmedianSyst[iCen] = (TH1F *)histDztPbPbCOLBTmedian[iCen]->Clone(Form("histDztPbPbCOLBTmedianSyst_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
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
  grDztNLOmedianpp = (TH1F *)fileNLO->Get(Form("grDztNLOmedian_pp"));
  grDztNLOmedianpp->SetLineWidth(8);
  grDztNLOmedianpp->SetLineColor(kRed - 4);
  grDztNLOmedianpp->SetLineStyle(10);
  // grDztNLOmedianpp->SetFillColorAlpha(kGray + 1, 0.40);

  TFile *fSystFile = new TFile(Form("%s/fAllSystFile%s%s%s%s.root", dirSyst.Data(), sMixed.Data(), shshBkg.Data(), sMirror.Data(), sPtAll.Data()));
  cout << dirSyst << endl;
  TH1F *hSystZt[nCen];
  TH1F *hsyst[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    hSystZt[iCen] = (TH1F *)fSystFile->Get(Form("hsystCen%d%d", cenBins[iCen], cenBins[iCen + 1]));
    PlotStyle(hSystZt[iCen], kMarkCen[iCen], 1, kColorMark[iCen], kColorMarkFill[iCen], "#it{z}_{T}", "1 / #it{N}^{ #it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d #it{z}_{T}", true);
    hsyst[iCen] = (TH1F *)hSystZt[iCen]->Clone(Form("hsyst%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    hsyst[iCen]->Divide(hZt_MC_Gen[iCen]);
    for (int ibin = 0; ibin < nAssoc; ibin++)
    {
      cout << cenBins[iCen] << "-" << cenBins[iCen + 1] << endl;
      cout << ibin << "-" << ibin + 1 << " systematics: " << hSystZt[iCen]->GetBinContent(ibin + 1) << " pm " << hSystZt[iCen]->GetBinError(ibin + 1) << endl;
    }
  }

  TCanvas *cPP_NLO = new TCanvas("cPP_NLO", "cPP_NLO", 800, 600);
  gPad->SetLogy();
  grDztNLOmedianpp->Draw("hist same c ");
  hZt_MC_Gen[0]->Draw("hist same c");
  cPP_NLO->Print(dirPlot + Form("/Checkpp_NLO.pdf"));
  Int_t nAssocNLO = 5;
  double assocZtNLO[] = {0.15, 0.2, 0.3, 0.4, 0.6, 1.0};
  TH1F *hPbPb_NLO[nCen];
  TH1F *hPbPbPYTHIA_NLO[nCen];
  TH1F *hPbPb_NLOSyst[nCen];

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    hPbPb_NLO[iCen] = new TH1F(Form("hPbPb%d_%d_NLO", cenBins[iCen], cenBins[iCen + 1]), Form("hPbPb%d_%d_NLO", cenBins[iCen], cenBins[iCen + 1]), nAssocNLO, assocZtNLO);
    hPbPb_NLOSyst[iCen] = new TH1F(Form("hPbPb%d_%d_NLOSyst", cenBins[iCen], cenBins[iCen + 1]), Form("hPbPb%d_%d_NLOSyst", cenBins[iCen], cenBins[iCen + 1]), nAssocNLO, assocZtNLO);

    hPbPbPYTHIA_NLO[iCen] = new TH1F(Form("hPbPbPYTHIA%d_%d_NLO", cenBins[iCen], cenBins[iCen + 1]), Form("hPbPbPYTHIA%d_%d_NLO", cenBins[iCen], cenBins[iCen + 1]), nAssocNLO, assocZtNLO);
    //hPbPbPYTHIA_NLOSyst[iCen] = new TH1F(Form("hPbPbPYTHIA%d_%d_NLOSyst", cenBins[iCen], cenBins[iCen + 1]), Form("hPbPb%d_%d_NLOSyst", cenBins[iCen], cenBins[iCen + 1]), nAssocNLO, assocZtNLO);
    for (int ibin = 0; ibin < nAssocNLO; ibin++)
    {
      cout << hPbPb_NLO[iCen]->GetBinCenter(ibin + 1) << "____" << hZtCent[iCen]->GetBinCenter(ibin + 2) << endl;
      cout << hPbPb_NLO[iCen]->GetBinContent(ibin + 1) << "____" << grDztNLOmedianpp->GetBinContent(ibin + 1) << " Ratio: " << hPbPb_NLO[iCen]->GetBinContent(ibin + 1) / grDztNLOmedianpp->GetBinContent(ibin + 1) << endl;

      hPbPb_NLO[iCen]->SetBinContent(ibin + 1, hZtCent[iCen]->GetBinContent(ibin + 2));
      hPbPb_NLO[iCen]->SetBinError(ibin + 1, hZtCent[iCen]->GetBinError(ibin + 2));
      hPbPbPYTHIA_NLO[iCen]->SetBinContent(ibin + 1, hZt_MC_Rec[iCen]->GetBinContent(ibin + 2));
      hPbPbPYTHIA_NLO[iCen]->SetBinError(ibin + 1, hZt_MC_Rec[iCen]->GetBinError(ibin + 2));
      hPbPb_NLOSyst[iCen]->SetBinContent(ibin + 1, hSystZt[iCen]->GetBinContent(ibin + 2));
      hPbPb_NLOSyst[iCen]->SetBinError(ibin + 1, hSystZt[iCen]->GetBinError(ibin + 2));
      cout << hPbPb_NLOSyst[iCen]->GetBinError(ibin + 1) << endl;
      cout << grDztNLOmedianpp->GetBinError(ibin + 1) << endl;
    }
    hPbPb_NLO[iCen]->Divide(grDztNLOmedianpp);
    hPbPbPYTHIA_NLO[0]->Divide(histDztPbPbCOLBTmedian[0]);
    hPbPb_NLOSyst[iCen]->Divide(grDztNLOmedianpp);

    PlotStyle(hPbPb_NLO[iCen], kMarkCen[iCen], 1, kColorMark[iCen], kColorMarkFill[iCen], "#it{z}_{T}", "Ratio", false);
    // PlotStyle(hPbPb_NLO[iCen], kMarkCen[iCen], 1, kColorMark[iCen], kColorMarkFill[iCen], "#it{z}_{T}", "#it{I}_{NLO} = Pb#font[122]{-}Pb/pQCD NLO", false);
    PlotStyle(hPbPb_NLOSyst[iCen], kMarkCen[iCen], 1, kColorMark[iCen], kColorMarkFill[iCen], "#it{z}_{T}", "Ratio", true);
    // PlotStyle(hPbPb_NLOSyst[iCen], kMarkCen[iCen], 1, kColorMark[iCen], kColorMarkFill[iCen], "#it{z}_{T}", "#it{I}_{NLO} = Pb#font[122]{-}Pb/pQCD NLO", true);
  }

  hPbPb_NLO[0]->GetXaxis()->SetRangeUser(0.15, 1.0);
  hPbPb_NLO[1]->GetXaxis()->SetRangeUser(0.15, 0.6);
  hPbPb_NLO[2]->GetXaxis()->SetRangeUser(0.15, 0.6);
  hPbPb_NLOSyst[0]->GetXaxis()->SetRangeUser(0.15, 1.0);
  hPbPb_NLOSyst[1]->GetXaxis()->SetRangeUser(0.15, 0.6);
  hPbPb_NLOSyst[2]->GetXaxis()->SetRangeUser(0.15, 0.6);

  hPbPbPYTHIA_NLO[0]->GetXaxis()->SetRangeUser(0.15, 1.0);
  hPbPbPYTHIA_NLO[1]->GetXaxis()->SetRangeUser(0.15, 0.6);
  hPbPbPYTHIA_NLO[2]->GetXaxis()->SetRangeUser(0.15, 0.6);

  TCanvas *cPbPbPYTHIA_NLORatio[nCen];
  TCanvas *cPbPb_NLORatio[nCen];
  TLatex *latPbPb_NLO[nCen];
  TLegend *legPbPb_NLOratio[nCen];
  TH1F *hGeneralRatio = new TH1F("hGeneralRatio", "hGeneralRatio", nZtBinThin, assocZtThinner);
  PlotStyle(hGeneralRatio, 20, 1, kWhite, kWhite, " #it{z}_{T} ", " #it{I}_{pQCD NLO} ", false);
  TGraph *lineNLO = DrawLine(lineNLO, 0, 0.5, 1.2, 0.5);
  TGraph *lineNLO1 = DrawLine(lineNLO1, 0, 1, 1.2, 1);
  TLegend *legRatioNLO = LegStd(legRatioNLO, 0.22, 0.86, 0.42, 0.94);
  // legRatioNLO->SetHeader("#it{I}_{pQCD NLO} = #frac{Pb#font[122]{-}Pb}{pQCD NLO pp} ");
  legRatioNLO->SetTextSize(0.03);
  legRatioNLO->SetHeader("#it{I}_{pQCD NLO} = #frac{Pb#font[122]{-}Pb}{pQCD NLO pp} ");
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    cPbPb_NLORatio[iCen] = canvasStd(Form("cPbPb_NLORatioCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), 1, 1);
    // new TCanvas(Form("cPbPb_NLORatioCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("cPbPb_NLORatioCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), 800, 600);
    cPbPb_NLORatio[iCen]->cd(iCen + 1);
    hGeneralRatio->GetYaxis()->SetRangeUser(-0.85, 2.2);
    hGeneralRatio->SetTitle(" ");

    hGeneralRatio->GetYaxis()->SetTitleSize(0.055);
    hGeneralRatio->GetYaxis()->SetTitleOffset(1.);
    hGeneralRatio->GetYaxis()->SetLabelSize(0.045);
    hGeneralRatio->GetXaxis()->SetRangeUser(0.05, 1.05);
    hGeneralRatio->GetXaxis()->SetLabelSize(0.045);
    hGeneralRatio->GetXaxis()->SetTitleSize(0.055);
    hGeneralRatio->Draw("histsame");
    hPbPb_NLOSyst[iCen]->Draw("E2Psame ");
    hPbPb_NLO[iCen]->Draw("EPX0same");
    grIaaNLOmedian[iCen]->Draw("pl3 same");

    if (iCen == 0)
    {
      legPbPb_NLOratio[iCen] = LegStd(legPbPb_NLOratio[iCen], 0.14, 0.16, 0.36, 0.40);
      legPbPb_NLOratio[iCen]->AddEntry(grDztNLOmedian[0], "pQCD NLO,", "lf");
      legPbPb_NLOratio[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      legPbPb_NLOratio[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
      histIaaCOLBTmedianSyst[0]->Draw("X0sameE3 ");
      histIaaCOLBTmedian[0]->SetFillStyle(0);
      histIaaCOLBTmedian[0]->Draw("histsameL ");
      legPbPb_NLOratio[iCen]->AddEntry(histIaaCOLBTmedianSyst[0], "CoLBT-hydro,", "lf");
      legPbPb_NLOratio[iCen]->AddEntry((TObject *)0, "X. N. Wang et al.", "");
    }
    else if (iCen == 1)
    {
      legPbPb_NLOratio[iCen] = LegStd(legPbPb_NLOratio[iCen], 0.14, 0.16, 0.36, 0.40);
      legPbPb_NLOratio[iCen]->AddEntry(grDztNLOmedian[1], "pQCD NLO,", "lf");
      legPbPb_NLOratio[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      legPbPb_NLOratio[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
      histIaaCOLBTmedianSyst[1]->Draw("X0sameE3 ");
      histIaaCOLBTmedian[1]->SetFillStyle(0);
      histIaaCOLBTmedian[1]->Draw("histsameL ");
      legPbPb_NLOratio[iCen]->AddEntry(histIaaCOLBTmedianSyst[1], "CoLBT-hydro,", "lf");
      legPbPb_NLOratio[iCen]->AddEntry((TObject *)0, "X. N. Wang et al.", "");
    }
    else if (iCen == 2)
    {
      legPbPb_NLOratio[iCen] = LegStd(legPbPb_NLOratio[iCen], 0.14, 0.24, 0.36, 0.40);
      legPbPb_NLOratio[iCen]->AddEntry(grDztNLOmedian[2], "pQCD NLO,", "lf");
      legPbPb_NLOratio[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      legPbPb_NLOratio[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
    }
    latPbPb_NLO[iCen] = LatexStdISO(latPbPb_NLO[iCen], 0.440, 0.86, 0.04, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true);
    legPbPb_NLOratio[iCen]->Draw("same");
    lineNLO->Draw("l");
    lineNLO1->Draw("l");
    legRatioNLO->Draw("same");
    cPbPb_NLORatio[iCen]->Print(dirPlot + Form("/I_NLOCen%d_%d.pdf", cenBins[iCen], cenBins[iCen + 1]));
    cPbPbPYTHIA_NLORatio[iCen] = canvasStd(Form("cPbPbPYTHIA_NLORatioCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), 1, 1);
    hPbPbPYTHIA_NLO[iCen]->Draw("HIST");
    cPbPbPYTHIA_NLORatio[iCen]->Print(dirPlot + Form("/RatioPYTHIANLOCen%d_%d.pdf", cenBins[iCen], cenBins[iCen + 1]));
  }

  TCanvas *cPbPb_NLORatioAll = canvasStd(Form("cPbPb_NLORatioAll"), 1, 1);
  TLatex *latPbPb_NLOAll;
  TLegend *legPbPb_NLOratioAll = LegStd(legPbPb_NLOratioAll, 0.54, 0.52, 0.92, 0.72);
  legPbPb_NLOratioAll->SetNColumns(2);
  legPbPb_NLOratioAll->SetTextSize(0.038);
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    // hGeneralRatio->GetYaxis()->SetNdivisions(505);
    hGeneralRatio->GetYaxis()->SetRangeUser(-0.01, 2.3);
    // hGeneralRatio->GetYaxis()->SetNdivisions(505);
    //  hGeneralRatio->GetXaxis()->SetRangeUser(0.15, 1.0);
    hGeneralRatio->SetTitle(" ");
    hGeneralRatio->GetYaxis()->CenterTitle(false);
    hGeneralRatio->GetYaxis()->SetRangeUser(-0.01, 2.3);
    hGeneralRatio->GetYaxis()->SetTitleSize(0.055);
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
    hPbPb_NLOSyst[iCen]->SetLineColor(kWhite);
    legPbPb_NLOratioAll->AddEntry(hPbPb_NLO[iCen], Form("%d#font[122]{-}%d%% stat.", cenBins[iCen], cenBins[iCen + 1]), "ep");
    legPbPb_NLOratioAll->AddEntry(hPbPb_NLOSyst[iCen], "syst. unc.", "f");
  }
  latPbPb_NLOAll = LatexStdISORatio(latPbPb_NLOAll, 0.48, 0.94, 0.040, cenBins[0], cenBins[1], ptMin, ptMax, false);
  legPbPb_NLOratioAll->Draw("same");
  legRatioNLO->Draw("same");
  lineNLO->Draw("l");
  lineNLO1->Draw("l");
  cPbPb_NLORatioAll->Print(dirPlot + Form("/I_NLOCenAll.pdf"));

  TH1F *hGeneral = new TH1F("hGeneral", "hGeneral", nZtBinThin, assocZtThinner);
  PlotStyle(hGeneral, 20, 1, kWhite, kWhite, "#it{z}_{T}", "1 / #it{N}^{ #it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d #it{z}_{T}", false);
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
    PlotStyle(hZtDiffCent[iCen], kMarkCen[iCen], 1, kColorMark[iCen], kColorMarkFill[iCen], "#it{z}_{T}", "1 / #it{N}^{ #it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d #it{z}_{T}", false);
    hZtDiffCent[iCen]->SetTitle(" ");
    hZtDiffCentSys[iCen]->SetTitle(" ");
    hZtDiffCent[iCen]->GetYaxis()->SetRangeUser(1e-3, 90);
    hZtDiffCentSys[iCen]->SetMarkerStyle(kMarkCen[iCen]);
    hZtDiffCentSys[iCen]->SetMarkerColor(kColorMark[iCen]);
    hZtDiffCentSys[iCen]->SetFillStyle(0);
    hGeneral->GetYaxis()->SetRangeUser(2 * 1e-4, 99);
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
    hGeneral->GetYaxis()->SetRangeUser(2 * 1e-4, 99);
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
    cDiffCent_Pythia[iCen] = new TCanvas(Form("cDiffCent_Pythia%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("cDiffCent_Pythia%d_%d", cenBins[iCen], cenBins[iCen + 1]), 600, 800);
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

    hZt_MC_Gen[iCen]->SetFillStyle(0);
    // hZt_MC_Gen[iCen]->SetLineStyle(0);
    hZt_MC_Gen[iCen]->SetLineWidth(10);
    hZt_MC_Gen[iCen]->SetLineColorAlpha(kBlack, 1);

    hSystZt[iCen]->SetLineColor(kWhite);
    legdiffcen[iCen]->AddEntry(hZtCent[iCen], "Pb#font[122]{-}Pb stat. unc. ", "ep");
    legdiffcen[iCen]->AddEntry(hSystZt[iCen], "Pb#font[122]{-}Pb syst. unc. ", "f");
    legdiffcen[iCen]->Draw("same");
    legdiffcenpp[iCen]->Draw("same");
    TLatex *latdiffCent1 = LatexStdISO(latdiffCent1, 0.320, 0.94, 0.04, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true);
    // cDiffCent_Pythia[iCen]->Print(dirPlot + Form("/withModelNoppztPbPb_pp_Cen%d_%d%s.pdf", cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
    cDiffCent_Pythia[iCen]->Print(dirPlot + Form("/ztPbPb_Cen%d_%d%s.pdf", cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
  }

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
  //  PlotStyle(hGeneralRatio, 20, 1, kWhite, kWhite,"#it{z}_{T}", "1 / #it{N}^{ #it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d #it{z}_{T}", false);
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
      legZTPbPbppSingleMod[iCen]->AddEntry(grDztNLOmedian[0], "pQCD NLO,", "lf");
      legZTPbPbppSingleMod[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      legZTPbPbppSingleMod[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
      legZTPbPbppSingleMod[iCen]->AddEntry(histDztPbPbCOLBTmedianSyst[0], "CoLBT-hydro,", "lf");
      legZTPbPbppSingleMod[iCen]->AddEntry((TObject *)0, "X. N. Wang et al.", "");
    }
    else if (iCen == 1)
    {
      legZTPbPbppSingleMod[iCen] = LegStd(legZTPbPbppSingleMod[iCen], 0.170, 0.02, 0.350, 0.305);
      legZTPbPbppSingleMod[iCen]->AddEntry(grDztNLOmedian[1], "pQCD NLO,", "lf");
      legZTPbPbppSingleMod[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      legZTPbPbppSingleMod[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
      legZTPbPbppSingleMod[iCen]->AddEntry(histDztPbPbCOLBTmedianSyst[1], "CoLBT-hydro,", "lf");
      legZTPbPbppSingleMod[iCen]->AddEntry((TObject *)0, "X. N. Wang et al.", "");
    }
    else if (iCen == 2)
    {
      legZTPbPbppSingleMod[iCen] = LegStd(legZTPbPbppSingleMod[iCen], 0.170, 0.02, 0.350, 0.195);
      legZTPbPbppSingleMod[iCen]->AddEntry(grDztNLOmedian[2], "pQCD NLO,", "lf");
      legZTPbPbppSingleMod[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      legZTPbPbppSingleMod[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
    }
    legZTPbPbppSingle[iCen]->Draw("same");
    legZTPbPbppSingleMod[iCen]->Draw("same");
    legRatioSingleMod[iCen] = LegStd(legRatioSingleMod[iCen], 0.40, 0.80, 0.80, 0.995);
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
    legNLO[iCen]->AddEntry(h3[iCen], "#frac{Pb#font[122]{-}Pb}{PYTHIA 8 pp}", "ep");
    legRatioSingleMod[iCen]->AddEntry(grIaaNLOmedian[iCen], "#it{I}_{AA}: pQCD NLO (X. N. Wang and M. Xie)", "lf");
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
    else if (iCen == 2)
      legRatioSingleMod[iCen]->SetY1(0.8975);

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
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    legNLOpQCD[iCen] = LegStd(legNLOpQCD[iCen], 0.16, 0.80, 0.35, 0.995);
    legNLOpQCD[iCen]->SetTextSize(0.06);
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
    // histDztPbPbCOLBTmedian[0]->SetMarkerStyle(0);
    // histDztPbPbCOLBTmedianSyst[0]->SetMarkerStyle(0);
    hSystZt[iCen]->Draw("E2P same");
    hZtCent[iCen]->Draw("EPX0same");
    grDztNLOmedianpp->SetLineWidth(5);
    grDztNLOmedianpp->Draw("histsameC");
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
    legZTPbPbppSingleNLOpQCD[iCen]->AddEntry(grDztNLOmedianpp, Form("pQCD NLO pp"), "l");
    legZTPbPbppSingleNLOpQCD[iCen]->AddEntry(hZtCent[iCen], "Pb#font[122]{-}Pb stat. unc.", "ep");
    legZTPbPbppSingleNLOpQCD[iCen]->AddEntry(hSystZt[iCen], "Pb#font[122]{-}Pb syst. unc.", "f");

    if (iCen == 0)
    {
      legZTPbPbppSingleModNLO[iCen] = LegStd(legZTPbPbppSingleModNLO[iCen], 0.170, 0.02, 0.350, 0.305);
      legZTPbPbppSingleModNLO[iCen]->AddEntry(grDztNLOmedian[0], "pQCD NLO,", "lf");
      legZTPbPbppSingleModNLO[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      legZTPbPbppSingleModNLO[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
      legZTPbPbppSingleModNLO[iCen]->AddEntry(histDztPbPbCOLBTmedianSyst[0], "CoLBT-hydro,", "lf");
      legZTPbPbppSingleModNLO[iCen]->AddEntry((TObject *)0, "X. N. Wang et al.", "");
    }
    else if (iCen == 1)
    {
      legZTPbPbppSingleModNLO[iCen] = LegStd(legZTPbPbppSingleModNLO[iCen], 0.170, 0.02, 0.350, 0.305);
      legZTPbPbppSingleModNLO[iCen]->AddEntry(grDztNLOmedian[1], "pQCD NLO,", "lf");
      legZTPbPbppSingleModNLO[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      legZTPbPbppSingleModNLO[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
      legZTPbPbppSingleModNLO[iCen]->AddEntry(histDztPbPbCOLBTmedianSyst[1], "CoLBT-hydro,", "lf");
      legZTPbPbppSingleModNLO[iCen]->AddEntry((TObject *)0, "X. N. Wang et al.", "");
    }
    else if (iCen == 2)
    {
      legZTPbPbppSingleModNLO[iCen] = LegStd(legZTPbPbppSingleModNLO[iCen], 0.170, 0.134, 0.350, 0.305);
      legZTPbPbppSingleModNLO[iCen]->AddEntry(grDztNLOmedian[2], "pQCD NLO,", "lf");
      legZTPbPbppSingleModNLO[iCen]->AddEntry((TObject *)0, "CT18A + EPPS21 nPDFs, KKP FFs,", "");
      legZTPbPbppSingleModNLO[iCen]->AddEntry((TObject *)0, "X. N. Wang and M. Xie", "");
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
    // hGeneralRatio->GetYaxis()->SetTitle("Pb#font[122]{-}Pb/pQCD NLO pp");
    hGeneralRatio->GetYaxis()->SetTitle("Ratio");
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
    legNLOpQCD[iCen]->AddEntry(h3[iCen], "#frac{Pb#font[122]{-}Pb}{pQCD NLO pp} ", "ep");
    if (iCen == 0)
    {
      histIaaCOLBTmedianSyst[0]->Draw("X0sameE3 ");
      histIaaCOLBTmedian[0]->SetFillStyle(0);
      histIaaCOLBTmedian[0]->Draw("histsameL ");
    }
    else if (iCen == 1)
    {
      histIaaCOLBTmedianSyst[1]->Draw("X0sameE3 ");
      histIaaCOLBTmedian[1]->SetFillStyle(0);
      histIaaCOLBTmedian[1]->Draw("histsameL ");
    }

    legNLOpQCD[iCen]->Draw("same");
    //  cPbPbpp->Update();
    //  pad2->Update();
    legRatioSingleMod[iCen]->Draw("same");
    // h3[iCen]->GetXaxis()->SetTitle("z_{#rm{T}} = p_{#rm{T}}^{had}}/p_{#rm{T}}^{#gamma}}");
    TGraph *linea = DrawLine(linea, 0, 0.5, 1.2, 0.5);
    linea->Draw("l");

    TGraph *lineb = DrawLine(lineb, 0, 1, 1.2, 1);
    lineb->Draw("l");
    cPbPbppRatioSingleNLOpQCD[iCen]->Print(dirPlot + Form("/SingleztPbPb_pp_Cen%d_%d%s_ppNLOpQCD.pdf", cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
  }

  // TCanvas *cAllZt = new TCanvas("cAllZt", "cAllZt", 800, 600);
  TCanvas *cAllZt = canvasStd("cAllZt", 1, 1);
  TLegend *legdiffcenZtPYTHIA = LegStd(legdiffcenZtPYTHIA, 0.14, 0.30, 0.4, 0.35);
  TLegend *legdiffcenZt = LegStd(legdiffcenZt, 0.14, 0.12, 0.5, 0.30);
  // cAllZt->cd();
  // cAllZt->SetTopMargin(0.015);
  // cAllZt->SetRightMargin(0.02);
  // cAllZt->SetLeftMargin(0.12);
  // cAllZt->SetBottomMargin(0.11);
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
  legdiffcenZt->SetNColumns(2);
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    legdiffcenZt->AddEntry(hZtCent[iCen], Form(" %d#font[122]{-}%d%% stat. ", cenBins[iCen], cenBins[iCen + 1]), "ep");
    legdiffcenZt->AddEntry(hSystZt[iCen], Form(" syst. unc"), "f");
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

  TCanvas *cAllZtNLOpQCD = canvasStd("cAllZtNLOpQCD", 1, 1);
  TLegend *legdiffcenZtNLOpQCD = LegStd(legdiffcenZtNLOpQCD, 0.14, 0.30, 0.4, 0.35);
  // cAllZtNLOpQCD->cd();
  // cAllZtNLOpQCD->SetTopMargin(0.015);
  // cAllZtNLOpQCD->SetRightMargin(0.02);
  // cAllZtNLOpQCD->SetLeftMargin(0.12);
  // cAllZtNLOpQCD->SetBottomMargin(0.11);
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
  grDztNLOmedianpp->Draw("HISTSAMEC ");
  legdiffcenZtNLOpQCD->AddEntry(grDztNLOmedianpp, "pQCD NLO pp ", "l");
  legdiffcenZtNLOpQCD->Draw("SAME");
  hSystZt[0]->Draw("samee2");
  hZtCent[0]->Draw("EP X0 same");
  hSystZt[2]->Draw("samee2");
  hZtCent[2]->Draw("EP X0 same");
  hSystZt[1]->Draw("samee2");
  hZtCent[1]->Draw("EP X0 same");
  legdiffcenZt->Draw("SAME");
  TLatex *latAliceZtOnlyNLOpQCD = LatexStdISORatio(latAliceZtOnlyNLOpQCD, 0.420, 0.92, 0.040, cenBins[0], cenBins[1], ptMin, ptMax, false);
  cAllZtNLOpQCD->Print(dirPlot + Form("/NLOpQCDZtAllCent030%s.pdf", sPtAll.Data()));

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

  // Zt STAR and PHENIX
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
  PlotStyle(hzTPHENIX, 20, 1.2, kAzure + 2, kWhite, "#it{z}_{T}", "1 / #it{N}^{ #it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d #it{z}_{T}", false);
  TH1F *hGeneralCMS = new TH1F("hGeneral", "hGeneral", 10, -0.05, 1.1);
  PlotStyle(hGeneralCMS, 71, 1, kWhite, kWhite, "#it{z}_{T}", "1 / #it{N}^{ #it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d #it{z}_{T}", false);

  TCanvas *cPHENIX = new TCanvas("cPHENIX", "cPHENIX", 2 * 800, 1 * 600);
  cPHENIX->Divide(2, 1);
  cPHENIX->cd(1);
  cPHENIX->cd(1)->SetTopMargin(0.015);
  cPHENIX->cd(1)->SetRightMargin(0.02);
  cPHENIX->cd(1)->SetLeftMargin(0.12);
  cPHENIX->cd(1)->SetBottomMargin(0.11);
  hGeneralCMS->SetTitle(" ");
  hGeneralCMS->GetYaxis()->SetTitleSize(0.045);
  hGeneralCMS->GetXaxis()->SetTitleSize(0.05);
  hGeneralCMS->GetYaxis()->SetLabelSize(0.04);
  hGeneralCMS->GetXaxis()->SetLabelSize(0.045);
  hGeneralCMS->GetYaxis()->SetLabelOffset(0.02);
  hGeneralCMS->GetXaxis()->SetLabelOffset(0.02);
  gPad->SetLogy();
  hGeneralCMS->GetYaxis()->SetRangeUser(1e-4, 99);
  // grSTARSys->GetXaxis()->SetRangeUser(0, 1.00);
  grSTARSys->GetYaxis()->SetTitle("1 / #it{N}^{ #it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d #it{z}_{T}");
  grSTARSys->GetXaxis()->SetTitle("#it{z}_{T}");
  grSTARSys->SetTitle(" ");
  hzTPHENIXSys->SetTitle(" ");
  hzTPHENIXSys->SetLineWidth(2);
  hzTPHENIXSys->SetLineColorAlpha(kAzure + 2, 1.00);
  hzTPHENIXSys->SetFillColorAlpha(kGray, 0.80);
  hzTPHENIXSys->SetFillStyle(0);
  hzTPHENIXSys->SetLineWidth(2);
  hGeneralCMS->Draw("hist");
  hzTPHENIXSys->Draw("samee2");
  hzTPHENIX->Draw("EP X0 same");
  hSystZt[0]->SetMarkerColor(kRed + 1);
  hSystZt[0]->SetFillColor(kRed + 1);
  hSystZt[0]->SetLineColor(kRed + 1);
  hSystZt[0]->SetLineWidth(2);
  hZtCent[0]->SetMarkerColor(kRed + 1);
  hZtCent[0]->SetLineColor(kRed + 1);
  hZtCent[0]->SetMarkerSize(1.2);
  hSystZt[0]->SetFillStyle(0);
  grSTARSys->SetLineColorAlpha(kTeal + 4, 1.00);
  grSTARSys->SetFillStyle(0);
  grSTARSys->SetLineWidth(2);
  grSTARSys->Draw(" p2 same");
  grSTAR->Draw("p SAME");
  hSystZt[0]->Draw("samee2");
  hZtCent[0]->Draw("EP X0 same");
  TLatex *latALICEPHENIX = LatexStdISORatio(latALICEPHENIX, 0.400, 0.92, 0.045, cenBins[0], cenBins[1], ptMin, ptMax, true);
  TLegend *legALICE1 = LegStd(legALICE1, 0.140, 0.16, 0.40, 0.40);
  legALICE1->SetTextSize(0.042);
  legALICE1->AddEntry(hZtCent[0], "ALICE, stat. unc.", "ep");
  legALICE1->AddEntry(grSTAR, "STAR, stat. unc. ", "ep");
  legALICE1->AddEntry(hzTPHENIX, "PHENIX, stat. unc. ", "ep");
  legALICE1->AddEntry(hSystZt[0], " syst. unc.", "f");
  legALICE1->Draw("same");
  cPHENIX->cd(2);
  TLegend *legSTAR1 = LegStd(legSTAR1, 0, 0.68, 0.10, 0.98);
  legSTAR1->SetHeader("STAR, Phys.Lett.B 760 (2016) 689-696");
  legSTAR1->AddEntry((TObject *)0, "#bf{0#font[122]{-}12%} Au#font[122]{-}Au, #sqrt{#it{s}_{NN}} = 200 GeV ", "");
  legSTAR1->AddEntry((TObject *)0, "|#Delta#it{#varphi}_{#it{#gamma}#font[122]{-}h} #font[122]{-} #it{#pi}| #leq 1.4 ", "");
  legSTAR1->AddEntry((TObject *)0, "12 < #it{p}_{T}^{ #it{#gamma}} < 20 GeV/#it{c} #otimes #it{p}_{T}^{ h} > 1.2 GeV/#it{c} ", "");
  legSTAR1->Draw("same");
  TLegend *legPHENIX = LegStd(legPHENIX, 0, 0.32, 0.10, 0.62);
  legPHENIX->SetHeader("PHENIX, PRL 111, 032301 (2013)");
  legPHENIX->AddEntry((TObject *)0, "#bf{0#font[122]{-}40%} Au#font[122]{-}Au, #sqrt{#it{s}_{NN}} = 200 GeV ", "");
  legPHENIX->AddEntry((TObject *)0, "|#Delta#it{#varphi}_{#it{#gamma}#font[122]{-}h} #font[122]{-} #it{#pi}| < #it{#pi}/2, |#it{y}| < 0.35 ", "");
  legPHENIX->AddEntry((TObject *)0, "5 < #it{p}_{T}^{ #it{#gamma}} < 9 GeV/#it{c} #otimes 0.5 < #it{p}_{T}^{ h} < 7 GeV/#it{c} ", "");
  legPHENIX->Draw("same");

  cPHENIX->Print(dirPlot + Form("/ALICESTARPHENIX%s.pdf", sPtAll.Data()));

  // Iaa STAR and PHENIX

  TH1F *hGeneralIaa = new TH1F("hGeneralIaa", "hGeneralIaa", 10, -0.05, 1.05);
  PlotStyle(hGeneralIaa, 20, 0, kWhite, kWhite, "#it{z}_{T}", "#it{I}_{pQCD NLO}, #it{I}_{AA}", false);
  // PlotStyle(hGeneralIaa, 20, 0, kWhite, kWhite, "#it{z}_{T}", "#it{I}_{PYTHIA}, #it{I}_{AA}", false);
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
  grIaaSTAR->SetMarkerStyle(20);
  grIaaSTAR->SetMarkerColor(kTeal + 4);
  grIaaSTAR->SetLineColor(kTeal + 4);
  grIaaSTAR->SetLineWidth(2);
  grIaaSTARSys->SetMarkerStyle(20);
  grIaaSTARSys->SetMarkerColor(kTeal + 4);
  grIaaSTARSys->SetLineColor(kTeal + 4);
  grIaaSTARSys->SetFillColor(kTeal + 4);
  grIaaSTARSys->SetLineWidth(2);
  grIaaSTARSys->SetMarkerSize(1.2);
  // grIaaSTARSys->SetFillStyle(3001);
  grIaaSTARSys->SetFillStyle(0);
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
  PlotStyle(hIaaPHENIX, 20, 1, kAzure + 2, kAzure + 2, "#it{z}_{T}", "#it{I}_{PYTHIA}, #it{I}_{AA}", false);
  PlotStyle(hIaaPHENIXSys, 20, 1, kAzure + 2, kAzure + 2, "#it{z}_{T}", "#it{I}_{PYTHIA}, #it{I}_{AA}", true);

  TCanvas *cIaaPHENIX = new TCanvas("cIaaPHENIX", "cIaaPHENIX", 2 * 800, 1 * 600);
  cIaaPHENIX->Divide(2, 1);
  cIaaPHENIX->cd(1);
  cIaaPHENIX->cd(1)->SetTopMargin(0.015);
  cIaaPHENIX->cd(1)->SetRightMargin(0.02);
  cIaaPHENIX->cd(1)->SetLeftMargin(0.12);
  cIaaPHENIX->cd(1)->SetBottomMargin(0.11);
  hGeneralIaa->SetTitle(" ");
  hGeneralIaa->GetYaxis()->SetTitleSize(0.052);
  hGeneralIaa->GetXaxis()->SetTitleSize(0.05);
  hGeneralIaa->GetYaxis()->SetLabelSize(0.04);
  hGeneralIaa->GetXaxis()->SetLabelSize(0.045);
  hGeneralIaa->GetYaxis()->SetTitleOffset(1.25);
  hGeneralIaa->GetYaxis()->SetLabelOffset(0.02);
  hGeneralIaa->GetXaxis()->SetLabelOffset(0.02);
  hGeneralIaa->GetYaxis()->SetRangeUser(-0.3, 3.2);
  hGeneralIaa->SetLineWidth(0);
  hGeneralIaa->Draw("histsame");
  hIaaPHENIXSys->SetFillStyle(0);
  hIaaPHENIXSys->SetLineWidth(2);
  hPbPb_NLO[0]->SetLineColor(kRed + 1);
  hPbPb_NLO[0]->SetMarkerColor(kRed + 1);
  hPbPb_NLO[0]->SetMarkerSize(1.3);
  hPbPb_NLOSyst[0]->SetFillColor(kRed + 1);
  hPbPb_NLOSyst[0]->SetLineColor(kRed + 1);
  hPbPb_NLOSyst[0]->SetFillStyle(0);
  hPbPb_NLOSyst[0]->SetLineWidth(2);
  hPbPb_NLOSyst[0]->SetMarkerSize(1.3);
  // h3[0]->SetLineColor(kRed + 1);
  // h3[0]->SetMarkerColor(kRed + 1);
  // h3[0]->SetMarkerSize(1.3);
  // hsyst[0]->SetLineWidth(2);
  // hsyst[0]->SetFillColor(kRed + 1);
  // hsyst[0]->SetLineColor(kRed + 1);
  // hsyst[0]->SetFillStyle(0);
  // hsyst[0]->SetMarkerSize(1.3);
  hIaaPHENIXSys->Draw("samee2");
  hIaaPHENIX->Draw("EP X0 same");
  grIaaSTARSys->Draw("p2 same");
  grIaaSTAR->Draw("p same");
  hPbPb_NLOSyst[0]->Draw("p samee2");
  hPbPb_NLO[0]->Draw("EP X0 same");
  // hsyst[0]->Draw("p samee2");
  // h3[0]->Draw("EP X0 same");
  cIaaPHENIX->cd(2);
  legPHENIX->Draw("same");
  legSTAR1->Draw("same");
  cIaaPHENIX->cd(1);
  TLatex *latALICEIaa = LatexStdISORatio(latALICEIaa, 0.400, 0.920, 0.045, cenBins[0], cenBins[1], ptMin, ptMax, true);
  TLegend *legALICE3 = LegStd(legALICE3, 0.60, 0.50, 0.84, 0.70);
  legALICE3->SetTextSize(0.042);
  legALICE3->AddEntry(hPbPb_NLO[0], "ALICE, stat. unc.", "ep");
  // legALICE3->AddEntry(h3[0], "ALICE, stat. unc.", "ep");
  legALICE3->AddEntry(grIaaSTAR, "STAR, stat. unc.", "ep");
  legALICE3->AddEntry(hIaaPHENIX, "PHENIX, stat. unc.", "ep");
  legALICE3->AddEntry(hPbPb_NLOSyst[0], " syst. unc.", "f");
  // legALICE3->AddEntry(hsyst[0], " syst. unc.", "f");
  legALICE3->Draw("same");
  TGraph *linePHENIX = DrawLine(linePHENIX, -0.05, 0.5, 1.05, 0.5);
  linePHENIX->Draw("l");
  TGraph *linePHENIX1 = DrawLine(linePHENIX1, -0.05, 1, 1.05, 1);
  linePHENIX1->Draw("l");

  cIaaPHENIX->Print(dirPlot + Form("/ppNLOpQCDIaaALICESTARPHENIX%s.pdf", sPtAll.Data()));

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
  PlotStyle(hzTCMS, 20, 1.2, kBlue + 1, kBlue + 1, "#it{z}_{T}", "1 / #it{N}^{ #it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d #it{z}_{T}", false);
  PlotStyle(hzTCMSSys, 20, 1.2, kBlue + 1, kBlue + 1, "#it{z}_{T}", "1 / #it{N}^{ #it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d #it{z}_{T}", true);

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
  PlotStyle(hIaaCMS, 20, 1.2, kBlue + 1, kBlue + 1, "#it{z}_{T}", "#it{I}_{PYTHIA}, #it{I}_{AA}", false);
  PlotStyle(hIaaCMSSys, 20, 1.2, kBlue + 1, kBlue + 1, "#it{z}_{T}", "#it{I}_{PYTHIA}, #it{I}_{AA}", true);

  /*TCanvas *cCMS = new TCanvas("cCMS", "cCMS", 2 * 800, 1 * 600);
  cCMS->Divide(2, 1);
  cCMS->cd(1);
  gPad->SetLogy();
  hzTCMSSys->SetTitle(" ");
  hzTCMSSys->GetYaxis()->SetTitle("1 / #it{N}^{ #it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d #it{z}_{T}");
  hzTCMSSys->GetXaxis()->SetTitle("#it{z}_{T}");
  hzTCMSSys->SetFillStyle(0);
  hGeneral->GetXaxis()->SetRangeUser(0, 1.05);
  hGeneral->GetYaxis()->SetRangeUser(5e-4, 100);
  hGeneral->Draw("hist");
  hSystZt[0]->Draw("samee2");
  hZtCent[0]->Draw("EP X0 same");
  hzTCMSSys->Draw("samee2");
  hzTCMS->Draw("EP X0 same");
  TLatex *latALICEcms = LatexStdISO(latALICEcms, 0.340, 0.840, 0.04, cenBins[0], cenBins[1], ptMin, ptMax, true);
  TLegend *legALICE2 = LegStd(legALICE2, 0.55, 0.60, 0.80, 0.750);
  legALICE2->AddEntry(hZtCent[0], " stat. unc.", "ep");
  legALICE2->AddEntry(hSystZt[0], " syst. unc.", "f");
  legALICE2->Draw("same");
  cCMS->cd(2);
  TLegend *legCMS = LegStd(legCMS, 0.05, 0.50, 0.60, 0.88);
  legCMS->SetHeader("CMS, Phys.Rev.Lett. 121 (2018) 712301, 2018");
  legCMS->AddEntry((TObject *)0, " 0#font[122]{-}10% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV ", "");
  legCMS->AddEntry((TObject *)0, " anti-k_{T} jet R = 0.3, #it{p}_{T}^{ jet} > 30 GeV/#it{c}, |#it{#eta}^{jet}| < 1.6 ", "");
  legCMS->AddEntry((TObject *)0, " |#it{#eta}^{ #it{#gamma}}| < 1.44 #it{p}_{T}^{ #it{#gamma}} > 60 GeV/#it{c} #otimes #it{p}_{T}^{ h} > 1 GeV/#it{c} ", "");
  legCMS->AddEntry(hzTCMS, " stat. unc. ", "ep");
  legCMS->AddEntry(hzTCMSSys, " syst. unc. ", "f");
  legCMS->Draw("same");

  cCMS->Print(dirPlot + Form("/ALICECMS%s.pdf", sPtAll.Data()));
*/
  // Iaa CMS

  /* TCanvas *cIaaCMS = new TCanvas("cIaaCMS", "cIaaCMS", 2 * 800, 1 * 600);
   cIaaCMS->Divide(2, 1);
   cIaaCMS->cd(1);
   hIaaCMSSys->SetFillStyle(0);
   hGeneralIaa->SetTitle(" ");
   hGeneralIaa->GetYaxis()->SetRangeUser(0, 3);
   hGeneralIaa->GetXaxis()->SetRangeUser(-0.1, 1.1);
   hGeneralIaa->GetYaxis()->SetTitleSize(0.045);
   hGeneralIaa->GetXaxis()->SetTitleSize(0.045);
   hGeneralIaa->GetYaxis()->SetLabelSize(0.04);
   hGeneralIaa->GetXaxis()->SetLabelSize(0.04);
   hGeneralIaa->SetLineWidth(0);
   hGeneralIaa->Draw("hist");

   hIaaCMSSys->Draw("samee2");
   hIaaCMS->Draw("EP X0 same");
   hsyst[0]->Draw("samee2");
   h3[0]->Draw("EP X0 same");
   // grSTAR->Draw("p same");
   // grSTARSys->Draw(" p2 same");
   TLatex *latALICEcms1 = LatexStdISO(latALICEcms1, 0.360, 0.84, 0.04, cenBins[0], cenBins[1], ptMin, ptMax, true);
   TLegend *legALICEIaa = LegStd(legALICEIaa, 0.50, 0.50, 0.70, 0.65);
   legALICEIaa->AddEntry(h3[0], " stat. unc.", "ep");
   legALICEIaa->AddEntry(hsyst[0], " syst. unc.", "f");
   legALICEIaa->Draw("same");
   cIaaCMS->cd(2);
   TLegend *legCMS1 = LegStd(legCMS1, 0.05, 0.50, 0.58, 0.88);
   legCMS1->SetHeader("CMS, Phys.Rev.Lett. 121 (2018) 712301, 2018");
   legCMS1->AddEntry((TObject *)0, " 0#font[122]{-}10% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV ", "");
   legCMS1->AddEntry((TObject *)0, " anti-k_{T} jet R = 0.3, #it{p}_{T}^{ jet} > 30 GeV/#it{c}, |#it{#eta}^{jet}| < 1.6 ", "");
   legCMS1->AddEntry((TObject *)0, " |#it{#eta}^{#gamma}| < 1.44 #it{p}_{T}^{ #it{#gamma}} > 60 GeV/#it{c} #otimes #it{p}_{T}^{ h} > 1 GeV/#it{c} ", "");
   legCMS1->AddEntry(hIaaCMS, " stat. unc. ", "ep");
   legCMS1->AddEntry(hIaaCMSSys, " syst. unc. ", "f");
   legCMS1->Draw("same");

   cIaaCMS->Print(dirPlot + Form("/IaaALICECMS%s.pdf", sPtAll.Data()));
 */
  // Z-hadron correlations

  // ATLAS

  // CMS

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
  PlotStyle(hzTCMS_Zhad, kFullCross, 1.2, kGreen + 3, kGreen + 3, "#it{z}_{T}", "1 / #it{N}^{ #it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d #it{z}_{T}", false);
  PlotStyle(hzTCMSSys_Zhad, kFullCross, 1.2, kGreen + 3, kGreen + 3, "#it{z}_{T}", "1 / #it{N}^{ #it{#gamma}} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d #it{z}_{T}", true);

  TH1F *hIaaCMS_Zhad = new TH1F("hIaaCMS_Zhad", "hIaaCMS_Zhad", 10, zTCMS_Zhad);
  TH1F *hIaaCMSSys_Zhad = new TH1F("hIaaCMSSys_Zhad", "hIaaCMSSys_Zhad", 10, zTCMS_Zhad);
  for (int ibin = 0; ibin < 10; ibin++)
  {
    hIaaCMS_Zhad->SetBinContent(ibin + 1, IaaCMS0_10_Zhad[9 - ibin]);
    hIaaCMS_Zhad->SetBinError(ibin + 1, IaaCMS_Stat0_10_Zhad[9 - ibin]);
    hIaaCMSSys_Zhad->SetBinContent(ibin + 1, IaaCMS0_10_Zhad[9 - ibin]);
    hIaaCMSSys_Zhad->SetBinError(ibin + 1, IaaCMS_Sys0_10_Zhad[9 - ibin]);
  }
  PlotStyle(hIaaCMS_Zhad, kFullCross, 1.2, kGreen + 3, kGreen + 3, "#it{z}_{T}", "#it{I}_{pQCD NLO}, #it{I}_{AA}", false);
  // PlotStyle(hGeneralCMSIaa, 71, 1, kBlack, "#it{z}_{T}", "#it{I}_{PYTHIA}, #it{I}_{AA}", false);
  PlotStyle(hIaaCMSSys_Zhad, kFullCross, 1.2, kGreen + 3, kGreen + 3, "#it{z}_{T}", "#it{I}_{PYTHIA}, #it{I}_{AA}", true);

  TCanvas *cCMS_Zhad = new TCanvas("cCMS_Zhad", "cCMS_Zhad", 2 * 800, 1 * 600);
  cCMS_Zhad->Divide(2, 1);
  cCMS_Zhad->cd(1);
  cCMS_Zhad->cd(1)->SetTopMargin(0.015);
  cCMS_Zhad->cd(1)->SetRightMargin(0.02);
  cCMS_Zhad->cd(1)->SetLeftMargin(0.12);
  cCMS_Zhad->cd(1)->SetBottomMargin(0.11);
  gPad->SetLogy();
  hGeneralCMS->SetTitle(" ");
  hGeneralCMS->GetYaxis()->SetTitleSize(0.045);
  hGeneralCMS->GetXaxis()->SetTitleSize(0.05);
  hGeneralCMS->GetYaxis()->SetLabelSize(0.04);
  hGeneralCMS->GetXaxis()->SetLabelSize(0.045);
  hGeneralCMS->GetYaxis()->SetLabelOffset(0.02);
  hGeneralCMS->GetXaxis()->SetLabelOffset(0.02);

  hGeneralCMS->GetYaxis()->SetRangeUser(5e-4, 100);
  hGeneralCMS->GetXaxis()->SetRangeUser(-0.02, 1.08);
  // hGeneral->Draw("hist");
  hGeneralCMS->SetLineWidth(0);
  hGeneralCMS->Draw("hist");
  hzTCMSSys_Zhad->SetFillStyle(0);
  hzTCMSSys_Zhad->SetLineWidth(2);
  hzTCMSSys_Zhad->Draw("samee2");
  hzTCMS_Zhad->Draw("EP X0 same");
  hzTCMSSys->SetFillStyle(0);
  hzTCMSSys->SetLineWidth(2);
  hzTCMSSys->Draw("samee2");
  hzTCMS->Draw("EP X0 same");
  hSystZt[0]->SetMarkerColor(kRed + 1);
  hSystZt[0]->SetFillColor(kRed + 1);
  hSystZt[0]->SetLineColor(kRed + 1);
  hSystZt[0]->SetLineWidth(2);
  hZtCent[0]->SetMarkerColor(kRed + 1);
  hZtCent[0]->SetLineColor(kRed + 1);
  hZtCent[0]->SetMarkerSize(1.2);
  hSystZt[0]->SetFillStyle(0);
  hSystZt[0]->Draw("samee2");
  hZtCent[0]->Draw("EP X0 same");
  TLatex *latALICEcms_Zhad = LatexStdISORatio(latALICEcms_Zhad, 0.380, 0.92, 0.045, cenBins[0], cenBins[1], ptMin, ptMax, true);
  TLegend *legALICE2_Zhad = LegStd(legALICE2_Zhad, 0.140, 0.16, 0.40, 0.40);
  legALICE2_Zhad->SetTextSize(0.045);
  legALICE2_Zhad->AddEntry(hZtCent[0], "ALICE, stat. unc. ", "ep");
  legALICE2_Zhad->AddEntry(hzTCMS, "CMS, #bf{#it{#gamma}#font[122]{-}jet}, stat. unc. ", "ep");
  legALICE2_Zhad->AddEntry(hzTCMS_Zhad, "CMS, #bf{#it{Z}#font[122]{-}hadron}, stat. unc. ", "ep");
  legALICE2_Zhad->AddEntry(hSystZt[0], " syst. unc. ", "f");
  legALICE2_Zhad->Draw("same");
  cCMS_Zhad->cd(2);
  cCMS_Zhad->cd(2)->SetRightMargin(0);
  TLegend *legCMS = LegStd(legCMS, 0, 0.68, 0.10, 0.98);
  legCMS->SetHeader("CMS, Phys.Rev.Lett. 121 (2018) 24, 242301, 2018");
  legCMS->AddEntry((TObject *)0, "#bf{#it{#gamma}#font[122]{-}jet}, #bf{0#font[122]{-}10%}", "");
  legCMS->AddEntry((TObject *)0, "anti-k_{T} jet R = 0.3, #it{p}_{T}^{ jet} > 30 GeV/#it{c}, |#it{#eta}^{jet}| < 1.6", "");
  legCMS->AddEntry((TObject *)0, "|#Delta#it{#varphi}_{#it{#gamma}#font[122]{-}jet}| > #frac{7}{8} #it{#pi}, |#it{#eta}^{ #it{#gamma}}| < 1.44, #it{p}_{T}^{ #it{#gamma}} > 60 GeV/#it{c} #otimes #it{p}_{T}^{ h} > 1 GeV/#it{c}", "");

  legCMS->Draw("same");
  TLegend *legCMSzt_Zhad = LegStd(legCMSzt_Zhad, 0, 0.42, 0.10, 0.62);
  legCMSzt_Zhad->SetHeader("CMS, Phys.Rev.Lett. 128 (2022) 12, 122301, 2022");
  legCMSzt_Zhad->AddEntry((TObject *)0, "#bf{#it{Z}#font[122]{-}hadron}, #bf{0#font[122]{-}30%} ", "");
  legCMSzt_Zhad->AddEntry((TObject *)0, "|#Delta#it{#varphi}_{#it{Z}#font[122]{-}h}| > #frac{7}{8} #it{#pi}, #it{p}_{T}^{ #it{Z}} > 30 GeV/#it{c} #otimes #it{p}_{T}^{ h} > 1 GeV/#it{c}", "");

  legCMSzt_Zhad->Draw("same");

  cCMS_Zhad->Print(dirPlot + Form("/Z_hadALICECMS%s.pdf", sPtAll.Data()));

  // Iaa CMS - Z-hadron
  TCanvas *cIaaCMS_Zhad = new TCanvas("cIaaCMS_Zhad", "cIaaCMS_Zhad", 2 * 800, 1 * 600);
  cIaaCMS_Zhad->Divide(2, 1);
  cIaaCMS_Zhad->cd(1);
  cIaaCMS_Zhad->cd(1)->SetTopMargin(0.015);
  cIaaCMS_Zhad->cd(1)->SetRightMargin(0.02);
  cIaaCMS_Zhad->cd(1)->SetLeftMargin(0.12);
  cIaaCMS_Zhad->cd(1)->SetBottomMargin(0.11);
  hIaaCMSSys_Zhad->SetFillStyle(0);
  hIaaCMSSys_Zhad->SetLineWidth(2);
  hGeneralIaa->SetTitle(" ");
  hGeneralIaa->GetYaxis()->SetRangeUser(0, 3);
  hGeneralIaa->GetXaxis()->SetRangeUser(-0.0, 1.1);
  hGeneralIaa->GetYaxis()->SetTitleSize(0.052);
  hGeneralIaa->GetYaxis()->SetTitleOffset(1.25);
  hGeneralIaa->GetXaxis()->SetTitleSize(0.045);
  hGeneralIaa->GetYaxis()->SetLabelSize(0.04);
  hGeneralIaa->GetXaxis()->SetLabelSize(0.04);
  hGeneralIaa->SetLineWidth(0);
  hGeneralIaa->Draw("hist");
  hIaaCMSSys->SetLineWidth(2);
  hIaaCMSSys->SetFillStyle(0);
  hIaaCMSSys->Draw("samee2");
  hIaaCMS->Draw("EP X0 same");
  hIaaCMSSys_Zhad->Draw("samee2");
  hIaaCMS_Zhad->Draw("EP X0 same");
  hPbPb_NLO[0]->SetLineColor(kRed + 1);
  hPbPb_NLO[0]->SetMarkerColor(kRed + 1);
  hPbPb_NLO[0]->SetMarkerSize(1.3);
  hPbPb_NLO[0]->SetLineWidth(2);
  hPbPb_NLOSyst[0]->SetLineWidth(2);
  hPbPb_NLOSyst[0]->SetFillColor(kRed + 1);
  hPbPb_NLOSyst[0]->SetLineColor(kRed + 1);
  hPbPb_NLOSyst[0]->SetMarkerSize(1.3);
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
  TLatex *latALICEcms2 = LatexStdISORatio(latALICEcms2, 0.380, 0.92, 0.045, cenBins[0], cenBins[1], ptMin, ptMax, true);
  TLegend *legALICEIaa1 = LegStd(legALICEIaa1, 0.50, 0.44, 0.70, 0.68);
  legALICEIaa1->SetTextSize(0.045);
  // legALICEIaa1->AddEntry(h3[0], "ALICE, stat. unc.", "ep");
  legALICEIaa1->AddEntry(hPbPb_NLO[0], "ALICE, stat. unc.", "ep");
  legALICEIaa1->AddEntry(hIaaCMS, "CMS, #bf{#it{#gamma}#font[122]{-}jet}, stat. unc. ", "ep");
  legALICEIaa1->AddEntry(hIaaCMS_Zhad, "CMS, #bf{#it{Z}#font[122]{-}hadron}, stat. unc. ", "ep");
  // legALICEIaa1->AddEntry(hsyst[0], " syst. unc.", "f");
  legALICEIaa1->AddEntry(hPbPb_NLOSyst[0], " syst. unc.", "f");
  legALICEIaa1->Draw("same");
  TGraph *lineCMS = DrawLine(lineCMS, -0.05, 0.5, 1.1, 0.5);
  lineCMS->Draw("l");
  TGraph *lineCMS1 = DrawLine(lineCMS1, -0.05, 1, 1.1, 1);
  lineCMS1->Draw("l");
  cIaaCMS_Zhad->cd(2);
  TLegend *legCMS1 = LegStd(legCMS1, 0, 0.68, 0.10, 0.98);
  legCMS1->SetHeader("CMS, Phys.Rev.Lett. 121 (2018) 242301, 2018");
  legCMS1->AddEntry((TObject *)0, "#bf{#it{#gamma}#font[122]{-}jet}, #bf{0#font[122]{-}10%} ", "");
  legCMS1->AddEntry((TObject *)0, "anti-k_{T} jet R = 0.3, #it{p}_{T}^{ jet} > 30 GeV/#it{c}, |#it{#eta}^{jet}| < 1.6 ", "");
  legCMS1->AddEntry((TObject *)0, "|#Delta#it{#varphi}_{#it{#gamma}#font[122]{-}jet}| > #frac{7}{8} #it{#pi}, |#it{#eta}^{ #it{#gamma}}| < 1.44 #it{p}_{T}^{ #it{#gamma}} > 60 GeV/#it{c} #otimes #it{p}_{T}^{ h} > 1 GeV/#it{c} ", "");
  legCMS1->Draw("same");
  TLegend *legCMS1_Zhad = LegStd(legCMS1_Zhad, 0, 0.42, 0.10, 0.62);
  legCMS1_Zhad->SetHeader("CMS, Phys.Rev.Lett. 128 (2022) 122301, 2022 ");
  legCMS1_Zhad->AddEntry((TObject *)0, "#bf{#it{Z}#font[122]{-}hadron}, #bf{0#font[122]{-}30%} ", "");
  legCMS1_Zhad->AddEntry((TObject *)0, "|#Delta#it{#varphi}_{#it{Z}#font[122]{-}h}| > #frac{7}{8} #it{#pi}, #it{p}_{T}^{ #it{Z}} > 30 GeV/#it{c} #otimes #it{p}_{T}^{ h} > 1 GeV/#it{c} ", "");
  legCMS1_Zhad->Draw("same");
  cIaaCMS_Zhad->Print(dirPlot + Form("/ppNLOpQCDIaa_ZhadALICECMS%s.pdf", sPtAll.Data()));
}

/*void PlotStyle(TH1F *hPlot, int kMarker, double kMarkerSize, int kColor, TString titleX, TString titleY)
{
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(1);
  gStyle->SetOptFit(1111111);
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleAlign(23);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadLeftMargin(0.15);
  // gStyle->SetTickX();
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  hPlot->SetMarkerStyle(kMarker);
  hPlot->SetMarkerSize(kMarkerSize);
  hPlot->SetMarkerColor(kColor);
  hPlot->SetLineColor(kColor);
  hPlot->SetFillColorAlpha(kColor, 0.250);
  hPlot->GetYaxis()->SetTitle(Form("%s", titleY.Data()));
  hPlot->GetXaxis()->SetTitle(Form("%s", titleX.Data()));
  hPlot->GetXaxis()->SetTickLength(0.015);
  hPlot->GetYaxis()->SetTickLength(0.02);

  // leg->SetFillColor(kWhite);
  // leg->SetLineColor(0);
}

TLatex *LatexStdISO(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax, bool bCen)
{
  lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetTextSize(0.04);
  lat->SetNDC();
  lat->DrawLatex(xpos, ypos, Form("#bf{ALICE preliminary}"));
  if (bCen)
    lat->DrawLatex(xpos, ypos - 0.06, Form("#bf{%d#font[122]{-}%d %%} Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, |#it{#eta}^{ #it{#gamma}}| < 0.67", cenMin, cenMax));
  else if (!bCen)
    lat->DrawLatex(xpos, ypos - 0.06, Form("Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, |#it{#eta}^{ #it{#gamma}}| < 0.67"));
  lat->DrawLatex(xpos, ypos - 2 * 0.06, Form("%2.0f < #it{p}_{T}^{ #it{#gamma}} < %2.0f GeV/#it{c} #otimes #it{p}_{T}^{ h} > 0.5 GeV/#it{c}", ptMin, ptMax));
  return lat;
}

TLegend *LegStd(TLegend *leg, double xpos1, double ypos1, double xpos2, double ypos2)
{
  leg = new TLegend(xpos1, ypos1, xpos2, ypos2);
  leg->SetFillColor(kWhite);
  leg->SetLineWidth(0);
  leg->SetTextSize(0.04);
  return leg;
}
*/