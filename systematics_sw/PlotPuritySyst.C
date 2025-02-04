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
#include "TF1.h"
#include <vector>
#include <algorithm>

using std::cout;
using std::endl;

Int_t nPtTrig = 9;
Float_t PtTrigger[] = {12, 14, 16, 18, 20, 25, 30, 40, 50, 60};
int nZtBin = 6;
double assocZt[] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 1.00};

// Define purity ranges directory file names
const int nPurUsed = 3;
// TString PurUsed[] = {"Results", "SystematicsPur11", "SystematicsPur09"};
TString PurUsed[] = {"", "SystPur11", "SystPur09"};

TFile *fPlot[nPurUsed];

Int_t kColor[] = {kAzure + 2, kPink + 6, kViolet + 6};
Int_t kStyle[] = {20, 71, 72};
Int_t kColorCen[] = {kAzure + 1, kAzure + 2, kAzure + 4, kAzure + 3};
Int_t kStyleCen[] = {21, 71, 20, 73};

void PlotStyle(TH1F *hPlot, int kMarker, double kMarkerSize, int kColor, TString titleX, TString titleY);
TLatex *LatexStd(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax, Bool_t bCen);
TLegend *LegStd(TLegend *leg, double xpos1, double ypos1, double xpos2, double ypos2);
TLatex *LatexStdISORatio(TLatex *lat, double xpos, double ypos, double texSize, int cenMin, int cenMax, float ptMin, float ptMax, bool bCen);

void PlotPuritySyst(Float_t ptMin = 18, Float_t ptMax = 40, TString sMixed = "Mixed", TString sDirRefFiles = "~/work/histogram/FromScratch/checkCode", bool Mirror = true,TString shshBkg = "0.40-1.00", bool b0_30 = true)
{

  Int_t nCen;
  std::vector<Int_t> cenBins;
  if (b0_30)
  {
    nCen = 3;
    cenBins.push_back(0);
    cenBins.push_back(30);
    cenBins.push_back(50);
    cenBins.push_back(90);
  }
  else if (!b0_30)
  {
    nCen = 4;
    cenBins.push_back(0);
    cenBins.push_back(10);
    cenBins.push_back(30);
    cenBins.push_back(50);
    cenBins.push_back(90);
  }

  TH1F *hZt[nPurUsed][nCen];
  TH1F *hPur[nPurUsed][nCen];
  TH1F *hZtUncertSyst[nCen];

  TH1F *hIcp[nPurUsed][nCen];
  TH1F *hIcpSyst[nCen];
  TH1F *hIcpUncertSyst[nCen];

  TString mirror;
  if (Mirror)
    mirror = "Mirror";
  else
    mirror = "NoMirror";

  TString sPtAll = Form("_Pt%2.0f_%2.0f", ptMin, ptMax);
  TString dirPlot; // define directory to save files

  if (b0_30)
  {
    dirPlot = Form("~/work/histogram/Systematics_checkCode0_30/Purity");
    gSystem->Exec(Form("mkdir %s", dirPlot.Data()));
  }
  if (!b0_30)
  {
    dirPlot = Form("~/work/histogram/Systematics_checkCode/Purity");
    gSystem->Exec(Form("mkdir %s", dirPlot.Data()));
  }

  TFile *fPurSyst = new TFile(Form("%s/fPurSyst%s%s%s.root", dirPlot.Data(), sMixed.Data(), shshBkg.Data(), sPtAll.Data()), "RECREATE");

  TLegend *legPurData[nCen];
  TCanvas *cPur[nCen];
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////// Get the different zT functions obtained with various purity ////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    for (int iPur = 0; iPur < nPurUsed; iPur++)
    {
      fPlot[iPur] = new TFile(Form("%s%s/fPlot%s_%s%s.root", sDirRefFiles.Data(), PurUsed[iPur].Data(), shshBkg.Data(), sCent.Data(), sPtAll.Data()));
      hZt[iPur][iCen] = (TH1F *)fPlot[iPur]->Get(Form("hZtEffCorrIso1Photon_%s%s", sCent.Data(), sPtAll.Data()));
      cout << hZt[iPur][iCen] << endl;
    }
  }
  ///////////////////////////////////////////////////////////////////////////
  ////////////////// Estimate Icp for various purity ///////////////////////
  /////////////////////////////////////////////////////////////////////////
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    for (int iPur = 0; iPur < nPurUsed; iPur++)
    {
      hIcp[iPur][iCen] = (TH1F *)hZt[iPur][iCen]->Clone(Form("hIcpIso1Photon_%s%s", sCent.Data(), sPtAll.Data()));
      if (b0_30)
        hIcp[iPur][iCen]->Divide(hZt[iPur][2]);
      else if (!b0_30)
        hIcp[iPur][iCen]->Divide(hZt[iPur][3]);
    }
  }

  //////////////////////////////////////////////////////////
  ////////////////// Plotting style ///////////////////////
  ////////////////////////////////////////////////////////

  TCanvas *cZt[nCen];
  TCanvas *cIcp[nCen];
  TLegend *legZTData[nCen];
  TLegend *legIcpData[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    TString PathPlot = {Form("%s/%s_%s%s", dirPlot.Data(), sCent.Data(), sMixed.Data(), mirror.Data())};
    gSystem->Exec(Form("mkdir %s", PathPlot.Data())); // define directory for all centralities

    cZt[iCen] = new TCanvas(Form("cZt%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("cZt%d_%d", cenBins[iCen], cenBins[iCen + 1]), 800, 600);
    cIcp[iCen] = new TCanvas(Form("cIcp%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("cIcp%d_%d", cenBins[iCen], cenBins[iCen + 1]), 800, 600);

    for (int iPur = 0; iPur < nPurUsed; iPur++)
    {
      PlotStyle(hZt[iPur][iCen], kStyle[iPur], 2, kColor[iPur], "#it{z}_{T}", "1/N^{trig}dN^{charg}/d#font[12]{z}_{T}");
      hZt[iPur][iCen]->SetDirectory(0);
      // gStyle->SetPadLeftMargin(0.24);
      hZt[iPur][iCen]->GetXaxis()->SetTitleSize(0.045);
      hZt[iPur][iCen]->GetYaxis()->SetTitleSize(0.045);
      hZt[iPur][iCen]->GetYaxis()->SetTitleOffset(0.72);
      hZt[iPur][iCen]->GetXaxis()->SetLabelSize(0.04);
      hZt[iPur][iCen]->GetYaxis()->SetLabelSize(0.04);
      cZt[iCen]->cd();
      gPad->SetLogy();
      hZt[iPur][iCen]->Draw("same");

      cIcp[iCen]->cd();
      PlotStyle(hIcp[iPur][iCen], kStyle[iPur], 2, kColor[iPur], "#it{z}_{T}", "#it{I}_{CP}");
      hIcp[iPur][iCen]->SetDirectory(0);
      hIcp[iPur][iCen]->GetYaxis()->SetRangeUser(0, 3);
      hIcp[iPur][iCen]->Draw("same");
    }
    legZTData[iCen] = LegStd(legZTData[iCen], 0.60, 0.50, 0.85, 0.70);
    legZTData[iCen]->AddEntry(hZt[0][iCen], "Nominal");
    legZTData[iCen]->AddEntry(hZt[1][iCen], "Nominal + #sigma_{syst}^{#it{P}}");
    legZTData[iCen]->AddEntry(hZt[2][iCen], "Nominal - #sigma_{syst}^{#it{P}}");
    cZt[iCen]->cd();
    TLatex *ALICEtex1 = LatexStd(ALICEtex1, 0.500, 0.84, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true);
    legZTData[iCen]->Draw("same");
    cZt[iCen]->Print(Form("%s/hZtPuritySystCen%d_%d.pdf", PathPlot.Data(), cenBins[iCen], cenBins[iCen + 1]));

    legIcpData[iCen] = LegStd(legIcpData[iCen], 0.120, 0.50, 0.55, 0.65);
    legIcpData[iCen]->AddEntry(hIcp[0][iCen], "Nominal");
    legIcpData[iCen]->AddEntry(hIcp[1][iCen], "Nominal + #sigma_{syst}^{#it{P}}");
    legIcpData[iCen]->AddEntry(hIcp[2][iCen], "Nominal - #sigma_{syst}^{#it{P}}");
    cIcp[iCen]->cd();
    legIcpData[iCen]->Draw("same");
    TLatex *ALICEtex10 = new TLatex();
    ALICEtex10->SetTextFont(42);
    ALICEtex10->SetTextSize(0.04);
    ALICEtex10->SetNDC();
    ALICEtex10->DrawLatex(0.15, 0.85, Form("#it{This Thesis}"));
    ALICEtex10->DrawLatex(0.15, 0.85 - 0.06, Form("#bf{%d-%d %% / 50-90 %%} Pb-Pb, #sqrt{s_{NN}} = 5.02 TeV", cenBins[iCen], cenBins[iCen + 1]));
    ALICEtex10->DrawLatex(0.15, 0.85 - 2 * 0.06, Form("|#it{#eta}^{ #it{#gamma}}| < 0.67, %2.0f < #it{p}_{#it{T}}^{#gamma} < %2.0f GeV/#it{c}", ptMin, ptMax));
    TLatex *ALICEtex1Icp = LatexStd(ALICEtex1Icp, 0.500, 0.84, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true);
    cIcp[iCen]->Print(Form("%s/hIcpPuritySystCen%d_%d.pdf", PathPlot.Data(), cenBins[iCen], cenBins[iCen + 1]));
  }
  ///////////////////////////////////////////////////////////////////////////
  ///////////// Estimate systematics (|fmax-fmin|/2)fref ///////////////////
  /////////////////////////////////////////////////////////////////////////
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    cout << "Systematics: " << endl;
    cout << hZt[0][iCen]->GetNbinsX() << endl;
    hZtUncertSyst[iCen] = new TH1F(Form("hZtUncertSystCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("hZtUncertSystCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), nZtBin, assocZt);
    hIcpUncertSyst[iCen] = new TH1F(Form("hIcpUncertSystCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("hIcpUncertSystCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), nZtBin, assocZt);

    for (int ibin = 0; ibin < (hZt[0][iCen]->GetNbinsX()); ibin++)
    {
      ////////////////////////// Debugging print //////////////////////////

      // cout << "Zt bin [" << assocZt[ibin] << "-" << assocZt[ibin + 1] << "] : ";
      // cout << "Pur 1: " << hZt[0][iCen]->GetBinContent(ibin + 1) << ", Error: " << hZt[0][iCen]->GetBinError(ibin + 1) << endl;
      // cout << "Pur + 10%%: " << hZt[1][iCen]->GetBinContent(ibin + 1) << ", Error: " << hZt[1][iCen]->GetBinError(ibin + 1) << endl;
      // cout << "Pur - 10%%: " << hZt[2][iCen]->GetBinContent(ibin + 1) << ", Error: " << hZt[2][iCen]->GetBinError(ibin + 1) << endl;
      // cout << ((hZt[1][iCen]->GetBinContent(ibin + 1) - hZt[2][iCen]->GetBinContent(ibin + 1)) / 2) / hZt[0][iCen]->GetBinContent(ibin + 1) << endl;

      hZtUncertSyst[iCen]->SetBinContent(ibin + 1, abs(((hZt[1][iCen]->GetBinContent(ibin + 1) - hZt[2][iCen]->GetBinContent(ibin + 1)) / 2) / hZt[0][iCen]->GetBinContent(ibin + 1)));
      hZtUncertSyst[iCen]->SetBinError(ibin + 1, hZt[0][iCen]->GetBinError(ibin + 1) / hZt[0][iCen]->GetBinContent(ibin + 1));
      hIcpUncertSyst[iCen]->SetBinContent(ibin + 1, abs(((hIcp[1][iCen]->GetBinContent(ibin + 1) - hIcp[2][iCen]->GetBinContent(ibin + 1)) / 2) / hIcp[0][iCen]->GetBinContent(ibin + 1)));
    }
  }
  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////// Compute the fit for polishing the trend for zT:   ///////////////
  ////////////////// pol0 for 0-10% and expo for 0-30%, 30-50%, 50-90%  //////////////
  ////////////////////////////////////////////////////////////////////////////////////
  TCanvas *cPurSyst[nCen];
  TH1F *hPurUncertFromFit[nCen];
  TLatex *parPol0[nCen];
  TF1 *fa0[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    TString PathPlot = {Form("%s/%s_%s%s", dirPlot.Data(), sCent.Data(), sMixed.Data(), mirror.Data())};
    hPurUncertFromFit[iCen] = new TH1F(Form("hPurUncertFromFit%s", sCent.Data()), Form("hPurUncertFromFit%s", sCent.Data()), nZtBin, assocZt);
    hZtUncertSyst[iCen]->SetDirectory(0);
    cPurSyst[iCen] = new TCanvas(Form("cZtPurSys%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("cZtPurSys%d_%d", cenBins[iCen], cenBins[iCen + 1]), 800, 600);
    PlotStyle(hZtUncertSyst[iCen], kStyleCen[iCen], 2, kColorCen[iCen], "#font[12]{z}_{T}", "Uncertainty %");
    cPurSyst[iCen]->cd();
    hZtUncertSyst[iCen]->Scale(100);
    hZtUncertSyst[iCen]->Draw("plhist");
    if (iCen == 0 && !b0_30) // 0-10%
    {
      fa0[iCen] = new TF1(Form("fitpol0Purity%d_%d", cenBins[iCen], cenBins[iCen + 1]), "pol0", 0.30, 1.00);
    }
    else
    {
      fa0[iCen] = new TF1(Form("fitpol0Purity%d_%d", cenBins[iCen], cenBins[iCen + 1]), "expo", 0.10, 0.8); // 0-30%, 30-50%, 50-90%
    }
    // gStyle->SetOptFit(1111);
    hZtUncertSyst[iCen]->GetYaxis()->SetRangeUser(0, 70);
    hZtUncertSyst[iCen]->Fit(Form("fitpol0Purity%d_%d", cenBins[iCen], cenBins[iCen + 1]), "R");
    parPol0[iCen] = new TLatex();
    parPol0[iCen]->SetTextSize(0.04);
    parPol0[iCen]->SetTextFont(42);
    parPol0[iCen]->SetNDC();
    fa0[iCen]->SetLineColor(kColor[iCen]);
    fa0[iCen]->SetLineStyle(10);
    fa0[iCen]->SetLineWidth(5);
    fa0[iCen]->Draw("same");
    parPol0[iCen]->DrawLatex(0.50, 0.8, Form("#chi^{2}/NDF: %f/%d", fa0[iCen]->GetChisquare(), fa0[iCen]->GetNDF()));
    parPol0[iCen]->DrawLatex(0.50, 0.75, Form("par0 = %f#pm%f", fa0[iCen]->GetParameter(0), fa0[iCen]->GetParError(0)));
    TLatex *ALICEtex3 = LatexStd(ALICEtex3, 0.120, 0.84, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true);
    cPurSyst[iCen]->Print(Form("%s/hUncertPuritySyst.pdf", PathPlot.Data()));
    for (int ibin = 0; ibin < nZtBin; ibin++)
    {
      hPurUncertFromFit[iCen]->SetBinContent(ibin + 1, fa0[iCen]->Eval(hPurUncertFromFit[iCen]->GetBinCenter(ibin + 1)));
    }
    PlotStyle(hPurUncertFromFit[iCen], kStyleCen[iCen], 2, kColorCen[iCen], "#font[12]{z}_{T}", "Uncertainty %");
    // new TCanvas();
    // hPurUncertFromFit[iCen]->SetDirectory(0);
    // hPurUncertFromFit[iCen]->Draw("");
    fPurSyst->cd();
    hZtUncertSyst[iCen]->Write();
    fa0[iCen]->Write();
    hPurUncertFromFit[iCen]->Write();
  }
  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////// Compute the fit for polishing the trend for Icp:   ///////////////
  ////////////////// pol0 for 0-10% and expo for 0-30%, 30-50%, 50-90%  //////////////
  ////////////////////////////////////////////////////////////////////////////////////
  TCanvas *cPurSystIcp[nCen];
  TH1F *hPurUncertFromFitIcp[nCen];
  TLatex *parPol0Icp[nCen];
  TF1 *fa0Icp[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    TString PathPlot = {Form("%s/%s_%s%s", dirPlot.Data(), sCent.Data(), sMixed.Data(), mirror.Data())};
    hPurUncertFromFitIcp[iCen] = new TH1F(Form("hPurUncertFromFitIcp%s", sCent.Data()), Form("hPurUncertFromFitIcp%s", sCent.Data()), nZtBin, assocZt);
    hIcpUncertSyst[iCen]->SetDirectory(0);
    cPurSystIcp[iCen] = new TCanvas(Form("cIcpPurSys%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("cIcpPurSys%d_%d", cenBins[iCen], cenBins[iCen + 1]), 800, 600);
    PlotStyle(hIcpUncertSyst[iCen], kStyleCen[iCen], 2, kColorCen[iCen], "#font[12]{z}_{T}", "Uncertainty %");
    cPurSystIcp[iCen]->cd();
    hIcpUncertSyst[iCen]->Scale(100);
    hIcpUncertSyst[iCen]->Draw("p");
    if (iCen == 0 && !b0_30)
    {
      fa0Icp[iCen] = new TF1(Form("fitpol1IcpPurity%d_%d", cenBins[iCen], cenBins[iCen + 1]), "pol0", 0.30, 1.00);
    }
    else
    {
      fa0Icp[iCen] = new TF1(Form("fitpol1IcpPurity%d_%d", cenBins[iCen], cenBins[iCen + 1]), "expo", 0.10, 0.6);
    }
    // gStyle->SetOptFit(1111);
    hIcpUncertSyst[iCen]->Fit(Form("fitpol1IcpPurity%d_%d", cenBins[iCen], cenBins[iCen + 1]), "R");
    parPol0Icp[iCen] = new TLatex();
    parPol0Icp[iCen]->SetTextSize(0.04);
    parPol0Icp[iCen]->SetTextFont(42);
    parPol0Icp[iCen]->SetNDC();
    fa0Icp[iCen]->SetLineColor(kColor[iCen]);
    fa0Icp[iCen]->SetLineStyle(10);
    fa0Icp[iCen]->SetLineWidth(5);
    fa0Icp[iCen]->Draw("same");
    parPol0Icp[iCen]->DrawLatex(0.60, 0.8, Form("#chi^{2}/NDF: %f/%d", fa0Icp[iCen]->GetChisquare(), fa0Icp[iCen]->GetNDF()));
    parPol0Icp[iCen]->DrawLatex(0.60, 0.75, Form("const = %f#pm%f", fa0Icp[iCen]->GetParameter(0), fa0Icp[iCen]->GetParError(0)));
    TLatex *ALICEtex3Icp = LatexStd(ALICEtex3Icp, 0.120, 0.84, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true);
    cPurSystIcp[iCen]->Print(Form("%s/hUncertPuritySystIcp.pdf", PathPlot.Data()));

    for (int ibin = 0; ibin < nZtBin; ibin++)
    {
      if (iCen == 0)
      {
        hPurUncertFromFitIcp[iCen]->SetBinContent(ibin + 1, fa0Icp[iCen]->Eval(hPurUncertFromFitIcp[iCen]->GetBinCenter(ibin + 1)));
      }
      else
      {
        hPurUncertFromFitIcp[iCen]->SetBinContent(ibin + 1, hIcpUncertSyst[iCen]->GetBinContent(ibin + 1));
      }
    }
    PlotStyle(hPurUncertFromFitIcp[iCen], kStyleCen[iCen], 2, kColorCen[iCen], "#font[12]{z}_{T}", "Uncertainty %");

    fPurSyst->cd();
    hIcpUncertSyst[iCen]->Write();
    fa0[iCen]->Write();
    hPurUncertFromFitIcp[iCen]->Write();
  }

  TLegend *legUncer = LegStd(legUncer, 0.580, 0.70, 0.7, 0.850);
  TLegend *legUncerFit = LegStd(legUncerFit, 0.70, 0.70, 0.88, 0.80);
  TCanvas *cPurSystAllCen = new TCanvas("cPurSystAllCen", "cPurSystAllCen", 800, 600);
  TCanvas *cPurSystAllCenNoFit = new TCanvas("cPurSystAllCenNoFit", "cPurSystAllCenNoFit", 800, 600);

  TLegend *legUncerIcp = LegStd(legUncerIcp, 0.55, 0.70, 0.85, 0.86);
  legUncerIcp->SetNColumns(2);
  TCanvas *cPurSystAllCenIcp = new TCanvas("cPurSystAllCenIcp", "cPurSystAllCenIcp", 800, 600);
  TCanvas *cPurSystAllCenNoFitIcp = new TCanvas("cPurSystAllCenNoFitIcp", "cPurSystAllCenNoFitIcp", 800, 600);
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    cPurSystAllCen->cd();
    hPurUncertFromFit[iCen]->GetYaxis()->SetRangeUser(0, 70);
    hPurUncertFromFit[iCen]->Draw("histsamepl");

    legUncer->AddEntry(hPurUncertFromFit[iCen], Form("%d-%d %%", cenBins[iCen], cenBins[iCen + 1]));

    cPurSystAllCenNoFit->cd();
    hZtUncertSyst[iCen]->GetYaxis()->SetRangeUser(0, 70);
    hZtUncertSyst[iCen]->Draw("histsamepl");
    legUncerFit->AddEntry(fa0[iCen], "", "l");
    fa0[iCen]->Draw("same");
  }
  cPurSystAllCen->cd();
  legUncer->Draw("same");

  TLatex *ALICEtex2 = LatexStd(ALICEtex2, 0.120, 0.84, cenBins[0], cenBins[1], ptMin, ptMax, false);
  cPurSystAllCenNoFit->cd();
  legUncer->Draw("same");
  legUncerFit->Draw("same");
  TLatex *ALICEtex4 = LatexStd(ALICEtex4, 0.120, 0.84, cenBins[0], cenBins[1], ptMin, ptMax, false);
  cPurSystAllCen->Print(Form("%s/UncertPuritySystAllCent%s.pdf", dirPlot.Data(), sPtAll.Data()));
  cPurSystAllCenNoFit->Print(Form("%s/UncertPuritySystAllCentNoFit%s.pdf", dirPlot.Data(), sPtAll.Data()));

  // Systematics Icp
  for (int iCen = 0; iCen < nCen - 1; iCen++)
  {

    cPurSystAllCenIcp->cd();
    hPurUncertFromFitIcp[iCen]->GetYaxis()->SetRangeUser(-0.1, 50);
    hPurUncertFromFitIcp[iCen]->Draw("histsamepl");

    legUncerIcp->AddEntry(hPurUncertFromFitIcp[iCen], Form("#bf{%d#font[122]{-}%d / 50#font[122]{-}90 %%}", cenBins[iCen], cenBins[iCen + 1]));
    legUncerIcp->AddEntry(fa0[iCen], "", "l");
    cPurSystAllCenNoFitIcp->cd();
    hIcpUncertSyst[iCen]->GetYaxis()->SetRangeUser(-0.1, 50);
    hIcpUncertSyst[iCen]->Draw("histsamepl");
    fa0Icp[iCen]->Draw("same");
  }

  cPurSystAllCenIcp->cd();
  legUncerIcp->Draw("same");
  TLatex *ALICEtex2Icp = LatexStd(ALICEtex2Icp, 0.120, 0.84, cenBins[0], cenBins[1], ptMin, ptMax, false);

  cPurSystAllCenNoFitIcp->cd();
  legUncerIcp->Draw("same");
  TLatex *ALICEtex4Icp = LatexStd(ALICEtex4Icp, 0.120, 0.84, cenBins[0], cenBins[1], ptMin, ptMax, false);

  cPurSystAllCenIcp->Print(Form("%s/UncertIcpPuritySystAllCent%s.pdf", dirPlot.Data(), sPtAll.Data()));
  cPurSystAllCenNoFitIcp->Print(Form("%s/UncertIcpPuritySystAllCentNoFit%s.pdf", dirPlot.Data(), sPtAll.Data()));

  fPurSyst->Close();
}

void PlotStyle(TH1F *hPlot, int kMarker, double kMarkerSize, int kColor, TString titleX, TString titleY)
{
  gStyle->SetTitleX(0.56);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetTickLength(0.02, "x");
  gStyle->SetTickLength(0.02, "y");
  // gStyle->SetOptFit(1111111);
  // gStyle->SetLabelSize(.5, "XY");
  gStyle->SetLineScalePS(1);
  gStyle->SetTitleAlign(23);
  // gStyle->SetPadRightMargin(0.03);
  // gStyle->SetPadLeftMargin(0.12);
  hPlot->SetMarkerStyle(kMarker);
  hPlot->SetMarkerSize(1.1);
  // hPlot->SetMarkerSize(2);
  hPlot->SetMarkerColor(kColor);
  hPlot->SetLineColor(kColor);
  hPlot->SetLineWidth(1);
  hPlot->GetYaxis()->SetTitle(Form("%s", titleY.Data()));
  hPlot->GetXaxis()->SetTitle(Form("%s", titleX.Data()));
  // hPlot->GetXaxis()->SetLabelSize(.05);
  // hPlot->GetYaxis()->SetLabelSize(.05);
  //  leg->SetFillColor(kWhite);
  //  leg->SetLineColor(0);
}
TLatex *LatexStd(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax, Bool_t bCen)
{
  lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetTextSize(0.04);
  lat->SetNDC();
  if (bCen)
  {
    lat->DrawLatex(xpos, ypos, Form("#it{This Thesis}"));
    lat->DrawLatex(xpos, ypos - 0.06, Form("%d#font[122]{-}%d %% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV", cenMin, cenMax));
    lat->DrawLatex(xpos, ypos - 2 * 0.06, Form("|#it{#eta}^{ #it{#gamma}}| < 0.67 , %2.0f < #it{p}_{T}^{#gamma} < %2.0f GeV/#it{c}", ptMin, ptMax));
  }
  else if (!bCen)
  {
    lat->DrawLatex(xpos, ypos, Form("#it{This Thesis}"));
    lat->DrawLatex(xpos, ypos - 0.06, Form("Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV"));
    lat->DrawLatex(xpos, ypos - 2 * 0.06, Form("|#it{#eta}^{ #it{#gamma}}| < 0.67 , %2.0f < #it{p}_{T}^{#gamma} < %2.0f GeV/#it{c}", ptMin, ptMax));
  }
  return lat;
}

TLegend *LegStd(TLegend *leg, double xpos1, double ypos1, double xpos2, double ypos2)
{
  leg = new TLegend(xpos1, ypos1, xpos2, ypos2);
  leg->SetFillColor(kWhite);
  leg->SetLineWidth(0);
  leg->SetTextSize(0.03);
  return leg;
}
TLatex *LatexStdISORatio(TLatex *lat, double xpos, double ypos, double texSize, int cenMin, int cenMax, float ptMin, float ptMax, bool bCen)
{
  lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetTextSize(texSize);
  lat->SetNDC();
  // lat->DrawLatex(xpos, ypos, Form("#font[42]{ALICE preliminary}"));
  lat->DrawLatex(xpos, ypos, Form("#it{This Thesis}"));
  if (bCen)
    lat->DrawLatex(xpos, ypos - 0.06, Form("#bf{%d#font[122]{-}%d%%} Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV ", cenMin, cenMax));
  else if (!bCen)
    lat->DrawLatex(xpos, ypos - 0.064, Form("Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV "));
  lat->DrawLatex(xpos, ypos - 2 * 0.064, Form("|#Delta#it{#varphi}_{#it{#gamma}#font[122]{-}h}| > #frac{3}{5} #it{#pi}, |#it{#eta}^{ #it{#gamma}}| < 0.67 "));
  lat->DrawLatex(xpos, ypos - 3 * 0.064, Form("%2.0f < #it{p}_{T}^{ #it{#gamma}} < %2.0f GeV/#it{c} #otimes #it{p}_{T}^{ h} > 0.5 GeV/#it{c} ", ptMin, ptMax));
  return lat;
}