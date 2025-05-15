
////////////////////////////////////////////////
///////////  PlotAllSystem.c  //////////////////
////////////////////////////////////////////////
//// Created by Carolina ARATA on 11/01/2022. //
////////////////////////////////////////////////

#include <TCanvas.h>
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
#include "TLegend.h"
#include <vector>
#include <algorithm>
#include <iterator>
#include "../Plotting.h"

const int nZtBin = 6;
double assocZt[] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 1.00};

int nZtBinThin = 9;
double assocZtThinner[] = {0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80, 1.00, 1.05};

const Int_t fillColors[] = {kGray + 1, kRed - 10, kBlue - 9, kGreen - 8, kMagenta - 9, kOrange - 9, kCyan - 8, kYellow - 7}; // for syst bands
const Int_t colors[] = {kBlack, kRed, kOrange + 7, kGreen + 3, kAzure + 1, kBlue + 1, kMagenta + 1};
const Int_t markers[] = {kFullCircle, kFullSquare, kFullCross, 24, kFullStar, 29, 47};

TFile *fNMixSyst;
void PrintSyst(TH1F *hSyst);
void PlotAllSystem(Float_t ptMin = 18, Float_t ptMax = 40, bool bMirror = true, TString Mixed = "Mixed", TString shshBkg = "0.40-1.00", TString dirRef = "Output_checkCode", bool b0_30 = true)
{

  Int_t nCen;
  std::vector<Int_t> cenBins;
  TString dirPlot; // Define directory where the root files and plot will be saved
  if (b0_30)
  {
    nCen = 3;
    // Int_t cenBins[] = {0, 10, 30, 50, 90};
    cenBins.push_back(0);
    cenBins.push_back(30);
    cenBins.push_back(50);
    cenBins.push_back(90);
    dirPlot = "Systematics_checkCode0_30";
  }
  else if (!b0_30)
  {
    nCen = 4;
    cenBins.push_back(0);
    cenBins.push_back(10);
    cenBins.push_back(30);
    cenBins.push_back(50);
    cenBins.push_back(90);
    dirPlot = "Systematics_checkCode";
  }

	TString processline = Form(".! mkdir -pv %s",dirPlot.Data()) ;
  gROOT->ProcessLine(processline.Data());
  
  TString sMirror = " ";
  if (bMirror)
    sMirror = "Mirror";

  TString sPtAll = Form("_Pt%2.0f_%2.0f", ptMin, ptMax); // Pt Range

  // Define output directory to save root files
  TFile *fSystFile = new TFile(Form("%s/fAllSystFile%s%s%s%s.root", dirPlot.Data(), Mixed.Data(), shshBkg.Data(), sMirror.Data(), sPtAll.Data()), "RECREATE");

  // Get zT distributions
  cout << "Get zT distributions" << endl;
  TFile *fResultsZt[nCen];
  TFile *fileNLO = new TFile("RootFiles/fileNLO.root"); // directory for theory

  // Get all the systematics from the other files
  cout << "Define all the directories to get the systematic uncertainties" << endl;
  // Purity
  TFile *fPurSyst = new TFile(Form("%s/Purity/fPurSyst%s%s%s.root", dirPlot.Data(), Mixed.Data(), shshBkg.Data(), sPtAll.Data()));
  cout << "Purity syst: " << fPurSyst << endl;
  // TFile *fUESyst = new TFile(Form("%s/fUESyst%s%s%s.root", dirPlot.Data(), Mixed.Data(), shshBkg.Data(), sPtAll.Data()));

  // ShSh
  TFile *fShShSyst = new TFile(Form("%s/ShShSyst/fShSyst%s%s%s.root", dirPlot.Data(), Mixed.Data(), shshBkg.Data(), sPtAll.Data()));
  cout << "ShSh syst: " << fShShSyst << endl;
  // Tracking efficiency
  TFile *fTrackEffSyst = new TFile(Form("%s/TrackIneff/fSystTrackIneff%s.root", dirPlot.Data(), sPtAll.Data()));
  // N centrality bins used for mixing
  TFile *fNMixSyst = new TFile(Form("%s/SystNCentXMix/fNMixCentSyst%s.root", dirPlot.Data(), sPtAll.Data()));
  // Residual UE
  TFile *fUEResidSyst = new TFile(Form("Systematics_checkCode/SystResidualUE/fUEResidSyst%s%s.root", shshBkg.Data(), sPtAll.Data()));
  cout << "Residual UE syst: " << fUEResidSyst << endl;

  // Define histograms
  double rangeMin[3] = {0.1, 0.1, 0.1};
  double rangeMax[3] = {1.00, 0.6, 0.6};
  TH1F *hUEUncert[nCen];
  TH1F *hPurUncert[nCen];
  TH1F *hCentMatchMixUncert[nCen];
  TH1F *hTrackIneffUncer[nCen];
  TH1F *hShShUncert[nCen];
  TH1F *hNMixCentUncert[nCen];
  TH1F *hUEresidUncert[nCen];

  TH1F *fZt[nCen];
  TH1F *fIcp[nCen];
  TH1F *hZtStatUncert[nCen];
  TH1F *hZtSystSumQuadr[nCen];

  gSystem->Exec(Form("mkdir %s/SystSh%s", dirPlot.Data(), shshBkg.Data()));

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////// Getter of all the systematics distribution, saved in previous declared histograms ///////////////////////////////////
  ///////// For the residual UE systematic, the files are all in ~/work/histogram/Systematics_checkCode/SystResidualUE/ ////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if (b0_30)
  {
    hUEresidUncert[0] = (TH1F *)fUEResidSyst->Get(Form("histSystErrFitNoUENarr0_30PtRangeFitTrend_Pt18_40"));
    //hUEresidUncert[1] = (TH1F *)fUEResidSyst->Get(Form("hErrSystNarrNoUE_Cen30_50_Pt18_40"));
    //hUEresidUncert[2] = (TH1F *)fUEResidSyst->Get(Form("hErrSystNarrNoUE_Cen50_90_Pt18_40"));
    hUEresidUncert[1] = (TH1F *)fUEResidSyst->Get(Form("hErrSystNarrNoUEFitTrend_Cen30_50_Pt18_40"));
    hUEresidUncert[2] = (TH1F *)fUEResidSyst->Get(Form("hErrSystNarrNoUEFitTrend_Cen50_90_Pt18_40"));
    
  }
  else if (!b0_30)
  {
    //hUEresidUncert[0] = (TH1F *)fUEResidSyst->Get(Form("hErrSystNarrNoUE_Cen0_10_Pt18_40"));
    //hUEresidUncert[1] = (TH1F *)fUEResidSyst->Get(Form("hErrSystNarrNoUE_Cen10_30_Pt18_40"));
    //hUEresidUncert[2] = (TH1F *)fUEResidSyst->Get(Form("hErrSystNarrNoUE_Cen30_50_Pt18_40"));
    //hUEresidUncert[3] = (TH1F *)fUEResidSyst->Get(Form("hErrSystNarrNoUE_Cen50_90_Pt18_40"));
    hUEresidUncert[0] = (TH1F *)fUEResidSyst->Get(Form("hErrSystNarrNoUEFitTrend_Cen0_10_Pt18_40"));
    hUEresidUncert[1] = (TH1F *)fUEResidSyst->Get(Form("hErrSystNarrNoUEFitTrend_Cen10_30_Pt18_40"));
    hUEresidUncert[2] = (TH1F *)fUEResidSyst->Get(Form("hErrSystNarrNoUEFitTrend_Cen30_50_Pt18_40"));
    hUEresidUncert[3] = (TH1F *)fUEResidSyst->Get(Form("hErrSystNarrNoUEFitTrend_Cen50_90_Pt18_40"));
  }

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    cout << sCent << endl;
    fResultsZt[iCen] = new TFile(Form("%s/fPlot%s_%s%s.root", dirRef.Data(), shshBkg.Data(), sCent.Data(), sPtAll.Data()));
    cout << fResultsZt[iCen] << endl;
    fZt[iCen] = (TH1F *)fResultsZt[iCen]->Get(Form("hZtEffCorrIso1Photon_%s%s", sCent.Data(), sPtAll.Data()));
    cout << fZt[iCen] << endl;
    fIcp[iCen] = (TH1F *)fZt[iCen]->Clone(Form("hIcp_%s", sCent.Data()));
    cout << fZt[iCen] << endl;
    // hUEUncert[iCen] = (TH1F *)fUESyst->Get(Form("hRelUncertUESumShCen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    // cout << hUEUncert[iCen] << endl;
    hPurUncert[iCen] = (TH1F *)fPurSyst->Get(Form("hPurUncertFromFitCen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    hShShUncert[iCen] = (TH1F *)fShShSyst->Get(Form("hShShUncertFromFitCen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    hTrackIneffUncer[iCen] = (TH1F *)fTrackEffSyst->Get(Form("hTrackIneffUncerFromFitCen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    cout << hPurUncert[iCen] << endl;
    cout << hShShUncert[iCen] << endl;
    cout << hTrackIneffUncer[iCen] << endl;
    hNMixCentUncert[iCen] = (TH1F *)fNMixSyst->Get(Form("hFromFitUncertNMix_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    cout << hNMixCentUncert[iCen] << endl;
    hZtStatUncert[iCen] = new TH1F(Form("hZTStatistCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("hZTStatistCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), nZtBin, assocZt);
    hZtSystSumQuadr[iCen] = new TH1F(Form("hZtSystSumQuadrCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("hhZtSystSumQuadrCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), nZtBin, assocZt);

    cout << "Estimate statistical uncertainty" << endl;
    for (int ibin = 0; ibin < nZtBin; ibin++)
    {
      hZtStatUncert[iCen]->SetBinContent(ibin + 1, fZt[iCen]->GetBinError(ibin + 1) / fZt[iCen]->GetBinContent(ibin + 1));
    }

    cout << "Define histograms style" << endl;

    PlotStyle(hShShUncert[iCen], markers[1], 1.1, colors[1], colors[1], "#it{z}_{T}", "Uncertainty %", false);
    // PlotStyle(hUEUncert[iCen], markers[2], 1.1, colors[2], colors[2], "#it{z}_{T}", "Uncertainty %", false);
    PlotStyle(hPurUncert[iCen], markers[3], 1.5, colors[3], colors[3], "#it{z}_{T}", "Uncertainty %", false);
    PlotStyle(hTrackIneffUncer[iCen], markers[4], 1.5, colors[4], colors[4], "#it{z}_{T}", "Uncertainty %", false);
    PlotStyle(hUEresidUncert[iCen], markers[5], 1.5, colors[5], colors[5], "#it{z}_{T}", "Uncertainty %", false);
    PlotStyle(hNMixCentUncert[iCen], markers[6], 1.5, colors[6], colors[6], "#it{z}_{T}", "Uncertainty %", false);
    PlotStyle(hZtStatUncert[iCen], 25, 1.5, kBlue + 1, 1, "#it{z}_{T}", "Uncertainty %", false);
    PlotStyle(fZt[iCen], 20, 1, kAzure - 3, kAzure - 3, "#it{z}_{T}", "1/N^{trig}dN^{charg}/d#it{z}_{T}", false);
  }

  cout << "Estimate total systematic uncertainty" << endl;
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    for (int ibin = 0; ibin < nZtBin; ibin++)
    {
      double binContShSh = hShShUncert[iCen]->GetBinContent(ibin + 1);
      // double binContUE = hUEUncert[iCen]->GetBinContent(ibin + 1);
      double binContUE = 0;
      double binContPur = hPurUncert[iCen]->GetBinContent(ibin + 1);
      double binContTrackEff = hTrackIneffUncer[iCen]->GetBinContent(ibin + 1);
      double binUEresid = hUEresidUncert[iCen]->GetBinContent(ibin + 1);
      double binContMixed = hNMixCentUncert[iCen]->GetBinContent(ibin + 1);
      double systSumQuad = TMath::Sqrt(binContShSh * binContShSh + binContUE * binContUE + binContPur * binContPur + binContTrackEff * binContTrackEff + binUEresid * binUEresid + binContMixed * binContMixed);
      hZtSystSumQuadr[iCen]->SetBinContent(ibin + 1, systSumQuad);
    }

    PlotStyle(hZtSystSumQuadr[iCen], 20, 1.1, kBlack, kBlack, "#it{z}_{T}", "Uncertainty %", false);
  }

  ////////////////////////////////////////
  //////// Plotting final results ///////////
  //////////////////////////////////////
  //TLegend *legUncert[nCen];
  TCanvas *cAllSyst = new TCanvas(Form("cAllSyst"), Form("cAllSyst"), 3 * 800, 1 * 600);
  cAllSyst->Divide(3, 1, 0.006, 0.001);

  TCanvas *cAllSyst_Stat = new TCanvas(Form("cAllSyst_Stat"), Form("cAllSyst_Stat"), 3 * 800, 1 * 600);
  cAllSyst_Stat->Divide(3, 1, 0.006, 0.001);
  //TLegend *legUncert_Stat[nCen];

  TLatex *lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetTextSize(0.04);
  lat->SetNDC();
  //TH1F *hGeneral = new TH1F("hGeneral", "hGeneral", nZtBinThin, assocZtThinner);
  TH1F *hGeneral = new TH1F("hGeneral", "hGeneral", 100, 0, 2);
  PlotStyle(hGeneral, 20, 1, kWhite, kWhite, "#it{z}_{T}", " D(#it{z}_{T}) uncert. (%) ", false);
  hGeneral->SetTitle(" ");
  hGeneral->GetYaxis()->SetTitleOffset(1.35);
  hGeneral->GetYaxis()->SetTitleSize(0.052);
  hGeneral->GetXaxis()->SetTitleSize(0.052);
  hGeneral->GetYaxis()->SetLabelSize(0.042);
  hGeneral->GetXaxis()->SetLabelSize(0.042);
  for (int iCen = 0; iCen < nCen; iCen++)
  {
//    legUncert[iCen] = LegStd(legUncert[iCen], 0.18, 0.430, 0.65, 0.66);
//    legUncert[iCen]->SetNColumns(2);
//    legUncert[iCen]->SetTextSize(0.035);
//    // legUncert[iCen]->AddEntry(hZtStatUncert[iCen], "Statistical errors", "ep");
//    legUncert[iCen]->AddEntry(hZtSystSumQuadr[iCen], "Total Syst.", "ep");
//    legUncert[iCen]->AddEntry(hShShUncert[iCen], "Bkg #it{#sigma}^{2}_{long, 5x5}", "ep");
//    legUncert[iCen]->AddEntry(hPurUncert[iCen], "Purity", "ep");
//    // legUncert[iCen]->AddEntry(hUEUncert[iCen], "UE estimation", "ep");
//    legUncert[iCen]->AddEntry(hUEresidUncert[iCen], "#it{#varepsilon}_{ME}", "ep");
//    legUncert[iCen]->AddEntry(hTrackIneffUncer[iCen], "#it{#varepsilon}_{Tracking}", "ep");
//    legUncert[iCen]->AddEntry(hNMixCentUncert[iCen], "ME centrality match", "ep");
    //////////////////////////////////////
    ////// Plot only systematics ////////
    /////////////////////////////////////
    cAllSyst->cd(iCen + 1);
    gPad->SetTickx();
    gPad->SetTicky();
    
    cAllSyst->cd(iCen + 1)->SetTopMargin(0.015);
    cAllSyst->cd(iCen + 1)->SetRightMargin(0.022);
    cAllSyst->cd(iCen + 1)->SetLeftMargin(0.1);
    cAllSyst->cd(iCen + 1)->SetBottomMargin(0.11);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadLeftMargin(0.20);
    gStyle->SetPadBottomMargin(0.15);
    hZtStatUncert[iCen]->Scale(100);

    hZtStatUncert[iCen]->GetYaxis()->SetNdivisions(511);
    hZtStatUncert[iCen]->GetYaxis()->SetDecimals();
    hGeneral->GetYaxis()->SetRangeUser(0, 75);
    hGeneral->SetTitle(" ");
    hGeneral->SetMinimum(-2);
    hGeneral->Draw("same");
    // hZtStatUncert[iCen]->GetYaxis()->SetRangeUser(0, 140);
    //  hZtStatUncert[iCen]->Draw("histsame pl");
    hZtSystSumQuadr[iCen]->GetXaxis()->SetRangeUser(rangeMin[iCen], rangeMax[iCen]);
    hPurUncert[iCen]->GetXaxis()->SetRangeUser(rangeMin[iCen], rangeMax[iCen]);
    hTrackIneffUncer[iCen]->GetXaxis()->SetRangeUser(rangeMin[iCen], rangeMax[iCen]);
    hUEresidUncert[iCen]->GetXaxis()->SetRangeUser(rangeMin[iCen], rangeMax[iCen]);
    hShShUncert[iCen]->GetXaxis()->SetRangeUser(rangeMin[iCen], rangeMax[iCen]);
    // hUEUncert[iCen]->GetXaxis()->SetRangeUser(rangeMin[iCen], rangeMax[iCen]);
    hNMixCentUncert[iCen]->GetXaxis()->SetRangeUser(rangeMin[iCen], rangeMax[iCen]);
    hZtSystSumQuadr[iCen]->Draw("histsame pl ");
    hPurUncert[iCen]->Draw("histsame pl ");
    hTrackIneffUncer[iCen]->Draw("histsame pl ");
    hUEresidUncert[iCen]->Draw("histsame pl ");
    hShShUncert[iCen]->Draw("histsame pl ");
    // hUEUncert[iCen]->Draw("histsame pl ");
    hNMixCentUncert[iCen]->Draw(" histsame pl ");
    // lat = LatexStdISO(lat, 0.180, 0.180, 0.045, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true);
    lat = new TLatex();
    lat->SetTextFont(42);
    lat->SetTextSize(0.05);
    lat->SetNDC();
    //lat->DrawLatex(0.45, 0.920, Form("#it{This Thesis}"));
    //lat->DrawLatex(0.4, 0.920 - 0.055, Form("#bf{%d#font[122]{-}%d%%} Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV ", cenBins[iCen], cenBins[iCen + 1]));
    lat->DrawLatex(0.4, 0.920 - 0.055, Form("#bf{%d#font[122]{-}%d%%}", cenBins[iCen], cenBins[iCen + 1]));
//    lat->DrawLatex(0.45, 0.920 - 2 * 0.055, Form("|#Delta#it{#varphi}_{#it{#gamma}#font[122]{-}h}| > #frac{3}{5} #it{#pi}, |#it{#eta}^{ #it{#gamma}}| < 0.67 "));
//    lat->DrawLatex(0.45, 0.920 - 3 * 0.058, Form("%2.0f < #it{p}_{T}^{ #it{#gamma}} < %2.0f GeV/#it{c} ", ptMin, ptMax));
//    lat->DrawLatex(0.45, 0.920 - 4 * 0.058, Form("#it{p}_{T}^{ h} > 0.5 GeV/#it{c} "));

    //////////////////////////////////////////////////////
    ////// Plot systematics + statistical errors ////////
    ////////////////////////////////////////////////////
    
    cAllSyst_Stat->cd(iCen + 1)->SetTopMargin(0.015);
    cAllSyst_Stat->cd(iCen + 1)->SetRightMargin(0.022);
    cAllSyst_Stat->cd(iCen + 1)->SetLeftMargin(0.1);
    cAllSyst_Stat->cd(iCen + 1)->SetBottomMargin(0.11);
//    legUncert_Stat[iCen] = LegStd(legUncert_Stat[iCen], 0.6, 0.45, 0.85, 0.56);
//    legUncert_Stat[iCen]->SetNColumns(2);
//    legUncert_Stat[iCen]->SetTextSize(0.035);
//    legUncert_Stat[iCen]->AddEntry(hZtStatUncert[iCen], "Statistical errors", "ep");
//    legUncert_Stat[iCen]->AddEntry(hZtSystSumQuadr[iCen], "Total Syst.", "ep");
    hGeneral->GetYaxis()->SetRangeUser(0, 43);
    hGeneral->GetXaxis()->SetRangeUser(0.05, 0.85);
    hGeneral->GetYaxis()->SetTitleSize(0.05);
    hGeneral->GetXaxis()->SetTitleSize(0.05);
    hGeneral->GetYaxis()->SetLabelSize(0.05);
    hGeneral->GetXaxis()->SetLabelSize(0.05);
    hGeneral->GetXaxis()->SetTitleOffset(1.0);
    hGeneral->GetYaxis()->SetTitleOffset(1.);
    
    
    hGeneral->Draw("same");
    hZtSystSumQuadr[iCen]->GetXaxis()->SetRangeUser(rangeMin[iCen], rangeMax[iCen]);
    hZtStatUncert[iCen]->GetXaxis()->SetRangeUser(rangeMin[iCen], rangeMax[iCen]);
    hZtStatUncert[iCen]->GetYaxis()->SetRangeUser(0, 50);
    hZtStatUncert[iCen]->SetTitle("");
    // hGeneral->Draw("same");
    hZtStatUncert[iCen]->Draw("histsame pl");
    hZtSystSumQuadr[iCen]->Draw("histsame pl ");
    TLatex *lat_Syst_Stat = LatexStdISO(lat_Syst_Stat, 0.55, 0.920, 0.035, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true);
    lat_Syst_Stat->SetTextSize(0.05);
  }

  //cAllSyst->cd(2);
  //legUncert[0]->Draw("same");
  cAllSyst->Print(Form("%s/SystSh%s/%s_AllSystematicsNoueVis%s.pdf", dirPlot.Data(), shshBkg.Data(), Mixed.Data(), sPtAll.Data()));

  cAllSyst_Stat->cd(2);
  //legUncert_Stat[0]->Draw("same");
  cAllSyst_Stat->Print(Form("%s/SystSh%s/%s_AllSystematicsStatNoueVis%s.pdf", dirPlot.Data(), shshBkg.Data(), Mixed.Data(), sPtAll.Data()));

  //////////////////////////////////////
  ///////// Debugging printing ////////
  /////////////////////////////////////
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    cout << "Centrality: " << cenBins[iCen] << " - " << cenBins[iCen + 1] << endl;
    cout << "Tot syst: " << endl;
    PrintSyst(hZtSystSumQuadr[iCen]);
    cout << "Tracking Efficiecy " << endl;
    PrintSyst(hTrackIneffUncer[iCen]);
    cout << "hPurity " << endl;
    PrintSyst(hPurUncert[iCen]);
    // cout << "UE statistical " << endl;
    // PrintSyst(hUEUncert[iCen]);
    cout << "UE MIXED " << endl;
    PrintSyst(hUEresidUncert[iCen]);
    cout << "ShSh " << endl;
    PrintSyst(hShShUncert[iCen]);
    cout << "nMixCent " << endl;
    PrintSyst(hNMixCentUncert[iCen]);
  }

  /////////////////////////////////////////
  ////// Plot zT with systematics ////////
  ///////////////////////////////////////

  TCanvas *cOverlap[nCen];
  TH1F *hfZtSyst[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    cOverlap[iCen] = new TCanvas(Form("cOverlapCen%d%d", cenBins[iCen], cenBins[iCen + 1]), Form("cOverlapCen%d%d", cenBins[iCen], cenBins[iCen + 1]), 800, 600);
    hfZtSyst[iCen] = (TH1F *)fZt[iCen]->Clone(Form("hsystCen%d%d", cenBins[iCen], cenBins[iCen + 1]));
    for (int ibin = 0; ibin < nZtBin; ibin++)
    {
      hfZtSyst[iCen]->SetBinError(ibin + 1, ((hZtSystSumQuadr[iCen]->GetBinContent(ibin + 1) / 100.) * fZt[iCen]->GetBinContent(ibin + 1)));
      cout << "error" << endl;
      cout << hfZtSyst[iCen]->GetBinContent(ibin + 1) << " " << fZt[iCen]->GetBinContent(ibin + 1) << endl;
    }
    fZt[iCen]->SetDirectory(0);
    hfZtSyst[iCen]->SetDirectory(0);
    gPad->SetLogy();
    hfZtSyst[iCen]->GetYaxis()->SetRangeUser(1e-3, 50);
    hfZtSyst[iCen]->SetFillColor(kBlue - 10);
    hfZtSyst[iCen]->Draw("samee2");
    fZt[iCen]->Draw("same");
    TLatex *ALICEtex1 = LatexStdISO(ALICEtex1, 0.50, 0.84, 0.04, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true);
    cOverlap[iCen]->Print(Form("%s/SystSh%s/zTwithSystCent%d_%d%s.pdf", dirPlot.Data(), shshBkg.Data(), cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));

    fSystFile->cd();
    fZt[iCen]->Write();
    hfZtSyst[iCen]->Write();
  }

  /////////////////////////////////////////////////////////////////////////
  //////////////         Compute Icp systematics       ///////////////////
  //////////////      - compute UE residual for Icp   ///////////////////
  ////////////// - get the other Icp systematics from the files ////////
  /////////////////////////////////////////////////////////////////////

  TH1F *IcpResidUESyst[nCen];
//  TFile *fzTFinal[nCen];
//  TH1F *hzTFinal[nCen];
//  for (int iCen = 0; iCen < nCen; iCen++)
//  {
//    fzTFinal[iCen] = TFile::Open(Form("Output_checkCode/fPlot0.40-1.00_Cen%d_%d_Pt18_40.root",cenBins[iCen],cenBins[iCen+1]));
//    printf("icen %d file %p\n",iCen, fzTFinal[iCen]);
//    hzTFinal[iCen] = (TH1F*) fzTFinal[iCen]
//    ->Get(Form("hZtEffCorrIso1Photon_Cen%d_%d_Pt18_40",cenBins[iCen],cenBins[iCen+1]));
//    printf("icen %d histo %p\n",iCen,hzTFinal[iCen]);
//  }
  for (int iCen = 0; iCen < nCen-1; iCen++)
    {
    IcpResidUESyst[iCen] = (TH1F *)hUEresidUncert[iCen]->Clone(Form("IcpResidUESyst_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    for (int ibin = 0; ibin < nZtBin; ibin++)
    {
      double errC = hUEresidUncert[iCen]->GetBinContent(ibin + 1);
      double errP = hUEresidUncert[nCen - 1]->GetBinContent(ibin + 1);
      double errorCP = sqrt(errC * errC + errP * errP);
      
//      double valC = hzTFinal[iCen]->GetBinContent(ibin + 1);
//      double valP = hzTFinal[nCen - 1]->GetBinContent(ibin + 1);
//      errC = errC / 100. * valC;
//      errP = errP / 100. * valP;
//      double errorCPfrac = sqrt( (errC * errC) / (valP*valP) +
//                             (valC * valC * errP * errP) / (valP * valP * valP * valP));
//
//      errorCPfrac *= 100 / (valC/valP);
//
//      printf("ICP: icen %d, ibin %d, val C %1.2f (%1.2f), val P %1.2f (%1.2f), val CP %1.2f errSum %1.2f err %1.2f\n", iCen,ibin,valC,errP,valP,errP, valC/valP, errorCP, errorCPfrac);
      
      IcpResidUESyst[iCen]->SetBinContent(ibin + 1, errorCP);
    }

    IcpResidUESyst[iCen]->SetDirectory(0);
    PlotStyle(IcpResidUESyst[iCen], 20, 1, kAzure + 3, kAzure + 3, "#font[12]{z_{T}}", "Uncertainty %", false);
  }
  
  TCanvas *cIcpResidUESyst = new TCanvas("cIcpResidUESyst", "cIcpResidUESyst", (nCen - 1) * 800, 1 * 600);
  cIcpResidUESyst->Divide((nCen - 1), 1);
  for (int iCen = 0; iCen < nCen - 1; iCen++)
  {
    cIcpResidUESyst->cd(iCen + 1);
    IcpResidUESyst[iCen]->Draw();
    //TLatex *latexIcp = LatexStdIcp(latexIcp, 0.15, 0.84, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax);
  }
  cIcpResidUESyst->Print(dirPlot + Form("/SystSh%s/SystResErrUEIcp%s.pdf", shshBkg.Data(), sPtAll.Data()));

  TH1F *hUEUncertIcp[nCen];
  TH1F *hPurUncertIcp[nCen];
  TH1F *hShShUncertIcp[nCen];
  TH1F *hNMixCentUncertIcp[nCen];
  TH1F *hZtSystSumQuadrIcp[nCen];
  TH1F *hIcpSyst[nCen];
  for (int iCen = 0; iCen < nCen - 1; iCen++)
  {
    TString sCent = Form("Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    fIcp[iCen]->Divide(fZt[nCen - 1]);
    // hUEUncertIcp[iCen] = (TH1F *)fUESyst->Get(Form("SystIcp_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    // cout << hUEUncertIcp[iCen] << endl; //
    //hPurUncertIcp[iCen] = (TH1F *)fPurSyst->Get(Form("hIcpUncertSystCen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    hPurUncertIcp[iCen] = (TH1F *)fPurSyst->Get(Form("hPurUncertFromFitIcpCen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    cout << hPurUncertIcp[iCen] << endl;
    hShShUncertIcp[iCen] = (TH1F *)fShShSyst->Get(Form("hSystUncertIcp_ShShFromFitCen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    cout << hShShUncertIcp[iCen] << endl;
    //hNMixCentUncertIcp[iCen] = (TH1F *)fNMixSyst->Get(Form("hUncertIcpNcentMix_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1])); //
    hNMixCentUncertIcp[iCen] = (TH1F *)fNMixSyst->Get(Form("hFromFitUncert_Icp_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1])); //
    cout << hNMixCentUncertIcp[iCen] << endl;
    //  hZtStatUncert[iCen] = new TH1F(Form("hZTStatistCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("hZTStatistCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), nZtBin, assocZt);
    hZtSystSumQuadrIcp[iCen] = new TH1F(Form("hZtSystSumQuadrIcp_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("hhZtSystSumQuadrIcp_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), nZtBin, assocZt);
    hIcpSyst[iCen] = new TH1F(Form("hIcpSyst_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("hhIcpSyst_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), nZtBin, assocZt);
    /////////////////////////////////////////////////////////////////////////
    //////////////         Compute quadratic sum        ////////////////////
    ///////////////////////////////////////////////////////////////////////
    for (int ibin = 0; ibin < nZtBin; ibin++)
    {
      double binContShShIcp = hShShUncertIcp[iCen]->GetBinContent(ibin + 1);
      double binContPurIcp = hPurUncertIcp[iCen]->GetBinContent(ibin + 1);
      // double binContUEIcp = hUEUncertIcp[iCen]->GetBinContent(ibin + 1);
      double binContUEIcp = 0;
      double binContNMixCenIcp = hNMixCentUncertIcp[iCen]->GetBinContent(ibin + 1);
      double binResidUE = IcpResidUESyst[iCen]->GetBinContent(ibin + 1);
      double systQuadIcp = sqrt(binContShShIcp * binContShShIcp + binContPurIcp * binContPurIcp + binContUEIcp * binContUEIcp + binContNMixCenIcp * binContNMixCenIcp + binResidUE * binResidUE);
      hZtSystSumQuadrIcp[iCen]->SetBinContent(ibin + 1, systQuadIcp);
      hIcpSyst[iCen]->SetBinContent(ibin + 1, fIcp[iCen]->GetBinContent(ibin + 1));
      hIcpSyst[iCen]->SetBinError(ibin + 1, abs((systQuadIcp / 100.) * fIcp[iCen]->GetBinContent(ibin + 1)));
      cout << "ErrSyst: " << systQuadIcp << " Icp: " << fIcp[iCen]->GetBinContent(ibin + 1) << endl;
    }

    
    PlotStyle(hShShUncertIcp[iCen], markers[1], 1.1, colors[1], colors[1], "#it{z}_{T}", "Uncertainty %", false);
    // PlotStyle(hUEUncertIcp[iCen], markers[2], 1.1, colors[2], colors[2], "#it{z}_{T}", "Uncertainty %", false);
    PlotStyle(hPurUncertIcp[iCen], markers[3], 1.5, colors[3], colors[3], "#it{z}_{T}", "Uncertainty %", false);
    PlotStyle(IcpResidUESyst[iCen], markers[5], 1.5, colors[5], colors[5], "#it{z}_{T}", "Uncertainty %", false);
    PlotStyle(hNMixCentUncertIcp[iCen], markers[6], 1.5, colors[6], colors[6], "#it{z}_{T}", "Uncertainty %", false);
    PlotStyle(hZtSystSumQuadrIcp[iCen], markers[0], 1.1, colors[0], colors[0], "#it{z}_{T}", "Uncertainty %", false);
    // PlotStyle(hZtStatUncertIcp[iCen], 20, 2, 1, "#it{z}_{T}", "Uncertainty %", false);
    hShShUncertIcp[iCen]->SetDirectory(0);
    // hUEUncertIcp[iCen]->SetDirectory(0);
    hPurUncertIcp[iCen]->SetDirectory(0);
    hNMixCentUncertIcp[iCen]->SetDirectory(0);
    hZtSystSumQuadrIcp[iCen]->SetDirectory(0);
  }

  //////////////////////////////////////
  ///////// Debugging printing ////////
  /////////////////////////////////////

  cout << "Icp systematics: " << endl;
  for (int iCen = 0; iCen < nCen - 1; iCen++)
  {
    cout << "Centrality: " << cenBins[iCen] << " - " << cenBins[iCen + 1] << endl;
    cout << "ICP Tot syst: " << endl;
    PrintSyst(hZtSystSumQuadrIcp[iCen]);
    cout << "ICP hPurity " << endl;
    PrintSyst(hPurUncertIcp[iCen]);
    cout << "ICP UE statistical " << endl;
    // PrintSyst(hUEUncertIcp[iCen]);
    cout << "ICP UE MIXED " << endl;
    PrintSyst(IcpResidUESyst[iCen]);
    cout << "ICP ShSh " << endl;
    PrintSyst(hShShUncertIcp[iCen]);
    cout << "ICP nMixCent " << endl;
    PrintSyst(hNMixCentUncertIcp[iCen]);
  }
  /////////////////////////////////////////////////
  ///////// Plotting total Icp systematic ////////
  ///////////////////////////////////////////////
  TCanvas *cSystIcp = new TCanvas(Form("cSystIcp"), Form("cSystIcp"), 3 * 800, 1 * 600);
  cSystIcp->Divide(3, 1);
  //TLegend *legSystIcp[nCen];
  for (int iCen = 0; iCen < nCen - 1; iCen++)
  {
//    legSystIcp[iCen] = LegStd(legSystIcp[iCen], 0.14, 0.48, 0.46, 0.70);
//    legSystIcp[iCen]->SetNColumns(2);
//    legSystIcp[iCen]->AddEntry(hZtSystSumQuadrIcp[iCen], "Total Syst.", "ep");
//    legSystIcp[iCen]->AddEntry(hShShUncertIcp[iCen], "Bkg #it{#sigma}^{2}_{long, 5x5}", "ep");
//    legSystIcp[iCen]->AddEntry(hPurUncertIcp[iCen], "Purity", "ep");
//    // legSystIcp[iCen]->AddEntry(hUEUncertIcp[iCen], "UE estimation", "ep");
//    legSystIcp[iCen]->AddEntry(IcpResidUESyst[iCen], "#it{#varepsilon}_{ME}", "ep");
//    legSystIcp[iCen]->AddEntry(hNMixCentUncertIcp[iCen], "ME centrality match", "ep");
    cSystIcp->cd(iCen + 1);
    cSystIcp->cd(iCen + 1)->SetTopMargin(0.015);
    cSystIcp->cd(iCen + 1)->SetRightMargin(0.022);
    cSystIcp->cd(iCen + 1)->SetLeftMargin(0.1);
    cSystIcp->cd(iCen + 1)->SetBottomMargin(0.11);
    gStyle->SetHistFillStyle(0);
    // gStyle->SetPalette(kCMYK);
    hZtSystSumQuadrIcp[iCen]->GetYaxis()->SetRangeUser(0, 80);
    hZtSystSumQuadrIcp[iCen]->GetXaxis()->SetRangeUser(0.05, 0.6);
    hPurUncertIcp[iCen]->GetXaxis()->SetRangeUser(0.05, 0.6);
    IcpResidUESyst[iCen]->GetXaxis()->SetRangeUser(0.05, 0.6);
    hShShUncertIcp[iCen]->GetXaxis()->SetRangeUser(0.05, 0.6);
    // hUEUncertIcp[iCen]->GetXaxis()->SetRangeUser(0.05, 0.6);
    hNMixCentUncertIcp[iCen]->GetXaxis()->SetRangeUser(0.05, 0.6);
    hGeneral->GetXaxis()->SetRangeUser(0.05, 0.65);
    hGeneral->GetYaxis()->SetRangeUser(0, 44);
    hGeneral->GetYaxis()->SetTitle(" #it{I}_{CP} uncert. (%) ");
    hGeneral->GetYaxis()->SetTitleSize(0.05);
    hGeneral->GetXaxis()->SetTitleSize(0.05);
    hGeneral->GetYaxis()->SetLabelSize(0.05);
    hGeneral->GetXaxis()->SetLabelSize(0.05);
    hGeneral->GetXaxis()->SetTitleOffset(1.1);
    hGeneral->GetYaxis()->SetTitleOffset(1.);
    hGeneral->Draw("same");
    hZtSystSumQuadrIcp[iCen]->Draw("same hist pl ");
    hShShUncertIcp[iCen]->Draw("same hist pl ");
    // hUEUncertIcp[iCen]->Draw("same hist pl ");
    hPurUncertIcp[iCen]->Draw("same hist pl ");
    IcpResidUESyst[iCen]->Draw("same hist pl ");
    hNMixCentUncertIcp[iCen]->Draw("same hist pl ");
    TLatex *latIcp = new TLatex();
    latIcp->SetTextFont(42);
    latIcp->SetTextSize(0.05);
    latIcp->SetNDC();
    latIcp->DrawLatex(0.5, 0.84, Form("#bf{%d#font[122]{-}%d%% / 50#font[122]{-}90%%}", cenBins[iCen], cenBins[iCen + 1]));
    //TLatex *ALICEtexIcp1 = LatexStdIcp(ALICEtexIcp1, 0.46, 0.92, cenBins[0], cenBins[1], ptMin, ptMax);
    //legSystIcp[iCen]->Draw("same");
  }
  
  cSystIcp->cd(3);
  
  TLatex *ALICEtexIcp1 = LatexStdIcp(ALICEtexIcp1, 0., 0.85, cenBins[0], cenBins[1], ptMin, ptMax);
  ALICEtexIcp1->SetTextSize(0.05);
  
  TLegend * legSystIcp = LegStd(legSystIcp, 0., 0.2, 0.7, 0.50);
  legSystIcp->SetTextSize(0.05);
  legSystIcp->SetNColumns(2);
  legSystIcp->AddEntry(hZtSystSumQuadrIcp[0], "Total Syst.", "pl");
  legSystIcp->AddEntry(hShShUncertIcp[0], "Bkg #it{#sigma}^{2}_{long, 5x5}", "pl");
  legSystIcp->AddEntry(hPurUncertIcp[0], "Purity", "pl");
  // legSystIcp->AddEntry(hUEUncertIcp[iCen], "UE estimation", "pl");
  legSystIcp->AddEntry(hTrackIneffUncer[0], "#it{#varepsilon}_{Tracking}", "pl");

  legSystIcp->AddEntry(IcpResidUESyst[0], "#it{#varepsilon}_{ME}", "pl");
  legSystIcp->AddEntry(hNMixCentUncertIcp[0], "ME centrality match", "pl");
  legSystIcp->Draw("same");
  cSystIcp->Print(Form("%s/SystSh%s/%s_cSystIcp_AllSystematicsCen%s.pdf", dirPlot.Data(), shshBkg.Data(), Mixed.Data(), sPtAll.Data()));

  int kMarkerColIcp[] = {kAzure + 2, kOrange + 8, kCyan - 2};
  int kMarkerStyleIcp[] = {21, 20, 24};
  // TCanvas *cIcp = new TCanvas("Icp", "Icp", 800, 600);
  TCanvas *cIcp = canvasStd("Icp", 1, 1);
  TLegend *legIcp = LegStd(legIcp, 0.14, 0.72, 0.360, 0.96);
  
  cIcp->cd(1)->SetTopMargin(0.015);
  cIcp->cd(1)->SetRightMargin(0.015);
  cIcp->cd(1)->SetLeftMargin(0.12);
  cIcp->cd(1)->SetBottomMargin(0.13);
  
  // hIcpSyst[0]->SetFillStyle(0);
  hIcpSyst[0]->SetFillColorAlpha(kAzure + 5, 0.30);
  hIcpSyst[0]->SetLineColorAlpha(kOrange + 7, 1.00);
  hIcpSyst[0]->GetXaxis()->SetRangeUser(0.0, 0.6);
  hIcpSyst[1]->GetXaxis()->SetRangeUser(0.0, 0.6);
  hIcpSyst[1]->SetLineColorAlpha(kMarkerColIcp[1], 1.00);
  hIcpSyst[1]->SetFillColorAlpha(kMarkerColIcp[1], 0.30);
  TH1F *hGeneralIcp = new TH1F("hGeneralIcp", "hGeneralIcp", 10, 0, 0.7);
  hGeneralIcp->SetDirectory(0);
  hGeneralIcp->Draw("histsame");
  PlotStyle(hGeneralIcp, 20, 1, kWhite, kWhite, "#it{z}_{T}", "#it{I}_{CP}", false);
  for (int iCen = 0; iCen < nCen - 1; iCen++)
  {
    hIcpSyst[iCen]->SetDirectory(0);
    PlotStyle(hIcpSyst[iCen], kMarkerStyleIcp[iCen], 1, kMarkerColIcp[iCen], kMarkerColIcp[iCen], "#it{z}_{T}", "#it{I}_{CP}", false);
    PlotStyle(fIcp[iCen], kMarkerStyleIcp[iCen], 1, kMarkerColIcp[iCen], kMarkerColIcp[iCen], "#it{z}_{T}", "#it{I}_{CP}", false);
    //hGeneralIcp->GetYaxis()->SetRangeUser(-0.1, 2.25);
    hGeneralIcp->GetYaxis()->SetRangeUser(-0.075, 1.95);
    hGeneralIcp->GetXaxis()->SetRangeUser(0.05, 0.65);
    hGeneralIcp->SetTitle(" ");
    hGeneralIcp->GetYaxis()->CenterTitle(false);
    hGeneralIcp->GetYaxis()->SetTitleSize(0.05);
    hGeneralIcp->GetYaxis()->SetTitleOffset(1.1);
    hGeneralIcp->GetYaxis()->SetLabelSize(0.04);
    hGeneralIcp->GetXaxis()->SetRangeUser(0.05, 1.05);
    hGeneralIcp->GetXaxis()->SetLabelSize(0.04);
    hGeneralIcp->GetXaxis()->SetTitleSize(0.05);
    hGeneralIcp->SetLineWidth(0);
    hIcpSyst[iCen]->SetLineWidth(0);
    hIcpSyst[iCen]->Draw("samee2");
    fIcp[iCen]->Draw("EP X0same");
    // grIcpNLOmedian[iCen]->Draw("pl3 same");
  }
  
  legIcp->SetTextSize(0.04);
  legIcp->AddEntry(fIcp[0], Form("0#font[122]{-}30%%  / 50#font[122]{-}90%% stat. unc."), "ep");
  legIcp->AddEntry(hIcpSyst[0], Form("syst. unc."), "f");
  legIcp->AddEntry(fIcp[1], Form("30#font[122]{-}50%% / 50#font[122]{-}90%% stat. unc."), "ep");
  legIcp->AddEntry(hIcpSyst[1], Form("syst. unc."), "f");
  legIcp->Draw("same");

  TGraph *line = DrawLine(line, 0, 1, 1.2, 1);
  line->Draw("same");
  TGraph *line1 = DrawLine(line1, 0, 0.5, 1.2, 0.5);
  line1->Draw("same");
  // TLatex *ALICEtexIcp2 = LatexStdIcp(ALICEtexIcp2, 0.52, 0.92, cenBins[0], cenBins[1], ptMin, ptMax);
  TLatex *ALICEtexIcp2 = new TLatex();
  ALICEtexIcp2->SetTextFont(42);
  ALICEtexIcp2->SetTextSize(0.04);
  ALICEtexIcp2->SetNDC();
  //ALICEtexIcp2->DrawLatex(0.56, 0.86, Form("#it{This Thesis}"));
  ALICEtexIcp2->DrawLatex(0.6, 0.92, Form("ALICE"));
  ALICEtexIcp2->DrawLatex(0.6, 0.92 - 0.055, Form("Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV "));
  ALICEtexIcp2->DrawLatex(0.6, 0.92 - 2 * 0.06, Form("|#Delta#it{#varphi}_{#it{#gamma}#font[122]{-}h}| > #frac{3}{5} #it{#pi}, |#it{#eta}^{ #it{#gamma}}| < 0.67 "));
  ALICEtexIcp2->DrawLatex(0.6, 0.92 - 3 * 0.06, Form("%2.0f < #it{p}_{T}^{ #it{#gamma}} < %2.0f GeV/#it{c} ", ptMin, ptMax));
  ALICEtexIcp2->DrawLatex(0.6, 0.92 - 4 * 0.06, Form("#it{p}_{T}^{ h} > 0.5 GeV/#it{c} "));
  cIcp->Print(Form("%s/SystSh%s/%s_cSystIcp%s.pdf", dirPlot.Data(), shshBkg.Data(), Mixed.Data(), sPtAll.Data()));

  fSystFile->Close();
  
  TFile * fileICP = new TFile(Form("%s/fIcpWithSyst%s%s%s%s.root", dirPlot.Data(), Mixed.Data(), shshBkg.Data(), sMirror.Data(), sPtAll.Data()), "RECREATE");
  printf("Write Icp in: %s\n",fileICP->GetName());
  
  
  for(Int_t iCen = 0; iCen < nCen-1; iCen++)
  {
    hIcpSyst[iCen]->Write();
    fIcp[iCen]->Write();
  }
  
  fileICP->Close();
}

void PrintSyst(TH1F *hSyst)
{
  for (int ibin = 0; ibin < nZtBin; ibin++)
  {
    cout << "Error: " << hSyst->GetBinContent(ibin + 1) << "%" << endl;
  }
}

