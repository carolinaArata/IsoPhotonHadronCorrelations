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
#include <vector>
#include <algorithm>
#include <iterator>

using std::cout;
using std::endl;


int nZtBin = 6;
double assocZt[] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 1.00};

Int_t nPtTrig = 8;
Float_t PtTrigger[] = {14, 16, 18, 20, 25, 30, 40, 50, 60};

Int_t const nIso = 2;
Int_t const nShSh = 2;
Int_t nCenBin = 3;
Int_t kColor[] = {2, 4, kOrange + 7, kAzure + 3, kPink - 4, kViolet - 7, kBlue + 2, kTeal - 6};
Bool_t deltaPhiPlot = kFALSE;

void ZtFunction(TH1F *hDeltaPhi, TH1F *hZT, int bin);
void PlotStyle(TH1F *hPlot, int kMarker, double kMarkerSize, int kColor, TString titleX, TString titleY);
TLatex *LatexStd(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax);

void TrackIneffSyst(Float_t ptMin = 18, Float_t ptMax = 40, bool Mirror = true, TString fName = "fPlot", TString shshBkg = "0.40-1.00", TString dirRefData = "~/work/histogram/FromScratch/checkCode", TString dirFiles = "~/work/histogram/FromScratch/checkCodeTrackEff", bool b0_30 = true)
{

  TString sPtAll = Form("_Pt%2.0f_%2.0f", ptMin, ptMax); // define pT all range

  // define centrality bins
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

  // define directory for the plots and the root file where I save the systematics
  TString dirPlot;
  TFile *fSystTrackIneff;
  if (b0_30)
  {
    dirPlot = "~/work/histogram/Systematics_checkCode0_30/TrackIneff";
    gSystem->Exec(Form("mkdir %s", dirPlot.Data()));
    fSystTrackIneff = new TFile(Form("%s/fSystTrackIneff%s.root", dirPlot.Data(), sPtAll.Data()), "RECREATE");
  }
  else if (!b0_30)
  {
    dirPlot = "~/work/histogram/Systematics_checkCode/TrackIneff";
    gSystem->Exec(Form("mkdir %s", dirPlot.Data()));
    fSystTrackIneff = new TFile(Form("%s/fSystTrackIneff%s.root", dirPlot.Data(), sPtAll.Data()), "RECREATE");
  }

  // ZT Distribution plots
  TH1F *hZtPlot[nCen];
  TH1F *hZtPlot_MC_Gen[nCen];
  TH1F *hZtPlot_MC_Rec[nCen];
  TH1F *hZtPlot_MC_Ratio[nCen];
  TH1F *hZtPlotTrackIneff_MC_Gen[nCen];
  TH1F *hZtPlotTrackIneff_MC_Rec[nCen];
  TH1F *hZtPlotTrackIneff_MC_Ratio[nCen];

  TFile *fPlot[nCen];
  TFile *fPlotTrackIneff[nCen];

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    fPlot[iCen] = TFile::Open(Form("%s/%s%s_%s%s.root", dirRefData.Data(), fName.Data(), shshBkg.Data(), sCent.Data(), sPtAll.Data()));         // zT function from the std MC
    fPlotTrackIneff[iCen] = TFile::Open(Form("%s/%s%s_%s%s.root", dirFiles.Data(), fName.Data(), shshBkg.Data(), sCent.Data(), sPtAll.Data())); // zT function from the MC with TrackIneff

    /////////////////////////////////////////////////////////////////////////////////////////
    ///////// Getter MC Zt distribution from Std MC and from MC with TrackIneff ////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    hZtPlot_MC_Gen[iCen] = (TH1F *)fPlot[iCen]->Get(Form("hZtMCGenIso1Photon_%s_Pt%2.0f_%2.0f", sCent.Data(), ptMin, ptMax));
    hZtPlot_MC_Rec[iCen] = (TH1F *)fPlot[iCen]->Get(Form("hZtMCRecIso1Photon_%s_Pt%2.0f_%2.0f", sCent.Data(), ptMin, ptMax));
    hZtPlot_MC_Ratio[iCen] = (TH1F *)fPlot[iCen]->Get(Form("hRatioEffCorrIso1Photon"));

    hZtPlotTrackIneff_MC_Gen[iCen] = (TH1F *)fPlotTrackIneff[iCen]->Get(Form("hZtMCGenIso1Photon_%s_Pt%2.0f_%2.0f", sCent.Data(), ptMin, ptMax));
    hZtPlotTrackIneff_MC_Rec[iCen] = (TH1F *)fPlotTrackIneff[iCen]->Get(Form("hZtMCRecIso1Photon_%s_Pt%2.0f_%2.0f", sCent.Data(), ptMin, ptMax));
    hZtPlotTrackIneff_MC_Ratio[iCen] = (TH1F *)fPlotTrackIneff[iCen]->Get(Form("hRatioEffCorrIso1Photon"));
    /////////// Plot style //////////////
    PlotStyle(hZtPlot_MC_Gen[iCen], 21, 1, kAzure + 7, "#font[12]{z}_{T}", "1/N^{trig}dN^{charg}/d#font[12]{z}_{T}");
    PlotStyle(hZtPlot_MC_Rec[iCen], 21, 1, kOrange - 3, "#font[12]{z}_{T}", "1/N^{trig}dN^{charg}/d#font[12]{z}_{T}");
    PlotStyle(hZtPlot_MC_Ratio[iCen], 21, 1, kCyan + 2, "#font[12]{z}_{T}", "1/N^{trig}dN^{charg}/d#font[12]{z}_{T}");
    PlotStyle(hZtPlotTrackIneff_MC_Gen[iCen], 89, 1, kAzure + 3, "#font[12]{z}_{T}", "1/N^{trig}dN^{charg}/d#font[12]{z}_{T}");
    PlotStyle(hZtPlotTrackIneff_MC_Rec[iCen], 89, 1, kOrange + 10, "#font[12]{z}_{T}", "1/N^{trig}dN^{charg}/d#font[12]{z}_{T}");
    PlotStyle(hZtPlotTrackIneff_MC_Ratio[iCen], 24, 1, kCyan + 2, "#font[12]{z}_{T}", "1/N^{trig}dN^{charg}/d#font[12]{z}_{T}");
  }
  /////////////////////////////////////////////////////////////////////////////////////////
  //// plot Zt distributions from Std MC and from MC with TrackIneff for gen and rec /////
  ///////////////////////////////////////////////////////////////////////////////////////
  TCanvas *cZtOverlap_MC_Gen = new TCanvas(Form("cZtOverlap_MC_GenCen"), Form("cZtOverlap_MC_GenCen"), 2 * 800, 2 * 600);
  cZtOverlap_MC_Gen->Divide(2, 2);
  TLegend *legZtOverlap_MC_Gen[nCen];
  TLatex *ALICEtex[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    cZtOverlap_MC_Gen->cd(iCen + 1);
    gPad->SetLogy();
    hZtPlot_MC_Gen[iCen]->Draw();
    hZtPlotTrackIneff_MC_Gen[iCen]->Draw("same");
    legZtOverlap_MC_Gen[iCen] = new TLegend(0.7, 0.50, 0.95, 0.60);
    legZtOverlap_MC_Gen[iCen]->SetFillColor(kWhite);
    legZtOverlap_MC_Gen[iCen]->SetLineWidth(0);
    legZtOverlap_MC_Gen[iCen]->AddEntry(hZtPlot_MC_Gen[iCen], "Nominal", "lp");
    legZtOverlap_MC_Gen[iCen]->AddEntry(hZtPlotTrackIneff_MC_Gen[iCen], "zT x TrackIneff", "lp");
    legZtOverlap_MC_Gen[iCen]->Draw("same");
    ALICEtex[iCen] = LatexStd(ALICEtex[iCen], 0.550, 0.80, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax);
    cZtOverlap_MC_Gen->Update();
  }
  cZtOverlap_MC_Gen->Print(Form("%s/zTOverlap_MCGen_TrackIneff%s.pdf", dirPlot.Data(), sPtAll.Data()));
  TCanvas *cZtOverlap_MC_Rec = new TCanvas(Form("cZtOverlap_MC_Rec"), Form("cZtOverlap_MC_Rec"), 2 * 800, 2 * 600);
  TLegend *legZtOverlap_MC_Rec[nCen];
  TLatex *ALICEtex1[nCen];
  cZtOverlap_MC_Rec->Divide(2, 2);
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    cZtOverlap_MC_Rec->cd(iCen + 1);
    gPad->SetLogy();
    hZtPlot_MC_Rec[iCen]->Draw("SAME");
    hZtPlotTrackIneff_MC_Rec[iCen]->Draw("same");
    legZtOverlap_MC_Rec[iCen] = new TLegend(0.7, 0.50, 0.95, 0.60);
    legZtOverlap_MC_Rec[iCen]->SetFillColor(kWhite);
    legZtOverlap_MC_Rec[iCen]->SetLineWidth(0);
    legZtOverlap_MC_Rec[iCen]->AddEntry(hZtPlot_MC_Rec[iCen], "Nominal", "lp");
    legZtOverlap_MC_Rec[iCen]->AddEntry(hZtPlotTrackIneff_MC_Rec[iCen], "zT x TrackIneff", "lp");
    legZtOverlap_MC_Rec[iCen]->Draw("same");
    ALICEtex1[iCen] = LatexStd(ALICEtex1[iCen], 0.550, 0.80, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax);
    cZtOverlap_MC_Rec->Update();
  }
  cZtOverlap_MC_Rec->Print(Form("%s/zTOverlap_MCRec_TrackIneff%s.pdf", dirPlot.Data(), sPtAll.Data()));

  TH1F *hRatio_TrackIneff_Nom_Rec[nCen];
  TCanvas *cRatio_TrackIneff_Nom_Rec = new TCanvas("cRatio_TrackIneff_Nom_Rec", "cRatio_TrackIneff_Nom_Rec", 2 * 800, 2 * 600);
  cRatio_TrackIneff_Nom_Rec->Divide(2, 2);
  TLatex *ALICEtex2[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    hRatio_TrackIneff_Nom_Rec[iCen] = (TH1F *)hZtPlotTrackIneff_MC_Rec[iCen]->Clone(Form("hZtPlotTrackIneff_MC_RecCen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    cRatio_TrackIneff_Nom_Rec->cd(iCen + 1);
    hRatio_TrackIneff_Nom_Rec[iCen]->Divide(hZtPlot_MC_Rec[iCen]);
    hRatio_TrackIneff_Nom_Rec[iCen]->SetDirectory(0);
    hRatio_TrackIneff_Nom_Rec[iCen]->SetTitle("");
    hRatio_TrackIneff_Nom_Rec[iCen]->GetYaxis()->SetTitle(" TrackIneff / Nominal");
    hRatio_TrackIneff_Nom_Rec[iCen]->Draw();
    ALICEtex2[iCen] = LatexStd(ALICEtex2[iCen], 0.200, 0.82, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax);
    cRatio_TrackIneff_Nom_Rec->Update();
  }
  cRatio_TrackIneff_Nom_Rec->cd();
  TPad *padtitle = new TPad("padtitle", "padtitle", 0.3, 0.95, 0.7, 0.99);
  padtitle->Draw();
  padtitle->cd();
  padtitle->SetFillStyle(0);

  auto tex = new TLatex(0.5, 0.5, "#bf{#bf{Reconstructed}}");
  tex->SetTextAlign(22);
  tex->SetTextFont(42);
  tex->SetTextSize(0.6);
  tex->Draw();
  cRatio_TrackIneff_Nom_Rec->Update();
  cRatio_TrackIneff_Nom_Rec->Print(Form("%s/Ratio_Rec%s.pdf", dirPlot.Data(), sPtAll.Data()));

  TH1F *hRatio_TrackIneff_Nom_Gen[nCen];
  TCanvas *cRatio_TrackIneff_Nom_Gen = new TCanvas("cRatio_TrackIneff_Nom_Gen", "cRatio_TrackIneff_Nom_Gen", 2 * 800, 2 * 600);
  cRatio_TrackIneff_Nom_Gen->Divide(2, 2);
  TLatex *ALICEtex3[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    hRatio_TrackIneff_Nom_Gen[iCen] = (TH1F *)hZtPlotTrackIneff_MC_Gen[iCen]->Clone(Form("hZtPlotTrackIneff_MC_GenCen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    hRatio_TrackIneff_Nom_Gen[iCen]->Divide(hZtPlot_MC_Gen[iCen]);
    hRatio_TrackIneff_Nom_Gen[iCen]->SetDirectory(0);
    hRatio_TrackIneff_Nom_Gen[iCen]->SetTitle(Form("hZt Gen Cen %d-%d %%", cenBins[iCen], cenBins[iCen + 1]));
    hRatio_TrackIneff_Nom_Gen[iCen]->GetYaxis()->SetTitle(" TrackIneff / Nominal");
    hRatio_TrackIneff_Nom_Gen[iCen]->SetTitle("");
    cRatio_TrackIneff_Nom_Gen->cd(iCen + 1);
    hRatio_TrackIneff_Nom_Gen[iCen]->Draw();
    ALICEtex3[iCen] = LatexStd(ALICEtex3[iCen], 0.200, 0.82, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax);
    cRatio_TrackIneff_Nom_Gen->Update();
  }
  cRatio_TrackIneff_Nom_Gen->cd();
  TPad *padtitle1 = new TPad("padtitle1", "padtitle1", 0.3, 0.95, 0.7, 0.99);
  padtitle1->Draw();
  padtitle1->cd();
  padtitle1->SetFillStyle(0);

  auto tex1 = new TLatex(0.5, 0.5, "Generated");
  tex1->SetTextAlign(22);
  tex1->SetTextFont(42);
  tex1->SetTextSize(0.6);
  tex1->Draw();
  cRatio_TrackIneff_Nom_Gen->Update();
  cRatio_TrackIneff_Nom_Gen->Print(Form("%s/Ratio_Gen%s.pdf", dirPlot.Data(), sPtAll.Data()));

  TCanvas *cZtOverlap_MC_Ratio[nCen];
  TLegend *legZtOverlap_MC_Ratio[nCen];
  TLatex *ALICEtex4[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    cZtOverlap_MC_Ratio[iCen] = new TCanvas(Form("cZtOverlap_MC_RatioCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("cZtOverlap_MC_RatioCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), 800, 600);
    // gPad->SetLogy();
    hZtPlotTrackIneff_MC_Ratio[iCen]->Draw("same");
    hZtPlot_MC_Ratio[iCen]->Draw("same");
    legZtOverlap_MC_Ratio[iCen] = new TLegend(0.7, 0.50, 0.95, 0.60);
    legZtOverlap_MC_Ratio[iCen]->SetFillColor(kWhite);
    legZtOverlap_MC_Ratio[iCen]->SetLineWidth(0);
    legZtOverlap_MC_Ratio[iCen]->AddEntry(hZtPlot_MC_Ratio[iCen], "Nominal", "lp");
    legZtOverlap_MC_Ratio[iCen]->AddEntry(hZtPlotTrackIneff_MC_Ratio[iCen], "zT x TrackIneff", "lp");
    legZtOverlap_MC_Ratio[iCen]->Draw("same");
    ALICEtex4[iCen] = LatexStd(ALICEtex4[iCen], 0.550, 0.80, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax);
    cZtOverlap_MC_Ratio[iCen]->Update();
    cZtOverlap_MC_Ratio[iCen]->Print(Form("%s/zTOverlapMC_RatioTrackIneffCen%d_%d%s.pdf", dirPlot.Data(), cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
  }

  TCanvas *cSystTrackIneff[nCen];
  TF1 *hfitpol[nCen];
  TH1F *hRatio_TrackInef_Nom[nCen];
  TLatex *ALICEtex5[nCen];
  TCanvas *cRatio_TrackInef_Nom[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    cRatio_TrackInef_Nom[iCen] = new TCanvas("cRatio_TrackInef_Nom", "cRatio_TrackInef_Nom", 2 * 800, 2 * 600);
    hRatio_TrackInef_Nom[iCen] = (TH1F *)hZtPlotTrackIneff_MC_Ratio[iCen]->Clone(Form("hRatio_TrackInef_NomCen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
    hRatio_TrackInef_Nom[iCen]->Add(hZtPlot_MC_Ratio[iCen], -1);
    hRatio_TrackInef_Nom[iCen]->Divide(hZtPlot_MC_Ratio[iCen]);
    hRatio_TrackInef_Nom[iCen]->SetDirectory(0);
    hRatio_TrackInef_Nom[iCen]->SetMarkerStyle(21);
    hRatio_TrackInef_Nom[iCen]->SetMarkerSize(2);
    hRatio_TrackInef_Nom[iCen]->Scale(100);
    // cRatio_TrackInef_Nom->cd(iCen + 1);
    hRatio_TrackInef_Nom[iCen]->SetTitle("");
    hRatio_TrackInef_Nom[iCen]->GetYaxis()->SetRangeUser(-0.1, 8);
    hRatio_TrackInef_Nom[iCen]->GetYaxis()->SetTitle("Uncertainty %");
    hRatio_TrackInef_Nom[iCen]->Draw();
    hfitpol[iCen] = new TF1(Form("fitpol0Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), "expo", 0, 1.00);
    // hRatio_TrackInef_Nom[iCen]->Fit(Form("fitpol0Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), "");
    hRatio_TrackInef_Nom[iCen]->Draw("histPL");
    // hfitpol[iCen]->Draw("");
    ALICEtex5[iCen] = LatexStd(ALICEtex5[iCen], 0.5, 0.84, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax);
    cRatio_TrackInef_Nom[iCen]->Print(Form("%s/Ratio%sCen%d_%d.pdf", dirPlot.Data(), sPtAll.Data(), cenBins[iCen], cenBins[iCen + 1]));
  }

  TH1F *hSystTrackIneff[nCen];
  TH1F *hTrackIneffUncerFromFit[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    hTrackIneffUncerFromFit[iCen] = new TH1F(Form("hTrackIneffUncerFromFitCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("hTrackIneffUncerFromFitCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), nZtBin, assocZt);
    for (int ibin = 0; ibin < nZtBin; ibin++)
    {
      hTrackIneffUncerFromFit[iCen]->SetBinContent(ibin + 1, (hRatio_TrackInef_Nom[iCen]->GetBinContent(ibin + 1)));
    }
  }
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    fSystTrackIneff->cd();
    hfitpol[iCen]->Write();
    hTrackIneffUncerFromFit[iCen]->Write();
  }

  // fPlot->Close();
  // fPlotTrackIneff->Close();
}

void PlotStyle(TH1F *hPlot, int kMarker, double kMarkerSize, int kColor, TString titleX, TString titleY)
{
  gStyle->SetTitleX(0.56);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111111);
  // gStyle->SetLabelSize(.5, "XY");
  gStyle->SetLineScalePS(1);
  gStyle->SetTitleAlign(23);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadLeftMargin(0.12);
  hPlot->SetMarkerStyle(kMarker);
  hPlot->SetMarkerSize(kMarkerSize);
  // hPlot->SetMarkerSize(2);
  hPlot->SetMarkerColor(kColor);
  hPlot->SetLineColor(kColor);
  // hPlot->SetLineWidth(2);
  hPlot->GetYaxis()->SetTitle(Form("%s", titleY.Data()));
  hPlot->GetXaxis()->SetTitle(Form("%s", titleX.Data()));
  // hPlot->GetXaxis()->SetLabelSize(.05);
  // hPlot->GetYaxis()->SetLabelSize(.05);
  //  leg->SetFillColor(kWhite);
  //  leg->SetLineColor(0);
}

TLatex *LatexStd(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax)
{
  lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetTextSize(0.04);
  lat->SetNDC();
  lat->DrawLatex(xpos, ypos, Form("#it{This Thesis}"));
  lat->DrawLatex(xpos, ypos - 0.06, Form("#bf{%d-%d %%} Pb#font[122]{-}Pb, #sqrt{s_{NN}} = 5.02 TeV", cenMin, cenMax));
  lat->DrawLatex(xpos, ypos - 2 * 0.06, Form("|#it{#eta}^{ #it{#gamma}}| < 0.67 , %2.0f < #it{p}_{T}^{#gamma} < %2.0f GeV/#it{c}", ptMin, ptMax));

  return lat;
}
