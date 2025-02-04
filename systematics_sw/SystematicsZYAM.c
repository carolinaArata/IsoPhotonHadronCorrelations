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

// Int_t nCen = 3;
// Int_t cenBins[] = {0, 30, 50, 90};
//  Int_t cenBins[] = {30, 50, 90};

// int nZtBin = 10;
// double assocZt[] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20};
int nZtBin = 6;
double assocZt[] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 1.00};

int npt = 11;
float ptTrig[] = {14, 16, 18, 20, 25, 30, 35, 40, 50, 60, 80};

Int_t const nIso = 2;
Int_t const nShSh = 1;

// Int_t kColorData[2][2] = {{}{}};
Int_t kMarkStyle[2] = {20, 24};
Int_t kMarkStyleIso_NotIso[2] = {21, 48};

//{2, 4, kOrange + 7, kAzure + 3, kPink - 4, kViolet - 7, kBlue + 2, kTeal - 6};
Int_t kColorMC[] = {2, 4, kOrange + 7, kAzure + 3, kPink - 4, kViolet - 7, kBlue + 2, kTeal - 6};
Int_t kColorIsoREC[] = {kAzure + 2, kOrange - 3};
Int_t kColorShSh[] = {kAzure + 2, kViolet};
TString sNamePtTrigGen[] = {"Photon", "Pi0"};
TLatex *LatexStd(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax, Bool_t bCen);
void PlotStyle(TH1F *hPlot, int kMarker, double kMarkerSize, int kColor, TString titleX, TString titleY);
TLatex *LatexDPhi(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax);
TLegend *LegStd(TLegend *leg, double xpos1, double ypos1, double xpos2, double ypos2);
TH1F *SumPtBinXzt(TH1F *hTrigSame, Float_t PtTrigger[npt], int index1, int index2, TH1F *hzTbin[npt], TH1F *hzTbinAll, TH1F *hPur, TF1 *fPur, double systPur, Bool_t bData);
void SystematicsZYAM(Float_t ptMin = 18, Float_t ptMax = 40, bool Mirror = true, TString sMixed = "Mixed", TString shshBkg = "0.40-1.00", TString dirFileResults = "~/work/histogram/FromScratch/checkCode", bool b0_30 = true)
{

  Int_t nCen;
  std::vector<Int_t> cenBins;
  TString dirPlot;
  if (b0_30)
  {
    nCen = 3;
    // Int_t cenBins[] = {0, 10, 30, 50, 90};
    cenBins.push_back(0);
    cenBins.push_back(30);
    cenBins.push_back(50);
    cenBins.push_back(90);
    dirPlot = "~/work/histogram/Systematics_checkCode0_30/crossZYAM";
  }
  else if (!b0_30)
  {
    nCen = 4;
    cenBins.push_back(0);
    cenBins.push_back(10);
    cenBins.push_back(30);
    cenBins.push_back(50);
    cenBins.push_back(90);
    dirPlot = "~/work/histogram/Systematics_checkCode/crossZYAM";
  }

  TString mirror;
  if (Mirror)
    mirror = "Mirror";
  else
    mirror = "NoMirror";

  TString shshString[2] = {"0.10-0.30", shshBkg};
  TString sPtAll = Form("_Pt%2.0f_%2.0f", ptMin, ptMax);

  gSystem->Exec(Form("mkdir %s", dirPlot.Data()));

  // index pt start and stop
  int nsize = sizeof(ptTrig) / sizeof(ptTrig[0]);
  auto itr1 = find(ptTrig, ptTrig + nsize, ptMin);
  auto itr2 = find(ptTrig, ptTrig + nsize, ptMax);
  int index1 = distance(ptTrig, itr1);
  int index2 = distance(ptTrig, itr2);
  int nPtTrig = index2 - index1;
  cout << index1 << "___" << index2 << ", " << nPtTrig << endl;
  
  TH1F *histPur[nCen];
  TH1F *histPurStat[nCen];
  TF1 *funcPur[nCen];

  TFile *fPlotMix[nCen];
  TFile *fPlotZYAM[nCen];
  TH1F *hZtMix[nCen];
  TH1F *hZtZYAM[nCen];
  TH1F *hZtUncertSyst[nCen];
  TH1F *hdPhiSamMCRecNoUE[nCen][nIso][nShSh][nZtBin][nPtTrig];
  TH1F *hdPhiSamMCRec[nCen][nIso][nShSh][nZtBin][nPtTrig];
  TH1F *hdPhiMixMCRec[nCen][nIso][nShSh][nZtBin][nPtTrig];
  TH1F *hdPhiSamMCRecNoUEZYAM[nCen][nIso][nShSh][nZtBin][nPtTrig];
  TH1F *hdPhiSamMCRecNoUEDiff[nCen][nZtBin][nPtTrig];
  TH1F *hZtMixPtBin[nCen][nPtTrig];
  TH1F *hZtZYAMPtBin[nCen][nPtTrig];

  TF1 *fitSamNOUE[nCen][nIso][nShSh][nZtBin][nPtTrig];
  TH1F *histFitNoUE[nCen][nPtTrig];
  TH1F *histSystErrFitNoUE[nCen][nPtTrig];
  TH1F *histSystErrFitNoUEXPtBin[nCen][nPtTrig];
  TH1F *histSystErrFitNoUEXPtAll[nCen];
  TH1F *hTriggerSamMCRec[nCen];
  double fitSystErrNoUE[nCen][nZtBin][nPtTrig];

  TFile *fZYAMSyst = new TFile(Form("%s/fZYAMSyst%s%s%s.root", dirPlot.Data(), sMixed.Data(), shshBkg.Data(), sPtAll.Data()), "RECREATE");
  TFile *fileMC = TFile::Open(Form("~/work/histogram/MCPtAssoc500/MC_GJ_0_90.root"));
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    // Purity
    TFile *fPurity = new TFile("~/work/histogram/IsoPhotonHadronCorrelations/Purity_IsoSig1.5_M02Sig0.10-0.30_IsoBkg_4.0_25.0_M02Bkg0.40_2.00_LHC15o_18qr_L1MB.root");
    cout << iCen << " " << cenBins[iCen] << "-" << cenBins[iCen + 1] << endl;
    histPur[iCen] = (TH1F *)fPurity->Get(Form("Purity_Cen%d_R0.2_Sys", iCen));
    histPurStat[iCen] = (TH1F *)fPurity->Get(Form("Purity_Cen%d_R0.2", iCen));
    funcPur[iCen] = histPur[iCen]->GetFunction("purityFitCombinedSigmoid");
    // Results
    cout << "Open files with results from Mixed and from ZYAM Monte Carlo" << endl;
    fPlotMix[iCen] = new TFile(Form("%s/fPlot%s%s%s.root", dirFileResults.Data(), shshBkg.Data(), sCent.Data(), sPtAll.Data()));
    fPlotZYAM[iCen] = new TFile(Form("%sZYAM/fPlot%s%s%s.root", dirFileResults.Data(), shshBkg.Data(), sCent.Data(), sPtAll.Data()));
    hTriggerSamMCRec[iCen] = (TH1F *)fPlotMix[iCen]->Get("AnaPhotonHadronCorr_Iso1_ShSh0.10-0.30" + sCent + "_hPtTrigger_MCPhoton");
    cout << hTriggerSamMCRec[iCen] << endl;

    for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
    {
      TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
      histFitNoUE[iCen][iptTr] = new TH1F("histFitNoUE" + sCent + sPtTrig, "histFitNoUE" + sCent + sPtTrig, nZtBin, assocZt);
      histSystErrFitNoUE[iCen][iptTr] = new TH1F("histSystErrFitNoUE" + sCent + sPtTrig, "histSystErrFitNoUE" + sCent + sPtTrig, nZtBin, assocZt);
      PlotStyle(histSystErrFitNoUE[iCen][iptTr], 20, 1.5, kBlue, "z_{T}", "(fit/away)%");
      for (Int_t izt = 0; izt < nZtBin; izt++)
      {
        TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);
        for (int iso = 1; iso < nIso; iso++)
        {
          TString sIso = Form("Iso%d", iso);
          for (int iSh = 0; iSh < nShSh; iSh++)
          {
            TString sShSh = Form("_ShSh%s", shshString[iSh].Data());
            cout << " Get azimuthal distributions " << endl;
            hdPhiSamMCRecNoUE[iCen][iso][iSh][izt][iptTr] = (TH1F *)fPlotMix[iCen]->Get("hdPhiSameMCRecNoUE" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh]);
            cout << hdPhiSamMCRecNoUE[iCen][iso][iSh][izt][iptTr] << endl;

            hdPhiMixMCRec[iCen][iso][iSh][izt][iptTr] = (TH1F *)fPlotMix[iCen]->Get("hdPhiMixMCRecMirror" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh]);
            cout << hdPhiMixMCRec[iCen][iso][iSh][izt][iptTr] << endl;
            hdPhiSamMCRec[iCen][iso][iSh][izt][iptTr] = (TH1F *)fPlotMix[iCen]->Get("hdPhiSameMCRecMirror" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh]);
            cout << hdPhiSamMCRec[iCen][iso][iSh][izt][iptTr] << endl;
            hdPhiSamMCRecNoUEZYAM[iCen][iso][iSh][izt][iptTr] = (TH1F *)fPlotZYAM[iCen]->Get("hdPhiSameMCRecNoUE" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh]);
            PlotStyle(hdPhiSamMCRecNoUE[iCen][iso][iSh][izt][iptTr], kMarkStyle[iSh], 1, kBlue, "#Delta#varphi (rad)", "1/N^{trig}dN^{charg}/d#Delta#varphi");
            PlotStyle(hdPhiSamMCRec[iCen][iso][iSh][izt][iptTr], kMarkStyle[iSh], 1, kBlack, "#Delta#varphi (rad)", "1/N^{trig}dN^{charg}/d#Delta#varphi");
            PlotStyle(hdPhiMixMCRec[iCen][iso][iSh][izt][iptTr], kMarkStyle[iSh], 1, kPink, "#Delta#varphi (rad)", "1/N^{trig}dN^{charg}/d#Delta#varphi");
            PlotStyle(hdPhiSamMCRecNoUEZYAM[iCen][iso][iSh][izt][iptTr], 25, 1, kRed, "#Delta#varphi (rad)", "1/N^{trig}dN^{charg}/d#Delta#varphi");

            hZtMixPtBin[iCen][iptTr] = (TH1F *)fPlotMix[iCen]->Get(Form("hZtPtBinMCRec_%s%s%s", sIso.Data(), sPtTrig.Data(), sNamePtTrigGen[iSh].Data()));
            hZtZYAMPtBin[iCen][iptTr] = (TH1F *)fPlotZYAM[iCen]->Get(Form("hZtPtBinMCRec_%s%s%s", sIso.Data(), sPtTrig.Data(), sNamePtTrigGen[iSh].Data()));
            PlotStyle(hZtMixPtBin[iCen][iptTr], 20, 2, kBlue, "z_{T}", "1/N^{trig}dN^{charg}/dz_{T}");
            PlotStyle(hZtZYAMPtBin[iCen][iptTr], 25, 2, kRed, "z_{T}", "1/N^{trig}dN^{charg}/dz_{T}");

            fitSamNOUE[iCen][iso][iSh][izt][iptTr] = new TF1("fitSamNOUE" + sCent + sPtTrig, "pol0", 3 * (TMath::Pi()) / 10, TMath::Pi() / 2, "R");
            hdPhiSamMCRecNoUE[iCen][iso][iSh][izt][iptTr]->Fit("fitSamNOUE" + sCent + sPtTrig, "R");
          }
        }
        histFitNoUE[iCen][iptTr]->SetBinContent(izt + 1, ((fitSamNOUE[iCen][1][0][izt][iptTr]->GetParameter(0))));
        histFitNoUE[iCen][iptTr]->SetBinError(izt + 1, ((fitSamNOUE[iCen][1][0][izt][iptTr]->GetParError(0))));
        int binMin = hdPhiSamMCRecNoUE[iCen][1][0][izt][iptTr]->FindBin(TMath::Pi() * 3 / 5);
        int binMax = hdPhiSamMCRecNoUE[iCen][1][0][izt][iptTr]->FindBin(TMath::Pi() - 0.0001);
        double integ = hdPhiSamMCRecNoUE[iCen][1][0][izt][iptTr]->Integral(binMin, binMax);

        if (integ > 0)
        {
          fitSystErrNoUE[iCen][izt][iptTr] = (abs(fitSamNOUE[iCen][1][0][izt][iptTr]->GetParameter(0)) * 4 / 2) / (integ - abs(fitSamNOUE[iCen][1][0][izt][iptTr]->GetParameter(0)) * 4 / 2);
          histSystErrFitNoUE[iCen][iptTr]->SetBinContent(izt + 1, (fitSystErrNoUE[iCen][izt][iptTr]) * 100);
          histSystErrFitNoUE[iCen][iptTr]->SetBinError(izt + 1, (fitSamNOUE[iCen][1][0][izt][iptTr]->GetParError(0) / integ) * 100);
        }

        hdPhiSamMCRecNoUEDiff[iCen][izt][iptTr] = (TH1F *)hdPhiSamMCRecNoUE[iCen][1][0][izt][iptTr]->Clone("hdPhiSameMCRecNoUEDiff" + sZtBin + sPtTrig);
        hdPhiSamMCRecNoUEDiff[iCen][izt][iptTr]->Add(hdPhiSamMCRecNoUEZYAM[iCen][1][0][izt][iptTr], -1);
      }
      fZYAMSyst->cd();
      histSystErrFitNoUE[iCen][iptTr]->Write();
      histSystErrFitNoUEXPtBin[iCen][iptTr] = (TH1F *)histSystErrFitNoUE[iCen][iptTr]->Clone("histSystErrFitNoUEXPtBin" + sCent + sPtTrig);
    }
    histSystErrFitNoUEXPtAll[iCen] = new TH1F("histSystErrFitNoUEXPtAll" + sCent, "histSystErrFitNoUEXPtAll" + sCent, nZtBin, assocZt);
    histSystErrFitNoUEXPtAll[iCen] = SumPtBinXzt(hTriggerSamMCRec[iCen], ptTrig, index1, index2, histSystErrFitNoUEXPtBin[iCen], histSystErrFitNoUEXPtAll[iCen], histPur[iCen], funcPur[iCen], 1, false);

    hZtMix[iCen] = (TH1F *)fPlotMix[iCen]->Get(Form("hZtMCRecIso1Photon%s%s", sCent.Data(), sPtAll.Data()));
    hZtZYAM[iCen] = (TH1F *)fPlotZYAM[iCen]->Get(Form("hZtMCRecIso1Photon%s%s", sCent.Data(), sPtAll.Data()));

    cout << hZtMix[iCen] << endl;
    cout << hZtZYAM[iCen] << endl;
    PlotStyle(hZtMix[iCen], 20, 2, kBlue, "z_{T}", "1/N^{trig}dN^{charg}/dz_{T}");
    PlotStyle(hZtZYAM[iCen], 25, 2, kRed, "z_{T}", "1/N^{trig}dN^{charg}/dz_{T}");
    fZYAMSyst->cd();
    histSystErrFitNoUEXPtAll[iCen]->Write();
  }

  TCanvas *cIsoClustSameNoUE_ZYAMMC[nCen][nPtTrig];
  TLegend *legIsoClustSameNoUE_ZYAMMC[nCen][nPtTrig];
  TCanvas *cIsoClustSameNoUE_ZYAMMCDiff[nCen][nPtTrig];
  TLegend *legIsoClustSameNoUE_ZYAMMCDiff[nCen][nPtTrig];
  // TF1 *fitNoUEDiff[nCen][nZtBin][nPtTrig];
  TLatex *latDphi[nCen];
  TCanvas *cIsoSameMix[nCen][nPtTrig];
  TCanvas *cIsoZtPtBinZYAM_Mix[nCen];
  TCanvas *chistErrFit[nCen];

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    TString sdirPlotXCent = Form("%s/Cen%d_%d", dirPlot.Data(), cenBins[iCen], cenBins[iCen + 1]);
    gSystem->Exec(Form("mkdir %s", sdirPlotXCent.Data()));
    cIsoZtPtBinZYAM_Mix[iCen] = new TCanvas("cIsoZtPtBinZYAM_Mix" + sCent, "cIsoZtPtBinZYAM_Mix" + sCent, 4 * 800, 2 * 600);
    cIsoZtPtBinZYAM_Mix[iCen]->Divide(4, 2);
    chistErrFit[iCen] = new TCanvas("chistErrFit" + sCent, "chistErrFit" + sCent, 3 * 800, 2 * 600);
    chistErrFit[iCen]->Divide(3, 2);
    for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
    {
      TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);

      cIsoZtPtBinZYAM_Mix[iCen]->cd(iptTr + 1);
      hZtMixPtBin[iCen][iptTr]->Draw("same");
      hZtZYAMPtBin[iCen][iptTr]->Draw("same");

      chistErrFit[iCen]->cd(iptTr + 1);
      histSystErrFitNoUE[iCen][iptTr]->SetTitle(Form("%2.0f < #it{p}_{T}^{tr} < %2.0f GeV/#it{c}", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]));
      histSystErrFitNoUE[iCen][iptTr]->Draw();

      cIsoClustSameNoUE_ZYAMMC[iCen][iptTr] = new TCanvas("cIsoClustSameNoUE_ZYAMMC" + sCent + sPtTrig, "cIsoClustSameNoUE_ZYAMMC" + sCent + sPtTrig, 4 * 800, 1 * 600);
      cIsoClustSameNoUE_ZYAMMC[iCen][iptTr]->Divide(4, 1);
      legIsoClustSameNoUE_ZYAMMC[iCen][iptTr] = LegStd(legIsoClustSameNoUE_ZYAMMC[iCen][iptTr], 0.12, 0.20, 0.12, 0.56);
      cIsoClustSameNoUE_ZYAMMCDiff[iCen][iptTr] = new TCanvas("cIsoClustSameNoUE_ZYAMMCDiff" + sCent + sPtTrig, "cIsoClustSameNoUE_ZYAMMCDiff" + sCent + sPtTrig, 4 * 800, 1 * 600);
      cIsoClustSameNoUE_ZYAMMCDiff[iCen][iptTr]->Divide(4, 1);
      legIsoClustSameNoUE_ZYAMMCDiff[iCen][iptTr] = LegStd(legIsoClustSameNoUE_ZYAMMCDiff[iCen][iptTr], 0.12, 0.24, 0.12, 0.40);
      cIsoSameMix[iCen][iptTr] = new TCanvas("cIsoSameMix" + sCent + sPtTrig, "cIsoSameMix" + sCent + sPtTrig, 4 * 800, 2 * 600);
      cIsoSameMix[iCen][iptTr]->Divide(4, 2);
      for (Int_t izt = 0; izt < nZtBin; izt++)
      {
        TString sTitle = Form("%2.0f < #it{p}_{T}^{tr} < %2.0f GeV/#it{c}, %2.2f < #it{z}_{T}^{as} < %2.2f ",
                              ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], assocZt[izt], assocZt[izt + 1]);
        // gStyle->SetPadRightMargin(0.018);
        // gStyle->SetPadLeftMargin(0.12);
        TGaxis::SetMaxDigits(2);
        cIsoClustSameNoUE_ZYAMMC[iCen][iptTr]->cd(izt + 1)->SetGridx();
        hdPhiSamMCRecNoUE[iCen][1][0][izt][iptTr]->SetTitle(sTitle);
        hdPhiSamMCRecNoUEZYAM[iCen][1][0][izt][iptTr]->SetTitle(sTitle);
        hdPhiSamMCRecNoUEZYAM[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(0.8 * hdPhiSamMCRecNoUEZYAM[iCen][1][0][izt][iptTr]->GetMinimum(), 1.2 * hdPhiSamMCRecNoUE[iCen][1][0][izt][iptTr]->GetMaximum());
        hdPhiSamMCRecNoUEZYAM[iCen][1][0][izt][iptTr]->Draw("same");
        hdPhiSamMCRecNoUE[iCen][1][0][izt][iptTr]->Draw("same");
        hdPhiSamMCRecNoUEZYAM[iCen][1][0][izt][iptTr]->Draw("same");
        // fitNoUEDiff[iCen][izt][iptTr] = new TF1("fitNoUEDiff", "pol0", 0.97, TMath::Pi() / 2, "R");
        cIsoClustSameNoUE_ZYAMMCDiff[iCen][iptTr]->cd(izt + 1)->SetGridx();
        hdPhiSamMCRecNoUEDiff[iCen][izt][iptTr]->Draw("same");
        //  hdPhiSamNoUEDiff[iCen][izt][iptTr]->SetMinimum(-1e-1);
        // hdPhiSamNoUEDiff[iCen][izt][iptTr]->Draw("same");
        // fitNoUEDiff[iCen][izt][iptTr]->Draw("same");
        cIsoSameMix[iCen][iptTr]->cd(izt + 1)->SetGridx();
        hdPhiMixMCRec[iCen][1][0][izt][iptTr]->SetTitle(sTitle);
        hdPhiMixMCRec[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(0.8 * hdPhiMixMCRec[iCen][1][0][izt][iptTr]->GetMinimum(), 1.2 * hdPhiMixMCRec[iCen][1][0][izt][iptTr]->GetMaximum());
        hdPhiMixMCRec[iCen][1][0][izt][iptTr]->Draw("same");
        hdPhiSamMCRec[iCen][1][0][izt][iptTr]->Draw("same");
      }

      cIsoClustSameNoUE_ZYAMMC[iCen][iptTr]->Print(sdirPlotXCent + "/cIsoClustSameNoUE_ZYAMMC" + sCent + sPtTrig + ".pdf");
      cIsoClustSameNoUE_ZYAMMCDiff[iCen][iptTr]->Print(sdirPlotXCent + "/cIsoClustSameNoUE_ZYAMMCDiff" + sCent + sPtTrig + ".pdf");
      cIsoSameMix[iCen][iptTr]->Print(sdirPlotXCent + "/cIsoSameMix" + sCent + sPtTrig + ".pdf");
    }
    cIsoZtPtBinZYAM_Mix[iCen]->Print(sdirPlotXCent + "/cIsoZtPtBinZYAM_Mix" + sCent + ".pdf");
    chistErrFit[iCen]->Print(sdirPlotXCent + "/cHistSystErrFit" + sCent + ".pdf");
  }

  TCanvas *cUncert = new TCanvas("cUncert", "cUncert", 2 * 800, 2 * 600);
  TCanvas *cZtMix_ZYAM = new TCanvas("cZtMix_ZYAM", "cZtMix_ZYAM", 2 * 800, 2 * 600);
  cUncert->Divide(2, 2);
  cZtMix_ZYAM->Divide(2, 2);
  TLatex *lat1[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    cZtMix_ZYAM->cd(iCen + 1);
    gPad->SetLogy();
    hZtMix[iCen]->SetTitle(" ");
    hZtMix[iCen]->Draw("same");
    hZtZYAM[iCen]->Draw("same");
    lat1[iCen] = LatexStd(lat1[iCen], 0.48, 0.8, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true);
    hZtUncertSyst[iCen] = new TH1F("hZYAMCrossCheck" + sCent, "hZYAMCrossCheck" + sCent, nZtBin, assocZt);
    for (int ibin = 1; ibin <= nZtBin; ibin++)
    {
      double bincont = abs((hZtMix[iCen]->GetBinContent(ibin) - hZtZYAM[iCen]->GetBinContent(ibin)) / hZtMix[iCen]->GetBinContent(ibin));
      hZtUncertSyst[iCen]->SetBinContent(ibin, bincont);
      hZtUncertSyst[iCen]->SetBinError(ibin, hZtMix[iCen]->GetBinError(ibin) / hZtMix[iCen]->GetBinContent(ibin));
    }
    hZtUncertSyst[iCen]->Scale(100);
    cUncert->cd(iCen + 1);
    hZtUncertSyst[iCen]->SetMarkerStyle(20);
    hZtUncertSyst[iCen]->SetMarkerColor(kAzure + 2);

    hZtUncertSyst[iCen]->SetTitle("");
    hZtUncertSyst[iCen]->SetStats(0);

    hZtUncertSyst[iCen]->Draw();
    TLatex *latStd = LatexStd(latStd, 0.220, 0.84, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true);
    fZYAMSyst->cd();
    hZtUncertSyst[iCen]->Write();
  }
  cUncert->Print(Form("%s/SystZtMix_ZYAM.pdf", dirPlot.Data()));
  cZtMix_ZYAM->Print(Form("%s/MCRecZtMix_ZYAM.pdf", dirPlot.Data()));

  TCanvas *cErrSystNoUE = new TCanvas("cErrSystNoUE", "cErrSystNoUE", 2 * 800, 2 * 600);
  cErrSystNoUE->Divide(2, 2);
  TLatex *latexGlob[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    // TGaxis::SetMaxDigits(3);
    PlotStyle(histSystErrFitNoUEXPtAll[iCen], 20, 1, kAzure + 3, "#font[12]{z_{T}}", "Uncertainty %");
    cErrSystNoUE->cd(iCen + 1);
    TGaxis::SetMaxDigits(4);
    histSystErrFitNoUEXPtAll[iCen]->SetTitle(" ");
    histSystErrFitNoUEXPtAll[iCen]->GetYaxis()->SetRangeUser(0, 100);
    histSystErrFitNoUEXPtAll[iCen]->Draw();
    latexGlob[iCen] = LatexStd(latexGlob[iCen], 0.15, 0.84, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true);
  }
  cErrSystNoUE->Print(Form("%s/MCSystErrUECen.pdf", dirPlot.Data()));
  fZYAMSyst->Close();
}

void PlotStyle(TH1F *hPlot, int kMarker, double kMarkerSize, int kColor, TString titleX, TString titleY)
{
  gStyle->SetTitleX(0.52);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111111);

  // gStyle->SetLabelSize(.5, "XY");
  gStyle->SetLineScalePS(1);
  gStyle->SetTitleAlign(23);
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
TLatex *LatexDPhi(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax)
{
  lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetTextSize(0.06);
  lat->SetNDC();
  lat->DrawLatex(xpos, ypos, Form("#it{Work in progress}"));
  lat->DrawLatex(xpos, ypos - 0.10, Form("%d-%d %% Pb-Pb, #sqrt{s_{NN}} = 5.02 TeV", cenMin, cenMax));
  return lat;
}
TLatex *LatexStd(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax, bool sPttrig)
{
  lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetTextSize(0.04);
  lat->SetNDC();
  lat->DrawLatex(xpos, ypos, Form("#it{Work in progress}"));
  lat->DrawLatex(xpos, ypos - 0.06, Form("%d-%d %% Pb-Pb, #sqrt{s_{NN}} = 5.02 TeV", cenMin, cenMax));
  if (sPttrig)
    lat->DrawLatex(xpos, ypos - 2 * 0.06, Form("%2.0f < #it{p}_{#it{T}}^{trig} < %2.0f GeV/#it{c}", ptMin, ptMax));
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

TH1F *SumPtBinXzt(TH1F *hTrigSame, Float_t PtTrigger[npt], int index1, int index2, TH1F *hzTbin[npt], TH1F *hzTbinAll, TH1F *hPur, TF1 *fPur, double systPur, Bool_t bData)
{
  Float_t intPtAll = 0;
  cout << "(1/NintPur_tot) * Sum(pur * NInt * f(Zt))" << endl;
  for (int iptTrig = 0; iptTrig < index2 - index1; iptTrig++)
  {
    cout << PtTrigger[iptTrig + index1] << endl;
    Float_t nPtRec = hTrigSame->Integral(hTrigSame->FindBin(PtTrigger[iptTrig + index1]), hTrigSame->FindBin(PtTrigger[iptTrig + index1 + 1] - 0.0001));
    Float_t pur = 1;
    if (bData)
    {
      double syst = hPur->GetBinError(hPur->FindBin(PtTrigger[iptTrig + index1 + 1] - 0.0001)) / (hPur->GetBinContent(hPur->FindBin(PtTrigger[iptTrig + index1 + 1] - 0.0001)));
      pur = fPur->Eval(hPur->GetBinCenter(hPur->FindBin(PtTrigger[iptTrig + index1 + 1] - 0.0001)));
      if (systPur == 1.1)
      {
        cout << "Upper limit" << systPur << endl;
        pur = pur * (1 + syst);
      }
      if (systPur == 0.9)
      {
        cout << "Lower limit" << systPur << endl;
        pur = pur * (1 - syst);
      }
    }

    cout << pur << endl;
    cout << hzTbin[iptTrig] << endl;
    hzTbin[iptTrig]->Scale(nPtRec * pur);
    cout << nPtRec << endl;
    hzTbinAll->Add(hzTbin[iptTrig]);
    // cout << nPtRec << endl;
    intPtAll = intPtAll + nPtRec * pur;
    // cout << nPtRec << endl;
  }
  cout << intPtAll << endl;
  hzTbinAll->Scale(1 / intPtAll);

  return hzTbinAll;
}
