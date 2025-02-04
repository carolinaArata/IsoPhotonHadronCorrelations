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
#include "../Plotting.h"

using std::cout;
using std::endl;
// Int_t nCen = 3;
// Int_t cenBins[] = {0, 30, 50, 90};
//  Int_t cenBins[] = {10, 30, 50, 90};

// int nZtBin = 10;
// double assocZt[] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20};

int nZtBin = 6;
double assocZt[] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 1.00};

int npt = 13;
float ptTrig[] = {12, 14, 16, 18, 20, 25, 30, 35, 40, 50, 60, 80};

Int_t const nIso = 2;
Int_t const nShSh = 2;

// Int_t kColorData[2][2] = {{}{}};
Int_t kMarkStyle[2] = {20, 24};
Int_t kMarkStyleIso_NotIso[2] = {21, 48};

//{2, 4, kOrange + 7, kAzure + 3, kPink - 4, kViolet - 7, kBlue + 2, kTeal - 6};
Int_t kColorMC[] = {2, 4, kOrange + 7, kAzure + 3, kPink - 4, kViolet - 7, kBlue + 2, kTeal - 6};
Int_t kColorIsoREC[] = {kAzure + 2, kOrange - 3};
Int_t kColorShSh[] = {kRed, kBlue};

// void PlotStyle(TH1F *hPlot, int kMarker, double kMarkerSize, int kColor, TString titleX, TString titleY);
// TLatex *LatexStd(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax, bool sPttrig);
// TLatex *LatexDPhi(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax);
TLatex *LatexStdIcp(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax, bool sPttrig);
// TLegend *LegStd(TLegend *leg, double xpos1, double ypos1, double xpos2, double ypos2);
TH1F *SumPtBinXzt(TH1F *hTrigSame, Float_t PtTrigger[npt], int index1, int index2, TH1F *hzTbin[npt], TH1F *hzTbinAll, TH1F *hPur, TF1 *fPur, double systPur, Bool_t bData);
TH1F *SumPtBinXzt0_30(Double_t nTrig0_30[npt], Float_t PtTrigger[npt], int index1, int index2, TH1F *hzTbin[npt], TH1F *hzTbinAll, TH1F *hPur010, TF1 *fPur010, TH1F *hPur1030, TF1 *fPur1030, Double_t nTrig0_10[npt], Double_t nTrig10_30[npt], double systPur, Bool_t bData);

void UEforClusterCheck(float ptMin = 18, float ptMax = 40, TString shshBkg = "0.40-1.00", TString dirRefFiles = "~/work/histogram/FromScratch/checkCode", bool b0_30 = true)
{
  Int_t nCen;
  std::vector<Int_t> cenBins;
  TString dirPlot;

  nCen = 4;
  cenBins.push_back(0);
  cenBins.push_back(10);
  cenBins.push_back(30);
  cenBins.push_back(50);
  cenBins.push_back(90);
  dirPlot = "~/work/histogram/Systematics_checkCode";

  // Define strings

  TString shshString[2] = {"0.10-0.30", shshBkg};
  TString sPtAll = Form("_Pt%2.0f_%2.0f", ptMin, ptMax);
  TString sHistName = "AnaPhotonHadronCorr_";

  // index pt start and stop
  int nsize = sizeof(ptTrig) / sizeof(ptTrig[0]);
  auto itr1 = find(ptTrig, ptTrig + nsize, ptMin);
  auto itr2 = find(ptTrig, ptTrig + nsize, ptMax);
  int index1 = distance(ptTrig, itr1);
  int index2 = distance(ptTrig, itr2);
  int nPtTrig = index2 - index1;
  cout << index1 << "___" << index2 << ", " << nPtTrig << endl;
  gSystem->Exec(Form("mkdir %s/SystResidualUE", dirPlot.Data()));

  // Output file
  TFile *fUEResidSyst = new TFile(Form(dirPlot + "/SystResidualUE/fUEResidSyst%s%s.root", shshBkg.Data(), sPtAll.Data()), "RECREATE");

  TH1F *histPur[nCen];
  TH1F *histPurStat[nCen];
  TF1 *funcPur[nCen];

  TH1F *hSamPtTrigger[nCen][nShSh];
  TH1F *hdPhiSamNoUE[nCen][nShSh][nZtBin][nPtTrig];
  TH1F *hdPhiSamNoUEDiff[nCen][nZtBin][nPtTrig];
  TH1F *hdPhiIsoPhoton[nCen][nZtBin][nPtTrig];
  TH1F *hdPhiSamePi0Pur[nCen][nZtBin][nPtTrig];
  TFile *fPlot[nCen];

  TF1 *fitNoUE[nCen][nShSh][nZtBin][nPtTrig];
  TH1F *histFitValueNoUE[nCen][nShSh][nPtTrig];   // fit for narrow and wide
  TH1F *histSystErrFitNoUE[nCen][nShSh][nPtTrig]; // fit normalised to awayside narrow and wide
  TH1F *histSystErrFitNoUEXSumPt[nCen][nShSh][nPtTrig];
  TH1F *hErrSystWideNoUE[nCen];
  TH1F *hErrSystNarrNoUE[nCen];
  TH1F *hErrSystIsoPhoton[nCen];
  TH1F *hErrSystMediaNoUE[nCen];

  TH1F *histFitValueMediaNoUE[nCen][nPtTrig];
  TH1F *histSystErrFitValueMediaNoUE[nCen][nPtTrig];
  TH1F *histSystErrFitValueMediaNoUEXSumPt[nCen][nPtTrig];
  TF1 *fitIsoPhoton[nCen][nZtBin][nPtTrig];
  TH1F *histFitValueIsoPhoton[nCen][nPtTrig];
  TH1F *histSystErrFitIsoPhoton[nCen][nPtTrig];
  TH1F *histSystErrFitIsoPhotonXSumPt[nCen][nPtTrig];
  TF1 *fitNoUEDiff[nCen][nZtBin][nPtTrig];
  double integIsoPhoton[nCen][nZtBin][nPtTrig];
  double errIntegIsoPhoton[nCen][nZtBin][nPtTrig];
  double integSameNoUE[nCen][nZtBin][nShSh][nPtTrig];
  double errIntegSameNoUE[nCen][nZtBin][nShSh][nPtTrig];
  double fitSystErrNoUE[nCen][nZtBin][nShSh][nPtTrig];
  double errFitSystErrNoUE[nCen][nZtBin][nShSh][nPtTrig];
  double fitSystErrIsoPhoton[nCen][nZtBin][nPtTrig];
  double errFitSystErrIsoPhoton[nCen][nZtBin][nPtTrig];
  double fitSystErrFitValueMediaNoUE[nCen][nZtBin][nPtTrig];

  TFile *fZYAMSyst = new TFile(Form("%s/crossZYAM/fZYAMSystMixed%s%s.root", dirPlot.Data(), shshBkg.Data(), sPtAll.Data()));
  TH1F *histSystMC[nCen][nPtTrig];
  TH1F *histSystMCXPtAll[nCen];

  double numbTrig[nCen][nShSh][nPtTrig];
  double numbTrig0_30[nShSh][nPtTrig];

  // Getter Data plots + make fits

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    fPlot[iCen] = new TFile(dirRefFiles + "/fPlot" + shshBkg + sCent + sPtAll + ".root");
    cout << "Get purity for centrality: " << sCent << endl;
    TFile *fPurity = new TFile("~/work/histogram/IsoPhotonHadronCorrelations/Purity_IsoSig1.5_M02Sig0.10-0.30_IsoBkg_4.0_25.0_M02Bkg0.40_2.00_LHC15o_18qr_L1MB.root");
    histPur[iCen] = (TH1F *)fPurity->Get(Form("Purity_Cen%d_R0.2_Sys", iCen));
    histPurStat[iCen] = (TH1F *)fPurity->Get(Form("Purity_Cen%d_R0.2", iCen));
    funcPur[iCen] = histPur[iCen]->GetFunction("purityFitCombinedSigmoid");
    // DataPlot
    cout << " Get Data plots + Set plots style" << endl;
    for (int iSh = 0; iSh < nShSh; iSh++)
    {
      TString sShSh = Form("_ShSh%s", shshString[iSh].Data());
      hSamPtTrigger[iCen][iSh] = (TH1F *)fPlot[iCen]->Get("AnaPhotonHadronCorr_Iso1" + sShSh + sCent + "_hPtTrigger"); // Pt trigger distribution for iso narrow and wide
      cout << hSamPtTrigger[iCen][iSh] << endl;

      for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
      {
        cout << iptTr << endl;
        TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
        cout << sPtTrig << endl;
        histFitValueNoUE[iCen][iSh][iptTr] = new TH1F("histFitValueNoUE" + sCent + sShSh + sPtTrig, "histFitValueNoUE" + sCent + sShSh + sPtTrig, nZtBin, assocZt);
        histSystErrFitNoUE[iCen][iSh][iptTr] = new TH1F("histSystErrFitNoUE" + sCent + sShSh + sPtTrig, "histSystErrFitNoUE" + sCent + sShSh + sPtTrig, nZtBin, assocZt);
        for (Int_t izt = 0; izt < nZtBin; izt++)
        {
          TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);
          hdPhiSamNoUE[iCen][iSh][izt][iptTr] = (TH1F *)fPlot[iCen]->Get("hdPhiSameNoUEIso1" + sShSh + sZtBin + sPtTrig);
          PlotStyle(hdPhiSamNoUE[iCen][iSh][izt][iptTr], kMarkStyle[iSh], 1.5, kCyan + 2, 0, "#Delta#varphi (rad)", "1/N^{trig}dN^{charg}/d#Delta#varphi", false);
          cout << hdPhiSamNoUE[iCen][iSh][izt][iptTr] << endl;
          cout << "Compute fits UE" << endl;
          fitNoUE[iCen][iSh][izt][iptTr] = new TF1(("fitNoUE" + sShSh), "pol0", 3 * (TMath::Pi()) / 10, TMath::Pi() / 2, "R");
          fitNoUE[iCen][iSh][izt][iptTr]->SetLineColor(kColorShSh[iSh]);
          hdPhiSamNoUE[iCen][iSh][izt][iptTr]->Fit(("fitNoUE" + sShSh), "R");
          histFitValueNoUE[iCen][iSh][iptTr]->SetBinContent(izt + 1, fitNoUE[iCen][iSh][izt][iptTr]->GetParameter(0));
          histFitValueNoUE[iCen][iSh][iptTr]->SetBinError(izt + 1, fitNoUE[iCen][iSh][izt][iptTr]->GetParError(0));

          integSameNoUE[iCen][iSh][izt][iptTr] = hdPhiSamNoUE[iCen][iSh][izt][iptTr]->IntegralAndError(hdPhiSamNoUE[iCen][iSh][izt][iptTr]->FindBin(TMath::Pi() * 3 / 5), hdPhiSamNoUE[iCen][iSh][izt][iptTr]->FindBin(TMath::Pi() - 0.0001), errIntegSameNoUE[iCen][iSh][izt][iptTr]);

          if (integSameNoUE[iCen][iSh][izt][iptTr] > 0)
          {
            fitSystErrNoUE[iCen][iSh][izt][iptTr] = (abs(fitNoUE[iCen][iSh][izt][iptTr]->GetParameter(0)) * 4 / 2) / (integSameNoUE[iCen][iSh][izt][iptTr] - (abs(fitNoUE[iCen][iSh][izt][iptTr]->GetParameter(0)) * 4 / 2));
            histSystErrFitNoUE[iCen][iSh][iptTr]->SetBinContent(izt + 1, fitSystErrNoUE[iCen][iSh][izt][iptTr] * 100);
            histSystErrFitNoUE[iCen][iSh][iptTr]->SetBinError(izt + 1, (fitNoUE[iCen][iSh][izt][iptTr]->GetParError(0) / integSameNoUE[iCen][iSh][izt][iptTr]) * 100);
          }
        }
        PlotStyle(histFitValueNoUE[iCen][iSh][iptTr], kMarkStyle[iSh], 1.5, kBlack, 0, "z_{T}", "fit value", false);
        PlotStyle(histSystErrFitNoUE[iCen][iSh][iptTr], kMarkStyle[iSh], 1.5, kBlue, 0, "z_{T}", "(fit value / away side) %", false);
        histSystErrFitNoUEXSumPt[iCen][iSh][iptTr] = (TH1F *)histSystErrFitNoUE[iCen][iSh][iptTr]->Clone("histSystErrFitNoUEXSumPt" + sShSh + sCent + sPtTrig);

        fUEResidSyst->cd();
        histSystErrFitNoUE[iCen][iSh][iptTr]->Write();
      }
    }

    hErrSystNarrNoUE[iCen] = new TH1F(Form("hErrSystNarrNoUE%s%s", sCent.Data(), sPtAll.Data()), Form("hErrSystNarrNoUE%s%s", sCent.Data(), sPtAll.Data()), nZtBin, assocZt);
    hErrSystWideNoUE[iCen] = new TH1F(Form("hErrSystWideNoUE%s%s", sCent.Data(), sPtAll.Data()), Form("hErrSystWideNoUE%s%s", sCent.Data(), sPtAll.Data()), nZtBin, assocZt);

    cout << "NARROW" << endl;
    hErrSystNarrNoUE[iCen] = SumPtBinXzt(hSamPtTrigger[iCen][0], ptTrig, index1, index2, histSystErrFitNoUEXSumPt[iCen][0], hErrSystNarrNoUE[iCen], histPur[iCen], funcPur[iCen], 1, true);
    cout << "pippo" << endl;
    PlotStyle(hErrSystNarrNoUE[iCen], kMarkStyle[0], 1, kBlue, 0, "z_{T}", "(fit value / away side) %", false);
    cout << "WIDE" << endl;
    hErrSystWideNoUE[iCen] = SumPtBinXzt(hSamPtTrigger[iCen][1], ptTrig, index1, index2, histSystErrFitNoUEXSumPt[iCen][1], hErrSystWideNoUE[iCen], histPur[iCen], funcPur[iCen], 1, false);
    PlotStyle(hErrSystWideNoUE[iCen], kMarkStyle[1], 1, kBlue, 0, "z_{T}", "(fit value / away side) %", false);

    histSystMCXPtAll[iCen] = (TH1F *)fZYAMSyst->Get(Form("histSystErrFitNoUEXPtAll_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1])); // ZYAM from MC
    cout << histSystMCXPtAll[iCen] << endl;
    PlotStyle(histSystMCXPtAll[iCen], 21, 1, kOrange + 7, 0, "#font[12]{z_{T}}", "Uncertainty %", false);
    fUEResidSyst->cd();
    hErrSystNarrNoUE[iCen]->Write();
    hErrSystWideNoUE[iCen]->Write();
    histSystMCXPtAll[iCen]->Write();
    for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
    {
      TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
      // Fit histograms data
      histFitValueMediaNoUE[iCen][iptTr] = new TH1F("histFitValueMediaNoUE" + sCent + sPtTrig, "histFitValueMediaNoUE" + sCent + sPtTrig, nZtBin, assocZt);
      histSystErrFitValueMediaNoUE[iCen][iptTr] = new TH1F("histSystErrFitValueMediaNoUE" + sCent + sPtTrig, "histSystErrFitValueMediaNoUE" + sCent + sPtTrig, nZtBin, assocZt);
      histFitValueIsoPhoton[iCen][iptTr] = new TH1F("histFitValueIsoPhoton" + sCent + sPtTrig, "histFitValueIsoPhoton" + sCent + sPtTrig, nZtBin, assocZt);
      histSystErrFitIsoPhoton[iCen][iptTr] = new TH1F("histSystErrFitIsoPhoton" + sCent + sPtTrig, "histSystErrFitIsoPhoton" + sCent + sPtTrig, nZtBin, assocZt);
      // Get fit and syst erro histograms MC
      histSystMC[iCen][iptTr] = (TH1F *)fZYAMSyst->Get("histSystErrFitNoUE" + sCent + sPtTrig);
      cout << histSystMC[iCen][iptTr] << endl;
      PlotStyle(histSystMC[iCen][iptTr], 21, 1.5, kOrange + 7, 0, "#font[12]{z_{T}}", "Uncertainty %", false);
      for (Int_t izt = 0; izt < nZtBin; izt++)
      {
        TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);
        hdPhiSamePi0Pur[iCen][izt][iptTr] = (TH1F *)fPlot[iCen]->Get(Form("hdPhiSamePi0PurIso1%s_%s", sZtBin.Data(), sPtTrig.Data()));
        // DeltaPhi difference and Fit
        hdPhiSamNoUEDiff[iCen][izt][iptTr] = (TH1F *)hdPhiSamNoUE[iCen][0][izt][iptTr]->Clone("hdPhiSameNoUEDiff" + sZtBin + sPtTrig);
        hdPhiSamNoUEDiff[iCen][izt][iptTr]->Add(hdPhiSamePi0Pur[iCen][izt][iptTr], -1);
        fitNoUEDiff[iCen][izt][iptTr] = new TF1("fitNoUEDiff", "pol0", 3 * (TMath::Pi()) / 10, TMath::Pi() / 2, "R");
        hdPhiSamNoUEDiff[iCen][izt][iptTr]->Fit("fitNoUEDiff", "R");
        // DeltaPhi IsoPhoton and Fit
        hdPhiIsoPhoton[iCen][izt][iptTr] = (TH1F *)fPlot[iCen]->Get(Form("hdPhiIso1Photon%s_%s", sZtBin.Data(), sPtTrig.Data()));
        integIsoPhoton[iCen][izt][iptTr] = hdPhiIsoPhoton[iCen][izt][iptTr]->IntegralAndError(hdPhiIsoPhoton[iCen][izt][iptTr]->FindBin(TMath::Pi() * 3 / 5 + 0.00001), hdPhiIsoPhoton[iCen][izt][iptTr]->FindBin(TMath::Pi()), errIntegIsoPhoton[iCen][izt][iptTr]);
        fitIsoPhoton[iCen][izt][iptTr] = new TF1("fitIsoPhotonCorr", "pol0", 3 * (TMath::Pi()) / 10, TMath::Pi() / 2, "R");
        hdPhiIsoPhoton[iCen][izt][iptTr]->Fit("fitIsoPhotonCorr", "R");
        histFitValueIsoPhoton[iCen][iptTr]->SetBinContent(izt + 1, (fitIsoPhoton[iCen][izt][iptTr]->GetParameter(0)));
        histFitValueIsoPhoton[iCen][iptTr]->SetBinError(izt + 1, (fitIsoPhoton[iCen][izt][iptTr]->GetParError(0)));
        if (integIsoPhoton[iCen][izt][iptTr] > 0)
        {
          fitSystErrIsoPhoton[iCen][izt][iptTr] = ((abs(fitIsoPhoton[iCen][izt][iptTr]->GetParameter(0)) * 4 / 2) / (integIsoPhoton[iCen][izt][iptTr] - (abs(fitIsoPhoton[iCen][izt][iptTr]->GetParameter(0)) * 4 / 2)));
          histSystErrFitIsoPhoton[iCen][iptTr]->SetBinContent(izt + 1, (fitSystErrIsoPhoton[iCen][izt][iptTr]) * 100);
          errFitSystErrIsoPhoton[iCen][izt][iptTr] = sqrt((4 * integIsoPhoton[iCen][izt][iptTr] * integIsoPhoton[iCen][izt][iptTr] * (fitIsoPhoton[iCen][izt][iptTr]->GetParError(0)) * (fitIsoPhoton[iCen][izt][iptTr]->GetParError(0)) - 4 * (fitIsoPhoton[iCen][izt][iptTr]->GetParameter(0)) * (fitIsoPhoton[iCen][izt][iptTr]->GetParameter(0)) * errIntegIsoPhoton[iCen][izt][iptTr] * errIntegIsoPhoton[iCen][izt][iptTr]) / (pow((integIsoPhoton[iCen][izt][iptTr] - 2 * fitIsoPhoton[iCen][izt][iptTr]->GetParameter(0)), 4)));
          // histSystErrFitIsoPhoton[iCen][iptTr]->SetBinError(izt + 1, (fitIsoPhoton[iCen][izt][iptTr]->GetParError(0) / integIsoPhoton[iCen][izt][iptTr]) * 100);
          histSystErrFitIsoPhoton[iCen][iptTr]->SetBinError(izt + 1, (errFitSystErrIsoPhoton[iCen][izt][iptTr] / integIsoPhoton[iCen][izt][iptTr]) * 100);
        }

        // histogram with media of iso wide and iso narrow fit UE range
        double mediaFitBinContent = (histFitValueNoUE[iCen][0][iptTr]->GetBinContent(izt + 1) + histFitValueNoUE[iCen][1][iptTr]->GetBinContent(izt + 1)) / 2;
        cout << "zt: " << izt << " bin content: " << mediaFitBinContent << endl;
        double mediaFitBinError = sqrt(((histFitValueNoUE[iCen][0][iptTr]->GetBinError(izt + 1) * histFitValueNoUE[iCen][0][iptTr]->GetBinError(izt + 1)) + (histFitValueNoUE[iCen][1][iptTr]->GetBinError(izt + 1) * histFitValueNoUE[iCen][1][iptTr]->GetBinError(izt + 1))) / 2);
        histFitValueMediaNoUE[iCen][iptTr]->SetBinContent(izt + 1, mediaFitBinContent);
        histFitValueMediaNoUE[iCen][iptTr]->SetBinError(izt + 1, mediaFitBinError);
        if (integIsoPhoton[iCen][izt][iptTr] > 0)
        {
          fitSystErrFitValueMediaNoUE[iCen][izt][iptTr] = (mediaFitBinContent * 4 / 2) / (integIsoPhoton[iCen][izt][iptTr] - mediaFitBinContent * 4 / 2);
          histSystErrFitValueMediaNoUE[iCen][iptTr]->SetBinContent(izt + 1, 100 * fitSystErrFitValueMediaNoUE[iCen][izt][iptTr]);
        }
        // Plotting style
        PlotStyle(hdPhiSamePi0Pur[iCen][izt][iptTr], 24, 1.6, kCyan + 2, 0, "#Delta#varphi (rad)", "1/N^{trig}dN^{charg}/d#Delta#varphi", false);
        PlotStyle(hdPhiSamNoUEDiff[iCen][izt][iptTr], 20, 1.6, kBlue + 2, 0, "#Delta#varphi (rad)", "#Delta#varphi(cluster^{narrow}) - (1-P)#cdot#Delta#varphi(cluster^{wide})", false);
        PlotStyle(hdPhiIsoPhoton[iCen][izt][iptTr], 20, 1.6, kViolet + 2, 0, "#Delta#varphi (rad)", "1/N^{trig}dN^{charg}/d#Delta#varphi", false);
      }

      cout << "pippo" << endl;
      histSystErrFitIsoPhotonXSumPt[iCen][iptTr] = (TH1F *)histSystErrFitIsoPhoton[iCen][iptTr]->Clone("histSystErrFitIsoPhotonXSumPt" + sCent + sPtTrig);
      cout << histSystErrFitIsoPhotonXSumPt[iCen][iptTr] << endl;
      histSystErrFitValueMediaNoUEXSumPt[iCen][iptTr] = (TH1F *)histSystErrFitValueMediaNoUE[iCen][iptTr]->Clone("histSystErrFitValueMediaNoUEXSumPt" + sCent + sPtTrig);
      PlotStyle(histFitValueIsoPhoton[iCen][iptTr], 20, 1, kBlack, 0, "z_{#rm{T}}", "fit value", false);
      PlotStyle(histSystErrFitIsoPhoton[iCen][iptTr], 20, 1, kBlack, 0, "z_{#rm{T}}", "(fit value / away side) %", false);
      PlotStyle(histFitValueMediaNoUE[iCen][iptTr], 20, 1, kBlack, 0, "z_{#rm{T}}", "fit value", false);
      PlotStyle(histSystErrFitValueMediaNoUE[iCen][iptTr], 20, 1, kBlack, 0, "z_{#rm T}", "(fit value / away side) %", false);
      fUEResidSyst->cd();
      histSystErrFitValueMediaNoUE[iCen][iptTr]->Write();
    }
    hErrSystIsoPhoton[iCen] = new TH1F(Form("hErrSystIsoPhoton%s%s", sCent.Data(), sPtAll.Data()), Form("hErrSystIsoPhoton%s%s", sCent.Data(), sPtAll.Data()), nZtBin, assocZt);
    hErrSystIsoPhoton[iCen] = SumPtBinXzt(hSamPtTrigger[iCen][0], ptTrig, index1, index2, histSystErrFitIsoPhotonXSumPt[iCen], hErrSystIsoPhoton[iCen], histPur[iCen], funcPur[iCen], 1, true);
    PlotStyle(hErrSystIsoPhoton[iCen], 20, 1, kBlack, 0, "z_{#rm{T}}", "(fit value / away side) %", false);
    hErrSystMediaNoUE[iCen] = new TH1F(Form("hErrSystMediaNoUE%s%s", sCent.Data(), sPtAll.Data()), Form("hErrSystMediaNoUE%s%s", sCent.Data(), sPtAll.Data()), nZtBin, assocZt);
    hErrSystMediaNoUE[iCen] = SumPtBinXzt(hSamPtTrigger[iCen][0], ptTrig, index1, index2, histSystErrFitValueMediaNoUEXSumPt[iCen], hErrSystMediaNoUE[iCen], histPur[iCen], funcPur[iCen], 1, true);
    PlotStyle(hErrSystMediaNoUE[iCen], 20, 1, kViolet - 1, 0, "z_{#rm{T}}", "(fit value / away side) %", false);
  }

  // Plot DeltaPhi distributions and corresponding fits
  TCanvas *cIsoClustNoUE[nCen][nPtTrig];
  TLegend *legfitNoUE[nCen][nZtBin][nPtTrig];
  TCanvas *cNoUEDiff[nCen][nPtTrig];
  TLegend *legfitNoUEDiff[nCen][nZtBin][nPtTrig];
  TCanvas *cIsoPhoton[nCen][nPtTrig];
  TLegend *legfitIsoPhoton[nCen][nZtBin][nPtTrig];

  // Plot fits
  TCanvas *cOnlyFitIsoClustNoUE[nCen];
  TCanvas *cSystErrIsoClustNoUE[nCen];
  TCanvas *cOnlyFitIsoClustMediaNoUE[nCen];
  TCanvas *cSystErrIsoClustMediaNoUE[nCen];
  TCanvas *cOnlyFitIsoPhoton[nCen];
  TCanvas *cSystErrIsoPhoton[nCen];

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    TString sdirPlotXCent = Form("%s/SystResidualUE/Cen%d_%d", dirPlot.Data(), cenBins[iCen], cenBins[iCen + 1]);
    gSystem->Exec(Form("mkdir %s", sdirPlotXCent.Data()));

    cout << "Azimuthal plots and corresponding fits" << endl;
    for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
    {
      TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
      cIsoClustNoUE[iCen][iptTr] = canvasStd("cIsoClustNoUE" + sCent + sPtTrig, 3, 2);
      cNoUEDiff[iCen][iptTr] = canvasStd("cNoUEDiff" + sCent + sPtTrig, 3, 2);
      cIsoPhoton[iCen][iptTr] = canvasStd("cIsoPhoton" + sCent + sPtTrig, 3, 2);
      for (Int_t izt = 0; izt < nZtBin; izt++)
      {
        TString sTitle = Form(" %2.2f < #it{z}_{T}< %2.2f ", assocZt[izt], assocZt[izt + 1]);
        cIsoClustNoUE[iCen][iptTr]->cd(izt + 1);
        // cIsoClustNoUE[iCen][iptTr]->cd(izt + 1)->SetRightMargin(0.0001);
        TGaxis::SetMaxDigits(2);
        gStyle->SetPadRightMargin(0.05);
        gStyle->SetPadLeftMargin(0.18);
        gStyle->SetPadBottomMargin(0.15);
        gStyle->SetTitleX(0.56);
        hdPhiSamNoUE[iCen][0][izt][iptTr]->SetTitle(sTitle);
        hdPhiSamNoUE[iCen][0][izt][iptTr]->GetXaxis()->SetTitleSize(0.06);
        hdPhiSamNoUE[iCen][0][izt][iptTr]->GetYaxis()->SetTitleSize(0.06);
        hdPhiSamNoUE[iCen][0][izt][iptTr]->GetXaxis()->SetLabelSize(0.05);
        hdPhiSamNoUE[iCen][0][izt][iptTr]->GetYaxis()->SetLabelSize(0.05);
        hdPhiSamNoUE[iCen][0][izt][iptTr]->Draw("same");
        // hdPhiSamNoUE[iCen][1][izt][iptTr]->Draw("same");
        legfitNoUE[iCen][izt][iptTr] = LegStd(legfitNoUE[iCen][izt][iptTr], 0.25, 0.7, 0.6, 0.85);
        legfitNoUE[iCen][izt][iptTr]->SetTextSize(0.05);
        legfitNoUE[iCen][izt][iptTr]->AddEntry(fitNoUE[iCen][0][izt][iptTr], Form("clust^{iso}_{narr} p0 = %0.5f#pm%0.5f", fitNoUE[iCen][0][izt][iptTr]->GetParameter(0), fitNoUE[iCen][0][izt][iptTr]->GetParError(0)), "lp");
        // legfitNoUE[iCen][izt][iptTr]->AddEntry(fitNoUE[iCen][1][izt][iptTr], Form("clust^{iso}_{wide} p0 = %0.5f#pm%0.5f", fitNoUE[iCen][1][izt][iptTr]->GetParameter(0), fitNoUE[iCen][1][izt][iptTr]->GetParError(0)), "lp");
        legfitNoUE[iCen][izt][iptTr]->Draw("same");

        cNoUEDiff[iCen][iptTr]->cd(izt + 1);

        // TGaxis::SetMaxDigits(2);
        hdPhiSamNoUEDiff[iCen][izt][iptTr]->SetTitle(sTitle);
        hdPhiSamNoUEDiff[iCen][izt][iptTr]->Draw("same");
        legfitNoUEDiff[iCen][izt][iptTr] = LegStd(legfitNoUEDiff[iCen][izt][iptTr], 0.12, 0.75, 0.6, 0.85);
        legfitNoUEDiff[iCen][izt][iptTr]->AddEntry(fitNoUEDiff[iCen][izt][iptTr], Form("clust^{iso}_{narr} p0 = %f#pm%f", fitNoUEDiff[iCen][izt][iptTr]->GetParameter(0), fitNoUEDiff[iCen][izt][iptTr]->GetParError(0)), "lp");
        legfitNoUEDiff[iCen][izt][iptTr]->Draw("same");

        cIsoPhoton[iCen][iptTr]->cd(izt + 1);
        // TGaxis::SetMaxDigits(2);

        gStyle->SetPadRightMargin(0.05);
        gStyle->SetPadLeftMargin(0.18);
        gStyle->SetPadBottomMargin(0.15);
        gStyle->SetTitleX(0.56);
        hdPhiIsoPhoton[iCen][izt][iptTr]->SetTitle(sTitle);
        hdPhiIsoPhoton[iCen][izt][iptTr]->GetXaxis()->SetTitleSize(0.06);
        hdPhiIsoPhoton[iCen][izt][iptTr]->GetYaxis()->SetTitleSize(0.06);
        hdPhiIsoPhoton[iCen][izt][iptTr]->GetXaxis()->SetLabelSize(0.05);
        hdPhiIsoPhoton[iCen][izt][iptTr]->GetYaxis()->SetLabelSize(0.05);
        hdPhiIsoPhoton[iCen][izt][iptTr]->Draw("same");
        legfitIsoPhoton[iCen][izt][iptTr] = LegStd(legfitIsoPhoton[iCen][izt][iptTr], 0.25, 0.7, 0.6, 0.85);
        legfitIsoPhoton[iCen][izt][iptTr]->SetTextSize(0.05);
        legfitIsoPhoton[iCen][izt][iptTr]->AddEntry(fitIsoPhoton[iCen][izt][iptTr], Form("clust^{iso}_{narr} p0 = %f#pm%f", fitIsoPhoton[iCen][izt][iptTr]->GetParameter(0), fitIsoPhoton[iCen][izt][iptTr]->GetParError(0)), "lp");
        legfitIsoPhoton[iCen][izt][iptTr]->Draw("same");
      }
      cIsoClustNoUE[iCen][iptTr]->Print(dirPlot + Form("/SystResidualUE/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/cIsoClustNoUE" + sCent + sPtTrig + ".pdf");
      cNoUEDiff[iCen][iptTr]->Print(dirPlot + Form("/SystResidualUE/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/cNoUEDiff" + sCent + sPtTrig + ".pdf");
      cIsoPhoton[iCen][iptTr]->Print(dirPlot + Form("/SystResidualUE/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/cIsoPhoton" + sCent + sPtTrig + ".pdf");
    }

    cout << "Plots fit values" << endl;
    cOnlyFitIsoClustNoUE[iCen] = canvasStd("cOnlyFitIsoClustNoUE" + sCent, 4, 2);
    cSystErrIsoClustNoUE[iCen] = canvasStd("cSystErrIsoClustNoUE" + sCent, 4, 2);
    cOnlyFitIsoClustMediaNoUE[iCen] = canvasStd("cOnlyFitIsoClustMediaNoUE" + sCent, 4, 2);
    cSystErrIsoClustMediaNoUE[iCen] = canvasStd("cSystErrIsoClustMediaNoUE" + sCent, 4, 2);
    cOnlyFitIsoPhoton[iCen] = canvasStd("cOnlyFitIsoPhoton" + sCent, 4, 2);
    cSystErrIsoPhoton[iCen] = canvasStd("cSystErrIsoPhoton" + sCent, 4, 2);
    for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
    {
      TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
      cOnlyFitIsoClustNoUE[iCen]->cd(iptTr + 1);
      histFitValueNoUE[iCen][0][iptTr]->SetTitle(Form("%2.0f < #it{p}_{T}^{tr} < %2.0f GeV/#it{c} ", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]));
      histFitValueNoUE[iCen][0][iptTr]->Draw("same");
      histFitValueNoUE[iCen][1][iptTr]->Draw("same");
      cSystErrIsoClustNoUE[iCen]->cd(iptTr + 1);
      histSystErrFitNoUE[iCen][0][iptTr]->SetTitle(Form("%2.0f < #it{p}_{T}^{tr} < %2.0f GeV/#it{c} ", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]));
      histSystErrFitNoUE[iCen][0][iptTr]->Draw("same");
      histSystErrFitNoUE[iCen][1][iptTr]->Draw("same");
      histSystMC[iCen][iptTr]->Draw("same");

      cOnlyFitIsoClustMediaNoUE[iCen]->cd(iptTr + 1);
      TGaxis::SetMaxDigits(4);
      histFitValueMediaNoUE[iCen][iptTr]->SetTitle(Form("%2.0f < #it{p}_{T}^{tr} < %2.0f GeV/#it{c} ", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]));
      histFitValueMediaNoUE[iCen][iptTr]->Draw("same");

      cSystErrIsoClustMediaNoUE[iCen]->cd(iptTr + 1);
      histSystErrFitValueMediaNoUE[iCen][iptTr]->SetTitle(Form("%2.0f < #it{p}_{T}^{tr} < %2.0f GeV/#it{c} ", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]));
      histSystErrFitValueMediaNoUE[iCen][iptTr]->SetMaximum(50);
      histSystErrFitValueMediaNoUE[iCen][iptTr]->SetMinimum(-10);

      histSystErrFitValueMediaNoUE[iCen][iptTr]->Draw("same");
      histSystMC[iCen][iptTr]->Draw("same");

      cOnlyFitIsoPhoton[iCen]->cd(iptTr + 1);
      histFitValueIsoPhoton[iCen][iptTr]->SetTitle(Form("%2.0f < #it{p}_{T}^{tr} < %2.0f GeV/#it{c} ", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]));
      histFitValueIsoPhoton[iCen][iptTr]->Draw();
      cSystErrIsoPhoton[iCen]->cd(iptTr + 1);
      histSystErrFitIsoPhoton[iCen][iptTr]->SetTitle(Form("%2.0f < #it{p}_{T}^{tr} < %2.0f GeV/#it{c} ", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]));
      histSystErrFitIsoPhoton[iCen][iptTr]->Draw();
    }
    cOnlyFitIsoClustNoUE[iCen]->Print(dirPlot + Form("/SystResidualUE/IsoClustNoUEFitCen%d_%d%s.pdf", cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
    cSystErrIsoClustNoUE[iCen]->Print(dirPlot + Form("/SystResidualUE/IsoClustNoUESystErrCen%d_%d%s.pdf", cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
    cOnlyFitIsoClustMediaNoUE[iCen]->Print(dirPlot + Form("/SystResidualUE/IsoClustNoUEFitMediaCen%d_%d%s.pdf", cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
    cSystErrIsoClustMediaNoUE[iCen]->Print(dirPlot + Form("/SystResidualUE/IsoClustNoUESystErrMediaCen%d_%d%s.pdf", cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
    cOnlyFitIsoPhoton[iCen]->Print(dirPlot + Form("/SystResidualUE/IsoPhotonFitCen%d_%d%s.pdf", cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
    cSystErrIsoPhoton[iCen]->Print(dirPlot + Form("/SystResidualUE/IsoPhotonSystErrCen%d_%d%s.pdf", cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
  }

  TCanvas *cSystErrAllPt_NarrWideMC = canvasStd("cSystErrAllPt_NarrWideMC", 2, 2);
  TLatex *latSystErrAllCen[nCen];
  TLegend *legSystNarrWideMC[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    PlotStyle(hErrSystNarrNoUE[iCen], 20, 1.5, kAzure + 7, 0, "#it{z}_{T}", "Uncertainty %", false);
    hErrSystNarrNoUE[iCen]->SetTitle("");
    cSystErrAllPt_NarrWideMC->cd(iCen + 1);
    hErrSystNarrNoUE[iCen]->SetMaximum(100);
    hErrSystNarrNoUE[iCen]->SetMinimum(-0.10);
    hErrSystNarrNoUE[iCen]->Draw("same");
    hErrSystWideNoUE[iCen]->Draw("same");
    // histSystMCXPtAll[iCen]->Draw("same");
    latSystErrAllCen[iCen] = LatexStdSyst(latSystErrAllCen[iCen - 1], 0.450, 0.84, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true, " ");
    legSystNarrWideMC[iCen] = LegStd(legSystNarrWideMC[iCen], 0.55, 0.55, 0.8, 0.65);
    legSystNarrWideMC[iCen]->AddEntry(hErrSystNarrNoUE[iCen], "Data #it{N}^{iso}_{narrow}", "lp");
    legSystNarrWideMC[iCen]->AddEntry(hErrSystWideNoUE[iCen], "Data #it{N}^{iso}_{wide}", "lp");
    // legSystNarrWideMC[iCen]->AddEntry(histSystMCXPtAll[iCen], "MC GJ cluster^{iso}_{narr}", "lp");
    //  hErrSystIsoPhoton[iCen]->Draw("same");
    //  hErrSystMediaNoUE[iCen]->Draw("same");
  }
  cSystErrAllPt_NarrWideMC->cd(2);
  legSystNarrWideMC[0]->Draw("same");
  cSystErrAllPt_NarrWideMC->Print(dirPlot + Form("/SystResidualUE/cSystErrPtAll%s_NarrWide.pdf", sPtAll.Data()));

  TCanvas *cSystErrAllPt_NarrAllCen = canvasStd("cSystErrAllPt_NarrAllCen", 2, 2);
  TLatex *latSystErrAllCen_Narr[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    cSystErrAllPt_NarrAllCen->cd(iCen + 1);
    hErrSystNarrNoUE[iCen]->SetMaximum(100);
    hErrSystNarrNoUE[iCen]->SetMinimum(-0.10);
    hErrSystNarrNoUE[iCen]->Draw("same");
    latSystErrAllCen_Narr[iCen] = LatexStdSyst(latSystErrAllCen_Narr[iCen - 1], 0.450, 0.84, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true, "");
    // hErrSystWideNoUE[iCen]->Draw("same");
    // histSystMCXPtAll[iCen]->Draw("same");
    //   hErrSystIsoPhoton[iCen]->Draw("same");
    //  hErrSystMediaNoUE[iCen]->Draw("same");
  }
  cSystErrAllPt_NarrAllCen->Print(dirPlot + Form("/SystResidualUE/cSystErrPtAll%s_Narr.pdf", sPtAll.Data()));

  TCanvas *cSystErrAllPt_Narr[nCen];
  TLatex *latSystErrAllPt_Narr[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    cSystErrAllPt_Narr[iCen] = canvasStd(Form("cSystErrAllPt_Narr_Cen%d_%d%s.pdf", cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()), 1, 1);
    cSystErrAllPt_Narr[iCen]->cd();
    hErrSystNarrNoUE[iCen]->GetYaxis()->SetRangeUser(-0.10, 60);
    hErrSystNarrNoUE[iCen]->Draw("pl");
    latSystErrAllPt_Narr[iCen] = LatexStdSyst(latSystErrAllPt_Narr[iCen - 1], 0.450, 0.84, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, true, " ");
    cSystErrAllPt_Narr[iCen]->Print(dirPlot + Form("/SystResidualUE/SystErrNarrCen%d_%d.pdf", cenBins[iCen], cenBins[iCen + 1]));
  }

  TF1 *fitNoUENarrCl[nCen][nZtBin][nPtTrig];
  TH1F *histSystErrFitNoUENarrClust[nCen][nPtTrig];
  TH1F *histOnlyFitNoUENarrClust[nCen][nPtTrig];

  TF1 *fitNoUEWideCl[nCen][nZtBin][nPtTrig];
  TH1F *histSystErrFitNoUEWideClust[nCen][nPtTrig];
  TH1F *histSystErrFitNoUEWideClustXPtAll[nCen][nPtTrig];
  TH1F *histSystErrFitNoUENarrClustXPtAll[nCen][nPtTrig];
  TH1F *histOnlyFitNoUEWideClust[nCen][nPtTrig];
  TH1F *hTriggerSam[nCen][nShSh];

  TH1F *histSystErrFitNoUE0_10[nShSh][nPtTrig];
  TH1F *histSystErrFitNoUE10_30[nShSh][nPtTrig];
  TH1F *histSystErrFitNoUE0_30[nShSh][nPtTrig];
  TH1F *histSystErrFitNoUE0_30XPt[nShSh][nPtTrig];

  TH1F *histSystErrFitNoUENarr0_30PtRange = new TH1F(Form("histSystErrFitNoUENarr0_30PtRange%s", sPtAll.Data()), Form("histSystErrFitNoUENarr0_30PtRange%s", sPtAll.Data()), nZtBin, assocZt);
  TH1F *histSystErrFitNoUENarr0_30PtRangeFitTrend;
  TF1 *fitConst;
  TF1 *fitExpo;

  if (b0_30)
  {
    for (int iCen = 0; iCen < 2; iCen++)
    {
      cout << iCen << " " << cenBins[iCen] << "-" << cenBins[iCen + 1] << endl;

      TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
      cout << "Getter Pt Trig distrib centrality: " << iCen << endl;
      for (int iSh = 0; iSh < nShSh; iSh++)
      {
        TString sShSh = Form("_ShSh%s", shshString[iSh].Data());
        hTriggerSam[iCen][iSh] = (TH1F *)fPlot[iCen]->Get("AnaPhotonHadronCorr_Iso1" + sShSh + sCent + "_hPtTrigger");

        cout << fPlot[iCen]->GetName() << endl;
        for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
        {
          numbTrig[iCen][iSh][iptTr] = hTriggerSam[iCen][iSh]->Integral(hTriggerSam[iCen][iSh]->FindBin(ptTrig[index1 + iptTr]), hTriggerSam[iCen][iSh]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001));
        }
      }
    }
  }
  for (int iSh = 0; iSh < nShSh; iSh++)
  {
    for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
    {
      numbTrig0_30[iSh][iptTr] = numbTrig[0][iSh][iptTr] + numbTrig[1][iSh][iptTr];

      cout << "Pt: " << ptTrig[index1 + iptTr] << "Same 0-10\%: " << numbTrig[0][iSh][iptTr] << endl;
      cout << "Pt: " << ptTrig[index1 + iptTr] << "Same 10-30\%: " << numbTrig[1][iSh][iptTr] << endl;
      cout << "Pt: " << ptTrig[index1 + iptTr] << "Same 0-30\%: " << numbTrig0_30[iSh][iptTr] << endl;
    }
  }

  for (int iSh = 0; iSh < nShSh; iSh++)
  {
    TString sShSh = Form("_ShSh%s", shshString[iSh].Data());
    for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
    {
      TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
      histSystErrFitNoUE0_10[iSh][iptTr] = (TH1F *)fUEResidSyst->Get("histSystErrFitNoUE_Cen0_10" + sShSh + sPtTrig);
      cout << "Sh: " << iSh << " " << histSystErrFitNoUE0_10[iSh][iptTr] << endl;
      histSystErrFitNoUE10_30[iSh][iptTr] = (TH1F *)fUEResidSyst->Get("histSystErrFitNoUE_Cen10_30" + sShSh + sPtTrig);
      cout << "Sh: " << iSh << " " << histSystErrFitNoUE10_30[iSh][iptTr] << endl;
      histSystErrFitNoUE0_30[iSh][iptTr] = (TH1F *)histSystErrFitNoUE0_10[iSh][iptTr]->Clone("histSystErrFitNoUE_Cen0_30" + sShSh + sPtTrig);
      cout << "Sh: " << iSh << " " << histSystErrFitNoUE0_30[iSh][iptTr] << endl;
      histSystErrFitNoUE0_30[iSh][iptTr]->Scale(numbTrig[0][iSh][iptTr]);
      histSystErrFitNoUE10_30[iSh][iptTr]->Scale(numbTrig[1][iSh][iptTr]);
      histSystErrFitNoUE0_30[iSh][iptTr]->Add(histSystErrFitNoUE10_30[iSh][iptTr]);
      histSystErrFitNoUE0_30[iSh][iptTr]->Scale(1 / numbTrig0_30[iSh][iptTr]);
      // new TCanvas();
      // histSystErrFitNoUE0_30[iSh][iptTr]->Draw("");
    }
  }
  // cout << "pippo2" << endl;

  // TString sPtTrig = Form("_PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
  // histSystErrFitNoUE0_30XPt[0][iptTr] = (TH1F *)histSystErrFitNoUE0_30[0][iptTr]->Clone("histSystErrFitNoUEXPt_ShSh0.10_0.30"+ sPtTrig);
  histSystErrFitNoUENarr0_30PtRange = SumPtBinXzt0_30(numbTrig0_30[0], ptTrig, index1, index2, histSystErrFitNoUE0_30[0], histSystErrFitNoUENarr0_30PtRange, histPur[0], funcPur[0], histPur[1], funcPur[1], numbTrig[0][0], numbTrig[1][0], 1, true);

  fitExpo = new TF1(Form("fitExpo_Cen0_30"), "expo", 0.2, 0.8);
  fitExpo->SetLineColor(kGreen + 3);
  fitExpo->SetLineStyle(10);
  fitExpo->SetLineWidth(3);
  fitConst = new TF1(Form("fitConst_Cen0_30"), "pol0", 0.05, 0.2);
  histSystErrFitNoUENarr0_30PtRange->Fit("fitExpo_Cen0_30", "R");
  histSystErrFitNoUENarr0_30PtRange->Fit("fitConst_Cen0_30", "R");
  TCanvas *cSystNarr0_30 = new TCanvas("cSystNarr0_30", "cSystNarr0_30", 800, 600);
  PlotStyle(histSystErrFitNoUENarr0_30PtRange, 20, 1.5, kAzure + 7, 0, "#it{z}_{T}", "Uncertainty %", false);
  // gStyle->SetOptFit(00000);
  histSystErrFitNoUENarr0_30PtRange->SetTitle("");
  histSystErrFitNoUENarr0_30PtRange->Draw("same");
  fitExpo->Draw("same");
  fitConst->Draw("same");
  TLatex *latUncer = LatexStd(latUncer, 0.450, 0.84, 0, 30, ptMin, ptMax, true);
  TLegend *legUncert = LegStd(legUncert, 0.45, 0.55, 0.80, 0.68);
  legUncert->AddEntry(fitExpo, "expo", "lp");
  legUncert->AddEntry(fitConst, "const", "lp");
  legUncert->Draw("same");
  cSystNarr0_30->Print(dirPlot + Form("/SystResidualUE/SystErrNarrCen0_30.pdf"));
  histSystErrFitNoUENarr0_30PtRangeFitTrend = new TH1F(Form("histSystErrFitNoUENarr0_30PtRangeFitTrend%s", sPtAll.Data()), Form("histSystErrFitNoUENarr0_30PtRangeFitTrend%s", sPtAll.Data()), nZtBin, assocZt);
  for (int ibin = 0; ibin < nZtBin; ibin++)
  {
    if (ibin == 0 || ibin == 1)
    {
      histSystErrFitNoUENarr0_30PtRangeFitTrend->SetBinContent(ibin + 1, fitConst->Eval(histSystErrFitNoUENarr0_30PtRangeFitTrend->GetBinCenter(ibin + 1)));
    }
    else
    {
      histSystErrFitNoUENarr0_30PtRangeFitTrend->SetBinContent(ibin + 1, fitExpo->Eval(histSystErrFitNoUENarr0_30PtRangeFitTrend->GetBinCenter(ibin + 1)));
    }
  }
  TCanvas *cSystNarr0_30FitTrend = new TCanvas("cSystNarr0_30FitTrend", "cSystNarr0_30FitTrend", 800, 600);
  PlotStyle(histSystErrFitNoUENarr0_30PtRangeFitTrend, 20, 1.5, kAzure + 7, 0, "#it{z}_{T}", "Uncertainty %", false);
  histSystErrFitNoUENarr0_30PtRangeFitTrend->SetTitle("");
  histSystErrFitNoUENarr0_30PtRangeFitTrend->Draw("pl");
  TLatex *latUncerFit = LatexStd(latUncerFit, 0.450, 0.84, 0, 30, ptMin, ptMax, true);

  // fitExpo->Draw("same");
  // fitConst->Draw("same");
  cSystNarr0_30FitTrend->Print(dirPlot + Form("/SystResidualUE/SystErrNarrCen0_30FitTrend.pdf"));
  fUEResidSyst->cd();
  histSystErrFitNoUENarr0_30PtRange->Write();
  histSystErrFitNoUENarr0_30PtRangeFitTrend->Write();
}


TLatex *LatexStdIcp(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax, bool sPttrig)
{
  lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetTextSize(0.04);
  lat->SetNDC();
  lat->DrawLatex(xpos, ypos, Form("#it{Work in progress}"));
  lat->DrawLatex(xpos, ypos - 0.06, Form("%d-%d%% / 50-90%% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV", cenMin, cenMax));
  if (sPttrig)
    lat->DrawLatex(xpos, ypos - 2 * 0.06, Form("%2.0f < #it{p}_{#it{T}}^{trig} < %2.0f GeV/#it{c}", ptMin, ptMax));
  return lat;
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
      cout << "Pippoooooo" << endl;
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
    for (int ibin = 0; ibin < hzTbin[iptTrig]->GetNbinsX(); ibin++)
    {
      cout << "Attenzione_______ibin: " << ibin + 1 << ")" << hzTbin[iptTrig]->GetBinCenter(ibin + 1) << "; content: " << hzTbin[iptTrig]->GetBinContent(ibin + 1) << endl;
    }
    cout << nPtRec << endl;
    cout << "hzTbinAll" << hzTbinAll << endl;
    cout << "hztbin: " << hzTbin[iptTrig] << endl;
    cout << iptTrig << endl;
    hzTbinAll->Add(hzTbin[iptTrig]);
    // cout << nPtRec << endl;
    intPtAll = intPtAll + nPtRec * pur;
    // cout << nPtRec << endl;
  }
  cout << "intPtAll" << intPtAll << endl;
  hzTbinAll->Scale(1 / intPtAll);
  cout << "pippo11111" << endl;
  return hzTbinAll;
}

TH1F *SumPtBinXzt0_30(Double_t nTrig0_30[npt], Float_t PtTrigger[npt], int index1, int index2, TH1F *hzTbin[npt], TH1F *hzTbinAll, TH1F *hPur010, TF1 *fPur010, TH1F *hPur1030, TF1 *fPur1030, Double_t nTrig0_10[npt], Double_t nTrig10_30[npt], double systPur, Bool_t bData)
{
  Float_t intPtAll = 0;
  cout << "(1/NintPur_tot) * Sum(pur * NInt * f(Zt))" << endl;
  for (int iptTrig = 0; iptTrig < index2 - index1; iptTrig++)
  {
    cout << PtTrigger[iptTrig + index1] << endl;
    // Float_t nPtRec = hTrigSame->Integral(hTrigSame->FindBin(PtTrigger[iptTrig + index1]), hTrigSame->FindBin(PtTrigger[iptTrig + index1 + 1] - 0.0001));
    Float_t nPtRec = nTrig0_30[iptTrig];
    Float_t pur = 1;
    Float_t pur010 = 1;
    Float_t pur1030 = 1;
    if (bData)
    {
      double syst = hPur1030->GetBinError(hPur1030->FindBin(PtTrigger[iptTrig + index1 + 1] - 0.0001)) / (hPur1030->GetBinContent(hPur1030->FindBin(PtTrigger[iptTrig + index1 + 1] - 0.0001)));
      pur010 = fPur010->Eval(hPur010->GetBinCenter(hPur010->FindBin(PtTrigger[iptTrig + index1 + 1] - 0.0001)));
      pur1030 = fPur1030->Eval(hPur1030->GetBinCenter(hPur1030->FindBin(PtTrigger[iptTrig + index1 + 1] - 0.0001)));
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
    pur = (pur010 * nTrig0_10[iptTrig] + pur1030 * nTrig10_30[iptTrig]) / nPtRec;
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