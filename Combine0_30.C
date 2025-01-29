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
#include <vector>
#include <algorithm>
#include <iterator>

using std::cout;
using std::endl;
int const nCen = 2;
int cenBins[] = {0, 10, 30};

int nIso = 2;
int nShSh = 2;

int nZtBin = 6;
double assocZt[] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 1.00};

int npt = 12;
float ptTrig[] = {10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 50, 60, 80};

double phiMin = TMath::Pi() * 3 / 5;
double phiMax = TMath::Pi();
TFile *fileDataShpp = 0;
TFile *fileDataShStd = 0;
TFile *fileData = 0;
TFile *fileDataShBkg = 0;
// TFile *fileDataMix = 0;
// TFile *fileDataMix1 = 0;
TFile *fileMC = 0;
TH1F *histPur[nCen];
TH1F *histPurStat[nCen];
TF1 *funcPur[nCen];
void PlotStyle(TH1F *hPlot, int kMarker, double kMarkerSize, int kColor, TString titleX, TString titleY);
TLatex *LatexStd(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax, bool sPttrig);
TLatex *LatexDPhi(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax);
TLegend *LegStd(TLegend *leg, double xpos1, double ypos1, double xpos2, double ypos2);
void ZtFunction(TH1F *hDeltaPhi, TH1F *hZT, int bin);
static void ScaleBinBySize(TH1F *h);
TH1F *SumPtBinXzt(Double_t nTrig0_30[npt], Float_t PtTrigger[npt], int index1, int index2, TH1F *hzTbin[npt], TH1F *hzTbinAll, TH1F *hPur010, TF1 *fPur010, TH1F *hPur1030, TF1 *fPur1030, Double_t nTrig0_10[npt], Double_t nTrig10_30[npt], double systPur, Bool_t bData);
void fZYAM(TH1F *hSame, double rangeMin = 3 * (TMath::Pi()) / 10, double rangeMax = TMath::Pi() / 2);
double fZYAM_Mix(TH1F *hSame, TH1F *hMix);
////////////////////////////////////////////////////////////////////////////////////////////////
/////// This macro compute the combination between 0-10% and 10-30% centrality intervals ///////
////////////////////////////////////////////////////////////////////////////////////////////////

void Combine0_30(float ptMin = 18, float ptMax = 40, int iCen = 0, bool bMirror = true, TString shshBkg = "0.40-1.00", TString dirFiles = "~/work/histogram/FromScratch/checkCode", double systPur = 1, bool bZYAM = false, bool bPlot = true, TString dirPlot = "~/work/histogram/FromScratch/FigcheckCode")
{
  TString sHistName = "AnaPhotonHadronCorr_";
  TString shshString[2] = {"0.10-0.30", shshBkg};
  TString sPtAll = Form("_Pt%2.0f_%2.0f", ptMin, ptMax);

  // index pt start and stop
  int nsize = sizeof(ptTrig) / sizeof(ptTrig[0]);
  auto itr1 = find(ptTrig, ptTrig + nsize, ptMin);
  auto itr2 = find(ptTrig, ptTrig + nsize, ptMax);
  int index1 = distance(ptTrig, itr1);
  int index2 = distance(ptTrig, itr2);
  int nPtTrig = index2 - index1;
  cout << index1 << "___" << index2 << ", " << nPtTrig << endl;

  // histogram for the analysis
  //
  TH1F *hTriggerSam[nCen][nShSh]; // PtTrig distrib SAME
  TH1F *hTriggerMix[nCen][nShSh]; // PtTrig distrib MIX

  TH1F *hTriggerSam0_30[nShSh]; // PtTrig distrib SAME

  // SameEvent and MixedEvent Mirror
  TH1F *hdPhiSamMirror[nCen][nShSh][nZtBin][nPtTrig];
  TH1F *hdPhiMixMirror[nCen][nShSh][nZtBin][nPtTrig];

  TH1F *hdPhiSame0_30[nShSh][nZtBin][nPtTrig];
  TH1F *hdPhiMix0_30[nShSh][nZtBin][nPtTrig];
  // NO UE (Same-Mix)
  TH1F *hdPhiSamNoUE0_30[nShSh][nZtBin][nPtTrig]; // Same - Mix

  // Purity correction
  TH1F *hdPhiSamPur[nCen][nZtBin][nPtTrig]; // IsoClusterGamma wo UE, subtract (1-p)IsoClPi0 and divided by Purity
  TH1F *hdPhiSamPi0Pur[nCen][nZtBin][nPtTrig];
  TH1F *hdPhiSamPi0Pur0_30[nZtBin][nPtTrig];
  TH1F *hdPhiMixPur[nCen][nZtBin][nPtTrig];
  TH1F *hdPhiMixPi0Pur[nCen][nZtBin][nPtTrig];

  double numbTrig[nCen][nShSh][nPtTrig];
  double numbTrig0_30[nShSh][nPtTrig];
  double numbTrigMix[nCen][nShSh][nPtTrig];
  double numbTrigMix0_30[nShSh][nPtTrig];
  // DeltaPhi for IsoGamma
  TH1F *hdPhiIsoGamma[nCen][nZtBin][nPtTrig]; // Final IsoGamma azimuthal distribution
  TH1F *hdPhiIsoGamma0_30[nZtBin][nPtTrig];

  // Zt Distribution
  TH1F *hZtIsoGammaPtBin[nPtTrig]; // Zt distribution X every PtBin (only purity correction)
  TH1F *hZtIsoGamma;               // Zt distribution full range (only purity correction)

  // Purity
  TFile *fileData[nCen];
  TFile *fPurity = new TFile("~/work/histogram/FromScratch/Purity_IsoSig1.5_M02Sig0.10-0.30_IsoBkg_4.0_25.0_M02Bkg0.40_2.00_LHC15o_18qr_L1MB.root");
  TFile *fOutPut = new TFile(Form("%s/fPlot%s_Cen0_30%s.root", dirFiles.Data(), shshBkg.Data(), sPtAll.Data()), "RECREATE");
  cout << fOutPut->GetName() << endl;
  // Data Results
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    cout << iCen << " " << cenBins[iCen] << "-" << cenBins[iCen + 1] << endl;
    histPur[iCen] = (TH1F *)fPurity->Get(Form("Purity_Cen%d_R0.2_Sys", iCen));
    histPurStat[iCen] = (TH1F *)fPurity->Get(Form("Purity_Cen%d_R0.2", iCen));

    funcPur[iCen] = histPur[iCen]->GetFunction("purityFitCombinedSigmoid");

    bool bPlot = true;

    gSystem->Exec(Form("mkdir %s", dirFiles.Data()));

    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);

    fileData[iCen] = new TFile(Form("%s/fPlot%s%s%s.root", dirFiles.Data(), shshBkg.Data(), sCent.Data(), sPtAll.Data()));
    // fileMix[iCen] = new TFile(Form("~/work/histogram/DataSh100_AssocPt500/EMCAL_MB_%d_%d.root", cenBins[iCen], cenBins[iCen + 1]));

    cout << "Getter Pt Trig distrib centrality: " << iCen << endl;
    for (int iSh = 0; iSh < nShSh; iSh++)
    {
      TString sShSh = Form("_ShSh%s", shshString[iSh].Data());
      hTriggerSam[iCen][iSh] = (TH1F *)fileData[iCen]->Get("AnaPhotonHadronCorr_Iso1" + sShSh + sCent + "_hPtTrigger");
      hTriggerMix[iCen][iSh] = (TH1F *)fileData[iCen]->Get("AnaPhotonHadronCorr_Iso1" + sShSh + sCent + "_hPtTriggerMixed");

      cout << fileData[iCen]->GetName() << endl;
      for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
      {
        numbTrig[iCen][iSh][iptTr] = hTriggerSam[iCen][iSh]->Integral(hTriggerSam[iCen][iSh]->FindBin(ptTrig[index1 + iptTr]), hTriggerSam[iCen][iSh]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001));
        numbTrigMix[iCen][iSh][iptTr] = hTriggerMix[iCen][iSh]->Integral(hTriggerMix[iCen][iSh]->FindBin(ptTrig[index1 + iptTr]), hTriggerMix[iCen][iSh]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001));
      }
    }
  }

  for (int iSh = 0; iSh < nShSh; iSh++)
  {
    for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
    {
      numbTrig0_30[iSh][iptTr] = numbTrig[0][iSh][iptTr] + numbTrig[1][iSh][iptTr];
      numbTrigMix0_30[iSh][iptTr] = numbTrigMix[0][iSh][iptTr] + numbTrigMix[1][iSh][iptTr];

      cout << "Pt: " << ptTrig[index1 + iptTr] << "Same 0-10\%: " << numbTrig[0][iSh][iptTr] << endl;
      cout << "Pt: " << ptTrig[index1 + iptTr] << "Same 10-30\%: " << numbTrig[1][iSh][iptTr] << endl;
      cout << "Pt: " << ptTrig[index1 + iptTr] << "Same 0-30\%: " << numbTrig0_30[iSh][iptTr] << endl;
      cout << "Pt: " << ptTrig[index1 + iptTr] << "Mix 0-10\%: " << numbTrigMix[0][iSh][iptTr] << endl;
      cout << "Pt: " << ptTrig[index1 + iptTr] << "Mix 10-30\%: " << numbTrigMix[1][iSh][iptTr] << endl;
      cout << "Pt: " << ptTrig[index1 + iptTr] << "Mix 0-30\%: " << numbTrigMix0_30[iSh][iptTr] << endl;
    }
  }

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    cout << "Getter Same & Mixed Azimuthal distrib (both ShSh) centrality: " << iCen << endl;
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    for (int iSh = 0; iSh < nShSh; iSh++)
    {
      TString sShSh = Form("_ShSh%s", shshString[iSh].Data());
      for (Int_t izt = 0; izt < nZtBin; izt++)
      {
        TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);
        for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
        {
          TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
          hdPhiSamMirror[iCen][iSh][izt][iptTr] = (TH1F *)fileData[iCen]->Get("hdPhiSameMirrorIso1" + sShSh + sZtBin + sPtTrig);
          hdPhiMixMirror[iCen][iSh][izt][iptTr] = (TH1F *)fileData[iCen]->Get("hdPhiMixMirrorIso1" + sShSh + sZtBin + sPtTrig);
        }
      }
    }
  }
  for (int iSh = 0; iSh < nShSh; iSh++)
  {
    TString sCent030 = "_Cen0_30";
    TString sShSh = Form("_ShSh%s", shshString[iSh].Data());
    hTriggerSam0_30[iSh] = (TH1F *)hTriggerSam[0][iSh]->Clone("AnaPhotonHadronCorr_Iso1" + sShSh + sCent030 + "_hPtTrigger");
    hTriggerSam0_30[iSh]->Add(hTriggerSam[1][iSh]);
    fOutPut->cd();
    hTriggerSam0_30[iSh]->Write();
  }

  for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
  {
    TString sCent030 = "_Cen0_30";
    TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
    for (int iSh = 0; iSh < nShSh; iSh++)
    {
      TString sShSh = Form("_ShSh%s", shshString[iSh].Data());
      for (Int_t izt = 0; izt < nZtBin; izt++)
      {
        TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);
        hdPhiSame0_30[iSh][izt][iptTr] = (TH1F *)hdPhiSamMirror[0][iSh][izt][iptTr]->Clone("hdPhiSameMirrorIso1" + sShSh + sZtBin + sPtTrig);
        hdPhiSame0_30[iSh][izt][iptTr]->Scale(numbTrig[0][iSh][iptTr]);
        hdPhiSamMirror[1][iSh][izt][iptTr]->Scale(numbTrig[1][iSh][iptTr]);
        hdPhiSame0_30[iSh][izt][iptTr]->Add(hdPhiSamMirror[1][iSh][izt][iptTr]);
        hdPhiSame0_30[iSh][izt][iptTr]->Scale(1 / numbTrig0_30[iSh][iptTr]);

        hdPhiMix0_30[iSh][izt][iptTr] = (TH1F *)hdPhiMixMirror[0][iSh][izt][iptTr]->Clone("hdPhiMixMirrorIso1" + sShSh + sZtBin + sPtTrig);
        hdPhiMix0_30[iSh][izt][iptTr]->Scale(numbTrigMix[0][iSh][iptTr]);
        hdPhiMixMirror[1][iSh][izt][iptTr]->Scale(numbTrigMix[1][iSh][iptTr]);
        hdPhiMix0_30[iSh][izt][iptTr]->Add(hdPhiMixMirror[1][iSh][izt][iptTr]);
        hdPhiMix0_30[iSh][izt][iptTr]->Scale(1 / numbTrigMix0_30[iSh][iptTr]);

        hdPhiSamNoUE0_30[iSh][izt][iptTr] = (TH1F *)hdPhiSame0_30[iSh][izt][iptTr]->Clone("hdPhiSameNoUEIso1" + sShSh + sZtBin + sPtTrig);
        hdPhiSamNoUE0_30[iSh][izt][iptTr]->Add(hdPhiMix0_30[iSh][izt][iptTr], -1);

        fOutPut->cd();
        hdPhiSame0_30[iSh][izt][iptTr]->Write();
        hdPhiMix0_30[iSh][izt][iptTr]->Write();
        hdPhiSamNoUE0_30[iSh][izt][iptTr]->Write();
      }
    }
  }

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    cout << "Getter IsoGamma Azimuthal distrib centrality: " << iCen << endl;
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    for (Int_t izt = 0; izt < nZtBin; izt++)
    {
      TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);
      for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
      {
        TString sPtTrig = Form("_PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
        hdPhiSamPi0Pur[iCen][izt][iptTr] = (TH1F *)fileData[iCen]->Get("hdPhiSamePi0PurIso1" + sZtBin + sPtTrig);
        hdPhiIsoGamma[iCen][izt][iptTr] = (TH1F *)fileData[iCen]->Get("hdPhiIso1Photon" + sZtBin + sPtTrig);
      }
    }
  }
  for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
  {
    TString sCent030 = "_Cen0_30";
    TString sPtTrig = Form("_PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
    hZtIsoGammaPtBin[iptTr] = new TH1F(Form("hZtIsoGammaPtBin%s", sPtTrig.Data()), Form("hZtIsoGammaPtBin%s", sPtTrig.Data()), nZtBin, assocZt);
    for (Int_t izt = 0; izt < nZtBin; izt++)
    {
      TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);
      cout << "pippo" << endl;
      hdPhiSamPi0Pur0_30[izt][iptTr] = (TH1F *)hdPhiSamPi0Pur[0][izt][iptTr]->Clone("hdPhiSamePi0Pur" + sZtBin + sPtTrig);
      hdPhiSamPi0Pur0_30[izt][iptTr]->Scale(numbTrig[0][1][iptTr]);
      hdPhiSamPi0Pur[1][izt][iptTr]->Scale(numbTrig[1][1][iptTr]);
      hdPhiSamPi0Pur0_30[izt][iptTr]->Add(hdPhiSamPi0Pur[1][izt][iptTr]);
      hdPhiSamPi0Pur0_30[izt][iptTr]->Scale(1 / numbTrig0_30[1][iptTr]);
      cout << "pippo" << endl;
      cout << hdPhiIsoGamma[0][izt][iptTr] << endl;
      hdPhiIsoGamma0_30[izt][iptTr] = (TH1F *)hdPhiIsoGamma[0][izt][iptTr]->Clone("hdPhiIso1Photon" + sZtBin + sPtTrig);
      cout << "pippo" << endl;
      hdPhiIsoGamma0_30[izt][iptTr]->Scale(numbTrig[0][0][iptTr]);
      hdPhiIsoGamma[1][izt][iptTr]->Scale(numbTrig[1][0][iptTr]);
      hdPhiIsoGamma0_30[izt][iptTr]->Add(hdPhiIsoGamma[1][izt][iptTr]);
      hdPhiIsoGamma0_30[izt][iptTr]->Scale(1 / numbTrig0_30[0][iptTr]);

      ZtFunction(hdPhiIsoGamma0_30[izt][iptTr], hZtIsoGammaPtBin[iptTr], izt);
      // new TCanvas();
      // hdPhiIsoGamma0_30[izt][iptTr]->Draw();
      fOutPut->cd();
      hdPhiIsoGamma0_30[izt][iptTr]->Write();
      hdPhiSamPi0Pur0_30[izt][iptTr]->Write();
    }
    ScaleBinBySize(hZtIsoGammaPtBin[iptTr]);
    fOutPut->cd();
    hZtIsoGammaPtBin[iptTr]->Write();
  }
  cout << "pippo" << endl;
  hZtIsoGamma = new TH1F(Form("hZtIso1Photon_Cen0_30%s", sPtAll.Data()), Form("hZtIso1Photon_Cen0_30%s", sPtAll.Data()), nZtBin, assocZt);
  hZtIsoGamma = SumPtBinXzt(numbTrig0_30[0], ptTrig, index1, index2, hZtIsoGammaPtBin, hZtIsoGamma, histPur[0], funcPur[0], histPur[1], funcPur[1], numbTrig[0][0], numbTrig[1][0], systPur, true);
  // hZtIsoGamma = SumPtBinXzt(numbTrig0_30, ptTrig, index1, index2, hZtIsoGammaPtBin, hZtIsoGamma, histPur[1], funcPur[1], systPur, true);
  new TCanvas();
  gPad->SetLogy();
  hZtIsoGamma->Draw();
  fOutPut->cd();
  hZtIsoGamma->Write();
  
  // MC Results

  TH1F *hTriggerSamMCGen[nCen];
  TH1F *hTriggerSamMCRec[nCen];
  TH1F *hTriggerSamMCRec0_30;
  TH1F *hTriggerMixMCRec[nCen];
  TH1F *hTriggerMixMCRec0_30;

  double numbTrigMCGen[nCen][nPtTrig];
  double numbTrigMCRec[nCen][nPtTrig];
  double numbTrigMixMCRec[nCen][nPtTrig];
  double numbTrigMCGen0_30[nPtTrig];
  double numbTrigMCRec0_30[nPtTrig];
  double numbTrigMixMCRec0_30[nPtTrig];

  TH1F *hdPhiMCGenMirror[nCen][nZtBin][nPtTrig];
  TH1F *hdPhiSamMCRec[nCen][nZtBin][nPtTrig];
  TH1F *hdPhiMixMCRec[nCen][nZtBin][nPtTrig];
  TH1F *hdPhiSamMCRecNoUE[nCen][nZtBin][nPtTrig]; // Same - Mixed Reconstructed wo UE
  TH1F *hdPhiMCGenMirror0_30[nZtBin][nPtTrig];
  TH1F *hdPhiSamMCRecNoUE0_30[nZtBin][nPtTrig]; // Same - Mixed Reconstructed wo UE
  TH1F *hdPhiSamMCRec0_30[nZtBin][nPtTrig];     // Same - Mixed Reconstructed wo UE
  TH1F *hdPhiMixMCRec0_30[nZtBin][nPtTrig];     // Same - Mixed Reconstructed wo UE
  TH1F *hZtIsoGammaPtBinMCGen[nPtTrig];         // Zt distribution MC Generated X Ptbin
  TH1F *hZtIsoGammaPtBinMCRec[nPtTrig];         // Zt distribution MC Recostructed X Ptbin
  TH1F *hRatioEffCorrPtBin[nPtTrig];            // Efficiency correction Gen/Rec X Ptbin

  TH1F *hZtIsoGammaMCGen;
  TH1F *hZtIsoGammaMCRec;

  TFile *fileMC = TFile::Open(Form("~/work/histogram/MCPtAssoc500/MC_GJ_0_90.root"));
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    cout << " MC Getter Pt Trig distrib centrality: " << iCen << endl;
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    hTriggerSamMCGen[iCen] = (TH1F *)fileMC->Get("AnaPhotonHadronCorr_Iso1_ShSh0.10-0.30" + sCent + "_hMCPtTrigger_Photon");
    hTriggerSamMCRec[iCen] = (TH1F *)fileMC->Get("AnaPhotonHadronCorr_Iso1_ShSh0.10-0.30" + sCent + "_hPtTrigger_MCPhoton");
    hTriggerMixMCRec[iCen] = (TH1F *)fileMC->Get("AnaPhotonHadronCorr_Iso1_ShSh0.10-0.30" + sCent + "_hPtTriggerMixed");
    cout << fileMC->GetName() << endl;
    for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
    {
      numbTrigMCGen[iCen][iptTr] = hTriggerSamMCGen[iCen]->Integral(hTriggerSamMCGen[iCen]->FindBin(ptTrig[index1 + iptTr]), hTriggerSamMCGen[iCen]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001));
      numbTrigMCRec[iCen][iptTr] = hTriggerSamMCRec[iCen]->Integral(hTriggerSamMCRec[iCen]->FindBin(ptTrig[index1 + iptTr]), hTriggerSamMCRec[iCen]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001));
      numbTrigMixMCRec[iCen][iptTr] = hTriggerMixMCRec[iCen]->Integral(hTriggerMixMCRec[iCen]->FindBin(ptTrig[index1 + iptTr]), hTriggerMixMCRec[iCen]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001));
    }
  }

  TString sCent030 = "_Cen0_30";
  hTriggerSamMCRec0_30 = (TH1F *)hTriggerSamMCRec[0]->Clone("AnaPhotonHadronCorr_Iso1_ShSh0.10-0.30" + sCent030 + "_hPtTrigger_MCPhoton");
  hTriggerSamMCRec0_30->Add(hTriggerSamMCRec[1]);

  hTriggerMixMCRec0_30 = (TH1F *)hTriggerMixMCRec[0]->Clone("AnaPhotonHadronCorr_Iso1_ShSh0.10-0.30" + sCent030 + "_hPtTriggerMixed");
  hTriggerMixMCRec0_30->Add(hTriggerMixMCRec[1]);

  fOutPut->cd();
  hTriggerSamMCRec0_30->Write();
  hTriggerMixMCRec0_30->Write();

  for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
  {
    numbTrigMCGen0_30[iptTr] = numbTrigMCGen[0][iptTr] + numbTrigMCGen[1][iptTr];
    numbTrigMCRec0_30[iptTr] = numbTrigMCRec[0][iptTr] + numbTrigMCRec[1][iptTr];
    numbTrigMixMCRec0_30[iptTr] = numbTrigMixMCRec[0][iptTr] + numbTrigMixMCRec[1][iptTr];
  }

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    cout << "Getter MC Azimuthal distrib centrality: " << iCen << endl;
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    for (Int_t izt = 0; izt < nZtBin; izt++)
    {
      TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);
      for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
      {
        TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
        hdPhiMCGenMirror[iCen][izt][iptTr] = (TH1F *)fileData[iCen]->Get("hdPhiMCGenMirrorIso1_ShSh0.10-0.30" + sZtBin + sPtTrig + "Photon");
        hdPhiSamMCRecNoUE[iCen][izt][iptTr] = (TH1F *)fileData[iCen]->Get("hdPhiSameMCRecNoUEIso1_ShSh0.10-0.30" + sZtBin + sPtTrig + "Photon");
        hdPhiSamMCRec[iCen][izt][iptTr] = (TH1F *)fileData[iCen]->Get("hdPhiSameMCRecMirrorIso1_ShSh0.10-0.30" + sZtBin + sPtTrig + "Photon");
        hdPhiMixMCRec[iCen][izt][iptTr] = (TH1F *)fileData[iCen]->Get("hdPhiMixMCRecMirrorIso1_ShSh0.10-0.30" + sZtBin + sPtTrig + "Photon");
        // cout << fileData[iCen]->GetName() << endl;
        // cout << hdPhiMCGenMirror[0][izt][iptTr] << endl;
      }
    }
  }

  for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
  {

    TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
    for (Int_t izt = 0; izt < nZtBin; izt++)
    {
      TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);

      cout << hdPhiSamMCRec0_30[izt][iptTr] << endl;
      hdPhiSamMCRec0_30[izt][iptTr] = (TH1F *)hdPhiSamMCRec[0][izt][iptTr]->Clone("hdPhiSameMCRecMirrorIso1_ShSh0.10-0.30" + sZtBin + sPtTrig + "Photon");
      hdPhiSamMCRec0_30[izt][iptTr]->Scale(numbTrigMCRec[0][iptTr]);
      hdPhiSamMCRec[1][izt][iptTr]->Scale(numbTrigMCRec[1][iptTr]);
      hdPhiSamMCRec0_30[izt][iptTr]->Add(hdPhiSamMCRec[1][izt][iptTr]);
      hdPhiSamMCRec0_30[izt][iptTr]->Scale(1 / numbTrigMCRec0_30[iptTr]);

      cout << hdPhiMixMCRec0_30[izt][iptTr] << endl;
      hdPhiMixMCRec0_30[izt][iptTr] = (TH1F *)hdPhiMixMCRec[0][izt][iptTr]->Clone("hdPhiMixMCRecMirrorIso1_ShSh0.10-0.30" + sZtBin + sPtTrig + "Photon");
      hdPhiMixMCRec0_30[izt][iptTr]->Scale(numbTrigMixMCRec[0][iptTr]);
      hdPhiMixMCRec[1][izt][iptTr]->Scale(numbTrigMixMCRec[1][iptTr]);
      hdPhiMixMCRec0_30[izt][iptTr]->Add(hdPhiMixMCRec[1][izt][iptTr]);
      hdPhiMixMCRec0_30[izt][iptTr]->Scale(1 / numbTrigMixMCRec0_30[iptTr]);

      // new TCanvas();
      // hdPhiIsoGamma0_30[izt][iptTr]->Draw();
      fOutPut->cd();
      hdPhiSamMCRec0_30[izt][iptTr]->Write();
      hdPhiMixMCRec0_30[izt][iptTr]->Write();
    }
  }

  TH1F *hRatioEffCorr;
  TH1F *hZtEffCorr;
  for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
  {
    TString sCent030 = "_Cen0_30";
    TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
    hZtIsoGammaPtBinMCGen[iptTr] = new TH1F(Form("hZtPtBinMCGen_Iso1%sPhoton", sPtTrig.Data()), Form("hZtPtBinMCGen_Iso1%sPhoton", sPtTrig.Data()), nZtBin, assocZt);
    hZtIsoGammaPtBinMCRec[iptTr] = new TH1F(Form("hZtPtBinMCRec_Iso1%sPhoton", sPtTrig.Data()), Form("hZtPtBinMCRec_Iso1%sPhoton", sPtTrig.Data()), nZtBin, assocZt);
    for (Int_t izt = 0; izt < nZtBin; izt++)
    {
      TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);

      hdPhiMCGenMirror0_30[izt][iptTr] = (TH1F *)hdPhiMCGenMirror[0][izt][iptTr]->Clone("hdPhiMCGenMirrorIso1_ShSh0.10-0.30" + sZtBin + sPtTrig + "Photon");
      hdPhiMCGenMirror0_30[izt][iptTr]->Scale(numbTrigMCGen[0][iptTr]);
      hdPhiMCGenMirror[1][izt][iptTr]->Scale(numbTrigMCGen[1][iptTr]);
      hdPhiMCGenMirror0_30[izt][iptTr]->Add(hdPhiMCGenMirror[1][izt][iptTr]);
      hdPhiMCGenMirror0_30[izt][iptTr]->Scale(1 / numbTrigMCGen0_30[iptTr]);

      // cout << hdPhiSamMCRecNoUE[0][izt][iptTr] << endl;
      hdPhiSamMCRecNoUE0_30[izt][iptTr] = (TH1F *)hdPhiSamMCRecNoUE[0][izt][iptTr]->Clone("hdPhiSameMCRecNoUEIso1_ShSh0.10-0.30" + sZtBin + sPtTrig + "Photon");
      hdPhiSamMCRecNoUE0_30[izt][iptTr]->Scale(numbTrigMCRec[0][iptTr]);
      hdPhiSamMCRecNoUE[1][izt][iptTr]->Scale(numbTrigMCRec[1][iptTr]);
      hdPhiSamMCRecNoUE0_30[izt][iptTr]->Add(hdPhiSamMCRecNoUE[1][izt][iptTr]);
      hdPhiSamMCRecNoUE0_30[izt][iptTr]->Scale(1 / numbTrigMCRec0_30[iptTr]);

      ZtFunction(hdPhiMCGenMirror0_30[izt][iptTr], hZtIsoGammaPtBinMCGen[iptTr], izt);
      ZtFunction(hdPhiSamMCRecNoUE0_30[izt][iptTr], hZtIsoGammaPtBinMCRec[iptTr], izt);
      // new TCanvas();
      // hdPhiIsoGamma0_30[izt][iptTr]->Draw();
      fOutPut->cd();
      hdPhiMCGenMirror0_30[izt][iptTr]->Write();
      hdPhiSamMCRecNoUE0_30[izt][iptTr]->Write();
    }
    ScaleBinBySize(hZtIsoGammaPtBinMCGen[iptTr]);
    ScaleBinBySize(hZtIsoGammaPtBinMCRec[iptTr]);
    hRatioEffCorrPtBin[iptTr] = (TH1F *)hZtIsoGammaPtBinMCGen[iptTr]->Clone("hRatioEffCorrPtBin" + sPtTrig);
    hRatioEffCorrPtBin[iptTr]->Sumw2();
    hRatioEffCorrPtBin[iptTr]->Divide(hZtIsoGammaPtBinMCRec[iptTr]);
    fOutPut->cd();
    hZtIsoGammaPtBinMCGen[iptTr]->Write();
    hZtIsoGammaPtBinMCRec[iptTr]->Write();
    hRatioEffCorrPtBin[iptTr]->Write();
  }

  hZtIsoGammaMCGen = new TH1F(Form("hZtMCGenIso1Photon_Cen0_30%s", sPtAll.Data()), Form("hZtMCGenIso1Photon_Cen0_30%s", sPtAll.Data()), nZtBin, assocZt);
  hZtIsoGammaMCGen = SumPtBinXzt(numbTrigMCGen0_30, ptTrig, index1, index2, hZtIsoGammaPtBinMCGen, hZtIsoGammaMCGen, histPur[0], funcPur[0], histPur[1], funcPur[1], numbTrigMCGen[0], numbTrigMCGen[1], systPur, false);

  fOutPut->cd();
  hZtIsoGammaMCGen->Write();

  hZtIsoGammaMCRec = new TH1F(Form("hZtMCRecIso1Photon_Cen0_30%s", sPtAll.Data()), Form("hZtMCRecIso1Photon_Cen0_30%s", sPtAll.Data()), nZtBin, assocZt);
  hZtIsoGammaMCRec = SumPtBinXzt(numbTrigMCRec0_30, ptTrig, index1, index2, hZtIsoGammaPtBinMCRec, hZtIsoGammaMCRec, histPur[0], funcPur[0], histPur[1], funcPur[1], numbTrigMCRec[0], numbTrigMCRec[1], systPur, false);
  new TCanvas();
  gPad->SetLogy();
  hZtIsoGammaMCGen->Draw("same");
  hZtIsoGammaMCRec->Draw("same");

  hRatioEffCorr = (TH1F *)hZtIsoGammaMCGen->Clone("hRatioEffCorrIso1Photon");
  hRatioEffCorr->Divide(hZtIsoGammaMCRec);

  hZtEffCorr = (TH1F *)hZtIsoGamma->Clone(Form("hZtEffCorrIso1Photon_Cen0_30%s", sPtAll.Data()));
  hZtEffCorr->Multiply(hRatioEffCorr);

  new TCanvas();
  hZtEffCorr->Draw();
  fOutPut->cd();
  hZtIsoGammaMCRec->Write();
  hRatioEffCorr->Write();
  hZtEffCorr->Write();

  if (bPlot)
  {
    PlotStyle(hZtIsoGammaMCGen, 21, 1, kAzure + 7, "#font[12]{z}_{T}", "1/N^{trig}dN^{charg}/d#font[12]{z}_{T}");
    PlotStyle(hZtIsoGammaMCRec, 21, 1, kOrange + 7, "#font[12]{z}_{T}", "1/N^{trig}dN^{charg}/d#font[12]{z}_{T}");
    PlotStyle(hZtEffCorr, 20, 1, kCyan + 2, "#it{z}_{T}", "1 / #it{N}^{ trig} d^{3}#it{N} / d#Delta#it{#eta} d |#Delta#it{#varphi}| d#it{z}_{T}");

    gSystem->Exec(Form("mkdir %s", dirPlot.Data()));
    cout << "Directory: " << dirPlot.Data() << endl;
    TString sdirPlotXCent = Form("%s/Cen0_30", dirPlot.Data());
    gSystem->Exec(Form("mkdir %s", sdirPlotXCent.Data()));
    TString sCent030 = "_Cen0_30";
    TCanvas *cZtIsoGammaEffCorrGlob_MC;
    TCanvas *cZtIsoGammaEffCorrGlob;
    TLatex *latexGlob;

    TString sCent = Form("_Cen0_30");
    cZtIsoGammaEffCorrGlob_MC = new TCanvas("cZtIsoGammaEffCorrGlob_MC" + sCent030, "cZtIsoGammaEffCorrGlob_MC" + sCent030, 800, 600);
    cZtIsoGammaEffCorrGlob_MC->cd();
    gPad->SetLogy();
    hZtIsoGammaMCGen->SetTitle("");
    hZtIsoGammaMCGen->Draw("same");
    hZtIsoGammaMCRec->Draw("same");
    hZtEffCorr->Draw("same");
    latexGlob = LatexStd(latexGlob, 0.440, 0.84, 0, 30, ptTrig[index1], ptTrig[index2], true);
    cZtIsoGammaEffCorrGlob_MC->Print(dirPlot + Form("/Cen0_30") + "/ZtDistribution_Data_MC" + sCent030 + sPtAll + ".pdf");

    cZtIsoGammaEffCorrGlob = new TCanvas("cZtIsoGammaEffCorrGlob" + sCent030, "cZtIsoGammaEffCorrGlob" + sCent030, 800, 600);
    cZtIsoGammaEffCorrGlob->cd();
    gPad->SetLogy();
    hZtEffCorr->GetYaxis()->SetRangeUser(1e-3, 50);
    hZtEffCorr->SetTitle("");
    hZtEffCorr->Draw("same");
    latexGlob = LatexStd(latexGlob, 0.440, 0.84, 0, 30, ptTrig[index1], ptTrig[index2], true);
    cZtIsoGammaEffCorrGlob->Print(dirPlot + Form("/Cen0_30") + "/ZtDistribution_Data" + sCent030 + sPtAll + ".pdf");
  }
}

void ZtFunction(TH1F *hDeltaPhi, TH1F *hZT, int bin)
{

  cout << "Angle: " << phiMin << " " << phiMax << endl;
  double binPhiMin = hDeltaPhi->FindBin(phiMin);
  double binPhiMax = hDeltaPhi->FindBin(phiMax - 0.0001);
  // cout<<hDeltaPhi->GetNbinsX()<<"___"<<hDeltaPhi->GetBinCenter(hDeltaPhi->GetNbinsX())<<endl;
  // cout<<binPhiMin<<"__"<<binPhiMax<<endl;
  // cout<<"width: "<<hDeltaPhi->GetBinWidth(binPhiMin)<<endl;
  // cout << "Valori"<<hDeltaPhi->GetBinCenter(binPhiMin) << "   " << hDeltaPhi->GetBinCenter(binPhiMax) << endl;
  double intPhiErr;
  double intPhi = hDeltaPhi->IntegralAndError(binPhiMin, binPhiMax, intPhiErr);
  hZT->SetBinContent(bin + 1, intPhi);
  hZT->SetBinError(bin + 1, intPhiErr);
  cout << "_____Bin: " << assocZt[bin] << "-" << assocZt[bin + 1] << "Integral: " << intPhi << endl;
  if (intPhi < 0)
    cout << "_____Bin: " << assocZt[bin] << "-" << assocZt[bin + 1] << "Integral: " << intPhi << endl;
}

TH1F *SumPtBinXzt(Double_t nTrig0_30[npt], Float_t PtTrigger[npt], int index1, int index2, TH1F *hzTbin[npt], TH1F *hzTbinAll, TH1F *hPur010, TF1 *fPur010, TH1F *hPur1030, TF1 *fPur1030, Double_t nTrig0_10[npt], Double_t nTrig10_30[npt], double systPur, Bool_t bData)
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

static void ScaleBinBySize(TH1F *h)
{
  for (Int_t ibin = 1; ibin <= h->GetNbinsX(); ibin++)
  {
    Double_t width = h->GetBinWidth(ibin);
    Double_t content = h->GetBinContent(ibin);
    Double_t error = h->GetBinError(ibin);
    // printf("bin %d, width %f, content %e\n",ibin,width,content);
    // cout<<h->GetNbinsX()<<endl;
    h->SetBinContent(ibin, content / width);
    h->SetBinError(ibin, error / width);
  }
}

void PlotStyle(TH1F *hPlot, int kMarker, double kMarkerSize, int kColor, TString titleX, TString titleY)
{
  gStyle->SetTitleX(0.52);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0000);

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
  lat->DrawLatex(xpos, ypos - 0.10, Form("%d-%d %% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV", cenMin, cenMax));
  return lat;
}
TLatex *LatexStd(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax, bool sPttrig)
{
  lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetTextSize(0.04);
  lat->SetNDC();
  lat->DrawLatex(xpos, ypos, Form("#it{This Thesis}"));
  lat->DrawLatex(xpos, ypos - 0.06, Form("#bf{%d#font[122]{-}%d %%} Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV,", cenMin, cenMax));
  if (sPttrig)
    lat->DrawLatex(xpos, ypos - 2 * 0.06, Form(" |#it{#eta}^{ trig}| < 0.67, %2.0f < #it{p}_{#it{T}}^{ trig} < %2.0f GeV/#it{c}", ptMin, ptMax));
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