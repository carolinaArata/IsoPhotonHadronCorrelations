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
#include "../Plotting.h"

using std::cout;
using std::endl;

int const nCen = 4;
int cenBins[] = {0, 10, 30, 50, 90};
// int cenBins[] = {30, 50, 90};

int nIso = 2;
int const nShSh = 2;
int nZtBinThin = 10;
double assocZtThinner[] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20};
// double assocZt[] = { 0.20, 0.30, 0.40};
int nZtBin = 6;
double assocZt[] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 1.00};
int npt = 12;
float ptTrig[] = {5, 8, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 50, 60, 80};
// double phiMin = TMath::DegToRad() * 120;
// double phiMax = TMath::DegToRad() * 240;
// double phiMin = TMath::Pi() * 3 / 5; default
// double phiMax = TMath::Pi(); default
Int_t kMarkStyleIso_NotIso[2] = {21, 93};
Int_t kColorMCRec[2] = {kAzure + 3, kPink - 4};

//{2, 4, kOrange + 7, kAzure + 3, kPink - 4, kViolet - 7, kBlue + 2, kTeal - 6};
Int_t kColorMC[2][3] = {{2, 4, kOrange + 7}, {kPink - 4, kViolet - 7, kBlue + 2}};
TFile *fileDataShpp = 0;
TFile *fileDataShStd = 0;
TFile *fileData = 0;
TFile *fileDataShBkg = 0;
// TFile *fileDataMix = 0;
TFile *fileMC;
TH1F *histPur[nCen];
TH1F *histPurStat[nCen];
TF1 *funcPur[nCen];
void Mirroring(TH1F *hMir, TH1F *hMirXtrue);
// void PlotStyle(TH1F *hPlot, int kMarker, int kMarkerSize, int kColor, TString titleX, TString titleY);
void ZtFunction(TH1F *hDeltaPhi, TH1F *hZT, int bin, double phiMin, double phiMax);
static void ScaleBinBySize(TH1F *h);
TH1F *SumPtBinXzt(TH1F *hTrigSame, Float_t PtTrigger[npt], int index1, int index2, TH1F *hzTbin[npt], TH1F *hzTbinAll, TH1F *hPur, TF1 *fPur, double systPur, Bool_t bData);
void fZYAM(TH1F *hSame, double rangeMin = 3 * (TMath::Pi()) / 10, double rangeMax = TMath::Pi() / 2);
double fZYAM_Mix(TH1F *hSame, TH1F *hMix);

void Exec(float ptMin = 18, float ptMax = 20, int iCen = 0, bool bMirror = true, TString shshBkg = "0.40-2.00", TString dirFiles = "ResultsStudyMC", double systPur = 1, bool bZYAM = false, double phiMin = TMath::Pi() * 3 / 5, double phiMax = TMath::Pi(), bool GJppOn = true)
{
  TString sHistName = "AnaPhotonHadronCorr_";
  TString shshString[2] = {"0.10-0.30", shshBkg};

  int const nGenType = 2;
  TString sTypeMCGenType[nGenType] = {"Pi0", "Pi0Decay"};
  // int const nGenType = 2;
  // TString sTypeMCGenType[nGenType] = {"Pi0", "Pi0Decay"};

  // index pt start and stop
  int nsize = sizeof(ptTrig) / sizeof(ptTrig[0]);
  auto itr1 = find(ptTrig, ptTrig + nsize, ptMin);
  auto itr2 = find(ptTrig, ptTrig + nsize, ptMax);
  int index1 = distance(ptTrig, itr1);
  int index2 = distance(ptTrig, itr2);
  int nPtTrig = index2 - index1;
  cout << index1 << "___" << index2 << ", " << nPtTrig << endl;

  cout << iCen << " " << cenBins[iCen] << "-" << cenBins[iCen + 1] << endl;
  // Purity

  TFile *fPurity = new TFile("RootFiles/Purity.root");
  histPur[iCen] = (TH1F *)fPurity->Get(Form("Purity_Cen%d_R0.2_Sys", iCen));
  histPurStat[iCen] = (TH1F *)fPurity->Get(Form("Purity_Cen%d_R0.2", iCen));

  funcPur[iCen] = histPur[iCen]->GetFunction("purityFitCombinedSigmoid");

  bool bPlot = true;  
  TString processline = Form(".! mkdir -pv %s", dirFiles.Data()) ;
  gROOT->ProcessLine(processline.Data());
  
  cout << iCen << " " << cenBins[iCen] << "-" << cenBins[iCen + 1] << endl;
  TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
  TString sPtAll = Form("_Pt%2.0f_%2.0f", ptMin, ptMax);

  TFile *fOutPut = new TFile(Form("%s/fPlot%s%s%s.root", dirFiles.Data(), shshBkg.Data(), sCent.Data(), sPtAll.Data()), "RECREATE");
  cout << fOutPut->GetName() << endl;

  if (bZYAM)
    cout << "ZYAM UE subtraction has been chosen" << endl;

  // MC_GJ and JJlow
  TH1F *hTriggerMCGen[nIso][nShSh][nGenType];
  TH1F *hTriggerSamMCRec[nIso][nShSh][nGenType];
  TH1F *hTriggerMixMCRec[nIso][nShSh];

  TH2F *h2DdPhidEtaMCGen[nIso][nShSh][nZtBin][nGenType][nPtTrig];
  TH2F *h2DdPhidEtaSamMCRec[nIso][nShSh][nGenType][nZtBin][nPtTrig];
  TH2F *h2DdPhidEtaMixMCRec[nIso][nShSh][nZtBin][nPtTrig];

  TH1F *hdPhiMCGen[nIso][nShSh][nGenType][nZtBin][nPtTrig];         // Generated level
  TH1F *hdPhiSamMCRec[nIso][nShSh][nGenType][nZtBin][nPtTrig];      // Same Reconstructed w/ UE
  TH1F *hdPhiMCGenMirrorUE[nIso][nShSh][nGenType][nZtBin][nPtTrig]; // Same Generated with UE
  TH1F *hdPhiMixMCRec[nIso][nShSh][nZtBin][nPtTrig];
  TH1F *hdPhiMCGenMirror[nIso][nShSh][nGenType][nZtBin][nPtTrig]; // Same Generated w/ UE
  TH1F *hdPhiMixMCRecMirror[nIso][nShSh][nZtBin][nPtTrig];
  TH1F *hdPhiSamMCRecMirror[nIso][nShSh][nGenType][nZtBin][nPtTrig];
  TH1F *hdPhiSamMCRecNoUE[nIso][nShSh][nGenType][nZtBin][nPtTrig];     // Same - Mixed Reconstructed wo UE
  TH1F *hdPhiSamMCRecNoUEZYAM[nIso][nShSh][nGenType][nZtBin][nPtTrig]; // ZYAM Reconstructed wo UE
  TH1F *hdPhiIsoGammaMCGen[nGenType][nZtBin][nPtTrig];

  TH1F *hdPhiIsoGammaMCRec[nGenType][nZtBin][nPtTrig];
  TH1F *hdPhiIsoGammaMCRecZYAM[nGenType][nZtBin][nPtTrig];
  TH1F *hZtIsoGammaPtBinMCGen[nIso][nShSh][nGenType][nPtTrig]; // Zt distribution MC Generated X Ptbin
  TH1F *hZtIsoGammaPtBinMCRec[nIso][nShSh][nGenType][nPtTrig]; // Zt distribution MC Recostructed X Ptbin
  TH1F *hZtIsoGammaPtBinMCRecZYAM[nGenType][nPtTrig];          // Zt distribution MC Recostructed X Ptbin
  TH1F *hRatioEffCorrPtBin[nGenType][nPtTrig];                 // Efficiency correction Gen/Rec X Ptbin
  TH1F *hRatioEffCorrPtBinZYAM[nGenType][nPtTrig];             // Efficiency correction Gen/Rec X Ptbin ZYAM

  TH1F *hZtIsoGammaMCGen[nIso][nShSh][nGenType];

  TH1F *hZtIsoGammaMCRec[nIso][nShSh][nGenType];
  TH1F *hZtIsoGammaMCRecZYAM;
  TH1F *hRatioEffCorr;
  TH1F *hRatioEffCorrZYAM;
  TH1F *hZtMixvsZYAM;

  //          h3Sam = (TH3F *)fileData->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt], assocZtThinner[izt + 1]));
  // h3Mix = (TH3F *)fileData->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt], assocZtThinner[izt + 1]));

  for (int iso = 0; iso < nIso; iso++)
  {
    TString sIso = Form("Iso%d", iso);

    for (int iSh = 0; iSh < nShSh; iSh++)
    {
      TString sShSh = Form("_ShSh%s", shshString[iSh].Data());
      for (int iGenType = 0; iGenType < nGenType; iGenType++)
      {
        hTriggerMCGen[iso][iSh][iGenType] = (TH1F *)fileMC->Get(sHistName + sIso + sShSh + sCent + Form("_hMCPtTrigger_%s", sTypeMCGenType[iGenType].Data()));
        cout << hTriggerMCGen[iso][iSh][iGenType] << endl;
        hTriggerSamMCRec[iso][iSh][iGenType] = (TH1F *)fileMC->Get(sHistName + sIso + sShSh + sCent + Form("_hPtTrigger_MC%s", sTypeMCGenType[iGenType].Data()));
        cout << hTriggerSamMCRec[iso][iSh][iGenType] << endl;
        fOutPut->cd();
        hTriggerSamMCRec[iso][iSh][iGenType]->Write();
        hTriggerMCGen[iso][iSh][iGenType]->Write();
      }
      hTriggerMixMCRec[iso][iSh] = (TH1F *)fileMC->Get(sHistName + sIso + sShSh + sCent + "_hPtTriggerMixed");

      for (Int_t izt = 0; izt < nZtBin; izt++)
      {
        TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);
        TH3F *h3MCGen[nShSh][nGenType];
        TH3F *h3MCGenNext[nShSh][nGenType];
        TH3F *h3SamMCRec[nShSh][nGenType];
        TH3F *h3SamMCRecNext[nShSh][nGenType];
        TH3F *h3MixMCRec[nShSh];
        TH3F *h3MixMCRecNext[nShSh];
        if (izt <= 3)
        {
          for (int iGenType = 0; iGenType < nGenType; iGenType++)
          {
            // generated
            h3MCGen[iSh][iGenType] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hMCDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_%s", assocZt[izt], assocZt[izt + 1], sTypeMCGenType[iGenType].Data()));
            // reconstructed
            h3SamMCRec[iSh][iGenType] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_RecoMC%s", assocZtThinner[izt], assocZtThinner[izt + 1], sTypeMCGenType[iGenType].Data()));
          }
          h3MixMCRec[iSh] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt], assocZtThinner[izt + 1]));
        }
        else if (izt == 4)
        {
          for (int iGenType = 0; iGenType < nGenType; iGenType++)
          {
            // generated
            h3MCGen[iSh][iGenType] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hMCDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_%s", assocZtThinner[izt], assocZtThinner[izt + 1], sTypeMCGenType[iGenType].Data()));
            h3MCGenNext[iSh][iGenType] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hMCDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_%s", assocZtThinner[izt + 1], assocZtThinner[izt + 2], sTypeMCGenType[iGenType].Data()));
            h3MCGen[iSh][iGenType]->Add(h3MCGenNext[iSh][iGenType]);
            // reconstructed
            h3SamMCRec[iSh][iGenType] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_RecoMC%s", assocZtThinner[izt], assocZtThinner[izt + 1], sTypeMCGenType[iGenType].Data()));
            h3SamMCRecNext[iSh][iGenType] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_RecoMC%s", assocZtThinner[izt + 1], assocZtThinner[izt + 2], sTypeMCGenType[iGenType].Data()));
            h3SamMCRec[iSh][iGenType]->Add(h3SamMCRecNext[iSh][iGenType]);
          }

          h3MixMCRec[iSh] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt], assocZtThinner[izt + 1]));
          h3MixMCRecNext[iSh] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 1], assocZtThinner[izt + 2]));
          h3MixMCRec[iSh]->Add(h3MixMCRecNext[iSh]);
        }
        else if (izt == 5)
        {
          for (int iGenType = 0; iGenType < nGenType; iGenType++)
          {
            // generated
            h3MCGen[iSh][iGenType] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hMCDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_%s", assocZtThinner[izt + 1], assocZtThinner[izt + 2], sTypeMCGenType[iGenType].Data()));
            h3MCGenNext[iSh][iGenType] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hMCDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_%s", assocZtThinner[izt + 2], assocZtThinner[izt + 3], sTypeMCGenType[iGenType].Data()));
            h3MCGen[iSh][iGenType]->Add(h3MCGenNext[iSh][iGenType]);

            // reconstructed
            h3SamMCRec[iSh][iGenType] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_RecoMC%s", assocZtThinner[izt + 1], assocZtThinner[izt + 2], sTypeMCGenType[iGenType].Data()));
            h3SamMCRecNext[iSh][iGenType] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_RecoMC%s", assocZtThinner[izt + 2], assocZtThinner[izt + 3], sTypeMCGenType[iGenType].Data()));
            h3SamMCRec[iSh][iGenType]->Add(h3SamMCRecNext[iSh][iGenType]);
          }

          h3MixMCRec[iSh] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 1], assocZtThinner[izt + 2]));
          h3MixMCRecNext[iSh] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 2], assocZtThinner[izt + 3]));
          h3MixMCRec[iSh]->Add(h3MixMCRecNext[iSh]);
        }
        else if (izt == 6)
        {
          for (int iGenType = 0; iGenType < nGenType; iGenType++)
          {
            h3MCGen[iSh][iGenType] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hMCDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_%s", assocZtThinner[izt + 2], assocZtThinner[izt + 3], sTypeMCGenType[iGenType].Data()));
            h3MCGenNext[iSh][iGenType] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hMCDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_%s", assocZtThinner[izt + 3], assocZtThinner[izt + 4], sTypeMCGenType[iGenType].Data()));
            h3MCGen[iSh][iGenType]->Add(h3MCGenNext[iSh][iGenType]);

            h3SamMCRec[iSh][iGenType] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_RecoMC%s", assocZtThinner[izt + 2], assocZtThinner[izt + 3], sTypeMCGenType[iGenType].Data()));
            h3SamMCRecNext[iSh][iGenType] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_RecoMC%s", assocZtThinner[izt + 3], assocZtThinner[izt + 4], sTypeMCGenType[iGenType].Data()));
            h3SamMCRec[iSh][iGenType]->Add(h3SamMCRecNext[iSh][iGenType]);
          }
          h3MixMCRec[iSh] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 2], assocZtThinner[izt + 3]));
          h3MixMCRecNext[iSh] = (TH3F *)fileMC->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 3], assocZtThinner[izt + 4]));
          h3MixMCRec[iSh]->Add(h3MixMCRecNext[iSh]);
        }
        int nbinsX = 0; //to be checked
        for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
        {
          cout << "pippo" << endl;
          TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
          for (int iGenType = 0; iGenType < nGenType; iGenType++)
          {
            // Generated
            // 3D Gen
            h3MCGen[iSh][iGenType]->SetAxisRange(ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1] - 0.0001, "X"); 
            // 2D Gen
            h2DdPhidEtaMCGen[iso][iSh][iGenType][izt][iptTr] = (TH2F *)h3MCGen[iSh][iGenType]->Project3D("zy");
            h2DdPhidEtaMCGen[iso][iSh][iGenType][izt][iptTr]->SetName(Form("h2DdPhidEtaMCGen%s", sTypeMCGenType[iGenType].Data()) + sIso + sShSh + sZtBin + sPtTrig);
            // 1D Gen
            hdPhiMCGen[iso][iSh][iGenType][izt][iptTr] = (TH1F *)h2DdPhidEtaMCGen[iso][iSh][iGenType][izt][iptTr]->ProjectionX(Form("hdPhiMCGen%s", sTypeMCGenType[iGenType].Data()) + sIso + sShSh + sZtBin + sPtTrig);
            // Normalisation Gen
            hdPhiMCGen[iso][iSh][iGenType][izt][iptTr]->Scale(1. / hTriggerMCGen[iso][iSh][iGenType]->Integral(hTriggerMCGen[iso][iSh][iGenType]->FindBin(ptTrig[index1 + iptTr]), hTriggerMCGen[iso][iSh][iGenType]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001)));
            // Reconstructed
            // 3D Rec
            h3SamMCRec[iSh][iGenType]->SetAxisRange(ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1] - 0.0001, "X");
            // 2D Rec
            h2DdPhidEtaSamMCRec[iso][iSh][iGenType][izt][iptTr] = (TH2F *)h3SamMCRec[iSh][iGenType]->Project3D("zy");
            h2DdPhidEtaSamMCRec[iso][iSh][iGenType][izt][iptTr]->SetName(Form("h2DdPhidEtaSamMCRec%s", sTypeMCGenType[iGenType].Data()) + sIso + sShSh + sZtBin + sPtTrig);
            // 1D Rec
            hdPhiSamMCRec[iso][iSh][iGenType][izt][iptTr] = (TH1F *)h2DdPhidEtaSamMCRec[iso][iSh][iGenType][izt][iptTr]->ProjectionX(Form("hdPhiSameMCRec%s", sTypeMCGenType[iGenType].Data()) + sIso + sShSh + sZtBin + sPtTrig);
            // Normalisation Normalisation
            hdPhiSamMCRec[iso][iSh][iGenType][izt][iptTr]->Scale(1. / hTriggerSamMCRec[iso][iSh][iGenType]->Integral(hTriggerSamMCRec[iso][iSh][iGenType]->FindBin(ptTrig[index1 + iptTr]), hTriggerSamMCRec[iso][iSh][iGenType]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001)));
            // Mirroring
            nbinsX = hdPhiSamMCRec[iso][iSh][iGenType][izt][iptTr]->GetNbinsX() / 2;
            hdPhiMCGenMirror[iso][iSh][iGenType][izt][iptTr] = new TH1F(Form("hdPhiMCGenMirror%s", sTypeMCGenType[iGenType].Data()) + sIso + sShSh + sZtBin + sPtTrig, Form("hdPhiMCGenMirror%s", sTypeMCGenType[iGenType].Data()) + sIso + sShSh + sZtBin + sPtTrig, nbinsX, 0, TMath::Pi());
            hdPhiSamMCRecMirror[iso][iSh][iGenType][izt][iptTr] = new TH1F(Form("hdPhiSameMCRecMirror%s", sTypeMCGenType[iGenType].Data()) + sIso + sShSh + sZtBin + sPtTrig, Form("hdPhiSameMCRecMirror%s", sTypeMCGenType[iGenType].Data()) + sIso + sShSh + sZtBin + sPtTrig, nbinsX, 0, TMath::Pi());
          cout<<"nbinsX: "<<nbinsX<<endl;
          }
          h3MixMCRec[iSh]->SetAxisRange(ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1] - 0.0001, "X");
          h2DdPhidEtaMixMCRec[iso][iSh][izt][iptTr] = (TH2F *)h3MixMCRec[iSh]->Project3D("zy");
          h2DdPhidEtaMixMCRec[iso][iSh][izt][iptTr]->SetName("h2DdPhidEtaMixMCRec" + sIso + sShSh + sZtBin + sPtTrig);
          hdPhiMixMCRec[iso][iSh][izt][iptTr] = (TH1F *)h2DdPhidEtaMixMCRec[iso][iSh][izt][iptTr]->ProjectionX("hdPhiMixMCRec" + sIso + sShSh + sZtBin + sPtTrig);
          hdPhiMixMCRec[iso][iSh][izt][iptTr]->Scale(1. / hTriggerMixMCRec[iso][iSh]->Integral(hTriggerMixMCRec[iso][iSh]->FindBin(ptTrig[index1 + iptTr]), hTriggerMixMCRec[iso][iSh]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001)));
          nbinsX = hdPhiMixMCRec[iso][iSh][izt][iptTr]->GetNbinsX() / 2;
          hdPhiMixMCRecMirror[iso][iSh][izt][iptTr] = new TH1F("hdPhiMixMCRecMirror" + sIso + sShSh + sZtBin + sPtTrig, "hdPhiMixMCRecMirror" + sIso + sShSh + sZtBin + sPtTrig, nbinsX, 0, TMath::Pi());
cout<<"nbinsXMix: "<<nbinsX<<endl;
          if (bMirror)
          {
            Mirroring(hdPhiMixMCRec[iso][iSh][izt][iptTr], hdPhiMixMCRecMirror[iso][iSh][izt][iptTr]);
            hdPhiMixMCRecMirror[iso][iSh][izt][iptTr]->Rebin(5);
            for (int iGenType = 0; iGenType < nGenType; iGenType++)
            {
              Mirroring(hdPhiMCGen[iso][iSh][iGenType][izt][iptTr], hdPhiMCGenMirror[iso][iSh][iGenType][izt][iptTr]);
              hdPhiMCGenMirror[iso][iSh][iGenType][izt][iptTr]->Rebin(5);
              hdPhiMCGenMirrorUE[iso][iSh][iGenType][izt][iptTr] = (TH1F *)hdPhiMCGenMirror[iso][iSh][iGenType][izt][iptTr]->Clone(Form("hdPhiMCGenMirrorUE%s", sTypeMCGenType[iGenType].Data()) + sIso + sShSh + sZtBin + sPtTrig);
              fZYAM(hdPhiMCGenMirror[iso][iSh][iGenType][izt][iptTr]);
              cout << hdPhiMCGen[iso][iSh][iGenType][izt][iptTr] << endl;

              Mirroring(hdPhiSamMCRec[iso][iSh][iGenType][izt][iptTr], hdPhiSamMCRecMirror[iso][iSh][iGenType][izt][iptTr]);
              hdPhiSamMCRecMirror[iso][iSh][iGenType][izt][iptTr]->Rebin(5);

              hdPhiSamMCRecNoUE[iso][iSh][iGenType][izt][iptTr] = (TH1F *)hdPhiSamMCRecMirror[iso][iSh][iGenType][izt][iptTr]->Clone(Form("hdPhiSameMCRecNoUE%s", sTypeMCGenType[iGenType].Data()) + sIso + sShSh + sZtBin + sPtTrig);
              hdPhiSamMCRecNoUE[iso][iSh][iGenType][izt][iptTr]->Sumw2();

              // double scaleFactMC = fZYAM_Mix(hdPhiSamMCRecMirror[iso][iSh][izt][iptTr], hdPhiMixMCRecMirror[iso][iSh][izt][iptTr]);
              if (bZYAM)
              {
                cout << "Mirror - ZYAM" << endl;
                cout << "ZYAM" << endl;
                fZYAM(hdPhiSamMCRecNoUE[iso][iSh][iGenType][izt][iptTr]);
              }
              else
              {
                cout << "Mirror - MIXED EVENT" << endl;
                cout << hdPhiSamMCRecNoUE[iso][iSh][iGenType][izt][iptTr] << endl;
                hdPhiSamMCRecNoUE[iso][iSh][iGenType][izt][iptTr]->Add(hdPhiMixMCRecMirror[iso][iSh][izt][iptTr], -1);
              }
            }
          }
          else
          {
            for (int iGenType = 0; iGenType < nGenType; iGenType++)
            {
              hdPhiMCGen[iso][iSh][iGenType][izt][iptTr]->Rebin(5);
              fZYAM(hdPhiMCGen[iso][iSh][iGenType][izt][iptTr]);
              cout << hdPhiMCGen[iso][iSh][iGenType][izt][iptTr] << endl;
              hdPhiSamMCRecNoUE[iso][iSh][iGenType][izt][iptTr] = (TH1F *)hdPhiSamMCRec[iso][iSh][iGenType][izt][iptTr]->Clone(Form("hdPhiSameMCRecNoUE%s", sTypeMCGenType[iGenType].Data()) + sIso + sShSh + sZtBin + sPtTrig);
              hdPhiSamMCRecNoUE[iso][iSh][iGenType][izt][iptTr]->Sumw2();
              // double scaleFactMC = fZYAM_Mix(hdPhiSamMCRecMirror[iso][iSh][izt][iptTr], hdPhiMixMCRecMirror[iso][iSh][izt][iptTr]);
              if (bZYAM)
              {
                cout << "NO Mirroring - ZYAM" << endl;
                cout << "ZYAM" << endl;
                fZYAM(hdPhiSamMCRecNoUE[iso][iSh][iGenType][izt][iptTr]);
              }
              else
              {
                cout << "NO Mirroring - MIXED EVENT" << endl;
                cout << hdPhiSamMCRecNoUE[iso][iSh][iGenType][izt][iptTr] << endl;
                hdPhiSamMCRecNoUE[iso][iSh][iGenType][izt][iptTr]->Add(hdPhiMixMCRec[iso][iSh][izt][iptTr], -1);
              }
            }
          }
          cout << "CiaoCiao" << endl;
          fOutPut->cd();
          for (int iGenType = 0; iGenType < nGenType; iGenType++)
          {
            hdPhiMCGen[iso][iSh][iGenType][izt][iptTr]->Write();
            hdPhiMCGenMirrorUE[iso][iSh][iGenType][izt][iptTr]->Write();
            hdPhiMCGenMirror[iso][iSh][iGenType][izt][iptTr]->Write();
            hdPhiSamMCRec[iso][iSh][iGenType][izt][iptTr]->Write();
            hdPhiSamMCRecMirror[iso][iSh][iGenType][izt][iptTr]->Write();
            hdPhiSamMCRecNoUE[iso][iSh][iGenType][izt][iptTr]->Write();
          }
          hdPhiMixMCRec[iso][iSh][izt][iptTr]->Write();
          hdPhiMixMCRecMirror[iso][iSh][izt][iptTr]->Write();
          cout << "CiaoCiao1" << endl;
        }
        cout << "CiaoCiao2" << endl;
      }
      cout << "CiaoCiao3" << endl;
    }
    cout << "CiaoCiao4" << endl;
  }

  for (int iSh = 0; iSh < nShSh; iSh++)
  {
    TString sShSh = Form("_ShSh%s", shshString[iSh].Data());
    for (int iso = 0; iso < nIso; iso++)
    {
      for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
      {
        cout << nPtTrig << endl;
        cout << index2 - index1 << endl;
        TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
        for (int iGenType = 0; iGenType < nGenType; iGenType++)
        {
          hZtIsoGammaPtBinMCGen[iso][iSh][iGenType][iptTr] = new TH1F(Form("hZtIso%dGammaPtBinMCGen%s_%s", iso, sTypeMCGenType[iGenType].Data(), sPtTrig.Data()) + sShSh, Form("hZtIso%dGammaPtBinMCGen%s_%s", iso, sTypeMCGenType[iGenType].Data(), sPtTrig.Data()) + sShSh, nZtBin, assocZt);

          hZtIsoGammaPtBinMCRec[iso][iSh][iGenType][iptTr] = new TH1F(Form("hZtIso%dGammaPtBinMCRec%s_%s", iso, sTypeMCGenType[iGenType].Data(), sPtTrig.Data()) + sShSh, Form("hZtIso%dGammaPtBinMCRec%s_%s", iso, sTypeMCGenType[iGenType].Data(), sPtTrig.Data()) + sShSh, nZtBin, assocZt);

          for (Int_t izt = 0; izt < nZtBin; izt++)
          {
            // IsoGammaMC
            if (!bMirror)
            {
              ZtFunction(hdPhiMCGen[iso][iSh][iGenType][izt][iptTr], hZtIsoGammaPtBinMCGen[iso][iSh][iGenType][iptTr], izt, phiMin, phiMax);
            }
            else if (bMirror)
            {
              ZtFunction(hdPhiMCGenMirror[iso][iSh][iGenType][izt][iptTr], hZtIsoGammaPtBinMCGen[iso][iSh][iGenType][iptTr], izt, phiMin, phiMax);
            }
            ZtFunction(hdPhiSamMCRecNoUE[iso][iSh][iGenType][izt][iptTr], hZtIsoGammaPtBinMCRec[iso][iSh][iGenType][iptTr], izt, phiMin, phiMax);
          }

          ScaleBinBySize(hZtIsoGammaPtBinMCGen[iso][iSh][iGenType][iptTr]);
          PlotStyle(hZtIsoGammaPtBinMCGen[iso][iSh][iGenType][iptTr], 20, 1, kBlue - 4, 0, "zt", "1/N dN/dzt", false);

          ScaleBinBySize(hZtIsoGammaPtBinMCRec[iso][iSh][iGenType][iptTr]);
          PlotStyle(hZtIsoGammaPtBinMCRec[iso][iSh][iGenType][iptTr], 21, 1, kOrange + 7, 0, "zt", "1/N dN/dzt", false);

          fOutPut->cd();
          hZtIsoGammaPtBinMCGen[iso][iSh][iGenType][iptTr]->Write();
          hZtIsoGammaPtBinMCRec[iso][iSh][iGenType][iptTr]->Write();
        }
      }
      for (int iGenType = 0; iGenType < nGenType; iGenType++)
      {
        hZtIsoGammaMCGen[iso][iSh][iGenType] = new TH1F(Form("hZtIso%dGammaMCGen%s%s", iso, sTypeMCGenType[iGenType].Data(), sPtAll.Data()) + sShSh, Form("hZtIso%dGammaMCGen%s%s", iso, sTypeMCGenType[iGenType].Data(), sPtAll.Data()) + sShSh, nZtBin, assocZt);
        hZtIsoGammaMCGen[iso][iSh][iGenType] = SumPtBinXzt(hTriggerMCGen[iso][iSh][iGenType], ptTrig, index1, index2, hZtIsoGammaPtBinMCGen[iso][iSh][iGenType], hZtIsoGammaMCGen[iso][iSh][iGenType], histPur[iCen], funcPur[iCen], systPur, false);

        hZtIsoGammaMCRec[iso][iSh][iGenType] = new TH1F(Form("hZtIso%dGammaMCRec%s%s", iso, sTypeMCGenType[iGenType].Data(), sPtAll.Data()) + sShSh, Form("hZtIso%dGammaMCRec%s%s", iso, sTypeMCGenType[iGenType].Data(), sPtAll.Data()) + sShSh, nZtBin, assocZt);
        hZtIsoGammaMCRec[iso][iSh][iGenType] = SumPtBinXzt(hTriggerSamMCRec[iso][iSh][iGenType], ptTrig, index1, index2, hZtIsoGammaPtBinMCRec[iso][iSh][iGenType], hZtIsoGammaMCRec[iso][iSh][iGenType], histPur[iCen], funcPur[iCen], systPur, false);

        fOutPut->cd();

        hZtIsoGammaMCGen[iso][iSh][iGenType]->Write();
        hZtIsoGammaMCRec[iso][iSh][iGenType]->Write();
      }
    }
  }

  /////////////////////////////////////////////////////////////////
  ///////////////////// Plotting section /////////////////////////
  ///////////////////////////////////////////////////////////////

  TLegend *dPhiGenleg[nPtTrig];
  for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
  {
    TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
    dPhiGenleg[iptTr] = LegStd(dPhiGenleg[iptTr], 0.5, 0.5, 0.8, 0.8);
    TCanvas *cSameNoUEGenMC = new TCanvas(Form("cSameNoUEGenMC_%s", sPtTrig.Data()), Form("cSameNoUEGenMC_%s", sPtTrig.Data()), 4 * 800, 2 * 600);
    cSameNoUEGenMC->Divide(4, 2);
    TCanvas *cSameNoUEMC = new TCanvas(Form("cSameNoUEMC_%s", sPtTrig.Data()), Form("cSameNoUEMC_%s", sPtTrig.Data()), 4 * 800, 2 * 600);
    cSameNoUEMC->Divide(4, 2);
    TCanvas *cMixMC = new TCanvas(Form("cMixMC_IsoNotIso_%s", sPtTrig.Data()), Form("cMixMC_IsoNotIso_%s", sPtTrig.Data()), 4 * 800, 2 * 600);
    cMixMC->Divide(4, 2);
    TCanvas *cSameMixMC = new TCanvas(Form("cSameMixMC_%s", sPtTrig.Data()), Form("cSameMixMC_%s", sPtTrig.Data()), 4 * 800, 2 * 600);
    cSameMixMC->Divide(4, 2);
    TCanvas *cSameNoUEMCMixvsZYAM = new TCanvas(Form("cSameNoUEMCMixvsZYAM_%s", sPtTrig.Data()), Form("cSameNoUEMCMixvsZYAM_%s", sPtTrig.Data()), 4 * 800, 2 * 600);
    cSameNoUEMCMixvsZYAM->Divide(4, 2);
    cout << "testttttt" << endl;
    for (Int_t izt = 0; izt < nZtBin; izt++)
    {
      TString sTitle = Form("%2.0f < #it{p}_{T}^{tr} < %2.0f GeV/#it{c}, %2.2f < #it{z}_{T}^{as} < %2.2f ", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], assocZt[izt], assocZt[izt + 1]);
      cSameNoUEGenMC->cd(izt + 1);
      for (int iso = 0; iso < nIso; iso++)
      {
        for (int iGenType = 0; iGenType < nGenType; iGenType++)
        {
          PlotStyle(hdPhiSamMCRecMirror[1][0][iGenType][izt][iptTr], 20, 2, kBlue + 1, 0, "#Delta#varphi (rad)", "1/N^{trig}dN/d#varphi", false);
          PlotStyle(hdPhiSamMCRecNoUE[1][0][iGenType][izt][iptTr], 25, 2, kOrange + 3, 0, "#Delta#varphi (rad)", "1/N^{trig}dN/d#varphi", false);

          hdPhiSamMCRecMirror[1][0][iGenType][izt][iptTr]->SetTitle(sTitle);
          hdPhiSamMCRecNoUE[1][0][iGenType][izt][iptTr]->SetTitle(sTitle);
          PlotStyle(hdPhiMCGenMirror[iso][0][iGenType][izt][iptTr], kMarkStyleIso_NotIso[iso], 2, kColorMC[iso][iGenType], 0, "#Delta#varphi (rad)", "1/N^{trig}dN/d#varphi", false);

          hdPhiMCGenMirror[0][0][0][izt][iptTr]->SetTitle(sTitle);
          hdPhiMCGenMirror[iso][0][iGenType][izt][iptTr]->Draw("same");
          // hdPhiMCGenMirror[iso][0][1][izt][iptTr]->Draw("same");
          // hdPhiMCGenMirror[iso][0][2][izt][iptTr]->Draw("same");
          cSameNoUEMC->cd(izt + 1);
          hdPhiSamMCRecMirror[1][0][iGenType][izt][iptTr]->GetYaxis()->SetRangeUser(-0.8 * hdPhiSamMCRecMirror[1][0][iGenType][izt][iptTr]->GetMinimum(), 1.1 * hdPhiSamMCRecMirror[1][0][iGenType][izt][iptTr]->GetMaximum());
          hdPhiSamMCRecMirror[1][0][iGenType][izt][iptTr]->Draw("same");
          hdPhiSamMCRecNoUE[1][0][iGenType][izt][iptTr]->Draw("same");

          cSameNoUEMCMixvsZYAM->cd(izt + 1);
          hdPhiSamMCRecNoUE[1][0][iGenType][izt][iptTr]->Draw("same");
        }
      }
    }
    cSameNoUEGenMC->cd(7);
    for (int iso = 0; iso < nIso; iso++)
    {
      for (int iGenType = 0; iGenType < nGenType; iGenType++)
      {
        dPhiGenleg[iptTr]->AddEntry(hdPhiMCGenMirror[iso][0][iGenType][0][iptTr], Form("Iso%d %s", iso, sTypeMCGenType[iGenType].Data()), "lep");
      }
    }
    dPhiGenleg[iptTr]->Draw("same");

    cSameNoUEGenMC->Print(Form("%s/cSameNoUEGeneratePhoton_Pi0_Pi0DecayMC_%s.pdf", dirFiles.Data(), sPtTrig.Data()));
    cSameNoUEMCMixvsZYAM->Print(Form("%s/cSameNoUEMCMixvsZYAM_%s.pdf", dirFiles.Data(), sPtTrig.Data()));
  }
  TLegend *zTleg = LegStd(zTleg, 0.6, 0.4, 0.8, 0.7);
  TCanvas *cZtIsoGamma = new TCanvas(Form("cZtIsoGamma"), Form("cZtIsoGamma"), 800, 600);
  gPad->SetLogy();
  /*
  for (int iso = 0; iso < nIso; iso++)
  {
    for (int iGenType = 0; iGenType < nGenType; iGenType++)
    {
      PlotStyle(hZtIsoGammaMCGen[iso][iSh][iGenType], kMarkStyleIso_NotIso[iso], 1, kColorMC[iso][iGenType], 0, "#it{z}_{#rm T}", "1/N^{trig}dN/d#it{z}_{T}", true);
      hZtIsoGammaMCGen[iso][iSh][iGenType]->Draw("same");
      zTleg->AddEntry(hZtIsoGammaMCGen[iso][iSh][iGenType], Form("Gen Iso%d ShSh %s , %s", iso, shshString[iSh].Data(), sTypeMCGenType[iGenType].Data()), "lep");
    PlotStyle(hZtIsoGammaMCRec[iso][iSh][iGenType], 20, 1, kColorMC[iso][iGenType], 0, "#it{z}_{#rm T}", "1/N^{trig}dN/d#it{z}_{T}", true);
    hZtIsoGammaMCRec[iso][iSh][iGenType]->Draw("same");
    zTleg->AddEntry(hZtIsoGammaMCRec[iso][iSh][iGenType], Form("Rec Iso%d ShSh %s, %s ", iso, shshString[iSh].Data(), sTypeMCGenType[iGenType].Data()), "lep");
    }
  }
*/

}

void IsoGammaHadronMCJJlow(float ptTrMin = 18, float ptTrMax = 40, TString sFileDirShSig = "RootFiles/DataSh100_AssocPt500", Bool_t bMirror = true, TString shshBkg = "0.40-2.00", TString dirFiles = "ResultsStudyMCJJpi0test", double systPur = 1, bool bZYAM = true, double phiMin = TMath::Pi() * 3 / 5, double phiMax = TMath::Pi(), bool GJppOn = true)
{
  TString tagFile[nCen];
  fileMC = TFile::Open(Form("RootFiles/Pi0Hadron/MCJJ_Pi0.root"));
  if (!fileMC)
    cout << "MC File doesn't exist" << endl;

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    tagFile[iCen] = Form("EMCAL_MB_%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    fileDataShStd = TFile::Open(Form("%s/%s.root", sFileDirShSig.Data(), tagFile[iCen].Data()));
    if (!fileDataShStd)
      cout << "Data File doesn't exist" << endl;

    cout << iCen << endl;
    Exec(ptTrMin, ptTrMax, iCen, bMirror, shshBkg, dirFiles, systPur, bZYAM, phiMin, phiMax, GJppOn);
  }
  fileDataShStd->Close();
  delete fileDataShStd;
  fileMC->Close();
  delete fileMC;
  // fileMC[1]->Close();
  // delete fileMC[1];
}

void Mirroring(TH1F *hMir, TH1F *hMirXtrue)
{
  int bin0 = hMir->FindBin(0 + 0.00001);
  int binPi = hMir->FindBin(TMath::Pi());
  int totLeft = hMir->GetNbinsX();
  /*cout << "bin0: " << bin0 << "Bin Content: " << hMir->GetBinContent(bin0) << endl;
  cout << "binPi: " << binPi << "Bin Content: " << hMir->GetBinContent(binPi) << endl;
  cout << "binALL: " << totLeft << endl;
  cout << hMir->GetBinCenter(hMir->FindBin(0)) << endl;
  cout << hMir->GetBinCenter(bin0) << endl;
  cout << hMir->GetBinCenter(binPi) << endl;
  cout << hMir->GetBinWidth(binPi) << endl;*/
  cout << hMir->GetNbinsX() << endl;

  for (int jbin = 0; jbin < bin0; jbin++) // calcolo bin(pi+j)+bin(pi-j) e j parte da 1
  {
    cout << "bin: " << bin0 + jbin << ", content1: " << hMir->GetBinContent(bin0 + jbin) << "bin: " << bin0 - jbin - 1 << ", bincontent2:" << hMir->GetBinContent(bin0 - jbin - 1) << endl;
    double binCont = (hMir->GetBinContent(bin0 + jbin)) + (hMir->GetBinContent(bin0 - jbin - 1));
    double binErr = sqrt((hMir->GetBinError(bin0 + jbin)) * (hMir->GetBinError(bin0 + jbin)) + (hMir->GetBinError(bin0 - jbin - 1)) * (hMir->GetBinError(bin0 - jbin - 1)));
    cout << "bin:" << (bin0 + jbin) << " content:" << binCont << endl;
    // hMir->SetBinContent((bin0 + jbin), binCont);
    // hMir->SetBinContent(bin0 - jbin - 1, 0);
    hMirXtrue->SetBinContent((jbin + 1), binCont);
    hMirXtrue->SetBinError((jbin + 1), binErr);
  }
  cout << "_________________" << endl;
  cout << "binPi: " << binPi << endl;
  for (int jbin = 0; jbin < (totLeft - binPi); jbin++)
  {
    cout << "bin: " << binPi - jbin << ", bincontent1:" << hMir->GetBinContent(binPi - jbin) << "bin: " << binPi + jbin + 1 << ", bin content2: " << hMir->GetBinContent(binPi + jbin + 1) << endl;
    double binCont = hMir->GetBinContent(binPi + jbin + 1) + hMir->GetBinContent(binPi - jbin);
    double binErr = sqrt((hMir->GetBinError(binPi + jbin + 1)) * (hMir->GetBinError(binPi + jbin + 1)) + (hMir->GetBinError(binPi - jbin)) * (hMir->GetBinError(binPi - jbin)));
    cout << "bin:" << (binPi - jbin) << " content:" << binCont << endl;
    // hMir->SetBinContent((binPi - jbin), binCont);
    // hMir->SetBinContent(binPi + jbin + 1, 0);
    hMirXtrue->SetBinContent((hMirXtrue->GetNbinsX() - jbin), binCont);
    hMirXtrue->SetBinError((hMirXtrue->GetNbinsX() - jbin), binErr);
  }
}
void ZtFunction(TH1F *hDeltaPhi, TH1F *hZT, int bin, double phiMin, double phiMax)
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
static void ScaleBinBySize(TH1F *h)
{
  for (Int_t ibin = 1; ibin <= h->GetNbinsX(); ibin++)
  {
    Double_t width = h->GetBinWidth(ibin);
    Double_t content = h->GetBinContent(ibin);
    Double_t error = h->GetBinError(ibin);
    // printf("bin %d, width %f, content %e\n",ibin,width,content);
    //  cout<<h->GetNbinsX()<<endl;
    h->SetBinContent(ibin, content / width);
    h->SetBinError(ibin, error / width);
  }
}

/*void PlotStyle(TH1F *hPlot, int kMarker, int kMarkerSize, int kColor, TString titleX, TString titleY)
{
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111111);
  // gStyle->SetPadTopMargin(0.10);
  hPlot->SetMarkerStyle(kMarker);
  hPlot->SetMarkerSize(kMarkerSize);
  hPlot->SetMarkerColor(kColor);
  hPlot->SetLineColor(kColor);
  hPlot->GetYaxis()->SetTitle(Form("%s", titleY.Data()));
  hPlot->GetXaxis()->SetTitle(Form("%s", titleX.Data()));
  // leg->SetFillColor(kWhite);
  // leg->SetLineColor(0);
}*/
TH1F *SumPtBinXzt(TH1F *hTrigSame, Float_t PtTrigger[npt], int index1, int index2, TH1F *hzTbin[npt], TH1F *hzTbinAll, TH1F *hPur, TF1 *fPur, double systPur, Bool_t bData)
{
  Float_t intPtAll = 0;
  // cout << "(1/NintPur_tot) * Sum(pur * NInt * f(Zt))" << endl;
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

void fZYAM(TH1F *hSame, double rangeMin = 3 * (TMath::Pi()) / 10, double rangeMax = TMath::Pi() / 2)
{
  int minBin = hSame->GetXaxis()->FindBin(rangeMin + 0.00001);
  int maxBin = hSame->GetXaxis()->FindBin(rangeMax);
  double IntSame = hSame->Integral(minBin, maxBin);
  cout << "bincontent min: " << hSame->GetBinCenter(minBin) << ", bincontent max: " << hSame->GetBinCenter(maxBin) << endl;
  printf("bin  %d %d \n", minBin, maxBin);
  double CorrFact = IntSame / (maxBin - minBin + 1);

  for (Int_t ibin = 1; ibin <= hSame->GetNbinsX(); ibin++)
  {
    double binCont = (hSame->GetBinContent(ibin) - CorrFact);
    // printf("Final Val %f, preval %f, int %f \n", binCont, hSame->GetBinContent(ibin),IntSame);
    hSame->SetBinContent(ibin, binCont);
  }
}

double fZYAM_Mix(TH1F *hSame, TH1F *hMix)
{
  double minBin = hSame->GetXaxis()->FindBin((3 * (TMath::Pi()) / 10) + 0.00001);
  double maxBin = hSame->GetXaxis()->FindBin(TMath::Pi() / 2);
  double IntSame = hSame->Integral(minBin, maxBin);
  printf("norm Same %f \n", IntSame);
  minBin = hMix->GetXaxis()->FindBin((3 * (TMath::Pi()) / 10) + 0.00001);
  maxBin = hMix->GetXaxis()->FindBin(TMath::Pi() / 2);
  double IntMix = hMix->Integral(minBin, maxBin);
  double CorrFact = 1;
  if (IntMix != 0 && IntSame != 0)
  {
    printf("norm Mix %f \n", IntMix);
    CorrFact = IntSame / IntMix;
    printf("CorrFactor %f, Name: %s \n", CorrFact, hMix->GetName());
    hMix->Scale(CorrFact);
  }
  return CorrFact;
}
