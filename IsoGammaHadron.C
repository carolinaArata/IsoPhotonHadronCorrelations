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

int const nCen = 4;
int cenBins[] = {0, 10, 30, 50, 90}; // definition of different centralities

int nIso = 2;        // isolation: not iso = 0; iso = 1
int const nShSh = 2; // shower shape ShSh: Signal = 0.10-0.30; Bkg = x
int nZtBinThin = 10; // all zT intervals
double assocZtThinner[] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20};

int nZtBin = 6; // larger zT intervals defined for QM
double assocZt[] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 1.00};
int npt = 12; // pT trigger intervals
float ptTrig[] = {10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 50, 60, 80};

TString sShShNCentMix;
bool systNMix = false;

TFile *fileData = 0;
TFile *fileDataShStd = 0;
TFile *fileDataShBkg = 0;
TFile *fileDataMix = 0;
TFile *fileMC[nShSh] = {0, 0};
TH1F *histPur[nCen];
TH1F *histPurStat[nCen];
TF1 *funcPur[nCen];

void Mirroring(TH1F *hMir, TH1F *hMirXtrue);                                          // function for mirroring
void PlotStyle(TH1F *hPlot, int kMarker, int kColor, TString titleX, TString titleY); // plot aesthetic
void ZtFunction(TH1F *hDeltaPhi, TH1F *hZT, int bin, double phiMin, double phiMax);   // calculation of zT function
static void ScaleBinBySize(TH1F *h);
TH1F *SumPtBinXzt(TH1F *hTrigSame, Float_t PtTrigger[npt], int index1, int index2, TH1F *hzTbin[npt], TH1F *hzTbinAll, TH1F *hPur, TF1 *fPur, double systPur, Bool_t bData);
void fZYAM(TH1F *hSame, double rangeMin = 3 * (TMath::Pi()) / 10, double rangeMax = TMath::Pi() / 2);
double fZYAM_Mix(TH1F *hSame, TH1F *hMix);

void Exec(float ptMin = 18, float ptMax = 20, int iCen = 0, bool bMirror = true, TString shshBkg = "0.40-1.00", TString dirFiles = "ResultsProvaUnisci", double systPur = 1, bool bZYAM = false, bool bPlot = false, double phiMin = TMath::Pi() * 3 / 5, double phiMax = TMath::Pi(), bool systShSh = false, bool systNMix = false)
{
  TString sHistName = "AnaPhotonHadronCorr_";
  TString shshString[2] = {"0.10-0.30", shshBkg};
  TString shshStringMC[2] = {"0.10-0.30", "0.10-0.30"};

  // index pT trigger start and stop
  int nsize = sizeof(ptTrig) / sizeof(ptTrig[0]);
  auto itr1 = find(ptTrig, ptTrig + nsize, ptMin);
  auto itr2 = find(ptTrig, ptTrig + nsize, ptMax);
  int index1 = distance(ptTrig, itr1);
  int index2 = distance(ptTrig, itr2);
  int nPtTrig = index2 - index1;
  cout << index1 << "___" << index2 << ", " << nPtTrig << endl;

  /////////////////////////////////////////////////////////////
  ////// Definition of histograms used in the analysis ///////
  ///////////////////////////////////////////////////////////

  TH1F *hTriggerSam[nIso][nShSh];                     // PtTrig distrib same
  TH1F *hTriggerMix[nIso][nShSh];                     // PtTrig distrib mixed
  TH2F *h2DdPhidEtaSam[nIso][nShSh][nZtBin][nPtTrig]; // 2D correlation distributions deltaPhi and deltaEta same
  TH2F *h2DdPhidEtaMix[nIso][nShSh][nZtBin][nPtTrig]; // 2D correlation distributions deltaPhi and deltaEta mix
  TH1F *hdPhiSam[nIso][nShSh][nZtBin][nPtTrig];       // 1D correlation distributions deltaPhi same
  TH1F *hdPhiMix[nIso][nShSh][nZtBin][nPtTrig];       // 1D correlation distributions deltaPhi mix
  // SameEvent and MixedEvent Mirror
  TH1F *hdPhiSamMirror[nIso][nShSh][nZtBin][nPtTrig]; // 1D correlation distributions deltaPhi same Mirrored
  TH1F *hdPhiMixMirror[nIso][nShSh][nZtBin][nPtTrig]; // 1D correlation distributions deltaPhi mix Mirrored

  double scaleFact[nIso][nShSh][nZtBin][nPtTrig]; // scaling factor for mixed event
  // NO UE (Same-Mix)
  TH1F *hdPhiSamNoUE[nIso][nShSh][nZtBin][nPtTrig]; // Same - Mix for every pT trig bin
  TH1F *hdPhiSamNoUERatio[nIso][nZtBin][nPtTrig];
  TH1F *hdPhiSamNoUEPtAll[nIso][nShSh][nZtBin]; // Same - Mix for LARGE pT trig range
  // Purity correction
  TH1F *hdPhiSamPur[nIso][nZtBin][nPtTrig]; // IsoClusterGamma wo UE, subtract (1-p)IsoClPi0 and divided by Purity
  TH1F *hdPhiSamPi0Pur[nIso][nZtBin][nPtTrig];
  TH1F *hdPhiMixPur[nIso][nZtBin][nPtTrig];
  TH1F *hdPhiMixPi0Pur[nIso][nZtBin][nPtTrig];
  // DeltaPhi for IsoGamma
  TH1F *hdPhiPhoton[nIso][nZtBin][nPtTrig];      // iso/not iso Gamma correlation distributions
  TH1F *hdPhiPhotonPtBin[nIso][nZtBin][nPtTrig]; // iso/not iso Gamma correlation distributions for combining pT trig bin
  TH1F *hdPhiPhotonPtAll[nIso][nZtBin];          // iso/not iso Gamma correlation distributions for LARGE pT trig bin
  // DeltaPhi for IsoPi0
  TH1F *hdPhiPi0[nIso][nZtBin][nPtTrig];      // iso/not iso pi0 correlation distributions
  TH1F *hdPhiPi0PtBin[nIso][nZtBin][nPtTrig]; // iso/not iso pi0 correlation distributions for combining pT trig bin
  TH1F *hdPhiPi0PtAll[nIso][nZtBin];          // iso/not iso pi0 correlation distributions for LARGE pT trig bin
  // DeltaPhi for NotIsoPi0
  TH1F *hdPhiNotIsoPi0[nZtBin][nPtTrig];
  TH1F *hdPhiNotIsoPi0PtBin[nZtBin][nPtTrig];
  TH1F *hdPhiNotIsoPi0PtAll[nZtBin];
  // Zt Distribution
  TH1F *hZtPhotonPtBin[nIso][nPtTrig]; // Zt distribution iso/not iso Gamma X every PtBin (only purity correction)
  TH1F *hZtPhoton[nIso];               // Zt distribution iso/not iso Gamma full range (only purity correction)
  TH1F *hZtPi0PtBin[nIso][nPtTrig];    // Zt distribution iso/not iso Pi0 X every PtBin (only purity correction)
  TH1F *hZtPi0[nIso];                  // Zt distribution iso/not iso Pi0 full range (only purity correction)

  // Purity root file definition
  TFile *fPurity = new TFile("~/work/histogram/IsoPhotonHadronCorrelations/Purity_IsoSig1.5_M02Sig0.10-0.30_IsoBkg_4.0_25.0_M02Bkg0.40_2.00_LHC15o_18qr_L1MB.root");
  histPur[iCen] = (TH1F *)fPurity->Get(Form("Purity_Cen%d_R0.2_Sys", iCen));
  histPurStat[iCen] = (TH1F *)fPurity->Get(Form("Purity_Cen%d_R0.2", iCen));
  funcPur[iCen] = histPur[iCen]->GetFunction("purityFitCombinedSigmoid");

  // Define directory where you will save the root files and the root files for every centrality
  TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
  TString sPtAll = Form("_Pt%2.0f_%2.0f", ptMin, ptMax);

  gSystem->Exec(Form("mkdir %s", dirFiles.Data()));
  TFile *fOutPut = new TFile(Form("%s/fPlot%s%s%s.root", dirFiles.Data(), shshBkg.Data(), sCent.Data(), sPtAll.Data()), "RECREATE");
  cout << fOutPut->GetName() << endl;

  cout << "Centrality: " << iCen << " " << cenBins[iCen] << "-" << cenBins[iCen + 1] << endl;
  if (bZYAM)
    cout << "ZYAM UE subtraction has been chosen" << endl;

  cout << "Mixed Event UE subtraction has been chosen" << endl;
  for (int iso = 0; iso < nIso; iso++)
  {
    TString sIso = Form("Iso%d", iso);

    for (int iSh = 0; iSh < nShSh; iSh++)
    {
      TString sShSh = Form("_ShSh%s", shshString[iSh].Data());
      if (iSh == 0)
        fileData = fileDataShStd;
      else if (iSh == 1 && iso == 0 && systShSh == true) // this condition is necessary because the train for the multi shsh background were run for Isolation only
      {
        fileData = fileDataShStd;
        sShSh = Form("_ShSh0.40-1.00"); // use standard
        cout << fileData->GetName() << endl;
      }
      else if (iSh == 1 && iso == 1 && systShSh == true) // this condition is necessary because the train for the multi shsh background were run for Isolation only
      {
        fileData = fileDataShBkg;
        cout << fileData->GetName() << endl;
      }
      cout << "Getter ptTrig histograms for Same and Mixed:" << sIso << sShSh << endl;
      hTriggerSam[iso][iSh] = (TH1F *)fileData->Get(sHistName + sIso + sShSh + sCent + "_hPtTrigger"); // pT trig distribution X same
      cout << fileData->GetName() << endl;
      cout << fileData << endl;
      cout << hTriggerSam[iso][iSh] << endl;
      if (!systNMix)
        hTriggerMix[iso][iSh] = (TH1F *)fileData->Get(sHistName + sIso + sShSh + sCent + "_hPtTriggerMixed"); // pT trig distribution X mixed
      else if (iSh==0 && systNMix)
        hTriggerMix[iso][iSh] = (TH1F *)fileDataMix->Get(sHistName + sIso + sShSh + sCent + "_hPtTriggerMixed"); // pT trig distribution X mixed
      else if (iSh==1 && systNMix)
      {
        hTriggerMix[iso][iSh] = (TH1F *)fileDataMix->Get(sHistName + sIso + sShShNCentMix + sCent + "_hPtTriggerMixed"); // pT trig distribution X mixed
        cout<<"hellooooo"<<sShShNCentMix<<endl;
      }
      if (bPlot)
      {
        cout << "Plot pT trig distributions" << endl;
        TCanvas *cIsoPtTrig = new TCanvas(Form("Iso%dPtTrig" + sShSh + sCent, iso), Form("Iso%dPtTrig" + sShSh + sCent, iso), 800, 600);
        hTriggerSam[iso][iSh]->SetDirectory(0);
        hTriggerMix[iso][iSh]->SetDirectory(0);
        hTriggerMix[iso][iSh]->SetLineColor(kRed);
        hTriggerMix[iso][iSh]->Draw("same");
        hTriggerSam[iso][iSh]->Draw("same");
      }
      fOutPut->cd();
      hTriggerSam[iso][iSh]->Write();
      if (systNMix)
        hTriggerMix[iso][iSh]->SetName(sHistName + sIso + sShSh + sCent + "_hPtTriggerMixed");
      
      hTriggerMix[iso][iSh]->Write();

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////// Get TH3F correlation distributions for every zT bin obtained from the train and combined in the new zT bins ////////
      ////// The "Next" in the files name are used for combining the zT intervals.                                      ////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      cout << "Getter 3D azimuthal distributions: " << "iso = " << sIso << ", ShSh = " << sShSh << endl;
      for (Int_t izt = 0; izt < nZtBin; izt++)
      {
        TH3F *h3Sam;
        TH3F *h3SamNext;
        TH3F *h3SamNext1;
        TH3F *h3SamNext2;
        TH3F *h3Mix;
        TH3F *h3MixNext;
        TH3F *h3MixNext1;
        TH3F *h3MixNext2;
        if (izt <= 3)
        {
          h3Sam = (TH3F *)fileData->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt], assocZtThinner[izt + 1]));
          if (!systNMix)
            h3Mix = (TH3F *)fileData->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt], assocZtThinner[izt + 1]));
          else if (iSh==0 && systNMix)
            h3Mix = (TH3F *)fileDataMix->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt], assocZtThinner[izt + 1]));
          else if (iSh==1 && systNMix)
            h3Mix = (TH3F *)fileDataMix->Get(Form(sHistName + sIso + sShShNCentMix + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt], assocZtThinner[izt + 1]));
        }
        else if (izt == 4)
        {
          h3Sam = (TH3F *)fileData->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt], assocZtThinner[izt + 1]));
          h3SamNext = (TH3F *)fileData->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 1], assocZtThinner[izt + 2]));
          h3Sam->Add(h3SamNext);
          if (!systNMix)
          {
            h3Mix = (TH3F *)fileData->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt], assocZtThinner[izt + 1]));
            h3MixNext = (TH3F *)fileData->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 1], assocZtThinner[izt + 2]));
          }
          else if (iSh==0 && systNMix)
          {
            h3Mix = (TH3F *)fileDataMix->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt], assocZtThinner[izt + 1]));
            h3MixNext = (TH3F *)fileDataMix->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 1], assocZtThinner[izt + 2]));
          }
          else if (iSh==1 && systNMix)
          {
            h3Mix = (TH3F *)fileDataMix->Get(Form(sHistName + sIso + sShShNCentMix + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt], assocZtThinner[izt + 1]));
            h3MixNext = (TH3F *)fileDataMix->Get(Form(sHistName + sIso + sShShNCentMix + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 1], assocZtThinner[izt + 2]));
          }
          h3Mix->Add(h3MixNext);
        }
        else if (izt == 5)
        {
          h3Sam = (TH3F *)fileData->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 1], assocZtThinner[izt + 2]));
          h3SamNext = (TH3F *)fileData->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 2], assocZtThinner[izt + 3]));
          h3Sam->Add(h3SamNext);
          h3SamNext1 = (TH3F *)fileData->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 3], assocZtThinner[izt + 4]));
          h3Sam->Add(h3SamNext1);
          h3SamNext2 = (TH3F *)fileData->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 4], assocZtThinner[izt + 5]));
          h3Sam->Add(h3SamNext2);
          if (!systNMix)
          {
            h3Mix = (TH3F *)fileData->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 1], assocZtThinner[izt + 2]));
            h3MixNext = (TH3F *)fileData->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 2], assocZtThinner[izt + 3]));
            h3Mix->Add(h3MixNext);
            h3MixNext1 = (TH3F *)fileData->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 3], assocZtThinner[izt + 4]));
            h3Mix->Add(h3MixNext1);
            h3MixNext2 = (TH3F *)fileData->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 4], assocZtThinner[izt + 5]));
            h3Mix->Add(h3MixNext2);
          }
          else if (iSh==0 && systNMix)
          {
            h3Mix = (TH3F *)fileDataMix->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 1], assocZtThinner[izt + 2]));
            h3MixNext = (TH3F *)fileDataMix->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 2], assocZtThinner[izt + 3]));
            h3Mix->Add(h3MixNext);
            h3MixNext1 = (TH3F *)fileDataMix->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 3], assocZtThinner[izt + 4]));
            h3Mix->Add(h3MixNext1);
            h3MixNext2 = (TH3F *)fileDataMix->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 4], assocZtThinner[izt + 5]));
            h3Mix->Add(h3MixNext2);
          }
          else if (iSh==1 && systNMix)
          {
            h3Mix = (TH3F *)fileDataMix->Get(Form(sHistName + sIso + sShShNCentMix + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 1], assocZtThinner[izt + 2]));
            h3MixNext = (TH3F *)fileDataMix->Get(Form(sHistName + sIso + sShShNCentMix + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 2], assocZtThinner[izt + 3]));
            h3Mix->Add(h3MixNext);
            h3MixNext1 = (TH3F *)fileDataMix->Get(Form(sHistName + sIso + sShShNCentMix + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 3], assocZtThinner[izt + 4]));
            h3Mix->Add(h3MixNext1);
            h3MixNext2 = (TH3F *)fileDataMix->Get(Form(sHistName + sIso + sShShNCentMix + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 4], assocZtThinner[izt + 5]));
            h3Mix->Add(h3MixNext2);
          }
        }

        cout << "Define the pT trig range on the 3D distributions, projecting and getting 2D and 1D distributions" << endl;
        TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);
        for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
        {
          TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
          h3Sam->SetAxisRange(ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1] - 0.0001, "X");
          h3Mix->SetAxisRange(ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1] - 0.0001, "X");

          h2DdPhidEtaSam[iso][iSh][izt][iptTr] = (TH2F *)h3Sam->Project3D("zy"); // 2D correlation distributions deltaPhi and deltaEta same
          h2DdPhidEtaSam[iso][iSh][izt][iptTr]->SetName("h2DdPhidEtaSam" + sIso + sShSh + sZtBin + sPtTrig);
          h2DdPhidEtaMix[iso][iSh][izt][iptTr] = (TH2F *)h3Mix->Project3D("zy"); // 2D correlation distributions deltaPhi and deltaEta mix
          h2DdPhidEtaMix[iso][iSh][izt][iptTr]->SetName("h2DdPhidEtaMix" + sIso + sShSh + sZtBin + sPtTrig);

          hdPhiSam[iso][iSh][izt][iptTr] = (TH1F *)h2DdPhidEtaSam[iso][iSh][izt][iptTr]->ProjectionX("hdPhiSame" + sIso + sShSh + sZtBin + sPtTrig); // 1D correlation distributions deltaPhi and deltaEta same
          hdPhiMix[iso][iSh][izt][iptTr] = (TH1F *)h2DdPhidEtaMix[iso][iSh][izt][iptTr]->ProjectionX("hdPhiMix" + sIso + sShSh + sZtBin + sPtTrig);  // 1D correlation distributions deltaPhi and deltaEta mix

          // Normalization with the # of triggers: scale of the distributions for the pT trig integral
          cout << "Normalization with the # of triggers 1D distributions" << endl;
          hdPhiSam[iso][iSh][izt][iptTr]->Scale(1. / (hTriggerSam[iso][iSh]->Integral(hTriggerSam[iso][iSh]->FindBin(ptTrig[index1 + iptTr]), hTriggerSam[iso][iSh]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001))));
          hdPhiMix[iso][iSh][izt][iptTr]->Scale(1. / (hTriggerMix[iso][iSh]->Integral(hTriggerMix[iso][iSh]->FindBin(ptTrig[index1 + iptTr]), hTriggerMix[iso][iSh]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001))));
          cout << "Normalization with the # of triggers 2D distributions" << endl;
          h2DdPhidEtaSam[iso][iSh][izt][iptTr]->Scale(1. / hTriggerSam[iso][iSh]->Integral(hTriggerSam[iso][iSh]->FindBin(ptTrig[index1 + iptTr]), hTriggerSam[iso][iSh]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001)));
          h2DdPhidEtaMix[iso][iSh][izt][iptTr]->Scale(1. / hTriggerMix[iso][iSh]->Integral(hTriggerMix[iso][iSh]->FindBin(ptTrig[index1 + iptTr]), hTriggerMix[iso][iSh]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001)));

          // Definition of Same and Mix histograms for the mirroring
          cout << "Definition of Same and Mix histograms for the mirroring" << endl;
          hdPhiSamMirror[iso][iSh][izt][iptTr] = new TH1F("hdPhiSameMirror" + sIso + sShSh + sZtBin + sPtTrig, "hdPhiSameMirror" + sIso + sShSh + sZtBin + sPtTrig, hdPhiSam[iso][iSh][izt][iptTr]->GetNbinsX() / 2, 0, TMath::Pi());
          hdPhiMixMirror[iso][iSh][izt][iptTr] = new TH1F("hdPhiMixMirror" + sIso + sShSh + sZtBin + sPtTrig, "hdPhiMixMirror" + sIso + sShSh + sZtBin + sPtTrig, hdPhiMix[iso][iSh][izt][iptTr]->GetNbinsX() / 2, 0, TMath::Pi());
          if (bMirror) // histogram with mirroring
          {
            cout << "Mirroring" << endl;
            Mirroring(hdPhiSam[iso][iSh][izt][iptTr], hdPhiSamMirror[iso][iSh][izt][iptTr]);
            Mirroring(hdPhiMix[iso][iSh][izt][iptTr], hdPhiMixMirror[iso][iSh][izt][iptTr]);
            hdPhiSamMirror[iso][iSh][izt][iptTr]->Rebin(5);
            hdPhiMixMirror[iso][iSh][izt][iptTr]->Rebin(5);
            cout << "UE subtraction" << endl;
            hdPhiSamNoUE[iso][iSh][izt][iptTr] = (TH1F *)hdPhiSamMirror[iso][iSh][izt][iptTr]->Clone("hdPhiSameNoUE" + sIso + sShSh + sZtBin + sPtTrig);
            hdPhiSamNoUE[iso][iSh][izt][iptTr]->Sumw2();
            if (bZYAM)
            {
              cout << "Mirror and UE subtraction with ZYAM" << endl;
              fZYAM(hdPhiSamNoUE[iso][iSh][izt][iptTr]);
            }
            else
            {
              cout << "Mirror and UE subtraction with MIXED EVENT" << endl;
              hdPhiSamNoUE[iso][iSh][izt][iptTr]->Add(hdPhiMixMirror[iso][iSh][izt][iptTr], -1);
              cout << "Check errors before and after subtraction" << endl;
              for (int ibin = 1; ibin <= hdPhiSamNoUE[iso][iSh][izt][iptTr]->GetNbinsX(); ibin++)
              {
                cout << "Err before  Mix: " << hdPhiMixMirror[iso][iSh][izt][iptTr]->GetBinError(ibin) << " , Same: " << hdPhiSamMirror[iso][iSh][izt][iptTr]->GetBinError(ibin) << endl;
                cout << "Sum: " << sqrt(hdPhiMixMirror[iso][iSh][izt][iptTr]->GetBinError(ibin) * hdPhiMixMirror[iso][iSh][izt][iptTr]->GetBinError(ibin) + hdPhiSamMirror[iso][iSh][izt][iptTr]->GetBinError(ibin) * hdPhiSamMirror[iso][iSh][izt][iptTr]->GetBinError(ibin)) << endl;
                cout << "Err after subtraction: " << hdPhiSamNoUE[iso][iSh][izt][iptTr]->GetBinError(ibin) << endl;
              }

              // fZYAM(hdPhiSamNoUE[iso][iSh][izt][iptTr]); //Apply ZYAM after Mixed event
            }
          }
          else // histogram without mirroring
          {
            hdPhiMix[iso][iSh][izt][iptTr]->Rebin(5);
            hdPhiSam[iso][iSh][izt][iptTr]->Rebin(5);
            cout << "UE Subtraction" << endl;
            hdPhiSamNoUE[iso][iSh][izt][iptTr] = (TH1F *)hdPhiSam[iso][iSh][izt][iptTr]->Clone("hdPhiSameNoUE" + sIso + sShSh + sZtBin + sPtTrig);
            hdPhiSamNoUE[iso][iSh][izt][iptTr]->Sumw2();
            if (bZYAM)
            {
              cout << "ZYAM" << endl;
              fZYAM(hdPhiSamNoUE[iso][iSh][izt][iptTr]);
            }
            else
            {
              cout << "MIXED EVENT" << endl;
              hdPhiSamNoUE[iso][iSh][izt][iptTr]->Add(hdPhiMix[iso][iSh][izt][iptTr], -1);
            }
          }
          cout << "Save all the histogram on output file" << endl;
          fOutPut->cd();
          hdPhiSam[iso][iSh][izt][iptTr]->Write();
          hdPhiMix[iso][iSh][izt][iptTr]->Write();
          hdPhiMixMirror[iso][iSh][izt][iptTr]->Write();
          hdPhiSamMirror[iso][iSh][izt][iptTr]->Write();
          hdPhiSamNoUE[iso][iSh][izt][iptTr]->Write();
          h2DdPhidEtaSam[iso][iSh][izt][iptTr]->Write();
          h2DdPhidEtaMix[iso][iSh][izt][iptTr]->Write();
        }
      }
    }
  }
  ////////////////////////////////////////////////////////////////
  /////////////////// Purity correction///////////////////////////
  ////////////////////////////////////////////////////////////////
  cout << "Purity correction for photons" << endl;
  for (int iso = 0; iso < nIso; iso++)
  {
    TString sIso = Form("Iso%d", iso);
    for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
    {
      TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
      hZtPhotonPtBin[iso][iptTr] = new TH1F(Form("hZt%sPhotonPtBin_%s", sIso.Data(), sPtTrig.Data()), Form("hZt%sPhotonPtBin_%s", sIso.Data(), sPtTrig.Data()), nZtBin, assocZt);
      hZtPi0PtBin[iso][iptTr] = new TH1F(Form("hZt%sPi0PtBin_%s", sIso.Data(), sPtTrig.Data()), Form("hZt%sPi0PtBin_%s", sIso.Data(), sPtTrig.Data()), nZtBin, assocZt);

      for (Int_t izt = 0; izt < nZtBin; izt++)
      {
        double systValPur = histPur[iCen]->GetBinError(histPur[iCen]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001)) / (histPur[iCen]->GetBinContent(histPur[iCen]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001)));
        cout << "Pt: " << histPur[iCen]->GetBinCenter(histPur[iCen]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001)) << endl;
        double valPur = funcPur[iCen]->Eval(histPur[iCen]->GetBinCenter(histPur[iCen]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001)));
        double statValPur = histPurStat[iCen]->GetBinError(histPur[iCen]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001));
        if (systPur == 1)
          cout << "Standard Purity" << endl;

        if (systPur == 1.1)
        {
          cout << "Upper Limit Systematic Purity " << endl;
          valPur = valPur * (1 + systValPur);
        }
        if (systPur == 0.9)
        {
          cout << "Lower Limit Systematic Purity " << endl;
          valPur = valPur * (1 - systValPur);
        }

        cout << "Purity: " << valPur << endl;
        TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);

        //////////////////////////////////////////////////////////////////////
        ////// Definition of histograms at the purity correction step ///////
        ////////////////////////////////////////////////////////////////////
        hdPhiSamPur[iso][izt][iptTr] = (TH1F *)hdPhiSamNoUE[iso][0][izt][iptTr]->Clone(Form("hdPhiSamePur%s%s_%s", sIso.Data(), sZtBin.Data(), sPtTrig.Data()));
        hdPhiSamPi0Pur[iso][izt][iptTr] = (TH1F *)hdPhiSamNoUE[iso][1][izt][iptTr]->Clone(Form("hdPhiSamePi0Pur%s%s_%s", sIso.Data(), sZtBin.Data(), sPtTrig.Data())); // Distribution of Iso wide multiplied by 1-P
        hdPhiPi0[iso][izt][iptTr] = (TH1F *)hdPhiSamNoUE[iso][1][izt][iptTr]->Clone(Form("hdPhi%sPi0%s_%s", sIso.Data(), sZtBin.Data(), sPtTrig.Data()));              // Distribution of iso wide

        hdPhiSamNoUERatio[iso][izt][iptTr] = (TH1F *)hdPhiSamNoUE[iso][0][izt][iptTr]->Clone(Form("hdPhiSameNoUERatio%s%s_%s", sIso.Data(), sZtBin.Data(), sPtTrig.Data()));
        hdPhiSamNoUERatio[iso][izt][iptTr]->Divide(hdPhiSamNoUE[iso][1][izt][iptTr]);
        /////////////////////////////////////////////////////////////////////////
        ////// Purity correction procedure: [IsoNarrow - (1-P)IsoWide]/P ///////
        ///////////////////////////////////////////////////////////////////////
        cout << "Purity correction: [IsoNarrow - (1-Purity)IsoWide] / Purity" << endl;
        cout << "Multiply IsoWide for (1-Purity)" << endl;
        hdPhiSamPi0Pur[iso][izt][iptTr]->Scale((1 - valPur));
        cout << "IsoNarrow - (1-Purity)IsoWide" << endl;
        hdPhiSamPur[iso][izt][iptTr]->Sumw2();
        hdPhiSamPur[iso][izt][iptTr]->Add(hdPhiSamPi0Pur[iso][izt][iptTr], -1);
        cout << "Divide by Purity" << endl;
        hdPhiSamPur[iso][izt][iptTr]->Scale(1. / valPur);
        cout << "Error propagation for purity correction" << endl;
        for (int ibin = 1; ibin <= hdPhiSamPur[iso][izt][iptTr]->GetNbinsX(); ibin++)
        {
          cout << "Error before propagation: " << hdPhiSamPur[iso][izt][iptTr]->GetBinError(ibin) << endl;

          double errBkg = hdPhiSamNoUE[iso][1][izt][iptTr]->GetBinError(ibin);
          double BkgBin = hdPhiSamNoUE[iso][1][izt][iptTr]->GetBinContent(ibin);
          double errSig = hdPhiSamNoUE[iso][0][izt][iptTr]->GetBinError(ibin);
          double SigBin = hdPhiSamNoUE[iso][0][izt][iptTr]->GetBinContent(ibin);
          // cout<<"bin gamma "<<SigBin<<" binPi0 "<<BkgBin<< " purity "<< nPur<<endl;
          // cout<<"Bin "<<ibin<<" Err_gamma = "<<errSig<<", "<<"Err_pi0 = "<<errBkg<<"Err_Pur"<<nPurErr<<endl;
          double PropErr1 = TMath::Sqrt((errSig * errSig + errBkg * errBkg) / (valPur * valPur));
          double PropErr2 = TMath::Sqrt((statValPur * statValPur) * (SigBin * SigBin + BkgBin * BkgBin) / (valPur * valPur * valPur * valPur));
          double PropErr3 = TMath::Sqrt(errBkg * errBkg);
          // cout<<PropErr1<<" "<<PropErr2<<" "<<PropErr3<<endl;
          double PropErrAll = PropErr1 + PropErr2 + PropErr3;
          cout << "Error after propagation: " << PropErrAll << endl;
          hdPhiSamPur[iso][izt][iptTr]->SetBinError(ibin, PropErrAll);
        }

        hdPhiPhoton[iso][izt][iptTr] = (TH1F *)hdPhiSamPur[iso][izt][iptTr]->Clone(Form("hdPhi%sPhoton%s_%s", sIso.Data(), sZtBin.Data(), sPtTrig.Data())); // final result after purity correction
        // save in output root file all the results after purity correction
        fOutPut->cd();
        hdPhiSamPur[iso][izt][iptTr]->Write();
        hdPhiSamPi0Pur[iso][izt][iptTr]->Write();
        hdPhiPhoton[iso][izt][iptTr]->Write();
        hdPhiPi0[iso][izt][iptTr]->Write();
        hdPhiSamNoUERatio[iso][izt][iptTr]->Write();

        hdPhiPhotonPtBin[iso][izt][iptTr] = (TH1F *)hdPhiPhoton[iso][izt][iptTr]->Clone(Form("hdPhi%sPhotonPtBin%s_%s", sIso.Data(), sZtBin.Data(), sPtTrig.Data()));
        hdPhiPi0PtBin[iso][izt][iptTr] = (TH1F *)hdPhiPi0[iso][izt][iptTr]->Clone(Form("hdPhi%sPi0PtBin%s_%s", sIso.Data(), sZtBin.Data(), sPtTrig.Data()));

        cout << "Produce Zt function for every Pt bin" << endl;
        ZtFunction(hdPhiPhoton[iso][izt][iptTr], hZtPhotonPtBin[iso][iptTr], izt, phiMin, phiMax);
        ZtFunction(hdPhiPi0[iso][izt][iptTr], hZtPi0PtBin[iso][iptTr], izt, phiMin, phiMax);
      }
      cout << "Scale Zt distributions for the bin size" << endl;
      ScaleBinBySize(hZtPhotonPtBin[iso][iptTr]);
      ScaleBinBySize(hZtPi0PtBin[iso][iptTr]);

      fOutPut->cd();
      hZtPhotonPtBin[iso][iptTr]->Write();
      hZtPi0PtBin[iso][iptTr]->Write();
    }
    cout << "Azimuthal distribution on the large Pt Range " << endl;
    for (Int_t izt = 0; izt < nZtBin; izt++)
    {
      TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);

      for (int iSh = 0; iSh < nShSh; iSh++)
      {
        TString sShSh = Form("_ShSh%s", shshString[iSh].Data());
        hdPhiSamNoUEPtAll[iso][iSh][izt] = new TH1F("hdPhiSamNoUEPtAll" + sIso + sShSh + sCent + sZtBin + sPtAll, "hdPhiSamNoUEPtAll" + sIso + sShSh + sCent + sZtBin + sPtAll, hdPhiSamNoUE[iso][iSh][0][0]->GetNbinsX(), hdPhiSamNoUE[iso][iSh][0][0]->GetXaxis()->GetBinLowEdge(1), hdPhiSamNoUE[iso][iSh][0][0]->GetXaxis()->GetBinUpEdge(hdPhiSamNoUE[iso][iSh][0][0]->GetNbinsX()));
        hdPhiSamNoUEPtAll[iso][iSh][izt] = SumPtBinXzt(hTriggerSam[iso][iSh], ptTrig, index1, index2, hdPhiSamNoUE[iso][iSh][izt], hdPhiSamNoUEPtAll[iso][iSh][izt], histPur[iCen], funcPur[iCen], systPur, false);
      }

      hdPhiPhotonPtAll[iso][izt] = new TH1F(Form("hdPhi%sPhotonPtAll%s%s%s", sIso.Data(), sCent.Data(), sZtBin.Data(), sPtAll.Data()), Form("hdPhi%sPhotonPtAll%s%s%s", sIso.Data(), sCent.Data(), sZtBin.Data(), sPtAll.Data()), hdPhiPhotonPtBin[iso][0][0]->GetNbinsX(), hdPhiPhotonPtBin[iso][0][0]->GetXaxis()->GetBinLowEdge(1), hdPhiPhotonPtBin[iso][0][0]->GetXaxis()->GetBinUpEdge(hdPhiPhotonPtBin[iso][0][0]->GetNbinsX()));
      hdPhiPhotonPtAll[iso][izt] = SumPtBinXzt(hTriggerSam[iso][0], ptTrig, index1, index2, hdPhiPhotonPtBin[iso][izt], hdPhiPhotonPtAll[iso][izt], histPur[iCen], funcPur[iCen], systPur, true);

      hdPhiPi0PtAll[iso][izt] = new TH1F(Form("hdPhi%sPi0PtAll%s%s%s", sIso.Data(), sCent.Data(), sZtBin.Data(), sPtAll.Data()), Form("hdPhi%sPi0PtAll%s%s%s", sIso.Data(), sCent.Data(), sZtBin.Data(), sPtAll.Data()), hdPhiPi0PtBin[iso][0][0]->GetNbinsX(), hdPhiPi0PtBin[iso][0][0]->GetXaxis()->GetBinLowEdge(1), hdPhiPi0PtBin[iso][0][0]->GetXaxis()->GetBinUpEdge(hdPhiPi0PtBin[iso][0][0]->GetNbinsX()));
      hdPhiPi0PtAll[iso][izt] = SumPtBinXzt(hTriggerSam[iso][1], ptTrig, index1, index2, hdPhiPi0PtBin[iso][izt], hdPhiPi0PtAll[iso][izt], histPur[iCen], funcPur[iCen], systPur, false);
      fOutPut->cd();
      hdPhiPhotonPtAll[iso][izt]->Write();
      hdPhiPi0PtAll[iso][izt]->Write();
      for (int iSh = 0; iSh < nShSh; iSh++)
      {
        hdPhiSamNoUEPtAll[iso][iSh][izt]->Write();
      }
    }

    // Definition zT distribution for Photons and Pi0 on the large Pt
    cout << "Produce Zt function on the full Pt range" << endl;
    hZtPhoton[iso] = new TH1F(Form("hZt%sPhoton%s%s", sIso.Data(), sCent.Data(), sPtAll.Data()), Form("hZt%sPhoton%s%s", sIso.Data(), sCent.Data(), sPtAll.Data()), nZtBin, assocZt);
    hZtPhoton[iso] = SumPtBinXzt(hTriggerSam[iso][0], ptTrig, index1, index2, hZtPhotonPtBin[iso], hZtPhoton[iso], histPur[iCen], funcPur[iCen], systPur, true);

    // Definition zT distribution for Iso and Not Iso pi0
    hZtPi0[iso] = new TH1F(Form("hZt%sPi0%s%s", sIso.Data(), sCent.Data(), sPtAll.Data()), Form("hZt%sPi0%s%s", sIso.Data(), sCent.Data(), sPtAll.Data()), nZtBin, assocZt);
    hZtPi0[iso] = SumPtBinXzt(hTriggerSam[iso][1], ptTrig, index1, index2, hZtPi0PtBin[iso], hZtPi0[iso], histPur[iCen], funcPur[iCen], systPur, false);
    fOutPut->cd();
    hZtPhoton[iso]->Write();
    hZtPi0[iso]->Write();
  }
  // new TCanvas();
  // hZtPhoton->Draw();

  if (bPlot)
  {

    for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
    {
      TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
      TCanvas *cSame = new TCanvas(Form("cSame_%s" + sCent, sPtTrig.Data()), Form("cSame_%s" + sCent, sPtTrig.Data()), 5 * 800, 2 * 600);
      cSame->Divide(5, 2);
      TCanvas *cSamePi0 = new TCanvas(Form("cSamePi0_%s" + sCent, sPtTrig.Data()), Form("cSamePi0_%s" + sCent, sPtTrig.Data()), 5 * 800, 2 * 600);
      cSamePi0->Divide(5, 2);
      TCanvas *cMix = new TCanvas(Form("cMix_%s" + sCent, sPtTrig.Data()), Form("cMix_%s" + sCent, sPtTrig.Data()), 5 * 800, 2 * 600);
      cMix->Divide(5, 2);
      for (Int_t izt = 0; izt < nZtBin; izt++)
      {
        PlotStyle(hdPhiSam[1][0][izt][iptTr], 20, kBlack, "#Delta#varphi (rad)", "1/N^{trig}dN/d#varphi");
        PlotStyle(hdPhiSam[1][1][izt][iptTr], 24, kBlack, "#Delta#varphi (rad)", "1/N^{trig}dN/d#varphi");
        PlotStyle(hdPhiMix[1][0][izt][iptTr], 20, kRed, "#Delta#varphi (rad)", "1/N^{trig}dN/d#varphi");
        PlotStyle(hdPhiMix[1][1][izt][iptTr], 24, kRed, "#Delta#varphi (rad)", "1/N^{trig}dN/d#varphi");

        cSame->cd(izt + 1);
        hdPhiSam[1][0][izt][iptTr]->SetMaximum(0.25);
        hdPhiSam[1][0][izt][iptTr]->Draw("same");
        // hdPhiMix[1][0][izt][iptTr]->Draw("same");

        cSamePi0->cd(izt + 1);
        hdPhiSam[1][1][izt][iptTr]->SetMaximum(0.25);
        hdPhiSam[1][1][izt][iptTr]->Draw("same");
        // hdPhiMix[1][1][izt][iptTr]->Draw("same");

        cMix->cd(izt + 1);
        hdPhiMix[1][0][izt][iptTr]->Draw("same");
        hdPhiMix[1][1][izt][iptTr]->Draw("same");
      }

      TCanvas *cSamePur = new TCanvas(Form("cSamePur_%s" + sCent, sPtTrig.Data()), Form("cSamePur_%s" + sCent, sPtTrig.Data()), 5 * 800, 2 * 600);
      cSamePur->Divide(5, 2);
      for (Int_t izt = 0; izt < nZtBin; izt++)
      {
        PlotStyle(hdPhiPhoton[1][izt][iptTr], 25, kOrange + 7, "#Delta#varphi (rad)", "1/N^{trig}dN/d#varphi");
        PlotStyle(hdPhiSamNoUE[1][0][izt][iptTr], 20, kAzure + 2, "#Delta#varphi (rad)", "1/N^{trig}dN/d#varphi");
        PlotStyle(hdPhiSamPur[1][izt][iptTr], 20, kGreen, "#Delta#varphi (rad)", "1/N^{trig}dN/d#varphi");
        PlotStyle(hdPhiSamPi0Pur[1][izt][iptTr], 20, kPink + 2, "#Delta#varphi (rad)", "1/N^{trig}dN/d#varphi");
        cSamePur->cd(izt + 1);
        hdPhiSamNoUE[1][0][izt][iptTr]->SetMaximum(0.25);
        // hdPhiSamNoUE[1][0][izt][iptTr]->SetMinimum(0.9*hdPhiSamPi0Pur[izt][iptTr]->GetMinimum());
        hdPhiSamNoUE[1][0][izt][iptTr]->Draw("same");
        // hdPhiSamPur[izt][iptTr]->Draw("same");
        hdPhiSamPi0Pur[1][izt][iptTr]->Draw("same");
        hdPhiPhoton[1][izt][iptTr]->Draw("same");
      }
      TCanvas *cSameMixPur = new TCanvas(Form("cSameMixPur_%s" + sCent, sPtTrig.Data()), Form("cSameMixPur_%s" + sCent, sPtTrig.Data()), 5 * 800, 2 * 600);
      cSameMixPur->Divide(5, 2);
      TCanvas *cIsoGamma = new TCanvas(Form("cIsoGamma_%s" + sCent, sPtTrig.Data()), Form("cIsoGamma_%s" + sCent, sPtTrig.Data()), 5 * 800, 2 * 600);
      cIsoGamma->Divide(5, 2);
      for (Int_t izt = 0; izt < nZtBin; izt++)
      {
        TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);
        cSameMixPur->cd(izt + 1);
        hdPhiSamPur[1][izt][iptTr]->Draw("same");
        PlotStyle(hdPhiPhoton[1][izt][iptTr], 25, kOrange + 7, "#Delta#varphi (rad)", "1/N^{trig}dN/d#varphi");
        cIsoGamma->cd(izt + 1);
        hdPhiPhoton[1][izt][iptTr]->Draw();
        // Mirroring(hdPhiPhoton[izt]);
      }
    }
  }

  bool bPlotMC = true;
  if (bPlotMC)
  {
    // MC_GJ and JJlow
    TH1F *hTriggerMCGen[nIso][nShSh];
    TH1F *hTriggerSamMCRec[nIso][nShSh];
    TH1F *hTriggerMixMCRec[nIso][nShSh];
    TH2F *h2DdPhidEtaMCGen[nIso][nShSh][nZtBin][nPtTrig];
    TH2F *h2DdPhidEtaSamMCRec[nIso][nShSh][nZtBin][nPtTrig];
    TH2F *h2DdPhidEtaMixMCRec[nIso][nShSh][nZtBin][nPtTrig];
    TH1F *hdPhiMCGen[nIso][nShSh][nZtBin][nPtTrig];         // Generated level
    TH1F *hdPhiSamMCRec[nIso][nShSh][nZtBin][nPtTrig];      // Same Reconstructed with UE
    TH1F *hdPhiMixMCRec[nIso][nShSh][nZtBin][nPtTrig];      // Mixed Reconstructed with UE
    TH1F *hdPhiMCGenMirrorUE[nIso][nShSh][nZtBin][nPtTrig]; // Same Generated with UE
    TH1F *hdPhiMCGenMirror[nIso][nShSh][nZtBin][nPtTrig];   // Same Generated with UE
    TH1F *hdPhiSamMCRecMirror[nIso][nShSh][nZtBin][nPtTrig];
    TH1F *hdPhiMixMCRecMirror[nIso][nShSh][nZtBin][nPtTrig];
    TH1F *hdPhiSamMCRecNoUE[nIso][nShSh][nZtBin][nPtTrig]; // Same - Mixed Reconstructed wo UE
    TH1F *hdPhiPhotonMCGen[nIso][nShSh][nZtBin][nPtTrig];
    TH1F *hdPhiPhotonMCRec[nIso][nShSh][nZtBin][nPtTrig];
    TH1F *hZtPtBinMCGen[nIso][nShSh][nPtTrig];      // Zt distribution MC Generated X Ptbin
    TH1F *hZtPtBinMCRec[nIso][nShSh][nPtTrig];      // Zt distribution MC Recostructed X Ptbin
    TH1F *hRatioEffCorrPtBin[nIso][nShSh][nPtTrig]; // Efficiency correction Gen/Rec X Ptbin

    TH1F *hZtMCGen[nIso][nShSh];
    TH1F *hZtMCRec[nIso][nShSh];
    TH1F *hRatioEffCorr[nIso][nShSh];
    TH1F *hZtEffCorrPhoton[nIso];
    TH1F *hZtEffCorrPi0[nIso];
    TString sNamePtTrigGen[nShSh] = {"Photon", "Pi0"};
    TString sNamePtTrigRec[nShSh] = {"_MCPhoton", ""};

    //          h3Sam = (TH3F *)fileData->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt], assocZtThinner[izt + 1]));
    // h3Mix = (TH3F *)fileData->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt], assocZtThinner[izt + 1]));
    for (int iso = 0; iso < nIso; iso++)
    {
      TString sIso = Form("Iso%d", iso);

      for (int iSh = 0; iSh < 1; iSh++)
      {
        TCanvas *cIsoPtTrigMC = new TCanvas(Form("Iso%dPtTrigMC", iso), Form("Iso%dPtTrigMC", iso), 800, 600);
        TString sShSh = Form("_ShSh%s", shshStringMC[iSh].Data());
        cout << "pippo" << endl;
        cout << " MC Pt Trigger distributions for Same Generated and for Same and Mix Recostructed" << endl;
        cout << fileMC[iSh] << endl;
        cout << sHistName + sIso + sShSh + sCent << endl;
        hTriggerMCGen[iso][iSh] = (TH1F *)fileMC[iSh]->Get(sHistName + sIso + sShSh + sCent + Form("_hMCPtTrigger_%s", sNamePtTrigGen[iSh].Data()));
        cout << hTriggerMCGen[iso][iSh] << endl;
        hTriggerSamMCRec[iso][iSh] = (TH1F *)fileMC[iSh]->Get(sHistName + sIso + sShSh + sCent + Form("_hPtTrigger%s", sNamePtTrigRec[iSh].Data()));
        cout << hTriggerSamMCRec[iso][iSh] << endl;
        hTriggerMixMCRec[iso][iSh] = (TH1F *)fileMC[iSh]->Get(sHistName + sIso + sShSh + sCent + "_hPtTriggerMixed");
        hTriggerMixMCRec[iso][iSh]->SetName(sHistName + sIso + sShSh + sCent + "_hPtTriggerMixedMCRec"); // change name because otherwise same name Mix used in Data and double in the saving output file
        cout << hTriggerMixMCRec[iso][iSh] << endl;

        hTriggerMCGen[iso][iSh]->SetDirectory(0);
        hTriggerSamMCRec[iso][iSh]->SetDirectory(0);
        hTriggerMixMCRec[iso][iSh]->SetDirectory(0);
        hTriggerMixMCRec[iso][iSh]->SetMarkerColor(kRed);
        hTriggerMixMCRec[iso][iSh]->SetMarkerStyle(20);
        hTriggerSamMCRec[iso][iSh]->SetMarkerStyle(21);
        hTriggerMCGen[iso][iSh]->SetMarkerStyle(25);
        hTriggerMixMCRec[iso][iSh]->Draw("same");
        hTriggerSamMCRec[iso][iSh]->Draw("same");
        hTriggerMCGen[iso][iSh]->Draw("same");
        fOutPut->cd();
        hTriggerSamMCRec[iso][iSh]->Write();
        hTriggerMixMCRec[iso][iSh]->Write();
        for (Int_t izt = 0; izt < nZtBin; izt++)
        {
          TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);

          TH3F *h3MCGen;
          TH3F *h3MCGenNext;
          TH3F *h3SamMCRec;
          TH3F *h3SamMCRecNext;
          TH3F *h3MixMCRec;
          TH3F *h3MixMCRecNext;
          if (izt <= 3)
          {
            h3MCGen = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hMCDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_%s", assocZt[izt], assocZt[izt + 1], sNamePtTrigGen[iSh].Data()));
            h3SamMCRec = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt], assocZtThinner[izt + 1]));
            h3MixMCRec = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt], assocZtThinner[izt + 1]));
          }
          else if (izt == 4)
          {
            h3MCGen = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hMCDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_%s", assocZtThinner[izt], assocZtThinner[izt + 1], sNamePtTrigGen[iSh].Data()));
            h3MCGenNext = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hMCDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_%s", assocZtThinner[izt + 1], assocZtThinner[izt + 2], sNamePtTrigGen[iSh].Data()));
            h3MCGen->Add(h3MCGenNext);

            h3SamMCRec = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt], assocZtThinner[izt + 1]));
            h3SamMCRecNext = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 1], assocZtThinner[izt + 2]));
            h3SamMCRec->Add(h3SamMCRecNext);

            h3MixMCRec = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt], assocZtThinner[izt + 1]));
            h3MixMCRecNext = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 1], assocZtThinner[izt + 2]));
            h3MixMCRec->Add(h3MixMCRecNext);
          }
          else if (izt == 5)
          {
            h3MCGen = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hMCDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_%s", assocZtThinner[izt + 1], assocZtThinner[izt + 2], sNamePtTrigGen[iSh].Data()));
            h3MCGenNext = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hMCDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_%s", assocZtThinner[izt + 2], assocZtThinner[izt + 3], sNamePtTrigGen[iSh].Data()));
            h3MCGen->Add(h3MCGenNext);

            h3SamMCRec = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 1], assocZtThinner[izt + 2]));
            h3SamMCRecNext = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 2], assocZtThinner[izt + 3]));
            h3SamMCRec->Add(h3SamMCRecNext);

            h3MixMCRec = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 1], assocZtThinner[izt + 2]));
            h3MixMCRecNext = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 2], assocZtThinner[izt + 3]));
            h3MixMCRec->Add(h3MixMCRecNext);
          }
          else if (izt == 6)
          {
            h3MCGen = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hMCDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_%s", assocZtThinner[izt + 2], assocZtThinner[izt + 3], sNamePtTrigGen[iSh].Data()));
            h3MCGenNext = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hMCDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f_%s", assocZtThinner[izt + 3], assocZtThinner[izt + 4], sNamePtTrigGen[iSh].Data()));
            h3MCGen->Add(h3MCGenNext);

            h3SamMCRec = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 2], assocZtThinner[izt + 3]));
            h3SamMCRecNext = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 3], assocZtThinner[izt + 4]));
            h3SamMCRec->Add(h3SamMCRecNext);

            h3MixMCRec = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 2], assocZtThinner[izt + 3]));
            h3MixMCRecNext = (TH3F *)fileMC[iSh]->Get(Form(sHistName + sIso + sShSh + sCent + "_hMixDeltaPhiDeltaEtaChargedZTBin%1.2f_%1.2f", assocZtThinner[izt + 3], assocZtThinner[izt + 4]));
            h3MixMCRec->Add(h3MixMCRecNext);
          }
          for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
          {
            TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
            h3MCGen->SetAxisRange(ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1] - 0.0001, "X");
            h3SamMCRec->SetAxisRange(ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1] - 0.0001, "X");
            h3MixMCRec->SetAxisRange(ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1] - 0.0001, "X");

            h2DdPhidEtaMCGen[iso][iSh][izt][iptTr] = (TH2F *)h3MCGen->Project3D("zy");
            h2DdPhidEtaMCGen[iso][iSh][izt][iptTr]->SetName("h2DdPhidEtaMCGen" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh]);
            h2DdPhidEtaSamMCRec[iso][iSh][izt][iptTr] = (TH2F *)h3SamMCRec->Project3D("zy");
            h2DdPhidEtaSamMCRec[iso][iSh][izt][iptTr]->SetName("h2DdPhidEtaSamMCRec" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh]);
            h2DdPhidEtaMixMCRec[iso][iSh][izt][iptTr] = (TH2F *)h3MixMCRec->Project3D("zy");
            h2DdPhidEtaMixMCRec[iso][iSh][izt][iptTr]->SetName("h2DdPhidEtaMixMCRec" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh]);

            hdPhiMCGen[iso][iSh][izt][iptTr] = (TH1F *)h2DdPhidEtaMCGen[iso][iSh][izt][iptTr]->ProjectionX("hdPhiMCGen" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh]);
            hdPhiSamMCRec[iso][iSh][izt][iptTr] = (TH1F *)h2DdPhidEtaSamMCRec[iso][iSh][izt][iptTr]->ProjectionX("hdPhiSameMCRec" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh]);
            hdPhiMixMCRec[iso][iSh][izt][iptTr] = (TH1F *)h2DdPhidEtaMixMCRec[iso][iSh][izt][iptTr]->ProjectionX("hdPhiMixMCRec" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh]);
            hdPhiMixMCRec[iso][iSh][izt][iptTr]->SetMarkerColor(kRed);

            hdPhiMCGen[iso][iSh][izt][iptTr]->Scale(1. / hTriggerMCGen[iso][iSh]->Integral(hTriggerMCGen[iso][iSh]->FindBin(ptTrig[index1 + iptTr]), hTriggerMCGen[iso][iSh]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001)));
            hdPhiSamMCRec[iso][iSh][izt][iptTr]->Scale(1. / hTriggerSamMCRec[iso][iSh]->Integral(hTriggerSamMCRec[iso][iSh]->FindBin(ptTrig[index1 + iptTr]), hTriggerSamMCRec[iso][iSh]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001)));
            hdPhiMixMCRec[iso][iSh][izt][iptTr]->Scale(1. / hTriggerMixMCRec[iso][iSh]->Integral(hTriggerMixMCRec[iso][iSh]->FindBin(ptTrig[index1 + iptTr]), hTriggerMixMCRec[iso][iSh]->FindBin(ptTrig[index1 + iptTr + 1] - 0.0001)));

            cout << "Definition of MC histograms for the mirroring" << endl;
            int nbinsX = hdPhiMCGen[iso][iSh][izt][iptTr]->GetNbinsX() / 2;
            hdPhiMCGenMirror[iso][iSh][izt][iptTr] = new TH1F("hdPhiMCGenMirror" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh], "hdPhiMCGenMirror" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh], nbinsX, 0, TMath::Pi());
            hdPhiSamMCRecMirror[iso][iSh][izt][iptTr] = new TH1F("hdPhiSameMCRecMirror" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh], "hdPhiSameMCRecMirror" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh], nbinsX, 0, TMath::Pi());
            hdPhiMixMCRecMirror[iso][iSh][izt][iptTr] = new TH1F("hdPhiMixMCRecMirror" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh], "hdPhiMixMCRecMirror" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh], nbinsX, 0, TMath::Pi());

            if (bMirror)
            {
              Mirroring(hdPhiMCGen[iso][iSh][izt][iptTr], hdPhiMCGenMirror[iso][iSh][izt][iptTr]);
              Mirroring(hdPhiSamMCRec[iso][iSh][izt][iptTr], hdPhiSamMCRecMirror[iso][iSh][izt][iptTr]);
              Mirroring(hdPhiMixMCRec[iso][iSh][izt][iptTr], hdPhiMixMCRecMirror[iso][iSh][izt][iptTr]);

              hdPhiMCGenMirror[iso][iSh][izt][iptTr]->Rebin(5);
              hdPhiSamMCRecMirror[iso][iSh][izt][iptTr]->Rebin(5);
              hdPhiMixMCRecMirror[iso][iSh][izt][iptTr]->Rebin(5);

              cout << "MC UE subtraction" << endl;
              hdPhiMCGenMirrorUE[iso][iSh][izt][iptTr] = (TH1F *)hdPhiMCGenMirror[iso][iSh][izt][iptTr]->Clone("hdPhiMCGenMirrorUE" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh]);
              fZYAM(hdPhiMCGenMirror[iso][iSh][izt][iptTr]);

              hdPhiSamMCRecNoUE[iso][iSh][izt][iptTr] = (TH1F *)hdPhiSamMCRecMirror[iso][iSh][izt][iptTr]->Clone("hdPhiSameMCRecNoUE" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh]);
              hdPhiSamMCRecNoUE[iso][iSh][izt][iptTr]->Sumw2();
              // double scaleFactMC = fZYAM_Mix(hdPhiSamMCRecMirror[iso][iSh][izt][iptTr], hdPhiMixMCRecMirror[iso][iSh][izt][iptTr]);
              if (bZYAM)
              {
                cout << "Mirror - ZYAM" << endl;
                fZYAM(hdPhiSamMCRecNoUE[iso][iSh][izt][iptTr]);
              }
              else
              {
                cout << "Mirror - MIXED EVENT" << endl;
                cout << hdPhiSamMCRecNoUE[iso][iSh][izt][iptTr] << endl;
                hdPhiSamMCRecNoUE[iso][iSh][izt][iptTr]->Add(hdPhiMixMCRecMirror[iso][iSh][izt][iptTr], -1);
              }
              fOutPut->cd();
              hdPhiMCGenMirrorUE[iso][iSh][izt][iptTr]->Write();
              hdPhiMCGenMirror[iso][iSh][izt][iptTr]->Write();
              hdPhiSamMCRecMirror[iso][iSh][izt][iptTr]->Write();
              hdPhiMixMCRecMirror[iso][iSh][izt][iptTr]->Write();
              hdPhiSamMCRecNoUE[iso][iSh][izt][iptTr]->Write();
            }
            else
            {
              hdPhiMCGen[iso][iSh][izt][iptTr]->Rebin(5);
              hdPhiSamMCRec[iso][iSh][izt][iptTr]->Rebin(5);
              hdPhiMixMCRec[iso][iSh][izt][iptTr]->Rebin(5);

              cout << "MC UE subtraction" << endl;
              fZYAM(hdPhiMCGen[iso][iSh][izt][iptTr]);
              hdPhiSamMCRecNoUE[iso][iSh][izt][iptTr] = (TH1F *)hdPhiSamMCRec[iso][iSh][izt][iptTr]->Clone("hdPhiSameMCRecNoUE" + sIso + sShSh + sZtBin + sPtTrig + sNamePtTrigGen[iSh]);
              hdPhiSamMCRecNoUE[iso][iSh][izt][iptTr]->Sumw2();

              if (bZYAM)
              {
                cout << "ZYAM" << endl;
                fZYAM(hdPhiSamMCRecNoUE[iso][iSh][izt][iptTr]);
              }
              else
              {
                cout << "MIXED EVENT" << endl;
                hdPhiSamMCRecNoUE[iso][iSh][izt][iptTr]->Add(hdPhiMixMCRec[iso][iSh][izt][iptTr], -1);
              }
              fOutPut->cd();
              hdPhiMCGen[iso][iSh][izt][iptTr]->Write();
              hdPhiSamMCRec[iso][iSh][izt][iptTr]->Write();
              hdPhiMixMCRec[iso][iSh][izt][iptTr]->Write();
              hdPhiSamMCRecNoUE[iso][iSh][izt][iptTr]->Write();
            }
          }
        }
      }
    }
    cout << "MC Zt distributions for Pt bin" << endl;
    for (int iso = 0; iso < nIso; iso++)
    {
      TString sIso = Form("Iso%d", iso);
      for (int iSh = 0; iSh < 1; iSh++)
      {
        for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
        {
          cout << nPtTrig << endl;
          cout << index2 - index1 << endl;
          TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
          hZtPtBinMCGen[iso][iSh][iptTr] = new TH1F(Form("hZtPtBinMCGen_%s%s%s", sIso.Data(), sPtTrig.Data(), sNamePtTrigGen[iSh].Data()), Form("hZtPtBinMCGen_%s%s%s", sIso.Data(), sPtTrig.Data(), sNamePtTrigGen[iSh].Data()), nZtBin, assocZt);
          hZtPtBinMCRec[iso][iSh][iptTr] = new TH1F(Form("hZtPtBinMCRec_%s%s%s", sIso.Data(), sPtTrig.Data(), sNamePtTrigGen[iSh].Data()), Form("hZtPtBinMCRec_%s%s%s", sIso.Data(), sPtTrig.Data(), sNamePtTrigGen[iSh].Data()), nZtBin, assocZt);
          for (Int_t izt = 0; izt < nZtBin; izt++)
          {
            // IsoGammaMC
            if (!bMirror)
            {
              ZtFunction(hdPhiMCGen[iso][iSh][izt][iptTr], hZtPtBinMCGen[iso][iSh][iptTr], izt, phiMin, phiMax);
            }
            else if (bMirror)
            {
              ZtFunction(hdPhiMCGenMirror[iso][iSh][izt][iptTr], hZtPtBinMCGen[iso][iSh][iptTr], izt, phiMin, phiMax);
            }

            ZtFunction(hdPhiSamMCRecNoUE[iso][iSh][izt][iptTr], hZtPtBinMCRec[iso][iSh][iptTr], izt, phiMin, phiMax);
          }
          ScaleBinBySize(hZtPtBinMCGen[iso][iSh][iptTr]);
          ScaleBinBySize(hZtPtBinMCRec[iso][iSh][iptTr]);
          PlotStyle(hZtPtBinMCGen[iso][iSh][iptTr], 21, kBlue - 4, "zt", "1/N dN/dzt");
          PlotStyle(hZtPtBinMCRec[iso][iSh][iptTr], 21, kOrange + 7, "zt", "1/N dN/dzt");

          new TCanvas();
          hZtPtBinMCGen[iso][iSh][iptTr]->Draw("same");
          hZtPtBinMCRec[iso][iSh][iptTr]->Draw("same");

          hRatioEffCorrPtBin[iso][iSh][iptTr] = (TH1F *)hZtPtBinMCGen[iso][iSh][iptTr]->Clone("hRatioEffCorrPtBin" + sIso + sPtTrig + sNamePtTrigGen[iSh]);
          hRatioEffCorrPtBin[iso][iSh][iptTr]->Sumw2();
          hRatioEffCorrPtBin[iso][iSh][iptTr]->Divide(hZtPtBinMCRec[iso][iSh][iptTr]);
          fOutPut->cd();
          hZtPtBinMCGen[iso][iSh][iptTr]->Write();
          hZtPtBinMCRec[iso][iSh][iptTr]->Write();
          hRatioEffCorrPtBin[iso][iSh][iptTr]->Write();
        }
        hZtMCGen[iso][iSh] = new TH1F(Form("hZtMCGen%s%s%s%s", sIso.Data(), sNamePtTrigGen[iSh].Data(), sCent.Data(), sPtAll.Data()), Form("hZtMCGen%s%s%s%s", sIso.Data(), sNamePtTrigGen[iSh].Data(), sCent.Data(), sPtAll.Data()), nZtBin, assocZt);
        hZtMCGen[iso][iSh] = SumPtBinXzt(hTriggerMCGen[iso][iSh], ptTrig, index1, index2, hZtPtBinMCGen[iso][iSh], hZtMCGen[iso][iSh], histPur[iCen], funcPur[iCen], systPur, false);
        hZtMCRec[iso][iSh] = new TH1F(Form("hZtMCRec%s%s%s%s", sIso.Data(), sNamePtTrigGen[iSh].Data(), sCent.Data(), sPtAll.Data()), Form("hZtMCRec%s%s%s%s", sIso.Data(), sNamePtTrigGen[iSh].Data(), sCent.Data(), sPtAll.Data()), nZtBin, assocZt);
        hZtMCRec[iso][iSh] = SumPtBinXzt(hTriggerSamMCRec[iso][iSh], ptTrig, index1, index2, hZtPtBinMCRec[iso][iSh], hZtMCRec[iso][iSh], histPur[iCen], funcPur[iCen], systPur, false);
        hRatioEffCorr[iso][iSh] = (TH1F *)hZtMCGen[iso][iSh]->Clone("hRatioEffCorr" + sIso + sNamePtTrigGen[iSh]);
        hRatioEffCorr[iso][iSh]->Divide(hZtMCRec[iso][iSh]);

        fOutPut->cd();
        hZtMCGen[iso][iSh]->Write();
        hZtMCRec[iso][iSh]->Write();
        hRatioEffCorr[iso][iSh]->Write();
      }
      hZtEffCorrPhoton[iso] = (TH1F *)hZtPhoton[iso]->Clone(Form("hZtEffCorr%sPhoton%s%s", sIso.Data(), sCent.Data(), sPtAll.Data()));
      hZtEffCorrPhoton[iso]->Multiply(hRatioEffCorr[iso][0]);

      hZtEffCorrPi0[iso] = (TH1F *)hZtPi0[iso]->Clone(Form("hZtEffCorr%sPi0%s%s", sIso.Data(), sCent.Data(), sPtAll.Data()));
      // hZtEffCorrPi0[iso]->Multiply(hRatioEffCorr[iso][1]);
    }

    fOutPut->cd();
    for (int iso = 0; iso < nIso; iso++)
    {
      hZtEffCorrPhoton[iso]->Write();
      hZtEffCorrPi0[iso]->Write();
    }
    if (bPlotMC)
    {
      for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
      {
        TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
        TCanvas *cSameMC = new TCanvas(Form("cSameMC_%s", sPtTrig.Data()), Form("cSameMC_%s", sPtTrig.Data()), 5 * 800, 2 * 600);
        cSameMC->Divide(5, 2);
        TCanvas *cMixMC = new TCanvas(Form("cMixMC_IsoNotIso_%s", sPtTrig.Data()), Form("cMixMC_IsoNotIso_%s", sPtTrig.Data()), 5 * 800, 2 * 600);
        cMixMC->Divide(5, 2);
        for (Int_t izt = 0; izt < nZtBin; izt++)
        {
          PlotStyle(hdPhiMCGen[1][0][izt][iptTr], 20, kAzure + 3, "#Delta#varphi (rad)", "1/N^{trig}dN/d#varphi");
          PlotStyle(hdPhiSamMCRec[1][0][izt][iptTr], 20, kOrange + 7, "#Delta#varphi (rad)", "1/N^{trig}dN/d#varphi");
          PlotStyle(hdPhiSamMCRecNoUE[1][0][izt][iptTr], 25, kPink + 3, "#Delta#varphi (rad)", "1/N^{trig}dN/d#varphi");
          PlotStyle(hdPhiMixMCRec[1][0][izt][iptTr], 20, kGreen - 6, "#Delta#varphi (rad)", "1/N^{trig}dN/d#varphi");
          // PlotStyle(hdPhiMixMCRec[0][0][izt][iptTr], 25, kGreen - 6, "#Delta#varphi (rad)", "1/N^{trig}dN/d#varphi");

          cSameMC->cd(izt + 1);
          hdPhiSamMCRecNoUE[1][0][izt][iptTr]->Draw("same");
          hdPhiSamMCRec[1][0][izt][iptTr]->Draw("same");
          hdPhiMixMCRec[1][0][izt][iptTr]->Draw("same");
          // hdPhiMCGen[1][0][izt][iptTr]->Draw("same");

          cMixMC->cd(izt + 1);
          hdPhiMixMCRec[1][0][izt][iptTr]->Draw("same");
          // hdPhiMixMCRec[0][0][izt][iptTr]->Draw("same");
        }
      }
    }
    cout << "pippo" << endl;
    new TCanvas();
    cout << hZtMCGen[1][0] << endl; //[iso][shsh]
    cout << hZtPhoton[1] << endl;   //[iso]

    // hZtMCGen[1][1]->Draw("same"); //[iso][shsh]
    hZtPhoton[1]->Draw("same"); //[iso]
  }
}

void IsoGammaHadron(float ptTrMin = 18, float ptTrMax = 40, TString sFileDirShSig = "~/work/histogram/DataSh100_AssocPt500", Bool_t bMirror = true, TString shshBkg = "0.40-1.00", TString dirFiles = "~/work/histogram/FromScratch/checkCode", double systPur = 1, bool bZYAM = false, bool bPlot = true, double phiMin = TMath::Pi() * 3 / 5, double phiMax = TMath::Pi(), bool systShSh = false, bool systTrackIneff = false, bool systNMix18 = false, bool systNMix45 = false)
{
  ///////////////////////////////////////////////////////////////////
  /////// Define MC root files: one file for all centralities //////
  /////////////////////////////////////////////////////////////////
  if (!systTrackIneff)
    fileMC[0] = TFile::Open(Form("~/work/histogram/MC_GJShSh150/MC_GJSh150.root"));
  else if (systTrackIneff)
    fileMC[0] = TFile::Open(Form("~/work/histogram/MC_GJShSh150/MC_GJTrackInEff.root")); // MC root file for systematic on tracking efficiency

  fileMC[1] = TFile::Open(Form("~/work/histogram/MC_JJlow_Pi0/MC_JJlow_0_90.root"));

  // fileMC = TFile::Open(Form("~/work/histogram/MC_GJShSh150/MC_GJSh150.root")); //Systematic TrackEfficiency

  if (!fileMC[0] || !fileMC[1])
    cout << "MC File doesn't exist" << endl;
  ///////////////////////////////////////////////////////////////////
  ///////// Define data root files: one file per centrality ////////
  /////////////////////////////////////////////////////////////////
  TString tagFile[nCen];
  // comments for shshSyst and NCentBinMixSyst
  //  fileDataMix = TFile::Open(Form("~/work/histogram/NCentBinMix45/%s.root", tagFile.Data()));
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    tagFile[iCen] = Form("EMCAL_MB_%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    fileDataShStd = TFile::Open(Form("%s/%s.root", sFileDirShSig.Data(), tagFile[iCen].Data()));
    if (!fileDataShStd)
      cout << "Data File doesn't exist" << endl;

    if (systShSh)
    {
      cout << "ShSh Bkg Systematics: ON" << endl;
      if (shshBkg == "0.35-1.00")
      {
        fileDataShBkg = TFile::Open(Form("~/work/histogram/ShShSyst03510/%s.root", tagFile[iCen].Data())); // ShSh 0.35-1.00
      }
      else if (shshBkg == "0.40-2.00")
      {
        fileDataShBkg = TFile::Open(Form("~/work/histogram/ShShSystMultBkg/%s.root", tagFile[iCen].Data())); // ShSh 0.40-2.00
      }
      else if (shshBkg == "0.40-1.50")
      {
        fileDataShBkg = TFile::Open(Form("~/work/histogram/ShShSystMultBkg/%s.root", tagFile[iCen].Data())); // ShSh 0.40-1.50
      }
      cout << "ShSh bkg: " << shshBkg << endl;
    }

    if (systPur != 1)
      cout << "Purity Systematics: ON, purity value" << systPur << endl;

    if (systNMix18)
    {
      //fileDataMix = TFile::Open(Form("~/work/histogram/RootFiles/SystematicsNCentrBin/%s.root", tagFile[iCen].Data())); // estimated with shsh = 0.40-1.00
      // sShShNCentMix = "_ShSh0.40-1.00"
      fileDataMix = TFile::Open(Form("~/work/histogram/IsoPhotonHadronCorrelations/RootFiles/NCentBinMix18/EMCAL_MB_0_90.root")); // Old Files with shsh = 0.40-2.00
      sShShNCentMix = "_ShSh0.40-2.00";
      systNMix = true;
    }
    else if(systNMix45)
    {
      //fileDataMix = TFile::Open(Form("~/work/histogram/RootFiles/SystematicsNCentrBin/%s.root", tagFile[iCen].Data())); // estimated with shsh = 0.40-1.00
      fileDataMix = TFile::Open(Form("~/work/histogram/IsoPhotonHadronCorrelations/RootFiles/NCentBinMix45/EMCAL_MB_0_90.root")); // Old Files with shsh = 0.40-2.00
      sShShNCentMix = "_ShSh0.40-2.00";
      systNMix = true;
    }

    cout << iCen << endl;
    Exec(ptTrMin, ptTrMax, iCen, bMirror, shshBkg, dirFiles, systPur, bZYAM, bPlot, phiMin, phiMax, systShSh, systNMix);
  }

  cout << "Close root files" << endl;
  fileDataShStd->Close();
  delete fileDataShStd;
  fileMC[0]->Close();
  delete fileMC[0];
  fileMC[1]->Close();
  delete fileMC[1];
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

void PlotStyle(TH1F *hPlot, int kMarker, int kColor, TString titleX, TString titleY)
{
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111111);
  // gStyle->SetPadTopMargin(0.10);
  hPlot->SetMarkerStyle(kMarker);
  hPlot->SetMarkerColor(kColor);
  hPlot->SetLineColor(kColor);
  hPlot->GetYaxis()->SetTitle(Form("%s", titleY.Data()));
  hPlot->GetXaxis()->SetTitle(Form("%s", titleX.Data()));
  // leg->SetFillColor(kWhite);
  // leg->SetLineColor(0);
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

void fZYAM(TH1F *hSame, double rangeMin = 3 * (TMath::Pi()) / 10, double rangeMax = TMath::Pi() / 2)
{
  int minBin = hSame->GetXaxis()->FindBin(rangeMin + 0.00001);
  int maxBin = hSame->GetXaxis()->FindBin(rangeMax);
  double IntSame = hSame->Integral(minBin, maxBin);
  // cout << "bincontent min: " << hSame->GetBinCenter(minBin) << ", bincontent max: " << hSame->GetBinCenter(maxBin) << endl;
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
