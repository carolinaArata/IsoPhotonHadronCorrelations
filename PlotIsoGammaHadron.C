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
#include "Plotting.h"

using std::cout;
using std::endl;
Int_t nCen = 4;
Int_t cenBins[] = {0, 10, 30, 50, 90};
// Int_t cenBins[] = {50, 90};

// int nZtBin = 10;
// double assocZt[] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20};

int nZtBin = 6;
double assocZt[] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 1.00};

int npt = 13;
float ptTrig[] = {12, 14, 16, 18, 20, 25, 30, 35, 40, 50, 60, 80};

Int_t const nIso = 2;
Int_t const nShSh = 2;
int const nShShMC = 1;

// Int_t kColorData[2][2] = {{}{}};
Int_t kMarkStyle[2] = {24, 25};
Int_t kMarkStyleNoUE[2] = {20, 53};
Int_t kMarkStyleIso_NotIso[2] = {21, 48};
Int_t kMarkStyleShSh[2] = {20, 25};

//{2, 4, kOrange + 7, kAzure + 3, kPink - 4, kViolet - 7, kBlue + 2, kTeal - 6};
Int_t kColorMC[2][2] = {};
Int_t kColorIsoREC[] = {kAzure + 2, kOrange - 3};
Int_t kColorShSh[] = {kAzure + 2, kViolet - 2};

// void PlotStyle(TH1F *hPlot, int kMarker, double kMarkerSize, int kColor, TString titleX, TString titleY);
// TLatex *LatexStd(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax, bool sPttrig);
// TLatex *LatexDPhi(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax);
// TLegend *LegStd(TLegend *leg, double xpos1, double ypos1, double xpos2, double ypos2);

void PlotIsoGammaHadron(float ptMin = 18, float ptMax = 40, TString outDirPlot = "Output_FigcheckCode", TString shshBkg = "0.40-1.00", TString dirFiles = "Output_checkCode")
{

  TString shshString[2] = {"0.10-0.30", shshBkg};
  TString shshStringMC[2] = {"0.10-0.30", "0.10-0.30"};
  TString sGenPartMC[2] = {"Photon", "Pi0"};
  TString sPtAll = Form("_Pt%2.0f_%2.0f", ptMin, ptMax);

  // index pt start and stop
  int nsize = sizeof(ptTrig) / sizeof(ptTrig[0]);
  auto itr1 = find(ptTrig, ptTrig + nsize, ptMin);
  auto itr2 = find(ptTrig, ptTrig + nsize, ptMax);
  int index1 = distance(ptTrig, itr1);
  int index2 = distance(ptTrig, itr2);
  int nPtTrig = index2 - index1;
  cout << index1 << "___" << index2 << ", " << nPtTrig << endl;

  int xNumPad;
  int yNumPad;
  int legPad;
  TFile *fPlot[nCen];
  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////Definition of histograms to be plotted//////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  // Data
  TH1F *hTriggerSam[nCen][nIso][nShSh];
  TH1F *hTriggerMix[nCen][nIso][nShSh];
  TH1F *hdPhiSam[nCen][nIso][nShSh][nZtBin][nPtTrig];
  TH1F *hdPhiMix[nCen][nIso][nShSh][nZtBin][nPtTrig];
  TH1F *hdPhiSamNoUE[nCen][nIso][nShSh][nZtBin][nPtTrig];
  TH1F *hdPhiSamNoUERatio[nCen][nIso][nZtBin][nPtTrig];
  TH1F *hdPhiSamNoUEPtAll[nCen][nIso][nShSh][nZtBin];
  TH1F *hdPhiSamPi0Pur[nCen][nIso][nZtBin][nPtTrig];
  TH1F *hdPhiPhoton[nCen][nIso][nZtBin][nPtTrig];
  // Alternative plotting for purity correction
  TH1F *hdPhiSameNoUEScaledbyPur[nCen][nIso][nZtBin][nPtTrig];
  TH1F *hdPhiSamePi0ScaledbyPur[nCen][nIso][nZtBin][nPtTrig];
  TH1F *hdPhiPhotonScaledbyPur[nCen][nIso][nZtBin][nPtTrig];

  TH1F *hZtPhotonPtBin[nCen][nIso][nPtTrig];
  TH1F *hZtPhotonPtBinNearSide[nCen][nIso][nPtTrig];
  TH1F *hZtPhoton[nCen][nIso];
  TH1F *hZtPi0PtBin[nCen][nIso][nPtTrig];
  TH1F *hZtPi0[nCen][nIso];

  // MC
  TH1F *hTriggerMCGen[nCen][nIso][nShSh];
  TH1F *hTriggerSamMCRec[nCen][nIso][nShSh];
  TH1F *hTriggerMixMCRec[nCen][nIso][nShSh];
  TH1F *hdPhiMCGenUE[nCen][nIso][nShSh][nZtBin][nPtTrig];
  TH1F *hdPhiMCGen[nCen][nIso][nShSh][nZtBin][nPtTrig];
  TH1F *hdPhiSamMCRec[nCen][nIso][nShSh][nZtBin][nPtTrig];
  TH1F *hdPhiMixMCRec[nCen][nIso][nShSh][nZtBin][nPtTrig];
  TH1F *hdPhiSamMCRecNoUE[nCen][nIso][nShSh][nZtBin][nPtTrig];
  TH1F *hZtPtBinMCGen[nCen][nIso][nShSh][nPtTrig];
  TH1F *hZtPtBinMCRec[nCen][nIso][nShSh][nPtTrig];
  TH1F *hRatioEffCorrPtBin[nCen][nIso][nShSh][nPtTrig];

  TH1F *hZtMCGen[nCen][nIso][nShSh];
  TH1F *hZtMCRec[nCen][nIso][nShSh];
  TH1F *hRatioEffCorr[nCen][nIso][nShSh];
  TH1F *hZtEffCorr[nCen][nIso][nShSh];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    fPlot[iCen] = new TFile(dirFiles + "/fPlot" + shshBkg + sCent + sPtAll + ".root");
    // DataPlot
    cout << "Cen: " << cenBins[iCen] << "-" << cenBins[iCen + 1] << ": Get Data plots + Set plots style" << endl;
    for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
    {
      TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
      for (int iso = 0; iso < nIso; iso++)
      {
        TString sIso = Form("Iso%d", iso);
        for (Int_t izt = 0; izt < nZtBin; izt++)
        {
          TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);
          cout << "Get azimuthal distributions for Same and Mixed X [Cent][iso][shsh][zt]" << endl;
          for (int iSh = 0; iSh < nShSh; iSh++)
          {
            TString sShSh = Form("_ShSh%s", shshString[iSh].Data());
            hdPhiSam[iCen][iso][iSh][izt][iptTr] = (TH1F *)fPlot[iCen]->Get("hdPhiSameMirror" + sIso + sShSh + sZtBin + sPtTrig);
            // cout << hdPhiSam[iCen][iso][iSh][izt][iptTr] << endl;
            hdPhiMix[iCen][iso][iSh][izt][iptTr] = (TH1F *)fPlot[iCen]->Get("hdPhiMixMirror" + sIso + sShSh + sZtBin + sPtTrig);
            // cout << hdPhiMix[iCen][iso][iSh][izt][iptTr] << endl;
            hdPhiSamNoUE[iCen][iso][iSh][izt][iptTr] = (TH1F *)fPlot[iCen]->Get("hdPhiSameNoUE" + sIso + sShSh + sZtBin + sPtTrig);
            // cout << hdPhiSamNoUE[iCen][iso][iSh][izt][iptTr] << endl;
            PlotStyle(hdPhiSam[iCen][iso][iSh][izt][iptTr], kMarkStyle[iSh], 1.3, kBlack, kBlue + 2, "|#Delta#it{#varphi}| (rad)", "1 /#color[0]{.}#it{N}^{#color[0]{..}trig} d^{2}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}|", false);
            PlotStyle(hdPhiMix[iCen][iso][iSh][izt][iptTr], 21, 0, kRed - 4, kBlue + 2, "|#Delta#it{#varphi}| (rad)", "1 /#color[0]{.}#it{N}^{#color[0]{..}trig} d^{2}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}|", false);
            PlotStyle(hdPhiSamNoUE[iCen][iso][iSh][izt][iptTr], kMarkStyleNoUE[iSh], 1.5, kColorShSh[iSh], kBlue + 2, "|#Delta#it{#varphi}| (rad)", "1 /#color[0]{.}#it{N}^{#color[0]{..}trig} d^{2}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}|", false);
          }
          hdPhiSamNoUERatio[iCen][iso][izt][iptTr] = (TH1F *)fPlot[iCen]->Get(Form("hdPhiSameNoUERatio%s%s_%s", sIso.Data(), sZtBin.Data(), sPtTrig.Data())); // ratio between photon and pi0
          cout << "Get azimuthal distributions after purity correction" << endl;
          cout << hdPhiSamNoUERatio[iCen][iso][izt][iptTr] << endl;
          hdPhiSamPi0Pur[iCen][iso][izt][iptTr] = (TH1F *)fPlot[iCen]->Get(Form("hdPhiSamePi0Pur%s%s_%s", sIso.Data(), sZtBin.Data(), sPtTrig.Data())); //(1-P)IsoPi0Deltaphi
          cout << hdPhiSamPi0Pur[iCen][iso][izt][iptTr] << endl;
          hdPhiPhoton[iCen][iso][izt][iptTr] = (TH1F *)fPlot[iCen]->Get(Form("hdPhi%sPhoton%s_%s", sIso.Data(), sZtBin.Data(), sPtTrig.Data()));
          cout << hdPhiPhoton[iCen][iso][izt][iptTr] << endl;
          cout << "Set plots style" << endl;
          PlotStyle(hdPhiSamNoUERatio[iCen][iso][izt][iptTr], kMarkStyle[1], 2, kBlue + 2, kBlue + 2, "|#Delta#it{#varphi}| (rad)", "Ratio", false);
          PlotStyle(hdPhiSamPi0Pur[iCen][iso][izt][iptTr], 24, 1.5, kViolet + 5, kBlue + 2, "|#Delta#it{#varphi}| (rad)", "1 /#color[0]{.}#it{N}^{#color[0]{..}trig} d^{2}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}|", false);
          PlotStyle(hdPhiPhoton[iCen][iso][izt][iptTr], 21, 1.5, kOrange + 8, kOrange + 8, "|#Delta#it{#varphi}| (rad)", "1 /#color[0]{.}#it{N}^{#color[0]{..}trig} d^{2}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}|", false);
          
          cout<<"Alternative plotting for Purity correction"<<endl;
          hdPhiSameNoUEScaledbyPur[iCen][iso][izt][iptTr] = (TH1F *)fPlot[iCen]->Get(Form("hdPhiSameNoUEScaledbyPur%s%s_%s", sIso.Data(), sZtBin.Data(), sPtTrig.Data())); //(1-P)IsoPi0Deltaphi
          cout << hdPhiSameNoUEScaledbyPur[iCen][iso][izt][iptTr] << endl;
          hdPhiSamePi0ScaledbyPur[iCen][iso][izt][iptTr] = (TH1F *)fPlot[iCen]->Get(Form("hdPhiSamePi0ScaledbyPur%s%s_%s", sIso.Data(), sZtBin.Data(), sPtTrig.Data())); //(1-P)IsoPi0Deltaphi
          cout << hdPhiSamePi0ScaledbyPur[iCen][iso][izt][iptTr] << endl;
          hdPhiPhotonScaledbyPur[iCen][iso][izt][iptTr] = (TH1F *)fPlot[iCen]->Get(Form("hdPhi%sPhoton%s_%sScaledbyPur", sIso.Data(), sZtBin.Data(), sPtTrig.Data()));
          cout << hdPhiPhotonScaledbyPur[iCen][iso][izt][iptTr] << endl;
          cout << "Set plots style Alternative plotting" << endl;
          PlotStyle(hdPhiSameNoUEScaledbyPur[iCen][iso][izt][iptTr], kMarkStyle[1], 2, kBlue + 2, kBlue + 2, "|#Delta#it{#varphi}| (rad)", "Ratio", false);
          PlotStyle(hdPhiSamePi0ScaledbyPur[iCen][iso][izt][iptTr], 24, 1.5, kViolet + 5, kBlue + 2, "|#Delta#it{#varphi}| (rad)", "1 /#color[0]{.}#it{N}^{#color[0]{..}trig} d^{2}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}|", false);
          PlotStyle(hdPhiPhotonScaledbyPur[iCen][iso][izt][iptTr], 21, 1.5, kOrange + 8, kOrange + 8, "|#Delta#it{#varphi}| (rad)", "1 /#color[0]{.}#it{N}^{#color[0]{..}trig} d^{2}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}|", false);
          
        }
        cout << "Get Zt distributions for different Pt bins" << endl;
        hZtPhotonPtBin[iCen][iso][iptTr] = (TH1F *)fPlot[iCen]->Get(Form("hZt%sPhotonPtBin_%s", sIso.Data(), sPtTrig.Data()));
        PlotStyle(hZtPhotonPtBin[iCen][iso][iptTr], kMarkStyle[0], 1, kAzure + 7, kBlue + 2, "#it{z}_{T}", "1/N^{trig}dN^{charg}/d#it{z}_{T}", false);

        hZtPi0PtBin[iCen][iso][iptTr] = (TH1F *)fPlot[iCen]->Get(Form("hZt%sPi0PtBin_%s", sIso.Data(), sPtTrig.Data()));
        PlotStyle(hZtPi0PtBin[iCen][iso][iptTr], kMarkStyle[1], 1, kAzure + 7, kBlue + 2, "#it{z}_{T}", "1/N^{trig}dN^{charg}/d#it{z}_{T}", false);

        hZtPhotonPtBinNearSide[iCen][iso][iptTr] = (TH1F *)fPlot[iCen]->Get(Form("hZt%sPhotonPtBin_%sNearSide", sIso.Data(), sPtTrig.Data()));
        PlotStyle(hZtPhotonPtBinNearSide[iCen][iso][iptTr], kMarkStyle[0], 1, kAzure + 7, kBlue + 2, "#it{z}_{T}", "1/N^{trig}dN^{charg}/d#it{z}_{T}", false);
      }
    }

    for (Int_t izt = 0; izt < nZtBin; izt++)
    {
      TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);
      for (int iso = 1; iso < nIso; iso++)
      {
        TString sIso = Form("Iso%d", iso);
        for (int iSh = 0; iSh < nShSh; iSh++)
        {
          TString sShSh = Form("_ShSh%s", shshString[iSh].Data());
          hdPhiSamNoUEPtAll[iCen][iso][iSh][izt] = (TH1F *)fPlot[iCen]->Get(Form("hdPhiSamNoUEPtAll" + sIso + sShSh + "%s%s%s", sCent.Data(), sZtBin.Data(), sPtAll.Data()));
          PlotStyle(hdPhiSamNoUEPtAll[iCen][iso][iSh][izt], kMarkStyleNoUE[1], 2, kAzure + 2, kBlue + 2, "|#Delta#it{#varphi}| (rad)", "1 / #it{N}^{ trig} d^{2}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}|", false);
        }
      }
    }
    cout << "Get Zt distributions for all Pt range" << endl;
    for (int iso = 0; iso < nIso; iso++)
    {
      TString sIso = Form("Iso%d", iso);
      hZtPhoton[iCen][iso] = (TH1F *)fPlot[iCen]->Get(Form("hZt%sPhoton%s%s", sIso.Data(), sCent.Data(), sPtAll.Data()));
      hZtPi0[iCen][iso] = (TH1F *)fPlot[iCen]->Get(Form("hZt%sPi0%s%s", sIso.Data(), sCent.Data(), sPtAll.Data()));

      PlotStyle(hZtPhoton[iCen][iso], kMarkStyle[0], 1, kAzure + 7, kBlue + 2, "#it{z}_{T}", "1/N^{trig}dN^{charg}/d#it{z}_{T}", false);
      PlotStyle(hZtPi0[iCen][iso], kMarkStyle[1], 1, kAzure + 7, kBlue + 2, "#it{z}_{T}", "1/N^{trig}dN^{charg}/d#it{z}_{T}", false);
    }

    cout << "Get MonteCarlo plots + Set plots style" << endl;
    // MC Plot
    for (int iso = 0; iso < nIso; iso++)
    {
      TString sIso = Form("Iso%d", iso);
      for (int iSh = 0; iSh < nShShMC; iSh++)
      {
        TString sShSh = Form("_ShSh%s", shshStringMC[iSh].Data());
        for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
        {
          TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
          for (Int_t izt = 0; izt < nZtBin; izt++)
          {
            TString sZtBin = Form("ZTBin%1.2f_%1.2f", assocZt[izt], assocZt[izt + 1]);

            hdPhiMCGenUE[iCen][iso][iSh][izt][iptTr] = (TH1F *)fPlot[iCen]->Get("hdPhiMCGenMirrorUE" + sIso + sShSh + sZtBin + sPtTrig + sGenPartMC[iSh]);
            hdPhiMCGen[iCen][iso][iSh][izt][iptTr] = (TH1F *)fPlot[iCen]->Get("hdPhiMCGenMirror" + sIso + sShSh + sZtBin + sPtTrig + sGenPartMC[iSh]);
            hdPhiSamMCRec[iCen][iso][iSh][izt][iptTr] = (TH1F *)fPlot[iCen]->Get("hdPhiSameMCRecMirror" + sIso + sShSh + sZtBin + sPtTrig + sGenPartMC[iSh]);
            hdPhiMixMCRec[iCen][iso][iSh][izt][iptTr] = (TH1F *)fPlot[iCen]->Get("hdPhiMixMCRecMirror" + sIso + sShSh + sZtBin + sPtTrig + sGenPartMC[iSh]);
            hdPhiSamMCRecNoUE[iCen][iso][iSh][izt][iptTr] = (TH1F *)fPlot[iCen]->Get("hdPhiSameMCRecNoUE" + sIso + sShSh + sZtBin + sPtTrig + sGenPartMC[iSh]);

            PlotStyle(hdPhiMCGenUE[iCen][iso][iSh][izt][iptTr], 25, 2, kBlue + 2, kBlue + 2, "|#Delta#it{#varphi}| (rad)", "1 / #it{N}^{ trig} d^{2}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}|", false);
            PlotStyle(hdPhiMCGen[iCen][iso][iSh][izt][iptTr], 21, 2, kBlue + 2, kBlue + 2, "|#Delta#it{#varphi}| (rad)", "1 / #it{N}^{ trig} d^{2}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}|", false);
            PlotStyle(hdPhiSamMCRec[iCen][iso][iSh][izt][iptTr], 21, 2, kPink + 6, kBlue + 2, "|#Delta#it{#varphi}| (rad)", "1 / #it{N}^{ trig} d^{2}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}|", false);
            PlotStyle(hdPhiMixMCRec[iCen][iso][iSh][izt][iptTr], kMarkStyleIso_NotIso[iso], 2, kRed - 3, kBlue + 2, "|#Delta#it{#varphi}| (rad)", "1 / #it{N}^{ trig} d^{2}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}|", false);
            PlotStyle(hdPhiSamMCRecNoUE[iCen][iso][iSh][izt][iptTr], 21, 2, kOrange + 7, kBlue + 2, "|#Delta#it{#varphi}| (rad)", "1 / #it{N}^{ trig} d^{2}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}|", false);
          }

          hZtPtBinMCGen[iCen][iso][iSh][iptTr] = (TH1F *)fPlot[iCen]->Get(Form("hZtPtBinMCGen_%s%s%s", sIso.Data(), sPtTrig.Data(), sGenPartMC[iSh].Data()));
          hZtPtBinMCRec[iCen][iso][iSh][iptTr] = (TH1F *)fPlot[iCen]->Get(Form("hZtPtBinMCRec_%s%s%s", sIso.Data(), sPtTrig.Data(), sGenPartMC[iSh].Data()));
          hRatioEffCorrPtBin[iCen][iso][iSh][iptTr] = (TH1F *)fPlot[iCen]->Get(Form("hRatioEffCorrPtBin%s%s%s", sIso.Data(), sPtTrig.Data(), sGenPartMC[iSh].Data()));
          PlotStyle(hZtPtBinMCGen[iCen][iso][iSh][iptTr], kMarkStyleShSh[iSh], 1, kAzure + 7, kBlue + 2, "#it{z}_{T}", "1 / #it{N}^{ trig} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
          PlotStyle(hZtPtBinMCRec[iCen][iso][iSh][iptTr], kMarkStyleShSh[iSh], 1, kOrange + 7, kBlue + 2, "#it{z}_{T}", "1 / #it{N}^{ trig} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
          PlotStyle(hRatioEffCorrPtBin[iCen][iso][iSh][iptTr], kMarkStyleShSh[iSh], 1, kTeal - 7, kBlue + 2, "#it{z}_{T}", "1/N^{trig}dN^{charg}/d#it{z}_{T}", false);
        }

        hZtMCGen[iCen][iso][iSh] = (TH1F *)fPlot[iCen]->Get(Form("hZtMCGen%s%s%s%s", sIso.Data(), sGenPartMC[iSh].Data(), sCent.Data(), sPtAll.Data()));
        hZtMCRec[iCen][iso][iSh] = (TH1F *)fPlot[iCen]->Get(Form("hZtMCRec%s%s%s%s", sIso.Data(), sGenPartMC[iSh].Data(), sCent.Data(), sPtAll.Data()));
        hZtEffCorr[iCen][iso][iSh] = (TH1F *)fPlot[iCen]->Get(Form("hZtEffCorr%s%s%s%s", sIso.Data(), sGenPartMC[iSh].Data(), sCent.Data(), sPtAll.Data()));
        PlotStyle(hZtMCGen[iCen][iso][iSh], kMarkStyleShSh[iSh], 1, kAzure + 7, kBlue + 2, "#it{z}_{T}", "1 / #it{N}^{ trig} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
        PlotStyle(hZtMCRec[iCen][iso][iSh], kMarkStyleShSh[iSh], 1, kOrange + 7, kBlue + 2, "#it{z}_{T}", "1 / #it{N}^{ trig} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
        PlotStyle(hZtEffCorr[iCen][iso][iSh], kMarkStyleShSh[iSh], 1, kCyan + 2, kBlue + 2, "#it{z}_{T}", "1 / #it{N}^{ trig} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
      }
    }
  }

  cout << "Azimuthal distributions Data" << endl;
  gSystem->Exec(Form("mkdir %s", outDirPlot.Data()));

  TString processline = Form(".! mkdir -pv %s", outDirPlot.Data());
  gROOT->ProcessLine(processline.Data());

  cout << "Directory: " << outDirPlot.Data() << endl;

  TCanvas *cSame_MixIsoClust[nCen][nPtTrig];
  TLegend *legSame_MixIsoClust[nCen][nPtTrig];
  TCanvas *cSame_MixIsoPi0[nCen][nPtTrig];
  TLegend *legSame_MixIsoPi0[nCen][nPtTrig];
  TCanvas *cSameNoUE_pTBin_pTAllIsoClust[nCen];
  TLegend *legSameNoUE_pTBin_pTAllIsoClust[nCen];
  TCanvas *cMixIsoGamma_Clust[nCen][nPtTrig];
  TLegend *legMixIsoGamma_Clust[nCen][nPtTrig];
  TCanvas *cSame_Mix_NoUEIsoClust[nCen][nPtTrig];
  TCanvas *cSame_Mix_NoUEIsoClustSingleZt[nCen][nPtTrig];
  TLegend *legSame_Mix_NoUEIsoClust[nCen][nPtTrig];
  TCanvas *cIsoClust_Pi0NoUE[nCen][nPtTrig];
  TLegend *legIsoClust_Pi0NoUE[nCen][nPtTrig];
  TCanvas *cIsoClust_Pi0NoUERatio[nCen][nPtTrig];
  TLegend *legIsoClust_Pi0NoUERatio[nCen][nPtTrig];
  TCanvas *cIsoClust_Pi0Pur[nCen][nPtTrig];
  TCanvas *cIsoClust_Pi0Pur_SingleZt[nCen][nPtTrig];
  TLegend *legIsoClust_Pi0Pur[nCen][nPtTrig];
  TLatex *latDphi[nCen];

  TGraph *lineX0 = DrawLine(lineX0, 0, 0, TMath::Pi(), 0);

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    if (iCen == 0 || iCen == 1)
    {
      nZtBin = 6;
      xNumPad = 4;
      legPad = 7;
    }
    else if (iCen == 2 || iCen == 3)
    {
      nZtBin = 5;
      xNumPad = 3;
      legPad = 6;
    }

    for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
    {
      TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
      cSame_MixIsoClust[iCen][iptTr] = new TCanvas("cSame_MixIsoClust" + sCent + sPtTrig, "cSame_MixIsoClust" + sCent + sPtTrig, xNumPad * 800, 2 * 600);
      cSame_MixIsoClust[iCen][iptTr]->Divide(xNumPad, 2);
      legSame_MixIsoClust[iCen][iptTr] = LegStd(legSame_MixIsoClust[iCen][iptTr], 0.10, 0.11, 0.12, 0.46);
      cSame_MixIsoPi0[iCen][iptTr] = new TCanvas("cSame_MixIsoPi0" + sCent + sPtTrig, "cSame_MixIsoPi0" + sCent + sPtTrig, xNumPad * 800, 2 * 600);
      cSame_MixIsoPi0[iCen][iptTr]->Divide(xNumPad, 2);
      legSame_MixIsoPi0[iCen][iptTr] = LegStd(legSame_MixIsoPi0[iCen][iptTr], 0.10, 0.11, 0.12, 0.46);
      cMixIsoGamma_Clust[iCen][iptTr] = new TCanvas("cMixIsoGamma_Pi0" + sCent + sPtTrig, "ccMixIsoGamma_Pi0" + sCent + sPtTrig, xNumPad * 800, 2 * 600);
      cMixIsoGamma_Clust[iCen][iptTr]->Divide(xNumPad, 2);
      legMixIsoGamma_Clust[iCen][iptTr] = LegStd(legMixIsoGamma_Clust[iCen][iptTr], 0.10, 0.11, 0.12, 0.46);
      cSame_Mix_NoUEIsoClust[iCen][iptTr] = new TCanvas("cSame_Mix_NoUEIsoClust" + sCent + sPtTrig, "cSame_Mix_NoUEIsoClust" + sCent + sPtTrig, xNumPad * 800, 2 * 600);
      cSame_Mix_NoUEIsoClust[iCen][iptTr]->Divide(xNumPad, 2);

      cSame_Mix_NoUEIsoClustSingleZt[iCen][iptTr] = new TCanvas("cSame_Mix_NoUEIsoClustSingleZt" + sCent + sPtTrig, "cSame_Mix_NoUEIsoClustSingleZt" + sCent + sPtTrig, 2 * 800, 1 * 600);
      cSame_Mix_NoUEIsoClustSingleZt[iCen][iptTr]->Divide(2, 1);

      legSame_Mix_NoUEIsoClust[iCen][iptTr] = LegStd(legSame_Mix_NoUEIsoClust[iCen][iptTr], 0.10, 0.10, 0.12, 0.465);
      cIsoClust_Pi0NoUE[iCen][iptTr] = new TCanvas("cIsoClust_Pi0NoUE" + sCent + sPtTrig, "cIsoClust_Pi0NoUE" + sCent + sPtTrig, xNumPad * 800, 2 * 600);
      cIsoClust_Pi0NoUE[iCen][iptTr]->Divide(xNumPad, 2);
      legIsoClust_Pi0NoUE[iCen][iptTr] = LegStd(legIsoClust_Pi0NoUE[iCen][iptTr], 0.10, 0.13, 0.12, 0.44);
      cIsoClust_Pi0NoUERatio[iCen][iptTr] = new TCanvas("cIsoClust_Pi0NoUERatio" + sCent + sPtTrig, "cIsoClust_Pi0NoUERatio" + sCent + sPtTrig, xNumPad * 800, 2 * 600);
      cIsoClust_Pi0NoUERatio[iCen][iptTr]->Divide(xNumPad, 2);
      legIsoClust_Pi0NoUERatio[iCen][iptTr] = LegStd(legIsoClust_Pi0NoUERatio[iCen][iptTr], 0.12, 0.24, 0.12, 0.40);

      cIsoClust_Pi0Pur[iCen][iptTr] = new TCanvas("cIsoClust_Pi0Pur" + sCent + sPtTrig, "cIsoClust_Pi0Pur" + sCent + sPtTrig, xNumPad * 800, 2 * 600);
      cIsoClust_Pi0Pur[iCen][iptTr]->Divide(xNumPad, 2);
      legIsoClust_Pi0Pur[iCen][iptTr] = LegStd(legIsoClust_Pi0Pur[iCen][iptTr], 0.10, 0.08, 0.12, 0.435);

      cIsoClust_Pi0Pur_SingleZt[iCen][iptTr] = new TCanvas("cIsoClust_Pi0Pur_SingleZt" + sCent + sPtTrig, "cIsoClust_Pi0Pur_SingleZt" + sCent + sPtTrig, 2 * 800, 1 * 600);
      cIsoClust_Pi0Pur_SingleZt[iCen][iptTr]->Divide(2, 1);

      for (Int_t izt = 0; izt < nZtBin; izt++)
      {
        // TString sTitle = Form("%2.0f < #it{p}_{T}^{tr} < %2.0f GeV/#it{c}, %2.2f < #it{z}_{T}^{as} < %2.2f ", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], assocZt[izt], assocZt[izt + 1]);
        TString sTitle = Form(" %2.2f < #it{z}_{T} < %2.2f ", assocZt[izt], assocZt[izt + 1]);
        gStyle->SetPadTopMargin(0.06);
        // gStyle->SetPadRightMargin(0.05);
        gStyle->SetPadLeftMargin(0.15);
        gStyle->SetPadBottomMargin(0.15);
        // gStyle->SetTitleX(0.56);
        cSame_MixIsoClust[iCen][iptTr]->cd(izt + 1);
        TGaxis::SetMaxDigits(1);
        hdPhiSam[iCen][1][0][izt][iptTr]->SetTitle(sTitle);
        // hdPhiSam[iCen][1][0][izt][iptTr]->SetMinimum(0.9 * hdPhiMix[iCen][1][0][izt][iptTr]->GetMinimum());
        hdPhiSam[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(0.7 * hdPhiMix[iCen][1][0][izt][iptTr]->GetMinimum(), 1.25 * hdPhiSam[iCen][1][0][izt][iptTr]->GetMaximum());
        // hdPhiSam[iCen][1][0][izt][iptTr]->SetMaximum(1.2 * hdPhiSam[iCen][1][0][izt][iptTr]->GetMaximum());
        hdPhiSam[iCen][1][0][izt][iptTr]->Draw("sameE0");
        hdPhiMix[iCen][1][0][izt][iptTr]->Draw("same");

        cSame_MixIsoPi0[iCen][iptTr]->cd(izt + 1);

        hdPhiSam[iCen][1][1][izt][iptTr]->SetTitle(sTitle);
        hdPhiSamNoUE[iCen][1][1][izt][iptTr]->GetYaxis()->SetRangeUser(-0.2 * hdPhiSam[iCen][1][1][izt][iptTr]->GetMaximum(), 1.2 * hdPhiSam[iCen][1][1][izt][iptTr]->GetMaximum());
        // hdPhiSamNoUE[iCen][1][1][izt][iptTr]->SetMaximum(1.2 * hdPhiSam[iCen][1][1][izt][iptTr]->GetMaximum());
        // hdPhiSam[iCen][1][1][izt][iptTr]->SetMinimum(-5e-3);
        hdPhiSamNoUE[iCen][1][1][izt][iptTr]->Draw("same");
        hdPhiSam[iCen][1][1][izt][iptTr]->Draw("sameE0");
        hdPhiMix[iCen][1][1][izt][iptTr]->Draw("same");
        lineX0->Draw("same");

        cMixIsoGamma_Clust[iCen][iptTr]->cd(izt + 1);
        hdPhiMix[iCen][1][0][izt][iptTr]->Draw("same");
        hdPhiMix[iCen][1][1][izt][iptTr]->Draw("same");
        // line->Draw("same");

        cSame_Mix_NoUEIsoClust[iCen][iptTr]->cd(izt + 1);
        // hdPhiSamNoUE[iCen][1][0][izt][iptTr]->SetMinimum(-5e-3);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-0.2 * hdPhiSam[iCen][1][0][izt][iptTr]->GetMaximum(), 1.2 * hdPhiSam[iCen][1][0][izt][iptTr]->GetMaximum());
        TAxis* aNoUE = hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetXaxis();
        aNoUE->SetNdivisions(-505);
        // Replace first and last labels
        aNoUE->ChangeLabel(1, -1, -1, -1, -1, -1, "0");
        aNoUE->ChangeLabel(-1, -1, -1, -1, -1, -1, "#pi");
        
        // Replace the intermediate ones (2,3,4,5)
        aNoUE->ChangeLabel(2, -1, -1, -1, -1, -1, "#pi/5");
        aNoUE->ChangeLabel(3, -1, -1, -1, -1, -1, "2#pi/5");
        aNoUE->ChangeLabel(4, -1, -1, -1, -1, -1, "3#pi/5");
        aNoUE->ChangeLabel(5, -1, -1, -1, -1, -1, "4#pi/5");
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->Draw("same");
        hdPhiSam[iCen][1][0][izt][iptTr]->Draw("sameE0");
        hdPhiMix[iCen][1][0][izt][iptTr]->Draw("same");
        // cSame_Mix_NoUEIsoClust[iCen][iptTr]->cd(izt + 1)->Update();
        lineX0->Draw("same");

        cIsoClust_Pi0NoUE[iCen][iptTr]->cd(izt + 1);
        hdPhiSamNoUE[iCen][1][1][izt][iptTr]->SetTitle(sTitle);
        cout << "MAXIMUM: " << hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetMaximum() << "_____" << hdPhiSamNoUE[iCen][1][1][izt][iptTr]->GetMaximum() << endl;
        // hdPhiSamNoUE[iCen][1][1][izt][iptTr]->SetMaximum(1.2 * hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetMaximum());
        hdPhiSamNoUE[iCen][1][1][izt][iptTr]->Draw("same");
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->Draw("same");

        cIsoClust_Pi0NoUERatio[iCen][iptTr]->cd(izt + 1);
        hdPhiSamNoUERatio[iCen][1][izt][iptTr]->SetTitle(sTitle);
        hdPhiSamNoUERatio[iCen][1][izt][iptTr]->SetMaximum(2.5);
        hdPhiSamNoUERatio[iCen][1][izt][iptTr]->SetMinimum(-2.5);
        hdPhiSamNoUERatio[iCen][1][izt][iptTr]->Draw("same");
        // line->Draw("same");

        cIsoClust_Pi0Pur[iCen][iptTr]->cd(izt + 1);
        hdPhiPhoton[iCen][1][izt][iptTr]->SetTitle(sTitle);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->SetTitle(sTitle);
        // hdPhiPhoton[iCen][izt][iptTr]->GetYaxis()->SetRangeUser(-1*hdPhiSamPi0Pur[iCen][izt][iptTr]->GetMaximum(), hdPhiSamPi0Pur[iCen][izt][iptTr]->GetMaximum());
        //  cout << "MAXXXXX: " << hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetMaximum() << endl;

        // cout << (hdPhiPhoton[iCen][izt][iptTr]->GetMinimum() - hdPhiPhoton[iCen][izt][iptTr]->GetBinError(5)) << endl;
        // hdPhiSamNoUE[iCen][1][0][izt][iptTr]->SetMinimum(0.1 * (hdPhiPhoton[iCen][izt][iptTr]->GetMinimum() - hdPhiPhoton[iCen][izt][iptTr]->GetBinError(6)));
        if (iCen == 0 && izt == 1)
        {
          hdPhiPhoton[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-60 * 1e-3, 65 * 1e-3);
        }
        TAxis* aPhoton = hdPhiPhoton[iCen][1][izt][iptTr]->GetXaxis();
        aPhoton->SetNdivisions(-505);
        // Replace first and last labels
        aPhoton->ChangeLabel(1, -1, -1, -1, -1, -1, "0");
        aPhoton->ChangeLabel(-1, -1, -1, -1, -1, -1, "#pi");
        
        // Replace the intermediate ones (2,3,4,5)
        aPhoton->ChangeLabel(2, -1, -1, -1, -1, -1, "#pi/5");
        aPhoton->ChangeLabel(3, -1, -1, -1, -1, -1, "2#pi/5");
        aPhoton->ChangeLabel(4, -1, -1, -1, -1, -1, "3#pi/5");
        aPhoton->ChangeLabel(5, -1, -1, -1, -1, -1, "4#pi/5");

        hdPhiPhoton[iCen][1][izt][iptTr]->Draw("same");

        hdPhiSamPi0Pur[iCen][1][izt][iptTr]->Draw("same");
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->Draw("same");
        hdPhiPhoton[iCen][1][izt][iptTr]->Draw("same");
        lineX0->Draw("same");
      }
      cSame_MixIsoClust[iCen][iptTr]->cd(legPad);
      latDphi[iCen] = LatexDPhi(latDphi[iCen], 0.01, 0.84, cenBins[iCen], cenBins[iCen + 1], ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], true);
      latDphi[iCen]->DrawLatex(0.01, 0.84 - 3 * 0.10, "0.10 < #it{#sigma}^{2}_{long , 5x5} < 0.30, #it{p}_{T}^{iso} < 1.5 GeV/#it{c}, #it{R} = 0.2");
      legSame_MixIsoClust[iCen][iptTr]->AddEntry(hdPhiSam[iCen][1][0][0][0], " ", "lep");
      latDphi[iCen]->DrawLatex(0.15, 0.355, "Same Event");
      legSame_MixIsoClust[iCen][iptTr]->AddEntry(hdPhiMix[iCen][1][0][0][0], " ", "lep");
      latDphi[iCen]->DrawLatex(0.15, 0.205, "Mixed Event");
      legSame_MixIsoClust[iCen][iptTr]->Draw("same");
      cSame_MixIsoClust[iCen][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/Same_MixIsoClustGamma" + sCent + sPtTrig + ".pdf");

      cSame_MixIsoPi0[iCen][iptTr]->cd(legPad);
      latDphi[iCen] = LatexDPhi(latDphi[iCen], 0.08, 0.84, cenBins[iCen], cenBins[iCen + 1], ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], true);
      latDphi[iCen]->DrawLatex(0.08, 0.84 - 3 * 0.10, "cluster^{iso}_{wide}: 0.40 < #it{#sigma}^{2}_{long , 5x5} < 1.00");
      legSame_MixIsoPi0[iCen][iptTr]->AddEntry(hdPhiSam[iCen][1][1][0][0], " ", "lep");
      latDphi[iCen]->DrawLatex(0.15, 0.385, "Same Event");
      legSame_MixIsoPi0[iCen][iptTr]->AddEntry(hdPhiMix[iCen][1][1][0][0], " ", "lep");
      latDphi[iCen]->DrawLatex(0.15, 0.266, "Mixed Event");
      legSame_MixIsoPi0[iCen][iptTr]->AddEntry(hdPhiSamNoUE[iCen][1][1][0][0], " ", "lep");
      latDphi[iCen]->DrawLatex(0.15, 0.147, "Same Event #font[122]{-} Mixed Event");
      legSame_MixIsoPi0[iCen][iptTr]->Draw("same");
      cSame_MixIsoPi0[iCen][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/Same_MixIsoClustPi0" + sCent + sPtTrig + ".pdf");

      cMixIsoGamma_Clust[iCen][iptTr]->cd(legPad);
      legMixIsoGamma_Clust[iCen][iptTr]->AddEntry(hdPhiMix[iCen][1][0][0][0], " ", "lep");
      latDphi[iCen]->DrawLatex(0.35, 0.46, "Mixed Event cluster^{iso}_{narrow}");
      legMixIsoGamma_Clust[iCen][iptTr]->AddEntry(hdPhiMix[iCen][1][1][0][0], " ", "lep");
      latDphi[iCen]->DrawLatex(0.35, 0.34, "Mixed Event cluster^{iso}_{wide}");
      legMixIsoGamma_Clust[iCen][iptTr]->Draw("same");
      cMixIsoGamma_Clust[iCen][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/MixIsoGamma_Pi0" + sCent + sPtTrig + ".pdf");

      cSame_Mix_NoUEIsoClust[iCen][iptTr]->cd(legPad);
      latDphi[iCen] = LatexDPhi(latDphi[iCen], 0.08, 0.84, cenBins[iCen], cenBins[iCen + 1], ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], true);
      latDphi[iCen]->DrawLatex(0.08, 0.84 - 3 * 0.10, "cluster^{iso}_{narrow}: 0.10 < #it{#sigma}^{2}_{long , 5x5} < 0.30");
      legSame_Mix_NoUEIsoClust[iCen][iptTr]->AddEntry(hdPhiSam[iCen][1][0][0][0], " ", "lep");
      latDphi[iCen]->DrawLatex(0.15, 0.385, "Same Event");
      legSame_Mix_NoUEIsoClust[iCen][iptTr]->AddEntry(hdPhiMix[iCen][1][0][0][0], " ", "lep");
      latDphi[iCen]->DrawLatex(0.15, 0.266, "Mixed Event");
      legSame_Mix_NoUEIsoClust[iCen][iptTr]->AddEntry(hdPhiSamNoUE[iCen][1][0][0][0], " ", "lep");
      latDphi[iCen]->DrawLatex(0.15, 0.147, "Same Event - Mixed Event");
      legSame_Mix_NoUEIsoClust[iCen][iptTr]->Draw("same");
      cSame_Mix_NoUEIsoClust[iCen][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/Same_MixIsoClustGamma_NoUE" + sCent + sPtTrig + ".pdf");

      cIsoClust_Pi0NoUE[iCen][iptTr]->cd(legPad);
      latDphi[iCen] = LatexDPhi(latDphi[iCen], 0.08, 0.87, cenBins[iCen], cenBins[iCen + 1], ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], true);
      latDphi[iCen]->DrawLatex(0.08, 0.87 - 3 * 0.10, "cluster^{iso}_{narrow}: 0.10 < #it{#sigma}^{2}_{long , 5x5} < 0.30");
      latDphi[iCen]->DrawLatex(0.08, 0.87 - 4 * 0.10, "cluster^{iso}_{wide}: 0.40 < #it{#sigma}^{2}_{long , 5x5} < 1.00");
      legIsoClust_Pi0NoUE[iCen][iptTr]->AddEntry(hdPhiSamNoUE[iCen][1][0][0][0], " ", "lep");
      latDphi[iCen]->DrawLatex(0.15, 0.355, "cluster^{iso}_{narrow}");
      legIsoClust_Pi0NoUE[iCen][iptTr]->AddEntry(hdPhiSamNoUE[iCen][1][1][0][0], " ", "lep");
      latDphi[iCen]->DrawLatex(0.15, 0.205, "cluster^{iso}_{wide}");
      legIsoClust_Pi0NoUE[iCen][iptTr]->Draw("same");
      cIsoClust_Pi0NoUE[iCen][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/IsoClust_Pi0NoUE" + sCent + sPtTrig + ".pdf");

      cIsoClust_Pi0NoUERatio[iCen][iptTr]->cd(legPad);
      latDphi[iCen] = LatexDPhi(latDphi[iCen], 0.08, 0.87, cenBins[iCen], cenBins[iCen + 1], ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], true);
      latDphi[iCen]->DrawLatex(0.08, 0.87 - 3 * 0.10, "cluster^{iso}_{narrow}: 0.10 < #it{#sigma}^{2}_{long , 5x5} < 0.30");
      latDphi[iCen]->DrawLatex(0.08, 0.87 - 4 * 0.10, "cluster^{iso}_{wide}: 0.40 < #it{#sigma}^{2}_{long , 5x5} < 1.00");
      legIsoClust_Pi0NoUERatio[iCen][iptTr]->AddEntry(hdPhiSamNoUERatio[iCen][1][0][0], " ", "lep");
      latDphi[iCen]->DrawLatex(0.15, 0.305, "cluster^{iso}_{narrow}/cluster^{iso}_{wide}");
      legIsoClust_Pi0NoUERatio[iCen][iptTr]->Draw("same");
      cIsoClust_Pi0NoUERatio[iCen][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/IsoClust_Pi0NoUERatio" + sCent + sPtTrig + ".pdf");

      cIsoClust_Pi0Pur[iCen][iptTr]->cd(legPad);
      latDphi[iCen] = LatexDPhi(latDphi[iCen], 0.08, 0.87, cenBins[iCen], cenBins[iCen + 1], ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], true);
      // latDphi[iCen]->DrawLatex(0.08, 0.87 - 2 * 0.10, Form("%2.0f < #it{p}_{T}^{tr} < %2.0f GeV/#it{c}", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]));
      latDphi[iCen]->DrawLatex(0.08, 0.87 - 3 * 0.10, "cluster^{iso}_{narrow}: 0.10 < #it{#sigma}^{2}_{long , 5x5} < 0.30");
      latDphi[iCen]->DrawLatex(0.08, 0.87 - 4 * 0.10, "cluster^{iso}_{wide}: 0.40 < #it{#sigma}^{2}_{long , 5x5} < 1.00");
      legIsoClust_Pi0Pur[iCen][iptTr]->AddEntry(hdPhiSamNoUE[iCen][1][0][0][0], " ", "lep");
      latDphi[iCen]->DrawLatex(0.15, 0.36, "cluster^{iso}_{narrow}");
      legIsoClust_Pi0Pur[iCen][iptTr]->AddEntry(hdPhiSamPi0Pur[iCen][1][0][0], " ", "lep");
      latDphi[iCen]->DrawLatex(0.15, 0.24, "(1-#it{P}) #upoint cluster^{iso}_{wide}");
      legIsoClust_Pi0Pur[iCen][iptTr]->AddEntry(hdPhiPhoton[iCen][1][0][0], " ", "lep");
      latDphi[iCen]->DrawLatex(0.15, 0.12, "#it{#gamma}^{iso}");
      legIsoClust_Pi0Pur[iCen][iptTr]->Draw("same");
      cIsoClust_Pi0Pur[iCen][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/IsoGamma_SameNoUE_Pi0pur" + sCent + sPtTrig + ".pdf");

      //=========================================================================

      Int_t selZt = 2;
      cSame_Mix_NoUEIsoClustSingleZt[iCen][iptTr]->cd(1);
      // hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->SetMinimum(-5e-3);
      gStyle->SetPadRightMargin(0.01);
      gStyle->SetPadLeftMargin(0.150);
      hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->GetYaxis()->SetRangeUser(-0.04 * hdPhiSam[iCen][1][0][selZt][iptTr]->GetMaximum(), 1.0 * hdPhiSam[iCen][1][0][selZt][iptTr]->GetMaximum());
      TAxis* a = hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->GetXaxis();
      a->SetNdivisions(-505);
      // Replace first and last labels
      a->ChangeLabel(1, -1, -1, -1, -1, -1, "0");
      a->ChangeLabel(-1, -1, -1, -1, -1, -1, "#pi");

      // Replace the intermediate ones (2,3,4,5)
      a->ChangeLabel(2, -1, -1, -1, -1, -1, "#pi/5");
      a->ChangeLabel(3, -1, -1, -1, -1, -1, "2#pi/5");
      a->ChangeLabel(4, -1, -1, -1, -1, -1, "3#pi/5");
      a->ChangeLabel(5, -1, -1, -1, -1, -1, "4#pi/5");
      
      hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->Draw("same");
      TGaxis::SetMaxDigits(3);

      hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->GetXaxis()->SetLabelSize(0.06);
      hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->GetXaxis()->SetTitleSize(0.06);
      hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->GetYaxis()->SetLabelSize(0.06);
      hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->GetYaxis()->SetTitleSize(0.06);
      hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->GetYaxis()->SetLabelOffset(0.01);
      hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->GetYaxis()->SetTitleOffset(1.2);
      hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->GetXaxis()->SetTitleOffset(1.1);

      hdPhiSam[iCen][1][0][selZt][iptTr]->SetMarkerSize(1.5);
      // hdPhiMix[iCen][1][0][selZt][iptTr]->SetMarkerSize(0);
      hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->SetMarkerSize(1.5);

      hdPhiSam[iCen][1][0][selZt][iptTr]->SetLineColor(1);

      hdPhiSam[iCen][1][0][selZt][iptTr]->SetLineWidth(4);
      hdPhiMix[iCen][1][0][selZt][iptTr]->SetLineWidth(4);
      hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->SetLineWidth(4);

      hdPhiSam[iCen][1][0][selZt][iptTr]->Draw("sameE0");
      hdPhiMix[iCen][1][0][selZt][iptTr]->Draw("same");
      hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->SetTitle("");
      
      lineX0->Draw("same");

      //      cSame_Mix_NoUEIsoClustSingleZt[iCen][iptTr]->cd(2);
      //      // hdPhiSamNoUE[iCen][1][1][selZt][iptTr]->SetMinimum(-5e-3);
      //
      //      // Make both plots equivalent markers
      //      hdPhiSamNoUE[iCen][1][1][selZt][iptTr]->SetMarkerStyle(24);
      //      hdPhiSam[iCen][1][1][selZt][iptTr]->SetMarkerStyle(20);
      //
      //      hdPhiSamNoUE[iCen][1][1][selZt][iptTr]->GetXaxis()->SetLabelSize(0.06);
      //      hdPhiSamNoUE[iCen][1][1][selZt][iptTr]->GetXaxis()->SetTitleSize(0.06);
      //      hdPhiSamNoUE[iCen][1][1][selZt][iptTr]->GetYaxis()->SetLabelSize(0.06);
      //      hdPhiSamNoUE[iCen][1][1][selZt][iptTr]->GetYaxis()->SetTitleSize(0.06);
      //      hdPhiSamNoUE[iCen][1][1][selZt][iptTr]->GetYaxis()->SetLabelOffset(0.01);
      //      hdPhiSamNoUE[iCen][1][1][selZt][iptTr]->GetYaxis()->SetTitleOffset(1.25);
      //      hdPhiSamNoUE[iCen][1][1][selZt][iptTr]->GetXaxis()->SetTitleOffset(1.1);
      //
      //      // Make sure same axis as narrow
      //      hdPhiSamNoUE[iCen][1][1][selZt][iptTr]->GetYaxis()->SetRangeUser(-0.2 * hdPhiSam[iCen][1][0][selZt][iptTr]->GetMaximum(), 1.2 * hdPhiSam[iCen][1][0][selZt][iptTr]->GetMaximum());
      //
      ////      hdPhiSamNoUE[iCen][1][1][selZt][iptTr]->GetYaxis()->SetRangeUser(-0.45, 0.7);
      ////      hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->GetYaxis()->SetRangeUser(-0.45, 0.7);
      //
      //      hdPhiSamNoUE[iCen][1][1][selZt][iptTr]->Draw("same");
      //      hdPhiSam[iCen][1][1][selZt][iptTr]->Draw("sameE0");
      //      hdPhiMix[iCen][1][1][selZt][iptTr]->Draw("same");
      //      hdPhiSamNoUE[iCen][1][1][selZt][iptTr]->SetTitle("Wide: #bf{0.4 < #it{#sigma}^{2}_{long, 5x5} < 1}");
      //
      //      lineX0->Draw("same");

      cSame_Mix_NoUEIsoClustSingleZt[iCen][iptTr]->cd(2);

      TLegend *legZt = new TLegend(0.05, 0.4, 0.05, 0.95);
      legZt->SetTextSize(0.06);
      legZt->SetLineColorAlpha(0, 0);
      legZt->SetFillColorAlpha(0, 0);
      legZt->SetHeader(Form("ALICE, #bf{%d#font[122]{-}%d%%} Pb#font[122]{-}Pb, #sqrt{s_{NN}} = 5.02 TeV", cenBins[iCen], cenBins[iCen + 1]));
      legZt->AddEntry("", Form("Trigger: %1.0f <#color[0]{..}#it{p}_{T}^{trig} < %1.0f GeV/#it{c}, |#it{#eta}^{trig}| < 0.67 ", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]), "");
      legZt->AddEntry("", "0.1 <#color[0]{..}#it{#sigma}^{2}_{long, 5x5} < 0.3 ; #it{p}_{T}^{iso, ch} < 1.5 GeV/#it{c}, #it{R} = 0.2", "");
      // legZt->AddEntry("", "           #it{p}_{T}^{iso, ch} < 1.5 GeV/#it{c}, #it{R} = 0.2", "");

      legZt->AddEntry("", "Associated:#color[0]{..}#it{p}_{T}^{h} > 0.5 GeV/#it{c}, |#it{#eta}^{h}| < 0.9", "");
      legZt->AddEntry("", Form("%1.2f <#color[0]{..}#it{z}_{T} < %1.2f", assocZt[selZt], assocZt[selZt + 1]), "");
      legZt->Draw();

      TLegend *legZtM = new TLegend(0.05, 0.15, 0.5, 0.4);
      legZtM->SetTextSize(0.06);
      legZtM->SetLineColorAlpha(0, 0);
      legZtM->SetFillColorAlpha(0, 0);
      legZtM->AddEntry(hdPhiSam[iCen][1][0][selZt][0], "Same event", "lp");
      legZtM->AddEntry(hdPhiMix[iCen][1][0][selZt][0], "Mixed event", "lp");
      legZtM->AddEntry(hdPhiSamNoUE[iCen][1][0][selZt][0], "Same event - Mixed event", "lp");
      legZtM->Draw();

      TString zTTitle = Form("_Zt%2.2f_%2.2f", assocZt[selZt], assocZt[selZt + 1]);
      cSame_Mix_NoUEIsoClustSingleZt[iCen][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/Same_MixIsoClustGamma_NoUE" + sCent + sPtTrig + zTTitle + ".pdf");
      //=========================================================================

      cIsoClust_Pi0Pur_SingleZt[iCen][iptTr]->cd(1);
      hdPhiPhoton[iCen][1][selZt][iptTr]->SetTitle(""); //(sTitle);
      // hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->SetTitle(sTitle);
      //  hdPhiPhoton[iCen][izt][iptTr]->GetYaxis()->SetRangeUser(-1*hdPhiSamPi0Pur[iCen][izt][iptTr]->GetMaximum(), hdPhiSamPi0Pur[iCen][izt][iptTr]->GetMaximum());
      //   cout << "MAXXXXX: " << hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetMaximum() << endl;

      // cout << (hdPhiPhoton[iCen][izt][iptTr]->GetMinimum() - hdPhiPhoton[iCen][izt][iptTr]->GetBinError(5)) << endl;
      // hdPhiSamNoUE[iCen][1][0][izt][iptTr]->SetMinimum(0.1 * (hdPhiPhoton[iCen][izt][iptTr]->GetMinimum() - hdPhiPhoton[iCen][izt][iptTr]->GetBinError(6)));
      if (iCen == 0 && selZt == 1)
      {
        hdPhiPhoton[iCen][1][selZt][iptTr]->GetYaxis()->SetRangeUser(-60 * 1e-3, 65 * 1e-3);
      }

      // if (iCen == 2 && selZt == 1)
      //{
      //   hdPhiPhoton[iCen][1][selZt][iptTr]->GetYaxis()->SetRangeUser(-50 * 1e-3, 100*1e-3);
      // }

      hdPhiPhoton[iCen][1][selZt][iptTr]->GetXaxis()->SetLabelSize(0.06);
      hdPhiPhoton[iCen][1][selZt][iptTr]->GetXaxis()->SetTitleSize(0.06);
      hdPhiPhoton[iCen][1][selZt][iptTr]->GetYaxis()->SetLabelSize(0.06);
      hdPhiPhoton[iCen][1][selZt][iptTr]->GetYaxis()->SetTitleSize(0.06);
      hdPhiPhoton[iCen][1][selZt][iptTr]->GetYaxis()->SetLabelOffset(0.01);
      hdPhiPhoton[iCen][1][selZt][iptTr]->GetYaxis()->SetTitleOffset(1.2);
      hdPhiPhoton[iCen][1][selZt][iptTr]->GetXaxis()->SetTitleOffset(1.1);
      // TGaxis::SetMaxDigits(2);

      hdPhiSamPi0Pur[iCen][1][selZt][iptTr]->SetMarkerSize(1.5);
      // hdPhiSamPi0Pur[iCen][1][selZt][iptTr]->SetLineColor(kViolet-2);
      // hdPhiSamPi0Pur[iCen][1][selZt][iptTr]->SetMarkerColor(kViolet-2);
      //
      hdPhiSamPi0Pur[iCen][1][selZt][iptTr]->SetLineWidth(4);
      hdPhiPhoton[iCen][1][selZt][iptTr]->SetLineWidth(4);

      hdPhiPhoton[iCen][1][selZt][iptTr]->SetMarkerSize(1.5);
      hdPhiPhoton[iCen][1][selZt][iptTr]->Draw("same");
      hdPhiSamPi0Pur[iCen][1][selZt][iptTr]->Draw("same");
      // hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->SetMarkerStyle(24);
      hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->SetMarkerSize(1.5);
      hdPhiSamNoUE[iCen][1][0][selZt][iptTr]->Draw("same");
      lineX0->Draw("same");
      cIsoClust_Pi0Pur_SingleZt[iCen][iptTr]->cd(2);

      TLegend *legZt2 = new TLegend(0.05, 0.5, 0.05, 0.95);
      legZt2->SetTextSize(0.06);
      legZt2->SetLineColorAlpha(0, 0);
      legZt2->SetFillColorAlpha(0, 0);
      legZt2->SetHeader(Form("ALICE, #bf{%d#font[122]{-}%d%%} Pb#font[122]{-}Pb, #sqrt{s_{NN}} = 5.02 TeV", cenBins[iCen], cenBins[iCen + 1]));
      legZt2->AddEntry("", Form("Trigger: %1.0f#color[0]{.}<#color[0]{.}#it{p}_{T}^{trig} < %1.0f GeV/#it{c}, |#it{#eta}^{trig}| < 0.67, |#Delta#it{#varphi}| >#color[0]{..}#frac{3}{5}#color[0]{.}#it{#pi}", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]), "");
      // legZt2->AddEntry("","           Narrow: 0.1 < #it{#sigma}^{2}_{long, 5x5} < 0.3","");
      legZt2->AddEntry("", "#it{p}_{T}^{iso, ch} < 1.5 GeV/#it{c}, #it{R} = 0.2", "");

      legZt2->AddEntry("", "Associated:#color[0]{..}#it{p}_{T}^{h} > 0.5 GeV/#it{c}, |#it{#eta}^{h}| < 0.9", "");
      legZt2->AddEntry("", Form("%1.2f <#color[0]{..}#it{z}_{T} < %1.2f", assocZt[selZt], assocZt[selZt + 1]), "");
      legZt2->Draw();
      // legZt->Draw();

      TLegend *legZtP = new TLegend(0.05, 0.10, 0.5, 0.45);
      legZtP->SetTextSize(0.06);
      legZtP->SetLineColorAlpha(0, 0);
      legZtP->AddEntry(hdPhiSamNoUE[iCen][1][0][selZt][0], "#it{C}_{Narrow}", "lp");
      legZtP->AddEntry(hdPhiSamPi0Pur[iCen][1][selZt][0], "#it{C}_{Wide}#color[0]{.}#times#color[0]{.}(1#color[0]{.}#font[122]{-}#color[0]{.}#it{P})", "lp");
      // legZtP->AddEntry(hdPhiSamNoUE[iCen][1][0][selZt][0], "Same - Mixed event, Narrow 0.1 < #it{#sigma}^{2}_{long, 5x5} < 0.3", "lep");
      // legZtP->AddEntry(hdPhiSamPi0Pur[iCen][1][selZt][0] , "Same - Mixed event, Wide #times (1 #font[122]{-} #it{P}) 0.4 < #it{#sigma}^{2}_{long, 5x5} < 1.0", "lep");
      legZtP->AddEntry(hdPhiPhoton[iCen][1][selZt][0], "#it{C}_{#it{#gamma}^{iso}}", "lp");
      legZtP->Draw();

      cIsoClust_Pi0Pur_SingleZt[iCen][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/IsoGamma_SameNoUE_Pi0pur" + sCent + sPtTrig + zTTitle + ".pdf");
    }

    cSameNoUE_pTBin_pTAllIsoClust[iCen] = new TCanvas("cSameNoUE_pTBin_pTAllIsoClust" + sCent, "cSameNoUE_pTBin_pTAllIsoClust" + sCent, xNumPad * 800, 2 * 600);
    cSameNoUE_pTBin_pTAllIsoClust[iCen]->Divide(xNumPad, 2);
    legSameNoUE_pTBin_pTAllIsoClust[iCen] = LegStd(legSameNoUE_pTBin_pTAllIsoClust[iCen], 0.10, 0.11, 0.12, 0.46);
    for (Int_t izt = 0; izt < nZtBin; izt++)
    {
      // TString sTitle = Form("%2.0f < #it{p}_{T}^{tr} < %2.0f GeV/#it{c}, %2.2f < #it{z}_{T}^{as} < %2.2f ", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], assocZt[izt], assocZt[izt + 1]);
      TString sTitle = Form(" %2.2f < #it{z}_{T} < %2.2f ", assocZt[izt], assocZt[izt + 1]);
      // gStyle->SetPadRightMargin(0.02);
      gStyle->SetPadRightMargin(0.018);
      // gStyle->SetPadLeftMargin(0.20);
      gStyle->SetPadLeftMargin(0.15);
      // gStyle->SetPadBottomMargin(0.15);
      // gStyle->SetTitleX(0.56);
      cSameNoUE_pTBin_pTAllIsoClust[iCen]->cd(izt + 1);
      hdPhiSamNoUE[iCen][1][0][izt][0]->GetYaxis()->SetRangeUser(0.9 * hdPhiSamNoUE[iCen][1][0][izt][0]->GetMinimum(), 1.15 * hdPhiSamNoUEPtAll[iCen][1][0][izt]->GetMaximum());
      hdPhiSamNoUE[iCen][1][0][izt][0]->Draw("same");
      hdPhiSamNoUEPtAll[iCen][1][0][izt]->Draw("same");
    }
    cSameNoUE_pTBin_pTAllIsoClust[iCen]->cd(legPad);
    legSameNoUE_pTBin_pTAllIsoClust[iCen]->AddEntry(hdPhiSamNoUE[iCen][1][0][0][0], Form("%2.0f < p_{T}^{trig} < %2.0f", ptTrig[index1 + 0], ptTrig[index1 + 1]), "lep");
    legSameNoUE_pTBin_pTAllIsoClust[iCen]->AddEntry(hdPhiSamNoUEPtAll[iCen][1][0][0], Form("%2.0f < p_{T}^{trig} < %2.0f", ptMin, ptMax), "lep");
    legSameNoUE_pTBin_pTAllIsoClust[iCen]->Draw("same");
    cSameNoUE_pTBin_pTAllIsoClust[iCen]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + Form("/SameNoUE_pTBin%2.0f_%2.0f_pTAllIsoClust", ptTrig[index1 + 0], ptTrig[index1 + 1]) + sCent + ".pdf");
  }

  ////////////////////////////////////////////////////////////////////////////////////////////
  //////////////// Combination of figures for paper: DeltaPhi distributions /////////////////
  //////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *cIsoNarrSameMix_SameNoUEXpap[nCen][nPtTrig];
  TCanvas *cIsoNarrWidePurityIsoPhotXpap[nCen][nPtTrig];
  TCanvas *cIsoNarrWidePurityIsoPhotXpapAlternative[nCen][nPtTrig];
  TLegend *legZtTitleXpap[nZtBin];
  TLegend *legNarrSameMix_SameNoUEXpap[nCen][nPtTrig];
  TLegend *legIsoNarrWidePurityIsoPhotXpap[nCen][nPtTrig];
  TLegend *legIsoNarrWidePurityIsoPhotXpapAlternative[nCen][nPtTrig];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    if (iCen == 0 || iCen == 1)
    {
      nZtBin = 6;
      xNumPad = 2;
      yNumPad = 4;
      legPad = 7;
    }
    else if (iCen == 2 || iCen == 3)
    {
      nZtBin = 5;
      xNumPad = 2;
      yNumPad = 3;
      legPad = 6;
    }

    TLegend *legNarrSameMix_SameNoUEXpapVariab = new TLegend(0.05, 0.15, 0.5, 0.4);
    legNarrSameMix_SameNoUEXpapVariab->SetTextSize(0.06);
    legNarrSameMix_SameNoUEXpapVariab->SetLineColorAlpha(0, 0);
    legNarrSameMix_SameNoUEXpapVariab->SetFillColorAlpha(0, 0);
    legNarrSameMix_SameNoUEXpapVariab->AddEntry(hdPhiSam[iCen][1][0][1][0], "Same event", "lp");
    legNarrSameMix_SameNoUEXpapVariab->AddEntry(hdPhiMix[iCen][1][0][1][0], "Mixed event", "lp");
    legNarrSameMix_SameNoUEXpapVariab->AddEntry(hdPhiSamNoUE[iCen][1][0][1][0], "Same event - Mixed event", "lp");

    for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
    {
      TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
      cIsoNarrSameMix_SameNoUEXpap[iCen][iptTr] = new TCanvas("cIsoNarrSameMix_SameNoUEXpap" + sCent + sPtTrig, "cIsoNarrSameMix_SameNoUEXpap" + sCent + sPtTrig, xNumPad * 800, yNumPad * 600);
      cIsoNarrSameMix_SameNoUEXpap[iCen][iptTr]->Divide(xNumPad, yNumPad);
      legNarrSameMix_SameNoUEXpap[iCen][iptTr] = new TLegend(0.05, 0.4, 0.05, 0.95);
      legNarrSameMix_SameNoUEXpap[iCen][iptTr]->SetTextSize(0.06);
      legNarrSameMix_SameNoUEXpap[iCen][iptTr]->SetLineColorAlpha(0, 0);
      legNarrSameMix_SameNoUEXpap[iCen][iptTr]->SetFillColorAlpha(0, 0);
      legNarrSameMix_SameNoUEXpap[iCen][iptTr]->SetHeader(Form("ALICE, #bf{%d#font[122]{-}%d%%} Pb#font[122]{-}Pb, #sqrt{s_{NN}} = 5.02 TeV", cenBins[iCen], cenBins[iCen + 1]));
      legNarrSameMix_SameNoUEXpap[iCen][iptTr]->AddEntry("", Form("Trigger: %1.0f < #it{p}_{T}^{trig} < %1.0f GeV/#it{c}, |#it{#eta}^{trig}| < 0.67", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]), "");
      legNarrSameMix_SameNoUEXpap[iCen][iptTr]->AddEntry("", "0.1 <#color[0]{..}#it{#sigma}^{2}_{long, 5x5}#color[0]{..}< 0.3 ; #it{p}_{T}^{iso, ch} < 1.5 GeV/#it{c}, #it{R} = 0.2", "");
      // legNarrSameMix_SameNoUEXpap[iCen][iptTr]->AddEntry("", "           Isolated: #it{p}_{T}^{iso, ch} < 1.5 GeV/#it{c}, #it{R} = 0.2", "");
      legNarrSameMix_SameNoUEXpap[iCen][iptTr]->AddEntry("", "Associated: #it{p}_{T}^{h} > 0.5 GeV/#it{c}, |#it{#eta}^{h}| < 0.9", "");

      for (Int_t izt = 0; izt < nZtBin; izt++)
      {
        cIsoNarrSameMix_SameNoUEXpap[iCen][iptTr]->cd(izt + 1);
        // TString sTitle = Form("%2.0f < #it{p}_{T}^{tr} < %2.0f GeV/#it{c}, %2.2f < #it{z}_{T}^{as} < %2.2f ", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], assocZt[izt], assocZt[izt + 1]);
        legZtTitleXpap[izt] = LegStd(legZtTitleXpap[izt], 0.4, 0.85, 0.65, 0.90);
        legZtTitleXpap[izt]->SetTextSize(0.07);
        legZtTitleXpap[izt]->SetHeader(Form("%2.2f <#color[0]{..}#it{z}_{T}#color[0]{.}< %2.2f", assocZt[izt], assocZt[izt + 1]), "");
        hdPhiSam[iCen][1][0][izt][iptTr]->SetLineWidth(4);
        hdPhiSam[iCen][1][0][izt][iptTr]->SetMarkerSize(1.3);
        hdPhiMix[iCen][1][0][izt][iptTr]->SetLineWidth(4);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->SetMarkerSize(1.3);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->SetLineWidth(4);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetXaxis()->SetLabelSize(0.06);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetXaxis()->SetTitleSize(0.06);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetLabelSize(0.06);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetTitleSize(0.06);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetLabelOffset(0.01);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetTitleOffset(1.2);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetXaxis()->SetTitleOffset(1.1);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-0.2 * hdPhiSam[iCen][1][0][izt][iptTr]->GetMaximum(), 1.2 * hdPhiSam[iCen][1][0][izt][iptTr]->GetMaximum());
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->Draw("same");
        TGaxis::SetMaxDigits(3);
        hdPhiSam[iCen][1][0][izt][iptTr]->Draw("sameE0");
        hdPhiMix[iCen][1][0][izt][iptTr]->Draw("same");
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->SetTitle("");
        lineX0->Draw("same");
        legZtTitleXpap[izt]->Draw("same");
      }
      cIsoNarrSameMix_SameNoUEXpap[iCen][iptTr]->cd(legPad);
      legNarrSameMix_SameNoUEXpap[iCen][iptTr]->Draw("same");
      legNarrSameMix_SameNoUEXpapVariab->Draw("same");
      cIsoNarrSameMix_SameNoUEXpap[iCen][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/cIsoNarrSameMix_SameNoUEXpap" + sCent + sPtTrig + ".pdf");
    }

    TLegend *legIsoNarrWidePurityIsoPhotXpapVariab = new TLegend(0.05, 0.10, 0.5, 0.45);
    legIsoNarrWidePurityIsoPhotXpapVariab->SetTextSize(0.06);
    legIsoNarrWidePurityIsoPhotXpapVariab->SetLineColorAlpha(0, 0);
    legIsoNarrWidePurityIsoPhotXpapVariab->AddEntry(hdPhiSamNoUE[iCen][1][0][1][0], "#it{C}_{Narrow}", "lp");
    legIsoNarrWidePurityIsoPhotXpapVariab->AddEntry(hdPhiSamPi0Pur[iCen][1][1][0], "#it{C}_{Wide}#color[0]{.}#times#color[0]{.}(1#color[0]{..}#font[122]{-}#color[0]{.}#it{P})", "lp");
    // legIsoNarrWidePurityIsoPhotXpap->AddEntry(hdPhiSamNoUE[iCen][1][0][1][0], "Same - Mixed event, Narrow 0.1 < #it{#sigma}^{2}_{long, 5x5} < 0.3", "lep");
    // legIsoNarrWidePurityIsoPhotXpap->AddEntry(hdPhiSamPi0Pur[iCen][1][1][0] , "Same - Mixed event, Wide #times (1 #font[122]{-} #it{P}) 0.4 < #it{#sigma}^{2}_{long, 5x5} < 1.0", "lep");
    legIsoNarrWidePurityIsoPhotXpapVariab->AddEntry(hdPhiPhoton[iCen][1][1][0], "#it{C}_{#it{#gamma}^{iso}} = 1/#it{P}#color[0]{.}#times#color[0]{.}#it{C}_{Narrow}#color[0]{.}#font[122]{-}#color[0]{.}(1#color[0]{..}#font[122]{-}#color[0]{.}#it{P})/#it{P}#color[0]{.}#times#color[0]{.}#it{C}_{Wide}", "lp");

    for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
    {
      TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
      cIsoNarrWidePurityIsoPhotXpap[iCen][iptTr] = new TCanvas("cIsoNarrWidePurityIsoPhotXpap" + sCent + sPtTrig, "cIsoNarrWidePurityIsoPhotXpap" + sCent + sPtTrig, xNumPad * 800, yNumPad * 600);
      cIsoNarrWidePurityIsoPhotXpap[iCen][iptTr]->Divide(xNumPad, yNumPad);
      legIsoNarrWidePurityIsoPhotXpap[iCen][iptTr] = new TLegend(0.05, 0.45, 0.05, 0.95);
      legIsoNarrWidePurityIsoPhotXpap[iCen][iptTr]->SetTextSize(0.06);
      legIsoNarrWidePurityIsoPhotXpap[iCen][iptTr]->SetLineColorAlpha(0, 0);
      legIsoNarrWidePurityIsoPhotXpap[iCen][iptTr]->SetFillColorAlpha(0, 0);
      legIsoNarrWidePurityIsoPhotXpap[iCen][iptTr]->SetHeader(Form("ALICE, #bf{%d#font[122]{-}%d%%} Pb#font[122]{-}Pb, #sqrt{s_{NN}} = 5.02 TeV", cenBins[iCen], cenBins[iCen + 1]));
      legIsoNarrWidePurityIsoPhotXpap[iCen][iptTr]->AddEntry("", Form("Trigger: %1.0f#color[0]{.}<#color[0]{.}#it{p}_{T}^{trig} < %1.0f GeV/#it{c}, |#it{#eta}^{trig}| < 0.67, |#Delta#it{#varphi}| >#color[0]{..}#frac{3}{5}#color[0]{.}#it{#pi}", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]), "");
      //legIsoNarrWidePurityIsoPhotXpap[iCen][iptTr]->AddEntry("", "", "");
      legIsoNarrWidePurityIsoPhotXpap[iCen][iptTr]->AddEntry("", "#it{p}_{T}^{iso, ch} < 1.5 GeV/#it{c}, #it{R} = 0.2", "");
      legIsoNarrWidePurityIsoPhotXpap[iCen][iptTr]->AddEntry("", "Associated: #it{p}_{T}^{h} > 0.5 GeV/#it{c}, |#it{#eta}^{h}| < 0.9", "");

      for (Int_t izt = 0; izt < nZtBin; izt++)
      {
        cIsoNarrWidePurityIsoPhotXpap[iCen][iptTr]->cd(izt + 1);
        // TString sTitle = Form("%2.0f < #it{p}_{T}^{tr} < %2.0f GeV/#it{c}, %2.2f < #it{z}_{T}^{as} < %2.2f ", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], assocZt[izt], assocZt[izt + 1]);
        legZtTitleXpap[izt] = LegStd(legZtTitleXpap[izt], 0.4, 0.85, 0.65, 0.90);
        legZtTitleXpap[izt]->SetTextSize(0.07);
        legZtTitleXpap[izt]->SetHeader(Form("%2.2f <#color[0]{..}#it{z}_{T}#color[0]{.}< %2.2f", assocZt[izt], assocZt[izt + 1]), "");
        hdPhiSamPi0Pur[iCen][1][izt][iptTr]->SetLineWidth(3);
        hdPhiPhoton[iCen][1][izt][iptTr]->SetLineWidth(3);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->SetLineWidth(3);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetXaxis()->SetLabelSize(0.06);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetXaxis()->SetTitleSize(0.06);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetLabelSize(0.06);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetTitleSize(0.06);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetLabelOffset(0.01);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetTitleOffset(1.2);
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetXaxis()->SetTitleOffset(1.1);
        if(iCen==0 && izt==0 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-0.4, 0.8);
        if(iCen==0 && izt==1 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-0.1, 0.2);
        if(iCen==0 && izt==2 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-0.05, 0.1);
        if(iCen==0 && izt==3 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-15e-3, 20e-3);
        if(iCen==0 && izt==4 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-10e-3, 22e-3);
        if(iCen==0 && izt==5 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-10e-3, 20e-3);
        
        if(iCen==1 && izt==0 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-0.4, 0.6);
        if(iCen==1 && izt==1 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-0.1, 0.2);
        if(iCen==1 && izt==2 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-0.05, 0.1);
        if(iCen==1 && izt==3 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-15e-3, 35e-3);
        if(iCen==1 && izt==4 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-15e-3, 30e-3);
        if(iCen==1 && izt==5 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-15e-3, 20e-3);

        if(iCen==2 && izt==0 && iptTr==0)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-0.04, 0.4);
        if(iCen==2 && izt==1 && iptTr==0)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-0.05, 0.18);
        if(iCen==2 && izt==2 && iptTr==0)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-0.04, 0.18);
        if(iCen==2 && izt==3 && iptTr==0)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-20e-3, 52e-3);
        if(iCen==2 && izt==4 && iptTr==0)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-10e-3, 36e-3);

        if(iCen==2 && izt==0 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-0.04, 0.4);
        if(iCen==2 && izt==1 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-0.05, 0.1);
        if(iCen==2 && izt==2 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-0.04, 0.085);
        if(iCen==2 && izt==3 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-10e-3, 36e-3);
        if(iCen==2 && izt==4 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-10e-3, 26e-3);

        if(iCen==3 && izt==0 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-0.2, 0.4);
        if(iCen==3 && izt==1 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-0.15, 0.2);
        if(iCen==3 && izt==2 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-0.06, 0.250);
        if(iCen==3 && izt==3 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-35e-3, 120e-3);
        if(iCen==3 && izt==4 && iptTr==1)
          hdPhiSamNoUE[iCen][1][0][izt][iptTr]->GetYaxis()->SetRangeUser(-35e-3, 100e-3);
        
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->Draw("same");
        TGaxis::SetMaxDigits(3);
        hdPhiSamPi0Pur[iCen][1][izt][iptTr]->Draw("sameE0");
        hdPhiPhoton[iCen][1][izt][iptTr]->Draw("same");
        hdPhiSamNoUE[iCen][1][0][izt][iptTr]->SetTitle("");
        lineX0->Draw("same");
        legZtTitleXpap[izt]->Draw("same");
      }
      cIsoNarrWidePurityIsoPhotXpap[iCen][iptTr]->cd(legPad);
      legIsoNarrWidePurityIsoPhotXpap[iCen][iptTr]->Draw("same");
      legIsoNarrWidePurityIsoPhotXpapVariab->Draw("same");
      cIsoNarrWidePurityIsoPhotXpap[iCen][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/cIsoNarrWidePurityIsoPhotXpap" + sCent + sPtTrig + ".pdf");
    }
  

    cout<<"ALTERNATIVE WAY OF PLOTTING DeltaPhi Photon"<<endl;
    
    TLegend *legIsoNarrWidePurityIsoPhotXpapVariabAlternative = new TLegend(0.05, 0.10, 0.5, 0.45);
    legIsoNarrWidePurityIsoPhotXpapVariabAlternative->SetTextSize(0.06);
    legIsoNarrWidePurityIsoPhotXpapVariabAlternative->SetLineColorAlpha(0, 0);
    legIsoNarrWidePurityIsoPhotXpapVariabAlternative->AddEntry(hdPhiSamNoUE[iCen][1][0][1][0], "#it{C}_{Narrow}", "lp");
    legIsoNarrWidePurityIsoPhotXpapVariabAlternative->AddEntry(hdPhiSamPi0Pur[iCen][1][1][0], "#it{C}_{Wide}#color[0]{.}#times#color[0]{.}(1#color[0]{..}#font[122]{-}#color[0]{.}#it{P})", "lp");
    legIsoNarrWidePurityIsoPhotXpapVariabAlternative->AddEntry(hdPhiPhoton[iCen][1][1][0], "#it{C}_{#it{#gamma}^{iso}} = 1/#it{P}#color[0]{.}#times#color[0]{.}#it{C}_{Narrow}#color[0]{.}#font[122]{-}#color[0]{.}(1#color[0]{..}#font[122]{-}#color[0]{.}#it{P})/#it{P}#color[0]{.}#times#color[0]{.}#it{C}_{Wide}", "lp");

    for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
      {
        TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
        cIsoNarrWidePurityIsoPhotXpapAlternative[iCen][iptTr] = new TCanvas("cIsoNarrWidePurityIsoPhotXpapAlternative" + sCent + sPtTrig, "cIsoNarrWidePurityIsoPhotXpapAlternative" + sCent + sPtTrig, xNumPad * 800, yNumPad * 600);
        cIsoNarrWidePurityIsoPhotXpapAlternative[iCen][iptTr]->Divide(xNumPad, yNumPad);
        legIsoNarrWidePurityIsoPhotXpapAlternative[iCen][iptTr] = new TLegend(0.05, 0.45, 0.05, 0.95);
        legIsoNarrWidePurityIsoPhotXpapAlternative[iCen][iptTr]->SetTextSize(0.06);
        legIsoNarrWidePurityIsoPhotXpapAlternative[iCen][iptTr]->SetLineColorAlpha(0, 0);
        legIsoNarrWidePurityIsoPhotXpapAlternative[iCen][iptTr]->SetFillColorAlpha(0, 0);
        legIsoNarrWidePurityIsoPhotXpapAlternative[iCen][iptTr]->SetHeader(Form("ALICE, #bf{%d#font[122]{-}%d%%} Pb#font[122]{-}Pb, #sqrt{s_{NN}} = 5.02 TeV", cenBins[iCen], cenBins[iCen + 1]));
        legIsoNarrWidePurityIsoPhotXpapAlternative[iCen][iptTr]->AddEntry("", Form("Trigger: %1.0f#color[0]{.}<#color[0]{.}#it{p}_{T}^{trig} < %1.0f GeV/#it{c}, |#it{#eta}^{trig}| < 0.67, |#Delta#it{#varphi}| >#color[0]{..}#frac{3}{5}#color[0]{.}#it{#pi}", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]), "");
        //legIsoNarrWidePurityIsoPhotXpapAlternative[iCen][iptTr]->AddEntry("", "", "");
        legIsoNarrWidePurityIsoPhotXpapAlternative[iCen][iptTr]->AddEntry("", "#it{p}_{T}^{iso, ch} < 1.5 GeV/#it{c}, #it{R} = 0.2", "");
        legIsoNarrWidePurityIsoPhotXpapAlternative[iCen][iptTr]->AddEntry("", "Associated: #it{p}_{T}^{h} > 0.5 GeV/#it{c}, |#it{#eta}^{h}| < 0.9", "");

        for (Int_t izt = 0; izt < nZtBin; izt++)
        {
          cIsoNarrWidePurityIsoPhotXpapAlternative[iCen][iptTr]->cd(izt + 1);
          // TString sTitle = Form("%2.0f < #it{p}_{T}^{tr} < %2.0f GeV/#it{c}, %2.2f < #it{z}_{T}^{as} < %2.2f ", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], assocZt[izt], assocZt[izt + 1]);
          legZtTitleXpap[izt] = LegStd(legZtTitleXpap[izt], 0.4, 0.85, 0.65, 0.90);
          legZtTitleXpap[izt]->SetTextSize(0.07);
          legZtTitleXpap[izt]->SetHeader(Form("%2.2f <#color[0]{..}#it{z}_{T}#color[0]{.}< %2.2f", assocZt[izt], assocZt[izt + 1]), "");
          hdPhiSamePi0ScaledbyPur[iCen][1][izt][iptTr]->SetLineWidth(3);
          hdPhiPhotonScaledbyPur[iCen][1][izt][iptTr]->SetLineWidth(3);
          hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->SetLineWidth(3);
          hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetXaxis()->SetLabelSize(0.06);
          hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetXaxis()->SetTitleSize(0.06);
          hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetLabelSize(0.06);
          hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetTitleSize(0.06);
          hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetLabelOffset(0.01);
          hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetTitleOffset(1.2);
          hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetXaxis()->SetTitleOffset(1.1);
          if(iCen==0 && izt==0 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-0.4, 0.8);
          if(iCen==0 && izt==1 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-0.1, 0.2);
          if(iCen==0 && izt==2 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-0.05, 0.1);
          if(iCen==0 && izt==3 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-15e-3, 20e-3);
          if(iCen==0 && izt==4 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-10e-3, 22e-3);
          if(iCen==0 && izt==5 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-10e-3, 20e-3);
          
          if(iCen==1 && izt==0 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-0.4, 0.6);
          if(iCen==1 && izt==1 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-0.1, 0.2);
          if(iCen==1 && izt==2 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-0.05, 0.1);
          if(iCen==1 && izt==3 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-15e-3, 35e-3);
          if(iCen==1 && izt==4 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-15e-3, 30e-3);
          if(iCen==1 && izt==5 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-15e-3, 20e-3);

          if(iCen==2 && izt==0 && iptTr==0)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-0.04, 0.4);
          if(iCen==2 && izt==1 && iptTr==0)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-0.05, 0.18);
          if(iCen==2 && izt==2 && iptTr==0)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-0.04, 0.18);
          if(iCen==2 && izt==3 && iptTr==0)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-20e-3, 52e-3);
          if(iCen==2 && izt==4 && iptTr==0)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-10e-3, 36e-3);

          if(iCen==2 && izt==0 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-0.04, 0.4);
          if(iCen==2 && izt==1 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-0.05, 0.1);
          if(iCen==2 && izt==2 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-0.04, 0.085);
          if(iCen==2 && izt==3 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-10e-3, 36e-3);
          if(iCen==2 && izt==4 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-10e-3, 26e-3);

          if(iCen==3 && izt==0 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-0.2, 0.4);
          if(iCen==3 && izt==1 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-0.15, 0.2);
          if(iCen==3 && izt==2 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-0.06, 0.250);
          if(iCen==3 && izt==3 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-35e-3, 120e-3);
          if(iCen==3 && izt==4 && iptTr==1)
            hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->GetYaxis()->SetRangeUser(-35e-3, 100e-3);
          
          hdPhiSameNoUEScaledbyPur[iCen][1][izt][iptTr]->Draw("same");
          TGaxis::SetMaxDigits(3);
          hdPhiSamePi0ScaledbyPur[iCen][1][izt][iptTr]->Draw("sameE0");
          hdPhiPhotonScaledbyPur[iCen][1][izt][iptTr]->Draw("same");
          lineX0->Draw("same");
          legZtTitleXpap[izt]->Draw("same");
        }
        cIsoNarrWidePurityIsoPhotXpapAlternative[iCen][iptTr]->cd(legPad);
        legIsoNarrWidePurityIsoPhotXpapAlternative[iCen][iptTr]->Draw("same");
        legIsoNarrWidePurityIsoPhotXpapVariabAlternative->Draw("same");
        cIsoNarrWidePurityIsoPhotXpapAlternative[iCen][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/cIsoNarrWidePurityIsoPhotXpapAlternative" + sCent + sPtTrig + ".pdf");
      }
  }
  
  
  cout<<"---------------------------------------------"<<endl;

  cout << "Zt distributions Data" << endl;
  TCanvas *cZtPtBin[nCen][nIso][nPtTrig];
  TLatex *latexZtPtBin[nCen][nIso][nPtTrig];
  TCanvas *cZtGamma_Pi0PtBin[nCen][nIso][nPtTrig];
  TCanvas *cZtGamma_Pi0[nCen][nIso];
  TLatex *latexGamma_Pi0[nCen][nIso];

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    for (int iso = 0; iso < nIso; iso++)
    {
      TString sIso = Form("Iso%d", iso);
      cZtGamma_Pi0[iCen][iso] = new TCanvas("cZtGamma_Pi0" + sIso + sCent, "cZtGamma_Pi0" + sIso + sCent, 800, 600);
      cZtGamma_Pi0[iCen][iso]->cd();
      gPad->SetLogy();
      hZtPhoton[iCen][iso]->SetTitle("");
      hZtPhoton[iCen][iso]->GetXaxis()->SetTitleSize(0.040);
      hZtPhoton[iCen][iso]->GetYaxis()->SetTitleSize(0.040);
      hZtPhoton[iCen][iso]->GetXaxis()->SetLabelSize(0.030);
      hZtPhoton[iCen][iso]->GetYaxis()->SetLabelSize(0.030);
      hZtPhoton[iCen][iso]->Draw();
      // hZtPi0[iCen][iso]->Draw("same");
      LatexStd(latexGamma_Pi0[iCen][iso], 0.50, 0.84, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax, false);
      cZtGamma_Pi0[iCen][iso]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/ZtDistribution_Gamma_Pi0" + sIso + sCent + ".pdf");
      // hZtPhoton[iCen][1]->Draw();
      // hZtPi0[iCen][1]->Draw("same");
      cout << "pippo" << endl;
      for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
      {
        TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
        cZtPtBin[iCen][iso][iptTr] = new TCanvas("cZtGamma" + sIso + sCent + sPtTrig, "cZtGamma" + sIso + sCent + sPtTrig, 800, 600);
        cZtPtBin[iCen][iso][iptTr]->cd();
        gPad->SetLogy();
        hZtPhotonPtBin[iCen][iso][iptTr]->SetTitle("");
        hZtPhotonPtBin[iCen][iso][iptTr]->GetXaxis()->SetTitleSize(0.040);
        hZtPhotonPtBin[iCen][iso][iptTr]->GetYaxis()->SetTitleSize(0.040);
        hZtPhotonPtBin[iCen][iso][iptTr]->GetXaxis()->SetLabelSize(0.030);
        hZtPhotonPtBin[iCen][iso][iptTr]->GetYaxis()->SetLabelSize(0.030);
        hZtPhotonPtBin[iCen][iso][iptTr]->Draw();
        LatexStd(latexZtPtBin[iCen][iso][iptTr], 0.40, 0.8, cenBins[iCen], cenBins[iCen + 1], ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], true);
        cZtPtBin[iCen][iso][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/ZtDistribution_DataPurOnly" + sIso + sCent + sPtTrig + ".pdf");

        cZtGamma_Pi0PtBin[iCen][iso][iptTr] = new TCanvas("cZtGamma_Pi0" + sIso + sCent + sPtTrig, "cZtGamma_Pi0" + sIso + sCent + sPtTrig, 800, 600);
        cZtGamma_Pi0PtBin[iCen][iso][iptTr]->cd();
        gPad->SetLogy();
        hZtPhotonPtBin[iCen][iso][iptTr]->SetTitle("");
        hZtPhotonPtBin[iCen][iso][iptTr]->GetXaxis()->SetTitleSize(0.040);
        hZtPhotonPtBin[iCen][iso][iptTr]->GetYaxis()->SetTitleSize(0.040);
        hZtPhotonPtBin[iCen][iso][iptTr]->GetXaxis()->SetLabelSize(0.030);
        hZtPhotonPtBin[iCen][iso][iptTr]->GetYaxis()->SetLabelSize(0.030);
        hZtPhotonPtBin[iCen][iso][iptTr]->Draw();
        // hZtPi0PtBin[iCen][iso][iptTr]->Draw("same");

        LatexStd(latexZtPtBin[iCen][iso][iptTr], 0.40, 0.8, cenBins[iCen], cenBins[iCen + 1], ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], true);
        cZtGamma_Pi0PtBin[iCen][iso][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/ZtDistribution_Gamma_Pi0" + sIso + sCent + sPtTrig + ".pdf");
      }
    }
  }

  cout << "Azimuthal distributions MonteCarlo" << endl;
  TCanvas *cSame_IsoGen[nCen][nPtTrig];
  TLegend *legSame_IsoGen[nCen][nPtTrig];

  TCanvas *cMixIso_NotIsoRec[nCen][nPtTrig];
  TLegend *legMixIso_NotIsoRec[nCen][nPtTrig];

  TCanvas *cSame_MixIsoRec[nCen][nPtTrig];
  TLegend *legSame_MixIsoRec[nCen][nPtTrig];

  TCanvas *cSame_Mix_NoUEIsoRec[nCen][nPtTrig];
  TLegend *legSame_Mix_NoUEIsoRec[nCen][nPtTrig];
  TLatex *latDphiMC[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
    {
      TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);

      cSame_IsoGen[iCen][iptTr] = new TCanvas("cSame_IsoGen" + sCent + sPtTrig, "cSame_IsoGen" + sCent + sPtTrig, xNumPad * 800, 2 * 600);
      cSame_IsoGen[iCen][iptTr]->Divide(xNumPad, 2);
      legSame_IsoGen[iCen][iptTr] = LegStd(legSame_IsoGen[iCen][iptTr], 0.10, 0.275, 0.12, 0.51);

      cMixIso_NotIsoRec[iCen][iptTr] = new TCanvas("cMixIso_NotIsoRec" + sCent + sPtTrig, "cMixIso_NotIsoRec" + sCent + sPtTrig, xNumPad * 800, 2 * 600);
      cMixIso_NotIsoRec[iCen][iptTr]->Divide(xNumPad, 2);
      legMixIso_NotIsoRec[iCen][iptTr] = LegStd(legMixIso_NotIsoRec[iCen][iptTr], 0.10, 0.275, 0.12, 0.51);

      cSame_MixIsoRec[iCen][iptTr] = new TCanvas("cSame_MixIsoRec" + sCent + sPtTrig, "cSame_MixIsoRec" + sCent + sPtTrig, xNumPad * 800, 2 * 600);
      cSame_MixIsoRec[iCen][iptTr]->Divide(xNumPad, 2);
      legSame_MixIsoRec[iCen][iptTr] = LegStd(legSame_MixIsoRec[iCen][iptTr], 0.10, 0.275, 0.12, 0.51);

      cSame_Mix_NoUEIsoRec[iCen][iptTr] = new TCanvas("cSame_Mix_NoUEIsoRec" + sCent + sPtTrig, "cSame_Mix_NoUEIsoRec" + sCent + sPtTrig, xNumPad * 800, 2 * 600);
      cSame_Mix_NoUEIsoRec[iCen][iptTr]->Divide(xNumPad, 2);
      legSame_Mix_NoUEIsoRec[iCen][iptTr] = LegStd(legSame_Mix_NoUEIsoRec[iCen][iptTr], 0.10, 0.10, 0.12, 0.465);
      for (Int_t izt = 0; izt < nZtBin; izt++)
      {
        TString sTitle = Form("%2.0f < #it{p}_{T}^{tr} < %2.0f GeV/#it{c}, %2.2f < #it{z}_{T}^{as} < %2.2f ", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], assocZt[izt], assocZt[izt + 1]);
        gStyle->SetPadRightMargin(0.018);
        // gStyle->SetPadLeftMargin(0.12);
        TGaxis::SetMaxDigits(2);

        cSame_IsoGen[iCen][iptTr]->cd(izt + 1)->SetGridx();
        hdPhiMCGen[iCen][1][0][izt][iptTr]->SetTitle(sTitle);
        hdPhiMCGen[iCen][1][0][izt][iptTr]->Draw("same");
        hdPhiMCGenUE[iCen][1][0][izt][iptTr]->Draw("same");

        cMixIso_NotIsoRec[iCen][iptTr]->cd(izt + 1)->SetGridx();
        hdPhiMixMCRec[iCen][1][0][izt][iptTr]->SetTitle(sTitle);
        hdPhiMixMCRec[iCen][0][0][izt][iptTr]->SetMinimum(0.8 * hdPhiMixMCRec[iCen][1][0][izt][iptTr]->GetMinimum());
        hdPhiMixMCRec[iCen][0][0][izt][iptTr]->Draw("same");
        hdPhiMixMCRec[iCen][1][0][izt][iptTr]->Draw("same");

        cSame_MixIsoRec[iCen][iptTr]->cd(izt + 1)->SetGridx();
        hdPhiSamMCRec[iCen][1][0][izt][iptTr]->SetTitle(sTitle);
        hdPhiSamMCRec[iCen][1][0][izt][iptTr]->SetMinimum(0.8 * hdPhiMixMCRec[iCen][1][0][izt][iptTr]->GetMinimum());
        hdPhiSamMCRec[iCen][1][0][izt][iptTr]->Draw("same");
        hdPhiMixMCRec[iCen][1][0][izt][iptTr]->Draw("same");

        cSame_Mix_NoUEIsoRec[iCen][iptTr]->cd(izt + 1)->SetGridx();
        hdPhiSamMCRecNoUE[iCen][1][0][izt][iptTr]->SetMaximum(1.2 * hdPhiSamMCRec[iCen][1][0][izt][iptTr]->GetMaximum());
        hdPhiSamMCRec[iCen][1][0][izt][iptTr]->SetMinimum(0.5 * hdPhiSamMCRecNoUE[iCen][1][0][izt][iptTr]->GetMinimum());
        hdPhiSamMCRecNoUE[iCen][1][0][izt][iptTr]->SetTitle(sTitle);
        hdPhiSamMCRecNoUE[iCen][1][0][izt][iptTr]->Draw("same");
        hdPhiSamMCRec[iCen][1][0][izt][iptTr]->Draw("same");
        hdPhiMixMCRec[iCen][1][0][izt][iptTr]->Draw("same");
      }

      cSame_IsoGen[iCen][iptTr]->cd(legPad);
      latDphiMC[iCen] = LatexDPhi(latDphiMC[iCen], 0.08, 0.87, cenBins[iCen], cenBins[iCen + 1], ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], true);
      latDphiMC[iCen]->DrawLatex(0.08, 0.87 - 3 * 0.10, "cluster^{iso}_{narrow}: 0.10 < #it{#sigma}^{2}_{long , 5x5} < 0.30");
      legSame_IsoGen[iCen][iptTr]->AddEntry(hdPhiMCGenUE[iCen][1][0][0][0], " ", "lep");
      legSame_IsoGen[iCen][iptTr]->AddEntry(hdPhiMCGen[iCen][1][0][0][0], " ", "lep");
      legSame_IsoGen[iCen][iptTr]->Draw("same");
      latDphiMC[iCen]->DrawLatex(0.15, 0.425, "Same Event Gen");
      latDphiMC[iCen]->DrawLatex(0.15, 0.315, "Same Event Gen w/ UE");
      cSame_IsoGen[iCen][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/MC_SameIsoClustGamma_GEN" + sCent + sPtTrig + ".pdf");

      cSame_MixIsoRec[iCen][iptTr]->cd(legPad);
      latDphiMC[iCen] = LatexDPhi(latDphiMC[iCen], 0.08, 0.87, cenBins[iCen], cenBins[iCen + 1], ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], true);

      latDphiMC[iCen]->DrawLatex(0.08, 0.87 - 3 * 0.10, "cluster^{iso}_{narrow}: 0.10 < #it{#sigma}^{2}_{long , 5x5} < 0.30");
      legSame_MixIsoRec[iCen][iptTr]->AddEntry(hdPhiSamMCRec[iCen][1][0][0][0], " ", "lep");
      latDphiMC[iCen]->DrawLatex(0.15, 0.425, "Same Event Rec");
      legSame_MixIsoRec[iCen][iptTr]->AddEntry(hdPhiMixMCRec[iCen][1][0][0][0], " ", "lep");
      latDphiMC[iCen]->DrawLatex(0.15, 0.315, "Mixed Event Rec");
      legSame_MixIsoRec[iCen][iptTr]->Draw("same");
      cSame_MixIsoRec[iCen][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/MC_SameMixIsoClustGamma_REC" + sCent + sPtTrig + ".pdf");

      cMixIso_NotIsoRec[iCen][iptTr]->cd(legPad);
      latDphiMC[iCen] = LatexDPhi(latDphiMC[iCen], 0.08, 0.87, cenBins[iCen], cenBins[iCen + 1], ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], true);
      latDphiMC[iCen]->DrawLatex(0.08, 0.87 - 3 * 0.10, "cluster^{iso}_{narrow}: 0.10 < #it{#sigma}^{2}_{long , 5x5} < 0.30");
      legMixIso_NotIsoRec[iCen][iptTr]->AddEntry(hdPhiMixMCRec[iCen][1][0][0][0], " ", "lep");
      latDphiMC[iCen]->DrawLatex(0.15, 0.425, "Iso Mixed Event Rec");
      legMixIso_NotIsoRec[iCen][iptTr]->AddEntry(hdPhiMixMCRec[iCen][0][0][0][0], " ", "lep");
      latDphiMC[iCen]->DrawLatex(0.15, 0.315, "Not Iso Mixed Event Rec");
      legMixIso_NotIsoRec[iCen][iptTr]->Draw("same");
      cMixIso_NotIsoRec[iCen][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/MC_MixIso_NotIsoClustGamma_REC" + sCent + sPtTrig + ".pdf");

      cSame_Mix_NoUEIsoRec[iCen][iptTr]->cd(legPad);
      latDphiMC[iCen] = LatexDPhi(latDphiMC[iCen], 0.08, 0.87, cenBins[iCen], cenBins[iCen + 1], ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], true);
      latDphiMC[iCen]->DrawLatex(0.08, 0.87 - 3 * 0.10, "cluster^{iso}_{narrow}: 0.10 < #it{#sigma}^{2}_{long , 5x5} < 0.30");
      legSame_Mix_NoUEIsoRec[iCen][iptTr]->AddEntry(hdPhiSamMCRec[iCen][1][0][0][0], " ", "lep");
      latDphiMC[iCen]->DrawLatex(0.15, 0.385, "Same Event Rec");
      legSame_Mix_NoUEIsoRec[iCen][iptTr]->AddEntry(hdPhiMixMCRec[iCen][1][0][0][0], " ", "lep");
      latDphiMC[iCen]->DrawLatex(0.15, 0.266, "Mixed Event Rec");
      legSame_Mix_NoUEIsoRec[iCen][iptTr]->AddEntry(hdPhiSamMCRecNoUE[iCen][1][0][0][0], " ", "lep");
      latDphiMC[iCen]->DrawLatex(0.15, 0.147, "Same Event - Mixed Event Rec");
      legSame_Mix_NoUEIsoRec[iCen][iptTr]->Draw("same");
      cSame_Mix_NoUEIsoRec[iCen][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/MC_Same_SameNoUE_Mix_REC" + sCent + sPtTrig + ".pdf");
    }
  }

  cout << "Zt distributions MC" << endl;
  TCanvas *cZtMC_GenRec[nCen][nIso][nPtTrig];
  TLatex *latexMC[nCen][nIso][nPtTrig];
  TCanvas *cZtMCRec[nCen][nIso];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    for (int iso = 0; iso < nIso; iso++)
    {
      TString sIso = Form("Iso%d", iso);
      for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
      {
        TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
        cZtMC_GenRec[iCen][iso][iptTr] = new TCanvas("cZtMC_GenRec" + sIso + sCent + sPtTrig, "cZtMC_GenRec" + sIso + sCent + sPtTrig, 800, 600);
        cZtMC_GenRec[iCen][iso][iptTr]->cd();
        for (int iSh = 0; iSh < nShShMC; iSh++)
        {
          hZtPtBinMCGen[iCen][iso][iSh][iptTr]->SetTitle("");
          hZtPtBinMCGen[iCen][iso][iSh][iptTr]->GetXaxis()->SetTitleSize(0.040);
          hZtPtBinMCGen[iCen][iso][iSh][iptTr]->GetYaxis()->SetTitleSize(0.040);
          hZtPtBinMCGen[iCen][iso][iSh][iptTr]->GetXaxis()->SetLabelSize(0.030);
          hZtPtBinMCGen[iCen][iso][iSh][iptTr]->GetYaxis()->SetLabelSize(0.030);
          hZtPtBinMCGen[iCen][iso][iSh][iptTr]->Draw("same");
          hZtPtBinMCRec[iCen][iso][iSh][iptTr]->Draw("same");
        }

        latexMC[iCen][iso][iptTr] = LatexStd(latexMC[iCen][iso][iptTr], 0.50, 0.84, cenBins[iCen], cenBins[iCen + 1], ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], true);
        cZtMC_GenRec[iCen][iso][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/MC_ZtGen_Rec" + sIso + sCent + sPtTrig + ".pdf");
      }

      cZtMCRec[iCen][iso] = new TCanvas("cZtMCRec" + sIso + sCent, "cZtMCRec" + sIso + sCent, 800, 600);
      cZtMCRec[iCen][iso]->cd();
      for (int iSh = 0; iSh < nShShMC; iSh++)
      {
        for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
        {
          hZtPtBinMCRec[iCen][iso][iSh][iptTr]->Draw("same");
        }
      }
    }
  }
  cout << "pippo" << endl;
  TCanvas *cZtMC_GenRec_alphaCorr[nCen][nIso][nPtTrig];
  TPad *pad1[nCen][nIso][nPtTrig];
  TPad *pad2[nCen][nIso][nPtTrig];
  TLegend *legMCGen_Rec[nCen][nIso][nPtTrig];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    for (int iso = 0; iso < nIso; iso++)
    {
      TString sIso = Form("Iso%d", iso);
      for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
      {
        TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
        cZtMC_GenRec_alphaCorr[iCen][iso][iptTr] = new TCanvas("cZtMC_GenRec_alphaCorr" + sIso + sCent + sPtTrig, "cZtMC_GenRec_alphaCorr" + sIso + sCent + sPtTrig, 800, 600);
        legMCGen_Rec[iCen][iso][iptTr] = LegStd(legMCGen_Rec[iCen][iso][iptTr], 0.55, 0.50, 0.84, 0.65);
        cZtMC_GenRec_alphaCorr[iCen][iso][iptTr]->cd();
        pad1[iCen][iso][iptTr] = new TPad("pad1" + sIso + sCent + sPtTrig, "pad1" + sIso + sCent + sPtTrig, 0, 0.3, 1, 1.0);
        pad1[iCen][iso][iptTr]->SetBottomMargin(0); // Upper and lower plot are joined
        pad1[iCen][iso][iptTr]->SetLeftMargin(0.1);
        // pad1[iCen]->SetGridx();         // Vertical grid
        pad1[iCen][iso][iptTr]->Draw(); // Draw the upper pad: pad1
        pad1[iCen][iso][iptTr]->cd();
        gPad->SetLogy();
        for (int iSh = 0; iSh < nShShMC; iSh++)
        {
          hZtPtBinMCGen[iCen][iso][iSh][iptTr]->Draw("same");
          hZtPtBinMCRec[iCen][iso][iSh][iptTr]->Draw("same");
          legMCGen_Rec[iCen][iso][iptTr]->AddEntry(hZtPtBinMCGen[iCen][iso][iSh][iptTr], Form("PYTHIA8 Gen %s", sGenPartMC[iSh].Data()), "lep");
          legMCGen_Rec[iCen][iso][iptTr]->AddEntry(hZtPtBinMCRec[iCen][iso][iSh][iptTr], Form("PYTHIA8 Rec %s", sGenPartMC[iSh].Data()), "lep");
        }

        legMCGen_Rec[iCen][iso][iptTr]->Draw("same");
        latexMC[iCen][iso][iptTr] = LatexStd(latexMC[iCen][iso][iptTr], 0.50, 0.84, cenBins[iCen], cenBins[iCen + 1], ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], true);
        cZtMC_GenRec_alphaCorr[iCen][iso][iptTr]->cd(); // Go back to the main canvas before defining pad2
        pad2[iCen][iso][iptTr] = new TPad("pad2" + sIso + sCent + sPtTrig, "pad2" + sIso + sCent + sPtTrig, 0, 0.05, 1, 0.3);
        pad2[iCen][iso][iptTr]->SetTopMargin(0);
        pad2[iCen][iso][iptTr]->SetLeftMargin(0.1);
        pad2[iCen][iso][iptTr]->SetBottomMargin(0.25);
        // pad2[iCen][iso] ->SetGridx(); // vertical grid
        pad2[iCen][iso][iptTr]->Draw();
        pad2[iCen][iso][iptTr]->cd();

        for (int iSh = 0; iSh < nShShMC; iSh++)
        {
          hRatioEffCorrPtBin[iCen][iso][iSh][iptTr]->SetMinimum(0.6);
          hRatioEffCorrPtBin[iCen][iso][iSh][iptTr]->SetMaximum(2);
          hRatioEffCorrPtBin[iCen][iso][iSh][iptTr]->Draw("same");
          hRatioEffCorrPtBin[iCen][iso][iSh][iptTr]->SetTitle("");
          hRatioEffCorrPtBin[iCen][iso][iSh][iptTr]->GetYaxis()->SetNdivisions(505);
          hRatioEffCorrPtBin[iCen][iso][iSh][iptTr]->GetXaxis()->SetLabelSize(0.12);
          hRatioEffCorrPtBin[iCen][iso][iSh][iptTr]->GetXaxis()->SetTitleSize(0.12);
          hRatioEffCorrPtBin[iCen][iso][iSh][iptTr]->GetYaxis()->SetLabelSize(0.1);
          hRatioEffCorrPtBin[iCen][iso][iSh][iptTr]->GetYaxis()->SetTitleSize(0.1);
          hRatioEffCorrPtBin[iCen][iso][iSh][iptTr]->GetYaxis()->SetTitle("#alpha_{corr} = #frac{Gen}{Rec}");
          hRatioEffCorrPtBin[iCen][iso][iSh][iptTr]->GetYaxis()->CenterTitle(true);
          hRatioEffCorrPtBin[iCen][iso][iSh][iptTr]->GetYaxis()->SetTitleOffset(0.4);
        }

        cZtMC_GenRec_alphaCorr[iCen][iso][iptTr]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/MC_ZtGen_Rec_alphaCorr" + sIso + sCent + sPtTrig + ".pdf");
      }
    }
  }

  TCanvas *cZtIsoGammaEffCorrGlob_MC[nCen][nIso];
  TCanvas *cZtIsoGammaEffCorrGlob[nCen][nIso];
  TLatex *latexGlob[nCen][nIso];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    for (int iso = 0; iso < nIso; iso++)
    {
      TString sIso = Form("Iso%d", iso);
      cZtIsoGammaEffCorrGlob_MC[iCen][iso] = new TCanvas("cZtIsoGammaEffCorrGlob_MC" + sIso + sCent, "cZtIsoGammaEffCorrGlob_MC" + sIso + sCent, 800, 600);
      cZtIsoGammaEffCorrGlob_MC[iCen][iso]->cd();
      gPad->SetLogy();
      for (int iSh = 0; iSh < nShShMC; iSh++)
      {
        hZtMCGen[iCen][iso][iSh]->SetTitle("");
        hZtMCGen[iCen][iso][iSh]->Draw("same");
        hZtMCGen[iCen][iso][iSh]->GetXaxis()->SetTitleSize(0.040);
        hZtMCGen[iCen][iso][iSh]->GetYaxis()->SetTitleSize(0.040);
        hZtMCGen[iCen][iso][iSh]->GetXaxis()->SetLabelSize(0.030);
        hZtMCGen[iCen][iso][iSh]->GetYaxis()->SetLabelSize(0.030);
        hZtMCRec[iCen][iso][iSh]->Draw("same");
        hZtEffCorr[iCen][iso][iSh]->Draw("same");
      }

      latexGlob[iCen][iso] = LatexStd(latexGlob[iCen][iso], 0.50, 0.84, cenBins[iCen], cenBins[iCen + 1], ptTrig[index1], ptTrig[index2], true);
      cZtIsoGammaEffCorrGlob_MC[iCen][iso]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/ZtDistribution_Data_MC" + sIso + sCent + sPtAll + ".pdf");

      cZtIsoGammaEffCorrGlob[iCen][iso] = new TCanvas("cZtIsoGammaEffCorrGlob" + sIso + sCent, "cZtIsoGammaEffCorrGlob" + sIso + sCent, 800, 600);
      cZtIsoGammaEffCorrGlob[iCen][iso]->cd();
      gPad->SetLogy();

      for (int iSh = 0; iSh < nShShMC; iSh++)
      {
        hZtEffCorr[iCen][iso][iSh]->GetYaxis()->SetRangeUser(1e-3, 50);
        hZtEffCorr[iCen][iso][iSh]->SetTitle("");
        hZtEffCorr[iCen][iso][iSh]->GetXaxis()->SetTitleSize(0.040);
        hZtEffCorr[iCen][iso][iSh]->GetYaxis()->SetTitleSize(0.040);
        hZtEffCorr[iCen][iso][iSh]->GetXaxis()->SetLabelSize(0.030);
        hZtEffCorr[iCen][iso][iSh]->GetYaxis()->SetLabelSize(0.030);
        hZtEffCorr[iCen][iso][iSh]->Draw("same");
      }

      latexGlob[iCen][iso] = LatexStd(latexGlob[iCen][iso], 0.50, 0.84, cenBins[iCen], cenBins[iCen + 1], ptTrig[index1], ptTrig[index2], true);
      cZtIsoGammaEffCorrGlob[iCen][iso]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/ZtDistribution_Data" + sIso + sCent + sPtAll + ".pdf");
    }
  }

// Check over the near side

  cout << "Zt distributions Data Near side" << endl;
  TCanvas *cZtPtBinNearSide[nCen][nIso];
  TLatex *latexZtPtBinNearSide[nCen][nIso][nPtTrig];

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    for (int iso = 0; iso < nIso; iso++)
    {
      TString sIso = Form("Iso%d", iso);
      cZtPtBinNearSide[iCen][iso] = new TCanvas("cZtGammaNearSide" + sIso + sCent, "cZtGammaNearSide" + sIso + sCent, 2*800, 3*600);
      cZtPtBinNearSide[iCen][iso]->Divide(2,3);
      for (Int_t iptTr = 0; iptTr < nPtTrig; iptTr++)
      {
        TString sPtTrig = Form("PtTr%2.0f_%2.0f", ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1]);
        cZtPtBinNearSide[iCen][iso]->cd(iptTr+1);
        hZtPhotonPtBinNearSide[iCen][iso][iptTr]->SetTitle("");
        hZtPhotonPtBinNearSide[iCen][iso][iptTr]->GetXaxis()->SetTitleSize(0.040);
        hZtPhotonPtBinNearSide[iCen][iso][iptTr]->GetYaxis()->SetTitleSize(0.040);
        hZtPhotonPtBinNearSide[iCen][iso][iptTr]->GetXaxis()->SetLabelSize(0.030);
        hZtPhotonPtBinNearSide[iCen][iso][iptTr]->GetYaxis()->SetLabelSize(0.030);
        hZtPhotonPtBinNearSide[iCen][iso][iptTr]->Draw();
        LatexStd(latexZtPtBinNearSide[iCen][iso][iptTr], 0.40, 0.8, cenBins[iCen], cenBins[iCen + 1], ptTrig[index1 + iptTr], ptTrig[index1 + iptTr + 1], true);
      }
      cZtPtBinNearSide[iCen][iso]->Print(outDirPlot + Form("/Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]) + "/ZtDistribution_NEARSIDE" + sIso + sCent + ".pdf");
    }
  }
}
