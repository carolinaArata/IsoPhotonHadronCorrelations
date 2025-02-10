
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

using std::cout;
using std::endl;
// Int_t nCen = 4;
// Int_t cenBins[] = {0, 10, 30, 50, 90};
//  Int_t cenBins[] = { 50, 90};
Int_t kMarkCen[] = {47, 45, 33, 25};
Int_t kColorMark[] = {kCyan + 2, kAzure - 3, kViolet + 6, kCyan - 2};
Int_t nAssoc = 7;

int nZtBinThin = 9;
double assocZtThinner[] = {0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80, 1.00, 1.05};

void PlotStyle(TH1F *hPlot, int kMarker, double kMarkerSize, int kColor, TString titleX, TString titleY);
TLatex *LatexStd(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax, bool bCen);
TLegend *LegStd(TLegend *leg, double xpos1, double ypos1, double xpos2, double ypos2);

void PlotZtVsOtherExp(float ptMin = 18, float ptMax = 40, bool Mirror = true, TString sMixed = "Mixed", TString shshBkg = "0.40-1.00", TString dirPlot = "~/work/histogram/zTfunct_PbPb_ppZtMergedCopy", bool b0_30 = true)
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
    dirSyst = "~/work/histogram/Systematics0_30";
  }
  else if (!b0_30)
  {
    nCen = 4;
    cenBins.push_back(0);
    cenBins.push_back(10);
    cenBins.push_back(30);
    cenBins.push_back(50);
    cenBins.push_back(90);
    dirSyst = "~/work/histogram/Systematics";
  }

  TString shshString[2] = {"0.10-0.30", shshBkg};
  TString sPtAll = Form("_Pt%2.0f_%2.0f", ptMin, ptMax);

  // Getter zt distributions
  TFile *fPlot[nCen];
  TH1F *hZtCent[nCen];    // zt data
  TH1F *hZt_MC_Gen[nCen]; // MC Gen pp
  TH1F *hZt_MC_Rec[nCen]; // MC Rec pp

  for (int iCen = 0; iCen < nCen; iCen++)
  {
    TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
    cout << "Getter zt distributions: " << sCent << endl;
    fPlot[iCen] = new TFile(Form("~/work/histogram/FromScratch/ResultsNewMixMC_ZtMergedMore/fPlot%s%s%s.root", shshString[1].Data(), sCent.Data(), sPtAll.Data()));

    hZt_MC_Gen[iCen] = (TH1F *)fPlot[iCen]->Get(Form("hZtIsoGammaMCGen%s%s", sCent.Data(), sPtAll.Data()));
    hZt_MC_Rec[iCen] = (TH1F *)fPlot[iCen]->Get(Form("hZtIsoGammaMCRec%s%s", sCent.Data(), sPtAll.Data()));
    hZtCent[iCen] = (TH1F *)fPlot[iCen]->Get(Form("hZtEffCorr%s%s", sCent.Data(), sPtAll.Data()));
    PlotStyle(hZt_MC_Gen[iCen], 72, 1, kAzure - 6, "#font[12]{z}_{T}", "1 / #it{N}^{ #gamma} d^{3}#it{N} / d#font[12]{z}_{T}d|#Delta#varphi|d#Delta#eta");
    PlotStyle(hZt_MC_Rec[iCen], 21, 1, kOrange + 7, "#font[12]{z}_{T}", "1 / #it{N}^{ #gamma} d^{3}#it{N} / d#font[12]{z}_{T}d|#Delta#varphi|d#Delta#eta");
    PlotStyle(hZtCent[iCen], kMarkCen[iCen], 1, kColorMark[iCen], "#font[12]{z}_{T}", "1 / #it{N}^{ #gamma} d^{3}#it{N} / d#font[12]{z}_{T}d|#Delta#varphi|d#Delta#eta");
  }
  gSystem->Exec(Form("mkdir %s", dirPlot.Data()));
  cout << "Getter zt systematics : " << endl;
  TFile *fSystFile = new TFile(Form("%s/fSystFile%s%s%s%s.root", dirSyst.Data(), sMixed.Data(), shshBkg.Data(), sMirror.Data(), sPtAll.Data()));
  cout << dirSyst << endl;
  TH1F *hSystZt[nCen];
  for (int iCen = 0; iCen < nCen; iCen++)
  {
    hSystZt[iCen] = (TH1F *)fSystFile->Get(Form("hsystCen%d%d", cenBins[iCen], cenBins[iCen + 1]));
    PlotStyle(hSystZt[iCen], kMarkCen[iCen], 1, kColorMark[iCen], "#font[12]{z}_{T}", "1 / #it{N}^{ #gamma} d^{3}#it{N} / d#font[12]{z}_{T}d|#Delta#varphi|d#Delta#eta");
    for (int ibin = 0; ibin < nAssoc; ibin++)
    {
      cout << cenBins[iCen] << "-" << cenBins[iCen + 1] << endl;
      cout << ibin << "-" << ibin + 1 << " systematics: " << hSystZt[iCen]->GetBinContent(ibin + 1) << " pm " << hSystZt[iCen]->GetBinError(ibin + 1) << endl;
    }
  }

  
}
 