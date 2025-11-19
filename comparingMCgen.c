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
#include "Plotting.h"

using std::cout;
using std::endl;

int const nCen = 4;
int cenBins[] = {0, 10, 30, 50, 90}; // definition of different centralities

int nIso = 2;        // isolation: not iso = 0; iso = 1
int const nShSh = 2; // shower shape ShSh: Signal = 0.10-0.30; Bkg = x
int nZtBinThin = 10; // all zT intervals
double assocZtThinner[] = {-0.05, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20};

int nZtBin = 6; // larger zT intervals defined for QM
double assocZt[] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 1.00};
int npt = 12; // pT trigger intervals
float ptTrig[] = {10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 50, 60, 80};

TString sShShNCentMix;
bool systNMix = false;

void comparingMCgen(TString outDirPlot = "Output_FigXMCgenComparison", TString shshBkg = "0.40-1.00", TString dirFiles = "Output_checkCodeXppMC")
{

    TString shshString[2] = {"0.10-0.30", shshBkg};
    TString shshStringMC[2] = {"0.10-0.30", "0.10-0.30"};

    gSystem->Exec(Form("mkdir %s", outDirPlot.Data()));

    TString processline = Form(".! mkdir -pv %s", outDirPlot.Data());
    gROOT->ProcessLine(processline.Data());

    cout << "Directory: " << outDirPlot.Data() << endl;

    int xNumPad;
    int legPad;
    TFile *fPlot12_18[nCen];
    TFile *fPlot12_40[nCen];
    TFile *fPlot18_40[nCen];
    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////Definition of histograms to be plotted//////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////

    TH1F *hZtMCGen12_18[nCen];
    TH1F *hZtMCGen12_40[nCen];
    TH1F *hZtMCGen18_40[nCen];

    TH1F *hZtMCGen12_18_Ratio18_40[nCen];
    TH1F *hZtMCGen12_40_Ratio18_40[nCen];
    TH1F *hZtMCGen18_40_Ratio18_40[nCen];

    for (int iCen = 0; iCen < nCen; iCen++)
    {
        TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
        fPlot12_18[iCen] = new TFile(dirFiles + "/fPlot" + shshBkg + sCent + "_Pt12_18" + ".root");
        fPlot12_40[iCen] = new TFile(dirFiles + "/fPlot" + shshBkg + sCent + "_Pt12_40" + ".root");
        fPlot18_40[iCen] = new TFile("Output_checkCode/fPlot" + shshBkg + sCent + "_Pt18_40" + ".root");
        // DataPlot
        cout << "Cen: " << cenBins[iCen] << "-" << cenBins[iCen + 1] << ": Get Data plots + Set plots style" << endl;
        hZtMCGen12_18[iCen] = (TH1F *)fPlot12_18[iCen]->Get("hZtMCGenIso1Photon" + sCent + "_Pt12_18");
        hZtMCGen12_40[iCen] = (TH1F *)fPlot12_40[iCen]->Get("hZtMCGenIso1Photon" + sCent + "_Pt12_40");
        hZtMCGen18_40[iCen] = (TH1F *)fPlot18_40[iCen]->Get("hZtMCGenIso1Photon" + sCent + "_Pt18_40");

        hZtMCGen12_18_Ratio18_40[iCen] = (TH1F *)hZtMCGen12_18[iCen]->Clone("hZtMCGen12_18_Ratio18_40" + sCent);
        hZtMCGen12_40_Ratio18_40[iCen] = (TH1F *)hZtMCGen12_40[iCen]->Clone("hZtMCGen12_40_Ratio18_40" + sCent);
        hZtMCGen18_40_Ratio18_40[iCen] = (TH1F *)hZtMCGen18_40[iCen]->Clone("hZtMCGen18_40_Ratio18_40" + sCent);

        hZtMCGen12_18_Ratio18_40[iCen]->Divide(hZtMCGen18_40[iCen]);
        hZtMCGen12_40_Ratio18_40[iCen]->Divide(hZtMCGen18_40[iCen]);
        hZtMCGen18_40_Ratio18_40[iCen]->Divide(hZtMCGen18_40[iCen]);

        PlotStyle(hZtMCGen12_18[iCen], 25, 1, kBlack, kBlack + 1, "#it{z}_{T}", "1 / #it{N}^{ trig} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
        PlotStyle(hZtMCGen12_40[iCen], 21, 1, kBlue + 1, kBlue + 1, "#it{z}_{T}", "1 / #it{N}^{ trig} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
        PlotStyle(hZtMCGen18_40[iCen], 24, 1, kRed, kRed, "#it{z}_{T}", "1 / #it{N}^{ trig} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
        PlotStyle(hZtMCGen12_18_Ratio18_40[iCen], 25, 1, kBlack, kBlack + 1, "#it{z}_{T}", "Ratio", false);
        PlotStyle(hZtMCGen12_40_Ratio18_40[iCen], 21, 1, kBlue + 1, kBlue + 1, "#it{z}_{T}", "Ratio", false);
        PlotStyle(hZtMCGen18_40_Ratio18_40[iCen], 24, 1, kRed, kRed, "#it{z}_{T}", "Ratio", false);
    }

    TCanvas *cZtGenDist[nCen];
    TCanvas *cZtGenRatioto18_40[nCen];
    TLegend *legZtGenDist[nCen];
    TLegend *legZtGen[nCen];
    TLegend *legZtCond[nCen];
    TLegend *legZtGenRatioto18_40[nCen];
    TH1F *hGeneralRatio = new TH1F("hGeneralRatio", "", 105, 0., 1.05);
    PlotStyle(hGeneralRatio, 20, 1, kWhite, kWhite, " #it{z}_{T} ", "#it{D}(z_{T}) ratio", false);
    TGraph *lineX0 = DrawLine(lineX0, 0, 1, 1.1, 1);
    for (int iCen = 0; iCen < nCen; iCen++)
    {
        TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
        cZtGenDist[iCen] = canvasStdIaa("cZtGenDist" + sCent, 1, 1);
        gPad->SetLogy();
        hZtMCGen12_18[iCen]->GetYaxis()->SetRangeUser(1e-2, 20);
        hZtMCGen12_18[iCen]->SetTitle("");
        hZtMCGen12_40[iCen]->SetTitle("");
        hZtMCGen18_40[iCen]->SetTitle("");
        hZtMCGen12_18[iCen]->Draw("histpesame");
        hZtMCGen12_40[iCen]->Draw("histpesame");
        hZtMCGen18_40[iCen]->Draw("histpesame");

        legZtGen[iCen] = LegStd(legZtGen[iCen], 0.13, 0.75, 0.13, 0.96);
        legZtGen[iCen]->SetHeader("ALICE simulation");
        legZtGen[iCen]->AddEntry("","pp,#color[0]{.}#sqrt{s} = 5.02 TeV","");
        legZtGen[iCen]->AddEntry("","PYTHIA 8,#color[0]{.}#gamma#font[122]{-}jet","");
        legZtGen[iCen]->AddEntry("","#it{p}_{T}^{iso, ch} < 1.5 GeV/#it{c}, #it{R} = 0.2","");
        legZtGen[iCen]->Draw("same");
        legZtGenDist[iCen] = LegStd(legZtGenDist[iCen], 0.54, 0.65, 0.71, 0.85);
        legZtGenDist[iCen]->AddEntry(hZtMCGen12_18[iCen], "#gamma#font[122]{-}jet, 12#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}18 GeV/#it{c}", "lep");
        legZtGenDist[iCen]->AddEntry(hZtMCGen12_40[iCen], "#gamma#font[122]{-}jet, 12#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}", "lep");
        legZtGenDist[iCen]->AddEntry(hZtMCGen18_40[iCen], "#gamma#font[122]{-}jet, 18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}", "lep");
        legZtGenDist[iCen]->Draw("same");
        cZtGenDist[iCen]->Print(outDirPlot + "/cZtGenDist" + sCent + ".pdf");

        cZtGenRatioto18_40[iCen] = canvasStdIaa("cZtGenRatioto18_40" + sCent, 1, 1);
        hZtMCGen12_18_Ratio18_40[iCen]->GetYaxis()->SetRangeUser(0.7, 1.15);
        //hZtMCGen12_18_Ratio18_40[iCen]
        hZtMCGen12_18_Ratio18_40[iCen]->SetTitle("");
        hZtMCGen12_40_Ratio18_40[iCen]->SetTitle("");
        hZtMCGen18_40_Ratio18_40[iCen]->SetTitle("");
        hGeneralRatio->GetYaxis()->SetRangeUser(0.8, 1.1);
        hGeneralRatio->GetYaxis()->SetTitleSize(0.046);
        hGeneralRatio->GetYaxis()->SetTitleOffset(1.1);
        hGeneralRatio->GetYaxis()->SetLabelSize(0.04);
        hGeneralRatio->GetXaxis()->SetRangeUser(0.05, 1.05);
        hGeneralRatio->GetXaxis()->SetNdivisions(505);
        hGeneralRatio->GetXaxis()->SetLabelSize(0.04);
        hGeneralRatio->GetXaxis()->SetTitleSize(0.05);
        hGeneralRatio->Draw("same");
        hZtMCGen12_18_Ratio18_40[iCen]->Draw("histpesame");
        hZtMCGen12_40_Ratio18_40[iCen]->Draw("histpesame");
        //hZtMCGen18_40_Ratio18_40[iCen]->Draw("histpesame");
        legZtGenRatioto18_40[iCen] = LegStd(legZtGenRatioto18_40[iCen], 0.62, 0.18, 0.95, 0.45);
        legZtGenRatioto18_40[iCen]->AddEntry(hZtMCGen12_18_Ratio18_40[iCen], "#frac{12#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}18 GeV/#it{c}}{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}", "lp");
        legZtGenRatioto18_40[iCen]->AddEntry(hZtMCGen12_40_Ratio18_40[iCen], "#frac{12#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}", "lp");
        //legZtGenRatioto18_40[iCen]->AddEntry(hZtMCGen18_40_Ratio18_40[iCen], "#gamma#font[122]{-}jet, #frac{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}", "lp");
        lineX0->Draw("same");
        legZtGenRatioto18_40[iCen]->Draw("same");
        legZtGen[iCen]->Draw("same");
        cZtGenRatioto18_40[iCen]->Print(outDirPlot + "/cZtGenRatioto18_40" + sCent + ".pdf");
    }

    TFile *fileNLO = new TFile("RootFiles/fileNLO.root ");
    TH1F * grDztNLOmedianpp = (TH1F *)fileNLO->Get(Form("grDztNLOmedian_pp"));
    grDztNLOmedianpp->SetLineWidth(8);
    grDztNLOmedianpp->SetLineColor(kRed - 4);
    grDztNLOmedianpp->SetLineStyle(1);

    TFile *finputpp_pPbPaper_DztpPb = new TFile("TheoryCalculations/HEPData-ins1798523-v1-root.root");
    TDirectory *dirPaper_DztpPb = (TDirectory *)finputpp_pPbPaper_DztpPb->Get("Figure 5 Top Panel");

    TH1F *hPYTHIApaper = (TH1F *)dirPaper_DztpPb->Get("Hist1D_y3");
    TH1F *hPYTHIApaper_stat = (TH1F *)dirPaper_DztpPb->Get("Hist1D_y3_e1");
    for (int ibin = 0; ibin < hPYTHIApaper->GetNbinsX(); ibin++)
    {
        hPYTHIApaper->SetBinError(ibin + 1, hPYTHIApaper_stat->GetBinContent(ibin + 1));
    }
    TCanvas *cZtGenDist_PYTHIA = new TCanvas("cZtGenDist_PYTHIA", "cZtGenDist_PYTHIA", 800, 600);
    gPad->SetLogy();
    TH1F* hGeneral = new TH1F("hGeneral", "hGeneral", 100, 0,1);
    //hPYTHIApaper->GetXaxis()->SetRangeUser(0.,1.05);
    //hZtMCGen12_40[2]->GetXaxis()->SetRangeUser(0.,1.05);
    hGeneral->Draw("histpe");
    grDztNLOmedianpp->Draw("histl same");
    hZtMCGen12_40[2]->Draw("histpesame");
    hPYTHIApaper->Draw("histpesame");
    hZtMCGen18_40[2]->Draw("histpesame");
    TLegend *legZtGenDistPYTHIApaper = LegStd(legZtGenDistPYTHIApaper, 0.5, 0.65, 0.85, 0.85);
    legZtGenDistPYTHIApaper->AddEntry(hPYTHIApaper, "PYTHIA pp paper", "lep");
    legZtGenDistPYTHIApaper->AddEntry(hZtMCGen12_40[2], "MC Gen 12_40", "lep");
    legZtGenDistPYTHIApaper->AddEntry(hZtMCGen18_40[2], "MC Gen18_40", "lep");
    legZtGenDistPYTHIApaper->AddEntry(grDztNLOmedianpp, "NLO pQCD", "lep");
    legZtGenDistPYTHIApaper->Draw("same");
    cZtGenDist_PYTHIA->Print(outDirPlot + "/cZtGenDist_PYTHIA.pdf");
}