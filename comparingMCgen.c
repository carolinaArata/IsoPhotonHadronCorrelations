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

int const nCen = 1;
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

void comparingMCgen(TString shshBkg = "0.40-1.00", bool bppMC = true)
{

    TString dirFiles[2] = {"ResultsStudyMCppR02", "ResultsStudyMCppR04"};
    TString outDirPlot[2] = {"Output_FigXMCgenComparisonR02", "Output_FigXMCgenComparisonR04"};
    TString shshString[2] = {"0.10-0.30", shshBkg};
    TString shshStringMC[2] = {"0.10-0.30", "0.10-0.30"};
    TString stringIso[2] = {"#it{R} = 0.2", "#it{R} = 0.4"};

    gSystem->Exec(Form("mkdir %s", outDirPlot[0].Data()));
    gSystem->Exec(Form("mkdir %s", outDirPlot[1].Data()));

    TString processline = Form(".! mkdir -pv %s", outDirPlot[0].Data());
    gROOT->ProcessLine(processline.Data());
    processline = Form(".! mkdir -pv %s", outDirPlot[1].Data());
    gROOT->ProcessLine(processline.Data());

    cout << "Directory: " << outDirPlot[0].Data() << endl;
    cout << "Directory: " << outDirPlot[1].Data() << endl;

    int xNumPad;
    int legPad;
    TFile *fPlot12_18[nCen][2];
    TFile *fPlot12_40[nCen][2];
    TFile *fPlot18_40[nCen][2];
    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////Definition of histograms to be plotted//////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////

    TH1F *hZtMCGen12_18[nCen][2];
    TH1F *hZtMCGen12_40[nCen][2];
    TH1F *hZtMCGen18_40[nCen][2];

    TH1F *hZtMCGen12_18_Ratio18_40[nCen][2];
    TH1F *hZtMCGen12_40_Ratio18_40[nCen][2];
    TH1F *hZtMCGen18_40_Ratio18_40[nCen][2];

    for (int iCen = 0; iCen < nCen; iCen++)
    {
        for (int iFile = 0; iFile < 2; iFile++)
        {
            TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
            fPlot12_18[iCen][iFile] = new TFile(dirFiles[iFile] + "/fPlot" + shshBkg + sCent + "_Pt12_18" + ".root");
            fPlot12_40[iCen][iFile] = new TFile(dirFiles[iFile] + "/fPlot" + shshBkg + sCent + "_Pt12_40" + ".root");
            fPlot18_40[iCen][iFile] = new TFile(dirFiles[iFile] + "/fPlot" + shshBkg + sCent + "_Pt18_40" + ".root");
            // DataPlot
            cout << "Cen: " << cenBins[iCen] << "-" << cenBins[iCen + 1] << ": Get Data plots + Set plots style" << endl;
            if (bppMC)
                sCent = ("");
            hZtMCGen12_18[iCen][iFile] = (TH1F *)fPlot12_18[iCen][iFile]->Get("hZtMCGenIso1Photon" + sCent + "_Pt12_18");
            cout << hZtMCGen12_18[iCen][iFile] << endl;
            hZtMCGen12_40[iCen][iFile] = (TH1F *)fPlot12_40[iCen][iFile]->Get("hZtMCGenIso1Photon" + sCent + "_Pt12_40");
            cout << hZtMCGen12_40[iCen][iFile] << endl;
            hZtMCGen18_40[iCen][iFile] = (TH1F *)fPlot18_40[iCen][iFile]->Get("hZtMCGenIso1Photon" + sCent + "_Pt18_40");
            cout << hZtMCGen18_40[iCen][iFile] << endl;

            hZtMCGen12_18_Ratio18_40[iCen][iFile] = (TH1F *)hZtMCGen12_18[iCen][iFile]->Clone("hZtMCGen12_18_Ratio18_40" + sCent);
            hZtMCGen12_40_Ratio18_40[iCen][iFile] = (TH1F *)hZtMCGen12_40[iCen][iFile]->Clone("hZtMCGen12_40_Ratio18_40" + sCent);
            hZtMCGen18_40_Ratio18_40[iCen][iFile] = (TH1F *)hZtMCGen18_40[iCen][iFile]->Clone("hZtMCGen18_40_Ratio18_40" + sCent);

            hZtMCGen12_18_Ratio18_40[iCen][iFile]->Divide(hZtMCGen18_40[iCen][iFile]);
            hZtMCGen12_40_Ratio18_40[iCen][iFile]->Divide(hZtMCGen18_40[iCen][iFile]);
            hZtMCGen18_40_Ratio18_40[iCen][iFile]->Divide(hZtMCGen18_40[iCen][iFile]);

            PlotStyle(hZtMCGen12_18[iCen][iFile], 25, 1, kBlack, kBlack + 1, "#it{z}_{T}", "1 / #it{N}^{ trig} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
            PlotStyle(hZtMCGen12_40[iCen][iFile], 21, 1, kBlue + 1, kBlue + 1, "#it{z}_{T}", "1 / #it{N}^{ trig} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
            PlotStyle(hZtMCGen18_40[iCen][iFile], 24, 1, kRed, kRed, "#it{z}_{T}", "1 / #it{N}^{ trig} d^{3}#it{N} / d#Delta#it{#eta} d|#Delta#it{#varphi}| d#it{z}_{T}", false);
            PlotStyle(hZtMCGen12_18_Ratio18_40[iCen][iFile], 25, 1, kBlack, kBlack + 1, "#it{z}_{T}", "Ratio", false);
            PlotStyle(hZtMCGen12_40_Ratio18_40[iCen][iFile], 21, 1, kBlue + 1, kBlue + 1, "#it{z}_{T}", "Ratio", false);
            PlotStyle(hZtMCGen18_40_Ratio18_40[iCen][iFile], 24, 1, kRed, kRed, "#it{z}_{T}", "Ratio", false);
        }
    }
    TCanvas *cZtGenDist[nCen][2];
    TCanvas *cZtGenRatioto18_40[nCen][2];
    TLegend *legZtGenDist[nCen][2];
    TLegend *legZtGen[nCen][2];
    TLegend *legZtCond[nCen][2];
    TLegend *legZtGenRatioto18_40[nCen][2];
    TH1F *hGeneralRatio = new TH1F("hGeneralRatio", "", 105, 0., 1.05);
    PlotStyle(hGeneralRatio, 20, 1, kWhite, kWhite, " #it{z}_{T} ", "D(z_{T}) ratio", false);
    TGraph *lineX0 = DrawLine(lineX0, 0, 1, 1.1, 1);
    for (int iCen = 0; iCen < nCen; iCen++)
    {
        for (int iFile = 0; iFile < 2; iFile++)
        {
            TString sCent = Form("_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
            cZtGenDist[iCen][iFile] = canvasStdIaa("cZtGenDist" + sCent + iFile, 1, 1);
            gPad->SetLogy();
            hZtMCGen12_18[iCen][iFile]->GetYaxis()->SetRangeUser(1e-2, 20);
            hZtMCGen12_18[iCen][iFile]->SetTitle("");
            hZtMCGen12_40[iCen][iFile]->SetTitle("");
            hZtMCGen18_40[iCen][iFile]->SetTitle("");
            hZtMCGen12_18[iCen][iFile]->Draw("histpesame");
            hZtMCGen12_40[iCen][iFile]->Draw("histpesame");
            hZtMCGen18_40[iCen][iFile]->Draw("histpesame");

            legZtGen[iCen][iFile] = LegStd(legZtGen[iCen][iFile], 0.13, 0.75, 0.13, 0.96);
            legZtGen[iCen][iFile]->SetHeader("ALICE simulation");
            legZtGen[iCen][iFile]->AddEntry("", "pp,#color[0]{.}#sqrt{s} = 5.02 TeV", "");
            legZtGen[iCen][iFile]->AddEntry("", "PYTHIA 8,#color[0]{.}#gamma#font[122]{-}jet", "");
            legZtGen[iCen][iFile]->AddEntry("", "#it{p}_{T}^{iso, ch} < 1.5 GeV/#it{c}, " + stringIso[iFile], "");
            legZtGen[iCen][iFile]->Draw("same");
            legZtGenDist[iCen][iFile] = LegStd(legZtGenDist[iCen][iFile], 0.54, 0.65, 0.71, 0.85);
            legZtGenDist[iCen][iFile]->AddEntry(hZtMCGen12_18[iCen][iFile], "#gamma#font[122]{-}jet, 12#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}18 GeV/#it{c}", "lep");
            legZtGenDist[iCen][iFile]->AddEntry(hZtMCGen12_40[iCen][iFile], "#gamma#font[122]{-}jet, 12#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}", "lep");
            legZtGenDist[iCen][iFile]->AddEntry(hZtMCGen18_40[iCen][iFile], "#gamma#font[122]{-}jet, 18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}", "lep");
            legZtGenDist[iCen][iFile]->Draw("same");
            cZtGenDist[iCen][iFile]->Print(outDirPlot[iFile] + "/cZtGenDist" + sCent + ".pdf");

            cZtGenRatioto18_40[iCen][iFile] = canvasStdIaa("cZtGenRatioto18_40" + sCent + iFile, 1, 1);
            hZtMCGen12_18_Ratio18_40[iCen][iFile]->GetYaxis()->SetRangeUser(0.7, 1.15);
            // hZtMCGen12_18_Ratio18_40[iCen][iFile]
            hZtMCGen12_18_Ratio18_40[iCen][iFile]->SetTitle("");
            hZtMCGen12_40_Ratio18_40[iCen][iFile]->SetTitle("");
            hZtMCGen18_40_Ratio18_40[iCen][iFile]->SetTitle("");
            hGeneralRatio->GetYaxis()->SetRangeUser(0.85, 1.1);
            hGeneralRatio->GetYaxis()->SetTitleSize(0.046);
            hGeneralRatio->GetYaxis()->SetTitleOffset(1.1);
            hGeneralRatio->GetYaxis()->SetLabelSize(0.04);
            hGeneralRatio->GetXaxis()->SetRangeUser(0.05, 1.05);
            hGeneralRatio->GetXaxis()->SetNdivisions(505);
            hGeneralRatio->GetXaxis()->SetLabelSize(0.04);
            hGeneralRatio->GetXaxis()->SetTitleSize(0.05);
            hGeneralRatio->Draw("same");
            hZtMCGen12_18_Ratio18_40[iCen][iFile]->Draw("histpesame");
            hZtMCGen12_40_Ratio18_40[iCen][iFile]->Draw("histpesame");
            // hZtMCGen18_40_Ratio18_40[iCen][iFile]->Draw("histpesame");
            legZtGenRatioto18_40[iCen][iFile] = LegStd(legZtGenRatioto18_40[iCen][iFile], 0.62, 0.18, 0.95, 0.45);
            legZtGenRatioto18_40[iCen][iFile]->AddEntry(hZtMCGen12_18_Ratio18_40[iCen][iFile], "#frac{12#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}18 GeV/#it{c}}{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}", "lp");
            legZtGenRatioto18_40[iCen][iFile]->AddEntry(hZtMCGen12_40_Ratio18_40[iCen][iFile], "#frac{12#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}", "lp");
            // legZtGenRatioto18_40[iCen][iFile]->AddEntry(hZtMCGen18_40_Ratio18_40[iCen][iFile], "#gamma#font[122]{-}jet, #frac{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}", "lp");
            lineX0->Draw("same");
            legZtGenRatioto18_40[iCen][iFile]->Draw("same");
            legZtGen[iCen][iFile]->Draw("same");
            cZtGenRatioto18_40[iCen][iFile]->Print(outDirPlot[iFile] + "/cZtGenRatioto18_40" + sCent + ".pdf");
        }
    }

    TCanvas *cZtGenRatioPbPb_pp = canvasStdIaa("cZtGenRatioPbPb_pp", 1, 1);
    TLegend *legZtGenRatioPbPb_pp = LegStd(legZtGenRatioPbPb_pp, 0.45, 0.18, 0.86, 0.45);
    TH1F* hRatioPbPb_pp = (TH1F*)hZtMCGen12_40[0][0]->Clone("hRatioPbPb_pp");
    hRatioPbPb_pp->Divide(hZtMCGen18_40[0][0]);
    hRatioPbPb_pp->GetYaxis()->SetRangeUser(0.7, 1.15);
    // hRatioPbPb_pp
    hRatioPbPb_pp->SetTitle("");
    hGeneralRatio->GetYaxis()->SetRangeUser(0.85, 1.1);
    hGeneralRatio->GetYaxis()->SetTitleSize(0.046);
    hGeneralRatio->GetYaxis()->SetTitleOffset(1.1);
    hGeneralRatio->GetYaxis()->SetLabelSize(0.04);
    hGeneralRatio->GetXaxis()->SetRangeUser(0.05, 1.05);
    hGeneralRatio->GetXaxis()->SetNdivisions(505);
    hGeneralRatio->GetXaxis()->SetLabelSize(0.04);
    hGeneralRatio->GetXaxis()->SetTitleSize(0.05);
    hGeneralRatio->Draw("same");
    hRatioPbPb_pp->Draw("histpesame");
    TLegend *legZtGenRatio_RTitle = LegStd(legZtGenRatio_RTitle, 0.13, 0.75, 0.13, 0.96);
    legZtGenRatio_RTitle->SetHeader("ALICE simulation");
    legZtGenRatio_RTitle->AddEntry("", "pp,#color[0]{.}#sqrt{s} = 5.02 TeV", "");
    legZtGenRatio_RTitle->AddEntry("", "PYTHIA 8,#color[0]{.}#gamma#font[122]{-}jet", "");
    legZtGenRatio_RTitle->AddEntry("", "#it{p}_{T}^{iso, ch} < 1.5 GeV/#it{c}" , "");
    legZtGenRatioPbPb_pp->AddEntry(hRatioPbPb_pp, "#frac{D(#it{z}_{T},#color[0]{.}12#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c})}{D(#it{z}_{T},#color[0]{.}18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c})}", "lp");
                                            
    // hZtMCGen18_40_Ratio18_40[iCen][iFile]->Draw("histpesame");
    //legZtGenRatioto18_40[iCen][iFile] = LegStd(legZtGenRatioto18_40[iCen][iFile], 0.62, 0.18, 0.95, 0.45);
    //legZtGenRatioto18_40[iCen][iFile]->AddEntry(hRatioPbPb_pp, "#frac{12#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}18 GeV/#it{c}}{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}", "lp");
    //legZtGenRatioto18_40[iCen][iFile]->AddEntry(hZtMCGen12_40_Ratio18_40[iCen][iFile], "#frac{12#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}", "lp");
    // legZtGenRatioto18_40[iCen][iFile]->AddEntry(hZtMCGen18_40_Ratio18_40[iCen][iFile], "#gamma#font[122]{-}jet, #frac{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}", "lp");
    lineX0->Draw("same");
    legZtGenRatioPbPb_pp->Draw("same");
    legZtGenRatio_RTitle->Draw("same");
    //legZtGen[0][0]->Draw("same");
    cZtGenRatioPbPb_pp->Print(outDirPlot[0] + "/cZtGenRatiotoPbPb_pp.pdf");

    TCanvas *cZtGenRatioPbPb_ppDiffR = canvasStdIaa("cZtGenRatioPbPb_ppDiffR", 1, 1);
    TLegend *legZtGenRatioPbPb_ppDiffR = LegStd(legZtGenRatioPbPb_ppDiffR, 0.45, 0.18, 0.86, 0.45);
    TH1F* hRatioPbPb_ppDiffR = (TH1F*)hZtMCGen12_40[0][1]->Clone("hRatioPbPb_pp");
    hRatioPbPb_ppDiffR->Divide(hZtMCGen18_40[0][0]);
    hRatioPbPb_ppDiffR->GetYaxis()->SetRangeUser(0.7, 1.15);
    // hRatioPbPb_pp
    hRatioPbPb_pp->SetTitle("");
    hGeneralRatio->GetYaxis()->SetRangeUser(0.85, 1.1);
    hGeneralRatio->GetYaxis()->SetTitleSize(0.046);
    hGeneralRatio->GetYaxis()->SetTitleOffset(1.1);
    hGeneralRatio->GetYaxis()->SetLabelSize(0.04);
    hGeneralRatio->GetXaxis()->SetRangeUser(0.05, 1.05);
    hGeneralRatio->GetXaxis()->SetNdivisions(505);
    hGeneralRatio->GetXaxis()->SetLabelSize(0.04);
    hGeneralRatio->GetXaxis()->SetTitleSize(0.05);
    hGeneralRatio->Draw("same");
    hRatioPbPb_ppDiffR->Draw("histpesame");
    TLegend *legZtGenRatio_RTitleDiffR = LegStd(legZtGenRatio_RTitleDiffR, 0.13, 0.75, 0.13, 0.96);
    legZtGenRatio_RTitleDiffR->SetHeader("ALICE simulation");
    legZtGenRatio_RTitleDiffR->AddEntry("", "pp,#color[0]{.}#sqrt{s} = 5.02 TeV", "");
    legZtGenRatio_RTitleDiffR->AddEntry("", "PYTHIA 8,#color[0]{.}#gamma#font[122]{-}jet", "");
    legZtGenRatio_RTitleDiffR->AddEntry("", "#it{p}_{T}^{iso, ch} < 1.5 GeV/#it{c}" , "");
    legZtGenRatioPbPb_ppDiffR->AddEntry(hRatioPbPb_pp, "#frac{D(#it{z}_{T},#color[0]{.}12#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c},#color[0]{.}#it{R} = 0.4)}{D(#it{z}_{T},#color[0]{.}18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c},#color[0]{.}#it{R} = 0.2)}", "lp");
                                            
    // hZtMCGen18_40_Ratio18_40[iCen][iFile]->Draw("histpesame");
    //legZtGenRatioto18_40[iCen][iFile] = LegStd(legZtGenRatioto18_40[iCen][iFile], 0.62, 0.18, 0.95, 0.45);
    //legZtGenRatioto18_40[iCen][iFile]->AddEntry(hRatioPbPb_pp, "#frac{12#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}18 GeV/#it{c}}{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}", "lp");
    //legZtGenRatioto18_40[iCen][iFile]->AddEntry(hZtMCGen12_40_Ratio18_40[iCen][iFile], "#frac{12#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}", "lp");
    // legZtGenRatioto18_40[iCen][iFile]->AddEntry(hZtMCGen18_40_Ratio18_40[iCen][iFile], "#gamma#font[122]{-}jet, #frac{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}", "lp");
    lineX0->Draw("same");
    legZtGenRatioPbPb_ppDiffR->Draw("same");
    legZtGenRatio_RTitleDiffR->Draw("same");
    //legZtGen[0][0]->Draw("same");
    cZtGenRatioPbPb_ppDiffR->Print(outDirPlot[0] + "/cZtGenRatiotoPbPb_ppDiffR.pdf");



    TCanvas *cZtGenRatio_R = canvasStdIaa("cZtGenRatio_R", 1, 1);
    TLegend *legZtGenRatio_R = LegStd(legZtGenRatio_R, 0.50, 0.20, 0.90, 0.55);
    legZtGenRatio_R->SetTextSize(0.036);
    TH1F* hRatio_R = (TH1F*)hZtMCGen18_40[0][1]->Clone("hRatio_R");
    hRatio_R->Divide(hZtMCGen18_40[0][0]);
    hRatio_R->GetYaxis()->SetRangeUser(0.7, 1.15);
    // hRatio_R
    hRatio_R->SetTitle("");
    hGeneralRatio->GetYaxis()->SetRangeUser(0.85, 1.1);
    hGeneralRatio->GetYaxis()->SetTitleSize(0.046);
    hGeneralRatio->GetYaxis()->SetTitleOffset(1.1);
    hGeneralRatio->GetYaxis()->SetLabelSize(0.04);
    hGeneralRatio->GetXaxis()->SetRangeUser(0.05, 1.05);
    hGeneralRatio->GetXaxis()->SetNdivisions(505);
    hGeneralRatio->GetXaxis()->SetLabelSize(0.04);
    hGeneralRatio->GetXaxis()->SetTitleSize(0.05);
    hGeneralRatio->Draw("same");
    hRatio_R->Draw("histpesame");
    hRatioPbPb_ppDiffR->Draw("histpesame");
    hRatioPbPb_pp->SetMarkerStyle(25);
    hRatioPbPb_pp->Draw("histpesame");
    //legZtGenRatio_R->AddEntry(hRatio_R, "#frac{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}#color[0]{.}}{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}#color[0]{.}}, #it{R}=0.4/0.2", "lp");
    legZtGenRatio_R->AddEntry(hRatio_R,"#frac{D(#it{z}_{T},#color[0]{.}18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c},#color[0]{.}#bf{#it{R} = 0.4})}{D(#it{z}_{T},#color[0]{.}18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c},#color[0]{.}#bf{#it{R} = 0.2})}", "lp");
    legZtGenRatio_R->AddEntry(hRatioPbPb_pp, "#frac{D(#it{z}_{T},#color[0]{.}#bf{12#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40} GeV/#it{c},#color[0]{.}#it{R} = 0.2)}{D(#it{z}_{T},#color[0]{.}#bf{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40} GeV/#it{c},#color[0]{.}#it{R} = 0.2)}", "lp");
    legZtGenRatio_R->AddEntry(hRatioPbPb_ppDiffR, "#frac{D(#it{z}_{T},#color[0]{.}#bf{12#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40} GeV/#it{c},#color[0]{.}#bf{#it{R} = 0.4})}{D(#it{z}_{T},#color[0]{.}#bf{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40} GeV/#it{c},#color[0]{.}#bf{#it{R} = 0.2})}", "lp");
                                            
    // hZtMCGen18_40_Ratio18_40[iCen][iFile]->Draw("histpesame");
    //legZtGenRatioto18_40[iCen][iFile] = LegStd(legZtGenRatioto18_40[iCen][iFile], 0.62, 0.18, 0.95, 0.45);
    //legZtGenRatioto18_40[iCen][iFile]->AddEntry(hRatio_R, "#frac{12#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}18 GeV/#it{c}}{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}", "lp");
    //legZtGenRatioto18_40[iCen][iFile]->AddEntry(hZtMCGen12_40_Ratio18_40[iCen][iFile], "#frac{12#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}", "lp");
    // legZtGenRatioto18_40[iCen][iFile]->AddEntry(hZtMCGen18_40_Ratio18_40[iCen][iFile], "#gamma#font[122]{-}jet, #frac{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}{18#color[0]{.}<#color[0]{.}#it{p}_{T}^{#gamma}#color[0]{.}<#color[0]{.}40 GeV/#it{c}}", "lp");
    lineX0->Draw("same");
    legZtGenRatio_R->Draw("same");
    legZtGenRatio_RTitle->Draw("same");
    cZtGenRatio_R->Print(outDirPlot[0] + "/Ratio_DiffpTtrig_DiffR.pdf");

    /*
    TFile *fileNLO = new TFile("RootFiles/fileNLO.root ");
    TH1F *grDztNLOmedianpp = (TH1F *)fileNLO->Get(Form("grDztNLOmedian_pp"));
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
    TH1F *hGeneral = new TH1F("hGeneral", "hGeneral", 100, 0, 1);
    // hPYTHIApaper->GetXaxis()->SetRangeUser(0.,1.05);
    // hZtMCGen12_40[2]->GetXaxis()->SetRangeUser(0.,1.05);
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
    cZtGenDist_PYTHIA->Print(outDirPlot[0] + "/cZtGenDist_PYTHIA.pdf");
    */
}