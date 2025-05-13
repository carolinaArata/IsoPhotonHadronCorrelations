#include <TCanvas.h>
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TStyle.h"
#include "stdio.h"
#include <iostream>
#include "TGraph.h"
#include "TSystem.h"
#include "TLegend.h"
#include <vector>
#include <algorithm>
#include <iterator>

using std::cout;
using std::endl;

Int_t nPtTrig = 11;
Float_t PtTrigger[] = {12, 14, 16, 18, 20, 25, 30, 40, 50, 60};
int nZtBins = 6;
double assocZt[] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 1.00};
Double_t range_UE[2] = {3 * (TMath::Pi()) / 10, TMath::Pi() / 2};

const int nCentMixUsed = 3;
TString CentMixUsed[] = {"Output_checkCode", "Output_checkCodeSystNMix18", "Output_checkCodeSystNMix45"};
// TString CentMixUsed[] = {"NMix9_ZtMergedMore", "NMix45_ZtMergedMore"}; files with ShSh Bkg 0.40-2.00
// TString CentMixUsed[] = {"~/work/histogram/FromScratch/checkCode", "~/work/histogram/FromScratch/checkCodeSystNCentrMix"};
TString shshLeg[] = {"0.10-0.30", "0.40-1.00"};
TFile *fPlot[nCentMixUsed];

Int_t kColor[] = {kRed - 7, kAzure - 2, kAzure + 3};
Int_t kStyle[] = {24, 25, 46};
Int_t kStyleSame[] = {20, 25, 47};

void PlotStyle(TH1F *hPlot, int kMarker, double kMarkerSize, int kColor, TString titleX, TString titleY);
double fZYAM_Mix(TH1F *hSame, TH1F *hMix);
TLatex *LatexStd(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax);
void PlotNCentMix(Float_t ptMin = 18, Float_t ptMax = 40, bool Mirror = true, TString shshBkg = "0.40-1.00", bool b0_30 = true)
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
		dirPlot = "Systematics_checkCode0_30/SystNCentXMix";
	}
	else if (!b0_30)
	{
		nCen = 4;
		cenBins.push_back(0);
		cenBins.push_back(10);
		cenBins.push_back(30);
		cenBins.push_back(50);
		cenBins.push_back(90);
		dirPlot = "Systematics_checkCode/SystNCentXMix";
	}

	TString processline = Form(".! mkdir -pv %s",dirPlot.Data()) ;
        gROOT->ProcessLine(processline.Data());

	TString sPtAll = Form("_Pt%2.0f_%2.0f", ptMin, ptMax);
	TFile *fNMixCentSyst = new TFile(Form("%s/fNMixCentSyst%s.root", dirPlot.Data(), sPtAll.Data()), "RECREATE");
	
	TH1F *hDPhiMixedPreZYAM[nCentMixUsed][nCen][nPtTrig][nZtBins];
	TH1F *hDPhiMixed[nCentMixUsed][nCen][nPtTrig][nZtBins];
	TH1F *hDPhiSameNoUE[nCentMixUsed][nCen][nPtTrig][nZtBins];
	TH1F *hDPhiSame[nCentMixUsed][nCen][nPtTrig][nZtBins];
	TH1F *hZt[nCentMixUsed][nCen];
	TH1F *hNMixCentUncer[nCen];
	TH1F *hDiv[nCentMixUsed][nCen];


	TString mirror;
	if (Mirror)
		mirror = "Mirror";
	else
		mirror = "NoMirror";

	TString direct = "Mixed";

	TF1 *fa0[nCen];
	TF1 *fConst[nCen];
	TLatex *ALICEtex[nCen];

	for (int iCen = 0; iCen < nCen; iCen++)
	{
		TString sCent = Form("Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);

		hNMixCentUncer[iCen] = new TH1F(Form("hZtNMixCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("hZtNMixCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), nZtBins, assocZt);
		for (int iMix = 0; iMix < nCentMixUsed; iMix++)
		{
			fPlot[iMix] = new TFile(Form("%s/fPlot%s_%s_Pt%2.0f_%2.0f.root", CentMixUsed[iMix].Data(), shshBkg.Data(), sCent.Data(), ptMin, ptMax));
			cout << fPlot[iMix] << endl;
			hZt[iMix][iCen] = (TH1F *)fPlot[iMix]->Get(Form("hZtEffCorrIso1Photon_%s%s", sCent.Data(), sPtAll.Data()));
			cout << hZt[iMix][iCen] << endl;
			hZt[iMix][iCen]->SetDirectory(0);

			PlotStyle(hZt[iMix][iCen], kStyleSame[iMix], 2, kColor[iMix], "#it{z}_{T}", "1/N^{trig}dN^{charg}/d#font[12]{z}_{T}");
			cout << "NCentMix: " << CentMixUsed[iMix] << endl;
		}
	}

	TCanvas *cZtNCentMix = new TCanvas("cZtNCentMix", "cZtNCentMix", 3 * 800, 1 * 600);
	cZtNCentMix->Divide(3, 1);
	TLegend *legZTData[nCen];

	TCanvas *cZtRatio = new TCanvas("cZtRatio", "cZtRatio", 3 * 800, 1 * 600);
	cZtRatio->Divide(3, 1);
	TLegend *legZTRatioNMix[nCen];
	TF1 *hfitRatio[nCentMixUsed][nCen];

	TCanvas *cNMixUncert = new TCanvas("cNMixUncert", "cNMixUncert", 3 * 800, 1 * 600);
	cNMixUncert->Divide(3, 1);
	double nMixBinSyst[2][nCen][nZtBins];
	TH1F *hFitUncert[nCen];
	for (int iCen = 0; iCen < nCen; iCen++)
	{
		TString sCent = Form("Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
		for (int iMix = 0; iMix < nCentMixUsed; iMix++)
		{
			cZtNCentMix->cd(iCen + 1);
			gPad->SetLogy();
			gStyle->SetPadRightMargin(0.02);
			gStyle->SetPadLeftMargin(0.35);
			gStyle->SetPadBottomMargin(0.15);
			hZt[iMix][iCen]->GetXaxis()->SetTitleSize(0.042);
			hZt[iMix][iCen]->GetYaxis()->SetTitleSize(0.042);
			hZt[iMix][iCen]->GetXaxis()->SetLabelSize(0.04);
			hZt[iMix][iCen]->GetYaxis()->SetLabelSize(0.04);
			hZt[iMix][iCen]->SetTitle("");
			hZt[iMix][iCen]->Draw("same");
			ALICEtex[iCen] = LatexStd(ALICEtex[iCen], 0.420, 0.84, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax);
		}
		legZTData[iCen] = new TLegend(0.65, 0.55, 0.85, 0.68);
		legZTData[iCen]->SetFillColor(kWhite);
		legZTData[iCen]->SetLineWidth(0);
		legZTData[iCen]->AddEntry(hZt[0][iCen], "Nominal");
		legZTData[iCen]->AddEntry(hZt[1][iCen], "N Mix = 18");
		legZTData[iCen]->AddEntry(hZt[2][iCen], "N Mix = 36");
		legZTData[iCen]->Draw("same");
		for (int iSysMix = 0; iSysMix < nCentMixUsed-1; iSysMix++)
		{
			hfitRatio[iSysMix][iCen] = new TF1(Form("hfitRatio%d_%d", cenBins[iCen], cenBins[iCen + 1]), "pol0", 0.2, 0.8);
			hDiv[iSysMix][iCen] = (TH1F *)hZt[iSysMix+1][iCen]->Clone(Form("Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
			hDiv[iSysMix][iCen]->SetDirectory(0);
			for (int ibin = 0; ibin < nZtBins; ibin++)
			{
				hDiv[iSysMix][iCen]->SetBinError(ibin + 1, 0);
			}
			
			PlotStyle(hDiv[iSysMix][iCen],kStyleSame[iSysMix], 2, kColor[iSysMix], "#it{z}_{T}", "NMixCentBin18-36 / Nomin");
			hDiv[iSysMix][iCen]->Divide(hZt[0][iCen]);
			cZtRatio->cd(iCen + 1);
			gStyle->SetPadRightMargin(0.02);
			gStyle->SetPadLeftMargin(0.35);
			gStyle->SetPadBottomMargin(0.15);
			hDiv[iSysMix][iCen]->GetXaxis()->SetTitleSize(0.042);
			hDiv[iSysMix][iCen]->GetYaxis()->SetTitleSize(0.042);
			hDiv[iSysMix][iCen]->GetXaxis()->SetLabelSize(0.04);
			hDiv[iSysMix][iCen]->GetYaxis()->SetLabelSize(0.04);
			hDiv[iSysMix][iCen]->Draw("plhist same");
			hDiv[iSysMix][iCen]->GetYaxis()->SetRangeUser(-0.5, 3);
			hDiv[iSysMix][iCen]->Fit(Form("hfitRatio%d_%d", cenBins[iCen], cenBins[iCen + 1]), "R0Q");
			ALICEtex[iCen] = LatexStd(ALICEtex[iCen], 0.160, 0.84, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax);
			hfitRatio[iSysMix][iCen]->Draw("same");
		}
		legZTRatioNMix[iCen] = new TLegend(0.65, 0.55, 0.85, 0.68);
		legZTRatioNMix[iCen]->SetFillColor(kWhite);
		legZTRatioNMix[iCen]->SetLineWidth(0);
		legZTRatioNMix[iCen]->AddEntry(hDiv[0][iCen], "N Mix = 18 / Nominal");
		legZTRatioNMix[iCen]->AddEntry(hDiv[1][iCen], "N Mix = 36 / Nominal");
		legZTRatioNMix[iCen]->Draw("same");
		//for (int ibin = 0; ibin < nZtBins; ibin++)
		//{
		//	double diff = 0;
		//	diff = abs(hZt[2][iCen]->GetBinContent(ibin + 1) - hZt[1][iCen]->GetBinContent(ibin + 1)) / abs(hZt[0][iCen]->GetBinContent(ibin + 1));
		//	hNMixCentUncer[iCen]->SetBinContent(ibin + 1, diff);
		//	hNMixCentUncer[iCen]->SetBinError(ibin + 1, hZt[0][iCen]->GetBinError(ibin + 1) / abs((hZt[0][iCen]->GetBinContent(ibin + 1))));
		//}

		for (int ibin = 0; ibin < nZtBins; ibin++)
		{
			nMixBinSyst[0][iCen][ibin] = abs(hZt[1][iCen]->GetBinContent(ibin + 1) - hZt[0][iCen]->GetBinContent(ibin + 1)) / abs(hZt[0][iCen]->GetBinContent(ibin + 1));
			nMixBinSyst[1][iCen][ibin] = abs(hZt[2][iCen]->GetBinContent(ibin + 1) - hZt[0][iCen]->GetBinContent(ibin + 1)) / abs(hZt[0][iCen]->GetBinContent(ibin + 1));
                        hNMixCentUncer[iCen]->SetBinContent(ibin + 1, (nMixBinSyst[0][iCen][ibin]+nMixBinSyst[1][iCen][ibin])/ 2);
			hNMixCentUncer[iCen]->SetBinError(ibin + 1, hZt[0][iCen]->GetBinError(ibin + 1) / abs((hZt[0][iCen]->GetBinContent(ibin + 1))));
		}


		PlotStyle(hNMixCentUncer[iCen], 20, 2, kBlue + 2, "#it{z}_{T}", "Uncertainty %");
		cNMixUncert->cd(iCen + 1);
		hNMixCentUncer[iCen]->Scale(100);
		hNMixCentUncer[iCen]->GetYaxis()->SetRangeUser(0, 50);
		gStyle->SetPadRightMargin(0.02);
		gStyle->SetPadLeftMargin(0.35);
		gStyle->SetPadBottomMargin(0.15);
		hNMixCentUncer[iCen]->GetXaxis()->SetTitleSize(0.042);
		hNMixCentUncer[iCen]->GetYaxis()->SetTitleSize(0.042);
		hNMixCentUncer[iCen]->GetXaxis()->SetLabelSize(0.04);
		hNMixCentUncer[iCen]->GetYaxis()->SetLabelSize(0.04);
		hNMixCentUncer[iCen]->SetMarkerColor(kBlue + 2);
		hNMixCentUncer[iCen]->SetTitle(" ");
                hNMixCentUncer[iCen]->Draw("histplsame");
		ALICEtex[iCen] = LatexStd(ALICEtex[iCen], 0.160, 0.84, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax);
		if (iCen == 2 || iCen == 3)
    {
      //fa0[iCen] = new TF1(Form("NMixfit0%d_%d", cenBins[iCen], cenBins[iCen + 1]), "expo", 0.1, 0.8);
      fa0[iCen] = new TF1(Form("NMixfit0%d_%d", cenBins[iCen], cenBins[iCen + 1]), "expo", 0.15, 0.59);
      fConst[iCen] = new TF1(Form("ConstNMixfit0%d_%d", cenBins[iCen], cenBins[iCen + 1]), "pol0", 0.2, 0.8);
    }
    if (iCen == 1)
    {
      //fa0[iCen] = new TF1(Form("NMixfit0%d_%d", cenBins[iCen], cenBins[iCen + 1]), "expo", 0.1, 0.8);
      fa0[iCen] = new TF1(Form("NMixfit0%d_%d", cenBins[iCen], cenBins[iCen + 1]), "expo", 0.25, 0.8);
      fConst[iCen] = new TF1(Form("ConstNMixfit0%d_%d", cenBins[iCen], cenBins[iCen + 1]), "pol0", 0.3, 0.6);
    }
    if (iCen == 0)
    {
      //fa0[iCen] = new TF1(Form("NMixfit0%d_%d", cenBins[iCen], cenBins[iCen + 1]), "expo", 0.2, 0.8);
      fa0[iCen] = new TF1(Form("NMixfit0%d_%d", cenBins[iCen], cenBins[iCen + 1]), "[0]*exp([1]*x)+[2]", 0.3, 0.8);
      fConst[iCen] = new TF1(Form("ConstNMixfit0%d_%d", cenBins[iCen], cenBins[iCen + 1]), "pol0", 0.3, 0.8);
      hNMixCentUncer[iCen]->Fit(Form("ConstNMixfit0%d_%d", cenBins[iCen], cenBins[iCen + 1]), "R");
      // fConst->SetLineColor(kBlue);
    }
   // gStyle->SetOptFit(1111);

    hNMixCentUncer[iCen]->Fit(Form("NMixfit0%d_%d", cenBins[iCen], cenBins[iCen + 1]), "R");

    fa0[iCen]->SetLineColor(kRed + 1);
    fa0[iCen]->SetLineStyle(1);
    fa0[iCen]->SetLineWidth(3);
    fa0[iCen]->Draw("same");
    hNMixCentUncer[iCen]->Fit(Form("ConstNMixfit0%d_%d", cenBins[iCen], cenBins[iCen + 1]), "R");
    fConst[iCen]->SetLineColor(kGreen + 3);
    fConst[iCen]->SetLineWidth(4);
    fConst[iCen]->Draw("same");

    hFitUncert[iCen] = new TH1F(Form("hFromFitUncertNMix_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("hFromFitUncertNMix_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), nZtBins, assocZt);
    for (int ibin = 0; ibin < nZtBins; ibin++)
    {
      if ( iCen == 1 )
      {
        hFitUncert[iCen]->SetBinContent(ibin + 1, fConst[iCen]->Eval(hFitUncert[iCen]->GetBinCenter(ibin + 1)));
        // Average Central and peripheral
        
        //        hFitUncert[iCen]->SetBinContent(ibin + 1,
        //                                        (fa0[0]->Eval(hFitUncert[iCen]->GetBinCenter(ibin + 1)) +
        //                                         fa0[2]->Eval(hFitUncert[iCen]->GetBinCenter(ibin + 1)))/2);
      }
      else
      {
        //hFitUncert[iCen]->SetBinContent(ibin + 1, fa0[iCen]->Eval(hFitUncert[iCen]->GetBinCenter(ibin + 1)));
        hFitUncert[iCen]->SetBinContent(ibin + 1, fConst[iCen]->Eval(hFitUncert[iCen]->GetBinCenter(ibin + 1)));
      }
//			if (iCen == 1 || iCen == 2)
//      {
//        //hFitUncert[iCen]->SetBinContent(ibin + 1, fa0[iCen]->Eval(hFitUncert[iCen]->GetBinCenter(ibin + 1)));
//        hFitUncert[iCen]->SetBinContent(ibin + 1, fConst[iCen]->Eval(hFitUncert[iCen]->GetBinCenter(ibin + 1)));
//      }
//			else if (iCen == 0)
//			{
//				if (ibin == 0 || ibin == 1 || ibin == 2)
//				{
//					hFitUncert[iCen]->SetBinContent(ibin + 1, fConst[iCen]->Eval(hFitUncert[iCen]->GetBinCenter(ibin + 1)));
//				}
//				else
//				{
//					//hFitUncert[iCen]->SetBinContent(ibin + 1, fa0[iCen]->Eval(hFitUncert[iCen]->GetBinCenter(ibin + 1)));
//					hFitUncert[iCen]->SetBinContent(ibin + 1, fConst[iCen]->Eval(hFitUncert[iCen]->GetBinCenter(ibin + 1)));
//				}
//			}
		}
		fNMixCentSyst->cd();
		fConst[iCen]->Write();
		fa0[iCen]->Write();
		hNMixCentUncer[iCen]->Write();
		hFitUncert[iCen]->Write();
	}
	cZtNCentMix->Print(Form("%s/ZtDistrib%s.pdf", dirPlot.Data(), sPtAll.Data()));
	cZtRatio->Print(Form("%s/ZtRatio%s.pdf", dirPlot.Data(), sPtAll.Data()));
	cNMixUncert->Print(Form("%s/NmixCentrMatch%s.pdf", dirPlot.Data(), sPtAll.Data()));

	// Systematics Icp
	TH1F *Icp[nCentMixUsed][nCen];
	for (int iCen = 0; iCen < nCen; iCen++)
	{
		Icp[0][iCen] = (TH1F *)hZt[0][iCen]->Clone(Form("IcpIso1Photon_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
		Icp[0][iCen]->Divide(hZt[0][nCen - 1]);
		PlotStyle(Icp[0][iCen], kStyleSame[0], 2, kColor[0], "#it{z}_{T}", "#it{I}_{CP}");
		Icp[1][iCen] = (TH1F *)hZt[1][iCen]->Clone(Form("IcpIso1PhotonNcenMix18_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
		Icp[1][iCen]->Divide(hZt[1][nCen - 1]);
		PlotStyle(Icp[1][iCen], kStyleSame[1], 2, kColor[1], "#it{z}_{T}", "#it{I}_{CP}");
		Icp[2][iCen] = (TH1F *)hZt[2][iCen]->Clone(Form("IcpIso1PhotonNcenMix45_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]));
		Icp[2][iCen]->Divide(hZt[2][nCen - 1]);
		PlotStyle(Icp[2][iCen], kStyleSame[2], 2, kColor[2], "#it{z}_{T}", "#it{I}_{CP}");
	}

	TCanvas *cIcpSyst = new TCanvas("cIcpSyst", "cIcpSyst", (nCen - 1) * 800, 1 * 600);
	cIcpSyst->Divide((nCen - 1), 1);
	TLegend *legZTDataIcp[(nCen - 1)];
	TLatex *lat1[(nCen - 1)];
        TF1 * fConstIcp[(nCen - 1)] ;
  
	for (int iCen = 0; iCen < (nCen - 1); iCen++)
	{
		for (int iMix = 0; iMix < nCentMixUsed; iMix++)
		{
			cIcpSyst->cd(iCen + 1);
			Icp[iMix][iCen]->GetYaxis()->SetRangeUser(-0.5, 3);
			Icp[iMix][iCen]->Draw("same");
		}
		lat1[iCen] = new TLatex();
		lat1[iCen]->SetTextFont(42);
		lat1[iCen]->SetTextSize(0.04);
		lat1[iCen]->SetNDC();
		lat1[iCen]->DrawLatex(0.15, 0.85, Form("ALICE #it{Work in progress}"));
		lat1[iCen]->DrawLatex(0.15, 0.85 - 0.06, Form("#bf{%d-%d %% / 50-90 %%} Pb#font[122]{-}Pb, #sqrt{s_{NN}} = 5.02 TeV", cenBins[iCen], cenBins[iCen + 1]));
		lat1[iCen]->DrawLatex(0.15, 0.85 - 2 * 0.06, Form("|#it{#eta}^{ #it{#gamma}}| < 0.67,%2.0f < #it{p}_{#it{T}}^{#gamma} < %2.0f GeV/#it{c}", ptMin, ptMax));
		legZTDataIcp[iCen] = new TLegend(0.70, 0.700, 0.86, 0.82);
		legZTDataIcp[iCen]->SetFillColor(kWhite);
		legZTDataIcp[iCen]->SetLineWidth(0);
		legZTDataIcp[iCen]->AddEntry(hZt[0][iCen], "N Mix = 9");
		legZTDataIcp[iCen]->AddEntry(hZt[1][iCen], "N Mix = 18");
		legZTDataIcp[iCen]->AddEntry(hZt[2][iCen], "N Mix = 36");
		legZTDataIcp[iCen]->Draw("same");
	}
	cIcpSyst->Print(Form("%s/IcpDistrib%s.pdf", dirPlot.Data(), sPtAll.Data()));

	TH1F *hUncertIcpNcentMix[nCen];
	TH1F *hfitIcpUncert[nCen];
	for (int iCen = 0; iCen < nCen; iCen++)
	{
	    hUncertIcpNcentMix[iCen] = new TH1F(Form("hUncertIcpNcentMix_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("hUncertIcpNcentMix_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), nZtBins, 
assocZt);
  
	    for (int ibin = 0; ibin < nZtBins; ibin++)
	    {
	      double binCont = 100 * (abs(Icp[2][iCen]->GetBinContent(ibin + 1) - Icp[1][iCen]->GetBinContent(ibin + 1)) / abs(Icp[0][iCen]->GetBinContent(ibin + 1)));
	      hUncertIcpNcentMix[iCen]->SetBinContent(ibin + 1, binCont);
	    }
	    hfitIcpUncert[iCen] = new TH1F(Form("hFromFitUncert_Icp_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]),
                                   Form("hFromFitUncert_Icp_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), nZtBins, assocZt);
	    fNMixCentSyst->cd();
    
            fConstIcp[iCen] = new TF1(Form("ConstIcpNMixfit0%d_%d", cenBins[iCen], cenBins[iCen + 1]), "pol0", 0.2, 0.55);
            hUncertIcpNcentMix[iCen]->Fit(Form("ConstIcpNMixfit0%d_%d", cenBins[iCen], cenBins[iCen + 1]),"R");
    
            for (int ibin = 0; ibin < nZtBins; ibin++)
            {
              hfitIcpUncert[iCen]->SetBinContent(ibin + 1, fConstIcp[iCen]->Eval(hUncertIcpNcentMix[iCen]->GetBinCenter(ibin + 1)));
            }
            hUncertIcpNcentMix[iCen]->Write();
            hfitIcpUncert[iCen]->Write();
        }

	TCanvas *cIcpSystUncert = new TCanvas("cIcpSystUncert", "cIcpSystUncert", (nCen - 1) * 800, 1 * 600);
	cIcpSystUncert->Divide((nCen - 1), 1);

	for (int iCen = 0; iCen < nCen - 1; iCen++)
	{
		PlotStyle(hUncertIcpNcentMix[iCen], 20, 1, kBlue + 3, "#it{z}_{T}", "Uncertainty %");
		hUncertIcpNcentMix[iCen]->GetYaxis()->SetRangeUser(0, 50);
		cIcpSystUncert->cd(iCen + 1);
		hUncertIcpNcentMix[iCen]->Fit(Form("hfitIcpUncert_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), "R");
		hUncertIcpNcentMix[iCen]->SetTitle("");
		hUncertIcpNcentMix[iCen]->Draw("hist pl");
    
                fConstIcp[iCen]->SetLineColor(kGreen+3);
                fConstIcp[iCen]->SetLineWidth(4);
                fConstIcp[iCen]->Draw("same");
		TLatex *lat = new TLatex();
		lat->SetTextFont(42);
		lat->SetTextSize(0.04);
		lat->SetNDC();
		lat->DrawLatex(0.15, 0.85, Form("ALICE #it{Work in progress}"));
		lat->DrawLatex(0.15, 0.85 - 0.06, Form("#bf{%d-%d %% / 50-90 %%} Pb-Pb, #sqrt{s_{NN}} = 5.02 TeV", cenBins[iCen], cenBins[iCen + 1]));
		lat->DrawLatex(0.15, 0.85 - 2 * 0.06, Form("|#it{#eta}^{ #it{#gamma}}| < 0.67, %2.0f < #it{p}_{#it{T}}^{#gamma} < %2.0f GeV/#it{c}", ptMin, ptMax));
	}
	cIcpSystUncert->Print(Form("%s/IcpUncertDistrib%s.pdf", dirPlot.Data(), sPtAll.Data()));
}

void PlotStyle(TH1F *hPlot, int kMarker, double kMarkerSize, int kColor, TString titleX, TString titleY)
{
	gStyle->SetTitleX(0.56);
	gStyle->SetOptTitle(1);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(00000);
	// gStyle->SetLabelSize(.5, "XY");
	gStyle->SetLineScalePS(1);
	gStyle->SetTitleAlign(23);
	gStyle->SetPadRightMargin(0.03);
	gStyle->SetPadLeftMargin(0.12);
	hPlot->SetMarkerStyle(kMarker);
	// hPlot->SetMarkerSize(1.1);
	hPlot->SetMarkerSize(1.5);
	hPlot->SetMarkerColor(kColor);
	hPlot->SetLineColor(kColor);
	hPlot->SetLineWidth(1);
	hPlot->GetYaxis()->SetTitle(Form("%s", titleY.Data()));
	hPlot->GetXaxis()->SetTitle(Form("%s", titleX.Data()));
	// hPlot->GetXaxis()->SetLabelSize(.05);
	// hPlot->GetYaxis()->SetLabelSize(.05);
	//  leg->SetFillColor(kWhite);
	//  leg->SetLineColor(0);
}
double fZYAM_Mix(TH1F *hSame, TH1F *hMix)
{
	double minBin = hSame->GetXaxis()->FindBin(range_UE[0] + 0.00001);
	double maxBin = hSame->GetXaxis()->FindBin(range_UE[1]);
	double IntSame = hSame->Integral(minBin, maxBin);
	// printf("norm Same %f \n", IntSame);
	minBin = hMix->GetXaxis()->FindBin(range_UE[0] + 0.00001);
	maxBin = hMix->GetXaxis()->FindBin(range_UE[1]);
	double IntMix = hMix->Integral(minBin, maxBin);
	// printf("norm Mix %f \n", IntMix);

	double CorrFact = IntSame / IntMix;
	// printf("CorrFactor %f \n", CorrFact);
	hMix->Scale(CorrFact);
	return CorrFact;
}

TLatex *LatexStd(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax)
{
	lat = new TLatex();
	lat->SetTextFont(42);
	lat->SetTextSize(0.04);
	lat->SetNDC();
	lat->DrawLatex(xpos, ypos, Form("ALICE #it{Work in progress}"));
	lat->DrawLatex(xpos, ypos - 0.06, Form("%d#font[122]{-}%d %% Pb#font[122]{-}Pb, #sqrt{s_{NN}} = 5.02 TeV", cenMin, cenMax));
	lat->DrawLatex(xpos, ypos - 2 * 0.06, Form("|#it{#eta}^{ #it{#gamma}}| < 0.67 , %2.0f < #it{p}_{T}^{#gamma} < %2.0f GeV/#it{c}", ptMin, ptMax));

	return lat;
}
