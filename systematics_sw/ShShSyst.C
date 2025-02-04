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

const Int_t nPtTrig = 9;
Float_t PtTrigger[] = {12, 14, 16, 18, 20, 25, 30, 40, 50, 60};

int nZtBin = 6;
double assocZt[] = {0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 1.00};

// Definition of various shsh ranges
int const nShSh = 5;
TString shshString[] = {"", "04_10PtDep", "04_15PtDep", "05_20PtDep", "06_20PtDep"};
TString shshLeg[] = {"0.40-1.00", "0.40-1.50", "0.40-2.00", "0.35-1.00", "0.50-2.00", "0.60-2.00"};
Double_t shMin[] = {0.40, 0.40, 0.40, 0.50, 0.6};
Double_t shMax[] = {2.00, 1.00, 1.50, 2.00, 2.00};

TFile *fPlot[6][4];

Int_t kColor[] = {2, kCyan + 2, kOrange + 7, kAzure + 3, kTeal - 6};
Int_t kStyle[] = {20, 20, 25, 47, 103};

void PlotStyle(TH1F *hPlot, int kMarker, double kMarkerSize, int kColor, TString titleX, TString titleY);

TLatex *LatexStd(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax);
TLegend *LegStd(TLegend *leg, double xpos1, double ypos1, double xpos2, double ypos2);
void ShShSyst(Float_t ptMin = 18, Float_t ptMax = 40, TString Mixed = "Mixed", bool Mirror = true, TString dirFiles = "~/work/histogram/FromScratch/checkCodeSystShSh", TString shshBkg = "0.40-1.00", TString dirRef = "~/work/histogram/FromScratch/checkCode", bool b0_30 = false)
{

	Int_t nCen;
	std::vector<Int_t> cenBins;
	if (b0_30)
	{
		nCen = 3;
		// Int_t cenBins[] = {0, 10, 30, 50, 90};
		cenBins.push_back(0);
		cenBins.push_back(30);
		cenBins.push_back(50);
		cenBins.push_back(90);
	}
	else if (!b0_30)
	{
		nCen = 4;
		cenBins.push_back(0);
		cenBins.push_back(10);
		cenBins.push_back(30);
		cenBins.push_back(50);
		cenBins.push_back(90);
	}

	TString mirror;
	if (Mirror)
		mirror = "Mirror";
	else
		mirror = "NoMirror";

	TString dirPlot;
	if (b0_30)
	{
		dirPlot = Form("~/work/histogram/Systematics_checkCode0_30/ShShSyst");
		gSystem->Exec(Form("mkdir %s", dirPlot.Data()));
	}
	else if (!b0_30)
	{
		dirPlot = Form("~/work/histogram/Systematics_checkCode/ShShSyst");
		gSystem->Exec(Form("mkdir %s", dirPlot.Data()));
	}
	/*TFile *fIn_Pi0PtDep = new TFile("~/work/histogram/DataShShMax150/EMCAL_MB_0_90.root");
	TFile *fIn_Pi0MultiShBkg = new TFile("~/work/histogram/ShShBkgPtDep/EMCAL_MB_0_90.root");

	TH1F *hPtTrigPi0[nShSh][nCen];
	TH1F *hPtTrigPi0Ratio[nShSh][nCen];

	TH1F *h1AnalyzNEvents[2];

	h1AnalyzNEvents[0] = (TH1F *)fIn_Pi0PtDep->Get("hNEvents");
	h1AnalyzNEvents[1] = (TH1F *)fIn_Pi0MultiShBkg->Get("hNEvents");

	double NumbAnalyzEv = h1AnalyzNEvents[1]->GetBinContent(1);
		h1AnalyzNEvents[0] = (TH1F *)fIn_Pi0PtDep->Get("hNEvents");
	h1AnalyzNEvents[1] = (TH1F *)fIn_Pi0MultiShBkg->Get("hNEvents");

	double NumbAnalyzEv[2];
	NumbAnalyzEv[0] = h1AnalyzNEvents[0]->GetBinContent(1);
	NumbAnalyzEv[1] = h1AnalyzNEvents[1]->GetBinContent(1);

	// Getter Plots
	for (int iCen = 0; iCen < nCen; iCen++)
	{
		hPtTrigPi0[0][iCen] = (TH1F *)fIn_Pi0PtDep->Get(Form("AnaPhotonHadronCorr_Iso1_ShSh%s_Cen%d_%d_hPtTrigger", shshLeg[0].Data(), cenBins[iCen], cenBins[iCen + 1]));
		hPtTrigPi0[1][iCen] = (TH1F *)fIn_Pi0MultiShBkg->Get(Form("AnaPhotonHadronCorr_Iso1_ShSh%s_Cen%d_%d_hPtTrigger", shshLeg[1].Data(), cenBins[iCen], cenBins[iCen + 1]));
		hPtTrigPi0[2][iCen] = (TH1F *)fIn_Pi0MultiShBkg->Get(Form("AnaPhotonHadronCorr_Iso1_ShSh%s_Cen%d_%d_hPtTrigger", shshLeg[2].Data(), cenBins[iCen], cenBins[iCen + 1]));
		hPtTrigPi0[3][iCen] = (TH1F *)fIn_Pi0MultiShBkg->Get(Form("AnaPhotonHadronCorr_Iso1_ShSh%s_Cen%d_%d_hPtTrigger", shshLeg[3].Data(), cenBins[iCen], cenBins[iCen + 1]));
		hPtTrigPi0[4][iCen] = (TH1F *)fIn_Pi0MultiShBkg->Get(Form("AnaPhotonHadronCorr_Iso1_ShSh%s_Cen%d_%d_hPtTrigger", shshLeg[4].Data(), cenBins[iCen], cenBins[iCen + 1]));

		hPtTrigPi0[0][iCen]->Scale(1 / NumbAnalyzEv[0]);
		hPtTrigPi0[1][iCen]->Scale(1 / NumbAnalyzEv[1]);
		hPtTrigPi0[2][iCen]->Scale(1 / NumbAnalyzEv[1]);
		hPtTrigPi0[3][iCen]->Scale(1 / NumbAnalyzEv[1]);
		hPtTrigPi0[4][iCen]->Scale(1 / NumbAnalyzEv[1]);

		for (int iSh = 0; iSh < nShSh; iSh++)
		{
			hPtTrigPi0[iSh][iCen]->SetDirectory(0);
			PlotStyle(hPtTrigPi0[iSh][iCen], kStyle[iSh], 2, kColor[iSh], "p_{#it{T}}^{trig} (GeV/#it{c})", "Entries");
			hPtTrigPi0Ratio[iSh][iCen] = (TH1F *)hPtTrigPi0[iSh][iCen]->Clone(Form("hPtTrig%s_Cent%d_%d", shshLeg[iSh].Data(), cenBins[iCen], cenBins[iCen + 1]));
		}

		h3PtM02SumPtCone_Cent[iCen] = (TH3F *)fIn_Pi0MultiShBkg->Get(Form("AnaIsolPhoton_hPtM02SumPtCone_Cent%d", iCen));
		// h3PtM02SumPtCone_Cent[iCen]->GetZaxis()->SetRange(h3PtM02SumPtCone_Cent[iCen]->GetZaxis()->FindBin(2.0001), h3PtM02SumPtCone_Cent[iCen]->GetZaxis()->GetNbins());
		h2PtM02Sum[iCen] = (TH2F *)h3PtM02SumPtCone_Cent[iCen]->Project3D("xy");
		h2PtM02Sum[iCen]->SetDirectory(0);
	}

	// Plots

	TLegend *legPtTrigPi0Ratio[nCen];
	for (int iCen = 0; iCen < nCen; iCen++)
	{

		TCanvas *cPtTrigPi0 = new TCanvas(Form("cPtTrigPi0"), Form("cPtTrigPi0"), 800, 600);
		TCanvas *cPtTrigPi0Ratio = new TCanvas(Form("cPtTrigPi0Ratio"), Form("cPtTrigPi0Ratio"), 800, 600);
		for (int iSh = 0; iSh < nShSh; iSh++)
		{
			cPtTrigPi0->cd();
			hPtTrigPi0[iSh][iCen]->Draw("same");

			hPtTrigPi0Ratio[iSh][iCen]->Divide(hPtTrigPi0[0][iCen]);
			cPtTrigPi0Ratio->cd();

			hPtTrigPi0Ratio[iSh][iCen]->SetLineWidth(3);
			hPtTrigPi0Ratio[iSh][iCen]->SetMaximum(4);
			hPtTrigPi0Ratio[iSh][iCen]->SetMinimum(-3);
			hPtTrigPi0Ratio[iSh][iCen]->Draw("same");
		}
		legPtTrigPi0Ratio[iCen] = new TLegend(0.65, 0.60, 0.85, 0.80);
		legPtTrigPi0Ratio[iCen]->SetFillColor(kWhite);
		legPtTrigPi0Ratio[iCen]->SetLineWidth(0);
		legPtTrigPi0Ratio[iCen]->AddEntry(hPtTrigPi0Ratio[0][iCen], "Nominal ShShBkg 0.40-2.00", "lp");
		for (int iSh = 1; iSh < nShSh; iSh++)
		{
			legPtTrigPi0Ratio[iCen]->AddEntry(hPtTrigPi0Ratio[iSh][iCen], Form("ShShBkg %s", shshLeg[iSh].Data()), "lp");
		}
		legPtTrigPi0Ratio[iCen]->Draw("same");
	}


	*/

	int nShShSyst = 4;
	TH1F *hZt[nShShSyst][nCen];
	TH1F *hZtSyst[nShShSyst][nCen];

	TH1F *hZtIcp[nShShSyst][nCen];
	TH1F *hZtSystIcp[nShShSyst][nCen];
	TString sPtAll = Form("_Pt%2.0f_%2.0f", ptMin, ptMax);

	TFile *fShSyst = new TFile(Form("%s/fShSyst%s%s%s.root", dirPlot.Data(), Mixed.Data(), shshBkg.Data(), sPtAll.Data()), "RECREATE");
	for (int iCen = 0; iCen < nCen; iCen++)
	{
		TString sCent = Form("Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
		fPlot[0][iCen] = new TFile(Form("%s/fPlot%s_%s%s.root", dirRef.Data(), shshLeg[0].Data(), sCent.Data(), sPtAll.Data()));
		hZt[0][iCen] = (TH1F *)fPlot[0][iCen]->Get(Form("hZtEffCorrIso1Photon_%s%s", sCent.Data(), sPtAll.Data()));
		cout << hZt[0][iCen] << endl;

		for (int iSh = 1; iSh < nShShSyst; iSh++)
		{
			fPlot[iSh][iCen] = new TFile(Form("%s/fPlot%s_%s%s.root", dirFiles.Data(), shshLeg[iSh].Data(), sCent.Data(), sPtAll.Data()));
			hZt[iSh][iCen] = (TH1F *)fPlot[iSh][iCen]->Get(Form("hZtEffCorrIso1Photon_%s%s", sCent.Data(), sPtAll.Data()));
			cout << hZt[iSh][iCen] << endl;
		}
		for (int iSh = 0; iSh < nShShSyst; iSh++)
		{
			hZt[iSh][iCen]->SetDirectory(0);
			PlotStyle(hZt[iSh][iCen], kStyle[iSh], 2, kColor[iSh], "#it{z}_{T}", "1/N^{trig}dN^{charg}/d#it{z}_{T}");
		}
	}

	for (int iCen = 0; iCen < nCen; iCen++)
	{
		TString sCent = Form("Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]);
		for (int iSh = 0; iSh < nShShSyst; iSh++)
		{
			hZtIcp[iSh][iCen] = (TH1F *)hZt[iSh][iCen]->Clone(Form("hIcpIso1Photon%s_Sh%s", sCent.Data(), shshLeg[iSh].Data()));
			if (b0_30)
			{
				hZtIcp[iSh][iCen]->Divide(hZt[iSh][2]);
			}
			else if (!b0_30)
			{
				hZtIcp[iSh][iCen]->Divide(hZt[iSh][3]);
			}
			hZtIcp[iSh][iCen]->SetDirectory(0);
			PlotStyle(hZtIcp[iSh][iCen], kStyle[iSh], 2, kColor[iSh], "#it{z}_{T}", "#it{I}_{CP}");
		}
	}

	// Plot Zt distributions for various ShSh ranges
	TCanvas *cDiffSh[nCen];
	TLegend *legZTData[nCen];
	int const nBins = hZt[0][0]->GetXaxis()->GetNbins();
	TLatex *ALICEtex[nCen];
	for (int iCen = 0; iCen < nCen; iCen++)
	{
		cDiffSh[iCen] = new TCanvas(Form("cDiffSh%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("cDiffSh%d_%d", cenBins[iCen], cenBins[iCen + 1]), 800, 600);
		cDiffSh[iCen]->cd();
		gPad->SetLogy();
		legZTData[iCen] = new TLegend(0.48, 0.480, 0.89, 0.70);
		legZTData[iCen]->SetFillColor(kWhite);
		legZTData[iCen]->SetLineWidth(0);
		legZTData[iCen]->SetTextSize(0.04);
		for (int iSyst = 0; iSyst < 4; iSyst++)
		{
			hZt[iSyst][iCen]->SetTitle("");
			hZt[iSyst][iCen]->Draw("Same");
			// hZt[1][iCen]->SetTitle("");
			// hZt[1][iCen]->Draw("Same");
			// hZt[2][iCen]->SetTitle("");
			// hZt[2][iCen]->Draw("Same");
			// hZt[3][iCen]->SetTitle("");
			// hZt[3][iCen]->Draw("Same");
		}
		legZTData[iCen]->AddEntry(hZt[0][iCen], "Nom. Bkg. #sigma^{2}_{long, 5#times5}: 0.40-1.00");
		legZTData[iCen]->AddEntry(hZt[1][iCen], Form("Bkg. #sigma^{2}_{long, 5#times5}: %s", shshLeg[1].Data()));
		legZTData[iCen]->AddEntry(hZt[2][iCen], Form("Bkg. #sigma^{2}_{long, 5#times5}: %s", shshLeg[2].Data()));
		legZTData[iCen]->AddEntry(hZt[3][iCen], Form("Bkg. #sigma^{2}_{long, 5#times5}: %s", shshLeg[3].Data()));
		legZTData[iCen]->Draw("same");
		ALICEtex[iCen] = LatexStd(ALICEtex[iCen], 0.50, 0.86, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax);
		cDiffSh[iCen]->Print(Form("%s/ZtDistribXdiffShCen%d_%d%s.pdf", dirPlot.Data(), cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
	}

	////////////////////////////////////////////////////////////////////////////////////
	//////////////////Estimate the systematic for the ShSh variation////////////////////
	////////////////////////////////////////////////////////////////////////////////////
	double shSyst[3][nCen][nBins];
	TH1F *hZtSystem[3][nCen];
	TH1F *hZtSystFinal[nCen];
	TLegend *legZtSyst[nCen];
	for (int iCen = 0; iCen < nCen; iCen++)
	{
		cout << "Cent: " << cenBins[iCen] << "-" << cenBins[iCen + 1] << " Systematics: " << endl;
		legZtSyst[iCen] = new TLegend(0.15, 0.42, 0.45, 0.70);
		legZtSyst[iCen]->SetFillColor(kWhite);
		legZtSyst[iCen]->SetLineWidth(0);
		legZtSyst[iCen]->SetTextSize(0.04);
		hZtSystem[2][iCen] = new TH1F(Form("hZtSystwrt035_15Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), "hZtSystwrt035_15", nZtBin, assocZt);
		hZtSystem[1][iCen] = new TH1F(Form("hZtSystwrt040_20Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), "hZtSystwrt040_20", nZtBin, assocZt);
		hZtSystem[0][iCen] = new TH1F(Form("hZtSystwrt040_15Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), "hZtSystwrt040_15", nZtBin, assocZt);
		hZtSystFinal[iCen] = new TH1F("hZtSystFinal", "hZtSystFinal", nZtBin, assocZt);
		for (int ibin = 0; ibin < nBins; ibin++)
		{
			shSyst[0][iCen][ibin] = abs(hZt[0][iCen]->GetBinContent(ibin + 1) - hZt[2][iCen]->GetBinContent(ibin + 1));
			shSyst[0][iCen][ibin] = abs(shSyst[0][iCen][ibin] / hZt[0][iCen]->GetBinContent(ibin + 1));
			hZtSystem[0][iCen]->SetBinContent(ibin + 1, shSyst[0][iCen][ibin]);
			// hZtSystem[iCen][0]->SetBinError(ibin + 1, hZt[2][iCen]->GetBinError(ibin + 1) / hZt[2][iCen]->GetBinContent(ibin + 1));

			cout << "Zt bin [" << assocZt[ibin] << "-" << assocZt[ibin + 1] << "] : ";
			cout << "Pur 1: " << hZt[0][iCen]->GetBinContent(ibin + 1) << ", Error: " << hZt[0][iCen]->GetBinError(ibin + 1) / hZt[0][iCen]->GetBinContent(ibin + 1) << endl;
			cout << "Zt bin: " << assocZt[ibin] << "-" << assocZt[ibin + 1] << ":\t";
			cout << "Sh0.4-2.00: " << shSyst[0][iCen][ibin] << ";\t";

			shSyst[1][iCen][ibin] = abs(hZt[0][iCen]->GetBinContent(ibin + 1) - hZt[1][iCen]->GetBinContent(ibin + 1));
			shSyst[1][iCen][ibin] = abs(shSyst[1][iCen][ibin] / hZt[0][iCen]->GetBinContent(ibin + 1));
			hZtSystem[1][iCen]->SetBinContent(ibin + 1, shSyst[1][iCen][ibin]);

			shSyst[2][iCen][ibin] = abs(hZt[0][iCen]->GetBinContent(ibin + 1) - hZt[3][iCen]->GetBinContent(ibin + 1));
			shSyst[2][iCen][ibin] = abs(shSyst[2][iCen][ibin] / hZt[0][iCen]->GetBinContent(ibin + 1));
			hZtSystem[2][iCen]->SetBinContent(ibin + 1, shSyst[2][iCen][ibin]);

			cout << "Sh 0.4-1.00: " << shSyst[1][iCen][ibin] << ";\t";

			cout << "Systematics: " << TMath::Sqrt((shSyst[0][iCen][ibin] * shSyst[0][iCen][ibin] + shSyst[1][iCen][ibin] * shSyst[1][iCen][ibin]) / 2) << "_____" << ((shSyst[0][iCen][ibin] + shSyst[1][iCen][ibin]) / 2) << endl;
			hZtSystFinal[iCen]->SetBinContent(ibin + 1, (shSyst[0][iCen][ibin] + shSyst[1][iCen][ibin] + shSyst[3][iCen][ibin]) / 3);
		}
	}

	///////////////////////////////////////////////////////////////////////////////////
	////////////////// Plot the systematic for the ShSh variation /////////////////////
	////////////////// Compute the fit for polishing the trend   //////////////////////
	///////////////////////////////////////////////////////////////////////////////////

	TCanvas *cSystZt[nCen];
	TLatex *ALICEtex1[nCen];
	int kColorSyst[3] = {kCyan + 2, kOrange + 6, kAzure + 3};

	TF1 *fa1[nCen];
	TLatex *parPol1[nCen];

	for (int iCen = 0; iCen < nCen; iCen++)
	{
		cSystZt[iCen] = new TCanvas(Form("cSystZt%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("cSystZt%d_%d", cenBins[iCen], cenBins[iCen + 1]), 800, 600);
		cSystZt[iCen]->cd();
		hZtSystem[0][iCen]->SetTitle(Form("Cent %d-%d", cenBins[iCen], cenBins[iCen + 1]));
		for (int iSyst = 0; iSyst < 3; iSyst++)
		{
			hZtSystem[iSyst][iCen]->Scale(100);
			hZtSystem[iSyst][iCen]->SetDirectory(0);
			PlotStyle(hZtSystem[iSyst][iCen], 33, 2, kColorSyst[iSyst], "z_{T}", "Uncertainty %");
			hZtSystem[iSyst][iCen]->GetYaxis()->SetRangeUser(-10, 110);
			hZtSystem[iSyst][iCen]->Draw("hist same p ");
		}
		hZtSystFinal[iCen]->Scale(100);
		hZtSystFinal[iCen]->SetDirectory(0);
		PlotStyle(hZtSystFinal[iCen], 20, 2, kPink + 5, "z_{T}", "Uncertainty %");
		hZtSystFinal[iCen]->Draw("hist same p ");
		legZtSyst[iCen]->AddEntry(hZtSystem[0][iCen], "Bkg. #sigma^{2}_{long, 5#times5}: 0.40-1.50", "lp");
		legZtSyst[iCen]->AddEntry(hZtSystem[1][iCen], "Bkg. #sigma^{2}_{long, 5#times5}: 0.40-2.00", "lp");
		legZtSyst[iCen]->AddEntry(hZtSystem[2][iCen], "Bkg. #sigma^{2}_{long, 5#times5}: 0.35-1.00", "lp");
		legZtSyst[iCen]->AddEntry(hZtSystFinal[iCen], "Mean", "lp");

		ALICEtex1[iCen] = LatexStd(ALICEtex1[iCen], 0.150, 0.86, cenBins[iCen], cenBins[iCen + 1], ptMin, ptMax);

		////////////////////////////////////////////////////////////////////////////////////
		////////////////// Compute the fit for polishing the trend:   //////////////////////
		////////////////// pol0 for 0-10% and expo for 0-30%, 30-50%, 50-90%  //////////////
		////////////////////////////////////////////////////////////////////////////////////

		if (iCen == 0 && !b0_30)
		{
			fa1[iCen] = new TF1(Form("fit%d_%d", cenBins[iCen], cenBins[iCen + 1]), "pol0", 0.15, 0.60);
		}
		else
		{
			fa1[iCen] = new TF1(Form("fit%d_%d", cenBins[iCen], cenBins[iCen + 1]), "expo", 0.15, 0.60);
		}
		gStyle->SetOptStat(0);
		gStyle->SetOptFit(111);

		parPol1[iCen] = new TLatex();
		parPol1[iCen]->SetTextSize(0.04);
		parPol1[iCen]->SetTextFont(42);
		parPol1[iCen]->SetNDC();
		hZtSystFinal[iCen]->Fit(Form("fit%d_%d", cenBins[iCen], cenBins[iCen + 1]), "R");

		fa1[iCen]->SetLineColor(kAzure + 7);
		fa1[iCen]->SetLineStyle(10);
		fa1[iCen]->SetLineWidth(6);
		legZtSyst[iCen]->AddEntry(fa1[iCen], "expo ", "l");
		fa1[iCen]->Draw("same");
		// parPol1[iCen]->DrawLatex(0.18, 0.44, Form("expo = %2.2f #pm %2.2f", fa1[iCen]->GetParameter(0), fa1[iCen]->GetParError(0)));
		//  parPol1[iCen]->DrawLatex(0.18, 0.40, Form("#chi^{2}/NDF: %2.2f/%d", fa1[iCen]->GetChisquare(), fa1[iCen]->GetNDF()));
		legZtSyst[iCen]->Draw("same");
		// parPol1[iCen]->DrawLatex(0.15, 0.70, Form("par1 %f#pm%f", fa1[iCen]->GetParameter(1), fa1[iCen]->GetParError(1)));
		cSystZt[iCen]->Print(Form("%s/ShShSystCen%d_%d%s.pdf", dirPlot.Data(), cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
	}
	/////////////////////////////////////////////////////////////////////////////////////
	/////////////////////// Filling an histogram with the various  //////////////////////
	////////////////// systematic uncertainties obtained from the fit ///////////////////
	/////////////////////////////////////////////////////////////////////////////////////
	TH1F *hShShUncert[nCen];
	for (int iCen = 0; iCen < nCen; iCen++)
	{
		hShShUncert[iCen] = new TH1F(Form("hShShUncertFromFitCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("hShShUncertFromFitCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), nZtBin, assocZt);
		for (int ibin = 0; ibin < nZtBin; ibin++)
		{
			hShShUncert[iCen]->SetBinContent(ibin + 1, fa1[iCen]->Eval((hShShUncert[iCen]->GetBinCenter(ibin + 1))));
		}
		fShSyst->cd();
		fa1[iCen]->Write();
		hZtSystFinal[iCen]->Write();
		hShShUncert[iCen]->Write();
	}

	///////////////////////////////////////////////////////////////////////////
	/////////////////////// Plot Icp for different ShSh  //////////////////////
	///////////////////////////////////////////////////////////////////////////

	TCanvas *cDiffShIcp[nCen];
	TLegend *legIcp[nCen];
	TLatex *ALICEtexIcp[nCen];
	for (int iCen = 0; iCen < nCen; iCen++)
	{
		legIcp[iCen] = LegStd(legIcp[iCen], 0.15, 0.45, 0.5, 0.70);
		legIcp[iCen]->SetTextSize(0.03);
		cDiffShIcp[iCen] = new TCanvas(Form("cDiffShIcpCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("cDiffShIcpCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), 800, 600);
		cDiffShIcp[iCen]->cd();
		gPad->SetLogy();
		hZtIcp[0][iCen]->GetYaxis()->SetRangeUser(0.15, 40);
		for (int iSyst = 0; iSyst < 4; iSyst++)
		{
			hZtIcp[iSyst][iCen]->SetTitle("");
			hZtIcp[iSyst][iCen]->Draw("Same");
		}

		ALICEtexIcp[iCen] = new TLatex();
		ALICEtexIcp[iCen]->SetTextFont(42);
		ALICEtexIcp[iCen]->SetTextSize(0.04);
		ALICEtexIcp[iCen]->SetNDC();
		ALICEtexIcp[iCen]->DrawLatex(0.15, 0.85, Form("#it{This Thesis}"));
		ALICEtexIcp[iCen]->DrawLatex(0.15, 0.85 - 0.06, Form("#bf{%d-%d %% / 50-90 %%} Pb-Pb, #sqrt{s_{NN}} = 5.02 TeV", cenBins[iCen], cenBins[iCen + 1]));
		ALICEtexIcp[iCen]->DrawLatex(0.15, 0.85 - 2 * 0.06, Form("|#it{#eta}^{ #it{#gamma}}| < 0.67, %2.0f < #it{p}_{#it{T}}^{#gamma} < %2.0f GeV/#it{c}", ptMin, ptMax));
		legIcp[iCen]->AddEntry(hZtIcp[0][iCen], "Nom. Bkg. #sigma^{2}_{long, 5#times5}: 0.40-1.00");
		for (int iSyst = 0; iSyst < 3; iSyst++)
		{
			legIcp[iCen]->AddEntry(hZtIcp[iSyst + 1][iCen], Form("Bkg. #sigma^{2}_{long, 5#times5}: %s", shshLeg[iSyst + 1].Data()));
		}
		legIcp[iCen]->Draw("same");
		cDiffShIcp[iCen]->Print(Form("%s/IcpDistribXdiffShCen%d_%d%s.pdf", dirPlot.Data(), cenBins[iCen], cenBins[iCen + 1], sPtAll.Data()));
	}
	///////////////////////////////////////////////////////////////////////////
	/////////////////////// Estimate Icp systematics //////////////////////////
	///////////////////////////////////////////////////////////////////////////
	TH1F *hSystIcp[3][nCen];
	for (int iCen = 0; iCen < nCen; iCen++)
	{
		hSystIcp[0][iCen] = new TH1F(Form("hSystIcpSh%s_Cen%d_%d", shshLeg[1].Data(), cenBins[iCen], cenBins[iCen + 1]), Form("hSystIcpSh%s_Cen%d_%d", shshLeg[1].Data(), cenBins[iCen], cenBins[iCen + 1]), nZtBin, assocZt);
		hSystIcp[1][iCen] = new TH1F(Form("hSystIcpSh%s_Cen%d_%d", shshLeg[2].Data(), cenBins[iCen], cenBins[iCen + 1]), Form("hSystIcpSh%s_Cen%d_%d", shshLeg[2].Data(), cenBins[iCen], cenBins[iCen + 1]), nZtBin, assocZt);
		hSystIcp[2][iCen] = new TH1F(Form("hSystIcpSh%s_Cen%d_%d", shshLeg[3].Data(), cenBins[iCen], cenBins[iCen + 1]), Form("hSystIcpSh%s_Cen%d_%d", shshLeg[3].Data(), cenBins[iCen], cenBins[iCen + 1]), nZtBin, assocZt);
		for (int ibin = 0; ibin < nZtBin; ibin++)
		{
			hSystIcp[0][iCen]->SetBinContent(ibin + 1, abs(hZtIcp[0][iCen]->GetBinContent(ibin + 1) - hZtIcp[1][iCen]->GetBinContent(ibin + 1)) / hZtIcp[0][iCen]->GetBinContent(ibin + 1));
			hSystIcp[1][iCen]->SetBinContent(ibin + 1, abs(hZtIcp[0][iCen]->GetBinContent(ibin + 1) - hZtIcp[2][iCen]->GetBinContent(ibin + 1)) / hZtIcp[0][iCen]->GetBinContent(ibin + 1));
			hSystIcp[2][iCen]->SetBinContent(ibin + 1, abs(hZtIcp[0][iCen]->GetBinContent(ibin + 1) - hZtIcp[3][iCen]->GetBinContent(ibin + 1)) / hZtIcp[0][iCen]->GetBinContent(ibin + 1));
		}
		hSystIcp[0][iCen]->Scale(100);
		hSystIcp[1][iCen]->Scale(100);
		hSystIcp[2][iCen]->Scale(100);
	}
	//////////////////////////////////////////////////////////////////////////////////
	/////////////////////// Evaluate the systematic average and fit //////////////////
	/////////////////////////////////////////////////////////////////////////////////
	TH1F *hSystMeanIcp[nCen];
	TF1 *fa0Icp[nCen];
	TH1F *hSystUncertIcp_ShSh[nCen];
	for (int iCen = 0; iCen < nCen; iCen++)
	{
		hSystMeanIcp[iCen] = new TH1F(Form("hSystMeanIcp_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("hSystMeanIcp_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), nZtBin, assocZt);
		for (int ibin = 0; ibin < nZtBin; ibin++)
		{
			double binContMean0 = hSystIcp[0][iCen]->GetBinContent(ibin + 1);
			double binContMean1 = hSystIcp[1][iCen]->GetBinContent(ibin + 1);
			double binContMean2 = hSystIcp[2][iCen]->GetBinContent(ibin + 1);
			hSystMeanIcp[iCen]->SetBinContent(ibin + 1, (binContMean0 + binContMean1 + binContMean2) / 3);
		}
		fa0Icp[iCen] = new TF1(Form("fa0Icp_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), "expo", 0.15, 0.60);
		hSystMeanIcp[iCen]->Fit(Form("fa0Icp_Cen%d_%d", cenBins[iCen], cenBins[iCen + 1]), "R");

		hSystUncertIcp_ShSh[iCen] = new TH1F(Form("hSystUncertIcp_ShShFromFitCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), Form("hSystUncertIcp_ShShFromFitCen%d_%d", cenBins[iCen], cenBins[iCen + 1]), nZtBin, assocZt);
		for (int ibin = 0; ibin < nZtBin; ibin++)
		{
			hSystUncertIcp_ShSh[iCen]->SetBinContent(ibin + 1, fa0Icp[iCen]->Eval((hSystMeanIcp[iCen]->GetBinCenter(ibin + 1))));
		}
		fShSyst->cd();
		hSystMeanIcp[iCen]->Write();
		hSystUncertIcp_ShSh[iCen]->Write();
	}
	TCanvas *cSystIcp = new TCanvas("cSystIcp", "cSystIcp", 2 * 800, 1 * 600);
	cSystIcp->Divide(2, 1);
	TLatex *ALICEtexIcpSyst[nCen - 1];
	TLatex *parPol0Icp[nCen - 1];
	TLegend *legShShIcp[nCen - 1];
	for (int iCen = 0; iCen < nCen - 1; iCen++)
	{
		legShShIcp[iCen] = LegStd(legShShIcp[iCen], 0.15, 0.45, 0.5, 0.70);
		legShShIcp[iCen]->SetTextSize(0.03);
		PlotStyle(hSystIcp[0][iCen], 33, 2, kCyan + 2, "#it{z}_{T}", "Uncertainty %");
		PlotStyle(hSystIcp[1][iCen], 34, 2, kOrange + 6, "#it{z}_{T}", "Uncertainty %");
		PlotStyle(hSystIcp[2][iCen], 22, 2, kAzure + 3, "#it{z}_{T}", "Uncertainty %");
		PlotStyle(hSystMeanIcp[iCen], 20, 2, kPink + 5, "#it{z}_{T}", "Uncertainty %");
		fa0Icp[iCen]->SetLineStyle(10);
		fa0Icp[iCen]->SetLineStyle(10);
		fa0Icp[iCen]->SetLineColor(kAzure + 3);
		cSystIcp->cd(iCen + 1);
		hSystIcp[0][iCen]->SetTitle(" ");
		hSystIcp[0][iCen]->GetYaxis()->SetRangeUser(0, 100);
		hSystIcp[0][iCen]->Draw("same hist p");
		hSystIcp[1][iCen]->Draw("same hist p");
		hSystIcp[2][iCen]->Draw("same hist p");
		hSystMeanIcp[iCen]->Draw("same hist p");
		fa0Icp[iCen]->Draw("same hist pc");
		legShShIcp[iCen]->AddEntry(hSystIcp[0][iCen], Form("Bkg. #sigma^{2}_{long, 5#times5}: %s", shshLeg[1].Data()));
		legShShIcp[iCen]->AddEntry(hSystIcp[1][iCen], Form("Bkg. #sigma^{2}_{long, 5#times5}: %s", shshLeg[2].Data()));
		legShShIcp[iCen]->AddEntry(hSystIcp[2][iCen], Form("Bkg. #sigma^{2}_{long, 5#times5}: %s", shshLeg[3].Data()));
		legShShIcp[iCen]->AddEntry(hSystMeanIcp[iCen], Form("Mean"));
		legShShIcp[iCen]->AddEntry(fa0Icp[iCen], Form("expo"));
		legShShIcp[iCen]->Draw("same");

		ALICEtexIcpSyst[iCen] = new TLatex();
		ALICEtexIcpSyst[iCen]->SetTextFont(42);
		ALICEtexIcpSyst[iCen]->SetTextSize(0.04);
		ALICEtexIcpSyst[iCen]->SetNDC();
		ALICEtexIcpSyst[iCen]->DrawLatex(0.15, 0.85, Form("#it{This Thesis}"));
		ALICEtexIcpSyst[iCen]->DrawLatex(0.15, 0.85 - 0.06, Form("#bf{%d-%d %% / 50-90 %%} Pb-Pb, #sqrt{s_{NN}} = 5.02 TeV", cenBins[iCen], cenBins[iCen + 1]));
		ALICEtexIcpSyst[iCen]->DrawLatex(0.15, 0.85 - 2 * 0.06, Form("|#it{#eta}^{ #it{#gamma}}| < 0.67, %2.0f < #it{p}_{#it{T}}^{#gamma} < %2.0f GeV/#it{c}", ptMin, ptMax));
	}
	cSystIcp->Print(Form("%s/hSystIcp%s.pdf", dirPlot.Data(), sPtAll.Data()));
}

void PlotStyle(TH1F *hPlot, int kMarker, double kMarkerSize, int kColor, TString titleX, TString titleY)
{
	gStyle->SetTitleX(0.56);
	// gStyle->SetOptTitle(1);
	gStyle->SetOptStat(0);
	// gStyle->SetOptStat(0);
	// gStyle->SetOptFit(1111111);
	//  gStyle->SetLabelSize(.5, "XY");
	gStyle->SetLineScalePS(1);
	gStyle->SetTitleAlign(23);
	gStyle->SetPadRightMargin(0.03);
	gStyle->SetPadLeftMargin(0.12);
	hPlot->SetMarkerStyle(kMarker);
	hPlot->SetMarkerSize(1.1);
	// hPlot->SetMarkerSize(2);
	hPlot->SetMarkerColor(kColor);
	hPlot->SetLineColor(kColor);
	hPlot->SetLineWidth(3);
	hPlot->GetYaxis()->SetTitle(Form("%s", titleY.Data()));
	hPlot->GetXaxis()->SetTitle(Form("%s", titleX.Data()));
	// hPlot->GetXaxis()->SetLabelSize(.05);
	// hPlot->GetYaxis()->SetLabelSize(.05);
	//  leg->SetFillColor(kWhite);
	//  leg->SetLineColor(0);
}

TLatex *LatexStd(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax)
{
	lat = new TLatex();
	lat->SetTextFont(42);
	lat->SetTextSize(0.04);
	lat->SetNDC();

	lat->DrawLatex(xpos, ypos, Form("#it{This Thesis}"));
	lat->DrawLatex(xpos, ypos - 0.06, Form("%d-%d %% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV", cenMin, cenMax));
	lat->DrawLatex(xpos, ypos - 2 * 0.06, Form("|#it{#eta}^{ #it{#gamma}}| < 0.67, %2.0f < #it{p}_{T}^{#gamma} < %2.0f GeV/#it{c}", ptMin, ptMax));

	return lat;
}
TLegend *LegStd(TLegend *leg, double xpos1, double ypos1, double xpos2, double ypos2)
{
	leg = new TLegend(xpos1, ypos1, xpos2, ypos2);
	leg->SetFillColor(kWhite);
	leg->SetLineWidth(0);
	leg->SetTextSize(0.045);
	return leg;
}