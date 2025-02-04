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

void PlotStyle(TH1F *hPlot, int kMarker, double kMarkerSize, int kColor, int kColorFill, TString titleX, TString titleY, bool onFill)
{
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(1);
  // gStyle->SetOptFit(1111111);
gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleAlign(23);

  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  hPlot->SetMarkerStyle(kMarker);
  hPlot->SetMarkerSize(kMarkerSize);
  hPlot->SetMarkerColor(kColor);
  hPlot->SetLineColor(kColor);
  hPlot->SetLineWidth(2);
  if (onFill)
  {
    hPlot->SetFillColorAlpha(kColorFill, 0.30);
  }

  hPlot->GetYaxis()->SetTitle(Form("%s", titleY.Data()));
  hPlot->GetXaxis()->SetTitle(Form("%s", titleX.Data()));
  //hPlot->GetXaxis()->SetTitleSize(0.06);
  //hPlot->GetYaxis()->SetTitleSize(0.06);
  //hPlot->GetXaxis()->SetLabelSize(0.05);
  //hPlot->GetYaxis()->SetLabelSize(0.05);
  //hPlot->GetYaxis()->SetLabelOffset(0.025);
  hPlot->GetXaxis()->SetLabelOffset(0.02);
  hPlot->GetXaxis()->SetTickLength(0.015);
  hPlot->GetYaxis()->SetTickLength(0.02);

  // leg->SetFillColor(kWhite);
  // leg->SetLineColor(0);
}

TLatex *LatexStdphi(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax, bool bCen, TString sTitle = "#font[42]{ALICE preliminary}")
{
  lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetTextSize(0.04);
  lat->SetNDC();
  //lat->DrawLatex(xpos, ypos + 0.06, Form("%s", sTitle.Data()));
  lat->DrawLatex(xpos, ypos, "#it{This Thesis}");
  if (bCen)
    lat->DrawLatex(xpos, ypos - 0.08, Form("#bf{%d#font[122]{-}%d%%} Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, |#it{#eta}^{ trig}| < 0.67", cenMin, cenMax));
  else if (!bCen)
    lat->DrawLatex(xpos, ypos - 0.08, Form("Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, |#it{#eta}^{ trig}| < 0.67"));
  lat->DrawLatex(xpos, ypos - 2 * 0.08, Form("%2.0f < #it{p}_{T}^{ trig} < %2.0f GeV/#it{c} #otimes #it{p}_{T}^{ h} > 0.5 GeV/#it{c}", ptMin, ptMax));
  return lat;
}

TLatex *LatexStd(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax, bool bCen, TString sTitle = "#font[42]{ALICE preliminary}")
{
  lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetTextSize(0.04);
  lat->SetNDC();
  //lat->DrawLatex(xpos, ypos + 0.06, Form("%s", sTitle.Data()));
  lat->DrawLatex(xpos, ypos, "#it{This Thesis}");
  if (bCen)
    lat->DrawLatex(xpos, ypos - 0.06, Form("#bf{%d#font[122]{-}%d%%} Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, |#it{#eta}^{ trig}| < 0.67", cenMin, cenMax));
  else if (!bCen)
    lat->DrawLatex(xpos, ypos - 0.06, Form("Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, |#it{#eta}^{ trig}| < 0.67"));
  lat->DrawLatex(xpos, ypos - 2 * 0.06, Form("%2.0f < #it{p}_{T}^{ trig} < %2.0f GeV/#it{c} #otimes #it{p}_{T}^{ h} > 0.5 GeV/#it{c}", ptMin, ptMax));
  return lat;
}


TLatex *LatexStdSyst(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax, bool bCen, TString sTitle = "#font[42]{ALICE preliminary}")
{
  lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetTextSize(0.04);
  lat->SetNDC();
  //lat->DrawLatex(xpos, ypos + 0.06, Form("%s", sTitle.Data()));
  lat->DrawLatex(xpos, ypos, "#it{This Thesis}");
  if (bCen)
    lat->DrawLatex(xpos, ypos - 0.06, Form("#bf{%d#font[122]{-}%d%%} Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, ", cenMin, cenMax));
  else if (!bCen)
    lat->DrawLatex(xpos, ypos - 0.06, Form("Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, |#it{#eta}^{ trig}| < 0.67"));
  lat->DrawLatex(xpos, ypos - 2 * 0.06, Form("|#it{#eta}^{ trig}| < 0.67, %2.0f < #it{p}_{T}^{ trig} < %2.0f GeV/#it{c}", ptMin, ptMax));
  return lat;
}


TLatex *LatexStdISO(TLatex *lat, double xpos, double ypos, double texSize, int cenMin, int cenMax, float ptMin, float ptMax, bool bCen)
{
  lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetTextSize(texSize);
  lat->SetNDC();
  //lat->DrawLatex(xpos, ypos, Form("#font[42]{ALICE preliminary}"));
    lat->DrawLatex(xpos, ypos, Form("#it{This Thesis}"));
  if (bCen)
    lat->DrawLatex(xpos, ypos - 0.045, Form("#bf{%d#font[122]{-}%d%%} Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV ", cenMin, cenMax));
  else if (!bCen)
    lat->DrawLatex(xpos, ypos - 0.045, Form("Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV "));
  lat->DrawLatex(xpos, ypos - 2 * 0.048, Form("|#Delta#it{#varphi}_{#it{#gamma}#font[122]{-}h}| > #frac{3}{5} #it{#pi}, |#it{#eta}^{ #it{#gamma}}| < 0.67 "));
  lat->DrawLatex(xpos, ypos - 3 * 0.048, Form("%2.0f < #it{p}_{T}^{ #it{#gamma}} < %2.0f GeV/#it{c} #otimes #it{p}_{T}^{ h} > 0.5 GeV/#it{c} ", ptMin, ptMax));
  return lat;
}

TLatex *LatexStdISORatio(TLatex *lat, double xpos, double ypos, double texSize, int cenMin, int cenMax, float ptMin, float ptMax, bool bCen)
{
  lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetTextSize(texSize);
  lat->SetNDC();
  //lat->DrawLatex(xpos, ypos, Form("#font[42]{ALICE preliminary}"));
    lat->DrawLatex(xpos, ypos, Form("#it{This Thesis}"));
  if (bCen)
    lat->DrawLatex(xpos, ypos - 0.06, Form("#bf{%d#font[122]{-}%d%%} Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV ", cenMin, cenMax));
  else if (!bCen)
    lat->DrawLatex(xpos, ypos - 0.064, Form("Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV "));
  lat->DrawLatex(xpos, ypos - 2 * 0.064, Form("|#Delta#it{#varphi}_{#it{#gamma}#font[122]{-}h}| > #frac{3}{5} #it{#pi}, |#it{#eta}^{ #it{#gamma}}| < 0.67 "));
  lat->DrawLatex(xpos, ypos - 3 * 0.064, Form("%2.0f < #it{p}_{T}^{ #it{#gamma}} < %2.0f GeV/#it{c} #otimes #it{p}_{T}^{ h} > 0.5 GeV/#it{c} ", ptMin, ptMax));
  return lat;
}

TLegend *LegStd(TLegend *leg, double xpos1, double ypos1, double xpos2, double ypos2)
{
  leg = new TLegend(xpos1, ypos1, xpos2, ypos2);
  leg->SetFillColor(kWhite);
  leg->SetLineWidth(0);
  leg->SetTextSize(0.04);
  return leg;
}

TLatex *LatexDPhi(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax, bool bCent)
{
  lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetTextSize(0.065);
  lat->SetNDC();
  //lat->DrawLatex(xpos, ypos + 0.10, Form("#font[42]{ALICE preliminary}"));
  lat->DrawLatex(xpos, ypos, Form("#it{This Thesis}"));
  if(bCent)
    lat->DrawLatex(xpos, ypos - 0.10, Form("#bf{%d#font[122]{-}%d%%} Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, |#it{#eta}^{ trig}| < 0.67", cenMin, cenMax));
  if(!bCent)
    lat->DrawLatex(xpos, ypos - 0.10, Form("Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, |#it{#eta}^{ trig}| < 0.67"));
  lat->DrawLatex(xpos, ypos - 2 * 0.10, Form("%2.0f < #it{p}_{T}^{ trig} < %2.0f GeV/#it{c} #otimes #it{p}_{T}^{h} > 0.5 GeV/#it{c}", ptMin, ptMax));
  return lat;
}
TLatex *LatexDPhiNopt(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax, bool bCent)
{
  lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetTextSize(0.065);
  lat->SetNDC();
  //lat->DrawLatex(xpos, ypos + 0.10, Form("#font[42]{ALICE preliminary}"));
  lat->DrawLatex(xpos, ypos, Form("#it{This Thesis}"));
  if(bCent)
    lat->DrawLatex(xpos, ypos - 0.10, Form("#bf{%d#font[122]{-}%d%%} Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, |#it{#eta}^{ trig}| < 0.67", cenMin, cenMax));
  if(!bCent)
    lat->DrawLatex(xpos, ypos - 0.10, Form("Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV, |#it{#eta}^{ trig}| < 0.67"));
  lat->DrawLatex(xpos, ypos - 2 * 0.10, Form("#it{p}_{T}^{ h} > 0.5 GeV/#it{c}"));
  return lat;
}

TGraph *DrawLine(TGraph *line, double x1, double y1, double x2, double y2)
{
  line = new TGraph(2);
  line->SetPoint(0, x1, y1);
  line->SetPoint(1, x2, y2);

  line->SetLineStyle(3);
  line->SetLineColor(kBlack);
  line->SetLineWidth(3);

  return line;
}

TLatex *LatexStdIcp(TLatex *lat, double xpos, double ypos, int cenMin, int cenMax, float ptMin, float ptMax)
{
  lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetTextSize(0.04);
  lat->SetNDC();
  //lat->DrawLatex(xpos, ypos, Form("#font[42]{ALICE preliminary}"));
  lat->DrawLatex(xpos, ypos, Form("#it{This Thesis}"));
  lat->DrawLatex(xpos, ypos - 0.055, Form("Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV "));
  lat->DrawLatex(xpos, ypos - 2 * 0.06, Form("|#Delta#it{#varphi}_{#it{#gamma}#font[122]{-}h}| > #frac{3}{5} #it{#pi}, |#it{#eta}^{ #it{#gamma}}| < 0.67 "));
  lat->DrawLatex(xpos, ypos - 3 * 0.06, Form("%2.0f < #it{p}_{T}^{ #it{#gamma}} < %2.0f GeV/#it{c} #otimes #it{p}_{T}^{ h} > 0.5 GeV/#it{c}", ptMin, ptMax));
  return lat;
}

TCanvas *canvasStd(TString name, int xPad, int yPad)
{
  TCanvas *canvas = new TCanvas(name, name, xPad * 800, yPad * 600);
  canvas->Divide(xPad, yPad);
  for (int iPad = 0; iPad < xPad * yPad; iPad++)
  {
    canvas->cd(iPad + 1);
    // canvas->GetPad(iPad + 1)->SetTopMargin(0.015);
    // canvas->GetPad(iPad + 1)->SetRightMargin(0.05);
    // canvas->GetPad(iPad + 1)->SetLeftMargin(0.2);
    // canvas->GetPad(iPad + 1)->SetBottomMargin(0.11);
    //gStyle->SetPadRightMargin(0.05);
    //gStyle->SetPadLeftMargin(0.20);
    //gStyle->SetPadBottomMargin(0.15);
    //gStyle->SetTitleX(0.56);
  }

  return canvas;
}
//-----------------------------------------------------------------------------
/// Divide 2 TGraphErrors
///
/// \return TGraphError result of division
/// \param  gNum TGraphError numerator
/// \param  gDen TGraphError denominator
//-----------------------------------------------------------------------------
static TGraphAsymmErrors * DivideGraphs(TGraphAsymmErrors* gNum, TGraphAsymmErrors *gDen)
{
  if ( !gDen || !gNum ) 
  {
    printf("Graph num %p or graph den %p not available\n",gNum,gDen);
    return 0x0;
  }
  const Int_t nBins = gNum->GetN();
  if ( nBins != gDen->GetN() )
  {
    printf("Cannot divide %s with %d bins and %s with %d bins!\n",
           gNum->GetName(),nBins,gDen->GetName(),gDen->GetN());
    return 0x0;
  }
  
  Double_t ratio   [nBins];
  Double_t ratioErr[nBins]; 
  Double_t x       [nBins];
  Double_t xErr    [nBins];
  for (Int_t ibin = 0; ibin < nBins; ibin++) 
  {
    Double_t num    =  gNum->GetY ()[ibin];
    Double_t den    =  gDen->GetY ()[ibin];
    Double_t numErr =  gNum->GetEY()[ibin];
    Double_t denErr =  gDen->GetEY()[ibin];
    
    x   [ibin]      =  gNum->GetX ()[ibin];
    xErr[ibin]      =  gNum->GetEX()[ibin];
    
    if ( num == 0 || den == 0 ) 
    {
      ratio   [ibin] = 0; 
      ratioErr[ibin] = 0; 
      continue;
    }
    
    ratio   [ibin] = num / den ;
    //ratioErr[ibin] =  GetFractionError(num,den,numErr,denErr);  
    
    //    printf("bin %d, x %f (%f) num %f (%f), den %f (%f), ratio %f (%f) \n",
    //           ibin,x[ibin],xErr[ibin],num,numErr,den,denErr,ratio[ibin],ratioErr[ibin]);
  } // do the ratio to sum
  
  return new TGraphAsymmErrors(nBins,x,ratio,xErr,ratioErr);
}

