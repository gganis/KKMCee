//////////////////////////////////////////////////////////////////////
//    make PlotFoam0-run
//////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
using namespace std;

#include <math.h>
#include <TLorentzVector.h>
#include <TLine.h>
#include <TArrow.h>
#include <TLatex.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TH2.h"
#include "TGaxis.h"
#include "TApplication.h"
#include "TMarker.h"
#include "TObjString.h"
#include "TFile.h"

#include "HisReMake.h"
#include "HisNorm.h"
#include "KKplot.h"

//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT
//=============================================================================
////  *** KKMC
TFile DiskFileA("../workKKMC/histo.root");
// Sept. 2017 runs
//TFile DiskFileA("../workKKMC/histo.root_95GeV_55M"); // It is 550M!!!
// August2017 runs
//TFile DiskFileA("../workKKMC/histo.root_10GeV_1G"); //
//TFile DiskFileA("../workKKMC/histo.root_88GeV_2.1G"); //
//TFile DiskFileA("../workKKMC/histo.root_95GeV_16G");
//TFile DiskFileA("../workKKMC/histo.root_91GeV_9G"); ????
// July2017 runs
//TFile DiskFileA("../workKKMC/histo.root_91GeV_6G"); //
//TFile DiskFileA("../workKKMC/histo.root_88GeV_4G"); //
//TFile DiskFileA("../workKKMC/histo.root_10GeV_5.7G"); //
//TFile DiskFileA("../workKKMC/histo.root_95GeV.4G");   //

////  *** FOAM
TFile DiskFileF("../workFOAM/histo.root"); // current
// Sept. 2017 runs
//TFile DiskFileF("../workFOAM/histo.root_95GeV_14G");
// August2017 runs
//TFile DiskFileF("../workFOAM/histo.root_10GeV_37G_vmax0.2");
//TFile DiskFileF("../workFOAM/histo.root_88GeV_16G");
//TFile DiskFileF("../workFOAM/histo.root_91GeV_45G");
//TFile DiskFileF("../workFOAM/histo.root_95GeV_10G");
//TFile DiskFileF("../workFOAM/histo.root_10GeV_32G");
//
//TFile DiskFileF("../workFOAM/histo.root_91GeV_35G");
//TFile DiskFileF("../workFOAM/histo.root_88GeV_32G");
//TFile DiskFileF("../workFOAM/histo.root_10GeV_15G");
//TFile DiskFileF("../workFOAM/histo.root_95GeV_15G");

//  ****** older FOAM and KKMC progs ******
//TFile DiskFileF("../workFoam0/rmain.root");
//TFile DiskFileF("../workFoam0rmain_95GeV_64M.root");
//
//TString XparFile="../workKKMC/workKKMC_95GeV.input"; // obsolete?
//
//TFile DiskFileA("../workAFB/rmain.root");
//TFile DiskFileA("../workAFB/rmain_95GeV.root"); // 100M new
//TFile DiskFileA("../workAFB/rmain.root_189GeV_100M"); // obsolete

TFile DiskFileB("RhoSemi.root","RECREATE","histograms");

// Interface to KKabox and some extra plotting facilities


///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
double gCMSene, gNevTot, gNevTot2; // from KKMC and KKfoam MC runs (histograms)
char   gTextEne[100], gTextNev[100], gTextNev2[100];
//
//int    gNbMax=45;         // gCosTheta = 45/50=0.90
//double gCosTheta=0.90;    // to be synchronized with gNbMax
//
int    gNbMax=50;         // gCosTheta = 45/50=0.90
double gCosTheta=1.00;    // to be synchronized with gNbMax
//
float  gXcanv = 50, gYcanv = 50;
KKplot LibSem("KKplot");
///////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////
void PlotSame(TH1D *HST, double &ycapt, Int_t kolor, TString opis)
{
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  HST->SetLineColor(kolor);
  HST->DrawCopy("hsame");      // Magenta
  CaptT->SetTextColor(kolor);
  ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt, opis);
}// PlotSame


///////////////////////////////////////////////////////////////////////////////////
void PlotSame2(TH1D *HST, double &ycapt, Int_t kolor, double xx,  TString label,  TString opis)
{
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  HST->SetLineColor(kolor);
  HST->DrawCopy("hsame");      // Magenta
  CaptT->SetTextColor(kolor);
  ycapt += -0.04;
  double xcapt = 0.40;
  CaptT->DrawLatex(xcapt,ycapt, opis);
  CaptT->DrawLatex(xcapt-0.05,ycapt, label);
  //
  TLatex *CaptS = new TLatex();
  CaptS->SetTextSize(0.040);
  CaptS->SetTextAlign(21);
  CaptS->SetTextColor(kolor);
  int ib = HST->FindBin(xx);
  double yy= HST->GetBinContent(ib);
  CaptS->DrawLatex(xx,yy,label);
}// PlotSame2



void TestNorm(){
// testing/debugging normalization
cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
cout<<"<<<<<<<<<<<<<<<<<<<<<<< TestNorm >>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;

TH1D *HST_FOAM_NORMA3 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA3");
TH1D *HST_tx_Ceex2n   = (TH1D*)DiskFileF.Get("HST_tx_Ceex2n");  // FOAM testing norm.

//int  NevTot = 4.5e6;
int NevTot  = HST_FOAM_NORMA3->GetEntries();
//NevTot = 4.5e6;
double Xsav = HST_FOAM_NORMA3->GetBinContent(0)/NevTot; // NANOBARNS
cout<<"TestNorm: Xsav = "<<Xsav<<"  NevTot =  "<< NevTot<<endl;
// first 3 bins only for vmax=0.006
double SWT =0;
double SSWT =0;
//for( int i=0; i <= 3; i++){ // 3 bins v<0.006
//for( int i=0; i <= 3; i++){ // 100 bins v<0.200
for( int i=0; i <= 3; i++){ // 100 bins v<0.020
  SWT  += HST_tx_Ceex2n->GetBinContent(i);
  SSWT += sqr(HST_tx_Ceex2n->GetBinError(i));
}//for

double AWT = SWT/NevTot;
double Sigma = sqrt(SSWT/NevTot-sqr(AWT));
double Error = Sigma/sqrt(NevTot);
double Errel = Error/AWT;

double Sigm0 = sqrt(SSWT/NevTot);       // approximation used in ROOT
double Erre0 = Sigm0/sqrt(NevTot)/AWT;  // approximation used in ROOT

cout<<" TestNorm:  AWT="<< AWT << endl;
cout<<" TestNorm:  Errel="<< Errel << endl;
cout<<" TestNorm:  Erre0="<< Erre0 << " ROOT approx. " <<endl;
cout<<" TestNorm:  ratio="<< Erre0/Errel<<endl;
cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;

}//TestNorm

///////////////////////////////////////////////////////////////////////////////////
void KKsemMakeHisto(){
  // Here we produce semianalytical plots using KKsem program, No plotting
  //------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ KKsem MakeHisto  BEGIN ============================"<<endl;
  //
  long KF=13; // muon
  long KeyDis, KeyFob;
  char chak[5];
  KeyFob=    0; // With EW (BornV_Dizet) With integration OK!
//------------------------------------------------------------------------
  TH1D *hstVtemplate = (TH1D*)DiskFileB.Get("HTot2_vTcPR_Ceex2");
//------------------------------------------------------------------------
//   MuMu  Sigma(vmax) and AFB(vmax) with ulimited c=cos(theta)
//
  KeyDis = 302302;        // ISR*FSR O(alf2)
  sprintf(chak,"XCHI2");  // ISR*FSR Mff
  TH1D *vcum_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("vcum_ISR2_FSR2");
  TH1D *afbv_ISR2_FSR2 =(TH1D*)hstVtemplate->Clone("afbv_ISR2_FSR2");
  LibSem.VVmake( vcum_ISR2_FSR2, afbv_ISR2_FSR2, KF, chak, KeyDis, KeyFob, gCosTheta);
  //
  KeyDis = 300300;        // ISR*FSR O(alf2)
  sprintf(chak,"XCHI2");  // ISR*FSR Mff
  TH1D *vcum_ISR0_FSR0 =(TH1D*)hstVtemplate->Clone("vcum_ISR0_FSR0");
  TH1D *afbv_ISR0_FSR0 =(TH1D*)hstVtemplate->Clone("afbv_ISR0_FSR0");
  LibSem.VVmake( vcum_ISR0_FSR0, afbv_ISR0_FSR0, KF, chak, KeyDis, KeyFob, gCosTheta);
  //
  cout<<"================ KKsem MakeHisto ENDs ============================="<<endl;
  cout<<"==================================================================="<<endl;
//------------------------------------------------------------------------
//------------------------------------------------------------------------
}//  KKsemMakeHisto


///////////////////////////////////////////////////////////////////////////////////
//  !!!!! OBSOLETE !!!!!
void FigVdist()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigVdist =========================== "<<endl;
 //
  TH1D *Hpro_vT_Ceex2n    = (TH1D*)DiskFileB.Get("Hpro_vT_Ceex2n");     // KKMC dsigma/dv IFI off, from scat.
  TH1D *HTot_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HTot_vTcPR_Ceex2n");  // KKMC sigma(vmax) from scat.
  //
  TH1D *HST_xx_Ceex2n     = (TH1D*)DiskFileF.Get("HST_xx_Ceex2n");      // Foam dsigma/d(v) direct
  TH1D *HST_xmax_Ceex2n   = (TH1D*)DiskFileB.Get("HST_xmax_Ceex2n");    // Foam sigma(vmax) direct
    //
  TH1D *Htot_xmax_Ceex2n  = (TH1D*)DiskFileB.Get("Htot_xmax_Ceex2n");    //Foam3 sigma(vmax) scatt.ISR+FSR
  TH1D *Htot_xmax_Ceex2   = (TH1D*)DiskFileB.Get("Htot_xmax_Ceex2");     //Foam5 sigma(vmax) scatt.ISR+FSR+IFI
  //
//  TH1D *vcum_ISR2_FSR2    = (TH1D*)DiskFileB.Get("vcum_ISR2_FSR2");     // KKsem  sigma(vmax)
//  TH1D *vdis_ISR2_FSR2    = (TH1D*)DiskFileB.Get("vdis_ISR2_FSR2");     // KKsem  dsigma/d(v)
  //
  //*****************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigVdist = new TCanvas("cFigVdist","FigVdist", gXcanv, gYcanv,    1200, 800);
  //                           Name    Title               xoff,   yoff,    WidPix,HeiPix
  gXcanv += 100; gYcanv += 50;
  cFigVdist->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigVdist->Divide( 2,  2);
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  //====================plot1========================
  //                sigma(vmax)
  cFigVdist->cd(1);
  //gPad->SetLogy(); // !!!!!!
  // MC v-true direct
  TH1D *Hst1 = HTot_vTcPR_Ceex2n;  //  KKMC sigma(vmax) from scat.
  //
  Hst1->SetStats(0);
  Hst1->SetTitle(0);
  //Hst1->SetMinimum(1e-3*Hst1->GetMaximum());
  Hst1->DrawCopy("h");
  //
  HST_xmax_Ceex2n->SetLineColor(kRed);     // red
  HST_xmax_Ceex2n->DrawCopy("hsame");      // Foam sigma(vmax) direct.
  //
  Htot_xmax_Ceex2n->SetLineColor(kBlue);   // blue ISR+FSR
  Htot_xmax_Ceex2n->DrawCopy("hsame");     // Foam sigma(vmax) scatt.
  //
  Htot_xmax_Ceex2->SetLineColor(kGreen);   // green ISR+FSR+IFI
  Htot_xmax_Ceex2->DrawCopy("hsame");      // Foam sigma(vmax) scatt.
  //
  //vcum_ISR2_FSR2->SetLineColor(kBlue);   // blue
  //vcum_ISR2_FSR2->DrawCopy("hsame");     // KKsem sigma(vmax)
  //
  CaptT->DrawLatex(0.02,0.95, "#sigma(v_{max}) ISR+FSR, Black KKMC_CEEX2, Blue FOAM");
  CaptT->DrawLatex(0.60,0.75,gTextEne);

  //====================plot2========================
  cFigVdist->cd(2);
  TH1D *Hst1_ratio =(TH1D*)Hst1->Clone("Hst1_ratio");
  //Hst1_ratio->Divide(vcum_ISR2_FSR2);   // divide by KKsem
  //Hst1_ratio->Divide(HST_xmax_Ceex2n);  // divide by Foam direct
  Hst1_ratio->Divide(Htot_xmax_Ceex2n);   // divide by Foam scatt.
  //Hst1_ratio->SetMinimum(0.95);
  //Hst1_ratio->SetMaximum(1.10);
  Hst1_ratio->SetLineColor(kBlue);
  Hst1_ratio->DrawCopy("h");
  //
  CaptT->DrawLatex(0.02,0.95,"#sigma(v_{max}) ISR+FSR, Ratio KKMC/FOAM");
  //====================plot3========================
  //                 dsigma/d(v)
  cFigVdist->cd(3);
  gPad->SetLogy(); // !!!!!!
  TH1D *Hst3 = Hpro_vT_Ceex2n;           // KKMC dsigma/dv from scat.
  Hst3->SetStats(0);
  Hst3->SetTitle(0);
  Hst3->SetLineColor(kRed); // red
  Hst3->DrawCopy("h");
  //
  // KKsem ISR+FSR
  //vdis_ISR2_FSR2->SetLineColor(kMagenta); // magenta
  //vdis_ISR2_FSR2->DrawCopy("hsame");      // KKsem dsigma/d(v) ISR+FSR

  HST_xx_Ceex2n->SetLineColor(kBlue);     // blue
  HST_xx_Ceex2n->DrawCopy("hsame");       // Foam dsigma/d(v)

  CaptT->DrawLatex(0.02,0.95,"d#sigma/dv(ISR+FSR),  Red KKMC_CEEX2, Blue FOAM");
  //====================plot4========================
  cFigVdist->cd(4);

  TH1D *Hst3_ratio =(TH1D*)Hst3->Clone("Hst3_ratio");
  //Hst3_ratio->Divide(vdis_ISR2_FSR2);
  Hst3_ratio->Divide(HST_xx_Ceex2n);  // direct

  Hst3_ratio->SetStats(0);
  Hst3_ratio->SetTitle(0);
//  Hst3_ratio->SetMinimum(0.85);
//  Hst3_ratio->SetMaximum(1.15);
//  Hst3_ratio->SetLineColor(kRed);
  Hst3_ratio->DrawCopy("h");  // black

  CaptT->DrawLatex(0.02,0.95,"d#sigma/dv(ISR+FSR ); Ratio KKMC/FOAM");
  CaptT->DrawLatex(0.60,0.75,gTextEne);

  //----------------------------
  cFigVdist->cd();
  //================================================
}//FigVdist



///////////////////////////////////////////////////////////////////////////////////
void FigAfb2()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigAfb2 =========================== "<<endl;

//  TH1D *HAfb2_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_Ceex2n");  //
//  TH1D *HAfb2_vTcPR_Ceex2  = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_Ceex2");  //
  //
  TH1D *HAfb2_vTcPL_Ceex2n = (TH1D*)DiskFileB.Get("HAfb2_vTcPL_Ceex2n");  //
  TH1D *HAfb2_vTcPL_Ceex2  = (TH1D*)DiskFileB.Get("HAfb2_vTcPL_Ceex2");  //
  //
  TH1D *Hafb2_xmax_Ceex2n  = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex2n");  // FOAM scatt.
  TH1D *Hafb2_xmax_Ceex2   = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex2");   // FOAM scatt.

  TH1D *afbv_ISR2_FSR2    = (TH1D*)DiskFileB.Get("afbv_ISR2_FSR2");    // KKsem


  TH1D *HST_PLBZ2 =(TH1D*)HAfb2_vTcPL_Ceex2->Clone("HST_PLBZ2");
  LibSem.Ord1fill(HST_PLBZ2,101);  // PRD41 IFIoff
  HST_PLBZ2->SetLineColor(kRed);
  //
  TH1D *HST_IFI4 =(TH1D*)HAfb2_vTcPL_Ceex2->Clone("HST_IFI4");
  LibSem.Ord1fill(HST_IFI4,105);   // PLB219 hard IFI
  HST_IFI4->SetLineColor(kCyan);
//
  //*****************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigAfb2 = new TCanvas("cFigAfb2","FigAfb2", gXcanv, gYcanv,   1200, 550);
  //                                 Name    Title      xoff,yoff, WidPix,HeiPix
  gXcanv += 100; gYcanv += 50;
  cFigAfb2->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigAfb2->Divide( 2,  0);
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  //====================plot1========================
  //                AFB(vmax)
  cFigAfb2->cd(1);
  //gPad->SetLogy(); // !!!!!!
  //TH1D *Hst1 = HAfb2_vTcPR_Ceex2n;       // KKMC AFB(vmax) from scat. IFI off
  //TH1D *Hst2 = HAfb2_vTcPR_Ceex2;        // KKMC AFB(vmax) from scat. IFI on
  TH1D *Hst1 = HAfb2_vTcPL_Ceex2n;         // KKMC AFB(vmax) from scat. IFI off
  TH1D *Hst2 = HAfb2_vTcPL_Ceex2;          // KKMC AFB(vmax) from scat. IFI on
  //
  Hst2->SetStats(0);
  Hst2->SetTitle(0);
  Hst2->SetLineColor(kMagenta);            // magenta

  Hst2->SetMinimum(-1);  //
  Hst2->SetMaximum( 1);  //

  if( fabs(gCMSene-94e0) <1.0 ) { Hst2->SetMinimum( 0.15); Hst2->SetMaximum( 0.35);}
  if( fabs(gCMSene-10e0) <1.0 ) { Hst2->SetMinimum(-0.02); Hst2->SetMaximum( 0.08);}

  Hst2->GetXaxis()->SetTitle("v_{max}");
  Hst2->DrawCopy("h");                     // KKMC AFB(vmax) from scat. IFI on

  CaptT->DrawLatex(0.06,0.95, "A_{FB}(v_{max}) ");
  double ycapt = 0.90; // starting value, to be decremented below
  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextEne);

  PlotSame2(Hst2,             ycapt, kMagenta,   0.015, "(a)", "KKMC.2   IFIon ");
  PlotSame2(Hafb2_xmax_Ceex2, ycapt, kGreen,     0.025, "(b)", "Foam5.2  IFIon ");
  PlotSame2(afbv_ISR2_FSR2,   ycapt, kRed,       0.120, "(c)", "KKsem.2  IFIoff ");
  PlotSame2(Hst1,             ycapt, kBlack,     0.135, "(d)", "KKMC.2   IFIoff ");
  PlotSame2(Hafb2_xmax_Ceex2n,ycapt, kBlue,      0.150, "(e)", "Foam3.2  IFIoff ");
  PlotSame2(HST_PLBZ2,        ycapt, kCyan,      0.100, "(f)", "PRD43.1  IFIoff ");

  //====================plot2========================
  cFigAfb2->cd(2);
  TH1D *Hst21_diff   = HstDiff("Hst21_diff",    HAfb2_vTcPL_Ceex2,  HAfb2_vTcPL_Ceex2n,  kBlack);
  TH1D *HST21_diff   = HstDiff("HST21_diff",    Hafb2_xmax_Ceex2,   Hafb2_xmax_Ceex2n,   kMagenta);
//  TH1D *HstKF_diff   = HstDiff("HstKF_diff",    HAfb2_vTcPL_Ceex2,  Hafb2_xmax_Ceex2,    kGreen);
  TH1D *HstPL_diff   = HstDiff("HstPL_diff",    HAfb2_vTcPL_Ceex2,  Hafb2_xmax_Ceex2,    kRed);
  TH1D *HstKFn_diff  = HstDiff("HstKFn_diff",   HAfb2_vTcPL_Ceex2n,  Hafb2_xmax_Ceex2n,  kBlue);

  Hst21_diff->SetMinimum(-0.004); Hst21_diff->SetMaximum( 0.004);  // zoom
  Hst21_diff->GetXaxis()->SetTitle("v_{max}");
  Hst21_diff->DrawCopy("h");

  CaptT->DrawLatex(0.06,0.95, "A_{FB}(v_{max}) ");
  ycapt = 0.99; // starting value, to be decremented below
  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextEne);

  PlotSame2(Hst21_diff,   ycapt, kBlack,     0.120, "(a)", "KKMC.2,   IFIon - IFIoff ");
  PlotSame2(HST21_diff,   ycapt, kMagenta,   0.140, "(b)", "Foam5.2   IFIon - IFIoff ");
  PlotSame2(HstPL_diff,   ycapt, kRed,       0.100, "(c)", "KKMC.2_IFIon  - Foam5.2_IFIon ");
  PlotSame2(HstKFn_diff,  ycapt, kBlue,      0.160, "(d)", "KKMC.2_IFIoff - Foam3.2_IFIoff ");
  PlotSame2(HST_IFI4,     ycapt, kCyan,      0.020, "(e)", "PLB219.1 IFI hard part ");


// zero line
  TH1D *hZero = (TH1D*)Hst1->Clone("hZero");  // zero line
  for(int i=1; i <= hZero->GetNbinsX() ; i++) { hZero->SetBinContent(i, 0); hZero->SetBinError(i, 0);}
  hZero->DrawCopy("hsame");

  cFigAfb2->cd();
  //================================================
}//FigAfb2



///////////////////////////////////////////////////////////////////////////////////
void FigTech0()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigTech0 IFI off =========================== "<<endl;

// sigma(vmax)
  TH1D *HTot2_vTcPL_Ceex0n = (TH1D*)DiskFileB.Get("HTot2_vTcPL_Ceex0n"); // KKMC sigma(vmax) from scat.
  TH1D *HAfb2_vTcPL_Ceex0n = (TH1D*)DiskFileB.Get("HAfb2_vTcPL_Ceex0n"); // KKMC AFB(vmax)   from scat.

  TH1D *vcum_ISR0_FSR0     = (TH1D*)DiskFileB.Get("vcum_ISR0_FSR0");     // KKsem
  TH1D *afbv_ISR0_FSR0     = (TH1D*)DiskFileB.Get("afbv_ISR0_FSR0");     // KKsem

  TH1D *Htot2_xmax_Ceex0n  = (TH1D*)DiskFileB.Get("Htot2_xmax_Ceex0n");  // FOAM3 scatt. GPS Born
  TH1D *Hafb2_xmax_Ceex0n  = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex0n");  // FOAM3 scatt. GPS Born

  TH1D *Htot2_xmax_EEX0    = (TH1D*)DiskFileB.Get("Htot2_xmax_EEX0");  // FOAM3 scatt. EEX Born
  TH1D *Hafb2_xmax_EEX0    = (TH1D*)DiskFileB.Get("Hafb2_xmax_EEX0");  // FOAM3 scatt. GPS Born

  TH1D *Htot2_xmax_Ceex0  = (TH1D*)DiskFileB.Get("Htot2_xmax_Ceex0");  // FOAM5 scatt. GPS Born
  TH1D *Hafb2_xmax_Ceex0  = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex0");  // FOAM5 scatt. GPS Born

  TH1D *hOne2 = (TH1D*)vcum_ISR0_FSR0->Clone("hOne2");  // unity line
  for(int i=1; i <= hOne2->GetNbinsX() ; i++) { hOne2->SetBinContent(i, 1); hOne2->SetBinError(i, 0);}
  hOne2->SetLineColor(kBlack);

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  //*****************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigTech0 = new TCanvas("cFigTech0","FigTech0", gXcanv, gYcanv,   1200, 550);
  //                                 Name    Title      xoff,    yoff, WidPix,HeiPix
  gXcanv += 100; gYcanv += 50;
  cFigTech0->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigTech0->Divide( 2,  0);
  //====================plot1========================
  cFigTech0->cd(1);
  TH1D *HstTech0_ratio  = HstRatio("HstTech0_ratio",   HTot2_vTcPL_Ceex0n, vcum_ISR0_FSR0, kBlack);
  TH1D *HstTech0_ratio2 = HstRatio("HstTech0_ratio2",  Htot2_xmax_Ceex0n,  vcum_ISR0_FSR0, kGreen);
  TH1D *HstTech0_ratio3 = HstRatio("HstTech0_ratio3",  Htot2_xmax_Ceex0,   vcum_ISR0_FSR0, kGreen);
  TH1D *HstTech0_ratio4 = HstRatio("HstTech0_ratio4",  Htot2_xmax_EEX0,    vcum_ISR0_FSR0, kRed);

  HstTech0_ratio->SetStats(0);
  HstTech0_ratio->SetTitle(0);
  HstTech0_ratio->SetMinimum(1 -0.0007); HstTech0_ratio->SetMaximum(1 +0.0007);  // zoom
//  HstTech0_ratio->SetMinimum(1 -0.007); HstTech0_ratio->SetMaximum(1 +0.007);  // zoom
  HstTech0_ratio->GetXaxis()->SetTitle("v_{max}");
  HstTech0_ratio->DrawCopy("h");

  double ycapt = 0.40; // starting value, to be decremented below
  PlotSame2(HstTech0_ratio,  ycapt, kBlack,    0.03, "(a)", "KKMC.0/KKsem.0  IFIoff");
  PlotSame2(HstTech0_ratio2, ycapt, kGreen,    0.06, "(b)", "FOAM3.0_gps/KKsem.0 IFIoff");
  PlotSame2(HstTech0_ratio4, ycapt, kMagenta,  0.09, "(c)", "FOAM3.0_eex/KKsem.0 IFIoff");
//  PlotSame2(HstTech0_ratio3, ycapt, kBlue,     0.16, "(d)", "FOAM5.0/KKsem.0 IFIon");

  hOne2->DrawCopy("hsame");
  CaptT->DrawLatex(0.00,0.96,"#sigma^{IFIoff}(v_{max}) ");

  //====================plot2========================
  cFigTech0->cd(2);
  TH1D *HstTech0_diff2  = HstDiff("HstTech0_diff2",   HAfb2_vTcPL_Ceex0n, afbv_ISR0_FSR0, kBlack);
  TH1D *HstTech0_diff1  = HstDiff("HstTech0_diff1",   Hafb2_xmax_Ceex0n,  afbv_ISR0_FSR0, kGreen);
  TH1D *HstTech0_diff3  = HstDiff("HstTech0_diff3",   Hafb2_xmax_EEX0,    afbv_ISR0_FSR0, kMagenta);

  TH1D *HstTech0_diff= HstTech0_diff2;
  HstTech0_diff->SetStats(0); HstTech0_diff->SetTitle(0);
  HstTech0_diff->SetMinimum(-0.0004); HstTech0_diff->SetMaximum( 0.0004);  // zoom

  HstTech0_diff->GetXaxis()->SetTitle("v_{max}");
  HstTech0_diff->DrawCopy("h");

  ycapt = 0.90; // starting value, to be decremented below
  PlotSame2(HstTech0_diff2,   ycapt, kBlack,    0.140, "(a)", "KKMC.0  - KKsem.0 IFIoff ");
  PlotSame2(HstTech0_diff1,   ycapt, kGreen,    0.180, "(b)", "Foam3.0_gps - KKsem.0 IFIoff ");
  PlotSame2(HstTech0_diff3,   ycapt, kMagenta,  0.080, "(c)", "Foam3.0_eex - KKsem.0 IFIoff ");

  TH1D *hZero2 = (TH1D*)HstTech0_diff->Clone("hZero2");  // unity line
  for(int i=1; i <= hZero2->GetNbinsX() ; i++) { hZero2->SetBinContent(i, 0); hZero2->SetBinError(i, 0);}
  hZero2->SetLineColor(kRed);
  hZero2->DrawCopy("hsame");

  CaptT->DrawLatex(0.00,0.96,"A_{FB}^{IFIoff}(v_{max})");
  ycapt =0.30;
  CaptT->DrawLatex(0.40,ycapt,gTextEne);  ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextNev);  ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextNev2); ycapt += -0.04;

  //================================================
}//FigTech0


///////////////////////////////////////////////////////////////////////////////////
void FigTech()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigTech IFI off =========================== "<<endl;

// sigma(vmax)
  TH1D *HTot2_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HTot2_vTcPR_Ceex2n"); // KKMC sigma(vmax) from scat.
  TH1D *HTot2_vTcPR_EEX2   = (TH1D*)DiskFileB.Get("HTot2_vTcPR_EEX2");   // KKMC sigma(vmax) from scat.
  TH1D *Htot2_xmax_Ceex2n  = (TH1D*)DiskFileB.Get("Htot2_xmax_Ceex2n");  // FOAM3 scatt. GPS Born
  TH1D *Htot2_xmax_EEX2    = (TH1D*)DiskFileB.Get("Htot2_xmax_EEX2");    // FOAM3 scatt. EEX Born
  TH1D *vcum_ISR2_FSR2     = (TH1D*)DiskFileB.Get("vcum_ISR2_FSR2");     // KKsem
// AFB(vmax)
  TH1D *HAfb2_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_Ceex2n"); // KKMC AFB(vmax) from scat
  TH1D *HAfb2_vTcPR_EEX2   = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_EEX2");   // KKMC AFB(vmax) from scat
  TH1D *Hafb2_xmax_Ceex2n  = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex2n");  // AFB FOAM3 scatt. GPS Born
  TH1D *Hafb2_xmax_EEX2    = (TH1D*)DiskFileB.Get("Hafb2_xmax_EEX2");    // AFB FOAM3 scatt. EEX Born
  TH1D *afbv_ISR2_FSR2     = (TH1D*)DiskFileB.Get("afbv_ISR2_FSR2");     // AFB KKsem

  TH1D *HST_txmax_Ceex2n   = (TH1D*)DiskFileB.Get("HST_txmax_Ceex2n");   // FOAM direct

  /*
  // correcting missing phase space beta factor in GPS Born
  double Mmu = 0.105, bin,vv, beta;
  TH1D *HST_bad = Htot2_xmax_EEX2;
  int Nbin    = HST_bad->GetNbinsX();
  double vmax = HST_bad->GetXaxis()->GetXmax();
  for(int i=1; i <= Nbin ; i++) {
	  bin= HST_bad->GetBinContent(i);
	  vv = (i*vmax)/Nbin;
	  beta = sqrt(1-4*sqr(Mmu/gCMSene)/(1-vv) );
	  cout<< " vv, beta ="<< vv << "   "<<beta<<endl;
	  HST_bad->SetBinContent(i, bin*beta);
  }*/

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  //*****************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigTech = new TCanvas("cFigTech","FigTech", gXcanv, gYcanv,   1200, 550);
  //                                 Name    Title      xoff,    yoff, WidPix,HeiPix
  gXcanv += 100; gYcanv += 50;
  cFigTech->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigTech->Divide( 2,  0);
  //====================plot1========================
  cFigTech->cd(1);
  TH1D *HstTech_ratio0  = HstRatio("HstTech_ratio0",  Htot2_xmax_EEX2,    vcum_ISR2_FSR2, kGreen);
  TH1D *HstTech_ratio1  = HstRatio("HstTech_ratio1",  Htot2_xmax_Ceex2n,  vcum_ISR2_FSR2, kGreen);
  TH1D *HstTech_ratio3  = HstRatio("HstTech_ratio3",  HST_txmax_Ceex2n,   vcum_ISR2_FSR2, kRed);
  TH1D *HstTech_ratio2  = HstRatio("HstTech_ratio2",  HTot2_vTcPR_Ceex2n, vcum_ISR2_FSR2, kBlack);
  TH1D *HstTech_ratio4  = HstRatio("HstTech_ratio4",  HTot2_vTcPR_EEX2,   vcum_ISR2_FSR2, kMagenta);

  TH1D *HstTech_ratio = HstTech_ratio0;   // FoamEEX2/KKsem IFIoff magenta

  HstTech_ratio->SetStats(0);
  HstTech_ratio->SetTitle(0);
  HstTech_ratio->SetMinimum(1 -0.0007); HstTech_ratio->SetMaximum(1 +0.0007);  // zoom
  HstTech_ratio->GetXaxis()->SetTitle("v_{max}");
  HstTech_ratio->DrawCopy("h");

  HstTech_ratio1->DrawCopy("hsame");      // FoamGPS/KKsem IFIoff   green
  //HstTech_ratio3->DrawCopy("hsame");    // testing norm. Foam/KKsem
  HstTech_ratio4->DrawCopy("hsame");      // KKMCeex/KKsem IFIoff   magenta
  HstTech_ratio2->DrawCopy("hsame");      // KKMCceexn/KKsem IFIoff black

  double ycapt = 0.40; // starting value, to be decremented below
  PlotSame2(HstTech_ratio4, ycapt, kBlue,    0.12, "(a)", "KKMC_EEX.2/KKsem.2   IFIoff ");
  PlotSame2(HstTech_ratio2, ycapt, kBlack,   0.16, "(b)", "KKMC_CEEX.2/KKsem.2  IFIoff");
  PlotSame2(HstTech_ratio0, ycapt, kMagenta, 0.04, "(c)", "Foam3.2_EEX/KKsem.2  IFIoff ");
  PlotSame2(HstTech_ratio1, ycapt, kGreen,   0.08, "(d)", "Foam3.2_GPS/KKsem.2  IFIoff ");

  TH1D *hOne = (TH1D*)HstTech_ratio->Clone("hOne");  // unity line
  for(int i=1; i <= hOne->GetNbinsX() ; i++) { hOne->SetBinContent(i, 1); hOne->SetBinError(i, 0);}
  hOne->SetLineColor(kBlack);
  hOne->DrawCopy("hsame");

  CaptT->DrawLatex(0.00,0.96,"#sigma^{IFIoff}(v_{max}) ");

  //====================plot2========================
  cFigTech->cd(2);
  TH1D *HstTech_diff0  = HstDiff("HstTech_diff0",   Hafb2_xmax_EEX2,    afbv_ISR2_FSR2, kMagenta);
  TH1D *HstTech_diff1  = HstDiff("HstTech_diff1",   Hafb2_xmax_Ceex2n,  afbv_ISR2_FSR2, kGreen);
  TH1D *HstTech_diff2  = HstDiff("HstTech_diff2",   HAfb2_vTcPR_Ceex2n, afbv_ISR2_FSR2, kBlack);
  TH1D *HstTech_diff4  = HstDiff("HstTech_diff4",   HAfb2_vTcPR_EEX2,   afbv_ISR2_FSR2, kBlue);

  TH1D *HstTech_diff= HstTech_diff4;
  HstTech_diff->SetStats(0); HstTech_diff->SetTitle(0);
  HstTech_diff->SetMinimum(-0.0004); HstTech_diff->SetMaximum( 0.0004);  // zoom

  HstTech_diff->GetXaxis()->SetTitle("v_{max}");
  HstTech_diff->DrawCopy("h");

  ycapt = 0.90; // starting value, to be decremented below
  PlotSame2(HstTech_diff4,   ycapt, kBlue,     0.120, "(a)", "KKMC_EEX.2  -KKsem.2 IFIoff ");
  PlotSame2(HstTech_diff2,   ycapt, kBlack,    0.140, "(b)", "KKMC_CEEX.2 -KKsem.2 IFIoff ");
  PlotSame2(HstTech_diff0,   ycapt, kMagenta,  0.160, "(c)", "Foam3_EEX.2 -KKsem.2 IFIoff ");
  PlotSame2(HstTech_diff1,   ycapt, kGreen,    0.180, "(d)", "Foam3.2_GPS -KKsem.2 IFIoff ");

  TH1D *hZero = (TH1D*)HstTech_diff->Clone("hZero");  // unity line
  for(int i=1; i <= hZero->GetNbinsX() ; i++) { hZero->SetBinContent(i, 0); hZero->SetBinError(i, 0);}
  hZero->SetLineColor(kRed);
  hZero->DrawCopy("hsame");

  CaptT->DrawLatex(0.00,0.96,"A_{FB}^{IFIoff}(v_{max})");
  ycapt =0.30;
  CaptT->DrawLatex(0.40,ycapt,gTextEne);  ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextNev);  ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextNev2); ycapt += -0.04;

  cFigTech->cd();
  //================================================

}//FigTech


///////////////////////////////////////////////////////////////////////////////////
void FigInfo()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigInfo =========================== "<<endl;

  TH1D *hst_weight3  = (TH1D*)DiskFileF.Get("hst_weight3"); // Foam3
  TH1D *hst_weight5  = (TH1D*)DiskFileF.Get("hst_weight5"); // Foam5
 //
  //*****************************************************************************
  ///////////////////////////////////////////////////////////////////////////////
  TCanvas *cFigInfo = new TCanvas("cFigInfo","FigInfo ", gXcanv, gXcanv,   1200, 550);
  //                            Name    Title                     xoff,yoff, WidPix,HeiPix
  gXcanv += 100; gYcanv += 50;
  cFigInfo->SetFillColor(10);
  ////////////////////////////////////////////////////////////////////////////////
  cFigInfo->Divide( 2,  0);
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  //====================plot1========================
  //                sigma(vmax)
  cFigInfo->cd(1);
  //gPad->SetLogy(); // !!!!!!
  // MC v-true direct
  TH1D *Hst1 = hst_weight3;  //  weight of Foam
  //
  Hst1->SetStats(0);
  Hst1->SetTitle(0);
  Hst1->DrawCopy("h");

  CaptT->DrawLatex(0.02,0.95, " Weight distribution Foam");
  //====================plot2========================
  cFigInfo->cd(2);

  Hst1 = hst_weight5;
  Hst1->SetTitle(0);
  Hst1->DrawCopy("h");

  cFigInfo->cd();
  //================================================
}//FigInfo


///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  /////////////////////////////////////////////////////////
  LibSem.Initialize(DiskFileA);  // for non-farm case
  /////////////////////////////////////////////////////////
  // Reading directly KKMC input (farming)
  /*
  double xpar[10001];
  int jmax = 10000;
  LibSem.ReaData("../../.KK2f_defaults",     jmax, xpar);  // numbering as in input!!!
  char dname[100];  sprintf(dname,XparFile);
  LibSem.ReaData(dname, -jmax, xpar);  // jmax<0 for append mode
  double CMSene  = xpar[ 1];
  cout<< "////// Main: CMSene = "<<  CMSene  <<endl;
  LibSem.Initialize(xpar);  // for non-farm case
  */
  int Nodes, Nodes2;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");
  Nodes    = HST_KKMC_NORMA->GetBinContent(511);       // No of farm nodes (trick)
  gCMSene  = HST_KKMC_NORMA->GetBinContent(1)/Nodes;   // CMSene=xpar(1), farn adjusted
  gNevTot  = HST_KKMC_NORMA->GetEntries();             // MC statistics from KKMC
  sprintf(gTextEne,"#sqrt{s} =%4.2fGeV", gCMSene);
  sprintf(gTextNev,"KKMC:%10.2e events", gNevTot);
  //
  TH1D *HST_FOAM_NORMA3 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA3");
  Nodes2   =  HST_FOAM_NORMA3->GetBinContent(511);    // No of farm nodes (trick)
  double  CMSeneF  = HST_FOAM_NORMA3->GetBinContent(1)/Nodes2; // CMSene=xpar(1)
  if( fabs(gCMSene/CMSeneF-1) >1e-4 ){
	  cout<<" +++++ Wrong input files !!!! KKMC "<< gCMSene <<"GeV and  FOAM "<< CMSeneF<<"GeV"<<endl;
	  exit(19);
  }
  gNevTot2  = HST_FOAM_NORMA3->GetEntries();       // MC statistics from KKMC
  sprintf(gTextNev2,"FOAM:%10.2e events", gNevTot2);

//////////////////////////////////////////////////////////////////////////
// ========= Preparing plots ==========
  DiskFileB.cd();
  TestNorm();          // special test of normalizaion
  HisReMakeKKMC(  &DiskFileA);   // reprocessing MC histos from KKC
  HisReMakeFoam35(&DiskFileF);   // reprocessing MC histos from Foam

  KKsemMakeHisto();    // prepare histos from KKsem
//========== PLOTTING ==========
// vmax=1, sigma(v) and sigma(vmax) KKMC/Foam
//  FigVdist(); // OBSOLETE
  FigInfo();
  FigAfb2();  // vmax =0.2
  FigTech();  // vmax =0.2
  FigTech0();  // vmax =0.2
// weight distribution
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //cout<<"------------------------------A.ls----------------------------------"<<endl;
  //DiskFileA.ls();
  //cout<<"------------------------------F.ls----------------------------------"<<endl;
  //DiskFileF.ls();
  cout<<"------------------------A.GetListOfKeys-----------------------------"<<endl;
  DiskFileA.GetListOfKeys()->Print();
  cout<<"------------------------F.GetListOfKeys-----------------------------"<<endl;
  DiskFileF.GetListOfKeys()->Print();
  //
  cout<< "CMSene[GeV] = "<< gCMSene<< endl;
  cout<< "KKMC: No. of farm nodes="<< Nodes  << "  Tot no. of events = "<<gNevTot<< endl;
  cout<< "FOAM: No. of farm nodes="<< Nodes2 << "  Tot no. of events = "<<gNevTot2<<endl;
  //cout<<"------------------------------end---------------------------------"<<endl;
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}


