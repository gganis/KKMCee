//////////////////////////////////////////////////////////////////////
//    make PlotFoam2-run
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

#include "HisNorm.h"
#include "KKplot.h"

//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT
//=============================================================================
////  *** KKMC
//TFile DiskFileA("../workKKMC/histo.root");
// August2017 runs
//TFile DiskFileA("../workKKMC/histo.root_10GeV_5.7G"); //
//TFile DiskFileA("../workKKMC/histo.root_88GeV_2.1G"); //
TFile DiskFileA("../workKKMC/histo.root_95GeV_16G");
//TFile DiskFileA("../workKKMC/histo.root_91GeV_9G"); ///????

// July2017 runs obsolete
//TFile DiskFileA("../workKKMC/histo.root"); // current
//TFile DiskFileA("../workKKMC/histo.root_91GeV_6G"); //
//TFile DiskFileA("../workKKMC/histo.root_88GeV_4G"); //
//TFile DiskFileA("../workKKMC/histo.root_10GeV_5.7G"); //
//TFile DiskFileA("../workKKMC/histo.root_95GeV.4G");   //

////  *** FOAM 5dim
//TFile DiskFileF("../workFOAM/histo.root"); // current
//TFile DiskFileF("../workFOAM/histo.root_10GeV_37G_vmax0.2");
//TFile DiskFileF("../workFOAM/histo.root_88GeV_16G");
//TFile DiskFileF("../workFOAM/histo.root_91GeV_45G");
TFile DiskFileF("../workFOAM/histo.root_95GeV_10G");
//TFile DiskFileF("../workFOAM/histo.root_10GeV_32G");


TFile DiskFileB("RhoSemi.root","RECREATE","histograms");

///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
double gCMSene, gNevTot, gNevTot2; // from KKMC and KKfoam MC runs (histograms)
char   gTextEne[100], gTextNev[100], gTextNev2[100];

//
int    gNbMax=50;         // gCosTheta = 45/50=0.90
double gCosTheta=1.00;    // to be synchronized with gNbMax
//
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
  //double yy=ycapt;
  CaptT->DrawLatex(0.40,ycapt, opis);
}// PlotSame




///////////////////////////////////////////////////////////////////////////////////
void ReMakeMChisto(){
	//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeMChisto  BEGIN  ============================"<<endl;
//////////////////////////////////////////////////////////////////
  cout<<"  Renormalizing  and reprocessing histograms from FOAM"<<endl;

  TH1D *HST_FOAM_NORMA3 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA3");
  TH1D *HST_FOAM_NORMA5 = (TH1D*)DiskFileF.Get("HST_FOAM_NORMA5");

  TH2D *SCT_xc_Ceex2n = (TH2D*)DiskFileF.Get("SCT_xc_Ceex2n");  // FOAM small range x<0.20
  TH2D *SCT_xc_EEX2   = (TH2D*)DiskFileF.Get("SCT_xc_EEX2");    // FOAM small range x<0.20
  TH2D *SCT_xc_Ceex2  = (TH2D*)DiskFileF.Get("SCT_xc_Ceex2");   // FOAM small range x<0.20

  // KKFoam1: sigma(vmax) and AFB(vmax) from scat. vmax<0.2, 100 bins in ctheta
  // gNbMax=45;           // cosThetaMax = 45/50=0.90 Now global variable
  TH1D                 *Htot2_xmax_Ceex2n, *Hafb2_xmax_Ceex2n;
  ProjV( SCT_xc_Ceex2n, Htot2_xmax_Ceex2n,  Hafb2_xmax_Ceex2n, gNbMax);  //!!!!
  Htot2_xmax_Ceex2n->SetName("Htot2_xmax_Ceex2n");
  Hafb2_xmax_Ceex2n->SetName("Hafb2_xmax_Ceex2n");
  //
  TH1D               *Htot2_xmax_EEX2, *Hafb2_xmax_EEX2;
  ProjV( SCT_xc_EEX2, Htot2_xmax_EEX2,  Hafb2_xmax_EEX2, gNbMax);  //!!!!
  Htot2_xmax_EEX2->SetName("Htot2_xmax_EEX2");
  Hafb2_xmax_EEX2->SetName("Hafb2_xmax_EEX2");
  //
  TH1D                *Htot2_xmax_Ceex2, *Hafb2_xmax_Ceex2;
  ProjV( SCT_xc_Ceex2, Htot2_xmax_Ceex2,  Hafb2_xmax_Ceex2, gNbMax);  //!!!!
  Htot2_xmax_Ceex2->SetName("Htot2_xmax_Ceex2");
  Hafb2_xmax_Ceex2->SetName("Hafb2_xmax_Ceex2");

  cout<<"================ ReMakeMChisto ENDs  ============================="<<endl;
  cout<<"==================================================================="<<endl;
}//RemakeMChisto





///////////////////////////////////////////////////////////////////////////////////
void ReMakeMChisto2(){
	//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeMChisto2  BEGIN ============================"<<endl;

  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA.Get("HST_KKMC_NORMA");

  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPR_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPR_Ceex2n") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPR_EEX2") );

  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2") );

////////////////////////////////////////////////////////////////////
// Pure KKMC reprocessing part
// from bigger scattergram and restricted vmax<0.2
  //////////////////////////////////////////////////////////////////
   // Wide range, vmax<1.
    TH2D *sct_vTcPR_Ceex2  = (TH2D*)DiskFileA.Get("sct_vTcPR_Ceex2");
    TH2D *sct_vTcPR_Ceex2n = (TH2D*)DiskFileA.Get("sct_vTcPR_Ceex2n");
    TH2D *sct_vTcPR_EEX2   = (TH2D*)DiskFileA.Get("sct_vTcPR_EEX2");
    TH2D *sct_vTcPL_Ceex2  = (TH2D*)DiskFileA.Get("sct_vTcPL_Ceex2");
    cout<<"ReMakeMChisto2 [2]"<<endl;
    ///****************************************************************************************
    ///****************************************************************************************
    /// Distributions of v=vTrue<vmax<0.20, c=cos(theta) with 100 bins
    //gNbMax=45;         // cosThetaMax = 45/50=0.90 Now global variable
    // IFI on
    TH1D                    *HTot2_vTcPR_Ceex2, *HAfb2_vTcPR_Ceex2;
    ProjV( sct_vTcPR_Ceex2,  HTot2_vTcPR_Ceex2,  HAfb2_vTcPR_Ceex2, gNbMax);  //!!!!
    HTot2_vTcPR_Ceex2->SetName("HTot2_vTcPR_Ceex2");
    HAfb2_vTcPR_Ceex2->SetName("HAfb2_vTcPR_Ceex2");
    // IFI off
    TH1D                    *HTot2_vTcPR_Ceex2n, *HAfb2_vTcPR_Ceex2n;
    ProjV( sct_vTcPR_Ceex2n, HTot2_vTcPR_Ceex2n,  HAfb2_vTcPR_Ceex2n, gNbMax);  //!!!!
    HTot2_vTcPR_Ceex2n->SetName("HTot2_vTcPR_Ceex2n");
    HAfb2_vTcPR_Ceex2n->SetName("HAfb2_vTcPR_Ceex2n");
    // IFI off
    TH1D                    *HTot2_vTcPR_EEX2, *HAfb2_vTcPR_EEX2;
    ProjV( sct_vTcPR_EEX2, HTot2_vTcPR_EEX2,  HAfb2_vTcPR_EEX2, gNbMax);  //!!!!
    HTot2_vTcPR_EEX2->SetName("HTot2_vTcPR_EEX2");
    HAfb2_vTcPR_EEX2->SetName("HAfb2_vTcPR_EEX2");
    // IFI on
    TH1D                    *HTot2_vTcPL_Ceex2, *HAfb2_vTcPL_Ceex2;
    ProjV( sct_vTcPL_Ceex2,  HTot2_vTcPL_Ceex2,  HAfb2_vTcPL_Ceex2, gNbMax);  //!!!!
    HTot2_vTcPL_Ceex2->SetName("HTot2_vTcPL_Ceex2");
    HAfb2_vTcPL_Ceex2->SetName("HAfb2_vTcPL_Ceex2");


    ///****************************************************************************************

  cout<<"================ ReMakeMChisto2 ENDs  ============================="<<endl;
  cout<<"==================================================================="<<endl;
}//RemakeMChisto2





///////////////////////////////////////////////////////////////////////////////////
void FigAfb3()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigAfb3 =========================== "<<endl;

  TH1D *HAfb2_vTcPR_Ceex2n = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_Ceex2n"); // KKMC[PR]
  TH1D *HAfb2_vTcPR_Ceex2  = (TH1D*)DiskFileB.Get("HAfb2_vTcPR_Ceex2");  // KKMC[PR]
  //
  TH1D *HAfb2_vTcPL_Ceex2  = (TH1D*)DiskFileB.Get("HAfb2_vTcPL_Ceex2");  // KKMC[PL]

  TH1D *Hafb2_xmax_Ceex2n  = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex2n");  // KKFoam scat.
  TH1D *Hafb2_xmax_Ceex2   = (TH1D*)DiskFileB.Get("Hafb2_xmax_Ceex2");   // KKFoam scat.

  TH1D *HST_PLBhard =(TH1D*)HAfb2_vTcPL_Ceex2->Clone("HST_PLBhard");
  LibSem.Ord1fill(HST_PLBhard,105);                           // IFI only, hard part, PLB219
  HST_PLBhard->SetLineColor(kCyan);

  // This HST_PL is temporrary insert
  double alfinv  = 137.035989;
  double alfpi   = 1/alfinv/3.1415926535897932;
  TH1D *HST_PL =(TH1D*)HAfb2_vTcPL_Ceex2->Clone("HST_PL");
  HST_PL->SetLineColor(kYellow);
  int Nbin    = HST_PL->GetNbinsX();
  double vmax = HST_PL->GetXaxis()->GetXmax();
  for(int i=1; i <= Nbin ; i++) {
	  double vv = (i*vmax)/Nbin;
	  // A_FB from PLB219,p103, pure gamma exch.
	  double afb = 3e0/2e0* alfpi*( 3*vv+log(1-vv/2) ); // only gamma
      double KFi=11; int KFf=13;
	  HST_PL->SetBinContent(i, afb);
	  HST_PL->SetBinError(i, 0);
  }// i

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
//
  //*****************************************************************************
  TCanvas *cFigAfb3a = new TCanvas("cFigAfb3a","FigAfb3a", 70, 70,   600, 600);
  //                                 Name    Title      xoff,yoff, WidPix,HeiPix
  cFigAfb3a->SetFillColor(10);
  cFigAfb3a->cd();

  TH1D *Hst21_diff   = HstDiff("Hst21_diff",    HAfb2_vTcPR_Ceex2,  HAfb2_vTcPR_Ceex2n, kBlack);
  TH1D *HST21_diff   = HstDiff("HST21_diff",    Hafb2_xmax_Ceex2,   Hafb2_xmax_Ceex2n,  kMagenta);
  TH1D *HstKF_diff   = HstDiff("HstKF_diff",    HAfb2_vTcPR_Ceex2,  Hafb2_xmax_Ceex2,   kGreen);
  TH1D *HstKFn_diff  = HstDiff("HstKFn_diff",   HAfb2_vTcPR_Ceex2n, Hafb2_xmax_Ceex2n,  kBlue);
  TH1D *HstPL_diff   = HstDiff("HstPL_diff",    HAfb2_vTcPL_Ceex2,  Hafb2_xmax_Ceex2,   kRed);

  Hst21_diff->SetStats(0);
  Hst21_diff->SetTitle(0);
  if( fabs(gCMSene -95.0) < 1.0) { Hst21_diff->SetMinimum(-0.004);  Hst21_diff->SetMaximum( 0.006);}  // 95GeV
  if( fabs(gCMSene -88.0) < 1.0) { Hst21_diff->SetMinimum(-0.002);  Hst21_diff->SetMaximum( 0.008);}  // 88GeV
  if( fabs(gCMSene -91.0) < 1.0) { Hst21_diff->SetMinimum(-0.0005); Hst21_diff->SetMaximum( 0.0035);} // 91GeV
  if( fabs(gCMSene -10.0) < 1.0) { Hst21_diff->SetMinimum(-0.002);  Hst21_diff->SetMaximum( 0.006);}  // 10GeV
  Hst21_diff->GetXaxis()->SetTitle("v_{max}");
  Hst21_diff->DrawCopy("h");

  double ycapt = 0.90; // starting value, to be decremented below
  CaptT->SetTextColor(kBlack); ycapt += -0.04;
  CaptT->DrawLatex(0.40,ycapt,gTextEne);
  CaptT->DrawLatex(0.06,0.95, "A_{FB}(v_{max}) ");

  PlotSame(Hst21_diff,      ycapt, kBlack,    "KKMC[vPR]: IFIon-IFIoff ");
  PlotSame(HST21_diff,      ycapt, kMagenta,  "KKFoam:    IFIon-IFIoff ");
  PlotSame(HstKFn_diff,     ycapt, kBlue,     "KKMC[vPR]-KKFoam:  IFIoff");
  PlotSame(HstKF_diff,      ycapt, kGreen,    "KKMC[vPR]-KKFoam:  IFIon");
  PlotSame(HstPL_diff,      ycapt, kRed,      "KKMC[PL] -KKFoam:  IFIon");
  PlotSame(HST_PLBhard,     ycapt, kCyan,     "PLB219, IFI only,  hard");

  // extras !!!??? ???
  HST_PL->DrawCopy("hsame");   // !!!??? ??? dubious agreement

// zero line
  TH1D *hZero = (TH1D*)HAfb2_vTcPR_Ceex2n->Clone("hZero");  // zero line
  for(int i=1; i <= hZero->GetNbinsX() ; i++) { hZero->SetBinContent(i, 0); hZero->SetBinError(i, 0);}

  hZero->DrawCopy("hsame");

  //*****************************************************************************
  TCanvas *cFigAfb3b = new TCanvas("cFigAfb3b","FigAfb3b", 170, 100,   600, 600);
  //                                 Name    Title      xoff,yoff, WidPix,HeiPix
  cFigAfb3b->SetFillColor(10);
  cFigAfb3b->cd();
  //=======================================================
  TH1D *HstPL2_diff = HstDiff("HstPL2_diff",  HAfb2_vTcPR_Ceex2,  HAfb2_vTcPL_Ceex2, kBlack);

  HstPL2_diff->SetMinimum(-0.0004);  HstPL2_diff->SetMaximum( 0.0004);  // zoom
  HstPL2_diff->SetStats(0);
  HstPL2_diff->SetTitle(0);
  HstPL2_diff->GetXaxis()->SetTitle("v_{max}");
  HstPL2_diff->DrawCopy("h");

  hZero->DrawCopy("hsame");
  CaptT->DrawLatex(0.22,0.95,"KKMC: A_{FB}[#theta_{PRD}]-A_{FB}[#theta_{PL}] ");
//  CaptT->DrawLatex(0.22,0.95,"A_{FB}^{#bullet}-A*_{FB} ");
  CaptT->DrawLatex(0.50,0.75,gTextEne);

  cFigAfb3b->cd();
  //================================================
}//FigAfb3




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
  //
  cout<< "CMSene[GeV] = "<< gCMSene<< endl;
  cout<< "KKMC: No. of farm nodes="<< Nodes  << "  Tot no. of events = "<<gNevTot<< endl;
  cout<< "FOAM: No. of farm nodes="<< Nodes2 << "  Tot no. of events = "<<gNevTot2<<endl;
//////////////////////////////////////////////////////////////////////////
// ========= Preparing plots ==========
  DiskFileB.cd();
  ReMakeMChisto();     // reprocessing MC histos from KKC and Foam
  ReMakeMChisto2();    // reprocessing MC histos from KKC and Foam
//========== PLOTTING ==========
  FigAfb3();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //cout<<"------------------------------A.ls----------------------------------"<<endl;
  //DiskFileA.ls();
  //cout<<"------------------------------F.ls----------------------------------"<<endl;
  //DiskFileF.ls();
  //cout<<"------------------------A.GetListOfKeys-----------------------------"<<endl;
  //DiskFileA.GetListOfKeys()->Print();
  cout<<"------------------------F.GetListOfKeys-----------------------------"<<endl;
  DiskFileF.GetListOfKeys()->Print();
  //cout<<"------------------------------end---------------------------------"<<endl;
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}


