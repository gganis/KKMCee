//    make PlotAFB2-run
//    Plot for FCC week 2016 in Rome and "Beyond precision..." paper

// This is only AFB study for several energies,
// renomalizing scattergrams is no necesssary!

#include <iostream>
#include <fstream>
#include <iomanip>
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

#include "KKplot.h"
#include "HisNorm.h"

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
//=============================================================================
//  ROOT  ROOT ROOT   ROOT  ROOT  ROOT  ROOT  ROOT  ROOT  ROOT   ROOT   ROOT 
//=============================================================================
// Latest from /workKKMC
//TFile DiskFileA88("../workKKMC/histo.root_88GeV.new");  // jan.2018
//TFile DiskFileA95("../workKKMC/histo.root_95GeV.new");  // jan.2018

//
TFile DiskFileA88("../workKKMC/histo.root_88GeV_10G");  // jan.2018
//
TFile DiskFileA95("../workKKMC/histo.root_95GeV_26G");   // oct.2017
////TFile DiskFileA88("../workKKMC/histo.root_88GeV_2.5G");  // oct.2017 OBSOLETE
TFile DiskFileA91("../workKKMC/histo.root_91GeV_3.5G");  // oct.2017
TFile DiskFileA10("../workKKMC/histo.root_10GeV_5.8G");  // oct.2017
//
//TFile DiskFileA88("../workKKMC/histo.root"); // actual
//
// Latest from /workFOAM
TFile DiskFileF95("../workFOAM/histo.root_95GeV_23G");    // Dec.  2017 runs
//TFile DiskFileF95("../workFOAM/histo.root_95GeV_57G");  // Sept. 2017 runs
//
TFile DiskFileF88("../workFOAM/histo.root_88GeV_22G");    // Dec.  2017 runs
//TFile DiskFileF88("../workFOAM/histo.root_88GeV_15G");  // Sept. 2017 runs


////////////////////////////////////////////////////////////////
// Archive from /workAFB, as for Rome
//TFile DiskFileA95("../workAFB/rmain.root_95GeV_100M");
//TFile DiskFileA88("../workAFB/rmain.root_88GeV_100M");
//TFile DiskFileA91("../workAFB/rmain.root_91GeV_48M");
//TFile DiskFileA10("../workAFB/rmain.root_10GeV_30M");
/// Current
//TFile DiskFileA95("../workAFB/rmain.root");
//TFile DiskFileA88("../test0/rmain.root");
TFile DiskFileB("RhoAFB.root","RECREATE","histograms");
///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
int    kGold=kOrange-3, kBrune=46, kPine=kGreen+3;
//
float  gXcanv = 50, gYcanv = 50;
//
//int    gTogEne = 1;   // 10 GeV and MZ included
int    gTogEne = 0; // 10 GeV and MZ exluded

//Double_t sqr( const Double_t x ){ return x*x;};


///////////////////////////////////////////////////////////////////////////////////
void PlotSame2(TH1D *HST, double &ycapt, Int_t kolor, double xx,  TString label,  TString opis)
{
  TLatex *CaptT = new TLatex();
  CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.035);
  HST->SetLineColor(kolor);
  HST->DrawCopy("hsame");      // Magenta
  CaptT->SetTextColor(kolor);
//  ycapt += -0.050;
  ycapt += -0.045;
  double xcapt = 0.50;
  CaptT->DrawLatex(xcapt,ycapt, opis);
  CaptT->DrawLatex(xcapt-0.07,ycapt, label);
  //
  TLatex *CaptS = new TLatex();
  CaptS->SetTextSize(0.040);
  CaptS->SetTextAlign(21);
  CaptS->SetTextColor(kolor);
  int ib = HST->FindBin(xx);
  double yy= HST->GetBinContent(ib);
  CaptS->DrawLatex(xx,yy,label);
}// PlotSame2




///////////////////////////////////////////////////////////////////////////////////
void HistNormalize(){
  //
  cout<<"----------------------------- HistNormalize ------------------------------------"<<endl;
  DiskFileA95.ls("");
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA95.Get("HST_KKMC_NORMA");
  //
  // renomalizing scattergrams is not necesssary!
  //  BIG scatergrams
  //HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA95.Get("sct_vAcPR_Ceex2") );
  //HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA95.Get("sct_vAcPR_Ceex2n") );
  //HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA95.Get("sct_vTcPL_Ceex2") );
  //HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA95.Get("sct_vKcPL_Ceex2") );
}

///////////////////////////////////////////////////////////////////////////////////
void ReMakeMChisto(){
	// Here we produce semianalytical plots using KKsem program, No plotting
	// also some MC histos are preprocessed
	//------------------------------------------------------------------------
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeMChisto  BEGIN  ============================"<<endl;
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA95.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR

  //****************************************************************************************
  // Pure MC reprocessing part
  //
  TH2D *sct9_vAcPR_Ceex2  = (TH2D*)DiskFileA95.Get("sct_vAcPR_Ceex2");
  TH2D *sct9_vAcPR_Ceex2n = (TH2D*)DiskFileA95.Get("sct_vAcPR_Ceex2n");
  //
  TH2D *sct8_vAcPR_Ceex2  = (TH2D*)DiskFileA88.Get("sct_vAcPR_Ceex2");
  TH2D *sct8_vAcPR_Ceex2n = (TH2D*)DiskFileA88.Get("sct_vAcPR_Ceex2n");
  //
  TH2D *sctZ_vAcPR_Ceex2  = (TH2D*)DiskFileA91.Get("sct_vAcPR_Ceex2");
  TH2D *sctZ_vAcPR_Ceex2n = (TH2D*)DiskFileA91.Get("sct_vAcPR_Ceex2n");
  ///
  TH2D *sct1_vAcPR_Ceex2  = (TH2D*)DiskFileA10.Get("sct_vAcPR_Ceex2");
  TH2D *sct1_vAcPR_Ceex2n = (TH2D*)DiskFileA10.Get("sct_vAcPR_Ceex2n");
  //
  /*
  //[[[[[[[[[[[[[[[[[[[
  TH2D *sct9_vTcPL_Ceex0  = (TH2D*)DiskFileA95.Get("sct_vTcPL_Ceex0");
  TH2D *sct9_vTcPL_Ceex0n = (TH2D*)DiskFileA95.Get("sct_vTcPL_Ceex0n");
  TH2D *sct9_vTcPL_Ceex1  = (TH2D*)DiskFileA95.Get("sct_vTcPL_Ceex1");
  TH2D *sct9_vTcPL_Ceex1n = (TH2D*)DiskFileA95.Get("sct_vTcPL_Ceex1n");
  //]]]]]]]]]]]]]]]]]]]
  */
  TH2D *sct9_vTcPL_Ceex2  = (TH2D*)DiskFileA95.Get("sct_vTcPL_Ceex2");
  TH2D *sct9_vTcPL_Ceex2n = (TH2D*)DiskFileA95.Get("sct_vTcPL_Ceex2n");
  /*
  //[[[[[[[[[[[[[[[[[[[
  TH2D *sct8_vTcPL_Ceex0  = (TH2D*)DiskFileA88.Get("sct_vTcPL_Ceex0");
  TH2D *sct8_vTcPL_Ceex0n = (TH2D*)DiskFileA88.Get("sct_vTcPL_Ceex0n");
  TH2D *sct8_vTcPL_Ceex1  = (TH2D*)DiskFileA88.Get("sct_vTcPL_Ceex1");
  TH2D *sct8_vTcPL_Ceex1n = (TH2D*)DiskFileA88.Get("sct_vTcPL_Ceex1n");
  //]]]]]]]]]]]]]]]]]]]
  */
  TH2D *sct8_vTcPL_Ceex2  = (TH2D*)DiskFileA88.Get("sct_vTcPL_Ceex2");
  TH2D *sct8_vTcPL_Ceex2n = (TH2D*)DiskFileA88.Get("sct_vTcPL_Ceex2n");
  //
  TH2D *sctZ_vTcPL_Ceex2  = (TH2D*)DiskFileA91.Get("sct_vTcPL_Ceex2");
  TH2D *sctZ_vTcPL_Ceex2n = (TH2D*)DiskFileA91.Get("sct_vTcPL_Ceex2n");
  ///
  TH2D *sct1_vTcPL_Ceex2  = (TH2D*)DiskFileA10.Get("sct_vTcPL_Ceex2");
  TH2D *sct1_vTcPL_Ceex2n = (TH2D*)DiskFileA10.Get("sct_vTcPL_Ceex2n");

  // ****************************************************************************************
  /// Distributions of v with limited c=cos(theta)
  //  without cutoff on c=cos(thetaPRD)
  int nbMax=0;   // cosThetaMax = 1.0, no cut
  //nbMax=50;      // cosThetaMax = 50/50=1.00
  nbMax=45;      // cosThetaMax = 45/50=0.90
  // ---------------------- 95GeV ----------------------------------
  TH1D  *Hsig9_vAcPR_Ceex2  = HstProjV("Hsig9_vAcPR_Ceex2", sct9_vAcPR_Ceex2, nbMax);
  TH1D  *Hafb9_vAcPR_Ceex2  = HstProjA("Hafb9_vAcPR_Ceex2", sct9_vAcPR_Ceex2, nbMax);
  //******** IFI off
  TH1D  *Hsig9_vAcPR_Ceex2n  = HstProjV("Hsig9_vAcPR_Ceex2n", sct9_vAcPR_Ceex2n, nbMax);
  TH1D  *Hafb9_vAcPR_Ceex2n  = HstProjA("Hafb9_vAcPR_Ceex2n", sct9_vAcPR_Ceex2n, nbMax);
  //
  // ---------------------- 88GeV ----------------------------------
  TH1D  *Hsig8_vAcPR_Ceex2  = HstProjV("Hsig8_vAcPR_Ceex2", sct8_vAcPR_Ceex2, nbMax);
  TH1D  *Hafb8_vAcPR_Ceex2  = HstProjA("Hafb8_vAcPR_Ceex2", sct8_vAcPR_Ceex2, nbMax);
  Hafb8_vAcPR_Ceex2->Scale(-1.0);
  //******** IFI off
  TH1D  *Hsig8_vAcPR_Ceex2n  = HstProjV("Hsig8_vAcPR_Ceex2n", sct8_vAcPR_Ceex2n, nbMax);
  TH1D  *Hafb8_vAcPR_Ceex2n  = HstProjA("Hafb8_vAcPR_Ceex2n", sct8_vAcPR_Ceex2n, nbMax);
  Hafb8_vAcPR_Ceex2n->Scale(-1.0);
  //
  // ---------------------- 91.2GeV ----------------------------------
  TH1D  *HsigZ_vAcPR_Ceex2  = HstProjV("HsigZ_vAcPR_Ceex2", sctZ_vAcPR_Ceex2, nbMax);
  TH1D  *HafbZ_vAcPR_Ceex2  = HstProjA("HafbZ_vAcPR_Ceex2", sctZ_vAcPR_Ceex2, nbMax);
  //******** IFI off
  TH1D  *HsigZ_vAcPR_Ceex2n  = HstProjV("HsigZ_vAcPR_Ceex2n", sctZ_vAcPR_Ceex2n, nbMax);
  TH1D  *HafbZ_vAcPR_Ceex2n  = HstProjA("HafbZ_vAcPR_Ceex2n", sctZ_vAcPR_Ceex2n, nbMax);
  // ---------------------- 10GeV ----------------------------------
  TH1D  *Hsig1_vAcPR_Ceex2  = HstProjV("Hsig1_vAcPR_Ceex2", sct1_vAcPR_Ceex2, nbMax);
  TH1D  *Hafb1_vAcPR_Ceex2  = HstProjA("Hafb1_vAcPR_Ceex2", sct1_vAcPR_Ceex2, nbMax);
  //******** IFI off
  TH1D  *Hsig1_vAcPR_Ceex2n  = HstProjV("Hsig1_vAcPR_Ceex2n", sct1_vAcPR_Ceex2n, nbMax);
  TH1D  *Hafb1_vAcPR_Ceex2n  = HstProjA("Hafb1_vAcPR_Ceex2n", sct1_vAcPR_Ceex2n, nbMax);

  // Warning! scattergrams for vTcPL noIFI were missing in Rome version
  // ---------------------- 95GeV ----------------------------------
/*
  //[[[[[[[[
  //******** IFI on
  TH1D  *Hsig9_vTcPL_Ceex0  = HstProjV("Hsig9_vTcPL_Ceex0", sct9_vTcPL_Ceex0, nbMax);
  TH1D  *Hafb9_vTcPL_Ceex0  = HstProjA("Hafb9_vTcPL_Ceex0", sct9_vTcPL_Ceex0, nbMax);
  //******** IFI off
  TH1D  *Hsig9_vTcPL_Ceex0n  = HstProjV("Hsig9_vTcPL_Ceex0n", sct9_vTcPL_Ceex0n, nbMax);
  TH1D  *Hafb9_vTcPL_Ceex0n  = HstProjA("Hafb9_vTcPL_Ceex0n", sct9_vTcPL_Ceex0n, nbMax);
  //******** IFI on
  TH1D  *Hsig9_vTcPL_Ceex1  = HstProjV("Hsig9_vTcPL_Ceex1", sct9_vTcPL_Ceex1, nbMax);
  TH1D  *Hafb9_vTcPL_Ceex1  = HstProjA("Hafb9_vTcPL_Ceex1", sct9_vTcPL_Ceex1, nbMax);
  //******** IFI off
  TH1D  *Hsig9_vTcPL_Ceex1n  = HstProjV("Hsig9_vTcPL_Ceex1n", sct9_vTcPL_Ceex1n, nbMax);
  TH1D  *Hafb9_vTcPL_Ceex1n  = HstProjA("Hafb9_vTcPL_Ceex1n", sct9_vTcPL_Ceex1n, nbMax);
  //]]]]]]]]
 */
  //******** IFI on
  TH1D  *Hsig9_vTcPL_Ceex2  = HstProjV("Hsig9_vTcPL_Ceex2", sct9_vTcPL_Ceex2, nbMax);
  TH1D  *Hafb9_vTcPL_Ceex2  = HstProjA("Hafb9_vTcPL_Ceex2", sct9_vTcPL_Ceex2, nbMax);
  //******** IFI off
  TH1D  *Hsig9_vTcPL_Ceex2n  = HstProjV("Hsig9_vTcPL_Ceex2n", sct9_vTcPL_Ceex2n, nbMax);
  TH1D  *Hafb9_vTcPL_Ceex2n  = HstProjA("Hafb9_vTcPL_Ceex2n", sct9_vTcPL_Ceex2n, nbMax);
  //
  // ---------------------- 88GeV ----------------------------------
  //[[[[[[[
  /*
  //******** IFI on
  TH1D  *Hsig8_vTcPL_Ceex0  = HstProjV("Hsig8_vTcPL_Ceex0", sct8_vTcPL_Ceex0, nbMax);
  TH1D  *Hafb8_vTcPL_Ceex0  = HstProjA("Hafb8_vTcPL_Ceex0", sct8_vTcPL_Ceex0, nbMax);
  Hafb8_vTcPL_Ceex0->Scale(-1.0);
  //******** IFI off
  TH1D  *Hsig8_vTcPL_Ceex0n  = HstProjV("Hsig8_vTcPL_Ceex0n", sct8_vTcPL_Ceex0n, nbMax);
  TH1D  *Hafb8_vTcPL_Ceex0n  = HstProjA("Hafb8_vTcPL_Ceex0n", sct8_vTcPL_Ceex0n, nbMax);
  Hafb8_vTcPL_Ceex0n->Scale(-1.0);
  //******** IFI on
  TH1D  *Hsig8_vTcPL_Ceex1  = HstProjV("Hsig8_vTcPL_Ceex1", sct8_vTcPL_Ceex1, nbMax);
  TH1D  *Hafb8_vTcPL_Ceex1  = HstProjA("Hafb8_vTcPL_Ceex1", sct8_vTcPL_Ceex1, nbMax);
  Hafb8_vTcPL_Ceex1->Scale(-1.0);
  //******** IFI off
  TH1D  *Hsig8_vTcPL_Ceex1n  = HstProjV("Hsig8_vTcPL_Ceex1n", sct8_vTcPL_Ceex1n, nbMax);
  TH1D  *Hafb8_vTcPL_Ceex1n  = HstProjA("Hafb8_vTcPL_Ceex1n", sct8_vTcPL_Ceex1n, nbMax);
  Hafb8_vTcPL_Ceex1n->Scale(-1.0);
  //]]]]]]]
   */
  //******** IFI on
  TH1D  *Hsig8_vTcPL_Ceex2  = HstProjV("Hsig8_vTcPL_Ceex2", sct8_vTcPL_Ceex2, nbMax);
  TH1D  *Hafb8_vTcPL_Ceex2  = HstProjA("Hafb8_vTcPL_Ceex2", sct8_vTcPL_Ceex2, nbMax);
  Hafb8_vTcPL_Ceex2->Scale(-1.0);
  //******** IFI off
  TH1D  *Hsig8_vTcPL_Ceex2n  = HstProjV("Hsig8_vTcPL_Ceex2n", sct8_vTcPL_Ceex2n, nbMax);
  TH1D  *Hafb8_vTcPL_Ceex2n  = HstProjA("Hafb8_vTcPL_Ceex2n", sct8_vTcPL_Ceex2n, nbMax);
  Hafb8_vTcPL_Ceex2n->Scale(-1.0);
  //
  // ---------------------- 91.2GeV ----------------------------------
  TH1D  *HsigZ_vTcPL_Ceex2  = HstProjV("HsigZ_vTcPL_Ceex2", sctZ_vTcPL_Ceex2, nbMax);
  TH1D  *HafbZ_vTcPL_Ceex2  = HstProjA("HafbZ_vTcPL_Ceex2", sctZ_vTcPL_Ceex2, nbMax);
  //******** IFI off
  TH1D  *HsigZ_vTcPL_Ceex2n  = HstProjV("HsigZ_vTcPL_Ceex2n", sctZ_vTcPL_Ceex2n, nbMax);
  TH1D  *HafbZ_vTcPL_Ceex2n  = HstProjA("HafbZ_vTcPL_Ceex2n", sctZ_vTcPL_Ceex2n, nbMax);
  // ---------------------- 10GeV ----------------------------------
  TH1D  *Hsig1_vTcPL_Ceex2  = HstProjV("Hsig1_vTcPL_Ceex2", sct1_vTcPL_Ceex2, nbMax);
  TH1D  *Hafb1_vTcPL_Ceex2  = HstProjA("Hafb1_vTcPL_Ceex2", sct1_vTcPL_Ceex2, nbMax);
  //******** IFI off
  TH1D  *Hsig1_vTcPL_Ceex2n  = HstProjV("Hsig1_vTcPL_Ceex2n", sct1_vTcPL_Ceex2n, nbMax);
  TH1D  *Hafb1_vTcPL_Ceex2n  = HstProjA("Hafb1_vTcPL_Ceex2n", sct1_vTcPL_Ceex2n, nbMax);

//======================================================================================
// FOAM corner
  int    NbMax   =0;          // for 100bins, default=0 for gCosTheta = 1.00
  NbMax=45;      // cosThetaMax = 45/50=0.90
  //////////////  95GeV /////////////////
  TH1D *HST_FOAM_NORMA3_95 = (TH1D*)DiskFileF95.Get("HST_FOAM_NORMA3");
  TH1D *HST_FOAM_NORMA5_95 = (TH1D*)DiskFileF95.Get("HST_FOAM_NORMA5");

  TH2D *SCT9_xc_Ceex2n = (TH2D*)DiskFileF95.Get("SCT_xc_Ceex2n");  // FOAM small range x<0.20
  TH2D *SCT9_xc_Ceex2  = (TH2D*)DiskFileF95.Get("SCT_xc_Ceex2");   // FOAM small range x<0.20

  HisNorm2(HST_FOAM_NORMA3_95, SCT9_xc_Ceex2n );  // normalizing
  HisNorm2(HST_FOAM_NORMA5_95, SCT9_xc_Ceex2 );   // normalizing

  // Foam3
  TH1D  *HHtot9_xmax_Ceex2n = HstProjV("HHtot9_xmax_Ceex2n",SCT9_xc_Ceex2n,NbMax);
  TH1D  *HHafb9_xmax_Ceex2n = HstProjA("HHafb9_xmax_Ceex2n",SCT9_xc_Ceex2n,NbMax);
  // Foam5
  TH1D  *HHtot9_xmax_Ceex2  = HstProjV("HHtot9_xmax_Ceex2",SCT9_xc_Ceex2,NbMax);
  TH1D  *HHafb9_xmax_Ceex2  = HstProjA("HHafb9_xmax_Ceex2",SCT9_xc_Ceex2,NbMax);

  //////////////  88GeV /////////////////
  TH1D *HST_FOAM_NORMA3_88 = (TH1D*)DiskFileF88.Get("HST_FOAM_NORMA3");
  TH1D *HST_FOAM_NORMA5_88 = (TH1D*)DiskFileF88.Get("HST_FOAM_NORMA5");

  TH2D *SCT8_xc_Ceex2n = (TH2D*)DiskFileF88.Get("SCT_xc_Ceex2n");  // FOAM small range x<0.20
  TH2D *SCT8_xc_Ceex2  = (TH2D*)DiskFileF88.Get("SCT_xc_Ceex2");   // FOAM small range x<0.20

  HisNorm2(HST_FOAM_NORMA3_88, SCT8_xc_Ceex2n );  // normalizing
  HisNorm2(HST_FOAM_NORMA5_88, SCT8_xc_Ceex2 );   // normalizing

  // Foam3
  TH1D  *HHtot8_xmax_Ceex2n = HstProjV("HHtot8_xmax_Ceex2n",SCT8_xc_Ceex2n,NbMax);
  TH1D  *HHafb8_xmax_Ceex2n = HstProjA("HHafb8_xmax_Ceex2n",SCT8_xc_Ceex2n,NbMax);
  HHafb8_xmax_Ceex2n->Scale(-1.0);
  // Foam5
  TH1D  *HHtot8_xmax_Ceex2  = HstProjV("HHtot8_xmax_Ceex2",SCT8_xc_Ceex2,NbMax);
  TH1D  *HHafb8_xmax_Ceex2  = HstProjA("HHafb8_xmax_Ceex2",SCT8_xc_Ceex2,NbMax);
  HHafb8_xmax_Ceex2->Scale(-1.0);


}//ReMakeMChisto


///////////////////////////////////////////////////////////////////////////////////
void FigTempl()
{
//------------------------------------------------------------------------
  cout<<" ========================= FigTempl =========================== "<<endl;
  // renormalize histograms in nanobarns
  Double_t CMSene;
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA95.Get("HST_KKMC_NORMA");
  CMSene  = HST_KKMC_NORMA->GetBinContent(1); // CMSene=xpar(1) stored in NGeISR
  //
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cTempl = new TCanvas("cTempl","cTempl", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cTempl->SetFillColor(10);
  cTempl->Divide( 2,  0);
  //cTempl->Divide( 2,  0,     0.0,     0.0,   10);
  //              nx, ny, xmargin, ymargin, color
  //////////////////////////////////////////////
  cTempl->cd(1);
  CaptT->DrawLatex(0.12,0.95,"A_{FB}(v_{max}), ????");
  //-------------------------------------
  cTempl->cd(2);
  //
  cTempl->cd();
//
}// FigTempl



///////////////////////////////////////////////////////////////////////////////////
void AfbIFIvA2()
{
//------------------------------------------------------------------------
  cout<<" ========================= AfbIFIvA2 =========================== "<<endl;
  ////////////////////////////////////////////////////////////////////////////////
  TH1D *Hafb9_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("Hafb9_vAcPR_Ceex2");
  TH1D *Hafb9_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb9_vAcPR_Ceex2n");
  //
  TH1D *Hafb8_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("Hafb8_vAcPR_Ceex2");
  TH1D *Hafb8_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb8_vAcPR_Ceex2n");
  //
  TH1D *HafbZ_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("HafbZ_vAcPR_Ceex2");
  TH1D *HafbZ_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("HafbZ_vAcPR_Ceex2n");
  //
  TH1D *Hafb1_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("Hafb1_vAcPR_Ceex2");
  TH1D *Hafb1_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb1_vAcPR_Ceex2n");
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cAfbIFIvA2 = new TCanvas("cAfbIFIvA2","cAfbIFIvA2", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cAfbIFIvA2->SetFillColor(10);
  cAfbIFIvA2->Divide( 2,  1);
//*****************************************************************************
  cAfbIFIvA2->cd(1);
  Hafb9_vAcPR_Ceex2->SetTitle(0);
  Hafb9_vAcPR_Ceex2->SetStats(0);
  Hafb9_vAcPR_Ceex2->GetXaxis()->SetTitle("v_{max}");
  //Hafb9_vAcPR_Ceex2->GetYaxis()->SetTitle("A_{FB}(v_{max})");
  Hafb9_vAcPR_Ceex2->SetLineColor(kBlue);
  Hafb9_vAcPR_Ceex2->SetMaximum( 0.33);
  Hafb9_vAcPR_Ceex2->SetMinimum( 0.15);
  Hafb9_vAcPR_Ceex2->DrawCopy("h");
  //
  Hafb9_vAcPR_Ceex2n->SetLineColor(kBlack);
  Hafb9_vAcPR_Ceex2n->DrawCopy("hsame");
  //
  Hafb8_vAcPR_Ceex2->SetLineColor(kBlue);
  Hafb8_vAcPR_Ceex2->DrawCopy("hsame");
  Hafb8_vAcPR_Ceex2n->SetLineColor(kBlack);
  Hafb8_vAcPR_Ceex2n->DrawCopy("hsame");

  CaptT->DrawLatex(0.02,0.95," Black=IFIoff,  Blue=IFIon, v=v_{ALEPH}");
  CaptT->DrawLatex(0.50,0.20,"  A_{FB}(v_{max}), #sqrt{s}=94.3GeV ");
  CaptT->DrawLatex(0.50,0.83," -A_{FB}(v_{max}), #sqrt{s}=87.9GeV ");
  //*****************************************************************************
  cAfbIFIvA2->cd(2);
  TH1D *hZero = (TH1D*)Hafb8_vAcPR_Ceex2n->Clone("hZero");  // zero line
  for(int i=1; i <= hZero->GetNbinsX() ; i++) { hZero->SetBinContent(i, 0); hZero->SetBinError(i, 0);}

  TH1D *Hafb9_vAcPR_IFIdiff = HstDiff("Hafb9_vAcPR_IFIdiff",  Hafb9_vAcPR_Ceex2, Hafb9_vAcPR_Ceex2n,  kBlack);
  TH1D *Hafb8_vAcPR_IFIdiff = HstDiff("Hafb8_vAcPR_IFIdiff",  Hafb8_vAcPR_Ceex2, Hafb8_vAcPR_Ceex2n,  kBlue);
  Hafb8_vAcPR_IFIdiff->Scale(-1.0); // undoing sign change
  TH1D *HafbZ_vAcPR_IFIdiff = HstDiff("HafbZ_vAcPR_IFIdiff",  HafbZ_vAcPR_Ceex2, HafbZ_vAcPR_Ceex2n,  kMagenta);
  TH1D *Hafb1_vAcPR_IFIdiff = HstDiff("Hafb1_vAcPR_IFIdiff",  Hafb1_vAcPR_Ceex2, Hafb1_vAcPR_Ceex2n,  kPine);
  TH1D *HDifSum             = HstDiff("HDifSum",            Hafb9_vAcPR_IFIdiff, Hafb8_vAcPR_IFIdiff,  kRed);
  //HDifSum->SetLineWidth(2);

  TH1D *Ddiff = Hafb9_vAcPR_IFIdiff;
//  TH1D *Ddiff = Hafb1_vAcPR_IFIdiff;
  Ddiff->SetTitle(0);
  Ddiff->SetStats(0);
  //Ddiff->GetYaxis()->SetTitle("#Delta A^{IFI}_{FB}(v_{max})");

  Ddiff->SetMaximum( 0.05); Ddiff->SetMinimum(-0.03);
  Ddiff->DrawCopy("h");

  CaptT->SetTextColor(kBlack);
  CaptT->DrawLatex(0.01, 0.95, " A^{IFIon}_{FB}(v_{max}) - A^{IFIoff}_{FB}(v_{max}) ");
  double ycapt =0.33;
  PlotSame2(Hafb9_vAcPR_IFIdiff, ycapt, kBlack,   0.020, "(a)", "#sqrt{s}=94.3GeV");
  PlotSame2(Hafb8_vAcPR_IFIdiff, ycapt, kBlue,    0.045, "(b)", "#sqrt{s}=87.9GeV");
  if( gTogEne ){
  PlotSame2(HafbZ_vAcPR_IFIdiff, ycapt, kMagenta, 0.060, "(c)", "#sqrt{s}=M_{Z}");
  PlotSame2(Hafb1_vAcPR_IFIdiff, ycapt, kPine,   0.090, "(d)", "#sqrt{s}=10GeV");
  }
  PlotSame2(HDifSum,             ycapt, kRed,     0.030, "(e)", "= (a) - (b) ");
  //
  hZero->DrawCopy("hsame");
  cAfbIFIvA2->cd();
//
}// AfbIFIvA2


///////////////////////////////////////////////////////////////////////////////////
void AfbIFIvA1()
{
//------------------------------------------------------------------------
  cout<<" ========================= AfbIFIvA1 =========================== "<<endl;
  ////////////////////////////////////////////////////////////////////////////////
  TH1D *Hafb9_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("Hafb9_vAcPR_Ceex2");
  TH1D *Hafb9_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb9_vAcPR_Ceex2n");
  //
  TH1D *Hafb8_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("Hafb8_vAcPR_Ceex2");
  TH1D *Hafb8_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb8_vAcPR_Ceex2n");
  //
  TH1D *HafbZ_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("HafbZ_vAcPR_Ceex2");
  TH1D *HafbZ_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("HafbZ_vAcPR_Ceex2n");
  //
  TH1D *Hafb1_vAcPR_Ceex2      = (TH1D*)DiskFileB.Get("Hafb1_vAcPR_Ceex2");
  TH1D *Hafb1_vAcPR_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb1_vAcPR_Ceex2n");

  TH1D *Hafb9_vAcPR_IFIdiff = HstDiff("Hafb9_vAcPR_IFIdiff",  Hafb9_vAcPR_Ceex2, Hafb9_vAcPR_Ceex2n,  kBlack);
  TH1D *Hafb8_vAcPR_IFIdiff = HstDiff("Hafb8_vAcPR_IFIdiff",  Hafb8_vAcPR_Ceex2, Hafb8_vAcPR_Ceex2n,  kBlue);
  Hafb8_vAcPR_IFIdiff->Scale(-1.0); // undoing sign change
  TH1D *HafbZ_vAcPR_IFIdiff = HstDiff("HafbZ_vAcPR_IFIdiff",  HafbZ_vAcPR_Ceex2, HafbZ_vAcPR_Ceex2n,  kMagenta);
  TH1D *Hafb1_vAcPR_IFIdiff = HstDiff("Hafb1_vAcPR_IFIdiff",  Hafb1_vAcPR_Ceex2, Hafb1_vAcPR_Ceex2n,  kPine);

  TH1D *hZero = (TH1D*)Hafb8_vAcPR_Ceex2n->Clone("hZero");  // zero line
  for(int i=1; i <= hZero->GetNbinsX() ; i++) { hZero->SetBinContent(i, 0); hZero->SetBinError(i, 0);}

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cAfbIFIvA1 = new TCanvas("cAfbIFIvA1","cAfbIFIvA1", gXcanv,  gYcanv,   600,  600);
  //                                   Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cAfbIFIvA1->SetFillColor(10);

  TH1D *Ddiff = Hafb9_vAcPR_IFIdiff;
  Ddiff->SetTitle(0);
  Ddiff->SetStats(0);

  Ddiff->SetMaximum( 0.05); Ddiff->SetMinimum(-0.03);
  Ddiff->DrawCopy("h");

  CaptT->SetTextColor(kBlack);
  CaptT->DrawLatex(0.01, 0.95, " A^{IFIon}_{FB}(v_{max}) - A^{IFIoff}_{FB}(v_{max}) ");
  double ycapt =0.33;
  PlotSame2(Hafb9_vAcPR_IFIdiff, ycapt, kBlack,   0.020, "(a)", "#sqrt{s}=94.3GeV");
  PlotSame2(Hafb8_vAcPR_IFIdiff, ycapt, kBlue,    0.045, "(b)", "#sqrt{s}=87.9GeV");
  if( gTogEne ){
  PlotSame2(HafbZ_vAcPR_IFIdiff, ycapt, kMagenta, 0.060, "(c)", "#sqrt{s}=M_{Z}");
  PlotSame2(Hafb1_vAcPR_IFIdiff, ycapt, kPine,   0.090, "(d)", "#sqrt{s}=10GeV");
  }
  hZero->DrawCopy("hsame");

  cAfbIFIvA1->SaveAs("cAfbIFIvA1.pdf");
  //
  }// AfbIFIvA1


///////////////////////////////////////////////////////////////////////////////////
void AfbIFIvT1()
{
//------------------------------------------------------------------------
  cout<<" ========================= AfbIFIvT1 =========================== "<<endl;
  ////////////////////////////////////////////////////////////////////////////////
  TH1D *Hafb9_vTcPL_Ceex2      = (TH1D*)DiskFileB.Get("Hafb9_vTcPL_Ceex2");
  TH1D *Hafb9_vTcPL_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb9_vTcPL_Ceex2n");
  //
  TH1D *Hafb8_vTcPL_Ceex2      = (TH1D*)DiskFileB.Get("Hafb8_vTcPL_Ceex2");
  TH1D *Hafb8_vTcPL_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb8_vTcPL_Ceex2n");
  //
  TH1D *HafbZ_vTcPL_Ceex2      = (TH1D*)DiskFileB.Get("HafbZ_vTcPL_Ceex2");
  TH1D *HafbZ_vTcPL_Ceex2n     = (TH1D*)DiskFileB.Get("HafbZ_vTcPL_Ceex2n");
  //
  TH1D *Hafb1_vTcPL_Ceex2      = (TH1D*)DiskFileB.Get("Hafb1_vTcPL_Ceex2");
  TH1D *Hafb1_vTcPL_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb1_vTcPL_Ceex2n");
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cAfbIFIvT1 = new TCanvas("cAfbIFIvT1","cAfbIFIvT1", gXcanv,  gYcanv,   600,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cAfbIFIvT1->SetFillColor(10);
//*****************************************************************************
  TH1D *hZero = (TH1D*)Hafb8_vTcPL_Ceex2n->Clone("hZero");  // zero line
  for(int i=1; i <= hZero->GetNbinsX() ; i++) { hZero->SetBinContent(i, 0); hZero->SetBinError(i, 0);}

  TH1D *Hafb9_vTcPL_IFIdiff = HstDiff("Hafb9_vTcPL_IFIdiff",  Hafb9_vTcPL_Ceex2, Hafb9_vTcPL_Ceex2n,  kBlack);
  TH1D *Hafb8_vTcPL_IFIdiff = HstDiff("Hafb8_vTcPL_IFIdiff",  Hafb8_vTcPL_Ceex2, Hafb8_vTcPL_Ceex2n,  kBlue);
  Hafb8_vTcPL_IFIdiff->Scale(-1.0); // undoing sign change
  TH1D *HafbZ_vTcPL_IFIdiff = HstDiff("HafbZ_vTcPL_IFIdiff",  HafbZ_vTcPL_Ceex2, HafbZ_vTcPL_Ceex2n,  kMagenta);
  TH1D *Hafb1_vTcPL_IFIdiff = HstDiff("Hafb1_vTcPL_IFIdiff",  Hafb1_vTcPL_Ceex2, Hafb1_vTcPL_Ceex2n,  kPine);
  TH1D *HDifPat_vTcPL       = HstDiff("HDifPat_vTcPL",        Hafb9_vTcPL_IFIdiff, Hafb8_vTcPL_IFIdiff,  kRed);
  //HDifPat_vTcPL->SetLineWidth(2);

  TH1D *Ddiff = Hafb9_vTcPL_IFIdiff;
  Ddiff->SetTitle(0);
  Ddiff->SetStats(0);
  //Ddiff->GetYaxis()->SetTitle("#Delta A^{IFI}_{FB}(v_{max})");

  Ddiff->SetMaximum( 0.05); Ddiff->SetMinimum(-0.03);
  Ddiff->DrawCopy("h");

  CaptT->SetTextColor(kBlack);
  CaptT->DrawLatex(0.01, 0.95, " A^{IFI}_{FB}(v_{max}) = A^{IFIon}_{FB}(v_{max}) - A^{IFIoff}_{FB}(v_{max}) ");
  double ycapt =0.33;
  PlotSame2(Hafb9_vTcPL_IFIdiff, ycapt, kBlack,   0.040, "(a)", "#sqrt{s}=94.3GeV");
  PlotSame2(Hafb8_vTcPL_IFIdiff, ycapt, kBlue,    0.065, "(b)", "#sqrt{s}=87.9GeV");
  if( gTogEne ){
  PlotSame2(HafbZ_vTcPL_IFIdiff, ycapt, kMagenta, 0.060, "(c)", "#sqrt{s}=M_{Z}");
  PlotSame2(Hafb1_vTcPL_IFIdiff, ycapt, kPine,    0.090, "(d)", "#sqrt{s}=10GeV");
  }
  //PlotSame2(HDifPat_vTcPL,             ycapt, kRed,     0.030, "(e)", "= (a) - (b) ");
  //
  hZero->DrawCopy("hsame");
  cAfbIFIvT1->cd();
  //
  cAfbIFIvT1->SaveAs("cAfbIFIvT1.pdf");
//
}// AfbIFIvT1



///////////////////////////////////////////////////////////////////////////////////
void AfbIFIvT2()
{
//------------------------------------------------------------------------
  cout<<" ========================= AfbIFIvT2 =========================== "<<endl;
  ////////////////////////////////////////////////////////////////////////////////
  TH1D *Hafb9_vTcPL_Ceex2      = (TH1D*)DiskFileB.Get("Hafb9_vTcPL_Ceex2");
  TH1D *Hafb9_vTcPL_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb9_vTcPL_Ceex2n");
  //
  TH1D *Hafb8_vTcPL_Ceex2      = (TH1D*)DiskFileB.Get("Hafb8_vTcPL_Ceex2");
  TH1D *Hafb8_vTcPL_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb8_vTcPL_Ceex2n");
  //
  TH1D *HafbZ_vTcPL_Ceex2      = (TH1D*)DiskFileB.Get("HafbZ_vTcPL_Ceex2");
  TH1D *HafbZ_vTcPL_Ceex2n     = (TH1D*)DiskFileB.Get("HafbZ_vTcPL_Ceex2n");
  //
  TH1D *Hafb1_vTcPL_Ceex2      = (TH1D*)DiskFileB.Get("Hafb1_vTcPL_Ceex2");
  TH1D *Hafb1_vTcPL_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb1_vTcPL_Ceex2n");

  TH1D *Hafb9m8_IFI_ceex2   = HstAddi("Hafb9m8_IFI_ceex2",  Hafb9_vTcPL_Ceex2,  Hafb8_vTcPL_Ceex2,  kBlack);
  TH1D *Hafb9m8_IFI_ceex2n  = HstAddi("Hafb9m8_IFI_ceex2n", Hafb9_vTcPL_Ceex2n, Hafb8_vTcPL_Ceex2n, kBlack);
  Hafb9m8_IFI_ceex2->Scale(0.5);
  Hafb9m8_IFI_ceex2n->Scale(0.5);

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cAfbIFIvT2 = new TCanvas("cAfbIFIvT2","cAfbIFIvT2", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cAfbIFIvT2->SetFillColor(10);
  cAfbIFIvT2->Divide( 2,  1);
//*****************************************************************************
  cAfbIFIvT2->cd(1);
  TH1D *HST = Hafb9_vTcPL_Ceex2;
  HST->SetTitle(0); HST->SetStats(0);
  HST->GetXaxis()->SetTitle("v_{max}");
  //
  HST->SetMaximum( 0.32); HST->SetMinimum( 0.15);
  HST->DrawCopy("h");

  CaptT->DrawLatex(0.02,0.95," KKMC: A_{FB}(v_{max}), |cos(#theta|<0.9");
  CaptT->DrawLatex(0.50,0.78," #sqrt{s_{-}}=87.9GeV ");
  CaptT->DrawLatex(0.50,0.22," #sqrt{s_{+}}=94.3GeV ");

  double ycapt =0.70;
  PlotSame2(Hafb9_vTcPL_Ceex2, ycapt, kBlue,  0.010, "[a+]", "= A_{FB}^{IFIon}(v_{max},s_{+})");
  PlotSame2(Hafb8_vTcPL_Ceex2, ycapt, kBlue,  0.050, "[a-]", "=-A_{FB}^{IFIon}(v_{max},s_{-})");
  //
  PlotSame2(Hafb9_vTcPL_Ceex2n,ycapt, kBlack, 0.010, "[b+]", "= A_{FB}^{IFIoff}(v_{max},s_{+})");
  PlotSame2(Hafb8_vTcPL_Ceex2n,ycapt, kBlack, 0.050, "[b-]", "=-A_{FB}^{IFIoff}(v_{max},s_{-})");
  //
  ycapt = 0.40;
  PlotSame2(Hafb9m8_IFI_ceex2,ycapt,  kBlue,   0.065, "(A)", "=([a+]+[a-])/2, IFI on");
  PlotSame2(Hafb9m8_IFI_ceex2n,ycapt, kBlack,  0.055, "(B)", "=([b+]+[b-])/2, IFI off");
  //*****************************************************************************
  cAfbIFIvT2->cd(2);
  TH1D *hZero = (TH1D*)Hafb8_vTcPL_Ceex2n->Clone("hZero");  // zero line
  for(int i=1; i <= hZero->GetNbinsX() ; i++) { hZero->SetBinContent(i, 0); hZero->SetBinError(i, 0);}

  TH1D *Hafb9_vTcPL_IFIdiff = HstDiff("Hafb9_vTcPL_IFIdiff",  Hafb9_vTcPL_Ceex2, Hafb9_vTcPL_Ceex2n,  kBlack);
  TH1D *Hafb8_vTcPL_IFIdiff = HstDiff("Hafb8_vTcPL_IFIdiff",  Hafb8_vTcPL_Ceex2, Hafb8_vTcPL_Ceex2n,  kBlue);
  Hafb8_vTcPL_IFIdiff->Scale(-1.0); // undoing sign change
  TH1D *HafbZ_vTcPL_IFIdiff = HstDiff("HafbZ_vTcPL_IFIdiff",  HafbZ_vTcPL_Ceex2, HafbZ_vTcPL_Ceex2n,  kMagenta);
  TH1D *Hafb1_vTcPL_IFIdiff = HstDiff("Hafb1_vTcPL_IFIdiff",  Hafb1_vTcPL_Ceex2, Hafb1_vTcPL_Ceex2n,  kPine);
  TH1D *HDifPat_vTcPL       = HstDiff("HDifPat_vTcPL",        Hafb9_vTcPL_IFIdiff, Hafb8_vTcPL_IFIdiff,  kRed);
  //HDifPat_vTcPL->SetLineWidth(2);

  TH1D *Ddiff = Hafb9_vTcPL_IFIdiff;
  Ddiff->SetTitle(0);
  Ddiff->SetStats(0);
  //Ddiff->GetYaxis()->SetTitle("#Delta A^{IFI}_{FB}(v_{max})");

  Ddiff->SetMaximum( 0.05); Ddiff->SetMinimum(-0.03);
  Ddiff->DrawCopy("h");

  CaptT->SetTextColor(kBlack);
  CaptT->DrawLatex(0.01, 0.95,
    "A^{IFI}_{FB}(v_{max}) = A^{IFIon}_{FB}(v_{max}) - A^{IFIoff}_{FB}(v_{max}), |cos(#theta|<0.9");
  ycapt =0.33;
  PlotSame2(Hafb9_vTcPL_IFIdiff, ycapt, kBlack,   0.040, "(a)", "= [a+]-[b+],  #sqrt{s}=94.3GeV");
  PlotSame2(Hafb8_vTcPL_IFIdiff, ycapt, kBlue,    0.065, "(b)", "=-[a-]+[b-],  #sqrt{s}=87.9GeV");
  if( gTogEne ){
  PlotSame2(HafbZ_vTcPL_IFIdiff, ycapt, kMagenta, 0.060, "(c)", "#sqrt{s}=M_{Z}");
  PlotSame2(Hafb1_vTcPL_IFIdiff, ycapt, kPine,    0.090, "(d)", "#sqrt{s}=10GeV");
  }
  PlotSame2(HDifPat_vTcPL,       ycapt, kRed,     0.030, "(e)", "= (a) - (b) = (A) - (B)");
  //
  hZero->DrawCopy("hsame");
  cAfbIFIvT2->cd();
  //
  cAfbIFIvT2->SaveAs("cAfbIFIvT2.pdf");
//
}// AfbIFIvT2




///////////////////////////////////////////////////////////////////////////////////
void AfbIFI_Foam()
{
//------------------------------------------------------------------------
  cout<<" ========================= AfbIFI_Foam =========================== "<<endl;
  ////////////////////////////////////////////////////////////////////////////////
  TH1D *Hafb9_xmax_Ceex2      = (TH1D*)DiskFileB.Get("HHafb9_xmax_Ceex2");
  TH1D *Hafb9_xmax_Ceex2n     = (TH1D*)DiskFileB.Get("HHafb9_xmax_Ceex2n");
  //HHafb9_xmax_Ceex2n
  TH1D *Hafb8_xmax_Ceex2      = (TH1D*)DiskFileB.Get("HHafb8_xmax_Ceex2");
  TH1D *Hafb8_xmax_Ceex2n     = (TH1D*)DiskFileB.Get("HHafb8_xmax_Ceex2n");
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cAfbIFI_Foam = new TCanvas("cAfbIFI_Foam","cAfbIFI_Foam", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cAfbIFI_Foam->SetFillColor(10);
  cAfbIFI_Foam->Divide( 2,  1);
//*****************************************************************************
  cAfbIFI_Foam->cd(1);
  Hafb9_xmax_Ceex2->SetTitle(0);
  Hafb9_xmax_Ceex2->SetStats(0);
  Hafb9_xmax_Ceex2->GetXaxis()->SetTitle("v_{max}");
  Hafb9_xmax_Ceex2->DrawCopy("h");

  double ycapt =0.70;
  PlotSame2(Hafb9_xmax_Ceex2, ycapt, kBlue,   0.010, "(a)", "A_{FB}(v_{max},s_{+}), IFI on");
  PlotSame2(Hafb9_xmax_Ceex2n,ycapt, kBlack,  0.010, "(b)", "A_{FB}(v_{max},s_{+}), IFI off");
  //
  PlotSame2(Hafb8_xmax_Ceex2, ycapt, kBlue,   0.050, "(c)", "-A_{FB}(v_{max},s_{-}), IFI on");
  PlotSame2(Hafb8_xmax_Ceex2n,ycapt, kBlack,  0.050, "(d)", "-A_{FB}(v_{max},s_{-}), IFI off");


  CaptT->DrawLatex(0.02,0.95," KKFoam: A_{FB}(v_{max}), |cos(#theta|<0.9");
  CaptT->DrawLatex(0.50,0.25," #sqrt{s_{+}}=94.3GeV ");
  CaptT->DrawLatex(0.50,0.83," #sqrt{s_{-}}=87.9GeV ");
  //*****************************************************************************
  cAfbIFI_Foam->cd(2);
  TH1D *hZero = (TH1D*)Hafb8_xmax_Ceex2n->Clone("hZero");  // zero line
  for(int i=1; i <= hZero->GetNbinsX() ; i++) { hZero->SetBinContent(i, 0); hZero->SetBinError(i, 0);}

  TH1D *Hafb9_xmax_IFIdiff = HstDiff("Hafb9_xmax_IFIdiff",  Hafb9_xmax_Ceex2, Hafb9_xmax_Ceex2n,  kBlack);
  TH1D *Hafb8_xmax_IFIdiff = HstDiff("Hafb8_xmax_IFIdiff",  Hafb8_xmax_Ceex2, Hafb8_xmax_Ceex2n,  kBlue);
  Hafb8_xmax_IFIdiff->Scale(-1.0); // undoing sign change
  TH1D *HDifPat_xmax             = HstDiff("HDifPat_xmax",            Hafb9_xmax_IFIdiff, Hafb8_xmax_IFIdiff,  kRed);
  //HDifPat_xmax->SetLineWidth(2);

  TH1D *Ddiff = Hafb9_xmax_IFIdiff;
  Ddiff->SetTitle(0);
  Ddiff->SetStats(0);
  //Ddiff->GetYaxis()->SetTitle("#Delta A^{IFI}_{FB}(v_{max})");

  Ddiff->SetMaximum( 0.05); Ddiff->SetMinimum(-0.03);
  Ddiff->DrawCopy("h");

  CaptT->SetTextColor(kBlack);
  CaptT->DrawLatex(0.01, 0.95, " A^{IFI}_{FB}(v_{max})= A^{IFIon}_{FB}(v_{max}) - A^{IFIoff}_{FB}(v_{max}) ");
  ycapt =0.33;
  PlotSame2(Hafb9_xmax_IFIdiff, ycapt, kBlack,   0.040, "(a)", "#sqrt{s}=94.3GeV");
  PlotSame2(Hafb8_xmax_IFIdiff, ycapt, kBlue,    0.065, "(b)", "#sqrt{s}=87.9GeV");
  PlotSame2(HDifPat_xmax,       ycapt, kRed,     0.030, "(e)", "= (a) - (b) ");
  //
  hZero->DrawCopy("hsame");
  cAfbIFI_Foam->cd();
//
}// AfbIFI_Foam


///////////////////////////////////////////////////////////////////////////////////
void AfbDifPat()
{
//------------------------------------------------------------------------
  cout<<" ========================= AfbDifPat =========================== "<<endl;
  ////////////////////////////////////////////////////////////////////////////////
  TH1D *HDifPat_vTcPL     = (TH1D*)DiskFileB.Get("HDifPat_vTcPL");
  TH1D *HDifPat_xmax      = (TH1D*)DiskFileB.Get("HDifPat_xmax");
  TH1D *Hafb9m8_IFI_ceex2 = (TH1D*)DiskFileB.Get("Hafb9m8_IFI_ceex2");
  Hafb9m8_IFI_ceex2->Scale(2.0);

  TH1D *Hafb_IFI_diff_Pat = HstDiff("Hafb_IFI_diff_Pat", HDifPat_vTcPL, HDifPat_xmax,  kRed);
  //
  TH1D *Hafb_IFI_rat_Pat= HstRatio("Hafb_IFI_rat_Pat", Hafb_IFI_diff_Pat, Hafb9m8_IFI_ceex2, kBlack);
  //
  HDifPat_vTcPL->Scale(0.10);
  HDifPat_xmax->Scale(0.10);
  //
  TH1D *hZero7 = (TH1D*)HDifPat_xmax->Clone("hZero7");  // zero line
  for(int i=1; i <= hZero7->GetNbinsX() ; i++) { hZero7->SetBinContent(i, 0); hZero7->SetBinError(i, 0);}
  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);

  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cAfbDifPat = new TCanvas("cAfbDifPat","cAfbDifPat", gXcanv,  gYcanv,   600,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cAfbDifPat->SetFillColor(10);

  cAfbDifPat->cd();

  TH1D *HDif1 = HDifPat_vTcPL;
  HDif1->SetMaximum( 0.0006); HDif1->SetMinimum(-0.0012);

  HDif1->SetTitle(0); HDif1->SetStats(0);
  HDif1->GetXaxis()->SetTitle("v_{max}");

  HDif1->DrawCopy("h");

  double ycapt =0.90;
  PlotSame2(HDifPat_vTcPL,     ycapt, kBlack, 0.155, "(a)", "#Delta A_{FB}^{IFI}x10^{-1},  KKMC");
  PlotSame2(HDifPat_xmax,      ycapt, kBlue,  0.120, "(b)", "#Delta A_{FB}^{IFI}x10^{-1},  KKFoam");
  PlotSame2(Hafb_IFI_diff_Pat, ycapt, kRed,   0.010, "(c)", "#Delta A_{FB}^{IFI}(a) - #Delta A_{FB}^{IFI}(b)");
  //
  PlotSame2(Hafb_IFI_rat_Pat,  ycapt, kPine,  0.010, "(d)", "#Delta A_{FB}^{IFI}(c)/(A_{FB}(+)-A_{FB}(-))");

  HDifPat_xmax->DrawCopy("hsame");

  hZero7->SetLineColor(kBlack);
  hZero7->DrawCopy("hsame");

  CaptT->DrawLatex(0.12,0.95,"#Delta A_{FB}^{IFI}(v_{max}) = A_{FB}^{IFI}(v_{max},s_{+} ) - A_{FB}^{IFI}(v_{max},s_{-} )");
  //
  cAfbDifPat->cd();
  cAfbDifPat->SaveAs("cAfbDifPat.pdf");

//
}// AfbDifPat


///////////////////////////////////////////////////////////////////////////////////
void AfbIFI_KKmc2()
{
//------------------------------------------------------------------------
  cout<<" ========================= AfbIFI_KKmc2 =========================== "<<endl;
  ////////////////////////////////////////////////////////////////////////////////
  TH1D *Hafb9_vTcPL_Ceex2      = (TH1D*)DiskFileB.Get("Hafb9_vTcPL_Ceex2");
  TH1D *Hafb9_vTcPL_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb9_vTcPL_Ceex2n");
  //
  TH1D *Hafb8_vTcPL_Ceex2      = (TH1D*)DiskFileB.Get("Hafb8_vTcPL_Ceex2");
  TH1D *Hafb8_vTcPL_Ceex2n     = (TH1D*)DiskFileB.Get("Hafb8_vTcPL_Ceex2n");
  ////////
  TH1D *Hafb9_vTcPL_Ceex1      = (TH1D*)DiskFileB.Get("Hafb9_vTcPL_Ceex1");
  TH1D *Hafb9_vTcPL_Ceex1n     = (TH1D*)DiskFileB.Get("Hafb9_vTcPL_Ceex1n");
  //
  TH1D *Hafb8_vTcPL_Ceex1      = (TH1D*)DiskFileB.Get("Hafb8_vTcPL_Ceex1");

  TH1D *Hafb8_vTcPL_Ceex1n     = (TH1D*)DiskFileB.Get("Hafb8_vTcPL_Ceex1n");
  ////////
  TH1D *Hafb9_vTcPL_Ceex0      = (TH1D*)DiskFileB.Get("Hafb9_vTcPL_Ceex0");
  TH1D *Hafb9_vTcPL_Ceex0n     = (TH1D*)DiskFileB.Get("Hafb9_vTcPL_Ceex0n");
  //
  TH1D *Hafb8_vTcPL_Ceex0      = (TH1D*)DiskFileB.Get("Hafb8_vTcPL_Ceex0");
  TH1D *Hafb8_vTcPL_Ceex0n     = (TH1D*)DiskFileB.Get("Hafb8_vTcPL_Ceex0n");
  ////////
  TH1D *Hafb9_vTcPL_IFI2 = HstDiff("Hafb9_vTcPL_IFI2",  Hafb9_vTcPL_Ceex2, Hafb9_vTcPL_Ceex2n,  kBlack);
  TH1D *Hafb8_vTcPL_IFI2 = HstDiff("Hafb8_vTcPL_IFI2",  Hafb8_vTcPL_Ceex2, Hafb8_vTcPL_Ceex2n,  kBlue);
  Hafb8_vTcPL_IFI2->Scale(-1.0); // undoing sign change
  //
  TH1D *Hafb9_vTcPL_IFI1 = HstDiff("Hafb9_vTcPL_IFI1",  Hafb9_vTcPL_Ceex1, Hafb9_vTcPL_Ceex1n,  kBlack);
  TH1D *Hafb8_vTcPL_IFI1 = HstDiff("Hafb8_vTcPL_IFI1",  Hafb8_vTcPL_Ceex1, Hafb8_vTcPL_Ceex1n,  kBlue);
  Hafb8_vTcPL_IFI1->Scale(-1.0); // undoing sign change
  //
  TH1D *Hafb9_vTcPL_IFI0 = HstDiff("Hafb9_vTcPL_IFI0",  Hafb9_vTcPL_Ceex0, Hafb9_vTcPL_Ceex0n,  kBlack);
  TH1D *Hafb8_vTcPL_IFI0 = HstDiff("Hafb8_vTcPL_IFI0",  Hafb8_vTcPL_Ceex0, Hafb8_vTcPL_Ceex0n,  kBlue);
  Hafb8_vTcPL_IFI0->Scale(-1.0); // undoing sign change
  //
  TH1D *Hafb9m8_vTcPL_IFI2 = HstDiff("Hafb9m8_vTcPL_IFI2",  Hafb9_vTcPL_IFI2,  Hafb8_vTcPL_IFI2,  kBlack);
  TH1D *Hafb9m8_vTcPL_IFI1 = HstDiff("Hafb9m8_vTcPL_IFI1",  Hafb9_vTcPL_IFI1,  Hafb8_vTcPL_IFI1,  kBlack);
  TH1D *Hafb9m8_vTcPL_IFI0 = HstDiff("Hafb9m8_vTcPL_IFI0",  Hafb9_vTcPL_IFI0,  Hafb8_vTcPL_IFI0,  kBlack);

  TH1D *hZero7      = (TH1D*)DiskFileB.Get("hZero7");

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cAfbIFI_KKmc2 = new TCanvas("cAfbIFI_KKmc2","cAfbIFI_KKmc2", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cAfbIFI_KKmc2->SetFillColor(10);
  cAfbIFI_KKmc2->Divide( 2,  1);
//*****************************************************************************
  cAfbIFI_KKmc2->cd(1);

  TH1D *Ddiff = Hafb9_vTcPL_IFI2;
  Ddiff->SetTitle(0); Ddiff->SetStats(0);

  Ddiff->SetMaximum( 0.045); Ddiff->SetMinimum(-0.015);
  Ddiff->DrawCopy("h");

  CaptT->SetTextColor(kBlack);
  CaptT->DrawLatex(0.01, 0.95, "KKMC: A^{IFI}_{FB}(v_{max}) = A^{IFIon}_{FB}(v_{max}) - A^{IFIoff}_{FB}(v_{max}) ");
  double ycapt =0.80;
  PlotSame2(Hafb9_vTcPL_IFI2, ycapt, kBlack,   0.10, "(a2)", "O(#alpha^{2}) #sqrt{s_{+}}=94.3GeV");
  PlotSame2(Hafb9_vTcPL_IFI1, ycapt, kGold,    0.14, "(a1)", "O(#alpha^{1}) #sqrt{s_{+}}=94.3GeV");
  PlotSame2(Hafb9_vTcPL_IFI0, ycapt, kPine,    0.18, "(a0)", "O(#alpha^{0}) #sqrt{s_{+}}=94.3GeV");
  //
  ycapt += -0.04;
  PlotSame2(Hafb8_vTcPL_IFI2, ycapt, kBlack,   0.10, "(b2)", "O(#alpha^{2}) #sqrt{s_{-}}=87.9GeV");
  PlotSame2(Hafb8_vTcPL_IFI1, ycapt, kGold,    0.14, "(b1)", "O(#alpha^{1}) #sqrt{s_{-}}=87.9GeV");
  PlotSame2(Hafb8_vTcPL_IFI0, ycapt, kPine,    0.18, "(b0)", "O(#alpha^{0}) #sqrt{s_{-}}=87.9GeV");

  hZero7->DrawCopy("hsame");

  /////////////////////////////
  cAfbIFI_KKmc2->cd(2);

  Ddiff = Hafb9m8_vTcPL_IFI2;
  Ddiff->SetTitle(0); Ddiff->SetStats(0);
  Ddiff->SetMaximum( 0.005); Ddiff->SetMinimum(-0.012);
  Ddiff->DrawCopy("h");

  CaptT->DrawLatex(0.01, 0.95, "KKMC: A^{+-}_{FB}(v_{max}) = A^{IFI}_{FB}(v_{max},s_{+}) - A^{IFI}_{FB}(v_{max},s_{-}) ");
  ycapt =0.80;
  PlotSame2(Hafb9m8_vTcPL_IFI2, ycapt, kBlack,   0.10, "(2)", "O(#alpha^{2}), a2-b2 ");
  PlotSame2(Hafb9m8_vTcPL_IFI1, ycapt, kGold,    0.10, "(1)", "O(#alpha^{1}), a1-b1 ");
  PlotSame2(Hafb9m8_vTcPL_IFI0, ycapt, kPine,    0.10, "(0)", "O(#alpha^{0}), a0-b0 ");

  hZero7->DrawCopy("hsame");

  //
  }// AfbIFI_KKmc2


///////////////////////////////////////////////////////////////////////////////////
void AfbIFI_KKmc4()
{
//------------------------------------------------------------------------
  cout<<" ========================= AfbIFI_KKmc4 =========================== "<<endl;
  TH1D *Hafb9_vTcPL_IFI2      = (TH1D*)DiskFileB.Get("Hafb9_vTcPL_IFI2");
  TH1D *Hafb9_vTcPL_IFI1      = (TH1D*)DiskFileB.Get("Hafb9_vTcPL_IFI1");
  TH1D *Hafb9_vTcPL_IFI0      = (TH1D*)DiskFileB.Get("Hafb9_vTcPL_IFI0");
  //
  TH1D *Hafb8_vTcPL_IFI2      = (TH1D*)DiskFileB.Get("Hafb8_vTcPL_IFI2");
  TH1D *Hafb8_vTcPL_IFI1      = (TH1D*)DiskFileB.Get("Hafb8_vTcPL_IFI1");
  TH1D *Hafb8_vTcPL_IFI0      = (TH1D*)DiskFileB.Get("Hafb8_vTcPL_IFI0");
  //
  TH1D *Hafb95_IFI2m1 = HstDiff("Hafb95_IFI2m1",  Hafb9_vTcPL_IFI2,  Hafb9_vTcPL_IFI1,  kBlack);
  TH1D *Hafb95_IFI1m0 = HstDiff("Hafb95_IFI2m1",  Hafb9_vTcPL_IFI1,  Hafb9_vTcPL_IFI0,  kBlack);
  TH1D *Hafb95_IFI2m0 = HstDiff("Hafb95_IFI2m1",  Hafb9_vTcPL_IFI2,  Hafb9_vTcPL_IFI0,  kBlack);
  //
  TH1D *Hafb88_IFI2m1 = HstDiff("Hafb88_IFI2m1",  Hafb8_vTcPL_IFI2,  Hafb8_vTcPL_IFI1,  kBlack);
  TH1D *Hafb88_IFI1m0 = HstDiff("Hafb88_IFI1m0",  Hafb8_vTcPL_IFI1,  Hafb8_vTcPL_IFI0,  kBlack);
  TH1D *Hafb88_IFI2m0 = HstDiff("Hafb88_IFI2m0",  Hafb8_vTcPL_IFI2,  Hafb8_vTcPL_IFI0,  kBlack);

  TH1D *hZero7      = (TH1D*)DiskFileB.Get("hZero7");

  //////////////////////////////////////////////
  TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
  CaptT->SetTextSize(0.04);
  ////////////////////////////////////////////////////////////////////////////////
  TCanvas *cAfbIFI_KKmc4 = new TCanvas("cAfbIFI_KKmc4","cAfbIFI_KKmc4", gXcanv,  gYcanv,   1200,  600);
  //                            Name    Title            xoff,yoff, WidPix,HeiPix
  ////////////////////////////////////////////////////////////////////////////////
  gXcanv += 50; gYcanv += 50;
  cAfbIFI_KKmc4->SetFillColor(10);
  cAfbIFI_KKmc4->Divide( 2,  1);
//*****************************************************************************
  cAfbIFI_KKmc4->cd(1);

  TH1D *Ddiff = Hafb9_vTcPL_IFI2;
  Ddiff->SetTitle(0); Ddiff->SetStats(0);

  Ddiff->SetMaximum( 0.045); Ddiff->SetMinimum(-0.015);
  Ddiff->DrawCopy("h");

  CaptT->SetTextColor(kBlack);
  CaptT->DrawLatex(0.01, 0.95, "KKMC:  A^{IFI}_{FB}(v_{max}) = A^{IFIon}_{FB}(v_{max}) - A^{IFIoff}_{FB}(v_{max}) ");
  double ycapt =0.80;
  PlotSame2(Hafb9_vTcPL_IFI2, ycapt, kBlack,   0.10, "(a2)", "O(#alpha^{2}) #sqrt{s_{+}}=94.3GeV");
  PlotSame2(Hafb9_vTcPL_IFI1, ycapt, kGold,    0.14, "(a1)", "O(#alpha^{1}) #sqrt{s_{+}}=94.3GeV");
  PlotSame2(Hafb9_vTcPL_IFI0, ycapt, kPine,    0.18, "(a0)", "O(#alpha^{0}) #sqrt{s_{+}}=94.3GeV");
  //
  ycapt += -0.04;
  PlotSame2(Hafb8_vTcPL_IFI2, ycapt, kBlack,   0.10, "(b2)", "O(#alpha^{2}) #sqrt{s_{-}}=87.9GeV");
  PlotSame2(Hafb8_vTcPL_IFI1, ycapt, kGold,    0.14, "(b1)", "O(#alpha^{1}) #sqrt{s_{-}}=87.9GeV");
  PlotSame2(Hafb8_vTcPL_IFI0, ycapt, kPine,    0.18, "(b0)", "O(#alpha^{0}) #sqrt{s_{-}}=87.9GeV");

  hZero7->DrawCopy("hsame");
//*****************************************************************************
  cAfbIFI_KKmc4->cd(2);

  Ddiff = Hafb95_IFI2m1;
  Ddiff->SetTitle(0); Ddiff->SetStats(0);
  Ddiff->SetMaximum( 0.002); Ddiff->SetMinimum(-0.003);
//  Ddiff->SetMaximum( 0.008); Ddiff->SetMinimum(-0.008);
  Ddiff->DrawCopy("h");

  CaptT->DrawLatex(0.01, 0.95, "KKMC: A^{IFI}_{FB}(v_{max})");
  ycapt =0.50;
  PlotSame2(Hafb95_IFI2m1, ycapt, kBlack,   0.07, "(a+)", "O(#alpha^{2})-O(#alpha^{1}), #sqrt{s_{+}}=94.3GeV ");
  PlotSame2(Hafb95_IFI1m0, ycapt, kGold ,   0.09, "(b+)", "O(#alpha^{1})-O(#alpha^{0}), #sqrt{s_{+}}=94.3GeV ");
  //PlotSame2(Hafb95_IFI2m0, ycapt, kMagenta,    0.11, "(c+)", "O(#alpha^{2})-O(#alpha^{0}) ");
  //
  PlotSame2(Hafb88_IFI2m1, ycapt, kBlue,    0.15, "(a-)", "O(#alpha^{2})-O(#alpha^{1}), #sqrt{s_{-}}=87.9GeV");
  PlotSame2(Hafb88_IFI1m0, ycapt, kPine,    0.17, "(b-)", "O(#alpha^{1})-O(#alpha^{0}), #sqrt{s_{-}}=87.9GeV");
  //PlotSame2(Hafb88_IFI2m0, ycapt, kBrune,   0.19, "(c-)", "O(#alpha^{2})-O(#alpha^{0}) ");

  //
}// AfbIFI_KKmc4

///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  HistNormalize();     // Renormalization of MC histograms
  //KKsemMakeHisto();    // prepare histos for plotting
  ReMakeMChisto();     // reprocessing MC histos
  //========== PLOTTING ==========
  // Template empty canvas  with 2 figures
  //FigTempl();
  //AfbIFIvA1();
  //AfbIFIvA2();
  //
  AfbIFIvT1();
  AfbIFIvT2();
  //
  AfbIFI_Foam();
  //
  AfbDifPat();
  //
  //AfbIFI_KKmc2();
  //AfbIFI_KKmc4();
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileA95.ls();
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();
  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}
