// make PlotBorn-run
//////////////////////////////////////////////////////////
//line-shape i Afb, with and without box corrections
//    w funkcji sqrt(s), zakres 70 - 150 GeV
//    fixed costheta = 0.0, 0.3, 0.6, 0.9
//    procesy: ee-> mumu, uubar, ddbar
//////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
using namespace std;

#include <math.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TH2.h"
#include "TLatex.h"
#include "TApplication.h"
#include "TObjString.h"
#include "TFile.h"
//
#include "TLine.h"

#include "HisNorm.h"
#include "KKplot.h"

KKplot LibSem("kkplot");

//
TFile DiskFileB("RhoSemi.root","RECREATE","histograms");
FILE *DFile;

///////////////////////////////////////////////////////////////////////////////////
//              GLOBAL stuff
///////////////////////////////////////////////////////////////////////////////////
double m_xpar[10001];      // complete input of KKMC run
double gCMSene, gNevTot; // from KKMC run
char   gTextEne[100], gTextNev[100], gTextNev2[100];
int    kGold=kOrange-3, kBrune=46, kPine=kGreen+3;
//
float  gXcanv = 20, gYcanv = 20, gDcanv = 30;

#define SP21 setw(21)<<setprecision(13)
#define SP15 setw(15)<<setprecision(9)
#define SP10 setw(10)<<setprecision(5)

///////////////////////////////////////////////////////////////////////////////////

void Vdef(double v[4], const double v1, const double v2, const double v3, const double v4)
  { // define a 4-vector (avoids initialization warnings)
       v[0] = v1; v[1] = v2; v[2] = v3; v[3] = v4;
  }


void InitBorn(){
cout<<"==========================================================="<<endl;
cout<<"================ InitBorn BEGIN ============================"<<endl;

//////////////////////////////////////////////////////////////
//   Initialize MC generator and analysis programs          //
//////////////////////////////////////////////////////////////
const int jmax =10000;
// in the input below only AMH=125e0, AMTOP=173e0 are actualized
LibSem.ReaData("./KK2f_defaults_ERW",    jmax, m_xpar);  // AMH AMTOP actualized
// standard old input file
//LibSem.ReaData("../../.KK2f_defaults",     jmax, m_xpar);  // numbering as in input!!!
// User input data with AMH AMTOP aclualized
//LibSem.ReaData("../workKKMC/workKKMC_189GeV.input", -jmax, m_xpar);  // jmax<0 means no-zeroing
LibSem.ReaData("../workKKMC/PlotBorn_189GeV.input", -jmax, m_xpar);  // jmax<0 means no-zeroing
//*****
double ypar[jmax];
for(int j=0;j<jmax;j++) ypar[j]=m_xpar[j+1];    // ypar has c++ numbering
// KKMC and KKsem initialization
const char *output_file = "./kksem.output";
long stl2 = strlen(output_file);
long mout =16;
//kk2f_fort_open_(mout,output_file,stl2);
//kk2f_initialize_(ypar);
//kksem_initialize_(ypar);
LibSem.Initialize(m_xpar);

//************************************
DFile = fopen("TableBorn.txt","w");
//************************************

cout<<"================ InitBorn END   ============================"<<endl;

}//InitBorn


void TabBorn(){
cout<<"==========================================================="<<endl;
cout<<"================ TabBorn BEGIN ============================"<<endl;
/*
Zalaczam wersje ktora uzywam do produkowania lookup tables
dla EW poprawek:
  .KK2f_defaults
   BornV.h
Interesuja mnie bechmarki dla
   "doubly-deconvoluted" = no ISR/FSR
line-shape i Afb, with and without box corrections
    w funkcji sqrt(s), zakres 70 - 150 GeV
    fixed costheta = 0.0, 0.3, 0.6, 0.9
    procesy: ee-> mumu, uubar, ddbar
papier gdzie sa przyklady takich benchmarkow:
https://arxiv.org/pdf/hep-ph/9902452.pdf
pozdrawiam
Ela
*/
//  ************************************************************************
double Emin =  70;  //GeV
double Emax = 150;  //GeV
int    nPt  =  40;  // No of. points
//
double svar2, sigma,  CosTheta;
double dSig0, dSig3, dSig6, dSig9;
int m_KFini = 11;
//----------------------
int m_KFf   = 13;
fprintf(DFile," e+e- --> mu+ mu-,  results from BornV_Dizet \n");
//----------------------
//int m_KFf   = 1;
//fprintf(DFile," e+e- --> d dbar,   results from BornV_Dizet \n");
//int m_KFf   = 2;
//fprintf(DFile," e+e- --> u ubar,   results from BornV_Dizet \n");
//
fprintf(DFile," d(sigma)/d(cos_theta) [nb], cos_theta = 0.0, 0.3, 0.6, 0.9 \n");
double Ene;
//
for( int i=0; i<=nPt; i++ ){
   Ene = Emin +i*((Emax-Emin)/nPt);
   svar2 = Ene*Ene;
   CosTheta=0; bornv_interpogsw_(m_KFf,svar2, CosTheta);
   dSig0 = bornv_dizet_( 1, m_KFini, m_KFf, svar2, CosTheta, 0.0, 0.0, 0.0, 0.0);
   //
   CosTheta=0.3; bornv_interpogsw_(m_KFf,svar2, CosTheta);
   dSig3 = bornv_dizet_( 1, m_KFini, m_KFf, svar2, CosTheta, 0.0, 0.0, 0.0, 0.0);
   //
   CosTheta=0.6; bornv_interpogsw_(m_KFf,svar2, CosTheta);
   dSig6 = bornv_dizet_( 1, m_KFini, m_KFf, svar2, CosTheta, 0.0, 0.0, 0.0, 0.0);
   //
   CosTheta=0.9; bornv_interpogsw_(m_KFf,svar2, CosTheta);
   dSig9 = bornv_dizet_( 1, m_KFini, m_KFf, svar2, CosTheta, 0.0, 0.0, 0.0, 0.0);
   fprintf(DFile,"Ene= %10.5f  dSig_dCos= %12.7f  %12.7f  %12.7f  %12.7f  \n", Ene, dSig0, dSig3, dSig6, dSig9);
//   cout<<"Ene="<<Ene; cout<< "  Sig_dCos= "<< dSig0 <<"  "<< dSig3 <<"  "<< dSig6 <<"  "<< dSig9<<endl;
}//i

cout<<"==========================================================="<<endl;

double KeyFob;
KeyFob=   10; // BornV_Dizet, with EW and without integration ???
KeyFob=  -11; // BornV_Simple, for KeyLib=0, NO EW, NO integration OK
KeyFob=    0; // With EW (BornV_Dizet) With integration OK!
kksem_setkeyfob_( KeyFob );
fprintf(DFile,"****************************************************************** \n");
fprintf(DFile," Born with EW (BornV_Dizet) With Gauss integration over cos(theta) \n");
double xBorn;
for( int i=0; i<=nPt; i++ ){
   Ene = Emin +i*((Emax-Emin)/nPt);
   svar2 = Ene*Ene;
   kksem_makeborn_( svar2, xBorn);
   //cout<< "*** Ene= "<<Ene<< "  xBorn [nb]= "<<xBorn<<endl;
   fprintf(DFile,"Ene= %10.5f  sigma [nb] = %12.7f   \n", Ene, xBorn);
}

//************************************
  fclose(DFile);
//************************************

cout<<"================ TabBorn END   ============================"<<endl;
}//TabBorn


void FigBorn3(){
cout<<"==========================================================="<<endl;
cout<<"================ FigBorn3 BEGIN ==========================="<<endl;
double Emin = 85;
double Emax = 95;
int    nPt  = 200;
double Ene;
TH1D *hst_sigma3 = new TH1D("hst_sigma3" ,  "sigma(CMSene) [nb]",    nPt, Emin, Emax);
hst_sigma3->Sumw2();

long KeyFob;
KeyFob=   10; // BornV_Dizet, with EW and without integration ???
KeyFob=  -11; // BornV_Simple, for KeyLib=0, NO EW, NO integration OK
KeyFob=    0; // BornV_Dizet, with EW and with integration
kksem_setkeyfob_( KeyFob );
double xBorn, svar2;
for(int ix=1; ix <= nPt; ix++){
   Ene = Emin +(ix-1)*((Emax-Emin)/nPt);
   svar2 = Ene*Ene;
   kksem_makeborn_( svar2, xBorn);
   hst_sigma3->SetBinContent(  ix, xBorn );
   hst_sigma3->SetBinError(    ix, 0.0 );
   cout<< "Ene= "<<Ene<< "  xBorn [nb]= "<<xBorn<<endl;
}//ix
//------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////
TCanvas *cFigBorn3 = new TCanvas("cFigBorn3","cFigBorn3: sigma(CMSene) ", gXcanv, gYcanv,    500,  500);
//                                      Name    Title                     xoff,     yoff, WidPix,HeiPix
gXcanv += 25, gYcanv += 25;
//
cFigBorn3->SetFillColor(10);
TH1D *Hst1 = hst_sigma3;
cFigBorn3->cd();
Hst1->SetStats(0);
//Hst1->SetTitle(0);
Hst1->SetMinimum(0);
Hst1->SetLineColor(kRed);
Hst1->DrawCopy("h");
cout<<"================ FigBorn3 END   ==========================="<<endl;
}//FigBorn3



void FigBorn1(){
cout<<"==========================================================="<<endl;
cout<<"================ FigBorn1 BEGIN ==========================="<<endl;
int    nPt  = 160;
double CosTheta = 0; // for some reason 0 no allowed!!!
CosTheta = 0.01;
double Emin =  70;
double Emax = 150;
//Emin=30;
//
//Emin=0.250;
//Emax = 10;
//
double CMSene;
//
double m_p1[4];         //!
double m_p2[4];         //!
double m_p3[4];         //!
double m_p4[4];         //!
//
TH1D *hst_sigEEX = new TH1D("hst_sigEEX" ,  "dSigma(CMSene,costhe)",    nPt, Emin, Emax);
hst_sigEEX->Sumw2();
TH1D *hst_sigGPS = new TH1D("hst_sigGPS" ,  "dSigma(CMSene,costhe)",    nPt, Emin, Emax);
hst_sigGPS->Sumw2();

long KeyFob;
KeyFob=   10; // BornV_Dizet, with EW and without integration ???
KeyFob=  -11; // BornV_Simple, for KeyLib=0, NO EW, NO integration OK
KeyFob=    0; // BornV_Dizet, with EW and with integration
kksem_setkeyfob_( KeyFob );
int KFf =13, KFini=11;
double xBorn, svar2;
for(int ix=1; ix <= nPt; ix++){
   CMSene = Emin +(ix-1)*((Emax-Emin)/nPt);
   svar2 = CMSene*CMSene;
   //
   bornv_interpogsw_(KFf,svar2, CosTheta);
   double dSig_EEX = bornv_dizet_( 1, KFini, KFf, svar2, CosTheta, 0.0, 0.0, 0.0, 0.0);
   //
   // =============== Sigm/dOmega from spin amplitudes ===============
   // Effective 4-momenta, KKMC convention: p={px,py,pz,E)
   double Ene   = CMSene/2;
   double   m_beam  = 0.510999e-3;
   double	m_fin   = 0.1056583;
   double Pmb  = sqrt( (Ene-m_beam)*(Ene+m_beam) ); // modulus
   Vdef(m_p1, 0, 0 , Pmb, Ene);  // beam
   Vdef(m_p2, 0, 0 ,-Pmb, Ene);  // beam
   double Pmf  =sqrt( (Ene-m_fin)*(Ene+m_fin) ); // modulus
   Vdef(m_p3, Pmf*sqrt(1-sqr(CosTheta)), 0 , Pmf*CosTheta,  Ene); // final
   Vdef(m_p4,-Pmf*sqrt(1-sqr(CosTheta)), 0 ,-Pmf*CosTheta,  Ene); // final
   double PX[4] = {0, 0, 0, 2*Ene};
   //***** pure Born of CEEX
   double dSig_GPS;
   gps_bornf_(KFini, KFf ,PX, CosTheta, m_p1,m_beam, m_p2, -m_beam,
                                        m_p3,m_fin,  m_p4, -m_fin,   dSig_GPS);
   /////////////////////////////////////////////////////////////////
   cout<< "CMSene= "<<CMSene<< "  dSig_EEX= "<<dSig_EEX;
   cout<<                      "  dSig_GPS= "<<dSig_GPS<<"  ratio="<< dSig_GPS/dSig_EEX  <<endl;

   hst_sigEEX->SetBinContent(  ix, dSig_EEX );
   hst_sigEEX->SetBinError(    ix, 0.0 );
   hst_sigGPS->SetBinContent(  ix, dSig_GPS );
   hst_sigGPS->SetBinError(    ix, 0.0 );
/*
   double corQED;
   bornv_getqedcor_( corQED ); // after calling bornv_interpogsw_ !
   double AlfRunInv = 137.04/corQED;
   cout<< "CMSene= "<<CMSene<< "  corQED= "<<corQED << "   AlfRunInv= "<< AlfRunInv <<  endl;
*/
}//ix

//------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////
TCanvas *cFigBorn1 = new TCanvas("cFigBorn1","cFigBorn1: sigma(CMSene) ", gXcanv, gYcanv, 900,  600);
//                                      Name    Title                  xoff,yoff, WidPix,HeiPix
gXcanv += 25, gYcanv += 25;
//
cFigBorn1->SetFillColor(10);
TPad*    UpperPad = new TPad("upperPad", "upperPad",
			       .005, .3525, .995, .995);
TPad*    LowerPad = new TPad("lowerPad", "lowerPad",
			       .005, .005, .995, .3575);
UpperPad->Draw();
LowerPad->Draw();
/////////////////////////////////////////////////////
UpperPad->cd();
TH1D *Hst1 = hst_sigGPS;
Hst1->SetStats(0);
//Hst1->SetTitle(0);
Hst1->SetMinimum(0);
Hst1->SetLineColor(kRed);       //GPS
Hst1->DrawCopy("h");
hst_sigEEX->SetLineColor(kBlue);
hst_sigEEX->DrawCopy("hsame");  //EEX

TLatex *CaptT = new TLatex(); CaptT->SetNDC(); // !!!
char TextCos[100]; sprintf(TextCos,"cos#theta =%6.2f ", CosTheta);
CaptT->DrawLatex(0.60,0.75,TextCos);

////////////////////////////////////////////////////
LowerPad->cd();

TH1D *dSigRatio =(TH1D*)hst_sigGPS->Clone("dSigRatio");
dSigRatio->Divide(dSigRatio, hst_sigEEX); // sigma(xmax), FoamEEX/KKsem IFIoff
dSigRatio->SetLineColor(kRed);                 // FoamEEX2/KKsem IFIoff

dSigRatio->GetYaxis()->SetLabelSize(0.08); // default is 0.04
dSigRatio->GetXaxis()->SetLabelSize(0.08); // default is 0.04

dSigRatio->SetMaximum(1+0.00003);
dSigRatio->SetMinimum(1-0.00003);
//dSigRatio->SetMaximum(1+0.0003);
//dSigRatio->SetMinimum(1-0.0003);
//dSigRatio->SetMaximum(1+0.3);
//dSigRatio->SetMinimum(1-0.3);
dSigRatio->SetStats(0);
dSigRatio->SetTitle(0);
dSigRatio->DrawCopy("h");


TH1D *hOne = (TH1D*)dSigRatio->Clone("hOne");  // unity line
for(int i=1; i <= hOne->GetNbinsX() ; i++) { hOne->SetBinContent(i, 1); hOne->SetBinError(i, 0);}
hOne->SetLineColor(kRed);
hOne->DrawCopy("hsame");


cout<<"================ FigBorn1 END   ==========================="<<endl;
}//FigBorn1


////////////////////////////////////////////////////////////////////////////////
// reproducing benchmarks of Elzbieta Richter-Was
////////////////////////////////////////////////////////////////////////////////
void FigAlfE1(){
cout<<"==========================================================="<<endl;
cout<<"================ FigAlfE1 BEGIN ==========================="<<endl;
int    nPt  = 160;
double CosTheta = 0; // for some reason 0 no allowed!!!
CosTheta = 0.01;
double Emin =  70;
double Emax = 150;
double MZ = 91.187e0;
//
TH1D *hst_AlfRun = new TH1D("hst_AlfRun" ,  " Alfa(sqrt(s))",    nPt, Emin, Emax);
hst_AlfRun->Sumw2();
TH1D *hst_AlfRat = new TH1D("hst_AlfRat" , "Alfa(s)/Alfa(0)=1/(2-GSW6)",    nPt, Emin, Emax);
hst_AlfRat->Sumw2();

int KFf =13, KFini=11;
double CMSene, xBorn, svar2, GSW6, AlfRatio, alfinv0;
for(int ix=1; ix <= nPt; ix++){
   CMSene = Emin +(ix-1)*((Emax-Emin)/nPt);
   svar2 = CMSene*CMSene;
   bornv_interpogsw_(KFf,svar2, CosTheta);
   bornv_getqedcor_( AlfRatio ); // must be after calling bornv_interpogsw_ !
   alfinv0 = 137.0359895e0;
   double AlfRun    = AlfRatio/alfinv0;
   double AlfRunInv = alfinv0/AlfRatio;
   cout<< "CMSene= "<<CMSene<< " AlfRatio= "<< AlfRatio << "   AlfRunInv= "<< AlfRunInv <<  endl;
   hst_AlfRat->SetBinContent(  ix, AlfRatio );
   hst_AlfRat->SetBinError(    ix, 0.0 );
   hst_AlfRun->SetBinContent(  ix, AlfRun );
   hst_AlfRun->SetBinError(    ix, 0.0 );
}//ix
//!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  TCanvas *cAlfE1 = new TCanvas("cAlfE1","cAlfE1", gXcanv,  gYcanv,  1200, 600);
  gXcanv += gDcanv; gYcanv += gDcanv;
//
  cAlfE1->SetFillColor(10);
  cAlfE1->Draw();
  cAlfE1->Divide(2, 0);
/////////////////////////////////////////////
  cAlfE1->cd(1);
  hst_AlfRat->SetStats(0);
  hst_AlfRat->SetMinimum(1.00); hst_AlfRat->SetMaximum(1.14);
  hst_AlfRat->DrawCopy("h");
  /////////////////////////////////////////////
  cAlfE1->cd(2);
  hst_AlfRun->SetStats(0);
  hst_AlfRun->SetMinimum(0.0077); hst_AlfRun->SetMaximum(0.0079);
  hst_AlfRun->DrawCopy("h");
  cAlfE1->cd();
//
  cAlfE1->SaveAs("cAlfE1.pdf");
cout<<"================ FigAlfE1 END   ==========================="<<endl;
}//FigAlfE1


////////////////////////////////////////////////////////////////////////////////
void FigAlfE2(){
cout<<"==========================================================="<<endl;
cout<<"================ FigAlfE2 BEGIN ==========================="<<endl;
int    nPt  = 100;
double CosTheta = 0.001; // for some reason 0 no allowed!!!
CosTheta = 0.01;
double MZ = 91.187e0;
double E1 =  87.90e0, E2 =  94.30e0;
double Emin   =  MZ-4.0,  Emax  =  MZ+4.0;
//
TH1D *hst_AlfRun = new TH1D("hst_AlfRun" ,  " Alfa(sqrt(s))",   nPt, Emin, Emax);
hst_AlfRun->Sumw2();
TH1D *hst_AlfInv = new TH1D("hst_AlfInv" , " 1/Alfa(sqrt(s)) ", nPt, Emin, Emax);
hst_AlfInv->Sumw2();
//
int KFf =13, KFini=11;
double CMSene, xBorn, svar2, GSW6, AlfRatio, AlfEst,xi;
double alfinv0 = 137.0359895e0;
double AlfRun;
// at Z peak
bornv_interpogsw_(KFf,MZ*MZ, CosTheta);
bornv_getqedcor_( AlfRatio ); // must be after calling bornv_interpogsw_ !
double AlfRunMZ    = AlfRatio/alfinv0;
// at E1
bornv_interpogsw_(KFf,E1*E1, CosTheta);
bornv_getqedcor_( AlfRatio ); // must be after calling bornv_interpogsw_ !
double AlfRun1    = AlfRatio/alfinv0;
// at E2
bornv_interpogsw_(KFf,E2*E2, CosTheta);
bornv_getqedcor_( AlfRatio ); // must be after calling bornv_interpogsw_ !
double AlfRun2    = AlfRatio/alfinv0;
double xres[1];
//
for(int ix=1; ix <= nPt; ix++){
   CMSene = Emin +(ix-0.5)*((Emax-Emin)/nPt);
   svar2 = CMSene*CMSene;
   bornv_interpogsw_(KFf,svar2, CosTheta);
   bornv_getqedcor_( AlfRatio ); // must be after calling bornv_interpogsw_ !
   AlfRun    = AlfRatio/alfinv0;
//[[[
//   double sig0,afb0; int KeyDist=0;
//   kksem_born_calc_(KFini,KFf,CMSene, AlfRun,xres);
//   sig0 = xres[1]; afb0 = xres[2];
//   cout<< ">>>>>>>>>>>>> From kksem_born_calc: sig0,afb0= "<<  sig0<<" ,   " << afb0 << endl;
//]]]
   //
   cout<< "CMSene= "<<CMSene<< " AlfRatio= "<< AlfRatio << "   AlfRunInv= "<< 1/AlfRun <<  endl;
   //
   xi = log(E1*E2/CMSene/CMSene)/log(E1/E2);       // interpolation based on RGE
   AlfEst = 2/( (1-xi)/AlfRun1 +(1+xi)/AlfRun2 );
   cout<< "Log. Estim: xi= "<<  xi << "  1/AlfEst ="<< 1/AlfEst <<"  Ratio ="<< AlfEst/AlfRun <<endl;
   xi = (CMSene*CMSene-E1*E1)/(E2*E2-E1*E1);       // interpolation in s=E*E
   xi = (CMSene-E1)/(E2-E1);                       // interpolation in E
   xi = (log(CMSene)-log(E1))/(log(E2)-(log(E1))); // interpolation in ln(E)
   AlfEst = (1-xi)*AlfRun1 +xi*AlfRun2;
   cout<< "Lin. Estim: xi= "<<  xi << "  1/AlfEst ="<< 1/AlfEst <<"  AlfEst/AlfRun-1 ="<< AlfEst/AlfRun-1 <<endl;
   //
   hst_AlfInv->SetBinContent(  ix, 1/AlfRun ); hst_AlfInv->SetBinError(    ix, 0.0 );
//   hst_AlfRun->SetBinContent(  ix, AlfRun );   hst_AlfRun->SetBinError(    ix, 0.0 );
   hst_AlfRun->SetBinContent(  ix, AlfEst/AlfRun-1 );   hst_AlfRun->SetBinError(    ix, 0.0 );
//   hst_AlfRun->SetBinContent(  ix, afb0 );   hst_AlfRun->SetBinError(    ix, 0.0 );
}//ix
double betaLN = (AlfRun2-AlfRun1)/AlfRunMZ/(2*log(E2/E1));
cout<< "Ln(s) interpol. betaLN= "<<  betaLN << endl;
double betaE  = (AlfRun2-AlfRun1)/AlfRunMZ/( (E2-E1)/MZ );
cout<< "interpol.in E,   betaE= "<<  betaE << endl;
double betaEE = (AlfRun2-AlfRun1)/AlfRunMZ/( (E2*E2-E1*E1)/MZ/MZ );
cout<< "interpol.in E*E, betaEE="<<  betaEE << endl;
double beta   = (AlfRun2-AlfRun1)/(2*log(E2/E1))/AlfRunMZ/AlfRunMZ;
cout<< "beta from alpha(s)= "<<  beta << endl;
double beta0  = ( 3*3*(1./9.) +3*2*(4./9.) +3)/(3*3.141594);
cout<< "O(alf1) RGE beta0 = "<<  beta0 << endl;
cout<< "zeta of Patrick:  log(E1*E2/MZ/MZ)/log(E1/E2) = "<<  log(E1*E2/MZ/MZ)/log(E1/E2) <<endl;
cout<< "ratio of factors  (E2*E2-MZ*MZ)/(MZ*MZ-E1*E1) = "<<  (E2*E2-MZ*MZ)/(MZ*MZ-E1*E1)<< endl;

char  TextAlfInv[100];
sprintf(TextAlfInv,"#alpha^{-1}_{QED}(M_{Z}) =%4.4f", 1/AlfRunMZ);

//////////////////////////////////////////////
  TLatex *CaptNDC = new TLatex();
  CaptNDC->SetNDC(); // !!!
  CaptNDC->SetTextSize(0.035);
//////////////////////////////////////////////
  TLatex *CaptT = new TLatex();
  CaptT->SetTextSize(0.035);
  //
  double alfmin=128.80,  alfmax = 129.00;
  TLine *lineMZ = new TLine(MZ, alfmin, MZ, alfmax);
  lineMZ->SetLineStyle(2);
  TLine *line88 = new TLine(E1, alfmin, E1, alfmax);
  line88->SetLineStyle(2);
  TLine *line94 = new TLine(E2, alfmin, E2, alfmax);
  line94->SetLineStyle(2);

//!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  TCanvas *cAlfE2 = new TCanvas("cAlfE2","cAlfE2", gXcanv,  gYcanv,  1200, 600);
  gXcanv += gDcanv; gYcanv += gDcanv;
//
  cAlfE2->SetFillColor(10);
  cAlfE2->Draw();
  cAlfE2->Divide(2, 0);
/////////////////////////////////////////////
  cAlfE2->cd(1);
  hst_AlfInv->SetStats(0);
  hst_AlfInv->SetTitle(0);
  hst_AlfInv->SetMinimum(alfmin);
  hst_AlfInv->SetMaximum(alfmax);
  hst_AlfInv->SetLineColor(kBlack);
  hst_AlfInv->DrawCopy("h");

  hst_AlfInv->Scale(1.0001);
  hst_AlfInv->SetLineColor(kRed);
  hst_AlfInv->DrawCopy("hsame");

  hst_AlfInv->Scale(1/1.0001/1.0001);
  hst_AlfInv->DrawCopy("hsame");

  lineMZ->Draw(); line88->Draw(); line94->Draw();

  CaptT->DrawLatex(MZ+0.1,    alfmax-0.01,"M_{Z}");
  CaptT->DrawLatex(E2 -0.6,alfmax-0.01,"94.3 GeV");
  CaptT->DrawLatex(E1 -0.5,alfmax-0.01,"87.9 GeV");

  CaptNDC->DrawLatex(0.03,0.95,"#alpha^{-1}_{QED}(s^{1/2})");
  CaptNDC->DrawLatex(0.50,0.02,"s^{1/2} [GeV]");

  CaptNDC->DrawLatex(0.20,0.20,TextAlfInv);

  /////////////////////////////////////////////
  cAlfE2->cd(2);
  hst_AlfRun->SetStats(0);
  hst_AlfRun->SetTitle(0);
  //hst_AlfRun->SetMinimum(0.0077); hst_AlfRun->SetMaximum(0.0079);
  hst_AlfRun->DrawCopy("h");

  CaptNDC->DrawLatex(0.03,0.95," Testing linear approximation ");
//
  cAlfE2->cd();
//
  cAlfE2->SaveAs("cAlfE2.pdf");
cout<<"================ FigAlfE2 END   ==========================="<<endl;
}//FigAlfE2



////////////////////////////////////////////////////////////////////////////////
void FigAlfE3(){
cout<<"   "<<endl;
cout<<"==========================================================="<<endl;
cout<<"================ FigAlfE3 BEGIN ==========================="<<endl;

double CosTheta = 0.001;    // for some reason 0 no allowed!!!
//
double MZ = 91.187e0;
double GammZ = 2.50072032;    // from KKdefaults!
double E1 =  87.90e0, E2 =  94.30e0;
double Emin   =  MZ-4.0,  Emax  =  MZ+4.0;
double GFermi = 1.166397e-5;  // one digit less in KKdefaults!
double Pi =3.1415926535897932;
//
//
int KFf =13, KFini=11;
double CMSene, xBorn, svar2, GSW6, AlfRatio, AlfEst,xi;
double alfinv0 = 137.0359895e0;
double AlfRun[3];
double AFBclc[3];
double xres[100];

for(int ix=0; ix <= 2; ix++){
  if( ix == 0) CMSene = MZ;
  if( ix == 1) CMSene = E1;
  if( ix == 2) CMSene = E2;
  svar2 = CMSene*CMSene;
  bornv_interpogsw_(KFf,svar2, CosTheta);
  bornv_getqedcor_( AlfRatio ); // must be after calling bornv_interpogsw_ !
  AlfRun[ix]    = AlfRatio/alfinv0;
  //
  cout<< " ix= "<< ix <<endl;
  cout<< "CMSene= "<<CMSene<< " AlfRatio= "<< AlfRatio << "   AlfRunInv= "<< 1/AlfRun[ix] <<  endl;

  double sig0,afb0;
  kksem_afb_calc_(0,KFini,KFf,CMSene,0.1,sig0);
  kksem_afb_calc_(1,KFini,KFf,CMSene,0.1,afb0);
  cout<< "+++ kksem_afb_calc:  sig0,afb0= "<<  sig0<<" ,   " << afb0 << endl;
  //
  kksem_born_calc_(KFini,KFf,CMSene, AlfRun[ix],xres);
  sig0= xres[0]; afb0=xres[1];
  cout<< "%%% kksem_born_cal:  sig0,afb0= "<<  sig0<<" ,   " << afb0 << endl;
  AFBclc[ix] = afb0;  // from ksem_born_cal
}//ix
double SinW2 =  m_xpar[503];
cout<< " SinW2 = m_xpar[503]=" << SP15 << m_xpar[503] << endl; // swsq=xpar[503]

cout<< "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
cout<< "*****  Checking quality of the interpolation of alfQED(E) in ln(E) " << endl;
xi = (log(MZ)-log(E1))/(log(E2)-(log(E1))); // interpolation in ln(E)
AlfEst = (1-xi)*AlfRun[1] +xi*AlfRun[2];
cout<< "xi= "<<  xi << "   1/AlfEst(MZ) ="<< 1/AlfEst <<"  AlfEst/AlfRun-1 ="<< AlfEst/AlfRun[0]-1 <<endl;
double chi2=log(E2/MZ), chi1=log(MZ/E1);
AlfEst = (chi2*AlfRun[1] +chi1*AlfRun[2])/(chi1+chi2);
cout<<"chi1,2/(1+2)= "<<chi1/(chi1+chi2) <<" "<< chi2/(chi1+chi2)<<"  AlfEst/AlfRun-1 ="<< AlfEst/AlfRun[0]-1 <<endl;

double svar = MZ*MZ;
double ga = -1.0/2.0;
double gv = ga*(1-4*SinW2);
double cG = AlfRun[0];
double cZ = MZ*MZ*GFermi/(2*sqrt(2)*Pi);
double BW = sqr(svar-MZ*MZ)+ sqr(svar*GammZ/MZ);
// coefficiets of (1+cos_theta**2)
double xGG =   cG*cG;
double xGZ = 2*cG*cZ *sqr(gv) *(svar-MZ*MZ)*svar/BW;
double xZZ =   cZ*cZ *sqr(ga*ga+gv*gv)*sqr(svar)/BW;
// coefficiets of 2*cos_theta
double yGZ = 2*cG*cZ*ga*ga    *(svar-MZ*MZ)*svar/BW;
double yZZ = cZ*cZ*4* ga*ga*gv*gv     *sqr(svar)/BW;
// complete AFB
double AFB = 3.0/4.0*(yGZ+yZZ)/(xGG+xGZ+xZZ);
cout<< "*****  variations on AFB " << endl;
cout<<" local AFB(MZ)= "<< AFB <<" from  kksem_afb_calc AFB="<< AFBclc[0] <<endl;
cout<<" PURE Z: AFB = 0.75*4*sqr(gv*ga)/sqr(sqr(gv)+sqr(ga)) = "<< 0.75*4*sqr(gv*ga)/sqr(sqr(gv)+sqr(ga)) <<endl;

//cout<<" 2*gv,2*ga = "<<  2*gv <<" "<< 2*ga <<" SinW2= " << SinW2<<endl;
//cout<<" cZ/4 ="<< cZ/4<<endl;
//cout<<" xGG, xZZ, yZZ= "<< xGG <<"  "<< xZZ <<"  "<< yZZ <<endl;
//cout<<" sqr(svar)/BW= "<< sqr(svar)/BW <<endl;
//cout<<" GammZ= " <<  GammZ <<endl;



kksem_born_calc_(KFini,KFf,MZ, AlfRun[0],xres);

cout<<"================ FigAlfE3 END   ==========================="<<endl;
}//FigAlfE3


///////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //++++++++++++++++++++++++++++++++++++++++
  TApplication theApp("theApp", &argc, argv);
  //++++++++++++++++++++++++++++++++++++++++
  //
  DiskFileB.cd();
  //
  InitBorn();
  TabBorn();
  //FigBorn1();
  //FigBorn3();
  //
  FigAlfE1();
  FigAlfE2();
  FigAlfE3();
  //
  //++++++++++++++++++++++++++++++++++++++++
  DiskFileB.ls();
  DiskFileB.Write();
  DiskFileB.Close();

  //++++++++++++++++++++++++++++++++++++++++
  theApp.Run();
  //++++++++++++++++++++++++++++++++++++++++
}

