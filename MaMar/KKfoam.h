//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class   KKfoam                                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
// KKfoam is multipurpose toolbox for KKMC testing.
//  1. Interfaces (wrappers) to KKMC and KKsem F77 subrograms
//  2. Integrand for Foam in semianalytical xcheck
//  3. A few routines for producing latex table out of histograms
//////////////////////////////////////////////////////////////////////////////

#ifndef KKFOAM_H
#define KKFOAM_H
#include<stdlib.h>
#include<stdio.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

#include "TRandom3.h"
#include "TFoam.h"
#include "TFoamIntegrand.h"

//------------------------------------------------------------
//  wrappers to f77 routines in KKMC and KKsem
extern "C" void kk2f_fort_open_( const int&, const char*, int);
extern "C" void kk2f_fort_close_(const int&);
//      SUBROUTINE KK2f_Initialize(xpar)
extern "C" void kk2f_initialize_(double xpar[]);
//-----------------------
//      SUBROUTINE BornV_SetKF(KFferm)
extern "C" void bornv_setkf_( const int& ); // set SINGLE Final State
//---------------------------------
//      DOUBLE PRECISION  FUNCTION BornV_Sig0nb(CMSene)
extern "C" double bornv_sig0nb_(const double&);
//      SUBROUTINE BornV_MakeGami(CMSene,gamiCR,gami,alfi)
extern "C" void bornv_makegami_(const double&, double&, double&, double&);
//      DOUBLE PRECISION  FUNCTION BornV_Simple(KFi,KFf,svar,costhe)
extern "C" double bornv_simple_( const int&,  const int&, const double&, const double&);
//------------------------------------------------------------
//      SUBROUTINE BornV_InterpoGSW(KFf,svar,CosThe)
extern "C" double bornv_interpogsw_( const int&,  const double&, const double&);
//      DOUBLE PRECISION FUNCTION BornV_Dizet(Mode,KFi,KFf,svar,CosThe,eps1,eps2,ta,tb)
extern "C" double bornv_dizet_(const int&, const int&, const int&,
		const double&, const double&, const double&, const double&, const double&, const double& );
//------------------------------------------------------------
//      SUBROUTINE GPS_BornF(KFi,KFf,PX,CosThe,p1,m1,p2,m2,p3,m3,p4,m4,Xborn)
extern "C" void gps_bornf_(const int&, const int&, double[], const double&,
		    double[], const double&, double[], const double&, double[], const double&, double[], const double&,
		    const double&);
//------------------------------------------------------------
//      SUBROUTINE GPS_BornFoam(Mode,KFi,KFf,CMSene,CosThe,Xborn)
extern "C" void gps_bornfoam_(const int&,   const int&,   const int&,
		                      const double&, const double&, const double&);
//      DOUBLE PRECISION  FUNCTION GPS_MakeRhoFoam(XNorm)
extern "C" double gps_makerhofoam_(const double&);

class KKfoam: public TFoamIntegrand{
// Interface and extensions to KKsem toolbox
 public:
    int       m_jmax;
    double    m_ypar[10000];   // input parameters of KKMC
//
 public:
 	double m_CMSene;
 	double m_Mmin;
 	double m_Mmax;
 	double m_vvmax;
 	//
 	double m_gnanob;
 	double m_pi;
 	double m_ceuler;
 	double m_alfinv;
 	double m_alfpi;
 	//
 	double m_beam;
 	double m_chini;
 	//
 	double m_fin;
 	double m_chfin;
 	int   m_KFini;  // electron
 	int   m_KFf;    // muon
 	//
 	int    m_kDim;
 	int    m_nCells;
 	int    m_nSampl;
 	int    m_KeyISR;
 	int    m_KeyFSR;
//
 	int    m_Mode;   // operation mode for Density
 	double m_del;
//******** MC EVENT ********
 	double m_CosTheta;
 	double m_vv;  // ISR
 	double m_uu;  // FSR
 	double m_r1;  // IFI
 	double m_r2;  // IFI
 	double m_xx;  // total
 	//
 	double m_Mka;
 	//
 	double m_p1[4];
 	double m_p2[4];
 	double m_p3[4];
 	double m_p4[4];
 	//
 	long  m_count;
//
//------ constructors destructors -------
 public:
  KKfoam(){;}
  ~KKfoam(){;}
  KKfoam(const char*);
public:
  // Interfaces to KKsem integration routines using Gauss method
  void Initialize( double ypar[10000]);
  void VVplot( TH1 *, int , char [], int, int );
  void Cplot(  TH1 *, int , char [], int, int, double, double);

  // Foam integrand
  double Fyfs( double );
  double gamISR( double );
  double gamFSR( double );
  double gamIFI( double );
  double Rho_isr(double, double );
  double Rho_fsr(double, double );
  double Rho_ifi(double, double , double );
  void MapIFI1( double, double, double, double &, double &);
  void MapIFI2( double, double, double, double &, double &);
  void MapIFI(  double, double, double, double &, double &);
  Double_t Density(int, Double_t*);
  Double_t Density3(int, Double_t*);
  Double_t Density5(int, Double_t*);

// auxiliary
  void Vdef(double[4], const double, const double, const double, const double);
////////////////////////////////////////////////////////////////////////////
};// KKfoam



#endif
