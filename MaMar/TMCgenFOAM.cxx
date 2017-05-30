#include "TMCgenFOAM.h"
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////
///     TMCgenFOAM class
/// This is class for axiliary exercises,  mainly integration with Monte Carlo

ClassImp(TMCgenFOAM);

TMCgenFOAM::TMCgenFOAM():
  TMCgen()
{
  /// This constructor is for ROOT streamers ONLY
  cout<< "----> TMCgenFOAM Default Constructor (for ROOT only) "<<endl;
  m_Foam3 = NULL;
}

///______________________________________________________________________________________
TMCgenFOAM::~TMCgenFOAM()
{
  //!Explicit destructor
  cout<< "----> TMCgenFOAM::TMCgenFOAM !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

///_____________________________________________________________
TMCgenFOAM::TMCgenFOAM(const char* Name):
  TMCgen(Name)
{
//! all defaults defined here can be changed by the user
//! before calling TMCgen::Initialize
  m_Foam3 = NULL;
///////////////////////////////////////////////////
// Physics
  m_gnanob  = 389.37966e3;
  m_pi      = 3.1415926535897932;
  m_ceuler  = 0.57721566;
  m_alfinv  = 137.035;
  m_alfpi   = 1/m_alfinv/m_pi;
  m_amel    = 0.510999e-3;
//
  m_beam    = 0.510999e-3;  // electron
  m_chini   = 1.0;          // electron
//
  m_fin     = 0.105;        // final ferm. muon mass
  m_chfin   = 1.0;          // final ferm. muon charge
//
  m_KFini   = 11;           // electron
  m_KFf     = 13;           // muon
//
  m_KeyISR  = 2;            // Type of ISR/QED switch, 0,1,2
  m_jmax    = 10000;        // length of xpar
///////////////////////////////////////////////////
/// Foam setup
  m_kDim    =    5;         // No. of dim. for Foam, =2,3 Machine energy spread OFF/ON
  m_nCells  = 2000;         // No. of cells, optional, default=2000
  m_nSampl  =  200;         // No. of MC evts/cell in exploration, default=200
  m_del     = 0.0001;       // limit for gamma*ln(eps) in IFI mapping
  m_Mode    = 5;
///////////////////////////////////////////////////
// debug
  m_count   =0;

cout<< "----> TMCgenFOAM::TMCgenFOAM USER Constructor "<<endl;
}///

///______________________________________________________________________________________
void TMCgenFOAM::Initialize(TRandom *RNgen, ofstream *OutFile, TH1D* h_NORMA)
{
  cout<< "----> TMCgenFOAM::Initialize, Entering "<<endl;
  ///	      SETTING UP RANDOM NUMBER GENERATOR
  TMCgen::Initialize(  RNgen, OutFile, h_NORMA);

/////////////////////////////////////////////////////////
//  m_NevGen=0;
  const int jmax =m_jmax;
  ReaData("../../.KK2f_defaults", jmax, m_xpar);  // numbering as in input!!!
  ReaData("./pro.input",         -jmax, m_xpar);  // jmax<0 means no-zeroing
  double ypar[jmax];
  for(int j=0;j<jmax;j++) ypar[j]=m_xpar[j+1];    // ypar has c++ numbering
  //
  //NevTot = (long)m_xpar[0];                       // NevTot hidden in xpar[0] !!!
  m_CMSene  = m_xpar[ 1];
  m_vvmax   = m_xpar[17];

  cout<<" TMCgen::Initialize: m_CMSene="<<m_CMSene<<endl;
  cout<<" TMCgen::Initialize: m_vvmax="<<m_vvmax<<endl;

  const char *output_file = "./kkmc.output";
  long stl2 = strlen(output_file);
  int mout = 16;
  kk2f_fort_open_(mout,output_file,stl2);
  kk2f_initialize_(ypar);

  /////////////////////////////////////////////////////////
  if(f_IsInitialized == 0)
  {
  /// ******  SETTING UP FOAM of base class  *****
  f_FoamI   = new TFOAM("FoamI");   // new instance of MC generator FOAM
  m_kDim    = 5;
  m_nCells  =  10000;
  m_nSampl  = 100000;
  f_FoamI->SetkDim(m_kDim);         // No. of dims. Obligatory!
  f_FoamI->SetnCells(m_nCells);     // No. of cells, optional, default=2000
  f_FoamI->SetnSampl(m_nSampl);     // No. of MC evts/cell in exploration, default=200
  f_FoamI->SetnBin(        16);     // No. of bins default 8
  f_FoamI->SetOptRej(0);            // wted events (=0), default wt=1 events (=1)
  m_Mode    = 5;
  f_FoamI->Initialize( f_RNgen, this);     // Initialize FOAM
  //////////////////////////////////////////////////////////////
  /// ******  SETTING UP additional FOAM of the user class *****
  m_Foam3   = new TFOAM("Foam3");   // new instance of MC generator FOAM
  m_kDim    = 3;
  m_Foam3->SetkDim(m_kDim);         // No. of dims. Obligatory!
  m_Foam3->SetnCells(m_nCells);     // No. of cells, optional, default=2000
  m_Foam3->SetnSampl(m_nSampl);     // No. of MC evts/cell in exploration, default=200
  m_Foam3->SetnBin(        16);     // No. of bins default 8
  m_Foam3->SetOptRej(0);            // wted events (=0), default wt=1 events (=1)
  m_Mode    = 3;
  m_count =0;
  m_Foam3->Initialize( f_RNgen, this);     // Initialize FOAM
  //////////////////////////////////////////////////////////////
  double errel;
  f_FoamI->GetIntNorm(m_Xnorm,errel);   // universal normalization
  m_Foam3->GetIntNorm(m_Xsav3,errel);
  //screen output
  BXOPE(*f_Out);
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"======   TMCgenFOAM::Initialize     ======");
  BXTXT(*f_Out,"========================================");
  BXTXT(*f_Out,"========================================");
  f_IsInitialized = 1;  /// re-initialization inhibited
  } else {
	  cout<< "----> TMCgenFOAM::Initialize, already initialized "<<endl;
  }
}/// Initialize

///______________________________________________________________________________________
void TMCgenFOAM::Generate()
{
  f_NevGen++;
  m_Mode = -5;
  f_FoamI->MakeEvent();         // Foam of base class
  m_WT   = f_FoamI->GetMCwt();  // get weight

  ///  Fill special normalization histogram f_TMCgen_NORMA
  f_TMCgen_NORMA->Fill(-1, m_Xnorm);    // New style
//  f_TMCgen_NORMA->Fill(0.5, m_Xnorm);   // 1-st bin = Normal*Nevtot
//  f_TMCgen_NORMA->Fill(1.5, 1);         // 2-nd bin = Nevtot

}//! Generate

///______________________________________________________________________________________
void TMCgenFOAM::Finalize()
{
  TMCgen::Finalize();
  ///   Finalize MC  run, final printouts, cleaning etc.
  BXOPE(*f_Out);
  BXTXT(*f_Out,"****************************************");
  BXTXT(*f_Out,"******     TMCgenFOAM::Finalize     ******");
  BXTXT(*f_Out,"****************************************");
  ///------------------------
  Double_t MCresult, MCerror;
  f_FoamI->GetIntegMC( MCresult, MCerror);  //! get MC integral, should be one
  cout << "**************************************************************"<<endl;
  cout << "**************** Foam::Finalize  ************************"<<endl;
  cout << "Directly from FOAM: MCresult= " << MCresult << " +- "<<MCerror <<endl;
  cout << "**************************************************************"<<endl;
  ///------------------------
}//!Finalize



///------------------------------------------------------------------------
double TMCgenFOAM::Fyfs(double gam){
	  return exp(-m_ceuler*gam)/TMath::Gamma(1+gam);
}

///------------------------------------------------------------------------
double TMCgenFOAM::gamISR( double svar){
	  return  sqr(m_chini)*2*m_alfpi*( log(svar/sqr(m_beam)) -1);
}

///------------------------------------------------------------------------
double TMCgenFOAM::gamFSR( double svar){
	  return  sqr(m_chfin)*2*m_alfpi*( log(svar/sqr(m_fin)) -1);
}

///------------------------------------------------------------------------
double TMCgenFOAM::gamIFI( double costhe){
	  return  m_chini*m_chfin*2*m_alfpi* log( (1-costhe)/(1+costhe));
}

///------------------------------------------------------------------------
void TMCgenFOAM::MapIFI1( double r, double gam, double eps, double &v, double &dJac){
// Maping for POSITIVE gam
// Input r in (0,1) is random number
// Returned v is distributed according to gam*v^{gam-1}
  if( fabs(gam*log(eps)) > m_del){
	  double eg = exp(gam*log(eps));
	  if( r< eg ){
		  v = 0;  dJac=1/eg;
	  } else {
		  v = exp((1/gam)*log(r)); // mapping
		  dJac = 1/(r*gam/v);      // jacobian
	  }
  } else {
	  double eg = 1+gam*log(eps);
	  if( r< eg ){
		  v = 0; dJac=1/eg;
	  } else {
		  v = exp(-(1/gam)*(1-r)); // mapping
		  dJac = 1/(gam/v);        // jacobian
	  }
  }
  if( v<0 || v>1) {
	  cout<<"STOP in TMCgenFOAM::MapIFI: +++ v = "<<v<<endl;
	  exit(11);
  }
}// MapIFI

///------------------------------------------------------------------------
void TMCgenFOAM::MapIFI2( double r, double gam, double eps, double &v, double &dJac){
// Maping for NEGATIVE gam
// Input r in (0,1) is random number
// Returned v is distributed according to gam*v^{gam-1}
// dJac is normalization (part of Jacobian) factor
  double r0, eg;
  if( fabs(gam*log(eps)) > m_del){
	  eg = exp(gam*log(eps));
	  r0 =eg/(2*eg-1);
	  if( r< r0 ){
		  v = 0; dJac= 1/r0;
	  } else {
		  v = exp( (1/gam)*log( 2*eg -(2*eg-1)*r ) ); // mapping
		  dJac = (2*eg-1)/(-gam/v*exp(gam*log(v)));   // jacobian
	  }
  } else {
	  eg = 1+gam*log(eps);
	  r0 = eg/(2*eg-1);
	  if( r< r0){
		  v = 0; dJac= 1/r0;
	  } else {
		  v = exp( (1/gam)*(1-r)*(2*eg-1) ); // mapping
		  dJac = (2*eg-1)/(-gam/v);          // jacobian
	  }
  }
  if( v<0 || v>1) {
	  cout<<"STOP in TMCgenFOAM::MapIFI2: +++ v = "<<v<<endl;
	  exit(11);
  }
}// MapIFI2

///------------------------------------------------------------------------
void TMCgenFOAM::MapIFI( double r, double gam, double eps, double &v, double &R){
//// mapping for IFI
if(gam > 0)
	MapIFI1( r, gam, eps, v, R);
else
    MapIFI2( r, gam, eps, v, R);
}// MapIFI


///------------------------------------------------------------------------
double TMCgenFOAM::Rho_isr(double svar, double vv){
/// ISR rho-function for ISR

  double alf1   = m_alfpi;
  double gami   = gamISR(svar);
  //gami = sqr(m_chini)*2*m_alfpi*( log(svar/sqr(m_beam)) -1);
///
  double gamfac = Fyfs(gami);
  double delb   = gami/4 +alf1*(-0.5  +sqr(m_pi)/3.0);
  double ffact  = gamfac*exp(delb);

  double rho,dels,delh;
  if(       m_KeyISR == 0){
/// zero   order exponentiated
	dels = 0;
	delh = 0;
	rho  = ffact*gami* exp( log(vv)*(gami-1) ) *(1 +dels +delh);
  }else if( m_KeyISR == 1){
/// first  order
	dels = gami/2;   /// NLO part =0 as for vector boson???
    delh = vv*(-1 +vv/2);
    rho = ffact*gami* exp( log(vv)*(gami-1) ) *(1 +dels +delh);
  }else if( m_KeyISR == 2){
/// second order without NLO part
    dels = gami/2 +sqr(gami)/8;
    delh = vv*(-1+vv/2.0)
          +gami*0.5*(-0.25*(4.0-6.0*vv+3.0*vv*vv)*log(1-vv)-vv);
    rho = ffact*gami* exp( log(vv)*(gami-1) ) *(1 +dels +delh);
  }else{
	  cout<<"+++++TMCgenFOAM::KKdistr: Wrong KeyISR = " << m_KeyISR<<endl;
	  exit(5);
  }
///
  return rho;
}//Rho_isr


///------------------------------------------------------------------------
double TMCgenFOAM::Rho_fsr(double svar, double uu){
/// ISR+FSR rho-function

  double alf1   = m_alfpi;
  double gamf   = gamFSR(svar*(1-uu));
///
  /////delb   = Chf2* alf1*(0.5d0*bilg -1d0  +pi**2/3d0)
  /////delb   = delb -betf/2 *dlog(1-uu)
  /////double gamfac = exp(-m_ceuler*gamf)/TMath::Gamma(1+gamf);
  double delb   = gamf/4 +alf1*(-0.5  +sqr(m_pi)/3.0)
		         -gamf/2 *log(1-uu);
  double ffact  = Fyfs(gamf)*exp(delb);

  double rho,dels,delh;
  if(       m_KeyISR == 0){
/// zero   order exponentiated
	dels = 0;
	delh = 0;
	rho  = ffact*gamf* exp( log(uu)*(gamf-1) ) *(1 +dels +delh);
  }else if( m_KeyISR == 1){
/// first  order
	dels = gamf/2;   /// NLO part =0 as for vector boson???
    delh = uu*(-1 +uu/2);
    rho = ffact*gamf* exp( log(uu)*(gamf-1) ) *(1 +dels +delh);
  }else if( m_KeyISR == 2){
/// second order without NLO part
    /////dels  = betf/2d0 +betf**2/8d0
    /////delh  = uu*(-1d0+uu/2d0)
    /////$        +betf*(-0.5d0*uu-0.25d0*uu*(-1d0+0.5d0*uu)*log(1d0-uu))
    dels = gamf/2 +sqr(gamf)/8;
    delh = uu*(-1+uu/2.0)
          +gamf*0.5*( -0.5*uu -0.25*uu*(-1.0 +0.5*uu)*log(1-uu));
    rho = ffact*gamf* exp( log(uu)*(gamf-1) ) *(1 +dels +delh);
  }else{
	  cout<<"+++++TMCgenFOAM::KKdistr: Wrong KeyISR = " << m_KeyISR<<endl;
	  exit(5);
  }
//
  //if(rho<0) rho=0;
  return rho;
}//Rho_fsr


///--------------------------------------------------------------
double TMCgenFOAM::Rho_ifi(double costhe, double uu, double eps){
/// ISR+FSR rho-function
  double rho, gami;
  gami   = gamIFI( costhe );
  if( fabs(gami*log(eps)) > m_del){
	  if( uu < eps){
		  rho = exp( log(eps)*gami );
	  } else {
		  rho = gami*exp( log(uu)*(gami-1) );
	  }
  } else {
	  if( uu < eps){
		  rho = 1+ gami*log(eps);
	  }  else {
		  rho = gami/uu;
	  }
  }
  rho *= Fyfs(gami);
  return rho;
}//Rho_ifi



///________________________________________________________________________
double TMCgenFOAM::Density(int nDim, double *Xarg){
	//
	if( abs(m_Mode) == 5 ){
	    return Density5(nDim, Xarg);
	} else if( abs(m_Mode) == 3 ){
		return Density3(nDim, Xarg);
	} else {
		cout<<" TMCgenFOAM::Density: wrong Mode ="<<m_Mode<<endl;
		exit(-9);
	}
}// Density

///////////////////////////////////////////////////////////////
double TMCgenFOAM::Density5(int nDim, double *Xarg)
{ // density distribution for Foam
	m_count++;  // counter for debug
	//
	Double_t Dist=1;
//	cout<< " Density5: Dist= "<< Dist <<endl;  ///%%%

	double svar = sqr(m_CMSene);
	double svarCum = svar;

// ******** mapping for ISR *******
	double gamiCR,gami,alfi;
	double CMSene1= sqrt(svar);
	bornv_makegami_( CMSene1, gamiCR,gami,alfi);   // from KKMC

	double R= Xarg[0];
	m_vv = exp(1.0/gami *log(R)) *m_vvmax; // mapping
	Dist *= m_vv/R/gami ;                  // Jacobian
	if( gami < 0 )      return 0.0;    // temporary fix
    if( m_vv < 1e-200 ) return 0.0;    // temporary fix
    // ISR photonic distribution
	double Rho2 = Rho_isr(svar,m_vv);  // remember take care of m_mbeam!!!
	Dist *= Rho2;
	svarCum *= (1-m_vv);
	double svar2 = svar*(1-m_vv);

// ******** mapping for FSR *******
    double rr= Xarg[1];
    double gamf   = gamFSR(svar2);
    m_uu = exp(1.0/gamf *log(rr));     // mapping
    Dist *= m_uu/rr/gamf;              // Jacobian
	if( gamf < 0 )      return 0.0;    // temporary fix
    if( m_uu < 1e-200 ) return 0.0;    // temporary fix
    // FSR photonic distribution
  	double Rho3 = Rho_fsr(svar2,m_uu);           // remember take care of m_mbeam!!!
  	Dist *= Rho3;
    svarCum *= (1-m_uu);

    // ******** mapping for polar angle *******
    double cmax = 0.999;
    m_CosTheta = cmax*( -1.0 + 2.0* Xarg[2] );
    Dist *= 2.0*cmax;

    // ******** mapping for IFI variaable *******
    double eps =1e-6;
    double gamint = gamIFI(m_CosTheta);
//
    double R1, R2;
    MapIFI( Xarg[3], gamint, eps, m_r1, R1);     // mapping
    MapIFI( Xarg[4], gamint, eps, m_r2, R2);     // mapping
    double RhoIFI1 = Rho_ifi( m_CosTheta, m_r1, eps);
    double RhoIFI2 = Rho_ifi( m_CosTheta, m_r2, eps);
//
//    m_r1 =0;
//    m_r2 =0;
    double WT1 = R1 *RhoIFI1;
    double WT2 = R2 *RhoIFI2;
    Dist *= WT1*WT2;
//
// ******* MC event *******
    double zz = (1-m_vv)*(1-m_uu)*(1-m_r1)*(1-m_r2);
    m_xx = 1-zz;
//    m_xx = m_vv + m_uu - m_vv*m_uu;  // numerically more stable
// effective masses
	m_Mka = sqrt(svar2);   // after ISR

// =============== Sigm/dOmega from spin amplitudes ===============
// Effective 4-momenta, KKMC convention: p={px,py,pz,E)
	double Ene = m_Mka/2;
	double Pmb  = sqrt( (Ene-m_beam)*(Ene+m_beam) ); // modulus
	Vdef(m_p1, 0, 0 , Pmb, Ene);  // beam
	Vdef(m_p2, 0, 0 ,-Pmb, Ene);  // beam
	double Pmf  =sqrt( (Ene-m_fin)*(Ene+m_fin) ); // modulus
	Vdef(m_p3, Pmf*sqrt(1-sqr(m_CosTheta)), 0 , Pmf*m_CosTheta,  Ene); // final
	Vdef(m_p4,-Pmf*sqrt(1-sqr(m_CosTheta)), 0 ,-Pmf*m_CosTheta,  Ene); // final
	double PX[4] = {0, 0, 0, 2*Ene};
	double dSigAngF,dSigAngF1,dSigAngF2, Misr1,Misr2;
	Misr1 = sqrt((1-m_vv)*(1-m_r1)*svar);
	Misr2 = sqrt((1-m_vv)*(1-m_r2)*svar);
	gps_bornfoam_( 0,m_KFini,m_KFf,Misr1,m_CosTheta,dSigAngF1);
	gps_bornfoam_( 1,m_KFini,m_KFf,Misr2,m_CosTheta,dSigAngF2);
    dSigAngF = gps_makerhofoam_(1.0);
//************ Debug*** Debug*** Debug*** Debug*** Debug ***********
    if( m_count <1 && fabs(svar/svar2-1)>0.20 ){  // debug
//    if( m_count <1000 ){  // debug
    	double Rat;
    	Rat = dSigAngF1/( dSigAngF2 );
    	cout<<" Density5 debug m_count= "<< m_count<< endl;
    	cout<<" dSigAngF1    = "<< dSigAngF1;
    	cout<<" dSigAngF2    = "<< dSigAngF2;
    	cout<<" svar/svar2 = "<< svar/svar2;
    	cout<<" Rat = "<<Rat<<endl;
    } //
//    if( m_count <10000 && m_r1 > 0 && m_r2 >0 ){  // debug
    if( m_count <1 && m_r1 > 0 && gamint <0 ){  // debug
    	cout<<" Density5 debug m_count= "<< m_count<< endl;
    	cout<<" m_r1= "<< m_r1 <<"  m_r2="<< m_r2<<"  m_xx="<< m_xx <<endl;
    	cout<<" m_CosTheta ="<< m_CosTheta <<" gamint= "<<gamint<<endl;
    	cout<<" WT1 ="<< WT1 <<"  WT2="<< WT2<<endl;
   }
//   **********  end debug **********
	double sig0nb = 4*m_pi* sqr(1/m_alfinv)/(3.0*svar2 )*m_gnanob;
	Dist *=  dSigAngF *3.0/8.0 *sig0nb;

//	if( Dist < 0 ){
//		cout<< " Density5: Dist= "<< Dist <<endl;  ///%%%
//	}

	if( svarCum < sqr(2*m_fin)) Dist = 1e-100;

	if(m_Mode > 0 ) Dist = fabs(Dist); // For initialization mode

	return Dist;
}// Density5



///////////////////////////////////////////////////////////////
Double_t TMCgenFOAM::Density3(int nDim, Double_t *Xarg)
{ // density distribution for Foam
	m_count++;  // counter for debug
	//
	Double_t Dist=1;

	double svar = sqr(m_CMSene);
	double svarCum = svar;

// ******** mapping for ISR *******
	double gamiCR,gami,alfi;
	double CMSene1= sqrt(svar);
	bornv_makegami_( CMSene1, gamiCR,gami,alfi);   // from KKMC
	//[[[[ debug
	//gami = gamISR(CMSene1);
	//]]]]
	// cout<<" CMSene1,gami= "<< CMSene1 <<"  "<< gami <<endl;
	double R= Xarg[0];
	m_vv = exp(1.0/gami *log(R)) *m_vvmax; // mapping
	Dist *= m_vv/R/gami ;                  // Jacobian
	if( gami < 0 )      return 0.0;    // temporary fix
    if( m_vv < 1e-200 ) return 0.0;    // temporary fix
    // ISR photonic distribution
	double Rho2 = Rho_isr(svar,m_vv);  // remember take care of m_mbeam!!!
	Dist *= Rho2;
	svarCum *= (1-m_vv);
	double svar2 = svar*(1-m_vv);

// ******** mapping for FSR *******
    double rr= Xarg[1];
    double gamf   = gamFSR(svar2);
    if( gamf <0 )       return 0.0;      // just in case
    m_uu = exp(1.0/gamf *log(rr));       // mapping
    if( m_uu < 1e-200 ) return 0.0;      // temporary fix
    Dist *= m_uu/rr/gamf;                // Jacobian
// FSR photonic distribution
  	double Rho3 = Rho_fsr(svar2,m_uu);   // remember take care of m_mbeam!!!
  	if( Rho3 <0 ) return 1e-100;
 	Dist *= Rho3;
    svarCum *= (1-m_uu);

    // ******** mapping for polar angle *******
    m_CosTheta = -1.0 + 2.0* Xarg[2];
    Dist *= 2.0;

    double zz = (1-m_vv)*(1-m_uu);
    m_xx = 1-zz;
    m_xx = m_vv + m_uu - m_vv*m_uu;  // numerically more stable

// ******** finally Born factor *******
    long KeyFob;
    KeyFob=   10; // BornV_Dizet, with EW and without integration ???
    KeyFob=  -11; // BornV_Simple, for KeyLib=0, NO EW, NO integration OK
    KeyFob=  -10; // KKsem_BornV, NO EW, NO integration OK!
    KeyFob= -100; // KKsem_BornV, NO EW, WITH integration, OK
    KeyFob=    0; // With EW (BornV_Dizet) With integration OK!
//  -----------------
//	kksem_setkeyfob_( KeyFob );

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//double xBorn;
	//Integrated Born from KKMC
	//kksem_makeborn_( svar2, xBorn);
	//Dist *= xBorn/2.0;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//    In BornV_Differential:
//    CALL BornV_InterpoGSW( ABS(KFf),  svar, CosThe)
//    Born= BornV_Dizet( 1,m_KFini,KFf, svar, CosThe, eps1,eps2,ta,tb)
//    Born = 4*pi*alfa**2/(3d0*svar )*BornY  *m_gnanob

//
	bornv_interpogsw_(m_KFf,svar2, m_CosTheta);
	double dSig_dCos = bornv_dizet_( 1, m_KFini, m_KFf, svar2, m_CosTheta, 0.0, 0.0, 0.0, 0.0);

// ******* effective masses *********
	m_Mka = sqrt(svar2);   // after ISR

// =============== Sigm/dOmega from spin amplitudes ===============
// Effective 4-momenta, KKMC convention: p={px,py,pz,E)
	double Ene = m_Mka/2;
	double Pmb  = sqrt( (Ene-m_beam)*(Ene+m_beam) ); // modulus
	Vdef(m_p1, 0, 0 , Pmb, Ene);  // beam
	Vdef(m_p2, 0, 0 ,-Pmb, Ene);  // beam
	double Pmf  =sqrt( (Ene-m_fin)*(Ene+m_fin) ); // modulus
	Vdef(m_p3, Pmf*sqrt(1-sqr(m_CosTheta)), 0 , Pmf*m_CosTheta,  Ene); // final
	Vdef(m_p4,-Pmf*sqrt(1-sqr(m_CosTheta)), 0 ,-Pmf*m_CosTheta,  Ene); // final
	double PX[4] = {0, 0, 0, 2*Ene};
//***** pure Born of CEEX
	double dSigAng;
    gps_bornf_(m_KFini, m_KFf ,PX, m_CosTheta, m_p1,m_beam, m_p2, -m_beam,
                                               m_p3,m_fin,  m_p4, -m_fin,   dSigAng);
//[[[[[[[
//    double dSigRef = bornv_dizet_( 1, m_KFini, m_KFf, svar2, 0.0 , 0.0, 0.0, 0.0, 0.0); // at cos(theta)=0
//************ Debug*** Debug*** Debug*** Debug*** Debug ***********
//    if( m_count <10 && m_xx>0.6 ){  // debug
//    if( m_count <100000 && fabs(dSigAng/dSig_dCos -1) >0.10 ){  // debug
//    if( m_count <10000 && fabs(dSigAng-dSig_dCos)/dSigRef >0.002 ){  // debug
//    	cout<<" ******************** Density3 debug m_count= "<< m_count<< endl;
//    	cout<<" (dSigAng-dSig_dCos)/ref  = "<< (dSigAng-dSig_dCos)/dSigRef ;
//   	  // Born+boxes, WARNING Z-box may be modified for KeyZet=2
//      double dSigAngF0,dSigAngF1;
//    	gps_bornfoam_( 0,m_KFini,m_KFf,m_Mka,m_CosTheta,dSigAngF0);
//    	cout<<" dSigAngF/dSig_dCos = "<< dSigAngF/dSig_dCos;
//    	gps_bornfoam_( 1,m_KFini,m_KFf,m_Mka,m_CosTheta,dSigAngF1);
//      double dSigAngFF = gps_makerhofoam_(1.0);
//      cout<<" // dSigAngFF: "<< (dSigAngFF-dSig_dCos)/dSigRef;
//      cout<<"    dSigAngF0: "<< (dSigAngF0-dSig_dCos)/dSigRef;
//  	  cout<<" m_CosTheta= "<< m_CosTheta;
//      cout<<" m_Mka= "<< m_Mka;
////    	cout<<" m_vv= "<< m_vv;
////    	cout<<" m_uu= "<< m_uu;
//      cout<<endl;
//    } // end debug **********
//
	double sig0nb = 4*m_pi* sqr(1/m_alfinv)/(3.0*svar2 )*m_gnanob;
//	Dist *=  dSig_dCos *3.0/8.0 *sig0nb;  // Born of EEX
	Dist *=  dSigAng   *3.0/8.0 *sig0nb;  // Born of CEEX

	if( svarCum < sqr(2*m_fin)) Dist = 1e-100;

	return Dist;
}// Density3


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                             UTILITIES                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void TMCgenFOAM::ReaData(const char *DiskFile, int imax, double xpar[])
//////////////////////////////////////////////////////////////
//    subprogram reading input data file and packing        //
//    entries into matrix xpar                              //
//    WARNING: input file cannot include empty lines        //
//    it cannot handle entries like 1d-30, has to be 1e-30! //
//////////////////////////////////////////////////////////////
{
  char trail[200];
  char ch1;
  int  foundB=0, foundE=0, line, indx;
  int  line_max =2000;
  double value;
  cout<<"============================ReaData=============================="<<endl;
  cout<<"===                     "<< DiskFile <<"               =========="<<endl;
  cout<<"================================================================="<<endl;
  ifstream InputFile;
  InputFile.open(DiskFile);
  for(indx=0;indx<imax; indx++) xpar[indx]=0.0;
  for(line=0;line<line_max; line++){
    InputFile.get(ch1);
    if( ch1 == 'B') foundB=1;
    InputFile.getline(trail,200);
    if(foundB) break;
  }
  for(line=0;line<line_max; line++){
    InputFile.get(ch1);
    if( ch1 == 'E'){
      foundE=1;
      InputFile.getline(trail,200);
      cout<<ch1<<trail<<"["<<line<<"]"<<endl;
      break;
    }
    if( ch1 == '*'){
      InputFile.getline(trail,200);
      cout<<ch1<<trail<<endl;
    }else{
      InputFile>>indx>>value;
      if(indx<0 || indx>abs(imax) ){
	cout<<" ++++++++ReaData: wrong indx = "<<indx<<endl;
	exit(0);
      }
      xpar[indx] = value;
      //xpar[indx-1] = value; // correction for fortran indexing in input file
      InputFile.getline(trail,200);
      cout<<ch1;
      cout<<setw(4)<<indx<<setw(15)<<value<<" ";
      cout<<trail<<endl;
    }
  }
  cout<<"================================================================="<<endl;
  InputFile.close();
}

void TMCgenFOAM::Vdef(double v[4], const double v1, const double v2, const double v3, const double v4)
  { // define a 4-vector (avoids initialization warnings)
       v[0] = v1; v[1] = v2; v[2] = v3; v[3] = v4;
  } 
