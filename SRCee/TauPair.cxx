///////////////////////////////////////////////////////////////////////////////
// Replaces f77 module TauPair
///////////////////////////////////////////////////////////////////////////////
#include "TauPair.h"
#include "KKpol/anomwt2-Z.h"

ClassImp(TauPair);

TauPair::TauPair()
{
  // This constructor is for ROOT streamers ONLY
  // all pointers has to be NULLed
  cout<< "----> TauPair Default Constructor (for ROOT only) "<<endl;
  m_Out      = NULL;
  DB         = NULL;
  m_Event    = NULL;
  m_Hvent    = NULL;
  m_GPS      = NULL;
  m_RNgen    = NULL;

}

///_____________________________________________________________
TauPair::TauPair(ofstream *OutFile)
{
  cout<< "----> TauPair USER Constructor "<<endl;
  m_Out      = OutFile;
  DB         = NULL;
  m_Event    = NULL;
  m_Hvent    = NULL;
  m_GPS      = NULL;
  m_RNgen    = NULL;
}//TauPair

///______________________________________________________________________________________
TauPair::~TauPair()
{
  //Explicit destructor
  cout<< "----> TauPair::TauPair !!!! DESTRUCTOR !!!! "<<endl;
}///destructor

double TauPair::sqr( const Double_t x ){ return x*x;};

///______________________________________________________________________________________
void TauPair::Initialize(double xpar[])
{
  cout  << "----> TauPair::Initialize, Entering "<<endl;
//=================================================================
// BX*** macros are in MCdev/BXFORMAT.h
  BXOPE(*m_Out);
  BXTXT(*m_Out,"========================================");
  BXTXT(*m_Out,"======    TauPair::Initialize     ======");
  BXTXT(*m_Out,"========================================");

  int ITAUXPAR=2000;
  m_IsInitialized = xpar[415-1];  // General mask for tau channel
// switches of tau+ tau- decay modes !!
  m_itdkRC        = xpar[ITAUXPAR+4-1];   // QED internal rad. in leptonic decays
  int Jak1        = xpar[ITAUXPAR+1-1];   // Decay Mask for first tau
  int Jak2        = xpar[ITAUXPAR+2-1];   // Decay Mask for second tau
  if( (Jak1 == -1) && (Jak2 == -1) ) m_IsInitialized = 0;

  m_KeyClone      = 1;       // dip-switch for cloning procedure, =1,2
  m_KeyClone      = 2;       // dip-switch for cloning procedure, =1,2
  BXTXT(*m_Out, " KK interface of Tauola                 ");
  BX1I( *m_Out, "  IsInit", m_IsInitialized, "xpar[415]       =");
  BX1I( *m_Out, "    Jak1",            Jak1, "xpar[2001]      =");
  BX1I( *m_Out, "    Jak2",            Jak2, "xpar[2002]      =");
  BX1I( *m_Out, "  itdkRC",        m_itdkRC, "xpar[2004]      =");
  BX1I( *m_Out, "KeyClone",      m_KeyClone, "Cloning proc.   =");
  BXCLO(*m_Out);

// Initialisation of tau decay package TAUOLA; ITAUXPAR is for indirect adressing.
  inietc_(&ITAUXPAR,xpar);
  if( m_IsInitialized == 0) {
     BXOPE(*m_Out);
     BXTXT(*m_Out, " !!!!! Tauola inhibited !!!!    ");
     BXCLO(*m_Out);
  } else {
// Initialisation of TAUOLA
    inimas_(&ITAUXPAR,xpar);
    initdk_(&ITAUXPAR,xpar);
    double xk0qed = 0.1;            // <=== It seems to be never used
    iniphy_(&xk0qed);
    int JAK =-1;
    dekay_(&JAK, m_HvecTau1);
////////////////////////////////////////////////
// Initialization of PHOTOS++
// KeyPhts =0 for off; =1 in non-leptonic; =2 in all decays
    double WTmax=4.0;
    if(DB->KeyPhts ==2 ){
      Photos::initialize();
      Photos::maxWtInterference(WTmax);
    } else if( DB->KeyPhts ==1){
// Suppressing Photos for leptonic decays
      Photos::initialize();
      Photos::maxWtInterference(WTmax);
//////////////////////////////////////////
// Flag selections below do not work properly.
// Leptonic tau decays are detected in the hepmc3 events
// and photos++ does not process them for KeyPhts=1
//      Photos::suppressAll();
//      Photos::forceBremForBranch(0, 15);
//      Photos::forceBremForBranch(0, -15);
//      Photos::suppressBremForDecay(3,  15,  16,  11, -12); //tau- => nutau,    e-, nuelbar
//      Photos::suppressBremForDecay(3, -15, -16, -11,  12); //tau+ => nutaubar, e+, nuel
//      Photos::suppressBremForDecay(3,  15,  16,  13, -14); //tau- => mu-
//      Photos::suppressBremForDecay(3, -15, -16, -13,  14); //tau+ => mu+
//////////////////////////////////////////
    }//KeyPhts
  }//IsInitialized
///////////////////////////////////////////////////
}// end if Initialize

///______________________________________________________________________________________
void TauPair::DecayInRest(){
  int J;
  if( m_IsInitialized != 0) {
    J=1; dekay_(&J,m_HvecTau1); // TAUOLA
    J=2; dekay_(&J,m_HvecTau2); // TAUOLA
  }
}//Make1

///______________________________________________________________________________________
void TauPair::RandRotor(){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   This routine is strongly interrelated with Tralor  !!!                        //
//                                                                                 //
//   Cloning tau decays by additional rotation tau decay products with respect     //
//   to frames  initially used in the decay simulation.                            //
//   This is perfectly legal because average spin weight is equal exactly one!!!   //
/////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////
//   Generation of random two independent Euler rotations                          //
/////////////////////////////////////////////////////////////////////////////////////
  double rrr[10];
  m_RNgen->RndmArray(3, rrr);
  m_alfa1  = 2.0*M_PI*rrr[2];        // azimuthal angle in (0,2*pi)
  m_beta1  = acos(2.0*rrr[0]-1.0);   // polar angle     in (0,  pi)
  m_gamma1 = 2.0*M_PI*rrr[1];        // azimuthal angle in (0,2*pi)
//------------------------------------------------
  m_RNgen->RndmArray(3, rrr);
  m_alfa2  = 2.0*M_PI*rrr[2];        // azimuthal angle in (0,2*pi)
  m_beta2  = acos(2.0*rrr[0]-1.0);   // polar angle     in (0,  pi)
  m_gamma2 = 2.0*M_PI*rrr[1];        // azimuthal angle in (0,2*pi)
//------------------------------------------------
  m_H1.SetPxPyPzE(m_HvecTau1[0],m_HvecTau1[1],m_HvecTau1[2],m_HvecTau1[3]);
  m_H2.SetPxPyPzE(m_HvecTau2[0],m_HvecTau2[1],m_HvecTau2[2],m_HvecTau2[3]);
  if(m_KeyClone == 1) {
/////////////////////////////////////////////////////////////////////////////////////
//   Cloning tau decay with help of  Euler rotations FIRST method                  //
    double Habs1 = (m_H1.Vect()).Mag();
    double Habs2 = (m_H2.Vect()).Mag();
    // Standart phi, theta for polarimeter fectors, phi in (0,2*pi), theta in (0,pi)
    m_phi1  =0.0; m_thet1 =0.0;
    m_phi2  =0.0; m_thet2 =0.0;
    if(Habs1 > 1e-5) { m_phi1  = m_H1.Phi(); m_thet1 = m_H1.Theta(); }
    if(Habs2 > 1e-5) { m_phi2  = m_H2.Phi(); m_thet2 = m_H2.Theta(); }
    m_H1.SetPxPyPzE(0.0, 0.0, Habs1, 1.0);
    m_H2.SetPxPyPzE(0.0, 0.0, Habs2, 1.0);
    m_Event->RotEul(m_beta1, m_gamma1, &m_H1);
    m_Event->RotEul(m_beta2, m_gamma2, &m_H2);
    for(int i=0; i<4;i++) m_HvClone1[i]=m_H1[i];
    for(int i=0; i<4;i++) m_HvClone2[i]=m_H2[i];
    //
  } else if(m_KeyClone == 2) {
/////////////////////////////////////////////////////////////////////////////////////
//   Cloning tau decay with help of  Euler rotations, SECOND method                //
     m_Event->RotEuler(m_alfa1, m_beta1, m_gamma1, &m_H1);
     m_Event->RotEuler(m_alfa2, m_beta2, m_gamma2, &m_H2);
     for(int i=0; i<4;i++) m_HvClone1[i]=m_H1[i];
     for(int i=0; i<4;i++) m_HvClone2[i]=m_H2[i];
 } else {
     (*m_Out)<< " ##### STOP in Taupair_Clone: wrong KeyClone= "<< m_KeyClone<< endl;
     cout <<    " ##### STOP in Taupair_Clone: wrong KeyClone= "<< m_KeyClone<< endl;
     exit(9);
  }
}//Clone


///______________________________________________________________________________________
void TauPair::ImprintSpin(){
//////////////////////////////////////////////////
//     introduces spin effects by rejection     //
//////////////////////////////////////////////////
  int loop=0;
  double rn,wt,wt0,wt1,wt2, wtmax=4.0;
e1099:
  loop=loop+1;
  RandRotor();   // Cloning tau decay by Euler rotation
  m_GPS->MakeRho2(m_HvClone1,m_HvClone2,wt0,wt1,wt2);
  wt = wt1;                         // why not wt2???
  rn = m_RNgen->Rndm();
  if (wt < wtmax*rn  && loop<100) goto e1099;
  // Save helicity information in the HepMC3 record as particle (Tau, Tau-) attributes
  if (m_Hvent) {
    std::vector<float> hel;
    for (auto p: m_Hvent->particles()){
      if (p->pid() == 15){ // Tau-
        hel.clear(); for(int i=0; i<3;i++) hel.push_back(m_HvClone2[i]);
        // Create attribute
        p->add_attribute("spin", std::make_shared<VectorFloatAttribute>(hel));
      } else if (p->pid() == -15){ // Tau+
        hel.clear(); for(int i=0; i<3;i++) hel.push_back(m_HvClone1[i]);
        // Create attribute
        p->add_attribute("spin", std::make_shared<VectorFloatAttribute>(hel));
      }
    }// end loop over particles
  } else {
     cout << "m_Hvent not defined! Cannot save helicity information" << endl;
  }

}//ImprintSpin

///______________________________________________________________________________________
void TauPair::TransExport(){
// Transforming decays from tau rest frame to LAB
// and appending event record (hepmc3) with tau decay particles
// Replacement for f77 taupair_make2_();
  int ih1,ih2;
  hepevt_getf_(   ih1);          // fermion is here
  hepevt_getfbar_(ih2);          // antifermion is here
  tauface_setfermpos_(ih1,ih2);  // set ffbar positions for /hepevt/ in Tauola
/////////// IMPORTANT!!! /////////////////////
// Inside fortran subroutine DEKAY of Tauola
// TauPair::Tralo4() and HepFace::FillHep3 are called!
// Interfaced from C++ into F77 through SRCee/globux.h
  int J;
  J=11;  dekay_(&J,m_HvClone1);
  J=12;  dekay_(&J,m_HvClone1);

}//Make2

/////////////////////////////////////////////////////////////////////////////////////
/// Run Photos
void TauPair::RunPhotosPP(){
// KeyPhts is flag for Photos c++
  if(DB->KeyPhts > 0) {

// test print before photos
  int buf= -m_Hvent->particles().size();
  int LimitPrint=50; // for debug only
  if(m_Event->m_EventCounter <= LimitPrint && DB->LevPri==3){
// test print before photos
    cout <<    "TauPair::RunPhotosPP:!!!!!!!!!!!!!!!!!!!!!!!!!! Before Photos !!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    (*m_Out) <<"TauPair::RunPhotosPP:!!!!!!!!!!!!!!!!!!!!!!!!!! Before Photos !!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    Print::listing(*m_Hvent);
  }//test print

  // check if tau has decayed leptonicaly (to avoid photon double counting)
  bool leptonicTauDecay = false;
  for (auto v : m_Hvent->vertices()){ // loop over vertices
    if (v->particles_in().size() == 1) { // search  for only decay of particles
        if (abs(v->particles_in()[0]->pid()) == 15){ // decay of tau+/-
          bool fermion = false;
          bool neutrino = false;
          for (auto p : v->particles_out()){ // check if leptonic decay
            if ( abs(p->pid()) == 11 || abs(p->pid()) == 13) fermion = true;
            if ( abs(p->pid()) == 12 || abs(p->pid()) == 14) neutrino = true;
            leptonicTauDecay = fermion & neutrino;
          }//if leptonic decay
        }//end tau decay
     }//end decay
   }//end verticles

//   if (!leptonicTauDecay)  m_TauGen->RunPhotosPP();   // Run PhotosPlusPlus for non leptonic dacays
//////////////////////////////////////////
// Beware: DB->KeyPhts ==2 requires m_itdkRC == 0
  if( DB->KeyPhts ==2 && m_itdkRC !=0 ){
    cout<<"TauPair::RunPhotosPP: +++STOP+++, KeyPhts ==2 && m_itdkRC ==1"<<endl;
    exit(33);
  }
// Process HEPMC3 event by PHOTOS++, Leptonic decays excluded,
// except the case of KeyPhts ==2 and m_itdkRC ==0
  if (!leptonicTauDecay || DB->KeyPhts ==2){
    PhotosHepMC3Event photosEvent(m_Hvent);
    photosEvent.process();
  }//
//////////////////////////////////////////

  // test print after photos
  buf += m_Hvent->particles().size();
  if(buf>0 && m_Event->m_EventCounter <= LimitPrint && DB->LevPri==3){
    cout<<      "TauPair::RunPhotosPP:!!!!!!!!!!!!!!!!!!!!!!!!!! After Photos !!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    (*m_Out) << "TauPair::RunPhotosPP:!!!!!!!!!!!!!!!!!!!!!!!!!! After Photos !!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    Print::listing(*m_Hvent);
    cout<<     ">>>>>>> TauPair::RunPhotosPP: ["<<m_Event->m_EventCounter<< "] PHOTOS++ added "<<buf<<" new photons !!!!!!"<<endl;
    (*m_Out) <<">>>>>>> TauPair::RunPhotosPP: ["<<m_Event->m_EventCounter<< "] PHOTOS++ added "<<buf<<" new photons !!!!!!"<<endl;
  }//test print


  }// if KeyPhts
}//end Run Photos

void TauPair::Finalize(){
////////////////////////////////
// Final printout from Tauola
////////////////////////////////
//  CALL DEKAY(100,HvecDummy)
//
    int JAK =100;
    dekay_(&JAK, m_HvecTau1);
//
}//Finalize

void TauPair::Tralo4(int Kto, float P[], float Q[], float &AM){
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//   This routine is strongly interrelated with Taupair::Clone!                    //
//                                                                                 //
//  SUBSITUTE OF TRALO4                                                            //
//  TRALO4 is called in TAUOLA --> /hepevt/ interface to boost from tau+-          //
//  restframe to lab. It includes rotations in tau rest frame due to spin effect   //
//  implementation                                                                 //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////
// locals
double  Pd[4];
//* ------------------------------------------------------------
AM = sqrt(abs( P[3]*P[3] -P[2]*P[2] -P[1]*P[1] -P[0]*P[0] ));  // Mass
//
for(int k=0; k<4; k++) Pd[k]=P[k]; // from REAL to DOUBLE PRECISION
//
if(m_KeyClone == 1) {
   m_PP.SetPxPyPzE(Pd[0],Pd[1],Pd[2],Pd[3]);
   if(   Kto == 1) {
      m_Event->RotEulInv( m_thet1, m_phi1,   &m_PP);
      m_Event->RotEul(    m_beta1, m_gamma1, &m_PP);
   }else if(Kto == 2) {
      m_Event->RotEulInv( m_thet2, m_phi2,   &m_PP);
      m_Event->RotEul(    m_beta2, m_gamma2, &m_PP);
   } else {
     (*m_Out)<<"###### STOP in TRALO4: Wrong Kto = "<<Kto<<endl;
     cout    <<"###### STOP in TRALO4: Wrong Kto = "<<Kto<<endl; exit(9);
   }
   for(int i=0; i<4;i++) Pd[i]=m_PP[i];
} else if(m_KeyClone == 2) {
   m_PP.SetPxPyPzE(Pd[0],Pd[1],Pd[2],Pd[3]);
   if(     Kto == 1) {
     m_Event->RotEuler(m_alfa1, m_beta1, m_gamma1, &m_PP);
   } else if( Kto == 2) {
     m_Event->RotEuler(m_alfa2, m_beta2, m_gamma2, &m_PP);
   } else {
   (*m_Out)<<"###### STOP in TauPair::Tralo4: Wrong Kto = "<<Kto<<endl;
   cout    <<"###### STOP in TauPair::Tralo4: Wrong Kto = "<<Kto<<endl; exit(9);
   }
   for(int i=0; i<4;i++) Pd[i]=m_PP[i];
} else {
   (*m_Out)<<"##### STOP in Taupair_Tralo4: wrong KeyClone="<<m_KeyClone<<endl;
   cout    <<"##### STOP in Taupair_Tralo4: wrong KeyClone="<<m_KeyClone<<endl; exit(9);
}
m_GPS->TralorDoIt(Kto,Pd,Pd);
// Translation from DOUBLE PRECISION  to REAL
for(int k=0; k<4; k++) P[k]=Pd[k];
//----------------------------------------------
}//Tralo4

void TauPair::FillBornVEssentials() {

   KKdbase *DD = m_GPS->DB;
   KKdizet *DZ = m_GPS->m_DZ;
   // Prepare the vectir information
   double GSW[KKdizet::m_poinG];
   for (int i = 1; i < KKdizet::m_poinG; i++) {
      GSW[i] = DZ->m_GSW[i].real() ;
   }
   int keyelw = DB->KeyElw;
   double swsq = DZ->D_swsq;
   double alfinv = DB->AlfinvZ;
   double mz = DZ->m_MZ;
   double gamz = DZ->D_GamZ;
   fillbornv_(&keyelw, &swsq, &alfinv, &mz, &gamz, GSW);
   return;
}

void TauPair::GetPolarizationInfo(double &wtME, double &wtSPIN, double &wtSPIN0) {
/*
   // Fill the common block
   int out = 6;
   int idyfs = 0;
   int ifphot = 0;
   filltaupair_(m_HvecTau1, m_HvecTau2, m_HvClone1, m_HvClone2,
                &m_beta1, &m_alfa1, &m_gamma1, &m_beta2, &m_alfa2, &m_gamma2,
                &m_phi1, &m_thet1, &m_phi2, &m_thet2, &out, &m_IsInitialized, &idyfs, &m_KeyClone, &ifphot);
   
   // Call the main function
   int iqed = 0;
   double Ar0(0), Ai0(0), Br0(0), Bi0(0);
   double wtME, wtSPIN, wtSPIN0;
   anomwt_(&iqed, &Ar0, &Ai0, &Br0, &Bi0, &wtME, &wtSPIN, &wtSPIN0);
*/
   // Main vectors
  // double h1[4], h2[4], hr1[4], hr2[4], p1[4], p2[4], am;
   double h1[4], h2[4], p1[4], p2[4];
   for(int i=0; i<4; i++) {
      h1[i] = m_HvClone1[i];
      h2[i] = m_HvClone2[i];
    //  hr1[i] = m_HvClone1[i];
    //  hr2[i] = m_HvClone2[i];
      p1[i] = 0;
      p2[i] = 0;
   }
   // Tau mass
   double massTau = 1.77705;
   p1[3] = massTau;
   p2[3] = massTau;

   // Rotations
   m_GPS->TralorPrepare(1);
   m_GPS->TralorPrepare(2);
   
   m_GPS->TralorDoIt(1, p1, p1);
   m_GPS->TralorDoIt(2, p2, p2);
   m_GPS->TralorDoIt(1, h1, h1);
   m_GPS->TralorDoIt(2, h2, h2);
   //Tralo4(1, hr1, hr1, am);
   //Tralo4(1, hr2, hr2, am);

   double qq[4], pb1[4], pb2[4];
   for(int i=0; i<4; i++) {
      qq[i] = p1[i] + p2[i];
      pb1[i] = 0;
      pb2[i] = 0;
   }
   pb1[3]= 1.0;
   pb1[2]= 1.0;
   pb2[3]= 1.0;
   pb2[2]=-1.0;
      
// We go to restframe of tau pair to work against bremstrahlungs
      
   KKceex::BostDQ(1, qq, h1, h1);
   KKceex::BostDQ(1, qq, h2, h2);
   KKceex::BostDQ(1, qq, p1, p1);
   KKceex::BostDQ(1, qq, p2, p2);
   KKceex::BostDQ(1, qq, pb1, pb1);
   KKceex::BostDQ(1, qq, pb2, pb2);
      
// We eliminate y-components of p1,p2
   double fi = KKceex::AngFi(p2[0],p2[1]);

   KKceex::RotoD3(-fi, h1, h1);
   KKceex::RotoD3(-fi, h2, h2);
   KKceex::RotoD3(-fi, p1, p1);
   KKceex::RotoD3(-fi, p2, p2);
   KKceex::RotoD3(-fi, pb1, pb1);
   KKceex::RotoD3(-fi, pb2, pb2);

// Set taus along z direction      
   double thet = KKceex::AngXY(p2[2], p2[0]);

   KKceex::RotoD2(-thet, h1, h1);
   KKceex::RotoD2(-thet, h2, h2);
   KKceex::RotoD2(-thet, p1, p1);
   KKceex::RotoD2(-thet, p2, p2);
   KKceex::RotoD2(-thet, pb1, pb1);
   KKceex::RotoD2(-thet, pb2, pb2);

// Set  beam difference to reaction plane  x-z
   double pbb[4];
   for(int i=0; i<4; i++) {
//      pbb[i] = pb1[i] - pb2[i]; // Inconsistece!
      pbb[i] = pb2[i] - pb1[i];
   }
      
   double fi1 = KKceex::AngFi(pbb[0],pbb[1]);

   KKceex::RotoD3(-fi1, h1, h1);
   KKceex::RotoD3(-fi1, h2, h2);
   KKceex::RotoD3(-fi1, p1, p1);
   KKceex::RotoD3(-fi1, p2, p2);
   KKceex::RotoD3(-fi1, pb1, pb1);
   KKceex::RotoD3(-fi1, pb2, pb2);

   KKceex::BostDQ(1, p1, h1, h1);
   KKceex::BostDQ(1, p1, h2, h2);

   double E = p1[3];
   double theta = TMath::ACos( -pb1[2] / sqrt(pb1[0]*pb1[0] + pb1[1]*pb1[1] + pb1[2]*pb1[2]));  // minus necessary for KKMC frames

// Fill BornV.h essentials (needed by dipolqqrijradcor_)
   FillBornVEssentials();

   int iqed(0), channel(1);  
   double ReA(0), ImA(0), ReB(0), ImB(0), ReX(0), ImX(0), ReY(0), ImY(0);
   double R0[4][4];

   dipolqqrijradcor_(&iqed, &E, &theta, &ReA, &ImA, &ReB, &ImB, &ReX, &ImX, &ReY, &ImY, R0, &channel);

   iqed = 1;
   double Ar0(0), Ai0(0), Br0(0), Bi0(0);
   ReA = Ar0; 
   ImA = Ai0;
   ReB = Br0;
   ImB = Bi0;
   ReX = Ar0;
   ImX = Ai0;
   ReY = Br0;
   ImY = Bi0;

   double R[4][4];
   dipolqqrijradcor_(&iqed, &E, &theta, &ReA, &ImA, &ReB, &ImB, &ReX, &ImX, &ReY, &ImY, R, &channel);

// Prepare outputs
   wtME = R[3][3] / R0[3][3];
   wtSPIN = 0;
   wtSPIN0 = 0;
   double sign[4] = {-1, 1, -1, 1}; // sign() introduces overall rotation around by pi around  y axis. Necessary for KKMC frames. 
   for(int i=0; i<4; i++) {
      for(int j=0; j<4; j++) {
          wtSPIN  = wtSPIN  + R [i][j] / R[3][3] * h1[i] * h2[j] * sign[i] * sign[j];
          wtSPIN0 = wtSPIN0 + R0[i][j] / R0[3][3]* h1[i] * h2[j] * sign[i] * sign[j];
      }
   }
   // Renormalize
   wtSPIN = wtSPIN / wtSPIN0;

  return;
}
