/////////////////////////////////////////////////////////////////////
#include "HisReMake.h"


void HisReMakeKKMC(TFile *DiskFileA){
  cout<<"==================================================================="<<endl;
  cout<<"================ HisReMakeKKMC  BEGIN  ============================"<<endl;
  //
  int NbMax=0;
  //DiskFileA->ls("");
  TH1D *HST_KKMC_NORMA = (TH1D*)DiskFileA->Get("HST_KKMC_NORMA");
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_vTrueMain") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_vTrueCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_vAlepCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_vXGenCeex2") );
  //
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_Cost1Ceex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_CosPLCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_CosPRCeex2") );
  HisNorm1(HST_KKMC_NORMA, (TH1D*)DiskFileA->Get("hst_CosPREex2") );
  //
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sca_vTcPR_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sca_vTcPR_Ceex2n") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sca_vTcPR_Eex2") );
  //
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sca_vTcPL_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sca_vTcPL_Ceex2n") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sca_vTcPL_Eex2") );
  //
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sct_vTcPR_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sct_vTcPR_Ceex2n") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sct_vTcPR_EEX2") );
  //
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sct_vTcPL_Ceex2") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sct_vTcPL_Ceex2n") );

  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sct_vTcPL_Ceex2") );

  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sct_vTcPL_Ceex0") );
  HisNorm2(HST_KKMC_NORMA, (TH2D*)DiskFileA->Get("sct_vTcPL_Ceex0n") );
  //

//]]]]]]]]]]]]]]]]
////////////////////////////////////////////////////////////////////
// Pure KKMC reprocessing part
// from bigger scattergram and restricted vmax<0.2
//////////////////////////////////////////////////////////////////
  cout<<"  Renormalizing  and reprocessing histograms from KKMC"<<endl;
  //****************************************************************************************
  // Distributions of v=vTrue<vmax<0.20, c=cos(theta) with 100 bins
  //****************************************************************************************
  TH2D *sct_vTcPR_Ceex2  = (TH2D*)DiskFileA->Get("sct_vTcPR_Ceex2");
  TH2D *sct_vTcPR_Ceex2n = (TH2D*)DiskFileA->Get("sct_vTcPR_Ceex2n");
  TH2D *sct_vTcPR_EEX2   = (TH2D*)DiskFileA->Get("sct_vTcPR_EEX2");
  TH2D *sct_vTcPL_Ceex2  = (TH2D*)DiskFileA->Get("sct_vTcPL_Ceex2");
  TH2D *sct_vTcPL_Ceex2n = (TH2D*)DiskFileA->Get("sct_vTcPL_Ceex2n");
  TH2D *sct_vTcPL_Ceex0  = (TH2D*)DiskFileA->Get("sct_vTcPL_Ceex0");
  TH2D *sct_vTcPL_Ceex0n = (TH2D*)DiskFileA->Get("sct_vTcPL_Ceex0n");
  //NbMax=45;         // cosThetaMax = 45/50=0.90 Now global variable
  // KKMC IFI on
  TH1D  *HTot2_vTcPR_EEX2   = HstProjV("HTot2_vTcPR_EEX2",sct_vTcPR_EEX2,NbMax);
  TH1D  *HAfb2_vTcPR_EEX2   = HstProjA("HAfb2_vTcPR_EEX2",sct_vTcPR_EEX2,NbMax);
  // KKMC IFI off
  TH1D  *HTot2_vTcPR_Ceex2  = HstProjV("HTot2_vTcPR_Ceex2",sct_vTcPR_Ceex2,NbMax);
  TH1D  *HAfb2_vTcPR_Ceex2  = HstProjA("HAfb2_vTcPR_Ceex2",sct_vTcPR_Ceex2,NbMax);
  // KKMC IFI off
  TH1D  *HTot2_vTcPR_Ceex2n = HstProjV("HTot2_vTcPR_Ceex2n",sct_vTcPR_Ceex2n,NbMax);
  TH1D  *HAfb2_vTcPR_Ceex2n = HstProjA("HAfb2_vTcPR_Ceex2n",sct_vTcPR_Ceex2n,NbMax);
  // KKMC IFI off
  TH1D  *HTot2_vTcPL_Ceex2  = HstProjV("HTot2_vTcPL_Ceex2",sct_vTcPL_Ceex2,NbMax);
  TH1D  *HAfb2_vTcPL_Ceex2  = HstProjA("HAfb2_vTcPL_Ceex2",sct_vTcPL_Ceex2,NbMax);
  // KKMC IFI off
  TH1D  *HTot2_vTcPL_Ceex2n = HstProjV("HTot2_vTcPL_Ceex2n",sct_vTcPL_Ceex2n,NbMax);
  TH1D  *HAfb2_vTcPL_Ceex2n = HstProjA("HAfb2_vTcPL_Ceex2n",sct_vTcPL_Ceex2n,NbMax);

  // KKMC IFI on
  TH1D  *HTot2_vTcPL_Ceex0  = HstProjV("HTot2_vTcPL_Ceex0",sct_vTcPL_Ceex0,NbMax);
  TH1D  *HAfb2_vTcPL_Ceex0  = HstProjA("HAfb2_vTcPL_Ceex0",sct_vTcPL_Ceex0,NbMax);
  // KKMC IFI on
  TH1D  *HTot2_vTcPL_Ceex0n = HstProjV("HTot2_vTcPL_Ceex0n",sct_vTcPL_Ceex0n,NbMax);
  TH1D  *HAfb2_vTcPL_Ceex0n = HstProjA("HAfb2_vTcPL_Ceex0n",sct_vTcPL_Ceex0n,NbMax);

///****************************************************************************************
// Pure KKMC reprocessing part
  int nbMax;
  cout<<"  Renormalizing  and reprocessing histograms from KKMC"<<endl;
  // Wide range, vmax<1.
  TH2D *sca_vTcPR_Eex2   = (TH2D*)DiskFileA->Get("sca_vTcPR_Eex2");
  TH2D *sca_vTcPR_Ceex2  = (TH2D*)DiskFileA->Get("sca_vTcPR_Ceex2");
  TH2D *sca_vTcPR_Ceex2n = (TH2D*)DiskFileA->Get("sca_vTcPR_Ceex2n");
//****************************************************************************************
// Distributions of v=vTrue<1.0 unlimited c=cos(theta), 50 bins
//****************************************************************************************
//    nbMax=0;            // cosThetaMax = 1.0
// KKMC IFI on
  TH1D  *HTot_vTcPR_Eex2    = HstProjV("HTot_vTcPR_Eex2",sca_vTcPR_Eex2,NbMax);
  TH1D  *HAfb_vTcPR_Eex2    = HstProjA("HAfb_vTcPR_Eex2",sca_vTcPR_Eex2,NbMax);
// KKMC IFI off
  TH1D  *HTot_vTcPR_Ceex2  = HstProjV("HTot_vTcPR_Ceex2",sca_vTcPR_Ceex2,NbMax);
  TH1D  *HAfb_vTcPR_Ceex2  = HstProjA("HAfb_vTcPR_Ceex2",sca_vTcPR_Ceex2,NbMax);
// KKMC IFI off
  TH1D  *HTot_vTcPR_Ceex2n = HstProjV("HTot_vTcPR_Ceex2n",sca_vTcPR_Ceex2n,NbMax);
  TH1D  *HAfb_vTcPR_Ceex2n = HstProjA("HAfb_vTcPR_Ceex2n",sca_vTcPR_Ceex2n,NbMax);
 ///****************************************************************************************
  // More Wide range, vmax<1.
  TH2D *sca_vTcPL_Eex2   = (TH2D*)DiskFileA->Get("sca_vTcPL_Eex2");
  TH2D *sca_vTcPL_Ceex2  = (TH2D*)DiskFileA->Get("sca_vTcPL_Ceex2");
  TH2D *sca_vTcPL_Ceex2n = (TH2D*)DiskFileA->Get("sca_vTcPL_Ceex2n");
  //
  TH1D  *HTot_vTcPL_Ceex2 = HstProjV("HTot_vTcPL_Ceex2",sca_vTcPL_Ceex2,NbMax);
  TH1D  *HAfb_vTcPL_Ceex2 = HstProjA("HAfb_vTcPL_Ceex2",sca_vTcPL_Ceex2,NbMax);
  //
  TH1D  *HTot_vTcPL_Ceex2n = HstProjV("HTot_vTcPL_Ceex2n",sca_vTcPL_Ceex2n,NbMax);
  TH1D  *HAfb_vTcPL_Ceex2n = HstProjA("HAfb_vTcPL_Ceex2n",sca_vTcPL_Ceex2n,NbMax);
  //****************************************************************************************
  //  dsigma/dv unlimited cos(theta)
  //****************************************************************************************
  TH1D *Hpro_vT_Ceex2;
  ProjX1(sca_vTcPR_Ceex2, Hpro_vT_Ceex2);
  Hpro_vT_Ceex2->SetName("Hpro_vT_Ceex2");
  //  dsigma/dv unlimited cos(theta)
  TH1D *Hpro_vT_Ceex2n;
  ProjX1(sca_vTcPR_Ceex2n, Hpro_vT_Ceex2n);
  Hpro_vT_Ceex2n->SetName("Hpro_vT_Ceex2n");

cout<<"================ HisReMakeKKMC ENDs  ============================="<<endl;
cout<<"==================================================================="<<endl;

}//HisReMakeKKMC

///////////////////////////////////////////////////////////////////////////////////
void HisReMakeFoam35(TFile *DiskFileF){
  //------------------------------------------------------------------------
  int NbMax=0;
  //
  cout<<"==================================================================="<<endl;
  cout<<"================ ReMakeFoam35  BEGIN   ============================"<<endl;
//////////////////////////////////////////////////////////////////
  cout<<"  Renormalizing  and reprocessing histograms from FOAM"<<endl;

  TH1D *HST_FOAM_NORMA3 = (TH1D*)DiskFileF->Get("HST_FOAM_NORMA3");
  TH1D *HST_FOAM_NORMA5 = (TH1D*)DiskFileF->Get("HST_FOAM_NORMA5");

  TH1D *HST_xx_Ceex2n = (TH1D*)DiskFileF->Get("HST_xx_Ceex2n");  // FOAM
  TH1D *HST_tx_Ceex2n = (TH1D*)DiskFileF->Get("HST_tx_Ceex2n");  // FOAM testing norm.

  TH2D *SCA_xc_Ceex2n = (TH2D*)DiskFileF->Get("SCA_xc_Ceex2n");  // FOAM big   range x<0.99
  TH2D *SCT_xc_Ceex2n = (TH2D*)DiskFileF->Get("SCT_xc_Ceex2n");  // FOAM small range x<0.20

  TH2D *SCT_xc_EEX2   = (TH2D*)DiskFileF->Get("SCT_xc_EEX2");    // FOAM small range x<0.20
  TH2D *SCT_xc_EEX0   = (TH2D*)DiskFileF->Get("SCT_xc_EEX0");    // FOAM small range x<0.20

  TH2D *SCA_xc_Ceex2  = (TH2D*)DiskFileF->Get("SCA_xc_Ceex2");   // FOAM big   range x<0.99
  TH2D *SCT_xc_Ceex2  = (TH2D*)DiskFileF->Get("SCT_xc_Ceex2");   // FOAM small range x<0.20


  TH2D *SCA_xc_Ceex0  = (TH2D*)DiskFileF->Get("SCA_xc_Ceex0");   // FOAM big   range x<0.99
  TH2D *SCT_xc_Ceex0  = (TH2D*)DiskFileF->Get("SCT_xc_Ceex0");   // FOAM small range x<0.20

  TH2D *SCA_xc_Ceex0n  = (TH2D*)DiskFileF->Get("SCA_xc_Ceex0n");   // FOAM big   range x<0.99
  TH2D *SCT_xc_Ceex0n  = (TH2D*)DiskFileF->Get("SCT_xc_Ceex0n");   // FOAM small range x<0.20

////////////////////////////////////////////////////////////////
  // sigma(vmax) direct histogramming
  HisNorm1(HST_FOAM_NORMA3, HST_xx_Ceex2n );  // normalizing
  HisNorm1(HST_FOAM_NORMA3, HST_tx_Ceex2n );  // normalizing testing norm.
  //
  TH1D *HST_xmax_Ceex2n   = HstCumul("HST_xmax_Ceex2n",   HST_xx_Ceex2n);
  TH1D *HST_txmax_Ceex2n  = HstCumul("HST_txmax_Ceex2n",  HST_tx_Ceex2n);

////////////////////////////////////////////////////////////////
  // FOAM3
  HisNorm2(HST_FOAM_NORMA3, SCT_xc_EEX2 );    // normalizing
  HisNorm2(HST_FOAM_NORMA3, SCT_xc_EEX0 );    // normalizing
  // FOAM3
  HisNorm2(HST_FOAM_NORMA3, SCA_xc_Ceex2n );  // normalizing
  HisNorm2(HST_FOAM_NORMA3, SCT_xc_Ceex2n );  // normalizing
// FOAM3
  HisNorm2(HST_FOAM_NORMA3, SCA_xc_Ceex0n );  // normalizing
  HisNorm2(HST_FOAM_NORMA3, SCT_xc_Ceex0n );  // normalizing
// FOAM5
  HisNorm2(HST_FOAM_NORMA5, SCA_xc_Ceex2 );   // normalizing
  HisNorm2(HST_FOAM_NORMA5, SCT_xc_Ceex2 );   // normalizing
// FOAM5
  HisNorm2(HST_FOAM_NORMA5, SCA_xc_Ceex0 );   // normalizing
  HisNorm2(HST_FOAM_NORMA5, SCT_xc_Ceex0 );   // normalizing

///////////////////////////////////////////////////////////////
  // sigma(vmax) and AFB(vmax) from KKfoam, vmax<1, 50 bins in costheta
  int nbMax=0;   // <--- CosThetaMax = 1.0
  // Foam3
  TH1D  *Htot_xmax_Ceex2n = HstProjV("Htot_xmax_Ceex2n",SCA_xc_Ceex2n,NbMax);
  TH1D  *Hafb_xmax_Ceex2n = HstProjA("Hafb_xmax_Ceex2n",SCA_xc_Ceex2n,NbMax);
  // Foam5
  TH1D  *Htot_xmax_Ceex2 = HstProjV("Htot_xmax_Ceex2",SCA_xc_Ceex2,NbMax);
  TH1D  *Hafb_xmax_Ceex2 = HstProjA("Hafb_xmax_Ceex2",SCA_xc_Ceex2,NbMax);

// sigma(vmax) and AFB(vmax) from KKfoam scat. vmax<0.2, 100 bins in ctheta
//NbMax=45;           // cosThetaMax = 45/50=0.90 Now global variable
  // Foam3
  TH1D  *Htot2_xmax_EEX2 = HstProjV("Htot2_xmax_EEX2",SCT_xc_EEX2,NbMax);
  TH1D  *Hafb2_xmax_EEX2 = HstProjA("Hafb2_xmax_EEX2",SCT_xc_EEX2,NbMax);
  //
  TH1D  *Htot2_xmax_EEX0 = HstProjV("Htot2_xmax_EEX0",SCT_xc_EEX0,NbMax);
  TH1D  *Hafb2_xmax_EEX0 = HstProjA("Hafb2_xmax_EEX0",SCT_xc_EEX0,NbMax);
  // Foam3
  TH1D  *Htot2_xmax_Ceex0n = HstProjV("Htot2_xmax_Ceex0n",SCT_xc_Ceex0n,NbMax);
  TH1D  *Hafb2_xmax_Ceex0n = HstProjA("Hafb2_xmax_Ceex0n",SCT_xc_Ceex0n,NbMax);
  // Foam3
   TH1D  *Htot2_xmax_Ceex2n = HstProjV("Htot2_xmax_Ceex2n",SCT_xc_Ceex2n,NbMax);
   TH1D  *Hafb2_xmax_Ceex2n = HstProjA("Hafb2_xmax_Ceex2n",SCT_xc_Ceex2n,NbMax);
  // Foam5
  TH1D  *Htot2_xmax_Ceex0 = HstProjV("Htot2_xmax_Ceex0",SCT_xc_Ceex0,NbMax);
  TH1D  *Hafb2_xmax_Ceex0 = HstProjA("Hafb2_xmax_Ceex0",SCT_xc_Ceex0,NbMax);
  // Foam5
  TH1D  *Htot2_xmax_Ceex2 = HstProjV("Htot2_xmax_Ceex2",SCT_xc_Ceex2,NbMax);
  TH1D  *Hafb2_xmax_Ceex2 = HstProjA("Hafb2_xmax_Ceex2",SCT_xc_Ceex2,NbMax);

  cout<<"================ ReMakeFoam35 ENDs  ============================="<<endl;
  cout<<"==================================================================="<<endl;
}//ReMakeFoam35



