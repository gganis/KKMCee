*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//                                                                                 //
*//                         Pseudo-CLASS  GPS                                       //
*//           ***     will be  Pseudo-CLASS  CEEX    ***                            //
*//                                                                                 //
*//             Purpose:  Calculate CEEX spin amplitudes                            //
*//                                                                                 //
*//    WARNIG:  !!!!!!!!!!!!!!!!!!!!!!![[[[[[[[]]]]]]]!!!!!!!!!!!!!!!               //
*//    Function BVR_SBvirt modified temporarily to avoid problem with pi**2/beta    //
*//    Not important for anything in practice                                       //
*//    WARNIG:  !!!!!!!!!!!!!!!!!!!!!!![[[[[[[[]]]]]]]!!!!!!!!!!!!!!!               //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////

      SUBROUTINE GPS_Initialize
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//    Class initialization                                                         //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
*
      INCLUDE 'BXformat.h'
      INCLUDE 'GPS.h'
      INTEGER j1,j2,k,KeyPia
      INTEGER init
      DOUBLE PRECISION    alfa

      SAVE    init
      DATA init/0/
*------------------------------------
      IF(init .EQ. 1) RETURN
      init = 1
*------------------------------------
      m_out    = 16
*
*/////////////////////////////////////////////////////////////////////////////////////
*//                        Import from   BornV                                      //
*/////////////////////////////////////////////////////////////////////////////////////
* EW parameters
      CALL BornV_GetSwsq(  m_Sw2 )
      CALL BornV_GetGmu(   m_Gmu)
      CALL BornV_GetMZ(    m_MZ)
      CALL BornV_GetGammZ( m_GammZ)
      CALL BornV_GetAlfInv(m_AlfInv)
      alfa   = 1d0/m_AlfInv
      m_e_QED  = DSQRT( 4d0*m_pi*alfa)
      m_Alfpi  = alfa/m_pi
      CALL BornV_GetKeyElw(m_KeyElw)
      CALL BornV_GetKeyZet(m_KeyZet)
*
      CALL KK2f_GetKeyISR(m_KeyISR)
      CALL KK2f_GetKeyFSR(m_KeyFSR)
      CALL KK2f_GetKeyINT(m_KeyINT)
      CALL KK2f_GetKeyGPS(m_KeyGPS)
      CALL KK2f_GetVcut(  m_Vcut)
      m_HasFSR = m_KeyFSR ! This might be redefined later!!!
*
*/////////////////////////////////////////////////////////////////////////////////////
      WRITE(m_out,bxope)
      WRITE(m_out,bxtxt) '  GPS   Initializator                '
      WRITE(m_out,bxl1f) m_MZ    ,   'Z mass     [GeV]   ','MZ    ','a1'
      WRITE(m_out,bxl1f) m_GammZ,    'Z width    [GeV]   ','GammZ ','a2'
      WRITE(m_out,bxl1f) m_Sw2,      'sin(theta_w)**2    ','Sw2   ','a3'
      WRITE(m_out,bxl1f) m_AlfInv,   '1/alfa_QED  at  Q=0','AlfInv','a4'
      WRITE(m_out,bxtxt) 'Test switches:                         '
      WRITE(m_out,bxl1i) m_KeyZet,   'Z on/off   switch  ','KeyZet','a5'
      WRITE(m_out,bxl1i) m_KeyElw,   'Electroweak lib.   ','KeyElw','a6'
      WRITE(m_out,bxl1i) m_KeyGPS,   'CEEX level         ','KeyGPS','a7'
      WRITE(m_out,bxl1i) m_KeyISR,   'ISR emission       ','KeyISR','a8'
      WRITE(m_out,bxl1i) m_KeyFSR,   'FSR emission       ','KeyFSR','a9'
      WRITE(m_out,bxl1i) m_KeyINT,   'ISR*FSR interferenc','KeyINT','a10'
      WRITE(m_out,bxclo)
*/////////////////////////////////////////////////////////////////////////////////////
* Key for switching on/off the use of m_b, can be reset with GPS_SetKeyArb
      m_KeyArb = 0  ! default zero value, m_b=Xi is assumed
*
      WRITE(m_out,bxope)
      WRITE(m_out,bxtxt) '     Initialization of GPS class      '
      WRITE(m_out,bxclo)
********************************
*     Define Pauli matrices
      DO k = 0,3
         DO j1 = 1,2
            DO j2 = 1,2
               m_Pauli( k,j1,j2) = DCMPLX(0d0,0d0)
            ENDDO
         ENDDO
      ENDDO
* Sigma0
      m_Pauli( 0,1,1) = DCMPLX( 1d0, 0d0)
      m_Pauli( 0,2,2) = DCMPLX( 1d0, 0d0)
* SigmaX
      m_Pauli( 1,1,2) = DCMPLX( 1d0, 0d0)
      m_Pauli( 1,2,1) = DCMPLX( 1d0, 0d0)
* SigmaY
      m_Pauli( 2,1,2) = DCMPLX( 0d0,-1d0)
      m_Pauli( 2,2,1) = DCMPLX( 0d0, 1d0)
* SigmaZ
      m_Pauli( 3,1,1) = DCMPLX( 1d0, 0d0)
      m_Pauli( 3,2,2) = DCMPLX(-1d0, 0d0)
*
* The other notation for 4-vector index
      DO k = 1,3
         DO j1 = 1,2
            DO j2 = 1,2
               m_Pauli4( k,j1,j2) = m_Pauli( k,j1,j2)
            ENDDO
         ENDDO
      ENDDO
      DO j1 = 1,2
         DO j2 = 1,2
            m_Pauli4( 4,j1,j2) = m_Pauli( 0,j1,j2)
         ENDDO
      ENDDO
********************************
*     Define GPS vectors in CMS
*     Note that they define uniquely the inner products
*     and spin quantization axises in fermion rest frame
      DO k = 1,4
         m_Xi(k)  = 0d0
         m_Eta(k) = 0d0
      ENDDO
      m_Xi( 4)  =  1d0
      m_Xi( 1)  =  1d0
      m_Eta(2)  =  1d0
* axial vectors (arbitrary lightlike ) for photon polarization
      m_b1( 1)   =  0.8723d0
      m_b1( 2)   = -0.7683d0
      m_b1( 3)   =  0.3348d0
      m_b1( 4)   =  DSQRT(m_b1( 1)**2 + m_b1( 2)**2 +m_b1( 3)**2)
* another random setting
      m_b2( 1)   = -0.78833d0
      m_b2( 2)   =  0.34788d0
      m_b2( 3)   = -0.33282d0
      m_b2( 4)   =  DSQRT(m_b2( 1)**2 + m_b2( 2)**2 +m_b2( 3)**2)
* this setting is very close to b=Xi, for special tests
      m_b3( 1)   =  1d0
      m_b3( 2)   =  1d-7
      m_b3( 3)   =  0d0
      m_b3( 4)   =  DSQRT(m_b3( 1)**2 + m_b3( 2)**2 +m_b3( 3)**2)
*
* default setting to m_b1, can be changed with CALL GPS_Setb2,3
      DO k=1,4
         m_b( k) = m_b1( k)
      ENDDO
*
      CALL KarFin_GetKeyPia(KeyPia)
      IF( (m_KeyFSR .EQ. 1) .AND.  (KeyPia .EQ. 0) ) THEN
         WRITE(*,*) ' +++++ STOP in GPS: you cannot have KeyPia=0 for GPS and FSR '
         STOP
      ENDIF
*
      END                       !!!end of GPS_Initialize!!!


      SUBROUTINE GPS_Make
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   Spin amplitudes O(alf1) and O(alf0) for exponentiation.                       //
*//   Helicities of emitted photons are choosen randomly!!!                         //
*//                                                                                 //
*//   INPUT: is transfered from  BornV and KK2f using Getters                       //
*//                                                                                 //
*//   OUTPUT:                                                                       //
*//   m_AmpExpo0 =  O(alf0) spin amplitudes                                         //
*//   m_AmpExpo1 =  O(alf1) spin amplitudes                                         //
*//   m_AmpExpo2 =  O(alf2) spin amplitudes                                         //
*//   m_AmpExpo2p=  O(alf2) spin amplitudes with virtual pairs                      //
*//   Wt0     =  O(alf0) exponentiation weight                                      //
*//   Wt1     =  O(alf1) exponentiation weight                                      //
*//                                                                                 //
*//   COMMON working space:                                                         //
*//   m_AmpExpo*  is working space, used by HiniPlus, HfinPlus, HfinMinus           //
*//   m_AmpBorn   is working space, used by HfinMinus                               //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'GPS.h'
*
      INTEGER               KFi,KFf,Nphot
      DOUBLE PRECISION      p1(4),p2(4),p3(4),p4(4),Phot(m_phmax,4)
*
      DOUBLE PRECISION      ph(4),ph1(4),ph2(4),PX(4),PP(4),QQ(4)
      INTEGER               i,j,k,l,n,j1,j2,last,loop,loop2,Hel,Hel1,Hel2
      DOUBLE PRECISION      ChaIni,ChaFin
      DOUBLE PRECISION      svar,svarX,svarX1,svarQ,Ene,betaf
      DOUBLE COMPLEX        GPS_soft,GPS_softb
      DOUBLE COMPLEX        Sini(2,m_phmax),Sfin(2,m_phmax)
      DOUBLE COMPLEX        CKine,sProd,Sactu,Sactu1,Sactu2,Cfact0,Cfact2,SactuA,SactuB
      DOUBLE PRECISION      Fleps,Massf,Mbeam,m1,m2,m3,m4,mph
      DOUBLE PRECISION      XborSum,XboxSum,CrudSum,Exp0Sum,Exp1Sum
      DOUBLE PRECISION      Xborn,Xboxy
      DOUBLE PRECISION      DistCru,BornCru,fLLux
      DOUBLE PRECISION      CrudNorm, ExpoNorm
      DOUBLE PRECISION      RhoCru3
      DOUBLE PRECISION      alfQED,alfpini,alfpfin, alfpmix, Emin, MasPhot
      DOUBLE PRECISION      BVR_SForFac,BVR_TForFac
      DOUBLE PRECISION      Yisr,Yfsr,Yint
      DOUBLE PRECISION      YFSkonIni, YFSkonFin, YFS_IRini, YFS_IRfin, YFS_isr, YFS_fsr
      DOUBLE PRECISION      BornV_GetMass, BornV_GetCharge, BornV_Simple,BornV_Differential
      DOUBLE PRECISION      Wt0,Wt1
      DOUBLE PRECISION      dummy
      INTEGER               BornV_GetColor, NCf
      INTEGER               iSaveCPU
*-----------------------------------------------
      DOUBLE PRECISION   Icont
      SAVE      Icont
      DATA      Icont /0d0/
*-----------------------------------------------
      CALL GPS_Initialize
*
      CALL KK2f_GetKFini(   KFi)     ! Normaly for beam KFi=11 is electron
      CALL KarLud_GetKFfin( KFf)     ! Actual KFcode of the final fermion
      Mbeam  =  BornV_GetMass(   KFi)
      Massf  =  BornV_GetMass(   KFf)
      ChaIni =  BornV_GetCharge( KFi)
      ChaFin =  BornV_GetCharge( KFf)
      NCf    =  BornV_GetColor(  KFf)
* check for existence of FSR
      CALL KarFin_GetIsFSR(m_HasFSR)
*

      Fleps =  1d-100
      m1  = Mbeam
      m2  = Mbeam
      m3  = Massf
      m4  = Massf
      mph = Fleps
      CALL KK2f_GetBeams(    p1,p2)
      CALL KK2f_GetFermions( p3,p4)
      CALL KK2f_GetPhotAll(Nphot,Phot) ! ordered in energy
      CALL KK2f_GetEmin( Emin)
      CALL KK2f_GetMasPhot(MasPhot)
*-------------------------------------------------------------
      DO k=1,4
         PP(k) = p1(k)+p2(k)
         QQ(k) = p3(k)+p4(k)
      ENDDO
      svar  = PP(4)**2 - PP(3)**2 - PP(2)**2 - PP(1)**2
      svarQ = QQ(4)**2 - QQ(3)**2 - QQ(2)**2 - QQ(1)**2
      Ene   = SQRT(svar)/2d0
*-------------------------------------------------------------
* Overall normalization factors
      CrudNorm  =  1d0
      ExpoNorm  =  2d0/(4d0*m_pi)*NCf  ! it is still quasi-empirical...
*////////////////////////////////////////////////////////////////////////////////////
*//                         YFS  FormFactors                                       //
*//  Note that FSR formfactor below cannot be used for KeyPia=0, Emin is in CMS!!! //
*////////////////////////////////////////////////////////////////////////////////////
      alfQED   = m_e_QED**2/(4d0*m_pi)
      alfpini  = m_Alfpi*ChaIni**2
      alfpfin  = m_Alfpi*ChaFin**2
      alfpmix  = m_Alfpi*ChaFin*ChaIni
      IF( m_KeyISR .NE. 0) THEN
         CALL  BornV_GetYFS_IR( YFS_IRini )
         CALL  BornV_GetYFSkon( YFSkonIni )
****>>   YFS_isr =  YFS_IRini*YFSkonIni
         YFS_isr =  YFS_IRini       !!!<<**  YFSkon not in WtCrud (historical reasons)
         CrudNorm  = CrudNorm *YFS_isr
         Yisr= BVR_SForFac( alfpini, p1,Mbeam, p2,Mbeam, Emin, MasPhot)
         ExpoNorm  = ExpoNorm *Yisr
      ENDIF
      IF( m_HasFSR .NE. 0) THEN
         CALL KarFin_GetYFS_IR( YFS_IRfin )
         CALL KarFin_GetYFSkon( YFSkonFin )
****>>   YFS_fsr =  YFS_IRfin*YFSkonFin
         YFS_fsr =  YFS_IRfin       !!!<<** YFSkon not in WtCrud (historical reasons)
         CrudNorm  = CrudNorm *YFS_fsr
         Yfsr= BVR_SForFac( alfpfin, p3,Massf, p4,Massf, Emin, MasPhot)
         ExpoNorm  = ExpoNorm *Yfsr
      ENDIF
* Remember Yint depends on Emin and provides angular asymmetry (MasPhot is dummy)
      IF(  m_KeyINT .NE. 0 .AND. m_KeyISR .NE. 0 .AND.  m_HasFSR .NE. 0  ) THEN
         Yint= BVR_TForFac( alfpmix, p1,Mbeam, p3,Massf, Emin, MasPhot)
     $        *BVR_TForFac( alfpmix, p2,Mbeam, p4,Massf, Emin, MasPhot)
     $        *BVR_TForFac(-alfpmix, p1,Mbeam, p4,Massf, Emin, MasPhot)
     $        *BVR_TForFac(-alfpmix, p2,Mbeam, p3,Massf, Emin, MasPhot)
         ExpoNorm  = ExpoNorm *Yint
      ENDIF
*//////////////////////////////////////////////////////////////////
*//                       S-factors                              //
*//////////////////////////////////////////////////////////////////
* List of soft-factors, note that we calculate them for helicity=+1
* The other one helicity=-1 is just minus complex conjugate!
* Sig = 3-2*Hel, Hel=1,2 --> Sig=1,-1
      DO j=1,nphot
         CrudNorm = CrudNorm *1d0/(2d0*m_pi)**3    !!<-- photon phase-space factor
         ExpoNorm = ExpoNorm *1d0/(2d0*m_pi)**3    !!<-- photon phase-space factor
         DO k=1,4
            ph(k) = Phot(j,k)
         ENDDO
         IF( m_KeyArb .EQ. 0 ) THEN
            Sini(1,j)  =  DCMPLX(ChaIni*m_e_QED) *GPS_soft(  1,ph,p1,p2)
            Sfin(1,j)  = -DCMPLX(ChaFin*m_e_QED) *GPS_soft(  1,ph,p3,p4)
         ELSE
            Sini(1,j)  =  DCMPLX(ChaIni*m_e_QED) *GPS_softb( 1,ph,p1,m1,p2,m2)
            Sfin(1,j)  = -DCMPLX(ChaFin*m_e_QED) *GPS_softb( 1,ph,p3,m3,p4,m4)
         ENDIF
         Sini(2,j) = -DCONJG(Sini(1,j))
         Sfin(2,j) = -DCONJG(Sfin(1,j))
      ENDDO
*//////////////////////////////////////////////////////////////////
*//             Define (randomly) photon helicities              //
*//////////////////////////////////////////////////////////////////
* Photon gelicities are now generated in the main module KK2f and imported here
* This is better in case of multiple calls of GPS for the single event
      CALL KK2f_GetPhel(m_Phel)
***   CALL GPS_PhelRandom(nphot)
***   WRITE(*,'(a,20i2)') 'Phel=   ',(m_Phel(j),j=1,nphot)
      CALL GPS_BornZero(m_AmpExpo0)
      CALL GPS_BornZero(m_AmpExpo1)
      CALL GPS_BornZero(m_AmpExpo2)
      CALL GPS_BornZero(m_AmpExpo2p)
      CrudSum = 0d0
      XborSum = 0d0
      XboxSum = 0d0
      iSaveCPU = 0
*//////////////////////////////////////////////////////////////////////
*//  *************************************************************   //
*//              LOOP OVER PARTITIONS STARTS HERE                    //
*//  Initialize loop over partitions, m_isr(j)=1,0 denotes isr,fsr   //
*//  *************************************************************   //
*//////////////////////////////////////////////////////////////////////
      CALL GPS_PartitionStart(nphot,last)
      DO loop=1,10000000
*        Calculate reduced 4-momentum to be used in gamma/Z propagator
         DO k=1,4
            PX(k) = p1(k)+p2(k)
         ENDDO
*   /////////////////////////////////////////////////////////
*   //         ===============================             //
*   //         Soft-part, soft-factors, m-zero             //
*   //         ===============================             //
*   /////////////////////////////////////////////////////////
         sProd=DCMPLX(1d0,0d0)
         DO j=1,nphot
            Hel  = m_Phel(j)+1
            Sactu  = m_isr(j)*Sini(Hel,j) + (1-m_isr(j))*Sfin(Hel,j)
            sProd  = sProd*Sactu
            DO k=1,4
               PX(k) = PX(k) -m_isr(j)*Phot(j,k)
            ENDDO
         ENDDO
         svarX = PX(4)**2 - PX(3)**2 - PX(2)**2 - PX(1)**2
         Cfact0 = sProd  *(svarX/svarQ)
         CALL GPS_BornPlus(iSaveCPU,KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4, Cfact0,Xborn,Xboxy)
*        ==================================================================================
*****    IF(nphot.EQ.0) CALL GPS_BPrint(6,' AmpBoxy ',m_AmpBoxy)
*****    IF(nphot.EQ.0) WRITE(*,'(a,5g20.10)') '#######>>> Boxy/Born = ',Xborn
*****    CALL GPS_BPrint(6,' AmpExp0',m_AmpExpo0)
*****    CALL GPS_BPrint(6,' AmpBoxy',m_AmpBoxy)
*****    WRITE(*,'(a,5g20.10)') '#######>>> Boxy/Born = ',Xborn
*  //////////////////////////////////////////////////////////////////////////////////////
*  //    Crude distribution as in low level Monte Carlo, see also KK2f_DsigOverDtau    //
*  //    We DO NOT use BornCru memorized in KK2f in order to be able to recalculate    //
*  //    Weights with another input parameters like MZ etc.                            //
*  //    This trick is susspended, have to think if it could be restored?              //
*  //////////////////////////////////////////////////////////////////////////////////////
*****    BornCru = 4d0/3d0*BornV_Simple( KFi,KFf,svarX, 0d0  )
         BornCru = 4d0/3d0*BornV_Differential(0,KFf,svarX,0d0,0d0,0d0,0d0,0d0)
         BornCru = BornCru *(svar/svarX)                                  !!<-- Born(svar)*svar
         fLLux   = svarX/svarQ                                            !!<-- extra LL factor
         betaf = DSQRT( 1d0 - 4*Massf**2/svarQ )                          !!<-- 2-body phase spase
         DistCru = BornCru/(4d0*m_pi) *fLLux *2d0/betaf *CDABS(sProd)**2  !!<-- CRUDE
         CrudSum =       CrudSum  + DistCru
         XborSum =       XborSum  + XBorn
         XboxSum =       XboxSum  + XBoxy
***      WRITE(*,'(a,5g20.10)') ' vv, Xborn/DistCru ', 1-svarX/svar, Xborn/DistCru
*   /////////////////////////////////////////////////////////////////
*   //         1-photon IR-finite contribution     ms-one           //
*   /////////////////////////////////////////////////////////////////
         iSaveCPU = 1
         IF(nphot .EQ. 0) GOTO 300
*        =========================
c[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      CALL GPS_BornZero(m_AmpExpo2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
         DO j1=1,nphot
            DO k=1,4
               ph1(k) = Phot(j1,k)
            ENDDO
            IF( ph1(4)/Ene .GT. m_Vcut(1) ) THEN ! accept 1 hard only
               Hel1  = m_Phel(j1)+1
               Sactu1  = m_isr(j1)*Sini(Hel1,j1) + (1-m_isr(j1))*Sfin(Hel1,j1)
               SvarX1  = (QQ(4)+ph1(4))**2-(QQ(3)+ph1(3))**2-(QQ(2)+ph1(2))**2-(QQ(1)+ph1(1))**2 !
               CKine   = (svarX1/svarQ)
               IF( m_isr(j1) .EQ. 1) THEN
                  CALL GPS_HiniPlus(KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,ph1,mph,Hel1,Sactu1,sProd) !
                  CALL GPS_HiniPlusW(1,KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,ph1,mph,Hel1,Sactu1,sProd) ! One calculates part due to W for neutrinos for 1 subtract for -1, fixed transfer approx for 0.
C[[[[[[[[[[
C                  IF (NPHOT.EQ.2.and.phot(1,4)/ene.gt.m_Vcut(2).and.phot(2,4)/ene.gt.m_Vcut(2)) then
C???                   CALL GPS_HiniPlusW(1,KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,ph1,mph,Hel1,Sactu1,sProd) ! check???
C                  ENDIF
C]]]]]]]]]]
               ELSE
                  CALL GPS_HfinPlus(KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,ph1,mph,Hel1,Sactu1,sProd,CKine) !
               ENDIF
            ENDIF
         ENDDO
*   /////////////////////////////////////////////////////////////////
*   //         2-photon IR-finite contribution   ms-two            //
*   /////////////////////////////////////////////////////////////////
*   O(alf2) correction to 2 photons in ISR or FSR, 
*   implicit Bose Einstein symmetrization inside HiiPlus, HffPlus
         DO j1=1,nphot
            DO j2=j1+1,nphot
               DO k=1,4
                  ph1(k) = Phot(j1,k)
                  ph2(k) = Phot(j2,k)
               ENDDO
               IF( (ph1(4)/Ene .GT. m_Vcut(2)) .AND. (ph2(4)/Ene .GT. m_Vcut(2)) ) THEN ! accept 2 hard only
                  Hel1  = m_Phel(j1)+1
                  Hel2  = m_Phel(j2)+1
                  Sactu2  = ( m_isr(j1)*Sini(Hel1,j1) + (1-m_isr(j1))*Sfin(Hel1,j1) ) !
     $                     *( m_isr(j2)*Sini(Hel2,j2) + (1-m_isr(j2))*Sfin(Hel2,j2) ) !
                  Cfact2 = sProd/Sactu2
                  IF(     (m_isr(j1) .EQ. 1) .AND. (m_isr(j2) .EQ. 1) ) THEN ! ini-ini
                     CALL GPS_HiiPlus(Cfact2,KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,Hel1,ph1,Hel2,ph2,mph,m_AmpExpo2) !
                     CALL GPS_HiiPlus(Cfact2,KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,Hel1,ph1,Hel2,ph2,mph,m_AmpExpo2p) !

                     SactuA  = ( m_isr(j1)*Sini(Hel1,j1) + (1-m_isr(j1))*Sfin(Hel1,j1) ) !
                     SactuB  = ( m_isr(j2)*Sini(Hel2,j2) + (1-m_isr(j2))*Sfin(Hel2,j2) ) !
                     CALL GPS_HiniPlusW(-1,KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,ph1,mph,Hel1,SactuA,sProd) !
                     CALL GPS_HiniPlusW(-1,KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,ph2,mph,Hel2,SactuB,sProd) !

                     CALL GPS_HiiPlusW(0,Cfact2,KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,Hel1,ph1,Hel2,ph2,mph,m_AmpExpo2)
c[[[???                     CALL GPS_HiiPlusW(0,Cfact2,KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,Hel1,ph1,Hel2,ph2,mph,m_AmpExpo2p)
                  ELSEIF( (m_isr(j1) .EQ. 0) .AND. (m_isr(j2) .EQ. 0) ) THEN ! fin-fin
                     CALL GPS_HffPlus(Cfact2,KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,Hel1,ph1,Hel2,ph2,mph,m_AmpExpo2) !
                     CALL GPS_HffPlus(Cfact2,KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,Hel1,ph1,Hel2,ph2,mph,m_AmpExpo2p) !
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
*   O(alf2) correction to 1 photon ISR and 1 photon FSR, explicit Bose-Einstein symmetrisation
         DO j1=1,nphot
            DO j2=1,nphot
               IF( j1 .NE. j2) THEN
                  DO k=1,4
                     ph1(k) = Phot(j1,k)
                     ph2(k) = Phot(j2,k)
                  ENDDO
                  IF( (ph1(4)/Ene .GT. m_Vcut(2)) .AND. (ph2(4)/Ene .GT. m_Vcut(2)) ) THEN ! accept 2 hard only
                     Hel1  = m_Phel(j1)+1
                     Hel2  = m_Phel(j2)+1
                     Sactu2  = ( m_isr(j1)*Sini(Hel1,j1) + (1-m_isr(j1))*Sfin(Hel1,j1) ) !
     $                        *( m_isr(j2)*Sini(Hel2,j2) + (1-m_isr(j2))*Sfin(Hel2,j2) ) !
                     Cfact2 = sProd/Sactu2
                     IF(     (m_isr(j1) .EQ. 1) .AND. (m_isr(j2) .EQ. 0) ) THEN ! ini-fin
                        CALL GPS_HifPlus(Cfact2,KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,Hel1,ph1,Hel2,ph2,mph,m_AmpExpo2) !
                        CALL GPS_HifPlus(Cfact2,KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,Hel1,ph1,Hel2,ph2,mph,m_AmpExpo2p) !
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
*-------------------
*****    WRITE(*,'(f20.10,a,20i2)') sqrt(svarX),' >>> m_isr(j)= ',(m_isr(j),j=1,nphot)
*/////////////////////////////////////////////////////////
*//   Update m_isr, check if it is the last partition   //
*/////////////////////////////////////////////////////////
         IF(last .EQ. 1) GOTO 300
         CALL GPS_PartitionPlus(nphot,last)
         IF(last .EQ. 2) GOTO 300
      ENDDO
      write(*,*) '########### INCORRECT EXIT from loop over partitions'
      STOP
 300  CONTINUE
*//////////////////////////////////////////////////////////////////////
*//  *************************************************************   //
*//              LOOP OVER PARTITIONS ENDS HERE                      //
*//  *************************************************************   //
*//////////////////////////////////////////////////////////////////////

 
      m_RhoCrud = CrudSum *CrudNorm   !<-- Crude (unpolarized for the time being)
      CALL GPS_MakeRho(ExpoNorm)      !<-- Defines m_RhoExp0, m_RhoExp1, m_RhoExp2

*/////////////////////////////////////////////////////////////////////////////////
*//                     Assigning weight list                                   //
*/////////////////////////////////////////////////////////////////////////////////
      IF( m_KeyINT .EQ. 0 ) THEN
         m_WtSet( 51) =   m_RhoExp0 /m_RhoCrud    !!! Interference OFF
         m_WtSet( 52) =   m_RhoExp1 /m_RhoCrud
         m_WtSet( 53) =   m_RhoExp2 /m_RhoCrud
         m_WtSet( 63) =   m_RhoExp2p/m_RhoCrud    !!! Pairs ON
         m_WtBest     =   m_WtSet( 53)
c[[[[[[[[[[[[[[[[[[[[
c      write(m_out,*) '>>>GPS_Make:  CrudSum =',CrudSum,'   CrudNorm=',CrudNorm
c      write(m_out,*) '>>>GPS_Make:  m_RhoCrud =',m_RhoCrud
c      write(m_out,*) '>>>GPS_Make:  m_WtSet(51-53)=', m_WtSet(51), m_WtSet(52), m_WtSet(53)
c]]]]]]]]]]]]]]]]]]]]
      ELSE
         m_WtSet( 1)  =   m_RhoExp0 /m_RhoCrud    !!! Interference ON
         m_WtSet( 2)  =   m_RhoExp1 /m_RhoCrud
         m_WtSet( 3)  =   m_RhoExp2 /m_RhoCrud
         m_WtSet(13)  =   m_RhoExp2p/m_RhoCrud    !!! Pairs ON
         m_WtBest     =   m_WtSet( 3)
c[[[[[[[[[[[[[[[[[[[[
c      write(m_out,*) '>>>GPS_Make:  m_WtSet( 1-3 )=', m_WtSet( 1), m_WtSet( 2), m_WtSet( 3)
c]]]]]]]]]]]]]]]]]]]]
      ENDIF
*/////////////////////////////////////////////////////////////////////////////////
*//                     X-Checks   Printouts                                    //
*/////////////////////////////////////////////////////////////////////////////////
* debug variables for tests
      m_debg( 3) = XBorSum
      m_debg( 4) = XBoxSum
*
      IF(Icont .GE. 100 )  RETURN
***      WRITE(*,'(3a)') '**************************************************'
***     $                            ,' !!!GPS_Make!!! ',
***     $                '**************************************************'
***      WRITE( *,'(a,i10,6g20.12)') 'GPS_Make::  nphot, vq, wt0,wt1,wt2',
***     $     nphot, 1-svarQ/svar, m_WtSet( 1),m_WtSet( 2),m_WtSet( 3)
***      CALL GPS_BPrint(6,'m_AmpEx0',m_AmpExpo0)
***      CALL GPS_BPrint(6,'m_AmpEx1',m_AmpExpo1)
***      CALL GPS_BPrint(6,'m_AmpEx2',m_AmpExpo2)
***   IF(Wt1 .LT. 10d0 )  RETURN
***   IF( (1-svarQ/svar) .LT. 0.1d0 ) RETURN
***   IF( nphot.NE.0 ) RETURN
***   Wt0 = Exp0Sum/XborSum/2**nphot
***   Wt1 = Exp1Sum/XborSum/2**nphot
***   CALL KK2f_Print1(6)
* X-check with DsigOverDtau makes sense ONLY for KeyINT=0
***   CALL  KK2f_DsigOverDtau(6,RhoCru3)
***   WRITE(*,'(a,5g20.12)') '***  nphot, vq, wt0,wt1', nphot, 1-svarQ/svar,wt0,wt1
***   WRITE(*,'(a,5g20.12)') '*** CDABS(sProd)**2 = ', CDABS(sProd)**2 *CrudNorm
***   WRITE(*,'(a,5g20.12)') '*** m_RhoCrud/RhoCru3 = ', m_RhoCrud/RhoCru3, m_RhoCrud

      Icont = Icont +1
      END                       !!!end of GPS_Make!!!


      DOUBLE PRECISION  FUNCTION GPS_MakeRhoFoam(XNorm)
*//////////////////////////////////////////////////////////////////////////////////
*//                                                                              //
*//   Calculate differential distributions (normalized to LIPS) from spin ampl.  //
*//                UNPOLARIZED final fermions                                    //
*//                                                                              //
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'GPS.h'
      DOUBLE PRECISION   XNorm
      DOUBLE PRECISION   Sum0
      INTEGER            j1,j2,j3,j4

      Sum0 = 0d0
      DO j1 = 1,2
         DO j2 = 1,2
            DO j3 = 1,2
               DO j4 = 1,2
                   Sum0 = Sum0 +DREAL( m_AmpBorn2(j1,j2,j3,j4)*DCONJG(m_AmpBorn1(j1,j2,j3,j4)) )
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      GPS_MakeRhoFoam = Sum0*XNorm
      END

      SUBROUTINE GPS_MakeRho(ExpoNorm)
*//////////////////////////////////////////////////////////////////////////////////
*//                                                                              //
*//   Calculate differential distributions (normalized to LIPS) from spin ampl.  //
*//                UNPOLARIZED final fermions                                    //
*//                                                                              //
*//   To be done:                                                                //
*//   One needs to put wigner rotation for initial polarizations somewhere       //
*//   either in setter called in KK2f, or in flight, for every event             //
*//                                                                              //
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'GPS.h'
* params
      DOUBLE PRECISION      ExpoNorm
* locals
      INTEGER               j1,j2,j3,j4,i1,i2,i3,i4
      DOUBLE PRECISION      Sum0, Sum1, Sum2, Sum2p
      DOUBLE PRECISION      polar1,  polar2
      DOUBLE COMPLEX        SDMprod, Tensor0, Tensor1, Tensor2, Tensor2p
*----------------------------------
* Normalized factor to LIPS for CEEX amplitudes is memorized for further use
      m_ExpoNorm = ExpoNorm
*
      polar1 = SQRT(ABS(m_PolBeam1(3)**2 +m_PolBeam1(2)**2 +m_PolBeam1(1)**2)) !
      polar2 = SQRT(ABS(m_PolBeam2(3)**2 +m_PolBeam2(2)**2 +m_PolBeam2(1)**2)) !
      Sum0 = 0d0
      Sum1 = 0d0
      Sum2 = 0d0
      Sum2p = 0d0
      IF( (polar1+polar2) .LT. 1d-6) THEN
*        The case with UNPOLARIZED beams and UNPOLARIZED final fermions
         DO j1 = 1,2
            DO j2 = 1,2
               DO j3 = 1,2
                  DO j4 = 1,2
                     Sum0 = Sum0 +m_AmpExpo0(j1,j2,j3,j4)*DCONJG(m_AmpExpo0(j1,j2,j3,j4)) !
                     Sum1 = Sum1 +m_AmpExpo1(j1,j2,j3,j4)*DCONJG(m_AmpExpo1(j1,j2,j3,j4)) !
                     Sum2 = Sum2 +m_AmpExpo2(j1,j2,j3,j4)*DCONJG(m_AmpExpo2(j1,j2,j3,j4)) !
                     Sum2p= Sum2p +m_AmpExpo2p(j1,j2,j3,j4)*DCONJG(m_AmpExpo2p(j1,j2,j3,j4)) !
                  ENDDO                  
               ENDDO
            ENDDO
         ENDDO
      ELSE
*        The case with POLARIZED beams and UNPOLARIZED final fermions
         DO j1=1,2
         DO i1=1,2
            DO j2=1,2
            DO i2=1,2
               DO j3=1,2
                  DO j4=1,2
                     SDMprod = m_SDMat1(i1,j1)*m_SDMat2(i2,j2)
                     Tensor0 = m_AmpExpo0(i1,i2,j3,j4)*DCONJG(m_AmpExpo0(j1,j2,j3,j4)) !
                     Tensor1 = m_AmpExpo1(i1,i2,j3,j4)*DCONJG(m_AmpExpo1(j1,j2,j3,j4)) !
                     Tensor2 = m_AmpExpo2(i1,i2,j3,j4)*DCONJG(m_AmpExpo2(j1,j2,j3,j4)) !
                     Tensor2p = m_AmpExpo2p(i1,i2,j3,j4)*DCONJG(m_AmpExpo2p(j1,j2,j3,j4)) !
                     Sum0 = Sum0 + Tensor0 *SDMprod
                     Sum1 = Sum1 + Tensor1 *SDMprod
                     Sum2 = Sum2 + Tensor2 *SDMprod
                     Sum2p = Sum2p + Tensor2p *SDMprod
                  ENDDO
               ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDIF
      m_RhoExp0 = Sum0 *ExpoNorm
      m_RhoExp1 = Sum1 *ExpoNorm
      m_RhoExp2 = Sum2 *ExpoNorm
      m_RhoExp2p = Sum2p *ExpoNorm
c[[[[[[[[[[[[[[[[[[[[[[[[[[[ *debug*
      write(16,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
c      write(16,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@GPS_MakeRho@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
c      write(16,*) 'ExpoNorm=',ExpoNorm
c      write(16,'(a,i1,a,i1,a,8g22.11)') (('m_AmpExpo0(',j1,',',j2,',*,*)=',((m_AmpExpo0(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c      write(16,*) '---------------------------------------------------------------------------------------------------------------'
c      write(16,'(a,i1,a,i1,a,8g22.11)') (('m_AmpExpo1(',j1,',',j2,',*,*)=',((m_AmpExpo1(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c      write(16,*) '---------------------------------------------------------------------------------------------------------------'
c      write(16,'(a,i1,a,i1,a,8g22.11)') (('m_AmpExpo2(',j1,',',j2,',*,*)=',((m_AmpExpo2(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
      write(16,*) '@@@@@@@@@ m_RhoExp0 =',m_RhoExp0,'  m_RhoExp1 =',m_RhoExp1,'  m_RhoExp2 =',m_RhoExp2
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
      END                       !!!GPS_MakeRho!!!

      SUBROUTINE GPS_MakeRho2(wt0,wt1,wt2)
*//////////////////////////////////////////////////////////////////////////////////
*//                                                                              //
*//   Used in Taupair_ImprintSpin                                                //
*//                                                                              //
*//   Calculate differential distributions (normalized to LIPS) from spin ampl.  //
*//                  POLARIZED final fermions                                    //
*//                                                                              //
*//                                                                              //
*//////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'GPS.h'
* params
      DOUBLE PRECISION   wt0,wt1,wt2
* locals
      INTEGER      j1,j2,j3,j4,i1,i2,i3,i4
      DOUBLE PRECISION        Sum0, Sum1, Sum2
      DOUBLE PRECISION        polar1,  polar2
      DOUBLE COMPLEX          Tensor0, Tensor1, Tensor2, SDMprod
      DOUBLE PRECISION        Rho0, Rho1, Rho2
*----------------------------------
      polar1 = SQRT(ABS(m_PolBeam1(3)**2 +m_PolBeam1(2)**2 +m_PolBeam1(1)**2))
      polar2 = SQRT(ABS(m_PolBeam2(3)**2 +m_PolBeam2(2)**2 +m_PolBeam2(1)**2))
      Sum0 = 0d0
      Sum1 = 0d0
      Sum2 = 0d0
      IF( (polar1+polar2) .LT. 1d-6) THEN
*        The case with UNPOLARIZED beams and POLARIZED final fermions
         DO j1 = 1,2
            DO j2 = 1,2
               DO j3 = 1,2
               DO i3 = 1,2
                     DO j4 = 1,2
                     DO i4 = 1,2
                        SDMprod = m_SDMat3(j3,i3)*m_SDMat4(j4,i4)
                        Tensor0 = m_AmpExpo0(j1,j2,i3,i4)*DCONJG(m_AmpExpo0(j1,j2,j3,j4))
                        Tensor1 = m_AmpExpo1(j1,j2,i3,i4)*DCONJG(m_AmpExpo1(j1,j2,j3,j4))
                        Tensor2 = m_AmpExpo2(j1,j2,i3,i4)*DCONJG(m_AmpExpo2(j1,j2,j3,j4))
                        Sum0 = Sum0 + Tensor0 *SDMprod
                        Sum1 = Sum1 + Tensor1 *SDMprod
                        Sum2 = Sum2 + Tensor2 *SDMprod
                     ENDDO                  
                     ENDDO
               ENDDO
               ENDDO
            ENDDO
         ENDDO
      ELSE
*        The case with POLARIZED beams and POLARIZED final fermions
         DO j1=1,2
         DO i1=1,2
            DO j2=1,2
            DO i2=1,2
               DO j3=1,2
               DO i3=1,2
                  DO j4=1,2
                  DO i4=1,2
                     SDMprod = m_SDMat1(i1,j1)*m_SDMat2(i2,j2)
     $                        *m_SDMat3(j3,i3)*m_SDMat4(j4,i4)
                     Tensor0 = m_AmpExpo0(i1,i2,i3,i4)*DCONJG(m_AmpExpo0(j1,j2,j3,j4))
                     Tensor1 = m_AmpExpo1(i1,i2,i3,i4)*DCONJG(m_AmpExpo1(j1,j2,j3,j4))
                     Tensor2 = m_AmpExpo2(i1,i2,i3,i4)*DCONJG(m_AmpExpo2(j1,j2,j3,j4))
                     Sum0 = Sum0 + Tensor0 *SDMprod
                     Sum1 = Sum1 + Tensor1 *SDMprod
                     Sum2 = Sum2 + Tensor2 *SDMprod
                  ENDDO
                  ENDDO
               ENDDO
               ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDIF
* distributions with polarized final fermions
      Rho0 = Sum0 *m_ExpoNorm
      Rho1 = Sum1 *m_ExpoNorm
      Rho2 = Sum2 *m_ExpoNorm
* Spin weight for  polarized final fermions
      wt0 = Rho0/m_RhoExp0
      wt1 = Rho1/m_RhoExp1
      wt2 = Rho2/m_RhoExp2
      END                       !!!GPS_MakeRho2!!!



      SUBROUTINE GPS_BornPlus(Mode,KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,Cfac,Xborn,Xboxy)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   Two-stroke version of GPS_Born, optimized for summation over partitions       //
*//   Virtual corrections (boxes vertices) to ms-one are also here!!!               //
*//   Warning! Input masses are TRUE (not signed as in GPS_Born)                    //
*//                                                                                 //
*//   Born spin amplitudes calculated with spinor methods.                          //
*//   Mass of the final fermion kept exactly.                                       //
*//                                                                                 //
*//   Input:                                                                        //
*//   KFi, Kff = beam and final fermion flavour codes (to define charges)           //
*//   PX       = s-chanel momentum for gamma and Z propagators (not for spinors)    //
*//   pi,mi    are for spinors, not for gamma and Z propagators                     //
*//   p1,m1    =fermion momentum and mass (beam)                                    //
*//   p2,m2    =fermion momentum and mass (beam)                                    //
*//   p3,m3    =fermion momentum and mass final state                               //
*//   p4,m4    =fermion momentum and mass final state                               //
*//                                                                                 //
*//   Output:                                                                       //
*//   Born     = spin summed squared amplitudes                                     //
*//                                                                                 //
*//   Common working space:                                                         //
*//   m_AmpBorn   is working space, used by HfinMinus                               //
*//   m_AmpExpo*  is working space, used by HiniPlus, HfinPlus, HfinMinus           //
*//                                                                                 //
*//   Notes:                                                                        //
*//   Electron mass neglected in spinors, this is why we may use Chisholm!          //
*//   Final fermion mass kept exactly.                                              //
*//   Gamma and Z in s-chanel.                                                      //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'GPS.h'
      SAVE                      !!! <-- necessary !!!
*
      INTEGER    Mode,KFi,KFf,Level
      DOUBLE PRECISION      PX(4),p1(4),p2(4),p3(4),p4(4)
      DOUBLE PRECISION      m1,m2,m3,m4,Xborn,Xboxy
*
      DOUBLE PRECISION      MasPhot
      DOUBLE PRECISION      BornSum,BoxySum
      DOUBLE PRECISION      SvarP, SvarQ, SvarX, Svar, s,t,u,s0,t0,u0
      LOGICAL               IFONE
      DOUBLE PRECISION      Fleps
*-----------------------------------------------------------------------------
      INTEGER         i,j,k,l
      INTEGER         j1,j2,j3,j4
      INTEGER         Hel1,Hel2,Hel3,Hel4
      DOUBLE COMPLEX  Cfac,AmpBorn,AmpBoxy
      DOUBLE COMPLEX  PropGam,PropZet
      DOUBLE COMPLEX  s31,s24,s14,s32
      DOUBLE COMPLEX  FFacTT(2),      FFacUU(2)
      DOUBLE COMPLEX  FFacTG(2),FFacTZ(2),      FFacUG(2),FFacUZ(2)
      DOUBLE COMPLEX  SpinoTT(2,2,2,2), SpinoUU(2,2,2,2)
      DOUBLE COMPLEX  BoxGG(2,2,2,2),   BoxGZ(2,2,2,2)
      DOUBLE COMPLEX  AmpBornW(2,2,2,2)
      DOUBLE COMPLEX  BVR_CBoxGG, BVR_CBoxGZ, BVR_IntIR, BVR_IntReson
      DOUBLE COMPLEX  Coef, IntIR
      DOUBLE COMPLEX  TT,UU
      DOUBLE COMPLEX  GPS_iProd1
      DOUBLE COMPLEX  GPS_iProd2
      DOUBLE PRECISION  F1iniPair,F1finPair
*-----------------------------------------------------------------------------
* Electroweak
      INTEGER           NCf,NCe
      DOUBLE PRECISION  T3e,Qe
      DOUBLE PRECISION  T3f,Qf
      DOUBLE COMPLEX    Ve,Ae,Vf,Af, VVCor, GamVPi, ZetVPi
      DOUBLE PRECISION  RsqV,RsqA ! QCD corrs.
      DOUBLE PRECISION  Svar9,CosThetD
*-----------------------------------------------------------------------------
      DOUBLE PRECISION        PP(4),QQ(4),dummy
c[[[[[[[[[[!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      DOUBLE PRECISION    icont
c      DATA       icont /0d0/
c]]]]]]]]]]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*-----------------------------------------------------------------------------
      Fleps =  1d-100
* Pure spinor part, does not depend on PX, Mode introduced to save CPU time
      IF(Mode .EQ. 0) THEN
         DO k=1,4
            PP(k) = p1(k)+p2(k)
            QQ(k) = p3(k)+p4(k)
         ENDDO
         SvarP  = PP(4)**2 - PP(3)**2 - PP(2)**2 - PP(1)**2
         SvarQ  = QQ(4)**2 - QQ(3)**2 - QQ(2)**2 - QQ(1)**2
*=============================================================
* Get charges, izospin, color
         CALL BornV_GetParticle(KFi, dummy, Qe,T3e,NCe)
         CALL BornV_GetParticle(KFf, dummy, Qf,T3f,NCf)
         CALL KK2f_GetMasPhot(MasPhot)
*=============================================================
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,*) '=================================GPS_BornPlus=========================================='
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
* Loop below correcponds to
****> CALL GPS_Born(KFi,KFf,PX, p1,Mbeam, p2,-Mbeam,  p3,Massf, p4,-Massf,m_AmpBorn) !!!!<****
         DO j1 = 1,2
            DO j2 = 1,2
               DO j3 = 1,2
                  DO j4 = 1,2
                     Hel1 = 3-2*j1
                     Hel2 = 3-2*j2
                     Hel3 = 3-2*j3
                     Hel4 = 3-2*j4
                     TT  = DCMPLX(0d0,0d0)
                     UU  = DCMPLX(0d0,0d0)
                     IF( Hel2 .EQ. -Hel1) THEN !!! <--helicity conservation imposed
                        s31 = GPS_iProd2(  Hel3, p3, m3,      Hel1, p1, Fleps) ! t
                        s24 = GPS_iProd2(  Hel2, p2,-Fleps,   Hel4, p4,-m4)    ! t1
                        s32 = GPS_iProd2(  Hel3, p3, m3,      Hel2, p2, Fleps) ! u1
                        s14 = GPS_iProd2(  Hel1, p1,-Fleps,   Hel4, p4,-m4)    ! u
                        TT  = s31*s24
                        UU  = s32*s14
                     ENDIF
                     SpinoTT(j1,j2,j3,j4) =  TT
                     SpinoUU(j1,j2,j3,j4) =  UU
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF                     !!! Mode.EQ.1
c[[[[[[[[[[[[[[[[[[[[[[[[[[[ *debug*
c      write(16,'(a,i1,a,i1,a,8f22.11)') (('SpinoTT(',j1,',',j2,',*,*)=',((SpinoTT(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c      write(16,'(a,i1,a,i1,a,8f22.11)') (('SpinoUU(',j1,',',j2,',*,*)=',((SpinoUU(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]

*////////////////////////////////////////////////////////////////////////////////////////////
*//                          Partition (PX) dependent part                                 //
*////////////////////////////////////////////////////////////////////////////////////////////
* Boxes: 
* IR-subtracted, helicity conservation imposed (small mass approx.)
      SvarX = PX(4)**2-PX(3)**2-PX(2)**2-PX(1)**2
      IF(SvarX .LE. (ABS(m3)+ABS(m4))**2 ) RETURN
      CALL KinLib_ThetaD(PX,p1,p2,p3,p4,Svar9,CosThetD)
      s =  SvarX
      t = -s*(1d0-CosThetD)/2d0
      u = -s*(1d0+CosThetD)/2d0
      Coef  = DCMPLX(m_Alfpi*Qe*Qf)
      IF(  m_KeyINT .NE. 0 .AND. m_KeyISR .NE. 0 .AND.  m_HasFSR .NE. 0  ) THEN
* Virtual 2*(B(t)-B(u)) Intereference IR part to be subtracted from boxes
         IntIR      = Coef*BVR_IntIR(MasPhot,s,t,u)                 !!<- asymetric in (t,u)
         m_IntReson = Coef*BVR_IntReson(MasPhot,m_MZ,m_GammZ,s,t,u) !!<- asymetric in (t,u)
         m_BoxGGtu  = Coef*( BVR_CBoxGG(MasPhot,             s,t,u)) -IntIR
         m_BoxGZtu  = Coef*( BVR_CBoxGZ(MasPhot,m_MZ,m_GammZ,s,t,u)) -IntIR
         m_BoxGGut  = Coef*(-BVR_CBoxGG(MasPhot,             s,u,t)) -IntIR
         m_BoxGZut  = Coef*(-BVR_CBoxGZ(MasPhot,m_MZ,m_GammZ,s,u,t)) -IntIR
* Exponentiate Resonance BigLogs according to Greco et al.
         IF( m_KeyInt .EQ. 2) THEN
            m_BoxGZtu = m_BoxGZtu -m_IntReson
            m_BoxGZut = m_BoxGZut -m_IntReson
         ENDIF
      ELSE
         m_IntReson = DCMPLX(0d0,0d0)
         m_BoxGGtu  = DCMPLX(0d0,0d0)
         m_BoxGGut  = DCMPLX(0d0,0d0)
         m_BoxGZtu  = DCMPLX(0d0,0d0)
         m_BoxGZut  = DCMPLX(0d0,0d0)
      ENDIF
c[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
c[[[[  switching off boxes by hand
c      m_BoxGZtu = DCMPLX(0d0,0d0)
c      m_BoxGZut = DCMPLX(0d0,0d0)
c      m_BoxGGtu = DCMPLX(0d0,0d0)
c      m_BoxGGut = DCMPLX(0d0,0d0)
c]]]]
c]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
      DO j1 = 1,2
         DO j2 = 1,2
            DO j3 = 1,2
               DO j4 = 1,2
                  Hel1 = 3-2*j1
                  Hel2 = 3-2*j2
                  Hel3 = 3-2*j3
                  Hel4 = 3-2*j4
                  BoxGG(j1,j2,j3,j4) =DCMPLX(0d0,0d0)
                  BoxGZ(j1,j2,j3,j4) =DCMPLX(0d0,0d0)
                  IF((Hel2 .EQ. -Hel1) .AND. (Hel4 .EQ. -Hel3)) THEN !!<--helicity conserv.
                     IF( Hel1*Hel3 .EQ. 1) THEN
                        BoxGG(j1,j2,j3,j4) = m_BoxGGtu
                        BoxGZ(j1,j2,j3,j4) = m_BoxGZtu
                     ELSE
                        BoxGG(j1,j2,j3,j4) = m_BoxGGut
                        BoxGZ(j1,j2,j3,j4) = m_BoxGZut
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
c[[[[[[[[[[[[[!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      icont=icont+1
c      IF(icont.LE.10) THEN
c         write(m_out,*) ' //////////////////////GPS///////////////////////////////////////'
cc         write(*,'(a,5g22.14)') 'CosThetD= ',CosThetD
cc         write(*,'(a,5g22.14)') 'Sw2= ',m_Sw2
c      ENDIF
c]]]]]]]]]]]]]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
* Propagators, with s-dependent width
      PropGam =    DCMPLX(  1d0/SvarX,  0d0)
      PropZet =    1d0/DCMPLX(SvarX-m_MZ**2, m_GammZ*SvarX/m_MZ)
* Possibility to switch off Z or gamma, etc.
      IF(m_KeyZet .LE. 0) PropZet =  DCMPLX(0d0)
      IF(m_KeyZet .EQ. 9) PropGam =  DCMPLX(0d0)
      IF(m_KeyZet .EQ.-1) PropZet =  1d0/DCMPLX(SvarX-m_MZ**2, m_GammZ*m_MZ)
* Exponentiate Resonance BigLogs according to Greco et al.
      IF(  m_KeyINT .EQ. 2 .AND. m_KeyISR .NE. 0 .AND.  m_HasFSR .NE. 0  ) THEN
         PropZet = PropZet * EXP(m_IntReson)
      ENDIF
*////////////////////////////////////////////////////////////////////////////////////////////
*//                        ElectroWeak Corrections                                         //    
*////////////////////////////////////////////////////////////////////////////////////////////
      CALL GPS_EWFFact(KFi,KFf,SvarX,CosThetD,Ve,Vf,Ae,Af,VVcor,GamVPi,ZetVPi,RsqV,RsqA) !
*////////////////////////////////////////////////////////////////////////////////////////////
*//     Primitives formfactor-type for construction of spin amplitudes                     //
*//     For boxes we need separately photon and Z parts                                    //
*//     (Ve -Hel1*Ae)*(Vf +Hel1*Af) is expanded because of double-vector f-factor          //
*////////////////////////////////////////////////////////////////////////////////////////////
      DO j1 = 1,2
         Hel1 = 3-2*j1
         FFacTG(j1) = PropGam*GamVPi *Qe*Qf *RsqV
         FFacTZ(j1) = PropZet*ZetVPi *(Ve*Vf*VVCor*RsqV -Hel1*Ae*Vf*RsqV +Hel1*Ve*Af*RsqA -Ae*Af*RsqA) !
         FFacUG(j1) = PropGam*GamVPi *Qe*Qf *RsqV
         FFacUZ(j1) = PropZet*ZetVPi *(Ve*Vf*VVCor*RsqV -Hel1*Ae*Vf*RsqV -Hel1*Ve*Af*RsqA +Ae*Af*RsqA) !
         FFacTT(j1) = FFacTG(j1)+FFacTZ(j1)
         FFacUU(j1) = FFacUG(j1)+FFacUZ(j1)
      ENDDO
c[[[[[[[[[[[[[[[[[[[[[[[[[[[*debug*
c      write(16,*) '///////////////////////////////////////////GPS_BornPlus///////////////////////////////////////////////////////'
c      write(16,'(a,4f22.11)') ' FFacTT=', (FFacTT(j),j=1,2)
c      write(16,'(a,4f22.11)') ' FFacUU=', (FFacUU(j),j=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
*///////////////////////////////////////////////////////////////////////////////////
*//       QED vertex  FFactor F1 minus B-vrtual (IR removed), exact mass terms    //
*///////////////////////////////////////////////////////////////////////////////////
      m_F1ini1 = DCMPLX(0d0)
      m_F1fin1 = DCMPLX(0d0)
      m_F1ini2 = DCMPLX(0d0)
      m_F1fin2 = DCMPLX(0d0)
      m_F1iniPair2 = DCMPLX(0d0)
      m_F1finPair2 = DCMPLX(0d0)
      IF( m_KeyISR .NE. 0) THEN
         CALL  BVR_MakeF1ini(SvarP,m1,m2,m_Alfpi,Qe,m_F1ini1,m_F1ini2)
         CALL BVR_MakeF2Pair(SvarP,m1,m_Alfpi,Qe,F1iniPair)
         m_F1iniPair2 = F1iniPair ! real to Complex
      ENDIF
      IF( m_HasFSR .NE. 0) THEN
         CALL  BVR_MakeF1fin(SvarQ,m3,m4,m_Alfpi,Qf,m_F1fin1,m_F1fin2)
         CALL BVR_MakeF2Pair(SvarQ,m3,m_Alfpi,Qf,F1finPair)
         m_F1finPair2 = F1finPair ! real to Complex
      ENDIF
c[[[[[[[
c      IF( DABS(SvarQ/SvarP-1d0) .LT. 1d-10 ) THEN
c      write(*,*) ' GPS:ABS(SvarQ/SvarP-1)', DABS(SvarQ/SvarP-1d0)
c      write(*,*) ' GPS: alf1 ',CDABS(1+ m_F1ini1)**2*CDABS(1+ m_F1fin1)**2
c      write(*,*) ' GPS: alf2 ',CDABS(1+ m_F1ini2)**2*CDABS(1+ m_F1fin2)**2
ccc      write(*,*) ' GPS: STOP'
ccc      STOP
c      ENDIF
c]]]]]]]

*///////////////////////////////////////////////////////////////////////////////////
*//                      Total result = Spinors*Formfactor                        //
*///////////////////////////////////////////////////////////////////////////////////
      Level=1  ! W amplitude calculated separately for  level=1, and masses are treated as in BornPlus

       IF ((p3(3)+p4(3))*p1(3).LT.0D0) THEN ! Where is dominant other photon?  It is instead of reduction procedure
         IFONE=.TRUE.      ! We assume all extra photons were emitted from p1 
        ELSE
         IFONE=.FALSE.     ! We assume all extra photons were emitted from p2
        ENDIF
C--------
         IF (IFONE) THEN
           s0=(p3(4)+p4(4))**2-(p3(3)+p4(3))**2-(p3(2)+p4(2))**2-(p3(1)+p4(1))**2
           t0=(p4(4)-p2(4))**2-(p4(3)-p2(3))**2-(p4(2)-p2(2))**2-(p4(1)-p2(1))**2
           u0=(p3(4)-p2(4))**2-(p3(3)-p2(3))**2-(p3(2)-p2(2))**2-(p3(1)-p2(1))**2
         ELSE
           s0=(p3(4)+p4(4))**2-(p3(3)+p4(3))**2-(p3(2)+p4(2))**2-(p3(1)+p4(1))**2
           t0=(p3(4)-p1(4))**2-(p3(3)-p1(3))**2-(p3(2)-p1(2))**2-(p3(1)-p1(1))**2
           u0=(p4(4)-p1(4))**2-(p4(3)-p1(3))**2-(p4(2)-p1(2))**2-(p4(1)-p1(1))**2
         ENDIF
      call GPS_BornWPlus(Mode,Level,KFi,KFf,s0,t0,u0,p1,m1,p2,m2,p3,m3,p4,m4,AmpBornW)
      BornSum = 0d0
      BoxySum  = 0d0

      DO j1 = 1,2
         DO j2 = 1,2
            DO j3 = 1,2
               DO j4 = 1,2
* Born,  Zero order
                  AmpBorn = SpinoTT(j1,j2,j3,j4)* FFacTT(j1)
     $                     +SpinoUU(j1,j2,j3,j4)* FFacUU(j1)
     $                     +AmpBornW(j1,j2,j3,j4)
* Boxes, First order
                  AmpBoxy = SpinoTT(j1,j2,j3,j4)* FFacTG(j1) *BoxGG(j1,j2,j3,j4)
     $                     +SpinoTT(j1,j2,j3,j4)* FFacTZ(j1) *BoxGZ(j1,j2,j3,j4)
     $                     +SpinoUU(j1,j2,j3,j4)* FFacUG(j1) *BoxGG(j1,j2,j3,j4)
     $                     +SpinoUU(j1,j2,j3,j4)* FFacUZ(j1) *BoxGZ(j1,j2,j3,j4)
* Store results
                  m_AmpBorn( j1,j2,j3,j4) = AmpBorn
                  m_AmpBoxy( j1,j2,j3,j4) = AmpBoxy*Cfac
                  m_AmpExpo0(j1,j2,j3,j4) = m_AmpExpo0(j1,j2,j3,j4) 
     $                 +Cfac*AmpBorn                               !!! pure Born
                  m_AmpExpo1(j1,j2,j3,j4) = m_AmpExpo1(j1,j2,j3,j4) 
     $                 +Cfac*AmpBorn*(1 +m_F1ini1)*(1 +m_F1fin1 )  !!! Born, O(alf1) FFactors
     $                 +Cfac*AmpBoxy                               !!! O(alf1) boxes
c[[[[[[[[[[[[[[[[[[[[[[[[[[[ *debug*
c[[[[[[test                  m_AmpExpo1(j1,j2,j3,j4) = dcmplx(0d0)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
*///////// new !!!!!!!! in construction !!!!!!! ////////////////////
                  m_AmpExpo2(j1,j2,j3,j4) = m_AmpExpo2(j1,j2,j3,j4) 
     $                 +Cfac*AmpBorn*(1 +m_F1ini2)*(1 +m_F1fin2 )  !!! Born, O(alf2) FFactors
     $                 +Cfac*AmpBoxy                               !!! O(alf1) boxes
                  m_AmpExpo2p(j1,j2,j3,j4) = m_AmpExpo2p(j1,j2,j3,j4) 
     $                 +Cfac*AmpBorn*(1 +m_F1ini2+m_F1iniPair2)*(1 +m_F1fin2+m_F1finPair2 ) ! +Pairs
     $                 +Cfac*AmpBoxy                               !!! O(alf1) boxes
*/////////////////////////////
                  BornSum = BornSum +Cfac*AmpBorn*DCONJG(Cfac*AmpBorn)
                  BoxySum = BoxySum +2*DREAL(Cfac*AmpBoxy*DCONJG(Cfac*AmpBorn))
               ENDDO                  
            ENDDO
         ENDDO
      ENDDO
c[[[[[[[[[[[[[[[[[[[[[[[[[[[ *debug*
c      write(16,*) '///////////////////////////////////////////GPS_BornPlus///////////////////////////////////////////////////////'
c      write(16,'(a,i1,a,i1,a,8f22.11)') (('m_AmpBorn(',j1,',',j2,',*,*)=',((m_AmpBorn(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c      write(16,*) '--------------------------------------------------------------------------------------------------------------'
c      write(16,'(a,i1,a,i1,a,8f22.11)') ((' AmpBornW(',j1,',',j2,',*,*)=',(( AmpBornW(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
      Xborn = BornSum
      Xboxy = BoxySum
* debug variables for tests
      m_debg( 1) = BornSum
      m_debg( 2) = BoxySum
      END                       !!!GPS_BornPlus!!!


      SUBROUTINE GPS_BornWPlus(Mode,Level,KFi,KFf,s,t,u,p1,m1,p2,m2,p3,m3,p4,m4,AmpBorn)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   correction  of W exchange to born                                             //
*//   Warning! Input masses are TRUE (not signed as in GPS_Born)                    //
*//                                                                                 //
*//   Born spin amplitudes calculated with spinor methods.                          //
*//   Mass of the final fermion kept exactly.                                       //
*//                                                                                 //
*//   Input:                                                                        //
*//   KFi, Kff = beam and final fermion flavour codes (to define charges)           //
*//   PX       = s-chanel momentum for W propagator (not for spinors)               //
*//   pi,mi    are for spinors, not for W propagator                                //
*//   p1,m1    =fermion momentum and mass (beam)                                    //
*//   p2,m2    =fermion momentum and mass (beam)                                    //
*//   p3,m3    =fermion momentum and mass final state                               //
*//   p4,m4    =fermion momentum and mass final state                               //
*//                                                                                 //
*//   Mode     time saving facility: not active                                     //
*//   Level    type of treatment of spinor products, probably irrelevant            //
*//   JakKoralZ=0 internal input for some old tests                                 //
*//                                                                                 //
*//   Output:                                                                       //
*//   AmpBorn     = spin amplitude                                                  //
*//                                                                                 //
*//                                                                                 //
*//   Notes:                                                                        //
*//   Electron mass neglected in spinors, this is why we may use Chisholm!          //
*//   Final fermion mass kept exactly.                                              //
*//   W in t-chanel.                                                                //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'GPS.h'
      SAVE                      !!! <-- necessary !!!
*
      INTEGER    Mode,KFi,KFf,Level,JakKoralZ
      DOUBLE PRECISION      p1(4),p2(4),p3(4),p4(4)
      DOUBLE PRECISION      m1,m2,m3,m4,Xborn
*
      DOUBLE PRECISION      BornSum
      DOUBLE PRECISION      s,t,u
      DOUBLE PRECISION      Fleps,m_MW, m_GammW
*-----------------------------------------------------------------------------
      INTEGER    i,j,k,l
      INTEGER    j1,j2,j3,j4
      INTEGER    Hel1,Hel2,Hel3,Hel4
      DOUBLE COMPLEX  Cfac,AmpBornW
      DOUBLE COMPLEX  PropW
      DOUBLE COMPLEX  s31,s24,s14,s32
      DOUBLE COMPLEX  FFacTT(2),      FFacUU(2)
      DOUBLE COMPLEX  FFacTG(2),FFacTZ(2),      FFacUG(2),FFacUZ(2)
      DOUBLE COMPLEX  SpinoTT(2,2,2,2), SpinoUU(2,2,2,2)
      DOUBLE COMPLEX  AmpBorn(2,2,2,2)
      DOUBLE COMPLEX  Coef, IntIR
      DOUBLE COMPLEX  TT,UU
      DOUBLE COMPLEX  GPS_iProd1
      DOUBLE COMPLEX  GPS_iProd2
*-----------------------------------------------------------------------------
* Electroweak
      INTEGER      NCf,NCe
      DOUBLE PRECISION        T3e,Qe
      DOUBLE PRECISION        T3f,Qf
      DOUBLE COMPLEX    Ve,Ae,Vf,Af, VVCor, WVPi
*-----------------------------------------------------------------------------
      DOUBLE PRECISION        dummy
c[[[[[[[[[[!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      DOUBLE PRECISION    icont
c      DATA       icont /0d0/
c]]]]]]]]]]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*-----------------------------------------------------------------------------
      IF (level.EQ.1) THEN
        DO j1 = 1,2
           DO j2 = 1,2
              DO j3 = 1,2
                 DO j4 = 1,2
                    AmpBorn( j1,j2,j3,j4) = 0D0
                 ENDDO                  
              ENDDO
           ENDDO
        ENDDO
        ENDIF
      IF (ABS(KFf) .NE. 12)  RETURN

      JakKoralZ=0  ! warning: defined in 2 places ! to fix a potential problem .eq.1 like KORALZ .eq.0 possiby OK
      Fleps =  1d-100

* Pure spinor part, does not depend on PX, Mode introduced to save CPU time
! zbw      IF(Mode .EQ. 0) THEN
*=============================================================
* Get charges, izospin, color
         CALL BornV_GetParticle(KFi, dummy, Qe,T3e,NCe)
         CALL BornV_GetParticle(KFf, dummy, Qf,T3f,NCf)
*=============================================================
* Loop below correcponds to
****> CALL GPS_Born(KFi,KFf,PX, p1,Mbeam, p2,-Mbeam,  p3,Massf, p4,-Massf,m_AmpBorn) !!!!<****
         DO j1 = 1,2
            DO j2 = 1,2
               DO j3 = 1,2
                  DO j4 = 1,2
                     Hel1 = 3-2*j1
                     Hel2 = 3-2*j2
                     Hel3 = 3-2*j3
                     Hel4 = 3-2*j4
                     TT  = DCMPLX(0d0,0d0)
                     UU  = DCMPLX(0d0,0d0)
                     IF( Hel2 .EQ. -Hel1) THEN !!! <--helicity conservation imposed
                       IF(level.eq.1) THEN
                        s31 = GPS_iProd2(  Hel3, p3, m3,      Hel1, p1, Fleps) ! t
                        s24 = GPS_iProd2(  Hel2, p2,-Fleps,   Hel4, p4,-m4)    ! t1
                        s32 = GPS_iProd2(  Hel3, p3, m3,      Hel2, p2, Fleps) ! u1
                        s14 = GPS_iProd2(  Hel1, p1,-Fleps,   Hel4, p4,-m4)    ! u
                       ELSE
                        s31 = GPS_iProd2(  Hel3, p3, m3,      Hel1, p1, m1)    ! t
                        s24 = GPS_iProd2(  Hel2, p2, m2,      Hel4, p4, m4)    ! t1
                        s32 = GPS_iProd2(  Hel3, p3, m3,      Hel2, p2,-m2)    ! u1
                        s14 = GPS_iProd2(  Hel1, p1,-m1,      Hel4, p4, m4)    ! u
                       ENDIF
                        TT  = s31*s24
                        UU  = s32*s14
                     ENDIF
                     SpinoTT(j1,j2,j3,j4) =  TT
                     SpinoUU(j1,j2,j3,j4) =  UU
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
! zbw      ENDIF                     !!! Mode.EQ.1
*////////////////////////////////////////////////////////////////////////////////////////////
*//                          Partition (PX) dependent part                                 //
*////////////////////////////////////////////////////////////////////////////////////////////

      IF(s .LE. (ABS(m3)+ABS(m4))**2 ) RETURN
c]]]]
*////////////////////////////////////////////////////////////////////////////////////////////
*//                        ElectroWeak Corrections in/and W propagator                     //    
*////////////////////////////////////////////////////////////////////////////////////////////
      Coef  =1.d0/2.d0/m_Sw2

      IF(JakKoralZ.eq.1) then
         CALL GPS_EWFFactW(KFi,KFf,s,u,PropW,WVPi) !W Propagator: it looks crazy in this case ....
      ELSE
        CALL GPS_EWFFactW(KFi,KFf,s,t,PropW,WVPi)
      ENDIF


*////////////////////////////////////////////////////////////////////////////////////////////
*//     Primitives formfactor-type for construction of spin amplitudes                     //
*////////////////////////////////////////////////////////////////////////////////////////////
      IF(JakKoralZ.eq.1) then ! switch for some tests
        Ve=-0.5  ! 
        Vf= 0.5  ! 
        Ae= 0.5 ! 
        Af= 0.5
      ELSE
        Ve= 0.5  ! 
        Vf= 0.5  ! 
        Ae= 0.5 ! 
        Af= 0.5
      ENDIF
      DO j1 = 1,2
         Hel1 = 3-2*j1
         FFacTZ(j1) = PropW*WVPi *(Ve*Vf -Hel1*Ae*Vf +Hel1*Ve*Af -Ae*Af)
         FFacUZ(j1) = PropW*WVPi *(Ve*Vf -Hel1*Ae*Vf -Hel1*Ve*Af +Ae*Af)
         FFacTT(j1) = FFacTZ(j1) * Coef
         FFacUU(j1) = FFacUZ(j1) * Coef
      ENDDO
*///////////////////////////////////////////////////////////////////////////////////
*//                      Total result = Spinors*Formfactor                        //
*///////////////////////////////////////////////////////////////////////////////////
      BornSum = 0d0
      DO j1 = 1,2
         DO j2 = 1,2
            DO j3 = 1,2
               DO j4 = 1,2
* Born,  Zero order
                  AmpBornW = SpinoTT(j1,j2,j3,j4)* FFacTT(j1)
     $                     +SpinoUU(j1,j2,j3,j4)* FFacUU(j1)
                  AmpBorn( j1,j2,j3,j4) = AmpBornW +AmpBorn( j1,j2,j3,j4)
*/////////////////////////////
                                     Cfac=1D0 ! Ad hoc, test may be destroyed !
                  BornSum = BornSum +Cfac*AmpBornW*DCONJG(Cfac*AmpBornW)
               ENDDO                  
            ENDDO
         ENDDO
      ENDDO
c[[[[[[[[[[[[[[[[[[[[[[[*debug*
c      write(16,*) '=================================GPS_BornWPlus=========================================='
c      write(16,'(a,i1,a,i1,a,8f22.11)') (( 'AmpBornW(',j1,',',j2,',*,*)=',(( AmpBorn(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c      write(16,*) 's,t,PropW,WVPi=',s,t,PropW,WVPi
c      write(16,'(a,4f22.11)') ' FFacTT=', (FFacTT(j),j=1,2)
c      write(16,'(a,4f22.11)') ' FFacUU=', (FFacUU(j),j=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]
!      write(*,*) ampbornw
!      stop
      Xborn = BornSum
* debug variables for tests
!      m_debg( 1) = BornSum
      END                       !!!GPS_BornWPlus!!!

      SUBROUTINE GPS_EWFFactW(KFi,KFf,s,t,PropW,WVPi)
*////////////////////////////////////////////////////////////////////////////////////////////
*//                        ElectroWeak Corrections                                         //    
*//                        empty prepared for nunu chanel                                  //
*////////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'GPS.h'
* Input
      INTEGER    KFi,KFf
      DOUBLE PRECISION      s,t,u
* Output
      DOUBLE COMPLEX    Ve,Vf,Ae,Af,VVcor,GamVPi,WVPi
* Local
      INTEGER            NCf,NCe,Jak
      DOUBLE PRECISION   T3e,Qe
      DOUBLE PRECISION   T3f,Qf,m_MW,m_GammW,m_MW_true
      DOUBLE COMPLEX     PropW,Coef,Coef1,DelW
      DOUBLE COMPLEX     RhoEW, VPgamma, CorEle, CorFin, CorEleFin, VVCef, CosThetD
      DOUBLE PRECISION   Deno,  dummy,row,Dum1,Dum2,Dum3,Dum4,Dum5,Dum6
      DOUBLE COMPLEX     GSW(100)
      DOUBLE PRECISION icont
      DATA icont /0d0/
      icont = icont +1
*===============================================================================
C we do it like for photon. I was unable to properly understand normalization
C this is half-guess half observation of the KK-KORALZ differences in code.
C it can be easily wrong. nothing like this tem is present in KORALZ, 
C nonetheless KORALZ results are 137/127 higher than KK for pure W exchange.
C in koralz this vacpol factor is installed at born step mode 1/0 ratio
C here not
      Jak= 1
      IF (Jak.eq.0) then ! old version for some tests only
        m_MW  =m_MZ*dsqrt(1d0-m_Sw2)
        m_GammW = 2.5d0
        PropW   =    1d0/DCMPLX(t-m_MW**2, m_GammW*m_MW)
        WVPi=DCMPLX(1.06794821,-0.0185004916)
        WVPi=DCMPLX(1.0,0.0)
      ELSE
        m_MW  =m_MZ*dsqrt(1d0-m_Sw2)
        IF ( m_KeyElw.GE.1) THEN
*         next will be mass fro common /cdzwsm/ from  DZface.f it is actually raw mass as above ...
          CALL BornV_GetMW(m_MW)
*         use of DZface forbiden, disk-mode does not work
*         CALL DZface_GetPrm( Dum1,Dum2,Dum3,m_MW,Dum4,Dum5,Dum6)
        ENDIF
        PropW   =    1d0/DCMPLX(t-m_MW**2,0D0 )
        WVPi=DCMPLX(1.06794821,-0.0185004916)
        WVPi=DCMPLX(1.0,0.0)
        Coef  =1.d0/2.d0/m_Sw2
c[[[[        Coef1=m_Gmu*m_MW**2* m_Alfinv/4D0 ! ERROR
        Coef1=m_Gmu*m_MW**2* m_Alfinv/(DSQRT(2d0)*m_pi)

        IF ( m_KeyElw.GE.1) THEN
        WVPi= WVPi*coef1/coef  ! effective coupling
           CosThetD=1+2*t/s 
! DIRTY trick to get around the problem. Anyway problem appears in non used
!           WVPi's where it is later overwritten. 
!           Physically it is ambiguity of order alpha**2
           IF (1D0.LT.(1+2*t/s)**2) THEN
             CosThetD=1D0/(1+2*t/s)
           ENDIF
           CALL BornV_InterpoGSW(KFf,S,CosThetD)
           CALL BornV_GetGSW(GSW)
*     use of rhocc forbiden, disk-mode does not work
c           CALL rhocc(s,-t,-s-t,-1D0,1D0,0D0,0D0,ROW)
c           IF(icont.LE.100) WRITE(*,*) '%%%%%%   ',GSW(5),ROW,'   ',s,t
c           IF (ABS(row-gsw(5)).GT.0.0002.AND.s.GT.10000) WRITE(*,*) GSW(5),ROW,'   ',s,t 
c           IF (ABS(row-1).gt.0.016) STOP

* this was added to GSW(5) in BornV_Dizet
           DelW= 1D0/m_AlfInv/m_pi/2*(-3D0/2*LOG(s/m_MW**2)+1D0/2*(LOG(-t/s))**2-4d0/6d0*m_pi**2+2D0) ! (b)
cccc       DelW= 1D0/m_AlfInv/m_pi/2*(-3D0/2*LOG(s/m_MW**2)+1D0/2*(LOG(-t/s))**2-1d0/6d0*m_pi**2+2D0) ! (a)
           WVPi= WVPi*(GSW(5)+DeLW)
        ENDIF
      ENDIF
c[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,*) '///////////////////////////////GPS_EWFFactW//////////////////////////////////////////////'
c      write(16,*) 'm_KeyElw,Jak=',m_KeyElw,Jak
c      write(16,*) 'm_MW,m_GammW, s,t=',m_MW,m_GammW, s,t
c      write(16,*) 'PropW,WVPi=',PropW,WVPi
c]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
      END


      SUBROUTINE GPS_Born(KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,AmpBorn)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   Used in construction of the hard Non-IR parts: in GPS_HiniPlus, GPS_HfinPlus  //
*//   It is esentialy a simplified clone of GPS_BornPlus                            //
*//                                                                                 //
*//   CAREFUL!!! p_i are sometimes be substituted for photons!!!                    //
*//                                                                                 //
*//   Mass of the final fermion kept exactly.                                       //
*//                                                                                 //
*//   Input:                                                                        //
*//   KFi, Kff = beam and final fermion flavour codes (to define charges)           //
*//   PX       = s-chanel momentum for gamma and Z propagators (not for spinors)    //
*//   pi,mi    are for spinors, not for gamma and Z propagators                     //
*//   p1,m1    =fermion momentum and mass (beam)                                    //
*//   p2,m2    =fermion momentum and mass (beam)                                    //
*//   p3,m3    =fermion momentum and mass final state                               //
*//   p4,m4    =fermion momentum and mass final state                               //
*//                                                                                 //
*//   Output:                                                                       //
*//   AmpBorn   = spin amplitudes                                                   //
*//                                                                                 //
*//   Notes:                                                                        //
*//   Electron mass neglected in spinors, this is why we may use Chisholm!          //
*//   Final fermion mass kept exactly.                                              //
*//   Gamma and Z in s-chanel.                                                      //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'GPS.h'
*
      INTEGER    KFi,KFf,Mode,Level
      DOUBLE PRECISION  PX(4),p1(4),p2(4),p3(4),p4(4)
      DOUBLE PRECISION  m1,m2,m3,m4
      DOUBLE COMPLEX        AmpBorn(2,2,2,2),AmpBornW(2,2,2,2)
      DOUBLE COMPLEX        Cfac
      DOUBLE PRECISION      Xborn
*
      DOUBLE PRECISION  T3e,Qe
      DOUBLE PRECISION  T3f,Qf
      DOUBLE COMPLEX    Ve,Vf,Ae,Af,VVcor,GamVPi,ZetVPi
      DOUBLE PRECISION  RsqV,RsqA ! QCD corrs.
      INTEGER           NCf,NCe
      DOUBLE PRECISION  svarX
*-----------------------------------------------------------------------------
      INTEGER           i,j,k,l
      INTEGER           j1,j2,j3,j4
      INTEGER           Hel1,Hel2,Hel3,Hel4
      DOUBLE COMPLEX    PropGam,PropZet
      DOUBLE COMPLEX    s31,s24,s14,s32
      DOUBLE COMPLEX    FFacTT(2),      FFacUU(2)
      DOUBLE COMPLEX    SpinoTT(2,2,2,2),SpinoUU(2,2,2,2)
      DOUBLE COMPLEX    TT,UU
      DOUBLE COMPLEX    GPS_iProd1
      DOUBLE COMPLEX    GPS_iProd2
      DOUBLE PRECISION  dummy
*-----------------------------------------------------------------------------
      CALL GPS_Initialize
*=============================================================
* Get charges, izospin, color
      CALL BornV_GetParticle(KFi, dummy, Qe,T3e,NCe)
      CALL BornV_GetParticle(KFf, dummy, Qf,T3f,NCf)
*=============================================================
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,*) '=================================GPS_Born=================================================='
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
      DO j1 = 1,2
         DO j2 = 1,2
            DO j3 = 1,2
               DO j4 = 1,2
                  Hel1 = 3-2*j1
                  Hel2 = 3-2*j2
                  Hel3 = 3-2*j3
                  Hel4 = 3-2*j4
                  TT  = DCMPLX(0d0,0d0)
                  UU  = DCMPLX(0d0,0d0)
                  IF( Hel2 .EQ. -Hel1) THEN   !!! <--helicity conservation imposed
                     s31 = GPS_iProd2(  Hel3, p3, m3,   Hel1, p1, m1) ! t
                     s24 = GPS_iProd2(  Hel2, p2, m2,   Hel4, p4, m4) ! t1
                     s32 = GPS_iProd2(  Hel3, p3, m3,   Hel2, p2,-m2) ! u1
                     s14 = GPS_iProd2(  Hel1, p1,-m1,   Hel4, p4, m4) ! u
                     TT  = s31*s24
                     UU  = s14*s32
                  ENDIF
                  SpinoTT(j1,j2,j3,j4) =  TT
                  SpinoUU(j1,j2,j3,j4) =  UU
               ENDDO
            ENDDO
         ENDDO
      ENDDO
c[[[[[[[[[[[[[[[[[[[[[[[[[[[ *debug*
c      write(16,'(a,i1,a,i1,a,8f22.11)') (('SpinoTT(',j1,',',j2,',*,*)=',((SpinoTT(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c      write(16,'(a,i1,a,i1,a,8f22.11)') (('SpinoUU(',j1,',',j2,',*,*)=',((SpinoUU(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
*////////////////////////////////////////////////////////////////////////////////////////////
*//                        ElectroWeak Corrections                                         //
*//   CosThetD = 0d0 is not mistake, we neglect miniscule second order effects only        //
*////////////////////////////////////////////////////////////////////////////////////////////
      svarX=PX(4)**2-PX(3)**2-PX(2)**2-PX(1)**2
      IF(svarX .LE. (ABS(m3)+ABS(m4))**2 ) RETURN
      CALL GPS_EWFFact(KFi,KFf,SvarX,0d0 ,Ve,Vf,Ae,Af,VVcor,GamVPi,ZetVPi,RsqV,RsqA) !
* Propagators, with s-dependent width
      PropGam =    DCMPLX(  1d0/svarX,  0d0)
      PropZet =    1d0/DCMPLX(svarX-m_MZ**2, m_GammZ*svarX/m_MZ)
* Possibility to switch off Z or gamma, etc.
      IF(m_KeyZet .LE. 0) PropZet =  DCMPLX(0d0)
      IF(m_KeyZet .EQ. 9) PropGam =  DCMPLX(0d0)
      IF(m_KeyZet .EQ.-1) PropZet =  1d0/DCMPLX(SvarX-m_MZ**2, m_GammZ*m_MZ)
* Exponentiate Resonance BigLogs according to Greco et al.
      IF(  m_KeyINT .EQ. 2 .AND. m_KeyISR .NE. 0 .AND.  m_HasFSR .NE. 0  ) THEN
         PropZet = PropZet * EXP(m_IntReson)
      ENDIF
*////////////////////////////////////////////////////////////////////////////////////////////
*//     Primitives formfactor-type for construction of spin amplitudes                     //
*//     (Ve -Hel1*Ae)*(Vf +Hel1*Af) is expanded because of double-vector f-factor          //
*////////////////////////////////////////////////////////////////////////////////////////////
      DO j1 = 1,2
         Hel1 = 3-2*j1
         FFacTT(j1) = PropGam*GamVPi *Qe*Qf *RsqV
     $               +PropZet*ZetVPi *(Ve*Vf*VVCor*RsqV -Hel1*Ae*Vf*RsqV +Hel1*Ve*Af*RsqA -Ae*Af*RsqA) !
         FFacUU(j1) = PropGam*GamVPi *Qe*Qf *RsqV
     $               +PropZet*ZetVPi *(Ve*Vf*VVCor*RsqV -Hel1*Ae*Vf*RsqV -Hel1*Ve*Af*RsqA +Ae*Af*RsqA) !
      ENDDO
c[[[[[[[[[[[[[[[[[[[[[[[[[[[*debug*
c      write(16,'(a,4f22.11)') ' FFacTT=', (FFacTT(j),j=1,2)
c      write(16,'(a,4f22.11)') ' FFacUU=', (FFacUU(j),j=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
*////////////////////////////////////////////////////////////////////////////////////////////
*//                     Total result = Spinors*Formfactor                                  //
*////////////////////////////////////////////////////////////////////////////////////////////
!>>>>      Mode=0
!>>>>      Level=1  ! t-dependent W-propagator for level=1
!>>>>      call GPS_BornWPlus(Mode,Level,KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,AmpBornW)
      DO j1 = 1,2
         DO j2 = 1,2
            DO j3 = 1,2
               DO j4 = 1,2
                  AmpBorn(j1,j2,j3,j4) =  SpinoTT(j1,j2,j3,j4)* FFacTT(j1)
     $                                   +SpinoUU(j1,j2,j3,j4)* FFacUU(j1)
!>>>>     $                                   +AmpBornW(j1,j2,j3,j4)
               ENDDO                  
            ENDDO
         ENDDO
      ENDDO
c[[[[[[[[[[[[[[[[[[[[[[[[[[[ *debug*
c      write(16,'(a,i1,a,i1,a,8f22.11)') (('AmpBorn(',j1,',',j2,',*,*)=',((AmpBorn(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
      END                       !!!GPS_Born!!!



      SUBROUTINE GPS_BornF(KFi,KFf,PX,CosThe,p1,m1,p2,m2,p3,m3,p4,m4,Xborn)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//                              Pure Born                                         //
*//                                                                                 //
*//   It is esentialy a clone of GPS_Born with elements of GPS_BornPlus             //
*//   To be used for integrand of FOAM                                              //
*//                                                                                 //
*//   Mass of the final fermion kept exactly.                                       //
*//                                                                                 //
*//   Input:                                                                        //
*//   KFi, Kff = beam and final fermion flavour codes (to define charges)           //
*//   PX       = s-chanel momentum for gamma and Z propagators (not for spinors)    //
*//   Costhe   = cos(theta) for EW formfactors
*//   pi,mi    are for spinors, not for gamma and Z propagators                     //
*//   p1,m1    =fermion momentum and mass (beam)                                    //
*//   p2,m2    =fermion momentum and mass (beam)                                    //
*//   p3,m3    =fermion momentum and mass final state                               //
*//   p4,m4    =fermion momentum and mass final state                               //
*//                                                                                 //
*//   Output:                                                                       //
*//   Xborn     = dsigma/dcosTheta                                                  //
*//   AmpBorn   = spin amplitudes    ???                                            //
*//                                                                                 //
*//   Notes:                                                                        //
*//   Electron mass neglected in spinors, this is why we may use Chisholm!          //
*//   Final fermion mass kept exactly.                                              //
*//   Gamma and Z in s-chanel.                                                      //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'GPS.h'
*
      INTEGER    KFi,KFf,Mode,Level
      DOUBLE PRECISION  PX(4),p1(4),p2(4),p3(4),p4(4)
      DOUBLE PRECISION  m1,m2,m3,m4, CosThe
      DOUBLE COMPLEX    AmpBorn(2,2,2,2)
      DOUBLE COMPLEX    Cfac   ! ???
      DOUBLE PRECISION  Xborn
*
      DOUBLE PRECISION  T3e,Qe
      DOUBLE PRECISION  T3f,Qf
      DOUBLE COMPLEX    Ve,Vf,Ae,Af,VVcor,GamVPi,ZetVPi
      DOUBLE PRECISION  RsqV,RsqA ! QCD corrs.
      INTEGER           NCf,NCe
      DOUBLE PRECISION  svarX
*-----------------------------------------------------------------------------
      INTEGER           i,j,k,l
      INTEGER           j1,j2,j3,j4
      INTEGER           Hel1,Hel2,Hel3,Hel4
      DOUBLE COMPLEX    PropGam,PropZet
      DOUBLE COMPLEX    s31,s24,s14,s32
      DOUBLE COMPLEX    FFacTT(2),      FFacUU(2)
      DOUBLE COMPLEX    SpinoTT(2,2,2,2),SpinoUU(2,2,2,2)
      DOUBLE COMPLEX    TT,UU
      DOUBLE COMPLEX    GPS_iProd1
      DOUBLE COMPLEX    GPS_iProd2
      DOUBLE PRECISION  dummy, BornSum
*-----------------------------------------------------------------------------
      CALL GPS_Initialize
*=============================================================
* Get charges, izospin, color
      CALL BornV_GetParticle(KFi, dummy, Qe,T3e,NCe)
      CALL BornV_GetParticle(KFf, dummy, Qf,T3f,NCf)
*=============================================================
      DO j1 = 1,2
         DO j2 = 1,2
            DO j3 = 1,2
               DO j4 = 1,2
                  Hel1 = 3-2*j1
                  Hel2 = 3-2*j2
                  Hel3 = 3-2*j3
                  Hel4 = 3-2*j4
                  TT  = DCMPLX(0d0,0d0)
                  UU  = DCMPLX(0d0,0d0)
                  IF( Hel2 .EQ. -Hel1) THEN   !!! <--helicity conservation imposed
                     s31 = GPS_iProd2(  Hel3, p3, m3,   Hel1, p1, m1) ! t
                     s24 = GPS_iProd2(  Hel2, p2, m2,   Hel4, p4, m4) ! t1
                     s32 = GPS_iProd2(  Hel3, p3, m3,   Hel2, p2,-m2) ! u1
                     s14 = GPS_iProd2(  Hel1, p1,-m1,   Hel4, p4, m4) ! u
                     TT  = s31*s24
                     UU  = s14*s32
                  ENDIF
                  SpinoTT(j1,j2,j3,j4) =  TT
                  SpinoUU(j1,j2,j3,j4) =  UU
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*////////////////////////////////////////////////////////////////////////////////////////////
*//                        ElectroWeak Corrections                                         //
*////////////////////////////////////////////////////////////////////////////////////////////
      svarX=PX(4)**2-PX(3)**2-PX(2)**2-PX(1)**2
      IF(svarX .LE. (ABS(m3)+ABS(m4))**2 ) RETURN
      CALL GPS_EWFFact(KFi,KFf,SvarX, CosThe ,Ve,Vf,Ae,Af,VVcor,GamVPi,ZetVPi,RsqV,RsqA) !
* Propagators, with s-dependent width
      PropGam =    DCMPLX(  1d0/svarX,  0d0)
      PropZet =    1d0/DCMPLX(svarX-m_MZ**2, m_GammZ*svarX/m_MZ)
* Possibility to switch off Z or gamma, etc.
      IF(m_KeyZet .LE. 0) PropZet =  DCMPLX(0d0)
      IF(m_KeyZet .EQ. 9) PropGam =  DCMPLX(0d0)
      IF(m_KeyZet .EQ.-1) PropZet =  1d0/DCMPLX(SvarX-m_MZ**2, m_GammZ*m_MZ)
*////////////////////////////////////////////////////////////////////////////////////////////
*//     Primitives formfactor-type for construction of spin amplitudes                     //
*//     (Ve -Hel1*Ae)*(Vf +Hel1*Af) is expanded because of double-vector f-factor          //
*////////////////////////////////////////////////////////////////////////////////////////////
      DO j1 = 1,2
         Hel1 = 3-2*j1
         FFacTT(j1) = PropGam*GamVPi *Qe*Qf *RsqV
     $               +PropZet*ZetVPi *(Ve*Vf*VVCor*RsqV -Hel1*Ae*Vf*RsqV +Hel1*Ve*Af*RsqA -Ae*Af*RsqA) !
         FFacUU(j1) = PropGam*GamVPi *Qe*Qf *RsqV
     $               +PropZet*ZetVPi *(Ve*Vf*VVCor*RsqV -Hel1*Ae*Vf*RsqV -Hel1*Ve*Af*RsqA +Ae*Af*RsqA) !
      ENDDO
*////////////////////////////////////////////////////////////////////////////////////////////
*//                     Total result = Spinors*Formfactor                                  //
*////////////////////////////////////////////////////////////////////////////////////////////
      BornSum = 0d0
      DO j1 = 1,2
         DO j2 = 1,2
            DO j3 = 1,2
               DO j4 = 1,2
                  AmpBorn(j1,j2,j3,j4) =  SpinoTT(j1,j2,j3,j4)* FFacTT(j1)
     $                                   +SpinoUU(j1,j2,j3,j4)* FFacUU(j1)
                  BornSum = BornSum +AmpBorn(j1,j2,j3,j4)*DCONJG(AmpBorn(j1,j2,j3,j4))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      Xborn = BornSum
      END                       !!!GPS_BornF!!!


      SUBROUTINE GPS_BornFoam(Mode,KFi,KFf,CMSene,CosThe,Yint)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//                         Born + Boxes                                            //
*//                                                                                 //
*//   It is esentialy a clone of GPS_Born with elements of GPS_BornPlus             //
*//   To be used for integrand of FOAM                                              //
*//                                                                                 //
*//   Masses of the fermion kept exactly.                                           //
*//                                                                                 //
*//   KFi, Kff = beam and final fermion flavour codes (to define charges)           //
*//   PX       = s-chanel momentum for gamma and Z propagators (not for spinors)    //
*//   Costhe   = cos(theta)                                                         //
*//   pi,mi    are for spinors, not for gamma and Z propagators                     //
*//   p1,m1    =fermion momentum and mass (beam)                                    //
*//   p2,m2    =fermion momentum and mass (beam)                                    //
*//   p3,m3    =fermion momentum and mass final state                               //
*//   p4,m4    =fermion momentum and mass final state                               //
*//                                                                                 //
*//   Output:                                                                       //
*//   Yint     = IFI part of YFS formfactor                                         //
*//   m_AmpBorn = Born spin amplitudes                                              //
*//                                                                                 //
*//   Notes:                                                                        //
*//   Electron mass neglected in spinors, this is why we may use Chisholm!          //
*//   Final fermion mass kept exactly.                                              //
*//   Gamma and Z in s-chanel.                                                      //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'GPS.h'
*
      INTEGER           KFi,KFf,Mode,Level
      DOUBLE PRECISION  PX(4),p1(4),p2(4),p3(4),p4(4)
      DOUBLE PRECISION  CMSene, CosThe, m1,m2,m3,m4,meps
      DOUBLE COMPLEX    AmpBorn
      DOUBLE PRECISION  Yint
*
      DOUBLE PRECISION  T3e,Qe
      DOUBLE PRECISION  T3f,Qf
      DOUBLE COMPLEX    Ve,Vf,Ae,Af,VVcor,GamVPi,ZetVPi
      DOUBLE PRECISION  RsqV,RsqA ! QCD corrs.
      INTEGER           NCf,NCe
      DOUBLE PRECISION  svarX
*-----------------------------------------------------------------------------
      INTEGER           i,j,k,l
      INTEGER           j1,j2,j3,j4
      INTEGER           Hel1,Hel2,Hel3,Hel4
      DOUBLE COMPLEX    PropGam,PropZet
      DOUBLE COMPLEX    s31,s24,s14,s32
      DOUBLE COMPLEX    FFacTT(2), FFacUU(2)
      DOUBLE COMPLEX    FFacTG(2), FFacTZ(2), FFacUG(2),FFacUZ(2)
      DOUBLE COMPLEX    SpinoTT(2,2,2,2),SpinoUU(2,2,2,2)
      DOUBLE COMPLEX    BoxGG(2,2,2,2),    BoxGZ(2,2,2,2)
      DOUBLE COMPLEX    TT,UU
      DOUBLE COMPLEX    GPS_iProd1
      DOUBLE COMPLEX    GPS_iProd2
      DOUBLE COMPLEX    AmpBoxy
      DOUBLE COMPLEX    BVR_CBoxGG, BVR_CBoxGZ, BVR_IntIR, BVR_IntReson
      DOUBLE COMPLEX    Coef, IntIR
      DOUBLE PRECISION  dummy, BornSum, Mini, Mfin, Ene, Pini, Pfin, CosThetD
      DOUBLE PRECISION  s,t,u, MasPhot, alfpmix, BVR_TForFac
c[[[[[[[[[[!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DOUBLE PRECISION    icont
      DATA       icont /0d0/
c]]]]]]]]]]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*-----------------------------------------------------------------------------
      CALL GPS_Initialize
*=============================================================
*=============================================================
* Get charges, izospin, color
      CALL BornV_GetParticle(KFi, Mini, Qe,T3e,NCe)
      CALL BornV_GetParticle(KFf, Mfin, Qf,T3f,NCf)
      CALL KK2f_GetMasPhot(MasPhot)
*      Mini =  0.510999e-3  ! electron
*      Mfin =  0.105        ! final ferm. muon
      Ene = CMSene/2
      DO j = 1,4
        P1(j) = 0d0
        P2(j) = 0d0
        P3(j) = 0d0
        P4(j) = 0d0
      ENDDO
      Pini  =  DSQRT( (Ene-Mini)*(Ene+Mini) )
      Pfin  =  DSQRT( (Ene-Mfin)*(Ene+Mfin) )
      P1(4) =  Ene
      P1(3) =  Pini
      P2(4) =  Ene
      P2(3) = -Pini
      P3(4) =  Ene
      P3(3) =  Pfin *CosThe
      P3(1) =  Pfin *DSQRT(1-CosThe**2)
      P4(4) =  Ene
      P4(3) = -P3(3)
      P4(1) = -P3(1)
      meps  = 1d-50
      m1 = Mini
      m2 =-Mini
      m3 = Mfin
      m4 =-Mfin
* YFS formfactor Real+Virtual
* Remember Yint depends on Emin and provides angular asymmetry (MasPhot is dummy)
      alfpmix  = m_Alfpi*Qe*Qf
      Yint  =  BVR_TForFac( alfpmix, p1,Mini, p3,Mfin, Ene, MasPhot)
     $        *BVR_TForFac( alfpmix, p2,Mini, p4,Mfin, Ene, MasPhot)
     $        *BVR_TForFac(-alfpmix, p1,Mini, p4,Mfin, Ene, MasPhot)
     $        *BVR_TForFac(-alfpmix, p2,Mini, p3,Mfin, Ene, MasPhot)
c[[[[[[[[[[[[[!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     icont=icont+1
c      IF(icont.LE.10) THEN
c         write(*,*) ' //////////////////////GPS_BornFoam///////////////////////////////////////'
c         write(*,*) 'GPS_BornFoam: Yint= ',Yint
c      ENDIF
c]]]]]]]]]]]]]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*=============================================================
* Basic spinor products
      DO j1 = 1,2
         DO j2 = 1,2
            DO j3 = 1,2
               DO j4 = 1,2
                  Hel1 = 3-2*j1
                  Hel2 = 3-2*j2
                  Hel3 = 3-2*j3
                  Hel4 = 3-2*j4
                  TT  = DCMPLX(0d0,0d0)
                  UU  = DCMPLX(0d0,0d0)
                  IF( Hel2 .EQ. -Hel1) THEN   !!! <--helicity conservation imposed
                     s31 = GPS_iProd2(  Hel3, p3, m3,   Hel1, p1, m1) ! t
                     s24 = GPS_iProd2(  Hel2, p2, m2,   Hel4, p4, m4) ! t1
                     s32 = GPS_iProd2(  Hel3, p3, m3,   Hel2, p2,-m2) ! u1
                     s14 = GPS_iProd2(  Hel1, p1,-m1,   Hel4, p4, m4) ! u
                     TT  = s31*s24
                     UU  = s14*s32
                  ENDIF
                  SpinoTT(j1,j2,j3,j4) =  TT
                  SpinoUU(j1,j2,j3,j4) =  UU
               ENDDO
            ENDDO
         ENDDO
      ENDDO
****************  Boxes:
* IR-subtracted, helicity conservation imposed (small mass approx.)
      svarX= CMSene**2
      IF(SvarX .LE. (ABS(m3)+ABS(m4))**2 ) RETURN
      s =  SvarX
      CosThetD = CosThe
      t = -s*(1d0-CosThetD)/2d0
      u = -s*(1d0+CosThetD)/2d0
      Coef  = DCMPLX(m_Alfpi*Qe*Qf)
* Virtual 2*(B(t)-B(u)) Intereference IR part to be subtracted from boxes
      IntIR      = Coef*BVR_IntIR(MasPhot,s,t,u)                 !!<- asymetric in (t,u)
      m_IntReson = Coef*BVR_IntReson(MasPhot,m_MZ,m_GammZ,s,t,u) !!<- asymetric in (t,u)
      m_BoxGGtu  = Coef*( BVR_CBoxGG(MasPhot,             s,t,u)) -IntIR
      m_BoxGZtu  = Coef*( BVR_CBoxGZ(MasPhot,m_MZ,m_GammZ,s,t,u)) -IntIR
      m_BoxGGut  = Coef*(-BVR_CBoxGG(MasPhot,             s,u,t)) -IntIR
      m_BoxGZut  = Coef*(-BVR_CBoxGZ(MasPhot,m_MZ,m_GammZ,s,u,t)) -IntIR
* Exponentiate Resonance BigLogs according to Greco et al.
      IF( m_KeyInt .EQ. 2) THEN
         m_BoxGZtu = m_BoxGZtu -m_IntReson
         m_BoxGZut = m_BoxGZut -m_IntReson
      ENDIF
      DO j1 = 1,2
         DO j2 = 1,2
            DO j3 = 1,2
               DO j4 = 1,2
                  Hel1 = 3-2*j1
                  Hel2 = 3-2*j2
                  Hel3 = 3-2*j3
                  Hel4 = 3-2*j4
                  BoxGG(j1,j2,j3,j4) =DCMPLX(0d0,0d0)
                  BoxGZ(j1,j2,j3,j4) =DCMPLX(0d0,0d0)
                  IF((Hel2 .EQ. -Hel1) .AND. (Hel4 .EQ. -Hel3)) THEN !!<--helicity conserv.
                     IF( Hel1*Hel3 .EQ. 1) THEN
                        BoxGG(j1,j2,j3,j4) = m_BoxGGtu
                        BoxGZ(j1,j2,j3,j4) = m_BoxGZtu
                     ELSE
                        BoxGG(j1,j2,j3,j4) = m_BoxGGut
                        BoxGZ(j1,j2,j3,j4) = m_BoxGZut
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
*////////////////////////////////////////////////////////////////////////////////////////////
*//                        ElectroWeak Corrections                                         //
*////////////////////////////////////////////////////////////////////////////////////////////
      svarX= CMSene**2
      IF(svarX .LE. (ABS(m3)+ABS(m4))**2 ) RETURN
      CALL GPS_EWFFact(KFi,KFf,SvarX, CosThe ,Ve,Vf,Ae,Af,VVcor,GamVPi,ZetVPi,RsqV,RsqA) !
* Propagators, with s-dependent width
      PropGam =    DCMPLX(  1d0/svarX,  0d0)
      PropZet =    1d0/DCMPLX(svarX-m_MZ**2, m_GammZ*svarX/m_MZ)
* Possibility to switch off Z or gamma, etc.
      IF(m_KeyZet .LE. 0) PropZet =  DCMPLX(0d0)
      IF(m_KeyZet .EQ. 9) PropGam =  DCMPLX(0d0)
      IF(m_KeyZet .EQ.-1) PropZet =  1d0/DCMPLX(SvarX-m_MZ**2, m_GammZ*m_MZ)
* Exponentiate Resonance BigLogs according to Greco et al.
      IF(  m_KeyINT .EQ. 2 .AND. MODE .GT. 1 ) THEN
          PropZet = PropZet * EXP(m_IntReson)
      ENDIF
*////////////////////////////////////////////////////////////////////////////////////////////
*//     Primitives formfactor-type for construction of spin amplitudes                     //
*//     For boxes we need separately photon and Z parts                                    //
*//     (Ve -Hel1*Ae)*(Vf +Hel1*Af) is expanded because of double-vector f-factor          //
*////////////////////////////////////////////////////////////////////////////////////////////
      DO j1 = 1,2
         Hel1 = 3-2*j1
         FFacTG(j1) = PropGam*GamVPi *Qe*Qf *RsqV
         FFacTZ(j1) = PropZet*ZetVPi *(Ve*Vf*VVCor*RsqV -Hel1*Ae*Vf*RsqV +Hel1*Ve*Af*RsqA -Ae*Af*RsqA) !
         FFacUG(j1) = PropGam*GamVPi *Qe*Qf *RsqV
         FFacUZ(j1) = PropZet*ZetVPi *(Ve*Vf*VVCor*RsqV -Hel1*Ae*Vf*RsqV -Hel1*Ve*Af*RsqA +Ae*Af*RsqA) !
         FFacTT(j1) = FFacTG(j1)+FFacTZ(j1)
         FFacUU(j1) = FFacUG(j1)+FFacUZ(j1)
      ENDDO
*////////////////////////////////////////////////////////////////////////////////////////////
*//                     Total result = Spinors*Formfactor                                  //
*////////////////////////////////////////////////////////////////////////////////////////////
      BornSum = 0d0
      DO j1 = 1,2
         DO j2 = 1,2
            DO j3 = 1,2
               DO j4 = 1,2
* Born,  Zero order
                  AmpBorn = SpinoTT(j1,j2,j3,j4)* FFacTT(j1)
     $                     +SpinoUU(j1,j2,j3,j4)* FFacUU(j1)
* Boxes, First order
                  AmpBoxy = SpinoTT(j1,j2,j3,j4)* FFacTG(j1) *BoxGG(j1,j2,j3,j4)
     $                     +SpinoTT(j1,j2,j3,j4)* FFacTZ(j1) *BoxGZ(j1,j2,j3,j4)
     $                     +SpinoUU(j1,j2,j3,j4)* FFacUG(j1) *BoxGG(j1,j2,j3,j4)
     $                     +SpinoUU(j1,j2,j3,j4)* FFacUZ(j1) *BoxGZ(j1,j2,j3,j4)
*
*////                  AmpBoxy = DCMPLX(0,0)
                  BornSum = BornSum +AmpBorn*DCONJG(AmpBorn)
                  IF(      Mode .EQ.  0 ) THEN
                    m_AmpBorn2(j1,j2,j3,j4) =AmpBorn
                  ELSE IF( Mode .EQ.  1 ) THEN
                    m_AmpBorn1(j1,j2,j3,j4) =AmpBorn
                  ELSE IF( Mode .EQ. 10 ) THEN
                    m_AmpBorn2(j1,j2,j3,j4) =AmpBorn
                  ELSE IF( Mode .EQ. 11 ) THEN
                    m_AmpBorn1(j1,j2,j3,j4) =AmpBorn
                  ELSE IF( Mode .EQ. 20 ) THEN
                    m_AmpBorn2(j1,j2,j3,j4) =AmpBorn+AmpBoxy
                  ELSE IF( Mode .EQ. 21 ) THEN
                    m_AmpBorn1(j1,j2,j3,j4) =AmpBorn+AmpBoxy
                  ELSE
                    WRITE(*,*) "+++ STOP in GPS_BornFoam, wrong Mode= ", Mode
                    STOP
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END                       !!!GPS_BornFoam!!!



      SUBROUTINE GPS_EWFFact(KFi,KFf,Svar,CosThetD,Ve,Vf,Ae,Af,VVcor,GamVPi,ZetVPi,RsqV,RsqA) !
*////////////////////////////////////////////////////////////////////////////////////////////
*//                        ElectroWeak Corrections                                         //    
*//  They are in Vector Couplings (multiplied by correcting f-factors)                     //
*//  Because of cost(theta) depenedence of WW boxes we need to define CosThetD variable    //
*////////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'GPS.h'
* Input
      INTEGER    KFi,KFf
      DOUBLE PRECISION      Svar,CosThetD
* Output
      DOUBLE COMPLEX    Ve,Vf,Ae,Af,VVcor,GamVPi,ZetVPi
      DOUBLE PRECISION  RsqV,RsqA ! QCD corrs.
* Local
      INTEGER           NCf,NCe
      DOUBLE PRECISION  T3e,Qe
      DOUBLE PRECISION  T3f,Qf
      DOUBLE COMPLEX    GSW(100)
      DOUBLE COMPLEX    RhoEW, VPgamma, CorEle, CorFin, CorEleFin, VVCef 
      DOUBLE PRECISION  Deno,  dummy
*===============================================================================
* Get charges, izospin, color
      CALL BornV_GetParticle(KFi, dummy, Qe,T3e,NCe)
      CALL BornV_GetParticle(KFf, dummy, Qf,T3f,NCf)
*
******IF( m_KeyElw .EQ. 0 .OR.  CosThetD .EQ. 0d0 ) THEN   !!! TEST TEST TEST
      IF( m_KeyElw .EQ. 0 ) THEN
* Vacuum polarization factors
         GamVPi = DCMPLX(1d0)
         ZetVPi = DCMPLX(1d0)
         VVCor  = DCMPLX(1d0)
* Couplings costants
         Deno   = DSQRT(16d0*m_Sw2*(1d0-m_Sw2))
         Ve     = (2*T3e -4*Qe*m_Sw2)/Deno
         Vf     = (2*T3f -4*Qf*m_Sw2)/Deno
         Ae     =  2*T3e             /Deno
         Af     =  2*T3f             /Deno
c[[[[[[[!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c         Ve     = 2*T3e             /Deno
c         Vf     = 2*T3e             /Deno
c         Ae     =  0d0
c         Af     =  0d0
c]]]]]]]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         RSQV=1d0
         RSQA=1d0
      ELSE
* Get EW form-factors
         CALL BornV_InterpoGSW(KFf,Svar,CosThetD)
         CALL BornV_GetGSW(GSW)
         CALL BornV_GetQCDcor2(KFf,RSQV,RSQA)
         RhoEW     = GSW(1)
         VPgamma   = GSW(6)
         CorEle    = GSW(2)
         CorFin    = GSW(3)
         CorEleFin = GSW(4)
* Vacuum polarization factors
         GamVPi = 1d0   /(2d0-VPgamma)
         ZetVPi = m_Gmu *m_MZ**2 *m_AlfInv /(DSQRT(2.d0)*8.d0*m_pi)
     $            *(m_Sw2*(1d0-m_Sw2)) *16d0
     $            *RhoEW
* Coupling costants times EW form-factors
         Deno   = DSQRT(16d0*m_Sw2*(1d0-m_Sw2))
         Ve     = (2*T3e -4*Qe*m_Sw2*CorEle)/Deno
         Vf     = (2*T3f -4*Qf*m_Sw2*CorFin)/Deno
         Ae     =  2*T3e             /Deno
         Af     =  2*T3f             /Deno
* Angle dependent double-vector extra-correction
         VVCef  = ( (2*T3e)      *(2*T3f) 
     $             -(4*Qe*m_Sw2) *(2*T3f)      *CorEle 
     $             -(4*Qf*m_Sw2) *(2*T3e)      *CorFin
     $             +(4*Qe*m_Sw2) *(4*Qf*m_Sw2) *CorEleFin )/Deno**2
         VVCor  = VVCef/(Ve*Vf)
* CosThetD = 1d0 is special
         IF( CosThetD .EQ. 0d0) VVCor  = DCMPLX(1d0)
      ENDIF
      END
      SUBROUTINE  GPS_KinExt(p1,m1,p2,m2,p3,m3,p4,m4,ph,mph,phv,p1v,p2v,p3v,p4v)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   kinematical extrapolation of complete to momenum conservation                 //
*//   pv3,pv4 are to replace p3 and p4 fulfilling that                              //
*//  used in reduction procedure for electron neutrino channel                      //
*//   p1,m1,p2,m2,p3,m3,p4,m4,ph,mph INPUT                                          //
*//   pv3,pv4                        OUTPUT                                         //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION      PX(4),p1(4),p2(4),p3(4),p4(4),ph(4),PFAT(4),PSUM(4),PTEST(4),p1v(4),p2v(4),p3v(4),p4v(4),phv(4)
      DOUBLE PRECISION      m1,m2,m3,m4,mph,F0,F1
      INTEGER k
       DO K=1,4
        PFAT(K)=p3(k)+p4(k)+ph(k)
        PTEST(K)=p1(k)+p2(k)-ph(k)
        psum(k)=p3(k)+p4(k)+ph(k)
        p1v(k)=p1(k)
        p2v(k)=p2(k)
        p3v(k)=p3(k)
        p4v(k)=p4(k)
        phv(k)=ph(k)
       ENDDO
!       return
       CALL KinLib_BostQ(1,Psum,ph,phv)
       CALL KinLib_BostQ(1,Psum,p3,p3v)
       CALL KinLib_BostQ(1,Psum,p4,p4v)
       CALL KinLib_BostQ(1,PFAT,ptest,ptest)
       CALL KinLib_BostQ(1,PFAT,psum,psum)
       F0=psum(4)/2
       F1=sqrt((F0**2-m3**2))
       p1v(4)=F0
       p2v(4)=F0
       p1v(3)=F1*p1(3)/abs(p1(3))
       p2v(3)=-p1v(3)
       DO K=1,2
         p1v(k)=0
         p2v(k)=0
       ENDDO


       CALL KinLib_BostQ(-1,PFAT,p1v,p1v)
       CALL KinLib_BostQ(-1,PFAT,p2v,p2v)
       CALL KinLib_BostQ(-1,PFAT,p3v,p3v)
       CALL KinLib_BostQ(-1,PFAT,p4v,p4v)
       CALL KinLib_BostQ(-1,PFAT,phv,phv)

       RETURN
C --   tests are below     
       write(*,*) '3 body kinematics missing momentum visible in last line'
       write(*,*) '--------------'
       write(*,*) p3
       write(*,*) p4
       do k=1,4
         write(*,*) p1(k)+p2(k),ph(k)+p3(k)+p4(k)
       enddo
       write(*,*) ' '
       write(*,*) ' after reparation (dirt moved under beams) '
       write(*,*) phv
       write(*,*) p1v
       write(*,*) p2v
       write(*,*) p3v
       write(*,*) p4v
       do k=1,4
         write(*,*) p1v(k)+p2v(k),phv(k)+p3v(k)+p4v(k)
       enddo
       stop
      END
      SUBROUTINE  GPS_KinExt2(p1,m1,p2,m2,p3,m3,p4,m4,ph1,ph2,mph,ph1v,ph2v,p1v,p2v,p3v,p4v)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   kinematical extrapolation of complete to momenum conservation                 //
*//   pv3,pv4 are to replace p3 and p4 fulfilling that                              //
*//  used in reduction procedure for electron neutrino channel                      //
*//   p1,m1,p2,m2,p3,m3,p4,m4,ph,mph INPUT                                          //
*//   pv3,pv4                        OUTPUT                                         //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION      PX(4),p1(4),p2(4),p3(4),p4(4),ph1(4),PFAT(4),PSUM(4),PTEST(4)
      DOUBLE PRECISION      p1v(4),p2v(4),p3v(4),p4v(4),ph1v(4),ph2v(4),ph2(4)
      DOUBLE PRECISION      m1,m2,m3,m4,mph,F0,F1
      INTEGER k
       DO K=1,4
        PFAT(K)=p3(k)+p4(k)+ph1(k)+ph2(k)
        PTEST(K)=p1(k)+p2(k)-ph1(k)-ph2(k)
        psum(k)=p3(k)+p4(k)+ph1(k)+ph2(k)
        p1v(k)=p1(k)
        p2v(k)=p2(k)
        ph1v(k)=ph1(k)
        ph2v(k)=ph2(k)
       ENDDO
       CALL KinLib_BostQ(1,Psum,ph1,ph1v)
       CALL KinLib_BostQ(1,Psum,ph2,ph2v)
       CALL KinLib_BostQ(1,Psum,p3,p3v)
       CALL KinLib_BostQ(1,Psum,p4,p4v)
       CALL KinLib_BostQ(1,PFAT,ptest,ptest)
       CALL KinLib_BostQ(1,PFAT,psum,psum)
       F0=psum(4)/2
       F1=sqrt((F0**2-m3**2))
       p1v(4)=F0
       p2v(4)=F0
       p1v(3)=F1*p1(3)/abs(p1(3))
       p2v(3)=-p1v(3)
       DO K=1,2
         p1v(k)=0
         p2v(k)=0
       ENDDO


       CALL KinLib_BostQ(-1,PFAT,p1v,p1v)
       CALL KinLib_BostQ(-1,PFAT,p2v,p2v)
       CALL KinLib_BostQ(-1,PFAT,p3v,p3v)
       CALL KinLib_BostQ(-1,PFAT,p4v,p4v)
       CALL KinLib_BostQ(-1,PFAT,ph1v,ph1v)
       CALL KinLib_BostQ(-1,PFAT,ph2v,ph2v)

       RETURN
C --   tests are below     
       write(*,*) '3 body kinematics missing momentum visible in last line'
       write(*,*) '--------------'
       write(*,*) p3
       write(*,*) p4
       do k=1,4
         write(*,*) p1(k)+p2(k),ph1(k)+ph2(k)+p3(k)+p4(k)
       enddo
       write(*,*) ' '
       write(*,*) ' after reparation (dirt moved under beams) '
       write(*,*) ph1v
       write(*,*) ph2v
       write(*,*) p1v
       write(*,*) p2v
       write(*,*) p3v
       write(*,*) p4v
       do k=1,4
         write(*,*) p1v(k)+p2v(k),ph1v(k)+ph2v(k)+p3v(k)+p4v(k)
       enddo
       stop
      END
      SUBROUTINE  GPS_KinExtB(p1,m1,p2,m2,p3,m3,p4,m4,ph,mph,phv,p1v,p2v,p3v,p4v)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   kinematical extrapolation of complete to momenum conservation                 //
*//   pv3,pv4 are to replace p3 and p4 fulfilling that                              //
*//  used in reduction procedure for electron neutrino channel                      //
*//   p1,m1,p2,m2,p3,m3,p4,m4,ph,mph INPUT                                          //
*//   pv3,pv4                        OUTPUT                                         //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE PRECISION      PX(4),p1(4),p2(4),p3(4),p4(4),ph(4),PFAT(4),PSUM(4),PTEST(4),p1v(4),p2v(4),p3v(4),p4v(4),phv(4)
      DOUBLE PRECISION      m1,m2,m3,m4,mph,F0,F1
      INTEGER k
       DO K=1,4
        PFAT(K)=p1(k)+p2(k)-ph(k)
        PTEST(K)=p1(k)+p2(k)-ph(k)
        psum(k)=p3(k)+p4(k) 
        p1v(k)=p1(k)
        p2v(k)=p2(k)
        phv(k)=ph(k)
       ENDDO
       CALL KinLib_BostQ(1,Psum,p3,p3v)
       CALL KinLib_BostQ(1,Psum,p4,p4v)
       CALL KinLib_BostQ(1,PFAT,ptest,ptest)
       CALL KinLib_BostQ(1,PFAT,psum,psum)
       F0=ptest(4)/(p3v(4)+p4v(4))
       F1=sqrt((ptest(4)**2-4*m3**2)/((p3v(4)+p4v(4))**2-4*m3**2))
       p3v(4)=p3v(4)*F0
       p4v(4)=p4v(4)*F0
       DO K=1,3
         p3v(k)=p3v(k)*F1
         p4v(k)=p4v(k)*F1
       ENDDO


       CALL KinLib_BostQ(-1,PFAT,p3v,p3v)
       CALL KinLib_BostQ(-1,PFAT,p4v,p4v)
       RETURN
C --   tests are below     
       write(*,*) '3 body kinematics missing momentum visible in last line'
       write(*,*) '--------------'
       write(*,*) p3
       write(*,*) p4
       do k=1,4
         write(*,*) p1(k)+p2(k),ph(k)+p3(k)+p4(k)
       enddo
       write(*,*) ' '
       write(*,*) ' after reparation (dirt moved under neutrionos) '
       write(*,*) p3v
       write(*,*) p4v
       do k=1,4
         write(*,*) p1(k)+p2(k),ph(k)+p3v(k)+p4v(k)
       enddo
       stop
      END
      SUBROUTINE GPS_HiniPlus(KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,ph,mph,Hel,Sactu,sProd)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   IR-finite part od 1-photon amplitudes for ISR  (equiv. to GPS_Hini)           //
*//   Photon helicity imported from the calling program                             //
*//                                                                                 //
*//   m_AmpExpo*  is working space                                                  //
*//   m_AmpBorn   is hidden INPUT                                                   //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'GPS.h'
*
      INTEGER               KFi,KFf
      DOUBLE PRECISION      PX(4),p1(4),p2(4),p3(4),p4(4),ph(4)
      DOUBLE PRECISION      m1,m2,m3,m4,mph
      DOUBLE COMPLEX        Sactu,sProd,GPS_Sof1,GPS_Sof1b
      INTEGER               Hel
      DOUBLE COMPLEX        AmpBornU(2,2,2,2)
      DOUBLE COMPLEX        AmpBornV(2,2,2,2),AmpBornX(2,2,2,2),AmpBornXm(2,2,2,2)
      DOUBLE COMPLEX        Csum1,Csum2,U(2,2),V(2,2),s1(2),s2(2)
      DOUBLE COMPLEX        AmpExpo1,AmpExpo2,AmpBorn
      INTEGER               j,j1,j2,j3,j4,k,Sig
      DOUBLE PRECISION      pr1,pr2,Fleps
      DOUBLE PRECISION      BornV_GetCharge,Qe
      DOUBLE COMPLEX        Vir1,Vir2
*----------------------------------------
      Fleps =  1d-100
      Qe =  BornV_GetCharge( KFi)

* Virtual corrections
      CALL  BVR_MakeVini(m_Alfpi,Qe,p1,m1,p2,m2,ph, Vir1,Vir2)
***   WRITE(*,*) 'ph/ene,Vir1,Vir2  =',ph(4)/p1(4),Vir1,Vir2

* ISR non-infrared two parts: (1) p1 -> photon, contracted with U-matrix
*                             (2) p2 -> photon, contracted with V-matrix
* Calculate Born spin amplitudes
cc      CALL GPS_Born     (    KFi,KFf,PX,       p1,Fleps,  p2,-Fleps,  p3,m3,   p4,-m4,AmpBornX) !!!!<****
cc      CALL GPS_Born     (    KFi,KFf,PX,       p1,m1,     p2,-m2,     p3,m3,   p4,-m4,AmpBornXm) !!!!<****
      CALL GPS_Born     (    KFi,KFf,PX,       ph,mph,    p2,-Fleps,  p3,m3,   p4,-m4,   AmpBornU)
      CALL GPS_Born(         KFi,KFf,PX,       p1,Fleps,  ph,-mph,    p3,m3,   p4,-m4,   AmpBornV)
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&GPS_HiniPlus&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
c      write(16,'(a,i1,a,i1,a,8g22.11)')(('AmpBornU(',j1,',',j2,',*,*)=',(( AmpBornU(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c      write(16,*) '---------------------------------------------------------------------------------------------'
c      write(16,'(a,i1,a,i1,a,8g22.11)')(('AmpBornV(',j1,',',j2,',*,*)=',(( AmpBornV(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]

* Fermion propagarotors
      pr1 = 1d0/(p1(4)*ph(4)-p1(3)*ph(3)-p1(2)*ph(2)-p1(1)*ph(1))/2d0
      pr2 =-1d0/(p2(4)*ph(4)-p2(3)*ph(3)-p2(2)*ph(2)-p2(1)*ph(1))/2d0
      Sig = 3-2*Hel
      IF( m_KeyArb .EQ. 0 ) THEN 
         CALL GPS_MakeU(ph,Sig,  ph,mph,  p1,m1,    U)
         CALL GPS_MakeV(ph,Sig,  p2,m2,   ph,mph,   V)
      ELSE
         CALL GPS_MakeUb(ph,Sig, ph,mph,  p1,m1,    U)
         CALL GPS_MakeVb(ph,Sig, p2,m2,   ph,mph,   V)
      ENDIF
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&GPS_HiniPlus&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
c      write(16,'(a,8f22.11)') 'U(*,*)= ',((U(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'V(*,*)= ',((V(j1,j2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
*--------------------
      IF( m_KeyArb  .EQ.  0 ) THEN
         s1(1)  = -GPS_Sof1( 1,ph,p1)
         s2(1)  =  GPS_Sof1( 1,ph,p2)
      ELSE
         s1(1)  = -GPS_Sof1b( 1,ph,p1,m1)
         s2(1)  =  GPS_Sof1b( 1,ph,p2,m2)
      ENDIF
      s1(2) = -DCONJG(s1(1))
      s2(2) = -DCONJG(s2(1))
*
c[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      CALL GPS_BornZero(m_AmpExpo1)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
      DO j1=1,2
         DO j2=1,2
            DO j3=1,2
               DO j4=1,2
                  Csum1=DCMPLX(0d0,0d0)
                  Csum2=DCMPLX(0d0,0d0)
                  DO j=1,2
                     Csum1=Csum1 +DCMPLX(Qe *m_e_QED) *U(j,j1)*pr1 *AmpBornU( j,j2,j3,j4)
                     Csum2=Csum2 +DCMPLX(Qe *m_e_QED) *V(j2,j)*pr2 *AmpBornV(j1, j,j3,j4)
                  ENDDO
                  AmpBorn  =  m_AmpBorn(j1,j2,j3,j4) ! this is possibly both small and wrong
                  AmpExpo1 =  sProd/Sactu*(Csum1+Csum2)
**/////////////// under construction
** (1+Vir1)*AmpBorn is already included in AmpExpo0 so we drop it to avoid double counting
** Note that remaining Vir2*AmpBorn is IR-finite because Vir2->0 in the IR limit
                  AmpExpo2 =  
     $                 sProd/Sactu*(Csum1+Csum2) *(1+Vir1+Vir2)*(1+m_F1fin1) ! non-IR sigle bremss. part
     $                       +sProd*AmpBorn*Vir2                             ! add virtual_non_IR*Born
**//////////////////////
                  m_AmpExpo1(j1,j2,j3,j4) =m_AmpExpo1(j1,j2,j3,j4) +AmpExpo1
                  m_AmpExpo2(j1,j2,j3,j4) =m_AmpExpo2(j1,j2,j3,j4) +AmpExpo2
                  m_AmpExpo2p(j1,j2,j3,j4)=m_AmpExpo2p(j1,j2,j3,j4)+AmpExpo2
               ENDDO
            ENDDO
         ENDDO
      ENDDO
c[[[[[[[[[[[[[[[[[[[[[[[*debug*
c      write(16,*) '+++++++++++++++++++++++++++++++++GPS_HiniPlus++++++++++++++++++++++++++++++++++++++++++'
c      write(16,*) 'Vir1,Vir2,m_F1fin1=',Vir1,Vir2,m_F1fin1
c      write(16,*) '----------------------------------------------------------------------------------------'
c      write(16,'(a,i1,a,i1,a,8g22.11)')(('m_AmpExpo2(',j1,',',j2,',*,*)=',(( m_AmpExpo2(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]
      END                       !!! GPS_HiniPlus

      SUBROUTINE GPS_HiniPlusW(Ikey,KFi,KFf,PX,p1,m1,p2,m2,p3,m3,p4,m4,ph,mph,Hel,Sactu,sProd)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   IR-finite part od 1-photon amplitudes for ISR W-exch only (ext. of GPS_Hini)  //
*//   Photon helicity imported from the calling program                             //
*//   Ibeta=1 normal mode of operation                                              //
*//   Ibeta=-1 removes action of Ibeta=1 (for m_AmpExpo2, m_AmpExpo2p)  part one    //
*//   Ibeta=-2 removes action of Ibeta=1 (for m_AmpExpo2, m_AmpExpo2p)  part two    //
*//   Ibeta=0 action at the `fixed transfer level'                                  //
*//                                                                                 //
*//   m_AmpExpo*  is working space                                                  //
*//   m_AmpBorn   is hidden INPUT (defunct)                                         //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'GPS.h'
*
      INTEGER    KFi,KFf    ,NUM,kiju,n,itest
      DOUBLE PRECISION      PX(4),p1(4),p2(4),p3(4),p4(4),ph(4),PTEST(4)
      DOUBLE PRECISION      phv(4),p1v(4),p2v(4),p3v(4),p4v(4)
      DOUBLE PRECISION      m1,m2,m3,m4,mph,CosThetD,s0,t0,u0,sa,ta,ua,sb,tb,ub
      DOUBLE COMPLEX        Sactu,sProd,GPS_Sof1,GPS_Sof1b,GPS_Sof1x,GPS_Sof1bx
      INTEGER               Hel
      DOUBLE COMPLEX        AmpBornU(2,2,2,2)
      DOUBLE COMPLEX        AmpBornV(2,2,2,2),AmpBornW(2,2,2,2),AmpBornWT(2,2,2,2)
      DOUBLE COMPLEX        Csum1,Csum2,Csum3,Csum4,U(2,2),V(2,2),UW(2,2),VW(2,2),s1v(2),s2v(2)
      DOUBLE COMPLEX        UWX(2,2),UWX0(2,2),VWX(2,2),VWX0(2,2)
      DOUBLE COMPLEX        Cnor,UW0(2,2),VW0(2,2),WC(2,2)
      DOUBLE COMPLEX        AmpExpo1,AmpExpo2,AmpBorn
      DOUBLE COMPLEX        PropW0,PropWa,PropWb,WVPi0,WVPia,WVPib,EpsDot(2)
      INTEGER               j,j1,j2,j3,j4,k,Sig,JakKoralZ,n1,n2,n3,n4
      DOUBLE PRECISION      pr1,pr2,pr1v,pr2v,pr1vx,pr2vx,Fleps,F,G,f0,g0,h,h0,R,R0
      DOUBLE PRECISION      BornV_GetCharge,Qe
      DOUBLE COMPLEX        Vir1,Vir2
      LOGICAL               IFONE
      INTEGER Ibeta,Ikey,kk,nn
      DOUBLE PRECISION ssum0,ssuma,ssumb
*----------------------------------------
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&GPS_HiniPlusW&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
c]]]]]]]]]]]]]]]]]]]]]]]]]]]

      IF (ABS(KFf).NE.12) RETURN
!      IF (IBETA.LT.0)  RETURN
      Ibeta=Ikey
      IF(Ikey.eq.-2) Ibeta=-1
      Fleps =  1d-100
      Cnor  =  1D0
!      Ibeta = 1   ! dip switch for t-transfers 1/0/-1 on/off/remove other wrong !!
      Qe =  BornV_GetCharge( KFi)

* Virtual corrections
      CALL  BVR_MakeVini(m_Alfpi,Qe,p1,m1,p2,m2,ph, Vir1,Vir2)
***   WRITE(*,*) 'ph/ene,Vir1,Vir2  =',ph(4)/p1(4),Vir1,Vir2

* ISR non-infrared two parts: (1) p1 -> photon, contracted with U-matrix
*                             (2) p2 -> photon, contracted with V-matrix
* Calculate Born spin amplitudes
****> CALL GPS_Born(KFi,KFf,PX, p1,Mbeam,  p2,-Mbeam,  p3,Massf,p4,-Massf,m_AmpBorn) !!!!<****

        CALL GPS_BornZero(AmpBornU)
        CALL GPS_BornZero(AmpBornV)
        CALL KinLib_ThetaD(PX,p1,p2,p3,p4,s0,CosThetD)
        t0 = -s0*(1d0-CosThetD)/2d0
        u0 = -s0*(1d0+CosThetD)/2d0
        IF ((p3(3)+p4(3))*p1(3).LT.0D0) THEN ! Where is dominant other photon?  It is instead of reduction procedure
         IFONE=.TRUE.      ! We assume all extra photons were emitted from p1 
        ELSE
         IFONE=.FALSE.     ! We assume all extra photons were emitted from p2
        ENDIF
C--------
         IF (IFONE) THEN
           s0=(p3(4)+p4(4))**2-(p3(3)+p4(3))**2-(p3(2)+p4(2))**2-(p3(1)+p4(1))**2
           t0=(p4(4)-p2(4))**2-(p4(3)-p2(3))**2-(p4(2)-p2(2))**2-(p4(1)-p2(1))**2
           u0=(p3(4)-p2(4))**2-(p3(3)-p2(3))**2-(p3(2)-p2(2))**2-(p3(1)-p2(1))**2
         ELSE
           s0=(p3(4)+p4(4))**2-(p3(3)+p4(3))**2-(p3(2)+p4(2))**2-(p3(1)+p4(1))**2
           t0=(p3(4)-p1(4))**2-(p3(3)-p1(3))**2-(p3(2)-p1(2))**2-(p3(1)-p1(1))**2
           u0=(p4(4)-p1(4))**2-(p4(3)-p1(3))**2-(p4(2)-p1(2))**2-(p4(1)-p1(1))**2
         ENDIF
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,*) 'GPS_HiniPlusW: IFONE= ', IFONE, '  s0= ',s0,'  t0= ',t0
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
        CALL GPS_BornZero(AmpBornW)
        CALL GPS_BornWPlus(1,0,KFi,KFf,s0,t0,u0, p1,Fleps,    p2,-Fleps,  p3,m3,   p4,-m4,   AmpBornW)
C --   call on test routine
        CALL GPS_BornWPlusT(KFi,KFf,s0,t0,u0, p1,Fleps,    p2,-Fleps,  p3,m3,   p4,-m4,   AmpBornW,AmpBornWT)
C--------
C reduction procedure .... It is normally not used.
C>>>        CALL GPS_KinExt(p1,m1,p2,m2,p3,m3,p4,m4,ph,mph,phv,p1v,p2v,p3v,p4v)
C ... instead 4-momenta are simply copied.
        DO k=1,4
         phv(k)=ph(k)
         p1v(k)=p1(k)
         p2v(k)=p2(k)
         p3v(k)=p3(k)
         p4v(k)=p4(k)
        ENDDO
        IF ((ph(3)+p3(3)+p4(3))*p1(3).LT.0D0) THEN ! Where is dominant other photon?  It is instead of reduction procedure
         IFONE=.TRUE.      ! We assume all extra photons were emitted from p1 
        ELSE
         IFONE=.FALSE.     ! We assume all extra photons were emitted from p2
        ENDIF
C--------
         IF (IFONE) THEN
           sa=(p3v(4)+p4v(4)+phv(4))**2-(p3v(3)+p4v(3)+phv(3))**2-(p3v(2)+p4v(2)+phv(2))**2-(p3v(1)+p4v(1)+phv(1))**2
           ta=(p4v(4)-p2v(4))**2-(p4v(3)-p2v(3))**2-(p4v(2)-p2v(2))**2-(p4v(1)-p2v(1))**2
           ua=(p3v(4)-p2v(4))**2-(p3v(3)-p2v(3))**2-(p3v(2)-p2v(2))**2-(p3v(1)-p2v(1))**2
         ELSE
           sa=(p3v(4)+p4v(4)+phv(4))**2-(p3v(3)+p4v(3)+phv(3))**2-(p3v(2)+p4v(2)+phv(2))**2-(p3v(1)+p4v(1)+phv(1))**2
           ta=(p3v(4)-p1v(4)+phv(4))**2-(p3v(3)-p1v(3)+phv(3))**2-(p3v(2)-p1v(2)+phv(2))**2-(p3v(1)-p1v(1)+phv(1))**2
           ua=(p4v(4)-p1v(4)+phv(4))**2-(p4v(3)-p1v(3)+phv(3))**2-(p4v(2)-p1v(2)+phv(2))**2-(p4v(1)-p1v(1)+phv(1))**2
         ENDIF
         IF (Ibeta.EQ.0) THEN
            sa=s0
            ta=t0
            ua=t0
         ENDIF
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,*) 'GPS_HiniPlusW:  IFONE= ', IFONE, '  sa= ',sa,'  ta= ',ta, '  Ibeta=',Ibeta
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
         CALL GPS_BornWPlus(1,0,KFi,KFf,sa,ta,ua, ph,mph,    p2,-Fleps,  p3,m3,   p4,-m4,   AmpBornU)

C--------
         IF (IFONE) THEN
           sb=(p3v(4)+p4v(4)+phv(4))**2-(p3v(3)+p4v(3)+phv(3))**2-(p3v(2)+p4v(2)+phv(2))**2-(p3v(1)+p4v(1)+phv(1))**2
           tb=(p4v(4)-p2v(4)+phv(4))**2-(p4v(3)-p2v(3)+phv(3))**2-(p4v(2)-p2v(2)+phv(2))**2-(p4v(1)-p2v(1)+phv(1))**2
           ub=(p3v(4)-p2v(4)+phv(4))**2-(p3v(3)-p2v(3)+phv(3))**2-(p3v(2)-p2v(2)+phv(2))**2-(p3v(1)-p2v(1)+phv(1))**2
         ELSE
           sb=(p3v(4)+p4v(4)+phv(4))**2-(p3v(3)+p4v(3)+phv(3))**2-(p3v(2)+p4v(2)+phv(2))**2-(p3v(1)+p4v(1)+phv(1))**2
           tb=(p3v(4)-p1v(4))**2-(p3v(3)-p1v(3))**2-(p3v(2)-p1v(2))**2-(p3v(1)-p1v(1))**2
           ub=(p4v(4)-p1v(4))**2-(p4v(3)-p1v(3))**2-(p4v(2)-p1v(2))**2-(p4v(1)-p1v(1))**2
         ENDIF
         IF (Ibeta.EQ.0) THEN
            sb=s0
            tb=t0
            ub=t0
         ENDIF
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,*) 'GPS_HiniPlusW: IFONE= ', IFONE, '   sb= ',sb,'  tb= ',tb
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
         CALL GPS_BornWPlus(1,0,KFi,KFf,sb,tb,ub, p1,Fleps,  ph,-mph,    p3,m3,   p4,-m4,   AmpBornV)
         JakKoralZ=0 ! warning: defined in 2 places
         IF(JakKoralZ.eq.1) then
*         W Propagator: it looks crazy in this case ....
          CALL GPS_EWFFactW(KFi,KFf,s0,u0,PropW0,WVPi0)
          CALL GPS_EWFFactW(KFi,KFf,sa,ua,PropWa,WVPia)
          CALL GPS_EWFFactW(KFi,KFf,sb,ub,PropWb,WVPib)
         ELSE
          CALL GPS_EWFFactW(KFi,KFf,s0,t0,PropW0,WVPi0)
          CALL GPS_EWFFactW(KFi,KFf,sa,ta,PropWa,WVPia)
          CALL GPS_EWFFactW(KFi,KFf,sb,tb,PropWb,WVPib)
         ENDIF
         WVPia=WVPi0 ! to keep gauge invariance we install t-transfer in formfactor at 0 order
         WVPib=WVPi0 ! to keep gauge invariance we install t-transfer in formfactor at 0 order
* Fermion propagarotors
      pr1 = 1d0/(p1(4)*ph(4)-p1(3)*ph(3)-p1(2)*ph(2)-p1(1)*ph(1))/2d0
      pr2 =-1d0/(p2(4)*ph(4)-p2(3)*ph(3)-p2(2)*ph(2)-p2(1)*ph(1))/2d0
      pr1v = 1d0/(p1v(4)*phv(4)-p1v(3)*phv(3)-p1v(2)*phv(2)-p1v(1)*phv(1))/2d0
      pr2v =-1d0/(p2v(4)*phv(4)-p2v(3)*phv(3)-p2v(2)*phv(2)-p2v(1)*phv(1))/2d0
      pr1vx= 1d0/(p1v(4)*phv(4))/2d0 * 1D4  ! arbitrary cut off set at 1D4
      pr2vx=-1d0/(p2v(4)*phv(4))/2d0 * 1D4
      Sig = 3-2*Hel
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,*) 'GPS_HiniPlusW: WVPi0= ', WVPi0, '   pr1= ',pr1,'  pr2= ',pr2
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
      IF( m_KeyArb .EQ. 0 ) THEN
         CALL GPS_MakeU(ph,Sig,  ph,mph,  p1,m1,    U)
         CALL GPS_MakeV(ph,Sig,  p2,m2,   ph,mph,   V)
         CALL GPS_MakeUW(Cnor,ph,Sig,  p3,m3,   p1,m1,    UW)
         CALL GPS_MakeVW(Cnor,ph,Sig,  p2,m2,   p4,m4,    VW)
      ELSE
         CALL GPS_MakeUb(ph,Sig, ph,mph,  p1,m1,    U)
         CALL GPS_MakeVb(ph,Sig, p2,m2,   ph,mph,   V)
         CALL GPS_MakeUWb(Cnor,ph,Sig, p3,m3,   p1,m1,    UW)
         CALL GPS_MakeVWb(Cnor,ph,Sig, p2,m2,   p4,m4,    VW)
      ENDIF

         CALL GPS_MakeUX(Cnor,ph,Fleps, p3,m3,   p1,m1,    UWX) ! v-a inside
         CALL GPS_MakeVX(Cnor,ph,Fleps, p2,m2,   p4,m4,    VWX) ! v-a inside
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,'(a,8f22.11)') 'U(*,*)= ',((U(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'V(*,*)= ',((V(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'UW(*,*)= ',((UW(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'VW(*,*)= ',((VW(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'UWX(*,*)= ',((UWX(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'VWX(*,*)= ',((VWX(j1,j2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]

         CALL  GPS_MakeUt(KFf,ph,Sig, p2,m2,   p4,m4, p3,m3,   p1,m1,UW,VW) ! test routine !!
         EpsDot(1)=0d0
*--------------------
      IF( m_KeyArb  .EQ.  0 ) THEN
         s1v(1)  = -GPS_Sof1( 1,phv,p1v)
         s2v(1)  =  GPS_Sof1( 1,phv,p2v)
C>>>           EpsDot(1)=0.5D0*(-GPS_Sof1x( 1,phv,p1v)+GPS_Sof1x( 1,phv,p3v)  ! 0.5 to compensate 2 from Sof1(b)x
C>>>     $                      +GPS_Sof1x( 1,phv,p2v)-GPS_Sof1x( 1,phv,p4v)) ! minis sign is in pr2/4

         IF (IFONE) THEN
           EpsDot(1)=(                                                    ! 0.5 to compensate 2 from Sof1(b)x
     $                      +GPS_Sof1x( 1,phv,p2v)-GPS_Sof1x( 1,phv,p4v)) ! minis sign is in pr2/4
         ELSE
           EpsDot(1)=(-GPS_Sof1x( 1,phv,p1v)+GPS_Sof1x( 1,phv,p3v)        ! 0.5 to compensate 2 from Sof1(b)x
     $                                                                  ) ! minis sign is in pr2/4
         ENDIF
      ELSE
         s1v(1)  = -GPS_Sof1b( 1,phv,p1v,m1)
         s2v(1)  =  GPS_Sof1b( 1,phv,p2v,m2)
C>>>           EpsDot(1)=0.5D0*(-GPS_Sof1bx( 1,phv,p1v,m1)+GPS_Sof1bx( 1,phv,p3v,m3)  ! 0.5 to compensate 2 from Sof1(b)x
C>>>     $                      +GPS_Sof1bx( 1,phv,p2v,m2)-GPS_Sof1bx( 1,phv,p4v,m4)) ! minis sign is in pr2/4

         IF (IFONE) THEN
           EpsDot(1)=(                                                      ! 0.5 to compensate 2 from Sof1(b)x
     $                +GPS_Sof1bx( 1,phv,p2v,m2)-GPS_Sof1bx( 1,phv,p4v,m4)) ! minis sign is in pr2/4
         ELSE
           EpsDot(1)=(-GPS_Sof1bx( 1,phv,p1v,m1)+GPS_Sof1bx( 1,phv,p3v,m3)  ! 0.5 to compensate 2 from Sof1(b)x
     $                                                                   ) ! minis sign is in pr2/4
         ENDIF
      ENDIF
      s1v(2) = -DCONJG(s1v(1))
      s2v(2) = -DCONJG(s2v(1))
      EpsDot(2)=-DCONJG(EpsDot(1))
*
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,'(a,4f22.11)') 's1v(*)= ',(s1v(j1),j1=1,2)
c      write(16,'(a,4f22.11)') 's2v(*)= ',(s2v(j1),j1=1,2)
c      write(16,'(a,4f22.11)') 'EpsDot(*)= ',(EpsDot(j1),j1=1,2)
c[[[[[[[[[[[[[[[[[[[[[[[*debug*
c      write(16,*) '=================================GPS_BornWPlus=========================================='
c      write(16,*) '----------------------------------------------------------------------------------------'
c      write(16,'(a,i1,a,i1,a,8g22.11)')(('m_AmpExpo1(',j1,',',j2,',*,*)=',(( m_AmpExpo1(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c      write(16,'(a,i1,a,i1,a,8g22.11)')(('AmpBornW(',j1,',',j2,',*,*)=',(( AmpBornW(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c      write(16,*) 'PropW0=',PropW0,'WVPi0=',WVPi0
c      write(16,*) 'PropWa=',PropWa,'WVPia=',WVPia
c      write(16,*) 'PropWb=',PropWb,'WVPib=',WVPib
c      write(16,*) 's1v(Hel)= ',s1v(Hel), 's2v(Hel)= ',s2v(Hel),'EpsDot(Hel)= ',EpsDot(Hel)
c      write(16,*) 'm_e_QED=',m_e_QED
c      write(16,*) 'Hel= ',Hel,'  Sig= ',Sig
c]]]]]]]]]]]]]]]]]]]]]]]

      DO j1=1,2
         DO j2=1,2
            DO j3=1,2
               DO j4=1,2
                  Csum1=DCMPLX(0d0,0d0)
                  Csum2=DCMPLX(0d0,0d0)
                  Csum3=DCMPLX(0d0,0d0)
                  Csum4=DCMPLX(0d0,0d0)
                  IF (Ibeta.EQ.0) THEN
                   DO j=1,2
                     Csum1=Csum1 +      DCMPLX(Qe *m_e_QED) *U(j,j1)*pr1 *AmpBornU( j,j2,j3,j4)
                     Csum2=Csum2 +      DCMPLX(Qe *m_e_QED) *V(j2,j)*pr2 *AmpBornV(j1, j,j3,j4)
                   ENDDO
                  ELSE
                   DO j=1,2
                     Csum1=Csum1 +Ibeta*DCMPLX(Qe *m_e_QED) *U(j,j1)*pr1 *AmpBornU( j,j2,j3,j4)
                     Csum2=Csum2 +Ibeta*DCMPLX(Qe *m_e_QED) *V(j2,j)*pr2 *AmpBornV(j1, j,j3,j4)
                   ENDDO
                  ENDIF 
                     Csum4=Csum4 +Ibeta*DCMPLX(    m_e_QED)*PropWa*WVPia*PropWb*WVPib/WVpi0 !denominator for gauge invariance see up
     $                            *2            ! from  feynman diagram
     $                            *(-0.5D0)     ! fixup originating from  test in  GPS_BornWPlusT
     $                            *(UW(j3,j1)*VWX(j2,j4)-VW(j2,j4)*UWX(j3,j1)) ! non-infrared part of emission from W

                     Csum3=Csum3 +Ibeta*DCMPLX(m_e_QED)*AmpBornW(j1,j2,j3,j4)/PropW0/WVPi0*(
     $                      Qe*s1v(Hel) *(PropWa*WVPia-PropW0*WVPi0)                   ! t-channel W-prop variation
     $                    + Qe*s2v(Hel) *(PropWb*WVPib-PropW0*WVPi0)                   ! t-channel W-prop variation
     $                    + EpsDot(Hel)* PropWa*WVPia*PropWb*WVPib/  WVPib             ! basically IR emission from W, reduction procedure used only here
     $                                                                               )
C      include 'GPS-t.h'   !!  printouts for tests

                  IF(Ikey.EQ.-1) THEN ! at present we still need to subtract single brem in case of double M.E. at two levels
!                  csum1=0D0
!                  csum2=0D0
!                  csum3=0D0
!                     csum4=0D0
                  ELSEIF(Ikey.EQ.-2) THEN
                   csum1=0D0
                   csum2=0D0
                   csum3=0D0
!                  csum4=0D0
                  ENDIF


                  AmpBorn  =  m_AmpBorn(j1,j2,j3,j4) ! not used
                  AmpExpo1 =  sProd/Sactu*(Csum1+Csum2+Csum3+Csum4)
**/////////////// under construction
** (1+Vir1)*AmpBorn is already included in AmpExpo0 so we drop it to avoid double counting
** Note that remaining Vir2*AmpBorn is IR-finite because Vir2->0 in the IR limit
                  AmpExpo2 =  
     $                 sProd/Sactu*(Csum1+Csum2+Csum3+Csum4) *(1+Vir1+Vir2)*(1+m_F1fin1) ! non-IR sigle bremss. part
**     $                       +sProd*AmpBorn*Vir2 *0                            ! add virtual_non_IR*Born
**//////////////////////
                  IF (Ibeta.NE.-1) THEN
                    m_AmpExpo1(j1,j2,j3,j4) =m_AmpExpo1(j1,j2,j3,j4) +AmpExpo1
                    m_AmpExpo2(j1,j2,j3,j4) =m_AmpExpo2(j1,j2,j3,j4) +AmpExpo2
                    m_AmpExpo2p(j1,j2,j3,j4)=m_AmpExpo2p(j1,j2,j3,j4)+AmpExpo2
                  ELSE
                    m_AmpExpo2(j1,j2,j3,j4) =m_AmpExpo2(j1,j2,j3,j4) +AmpExpo2/(1+Vir1+Vir2)/(1+m_F1fin1)
                    m_AmpExpo2p(j1,j2,j3,j4)=m_AmpExpo2p(j1,j2,j3,j4)+AmpExpo2/(1+Vir1+Vir2)/(1+m_F1fin1)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
c[[[[[[[[[[[[[[[[[[[[[[[*debug*
c      write(16,*) '=================================GPS_BornWPlus=========================================='
c      write(16,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
c      write(16,'(a,i1,a,i1,a,8g22.11)')(('m_AmpExpo1(',j1,',',j2,',*,*)=',(( m_AmpExpo1(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]
      END                       !!! GPS_HiniPlusW


      SUBROUTINE GPS_HfinPlus(KFi,KFf,PX, p1,m1,p2,m2,p3,m3,p4,m4,ph,mph,Hel,Sactu,sProd,CKine)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   IR-finite part od 1-photon amplitudes for FSR  (equiv. to GPS_HfinPlus)       //
*//   Photon helicity is give by the calling program                                //
*//                                                                                 //
*//   Missing contribution in FSR non-IR part due to  svarX/svarQ                   //
*//   Contribution -svarX/svarQ from HERE cancels exactly with svarX/svarQ in beta0 //
*//                                                                                 //
*//   m_AmpExpo*  is working space                                                  //
*//   m_AmpBorn   is hidden INPUT                                                   //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'GPS.h'
*
      INTEGER    KFi,KFf
      DOUBLE PRECISION      PX(4),p1(4),p2(4),p3(4),p4(4),ph(4)
      DOUBLE PRECISION      m1,m2,m3,m4,mph
*
      DOUBLE COMPLEX        Sactu,sProd,CKine
      INTEGER               Hel
      DOUBLE COMPLEX        AmpBornU(2,2,2,2)
      DOUBLE COMPLEX        AmpBornV(2,2,2,2)
      DOUBLE COMPLEX        Csum1,Csum2,U(2,2),V(2,2)
      INTEGER               j,j1,j2,j3,j4,k,Sig
      DOUBLE PRECISION      pr1,pr2,Fleps
      DOUBLE PRECISION      BornV_GetCharge,Qf
      DOUBLE COMPLEX        Vir1,Vir2
      DOUBLE COMPLEX        AmpExpo1,AmpExpo2,AmpBorn
*----------------------------------------
      Fleps =  1d-100
      Qf =  BornV_GetCharge( KFf)
* Virtual corrections
      CALL  BVR_MakeVfin(m_Alfpi,Qf,p3,m3,p4,m4,ph, Vir1,Vir2)
***      WRITE(*,*) 'ph/ene,Vir1,Vir2 =',ph(4)/p1(4),Vir1,Vir2
* FSR non-infrared two parts: (1) p1 -> photon, contracted with U-matrix
*                             (2) p2 -> photon, contracted with V-matrix
****> CALL GPS_Born(KFi,KFf,PX, p1,Mbeam, p2,-Mbeam,  p3,Massf, p4,-Massf,m_AmpBorn) !!!!<****
      CALL GPS_Born(KFi,KFf,PX, p1,Fleps, p2,-Fleps,  ph,mph,   p4,-m4,   AmpBornU)
      CALL GPS_Born(KFi,KFf,PX, p1,Fleps, p2,-Fleps,  p3,m3,    ph,-mph,  AmpBornV)
* Fermion propagarotors
      pr1 = 1d0/(p3(4)*ph(4)-p3(3)*ph(3)-p3(2)*ph(2)-p3(1)*ph(1))/2d0
      pr2 =-1d0/(p4(4)*ph(4)-p4(3)*ph(3)-p4(2)*ph(2)-p4(1)*ph(1))/2d0
      Sig = 3-2*Hel
      IF( m_KeyArb .EQ. 0 ) THEN
         CALL GPS_MakeU(ph,Sig,    p3,m3,  ph,mph,   U)
         CALL GPS_MakeV(ph,Sig,    ph,mph, p4,m4,    V)
      ELSE
         CALL GPS_MakeUb(ph,Sig,   p3,m3,  ph,mph,   U)
         CALL GPS_MakeVb(ph,Sig,   ph,mph, p4,m4,    V)
      ENDIF
      DO j1=1,2
         DO j2=1,2
            DO j3=1,2
               DO j4=1,2
                  Csum1=DCMPLX(0d0,0d0)
                  Csum2=DCMPLX(0d0,0d0)
                  DO j=1,2
                     Csum1=Csum1 +DCMPLX(Qf *m_e_QED) *U(j3,j)*pr1* AmpBornU(j1,j2, j,j4)
                     Csum2=Csum2 +DCMPLX(Qf *m_e_QED) *V(j,j4)*pr2* AmpBornV(j1,j2,j3, j)
                  ENDDO
                  AmpBorn  = m_AmpBorn(j1,j2,j3,j4)
*///// first order
                  AmpExpo1 =  
     $                 +sProd/Sactu*(Csum1+Csum2)             !! non-IR sigle bremss. part
     $                 -sProd*CKine*AmpBorn +sProd*AmpBorn    !! compensate for (svarX1/svarQ)
*///// second order
** (1+Vir1)*AmpBorn is already included in AmpExpo2 so we drop it to avoid double counting
** the remaining Vir2*AmpBorn is IR-finite because Vir2->0 in the IR limit
                  AmpExpo2 =  
     $                 +sProd/Sactu*(Csum1+Csum2)*(1+Vir1+Vir2)*(1+m_F1ini1) ! non-IR sigle bremss. part
     $                 +sProd*(-CKine*AmpBorn+AmpBorn)*(1+Vir1)*(1+m_F1ini1) ! compensate for (svarX1/svarQ)
     $                 +sProd*AmpBorn*Vir2                                   ! add virtual_non_IR*Born
*///////////////////////
                  m_AmpExpo1(j1,j2,j3,j4) =m_AmpExpo1(j1,j2,j3,j4) +AmpExpo1
                  m_AmpExpo2(j1,j2,j3,j4) =m_AmpExpo2(j1,j2,j3,j4) +AmpExpo2
                  m_AmpExpo2p(j1,j2,j3,j4)=m_AmpExpo2p(j1,j2,j3,j4)+AmpExpo2
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END                       !!! GPS_HfinPlus !!!






      SUBROUTINE GPS_HffPlus(CNorm,KFi,KFf,PX,pA,mA,pB,mB,pC,mC,pD,mD,Hel1,ph1,Hel2,ph2,mph,AmpWork) !
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   Genuine IR-finite non 1-photon amplitudes for FSR-FSR are added to AmpWork    //
*//   Photon helicity imported from the calling program.                            //
*//                                                                                 //
*//   For Y_IR=1 IR-finite 1-photon contributions are calculated here as well.      //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
*
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                                    FSR                                            //
*      //                                                                                   //
*      //                            1                  2                                   //
*      //                            |                  |                                   //
*      //             c              |                  |          d                        //
*      //    u  ------<------OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO----<----- v                 //
*      //                                     |                                             //
*      //                                     |X                                            //
*      //                                     |                                             //
*      //                                                                                   //
*      ///////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'GPS.h'
      INTEGER            KFi,KFf
      DOUBLE PRECISION   PX(4),pA(4),pB(4),pC(4),pD(4),ph1(4),ph2(4),ph(4)
      DOUBLE PRECISION   mA,mB,mC,mD,mph
      DOUBLE COMPLEX     CNorm,sProd
      INTEGER            Hel1, Hel2
      DOUBLE COMPLEX     BornAB1D(2,2,2,2), BornAB2D(2,2,2,2), BornABCD(2,2,2,2)
      DOUBLE COMPLEX     BornABC1(2,2,2,2), BornABC2(2,2,2,2)
      DOUBLE COMPLEX     BornAB12(2,2,2,2), BornAB21(2,2,2,2)
      DOUBLE COMPLEX     AmpWork(2,2,2,2)
      DOUBLE COMPLEX     Uc11(2,2),   V11d(2,2),   U122(2,2),  V221(2,2)
      DOUBLE COMPLEX     Uc22(2,2),   V22d(2,2),   U211(2,2),  V112(2,2)
      DOUBLE COMPLEX     U21c(2,2),   Vd12(2,2),   Uc12(2,2),  V21d(2,2)
      DOUBLE COMPLEX     U12c(2,2),   Vd21(2,2),   Uc21(2,2),  V12d(2,2)
      DOUBLE COMPLEX     U121(2,2),   V121(2,2),   U212(2,2),  V212(2,2)
      DOUBLE COMPLEX     Su1,Su2,Su3
      INTEGER            j,j1,j2,j3,j4,k,Sig,l
      DOUBLE PRECISION   Fprop1, Fprop2
      DOUBLE PRECISION   prC1, prD1, prC2, prD2, prC12, prD12
      DOUBLE PRECISION   BornV_GetCharge, ChaIni,ChaFin
      DOUBLE COMPLEX     GPS_Sof1,GPS_Sof1b
      DOUBLE COMPLEX     sC(2,2),sD(2,2)
      DOUBLE COMPLEX     gF
      INTEGER            Y_IR, N_IR
      DOUBLE PRECISION   PP12(4),PP1(4),PP2(4),QQ(4),SvarQ,SvarX12,SvarX1,SvarX2
*------------------------------------------------------------------------------------------
      Y_IR=1                  ! YES, IR included
      Y_IR=0                  ! No,  IR not included
      N_IR=1-Y_IR
*--------------------
      ChaFin =  BornV_GetCharge( KFf)
      gF = DCMPLX(ChaFin *m_e_QED)
      IF( m_KeyArb  .EQ.  0 ) THEN
         sC(1,1)  =  gF *GPS_Sof1( 1,ph1,pC)
         sC(2,1)  =  gF *GPS_Sof1( 1,ph2,pC)
         sD(1,1)  = -gF *GPS_Sof1( 1,ph1,pD)
         sD(2,1)  = -gF *GPS_Sof1( 1,ph2,pD)
      ELSE
         sC(1,1)  =  gF *GPS_Sof1b( 1,ph1,pC,mC)
         sC(2,1)  =  gF *GPS_Sof1b( 1,ph2,pC,mC)
         sD(1,1)  = -gF *GPS_Sof1b( 1,ph1,pD,mD)
         sD(2,1)  = -gF *GPS_Sof1b( 1,ph2,pD,mD)
      ENDIF
      sC(1,2) = -DCONJG(sC(1,1))
      sC(2,2) = -DCONJG(sC(2,1))
      sD(1,2) = -DCONJG(sD(1,1))
      sD(2,2) = -DCONJG(sD(2,1))
* Calculate Born spin amplitudes, also with substitutions
      CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  pC,   mC,   pD,   -mD, BornABCD) ! Standard
      CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  ph1, mph,   pD,   -mD, BornAB1D) ! C->1
      CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  pC,   mC,   ph1, -mph, BornABC1) ! D->1
      CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  ph2, mph,   pD,   -mD, BornAB2D) ! C->2
      CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  pC,   mC,   ph2, -mph, BornABC2) ! D->2
      CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  ph1, mph,   ph2, -mph, BornAB12) ! C->1,D->2
      CALL GPS_Born(KFi,KFf,PX, pA,mA,  pB,-mB,  ph2, mph,   ph1, -mph, BornAB21) ! C->2,D->1
      DO k=1,4
         PP12(k) = pC(k)+pD(k)+ph1(k)+ph2(k)
         PP1 (k) = pC(k)+pD(k)+ph1(k)
         PP2 (k) = pC(k)+pD(k)+ph2(k)
         QQ(k)   = pC(k)+pD(k)
      ENDDO
      svarX12 = PP12(4)**2 -PP12(3)**2 -PP12(2)**2 -PP12(1)**2
      svarX1  =  PP1(4)**2  -PP1(3)**2  -PP1(2)**2  -PP1(1)**2
      svarX2  =  PP2(4)**2  -PP2(3)**2  -PP2(2)**2  -PP2(1)**2
      svarQ   =   QQ(4)**2   -QQ(3)**2   -QQ(2)**2   -QQ(1)**2
* Fermion propagarotors 1
      prC1=  1d0/(pC(4)*ph1(4)-pC(3)*ph1(3)-pC(2)*ph1(2)-pC(1)*ph1(1))/2d0
      prD1= -1d0/(pD(4)*ph1(4)-pD(3)*ph1(3)-pD(2)*ph1(2)-pD(1)*ph1(1))/2d0
* Fermion propagarotors 2
      prC2=  1d0/(pC(4)*ph2(4)-pC(3)*ph2(3)-pC(2)*ph2(2)-pC(1)*ph2(1))/2d0
      prD2= -1d0/(pD(4)*ph2(4)-pD(3)*ph2(3)-pD(2)*ph2(2)-pD(1)*ph2(1))/2d0
* Double propagators
      prC12= 1d0/( pC(4)*ph1(4)- pC(3)*ph1(3)- pC(2)*ph1(2)- pC(1)*ph1(1)
     $           + pC(4)*ph2(4)- pC(3)*ph2(3)- pC(2)*ph2(2)- pC(1)*ph2(1)
     $           +ph1(4)*ph2(4)-ph1(3)*ph2(3)-ph1(2)*ph2(2)-ph1(1)*ph2(1))/2d0
      prD12=-1d0/( pD(4)*ph1(4)- pD(3)*ph1(3)- pD(2)*ph1(2)- pD(1)*ph1(1)
     $            +pD(4)*ph2(4)- pD(3)*ph2(3)- pD(2)*ph2(2)- pD(1)*ph2(1)
     $           +ph1(4)*ph2(4)-ph1(3)*ph2(3)-ph1(2)*ph2(2)-ph1(1)*ph2(1))/2d0
      Fprop1= (1d0/prC1+1d0/prC2)*prC12 -1d0
      Fprop2= (1d0/prD1+1d0/prD2)*prD12 -1d0
      Sig = 3-2*Hel1
      IF( m_KeyArb  .EQ.  0 ) THEN
* end of line
         CALL GPS_MatrV( gF, ph1,Sig,  ph1,mph, pD,mD,     V11d) ! <1|{1}|D>
         CALL GPS_MatrU( gF, ph1,Sig,  pC,mC,   ph1,mph,   Uc11) ! <C|[1]|1>
* false second
         CALL GPS_MatrV( gF, ph1,Sig,  ph2,mph, pD,mD,     V21d) ! <2|{1}|D>
         CALL GPS_MatrU( gF, ph1,Sig,  pC,mC,   ph2,mph,   Uc12) ! <C|[1]|2>
* reverse order
         CALL GPS_MatrV( gF, ph1,Sig,  pD,mD,   ph2,mph,   Vd12) ! <D|{1}|2>
         CALL GPS_MatrU( gF, ph1,Sig,  ph2,mph, pC,mC,     U21c) ! <2|[1]|C>
* xk-xk term case, ph2 first
         CALL GPS_MatrV( gF, ph1,Sig,  ph1,mph, ph2,mph,   V112) ! <1|{1}|2>
         CALL GPS_MatrU( gF, ph1,Sig,  ph2,mph, ph1,mph,   U211) ! <2|[1]|1>
* xk-xk term case ph2 first
         CALL GPS_MatrV( gF, ph1,Sig,  ph2,mph, ph2,mph,   V212) ! <2|{1}|2>
         CALL GPS_MatrU( gF, ph1,Sig,  ph2,mph, ph2,mph,   U212) ! <2|[1]|2>
      ELSE
         CALL GPS_MatrVb(gF, ph1,Sig,  ph1,mph, pD,mD,     V11d)
         CALL GPS_MatrUb(gF, ph1,Sig,  pC,mC,   ph1,mph,   Uc11)
* falSe second
         CALL GPS_MatrVb(gF, ph1,Sig,  ph2,mph, pD,mD,     V21d)
         CALL GPS_MatrUb(gF, ph1,Sig,  pC,mC,   ph2,mph,   Uc12)
* reverse order
         CALL GPS_MatrVb(gF, ph1,Sig,  pD,mD,   ph2,mph,   Vd12)
         CALL GPS_MatrUb(gF, ph1,Sig,  ph2,mph, pC,mC,     U21c)
* for the case when there was ph2 first xk-xk term
         CALL GPS_MatrVb(gF, ph1,Sig,  ph1,mph, ph2,mph,   V112)
         CALL GPS_MatrUb(gF, ph1,Sig,  ph2,mph, ph1,mph,   U211)
* for the case when there was ph2 first xk-xk term
         CALL GPS_MatrVb(gF, ph1,Sig,  ph2,mph, ph2,mph,   V212)
         CALL GPS_MatrUb(gF, ph1,Sig,  ph2,mph, ph2,mph,   U212)
      ENDIF
      Sig = 3-2*Hel2
      IF( m_KeyArb  .EQ.  0 ) THEN
         CALL GPS_MatrV( gF, ph2,Sig,  ph2,mph,  pD,mD,    V22d) ! <2|{2}|D>
         CALL GPS_MatrU( gF, ph2,Sig,  pC,mC,   ph2,mph,   Uc22) ! <C|[2]|2>
* falSe second
         CALL GPS_MatrV( gF, ph2,Sig,  ph1,mph,  pD,mD,    V12d) ! <1|{2}|D>
         CALL GPS_MatrU( gF, ph2,Sig,  pC,mC,   ph1,mph,   Uc21) ! <C|[2]|1>
* reverse order
         CALL GPS_MatrV( gF, ph2,Sig,  pD,mD,   ph1,mph,   Vd21) ! <D|{2}|1>
         CALL GPS_MatrU( gF, ph2,Sig,  ph1,mph, pC,mC,     U12c) ! <1|[2]|C>
* xk-xk term, ph1 first
         CALL GPS_MatrV( gF, ph2,Sig,  ph2,mph, ph1,mph,   V221) ! <2|{2}|1>
         CALL GPS_MatrU( gF, ph2,Sig,  ph1,mph, ph2,mph,   U122) ! <1|[2]|2>
* xk-xk term, ph1 first 
         CALL GPS_MatrV( gF, ph2,Sig,  ph1,mph, ph1,mph,   V121) ! <1|{2}|1>
         CALL GPS_MatrU( gF, ph2,Sig,  ph1,mph, ph1,mph,   U121) ! <1|[2]|1>
      ELSE
         CALL GPS_MatrVb(gF, ph2,Sig,  ph2,mph,  pD,mD,    V22d)
         CALL GPS_MatrUb(gF, ph2,Sig,  pC,mC,   ph2,mph,   Uc22)
* falSe second
         CALL GPS_MatrVb(gF, ph2,Sig,  ph1,mph,  pD,mD,    V12d)
         CALL GPS_MatrUb(gF, ph2,Sig,  pC,mC,   ph1,mph,   Uc21)
* reverse order
         CALL GPS_MatrVb(gF, ph2,Sig,  pD,mD,   ph1,mph,   Vd21)
         CALL GPS_MatrUb(gF, ph2,Sig,  ph1,mph, pC,mC,     U12c)
* for the case when there was ph1 first xk-xk term
         CALL GPS_MatrVb(gF, ph2,Sig,  ph2,mph, ph1,mph,   V221)
         CALL GPS_MatrUb(gF, ph2,Sig,  ph1,mph, ph2,mph,   U122)
* for the case when there was ph1 first xk-xk term
         CALL GPS_MatrVb(gF, ph2,Sig,  ph1,mph, ph1,mph,   V121)
         CALL GPS_MatrUb(gF, ph2,Sig,  ph1,mph, ph1,mph,   U121)
      ENDIF
      DO j1=1,2
         DO j2=1,2
            DO j3=1,2
               DO j4=1,2
                  Su1 = DCMPLX(0d0,0d0)
                  DO j=1,2
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                     1|            2|                                              //
*      //               c      |    c+m+1    |    c+m+1+2          -d                       //
*      //       u  -----<------S-----<-------U------<-------O-------<------ v               //
*      //                                                   |X                              //
*      ///////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +(sC(1,Hel1)*(prC12-prC2*N_IR))*Uc22(j3,j)*BornAB2D(j1,j2,j,j4) !<c|(1)c[2]2|X|d>
                     Su1=Su1  +sC(1,Hel1)* prC12            *Uc21(j3,j)*BornAB1D(j1,j2,j,j4) !<c|(1)c[2]1|X|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                     2|            1|                                              //
*      //               c      |    c+m+2    |    c+m+1+2          -d                       //
*      //       u  -----<------S-----<-------U------<-------O-------<------ v               //
*      //                                                   |X                              //
*      ///////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1+(sC(2,Hel2)*(prC12-prC1*N_IR))*Uc11(j3,j)*BornAB1D(j1,j2,j,j4) !<c|(2)c[1]1|X|d>
                     Su1=Su1 +sC(2,Hel2)* prC12            *Uc12(j3,j)*BornAB2D(j1,j2,j,j4) !<c|(2)c[1]2|X|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                                    |1             |2                              //
*      //               c          -d+m-1-2  |   -d+m-2     |       -d                      //
*      //       u  -----<------O-----<-------V-----<--------S--------<----- v               //
*      //                     X|                                                            //
*      ///////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +BornABC1(j1,j2,j3,j)*V11d(j,j4)*( sD(2,Hel2)*(prD12-prD1*N_IR))!<c|X|1{1}d(2)|d>
                     Su1=Su1 +BornABC2(j1,j2,j3,j)*V21d(j,j4)*  sD(2,Hel2)* prD12            !<c|X|2{1}d(2)|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                                    |2             |1                              //
*      //               c          -d+m-1-2  |   -d+m-1     |       -d                      //
*      //       u  -----<------O-----<-------V-----<--------S--------<----- v               //
*      //                     X|                                                            //
*      ///////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +BornABC2(j1,j2,j3,j)*V22d(j,j4)*( sD(1,Hel1)*(prD12-prD2*N_IR))!<c|X|2{2}d(1)|d>
                     Su1=Su1 +BornABC1(j1,j2,j3,j)*V12d(j,j4)  *sD(1,Hel1)* prD12            !<c|X|1{2}d(1)|d>
                  ENDDO
                  Su3 = DCMPLX(0d0,0d0)
                  DO j=1,2
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                     2|                            |1                              //
*      //               c      |    c+m+1         c+m+1     |      -d                       //
*      //       u  -----<------U-----<-------O------<-------S-------<------ v               //
*      //                                    |X                                             //
*      ///////////////////////////////////////////////////////////////////////////////////////
                     Su3=Su3 +Uc22(j3,j)*prC2  *BornAB2D(j1,j2,j,j4) *sD(1,Hel1) *Y_IR !<c|[2]c|X|1(1)|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                     1|                            |2                              //
*      //               c      |    c+m+1         c+m+2     |      -d                       //
*      //       u  -----<------U-----<-------O------<-------S-------<------ v               //
*      //                                    |X                                             //
*      ///////////////////////////////////////////////////////////////////////////////////////
                     Su3=Su3 +Uc11(j3,j)*prC1 *BornAB1D(j1,j2,j,j4)  *sD(2,Hel2) *Y_IR !<c|[1]c|X|2(2)|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                     2|                            |1                              //
*      //               c      |    c+m+2         c+m+1     |      -d                       //
*      //       u  -----<------S-----<-------O------<-------U-------<------ v               //
*      //                                    |X                                             //
*      ///////////////////////////////////////////////////////////////////////////////////////
                     Su3=Su3 +sC(2,Hel2) *BornABC1(j1,j2,j3,j) *prD1 *V11d(j,j4) *Y_IR !<c|(2)c|X|1{1}|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                     1|                            |2                              //
*      //               c      |    c+m+1         c+m+2     |      -d                       //
*      //       u  -----<------S-----<-------O------<-------U-------<------ v               //
*      //                                    |X                                             //
*      ///////////////////////////////////////////////////////////////////////////////////////
                     Su3=Su3 +sC(1,Hel1) *BornABC2(j1,j2,j3,j) *prD2 *V22d(j,j4) *Y_IR !<c|(1)c|X|2{2}|d>
                  ENDDO
                  Su2 = DCMPLX(0d0,0d0)
                  DO j=1,2
                     DO l=1,2
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                      |2                           |1                              //
*      //               c      |    c+m+2        -d+m-1     |       -d                      //
*      //       u  -----<------U-----<-------O-----<--------V--------<----- v               //
*      //                                    |X                                             //
*      ///////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Uc22( j3,l)*prC2 *BornAB21(j1,j2,l,j ) *V11d( j,j4)*prD1 !<c|[2]2|X|1{1}|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                      |1                           |2                              //
*      //               c      |    c+m+1        -d+m-2     |       -d                      //
*      //       u  -----<------U-----<-------O-----<--------V--------<----- v               //
*      //                                    |X                                             //
*      ///////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Uc11( j3,l)*prC1 *BornAB12(j1,j2,l,j ) *V22d( j,j4)*prD2 !<c|[1]1|X|2{2}|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                     1|            2|                                              //
*      //               c      |    c+m+1    |    c+m+1+2          -d                       //
*      //       u  -----<------U-----<-------O------<-------V-------<------ v               //
*      //                                                  X|                               //
*      ///////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Uc11( j3,l)*prC1  *U122(l,j)*prC12  *BornAB2D(j1,j2,j ,j4) ! <c|[1]1[2]2|X|d>
                        Su2=Su2 +Uc11( j3,l)*prC1  *U121(l,j)*prC12  *BornAB1D(j1,j2,j ,j4) ! <c|[1]1[2]1|X|d>
                        Su2=Su2 +Uc11( j3,l)*prC1  *U12c(l,j)*prC12  *BornABCD(j1,j2,j ,j4) ! <c|[1]1[2]c|X|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                     2|            1|                                              //
*      //               c      |    c+m+2    |    c+m+1+2          -d                       //
*      //       u  -----<------U-----<-------U------<-------O-------<------ v               //
*      //                                                  X|                               //
*      ///////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Uc22( j3,l)*prC2  *U211(l,j)*prC12  *BornAB1D(j1,j2,j ,j4) ! <c|[2]2[1]1|X|d>
                        Su2=Su2 +Uc22( j3,l)*prC2  *U212(l,j)*prC12  *BornAB2D(j1,j2,j ,j4) ! <c|[2]2[1]2|X|d>
                        Su2=Su2 +Uc22( j3,l)*prC2  *U21c(l,j)*prC12  *BornABCD(j1,j2,j ,j4) ! <c|[2]2[1]c|X|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                                    |2             |1                              //
*      //               c          -d+m-1-2  |   -d+m-1     |       -d                      //
*      //       u  -----<------O-----<-------V-----<--------V--------<----- v               //
*      //                      |X                                                           //
*      ///////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +BornABC2(j1,j2,j3,j ) *V221(j,l )*prD12 *V11d( l,j4)*prD1  ! <c|X|2{2}1{1}|d>
                        Su2=Su2 +BornABC1(j1,j2,j3,j ) *V121(j,l )*prD12 *V11d( l,j4)*prD1  ! <c|X|1{2}1{1}|d>
                        Su2=Su2 +BornABCD(j1,j2,j3,j ) *Vd21(j,l )*prD12 *V11d( l,j4)*prD1  ! <c|X|d{2}1{1}|d>
*      ///////////////////////////////////////////////////////////////////////////////////////
*      //                                    |1             |2                              //
*      //               c          -d+m-1-2  |   -d+m-2     |       -d                      //
*      //       u  -----<------O-----<-------V-----<--------V--------<----- v               //
*      //                      |X                                                           //
*      ///////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +BornABC1(j1,j2,j3,j ) *V112(j,l )*prD12 *V22d( l,j4)*prD2  ! <c|X|1{1}2{2}|d>
                        Su2=Su2 +BornABC2(j1,j2,j3,j ) *V212(j,l )*prD12 *V22d( l,j4)*prD2  ! <c|X|2{1}2{2}|d>
                        Su2=Su2 +BornABCD(j1,j2,j3,j ) *Vd12(j,l )*prD12 *V22d( l,j4)*prD2  ! <c|X|d{1}2{2}|d>
                     ENDDO
                  ENDDO
                  sProd = (sC(1,Hel1)+sD(1,Hel1)) *( sC(2,Hel2)+sD(2,Hel2)) !
                  AmpWork(j1,j2,j3,j4) = AmpWork(j1,j2,j3,j4)
     $                 +CNorm*( Su1 +Su2 +Su3)
     $                 +CNorm*BornABCD(j1,j2,j3,j4)*( sC(1,Hel1)*sC(2,Hel2)*Fprop1   !
     $                                               +sD(1,Hel1)*sD(2,Hel2)*Fprop2 ) !
     $                 +CNorm*BornABCD(j1,j2,j3,j4)*sProd *(1d0 -svarX12/svarQ)  *N_IR !
     $                 +CNorm*BornABCD(j1,j2,j3,j4)*sProd *( svarX1/svarQ -1d0)  *N_IR !
     $                 +CNorm*BornABCD(j1,j2,j3,j4)*sProd *( svarX2/svarQ -1d0)  *N_IR !
     $                 +CNorm*BornABCD(j1,j2,j3,j4)*sProd                        *Y_IR !
               ENDDO
            ENDDO
         ENDDO
      ENDDO


      END                       ! GPS_HffPlus



      SUBROUTINE GPS_HiiPlus(CNorm,KFi,KFf,PX,pA,mA,pB,mB,pC,mC,pD,mD,Hel1,ph1,Hel2,ph2,mph,AmpWork) !
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   Genuine IR-finite non 1-photon amplitudes for ISR-ISR are added to AmpWork    //
*//   That is for dip-switch Y_IR=0.                                                //
*//   For Y_IR=1 IR-finite 1-photon contributions are calculated here as well       //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
*//                                        |                                        //
*//                              1         |          2                             //
*//                              |         |X         |                             //
*//      _       -b              |         |          |          a                  //
*//      v  ------<------OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO----<----- u           //
*//                                                                                 //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'GPS.h'
      INTEGER           KFi,KFf
      DOUBLE PRECISION  PX(4),pA(4),pB(4),pC(4),pD(4),ph1(4),ph2(4),ph(4)
      DOUBLE PRECISION  mA,mB,mC,mD,mph
      DOUBLE COMPLEX    CNorm
      INTEGER           Hel1, Hel2
      DOUBLE COMPLEX    Born1BCD(2,2,2,2), Born2BCD(2,2,2,2), BornABCD(2,2,2,2)
      DOUBLE COMPLEX    BornA1CD(2,2,2,2), BornA2CD(2,2,2,2)
      DOUBLE COMPLEX    Born12CD(2,2,2,2), Born21CD(2,2,2,2)
      DOUBLE COMPLEX    AmpWork(2,2,2,2)
      DOUBLE COMPLEX    U11a(2,2),Vb11(2,2),U221(2,2),V122(2,2)
      DOUBLE COMPLEX    U22a(2,2),Vb22(2,2),U112(2,2),V211(2,2)
      DOUBLE COMPLEX    Ua12(2,2),V21b(2,2),U21a(2,2),Vb12(2,2)
      DOUBLE COMPLEX    Ua21(2,2),V12b(2,2),U12a(2,2),Vb21(2,2)
      DOUBLE COMPLEX    U121(2,2),V121(2,2),U212(2,2),V212(2,2)
      DOUBLE COMPLEX    Su1,Su2,sProd
      INTEGER           j,j1,j2,j3,j4,k,Sig,l
      DOUBLE PRECISION  prA1,prB1,prA2,prB2,prA12,prB12
      DOUBLE PRECISION  BornV_GetCharge,ChaIni
      DOUBLE PRECISION  Fprop1,Fprop2
      DOUBLE COMPLEX    gI
      DOUBLE COMPLEX    GPS_Sof1,GPS_Sof1b
      DOUBLE COMPLEX    sA(2,2),sB(2,2)
      INTEGER           Y_IR,N_IR
*------------------------------------------------------------------------------------------
      Y_IR=1                  ! YES, IR included
      Y_IR=0                  ! No,  IR not included
      N_IR=1-Y_IR
*--------------------
      ChaIni =  BornV_GetCharge( KFi)
      gI = DCMPLX(ChaIni*m_e_QED)
      IF( m_KeyArb  .EQ.  0 ) THEN
         sA(1,1)  = -gI*GPS_Sof1( 1,ph1,pA)
         sA(2,1)  = -gI*GPS_Sof1( 1,ph2,pA)
         sB(1,1)  =  gI*GPS_Sof1( 1,ph1,pB)
         sB(2,1)  =  gI*GPS_Sof1( 1,ph2,pB)
      ELSE
         sA(1,1)  = -gI*GPS_Sof1b( 1,ph1,pA,mA)
         sA(2,1)  = -gI*GPS_Sof1b( 1,ph2,pA,mA)
         sB(1,1)  =  gI*GPS_Sof1b( 1,ph1,pB,mB)
         sB(2,1)  =  gI*GPS_Sof1b( 1,ph2,pB,mB)
      ENDIF
      sA(1,2) = -DCONJG(sA(1,1))
      sA(2,2) = -DCONJG(sA(2,1))
      sB(1,2) = -DCONJG(sB(1,1))
      sB(2,2) = -DCONJG(sB(2,1))
* Calculate Born spin amplitudes
      CALL GPS_Born(KFi,KFf,PX, pA,mA,    pB,-mB,      pC,MC,   pD,-mD,   BornABCD) ! Standard
      CALL GPS_Born(KFi,KFf,PX, ph1,mph,  pB,-mB,      pC,mC,   pD,-mD,   Born1BCD) ! A->1
      CALL GPS_Born(KFi,KFf,PX, pA,mA,    ph1,-mph,    pC,mC,   pD,-mD,   BornA1CD) ! B->1
      CALL GPS_Born(KFi,KFf,PX, ph2,mph,  pB,-mB,      pC,mC,   pD,-mD,   Born2BCD) ! A->2
      CALL GPS_Born(KFi,KFf,PX, pA,mA,    ph2,-mph,    pC,mC,   pD,-mD,   BornA2CD) ! B->2
      CALL GPS_Born(KFi,KFf,PX, ph1,mph,  ph2,-mph,    pC,mC,   pD,-mD,   Born12CD) ! A->1,B->2
      CALL GPS_Born(KFi,KFf,PX, ph2,mph,  ph1,-mph,    pC,mC,   pD,-mD,   Born21CD) ! A->2,B->1
* Fermion propagarotors ini1
      prA1= 1d0/(pA(4)*ph1(4)-pA(3)*ph1(3)-pA(2)*ph1(2)-pA(1)*ph1(1))/2d0
      prB1=-1d0/(pB(4)*ph1(4)-pB(3)*ph1(3)-pB(2)*ph1(2)-pB(1)*ph1(1))/2d0
* Fermion propagarotors ini2
      prA2= 1d0/(pA(4)*ph2(4)-pA(3)*ph2(3)-pA(2)*ph2(2)-pA(1)*ph2(1))/2d0
      prB2=-1d0/(pB(4)*ph2(4)-pB(3)*ph2(3)-pB(2)*ph2(2)-pB(1)*ph2(1))/2d0
* DOUBLE propagators
      prA12= 1d0/( pA(4)*ph1(4)- pA(3)*ph1(3)- pA(2)*ph1(2)- pA(1)*ph1(1)
     $           + pA(4)*ph2(4)- pA(3)*ph2(3)- pA(2)*ph2(2)- pA(1)*ph2(1)
     $           -ph1(4)*ph2(4)+ph1(3)*ph2(3)+ph1(2)*ph2(2)+ph1(1)*ph2(1)
     $           )/2d0
      prB12=-1d0/( pB(4)*ph1(4)- pB(3)*ph1(3)- pB(2)*ph1(2)- pB(1)*ph1(1)
     $            +pB(4)*ph2(4)- pB(3)*ph2(3)- pB(2)*ph2(2)- pB(1)*ph2(1)
     $           -ph1(4)*ph2(4)+ph1(3)*ph2(3)+ph1(2)*ph2(2)+ph1(1)*ph2(1)
     $           )/2d0
      Fprop1=(1d0/prA1+1d0/prA2)*prA12-1d0
      Fprop2=(1d0/prB1+1d0/prB2)*prB12-1d0
      Sig = 3-2*Hel1
      IF( m_KeyArb  .EQ.  0 ) THEN
         CALL GPS_MatrU( gI, ph1,Sig,  ph1,mph, pA,mA,     U11a) ! <1|[1]|a>
         CALL GPS_MatrV( gI, ph1,Sig,  pB,mB,   ph1,mph,   Vb11) ! <b|{1}|1>
* falSe second
         CALL GPS_MatrU( gI, ph1,Sig,  ph2,mph, pA,mA,     U21a) ! <2|[1]|a>
         CALL GPS_MatrV( gI, ph1,Sig,  pB,mB,   ph2,mph,   Vb12) ! <b|{1}|2>
* reverse order
         CALL GPS_MatrU( gI, ph1,Sig,  pA,mA,   ph2,mph,   Ua12) ! <a|[1]|2>
         CALL GPS_MatrV( gI, ph1,Sig,  ph2,mph, pB,mB,     V21b) ! <2|{1}|b>
* for the case when there was ph2 first xk-xk term
         CALL GPS_MatrU( gI, ph1,Sig,  ph1,mph, ph2,mph,   U112) ! <1|[1]|2>
         CALL GPS_MatrV( gI, ph1,Sig,  ph2,mph, ph1,mph,   V211) ! <2|{1}|1>
* for the case when there was ph2 first xk-xk term
         CALL GPS_MatrU( gI, ph1,Sig,  ph2,mph, ph2,mph,   U212) ! <2|[1]|2>
         CALL GPS_MatrV( gI, ph1,Sig,  ph2,mph, ph2,mph,   V212) ! <2|{1}|2>
      ELSE
         CALL GPS_MatrUb(gI, ph1,Sig,  ph1,mph, pA,mA,     U11a)
         CALL GPS_MatrVb(gI, ph1,Sig,  pB,mB,   ph1,mph,   Vb11)
* falSe second
         CALL GPS_MatrUb(gI, ph1,Sig,  ph2,mph, pA,mA,     U21a)
         CALL GPS_MatrVb(gI, ph1,Sig,  pB,mB,   ph2,mph,   Vb12)
* reverse order
         CALL GPS_MatrUb(gI, ph1,Sig,  pA,mA,   ph2,mph,   Ua12)
         CALL GPS_MatrVb(gI, ph1,Sig,  ph2,mph, pB,mB,     V21b)
* for the case when there was ph2 first xk-xk term
         CALL GPS_MatrUb(gI, ph1,Sig,  ph1,mph, ph2,mph,   U112)
         CALL GPS_MatrVb(gI, ph1,Sig,  ph2,mph, ph1,mph,   V211)
* for the case when there was ph2 first xk-xk term
         CALL GPS_MatrUb(gI, ph1,Sig,  ph2,mph, ph2,mph,   U212)
         CALL GPS_MatrVb(gI, ph1,Sig,  ph2,mph, ph2,mph,   V212)
      ENDIF
      Sig = 3-2*Hel2
      IF( m_KeyArb  .EQ.  0 ) THEN
         CALL GPS_MatrU( gI, ph2,Sig,  ph2,mph,  pA,mA,    U22a) ! <2|[2]|a>
         CALL GPS_MatrV( gI, ph2,Sig,  pB,mB,   ph2,mph,   Vb22) ! <b|{2}|2>
* falSe second
         CALL GPS_MatrU( gI, ph2,Sig,  ph1,mph,  pA,mA,    U12a) ! <1|[2]|a>
         CALL GPS_MatrV( gI, ph2,Sig,  pB,mB,   ph1,mph,   Vb21) ! <b|{2}|1>
* reverse order
         CALL GPS_MatrU( gI, ph2,Sig,  pA,mA,   ph1,mph,   Ua21) ! <a|[2]|1>
         CALL GPS_MatrV( gI, ph2,Sig,  ph1,mph, pB,mB,     V12b) ! <1|{2}|b>
* for the case when there was ph1 first xk-xk term
         CALL GPS_MatrU( gI, ph2,Sig,  ph2,mph, ph1,mph,   U221) ! <2|[2]|1>
         CALL GPS_MatrV( gI, ph2,Sig,  ph1,mph, ph2,mph,   V122) ! <1|{2}|2>
c for the case when there was ph1 first xk-xk term
         CALL GPS_MatrU( gI, ph2,Sig,  ph1,mph, ph1,mph,   U121) ! <1|[2]|1>
         CALL GPS_MatrV( gI, ph2,Sig,  ph1,mph, ph1,mph,   V121) ! <1|{2}|1>
      ELSE
         CALL GPS_MatrUb(gI, ph2,Sig,  ph2,mph,  pA,mA,    U22a)
         CALL GPS_MatrVb(gI, ph2,Sig,  pB,mB,   ph2,mph,   Vb22)
c falSe second
         CALL GPS_MatrUb(gI, ph2,Sig,  ph1,mph,  pA,mA,    U12a)
         CALL GPS_MatrVb(gI, ph2,Sig,  pB,mB,   ph1,mph,   Vb21)
c reverse order
         CALL GPS_MatrUb(gI, ph2,Sig,  pA,mA,   ph1,mph,   Ua21)
         CALL GPS_MatrVb(gI, ph2,Sig,  ph1,mph, pB,mB,     V12b)
c for the case when there was ph1 first xk-xk term
         CALL GPS_MatrUb(gI, ph2,Sig,  ph2,mph, ph1,mph,   U221)
         CALL GPS_MatrVb(gI, ph2,Sig,  ph1,mph, ph2,mph,   V122)
c for the case when there was ph1 first xk-xk term
         CALL GPS_MatrUb(gI, ph2,Sig,  ph1,mph, ph1,mph,   U121)
         CALL GPS_MatrVb(gI, ph2,Sig,  ph1,mph, ph1,mph,   V121)
      ENDIF
      DO j1=1,2
         DO j2=1,2
            DO j3=1,2
               DO j4=1,2
                  Su1 = DCMPLX(0d0,0d0)
                  DO j=1,2
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    |                  1               2                         //
*        //                    |X                 |               |                         //
*        //      _       -b    |      b+m-1-2     |      a+m-2    |      a                  //
*        //      v  ------<----O--------<---------U--------<------S------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1+Born1BCD(j,j2,j3,j4) *U11a(j,j1)*(prA12-prA1*N_IR)*sA(2,Hel2) !<b|X|1[1]a(2)|a>
                     Su1=Su1+Born2BCD(j,j2,j3,j4) *U21a(j,j1)*   prA12         *sA(2,Hel2) !<b|X|2[1]a(2)|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    |                  2               1                         //
*        //                    |X                 |               |                         //
*        //      _       -b    |      b+m-1-2     |      a+m-1    |      a                  //
*        //      v  ------<----O--------<---------U-------<-------S------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1+Born2BCD(j,j2,j3,j4) *U22a(j,j1)*(prA12-prA2*N_IR)*sA(1,Hel1) !<b|X|2[2]a(1)|a>
                     Su1=Su1+Born1BCD(j,j2,j3,j4) *U12a(j,j1)* prA12           *sA(1,Hel1) !<b|X|1[2]a(1)|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    2                  1               |                         //
*        //                    |                  |               |X                        //
*        //      _       -b    |     -b+m+2       |    -b+m+1+2   |      a                  //
*        //      v  ------<----S--------<---------V--------<------O------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1+( sB(2,Hel2)*(-prB1*N_IR+prB12))*Vb11(j2,j)*BornA1CD(j1,j,j3,j4)!<b|(2)b[1]1|X|a>
                     Su1=Su1+( sB(2,Hel2))           *prB12  *Vb12(j2,j)*BornA2CD(j1,j,j3,j4)!<b|(2)b[1]2|X|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    1                  2               |                         //
*        //                    |                  |               |X                        //
*        //      _       -b    |     -b+m+1       |    -b+m+1+2   |      a                  //
*        //      v  ------<----S--------<---------V--------<------O------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +(sB(1,Hel1)*(-prB2*N_IR+prB12))*Vb22(j2,j)*BornA2CD(j1,j,j3,j4)!<b|(1)b[2]2|X|a>
                     Su1=Su1 +(sB(1,Hel1))           *prB12  *Vb21(j2,j)*BornA1CD(j1,j,j3,j4)!<b|(1)b[2]1|X|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    2                  |               1                         //
*        //                    |                  |X              |                         //
*        //      _       -b    |     -b+m+2       |     a+m-1     |      a                  //
*        //      v  ------<----S--------<---------O--------<------U------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +sB(2,Hel2)    *Born1BCD(j,j2,j3,j4) *prA1*U11a(j,j1) *Y_IR !<b|(2)2|X|1[1]|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    1                  |               2                         //
*        //                    |                  |X              |                         //
*        //      _       -b    |     -b+m+1       |     a+m-2     |      a                  //
*        //      v  ------<----S--------<---------O--------<------U------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +sB(1,Hel1)    *Born2BCD(j,j2,j3,j4) *prA2*U22a(j,j1) *Y_IR !<b|(1)1|X|2[2]|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    1                  |               2                         //
*        //                    |                  |X              |                         //
*        //      _       -b    |     -b+m+1       |     a+m-2     |      a                  //
*        //      v  ------<----V--------<---------O--------<------S------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +Vb11(j2,j)*prB1 *BornA1CD(j1,j,j3,j4) *sA(2,Hel2)    *Y_IR !<b|[1]1|X|a[2]|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    2                  |               1                         //
*        //                    |                  |X              |                         //
*        //      _       -b    |     -b+m+2       |     a+m-1     |      a                  //
*        //      v  ------<----V--------<---------O--------<------S------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +Vb22(j2,j)*prB2 *BornA2CD(j1,j,j3,j4) *sA(1,Hel1)    *Y_IR !<b|[2]2|X|a(1)|a>
                  ENDDO
                  Su2 = DCMPLX(0d0,0d0)
                  DO j=1,2
                     DO l=1,2
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    2                  |               1                         //
*        //                    |                  |X              |                         //
*        //      _       -b    |     -b+m+2       |     a+m-1     |      a                  //
*        //      v  ------<----*--------<---------O--------<------O------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2  +Vb22(j2,l)*prB2  *Born12CD(j,l,j3,j4 )  *U11a(j,j1)*prA1 ! <b|[2]2|X|1[1]|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    1                  |               2                         //
*        //                    |                  |X              |                         //
*        //      _       -b    |     -b+m+1       |     a+m-2     |      a                  //
*        //      v  ------<----*--------<---------O--------<------O------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2  +Vb11(j2,l)*prB1  *Born21CD(j,l,j3,j4 )  *U22a(j,j1)*prA2 ! <b|[1]1|X|2[2]|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    |                  2               1                         //
*        //                    |X                 |               |                         //
*        //      _       -b    |      b+m-1-2     |      a+m-1    |      a                  //
*        //      v  ------<----O--------<---------O--------<------*------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Born1BCD(j,j2,j3,j4) *U121(j,l)*prA12 *U11a(l,j1)*prA1  ! <b|X|1[2]1(1)|a>
                        Su2=Su2 +Born2BCD(j,j2,j3,j4) *U221(j,l)*prA12 *U11a(l,j1)*prA1  ! <b|X|2[2]1(1)|a>
                        Su2=Su2 -BornABCD(j,j2,j3,j4) *Ua21(j,l)*prA12 *U11a(l,j1)*prA1  ! <b|X|a[2]1(1)|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    |                  1               2                         //
*        //                    |X                 |               |                         //
*        //      _       -b    |      b+m-1-2     |      a+m-2    |      a                  //
*        //      v  ------<----O--------<---------O--------<------*------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Born2BCD(j,j2,j3,j4) *U212(j,l)*prA12  *U22a(l,j1)*prA2 ! <b|X|2[1]2(2)|a>
                        Su2=Su2 +Born1BCD(j,j2,j3,j4) *U112(j,l)*prA12  *U22a(l,j1)*prA2 ! <b|X|1[1]2(2)|a>
                        Su2=Su2 -BornABCD(j,j2,j3,j4) *Ua12(j,l)*prA12  *U22a(l,j1)*prA2 ! <b|X|a[1]2(2)|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    1                  2               |                         //
*        //                    |                  |               |X                        //
*        //      _       -b    |     -b+m+1       |    -b+m+1+2   |      a                  //
*        //      v  ------<----*--------<---------O--------<------O------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Vb11(j2,l)*prB1 *V121(l,j)*prB12 *BornA1CD(j1,j,j3,j4) ! <b|(1)1[2]1|X|a>
                        Su2=Su2 +Vb11(j2,l)*prB1 *V122(l,j)*prB12 *BornA2CD(j1,j,j3,j4) ! <b|(1)1[2]2|X|a>
                        Su2=Su2 -Vb11(j2,l)*prB1 *V12b(l,j)*prB12 *BornABCD(j1,j,j3,j4) ! <b|(1)1[2]b|X|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    2                  1               |                         //
*        //                    |                  |               |X                        //
*        //      _       -b    |     -b+m+2       |    -b+m+1+2   |      a                  //
*        //      v  ------<----*--------<---------O--------<------O------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Vb22(j2,l)*prB2 *V212(l,j)*prB12  *BornA2CD(j1,j,j3,j4) ! <b|(2)2[1]2|X|a>
                        Su2=Su2 +Vb22(j2,l)*prB2 *V211(l,j)*prB12  *BornA1CD(j1,j,j3,j4) ! <b|(2)2[1]1|X|a>
                        Su2=Su2 -Vb22(j2,l)*prB2 *V21b(l,j)*prB12  *BornABCD(j1,j,j3,j4) ! <b|(2)2[2]b|X|a>
                     ENDDO
                  ENDDO
                  sProd = ( sA(1,Hel1)+sB(1,Hel1))*( sA(2,Hel2)+sB(2,Hel2))
                  AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4)
     $                 +CNorm*( Su1+Su2 )
     $                 +CNorm*BornABCD(j1,j2,j3,j4)*( sA(1,Hel1)*sA(2,Hel2)*Fprop1  !
     $                                               +sB(1,Hel1)*sB(2,Hel2)*Fprop2) !
     $                 +CNorm*BornABCD(j1,j2,j3,j4)* sProd                   *Y_IR  !
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END                       ! GPS_HiiPlus    

      SUBROUTINE GPS_HifPlus(CNorm,KFi,KFf,PX,pA,mA,pB,mB,pC,mC,pD,mD,Hel1,ph1,Hel2,ph2,mph,AmpWork) !
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   Genuine IR-finite non 1-photon amplitudes for ISR-FSR are added to AmpWork    //
*//   1-st photon in Initial  state,  symmetrisation 1<-->2 is required!            //
*//                                                                                 //
*//   For Y_IR=1 IR-finite 1-photon contributions are calculated here as well       //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//                              1                2                                 //
*//                              |                |                                 //
*//                              |                |                                 //
*//               c              |    OOOOOOOOOOOOOOOOOO          d                 //
*//     u  -------<------------- | ---OOOOOOOOOOOOOOOOOO----------<----- v          //
*//                              |         |                                        //
*//                              |         |X                                       //
*//                              |         |                                        //
*//      _       -b          OOOOOOOOOOOOOOOOOOOOO                a                 //
*//      v  ------<----------OOOOOOOOOOOOOOOOOOOOO----------------<----- u          //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'GPS.h'
      INTEGER           KFi,KFf
      DOUBLE PRECISION  PX(4),pA(4),pB(4),pC(4),pD(4),ph1(4),ph2(4),ph(4)
      DOUBLE PRECISION  mA,mB,mC,mD,mph
      DOUBLE COMPLEX    CNorm
      INTEGER           Hel1,Hel2
      DOUBLE COMPLEX    Born1BCD(2,2,2,2),BornAB2D(2,2,2,2),BornABCD(2,2,2,2)
      DOUBLE COMPLEX    BornA1CD(2,2,2,2),BornABC2(2,2,2,2)
      DOUBLE COMPLEX    Born1B2D(2,2,2,2),Born1BC2(2,2,2,2)
      DOUBLE COMPLEX    BornA12D(2,2,2,2),BornA1C2(2,2,2,2)
      DOUBLE COMPLEX    AmpWork(2,2,2,2)
      DOUBLE COMPLEX    U11a(2,2),Vb11(2,2),Uc22(2,2),V22d(2,2)
      DOUBLE COMPLEX    Su1,Su2
      DOUBLE COMPLEX    GPS_soft,GPS_softb
      DOUBLE COMPLEX    Sini(2),Sfin(2)
      INTEGER           j,j1,j2,j3,j4,k,Sig,l
      DOUBLE PRECISION  prA1,prB1,prC2,prD2
      DOUBLE PRECISION  BornV_GetCharge,ChaIni,ChaFin
      DOUBLE COMPLEX    gI,gF,sProd
      INTEGER           Y_IR,N_IR
      DOUBLE PRECISION  PP2(4),QQ(4),SvarQ,SvarX2
*----------------------------------------
      Y_IR=1                  ! YES, IR included
      Y_IR=0                  ! No,  IR not included
      N_IR=1-Y_IR
*--------------------
      ChaIni =  BornV_GetCharge( KFi)
      ChaFin =  BornV_GetCharge( KFf)
      gI = DCMPLX(ChaIni *m_e_QED)
      gF = DCMPLX(ChaFin *m_e_QED)
      IF( m_KeyArb  .EQ.  0 ) THEN
         Sini(1)  =  gI *GPS_soft(  1,ph1,pA,pB)
         Sfin(1)  = -gF *GPS_soft(  1,ph2,pC,pD)
      ELSE
         Sini(1)  =  gI *GPS_softb( 1,ph1,pA,mA,pB,mB)
         Sfin(1)  = -gF *GPS_softb( 1,ph2,pC,mC,pD,mD)
      ENDIF
      Sini(2) = -DCONJG(Sini(1))
      Sfin(2) = -DCONJG(Sfin(1))
* Calculate Born spin amplitudes
      CALL GPS_Born(KFi,KFf,PX, pA,mA,    pB,-mB,    pC,MC,    pD,-md,   BornABCD) ! standard
      CALL GPS_Born(KFi,KFf,PX, ph1,mph,  pB,-mB,    pC,mC,    pD,-mD,   Born1BCD) ! A->1
      CALL GPS_Born(KFi,KFf,PX, pA,mA,    ph1,-mph,  pC,mC,    pD,-mD,   BornA1CD) ! B->1
      CALL GPS_Born(KFi,KFf,PX, pA,mA,    pB,-mB,    ph2,mph,  pD,-mD,   BornAB2D) ! C->2
      CALL GPS_Born(KFi,KFf,PX, pA,mA,    pB,-mB,    pC,mC,    ph2,-mph, BornABC2) ! D->2
      CALL GPS_Born(KFi,KFf,PX, ph1,mph,  pB,-mB,    ph2,mph,  pD,-mD,   Born1B2D) ! A->1,C->2
      CALL GPS_Born(KFi,KFf,PX, ph1,mph,  pB,-mB,    pC,mC,    ph2,-mph, Born1BC2) ! A->1,D->2
      CALL GPS_Born(KFi,KFf,PX, pA,mA,    ph1,-mph,  ph2,mph,  pD,-mD,   BornA12D) ! B->1,C->2
      CALL GPS_Born(KFi,KFf,PX, pA,mA,    ph1,-mph,  pC,mC,    ph2,-mph, BornA1C2) ! B->1,D->2
      DO k=1,4
         PP2 (k) = pC(k)+pD(k)+ph2(k)
         QQ(k)   = pC(k)+pD(k)
      ENDDO
      svarX2  =  PP2(4)**2  -PP2(3)**2  -PP2(2)**2  -PP2(1)**2
      svarQ   =   QQ(4)**2   -QQ(3)**2   -QQ(2)**2   -QQ(1)**2
* Fermion propagarotors ini
      prA1=  1d0/(pA(4)*ph1(4)-pA(3)*ph1(3)-pA(2)*ph1(2)-pA(1)*ph1(1))/2d0
      prB1= -1d0/(pB(4)*ph1(4)-pB(3)*ph1(3)-pB(2)*ph1(2)-pB(1)*ph1(1))/2d0
* Fermion propagarotors fin
      prC2=  1d0/(pC(4)*ph2(4)-pC(3)*ph2(3)-pC(2)*ph2(2)-pC(1)*ph2(1))/2d0
      prD2= -1d0/(pD(4)*ph2(4)-pD(3)*ph2(3)-pD(2)*ph2(2)-pD(1)*ph2(1))/2d0
      Sig = 3-2*Hel1
      IF( m_KeyArb  .EQ.  0 ) THEN
         CALL GPS_MatrU( gI, ph1,Sig,  ph1,mph,  pA,mA,    U11a) ! <1|[1]|a>
         CALL GPS_MatrV( gI, ph1,Sig,  pB,mB,   ph1,mph,   Vb11) ! <b|{1}|1>
      ELSE
         CALL GPS_MatrUb(gI, ph1,Sig,  ph1,mph,  pA,mA,    U11a)
         CALL GPS_MatrVb(gI, ph1,Sig,  pB,mB,   ph1,mph,   Vb11)
      ENDIF
      Sig = 3-2*Hel2
      IF( m_KeyArb  .EQ.  0 ) THEN
         CALL GPS_MatrU( gF, ph2,Sig,  pC,mC,   ph2,mph,   Uc22) ! <c|[2]|2>
         CALL GPS_MatrV( gF, ph2,Sig,  ph2,mph,  pD,mD,    V22d) ! <2|{2}|d>
      ELSE
         CALL GPS_MatrUb(gF, ph2,Sig,  pC,mC,   ph2,mph,   Uc22) ! 
         CALL GPS_MatrVb(gF, ph2,Sig, ph2,mph,  pD,mD,     V22d) ! 
      ENDIF
      DO j1=1,2
         DO j2=1,2
            DO j3=1,2
               DO j4=1,2
                  Su1= DCMPLX(0d0,0d0)
                  DO j=1,2
*      /////////////////////////////////////////////////////////////////////////////////////
*      //               c                |2                             d                 //
*      //      u  ------<-------SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS-----<----- v          //
*      //                                        |                                        //
*      //                                        |X             |1                        //
*      //      _       -b                        |     a+m-1    |       a                 //
*      //      v  ------<------------------------O--------------O-------<----- u          //
*      /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +Born1BCD( j,j2,j3,j4) *U11a(j,j1)*prA1 *Sfin(Hel2)*Y_IR !<b|X|1[1]|a><c|(+2)|X|(+2)|d>
*      /////////////////////////////////////////////////////////////////////////////////////
*      //               c                |2                            -d                 //
*      //      u  ------<-------SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS-----<----- v          //
*      //                                        |                                        //
*      //                          |1            |X                                       //
*      //      _       -b          |   -b+m+1    |                      a                 //
*      //      v  ------<----------O-------------O----------------------<----- u          //
*      /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +Vb11(j2,j)*prB1 *BornA1CD(j1,j,j3,j4) *Sfin(Hel2)*Y_IR !<b|[1]1|X|a><c|(+2)|X|(+2)|d>
*      /////////////////////////////////////////////////////////////////////////////////////
*      //                         |2                                                      //
*      //               c         |   c+m+2                            -d                 //
*      //      u  ------<---------O--------------O----------------------<----- v          //
*      //                                        |X                                       //
*      //      _       -b                        |       |1             a                 //
*      //      v  ------<------SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS-----<----- u          //
*      /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +Uc22(j3,j)*prC2 *BornAB2D(j1,j2,j,j4) *Sini(Hel1)*Y_IR !<b|(+1)|X|(+1)|a><c|[2]2|X|d>
*      /////////////////////////////////////////////////////////////////////////////////////
*      //                                                       |2                        //
*      //               c                            -d+m+2     |      -d                 //
*      //      u  ------<------------------------O----------------------<----- v          //
*      //                                        |X                                       //
*      //      _       -b                        |       |1             a                 //
*      //      v  ------<------SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS-----<----- u          //
*      /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +BornABC2(j1,j2,j3,j) *V22d(j,j4)*prD2  *Sini(Hel1)*Y_IR !<b|(+1)|X|(+1)|a><c|X|2[2]d>
                  ENDDO
                  Su2= DCMPLX(0d0,0d0)
                  DO j=1,2
                     DO l=1,2
*      /////////////////////////////////////////////////////////////////////////////////////
*      //                         |2                                                      //
*      //               c         |   c+m+2                             d                 //
*      //      u  ------<---------O--------------O----------------------<----- v          //
*      //                                        |                                        //
*      //                                        |X             |1                        //
*      //      _       -b                        |     a+m-1    |       a                 //
*      //      v  ------<------------------------O--------------O-------<----- u          //
*      /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Uc22(j3,l )*prC2 *Born1B2D(j,j2,l,j4) *U11a(j,j1)*prA1 !<b|X|1[1]|a><c|[2]2|X|d>
*      /////////////////////////////////////////////////////////////////////////////////////
*      //                                                       |2                        //
*      //               c                             -d+m-2    |       d                 //
*      //      u  ------<------------------------O--------------O-------<----- v          //
*      //                                        |                                        //
*      //                                        |X             |1                        //
*      //      _       -b                        |     a+m-1    |       a                 //
*      //      v  ------<------------------------O--------------O-------<----- u          //
*      /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Born1BC2(j,j2,j3,l) *U11a(j,j1)*prA1  *V22d(l,j4)*prD2 !<b|X|1[1]|a><c|X|2{2}|d>
*      /////////////////////////////////////////////////////////////////////////////////////
*      //                       |2                                                        //
*      //               c       |    c+m+2                              d                 //
*      //      u  ------<-------O----------------O----------------------<----- v          //
*      //                                        |                                        //
*      //                       |1               |X                                       //
*      //      _       -b       |   -b+m+1       |                      a                 //
*      //      v  ------<-------O----------------O----------------------<----- u          //
*      /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Vb11(j2,j)*prB1  *Uc22(j3,l )*prC2 *BornA12D(j1,j,l,j4)!<b|{1}1|X|a><c|[2]2|X|d>
*      /////////////////////////////////////////////////////////////////////////////////////
*      //                                                       |2                        //
*      //               c                             -d+m-2    |       d                 //
*      //      u  ------<------------------------O--------------O-------<----- v          //
*      //                                        |                                        //
*      //                       |1               |X                                       //
*      //      _       -b       |   -b+m+1       |                      a                 //
*      //      v  ------<-------O----------------O----------------------<----- u          //
*      /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Vb11(j2,j)*prB1  *BornA1C2(j1,j,j3,l) *V22d(l,j4)*prD2 !<b|{1}1|X|a><c|X|2[2]|d>
                     ENDDO
                  ENDDO
                  sProd = Sini(Hel1)* Sfin(Hel2)
                  AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4) 
     $                 +CNorm*( Su1+Su2 )
     $                 +CNorm*BornABCD(j1,j2,j3,j4)*sProd                      *Y_IR !
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END                       ! GPS_HifPlus



      SUBROUTINE GPS_HiiPlusW(Level,CNorm,KFi,KFf,PX,pAo,mA,pBo,mB,pCo,mC,pDo,mD,Hel1,ph1o,Hel2,ph2o,mph,AmpWork) !
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//   IR-semi-finite (only double IR-excluded) W-exchange amplitudes add to AmpWork //
*//   That is for dip-switch Y_IR=0, Y_IR1=1, other settings not checked.           //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
*//                                        |                                        //
*//                              1         |          2                             //
*//                              |         |X         |                             //
*//      _       -b              |         |          |          a                  //
*//      v  ------<------OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO----<----- u           //
*//                                                                                 //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'GPS.h'
      INTEGER           KFi,KFf
      DOUBLE PRECISION  pAo(4),pBo(4),pCo(4),pDo(4),ph1o(4),ph2o(4)
      DOUBLE PRECISION  PX(4),pA(4),pB(4),pC(4),pD(4),ph1(4),ph2(4),ph(4)
      DOUBLE PRECISION  mA,mB,mC,mD,mph
      DOUBLE COMPLEX    CNorm,Cnor
      INTEGER           Hel1, Hel2
      DOUBLE COMPLEX    Born1BCD(2,2,2,2), Born2BCD(2,2,2,2), BornABCD(2,2,2,2)
      DOUBLE COMPLEX    BornA1CD(2,2,2,2), BornA2CD(2,2,2,2)
      DOUBLE COMPLEX    Born12CD(2,2,2,2), Born21CD(2,2,2,2)
      DOUBLE COMPLEX    AmpWork(2,2,2,2)
      DOUBLE COMPLEX    U11a(2,2),Vb11(2,2),U221(2,2),V122(2,2)
      DOUBLE COMPLEX    U22a(2,2),Vb22(2,2),U112(2,2),V211(2,2)
      DOUBLE COMPLEX    Ua12(2,2),V21b(2,2),U21a(2,2),Vb12(2,2)
      DOUBLE COMPLEX    Ua21(2,2),V12b(2,2),U12a(2,2),Vb21(2,2)
      DOUBLE COMPLEX    U121(2,2),V121(2,2),U212(2,2),V212(2,2)
      DOUBLE COMPLEX    UCAW1(2,2),VBDW1(2,2),UC2W1(2,2),V2DW1(2,2)
      DOUBLE COMPLEX    UCAW2(2,2),VBDW2(2,2),UC1W2(2,2),V1DW2(2,2)
      DOUBLE COMPLEX    UCAWX1(2,2),VBDWX1(2,2),UC2WX1(2,2),V2DWX1(2,2),UC2WXA(2,2)
      DOUBLE COMPLEX    UCAWX2(2,2),VBDWX2(2,2),UC1WX2(2,2),V1DWX2(2,2),UC1WXA(2,2)
      DOUBLE COMPLEX    UCAWXB(2,2),VBDWXA(2,2)
      DOUBLE COMPLEX    UCAWXD(2,2),VBDWXC(2,2)
      DOUBLE COMPLEX    UC2WXB(2,2),V2DWXA(2,2),UC2WXD(2,2),V2DWXC(2,2),V2DWXB(2,2)
      DOUBLE COMPLEX    UC1WXB(2,2),V1DWXA(2,2),UC1WXD(2,2),V1DWXC(2,2),V1DWXB(2,2)
      DOUBLE COMPLEX    Su1,Su2,sProd,eps1(4),eps2(4),eps1D2
      INTEGER           j,j1,j2,j3,j4,k,Sig,l
      DOUBLE PRECISION  prA1,prB1,prA2,prB2,prA12,prB12
      DOUBLE COMPLEX    eps1pA,eps2pA,eps1pB,eps2pB,eps1pC,eps2pC,eps1pD,eps2pD,eps1p2,eps2p1,p2p1
      DOUBLE PRECISION  BornV_GetCharge,ChaIni
      DOUBLE PRECISION  s0,t0,u0,CosThetD
      DOUBLE PRECISION  s,t,u,s1,t1,u1,s2,t2,u2,s12,t12,u12,PSI(4),PSF(4)
      DOUBLE COMPLEX    PropW0,WVPi0,PropW,WVPi,PropW1,WVPi1,PropW2,WVPi2,PropW12,WVPi12,CPF0,CPF,CPF1,CPF2,CPF12
      DOUBLE COMPLEX    Fprop1,Fprop2,EpsDot1(2),EpsDot12(2),EpsDot2(2),EpsDot21(2),Fprop1B,Fprop2B
      DOUBLE COMPLEX    gI
      DOUBLE COMPLEX    GPS_Sof1,GPS_Sof1b,GPS_Sof1x,GPS_Sof1bx
      DOUBLE COMPLEX    sA(2,2),sB(2,2)
      INTEGER           Y_IR,Y_IR1,N_IR,I71,I72,I8,I9,I10,IA,IX,IY,IZ,IT,I9X,I9Y,I9Z,I9T,IVI,IV2,IV1
      INTEGER           I9s1,I9s2,I71b,I72b
      INTEGER           ICNT,nin,I9B,ICOL,ICOL1,NICOL,NICOL1,Level
      LOGICAL           IFONE
      DATA IX,IY,IZ,IT /0,0,0,0/
*------------------------------------------------------------------------------------------
      IF     (Level.EQ.0) THEN ! depending on Level some groups of terms are calculated. only Level=0 is used now.
        ICOL =1
        ICOL1=-1
        NICOL1=-1
        NICOL=-1
      ELSEIF (Level.EQ.1) THEN
        ICOL =0
        ICOL1=1
        NICOL1=0
        NICOL=0
      ELSEIF (Level.EQ.2) THEN
        ICOL =0
        ICOL1=0
        NICOL1=1
        NICOL=0
      ELSE
        ICOL =0
        ICOL1=-1
        NICOL1=-1 
        NICOL=-1
      ENDIF
      IF (ABS(KFf).NE.12) RETURN
!       technical switches for some tests. Must be all 1 or, better know what you do.
!       numerical stability not yet checked, problems like with single brem. possible.
        I71=ICOL1   ! single additional loop  1st photon ir-factor     >  ggchkok
        I72=ICOL1   ! single additional loop  2nd photon ir-factor     >  ggchkok
        I8= ICOL    ! double additional loop                           >  ggchkok 
        IA= ICOL    ! doble emission from single leg, but selfcanc     >  ggchkok
        I9= NICOL   ! no additional loop                               >  #######ggchkok
        I9B=NICOL   ! no additional loop                               >  #######ggchkok
        I9X=ICOL    !NICOL1    ! first finite upperline U, second for cancel      >  ggchkok
        I9Y=ICOL    !NICOL1    ! first finite upperline V, second for cancel      >  ggchkok
        I9Z=ICOL    !NICOL1    ! second finite upperline U, second for cancel     >  ggchkok
        I9T=ICOL    !NICOL1    ! second finite upperline V, second for cancel     >  ggchkok
        IVI=ICOL    !NICOL   ! simplified IV i.e. 'true double infrared'        >  ggchkok
        IV2=ICOL    !ICOL1 ! ICOL1   ! double photon rest (two from 1 leg) k1*k2 ....   >  ggchkok
        IV1=ICOL    !ICOL1 ! ICOL1   ! double photon rest (two from 1 leg) k1*k2 ....   >  ggchkok
        I10=NICOL   ! four boson (rest)                                >  #######ggchkok
        I9s1=ICOL1  ! first soft second gaugeinv nontrivial WWgamma   >  ggchkok
        I9s2=ICOL1  ! second soft first gaugeinv nontrivial WWgamma   >  ggchkok
        I71b=0      !(NICOL)  ! non-ir remnant of I71     >  outed to I9B (ggchkok)
        I72b=0      !(NICOL)  ! non-ir remnant of I72     >  outed to I9B (ggchkok)
      DO K=1,4
        PSF(K)=pCo(K)+pDo(K)+ph1o(K)+ph2o(K)
        PSI(K)=pAo(K)+pBo(K)
      ENDDO

*------------------------------------------------------------------------------------------
      Y_IR=1                   ! YES, IR included
      Y_IR=0                   ! No,  IR not included
      Y_IR1=1                  ! defunct, must be 1 use I9X-T YES, IR times beta 0 included
      N_IR=1-Y_IR1
      N_IR=1     ! #########################################################
*--------------------

      ChaIni =  BornV_GetCharge( KFi)
      gI = DCMPLX(ChaIni*m_e_QED)
        CALL KinLib_ThetaD(PX,pAo,pBo,pCo,pDo,s0,CosThetD)
        t0 = -s0*(1d0-CosThetD)/2d0
        u0 = -s0*(1d0+CosThetD)/2d0

        IF ((pCo(3)+pDo(3))*pAo(3).LT.0D0) THEN ! Where is dominant other photon?  It is instead of reduction procedure
         IFONE=.TRUE.      ! We assume all extra photons were emitted from p1 
        ELSE
         IFONE=.FALSE.     ! We assume all extra photons were emitted from p2
        ENDIF
C--------
         IF (IFONE) THEN
           s0=(pCo(4)+pDo(4))**2-(pCo(3)+pDo(3))**2-(pCo(2)+pDo(2))**2-(pCo(1)+pDo(1))**2
           t0=(pDo(4)-pBo(4))**2-(pDo(3)-pBo(3))**2-(pDo(2)-pBo(2))**2-(pDo(1)-pBo(1))**2
           u0=(pCo(4)-pBo(4))**2-(pCo(3)-pBo(3))**2-(pCo(2)-pBo(2))**2-(pCo(1)-pBo(1))**2
         ELSE
           s0=(pCo(4)+pDo(4))**2-(pCo(3)+pDo(3))**2-(pCo(2)+pDo(2))**2-(pCo(1)+pDo(1))**2
           t0=(pCo(4)-pAo(4))**2-(pCo(3)-pAo(3))**2-(pCo(2)-pAo(2))**2-(pCo(1)-pAo(1))**2
           u0=(pDo(4)-pAo(4))**2-(pDo(3)-pAo(3))**2-(pDo(2)-pAo(2))**2-(pDo(1)-pAo(1))**2
         ENDIF
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,*) '============================================================================================='
c      write(16,*) '=========================================GPS_HiiPlusW========================================'
c      write(16,*) 'GPS_HiniPlusW: IFONE= ', IFONE, '  s0= ',s0,'  t0= ',t0
c]]]]]]]]]]]]]]]]]]]]]]]]]]]

C--------
C reduction procedure (not used anymore) ....
!        CALL GPS_KinExt2(pAo,mA,pBo,mB,pCo,mC,pDo,mD,ph1o,ph2o,mph,ph1,ph2,pA,pB,pC,pD)
C kill reduction procedure at low levels ....
!         IF (Level.LE.1) THEN
           DO K=1,4
            ph1(K)=ph1o(k)
            ph2(K)=ph2o(k)
            pA(K)=pAo(k)
            pB(K)=pBo(k)
            pC(K)=pCo(k)
            pD(K)=pDo(k)
           ENDDO
!         ENDIF
        IF ((ph1(3)+ph2(3)+pC(3)+pD(3))*pA(3).LT.0D0) THEN
         IFONE=.TRUE.      ! We assume all extra photons were emitted from p1 
        ELSE
         IFONE=.FALSE.     ! We assume all extra photons were emitted from p2
        ENDIF

C--------
C transfers
        IF (IFONE) THEN        
        s=(pC(4)+pD(4)+ph1(4)+ph2(4))**2-(pC(3)+pD(3)+ph1(3)+ph2(3))**2
     $   -(pC(2)+pD(2)+ph1(2)+ph2(2))**2-(pC(1)+pD(1)+ph1(1)+ph2(1))**2
        t=(pD(4)-pB(4)+ph1(4)+ph2(4))**2-(pD(3)-pB(3)+ph1(3)+ph2(3))**2
     $   -(pD(2)-pB(2)+ph1(2)+ph2(2))**2-(pD(1)-pB(1)+ph1(1)+ph2(1))**2
        u=(pC(4)-pB(4)+ph1(4)+ph2(4))**2-(pC(3)-pB(3)+ph1(3)+ph2(3))**2
     $   -(pC(2)-pB(2)+ph1(2)+ph2(2))**2-(pC(1)-pB(1)+ph1(1)+ph2(1))**2

        s1=(pC(4)+pD(4)+ph1(4)+ph2(4))**2-(pC(3)+pD(3)+ph1(3)+ph2(3))**2
     $    -(pC(2)+pD(2)+ph1(2)+ph2(2))**2-(pC(1)+pD(1)+ph1(1)+ph2(1))**2
        t1=(pD(4)-pB(4)+ph2(4))**2-(pD(3)-pB(3)+ph2(3))**2-(pD(2)-pB(2)+ph2(2))**2-(pD(1)-pB(1)+ph2(1))**2
        u1=(pC(4)-pB(4)+ph2(4))**2-(pC(3)-pB(3)+ph2(3))**2-(pC(2)-pB(2)+ph2(2))**2-(pC(1)-pB(1)+ph2(1))**2


        s2=(pC(4)+pD(4)+ph1(4)+ph2(4))**2-(pC(3)+pD(3)+ph1(3)+ph2(3))**2
     $    -(pC(2)+pD(2)+ph1(2)+ph2(2))**2-(pC(1)+pD(1)+ph1(1)+ph2(1))**2
        t2=(pD(4)-pB(4)+ph1(4))**2-(pD(3)-pB(3)+ph1(3))**2-(pD(2)-pB(2)+ph1(2))**2-(pD(1)-pB(1)+ph1(1))**2
        u2=(pC(4)-pB(4)+ph1(4))**2-(pC(3)-pB(3)+ph1(3))**2-(pC(2)-pB(2)+ph1(2))**2-(pC(1)-pB(1)+ph1(1))**2


        s12=(pC(4)+pD(4)+ph1(4)+ph2(4))**2-(pC(3)+pD(3)+ph1(3)+ph2(3))**2
     $     -(pC(2)+pD(2)+ph1(2)+ph2(2))**2-(pC(1)+pD(1)+ph1(1)+ph2(1))**2
        t12=(pD(4)-pB(4))**2-(pD(3)-pB(3))**2-(pD(2)-pB(2))**2-(pD(1)-pB(1))**2
        u12=(pC(4)-pB(4))**2-(pC(3)-pB(3))**2-(pC(2)-pB(2))**2-(pC(1)-pB(1))**2

      ELSE
        s=(pC(4)+pD(4)+ph1(4)+ph2(4))**2-(pC(3)+pD(3)+ph1(3)+ph2(3))**2
     $   -(pC(2)+pD(2)+ph1(2)+ph2(2))**2-(pC(1)+pD(1)+ph1(1)+ph2(1))**2
        t=(pC(4)-pA(4))**2-(pC(3)-pA(3))**2-(pC(2)-pA(2))**2-(pC(1)-pA(1))**2
        u=(pD(4)-pA(4))**2-(pD(3)-pA(3))**2-(pD(2)-pA(2))**2-(pD(1)-pA(1))**2

        s1=(pC(4)+pD(4)+ph1(4)+ph2(4))**2-(pC(3)+pD(3)+ph1(3)+ph2(3))**2
     $    -(pC(2)+pD(2)+ph1(2)+ph2(2))**2-(pC(1)+pD(1)+ph1(1)+ph2(1))**2
        t1=(pC(4)-pA(4)+ph1(4))**2-(pC(3)-pA(3)+ph1(3))**2-(pC(2)-pA(2)+ph1(2))**2-(pC(1)-pA(1)+ph1(1))**2
        u1=(pD(4)-pA(4)+ph1(4))**2-(pD(3)-pA(3)+ph1(3))**2-(pD(2)-pA(2)+ph1(2))**2-(pD(1)-pA(1)+ph1(1))**2


        s2=(pC(4)+pD(4)+ph1(4)+ph2(4))**2-(pC(3)+pD(3)+ph1(3)+ph2(3))**2
     $    -(pC(2)+pD(2)+ph1(2)+ph2(2))**2-(pC(1)+pD(1)+ph1(1)+ph2(1))**2
        t2=(pC(4)-pA(4)+ph2(4))**2-(pC(3)-pA(3)+ph2(3))**2-(pC(2)-pA(2)+ph2(2))**2-(pC(1)-pA(1)+ph2(1))**2
        u2=(pD(4)-pA(4)+ph2(4))**2-(pD(3)-pA(3)+ph2(3))**2-(pD(2)-pA(2)+ph2(2))**2-(pD(1)-pA(1)+ph2(1))**2


        s12=(pC(4)+pD(4)+ph1(4)+ph2(4))**2-(pC(3)+pD(3)+ph1(3)+ph2(3))**2
     $     -(pC(2)+pD(2)+ph1(2)+ph2(2))**2-(pC(1)+pD(1)+ph1(1)+ph2(1))**2
        t12=(pC(4)-pA(4)+ph1(4)+ph2(4))**2-(pC(3)-pA(3)+ph1(3)+ph2(3))**2
     $     -(pC(2)-pA(2)+ph1(2)+ph2(2))**2-(pC(1)-pA(1)+ph1(1)+ph2(1))**2
        u12=(pD(4)-pA(4)+ph1(4)+ph2(4))**2-(pD(3)-pA(3)+ph1(3)+ph2(3))**2
     $     -(pD(2)-pA(2)+ph1(2)+ph2(2))**2-(pD(1)-pA(1)+ph1(1)+ph2(1))**2
      ENDIF
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,*) 'GPS_HiniPlusW:  s= ',s,'  t= ',t,' s1= ',s1,'  t1= ',t1,' s2= ',s2,'  t2= ',t2
c      write(16,*) 'GPS_HiniPlusW:  s12= ',s12,'  t12= ',t12
c]]]]]]]]]]]]]]]]]]]]]]]]]]]

C all W propagators ... some simplifications on vac pol set
          CALL GPS_EWFFactW(KFi,KFf,s0,t0,PropW0,WVPi0)
          CALL GPS_EWFFactW(KFi,KFf,s,t,PropW,WVPi)
          CALL GPS_EWFFactW(KFi,KFf,s1,t1,PropW1,WVPi1)
          CALL GPS_EWFFactW(KFi,KFf,s2,t2,PropW2,WVPi2)
          CALL GPS_EWFFactW(KFi,KFf,s12,t12,PropW12,WVPi12)
          WVPi  =1D0            !!!  WVPi0
          WVPi1 =1D0            !!!  WVPi0
          WVPi2 =1D0            !!!  WVPi0
          WVPi12=1D0            !!!  WVPi0
          CPF0 =PropW0 *1D0     !!!   *WVPi0
          CPF  =PropW  *WVPi
          CPF1 =PropW1 *WVPi1
          CPF2 =PropW2 *WVPi2
          CPF12=PropW12*WVPi12
c[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,*) 'PropW0=',PropW0, '   PropW=',PropW
c      write(16,*) 'PropW1=',PropW1,'  PropW2=',PropW2, ' PropW12=',PropW12
c      write(16,*) 'WVPi0=',WVPi0,'   WVPi=',WVPi
c      write(16,*) 'WVPi1=',WVPi1,'  WVPi2=',WVPi2,' WVPi12=',WVPi12
c]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
      IF( m_KeyArb  .EQ.  0 ) THEN
         sA(1,1)  = -gI*GPS_Sof1( 1,ph1,pA)
         sA(2,1)  = -gI*GPS_Sof1( 1,ph2,pA)
         sB(1,1)  =  gI*GPS_Sof1( 1,ph1,pB)
         sB(2,1)  =  gI*GPS_Sof1( 1,ph2,pB)
         IF (IFONE) THEN
           EpsDot1(1) = -gI*(
     $                       +GPS_Sof1x( 1,ph1,pB)-GPS_Sof1x( 1,ph1,pD)-GPS_Sof1x( 1,ph1,ph2))
           EpsDot12(1)= -gI*(
     $                       +GPS_Sof1x( 1,ph1,pB)-GPS_Sof1x( 1,ph1,pD))
           EpsDot2(1) = -gI*(
     $                       +GPS_Sof1x( 1,ph2,pB)-GPS_Sof1x( 1,ph2,pD)-GPS_Sof1x( 1,ph2,ph1))
           EpsDot21(1)= -gI*(
     $                       +GPS_Sof1x( 1,ph2,pB)-GPS_Sof1x( 1,ph2,pD))
         ELSE
           EpsDot1(1) = -gI*(-GPS_Sof1x( 1,ph1,pA)+GPS_Sof1x( 1,ph1,pC)  
     $                       )
           EpsDot12(1)= -gI*(-GPS_Sof1x( 1,ph1,pA)+GPS_Sof1x( 1,ph1,pC)+GPS_Sof1x( 1,ph1,ph2) 
     $                       )
           EpsDot2(1) = -gI*(-GPS_Sof1x( 1,ph2,pA)+GPS_Sof1x( 1,ph2,pC)  
     $                       )
           EpsDot21(1)= -gI*(-GPS_Sof1x( 1,ph2,pA)+GPS_Sof1x( 1,ph2,pC)+GPS_Sof1x( 1,ph2,ph1) 
     $                       )
         ENDIF
      ELSE
         sA(1,1)  = -gI*GPS_Sof1b( 1,ph1,pA,mA)
         sA(2,1)  = -gI*GPS_Sof1b( 1,ph2,pA,mA)
         sB(1,1)  =  gI*GPS_Sof1b( 1,ph1,pB,mB)
         sB(2,1)  =  gI*GPS_Sof1b( 1,ph2,pB,mB)
         IF (IFONE) THEN
           EpsDot1(1) = -gI*(
     $                       +GPS_Sof1bx( 1,ph1,pB,mB)-GPS_Sof1bx( 1,ph1,pD,mD)-GPS_Sof1bx( 1,ph1,ph2,mph))
           EpsDot12(1)= -gI*(
     $                       +GPS_Sof1bx( 1,ph1,pB,mB)-GPS_Sof1bx( 1,ph1,pD,mD))
           EpsDot2(1) = -gI*(
     $                       +GPS_Sof1bx( 1,ph2,pB,mB)-GPS_Sof1bx( 1,ph2,pD,mD)-GPS_Sof1bx( 1,ph2,ph1,mph))
           EpsDot21(1)= -gI*(
     $                       +GPS_Sof1bx( 1,ph2,pB,mB)-GPS_Sof1bx( 1,ph2,pD,mD))
         ELSE
           EpsDot1(1) = -gI*(-GPS_Sof1bx( 1,ph1,pA,mA)+GPS_Sof1bx( 1,ph1,pC,mC)  
     $                       )
           EpsDot12(1)= -gI*(-GPS_Sof1bx( 1,ph1,pA,mA)+GPS_Sof1bx( 1,ph1,pC,mC)+GPS_Sof1bx( 1,ph1,ph2,mph) 
     $                       )
           EpsDot2(1) = -gI*(-GPS_Sof1bx( 1,ph2,pA,mA)+GPS_Sof1bx( 1,ph2,pC,mC)  
     $                       )
           EpsDot21(1)= -gI*(-GPS_Sof1bx( 1,ph2,pA,mA)+GPS_Sof1bx( 1,ph2,pC,mC)+GPS_Sof1bx( 1,ph2,ph1,mph) 
     $                       )
         ENDIF
      ENDIF
      sA(1,2) = -DCONJG(sA(1,1))
      sA(2,2) = -DCONJG(sA(2,1))
      sB(1,2) = -DCONJG(sB(1,1))
      sB(2,2) = -DCONJG(sB(2,1))
      EpsDot1 (2) = -DCONJG(EpsDot1 (1))
      EpsDot12(2) = -DCONJG(EpsDot12(1))
      EpsDot2 (2) = -DCONJG(EpsDot2 (1))
      EpsDot21(2) = -DCONJG(EpsDot21(1))
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,'(a,8f22.11)') 'sA(*,*)= ',((sA(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'sB(*,*)= ',((sB(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,4f22.11)') ' EpsDot1(*,*)= ',(EpsDot1(j1),j1=1,2)
c      write(16,'(a,4f22.11)') ' EpsDot2(*,*)= ',(EpsDot2(j1),j1=1,2)
c      write(16,'(a,4f22.11)') 'EpsDot12(*,*)= ',(EpsDot12(j1),j1=1,2)
c      write(16,'(a,4f22.11)') 'EpsDot21(*,*)= ',(EpsDot21(j1),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
C photon polarization 4-vectors calculated explicitelly ...
      Sig = 3-2*Hel1
      CALL GPS_Make_eps(ph1,Sig,eps1)
      Sig = 3-2*Hel2
      CALL GPS_Make_eps(ph2,Sig,eps2)
      eps1D2=eps1(4)*eps2(4)-eps1(3)*eps2(3)-eps1(2)*eps2(2)-eps1(1)*eps2(1)
      eps1pA=eps1(4)*pA(4)-eps1(3)*pA(3)-eps1(2)*pA(2)-eps1(1)*pA(1)
      eps2pA=eps2(4)*pA(4)-eps2(3)*pA(3)-eps2(2)*pA(2)-eps2(1)*pA(1)
      eps1pB=eps1(4)*pB(4)-eps1(3)*pB(3)-eps1(2)*pB(2)-eps1(1)*pB(1)
      eps2pB=eps2(4)*pB(4)-eps2(3)*pB(3)-eps2(2)*pB(2)-eps2(1)*pB(1)
      eps1pC=eps1(4)*pC(4)-eps1(3)*pC(3)-eps1(2)*pC(2)-eps1(1)*pC(1)
      eps2pC=eps2(4)*pC(4)-eps2(3)*pC(3)-eps2(2)*pC(2)-eps2(1)*pC(1)
      eps1pD=eps1(4)*pD(4)-eps1(3)*pD(3)-eps1(2)*pD(2)-eps1(1)*pD(1)
      eps2pD=eps2(4)*pD(4)-eps2(3)*pD(3)-eps2(2)*pD(2)-eps2(1)*pD(1)
      eps1p2=eps1(4)*ph2(4)-eps1(3)*ph2(3)-eps1(2)*ph2(2)-eps1(1)*ph2(1)
      eps2p1=eps2(4)*ph1(4)-eps2(3)*ph1(3)-eps2(2)*ph1(2)-eps2(1)*ph1(1)
        p2p1=ph2(4)*ph1(4)-ph2(3)*ph1(3)-ph2(2)*ph1(2)-ph2(1)*ph1(1)
c[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,*) 'eps1pA=',eps1pA,' eps2pA=',eps2pA
c      write(16,*) 'eps1pB=',eps1pB,' eps2pB=',eps2pB
c      write(16,*) 'eps1pC=',eps1pC,' eps2pC=',eps2pC
c      write(16,*) 'eps1pD=',eps1pD,' eps2pD=',eps2pD
c      write(16,*) 'eps1p2=',eps1p2,' eps2p1=',eps2p1
c      write(16,*) 'eps1D2=',eps1D2,' p2p1=',p2p1
c]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
* Calculate Born spin amplitudes
      CALL GPS_BornZero(BornABCD)

      CALL GPS_BornZero(Born1BCD)
      CALL GPS_BornZero(BornA1CD)

      CALL GPS_BornZero(Born2BCD)
      CALL GPS_BornZero(BornA2CD)

      CALL GPS_BornZero(Born12CD)
      CALL GPS_BornZero(Born21CD)

      CALL GPS_BornWPlus(1,0,KFi,KFf,s0,t0,u0, pA,mA,    pB,-mB,      pC,MC,   pD,-mD,   BornABCD) ! Standard
      CALL GPS_BornWPlus(1,0,KFi,KFf,s0,t0,u0, ph1,mph,  pB,-mB,      pC,mC,   pD,-mD,   Born1BCD) ! A->1
      CALL GPS_BornWPlus(1,0,KFi,KFf,s0,t0,u0, pA,mA,    ph1,-mph,    pC,mC,   pD,-mD,   BornA1CD) ! B->1
      CALL GPS_BornWPlus(1,0,KFi,KFf,s0,t0,u0, ph2,mph,  pB,-mB,      pC,mC,   pD,-mD,   Born2BCD) ! A->2
      CALL GPS_BornWPlus(1,0,KFi,KFf,s0,t0,u0, pA,mA,    ph2,-mph,    pC,mC,   pD,-mD,   BornA2CD) ! B->2
      CALL GPS_BornWPlus(1,0,KFi,KFf,s0,t0,u0, ph1,mph,  ph2,-mph,    pC,mC,   pD,-mD,   Born12CD) ! A->1,B->2
      CALL GPS_BornWPlus(1,0,KFi,KFf,s0,t0,u0, ph2,mph,  ph1,-mph,    pC,mC,   pD,-mD,   Born21CD) ! A->2,B->1
* Calculate Born spin amplitudes
C      CALL GPS_Born(KFi,KFf,PX, pA,mA,    pB,-mB,      pC,MC,   pD,-mD,   BornABCD) ! Standard
C      CALL GPS_Born(KFi,KFf,PX, ph1,mph,  pB,-mB,      pC,mC,   pD,-mD,   Born1BCD) ! A->1
C      CALL GPS_Born(KFi,KFf,PX, pA,mA,    ph1,-mph,    pC,mC,   pD,-mD,   BornA1CD) ! B->1
C      CALL GPS_Born(KFi,KFf,PX, ph2,mph,  pB,-mB,      pC,mC,   pD,-mD,   Born2BCD) ! A->2
C      CALL GPS_Born(KFi,KFf,PX, pA,mA,    ph2,-mph,    pC,mC,   pD,-mD,   BornA2CD) ! B->2
C      CALL GPS_Born(KFi,KFf,PX, ph1,mph,  ph2,-mph,    pC,mC,   pD,-mD,   Born12CD) ! A->1,B->2
C      CALL GPS_Born(KFi,KFf,PX, ph2,mph,  ph1,-mph,    pC,mC,   pD,-mD,   Born21CD) ! A->2,B->1
c[[[[[[[[[[[[[[[[[[[[[[[*debug*
c      write(16,*) '----------------------------------------------------------------------------------------'
c      write(16,*) ' s0=',s0,' t0=',t0
c      write(16,'(a,i1,a,i1,a,8f22.11)') (( 'BornABCD(',j1,',',j2,',*,*)=  ',((BornABCD(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c      write(16,'(a,i1,a,i1,a,8f22.11)') (( 'Born1BCD(',j1,',',j2,',*,*)=',((Born1BCD(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c      write(16,'(a,i1,a,i1,a,8f22.11)') (( 'BornA1CD(',j1,',',j2,',*,*)=  ',((BornA1CD(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c      write(16,*) '----------------------------------------------------------------------------------------'
c      write(16,'(a,i1,a,i1,a,8f22.11)') (( 'Born2BCD(',j1,',',j2,',*,*)=',((Born2BCD(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c      write(16,'(a,i1,a,i1,a,8f22.11)') (( 'BornA2CD(',j1,',',j2,',*,*)=  ',((BornA2CD(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c      write(16,*) '----------------------------------------------------------------------------------------'
c      write(16,'(a,i1,a,i1,a,8f22.11)') (( 'Born12CD(',j1,',',j2,',*,*)=',((Born12CD(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c      write(16,'(a,i1,a,i1,a,8f22.11)') (( 'Born21CD(',j1,',',j2,',*,*)=  ',((Born21CD(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]

* Fermion propagarotors ini1
      prA1= 1d0/(pA(4)*ph1(4)-pA(3)*ph1(3)-pA(2)*ph1(2)-pA(1)*ph1(1))/2d0
      prB1=-1d0/(pB(4)*ph1(4)-pB(3)*ph1(3)-pB(2)*ph1(2)-pB(1)*ph1(1))/2d0
* Fermion propagarotors ini2
      prA2= 1d0/(pA(4)*ph2(4)-pA(3)*ph2(3)-pA(2)*ph2(2)-pA(1)*ph2(1))/2d0
      prB2=-1d0/(pB(4)*ph2(4)-pB(3)*ph2(3)-pB(2)*ph2(2)-pB(1)*ph2(1))/2d0
* DOUBLE propagators
      prA12= 1d0/( pA(4)*ph1(4)- pA(3)*ph1(3)- pA(2)*ph1(2)- pA(1)*ph1(1)
     $           + pA(4)*ph2(4)- pA(3)*ph2(3)- pA(2)*ph2(2)- pA(1)*ph2(1)
     $           -ph1(4)*ph2(4)+ph1(3)*ph2(3)+ph1(2)*ph2(2)+ph1(1)*ph2(1)
     $           )/2d0
      prB12=-1d0/( pB(4)*ph1(4)- pB(3)*ph1(3)- pB(2)*ph1(2)- pB(1)*ph1(1)
     $            +pB(4)*ph2(4)- pB(3)*ph2(3)- pB(2)*ph2(2)- pB(1)*ph2(1)
     $           -ph1(4)*ph2(4)+ph1(3)*ph2(3)+ph1(2)*ph2(2)+ph1(1)*ph2(1)
     $           )/2d0
      Fprop1 =(1d0/prA1+1d0/prA2)*prA12*CPF12/CPF0-1d0
      Fprop2 =(1d0/prB1+1d0/prB2)*prB12*CPF/CPF0-1d0
      Fprop1B=CPF12/CPF0-1d0
      Fprop2B=CPF/CPF0-1d0
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,*) 'prA1 = ',prA1,'   prB1 = ',prB1,' prA2= ',prA2,'  prB2= ',prB2
c      write(16,*) 'prA12= ',prA12,'  prB12= ',prB12
c      write(16,*) 'Fprop1=   ',Fprop1, '  Fprop2=  ',Fprop2
c      write(16,*) 'Fprop1B= ',Fprop1B,'  Fprop2B= ',Fprop2B
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
      Sig = 3-2*Hel1
      IF( m_KeyArb  .EQ.  0 ) THEN
         CALL GPS_MatrU( gI, ph1,Sig,  ph1,mph, pA,mA,     U11a) ! <1|[1]|a>
         CALL GPS_MatrV( gI, ph1,Sig,  pB,mB,   ph1,mph,   Vb11) ! <b|{1}|1>
* falSe second
         CALL GPS_MatrU( gI, ph1,Sig,  ph2,mph, pA,mA,     U21a) ! <2|[1]|a>
         CALL GPS_MatrV( gI, ph1,Sig,  pB,mB,   ph2,mph,   Vb12) ! <b|{1}|2>
* reverse order
         CALL GPS_MatrU( gI, ph1,Sig,  pA,mA,   ph2,mph,   Ua12) ! <a|[1]|2>
         CALL GPS_MatrV( gI, ph1,Sig,  ph2,mph, pB,mB,     V21b) ! <2|{1}|b>
* for the case when there was ph2 first xk-xk term
         CALL GPS_MatrU( gI, ph1,Sig,  ph1,mph, ph2,mph,   U112) ! <1|[1]|2>
         CALL GPS_MatrV( gI, ph1,Sig,  ph2,mph, ph1,mph,   V211) ! <2|{1}|1>
* for the case when there was ph2 first xk-xk term
         CALL GPS_MatrU( gI, ph1,Sig,  ph2,mph, ph2,mph,   U212) ! <2|[1]|2>
         CALL GPS_MatrV( gI, ph1,Sig,  ph2,mph, ph2,mph,   V212) ! <2|{1}|2>
      ELSE
         CALL GPS_MatrUb(gI, ph1,Sig,  ph1,mph, pA,mA,     U11a)
         CALL GPS_MatrVb(gI, ph1,Sig,  pB,mB,   ph1,mph,   Vb11)
* falSe second
         CALL GPS_MatrUb(gI, ph1,Sig,  ph2,mph, pA,mA,     U21a)
         CALL GPS_MatrVb(gI, ph1,Sig,  pB,mB,   ph2,mph,   Vb12)
* reverse order
         CALL GPS_MatrUb(gI, ph1,Sig,  pA,mA,   ph2,mph,   Ua12)
         CALL GPS_MatrVb(gI, ph1,Sig,  ph2,mph, pB,mB,     V21b)
* for the case when there was ph2 first xk-xk term
         CALL GPS_MatrUb(gI, ph1,Sig,  ph1,mph, ph2,mph,   U112)
         CALL GPS_MatrVb(gI, ph1,Sig,  ph2,mph, ph1,mph,   V211)
* for the case when there was ph2 first xk-xk term
         CALL GPS_MatrUb(gI, ph1,Sig,  ph2,mph, ph2,mph,   U212)
         CALL GPS_MatrVb(gI, ph1,Sig,  ph2,mph, ph2,mph,   V212)
      ENDIF
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,'(a,8f22.11)') 'U11a(*,*)= ',((U11a(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'Vb11(*,*)= ',((Vb11(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'U21a(*,*)= ',((U21a(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'Vb12(*,*)= ',((Vb12(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'Ua12(*,*)= ',((Ua12(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'V21b(*,*)= ',((V21b(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'U112(*,*)= ',((U112(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'V211(*,*)= ',((V211(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'U212(*,*)= ',((U212(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'V212(*,*)= ',((V212(j1,j2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
      Sig = 3-2*Hel2
      IF( m_KeyArb  .EQ.  0 ) THEN
         CALL GPS_MatrU( gI, ph2,Sig,  ph2,mph,  pA,mA,    U22a) ! <2|[2]|a>
         CALL GPS_MatrV( gI, ph2,Sig,  pB,mB,   ph2,mph,   Vb22) ! <b|{2}|2>
* falSe second
         CALL GPS_MatrU( gI, ph2,Sig,  ph1,mph,  pA,mA,    U12a) ! <1|[2]|a>
         CALL GPS_MatrV( gI, ph2,Sig,  pB,mB,   ph1,mph,   Vb21) ! <b|{2}|1>
* reverse order
         CALL GPS_MatrU( gI, ph2,Sig,  pA,mA,   ph1,mph,   Ua21) ! <a|[2]|1>
         CALL GPS_MatrV( gI, ph2,Sig,  ph1,mph, pB,mB,     V12b) ! <1|{2}|b>
* for the case when there was ph1 first xk-xk term
         CALL GPS_MatrU( gI, ph2,Sig,  ph2,mph, ph1,mph,   U221) ! <2|[2]|1>
         CALL GPS_MatrV( gI, ph2,Sig,  ph1,mph, ph2,mph,   V122) ! <1|{2}|2>
c for the case when there was ph1 first xk-xk term
         CALL GPS_MatrU( gI, ph2,Sig,  ph1,mph, ph1,mph,   U121) ! <1|[2]|1>
         CALL GPS_MatrV( gI, ph2,Sig,  ph1,mph, ph1,mph,   V121) ! <1|{2}|1>
      ELSE
         CALL GPS_MatrUb(gI, ph2,Sig,  ph2,mph,  pA,mA,    U22a)
         CALL GPS_MatrVb(gI, ph2,Sig,  pB,mB,   ph2,mph,   Vb22)
c falSe second
         CALL GPS_MatrUb(gI, ph2,Sig,  ph1,mph,  pA,mA,    U12a)
         CALL GPS_MatrVb(gI, ph2,Sig,  pB,mB,   ph1,mph,   Vb21)
c reverse order
         CALL GPS_MatrUb(gI, ph2,Sig,  pA,mA,   ph1,mph,   Ua21)
         CALL GPS_MatrVb(gI, ph2,Sig,  ph1,mph, pB,mB,     V12b)
c for the case when there was ph1 first xk-xk term
         CALL GPS_MatrUb(gI, ph2,Sig,  ph2,mph, ph1,mph,   U221)
         CALL GPS_MatrVb(gI, ph2,Sig,  ph1,mph, ph2,mph,   V122)
c for the case when there was ph1 first xk-xk term
         CALL GPS_MatrUb(gI, ph2,Sig,  ph1,mph, ph1,mph,   U121)
         CALL GPS_MatrVb(gI, ph2,Sig,  ph1,mph, ph1,mph,   V121)
      ENDIF
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,'(a,8f22.11)') 'U22a(*,*)= ',((U22a(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'Vb22(*,*)= ',((Vb22(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'U12a(*,*)= ',((U12a(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'Vb21(*,*)= ',((Vb21(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'Ua21(*,*)= ',((Ua21(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'V12b(*,*)= ',((V12b(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'U221(*,*)= ',((U221(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'V122(*,*)= ',((V122(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'U121(*,*)= ',((U121(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'V121(*,*)= ',((V121(j1,j2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
C fermion lines including v-a coupling contracted with photon polarization ... 
         Cnor=1D0
      IF( m_KeyArb .EQ. 0 ) THEN
         Sig = 3-2*Hel1
         CALL GPS_MakeUW(Cnor,ph1,Sig, pC,  mC,   pA,  mA,     UCAW1) ! v-a inside
         CALL GPS_MakeVW(Cnor,ph1,Sig, pB,  mB,   pD,  mD,     VBDW1) ! v-a inside
         CALL GPS_MakeUW(Cnor,ph1,Sig, pC,  mC,   ph2, mph,    UC2W1) ! v-a inside
         CALL GPS_MakeVW(Cnor,ph1,Sig, ph2,mph,   pD,  mD,     V2DW1) ! v-a inside
         Sig = 3-2*Hel2
         CALL GPS_MakeUW(Cnor,ph2,Sig, pC,  mC,   pA,   mA,    UCAW2) ! v-a inside
         CALL GPS_MakeVW(Cnor,ph2,Sig, pB,  mB,   pD,   mD,    VBDW2) ! v-a inside
         CALL GPS_MakeUW(Cnor,ph2,Sig, pC,  mC,   ph1, mph,    UC1W2) ! v-a inside
         CALL GPS_MakeVW(Cnor,ph2,Sig, ph1,mph,   pD,   mD,    V1DW2) ! v-a inside
      ELSE
         Sig = 3-2*Hel1
         CALL GPS_MakeUWb(Cnor,ph1,Sig, pC,  mC,   pA,  mA,     UCAW1) ! v-a inside
         CALL GPS_MakeVWb(Cnor,ph1,Sig, pB,  mB,   pD,  mD,     VBDW1) ! v-a inside
         CALL GPS_MakeUWb(Cnor,ph1,Sig, pC,  mC,   ph2, mph,    UC2W1) ! v-a inside
         CALL GPS_MakeVWb(Cnor,ph1,Sig, ph2,mph,   pD,  mD,     V2DW1) ! v-a inside
         Sig = 3-2*Hel2
         CALL GPS_MakeUWb(Cnor,ph2,Sig, pC,  mC,   pA,   mA,    UCAW2) ! v-a inside
         CALL GPS_MakeVWb(Cnor,ph2,Sig, pB,  mB,   pD,   mD,    VBDW2) ! v-a inside
         CALL GPS_MakeUWb(Cnor,ph2,Sig, pC,  mC,   ph1, mph,    UC1W2) ! v-a inside
         CALL GPS_MakeVWb(Cnor,ph2,Sig, ph1,mph,   pD,   mD,    V1DW2) ! v-a inside
      ENDIF
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,'(a,8f22.11)') 'UCAW1(*,*)= ',((UCAW1(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'VBDW1(*,*)= ',((VBDW1(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'UC2W1(*,*)= ',((UC2W1(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'V2DW1(*,*)= ',((V2DW1(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'UCAW2(*,*)= ',((UCAW2(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'VBDW2(*,*)= ',((VBDW2(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'UC1W2(*,*)= ',((UC1W2(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'V1DW2(*,*)= ',((V1DW2(j1,j2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
C fermion lines including v-a coupling contracted with incoming momentum ...

         CALL GPS_MakeUX(Cnor,ph1,mph, pC,mC,   pA,mA,    UCAWX1) ! v-a inside
         CALL GPS_MakeVX(Cnor,ph1,mph, pB,mB,   pD,mD,    VBDWX1) ! v-a inside
         CALL GPS_MakeUX(Cnor,ph2,mph, pC,mC,   pA,mA,    UCAWX2) ! v-a inside
         CALL GPS_MakeVX(Cnor,ph2,mph, pB,mB,   pD,mD,    VBDWX2) ! v-a inside

         CALL GPS_MakeUX(Cnor,pB,mB, pC,mC,   pA,mA,    UCAWXB) ! v-a inside
         CALL GPS_MakeVX(Cnor,pA,mA, pB,mB,   pD,mD,    VBDWXA) ! v-a inside
         CALL GPS_MakeUX(Cnor,pD,mD, pC,mC,   pA,mA,    UCAWXD) ! v-a inside
         CALL GPS_MakeVX(Cnor,pC,mC, pB,mB,   pD,mD,    VBDWXC) ! v-a inside
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,'(a,8f22.11)') 'UCAWX1(*,*)= ',((UCAWX1(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'VBDWX1(*,*)= ',((VBDWX1(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'UCAWX2(*,*)= ',((UCAWX2(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'VBDWX2(*,*)= ',((VBDWX2(j1,j2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,'(a,8f22.11)') 'UCAWXB(*,*)= ',((UCAWXB(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'VBDWXA(*,*)= ',((VBDWXA(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'UCAWXD(*,*)= ',((UCAWXD(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'VBDWXC(*,*)= ',((VBDWXC(j1,j2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]

         CALL GPS_MakeUX(Cnor,ph1,mph, pC, mC,   ph2,mph,    UC2WX1) ! v-a inside
         CALL GPS_MakeUX(Cnor,pA, mA,  pC, mC,   ph2,mph,    UC2WXA) ! v-a inside
         CALL GPS_MakeVX(Cnor,ph1,mph, ph2,mph,  pD, mD,     V2DWX1) ! v-a inside
         CALL GPS_MakeUX(Cnor,ph2,mph, pC, mC,   ph1,mph,    UC1WX2) ! v-a inside
         CALL GPS_MakeVX(Cnor,ph2,mph, ph1,mph,  pD, mD,     V1DWX2) ! v-a inside
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,'(a,8f22.11)') 'UC2WX1(*,*)= ',((UC2WX1(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'UC2WXA(*,*)= ',((UC2WXA(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'V2DWX1(*,*)= ',((V2DWX1(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'V1DWX2(*,*)= ',((V1DWX2(j1,j2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]
         CALL GPS_MakeUX(Cnor,pB,mB, pC, mC,   ph2,mph,    UC2WXB) ! v-a inside
         CALL GPS_MakeVX(Cnor,pA,mA, ph2,mph,  pD, mD,     V2DWXA) ! v-a inside
         CALL GPS_MakeVX(Cnor,pB,mB, ph2,mph,  pD, mD,     V2DWXB) ! v-a inside
         CALL GPS_MakeUX(Cnor,pB,mB, pC, mC,   ph1,mph,    UC1WXB) ! v-a inside
         CALL GPS_MakeUX(Cnor,pA,mA, pC, mC,   ph1,mph,    UC1WXA) ! v-a inside
         CALL GPS_MakeVX(Cnor,pA,mA, ph1,mph,  pD, mD,     V1DWXA) ! v-a inside
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,'(a,8f22.11)') 'UC2WXB(*,*)= ',((UC2WXB(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'V2DWXA(*,*)= ',((V2DWXA(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'V2DWXB(*,*)= ',((V2DWXB(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'UC1WXB(*,*)= ',((UC1WXB(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'UC1WXA(*,*)= ',((UC1WXA(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'V1DWXA(*,*)= ',((V1DWXA(j1,j2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]]

         CALL GPS_MakeUX(Cnor,pD,mD, pC, mC,   ph2,mph,    UC2WXD) ! v-a inside
         CALL GPS_MakeVX(Cnor,pC,mC, ph2,mph,  pD, mD,     V2DWXC) ! v-a inside
         CALL GPS_MakeUX(Cnor,pD,mD, pC, mC,   ph1,mph,    UC1WXD) ! v-a inside
         CALL GPS_MakeVX(Cnor,pC,mC, ph1,mph,  pD, mD,     V1DWXC) ! v-a inside
         CALL GPS_MakeVX(Cnor,pB,mB, ph1,mph,  pD, mD,     V1DWXB) ! v-a inside
c[[[[[[[[[[[[[[[[[[[[[[[[[[[
c      write(16,'(a,8f22.11)') 'UC2WXD(*,*)= ',((UC2WXD(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'V2DWXC(*,*)= ',((V2DWXC(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'UC1WXD(*,*)= ',((UC1WXD(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'V1DWXC(*,*)= ',((V1DWXC(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,8f22.11)') 'V1DWXB(*,*)= ',((V1DWXB(j1,j2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]]]]

c[[[[[[[[[[[[[[
c      CALL GPS_BornZero(AmpWork)
c]]]]]]]]]]]]]]
c[[[[[[[[[[[[[[
c      write(16,*) ' sB(2,Hel2)=',sB(2,Hel2)
c      write(16,*) '      prB12=',prB12
c      write(16,*) '   CPF,CPF0=',CPF,CPF0
c      write(16,'(a,8f22.11)') 'Vb12(*,*)= ',((Vb12(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,i1,a,i1,a,8f22.11)') (( 'BornA2CD(',j1,',',j2,',*,*)=  ',((BornA2CD(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c Su1=Su1 +EpsDot21(Hel2)  *Born1BCD(j,j2,j3,j4) *prA1*U11a(j,j1) *CPF1*CPF12/CPF0*I9X  !<b|(2)2|X|1[1]|a>
c      write(16,*) ' EpsDot21(Hel2)=',EpsDot21(Hel2)
c      write(16,*) '      prA1=',prA1
c      write(16,*) '   CPF1=',CPF1,' CPF0=',CPF0,' CPF12=',CPF12
c      write(16,'(a,8f22.11)') 'U11a(*,*)= ',((U11a(j1,j2),j2=1,2),j1=1,2)
c      write(16,'(a,i1,a,i1,a,8f22.11)') (( 'Born1BCD(',j1,',',j2,',*,*)=  ',((Born1BCD(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]
      DO j1=1,2
         DO j2=1,2
            DO j3=1,2
               DO j4=1,2
                  Su1 = DCMPLX(0d0,0d0)

                  DO j=1,2
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    |                  1               2                         //
*        //                    |X                 |               |                         //
*        //      _       -b    |      b+m-1-2     |      a+m-2    |      a                  //
*        //      v  ------<----O--------<---------U--------<------S------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1+Born1BCD(j,j2,j3,j4) *U11a(j,j1)*(prA12-prA1)*sA(2,Hel2)*CPF12/CPF0*IV2 !<b|X|1[1]a(2)|a>
                     Su1=Su1+Born1BCD(j,j2,j3,j4) *U11a(j,j1)*(      prA1)*sA(2,Hel2)*CPF12/CPF0*Y_IR1*I9X !<b|X|1[1]a(2)|a>
                     Su1=Su1+Born2BCD(j,j2,j3,j4) *U21a(j,j1)* prA12      *sA(2,Hel2)*CPF12/CPF0*IV2 !<b|X|2[1]a(2)|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    |                  2               1                         //
*        //                    |X                 |               |                         //
*        //      _       -b    |      b+m-1-2     |      a+m-1    |      a                  //
*        //      v  ------<----O--------<---------U-------<-------S------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1+Born2BCD(j,j2,j3,j4) *U22a(j,j1)*(prA12-prA2)*sA(1,Hel1) *CPF12/CPF0*IV2!<b|X|2[2]a(1)|a>
                     Su1=Su1+Born2BCD(j,j2,j3,j4) *U22a(j,j1)*(      prA2)*sA(1,Hel1) *CPF12/CPF0*Y_IR1*I9Z!<b|X|2[2]a(1)|a>
                     Su1=Su1+Born1BCD(j,j2,j3,j4) *U12a(j,j1)* prA12      *sA(1,Hel1) *CPF12/CPF0*IV2!<b|X|1[2]a(1)|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    2                  1               |                         //
*        //                    |                  |               |X                        //
*        //      _       -b    |     -b+m+2       |    -b+m+1+2   |      a                  //
*        //      v  ------<----S--------<---------V--------<------O------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1+( sB(2,Hel2)*(-prB1+prB12))*Vb11(j2,j)*BornA1CD(j1,j,j3,j4)*CPF/CPF0*IV1!<b|(2)b[1]1|X|a>
                     Su1=Su1+( sB(2,Hel2)*( prB1      ))*Vb11(j2,j)*BornA1CD(j1,j,j3,j4)*CPF/CPF0*Y_IR1*I9Y!<b|(2)b[1]1|X|a>
                     Su1=Su1+( sB(2,Hel2)*(      prB12))*Vb12(j2,j)*BornA2CD(j1,j,j3,j4)*CPF/CPF0*IV1!<b|(2)b[1]2|X|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    1                  2               |                         //
*        //                    |                  |               |X                        //
*        //      _       -b    |     -b+m+1       |    -b+m+1+2   |      a                  //
*        //      v  ------<----S--------<---------V--------<------O------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +(sB(1,Hel1)*(-prB2+prB12))*Vb22(j2,j)*BornA2CD(j1,j,j3,j4)*CPF/CPF0*IV1!<b|(1)b[2]2|X|a>
                     Su1=Su1 +(sB(1,Hel1)*( prB2      ))*Vb22(j2,j)*BornA2CD(j1,j,j3,j4)*CPF/CPF0*Y_IR1*I9T!<b|(1)b[2]2|X|a>
                     Su1=Su1 +(sB(1,Hel1)*(      prB12))*Vb21(j2,j)*BornA1CD(j1,j,j3,j4)*CPF/CPF0*IV1!<b|(1)b[2]1|X|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    2                  |               1                         //
*        //                    |                  |X              |                         //
*        //      _       -b    |     -b+m+2       |     a+m-1     |      a                  //
*        //      v  ------<----S--------<---------O--------<------U------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +sB(2,Hel2)    *Born1BCD(j,j2,j3,j4) *prA1*U11a(j,j1) *CPF1/CPF0*Y_IR1*I9X !<b|(2)2|X|1[1]|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    1                  |               2                         //
*        //                    |                  |X              |                         //
*        //      _       -b    |     -b+m+1       |     a+m-2     |      a                  //
*        //      v  ------<----S--------<---------O--------<------U------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +sB(1,Hel1)    *Born2BCD(j,j2,j3,j4) *prA2*U22a(j,j1) *CPF2/CPF0*Y_IR1*I9Z !<b|(1)1|X|2[2]|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    1                  |               2                         //
*        //                    |                  |X              |                         //
*        //      _       -b    |     -b+m+1       |     a+m-2     |      a                  //
*        //      v  ------<----V--------<---------O--------<------S------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +Vb11(j2,j)*prB1 *BornA1CD(j1,j,j3,j4) *sA(2,Hel2) *CPF2/CPF0   *Y_IR1*I9Y !<b|[1]1|X|a[2]|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    2                  |               1                         //
*        //                    |                  |X              |                         //
*        //      _       -b    |     -b+m+2       |     a+m-1     |      a                  //
*        //      v  ------<----V--------<---------O--------<------S------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                      Su1=Su1 +Vb22(j2,j)*prB2 *BornA2CD(j1,j,j3,j4) *sA(1,Hel1)*CPF1/CPF0      *Y_IR1*I9T !<b|[2]2|X|a(1)|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                                       E--2                                      //
*        //                                       |               1                         //
*        //                                       |X              |                         //
*        //      _       -b                       |     a+m-1     |      a                  //
*        //      v  ------<-------------<---------O--------<------U------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +EpsDot21(Hel2)  *Born1BCD(j,j2,j3,j4) *prA1*U11a(j,j1) *CPF1*CPF12/CPF0*I9X  !<b|(2)2|X|1[1]|a>
     $                       -0.5D0*DCMPLX(ChaIni*m_e_QED)*prA1*U11a(j,j1)*CPF1*CPF12*WVpi0                       *I71
     $                       *( UC1W2(j3,j )*(2*VBDWX2(j2,j4) )
     $                         -VBDW2(j2,j4)*(2*UC1WX2(j3,j ) )                        )
      Su1=Su1     -2*0.5D0*DCMPLX(ChaIni*m_e_QED)*prA1*CPF1*CPF12*WVpi0
     &                          *VBDWX2(j2,j4) *UC1W2( j3,j ) *U11a(j,j1)
      Su1=Su1     +2*0.5D0*DCMPLX(ChaIni*m_e_QED)*prA1*CPF1*CPF12*WVpi0
     &                          *VBDW2( j2,j4) *UC1WX2(j3,j ) *U11a(j,j1)
!!!     $                       -0.5D0*DCMPLX(ChaIni*m_e_QED)*prA1*U11a(j,j1)*CPF1*CPF12*WVpi0                       *I71b   ! (zeroed) included later
!!!     $                       *(
!!!     $                         -VBDW2(j2,j4)*(               -UC1WXA(j3,j) )           )
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                                       E--1                                      //
*        //                                       |               2                         //
*        //                                       |X              |                         //
*        //      _       -b                       |     a+m-2     |      a                  //
*        //      v  ------<-------------<---------O--------<------U------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +EpsDot12(Hel1)  *Born2BCD(j,j2,j3,j4) *prA2*U22a(j,j1) *CPF2*CPF12/CPF0*I9Z!<b|(1)1|X|2[2]|a>
     $                       -0.5D0*DCMPLX(ChaIni*m_e_QED)*prA2*U22a(j,j1)*CPF2*CPF12*WVpi0*I72
     $                       *( UC2W1(j3,j )*(2*VBDWX1(j2,j4) )
     $                         -VBDW1(j2,j4)*(2*UC2WX1(j3,j )              ) )
cccc      Su1=Su1       -2*0.5D0*DCMPLX(ChaIni*m_e_QED)*prA2*CPF2*CPF12*WVpi0*I72
cccc     &                       *VBDWX1(j2,j4) *UC2W1( j3,j) *U22a(j,j1)
cccc      Su1=Su1       +2*0.5D0*DCMPLX(ChaIni*m_e_QED)*prA2*CPF2*CPF12*WVpi0*I72
cccc     &                       *VBDW1(j2,j4)  *UC2WX1(j3,j) *U22a(j,j1)
!!!     $                       -0.5D0*DCMPLX(ChaIni*m_e_QED)*prA2*U22a(j,j1)*CPF2*CPF12*WVpi0*I72b  ! (zeroed) included later
!!!     $                       *(
!!!     $                         -VBDW1(j2,j4)*(               -UC2WXA(j3,j) ) )
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                                       E--2                                      //
*        //                    1                  |                                         //
*        //                    |                  |X                                        //
*        //      _       -b    |     -b+m+1       |                      a                  //
*        //      v  ------<----V--------<---------O--------<-------------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 +Vb11(j2,j)*prB1 *BornA1CD(j1,j,j3,j4) *EpsDot2(Hel2) *CPF *CPF2/CPF0*I9Y  !<b|[1]1|X|a[2]|a>
     $                       -0.5D0*DCMPLX(ChaIni*m_e_QED)*Vb11(j2,j)*prB1*CPF*CPF2*WVpi0*I71
     $                       *( UCAW2(j3,j1)*(2*V1DWX2(j ,j4)              )
     $                         -V1DW2(j ,j4)*(2*UCAWX2(j3,j1)) )
cccc      Su1=Su1  -2*0.5D0*DCMPLX(ChaIni*m_e_QED)*prB1*CPF*CPF2*WVpi0*I71
cccc     &                         *Vb11(j2,j)*V1DWX2(j ,j4)*UCAW2( j3,j1)
cccc      Su1=Su1  +2*0.5D0*DCMPLX(ChaIni*m_e_QED)*prB1*CPF*CPF2*WVpi0*I71
cccc     &                         *Vb11(j2,j)*V1DW2( j ,j4)*UCAWX2(j3,j1)
!!!     $                       -0.5D0*DCMPLX(ChaIni*m_e_QED)*Vb11(j2,j)*prB1*CPF*CPF2*WVpi0*I71b  ! (zeroed) included later
!!!     $                       *( UCAW2(j3,j1)*(               -V1DWXB(j,j4) )
!!!     $                                                         )
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                                       E--1                                      //
*        //                    2                  |                                         //
*        //                    |                  |X                                        //
*        //      _       -b    |     -b+m+2       |                      a                  //
*        //      v  ------<----V--------<---------O--------<-------------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                     Su1=Su1 !+Vb22(j2,j)*prB2 *BornA2CD(j1,j,j3,j4) *EpsDot1(Hel1)*CPF *CPF1/CPF0*I9T  !<b|[2]2|X|a(1)|a>
     $                       -0.5D0*DCMPLX(ChaIni*m_e_QED)*Vb22(j2,j)*prB2*CPF*CPF1*WVpi0*I72
     $                       *( UCAW1(j3,j1)*(2*V2DWX1(j ,j4)               )
     $                         -V2DW1(j ,j4)*(2*UCAWX1(j3,j1) ) )
cccc      Su1=Su1  -2*0.5D0*DCMPLX(ChaIni*m_e_QED) *prB2*CPF*CPF1*WVpi0  *I72
cccc     &                           *Vb22(j2,j) *V2DWX1(j,j4) *UCAW1( j3,j1)
cccc      Su1=Su1  +2*0.5D0*DCMPLX(ChaIni*m_e_QED) *prB2*CPF*CPF1*WVpi0  *I72
cccc     &                           *Vb22(j2,j) *V2DW1( j,j4) *UCAWX1(j3,j1)
!!!     $                       -0.5D0*DCMPLX(ChaIni*m_e_QED)*Vb22(j2,j)*prB2*CPF*CPF1*WVpi0*I72b  ! (zeroed) included later
!!!     $                       *( UCAW1(j3,j1)*(               -V2DWXB(j ,j4) )
!!!     $                                                          )
                  ENDDO

                  Su2 = DCMPLX(0d0,0d0)
                  DO j=1,2
                     DO l=1,2
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    2                  |               1                         //
*        //                    |                  |X              |                         //
*        //      _       -b    |     -b+m+2       |     a+m-1     |      a                  //
*        //      v  ------<----*--------<---------O--------<------O------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2  +Vb22(j2,l)*prB2  *Born12CD(j,l,j3,j4 )  *U11a(j,j1)*prA1*CPF1/CPF0*I8 ! <b|[2]2|X|1[1]|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    1                  |               2                         //
*        //                    |                  |X              |                         //
*        //      _       -b    |     -b+m+1       |     a+m-2     |      a                  //
*        //      v  ------<----*--------<---------O--------<------O------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2  +Vb11(j2,l)*prB1  *Born21CD(j,l,j3,j4 )  *U22a(j,j1)*prA2*CPF2/CPF0*I8 ! <b|[1]1|X|2[2]|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    |                  2               1                         //
*        //                    |X                 |               |                         //
*        //      _       -b    |      b+m-1-2     |      a+m-1    |      a                  //
*        //      v  ------<----O--------<---------O--------<------*------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Born1BCD(j,j2,j3,j4) *U121(j,l)*prA12 *U11a(l,j1)*prA1*CPF12/CPF0*IV2  ! <b|X|1[2]1(1)|a>
                        Su2=Su2 +Born2BCD(j,j2,j3,j4) *U221(j,l)*prA12 *U11a(l,j1)*prA1*CPF12/CPF0*IA  ! <b|X|2[2]1(1)|a>
                        Su2=Su2 -BornABCD(j,j2,j3,j4) *Ua21(j,l)*prA12 *U11a(l,j1)*prA1*CPF12/CPF0*IV2  ! <b|X|a[2]1(1)|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    |                  1               2                         //
*        //                    |X                 |               |                         //
*        //      _       -b    |      b+m-1-2     |      a+m-2    |      a                  //
*        //      v  ------<----O--------<---------O--------<------*------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Born2BCD(j,j2,j3,j4) *U212(j,l)*prA12  *U22a(l,j1)*prA2*CPF12/CPF0*IV2 ! <b|X|2[1]2(2)|a>
                        Su2=Su2 +Born1BCD(j,j2,j3,j4) *U112(j,l)*prA12  *U22a(l,j1)*prA2*CPF12/CPF0*IA ! <b|X|1[1]2(2)|a>
                        Su2=Su2 -BornABCD(j,j2,j3,j4) *Ua12(j,l)*prA12  *U22a(l,j1)*prA2*CPF12/CPF0*IV2 ! <b|X|a[1]2(2)|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    1                  2               |                         //
*        //                    |                  |               |X                        //
*        //      _       -b    |     -b+m+1       |    -b+m+1+2   |      a                  //
*        //      v  ------<----*--------<---------O--------<------O------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Vb11(j2,l)*prB1 *V121(l,j)*prB12 *BornA1CD(j1,j,j3,j4)*CPF/CPF0*IV1 ! <b|(1)1[2]1|X|a>
                        Su2=Su2 +Vb11(j2,l)*prB1 *V122(l,j)*prB12 *BornA2CD(j1,j,j3,j4)*CPF/CPF0*IA ! <b|(1)1[2]2|X|a>
                        Su2=Su2 -Vb11(j2,l)*prB1 *V12b(l,j)*prB12 *BornABCD(j1,j,j3,j4)*CPF/CPF0*IV1 ! <b|(1)1[2]b|X|a>
*        /////////////////////////////////////////////////////////////////////////////////////
*        //                    2                  1               |                         //
*        //                    |                  |               |X                        //
*        //      _       -b    |     -b+m+2       |    -b+m+1+2   |      a                  //
*        //      v  ------<----*--------<---------O--------<------O------<----- u           //
*        /////////////////////////////////////////////////////////////////////////////////////
                        Su2=Su2 +Vb22(j2,l)*prB2 *V212(l,j)*prB12  *BornA2CD(j1,j,j3,j4)*CPF/CPF0*IV1 ! <b|(2)2[1]2|X|a>
                        Su2=Su2 +Vb22(j2,l)*prB2 *V211(l,j)*prB12  *BornA1CD(j1,j,j3,j4)*CPF/CPF0*IA ! <b|(2)2[1]1|X|a>
                        Su2=Su2 -Vb22(j2,l)*prB2 *V21b(l,j)*prB12  *BornABCD(j1,j,j3,j4)*CPF/CPF0*IV1 ! <b|(2)2[2]b|X|a>
                     ENDDO
                  ENDDO
c===========================
      AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4)
     $                 +CNorm*( Su1+Su2 )
c===========================
      sProd = ( sA(1,Hel1)+sB(1,Hel1))*( sA(2,Hel2)+sB(2,Hel2))
      AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4) ! non-infrared part: Fprop1- Fprop1B) is prop to ph1*ph2
     $                 +CNorm*BornABCD(j1,j2,j3,j4)*( sA(1,Hel1)*sA(2,Hel2)*(Fprop1- Fprop1B)      *IV2 !
     $                                               +sB(1,Hel1)*sB(2,Hel2)*(Fprop2- Fprop2B)      *IV1
     $                                                                                      ) !
c===========================
      AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4) ! infrared-type part
     $                 +CNorm*BornABCD(j1,j2,j3,j4)*( sA(1,Hel1)*sA(2,Hel2)*Fprop1B                      *IVI !
     $                                               +sA(1,Hel1)*sB(2,Hel2)*(CPF1/CPF0-1d0)              *IVI
     $                                               +sB(1,Hel1)*sA(2,Hel2)*(CPF2/CPF0-1d0)              *IVI
     $                                               +sB(1,Hel1)*sB(2,Hel2)*Fprop2B                      *IVI
     $                                               +sA(1,Hel1)*CPF1*CPF12/CPF0*EpsDot21(Hel2)          *IVI! terms due to diag WWgamma
     $                                               +sA(2,Hel2)*CPF2*CPF12/CPF0*EpsDot12(Hel1)          *IVI
     $                                               +sB(1,Hel1)*CPF*CPF2/CPF0*EpsDot2(Hel2)             *IVI
     $                                               +sB(2,Hel2)*CPF*CPF1/CPF0*EpsDot1(Hel1)             *IVI
     $                                               +CPF*CPF1*CPF12/CPF0*EpsDot1(Hel1)*EpsDot21(Hel2)   *IVI
     $                                               +CPF*CPF2*CPF12/CPF0*EpsDot2(Hel2)*EpsDot12(Hel1)   *IVI
     $                                                                                      ) !
     $                 +CNorm*BornABCD(j1,j2,j3,j4)* sProd                   *Y_IR                        !
c===========================
! terms due to ph1 ph2 attached to W
      AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4) + CNorm*( 0.5D0) *DCMPLX(-ChaIni*m_e_QED)*CPF*CPF1*CPF12*WVpi0*(
ccc     $  ! absorbed to infraed-type      BornABCD(j1,j2,j3,j4)/CPF0/WVpi0*(2*eps1pC-2*eps1pA)*(2*eps2pB-2*eps2pD)* *IV ! 1
     $ +EpsDot1 (hel1)*VBDW2(j2,j4)*(               -2*UCAWX2(j3,j1)) *I9s1 ! 2
     $ +EpsDot1 (hel1)*VBDW2(j2,j4)*(-UCAWX1(j3,j1)                 ) *I9   ! 2
     $ +EpsDot1 (hel1)*2*VBDWX2(j2,j4) * UCAW2 (j3,j1)                *I9s1 ! 3
     $ -2*UCAWX1(j3,j1)*EpsDot21(hel2)* VBDW1 (j2,j4)                *I9s2 ! 4
     $ -UCAWX1(j3,j1)* VBDW2 (j2,j4) *(-EpsDot12(hel1)-2*eps1p2*DCMPLX(-ChaIni*m_e_QED))         *I9   ! 5
     $ + eps1D2*(-2*UCAWX1(j3,j1))*        2* VBDWX2(j2,j4)                *I9   *DCMPLX(-ChaIni*m_e_QED)! 6
     $ + UCAW1(j3,j1)*(2*VBDWX1(j2,j4)              )*EpsDot21(hel2)  *I9s2 ! 7
     $ + UCAW1(j3,j1)*(                VBDWX2(j2,j4))*EpsDot21(hel2)  *I9   ! 7
ccc!!!!     $ + UCAW1(j3,j1)*(t1-2*p2p1 )*VBDW2(j2,j4)                           *I9   *DCMPLX(-ChaIni*m_e_QED)! 8 (out)
     $ + UCAW1(j3,j1)*(1D0/CPF2-2*p2p1 )*VBDW2(j2,j4)                      *I9   *DCMPLX(-ChaIni*m_e_QED)! 8    'Higgs' contrib added
     $ + UCAW1(j3,j1)*VBDWX2(j2,j4)*(-EpsDot2(hel2)+2*eps2p1*DCMPLX(-ChaIni*m_e_QED))            *I9   ! 9
     $                                                                            )
c===========================
! terms due to ph2 ph1 attached to W
      AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4) + CNorm*( 0.5D0) *DCMPLX(-ChaIni*m_e_QED)*CPF*CPF2*CPF12*WVpi0*(
ccc     $   ! absorbed to infraed-type     BornABCD(j1,j2,j3,j4)/CPF0/WVpi0*(2*eps2pC-2*eps2pA)*(2*eps1pB-2*eps1pD)*IV ! 1
     $ +EpsDot2 (hel2)*VBDW1(j2,j4)*(               -2*UCAWX1(j3,j1)) *I9s2 ! 2
     $ +EpsDot2 (hel2)*VBDW1(j2,j4)*(-UCAWX2(j3,j1)                 ) *I9   ! 2
     $ +EpsDot2 (hel2)*2*VBDWX1(j2,j4) * UCAW1 (j3,j1)                *I9s2 ! 3
     $ - 2*UCAWX2(j3,j1)*EpsDot12(hel1)* VBDW2 (j2,j4)                *I9s1 ! 4
     $ -UCAWX2(j3,j1)* VBDW1 (j2,j4)*(-EpsDot21(hel2)-2*eps2p1*DCMPLX(-ChaIni*m_e_QED))           *I9  ! 5
     $ + eps1D2*(-2*UCAWX2(j3,j1))*        2* VBDWX1(j2,j4)                *I9  *DCMPLX(-ChaIni*m_e_QED) ! 6
     $ + UCAW2(j3,j1)*(2*VBDWX2(j2,j4)              )*EpsDot12(hel1)  *I9s1 ! 7
     $ + UCAW2(j3,j1)*(                VBDWX1(j2,j4))*EpsDot12(hel1)  *I9   ! 7
!!!!     $ + UCAW2(j3,j1)*(t2-2*p2p1)*VBDW1(j2,j4)                          *I9  *DCMPLX(-ChaIni*m_e_QED) ! 8 (out)
     $ + UCAW2(j3,j1)*(1D0/CPF1-2*p2p1)*VBDW1(j2,j4)                       *I9  *DCMPLX(-ChaIni*m_e_QED) ! 8   'Higgs' contrib added
     $ + UCAW2(j3,j1)*VBDWX1(j2,j4)*(-EpsDot1(hel1)+2*eps1p2*DCMPLX(-ChaIni*m_e_QED))            *I9   ! 9
     $                                                                            )
c===========================
      AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4)   ! terms due to 4 boson coupling
     $                 +CNorm*(1D0) *(m_e_QED**2)*CPF*CPF12*WVpi0*(        !!!+CNorm*(-1D0) *DCMPLX(0,-m_e_QED**2)*CPF*CPF12*WVpi0
     $                  -BornABCD(j1,j2,j3,j4)/CPF0/WVpi0*2*eps1D2          *IVI
     $                      +0.5D0 *UCAW1(j3,j1)*VBDW2(j2,j4)               *I10
     $                      +0.5D0 *UCAW2(j3,j1)*VBDW1(j2,j4)               *I10
     $                                                                )
c===========================
      AmpWork(j1,j2,j3,j4) =AmpWork(j1,j2,j3,j4)+ CNorm*(-0.5D0) *DCMPLX(ChaIni*m_e_QED)*WVpi0 *(  ! terms abcd due rem. single WWgam coupl (redistr)
     $  +( UCAW1(j3,j1)*( 2*VBDWX1(j2,j4) ) -VBDW1(j2,j4)*( 2*UCAWX1(j3,j1)             ) )*sA(2,Hel2)      *CPF2*CPF12 *I9s2      ! (a) plus below
     $  +( UCAW1(j3,j1)*( 2*VBDWX1(j2,j4)               )-VBDW1(j2,j4)*( 2*UCAWX1(j3,j1)) )*sB(2,Hel2)      *CPF *CPF1  *I9s2      ! (b) plus below
     $  +( UCAW2(j3,j1)*( 2*VBDWX2(j2,j4) ) -VBDW2(j2,j4)*( 2*UCAWX2(j3,j1)             ) )*sA(1,Hel1)      *CPF1*CPF12 *I9s1      ! (c) plus below
     $  +( UCAW2(j3,j1)*( 2*VBDWX2(j2,j4)               )-VBDW2(j2,j4)*( 2*UCAWX2(j3,j1)) )*sB(1,Hel1)      *CPF*CPF2   *I9s1      ! (d) plus below
     $                                      -VBDW1(j2,j4)*(                  UCAW2(j3,j1) )*DCMPLX( m_e_QED)*CPF2*CPF12 *I9B       ! (a) plus  I72b
     $  +( UCAW1(j3,j1)*(                 VBDW2(j2,j4)  )                                 )*DCMPLX(-m_e_QED)*CPF *CPF1  *I9B       ! (b) plus  I72b
     $                                      -VBDW2(j2,j4)*(                  UCAW1(j3,j1) )*DCMPLX( m_e_QED)*CPF1*CPF12 *I9B       ! (c) plus  I71b
     $  +( UCAW2(j3,j1)*(                 VBDW1(j2,j4) )                                  )*DCMPLX(-m_e_QED)*CPF*CPF2   *I9B       ! (d) plus  I71b
     $                                                                         )
               ENDDO
            ENDDO
         ENDDO
      ENDDO
c[[[[[[[[[[[[[[[[[[[[[[[*debug*
c      write(16,*) '////////////////////////////////////////////////////////////////////////////////////////'
c      write(16,*) '///////////////////////////////////////AmpWork//////////////////////////////////////////'
c      write(16,'(a,i1,a,i1,a,8g22.11)') (( 'AmpWork(',j1,',',j2,',*,*)=  ',((AmpWork(j1,j2,j3,j4),j4=1,2),j3=1,2),j2=1,2),j1=1,2)
c]]]]]]]]]]]]]]]]]]]]]]]
      END                       ! GPS_HiiPlusW


























