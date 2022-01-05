*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//                     Pseudo-CLASS  KarFin                                        //
*//                                                                                 //
*//   Purpose:                                                                      //
*//   Top level  Monte-Carlo event generator for FSR radiadion.                     //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////


      SUBROUTINE KarFin_Initialize(xpar_input)
*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//    Initialization of input and internal variables                               //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
      INCLUDE 'BXformat.h'
      DOUBLE PRECISION   xpar_input(*)
      DOUBLE PRECISION   vvmin,cmsene,delfac
*
      m_Nevgen = 0
      m_MarTot = 0
      m_out    = xpar_input(4)
      m_idyfs  = xpar_input(8)
      m_nmax   = xpar_input(19)/2
      
      m_KeyFSR = xpar_input(21)
      m_KeyPia = xpar_input(22)
      m_MltFSR = xpar_input(24)
      m_KeyQSR = xpar_input(29)
      m_alfinv = xpar_input(30)

      m_IsFSR  = m_KeyFSR ! initial value only

      vvmin     = xpar_input(16)
      cmsene    = xpar_input( 1)
      m_emin    = cmsene/2d0*vvmin
      m_MasPhot =1d-60 ! dummy photon mass
     
      delfac   = xpar_input(18)
      m_delta  = vvmin*delfac

      CALL KK2f_GetXenph(m_Xenph)

      m_WtMass = 1d0

      CALL GLK_Mbook(m_idyfs+60,'KarLud: phot. mult. raw $', 1, 1d0*m_nmax)
      CALL GLK_Mbook(m_idyfs+69,'KarLud: totat weight    $', 1, 2d0)

      IF(m_KeyFSR .EQ. 0) RETURN

      CALL GLK_Mbook(m_idyfs+64,'YFSfin: photon raw multiplicity$', 10, 0.2d0*m_nmax)
      CALL GLK_Mbook(m_idyfs+61,'YFSfin: wt1, kinematics $', 1, 2d0)
      CALL GLK_Mbook(m_idyfs+62,'YFSfin: wt2, jacobian   $', 1, 2d0)
      CALL GLK_Mbook(m_idyfs+63,'YFSfin: wt3, mass       $', 1, 2d0)
      CALL GLK_Mbook(m_idyfs+65,'YFSfin: piatek, wtctrl  $', 1, 2d0)
      CALL GLK_Mbook(m_idyfs+66,'YFSfin: piatek, wtrem   $', 1, 2d0)

      WRITE(m_out,bxope)
      WRITE(m_out,bxtxt) 'KarFin Initialize  START'
      WRITE(m_out,bxl1i) m_KeyFSR,  'FSR radiation on/off','KeyFSR','a1'
      WRITE(m_out,bxl1i) m_KeyQSR,  'radiation from quark','KeyQSR','a2'
      WRITE(m_out,bxl1i) m_KeyPia,  'removal    switch   ','KeyPia','a3'
      WRITE(m_out,bxl1g) delfac,    'infrared cut FACTOR ','delfac','a4'
      WRITE(m_out,bxl1g) m_delta,   'infrared cut itself ','delta ','a5'
      WRITE(m_out,bxl1g) m_emin ,   'EminCMS for removal ','[GeV] ','a6'
      WRITE(m_out,bxl1i) m_nmax,    'Max. photon mult.   ','nmax  ','a7'
      WRITE(m_out,bxtxt) 'KarFin Initialize  END  '
      WRITE(m_out,bxclo)

      END                       !!!KarFin_Initialize!!!


      SUBROUTINE KarFin_Finalize(mode)
*     ********************************
      IMPLICIT NONE
*
      INCLUDE 'BXformat.h'
      INCLUDE 'KarFin.h'
*-----------------------------------------------------------------------------
      INTEGER mode
      DOUBLE PRECISION   avmlt,ermlt,upmlt
      DOUBLE PRECISION    WTsup, AvUnd, AvOve, ROverf, AveWt, ERela
      INTEGER  Nevtot,Nevacc,Nevneg,Nevove,Nevzer
*
      DOUBLE PRECISION   awt69,dwt69
      DOUBLE PRECISION   evove,evneg,evacc,evzer,evtot,dumm2
*----------
      DOUBLE PRECISION   awt61,dwt61
      DOUBLE PRECISION   awt62,dwt62
      DOUBLE PRECISION   awt63,dwt63
      DOUBLE PRECISION   awt65,dwt65
      DOUBLE PRECISION   awt66,dwt66
      INTEGER ntot66,nzer66
*
      INTEGER Ntot, Nacc, Nneg, Nove, Nzer
*---------------------------------------------------------------------------
*
      IF(m_KeyFSR .EQ. 0) RETURN
* no printout for mode=1
      IF(mode .EQ. 1) RETURN
* no printout if there was no generation
      IF(m_Nevgen .EQ. 0) RETURN
*============================================================================
      CALL  GLK_Mprint(m_idyfs+64)
*
      WRITE(m_out,bxope)
      WRITE(m_out,bxtxt) 'KarFin Finalize    START'
      WRITE(m_out,bxl1i) m_NevGen,   ' generated events  ','nevgen','a2'

      CALL  GLK_MgetAve(m_idyfs+61,awt61,dwt61,WtSup)
      CALL  GLK_MgetAve(m_idyfs+62,awt62,dwt62,WtSup)
      CALL  GLK_MgetAve(m_idyfs+63,awt63,dwt63,WtSup)

      WRITE(m_out,bxl2f) awt61,dwt61,' kinematics, smin  ','wt1   ','a5'
      WRITE(m_out,bxl2f) awt62,dwt62,' jacobian          ','wt2   ','a6'
      WRITE(m_out,bxl2f) awt63,dwt63,' photon ang. dist. ','wt3   ','a7'
*
* Details on MASS WEIGHT REARRANGENMENT, and PHOTON REMOVAL
*
      CALL  GLK_MgetAve(m_idyfs+64,avmlt,ermlt,upmlt)
      CALL  GLK_MgetAve(m_idyfs+65,awt65,dwt65,WtSup)

      CALL GLK_MgetAll(m_idyfs+66,
     $     awt66,dwt66, WtSup, AvUnd, AvOve,
     $     ntot66,Nacc,Nneg,Nove,nzer66)
*
      WRITE(m_out,bxtxt) '    ON MASS WEIGHTS     '
      WRITE(m_out,bxl2f) awt66,dwt66,'removal wgt wtrem  ','     ',' b1'
      WRITE(m_out,bxl1i) ntot66,     'no. of raw events  ','     ',' b2'
      WRITE(m_out,bxl1i) nzer66,     'wt6=0      events  ','     ',' b3'
      WRITE(m_out,bxl2f) awt65,dwt65,'control wgt wctrl  ','     ',' b4'
      WRITE(m_out,bxl1i) m_MarTot,   'marked photons     ','MarTot','a5'
      WRITE(m_out,bxl1g) m_emin,     'emin               ','     ',' b6'
      WRITE(m_out,bxl1g) m_delta,    'delta              ','     ',' b7'
      WRITE(m_out,bxl1f) avmlt,      'raw ph. multipl.   ','     ',' b8'
      WRITE(m_out,bxl1f) upmlt,      'Highest phot. mult.','     ',' b9'
      WRITE(m_out,bxtxt) 'YFSfin Finalize    END  '
      WRITE(m_out,bxclo)
* test histograms
*[[      CALL GLK_Print(m_idyfs+20)
*[[      CALL gopera(m_idyfs+31,'/',m_idyfs+32,m_idyfs+33,1d0,1d0)
*[[      CALL GLK_Print(m_idyfs+33)
*=========================================================================
*-----------------------------------------------------------------------
*.........................output window a...............................
      CALL GLK_MgetAve(m_idyfs+60, avmlt, ermlt, upmlt)
      CALL GLK_MgetAve(m_idyfs+69,awt69,dwt69,WtSup)


* general information on weights
      WRITE(m_out,bxope)
      WRITE(m_out,bxtxt) '    KarFin Finalize     '
      WRITE(m_out,bxl1i) m_NevGen,   ' generated events  ','nevgen','a2'
      WRITE(m_out,bxl2f) awt69,dwt69,' general weight    ','wt    ','a1'
      WRITE(m_out,bxl1f) avmlt,      ' aver. ph. multi.  ','avmlt ','a3'
      WRITE(m_out,bxclo)
      END                       !!!KarFin_Finalize!!!


      SUBROUTINE KarFin_Make(PX,amfi1,amfi2,CharSq,wt)
*/////////////////////////////////////////////////////////////////////////////
* INPUT:
*     amfi1,2 = final masses, may change event per event (resonances)
*     CharSq  = charge squared
*     PX     = 4-momentum of the entire FSR system (fermions+photons)
*     KFfin     = final state fermion KF code
* OUTPUT:
*     m_qf1,2   = final state charged particle momenta
*     m_nphot   = FSR photon multiplicity
*     m_phot    = list of FSR photons
*     phsu    = total FSR photon 4-momentum
*     wt      = Crude mc weight
*/////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
      INCLUDE 'BXformat.h'
      DOUBLE PRECISION   CharSq,amfi1,amfi2
      DOUBLE PRECISION   PX(4)
      DOUBLE PRECISION   wt,WtFin,wtt
      INTEGER KFfin
      INTEGER k,j
*----------------------------------------------------------------------------
*
* Initialize momenta to zero
      DO k=1,4
         m_q1(k)   =0d0
         m_q2(k)   =0d0
         m_phsu(k) =0d0
      ENDDO
      m_nphot=0
      DO j=1,m_npmx
         DO k=1,4
            m_phot(j,k) =0d0
         ENDDO
      ENDDO
*
*     Generate photons and fermions in the rest frame of ferms Q=qf1+qf2
*     and transform to LAB system, 
*     PX=qf1+qf2+phsu is total momentum of the entire FSR system,
*     as seen from the LAB system.
*
      m_IsFSR = m_KeyFSR
      CALL KarLud_GetKFfin(KFfin)
* check for quarks KeyQSR flag
      IF( ABS(KFfin) .LT. 10) m_IsFSR = m_IsFSR*m_KeyQSR
* check for neutrinos
      IF( ABS(KFfin) .EQ. 12) m_IsFSR = 0 ! nu_el
      IF( ABS(KFfin) .EQ. 14) m_IsFSR = 0 ! nu_mu
      IF( ABS(KFfin) .EQ. 16) m_IsFSR = 0 ! nu_tau
*-----------------------------------------------------------------------
      IF(m_IsFSR .EQ. 1) THEN
         m_NevGen = m_NevGen+1  ! Only bremsstrahlung counts!!!
         CALL KarFin_YFSfin( PX,  amfi1, amfi2,  CharSq,  WtFin)
      ELSE
*-----------------------------------------------------------------------
*     No final state bremss, fermion momenta defined in Z frame
         CALL KinLib_phspc2( PX,amfi1,amfi2,m_q1,m_q2,wtt)
         WtFin = 1d0
         m_nphot  = 0
      ENDIF
*-----------------------------------------------------------------------
*     main weight
      wt=WtFin

      CALL GLK_Mfill(m_idyfs+69, wt,  1d0)
      CALL GLK_Mfill(m_idyfs+60, 1d0*m_nphot,  1d0)
* =============================================
      END                       !!!KarFin_Make!!!


      SUBROUTINE KarFin_YFSfin( PX,Mas1,Mas2,CharSq,WtFin)
*/////////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                         //
*//  Simulates final state bremsstrahlung out of pair of charged patricles.                 //
*//  (algorithm of 11-th March 1989, by S. Jadach)                                          //
*//  INPUT  :                                                                               //
*//      MltFSR  = normaly set to zero, otherwise enforces photon                           //
*//                multiplicity to be exactly equal MltFSR (special tests)                  //
*//      KeyPia  = photon removal switch                                                    //
*//      PX     = 4-momentum of FSR system as seen from LAB system                          //
*//      Mas1,2 = masses of the final charged pair                                          //
*//      delta   = lower energy bound (dimensionless)                                       //
*//      emin    = lower energy bound (GeV) in LAB system                                   //
*//      CharSq  = final state charge squared                                               //
*//      alfinv  = 1/alpha_QED                                                              //
*//  OUTPUT :                                                                               //
*//      qf1,2    = final fermion four momentum (GeV)                                       //
*//      nphot   = photon multiplicity                                                      //
*//      phot    = photon four momenta (GeV) in cms                                         //
*//      phsu    = sum of photon momenta                                                    //
*//      ygr,zet = Sudakov variables                                                        //
*//      WtFin   = MC weight at this stage                                                  //
*//      martot  = control variable, no of marked photons                                   //
*//  OTHER:                                                                                 //
*//      Mark    = marks on photons CLOSE to lower energy bound                             //
*//      qf1c,2c = ficticious final fermion four momenta for crude S-factor                 //
*//      wt1    = the weight - phase space limits for very hard phot.                       //
*//      wt2    = the weight - translation jacobian.                                        //
*//      wt3    = the weight - single final mass-weight                                     //
*//      wtm    = the list/matrix of mass-weights.                                          //
*//                                                                                         //
*/////////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
*
      DOUBLE PRECISION  pi
      PARAMETER(pi=3.1415926535897932d0)
*
      DOUBLE PRECISION   PX(4)

      DOUBLE PRECISION   WtMlist(100)
      DOUBLE PRECISION   rr(100),xk(100),cgx(100),sgx(100)
      DOUBLE PRECISION   dis0(100)
      REAL               rvec(100)
      INTEGER Mark(100)
*
      INTEGER i,j,k
      DOUBLE PRECISION  phi,cg,sg,eta1,eta2,uu,ul,del1,del2
      DOUBLE PRECISION  yy,xmk2,ener,smini
      DOUBLE PRECISION  wt1,wt2,wt3,WtFin
      DOUBLE PRECISION  betc,etc1,etc2,amd1
      DOUBLE PRECISION  average,Mas1,Mas2,svar,amfin,sprim,qmsene,amcru
      DOUBLE PRECISION  alf1,dis1,amc2,CharSq,xsum
      DOUBLE PRECISION  dl1,dl2,betn,xfact,gamf2,dist1,amd2
*-----------------------------------------------------------------------
      m_NevGen = m_NevGen + 1
*
      alf1 = 1d0/pi/m_alfinv
      wt1   =1d0
      wt2   =1d0
      wt3   =1d0
*-----------------------------------------------------------------------
      svar = PX(4)**2-PX(3)**2-PX(2)**2-PX(1)**2
      amfin = MIN(Mas1,Mas2)
      amc2  = 4.d0*amfin**2/svar  ! overvalued svar for crude
      betc  = sqrt(1d0-amc2)      ! overvalued svar for crude
      DO i=1,100
         rr(i)=0.d0
         Mark(i)=0
      ENDDO
      DO k=1,4
         m_phsu(k)=0.d0
      ENDDO
      DO i=1,100
         m_yfin(i)=0
         m_zfin(i)=0
         DO k=1,4
            m_phot(i,k)=0.d0
         ENDDO
      ENDDO
*-----------------------------------------------------------------------
*  generate photon multiplicity, average = average multiplicity (crude)
*-----------------------------------------------------------------------
      gamf2 = CharSq*alf1 *(1+betc**2)/betc *dlog((1d0+betc)**2/amc2)
      average = gamf2*dlog(1/m_delta)
      average = average *m_Xenph
    5 CONTINUE
      CALL KarFin_PoissGen(average,m_nmax,m_nphot,rr)
c[[[[[[[[[[[[[[[[[[[[
c      write(m_out,*) '//// KarFin_YFSfin:amc2=',amc2,' alf1=',alf1,'  m_delta=',m_delta
c      write(m_out,*) '//// KarFin_YFSfin:betc=',betc,' gamf2=',gamf2,'  CharSq=',CharSq
c      write(m_out,*) '//// KarFin_YFSfin: m_nphot=', m_nphot,'  average=',average, ' m_Xenph=',m_Xenph
c      write(m_out,*) '//// KarFin_YFSfin: ',('rr(',i,')=',rr(i),i=1,m_nphot)
c]]]]]]]]]]]]]]]]]]]]

** This is for tests of program at fixed multiplicity (advanc. users)
      IF((m_MltFSR .NE. 0) .AND. (m_nphot .NE. m_MltFSR)) GOTO 5
**
      IF(m_nphot .EQ. 0) THEN
         sprim=svar
      ELSE
*-----------------------------------------------------------------------
*     begin with photon energy
         xsum=0.d0
         DO  i=1,m_nphot
            xk(i)=m_delta**rr(i)
            IF(xk(i) .LT. sqrt(10.d0)*m_delta) Mark(i)=1
            xsum=xsum+xk(i)
         ENDDO
         IF(xsum .GE. 1.d0) GOTO 900
         xfact=1d0/(1.d0-xsum)
         DO i=1,m_nphot
            xk(i)=xk(i)*xfact
         ENDDO
         CALL PseuMar_MakeVec(rvec,m_nphot)
         DO i=1,m_nphot
*-----------------------------------------------------------------------
*     simplified photon angular distribution,
*     s'->s and m**2/(kp)**2 dropped
*     cg=cos(theta) and sg=sin(theta) memorized to avoid rounding err.
            CALL KarFin_AngBre(amc2,dl1,dl2,cg,sg,dis0(i),dis1)
*-----------------------------------------------------------------------
*     define photon momenta (in units of sqrt(s')/2 )
            phi=2.d0*pi*rvec(i)
            m_phot(i,1)=xk(i)*sg*cos(phi)
            m_phot(i,2)=xk(i)*sg*sin(phi)
            m_phot(i,3)=xk(i)*cg
            m_phot(i,4)=xk(i)
            DO k=1,4
               m_phsu(k)=m_phsu(k)+m_phot(i,k)
            ENDDO
            cgx(i)=cg
            sgx(i)=sg
         ENDDO
*-----------------------------------------------------------------------
*     determine rescaling factor and s', wt2 is dilatation jacobian
         xmk2 = m_phsu(4)**2-m_phsu(3)**2-m_phsu(2)**2-m_phsu(1)**2
         yy   = 1.d0/(1.d0 +m_phsu(4) +xmk2/4.d0 )
         wt2  = yy*(1.d0+m_phsu(4))
         sprim= svar*yy
*-----------------------------------------------------------------------
*     reject events with too hard photons
         smini= (Mas1+Mas2)**2
         IF(sprim .LT. smini) GOTO 900
*-----------------------------------------------------------------------
*     Recsale properly all photon momenta
*-----------------------------------------------------------------------
         ener = sqrt(sprim)/2.d0
         DO  k=1,4
            m_phsu(k)= m_phsu(k)*ener
            DO  i=1,m_nphot
               m_phot(i,k)=m_phot(i,k)*ener
            ENDDO
         ENDDO
      ENDIF ! m_nphot
*-----------------------------------------------------------------------
*     final fermion momenta
*-----------------------------------------------------------------------
      amcru  = amfin*SQRT(sprim/svar)
      qmsene = SQRT(sprim)
      ener   = qmsene/2d0
      CALL KinLib_givpair(qmsene,Mas1,Mas2,m_q1, m_q2 ,betn,eta1,eta2)! real
      CALL KinLib_givpair(qmsene,amcru,amcru,m_r1,m_r2,betc,etc1,etc2)! ghost
*-----------------------------------------------------------------------
*     Mass weight for theta distribution
*-----------------------------------------------------------------------
*     Mass weight compensates for s'->s and droping terms -m**2/(k.q)**2
*     Care is taken of machine rounding errors.
*     del1 and del2 RECALCULATED out of angles sgx(i),cgx(i)
*     with TRUE sprim, using EXACT formulas
      amd1 = (Mas1/ener)**2
      amd2 = (Mas2/ener)**2
      DO i=1,m_nphot
         IF( cgx(i) .GT. 0.d0 ) THEN
            del1 = amd1/(eta1+betn) +betn*sgx(i)**2/(1+cgx(i))
            del2 = eta2 +betn*cgx(i)
         ELSE
            del1 = eta1 -betn*cgx(i)
            del2 = amd2/(eta2+betn) +betn*sgx(i)**2/(1-cgx(i))
         ENDIF
         dist1=1d0/(del1*del2) 
     $        *(1d0 -(amd1+amd2)/4d0
     $                   -amd2/4d0*del1/del2 -amd1/4d0*del2/del1)
         WtMlist(i)= dist1/dis0(i)
         IF(WtMlist(i) .LT.  1.d-90) WtMlist(i)= 0.d0
***********
** dist1x below is exactly the same as dist1 but for small masses is 
** prone to severe rounding errors (in contrast to dist1 used above)
*         IF((1-sprim/svar) .gt. 0.01d0) THEN
*            qf1qf2= m_q1(4)*m_q2(4) -m_q1(3)*m_q2(3) -m_q1(2)*m_q2(2) -m_q1(1)*m_q2(1)
*            qf1k = m_q1(4)*m_phot(i,4)-m_q1(3)*m_phot(i,3) -m_q1(2)*m_phot(i,2)-m_q1(1)*m_phot(i,1)
*            qf2k = m_q2(4)*m_phot(i,4)-m_q2(3)*m_phot(i,3) -m_q2(2)*m_phot(i,2)-m_q2(1)*m_phot(i,1)
*            dist1x = 2*qf1qf2/qf1k/qf2k -Mas1**2/qf1k**2 -Mas2**2/qf2k**2
*            dist1x = dist1x/4d0*m_phot(i,4)**2
*            WRITE(*,'(a,5f20.10)') '===>: ',dist1x/dist1
*         ENDIF
***********
*     finaly define Sudakov variables
         m_yfin(i)=del1*xk(i)/2d0
         m_zfin(i)=del2*xk(i)/2d0
      ENDDO
*-----------------------------------------------------------------------
* Transform from rest frame of Q=qf1+qf2 down to CMS=Lab,
* through inetrmediate rest frame of PX=qf1+qf2+phsu.
*-----------------------------------------------------------------------
      CALL KarFin_Kinf1(PX,m_q1,m_q2,m_r1,m_r2,m_nphot,m_phot,m_phsu)
*-----------------------------------------------------------------------
* Calculate YFS formfactor (cut-off dependent part) and mass weights
* Optionally removing photons below emin from the list
      CALL KarFin_Piatek( Mas1,Mas2,CharSq,WtMlist, wt3)
*-----------------------------------------------------------------------
* Monitoring weights and other control variables,
* Non-essential for the MC generation itself.
      CALL GLK_Mfill(m_idyfs+64, 1d0*m_nphot  ,1d0)
*[[[      uu = 1d0-sprim/svar
*[[[      CALL GLK_Fil1(m_idyfs+31, uu  ,wctrl)
*[[[      CALL GLK_Fil1(m_idyfs+32, uu  ,  1d0)
* marked photons
      IF(m_nphot .GE. 1) THEN
         DO i=1,m_nphot
*[[[            ul= log10(m_phot(i,4)/m_emin)
*[[[            IF(Mark(i) .EQ. 1)   CALL GLK_Fil1(m_idyfs+20,   ul,1.d0)
            IF(Mark(i) .EQ. 1)   m_MarTot=m_MarTot+1
         ENDDO
      ENDIF
*-----------------------------------------------------------------------
* Final Monitoring weights
 1000 CONTINUE
      WtFin = wt1*wt2*wt3
      CALL GLK_Mfill(m_idyfs+61,     wt1  ,1d0)
      CALL GLK_Mfill(m_idyfs+62,     wt2  ,1d0)
      CALL GLK_Mfill(m_idyfs+63,     wt3  ,1d0)
*-----------------------------------------------------------------------
      RETURN
*-----------------------------------------------------------------------
* Event outside phase space (too hard photon)
 900  CONTINUE
      wt1   = 0.d0
      wt2   = 1.d0
      wt3   = 1.d0
      m_nphot = 0
      GOTO 1000   !!!! a litle bit clumsy jump
*
      END                       !!!KarFin_YFSfin!!!


      SUBROUTINE KarFin_Piatek( Mas1,Mas2,CharSq,WtMlist, Wt3)
*/////////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                         //
*// Written CERN, piatek, 22 sept. 1989  (S.J.)                                             //
*// This routine calculates YFS form-factor and optionaly                                   //
*// removes some of soft photons (below Emin in CMS), calculating compensating weight.      //
*// Note the action of this routine is not Loretnz invariant !!!!                           //
*// Input:                                                                                  //
*//     KeyPia   = 0, NO removal of photons below Emin                                      //
*//              = 1, with removal of photons below Emin                                    //
*//     Mas1,2  = fermion masses                                                           //
*//     delta    = infrared cut-off in generation (dimensionless)                           //
*//     CharSq   = FSR charge squared                                                       //
*//     alfinv   = 1/alpha_QED                                                              //
*//     qf1,2    = fermion momenta                                                          //
*//     qf1c,2c  = ghost fermion momenta in crude S-factor                                  //
*//     phsu     = sum of photon momenta                                                    //
*//     phot     = list of photon momenta                                                   //
*//     WtMlist  = list of mass-weights for all photons                                     //
*// OUTPUT:                                                                                 //
*//     WtRem    = mass-weight of removed photons, for tests                                //
*//              = 1,  for KeyPia=0                                                         //
*//     WtMas    = mass-weight for all photons, see comments below,                         //
*//                the MOST IMPORTANT OUTPUT of Piatek !!!                                  //
*//     WtCtrl    = control-weight for remowed photons, for tests                           //
*//                                                                                         //
*/////////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                         //
*// REMARK on r1, r2:                                                                       //
*// They are ghost-fermion momenta in the definition of the crude                           //
*// distribution, i.e. truncated S-factor is r1*r2/(k*r1)/(k*r2)                            //
*// as generated in the crude MC.                                                           //
*// In QMS frame 4-momenta r1, r2 have the same energy                                      //
*// and directions as q1, q2                                                                //
*// but different masses Mas1c, Mas2c and therefore longer 3-momenta.                       //
*// They are usefull because we may use the same analytical                                 //
*// expression for the exact Breal and  for 'crude Breal'                                   //
*//                                                                                         //
*/////////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
*
      DOUBLE PRECISION    Mas1,Mas2,CharSq,WtMlist(100),Wt3
*
      DOUBLE PRECISION  pi
      PARAMETER(pi=3.1415926535897932d0)
      INTEGER  i,j,k,nph
      DOUBLE PRECISION    QQ(4),PX(4)
      DOUBLE PRECISION    amc2,Mass,Delta1,Epsi1,wtm1,wtm2
      DOUBLE PRECISION    svarX,svarQ,QQk,q1q2,betc
      DOUBLE PRECISION    DelB,DelB2,DelB2u,DelVol,DelYFS
      DOUBLE PRECISION    BVR_A
      DOUBLE PRECISION    Eq1,Eq2,EminQ,EQQ,r1r2,Mas1c,Mas2c
      DOUBLE PRECISION    BtiXcru,BtiQcru,BtiQexa,BtiXexa
      DOUBLE PRECISION    BVR_Btildc,BVR_Btilda
      DOUBLE PRECISION    WtRem,WtCtrl
      DOUBLE PRECISION    VoluMC, YFS_IR, YFSkon
      DOUBLE PRECISION    alfpi,alfch,alfCR
*-----------------------------
      INTEGER  iCont
      DATA     iCont /0/
*-----------------------------
      WtCtrl =1d0
      WtRem  =1d0
      alfpi = 1d0/pi/m_AlfInv
      alfch = alfpi*CharSq
      alfCR = alfch* m_Xenph
      Mass = MIN(Mas1,Mas2)
      DO k=1,4
         PX(k) = m_q1(k)+m_q2(k) +m_phsu(k)
         QQ(k) = m_q1(k)+m_q2(k)
      ENDDO
      svarX = PX(4)**2-PX(3)**2 -PX(2)**2 -PX(1)**2
      svarQ = QQ(4)**2-QQ(3)**2 -QQ(2)**2 -QQ(1)**2
      EQQ   = 0.5d0*SQRT(svarQ)
      q1q2  =  m_q1(4)*m_q2(4) -m_q1(3)*m_q2(3) -m_q1(2)*m_q2(2) -m_q1(1)*m_q2(1)
      r1r2  =  m_r1(4)*m_r2(4) -m_r1(3)*m_r2(3) -m_r1(2)*m_r2(2) -m_r1(1)*m_r2(1)
      QQk   =  QQ(4)*m_phsu(4) -QQ(3)*m_phsu(3) -QQ(2)*m_phsu(2) -QQ(1)*m_phsu(1)
      Eq1   = (svarQ +Mas1**2 -Mas2**2)/(2*SQRT(svarQ))
      Eq2   = (svarQ -Mas1**2 +Mas2**2)/(2*SQRT(svarQ))
* Delta1 and Epsi1 are cutoffs located in YFS formfactor
* Note that xfact=(1+2*QQk/svarQ) is EXACTLY the same as in KarFin_YFSfin
* Delta1 = 2*Emin/sqrt(s'), where Emin is in QMS
      Delta1 = m_Delta*(1+ 2*QQk/svarQ)
      Epsi1  = DSQRT(m_Emin**2/m_q1(4)/m_q2(4))
      EminQ  = EQQ*Delta1
* The total phase space integral for crude x-section and YFS formfactor cut-off dependend part
* Note that delta is a lower limit on y-z-variables in crude generation
      amc2  =  4d0*Mass**2/svarX
      betc  =  DSQRT(1d0-amc2)
*     ===========================================================================================
      VoluMC =  2d0*alfCR *(1d0+betc**2)/(2d0*betc)* DLOG((1d0+betc)**2/amc2) *DLOG(1/m_Delta)
      YFS_IR = -2d0*alfch *(     q1q2 *BVR_A( q1q2, Mas1,Mas2) -1d0  )        *DLOG(1/Delta1)
* YFSkon is delegated/exported to QED3 (unused here).
      YFSkon =  1/4d0 *2*alfch*(DLOG(svarQ/Mass**2)-1) + alfch*( -.5d0  +pi**2/3d0) ! Mass<<sqrt(s)
*     ===========================================================================================
* Corrections necessary for photon remooval scenario (now default!)
      Mas1c  = Mass*SQRT(svarQ/svarX)
      Mas2c  = Mass*SQRT(svarQ/svarX)
      BtiXcru= BVR_Btildc(alfCR, r1r2, m_r1(4),m_r2(4), Mas1c,Mas2c, m_Emin, m_MasPhot) !crude
      BtiQcru= BVR_Btildc(alfCR, r1r2, EQQ,    EQQ,     Mas1c,Mas2c, EminQ,  m_MasPhot) !crude
      BtiXexa= BVR_Btilda(alfch, q1q2, m_q1(4),m_q2(4), Mas1, Mas2,  m_Emin, m_MasPhot) !exact
      BtiQexa= BVR_Btilda(alfch, q1q2, Eq1,    Eq2,     Mas1 ,Mas2,  EminQ,  m_MasPhot) !exact
      DelVol = BtiXcru -BtiQcru   ! positive
      DelYFS = BtiXexa -BtiQexa   ! positive
      DelB2  = -DelVol+DelYFS
*------ oldies ------
* Total QMS-CMS, result should be negative ( delta<<epsilon )
      DelB  =   VoluMC +YFS_IR
* Ultrarelativistic (small mass) old version is the following:
      DelB2u = -2*alfch*(DLOG(svarX/svarQ)+1) *dlog(Epsi1/Delta1)
* The average mass-weight for removed photon = exp(DelB2)
* It can be calculated analyticaly as a  ratio of YFS formfactors
* On the other hand, it is checked by MC, see control weight WtCtrl.
*********************************************************************************************
*      IF(iCont.LE.10 ) THEN
*         IF((1-svarQ/svarX) .GT. 0.1d0) THEN
*            iCont = iCont+1
*((( Approximate version of DelB2 without contr. with A4 terms for tests (surpisingly good!!!)
*      DelB2w = -2*alfch*(  (1d0+betc**2)/(2d0*betc)* DLOG((1d0+betc)**2/amc2) 
*     $                    -( q1q2*BVR_A(q1q2,Mas1,Mas2) -1d0 )
*     $                  ) *DLOG(Epsi1/Delta1)
*         WRITE(*,'(a,5f20.10)') 'piatek: ',1-svarQ/svarX,DelB2,DelB2u/DelB2,DelB2w/DelB2
*)))
*            VoluMC2 = 
*     $            BVR_Btildc(alfch, r1r2, EQQ,EQQ,  Mas1c,Mas2c, EQQ,          m_MasPhot)
*     $           -BVR_Btildc(alfch, r1r2, EQQ,EQQ,  Mas1c,Mas2c, EQQ*m_Delta,  m_MasPhot)
*            WRITE(*,'(a,5f20.10)') '###Piatek: ',1-svarQ/svarX, VoluMC2,VoluMC,VoluMC2/VoluMC
*         ENDIF
*      ENDIF
*********************************************************************************************
*--------------------------
      wtm1=1d0
      wtm2=1d0
* mass weight below and above Emin calculated separately
      DO i=1,m_nphot
         IF(m_phot(i,4) .LT. m_Emin) THEN
            wtm1=wtm1*WtMlist(i) /m_Xenph
            IF(wtm1 .LE. 1d-90) wtm1=0d0
         ELSE
            wtm2=wtm2*WtMlist(i) /m_Xenph
            IF(wtm2 .LE. 1d-90) wtm2=0d0
         ENDIF
      ENDDO
*------------------------------------------------------------------------------
* Control weight - its average should be =1 within statist. error!
      WtCtrl =wtm1*exp(-DelB2)
      IF(m_KeyPia .EQ. 0) THEN
         IF( ABS(DelB) .GT. 100d0 ) WRITE(*,*) '#### KarFin_Piatek: DelB= ',DelB
         WtRem    = 1d0
         m_WtMass = wtm1*wtm2      !!! <--removal OFF
         m_VoluMC = EXP( VoluMC)
         m_YFS_IR = EXP( YFS_IR)
         m_YFSkon = EXP( YFSkon)   !!! <--finite part of YFS formfactor
*        =====================================
      ELSE
* Optional removal of photons below Emin from the record
* in such a case WtMas includes exp(belb2)= <wt3> for removed ph.
         nph=m_nphot
         DO j=m_nphot,1,-1
            IF(m_phot(j,4) .LT. m_Emin) THEN
               DO i=j+1,nph
                  DO k=1,4
                     m_phot(i-1,k)=m_phot(i,k)
                  ENDDO
               ENDDO
               nph=nph-1
            ENDIF
         ENDDO
* Correction of Alex Read, probably obsolete because KarFin_Merge is also corrected
         DO j=nph+1,m_nphot
            DO k=1,4
               m_phot(j,k) = 0.d0
            ENDDO
         ENDDO
         m_nphot=nph
* Important to remember: EXP(YFS_IR +YFSkon)         = full YFS formfactor in QMS
*                        EXP(YFS_IR +YFSkon +DelYFS) = full YFS formfactor in CMS
         WtRem    = wtm1
         m_WtMass = wtm2           !!! <--removal ON
         m_VoluMC = EXP(VoluMC -DelVol)
         m_YFS_IR = EXP(YFS_IR +DelYFS)
* YFSkon is delegated/exported to QED3 (unused here).
         m_YFSkon = EXP(YFSkon)    !!! <--finite part of YFS formfactor in QMS frame!!!
*        ===============================================
      ENDIF
*     **********************************
      Wt3 = m_WtMass *m_VoluMC *m_YFS_IR
*     **********************************
* Monitoring
      CALL GLK_Mfill(m_idyfs+65,       WtCtrl ,1d0)
      CALL GLK_Mfill(m_idyfs+66,       WtRem  ,1d0)
*
      END                       !!! KarFin_Piatek !!!



      SUBROUTINE KarFin_PoissGen(average,nmax,mult,rr)
*/////////////////////////////////////////////////////////////////////////////////////////////
*//  Last corr. nov. 91                                                                     //
*//   This generates photon multipl. nphot according to poisson distr.                      //
*//   Input:  average = average multiplicity                                                //
*//           nmax  = maximum multiplicity                                                  //
*//   Output: mult = generated multiplicity                                                 //
*//           rr(1:100) list of ordered uniform random numbers,                             //
*//           a byproduct result, to be eventually used for some further                    //
*//           purpose (i.e.  generation of photon energies).                                //
*/////////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      DOUBLE   PRECISION average,  rr(*)
      INTEGER  nmax,mult
* locals
      REAL                 rvec(10)
      DOUBLE   PRECISION   sum,y
      INTEGER              nn,it,nfail
      DATA nfail/0/
*---------------------------------------------------------------------------------------------
      IF( average.LE.0 ) GOTO 900
 50   nn=0
      sum=0d0
      DO it=1,nmax
         CALL PseuMar_MakeVec(rvec,1)
         y= log(rvec(1))
         sum=sum+y
         nn=nn+1
         rr(nn)=sum/(-average)
         IF(sum .LT. -average) GOTO 130
      ENDDO
      nfail=nfail+1
      IF(nfail .GT. 10) GOTO 900
      GOTO 50
 130  mult=nn-1
      RETURN
 900  CONTINUE
      WRITE(*,*) ' STOP in KarFin_PoissGen: nmax,average= ',nmax,average
      STOP
      END

      SUBROUTINE KarFin_AngBre(am2,del1,del2,costhg,sinthg,dist0,dist1)
*/////////////////////////////////////////////////////////////////////////////////////////////
*//   This routine generates photon angular distribution                                    //
*//   in the rest frame of the fermion pair.                                                //
*//   The distribution is the S-factor without mass term,                                   //
*//   i.e. without terms 2p_1p_2/(kp_1)(kp_2)                                               //
*//   Fermion mass is treated exactly!                                                      //
*//   INPUT:                                                                                //
*//       am2 = 4*massf**2/s where massf is fermion mass                                    //
*//       and s is effective mass squared of the parent fermion-pair.                       //
*//   OUTPUT:                                                                               //
*//       del1= 1-beta*cos(theta)                                                           //
*//       del2= 1+beta*cos(theta)                                                           //
*//       costhg, sinthg, cos and sin of the photon                                         //
*//       angle with respect to fermions direction                                          //
*//       dist0 = distribution generated, without m**2/(kp)**2 terms                        //
*//       dist1 = distribution with m**2/(kp)**2 terms                                      //
*/////////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT  NONE
      DOUBLE PRECISION am2,del1,del2,costhg,sinthg,dist0,dist1
* Locals
      DOUBLE PRECISION  a,eps,beta
      REAL              rn(10)
*---------------------------------------------------------------------------------------------
      CALL PseuMar_MakeVec(rn,2)
      beta =SQRT(1.d0-am2)
      eps  =am2/(1.d0+beta)                     != 1-beta
      del1 =(2.d0-eps)*(eps/(2.d0-eps))**rn(1)  != 1-beta*costhg
      del2 =2.d0-del1                           != 1+beta*costhg
* calculation of sin and cos theta from internal variables
      costhg=(del2-del1)/(2*beta)               ! exact
      sinthg=SQRT(del1*del2-am2*costhg**2)      ! exact
* symmetrization
      IF(rn(2) .LE. 0.5d0) THEN
        a=del1
        del1=del2
        del2=a
        costhg= -costhg
      ENDIF
      dist0=1d0/(del1*del2)*(1d0 -am2/2d0)
      dist1=1d0/(del1*del2) 
     $     *(1d0 -am2/2d0 -am2/4d0*(del1/del2+del2/del1))
* totaly equivalent formula is the following
*     dist1=1d0/(del1*del2)   *beta*sinthg**2/(del1*del2)
      END


      SUBROUTINE KarFin_Kinf1(PX,q1,q2,q1c,q2c,nphot,phot,phsu)
*/////////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                         //
*// Transforms final fermions and photons from QMS through Z-frame to CMS.                  //
*// Random Euler rotation is applied in the intermediate PX frame (Z frame)                 //
*//                                                                                         //
*/////////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
*
      DOUBLE PRECISION   PX(4),phot(100,4),phsu(4),q1(4),q2(4),q1c(4),q2c(4)
      DOUBLE PRECISION   ph(4),qqk(4)
      INTEGER i,k,nphot

* Calculate total final state four-momentum ferms+phots
      DO k=1,4
         qqk(k)=q1(k)+q2(k)+phsu(k)
      ENDDO
* Transform fermions
      CALL KarFin_BostEul(-1,qqk,PX,q1,q1) ! <-- Initialize Euler angles!!!
      CALL KarFin_BostEul( 0,qqk,PX,q2,q2)
      CALL KarFin_BostEul( 0,qqk,PX,q1c,q1c)
      CALL KarFin_BostEul( 0,qqk,PX,q2c,q2c)
      CALL KarFin_BostEul( 0,qqk,PX,phsu,phsu)
* Transform photons
      DO i=1,nphot
         DO k=1,4
            ph(k)=phot(i,k)
         ENDDO
         CALL KarFin_BostEul( 0,qqk,PX,ph,ph)
         DO k=1,4
            phot(i,k)= ph(k)
         ENDDO
      ENDDO
      END

      SUBROUTINE KarFin_BostEul(mode,qqk,PX,pvec,qvec)
*/////////////////////////////////////////////////////////////////////////////////////////////
*//                                                                                         //
*// Three transformations:                                                                  //
*// (1) Boost from final fermions rest frame to ferms+phots rest frame (Z frame).           //
*// (2) Euler rotation erasing memory of fermion directions.                                //
*// (3) Boost to laboratory system.                                                         //
*// Note that this transformation generates 'flat 2-body Born' in PX frame.                 //
*// This is very easy to implement.                                                         //
*// Of course, more sophisticated angular distribution can be implemented.                  //
*// In such a case BostEul will be replaced with some other transformation.                 //
*// Photon removal procedure with Piatek will work for arbitrary transformation.            //
*//                                                                                         //
*// Note that the procedure is encapsulated functionaly by generating Euler                 //
*// angles localy, during first call for mode=-1.                                           //
*//                                                                                         //
*/////////////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INTEGER mode
      DOUBLE PRECISION  pvec(4),qvec(4),qqk(4),PX(4)
      DOUBLE PRECISION  pi
      PARAMETER(pi=3.1415926535897932d0)
      REAL              rvec(10)
      DOUBLE PRECISION  the,phi,cth
      SAVE              the,phi,cth
*
      IF(mode .EQ. -1) THEN
* Angles for Euler rotation
         CALL PseuMar_MakeVec(rvec,2)
         cth= 1.d0 -2.d0*rvec(1)
         the= acos(cth)
         phi= 2.d0*pi*rvec(2)
      ENDIF
* And the transformations
      CALL KinLib_BostQ(  1,qqk,pvec,qvec) ! Boost along qqk
      CALL KinLib_Rotor(2,3,the,qvec,qvec) ! Rotation y-z
      CALL KinLib_Rotor(1,2,phi,qvec,qvec) ! Rotation x-y
      CALL KinLib_BostQ( -1,PX,qvec,qvec)  ! Boost to CMS
      END


      SUBROUTINE KarFin_ZBoostAll(exe)
*///////////////////////////////////////////////////////////////////////////////
*//                                                                           //
*//   performs z-boost on all momenta of the event                            //
*//   this z-boost corresponds to beamstrahlung or beamspread                 //
*//   and is done at the very end of generation, after m.el. calculation      //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
      DOUBLE PRECISION  exe
      INTEGER           j,k
      DOUBLE PRECISION  ph(4)
*
      IF( exe.EQ. 1d0) RETURN
      CALL KinLib_Boost(3,exe,m_q1,m_q1)
      CALL KinLib_Boost(3,exe,m_q2,m_q2)
      CALL KinLib_Boost(3,exe,m_phsu,m_phsu)
      DO j=1,m_nphot
         DO k=1,4
            ph(k) = m_phot(j,k)
         ENDDO
         CALL KinLib_Boost(3,exe,ph,ph)
         DO k=1,4
            m_phot(j,k) = ph(k)
         ENDDO
      ENDDO
      END


      SUBROUTINE KarFin_GetIsFSR(IsFSR)
*///////////////////////////////////////////////////////////////////////////////
*//                                                                           //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
      INTEGER IsFSR
      IsFSR= m_IsFSR
      END

      SUBROUTINE KarFin_GetKeyPia(KeyPia)
*///////////////////////////////////////////////////////////////////////////////
*//                                                                           //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
      INTEGER KeyPia
      KeyPia= m_KeyPia
      END

      SUBROUTINE KarFin_GetNphot(nphot)     
*///////////////////////////////////////////////////////////////////////////////
*//                                                                           //
*//   Get photon multiplicity                                                 //
*//                                                                           //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
      INTEGER nphot
*
      nphot = m_nphot
      END

      SUBROUTINE KarFin_GetPhoton1(iphot,phot)
*///////////////////////////////////////////////////////////////////////////////
*//                                                                           //
*//   get i-th photon momentum                                                //
*//                                                                           //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
      INTEGER iphot
      DOUBLE PRECISION   phot(4)
      INTEGER k
*
      DO k=1,4
         phot(k) = m_phot(iphot,k)
      ENDDO
      END

      SUBROUTINE KarFin_GetSudakov1(iphot,y,z)
*///////////////////////////////////////////////////////////////////////////////
*//                                                                           //
*//   Get sudakovs of i-th photon                                             //
*//                                                                           //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
      INTEGER iphot
      DOUBLE PRECISION   y,z
*
      y = m_yfin(iphot)
      z = m_zfin(iphot)
      END


      SUBROUTINE KarFin_GetSudakov(nphot,yfin,zfin)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
*
      INTEGER nphot
      DOUBLE PRECISION   yfin(*),zfin(*)
      INTEGER i
*-------------
      nphot = m_nphot
      DO i=1,m_nphot
         yfin(i) = m_yfin(i)
         zfin(i) = m_zfin(i)
      ENDDO
      END

      SUBROUTINE KarFin_GetPhotons(nphot,sphot)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
*
      INTEGER nphot
      DOUBLE PRECISION   sphot(m_npmx,4)
      INTEGER j,k
*-------------
      nphot = m_nphot
      DO j=1,m_nphot
         DO k=1,4
            sphot(j,k) = m_phot(j,k)
         ENDDO
      ENDDO
      END

      SUBROUTINE KarFin_GetFermions(qf1,qf2)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
*
      DOUBLE PRECISION   qf1(4),qf2(4)
      INTEGER k
*---------------
      DO k=1,4
         qf1(k) = m_q1(k)
         qf2(k) = m_q2(k)
      ENDDO
      END

      SUBROUTINE KarFin_GetSvarQ(SvarQ)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
*
      DOUBLE PRECISION   SvarQ
*---------------
      SvarQ = ( m_q1(4)+m_q2(4) )**2 -( m_q1(3)+m_q2(3) )**2 
     $       -( m_q1(2)+m_q2(2) )**2 -( m_q1(1)+m_q2(1) )**2
      END

      SUBROUTINE KarFin_GetYFSkon(YFSkon)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//   Get finite part of YFS form-factor                                     //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
      DOUBLE PRECISION    YFSkon
*------------------
      YFSkon = m_YFSkon
      END

      SUBROUTINE KarFin_GetYFS_IR(YFS_IR)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//   Get finite part of YFS form-factor                                     //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
      DOUBLE PRECISION    YFS_IR
*------------------
      YFS_IR = m_YFS_IR
      END

      SUBROUTINE KarFin_WtMass(WtMass)
*///////////////////////////////////////////////////////////////////////////////
*//                                                                           //
*///////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
      DOUBLE PRECISION   WtMass
*
      WtMass = m_WtMass
      END

      SUBROUTINE KarFin_SetEmin(Emin)
*/////////////////////////////////////////////////////////////////////////////////////
*//   Photon minimum energy in LAB system                                           //
*/////////////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
      DOUBLE PRECISION   Emin
*
      m_Emin = Emin
      END

      SUBROUTINE KarFin_Print(iev,ie1,ie2)
*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*// Prints out four momenta of FINAL state                                   //
*// and the serial number of event iev on unit m_out                         //
*//                                                                          //
*//////////////////////////////////////////////////////////////////////////////
      IMPLICIT NONE
      INCLUDE 'KarFin.h'
      INTEGER  iev,ie1,ie2
      DOUBLE PRECISION    sphum(4)
      CHARACTER*8 txt
      DOUBLE PRECISION    sum(4),ams,amph,amf1,amf2
      INTEGER  i,k
*--------------------------------------------------------
      IF( (iev .GE. ie1) .AND. (iev .LE. ie2) ) THEN
         txt = '  KarFin '
         WRITE(m_out,*) 
     $        '=========== ',txt,' ======================>',iev
         amf1 = m_q1(4)**2-m_q1(3)**2-m_q1(2)**2-m_q1(1)**2
         amf1 = sqrt(abs(amf1))
         amf2 = m_q2(4)**2-m_q2(3)**2-m_q2(2)**2-m_q2(1)**2
         amf2 = sqrt(abs(amf2))
         WRITE(m_out,3100) 'qf1',(  m_q1(  k),k=1,4),amf1
         WRITE(m_out,3100) 'qf2',(  m_q2(  k),k=1,4),amf2
         DO i=1,m_nphot
            amph = m_phot(i,4)**2-m_phot(i,3)**2 -m_phot(i,2)**2-m_phot(i,1)**2
            amph = sqrt(abs(amph))
            WRITE(m_out,3100) 'pho',(m_phot(i,k),k=1,4),amph
         ENDDO
         DO k=1,4
            sum(k)=m_q1(k)+m_q2(k)
         ENDDO
         DO i=1,m_nphot
            DO k=1,4
               sum(k)=sum(k)+m_phot(i,k)
            ENDDO
         ENDDO
         ams = sum(4)**2-sum(3)**2-sum(2)**2-sum(1)**2
         ams = sqrt(abs(ams))
         WRITE(m_out,3100) 'sum',(  sum(  k),k=1,4),ams
      ENDIF
 3100 FORMAT(1x,a3,1x,5f20.14)
      END


*//////////////////////////////////////////////////////////////////////////////
*//                                                                          //
*//                      End of CLASS  KarFin                                //
*//////////////////////////////////////////////////////////////////////////////
