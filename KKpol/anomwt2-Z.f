c
c     Code provided by Z. Was, June 2024
c     except for subroutine fillborn, provided by G Ganis, CERN, Aug 2024
c
      
c       subroutine filltaupair(hvt1, hvt2, hvc1, hvc2, bet1, alf1, gam1, bet2, alf2, gam2, ph1, th1, ph2, th2,
c      +                       lout, isinit, idyfs, keyclone, ifphot)

c       DOUBLE PRECISION htv1(4), hvt2(4), hvc1(4), hvc2(4), bet1, alf1, gam1, bet2, alf2, gam2
c       DOUBLE PRECISION ph1, th1, ph2, th2
c       INTEGER lout, isinit, idyfs, keyclone, ifphot
c       INTEGER j

c       INCLUDE 'Taupair.h'
       
c       do j = 1, 4
c c    m_HvecTau1(j) = hvt1(j)
c c    m_HvecTau2(j) = hvt2(j)
c          m_HvClone1(j) = hvc1(j)
c          m_HvClone2(j) = hvc2(j)
c       enddo
c       m_beta1 = bet1
c       m_alfa1 = alf1
c       m_gamma1 = gam1
c       m_beta2 = bet2
c       m_alfa2 = alf2
c       m_gamma2 = gam2
c       m_phi1 = ph1
c       m_thet1 = th1
c       m_phi2 = ph2
c       m_thet2 = th2
c       m_out = lout
c       m_IsInitialized = isinit
c       m_idyfs = idyfs
c       m_KeyClone = keyclone
c       m_IFPHOT = ifphot

c       return
c       end

c       subroutine anomwt(iqed, Ar0, Ai0, Br0, Bi0, wtME,wtSPIN,wtSPIN0)
c C ####################################################################      
c       ! anomalous weight calculation, alternatively variants of R matrix
c       ! calculation old without Z exchange and new with, can be used.
c       ! electroweak loops included too
c       !
c       ! Frame orientation is refined with respect to old work
c       ! previously some symmetries of matrix elements were unintentionally
c       ! used. ZBW: May 4, 2023
c C ######################################################################

c       INCLUDE 'Taupair.h'
c c *
c        double precision H1(4),H2(4),P1(4),P2(4),QQ(4),PB1(4),PB2(4),pbb(4)
c        double precision R(4,4),R0(4,4),E, theta, Ar, Ai, Br, Bi,X
c        double precision ReA, ImA, ReB, ImB, ReX, ImX, ReY, ImY
c        double precision Ar0, Ai0, Br0, Bi0,sign(4)
c       REAL*8 wtME,wtSPIN,wtSPIN0
c       COMMON /HVEC/ H1,H2
c       REAL HR1(4),HR2(4),AM
c       Do k=1,4
c          H1(k)=m_HvClone1(k) ! differ from m_HvecTau1(k)
c          H2(k)=m_HvClone2(k)    !same as m_HvecTau2(k)
c          HR1(K)=H1(K)
c          HR2(K)=H2(K)
c          P1(k)=0.0
c          P2(k)=0.0
c       enddo
c       P1(4)=1.777
c       P2(4)=1.777

c       CALL GPS_tralorPrepare(p1,1)
c       CALL GPS_tralorPrepare(p2,2)
c       CALL GPS_TralorDoIt(1,P1,P1)
c       CALL GPS_TralorDoIt(2,P2,P2)
c       CALL GPS_TralorDoIt(1,H1,H1)
c       CALL GPS_TralorDoIt(2,H2,H2)
c c      CALL GPS_tralorPrepare(p1,1)
c c      CALL GPS_tralorPrepare(p2,2)
c       CALL Tralo4(1,HR1,HR1,AM)
c       CALL Tralo4(2,HR2,HR2,AM)

      
c       DO k=1,4
c          QQ(k)=P1(k)+P2(k)
c          PB1(K)=0.0
c          PB2(K)=0.0
c       enddo
c       PB1(4)= 1.0
c       PB2(4)= 1.0
c       PB1(3)= 1.0
c       PB2(3)=-1.0
      
c c     we go to restframe of tau pair to work against bremstrahlungs
      
c       call bostdq(1,qq,H1,H1)
c       call bostdq(1,qq,H2,H2)
c       call bostdq(1,qq,p1,p1)
c       call bostdq(1,qq,p2,p2)
c       call bostdq(1,qq,pb1,pb1)
c       call bostdq(1,qq,pb2,pb2)
      
c C     we eliminate y-components of p1,p2
c       FI=ANGFI(P2(1),P2(2))

c       call rotod3(-fi,h1,h1)
c       call rotod3(-fi,h2,h2)
c       call rotod3(-fi,p1,p1)
c       call rotod3(-fi,p2,p2)
c       call rotod3(-fi,pb1,pb1)
c       call rotod3(-fi,pb2,pb2)

c C set taus along z direction      
c       thet=angxy(p2(3),p2(1))

c       call rotod2(-thet,h1,h1)
c       call rotod2(-thet,h2,h2)
c       call rotod2(-thet,p1,p1)
c       call rotod2(-thet,p2,p2)
c       call rotod2(-thet,pb1,pb1)
c       call rotod2(-thet,pb2,pb2)

c C     set  beam difference to reaction plane  x-z
c       do k=1,4
c          PBB(k)=pb1(k)-pb2(k)
c           PBB(k)=pb2(k)-pb1(k)
c       enddo
      
c       FI1=ANGFI(PBB(1),PBB(2))

c       call rotod3(-fi1,pb1,pb1)
c       call rotod3(-fi1,pb2,pb2)
c       call rotod3(-fi1,h1,h1)
c       call rotod3(-fi1,h2,h2)
c       call rotod3(-fi1,p1,p1)
c       call rotod3(-fi1,p2,p2)
                            
c       call bostdq(1,p1,H1,H1)
c       call bostdq(1,p2,H2,H2)

c       E=p1(4)
c       theta=acos(-pb1(3)/sqrt(pb1(1)**2+pb1(2)**2+pb1(3)**2))  ! minus necessary for KKMC frames
      
c       ReA=0.0
c       ImA=0.0
c       ReB=0.0
c       ImB=0.0
c       ReX=0.0
c       ImX=0.0
c       ReY=0.0
c       ImY=0.0

c       call DipolQQRijRadCor (0,E, theta, ReA, ImA, ReB, ImB, ReX, ImX, ReY, ImY, R0,1)

c       iqed=1
c       ReA=Ar0!*0.0
c       ImA=Ai0!*0.0
c       ReB=Br0!*0.0           
c       ImB=Bi0!*0.0
c       ReX=Ar0!*0.0
c       ImX=Ai0!*0.0
c       ReY=Br0!*0.0
c       ImY=Bi0!*0.0

c       call DipolQQRijRadCor (iqed,E, theta, ReA, ImA, ReB, ImB, ReX, ImX, ReY, ImY, R,1)

c       wtME=R(4,4)/R0(4,4)
c       wtSPIN=0.0
c       wtSPIN0=0.0

c       ! sign() introduces overall rotation around by pi around  y axis. Necessary for KKMC frames. 
c       sign(1)=-1.0
c       sign(2)= 1.0
c       sign(3)=-1.0
c       sign(4)= 1.0
      
c       DO I=1,4
c          DO J=1,4
c             WTSPIN =WTSPIN +R (I,J)/R(4,4)*H1(I)*H2(J)*sign(I)*sign(j)
c             WTSPIN0=WTSPIN0+R0(I,J)/R0(4,4)*H1(I)*H2(J)*sign(I)*sign(j)
c          ENDDO
c       ENDDO
     
c       WTSPIN =WTSPIN/WTSPIN0

c       end

      subroutine fillbornv(keyelw, swsq, alfinv, mz, gammz, GSW)
c Routine provided by G Ganis, CERN, Aug 2024
      
      INCLUDE '../dizet/BornV.h'
      INTEGER keyelw, j
      DOUBLE PRECISION swsq, alfinv, mz, gammz, GSW(m_poinG)

      m_KeyElw = keyelw;
      m_swsq = swsq;
      m_AlfInv = alfinv;
      m_MZ = mz;
      m_GammZ = gammz;
      do j = 1, m_poinG;
         m_GSW(j) = GSW(J);
      enddo

      return
      end
	  
      SUBROUTINE DipolQQRijRadCor(iqed,Energy, theta, ReA0, ImA0, ReB, ImB, ReX, ImX, ReY, ImY, R, channel)
C     version of dipole routine with electroweak form factors included, 
C     If ported to different applications, beware of electroweak schemes, constant
C     and form-factors initialization
C     qed anomalous magnetic moment coded, but commented out      
C     (comments)  #####################################
C     1) code is adopted to run  with KKMC
C     2) and EW for channel=1  (initial state e, final state tau or mu.
C     3) internal activation of electroweak KKMC libs with IfGSW=KeyElw, call BornV_GetKeyElw(KeyElw)
C     4) note distinct conventions, papers:  e-Print: 2307.03526 versus Eur.Phys.J.C 79 (2019) 6, 480     
C        that is  about  photon/Z normalization too
C     ########################################################
      
c             PROCESS   fermion_i + fermionbar_i -> fermion_f- + fermionbar_f   
c  
c           IMPROVED BORN APPROXIMATION (IBA) with radiative corrections
c  Integer parameter channel = 1 for LEPTONS, 2 -for UP QUARK, 3 - for DOWN QUARK
c 
c  Photon and Z-boson exchanges are included. 
c  Calculates spin-correlation coefficients R(i,j) as functions of beam
c  energy: Energy = sqrt(s)/2 (in GeV), and scattaring angle theta.
c
c  Anomalous magnetic formfactor A = ReA +i*ImA, 
c  and electric dipole formfactor B = ReB + i*ImB are complex.
c  Also weak anomalous magnetic formfactor X = ReX + i*ImX and 
c  weak electric dipole formfactor Y = ReY + i*ImY are complex. 
c  
c  V is velocity of tau, gam is Lorentz factor, alpha is fine-structure constant,
c  G_F is Fermi constant of weak interaction.
c  RSM(i,j) - spin-correlation matrix in IBA without anomalous moments.
c  RDM(i,j) - spin-correlation matrix, linear in anomalous dipole moments A, B, X, Y.
c  R(i,j) - total spin-correlation matrix. 
c  Order of coefficients R(ij):  i = 1,2,3 correspond to S_x, S_y, S_z for tau-,
c                                j = 1,2,3 correspond to S'_x, S'_y, S'_z for tau+,
c                                i = j = 4 corresponds to 1 (no spin dependence).
c 
c   Functions below and their notation are from the paper: 
c   E. Richter-Was, Z. Was Eur. Phys.J. C (2019) 79, 480
c
c  indices i = 1, 2, 3 correspond to lepton, up quark, down quark, respectively 
c          Gamma_vp - vacuum polarization factor (function of s)
c          Rho_11, Rho_12, Rho_13 =rho_{lepton,i}(s,t) for i= L,Up,Down 	
c          K_1,K_2,K_3 = K_i(s,t)  for i= L, Up, Down
c          K_11,K_12,K_13  = K_{lepton,i}(s,t) for i = L,Up,Down
                                                                       
      Implicit none
  !    INCLUDE '../../basf2/KK2f/GPS.fi'
  !    INCLUDE '../../basf2/bornv/BornV.fi'
!      INCLUDE '../DZface/BornV.h'
      INCLUDE '../dizet/BornV.h'
cc not needed     external function K_1, K_2, K_3, K_11, K_12, K_13 
cc not needed     external function Rho_11,Rho_12,Rho_13,Gamma_vp
	  complex*16 A,B,X,Y, Pg, Pz, Den, vi, vf 
	  complex*16 K1,K2,K3,K11,K12,K13,Gammavp,Rho11,Rho12,Rho13
          DOUBLE COMPLEX     GSW(100)
          DOUBLE PRECISION Svar,CosThetD
       DOUBLE COMPLEX    RhoEW, VPgamma, CorEle, CorFin, CorEleFin, VVCef  	  
	  integer i, j, channel,iqed,KFf,IfGSW,KeyElw,IfPrint                  
      real*8 Energy,theta,ReA,ImA,ReA0,ImA0,ReB,ImB,ReX,ImX,ReY,ImY,gam,V
      real*8 ArQED, AiQED, Ar1, Ai1
	  real*8 e, f, m, Qi, Qf, ai, af, sw2, sw, cw2,cw, s, t, Fvivf  
      real*8 R(1:4, 1:4), RSM(1:4, 1:4), RDM(1:4, 1:4)
      real*8 PI/3.141592653589793238d0/,alpha/7.2973525693d-3/
	  real*8 Mz/91.1876d0/,Gz/2.4952d0/      ! Mass and width of Z in GeV
cc     real*8 v0/-0.03783d0/,a0/-0.50123d0/  ! vector and axial couplings of leptons
      real*8 mtau/1.77686d0/                 ! mass of tau in GeV

	  real*8 G_F/1.1663788d-5/                ! from PDG Fermi constant, GeV^{-2}
          real*8 G_mu/1.166389d-5/,Mw/80.353d0/ ! G_mu and M_W from Eur.Phys.J. C (2019) 79, 480
          
C     ==================================================================     
C     definition of form-factors as <in-line> functions. Complex form-factors
C     initialized from KKMC libraries will be used. Note that electroweak
C     scheme need to be reconsidered if function is used for other purposes.
      k3(s,t) = CorEle
      k1(s,t) = CorEle
      k13(s,t) = CorEleFin
      k2(s,t) = CorEle
      k12(s,t) = CorEleFin
      k11(s,t) = CorEleFin
      
      gammavp(s) = VPgamma
      rho11(s,t) = RhoEW
      rho12(s,t) = RhoEW
      rho13(s,t) = RhoEW
C     ===================================================================
      
      KFf=15                    ! 415-400
      Svar=4*Energy**2
      CosThetD=cos(theta)
!         RhoEW     = GSW(1)
!         VPgamma   = GSW(6)
!         CorEle    = GSW(2)
!         CorFin    = GSW(3)
!     CorEleFin = GSW(4)
      IfGSW=1                   ! switch to activate electroweak formfactor
C      call BornV_GetKeyElw(KeyElw)
C      IfGSW=KeyElw
      IfGSW = m_KeyElw  ! GG-Aug2024: from BornV.h
      IfPrint=0  ! switch to activate prints
      If(IfGSW.eq.1) then
c    CALL BornV_InterpoGSW(KFf,Svar,CosThetD)
c    CALL BornV_GetGSW(GSW)
c    RhoEW     = GSW(1)
c    VPgamma   = GSW(6)
c    VPgamma = 1.d0   /(2.d0-GSW(6))
c    CorEle    = GSW(2)
c    CorFin    = GSW(3)
c    CorEleFin = GSW(4)
         RhoEW     = m_GSW(1)     ! switch to activate electroweak formfactor 
         VPgamma   = m_GSW(6)
         VPgamma = 1.d0   /(2.d0 - m_GSW(6))
         CorEle    = m_GSW(2)
         CorFin    = m_GSW(3)
         CorEleFin = m_GSW(4)
!         sw2 = m_Sw2            ! 0.22351946d0     ! sin(theta_w)^2 from  Eur.Phys.J. C (2019) 79, 480
         sw2 = m_swsq
         alpha=1.d0/m_AlfInv
         Mz=m_MZ
         Gz=m_GammZ*Svar/Mz**2  ! running width
         if (ifPrint.eq.1) then
          write(*,*) ' m_KeyElw= ', m_KeyElw
          write(*,*) 'sw2,alpha,Mz,Gz,svar=',  sw2,alpha,Mz,Gz,svar
          write(*,*) 'G_mu, m_gmu= ', G_mu, m_gmu 
          write(*,*) 'RhoEW     = ',RhoEW,   m_GSW(1)
          write(*,*) 'VPgamma   = ',VPgamma, m_GSW(6)
          write(*,*) 'gsw(6) (7)= ', m_GSW(6), m_GSW(7)          
          write(*,*) 'CorEle    = ', CorEle, m_GSW(2)
          write(*,*) 'CorFin    = ',CorFin , m_GSW(3)
          write(*,*) 'CorEleFin = ',CorEleFin, m_GSW(4)
          write(*,*) 'CorEleFin -CorEle*CorFin = ',CorEleFin-CorEle*CorFin
!         stop
         endif
      else
         RhoEW     = 1.d0
         VPgamma   = 1.d0
         CorEle    = 1.d0
         CorFin    = 1.d0
         CorEleFin = 1.d0
         sw2 = 0.22351946d0     ! sin(theta_w)^2 from  Eur.Phys.J. C (2019) 79, 480
         
      endif
	  
	  V= dsqrt(1.d0 -(mtau/Energy)**2)       ! tau velocity

c     Contributions to ReA(s) and ImA(s) from QED in one loop
      if(iqed.eq.10) then    ! Warning: temporarily this contribution is blocked
       ArQED = -alpha*mtau**2 /(PI*V*4.d0*Energy**2)*dlog((1.d0+V)/(1.d0-V))
       AiQED = alpha *mtau**2 /(V *4.d0 *Energy**2)
      else
       ArQED = 0.d0   !switch for dipole QED part  
       AiQED = 0.d0
      endif
       
	  

      ReA = ArQED + ReA0                   ! QED + new physics
      ImA = AiQED + ImA0  	               ! QED + new physics 

      
	  m= mtau   
	  gam= Energy/mtau    	                 ! Lorentz factor of tau
      e= dsqrt(4.0d0 *PI *alpha)    	     ! electric charge
c       Below coupling f=g_w/(2*cos(theta_w))=e/(2*cos(theta_w)*sin(theta_w)), 
c       which is expressed via Fermi constant:  
      f = Mz *dsqrt(G_mu *dsqrt(2.d0))  
	  

      sw = dsqrt(sw2)        ! sin(theta_w) 
      cw2 = 1.d0 - sw2
      cw = dsqrt(cw2)	
  	  
c       Mandelstam variables s and t (mass of tau is included in t) 	  
      s = 4.d0 *Energy**2                        
      t = m**2 -0.5d0 *s *(1.d0 -V*dcos(theta))  

c       Denominator of the Z propagator with running width
      Den = s - Mz**2 + (0.d0, 1.d0) *Gz *Mz
c       Constant for effective photon like couplings correction coming from Z boson	  
      Fvivf = 4.d0 *(f**2/e**2) *sw**4 ! WARNING: the *sw**4 factor come from vector couplings numerators Fvivf 
c       Note: Fvivf can also be written as Fvivf= sw**2/cw**2, or through the constant G_mu:
c       Fvivf= sqrt(2) *G_mu *Mz**2 *sw**4 /(pi*alpha)  This agrees with definition below.   
 	  
c       FINAL FERMION IS TAU LEPTON (f=lepton) with couplings: 	  
      Qf= -1.d0                 
      vf= -0.03783d0  ! from PDG, used in previous code (v0 coupling)
      af= -0.50123d0  ! from PDG, used in previous code (a0 coupling)
      if(IfGSW.eq.1) then
        vf= -0.5d0 - 2.d0 *Qf *sw2 * K1(s,t) ! is function for lepto
        af= -0.5d0
       endif
          
c	  
c        COUPLINGS FOR INITIAL FERMIONS e u or d             
      IF (channel.EQ.1) THEN    ! LEPTONS       
	  Qi= -1.d0                                        
      	  vi= -0.03783d0  ! from PDG, used in the previous code
	  ai= -0.50123d0  ! from PDG, used in the previous code	  
c	  
c           Effective Z propagator  (reduces to 1/Den without rad. corr.):	  
          Pz = Rho11(s,t) /Den           
          if(IfGSW.eq.1) then
            vi= -0.5d0 - 2.d0 *Qi *sw2 *K1(s,t)       
            ai= -0.5d0
            if (ifPrint.eq.1) write(*,*) 'no EW-corr  Fvivf=',Fvivf,' Fvivf is for (v_if-v_i*v_f) implemented as add-up to photon'
            Fvivf=
     $            m_Gmu *m_MZ**2 *m_AlfInv /(DSQRT(2.d0)*8.d0*m_pi)
     $             *(sw2*(1.d0-sw2)) *16.d0
            Fvivf=Fvivf              *sw2/(1.d0-sw2) ! why this factor? Because  vi ai couplings in KKMC
                 ! are divided by deno= 4 sqrt(m_swsq*(1d0-m_swsq)) in addition multiplied by 2 and we are not using it here
                 ! Ve = (2*T3e -4*Qe*m_swsq)/deno and Ae = 2*T3e/deno. 
            if (ifPrint.eq.1) then
             write(*,*) 'with EW-corr  Fvivf=',Fvivf,' Fvivf is for (v_if-v_i*v_f) implemented as add-up to photon'
             stop
            endif
         endif
         
C          Effective photon propagator with v_{if}-v_i*v_f correction to Z exchange  
C          (NOTE: Pz reduces to 1/s if  rad. corr. are off): 
	 Pg = 1.d0/s *(Gammavp(s) + Fvivf *s/Den*
     $                              Rho11(s,t)*(K11(s,t)-K1(s,t)*K1(s,t)))
          
c	 
       ELSE IF (channel.EQ.2) THEN 	  ! UP QUARK  	   
        Qi= +2.d0/3.d0                                   
        vi= +0.266d0   ! from PDG, used in the previous code
        ai= +0.519d0   ! from PDG, used in the previous code
        if(IfGSW.eq.1) then 
         vi= +0.5d0 - 2.d0 *Qi *sw2 *K2(s,t)     
         ai= +0.5d0
         if (ifPrint.eq.1) write(*,*) 'no   EW-corr  Fvivf=',Fvivf,' Fvivf is for (v_if-v_i*v_f) implemented as add-up to photon'
         Fvivf=
     $   m_Gmu *m_MZ**2 *m_AlfInv /(DSQRT(2.d0)*8.d0*m_pi)
     $         *(sw2*(1.d0-sw2)) *16.d0
         Fvivf=Fvivf              *sw2/(1.d0-sw2) ! why this factor? Because  vi ai couplings in KKMC
            ! are divided by deno= 4 sqrt(m_swsq*(1d0-m_swsq)) in addition multiplied by 2. 
            ! Ve = (2*T3e -4*Qe*m_swsq)/deno and Ae = 2*T3e/deno. Note that Fvivf is used only for calculation of
            ! (vv_if-v_i*v_f)
         if (ifPrint.eq.1) then
           write(*,*) 'with EW-corr  Fvivf=',Fvivf,' Fvivf is for (v_if-v_i*v_f) implemented as add-up to photon'
           stop
          endif
        endif

c
c       Effective Z propagator (reduces to 1/Den without rad. corr.):
      Pz = Rho12(s,t) /Den 	  
c       Effective photon propagator with v_{if}-v_i*v_f correction to Z exchange  
C       (NOTE: Pz reduces to 1/s if  rad. corr. are off):
      Pg = 1.d0/s *(Gammavp(s) + Fvivf *s/Den*
     $                           Rho12(s,t)*(K12(s,t)-K1(s,t)*K2(s,t))) 
c	 
       ELSE IF (channel.EQ.3) THEN    ! DOWN QUARK 
        Qi= -1.d0/3.d0                                      
        vi= -0.38d0       ! from PDG, used in the previous code
        ai= -0.527d0      ! from PDG, used in the previous code
        if(IfGSW.eq.1) then
          vi= -0.5d0 - 2.d0 *Qi *sw2 *K3(s,t)      
          ai= -0.5d0
          if (ifPrint.eq.1) write(*,*) 'no   EW-corr  Fvivf=',Fvivf,' Fvivf is for (v_if-v_i*v_f) implemented as add-up to photon'
          Fvivf=
     $          m_Gmu *m_MZ**2 *m_AlfInv /(DSQRT(2.d0)*8.d0*m_pi)
     $                *(sw2*(1.d0-sw2)) *16.d0
          Fvivf=Fvivf              *sw2/(1.d0-sw2) ! why this factor? Because  vi ai couplings in KKMC
              ! are divided by deno= 4 sqrt(m_swsq*(1d0-m_swsq)) in addition multiplied by 2. 
              ! Ve = (2*T3e -4*Qe*m_swsq)/deno and Ae = 2*T3e/deno. Note that Fvivf is used only for calculation of
              ! (vv_if-v_i*v_f)
          if (ifPrint.eq.1) then
           write(*,*) 'with EW-corr  Fvivf=',Fvivf,' Fvivf is for (v_if-v_i*v_f) implemented as add-up to photon'
           stop
          endif
        endif
        
c		 
c         Effective Z propagator  (reduces to 1/Den without rad. corr.): 
	Pz = Rho13(s,t) /Den 
c         Effective photon propagator  (reduces to 1/s without rad. corr.): 
        Pg = 1.d0/s *(Gammavp(s) + Fvivf *s/Den*
     $                             Rho13(s,t) *(K13(s,t)-K1(s,t)*K3(s,t))) 
c
      ELSE                           ! OUTPUT WILL BE ZERO     
        Qi= 0.d0                                                
        vi= (0.d0, 0.d0) 
        ai= 0.d0 
        Pz= (0.d0, 0.d0)
        Pg= (0.d0, 0.d0) 	  
      END IF 	   

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c         Complex magnetic and electric form-factors 
      A = dcmplx(ReA, ImA)
      B = dcmplx(ReB, ImB)
      X = dcmplx(ReX, ImX)
      Y = dcmplx(ReY, ImY)

c  Contributions to A(s), X(s) from QED in one loop: NOT INCLUDED 
cc      ArQED = -alpha *mtau**2 /(PI *V *4.d0 *E**2) *dlog((1.d0+V)/(1.d0-V))
cc     AiQED = alpha *mtau**2 /(V *4.d0 *E**2)
cc      Ar1 = ArQED + Ar                     ! QED + new physics
cc      Ai1 = AiQED + Ai                     ! QED + new physics


 
c ----------------------------------------------------------------
c    CONTRIBUTIONS to R(i,j) IN STANDARD MODEL NO DIPOLE MOMENTS 
c ----------------------------------------------------------------
c 

      RSM(1,1)= 4.d0*gam**2*m**4*(e**2*(1.d0 + gam**2)*Qf*Qi*
     -     (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dconjg(Pg) + 
     -    f**2*dconjg(Pz)*
     -     (-(af**2*f**2*(-1.d0 + gam**2)*Pz*
     -          (ai**2 + vi*dconjg(vi))) + 
     -       (1.d0 + gam**2)*dconjg(vf)*
     -        (ai**2*f**2*Pz*vf + 
     -          (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dconjg(vi))))*
     -  dsin(theta)**2

      RSM(1,2)=  (0.d0, -4.d0)*af*f**2*gam**4*m**4*V*
     -  (e**2*Pz*Qf*Qi*vi*dconjg(Pg) + 
     -    dconjg(Pz)*(-(ai**2*f**2*Pz*vf) - 
     -       (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dconjg(vi) + 
     -       f**2*Pz*dconjg(vf)*(ai**2 + vi*dconjg(vi))))*
     -  dsin(theta)**2

      RSM(2,1)= RSM(1,2)

      RSM(2,2)=  4.d0*gam**2*(-1.d0 + gam**2)*m**4*
     -  (-(e**2*Qf*Qi*(e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*
     -       dconjg(Pg)) + 
     -    dconjg(Pz)*(af**2*f**4*Pz*
     -        (ai**2 + vi*dconjg(vi)) - 
     -       f**2*dconjg(vf)*
     -        (ai**2*f**2*Pz*vf + 
     -          (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dconjg(vi))))*
     -  dsin(theta)**2

      RSM(1,3)=   4.d0*gam**3*m**4*(e**2*Qf*Qi*dconjg(Pg)*
     -     (af*ai*f**2*Pz*V - 
     -       2.d0*gam**2*(-1.d0 + V**2)*
     -        (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dcos(theta)) + 
     -    f**2*dconjg(Pz)*
     -     (af*(ai*V*(e**2*Pg*Qf*Qi + 
     -             f**2*Pz*vf*(vi + dconjg(vi))) + 
     -          2.d0*af*f**2*Pz*(1.d0 + gam**2*(-1.d0 + V**2))*
     -           (ai**2 + vi*dconjg(vi))*dcos(theta)) + 
     -       dconjg(vf)*
     -        (ai*f**2*Pz*(af*V*vi - 
     -             2.d0*ai*gam**2*(-1.d0 + V**2)*vf*dcos(theta)) + 
     -          dconjg(vi)*
     -           (af*ai*f**2*Pz*V - 
     -             2.d0*gam**2*(-1.d0 + V**2)*
     -              (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dcos(theta))))
     -    )*dsin(theta)

      RSM(3,1)= RSM(1,3)

      RSM(2,3)=  (0.d0, -2.d0)*af*f**2*gam**3*m**4*V*
     -  (e**2*Pz*Qf*Qi*vi*dconjg(Pg) + 
     -    dconjg(Pz)*(-(ai**2*f**2*Pz*vf) - 
     -       (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dconjg(vi) + 
     -       f**2*Pz*dconjg(vf)*(ai**2 + vi*dconjg(vi))))*
     -  dsin(2.d0*theta)

      RSM(3,2)= RSM(2,3)
     
      RSM(3,3)=  2.d0*gam**2*m**4*(-(e**4*Pg*Qf**2*Qi**2*dconjg(Pg)*
     -       (1.d0 - 3.d0*gam**2 + 
     -         (-1.d0 + 3.d0*gam**2 + 4.d0*gam**4*(-1.d0 + V**2))*
     -          dcos(2.d0*theta))) - 
     -    e**2*f**2*Qf*Qi*(Pz*dconjg(Pg)*
     -        (-4.d0*af*ai*gam**2*V*dcos(theta) + 
     -          vf*vi*(1.d0 - 3.d0*gam**2 + 
     -          (-1.d0 + 3.d0*gam**2 + 4.d0*gam**4*(-1.d0 + V**2))*
     -              dcos(2.d0*theta))) + 
     -       Pg*dconjg(Pz)*
     -        (-4.d0*af*ai*gam**2*V*dcos(theta) + 
     -          dconjg(vf)*dconjg(vi)*
     -           (1.d0 - 3.d0*gam**2 + 
     -          (-1.d0 + 3.d0*gam**2 + 4.d0*gam**4*(-1.d0 + V**2))*
     -              dcos(2.d0*theta)))) + 
     -    f**4*Pz*dconjg(Pz)*
     -     (af*(4*ai*gam**2*V*vf*(vi + dconjg(vi))*
     -           dcos(theta) + 
     -          af*(ai**2 + vi*dconjg(vi))*
     -           (1.d0 + gam**2*(-1.d0 + 4.d0*V**2) + 
     -          (-1.d0 + 5.d0*gam**2 + 4.d0*gam**4*(-1.d0 + V**2))*
     -              dcos(2.d0*theta))) + 
     -       dconjg(vf)*
     -        (dconjg(vi)*
     -           (4.d0*af*ai*gam**2*V*dcos(theta) + 
     -             vf*vi*(-1.d0 + 3.d0*gam**2 + 
     -            (1.d0 - 3.d0*gam**2 - 4.d0*gam**4*(-1.d0 + V**2))*
     -                 dcos(2.d0*theta))) + 
     -          ai*(4.d0*af*gam**2*V*vi*dcos(theta) - 
     -             ai*vf*(1.d0 - 3.d0*gam**2 + 
     -            (-1.d0 + 3.d0*gam**2 + 4.d0*gam**4*(-1.d0 + V**2))*
     -                 dcos(2.d0*theta))))))    

      RSM(1,4)=  -4.d0*f**2*gam**3*m**4*
     -  (e**2*Pz*Qf*Qi*dconjg(Pg)*
     -     (2.d0*ai*vf + af*V*vi*dcos(theta)) + 
     -    dconjg(Pz)*(af*V*
     -        (ai**2*f**2*Pz*vf + 
     -          (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dconjg(vi))*
     -        dcos(theta) + 
     -       dconjg(vf)*
     -        (2.d0*ai*(e**2*Pg*Qf*Qi + 
     -             f**2*Pz*vf*(vi + dconjg(vi))) + 
     -          af*f**2*Pz*V*(ai**2 + vi*dconjg(vi))*
     -           dcos(theta))))*dsin(theta)

      RSM(4,1)= RSM(1,4)        

      RSM(2,4)=  (0.d0, 4.d0)*af*ai*f**2*gam**3*m**4*V*
     -  (e**2*Pz*Qf*Qi*dconjg(Pg) + 
     -    dconjg(Pz)*(-(e**2*Pg*Qf*Qi) - 
     -       f**2*Pz*vf*(vi + dconjg(vi)) + 
     -       f**2*Pz*dconjg(vf)*(vi + dconjg(vi))))*
     -  dsin(theta)

      RSM(4,2)= RSM(2,4)

      RSM(3,4)=  -4.d0*f**2*gam**2*m**4*
     -  (e**2*gam**2*Pz*Qf*Qi*dconjg(Pg)*
     -   (2.d0*ai*vf*dcos(theta) + af*V*vi*(1.d0 + dcos(theta)**2)) + 
     -    (dconjg(Pz)*(dconjg(vi)*
     -          (4.d0*af**2*ai*f**2*(-1.d0 + gam**2)*Pz*dcos(theta) + 
     -            4.d0*ai*f**2*gam**2*Pz*vf*dconjg(vf)*
     -             dcos(theta) + 
     -            af*gam**2*V*
     -             (e**2*Pg*Qf*Qi + 
     -               f**2*Pz*vi*(vf + dconjg(vf)))*
     -             (3.d0 + dcos(2.d0*theta))) + 
     -         ai*(gam**2*dconjg(vf)*
     -             (4.d0*(e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*
     -                dcos(theta) + 
     -               af*ai*f**2*Pz*V*(3.d0 + dcos(2.d0*theta))) + 
     -            af*f**2*Pz*
     -             (4.d0*af*(-1.d0 + gam**2)*vi*dcos(theta) + 
     -             ai*gam**2*V*vf*(3.d0 + dcos(2.d0*theta))))))/2.d0)

      RSM(4,3)= RSM(3,4) 


      RSM(4,4)=  2.d0*gam**2*m**4*(e**4*Pg*Qf**2*Qi**2*dconjg(Pg)*
     -     (1.d0 + 3.d0*gam**2 + (-1.d0 + gam**2)*dcos(2.d0*theta)) + 
     -    e**2*f**2*Qf*Qi*(Pz*dconjg(Pg)*
     -        (4.d0*af*ai*gam**2*V*dcos(theta) + 
     -  vf*vi*(1.d0 + 3.d0*gam**2 + (-1.d0 + gam**2)*dcos(2.d0*theta))
     -          ) + Pg*dconjg(Pz)*
     -        (4.d0*af*ai*gam**2*V*dcos(theta) + 
     -          dconjg(vf)*dconjg(vi)*
     -    (1.d0 + 3.d0*gam**2 + (-1.d0 + gam**2)*dcos(2.d0*theta)))) - 
     -    f**4*Pz*dconjg(Pz)*
     -     (-(af*(4.d0*ai*gam**2*V*vf*(vi + dconjg(vi))*
     -             dcos(theta) + 
     -            af*(-1.d0 + gam**2)*(ai**2 + vi*dconjg(vi))*
     -             (3.d0 + dcos(2.d0*theta)))) - 
     -       dconjg(vf)*
     -        (ai*(4.d0*af*gam**2*V*vi*dcos(theta) + 
     -             ai*vf*(1.d0 + 3.d0*gam**2 + 
     -                (-1.d0 + gam**2)*dcos(2.d0*theta))) + 
     -          dconjg(vi)*
     -           (4.d0*af*ai*gam**2*V*dcos(theta) + 
     -             vf*vi*(1.d0 + 3.d0*gam**2 + 
     -                (-1.d0 + gam**2)*dcos(2.d0*theta))))))


c -------------------------------------------------------------
c CONTRIBUTION TO R(i,j), LINEAR IN DIPOLE MOMENTS A, B, X, Y 
c -------------------------------------------------------------

      RDM(1,1)= 8.d0*gam**4*m**4*(e**2*Qf*Qi*
     -     (e**2*Pg*Qf*Qi*(A + dconjg(A)) + 
     -       f**2*Pz*vi*(X + vf*dconjg(A)))*dconjg(Pg) + 
     -    f**2*dconjg(Pz)*
     -     (dconjg(vf)*(ai**2*f**2*Pz*X + 
     -        (A*e**2*Pg*Qf*Qi + f**2*Pz*vi*X)*dconjg(vi)) +
     -         (ai**2*f**2*Pz*vf + 
     -          (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dconjg(vi))*
     -        dconjg(X)))*dsin(theta)**2

      RDM(1,2)=  4.d0*gam**4*m**4*V*(e**2*Qf*Qi*
     -     (e**2*Pg*Qf*Qi*(B + dconjg(B)) + 
     -       f**2*Pz*vi*(Y - (0.d0, 1.d0)*af*dconjg(A) + 
     -          vf*dconjg(B)))*dconjg(Pg) + 
     -    f**2*dconjg(Pz)*
     -     (dconjg(vf)*(ai**2*f**2*Pz*Y + 
     -          (B*e**2*Pg*Qf*Qi + f**2*Pz*vi*Y)*dconjg(vi)) +
     -         (0.d0, 1.d0)*af*(dconjg(vi)*
     -           (A*e**2*Pg*Qf*Qi + 
     -             f**2*Pz*vi*(X - dconjg(X))) + 
     -          ai**2*f**2*Pz*(X - dconjg(X))) + 
     -       (ai**2*f**2*Pz*vf + 
     -          (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dconjg(vi))*
     -        dconjg(Y)))*dsin(theta)**2

      RDM(2,1)=  -4.d0*gam**4*m**4*V*(e**2*Qf*Qi*
     -     (e**2*Pg*Qf*Qi*(B + dconjg(B)) + 
     -       f**2*Pz*vi*(Y + (0,1)*af*dconjg(A) + 
     -          vf*dconjg(B)))*dconjg(Pg) + 
     -    f**2*dconjg(Pz)*
     -     (dconjg(vf)*(ai**2*f**2*Pz*Y + 
     -          (B*e**2*Pg*Qf*Qi + f**2*Pz*vi*Y)*dconjg(vi)) -
     -         (0.d0 ,1.d0)*af*(dconjg(vi)*
     -           (A*e**2*Pg*Qf*Qi + 
     -             f**2*Pz*vi*(X - dconjg(X))) + 
     -          ai**2*f**2*Pz*(X - dconjg(X))) + 
     -       (ai**2*f**2*Pz*vf + 
     -          (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dconjg(vi))*
     -        dconjg(Y)))*dsin(theta)**2

      RDM(2,2)= 0.d0

      RDM(1,3)=  4.d0*gam**5*m**4*(-(e**2*Qf*Qi*dconjg(Pg)*
     -       ((0.d0, -1.d0)*ai*f**2*Pz*V*
     -          (Y - (0.d0, 1.d0)*af*dconjg(A) - vf*dconjg(B)) + 
     -         (e**2*Pg*Qf*Qi*(-2 + V**2)*(A + dconjg(A)) + 
     -            f**2*Pz*vi*
     -             ((-2.d0 + V**2)*(X + vf*dconjg(A)) + 
     -            (0.d0, 1.d0)*af*V**2*dconjg(B)))*dcos(theta))) + 
     -    f**2*dconjg(Pz)*
     -     (dconjg(vf)*(ai*
     -           ((0.d0, 1.d0)*V*(B*e**2*Pg*Qf*Qi + f**2*Pz*vi*Y) - 
     -             ai*f**2*Pz*(-2.d0 + V**2)*X*dcos(theta)) + 
     -          dconjg(vi)*
     -           ((0.d0, 1.d0)*ai*f**2*Pz*V*Y - 
     -             (-2.d0 + V**2)*(A*e**2*Pg*Qf*Qi + f**2*Pz*vi*X)*
     -              dcos(theta))) + 
     -       dconjg(vi)*
     -        ((0.d0,-1.d0)*ai*f**2*Pz*V*vf*dconjg(Y) - 
     -          (-2.d0 + V**2)*(e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*
     -           dconjg(X)*dcos(theta) + 
     -          af*V*(ai*f**2*Pz*(X + dconjg(X)) + 
     -             (0.d0, 1.d0)*V*
     -              (B*e**2*Pg*Qf*Qi + 
     -                f**2*Pz*vi*(Y - dconjg(Y)))*dcos(theta)))
     -       + ai*((0.d0, -1.d0)*V*(e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*
     -           dconjg(Y) + 
     -          ai*f**2*Pz*(2.d0 - V**2)*vf*dconjg(X)*
     -           dcos(theta) + 
     -          af*V*(A*e**2*Pg*Qf*Qi + 
     -             f**2*Pz*vi*(X + dconjg(X)) + 
     -             (0.d0, 1.d0)*ai*f**2*Pz*V*(Y - dconjg(Y))*
     -              dcos(theta)))))*dsin(theta)

      RDM(3,1)=  4.d0*gam**5*m**4*(e**2*Qf*Qi*dconjg(Pg)*
     -     (ai*f**2*Pz*V*(af*dconjg(A) - 
     -          (0.d0, 1.d0)*(Y - vf*dconjg(B))) - 
     -       (e**2*Pg*Qf*Qi*(-2.d0 + V**2)*(A + dconjg(A)) + 
     -          f**2*Pz*vi*
     -           ((-2.d0 + V**2)*(X + vf*dconjg(A)) - 
     -             (0.d0, 1.d0)*af*V**2*dconjg(B)))*dcos(theta)) + 
     -    f**2*dconjg(Pz)*
     -     (dconjg(vf)*(-(ai*
     -          ((0.d0, 1.d0)*V*(B*e**2*Pg*Qf*Qi + f**2*Pz*vi*Y) + 
     -               ai*f**2*Pz*(-2.d0 + V**2)*X*dcos(theta))) + 
     -          dconjg(vi)*
     -           ((0.d0, -1.d0)*ai*f**2*Pz*V*Y - 
     -             (-2.d0 + V**2)*(A*e**2*Pg*Qf*Qi + f**2*Pz*vi*X)*
     -              dcos(theta))) + 
     -       dconjg(vi)*
     -        ((0.d0, 1.d0)*ai*f**2*Pz*V*vf*dconjg(Y) - 
     -          (-2.d0 + V**2)*(e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*
     -           dconjg(X)*dcos(theta) + 
     -          af*V*(ai*f**2*Pz*(X + dconjg(X)) - 
     -             (0.d0, 1.d0)*V*
     -              (B*e**2*Pg*Qf*Qi + 
     -                f**2*Pz*vi*(Y - dconjg(Y)))*dcos(theta)))
     -         + ai*((0.d0,1.d0)*V*(e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*
     -           dconjg(Y) + 
     -          ai*f**2*Pz*(2.d0 - V**2)*vf*dconjg(X)*
     -           dcos(theta) + 
     -          af*V*(A*e**2*Pg*Qf*Qi + 
     -             f**2*Pz*vi*(X + dconjg(X)) - 
     -             (0.d0, 1.d0)*ai*f**2*Pz*V*(Y - dconjg(Y))*
     -              dcos(theta)))))*dsin(theta)

      RDM(2,3)=  (0.d0, 4.d0)*gam**5*m**4*V*
     -  (e**2*Qf*Qi*dconjg(Pg)*
     -     (ai*f**2*Pz*V*(X - vf*dconjg(A) + 
     -          (0.d0, 1.d0)*af*dconjg(B)) + 
     -       (0.d0, 1.d0)*(e**2*Pg*Qf*Qi*(B + dconjg(B)) + 
     -          f**2*Pz*vi*
     -           (Y + (0.d0, 1.d0)*af*dconjg(A) + vf*dconjg(B)))*
     -        dcos(theta)) + 
     -    f**2*dconjg(Pz)*
     -     (dconjg(vf)*(ai*
     -           (A*e**2*Pg*Qf*Qi*V + f**2*Pz*V*vi*X + 
     -             (0.d0, 1.d0)*ai*f**2*Pz*Y*dcos(theta)) + 
     -          dconjg(vi)*
     -           (ai*f**2*Pz*V*X + 
     -             (0.d0, 1.d0)*(B*e**2*Pg*Qf*Qi + f**2*Pz*vi*Y)*
     -              dcos(theta))) + 
     -       (0.d0, 1.d0)*(dconjg(vi)*
     -           ((0.d0, 1.d0)*ai*f**2*Pz*V*vf*dconjg(X) + 
     -             (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dconjg(Y)*
     -              dcos(theta) + 
     -             af*(ai*f**2*Pz*V*(Y + dconjg(Y)) - 
     -                (0.d0, 1.d0)*
     -                 (A*e**2*Pg*Qf*Qi + 
     -                   f**2*Pz*vi*(X - dconjg(X)))*
     -                 dcos(theta))) + 
     -          ai*((0.d0, 1.d0)*V*(e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*
     -              dconjg(X) + 
     -             ai*f**2*Pz*vf*dconjg(Y)*dcos(theta) + 
     -             af*(V*(B*e**2*Pg*Qf*Qi + 
     -                   f**2*Pz*vi*(Y + dconjg(Y))) - 
     -                (0.d0, 1.d0)*ai*f**2*Pz*(X - dconjg(X))*
     -                 dcos(theta))))))*dsin(theta)

      RDM(3,2)=  4.d0*gam**5*m**4*V*(e**2*Qf*Qi*dconjg(Pg)*
     -     (ai*f**2*Pz*V*((0.d0, 1.d0)*(X - vf*dconjg(A)) + 
     -          af*dconjg(B)) + 
     -       (e**2*Pg*Qf*Qi*(B + dconjg(B)) + 
     -          f**2*Pz*vi*
     -           (Y - (0.d0, 1.d0)*af*dconjg(A) + vf*dconjg(B)))*
     -        dcos(theta)) + 
     -    f**2*dconjg(Pz)*
     -     (dconjg(vf)*(ai*
     -         ((0.d0, 1.d0)*V*(A*e**2*Pg*Qf*Qi + f**2*Pz*vi*X) + 
     -             ai*f**2*Pz*Y*dcos(theta)) + 
     -          dconjg(vi)*
     -           ((0.d0, 1.d0)*ai*f**2*Pz*V*X + 
     -           (B*e**2*Pg*Qf*Qi + f**2*Pz*vi*Y)*dcos(theta))) +
     -         dconjg(vi)*
     -        ((0.d0, -1.d0)*ai*f**2*Pz*V*vf*dconjg(X) + 
     -          (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dconjg(Y)*
     -           dcos(theta) + 
     -          af*(ai*f**2*Pz*V*(Y + dconjg(Y)) + 
     -             (0.d0, 1.d0)*(A*e**2*Pg*Qf*Qi + 
     -             f**2*Pz*vi*(X - dconjg(X)))*dcos(theta)))
     -       + ai*((0.d0, -1.d0)*V*(e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*
     -           dconjg(X) + 
     -          ai*f**2*Pz*vf*dconjg(Y)*dcos(theta) + 
     -          af*(V*(B*e**2*Pg*Qf*Qi + 
     -                f**2*Pz*vi*(Y + dconjg(Y))) + 
     -          (0.d0, 1.d0)*ai*f**2*Pz*(X - dconjg(X))*dcos(theta)
     -             ))))*dsin(theta)

      RDM(3,3)=  -8.d0*gam**6*m**4*(-1.d0 + V**2)*dcos(theta)*
     -  (e**2*Qf*Qi*dconjg(Pg)*
     -     ((A*e**2*Pg*Qf*Qi + f**2*Pz*vi*X)*dcos(theta) + 
     -       dconjg(A)*(af*ai*f**2*Pz*V + 
     -          (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dcos(theta))) + 
     -    f**2*dconjg(Pz)*
     -     (ai*(af*V*(A*e**2*Pg*Qf*Qi + 
     -             f**2*Pz*vi*(X + dconjg(X))) + 
     -          ai*f**2*Pz*X*dconjg(vf)*dcos(theta) + 
     -          ai*f**2*Pz*vf*dconjg(X)*dcos(theta)) + 
     -       dconjg(vi)*
     -        (af*ai*f**2*Pz*V*(X + dconjg(X)) + 
     -          ((A*e**2*Pg*Qf*Qi + f**2*Pz*vi*X)*
     -              dconjg(vf) + 
     -             (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dconjg(X))*
     -           dcos(theta))))

      RDM(1,4)=  -4.d0*gam**3*m**4*(e**2*Qf*Qi*dconjg(Pg)*
     -     (ai*f**2*Pz*((1 + gam**2)*(X + vf*dconjg(A)) - 
     -          (0.d0, 1.d0)*af*(-1.d0 + gam**2)*dconjg(B)) + 
     -       (0.d0, 1.d0)*gam**2*V*
     -        (e**2*Pg*Qf*Qi*(B - dconjg(B)) + 
     -          f**2*Pz*vi*
     -         (Y - (0.d0, 1.d0)*af*dconjg(A) - vf*dconjg(B)))*
     -        dcos(theta)) + 
     -    f**2*dconjg(Pz)*
     -     (dconjg(vf)*(ai*
     -           ((1.d0 + gam**2)*
     -              (A*e**2*Pg*Qf*Qi + f**2*Pz*vi*X) + 
     -           (0.d0, 1.d0)*ai*f**2*gam**2*Pz*V*Y*dcos(theta)) + 
     -          dconjg(vi)*
     -           (ai*f**2*(1.d0 + gam**2)*Pz*X + 
     -             (0.d0, 1.d0)*gam**2*V*
     -              (B*e**2*Pg*Qf*Qi + f**2*Pz*vi*Y)*dcos(theta)))
     -         + ai*((1.d0 + gam**2)*
     -           (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dconjg(X) - 
     -          (0.d0, 1.d0)*ai*f**2*gam**2*Pz*V*vf*dconjg(Y)*
     -           dcos(theta) + 
     -          af*((0.d0, 1.d0)*(-1.d0 + gam**2)*
     -              (B*e**2*Pg*Qf*Qi + 
     -                f**2*Pz*vi*(Y - dconjg(Y))) + 
     -             ai*f**2*gam**2*Pz*V*(X + dconjg(X))*
     -              dcos(theta))) + 
     -       dconjg(vi)*
     -        (ai*f**2*(1.d0 + gam**2)*Pz*vf*dconjg(X) - 
     -        (0.d0, 1.d0)*gam**2*V*(e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*
     -           dconjg(Y)*dcos(theta) + 
     -          af*((0.d0, 1.d0)*ai*f**2*(-1.d0 + gam**2)*Pz*
     -              (Y - dconjg(Y)) + 
     -             gam**2*V*
     -              (A*e**2*Pg*Qf*Qi + 
     -                f**2*Pz*vi*(X + dconjg(X)))*dcos(theta)))
     -       ))*dsin(theta)

      RDM(4,1)=  -4.d0*gam**3*m**4*(e**2*Qf*Qi*dconjg(Pg)*
     -     (ai*f**2*Pz*((1.d0 + gam**2)*(X + vf*dconjg(A)) + 
     -          (0.d0, 1.d0)*af*(-1.d0 + gam**2)*dconjg(B)) - 
     -       (0.d0, 1.d0)*gam**2*V*
     -        (e**2*Pg*Qf*Qi*(B - dconjg(B)) + 
     -          f**2*Pz*vi*
     -           (Y + (0.d0, 1.d0)*af*dconjg(A) - vf*dconjg(B)))*
     -        dcos(theta)) + 
     -    f**2*dconjg(Pz)*
     -     (dconjg(vf)*(ai*
     -           ((1.d0 + gam**2)*
     -              (A*e**2*Pg*Qf*Qi + f**2*Pz*vi*X) - 
     -          (0.d0, 1.d0)*ai*f**2*gam**2*Pz*V*Y*dcos(theta)) + 
     -          dconjg(vi)*
     -           (ai*f**2*(1.d0 + gam**2)*Pz*X - 
     -             (0.d0, 1.d0)*gam**2*V*
     -              (B*e**2*Pg*Qf*Qi + f**2*Pz*vi*Y)*dcos(theta)))
     -         + ai*((1.d0 + gam**2)*
     -           (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dconjg(X) + 
     -          (0.d0, 1.d0)*ai*f**2*gam**2*Pz*V*vf*dconjg(Y)*
     -           dcos(theta) + 
     -          af*((0.d0, -1.d0)*(-1.d0 + gam**2)*
     -              (B*e**2*Pg*Qf*Qi + 
     -                f**2*Pz*vi*(Y - dconjg(Y))) + 
     -             ai*f**2*gam**2*Pz*V*(X + dconjg(X))*
     -              dcos(theta))) + 
     -       dconjg(vi)*
     -        (ai*f**2*(1.d0 + gam**2)*Pz*vf*dconjg(X) + 
     -        (0.d0, 1.d0)*gam**2*V*(e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*
     -           dconjg(Y)*dcos(theta) + 
     -          af*((0.d0, -1.d0)*ai*f**2*(-1.d0 + gam**2)*Pz*
     -              (Y - dconjg(Y)) + 
     -             gam**2*V*
     -              (A*e**2*Pg*Qf*Qi + 
     -                f**2*Pz*vi*(X + dconjg(X)))*dcos(theta)))
     -       ))*dsin(theta)

      RDM(2,4)=  4.d0*gam**3*m**4*(e**2*Qf*Qi*dconjg(Pg)*
     -     (ai*f**2*gam**2*Pz*V*
     -        (Y + (0.d0, 1.d0)*af*dconjg(A) + vf*dconjg(B)) - 
     -       (0.d0, 1.d0)*(-1.d0 + gam**2)*
     -        (e**2*Pg*Qf*Qi*(A - dconjg(A)) + 
     -          f**2*Pz*vi*
     -           (X - vf*dconjg(A) + (0.d0, 1.d0)*af*dconjg(B)))*
     -        dcos(theta)) + 
     -    f**2*dconjg(Pz)*
     -     (dconjg(vf)*(ai*
     -           (gam**2*V*(B*e**2*Pg*Qf*Qi + f**2*Pz*vi*Y) - 
     -        (0.d0, 1.d0)*ai*f**2*(-1.d0 + gam**2)*Pz*X*dcos(theta)) +
     -            dconjg(vi)*
     -           (ai*f**2*gam**2*Pz*V*Y - 
     -             (0.d0, 1.d0)*(-1.d0 + gam**2)*
     -              (A*e**2*Pg*Qf*Qi + f**2*Pz*vi*X)*dcos(theta)))
     -         + ai*(gam**2*V*(e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*
     -           dconjg(Y) + 
     -          (0.d0, 1.d0)*ai*f**2*(-1.d0 + gam**2)*Pz*vf*dconjg(X)*
     -           dcos(theta) + 
     -          af*((0.d0, -1.d0)*gam**2*V*
     -              (A*e**2*Pg*Qf*Qi + 
     -                f**2*Pz*vi*(X - dconjg(X))) + 
     -             ai*f**2*(-1.d0 + gam**2)*Pz*(Y + dconjg(Y))*
     -              dcos(theta))) + 
     -       dconjg(vi)*
     -        ((0.d0, 1.d0)*((0.d0, -1.d0)*ai*f**2*gam**2*Pz*V*vf*
     -              dconjg(Y) + 
     -             (-1.d0 + gam**2)*(e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*
     -              dconjg(X)*dcos(theta)) + 
     -          af*((0.d0, -1.d0)*ai*f**2*gam**2*Pz*V*
     -              (X - dconjg(X)) + 
     -             (-1.d0 + gam**2)*
     -              (B*e**2*Pg*Qf*Qi + 
     -                f**2*Pz*vi*(Y + dconjg(Y)))*dcos(theta)))
     -       ))*dsin(theta)

      RDM(4,2)=   4.d0*gam**3*m**4*(e**2*Qf*Qi*dconjg(Pg)*
     -     (-(ai*f**2*gam**2*Pz*V*
     -          (Y - (0.d0, 1.d0)*af*dconjg(A) + vf*dconjg(B))) - 
     -       (0.d0, 1.d0)*(-1.d0 + gam**2)*
     -        (e**2*Pg*Qf*Qi*(A - dconjg(A)) + 
     -          f**2*Pz*vi*
     -           (X - vf*dconjg(A) - (0.d0, 1.d0)*af*dconjg(B)))*
     -        dcos(theta)) - 
     -    f**2*dconjg(Pz)*
     -     (dconjg(vf)*(ai*
     -           (gam**2*V*(B*e**2*Pg*Qf*Qi + f**2*Pz*vi*Y) + 
     -      (0.d0, 1.d0)*ai*f**2*(-1.d0 + gam**2)*Pz*X*dcos(theta)) +
     -            dconjg(vi)*
     -           (ai*f**2*gam**2*Pz*V*Y + 
     -             (0.d0, 1.d0)*(-1.d0 + gam**2)*
     -              (A*e**2*Pg*Qf*Qi + f**2*Pz*vi*X)*dcos(theta)))
     -         + ai*(gam**2*V*(e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*
     -           dconjg(Y) - 
     -          (0.d0, 1.d0)*ai*f**2*(-1.d0 + gam**2)*Pz*vf*dconjg(X)*
     -           dcos(theta) + 
     -          af*((0.d0, 1.d0)*gam**2*V*
     -              (A*e**2*Pg*Qf*Qi + 
     -                f**2*Pz*vi*(X - dconjg(X))) + 
     -             ai*f**2*(-1.d0 + gam**2)*Pz*(Y + dconjg(Y))*
     -              dcos(theta))) + 
     -       dconjg(vi)*
     -        ((0.d0, -1.d0)*((0.d0, 1.d0)*ai*f**2*gam**2*Pz*V*vf*
     -              dconjg(Y) + 
     -             (-1.d0 + gam**2)*(e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*
     -              dconjg(X)*dcos(theta)) + 
     -          af*((0.d0, 1.d0)*ai*f**2*gam**2*Pz*V*
     -              (X - dconjg(X)) + 
     -             (-1.d0 + gam**2)*
     -              (B*e**2*Pg*Qf*Qi + 
     -                f**2*Pz*vi*(Y + dconjg(Y)))*dcos(theta)))
     -       ))*dsin(theta)

      RDM(3,4)=   -2.d0*gam**4*m**4*(f**2*dconjg(Pz)*
     -     (-(af*V*(ai**2*f**2*Pz*(X + dconjg(X)) + 
     -            dconjg(vi)*
     -             (A*e**2*Pg*Qf*Qi + 
     -               f**2*Pz*vi*(X + dconjg(X))))*
     -          (-2.d0 + gam**2*(-1.d0 + V**2) + 
     -            gam**2*(-1.d0 + V**2)*dcos(2.d0*theta))) + 
     -       (0.d0, 1.d0)*((0.d0, -4.d0)*ai*
     -           (e**2*Pg*Qf*Qi + 
     -             f**2*Pz*vf*(vi + dconjg(vi)))*dconjg(X)*
     -           dcos(theta) + 
     -          V*(ai**2*f**2*Pz*vf + 
     -             (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dconjg(vi))
     -            *dconjg(Y)*
     -           (2.d0 + gam**2*(-1.d0 + V**2) + 
     -             gam**2*(-1.d0 + V**2)*dcos(2.d0*theta))) + 
     -       dconjg(vf)*
     - ((0.d0, -1.d0)*ai*((0.d0, 4.d0)*(A*e**2*Pg*Qf*Qi + f**2*Pz*vi*X)*
     -              dcos(theta) + 
     -             ai*f**2*Pz*V*Y*
     -              (2.d0 + gam**2*(-1.d0 + V**2) + 
     -                gam**2*(-1.d0 + V**2)*dcos(2.d0*theta))) + 
     -          dconjg(vi)*
     -           (4.d0*ai*f**2*Pz*X*dcos(theta) - 
     -             (0.d0, 1.d0)*V*(B*e**2*Pg*Qf*Qi + f**2*Pz*vi*Y)*
     -              (2.d0 + gam**2*(-1.d0 + V**2) + 
     -                gam**2*(-1.d0 + V**2)*dcos(2.d0*theta))))) + 
     -    e**2*Qf*Qi*dconjg(Pg)*
     -     (4.d0*ai*f**2*Pz*(X + vf*dconjg(A))*dcos(theta) - 
     -       (0.d0, 1.d0)*V*(e**2*Pg*Qf*Qi*(B - dconjg(B))*
     -           (2.d0 + gam**2*(-1.d0 + V**2) + 
     -             gam**2*(-1.d0 + V**2)*dcos(2.d0*theta)) + 
     -          f**2*Pz*vi*
     -           ((0.d0, -1.d0)*af*dconjg(A)*
     -              (-2.d0 + gam**2*(-1.d0 + V**2) + 
     -                gam**2*(-1.d0 + V**2)*dcos(2.d0*theta)) + 
     -             (Y - vf*dconjg(B))*
     -              (2.d0 + gam**2*(-1.d0 + V**2) + 
     -                gam**2*(-1.d0 + V**2)*dcos(2.d0*theta))))))

      RDM(4,3)=  -2.d0*gam**4*m**4*(f**2*dconjg(Pz)*
     -     (-(af*V*(ai**2*f**2*Pz*(X + dconjg(X)) + 
     -            dconjg(vi)*
     -             (A*e**2*Pg*Qf*Qi + 
     -               f**2*Pz*vi*(X + dconjg(X))))*
     -          (-2.d0 + gam**2*(-1.d0 + V**2) + 
     -            gam**2*(-1.d0 + V**2)*dcos(2.d0*theta))) - 
     -       (0.d0, 1.d0)*((0.d0, 4.d0)*ai*
     -           (e**2*Pg*Qf*Qi + 
     -             f**2*Pz*vf*(vi + dconjg(vi)))*dconjg(X)*
     -           dcos(theta) + 
     -          V*(ai**2*f**2*Pz*vf + 
     -             (e**2*Pg*Qf*Qi + f**2*Pz*vf*vi)*dconjg(vi))
     -            *dconjg(Y)*
     -           (2.d0 + gam**2*(-1.d0 + V**2) + 
     -             gam**2*(-1.d0 + V**2)*dcos(2.d0*theta))) + 
     -       dconjg(vf)*
     - ((0.d0, 1.d0)*ai*((0.d0, -4.d0)*(A*e**2*Pg*Qf*Qi + f**2*Pz*vi*X)*
     -              dcos(theta) + 
     -             ai*f**2*Pz*V*Y*
     -              (2.d0 + gam**2*(-1.d0 + V**2) + 
     -                gam**2*(-1.d0 + V**2)*dcos(2.d0*theta))) + 
     -          dconjg(vi)*
     -           (4.d0*ai*f**2*Pz*X*dcos(theta) + 
     -             (0.d0, 1.d0)*V*(B*e**2*Pg*Qf*Qi + f**2*Pz*vi*Y)*
     -              (2.d0 + gam**2*(-1.d0 + V**2) + 
     -                gam**2*(-1.d0 + V**2)*dcos(2.d0*theta))))) + 
     -    e**2*Qf*Qi*dconjg(Pg)*
     -     (4.d0*ai*f**2*Pz*(X + vf*dconjg(A))*dcos(theta) + 
     -       (0.d0, 1.d0)*V*(e**2*Pg*Qf*Qi*(B - dconjg(B))*
     -           (2.d0 + gam**2*(-1.d0 + V**2) + 
     -             gam**2*(-1.d0 + V**2)*dcos(2.d0*theta)) + 
     -          f**2*Pz*vi*
     -           ((0.d0, 1.d0)*af*dconjg(A)*
     -              (-2.d0 + gam**2*(-1.d0 + V**2) + 
     -                gam**2*(-1.d0 + V**2)*dcos(2.d0*theta)) + 
     -             (Y - vf*dconjg(B))*
     -              (2.d0 + gam**2*(-1.d0 + V**2) + 
     -                gam**2*(-1.d0 + V**2)*dcos(2.d0*theta))))))

      RDM(4,4)=  8.d0*gam**4*m**4*(e**4*Pg*Qf**2*Qi**2*(A + dconjg(A))*
     -     dconjg(Pg) + 
     -    e**2*f**2*Qf*Qi*(Pz*vi*X*dconjg(Pg) + 
     -       A*Pg*dconjg(Pz)*dconjg(vf)*dconjg(vi) + 
     -       Pg*dconjg(Pz)*dconjg(vi)*dconjg(X) + 
     -       A*af*ai*Pg*V*dconjg(Pz)*dcos(theta) + 
     -       Pz*dconjg(A)*dconjg(Pg)*
     -        (vf*vi + af*ai*V*dcos(theta))) + 
     -    f**4*Pz*dconjg(Pz)*
     -     (ai*(ai*X*dconjg(vf) + ai*vf*dconjg(X) + 
     -          af*V*vi*X*dcos(theta) + 
     -          af*V*vi*dconjg(X)*dcos(theta)) + 
     -       dconjg(vi)*
     -        (vi*(X*dconjg(vf) + vf*dconjg(X)) + 
     -          af*ai*V*(X + dconjg(X))*dcos(theta))))

c-------------------------------------------------------------
c   TOTAL SPIN-CORRELATION MATRIX: STANDARD MODEL + DIPOLE MOMENTS
c-------------------------------------------------------------
      do i=1,4
       do j=1,4
         R(i,j)= RSM(i,j) + RDM(i,j)
       enddo
      enddo	   

      return
      end	






