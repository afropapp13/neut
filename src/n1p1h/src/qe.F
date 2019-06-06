******************************************************************************
*******         READ(4,*)DZZ,DAA,KLAVE,DNCXP,DNCA0P,DNCXN,DNCA0N
*******         READ($,*)qvalue,lepton
*******         READ(4,*) mn
*******         READ(4,*)ILIN
*******         IF (ILIN.EQ.1) THEN
*******          'LINHARD NORMAL'
*******         ELSEIF (ILIN.EQ.2) THEN
*******          'LINHARD RELATIVISTA'
*******         ELSEIF (ILIN.EQ.3) THEN
*******          'LINHARD F.SPECTRALES'
*******         ENDIF
******************************************************************************
*******  Units FMs and qvector goes in  z direction 

      SUBROUTINE ConstantsInitialization(xmagev,irpa)
      IMPLICIT REAL*8 (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      COMMON/datos/dpi,hbarc,GF0,DMNU,DMA
      COMMON/datos2/dmneutrino,dmlepton,dmi,dmf,
     f     coscabibbo,dmuon,dmelectron,dmtau,xuma

      INTEGER*4 irelativistic
      COMMON/rpavariables/fp0in,fp0ex,f,fstar,pilambda,
     f     cr0,rholambda,gp,xmpi,xmrho,irelativistic

      COMMON/testieta/ieta
      COMMON/enregistre/ipol
      COMMON/mavalue/xma      

      ipol=irpa

      if (ipol.ne.0 .and. ipol.ne.1 .and. ipol.ne.2 ) stop'wrong ipol'

      DPI=3.141592653589793D0
      hbarc = 197.3269602d0
      dmneutrino = 0.d0 / hbarc
      dmuon = 0.105658357*1000. / hbarc ! muon mass
      dmelectron = 0.000510998910*1000./ hbarc ! electron mass        
      dmtau = 1.77682*1000. / hbarc ! tau mass        
      dmnu = 939.d0 / hbarc
      coscabibbo=0.974d0
      GF0= 1.1664d-11*hbarc**2  ! units fm^2
      xma= xmagev*1.d3/ hbarc
      xuma = 931.494013 / hbarc

      PRINT*,"------------------------------------"
      PRINT*,"         CONSTANTS                  "
      PRINT*,"------------------------------------"
      PRINT*,"HBARC ",hbarc
      PRINT*,"Neutrino mass ",dmneutrino*hbarc
      PRINT*,"Muon  mass    ",dmuon*hbarc
      PRINT*,"Elec  mass    ",dmelectron*hbarc
      PRINT*,"Nucl  mass    ",dmnu*hbarc
      PRINT*,"Cabibbo       ",coscabibbo
      PRINT*,"GF0           ",GF0 
      PRINT*,"Axial mass    ",xma*hbarc
      PRINT*,"------------------------------------"

      RETURN 
      END


      SUBROUTINE rpaparameters(fp0inval,fp0exval,f0val,fstarval,
     f     pilambval,cr0val,rholambval,gpval,xmpival,xmrhoval,ir)
      IMPLICIT REAL*8 (A-H,O-T,V-Z)
      
      COMMON/datos/dpi,hbarc,GF0,DMNU,DMA
      
      INTEGER*4 irelativistic
      COMMON/rpavariables/fp0in,fp0ex,f,fstar,pilambda,
     f     cr0,rholambda,gp,xmpi,xmrho,irelativistic
      
      LOGICAL ISINITIALIZED
      COMMON/FLAGS/ISINITIALIZED
      
      if( .NOT.ISINITIALIZED ) THEN 
         PRINT*, " IT IS NOT INITIALIZED " 
         STOP 0
      ENDIF 
      
c     
c     RPA Parameters.       
c     
      fp0in = fp0inval
      fp0ex = fp0exval
      f = f0val
      fstar = fstarval
      pilambda =  pilambval
      cr0 = cr0val
      rholambda = rholambval
      gp = gpval
      xmpi  = xmpival/hbarc
      xmrho = xmrhoval/hbarc
      
      irelativistic = ir

      PRINT*,"---------------------------------------"
      PRINT*,"    Modified RPA Parameters            "
      PRINT*,"---------------------------------------"
      PRINT*," FP0 IN = ",fp0in
      PRINT*," FP0 EX = ", fp0ex
      PRINT*," F      = ",f
      PRINT*," F*     = ",fstar
      PRINT*," Lam_pi = ",pilambda
      PRINT*," Cr0    = ",cr0
      PRINT*," Lam_rho= ",rholambda
      PRINT*," gp     = ",gp 
      PRINT*," xmpi   = ",xmpi 
      PRINT*," xmrho  = ",xmrho 
      PRINT*," relat  = ",irelativistic
      PRINT*,"---------------------------------------"
      
      RETURN 
      end 

      SUBROUTINE Initialization(xmagev,irpa,rm,ianti,type,icontrlepton)
      IMPLICIT REAL*8 (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      DIMENSION drp(2000),df0(2000)
      
      LOGICAL  ISINITIALIZED
      COMMON/FLAGS/ISINITIALIZED 
      
      COMMON/NXRO/NXRO
      COMMON/VFPI/DROP,DRON
      COMMON/CTEROPN/DXP,DA0P,DRO0P,DXN,DA0N,DRO0N
      
      COMMON/NCXRO/NCXRO
      COMMON/CTEROPNC/DNCXP,DNCA0P,DRONC0P
      
      COMMON/ROFEROS/KLAVE
      COMMON/NUC/DZZ,DAA
      
      COMMON/datos/dpi,hbarc,GF0,DMNU,DMA
      COMMON/datos2/dmneutrino,dmlepton,dmi,dmf,
     f     coscabibbo,dmuon,dmelectron,dmtau,xuma

      INTEGER   NR
      DIMENSION dr(2000),vcd(2000)
      COMMON/vcfto/vcd,dr,nr
  
      COMMON/mavalue/xma
      
      COMMON/testieta/ieta
      COMMON/Rinit/RMAX
      COMMON/Qvalue/qvalue,qvaluecorr
      COMMON/enregistre/ipol
      COMMON/VERIF/ILIN,MN,MNR
      
      INTEGER*4 irelativistic
      COMMON/rpavariables/fp0in,fp0ex,f,fstar,pilambda,
     f     cr0,rholambda,gp,xmpi,xmrho,irelativistic
      
      ieta=ianti

      if(type.eq.12.)then  ! This is Carbon 
         DZZ=6.d0
         DAA=12.d0
         KLAVE=1
         DNCXP=1.692d0
         DNCA0P=1.082d0
         DNCXN=1.692d0
         DNCA0N=1.082d0
         if(ieta.eq.1)then
            qp=0.601
            qvalue=16.827+qp
         elseif(ieta.eq.-1)then
            qn=3.370
            qvalue=13.880+qn
         else
            stop'not a neutrino nor an antineutrino'
         endif    
      elseif(type.eq.16.)then  ! This is Oxygen 
         DZZ=8.d0
         DAA=16.d0
         KLAVE=1
         DNCXP=1.833d0
         DNCA0P=1.544d0
         DNCXN=1.833d0
         DNCA0N=1.544d0
         if(ieta.eq.1)then
            qp=-0.536
            qvalue=14.334+qp
         elseif(ieta.eq.-1)then
            qn=2.489
            qvalue=4.642359+qn
         else
            stop'not a neutrino nor an antineutrino'
         endif    

      elseif(type.eq.28.)then  ! This is Silicon
         DZZ=14.d0
         DAA=28.d0
         KLAVE=0
         DNCXP=2.93d0
         DNCA0P=0.569d0
         DNCXN=2.86d0
         DNCA0N=0.569d0
         if(ieta.eq.1)then
            qp=2.066
            qvalue=14.332+qp
         elseif(ieta.eq.-1)then
            qn=7.725
            qvalue=4.642+qn
         else
            stop'not a neutrino nor an antineutrino'
         endif    

      elseif(type.eq.40.)then  ! This is Argon. Not Calcium 
         DZZ=18.d0
         DAA=40.d0
         KLAVE=0
         DNCXP=3.47d0
         DNCA0P=0.569d0
         DNCXN=3.64d0
         DNCA0N=0.569d0
         if(ieta.eq.1)then
            qp=7.582
            qvalue=0.994+qp
         elseif(ieta.eq.-1)then
            qn=5.830
            qvalue=7.991+qn
         else
            stop'not a neutrino nor an antineutrino'
         endif    

      elseif(type.eq.56.)then  ! This is Iron 
       DZZ=26.d0
       DAA=56.d0
       KLAVE=0             ! Missing these values for Fe
       DNCXP=3.971d0
       DNCA0P=0.5935d0
       DNCXN=4.05d0
       DNCA0N=0.5935d0
      if(ieta.eq.1)then    ! These values are fine. 
       qp=5.849
       qvalue=3.045+qp
      elseif(ieta.eq.-1)then
       qn=8.257
       qvalue=4.2065+qn
      else
       stop'not a neutrino nor an antineutrino'
      endif 
 
      elseif(type.eq.27.)then  ! This is Aluminum 
       DZZ=13.d0
       DAA=27.d0
       KLAVE=0             ! Missing these values for Al
       DNCXP=3.05d0
       DNCA0P=0.535d0
       DNCXN=3.05d0
       DNCA0N=0.535d0
      if(ieta.eq.1)then    ! These values are fine. 
       qp=7.46296 
       qvalue=4.3+qp
      elseif(ieta.eq.-1)then
       qn=6.44339
       qvalue=3.1+qn
      else
      stop'not a neutrino nor an antineutrino'
      endif    
      
      elseif(type.eq.208.)then   ! This is Aluminum 
       DZZ=82.d0
       DAA=126.d0
       KLAVE=0             ! Missing these values for Al
       DNCXP=6.624d0
       DNCA0P=0.54d0
       DNCXN=6.89d0
       DNCA0N=0.54d0
      if(ieta.eq.1)then    ! These values are fine. 
       qp=3.706 
       qvalue=2.8784+qp
      elseif(ieta.eq.-1)then
       qn=3.789
       qvalue=4.998956+qn
      else
      stop'not a neutrino nor an antineutrino'
      endif    
                                 
      else
         stop'nucleus non implemented'
      endif
      
      DMA = DAA*xuma            ! nucleus mass.
c
c     read(4,*) qvalue,lepton
c      
      
      lepton=icontrlepton
      if (lepton.eq.1) then 
         dmlepton = dmuon 
      elseif (lepton.eq.0) then
         dmlepton = dmelectron 
      else
         stop'Wrong lepton identification (1 or 0)' 
      endif
      
      dmi = dmneutrino
      dmf = dmlepton
      
      mn=2
      mnr=10
      
c     
c     read(4,*)ieta
c     ieta = 1 
      
c     rewind(4)          
      
      
c********************************************************************
c********************************************************************
      ilin=2
c********************************************************************
c********************************************************************
      
      
      nxro=0                    ! density normalization 
      ncxro=0                   ! density normalization
      
      lc = 0                    ! 0 desconvoluciona cualquier otro numero no
      
      call convolucion (dncxp,dnca0p,dncxn,dnca0n,klave,lc,
     f     dxp,da0p,dxn,da0n)
      
      CALL XRO(0.D0)
      
      IF (KLAVE.EQ.0) THEN
         RMAX= MAX(DXP,DXN)+9.25*MAX(DA0P,DA0N)
      ELSEIF (KLAVE.EQ.1) THEN
         RMAX=DSQRT(20.D0)*MAX(DXN,DXP)
      ENDIF
      
c      rm=RMAX
      
      if (ipol.eq.1) then
c
c     Finite Coulomb Potencial calculation 
c        
         call dsg20r(0.d0,rmax,mnr,dr,nr)

         alpha = 1.d0/137.036d0
         
         do ir=1,nr
            r = dr(ir)
            
            CALL DSG20R (0.D0,R,5,DRP,NRP)
            do irp = 1,nrP
               rp=drp(irp) 
               df0 (irp) = (rp**2)*densq(rp)/r
            enddo                   
            CALL DRG20R (0.D0,R,5,df0,f1)
            
            CALL DSG20R (R,RMAX,5,DRP,NRP)
            do irp = 1,nrP
               rp=drp(irp)
               df0 (irp) = rp*densq(rp)
            enddo                  
            CALL DRG20R (R,RMAX,5,df0,f2)
            
            vcd(ir)=-alpha*ieta*4.d0*dpi*(f1+f2)
            
         enddo
         
      endif
      
      ISINITIALIZED = .TRUE.
      
      RETURN
      END
      


      SUBROUTINE SelectLepton(ilepton) 
      IMPLICIT REAL*8 (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      COMMON/datos2/dmneutrino,dmlepton,dmi,dmf,
     f     coscabibbo,dmuon,dmelectron,dmtau,xuma

      lepton=ilepton
      if(lepton.eq.2) then 
         dmlepton = dtau 
      else if (lepton.eq.1) then 
         dmlepton = dmuon 
      elseif (lepton.eq.0) then
         dmlepton = dmelectron 
      else
         stop 'Wrong lepton identification (1 or 0)' 
      endif

      dmi = dmneutrino
      
      dmf = dmlepton

      RETURN 
      END


      SUBROUTINE SetLeptonMass(ilepton,fmass) 
      IMPLICIT REAL*8 (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      COMMON/datos/dpi,hbarc,GF0,DMNU,DMA
      COMMON/datos2/dmneutrino,dmlepton,dmi,dmf,
     f     coscabibbo,dmuon,dmelectron,dmtau,xuma

      lepton=ilepton
      
      dmlepton = fmass/hbarc 
     
      dmi = dmneutrino
      
      dmf = dmlepton
      
      RETURN 
      END


      

      
      subroutine vcinterpolate(R,vc)
      IMPLICIT REAL*8 (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      
      COMMON/datos/dpi,hbarc,GF0,DMNU,DMA
      COMMON/datos2/dmneutrino,dmlepton,dmi,dmf,
     f     coscabibbo,dmuon,dmelectron,dmtau,xuma

      INTEGER   NR
      DIMENSION dr(2000),vcd(2000)
      COMMON/vcfto/vcd,dr,nr
 
      IMIN = 1
      IMAX = NR
c     
c Find radious with Newton's method and interpolate the potential that is precomputed. 
c     
      IF( R .LT.DR(IMIN) ) THEN
         IMAX = IMIN+1
      ELSE IF ( R .GT.DR(IMAX) ) THEN
         IMIN = IMAX-1
      ELSE 
         DO WHILE ( (IMAX-IMIN).GT.1 )
            
            IMID = (IMIN+IMAX)/2
            
            IF( R < DR(IMID) ) THEN
               IMAX = IMID
            ELSE
               IMIN = IMID
            ENDIF
            
         END DO
      ENDIF
      
      VC=((VCD(IMAX)-VCD(IMIN))/(DR(IMAX)-DR(IMIN))*
     f     (R-DR(IMIN))+VCD(IMIN))*(hbarc/1.d3)       
c     Here vc is calculated in GeV for external use(NRC=10.)
c     change is made to MeV/hbarc for internal use (NRC=1.and10.)in WSELF
      return
      end

      double precision function targetminmon(q0,dq,dmt,dmp)
      IMPLICIT DOUBLE PRECISION (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      COMMON/datos/dpi,hbarc,GF0,D,DMA

      ZLOC = q0*q0-dq*dq;
      DIFF = -dmp*dmp+dmt*dmt
      ALOC = 4.*ZLOC
      BLOC = 4.*dq*(ZLOC+DIFF)
      CLOC = -(ZLOC+DIFF)*(ZLOC+DIFF)+4.*dmt*dmt*q0*q0
c      BLOC = 4.*dq*(ZLOC)
c      CLOC = (ZLOC)*(ZLOC)-4.*dmt*dmt*q0*q0

      if (BLOC*BLOC-4.*ALOC*CLOC.gt.0) then
         Delta=DSQRT(BLOC*BLOC-4.*ALOC*CLOC)
      else if (abs(BLOC*BLOC-4.*ALOC*CLOC).le.1D-6) then
         Delta=0.
      else
         write(*,*) 'Weired BLOC*BLOC-4*ALOC*CLOC value'
         write(*,*) 'ALOC=',ALOC
         write(*,*) 'BLOC=',BLOC
         write(*,*) 'CLOC=',CLOC
         Delta=DSQRT(BLOC*BLOC-4.*ALOC*CLOC)
      endif

       Sol1=(-BLOC-Delta)/(2.*ALOC)
       Sol2=(-BLOC+Delta)/(2.*ALOC)
c
c      Sol3=(-BLOC-Delta)/(2.*ALOC)
c      Sol4=(-BLOC+Delta)/(2.*ALOC)
c      if(DABS(Sol1).ne.DABS(Sol4).or.DADp = PF1 BS(Sol2).ne.DABS(Sol3))then
c         stop'error: wrong solution computation'
c      endif
c

       targetminmon=MIN(DABS(Sol1),DABS(Sol2))

c      targetminmon = 0.

       return
       end

c************************************************************************
c************************************************************************
c************************************************************************
c************************************************************************
c************************************************************************

       double precision function SigthNacho(eingev,tmugev,cost) 

       IMPLICIT REAL*8 (A-H,O-T,V-Z)
       IMPLICIT COMPLEX*16 (U)
       DIMENSION WD(5)           ! Response functions 

       COMMON/NUC/DZZ,DAA

       COMMON/datos/dpi,hbarc,GF0,DMNU,DMA
       COMMON/datos2/dmneutrino,dmlepton,dmi,dmf,
     f     coscabibbo,dmuon,dmelectron,dmtau,xuma

       COMMON/mavalue/xma

       COMMON/testieta/ieta
       COMMON/Rinit/RMAX
       COMMON/Qvalue/qvalue,qvaluecorr
       COMMON/enregistre/ipol
       COMMON/VERIF/ILIN,MN,MNR

       LOGICAL ISINITIALIZED
       COMMON/FLAGS/ISINITIALIZED

       if( .NOT.ISINITIALIZED ) THEN 
          PRINT*, " IT IS NOT INITIALIZED " 
          STOP 0
       ENDIF 

       cosaux= cost
       theta=acos(cosaux)        !60.d0*dpi/180.d0

       eout = dmf+(tmugev*1.d3)/hbarc 
       POUT=DSQRT(EOUT**2-DMF**2)

       ein = eingev*1d3/hbarc
       PIN=DSQRT(EIN**2-DMI**2)

c       print*,qvalue

       Q0=EIN-eout-qvalue/hbarc

       if (q0.ge.0) then                  
          DQ=DSQRT(PIN**2+POUT**2-2.D0*DCOS(THETA)*PIN*POUT)

          call secciondif(q0,dq,theta,ein,pin,eout,pout,ilin,
     f        ieta,rmax,mn,mnr,wd,sig,ipol)
       else
          sig=0.d0
       endif

       if (ieta.eq.1) then 
          dnumber=daa-dzz
       elseif (ieta.eq.-1) then
          dnumber=dzz
       else
          stop'ieta wrong'
       endif

       SigthNacho = sig*2.d0*dpi*1.d18/hbarc/dnumber

c     if( SigthNacho.gt.10000 ) then 
c     PRINT*," LARGE ",DQ*hbarc,Q0*hbarc,cosaux,EIN,PIN,EOUT,POUT
c     endif
c-----------------------------------------------------------
c     UNITS 10^{-41} cm^2/GeV/num_(active nucleons)
c-----------------------------------------------------------
       RETURN
       END

***************************************************************
****************************************************************
****************************************************************
****************************************************************      

       double precision function SigthBruno(eingev,tmugev,cost,NRC,NPC, 
     f     R,dP,vc)
       IMPLICIT REAL*8 (A-H,O-T,V-Z)
       IMPLICIT COMPLEX*16 (U)
       DIMENSION WD(5)           ! Response functions

       COMMON/NUC/DZZ,DAA

       COMMON/datos/dpi,hbarc,GF0,DMNU,DMA
       COMMON/datos2/dmneutrino,dmlepton,dmi,dmf,
     f                coscabibbo,dmuon,dmelectron,dmtau,xuma

       COMMON/mavalue/xma
       COMMON/testieta/ieta
       COMMON/Rinit/RMAX
       COMMON/Qvalue/qvalue,qvaluecorr
       COMMON/VERIF/ILIN,MN,MNR
       COMMON/enregistre/ipol

       COMMON/ENERGY/Enu,Emu,Eq,En,Ep
       COMMON/ANGLE/theta,thetap
       COMMON/momentum/pm,dqt,dPt,P
       COMMON/coulomb/coul

       LOGICAL ISINITIALIZED
       COMMON/FLAGS/ISINITIALIZED

       if( .NOT.ISINITIALIZED ) THEN 
          PRINT*, " IT IS NOT INITIALIZED " 
          STOP 0
       ENDIF 

       sig = 0.0d0
       SigthBruno = 0.0d0 

       Precup = (dp*1.d3)/hbarc  ! to fermi 
       Rrecup = R                ! already fm units. 

       cosaux= cost
       theta=acos(cosaux)        ! 60.d0*dpi/180.d0      

       eout = dmf+(tmugev*1.d3)/hbarc
c     if( EOUT .lt. DMF ) RETURN
       POUT=DSQRT(EOUT**2-DMF**2)

c     SigthBruno = -2. 
       ein = eingev*1d3/hbarc
c     if( EIN .lt. DMI ) RETURN 
       PIN=DSQRT(EIN**2-DMI**2)

c      PRINT*," DMF ",DMF*hbarc

       Q0 = ein-eout-qvaluecorr*1000./hbarc
       Dmp = 0.93827208*1000./hbarc
       Dmn = 0.93956542*1000./hbarc

c       print*,"!!!!! ",Q0*hbarc,qvaluecorr*1000.

c     if (q0.gt.0) then                   
       DQ=DSQRT(PIN**2+POUT**2-2.D0*DCOS(THETA)*PIN*POUT)
       Q2=Q0**2-DQ**2

       factor = pout*eout*dma*(coscabibbo*GF0)**2/(dpi**2)
       cte = dmlepton**2/(eout*(eout+pout))

       sin12 = sin(theta/2)  
       cos12 = dsqrt(1.d0 - sin12**2)
       coseno = cos12**2-sin12**2             

       CALL WSELF(ILIN,Q0,dq,coseno,Q2,RMAX,MN,mnr,WD,ipol,eout,pin,
     f     ieta,NRC,NPC,Rrecup,Precup,CTHETAP,Panc,vc)

       if (abs(cthetap).ge.1.) then
          goto 10
       endif

       R=Rrecup
c     
c     Different target mass for neutrinos and antineutrinos. 
c     
       if( ieta.eq.1 ) then   
          Epqanc=DSQRT(dmp**2+Panc**2)
          Epanc=Epqanc-Q0
       else 
          Epqanc=DSQRT(dmn**2+Panc**2)
          Epanc=Epqanc-Q0
       endif 

ccccccccccccccoulomb potential cccccccccccccccccccccccccccc           
       vc=vc*hbarc
       
ccccccccccccccccccccccccenergycccccccccccccccccccccccccccccccccccccccccccccc
       Enu=ein*hbarc
       Emu=eout*hbarc
       Eq=q0*hbarc
       En=Epanc*hbarc
       Ep=Epqanc*hbarc
cccccccccccccccccccccccimpulsion and anglecccccccccccccccccccccccccccccccccc
       thetap=acos(cthetap)
       pm=Pout*hbarc
       dqt=DQ*hbarc
       dPt= (dp*1.d3)
       P=Panc*hbarc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     if( WD(1) .EQ. 0 ) return

      dsum1 = 2.d0*WD(1)*sin12**2 + WD(2)*cos12**2 
     f     -ieta*WD(3)*(ein+eout)*sin12**2/dma
      
      dsum2 = (wd(1)-wd(2)/2.d0)*coseno+0.5d0*ieta*wd(3)*
     f     ((eout+pout)/dma -(eout+ein)*coseno/dma)+
     f     0.5d0*wd(4)*(coseno*dmlepton**2/dma**2 + 
     f     2.d0*eout*(eout+pout)*sin12**2/dma**2) -
     f     0.5d0*wd(5)*(eout+pout)/dma
      
      SIG= factor*(dsum1+cte*dsum2)

c     if( NPC.EQ.10 ) print*," >>" ,factor,dsum1,cte,dsum2,sig
c     else
c     SigthBruno = -3.
c     return 
c     endif

 10   if( SIG .LE. 0 ) then 
         vc=0
ccccccccccccccccccccccccenergycccccccccccccccccccccccccccccccccccccccccccccc
         Enu=0
         Emu=0
         Eq=0
         En=0
         Ep=0
cccccccccccccccccccccccimpulsion and anglecccccccccccccccccccccccccccccccccc
         thetap=0
         pm=0
         dqt=0
         dPt=0
         P=0
c     SigthBruno = -4.
         return
c     
      endif
cccccccccccccccccRETURN proton momemtum
      dP=P/1d+3   
      coul=vc
c     
      if (ieta.eq.1) then 
         dnumber=daa-dzz 
      elseif (ieta.eq.-1) then
         dnumber=dzz
      else
         stop'ieta wrong'
      endif
      
      SigthBruno = sig*2.d0*dpi*1.d18/hbarc/dnumber
c-----------------------------------------------------------
c     UNITS 10^{-41} cm^2/GeV/num_(active nucleons)/fermi(R)/fermi(Fermi mom) 
c-----------------------------------------------------------

      RETURN
      END
      
      double precision function Sigthfree(eingev,tmugev,cost)
      IMPLICIT REAL*8 (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)

      COMMON/datos/dpi,hbarc,GF0,DMNU,DMA
      COMMON/datos2/dmneutrino,dmlepton,dmi,dmf,
     f     coscabibbo,dmuon,dmelectron,dmtau,xuma
        
      COMMON/mavalue/xma
      COMMON/testieta/ieta
      
      COMMON/ENERGY/Enu,Emu,Eq,En,Ep
      COMMON/ANGLE/theta,thetap
      COMMON/momentum/pm,dqt,dPt,P

      LOGICAL ISINITIALIZED
      COMMON/FLAGS/ISINITIALIZED
      
      if( .NOT.ISINITIALIZED ) THEN 
         PRINT*, " IT IS NOT INITIALIZED " 
         STOP 0
      ENDIF 
   
      if( ieta.eq.1 ) then 
         Sigthfree = 0
         return
      endif

      Dmp=0.93827208*1000./hbarc
      Dmn=0.93956542*1000./hbarc
      if(ieta.eq.1) dmnui=dmn
      if(ieta.eq.-1) dmnui=dmp
      if(ieta.eq.1) dmnuf=dmp
      if(ieta.eq.-1) dmnuf=dmn
c      dmnui=DMNU
c      dmnuf=DMNU
      eout = dmf+(tmugev*1.d3)/hbarc 
      POUT=DSQRT(EOUT**2-DMF**2)
      
      ein = eingev*1d3/hbarc
      PIN=DSQRT(EIN**2-DMI**2)
      
      Eex=Ein-Eout+dmnui

      cost=(dmnuf**2+Ein**2+Pout**2-Eex**2)/(2*Ein*Pout)
      cosaux= cost
      if( abs(cosaux).gt.1.D0 ) then 
         Sigthfree = 0
         return
      endif

      theta=acos(cosaux)

      Q0=ein-eout
      if (q0.ge.0.and.abs(cost).le.1.) then                        
         DQ=DSQRT(PIN**2+POUT**2-2.D0*DCOS(THETA)*PIN*POUT)
         Q2=Q0**2-DQ**2
         S=(2*ein+dmnui)*dmnui
         CALL WFREE(ILIN,Q0,dq,Q2,tensLA,ieta,S,dmnui)

         Q2test=dmnui**2+dmnuf**2-2*dmnui*Eex 
         IF( ABS(Q2-Q2test).GT.1.d-6 )THEN 
            PRINT*,"ERROR IN Q2 ",Q2,Q2test
         ENDIF 

         if(q2.gt.0) then
            Print*," Error in q2 ",q2
            sig=0.d0
         endif 

         factfree=(coscabibbo*GF0)**2/(8*dpi*(S-dmnui**2)**2)
         sig=tensLA*factfree
      else
         sig=0.d0
      endif

      Enu=ein*hbarc
      Emu=eout*hbarc
      Eq=Q0*hbarc
      if(ieta.eq.1) then
         Ep=Eex
         En=dmnui
      elseif(ieta.eq.-1)then
         En=Eex
         Ep=dmnui
      endif
      pm=Pout*hbarc
      dqt=DQ*hbarc
      dPt=0.d0
      thetap=0.d0
      
      aJacob=2.d0*dmnui

      Sigthfree = sig*aJacob/2.d0*1.d18/hbarc
cccccccccccc   dxsect/dq2*Jacobian/thetaangle*unit->10**41cm2/GeV     
      RETURN   
      END

***************************************************************
      SUBROUTINE comptmomentum(dnu)      
      IMPLICIT REAL*8 (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      
C      DATA DPI/3.141592653589793D0/

      COMMON/datos/dpi,hbarc,GF0,DMNU,DMA
      COMMON/datos2/dmneutrino,dmlepton,dmi,dmf,
     f     coscabibbo,dmuon,dmelectron,dmtau,xuma
      DIMENSION dt1(3),dt2(3),dt3(3),dnu(3)   
      DIMENSION y(3),clone(3)
      DIMENSION vPnu(4),vPmu(4),vPn(4),vPp(4),vq(4)
      COMMON/fourvectors/vPnu,vPmu,vq,vPn,vPp
      
      COMMON/ENERGY/Enu,Emu,Eq,En,Ep
      COMMON/ANGLE/theta,thetap
      COMMON/momentum/pm,dqt,dPt,P
      COMMON/coulomb/coul

#ifndef __GFORTRAN__
      INTEGER*4 IDUM
      REAL*4    RLU
      EXTERNAL  RLU
#endif

      DPI = 3.141592653589793D0
#ifdef __GFORTRAN__
      CALL RANDOM_NUMBER(A)
#else
      A = RLU(IDUM)
#endif      
      PHI=2.D0*DPI*A

      ct=dcos(theta)
      st=dsin(theta)
     
      cp=dcos(phi)
      sp=dsin(phi)
C
C  Reference system
C
      dnunorm = dsqrt(dnu(1)*dnu(1)+dnu(2)*dnu(2)+dnu(3)*dnu(3))

      dt3(1)  = dnu(1)/dnunorm
      dt3(2)  = dnu(2)/dnunorm
      dt3(3)  = dnu(3)/dnunorm
c
c If dt3 is parallel to y then use x as perpendicular. 
c
      if( dt3(2).eq.1 ) then   
         y(1) = 1. 
         y(2) = 0.
         y(3) = 0.
      else 
         y(1) = 0. 
         y(2) = 1. 
         y(3) = 0.
      endif

      dt1(1) =  y(2)*dt3(3)-y(3)*dt3(2)
      dt1(2) =  y(3)*dt3(1)-y(1)*dt3(3)
      dt1(3) =  y(1)*dt3(2)-y(2)*dt3(1)

      dt1norm = dsqrt(dt1(1)*dt1(1)+dt1(2)*dt1(2)+dt1(3)*dt1(3))

      dt1(1)  = dt1(1) / dt1norm 
      dt1(2)  = dt1(2) / dt1norm 
      dt1(3)  = dt1(3) / dt1norm 

      dt2(1) =  dt3(2)*dt1(3)-dt3(3)*dt1(2)
      dt2(2) =  dt3(3)*dt1(1)-dt3(1)*dt1(3)
      dt2(3) =  dt3(1)*dt1(2)-dt3(2)*dt1(1)
C
C   Energy 
C
      vPnu(1) = Enu
      vPmu(1) = Emu
      vq(1)   = Eq
      vPn(1)  = En
      vPp(1)  = Ep

      print*," 3 !!!!!!!!!! ",Eq*hbarc


c
C     coulomb potential modification
c
      Emuloc=Emu-coul
      
      if( Emuloc.lt.(DMF*hbarc) ) then
         Emuloc = DMF*hbarc
      endif
      
      pmuloc=DSQRT(Emuloc**2-(DMF*hbarc)**2)
C
C     3-momentum for the lepton.
C
      do i = 2 , 4 
         vPnu(i)= dt3(i-1)*Enu
         vPmu(i)=(dt3(i-1)*ct + dt1(i-1)*st*cp + dt2(i-1)*st*sp)*pm
         clone(i-1)=vPmu(i)*pmuloc/pm
c         vq(i)=vPnu(i)-vPmu(i)
         vq(i)=vPnu(i)-clone(i-1)
      enddo
C
C  Build the new reference system along q3 
C
      vqnorm = dsqrt(vq(2)*vq(2)+vq(3)*vq(3)+vq(4)*vq(4))

      dt3(1)  = vq(2) / vqnorm 
      dt3(2)  = vq(3) / vqnorm 
      dt3(3)  = vq(4) / vqnorm 

      if( dt3(2).eq.1 ) then   
         y(1) = 1 
         y(2) = 0 
         y(3) = 0
      else 
         y(1) = 0 
         y(2) = 1 
         y(3) = 0
      endif

      dt1(1) =  y(2)*dt3(3)-y(3)*dt3(2)
      dt1(2) =  y(3)*dt3(1)-y(1)*dt3(3)
      dt1(3) =  y(1)*dt3(2)-y(2)*dt3(1)

      dt1norm = dsqrt(dt1(1)*dt1(1)+dt1(2)*dt1(2)+dt1(3)*dt1(3))

      dt1(1)  = dt1(1) / dt1norm 
      dt1(2)  = dt1(2) / dt1norm 
      dt1(3)  = dt1(3) / dt1norm 

      dt2(1) =  dt3(2)*dt1(3)-dt3(3)*dt1(2)
      dt2(2) =  dt3(3)*dt1(1)-dt3(1)*dt1(3)
      dt2(3) =  dt3(1)*dt1(2)-dt3(2)*dt1(1)

#ifdef __GFORTRAN__
      CALL RANDOM_NUMBER(B)
#else
      B = RLU(IDUM)
#endif      
      PHIP=2.D0*DPI*B
      
      ctp=dcos(thetap)
      stp=dsin(thetap)

      cpp=dcos(phip)
      spp=dsin(phip)
C
C     3-momentum for the hadrons.
C
      do i = 2 , 4 
         vPn(i)=dPt*(dt3(i-1)*ctp + dt1(i-1)*stp*cpp + dt2(i-1)*stp*spp)
         vPp(i)= vq(i)+vPn(i)

      enddo

c      print*,En,dsqrt(En*En-vPn(2)*vPn(2)-vPn(3)*vPn(3)-vPn(4)*vPn(4))
c      print*,P-dsqrt(vPp(2)*vPp(2)+vPp(3)*vPp(3)+vPp(4)*vPp(4)),ctp


      RETURN
      END

******************************************************************
*****************************************************************

      SUBROUTINE secciondif(q0,dq,theta,ein,pin,eout,pout,ilin,
     f     ieta,rmax,mn,mnr,wd,sig,ipol) 
      
      IMPLICIT REAL*8 (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      
      DIMENSION wd(5)
      COMMON/datos/dpi,hbarc,GF0,DMNU,DMA
      COMMON/datos2/dmneutrino,dmlepton,dmi,dmf,
     f     coscabibbo,dmuon,dmelectron,dmtau,xuma
      
      Q2=Q0**2-DQ**2
      factor = pout*eout*dma*(coscabibbo*GF0)**2/(dpi**2)
      cte = dmlepton*dmlepton/(eout*(eout+pout))
      sin12 = sin(theta/2)  
      cos12 = dsqrt(1.d0 - sin12**2)
      coseno = cos12**2-sin12**2             
      
      CALL WSELF(ILIN,Q0,dq,coseno,Q2,RMAX,MN,mnr,WD,ipol,eout,pin,
     f     ieta,0,0,R,P,CTHETAP,Panc,vc)
      
      dsum1 = 2.d0*WD(1)*sin12*sin12 + WD(2)*cos12*cos12
     f     -ieta*WD(3)*(ein+eout)*sin12*sin12/dma

      dma2 = dma*dma
      
      dsum2 = (wd(1)-wd(2)/2.d0)*coseno+0.5d0*ieta*wd(3)*
     f     ((eout+pout-(eout+ein)*coseno)/dma)+
     f     0.5d0*wd(4)*(coseno*dmlepton**2/dma2 + 
     f     2.d0*eout*(eout+pout)*sin12**2/dma2) -
     f     0.5d0*wd(5)*(eout+pout)/dma
      
      
c      print*,factor,dsum1,dsum2,wd(1),wd(2),wd(3),wd(4),wd(5)

      SIG= factor*(dsum1+cte*dsum2)
      
      
      RETURN
      end



      
      SUBROUTINE WSELF(ILIN,Q0,dq,coseno,Q2,RMAX,MN,mnr,WD,ipol,
     f     eout,pin,ieta,NRC,NPC,Rrecup,Precup,CTHETAP,Panc,vc)
     
******* TODAS LAS UNIDADES EN FMS; Q0 Y DQ DEFINEN EL4-MOMENTO DEL FOTON
******* RMAX DEFINEL RANGO DE INTEGRACION Y MN EL NUMERO DE PTOS.
******* GAUSS. FINALMENTE WL Y WT SON LAS FUNCIONES DE ESTRUCTURA.
*******         IF (ILIN.EQ.1) THEN
*******          'LINHARD NORMAL'
*******         ELSEIF (ILIN.EQ.2) THEN
*******          'LINHARD RELATIVISTA'
*******         ELSEIF (ILIN.EQ.3) THEN
*******          'LINHARD F.SPECTRALES'
*******         ENDIF
*******         Unidades FMs y qvector va en la direccion z
 
      IMPLICIT DOUBLE PRECISION (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      complex*16 cunuc, delta_lind
      
      DIMENSION DFXX(2000),DFZZ(2000)
      DIMENSION DF0Z(2000),DFXY(2000),DF00(2000)
      DIMENSION tulin(0:3),rulin(0:3,0:3),WD(5)
      DIMENSION drp(2000),df0(2000)
      
      COMMON/VFPI/DROP,DRON
      COMMON/datos/dpi,hbarc,GF0,DMNU,DMA
c     COMMON/datos3/xmpi,gp
      COMMON/densidad/dro,dro0
      COMMON/cunuc/cunuc

      INTEGER   NR 
      DIMENSION dr(2000),vcd(2000)
      COMMON/vcfto/vcd,dr,nr

      COMMON/fsi/xlindfs
      
      COMMON/datos2/dmneutrino,dmlepton,dmi,dmf,
     f     coscabibbo,dmuon,dmelectron,dmtau,xuma
      COMMON/mavalue/xma
      
      INTEGER*4 irelativistic
      COMMON/rpavariables/fp0in,fp0ex,f,fstar,pilambda,
     f     cr0,rholambda,gp,xmpi,xmrho,irelativistic

#ifndef __GFORTRAN__
      INTEGER*4 IDUM
      REAL*4    RLU
      EXTERNAL  RLU
#endif

      if (ilin.ne.2) stop'not ready ilin.ne.2 '
 
      cte_r= -4.d0*dpi/(8.d0*dpi*dmnu*dmnu)

      if( NRC.EQ.0 ) THEN        
         NRI = NR
      ELSE 
         NRI = 1
      ENDIF
      
cccc  NRC control on loop or not 0->loop/other don t loop  

      DO IR=1,NRI

         iflag=0
c     q0old = q0
         dqold = dq
         q2old = q2
c         vc= 0.d0
c         qth = 0.d0
  
         if( NRC.eq.0 ) then    ! Loping 
            R = DR(IR)
            FACTINT = 1
         else if (NRC.eq.1)then ! Not looping 
#ifdef __GFORTRAN__
            CALL RANDOM_NUMBER(Z)
#else
            Z = RLU(IDUM)
#endif      
            R = Z*RMAX
            FACTINT = RMAX
         else if(NRC.eq.10 ) then ! Provided from outside the code. 
            R = Rrecup
            P = Precup
            FACTINT = 1
         endif

         CALL xro(R)
         
         dkfp=(3.D0*DPI*DPI*DROP)**(1.D0/3.D0)
         dkfn=(3.D0*DPI*DPI*DRON)**(1.D0/3.D0)
         
c     Coulombian potential correction 
         
         dmp=938.27208d0/hbarc
         dmn=939.56542d0/hbarc
         
         if (ipol.eq.1) then 
            
            if(NRC.eq.0)then
               vc = vcd(ir)
            elseif(NRC.eq.1)then

c               alpha = 1.d0 / 137.036d0               
c               CALL DSG20R (0.D0,R,5,DRP,NRP)
c              do irp = 1,nrP
c                  rp=drp(irp)
c                  df0 (irp) = (rp**2)*densq(rp)/r
c               enddo                      
c               CALL DRG20R (0.D0,R,5,df0,f1)
c               CALL DSG20R (R,RMAX,5,DRP,NRP)
c               do irp = 1,nrP
c                 rp=drp(irp)                      
c                  df0 (irp) = rp*densq(rp)
c               enddo                  
c               CALL DRG20R (R,RMAX,5,df0,f2)
c               vc=-alpha*ieta*4.d0*dpi*(f1+f2)
c               PRINT*," TESTING ",R,VC

c               IMIN = 1
c               IMAX = NR
      
               call vcinterpolate(R,VCtemp)  
    
               vc=VCtemp*1.d3/hbarc

            elseif(NRC.eq.10)then
               vc=vc*1.d3/hbarc
!     calculated before
            endif
            
c     epro = dsqrt(dmp**2+dkfp**2)              
c     eneu = dsqrt(dmn**2+dkfn**2)              
c     qth = ieta*(epro-eneu)

         else 
            vc = 0
         endif                  ! Endif ipol.eq.1 
        
c         print*," 2 !!!!!!!!!! ",Q0*hbarc,qvaluecorr*1000.
 
c     q0 = q0+qth      
         eoutlocal = eout-vc
         
         if ((eout-vc-dmf).lt.0) then
            iflag=1             ! asi vale cero la respuesta
            poutlocal = 0.d0
         else
            poutlocal = dsqrt(eoutlocal**2-dmf**2)
         endif

         pout = dsqrt(eout**2-dmf**2)
         dq =DSQRT(PIN**2+POUTlocal**2-2.D0*coseno*PIN*POUTlocal) 
         fema=eoutlocal*poutlocal/(eout*pout)
         q2=(q0**2-dq**2)

*     CONSTANTS
         
         UYI=(0.D0,1.D0)
         DMAGP=2.792847D0
         DMAGN=-1.913043D0
         FFAC=0.843D0*1.D3/hbarc
         xln=5.6d0

         if (ipol.eq.0) xln = 0.d0            
         
         GQ=1.D0/(1.D0-Q2/(FFAC*FFAC))**2
         tau = -q2/4.d0/dmnu**2
         
ccccccccccccccccccxma = 1.049*1.D3/hbarc
         ga = 1.267d0           !1.257d0 !!!
         
         f1v = 0.5d0*gq*(1.d0+tau*dmagp
     f        -dmagn*xln*tau**2/(1.d0+xln*tau))/(1.d0+tau)
         xmuf2v =0.5d0*gq*(-1.d0+dmagp
     f        -dmagn*(1.d0+xln*tau+tau)/(1.d0+xln*tau))/(1.d0+tau)
         
C     GAQ = ga /(1.d0-Q2/xma**2)**2
         GAQ =(-1.)*ELAXFF(-q2/((1.d3/ hbarc)**2),.false.,ieta)
         
         GPQ = 2.d0*dmnu*gaq/(xmpi**2-q2)
         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (ipol.gt.0) then
            C0POL=380.D0/hbarc
            VT=VTP(Q0,DQ)
            VL=VLP(Q0,DQ)
         endif
****
****
         if(NPC.eq.0) then
            if (ieta.eq.1) then   
               xlind=DUR(q0,DQ,DKFn,dkfp,DMN,DMP,t0,r00)
            elseif (ieta.eq.-1) then
               xlind=DUR(q0,DQ,DKFp,dkfn,DMP,DMN,t0,r00)
            endif

         elseif(NPC.eq.1)then
       
            if (ieta.eq.1) then   
               xlind=DUR1(q0,DQ,DKFn,dkfp,DMN,DMP,t0,r00)
            elseif (ieta.eq.-1) then
               xlind=DUR1(q0,DQ,DKFp,dkfn,DMP,DMN,t0,r00)
            endif

         elseif(NPC.eq.2.or.NPC.eq.10)then 
            if (ieta.eq.1) then
               xlind=DURT(q0,DQ,DKFn,dkfp,DMN,DMP,t0,
     f              r00,NPC,P,CTHETAP,Panc)
            elseif (ieta.eq.-1) then
               xlind=DURT(q0,DQ,DKFp,dkfn,DMP,DMN,t0,
     f              r00,NPC,P,CTHETAP,Panc)
            else
               print*,"ERROR ieta =", ieta , " invalid " 
            endif

         endif
         
         if (ipol.gt.0) then
            if (iflag.eq.1) then
               xlind=0.d0
               t0 = 0.d0
               r00 = 0.d0
            endif
         endif
         
         tulin(0) = t0
         tulin(1) = 0.d0
         tulin(2) = 0.d0
         tulin(3) = (0.5d0*q2 + q0*tulin(0))/dq 
         
         rulin(0,0) = r00
         rulin(0,1) = 0.d0
         rulin(0,2) = 0.d0
         rulin(0,3) = (0.5d0*q2*tulin(0) + q0*rulin(0,0))/dq 
         
         rulin(1,0) = rulin(0,1)
         rulin(2,0) = rulin(0,2)
         rulin(3,0) = rulin(0,3)
         
         trazaerre =r00-dmnu**2
         errezz = (q2**2+4.d0*r00*q0**2+4.d0*q2*q0*t0)/(4.d0*dq**2)
         gamma = (trazaerre-errezz)/2.d0
         delta = (-trazaerre+3.d0*errezz)/2.d0/dq**2
         
         rulin(1,1) = gamma
         rulin(1,2) = 0.d0
         rulin(1,3) = 0.d0
         rulin(2,2) = gamma
         rulin(2,3) = 0.d0
         rulin(3,3) = gamma+delta*dq*dq
        
         do i = 2,3
            do j=1,i-1
               rulin(i,j) = rulin(j,i)
            enddo
         enddo  
         
         
         if (ipol.gt.0) then
            if (iflag.eq.1) then
               xlind=0.d0
               t0 = 0.d0
               r00 = 0.d0
            endif           
         endif  
         
         
         fact=1.d0
         facl=1.d0
         f00=1.d0
         
         if (ipol.gt.0) then
            
            DKF=(1.5D0*DPI*DPI*DRO)**(1.D0/3.D0)
               
            irel= 1             ! relativista  
            
            
            if (abs(xlind).gt.1.d-10) then  
               unuc_j = ulinrel(q0,dq,dkf)     
               u_zero=q0+uyi*0.0d0
               udel=delta_lind(u_zero,dq,dro,dkf)
            else                ! no es que valga cero, pero como xlind es cero da igual
               unuc_j = 0.d0+uyi*0.d0
               udel = 0.d0 +uyi*0.d0
            endif

            cunuc= unuc_j
            ufuli = unuc_j+udel          
            DELTA=dreal(ufuli-cunuc)
            
            FACT=(1.D0/abs(1.d0-VT*UFULI))**2
            FACL=(1.D0/abs(1.d0-VL*UFULI))**2
            delfact=(1.d0-vt*delta)**2
            delfacl=(1.d0-vl*delta)**2

c     fact=fact*delfact  
c     facl=facl*delfacl  
            
            FPRIMA0=(DRO/dro0*fp0in+(1.D0-DRO/dro0)*fp0ex)*C0POL
            F00=ABS(1.D0/(1.D0-FPRIMA0*CUNUC))**2 

         
         endif
         
         call amunu_pol(axx,azz,a0z,a00,axy,tulin,rulin,q2,q0,dq,f1v,
     f        xmuf2v,gaq,gpq,fact,facl,f00,ipol)
         
         R2 = R*R

         if(NRC .EQ. 0)then
            DFXX(IR)=cte_r*R2*axx*xlind*fema
            DFZZ(IR)=cte_r*R2*azz*xlind*fema
            DF00(IR)=cte_r*R2*a00*xlind*fema
            DFXY(IR)=cte_r*R2*axy*xlind*fema
            DF0Z(IR)=cte_r*R2*a0Z*xlind*fema
            
c
c     q0=q0old
            q2=q2old
            dq =dqold
         elseif(NRC.EQ. 1.OR.NRC.EQ.10)then
            xlind=xlind*FACTINT
            RESXX=cte_r*R2*axx*xlind*fema
            RESZZ=cte_r*R2*azz*xlind*fema
            RES00=cte_r*R2*a00*xlind*fema
            RESXY=cte_r*R2*axy*xlind*fema
            RES0Z=cte_r*R2*a0Z*xlind*fema
c     q0=q0old
            q2=q2old
            dq =dqold
            Rrecup=R
            Precup=P
            
c     exit
         endif
         
      END DO                    ! End of the loop on the radius. 
      
      if(NRC .EQ. 0)then
         CALL DRG20R (0.D0,RMAX,MNR,DFXX,RESXX)
         CALL DRG20R (0.D0,RMAX,MNR,DFZZ,RESZZ)
         CALL DRG20R (0.D0,RMAX,MNR,DF00,RES00)
         CALL DRG20R (0.D0,RMAX,MNR,DFXY,RESXY)
         CALL DRG20R (0.D0,RMAX,MNR,DF0Z,RES0Z)
      endif
  
      wd(1) = resxx/2.d0/dma 
      wd(4) = dma*(reszz-resxx)/2.d0/dq**2
      wd(5) = (res0z-q0*(reszz-resxx)/dq)/dq
      wd(3) = resxy/dq
      wd(2) = (res00+resxx+(reszz-resxx)*(q0/dq)**2
     f     -2.d0*q0*res0z/dq)/2.d0/dma


c         print*," 4 !!!!!!!!!! ",Q0*hbarc,qvaluecorr*1000.
      
      RETURN
      END




      SUBROUTINE WFREE(ILIN,Q0,dq,Q2,tensLA,ieta,S,dmnui)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      
      COMMON/datos/dpi,hbarc,GF0,DMNU,DMA
c     COMMON/datos3/xmpi,gp
        
      COMMON/datos2/dmneutrino,dmlepton,dmi,dmf,
     f     coscabibbo,dmuon,dmelectron,dmtau,xuma
      COMMON/mavalue/xma
      
      ct_1=Q2-dmf**2
      ct_2=S-dmnui**2
      
*     CONSTANTS
      DMAGP=2.792847D0
      DMAGN=-1.913043D0
      FFAC=0.843D0*1.D3/hbarc
      xln=5.6d0  
      
      GQ=1.D0/(1.D0-Q2/(FFAC*FFAC))**2
      tau = -q2/4.d0/dmnui**2
ccccccccccccccccccxma = 1.049*1.D3/hbarc
      ga = 1.267d0              !1.2768d0!1.2742d0!!1.2755d0!1.267d0!1.257d0
      xmpi = 139.57d0 / hbarc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      f1v = 0.5d0*gq*(1.d0+tau*dmagp
     f     -dmagn*xln*tau**2/(1.d0+xln*tau))/(1.d0+tau)
      xmuf2v =0.5d0*gq*(-1.d0+dmagp
     f     -dmagn*(1.d0+xln*tau+tau)/(1.d0+xln*tau))/(1.d0+tau)
      
C      GAQ = ga /(1.d0-Q2/xma**2)**2
      GAQ =(-1.)*ELAXFF(-q2/((1.d3/ hbarc)**2),.false.,ieta)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      A1=8.d0*q2*((f1v+xmuf2v)**2+GAQ**2*(1.d0/4.d0-dmnui**2/q2))
      A2=32.d0*f1v**2-8.d0*(xmuf2v**2)*q2/dmnui**2+8.d0*GAQ**2
      A3=16.d0*GAQ*(f1v+xmuf2v)
      A4=-8.d0*q2/(dmnui**2)*(dmnui**2/q2+1.d0/4.d0)*xmuf2v**2
     f     -16.d0*f1v*xmuf2v
     f     -8.d0*dmnui**2*GAQ**2/(xmpi**2-q2)*(q2/(xmpi**2-q2)+2.d0) 
      
      V1=A1+s*A2/2.d0-ieta*q2*A3/2.d0-A4/2*dmf**2
      V2=(s-dmnui**2)*A2/2.d0-ieta*A3*q2
      
      tensLA=ct_1*V1+ct_2*V2
      
      RETURN
      END
      
      SUBROUTINE amunu_pol(axx,azz,a0z,a00,axy,tulin,rulin,q2,q0
     f     ,dq,f1v,xmuf2v,gaq,gpq,fact,facl,f00,ipol)
      
      
      IMPLICIT DOUBLE PRECISION (A-H,O-T,V-Z)
      DIMENSION tulin(0:3),rulin(0:3,0:3)
      COMMON/datos/dpi,hbarc,GF0,DMNU,DMA
      
        
      axx = 16.d0*f1v**2*(2.d0*rulin(1,1)-q2/2.d0)+ 
     f     2.d0*q2*xmuf2v**2*(-4.d0*fact-4.d0*rulin(1,1)/dmnu**2) +
     f     4.d0*gaq**2*(2.d0*rulin(1,1)-(q2/2.d0-2.d0*fact*dmnu**2))-
     f     16.d0*f1v*xmuf2v*fact*q2
      
      azz = 16.d0*f1v**2*(2.d0*rulin(3,3)+2.d0*dq*tulin(3)
     f     -q2/2.d0)+ 
     f     2.d0*q2*xmuf2v**2*(-4.d0-4.d0*rulin(3,3)/dmnu**2-
     f     4.d0*dq*tulin(3)/dmnu**2-dq**2*(4.d0/q2+1.d0/dmnu**2)) +
     f     4.d0*gaq**2*(2.d0*rulin(3,3)+2.d0*dq*tulin(3)
     f     -(q2/2.d0-2.d0*facl*dmnu**2))-
     f     (2.d0*facl*gpq**2*q2+8.d0*gaq*gpq*facl*dmnu)*dq**2-
     f     16.d0*f1v*xmuf2v*(q2+dq**2)
      
      
      a0z = 16.d0*f1v**2*((2.d0*rulin(0,3)+tulin(0)*dq)*f00
     f     +tulin(3)*q0)+
     f     2.d0*q2*xmuf2v**2*(-4.d0*rulin(0,3)/dmnu**2-
     f     2.d0*(dq*tulin(0)+q0*tulin(3))/dmnu**2
     f     -dq*q0*(4.d0/q2+1.d0/dmnu**2)) +
     f     4.d0*gaq**2*((2.d0*rulin(0,3)+dq*tulin(0))*facl+q0*tulin(3))-
     f     (2.d0*facl*gpq**2*q2+8.d0*gaq*gpq*facl*dmnu)*dq*q0-
c     c     f    16.d0*f1v*xmuf2v*f00*dq*q0
     f     16.d0*f1v*xmuf2v*dq*q0

c     c              axy = 16.d0*gaq*(xmuf2v+f1v)*fact*(-dq*tulin(0)
c     c     f         +q0*tulin(3))
      
      axy = 16.d0*gaq*(xmuf2v+f1v)*(-dq*tulin(0)*fact
     f     +q0*tulin(3))

      a00 = 16.d0*f1v**2*(2.d0*rulin(0,0)*f00+2.d0*q0*tulin(0)
     f     +q2/2.d0)+ 
     f     2.d0*q2*xmuf2v**2*(4.d0-4.d0*rulin(0,0)/dmnu**2-
     f     4.d0*q0*tulin(0)/dmnu**2-q0**2*(4.d0/q2+1.d0/dmnu**2)) +
     f     4.d0*gaq**2*(2.d0*rulin(0,0)+2.d0*q0*tulin(0)
     f     +(q2/2.d0-2.d0*dmnu**2))-
     f     (2.d0*facl*gpq**2*q2+8.d0*gaq*gpq*facl*dmnu)*q0**2-
     f     16.d0*f1v*xmuf2v*(-q2+q0**2)*f00
      
      RETURN 
      END
      
      
      
      SUBROUTINE XRO(DR)
      IMPLICIT REAL*8 (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      DIMENSION DID(2000),DFD(2000)
      COMMON/VFPI/DROP,DRON
      COMMON/ROFEROS/KLAVE
      COMMON/NXRO/NXRO
      COMMON/NUC/DZZ,DAA
      COMMON/DENSIDAD/dro,dro0
C************************************************************** 
C     DATA FOR PROTON AND NEUTRON DENSITY FUNCTION CALCULATION 
      COMMON/CTEROPN/DXP,DA0P,DRO0P,DXN,DA0N,DRO0N
C     DROPR...  proton density /dro0p
C     DRONR.... neutron density /dro0n
C     KLAVE=0 Fermi density . KLAVE=1 oscillator density 
C***************************************************************
      DATA DPI/3.141592653589793D0/
      
      IF(NXRO.EQ.0) THEN

         IF (KLAVE.EQ.0) THEN  ! FERMI 
C
C Compute the normalization constants for neutron and proton densities. 
C They are normalized to DZZ and DAA-DZZ. 
C
            DRMAX1P=5.D0*DXP
            CALL DSG20R(0.D0,DRMAX1P,10,DID,NPP)
            DO I=1,NPP
               DI=DID(I)
               DFD(I)=DI*DI/(DEXP((DI-DXP)/DA0P)+1.D0)
            END DO
            CALL DRG20R(0.D0,DRMAX1P,10,DFD,DRESP)
            DRO0P=DZZ/(DRESP*4.D0*DPI)
            
            DRMAX1N=5.D0*DXN
            CALL DSG20R(0.D0,DRMAX1N,10,DID,NPP)
            DO I=1,NPP
               DI=DID(I)
               DFD(I)=DI*DI/(DEXP((DI-DXN)/DA0N)+1.D0)
            END DO
            CALL DRG20R(0.D0,DRMAX1N,10,DFD,DRESN)
            DRO0N=(DAA-DZZ)/(DRESN*4.D0*DPI)
            
         ELSE                   ! oscillator
            DRMAX1P=8.D0*DXP
            CALL DSG20R(0.D0,DRMAX1P,10,DID,NPP)
            DO I=1,NPP
               DI=DID(I)
               DFD(I)=DI*DI*(1.D0+DA0P*(DI/DXP)**2)*DEXP(-(DI/DXP)**2)
            END DO
            CALL DRG20R(0.D0,DRMAX1P,10,DFD,DRESP)
            DRO0P=DZZ/(DRESP*4.D0*DPI)
            DRMAX1N=8.D0*DXN
            CALL DSG20R(0.D0,DRMAX1N,10,DID,NPP)
            DO I=1,NPP
               DI=DID(I)
               DFD(I)=DI*DI*(1.D0+DA0N*(DI/DXN)**2)*DEXP(-(DI/DXN)**2)
            END DO
            CALL DRG20R(0.D0,DRMAX1N,10,DFD,DRESN)
            DRO0N=(DAA-DZZ)/(DRESN*4.D0*DPI)
         ENDIF
        
         NXRO=1
         
      END IF
      
      IF(DR.GT.5.D0*DXN.AND.KLAVE.EQ.0)THEN
         DROPR=0.D0
         DRONR=0.D0
         DROP=0.D0
         DRON=0.D0
         dro  = drop+dron
         dro0 = dro0p+dro0n
         RETURN
      END IF
      
      IF(DR.GT.8.D0*DXN.AND.KLAVE.EQ.1)THEN
         DROPR=0.D0
         DRONR=0.D0
         DROP=0.D0
         DRON=0.D0
         dro  = drop+dron
         dro0 = dro0p+dro0n
         
         RETURN
      END IF
      
      
      IF(KLAVE.EQ.0)THEN
C     Fermi reduced density 
         DROPR=1.D0/(DEXP((DR-DXP)/DA0P)+1.D0)
         DRONR=1.D0/(DEXP((DR-DXN)/DA0N)+1.D0)
         DROP=DRO0P*DROPR
         DRON=DRO0N*DRONR
         
      ELSE
C     Oscillator reduced density 
         DROPR=(1.D0+DA0P*(DR/DXP)**2)*DEXP(-(DR/DXP)**2)
         DRONR=(1.D0+DA0N*(DR/DXN)**2)*DEXP(-(DR/DXN)**2)
         DROP=DRO0P*DROPR
         DRON=DRO0N*DRONR
      END IF
      dro  = drop+dron
      dro0 = dro0p+dro0n
      
      RETURN
      END
      
      
      double precision function densq(DR)
      IMPLICIT REAL*8 (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      DIMENSION DID(2000),DFD(2000)
      
      COMMON/ROFEROS/KLAVE
      COMMON/NCXRO/NCXRO
      COMMON/CTEROPNC/DXP,DA0P,DRO0P
      COMMON/NUC/DZZ,DAA
      
C     DROPR... proton density /dro0p
C     KLAVE=0 Fermi density . KLAVE=1 Oscillator density 
      DATA DPI/3.141592653589793D0/
      
      IF(NCXRO.EQ.0) THEN
         IF (KLAVE.EQ.0) THEN   ! Femi density 
C     
C     Proton density normalization constants normalized to DZZ.
C     
            DRMAX1P=5.D0*DXP
            CALL DSG20R(0.D0,DRMAX1P,10,DID,NPP)
            DO I=1,NPP
               DI=DID(I)
               DFD(I)=DI*DI/(DEXP((DI-DXP)/DA0P)+1.D0)
            END DO
            CALL DRG20R(0.D0,DRMAX1P,10,DFD,DRESP)
            DRO0P=DZZ/(DRESP*4.D0*DPI)
            
         ELSE                   ! oscillator density 
            DRMAX1P=8.D0*DXP
            CALL DSG20R(0.D0,DRMAX1P,10,DID,NPP)
            DO I=1,NPP
               DI=DID(I)
               DFD(I)=DI*DI*(1.D0+DA0P*(DI/DXP)**2)*DEXP(-(DI/DXP)**2)
            END DO
            CALL DRG20R(0.D0,DRMAX1P,10,DFD,DRESP)
            DRO0P=DZZ/(DRESP*4.D0*DPI)
            
         ENDIF
         
         NCXRO = 1
         
      END IF
      
      IF(DR.GT.5.D0*DXP.AND.KLAVE.EQ.0)THEN
         densq=0.D0
         RETURN
      END IF
      
      IF(DR.GT.8.D0*DXP.AND.KLAVE.EQ.1)THEN
         densq=0.D0
         RETURN
      END IF
      
      
      IF(KLAVE.EQ.0)THEN
C     Reduced Fermi density 
         DROPR=1.D0/(DEXP((DR-DXP)/DA0P)+1.D0)
         densq=DRO0P*DROPR
         
         
      ELSE
C     Reduced oscillation density 
         DROPR=(1.D0+DA0P*(DR/DXP)**2)*DEXP(-(DR/DXP)**2)
         densq=DRO0P*DROPR
         
      END IF
      
      RETURN
      END
      
      
      SUBROUTINE convolucion (dncxp,dnca0p,dncxn,dnca0n,klave,lc,
     f     dxp,da0p,dxn,da0n)
      IMPLICIT REAL*8 (A-H,O-T,V-Z)
      DATA DPI/3.141592653589793D0/
      
      PI2 = DPI*DPI
*
*     Parameters deconvolution 
*     
      if (lc.eq.0) then
         
         IF (KLAVE.EQ.0) THEN
            DXP=DNCXP+3.45D0*DNCXP/(15.*DNCXP**2+7*PI2*DNCA0P**2)
            DA0P=DSQRT((DNCXP**3+PI2*DXP*DNCA0P**2-DXP**3)/(DXP*PI2))
            DXN=DNCXN+3.45D0*DNCXN/(15.*DNCXN**2+7*PI2*DNCA0N**2)
            DA0N=DSQRT((DNCXN**3+PI2*DXN*DNCA0N**2-DXN**3)/(DXN*PI2))
         ELSE
            DXPP=0.46D0
*
*     DXPP IS EQUAL TO 2/3 OF THE AVERAGE QUADRATIC RADIOUS OF THE PROTOM IN FM2
*
            DXP=DSQRT(DNCXP**2-DXPP)
            DXXP=DNCA0P*(DNCXP**2)/((1.D0+1.5D0*DNCA0P)*DXP**2)
            DA0P=2.D0*DXXP/(2.D0-3.D0*DXXP)
            DXN=DSQRT(DNCXN**2-DXPP)
            DXXN=DNCA0N*(DNCXN**2)/((1.D0+1.5D0*DNCA0N)*DXN**2)
            DA0N=2.D0*DXXN/(2.D0-3.D0*DXXN)
            
         ENDIF
         
         
      else                      !lc.ne.0
         
         dxp=dncxp
         da0p=dnca0p
         dxn=dncxn
         da0n=dnca0n
         
      endif
      
      RETURN
      end
      
      
      

      FUNCTION DUR (Q0,DQ,PF1,PF2,DMNU,DMSO,t0,r00)
      IMPLICIT DOUBLE PRECISION (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)

c     IF (Q0.LT.0) THEN
c     WRITE(6,*)'EN LA LLAMADA A DUR,Q0 ES NEGATIVA'
c     DUR=0.
c     t0=0.d0
c     r00 =0.d0 
c     RETURN
c     ENDIF
      DPI=3.141592653589793D0
      CTE=-0.5D0*DMNU*DMNU/DPI/DQ
      Q2=Q0**2-DQ**2
      IF (Q2.gE.0) THEN
         DUR=0.
         t0=0.d0
         r00 =0.d0 
         RETURN
      endif
      EF1=DSQRT(DMNU**2+PF1**2)
      EF2=DSQRT(DMSO**2+PF2**2)
      IF ((EF2-Q0).GT.EF1) THEN
         DUR=0.
         t0=0.d0
         r00 =0.d0
       RETURN
      ENDIF
      A=0.5D0*(-Q0+DQ*SQRT(1.d0-(DMNU+DMSO)**2/Q2))
      IF (A.GT.EF1) THEN
         DUR=0.
         t0=0.d0
         r00 =0.d0
         RETURN
      ENDIF

      EINF=MAX(DMNU,EF2-Q0,a)
c      write(70,*)Q0,DQ,EF1,PF1
      t0  = 0.5d0*(ef1+einf)
      r00 = (ef1**2+einf**2+ef1*einf)/3.d0 
      DUR=(ef1-einf)*CTE
      RETURN
      END
      
      FUNCTION DUR1 (Q0,DQ,PF1,PF2,DMNU,dmso,t0,r00)
      IMPLICIT DOUBLE PRECISION (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      DIMENSION pn(100000)
c
ccc   PF1=KFn PF2=KFp
c
      DPI=3.141592653589793D0
      CTE=-DMNU*DMNU/(2.d0*DPI) 

      Q2=Q0**2-DQ**2
      IF (Q2.gE.0) THEN
         DUR=0.
         t0=0.d0
         r00 =0.d0 
         RETURN
      endif

c      DpminCalc=targetminmon(q0,dq,dmnu,dmso)
      DpminCalc = 0.

      if( DpminCalc.gt.PF1) then
         P=DpminCalc
         DUR1=0.d0
         t0=0.d0
         r00=0.d0
         RETURN 
      endif 

      DUR1=0.d0
      t0=0.d0
      r00=0.d0  

      CALL DSG20R(DpminCalc,PF1,500,pn,Np) 
      
      Deltap=(PF1-DpminCalc)/Np
      
      DO I=1,NP
         dP=pn(I)         
         
         Ep=DSQRT(dmnu**2+dP**2)     
c         
c     Epq = DSQRT(dmnu**2+pqNormsq)
c     delta(f(x))=delta(x)/abs(df/dx)
c     df=dp*dq*Epq
c
         CTHETA=0.5d0*((Ep+Q0)**2-(dmnu**2+dP**2+DQ**2))/(dP*DQ)
         pqNormsq=dP**2+DQ**2+2.d0*dP*DQ*CTHETA
         pqNorm=DSQRT(pqNormsq)
         Epq=Ep+Q0

         if((CTHETA.lt.-1).or.(CTHETA.gt.1))then
            DUR1=DUR1+0.d0
            t0=t0+0.d0
            r00=r00+0.d0
         else if(dP.gt.PF1)then
            DUR1=DUR1+0.d0
            t0=t0+0.d0
            r00=r00+0.d0
         else if(PF2.gt.pqNorm)then
            DUR1=DUR1+0.d0
            t0=t0+0.d0
            r00=r00+0.d0            
         else
            DUR1=DUR1+dP/Ep/dQ
            t0=t0+dP/dQ
            r00=r00+Ep*dP/dQ
         endif
      ENDDO
      
      if(DUR1.ne.0.d0)then
         t0=t0/DUR1
         r00=r00/DUR1
         DUR1=DUR1*CTE*Deltap
         RETURN
      else
         t0=t0
         r00=r00
         DUR1=DUR1*CTE*Deltap
         RETURN
      endif
      
      RETURN
      END


      FUNCTION DURT(Q0,DQ,PF1,PF2,DMNU,DMSO,t0,r00,NPC,P,CTHETAP,Panc)
      IMPLICIT DOUBLE PRECISION (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      COMMON/datos/dpi,hbarc,GF0,D,DMA

c      DIMENSION pn(100000)
c
ccc   PF1=KFn PF2=KFp
c
      DPI=3.141592653589793D0
      CTE=-DMNU*DMNU/(2.d0*DPI) 
     
      Q2=Q0**2-DQ**2
      IF (Q2.GE.0) THEN
         DURT=0.
         t0=0.d0
         r00 =0.d0 
         RETURN
      endif     
C
      if( NPC.EQ.2 ) THEN    

c     DpminCalc=targetminmon(q0,dq,dmnu,DMSO)
         DpminCalc=0.
         
c     DpminCalc = 0.    
         if( DpminCalc.gt.PF1 ) then ! Minimum momentum should be below the fermi level. 
            P=DpminCalc
            DURT=0.d0
            t0=0.d0
            r00=0.d0
            RETURN 
         endif          
     
#ifdef __GFORTRAN__
            CALL RANDOM_NUMBER(Z)
#else
            Z = RLU(IDUM)
#endif      
         Dp = DpminCalc+Z*(PF1-DpminCalc) ! Generate flat distribution between min and fermi level.
         FACTINT = (PF1-DpminCalc)
      ELSE IF ( NPC.EQ.10 ) THEN
         IF( P > PF1 ) THEN 
            DURT=0.d0
            t0=0.d0
            r00=0.d0
            PRINT*,"ERROR, wrong target momentum."
            RETURN
         ENDIF
         FACTINT = 1.
         dP = P
      ELSE 
         FACTINT = 1.
      ENDIF 
      
      EF1=DSQRT(DMNU**2+PF1**2)
      EF2=DSQRT(DMSO**2+PF2**2)
      IF((EF2-Q0).GT.EF1) THEN
         DURT=0.
         t0=0.d0
         r00 =0.d0 
         RETURN
      ENDIF
      A=0.5d0*(-Q0+DQ*SQRT(1.d0-(DMNU+DMSO)**2/Q2))      
      IF(A.GT.EF1) THEN
         DURT=0.
         t0=0.d0
         r00 =0.d0
c         PRINT*,"THIS HAPPENS",A*hbarc,EF1*hbarc
         RETURN
      ENDIF

C     ENKin = 0.5*Dp*dp/DMA 

      Ep=DSQRT(dmnu**2+dP**2)

c     Q0 = Q0-ENKin
      
      CTHETA=0.5d0*((Ep+Q0)*(Ep+Q0)-(dmso*dmso+dP*dP+DQ*DQ))/(dP*DQ)
      pqNormsq=dP**2+DQ**2+2.d0*dP*DQ*CTHETA
      pqNorm = DSQRT(pqNormsq)
      Epq=Ep+Q0
      
      if(CTHETA.lt.-1.or.ctheta.gt.1)then
         DURT=0.d0
         t0=0.d0
         r00=0.d0
      else if(dP.gt.PF1)then
         DURT=0.d0
         t0=0.d0
         r00=0.d0
      else if(PF2.gt.pqNorm)then
         DURT=0.d0
         t0=0.d0
         r00=0.d0
      else
         DURT=dP/Ep/dQ
         t0=dP/dQ
         r00=Ep*dP/dQ
c         if( dQ < 1d-1 ) Print*,dQ
      endif
 
      if(DURT.gt.0.d0)then
         Panc=pqNorm
         CTHETAP=CTHETA
         P=dP
         t0=t0/DURT
         r00=r00/DURT
         DURT=DURT*CTE*FACTINT
         RETURN
      else
         Panc=pqNorm
         CTHETAP=CTHETA
         P=dP
         DURT=0.d0
         t0=0.d0
         r00=0.d0
         RETURN
      endif
      
      END
c
c
      
      FUNCTION VTP (A0,A1)
      IMPLICIT REAL*8 (A-H,O-T,V-Z)
      COMMON/datos/dpi,hbarc,GF0,DMNU,DMA
      COMMON/mavalue/xma

      INTEGER*4 irelativistic
      COMMON/rpavariables/fp0in,fp0ex,f,fstar,pilambda,
     f     cr0,rholambda,gp,xmpi,xmrho,irelativistic
   
      DFM2=f/(XMPI*XMPI)
      A2=A0**2-A1**2
c      CRO=2.0
      XL2=(rholambda/hbarc)**2
      DMRO2=XMRHO*XMRHO
      AFORMA=cr0*((XL2-DMRO2)/(XL2-A2))**2

      if(irelativistic.eq.1)then
         VTP=DFM2*(-AFORMA*A2/(A2-DMRO2)+GP)
      else
         VTP=DFM2*(AFORMA*A1**2/(A2-DMRO2)+GP)
      endif
      RETURN
      END


      
c      FUNCTION VTP (A0,A1)
c      IMPLICIT REAL*8 (A-H,O-T,V-Z)
c      COMMON/datos/dpi,hbarc,GF0,DMNU,DMA
c      COMMON/datos3/xmpi,gp
c       
c      DFM2=0.08*4.D0*DPI/(XMPI*XMPI)
c      A2=A0**2-A1**2
c      XL2=(2500.d0/hbarc)**2
cc     CRO=2.0
c      DMRO2=(770.d0/hbarc)**2
c      AFORMA=CRO*((XL2-DMRO2)/(XL2-A2))**2
c      VTP=DFM2*(AFORMA*A1**2/(A2-DMRO2)+GP)
c      RETURN
c      END

      FUNCTION VLP (A0,A1)
      IMPLICIT REAL*8 (A-H,O-T,V-Z)
      COMMON/datos/dpi,hbarc,GF0,DMNU,DMA

      INTEGER*4 irelativistic
      COMMON/rpavariables/fp0in,fp0ex,f,fstar,pilambda,
     f     cr0,rholambda,gp,xmpi,xmrho,irelativistic

      DMPI2=xmpi*xmpi 
      DFM2= f/DMPI2
      A2=A0**2-A1**2
      XL2=(pilambda/hbarc)**2

      AFORMA=((XL2-DMPI2)/(XL2-A2))**2
      if(irelativistic.eq.1)then
         VLP=DFM2*(-AFORMA*A2/(A2-DMPI2)+GP)  
      else
         VLP=DFM2*(AFORMA*A1**2/(A2-DMPI2)+GP)
      endif
      RETURN
      END



c      FUNCTION VLP (A0,A1)
c      IMPLICIT REAL*8 (A-H,O-T,V-Z)
c      COMMON/datos/dpi,hbarc,GF0,DMNU,DMA
c      COMMON/datos3/xmpi,gp
c      COMMON/rpavariables/fp0in,fp0ex,f,fstar,pilambda,
c     f     cr0,rholambda,gp,xmpi,xmrho

c      DMPI2=xmpi*xmpi 
c      DFM2=0.08*4.D0*DPI/DMPI2
c      A2=A0**2-A1**2
c      XL2=(1200.d0/hbarc)**2      
c      AFORMA=((XL2-DMPI2)/(XL2-A2))**2
c      VLP=DFM2*(AFORMA*A1**2/(A2-DMPI2)+GP)
c      RETURN
c      END


*********************************************************************
*      FUNCIONES DE LINHARD
**********************************************************************
       SUBROUTINE ULIND(QZR,Q,XKF,CUFUN)
              implicit real*8 (a-b,d-h,o-z)
       IMPLICIT complex*16 (C)
       COMMON/cunuc/cunuc
       XMN=6.7348
       XMD=8.8315
       WRES=XMd-XMN
       FNS=1.
       FDS=4.52
       CYI=(0.,1.)
       PI=3.141592
       QA=Q/XKF
       QZ=ABS(QZR)
       QZA=ABS(QZ*XMN/XKF**2)
       IF(Q.EQ.0.) GO TO 15
       AM=QZA/QA-QA/2.
       AP=QZA/QA+QA/2.
       IF(ABS(1.-AM).LT.0.00001) AM=1.00001
       IF(ABS(1.-AP).LT.0.00001) AP=1.00001
       TERF=ABS((1.+AM)/(1.-AM))+1.e-15
       TERS=ABS((1.+AP)/(1.-AP))+1.e-15
       IF(QZA-0.00001) 16,16,18
 18      DCOCI=ABS(QA/QZA)
       IF(DCOCI-0.1) 15,15,16
 15      RUN=FNS*XMN*XKF/PI**2*2./3.*(QA/QZA)**2
       GO TO 17
 16      RUN=FNS*XMN*XKF/PI**2*(-1.+(1.-AM**2)/(2.
     F *QA)*log(TERF)-(1.-AP**2)/(2.*QA)*log(TERS))
 17      YMN=0.
       IF((QA**2/2.+QA).GE.QZA.AND.QZA.GE.(QA**2/2.-QA).AND.QA.
     F GE.2.)YMN=-2.*XMN*XKF*FNS/(4.*PI*QA)*(1.-AM**2)
       IF(QA.LT.2..AND.(QA+QA**2/2.).GE.QZA.AND.QZA.GE.(QA-QA**
     F 2/2.))YMN=-2.*XMN*XKF*FNS/(4.*PI*QA)*(1.-AM**2)
       IF(QA.LT.2..AND.0..LE.QZA.AND.QZA.LE.(QA-QA**2/2.))
     F YMN=-2.*XMN*XKF*FNS/(4.*PI*QA)*2.*QZA
       CUNUC=CMPLX(RUN,YMN)
       B=Q/XMD
       RO=XKF**3*2./(3.*PI**2)
       IF(Q.EQ.0.) GOTO 250
       FAC=(4./3.*XKF/(2.*PI))**2/B**3*FDS
 250      IF(ABS(QZ)-1.) 2,2,3
 2      A=(QZ-WRES-Q**2/(2.*XMD))/XKF
       AP=(-QZ-WRES-Q**2/(2.*XMD))/XKF
       IF(ABS(B/A)-0.1) 25,25,26
 25      TDIR=4./9.*FDS*RO/(A*XKF)
       GO TO 27
 26      TERF=ABS((A+B)/(A-B))+1.e-15
       TDIR=FAC*(B*A+(B**2-A**2)/2.*log(TERF))
 27      IF(ABS(B/AP)-0.1) 28,28,29
 28      TCROS=4./9.*FDS*RO/(AP*XKF)
       GO TO 30
 29      TERS=ABS((AP-B)/(AP+B))+1.e-15
       TCROS=FAC*(B*AP-(B**2-AP**2)/2.*log(TERS))
 30      CUDEL=TDIR+TCROS
       GO TO 10
 3      SQS=XMN+ABS(QZ)
       QFR=SQRT((QZ**2-1.)/(1.+QZ/XMN))
       GAMH=1./(3.*4.*PI)*FDS*XMN/SQS*QFR**3
       IF(QZ.GT.0.) CAC=(QZ-WRES+CYI*GAMH-Q**2/(2.*XMD))/XKF
       IF(QZ.GT.0.)CAPC=(-QZ-WRES-Q**2/(2.*XMD))/XKF
       IF(QZ.LE.0.) CAC=(QZ-WRES-Q**2/(2.*XMD))/XKF
       IF(QZ.LE.0.) CAPC=(-QZ-WRES+CYI*GAMH-Q**2/(2.*XMD))/XKF
       IF(abs(B/CAC)-0.1) 35,35,36
 35      CDIR=4./9.*FDS*RO/(CAC*XKF)
       GO TO 37
 36      CDIR=FAC*(B*CAC+(B**2-CAC**2)/2.*log((CAC+B)/(CAC-B)))
 37      IF(abs(B/CAPC)-0.1) 38,38,39
 38      CCROS=4./9.*FDS*RO/(CAPC*XKF)
       GO TO 40
 39      CCROS=FAC*(B*CAPC-(B**2-CAPC**2)/2.*log((CAPC-B)/(CAPC+B)))
 40      CUDEL=CDIR+CCROS
 10      CONTINUE
       CUFUN=CUNUC+CUDEL
       RETURN
       END





        function delta_lind (q_zero,q_mod,rho,k_fermi)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   complex Lindhard function for symmetric nuclear matter:
c                    from Appendix of
c                    E.Oset et al Phys. Rept. 188:79, 1990
c                    formula A.4 
c
c   input variables: 
c     q_zero [fm^-1] : Energy
c     q_mod  [fm^-1] : Momentum
c     rho    [fm^3]  : Nuclear density
c     k_fermi[fm^-1] : Fermi momentum
c
c   All variables are real*8
c
c   output variable: 
c     delta_lind [fm^-2]
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            ATTENTION!!!
c Only works properly for real q_zero,
c if q_zero has an imaginary part calculates the L. function
c assuming Gamma= 0.
c Therefore this SUBROUTINE provides two different functions
c depending on whether q_zero is real or not!!!!!!!!!!!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        IMPLICIT NONE
        real*8 q_mod,rho,k_fermi,m,rq_zero,gamma,s,srot,wr
        real*8 fdel_f,mpi,gammap,sp,srotp,qcm,qcmp,md,pi
        complex*16 delta_lind,q_zero,z,zp,i,pzeta,pzetap 
        real*8 fp0in,fp0ex,f,fstar,pilambda,
     f         cr0,rholambda,gp,xmpi,xmrho

c
c m = 939/197.3,md = 1232./179.3, mpi = 139./197.3 
c
        DATA md,m,i/6.2433d0,4.7592d0,(0.d0,1.d0)/
        DATA mpi,pi/0.7045d0,3.14159265358979323846264338d0/
 
        INTEGER*4 irelativistic
        COMMON/rpavariables/fp0in,fp0ex,f,fstar,pilambda,
     f         cr0,rholambda,gp,xmpi,xmrho,irelativistic


c
c  f*/f = 2.13 --> f*^2/4pi = .36
c  NOTE: fdel = fdel_f because f =1.0
c

        fdel_f = fstar
        
        wr = md-m
        
        
        gamma = 0.d0 
        gammap = 0.d0
        if(imag(q_zero)**2.lt.1.d-36)then
          rq_zero = dreal(q_zero)
          s = m**2+rq_zero**2-q_mod**2+
     &        2.d0*rq_zero *sqrt(m**2+3.d0/5.d0*k_fermi**2)
          if(s.gt.(m+mpi)**2)then
           srot = sqrt(s)
           qcm = sqrt(s**2+mpi**4+m**4-2.d0*(s*mpi**2+s*m**2+
     &    (mpi*m)**2)) /(2.d0*srot)
           gamma = 1.d0/3.d0 * 1.d0/(4.d0*pi) * fdel_f**2*
     &             qcm**3/srot*(m+sqrt(m**2+qcm**2))/mpi**2  
         endif          
          sp = m**2+rq_zero**2-q_mod**2-
     &        2.d0*rq_zero *sqrt(m**2+3.d0/5.d0*k_fermi**2)
          if(sp.gt.(m+mpi)**2)then
           srotp = sqrt(sp)
           qcmp=sqrt(sp**2+mpi**4+m**4-2.d0*(sp*mpi**2+sp*m**2+
     &                    (mpi*m)**2))/(2.d0*srotp)
           gammap = 1.d0/3.d0 * 1.d0/(4.d0*pi) * fdel_f**2*
     &             qcmp**3/srotp*(m+sqrt(m**2+qcmp**2))/mpi**2  
          endif          
        endif
              
        z=md/(q_mod*k_fermi)*(q_zero-q_mod**2/(2.d0*md)
     &                         -wr +i*gamma/2.d0)

        zp=md/(q_mod*k_fermi)*(-q_zero-q_mod**2/(2.d0*md)
     &                          -wr +i*gammap/2.d0)
c
c care with limit cases
c
        if(abs(z).gt.50.d0)then
            pzeta =  2.d0/(3.d0*z) +2.d0/(15.d0*z**3)
        else if(abs(z).lt.1.d-2)then
            pzeta =  2.d0*z -2.d0/3.d0*z**3 -i*pi/2.d0*(1.d0 - z**2) 
        else
            pzeta =  z +  (1.d0-z**2) * log((z+1.d0)/(z-1.d0))/2.d0
        endif
        
        if(abs(zp).gt.50.d0)then
            pzetap =  2.d0/(3.d0*zp) +2.d0/(15.d0*zp**3)
        else if(abs(zp).lt.1.d-2)then
            pzetap =  2.d0*zp -2.d0/3.d0*zp**3 -i*pi/2.d0*(1.d0 - zp**2) 
        else
            pzetap =  zp + (1.d0-zp**2) * log((zp+1.d0)/(zp-1.d0))/2.d0
        endif


       delta_lind = 2.d0/3.d0 * rho * md/(q_mod*k_fermi) * ( 
     &        pzeta +pzetap) * fdel_f **2    
     
        RETURN
        end

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c NUCLEON relativistic Lindhard Function
c Same normalization as ULIN
c Real part
c taken from Eur.Phys.J.A25:299-318,2005 (Barbaro et al)
c Eq. 61
c
c Im. part: Juan. 
c
c INPUT: Real*8 
c  q0:Energy   [fm]
c  qm: modulus 3mom [fm]
c  kf: Fermi mom [fm]
c
c OUTPUT: Complex*16 [fm]
c
c USES: RULINRELX, DUR_J
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      COMPLEX*16 FUNCTION ULINREL(q0,qm,kf)
      IMPLICIT REAL*8 (A-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      INTEGER  ieta
      COMMON/testieta/ieta
      
c     range validity
      if(q0.gt.qm)then
         write(6,*)' Rlinrel is not valid for q0 > q'
         write(6,*)' q0=',q0,' q=',qm
         stop
      endif
      if( ieta.eq.1. ) Then
         m = 939.6d0/197.3269602d0
      else
         m = 938.3d0/197.3269602d0
      endif 
      uyi=(0.d0,1.d0)
c     
      RLINREL=RulinrelX(q0,qm,kf)+RulinrelX(-q0,qm,kf)
      ximag= 2.d0*DUR_j(q0,qm,kf,kf,m)
      ULINREL=RLINREL+uyi*ximag 
      
c      PRINT*,"ULINREL ",XIMAG,UYi
      
      RETURN
      END
C     
      REAL*8 function RulinrelX(q0,qm,kf)
c     fm, only real part, only |q0|<|q|
      IMPLICIT REAL*8 (A-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      INTEGER  ieta
      COMMON/testieta/ieta
c
      if( ieta.eq.1. ) Then
         m = 939.6d0/197.3269602d0
      else
         m = 938.3d0/197.3269602d0
      endif 
      
      pi=3.1415926535d0
      pi2=pi**2
      uy=(0.d0,1.d0)
c     
      ef=sqrt(m**2+kf**2)
      q2=q0**2-qm**2
      ds=Sqrt(1.d0-(4.d0*m**2)/q2)
c     
      m2 = m**2 
      qm2 = qm**2
c
      L1=Log((kf + ef)/m)
      uL2=Log(Abs((ef + q0 - Sqrt(m2+(kf-qm)**2))/
     &     (ef + q0 - Sqrt(m2 + (kf + qm)**2)))* 
     &     Abs((ef + q0 + Sqrt(m2 + (kf - qm)**2))/
     &     (ef + q0 + Sqrt(m2 + (kf + qm)**2))))
      
c      uL3=Log(Abs(((2*kf + q0*ds)**2-qm**2)/
c     &     ((2*kf - q0*ds)**2-qm**2))) + 
c     &     Log(Abs(((kf-ef*ds)**2 - (4*m**4*qm**2)/q2**2)/
c     &     ((kf+ef*ds)**2 - (4*m**4*qm**2)/q2**2)))
c     
      uL3=Log(Abs(((2*kf + q0*ds)**2-qm2)/
     &     ((2*kf - q0*ds)**2-qm2))* 
     &     (Abs(((kf-ef*ds)**2 - (4*m2**2*qm2)/q2**2)/
     &     ((kf+ef*ds)**2 - (4*m2**2*qm2)/q2**2))))

c      if( uL3_1.ne.uL3 ) print*,"uL3  ",uL3,uL3_1

      RlinrelX = -L1/(16.d0*pi2)+uL2*(2.d0*ef+q0)/(32.d0*pi2*qm)
     &     -uL3*ds/(64.d0*Pi2)
c     factor definition
      RulinrelX=RlinrelX*16.*m2
      RETURN
      end
c     
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
c  JUAN REL.+norel LINDHARD, SLOW
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE unucjuanb(irel,Q0,DQ,DKF,unuc_jb)               
      IMPLICIT DOUBLE PRECISION (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      DIMENSION pd(4000), xmud(4000), xmud2(4000)
      DIMENSION fmu(4000), fp(4000)
      DIMENSION fmu2(4000)
c     irel=0 no relativista
c     irel=1 relativista 
      ener(irel,dmnu,pmon2)=dfloat(1-irel)*(dmnu+0.5d0*pmon2/dmnu) +
     f     dfloat(irel)*dsqrt(dmnu**2+pmon2)
      
      
      if (irel.ne.1 .and. irel.ne.0) stop'wrong irel'
      
      call dsg20r(0.d0, dkf, 20, pd,np)
      call dsg20r(-1.d0, 1.d0, 20, xmud2,nmu2)
      
      dkf2= dkf**2
      dq2=dq**2 
      uyi= (0.d0,1.d0)
      
      
      dmnu = 940.d0 / 197.3269602d0
      

      epsilon = (0.5d0/ 197.3269602d0)*dq*197.3269602d0/650.
      pi2= 3.14159265d0**2
      
      xnuc_j= 0.d0

c     contrib. q y -q
      
      do ij=1,2
         
         
         do ip=1,np
            p= pd(ip)
            p2= p**2
            ep=  ener(irel,dmnu,p2)              
            xmu0= (dkf2-p2-dq2)/(2.d0*p*dq)
            
            
            if (xmu0.le.1.d0) then
               
               xmu_min= max(-1.d0, xmu0)
               
               call dsg20r(xmu_min, 1.d0, 20, xmud,nmu)
               
               do imu=1,nmu
                  xmu = xmud(imu)
                  
                  pq2=p2+dq2+2.d0*p*dq*xmu  
                  
                  
                  epq= ener(irel,dmnu,pq2)
                  
                  if (ij.eq.1) uaux= q0+ep-epq+uyi*epsilon
                  if (ij.eq.2) uaux= -q0+ep-epq+uyi*epsilon
                  
                  xreal_uaux= dreal(1.d0/uaux)

                  if (irel.eq.1) fmu(imu) =xreal_uaux*dmnu/epq
                  if (irel.eq.0) fmu(imu) =xreal_uaux
c     enddo xmu
               enddo
               
               call DRG20R(xmu_min, 1.d0, 20,fmu,resmu)
               
            else
               resmu= 0.d0
            endif
            

            if (irel.eq.1) then !  nueva contrib
               if (q0.gt.2.d0*dmnu)stop'contribuc imag, de antipart'
               
               do imu=1,nmu2
                  xmubis = xmud2(imu)
                  
                  pq2=p2+dq2+2.d0*p*dq*xmubis  
                  epq= ener(irel,dmnu,pq2)
                  
                  if (ij.eq.1) uaux= q0+ep+epq-uyi*epsilon
                  if (ij.eq.2) uaux= -q0+ep+epq-uyi*epsilon
                  
                  xreal_uaux= dreal(1.d0/uaux)
                  
                  fmu2(imu) =xreal_uaux*dmnu/epq
                  
c enddo xmubis
               enddo
               
               
               call DRG20R(-1.d0, 1.d0, 20,fmu2,resmu2)
               
c nueva contri
            endif                
            
            if (irel.eq.1) fp(ip) = (resmu-resmu2)*p2*dmnu/ep 
            if (irel.eq.0) fp(ip) = resmu*p2
c     enddo p
         enddo  
         
         call DRG20R(0.d0, dkf, 20,fp,resp)
         
         
         xnuc_j= xnuc_j+resp/pi2
         
c     enddo ij
      enddo
      
      
      if (irel.eq.1) ximag= 2.d0*DUR_j(q0,DQ,DKF,dkf,DMNU)
      if (irel.eq.0) ximag= 2.d0*DUIM_j(q0,DQ,DKF,dkf,DMNU)
      
      unuc_jb=xnuc_j+uyi*ximag 
      
      RETURN
      end
c     
c     
      FUNCTION DUIM_j(W,Q,K1,K2,DMNU)
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z)
      DOUBLE PRECISION KF,L2,K1,K2,NUM
 
      PI=3.141592653589793D0
      L2=(DMNU*W/Q-Q/2.d0)**2
      NUM=K2**2-2*DMNU*W
        ef1 = dmnu + k1**2 / (2.d0*dmnu)         
      IF((L2.LE.NUM).AND.(NUM.LE.K1**2))THEN
        DUIM_j=-DMNU/(4.d0*PI*Q)*(K1**2-K2**2+2.d0*DMNU*W)
      ELSEIF((NUM.LE.L2).AND.(L2.LT.K1**2))THEN
        DUIM_j=-DMNU/(4.d0*PI*Q)*(K1**2-L2)
      ELSE
        DUIM_j=0.d0
       ENDIF
      RETURN
      END
c
c
      FUNCTION DUR_j (Q0,DQ,PF1,PF2,DMNU)
      IMPLICIT DOUBLE PRECISION (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      DIMENSION DE(2000),DFE(2000)

c      IF (Q0.LT.0) THEN
c       WRITE(6,*)'EN LA LLAMADA A DUR,Q0 ES NEGATIVA'
c       DUR=0.
c       RETURN
c      ENDIF
 
      DPI=3.141592653589793D0
      CTE=-0.5D0*DMNU*DMNU/DPI/DQ
      Q2=Q0**2-DQ**2
      IF (Q2.gE.0) THEN
         DUR_j=0.
       RETURN
      endif
      EF1=DSQRT(DMNU**2+PF1**2)
      EF2=DSQRT(DMNU**2+PF2**2)
      IF ((EF2-Q0).GT.EF1) THEN
         DUR_j=0.
C         PRINT*,"(EF2-Q0).GT.EF1 "
       RETURN
      ENDIF
      a=0.5d0*(-q0+dq*sqrt(1.d0-4.d0*dmnu**2/q2))
      IF (a.GT.EF1) THEN
C         PRINT*,"a.gt.EF1",a,EF1
       DUR_j=0.
       RETURN
      ENDIF
         
      EINF=MAX(DMNU,EF2-Q0,a)
      DUR_j=(ef1-einf)*CTE
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccc
c Juan Rel. anterior / additional approximation
cccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE unucjuan(irel,Q0,DQ,DKF,unuc_j)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-T,V-Z)
      IMPLICIT COMPLEX*16 (U)
      DIMENSION pd(4000), xmud(4000)
      DIMENSION fmu(4000), fp(4000)
c     irel=0 no relativista
c     irel=1 relativista 
      ener(irel,dmnu,pmon2)=dfloat(1-irel)*(dmnu+0.5d0*pmon2/dmnu) +
     f     dfloat(irel)*dsqrt(dmnu**2+pmon2)
      
      
      if (irel.ne.1 .and. irel.ne.0) stop'wrong irel'
      call dsg20r(0.d0, dkf, 20, pd,np)
      dkf2= dkf**2
      dq2=dq**2 
      uyi= (0.d0,1.d0)
      dmnu = 940.d0 / 197.3269602d0
      epsilon = (0.5d0/ 197.3269602d0)*dq*197.3269602d0/650.
      pi2= 3.14159265d0**2
      
      xnuc_j= 0.d0
c contrib. q y -q
      do ij=1,2
         do ip=1,np
            p= pd(ip)
            p2= p**2
            xmu0= (dkf2-p2-dq2)/(2.d0*p*dq)
            
            if (xmu0.le.1.d0) then  
               xmu_min= max(-1.d0, xmu0)
               
               call dsg20r(xmu_min, 1.d0, 20, xmud,nmu)
               do imu=1,nmu
                  xmu = xmud(imu)
                  
                  pq2=p2+dq2+2.d0*p*dq*xmu  
                  
                  ep=  ener(irel,dmnu,p2)              
                  epq= ener(irel,dmnu,pq2)
                  
                  if (ij.eq.1) uaux= q0+ep-epq+uyi*epsilon
                  if (ij.eq.2) uaux= -q0+ep-epq+uyi*epsilon
                  
                  xreal_uaux= dreal(1.d0/uaux)
                  
                  if (irel.eq.1) fmu(imu) =xreal_uaux*dmnu/epq
                  if (irel.eq.0) fmu(imu) =xreal_uaux
c     enddo xmu
               enddo
               
               call DRG20R(xmu_min, 1.d0, 20,fmu,resmu)
               
            else
               resmu= 0.d0
            endif

            if (irel.eq.1) fp(ip) = resmu*p2*dmnu/ep 
            if (irel.eq.0) fp(ip) = resmu*p2
c     enddo p
         enddo  
         call DRG20R(0.d0, dkf, 20,fp,resp)
         xnuc_j= xnuc_j+resp/pi2
c     enddo ij
      enddo
      if (irel.eq.1) ximag= 2.d0*DUR_j(q0,DQ,DKF,dkf,DMNU)
      if (irel.eq.0) ximag= 2.d0*DUIM_j(q0,DQ,DKF,dkf,DMNU)
      
      unuc_j=xnuc_j+uyi*ximag 
      
      RETURN 
      END
      
      




	
****************************************************************************
C    *******************************************************************GAU00010
C              GAUS INTEGRATION SUBROUTINES
C    *******************************************************************GAU00030
      SUBROUTINE SG20R(A,B,N,X,NP)                                      GAU00040
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               GAU00050
      DIMENSION Y(10),X(2000)                                           GAU00060
      DATA Y/.9931285991,.9639719272,.9122344282,.8391169718,           GAU00070
     F .7463319064,.6360536807,.5108670019,.3737060887,                 GAU00080
     F .2277858511,.0765265211/                                         GAU00090
      NP=20*N                                                           GAU00100
      DINT=(B-A)/FLOAT(N)                                               GAU00110
      DELT=DINT*0.5                                                     GAU00120
      ORIG=A-DELT                                                       GAU00130
      I1=-20                                                            GAU00140
      DO 1 I=1,N                                                        GAU00150
      ORIG=ORIG+DINT                                                    GAU00160
      DORIG=ORIG+ORIG                                                   GAU00170
      I1=I1+20                                                          GAU00180
      I2=I1+21                                                          GAU00190
      DO 2 J=1,10                                                       GAU00200
      J1=I1+J                                                           GAU00210
      J2=I2-J                                                           GAU00220
      X(J1)=ORIG-DELT*Y(J)                                              GAU00230
 2    X(J2)=DORIG-X(J1)                                                 GAU00240
 1    CONTINUE                                                          GAU00250
      RETURN                                                            GAU00260
      END                                                               GAU00270
      SUBROUTINE RG20R(A,B,N,CF,CRES)                                   GAU00280
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               GAU00290
      DIMENSION W(10),CF(2000)                                          GAU00300
      DATA W/.0176140071,.0406014298,.0626720483,.0832767415,           GAU00310
     F .1019301198,.1181945319,.1316886384,.1420961093,.1491729864,     GAU00320
     F .1527533871/                                                     GAU00330
      CR=(0.,0.)                                                        GAU00340
      I1=-20                                                            GAU00350
      DO 1 I=1,N                                                        GAU00360
      I1=I1+20                                                          GAU00370
      I2=I1+21                                                          GAU00380
      DO 2 J=1,10                                                       GAU00390
      J1=I1+J                                                           GAU00400
      J2=I2-J                                                           GAU00410
 2    CR=CR+W(J)*(CF(J1)+CF(J2))                                        GAU00420
 1    CONTINUE                                                          GAU00430
      CRES=CR*0.5*(B-A)/FLOAT(N)                                        GAU00440
      RETURN                                                            GAU00450
      END                                                               GAU00460
                                                                        GAU00470
                                                                        GAU00480
C    *******************************************************************GAU00490
C               SUBROUTINAS DE INTEGRACION NUMERICA POR GAUSS           GAU00500
C    *******************************************************************GAU00510
      SUBROUTINE DSG20R(A,B,N,X,NP)                                     GAU00520
      IMPLICIT REAL*8 (A-H,O-Z)                                         GAU00530
      DIMENSION Y(10),X(2000)                                           GAU00540
      DATA Y/.9931285991,.9639719272,.9122344282,.8391169718,           GAU00550
     F .7463319064,.6360536807,.5108670019,.3737060887,                 GAU00560
     F .2277858511,.0765265211/                                         GAU00570
      NP=20*N                                                           GAU00580
      DINT=(B-A)/DBLE(N)                                                GAU00590
      DELT=DINT*0.5D0                                                   GAU00600
      ORIG=A-DELT                                                       GAU00610
      I1=-20                                                            GAU00620
      DO 1 I=1,N                                                        GAU00630
      ORIG=ORIG+DINT                                                    GAU00640
      DORIG=ORIG+ORIG                                                   GAU00650
      I1=I1+20                                                          GAU00660
      I2=I1+21                                                          GAU00670
      DO 2 J=1,10                                                       GAU00680
      J1=I1+J                                                           GAU00690
      J2=I2-J                                                           GAU00700
      X(J1)=ORIG-DELT*Y(J)                                              GAU00710
 2    X(J2)=DORIG-X(J1)                                                 GAU00720
 1    CONTINUE                                                          GAU00730
      RETURN                                                            GAU00740
      END                                                               GAU00750
      SUBROUTINE DRG20R(A,B,N,CF,CRES)                                  GAU00760
      IMPLICIT REAL*8 (A-H,O-Z)                                         GAU00770
      DIMENSION W(10),CF(2000)                                          GAU00780
      DATA W/.0176140071,.0406014298,.0626720483,.0832767415,           GAU00790
     F .1019301198,.1181945319,.1316886384,.1420961093,.1491729864,     GAU00800
     F .1527533871/                                                     GAU00810
      CR=(0.D0,0.D0)                                                    GAU00820
      I1=-20                                                            GAU00830
      DO 1 I=1,N                                                        GAU00840
      I1=I1+20                                                          GAU00850
      I2=I1+21                                                          GAU00860
      DO 2 J=1,10                                                       GAU00870
      J1=I1+J                                                           GAU00880
      J2=I2-J                                                           GAU00890
 2    CR=CR+W(J)*(CF(J1)+CF(J2))                                        GAU00900
 1    CONTINUE                                                          GAU00910
      CRES=CR*0.5D0*(B-A)/DBLE(N)                                       GAU00920
      RETURN                                                            GAU00930
      END                                                               GAU00940


      SUBROUTINE CRG20R(A,B,N,CF,CRES)                                  GAU00760
      IMPLICIT REAL*8 (A-H,O-Z)
      complex*16 CF,cres,CR 
      DIMENSION W(10),CF(2000)                                          GAU00780
      DATA W/.0176140071,.0406014298,.0626720483,.0832767415,           GAU00790
     F .1019301198,.1181945319,.1316886384,.1420961093,.1491729864,     GAU00800
     F .1527533871/                                                     GAU00810
      CR=(0.D0,0.D0)                                                    GAU00820
      I1=-20                                                            GAU00830
      DO 1 I=1,N                                                        GAU00840
      I1=I1+20                                                          GAU00850
      I2=I1+21                                                          GAU00860
      DO 2 J=1,10                                                       GAU00870
      J1=I1+J                                                           GAU00880
      J2=I2-J                                                           GAU00890
 2    CR=CR+W(J)*(CF(J1)+CF(J2))                                        GAU00900
 1    CONTINUE                                                          GAU00910
      CRES=CR*0.5D0*(B-A)/DBLE(N)                                       GAU00920
      RETURN                                                            GAU00930
      END                                                               GAU00940
