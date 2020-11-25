      Program crscalc

      Implicit None

#include "neutmodel.h"
#include "rscons.h"

*     Lepton mass effects are taken into account
*     following the formalism by Ch.Berger and L.M.Sehgal
*     Phys. Rev. D 76, 113004 (2007)
*     Phys. Rev. D 77, 059901(E) (2008) !! Erratum !!

      integer ie,iflag,nq,iq,iw,nw,ipflag,intflag,inlep
      real wmin,dw,dq2,sigma,sigma3,qmax,smax,wwmax
      real siw,siw3,q2,e,w
      real dsi(2),dsi3(2)
      real dsi_dbg(2),dsi3_dbg(2)

      integer iparlop

      integer*4 iptbl(3)
      data iptbl /12,14,16/

      integer IEMAX /24/
      real etbl(24)

      DATA ETBL/
     $     0.370, 0.400, 0.500, 0.625, 0.750,
     $     0.875, 1.000, 1.125, 1.250, 1.375,
     $     1.500, 2.000, 2.500, 3.000, 3.500,
     $     4.750, 6.000, 10.00, 50.00, 100.0,
     $     200.0, 400.0, 650.0, 1000.0/

      character*15 fname_ref(3)
      character*100 fname
      character*50 argv1
      character*50 argv2

      data fname_ref/
     $     'spi_nue_xsec_ma',
     $     'spi_num_xsec_ma',
     $     'spi_nut_xsec_ma'/

      character*25 typesn(7),typesa(7)
      data typesn/
     $     "nu p -> mu- p pi+",
     $     "nu n -> mu- p pi0",
     $     "nu n -> mu- n pi+",
     $     "nu p -> nu- p pi0",
     $     "nu p -> nu- n pi+",
     $     "nu n -> nu- n pi0",
     $     "nu n -> nu- p pi-"/
      data typesa/
     $     "nubar p -> mu+   n pi-",
     $     "nubar p -> mu+   n pi0",
     $     "nubar p -> mu+   p pi-",
     $     "nubar n -> nubar n pi0",
     $     "nubar n -> nubar p pi-",
     $     "nubar p -> nubar p pi0",
     $     "nubar p -> nubar n pi+"/
      character*30 ffn(2)
      data ffn/
     $     "Rein-Seghal default",
     $     "Graczyk & Sobczyk default"/

      REAL*4 ERES,ELEP,PRES,PLEP,COS1
      REAL*4 XMLEP,XMLEP2

      real*4 XMASPI2, XMVSPI2

      integer i,j,k
      real Xsec(7,2,7,38), Xsec3(7,2,7,38), Enu(7,2,7,38)

      real W_max,s,eorg,q2org,worg,CosLep
      integer iW_max,Lambda(2)
      real Q2_min, Q2_max
      integer ieffect, iOrg
      data ieffect / 0 / ! Lepton mass effect on(0), off(1)
      real DW_real, DQ2_real

      if (iargc().lt.4) then
         print *, 'iargc is ', iargc()
         stop 'args: nu type, MAres, CA5(0), BGscl'
      endif

C     RS MA_Res 1.21, 1.11
C     G&S MA_Res Default 1: 0.95; 2: 1.15
      call getarg(2, argv2)
      read(argv2, '(F4.2)') XMASPI
      print *, 'XMA=', XMASPI


C     G&S CA5(0) Default 1: 1.2; 2: 0.89
      call getarg(3, argv2)
      read(argv2, '(F4.2)') RCA5ISPI
      print *, 'RCA5ISPI=', RCA5ISPI

C     Background scaling
      call getarg(4, argv2)
      read(argv2, '(F4.2)') RBGSCLSPI
      print *, 'RBGSCLSPI=', RBGSCLSPI


C     Form factor variables
      XMVSPI = 0.84
      XMASPI2 = XMASPI*XMASPI
      XMVSPI2 = XMVSPI*XMVSPI

C     Form factor models
C     Rein-Seghal: 0, Graczyk & Sobcyzk (default):1
      IFFSPI = 1

C     NRTYPESPI only affects RS model
C     Equations 12, 13 in Graczyk & Sobcyzk (Phys Rev D 77 053001)
C     Electroproduction: 0, Neutrino-production (default):1
      NRTYPESPI = 1

      print*, '==============================='
      if (ieffect.eq.0) then
         print*, '== LEPTON MASS EFFECT -> ON ==='
      elseif (ieffect.eq.1) then
         print*, '== LEPTON MASS EFFECT -> OFF =='
      else
         stop' please set LEPTONMASS'
      endif
      print*, '==============================='

      WMIN = XMN+XMPI
      WWMAX = 2.0

      DW = 0.02
      DQ2 = 0.05

      do iparlop=1,3
         do i=1,2
            do j=1,7
               do k=1,38
                  Enu(iparlop,i,j,k) = 0.
                  Xsec(iparlop,i,j,k) = 0.
                  Xsec3(iparlop,i,j,k) = 0.
               enddo
            enddo
         enddo
      enddo

      call getarg(1, argv1)
      read(argv1, '(I10)') IPARLOP
      print *, 'IPARLOP=', IPARLOP

      if ((IPARLOP.GT.3).or.(IPARLOP.lt.1)) then
         write(*,*) 'IPARLOP=', IPARLOP, ' invalid'
         goto 9999
      endif

      DO 110  IPFLAG = 0,1
         DO 105  IFLAG = 1,7

            INTFLAG = IFLAG+IPFLAG*10

C            DO 5000 IPARLOP=1,3
               iOrg = ieffect ! Consider lepton mass effecs(0) or not(1)
               IF (IPARLOP.eq.1)  XMLEP=XME
               IF (IPARLOP.eq.2)  XMLEP=XMMU
               IF (IPARLOP.eq.3)  XMLEP=XMTAU
               IF (IFLAG.gt.3)    then
                  iOrg = 1 ! Original R-S model
                  XMLEP=0.0E0
               ENDIF

               write(*,*)  "XMLEP=",XMLEP
               write(*,*) "WMAX=",WWMAX

               XMLEP2 = XMLEP**2

C              This skipped the NC contribution for numu and nutau.  
C              IF ((IFLAG.gt.3).and.(IPARLOP.gt.1)) GOTO 5000

C     Energy integration
               DO 95 IE = 1,IEMAX
C               DO 95 IE = 1,5
                  E = etbl(IE)
*     Lepton mass effects are negligible in higher energy
                  if (IPARLOP.eq.1.and.E.gt.1.5 ) iOrg = 1
                  if (IPARLOP.eq.2.and.E.gt.100.) iOrg = 1
                  if (IPARLOP.eq.3.and.E.gt.200.) iOrg = 1

                  SIGMA = 0.
                  SIGMA3 = 0.
                  SMAX = (2 * XMN * E) + XMN2
                  QMAX = (SMAX - XMN2)
                  W_max = MIN(SQRT(SMAX)-XMLEP,WWMAX)
                  iW_max = INT((W_max-WMIN)/dW) + 1

                  DO 70  IW = 1,iW_max
                     if (iW .eq. iW_max) then
                        W = 0.5 * (WMIN+DW*(FLOAT(IW)-1.) + W_max)
                        DW_real = 2. * (W_max - W)
                     else
                        W = WMIN+DW*(FLOAT(IW)-.5)
                        DW_real = DW
                     endif

                     SIW = 0.
                     SIW3 = 0.

                     call res_q2_limits(SMAX,W,XMLEP,Q2_min,Q2_max)

C                     q2_min = 0.
C                     q2_max = (( SMAX - XMN2 )*( SMAX - W*W ))/SMAX
                     NQ = int((Q2_max-Q2_min)/DQ2) + 1

                     DO 90 IQ = 1,NQ
                        if (IQ .eq. NQ) then
                           Q2 = -0.5*(Q2_min + DQ2*(float(IQ) - 1.) +
     $                          Q2_max)
                           DQ2_real = 2.*(Q2_max + Q2)
                        else
                           Q2 = -1.*(Q2_min + DQ2*(float(IQ) - 0.5))
                           DQ2_real = DQ2
                        endif


                        ERES=(XMN2+W**2-q2)/(2*XMN)
                        IF (ERES.lt.W) THEN
c     print*, ERES,W
                           GOTO 90
                        ENDIF

                        ELEP=E+XMN-ERES
                        IF (ELEP.lt.XMLEP) THEN
c     print*, Elep, XMlep
                           GOTO 90
                        ENDIF

                        PLEP=sqrt(ELEP**2-XMLEP**2)
                        PRES=sqrt(ERES**2-W**2)
                        CosLep=ABS(q2-XMLEP**2+2*E*ELEP)/(2*E*PLEP)
                        IF (CosLep.gt.1.) THEN
c     print*, 'coslep ',coslep
                           GOTO 90
                        ENDIF

                        COS1=ABS(2*E*XMN+XMN2-XMLEP**2-W**2-2*ELEP*ERES)
     $                       /(2*PLEP*PRES)

                        IF (COS1.gt.1.) THEN
c     print*, 'cos1 ',cos1
                           GOTO 90
                        ENDIF

                        if (IPFLAG.eq.0) then
                           Lambda(1) =  1
                           Lambda(2) = -1
                        else
                           Lambda(1) = -1
                           Lambda(2) =  1
                        endif
                        eorg  = e
                        q2org = q2
                        worg  = w

                        call rsdcrs(iNTflag,iOrg,Lambda(1),xmlep,
     $                       e,q2,w,dsi(1),dsi3(1))
                        e  = eorg
                        q2 = q2org
                        w  = worg
                        call rsdcrs(iNTflag,iOrg,Lambda(2),xmlep,
     $                       e,q2,w,dsi(2),dsi3(2))

                        if (dsi(1).lt.0. .or. dsi(2).lt.0.) stop

                        SIW  = SIW  + DQ2_real*(dsi(1) + dsi(2))
                        SIW3 = SIW3 + DQ2_real*(dsi3(1)+ dsi3(2))
 90                  CONTINUE

                     SIGMA  = SIGMA  + DW_real*SIW
                     SIGMA3 = SIGMA3 + DW_real*SIW3
 70               CONTINUE
                  WRITE(*,*) E,SIGMA,sigma3

                  Xsec( IPARLOP,IPFLAG+1,IFLAG,IE) = SIGMA
                  Xsec3(IPARLOP,IPFLAG+1,IFLAG,IE) = SIGMA3
                  Enu(  IPARLOP,IPFLAG+1,IFLAG,IE) = E


 95            CONTINUE

                  write(*,*)  "#",TYPESn(IFLAG)
                  write(*,*)  "#",TYPESa(IFLAG)

 5000       CONTINUE


 105     CONTINUE
 110  CONTINUE
      do I=iparlop,iparlop

         write(fname,'(A,F4.2,A,F4.2,A,F4.2,A)') 
     $     fname_ref(I), xmaspi, '_ca5I', RCA5ISPI, '_bgscl', 
     $        RBGSCLSPI, '.dat'
         write(*,*) 'Opening file ', fname

         open(60,FILE=FNAME,ERR=9999)

         write(60,'(A,A,A,F4.2,A,F4.2,A,F4.2,A,I1)') 
     $        '# 1pi :', ffn(IFFSPI+1), ':MA=', XMAspi, 
     $        ', CA5(0)=', RCA5ISPI, ', BGscale=', RBGSCLSPI, 
     $        ', NRTYPESPI=', NRTYPESPI

         do IE=1,IEMAX
            write(60,'(I4,I2,F8.2,7E13.6)')
     $           IPTBL(I),1,Enu(I,1,1,IE),
     %           Xsec(I,1,1,IE), Xsec(I,1,2,IE), Xsec(I,1,3,IE),
     %           Xsec(I,1,4,IE), Xsec(I,1,5,IE), Xsec(I,1,6,IE),
     $           Xsec(I,1,7,IE)
         enddo
         do IE=1,IEMAX
            write(60,'(I4,I2,F8.2,7E13.6)')
     $           IPTBL(I),2,Enu(I,1,1,IE),
     %           Xsec3(I,1,1,IE), Xsec3(I,1,2,IE), Xsec3(I,1,3,IE),
     %           Xsec3(I,1,4,IE), Xsec3(I,1,5,IE), Xsec3(I,1,6,IE),
     $           Xsec3(I,1,7,IE)
         enddo
         do IE=1,IEMAX
            write(60,'(I4,I2,F8.2,7E13.6)')
     $        -1*IPTBL(I),1,Enu(I,1,1,IE),
     $           Xsec(I,2,1,IE), Xsec(I,2,2,IE), Xsec(I,2,3,IE),
     $           Xsec(I,2,4,IE), Xsec(I,2,5,IE), Xsec(I,2,6,IE),
     $           Xsec(I,2,7,IE)
         enddo
         do IE=1,IEMAX
            write(60,'(I4,I2,F8.2,7E13.6)')
     $           -1*IPTBL(I),2,Enu(I,1,1,IE),
     $           Xsec3(I,2,1,IE), Xsec3(I,2,2,IE), Xsec3(I,2,3,IE),
     $           Xsec3(I,2,4,IE), Xsec3(I,2,5,IE), Xsec3(I,2,6,IE),
     $           Xsec3(I,2,7,IE)

         enddo
      enddo
      goto 1000

 999  write(*,*) "Error occured in rsdcrs."
      write(*,*) "intflag=",intflag
      write(*,*) "e    =",e
      write(*,*) "q2   =",q2
      write(*,*) "w    =",w
      write(*,*) "dsi  =",dsi(1),dsi(2)

 1000 continue
      STOP
 9999 write(*,*) "Error in opening file."
      END
