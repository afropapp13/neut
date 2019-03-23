C*****
C     
C     This Program calculates total cross-section for single
C     pion production by using Rein-Sehgal.
C     The output file will be used by other program.
C     
C     May-94 (Ver 0.00) By Y.Hayato
C     
C     
C***  Notice
C     Kind of neutrino is not considered.
C     Nuclear effect is not considered.
C     
      
      Implicit None

#include <rscons.h>
      
C      REAL  MN,MN2,MPI,MPI2,MV,MV2
C      REAL  MA,MA2,Z,PI,OMEG,XW,EPSI
      
      REAL*4 xmlep

C      COMMON /RSCONS/Z,PI,MN,MN2,MPI,MPI2,MV,MV2,MA,MA2,OMEG,XW,EPSI
      
      real junk
      
      integer ie,iflag,nq,iq,iw,nw,ipflag,intflag
      real wmin,wmax,dw,dq2,sigma,sigma3,qmax,smax
      real siw,siw3,q2,e,w
      real dsi,dsi3
      
      integer IEMAX / 20/
      real etbl(20)
      
      data etbl/0.2,0.4,0.5,0.625,0.75,0.875,1.00,1.125,1.25,
     $     1.375,1.50,2.00,2.50,3.00,3.50,4.75,6.00,10.0,50.0,
     $     100.0/
      
      character*10 fname1(7),fname2(7)
      data fname1/"nu.tcrs.aa","nu.tcrs.ab","nu.tcrs.ac",   
     $     "nu.tcrs.ad","nu.tcrs.ae","nu.tcrs.af",
     $     "nu.tcrs.ag"/
      data fname2/"an.tcrs.aa","an.tcrs.ab","an.tcrs.ac",   
     $     "an.tcrs.ad","an.tcrs.ae","an.tcrs.af",
     $     "an.tcrs.ag"/
      character*25 typesn(7),typesa(7)
      data typesn/"nu p -> mu- p pi+","nu n -> mu- p pi0",     
     $     "nu n -> mu- n pi+","nu p -> nu- p pi0",     
     $     "nu p -> nu- n pi+","nu n -> nu- n pi0",     
     $     "nu n -> nu- p pi-"/
      data typesa/"nubar p -> mu+   n pi-",
     $     "nubar p -> mu+   n pi0",
     $     "nubar p -> mu+   p pi-",
     $     "nubar n -> nubar n pi0",
     $     "nubar n -> nubar p pi-",
     $     "nubar p -> nubar p pi0",
     $     "nubar p -> nubar n pi+"/
C     
C     
C     Program begins from here.      
C     
C     
C      
C     set COMMON /RSCONS/
C
      call rsstcm
C      
C     Set Parameters
C     
      WMIN = XMN+XMPI
      WMAX = 1.4
      
      DW = .02
      DQ2 = 0.05
C      
C     Read Parameters
C      
      write(*,*) "WMAX=[",WMAX,"]"
      call rdpara(WMAX)
      write(*,*) "IEMAX=[",IEMAX,"]"
      junk=float(IEMAX)
      call rdpara(junk)
      IEMAX=int(junk)
      write(*,*) IEMAX
      if (IEMAX.eq.20) goto 12
      do 11 ie=1,IEMAX
         write(*,*) "ETBL(",ie,")=",ETBL(IE)
         call rdpara(ETBL(IE))
 11   continue
 12   continue
      write(*,*) "DW=[",DW,"]"
      call rdpara(DW)
      write(*,*) "DQ2=[",DQ2,"]"
      call rdpara(DQ2)
C
C     Main Loop
C
C        IPFLAG == 0 -> neutrino
C        IPFLAG == 1 -> anti neutrino
C      
      DO 110  IPFLAG = 0,1
         DO 105  IFLAG = 1,7
            
            if (IPFLAG.eq.0) then
               open(201,FILE=FNAME1(IFLAG),STATUS="NEW",ERR=9999)
            else
               open(201,FILE=FNAME2(IFLAG),STATUS="NEW",ERR=9999)
            endif
            
            INTFLAG = IFLAG+IPFLAG*10
            
            DO 95 IE = 1,IEMAX
               
               E = etbl(IE)
C               
C     NOW  E was FIXED
C               
               SIGMA = 0.
               SIGMA3 = 0.
C               
C     Calc Kinematics(1)
C
               SMAX = (2 * XMN * E) + XMN2
               QMAX = (SMAX - XMN2)
               
               NQ = int(QMAX / DQ2)
               
C
C     Q-integration
C               
               DO 90 IQ = 1,NQ
                  Q2 = -DQ2 * (float(IQ) - 0.5)
C                  
C     W-INTEGRATION
C                  
                  SIW = 0.
                  SIW3 = 0.
                  
                  NW = INT((WMAX - WMIN) / DW)
                  
                  DO 70  IW = 1,NW
                     
                     W = WMIN+DW*(FLOAT(IW)-.5)
C
C     Calc differential cross-section(dsigma/dq^2dW)
C                     
C                     call rsdcrs(iNTflag,e,q2,w,dsi,dsi3)
                     xmlep = 0.511 / 1000.
                     call rsdcrs(iNTflag,0,0,xmlep,
     $                    e,q2,w,dsi(1),dsi3(1))

                     if (dsi.lt.0.0) goto 999
                     if (dsi3.lt.0.0) goto 999
                     
                     SIW = SIW+DW*DSI
                     SIW3 = SIW3+DW*DSI3
                     
 70               CONTINUE
                  
                  SIGMA = SIGMA + DQ2 * SIW
                  SIGMA3 = SIGMA3 + DQ2 * SIW3
                  
 90            CONTINUE
               
               WRITE(201,*) E,SIGMA,sigma3
               
 95         CONTINUE
            if (IPFLAG.eq.0) then
               write(201,*) "#",TYPESn(IFLAG)
            else
               write(201,*) "#",TYPESa(IFLAG)
            endif
            close(201)
 105     CONTINUE
 110  CONTINUE
      
      goto 1000
      
C
C     Error trap(Differential Cross-section)
C      
 999  write(*,*) "Error occured in rsdcrs."
      write(*,*) "intflag=",intflag
      write(*,*) "e    =",e
      write(*,*) "q2   =",q2
      write(*,*) "w    =",w
      write(*,*) "dsi  =",dsi
      
 1000 continue
      STOP
C
C     Error Trap(file access)
C      
 9999 write(*,*) "Error in opening file."
      if (IPFLAG.eq.0) then
         write(*,*) "FILE Name=",FNAME1(IFLAG)
      else
         write(*,*) "FILE Name=",FNAME2(IFLAG)
      endif
      END
C
C     subroutine for reading parameter
C
      subroutine rdpara(para)
      real para
      real junk
      
      read(*,*) junk
      if (junk.ne.0.) para=junk
      end
