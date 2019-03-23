***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE RSSKVCT(IPAR,JMOD,E,DNEUT,IERR)
*     ---------------------------------------------------
*
*     ( purpose )
*       VECTOR GENERATION FOR SINGLE  K  PRODUCTION
*
*     ( input )
*       IPAR        : kind of incoming neutrino
*       JMOD        : interaction mode
*       E           : Energy of incident neutrino
*       dneut(3)    : direction of incoming neutrino
*     
*     ( output )
*       IERR        : ERROR CODE
*       nework.h
*
*         Output meson type is   K+:321,K-:-321,K0:311           
*                baryon type is  Lambda0 
*
*     ( creation date and author )
*       1993.Dec. ; for Super-Kamioka by Y.Hayato
*                   (Based on Rein and Sehgal's code)
*       1995.Sep. ; BUG fix(around calc. sinth,costh,sinphi,cosphi)
*       1995.Nov. ; BUG fix(set IORGNE,IFLGNE)
*       1996.Nov. ; BUG fix(Avoid Infinite LOOP,Check more carefully)
*       1998.Mar. ; J.Kameda for  K
*       
*     ( comment )
*       Reference:Rein and Sehgal,Ann. of Phys. 133(1981),79-153
*                 Rein           ,Z.Phys.C 35(1987),43-64
*                 Feynman et al. ,Phys.Rev.D3(1971),2706-
***********************************************************************


      IMPLICIT NONE

#include<nework.h>      

      REAL DNEUT(3),PNEU(3),DIRN(3)

c      INTEGER MAXNE
c      INTEGER MODENE,NUMNE,IPNE,IORGNE,IFLGNE,ICRNNE
C      REAL    PNE

      REAL*4  TMPEP,TMPEL,TMPEEE
      REAL*4  TMPSIG
      INTEGER*4 IMOD
      REAL*4  RSSKCRS
      EXTERNAL RSSKCRS

C      INTEGER IFREE,IPAR,JMOD,IERR,IEC,IDECAY,I,INO16,IWAVE
      INTEGER IPAR,JMOD,IERR,IEC,IDECAY,I,IWAVE
      INTEGER IFLAG
      REAL    E,PFABS,RTPLEP(3),RTPK(3),RTPNC(3),DX,DY,DZ
      REAL    PNCEPZ(3)
      REAL    COSTH,SINTH,COSPH,SINPH
#include<rscons.h>

#include "vcwork.h"      
#include "posinnuc.h"
C      REAL*4 XIN(3)
C      COMMON /SAVXIN/INO16,XIN 
 
C      COMMON/FREPRO/IFREE
C
      IERR=0
      IEC=0
      IDECAY=0
C
C -- SET IMOD AND LVECT
C
      NUMNE=5

C-CONVERT MODE
      IF (JMOD.eq.23)  IMOD=1
      IF (JMOD.eq.44)  IMOD=11
      IF (JMOD.eq.45)  IMOD=12
      IF (JMOD.eq.-23) IMOD=2
      IF (JMOD.eq.-44) IMOD=13
      IF (JMOD.eq.-45) IMOD=14

C
C -- SET NEUTRINO DIRECTION
C 
      DO 100 I=1,3
         PNE(I,1)=E*DNEUT(I)
         PNEU(I)=PNE(I,1)
  100 CONTINUE
      IPNE(1)=IPAR

C---

C      INO16=0
C      IF(IFREE.EQ.0) INO16=1
      
      PFABS=0.
C 1500 IF(INO16.EQ.1 .AND. ABS(IPAR).NE.16) 
C 1500 IF(IBOUND.EQ.1 .AND. ABS(IPAR).NE.16) 
 1500 IF(IBOUND.EQ.1)
     & CALL nefrmmom(PFABS,IWAVE)
C     & CALL RNFERM(PFABS,IWAVE)
      PFABS = PFABS/1000.

      CALL RNDIR(DIRN)
      PNE(1,2)=DIRN(1)*PFABS
      PNE(2,2)=DIRN(2)*PFABS
      PNE(3,2)=DIRN(3)*PFABS

      TMPEP=sqrt(XMN2+(PNE(1,2)**2+PNE(2,2)**2+PNE(3,2)**2))
      TMPEL=E

      TMPEEE=( (TMPEP+TMPEL)**2
     $        -( (PNE(1,2)+PNEU(1))**2+(PNE(2,2)+PNEU(2))**2
     $          +(PNE(3,2)+PNEU(3))**2)-XMN2)/2./XMN

C----------------------------------------------------------------
C      write(*,*) "E=",E,"/P(Nucl)=",PFABS
C      write(*,*) "ENEU(REST)=",TMPEEE
C      write(*,*) "S=",( (TMPEP+TMPEL)**2
C     $        -( (PNE(1,2)+PNEU(1))**2+(PNE(2,2)+PNEU(2))**2
C     $          +(PNE(3,2)+PNEU(3))**2))
C      write(*,*) "S=",(TMPEEE+XMN)**2-TMPEEE**2
C----------------------------------------------------------------

C-Check IF sigma < 1.e-43 then re-calculate fermi-momentum 
      TMPSIG=RSSKCRS(TMPEEE,IPAR,IMOD)
      IF (TMPSIG.le.1.e-5) GO TO 1500

C--- initialize/set some commons
      IORGNE(1)=0
      IORGNE(2)=0
C    lepton :  from initial lepton
      IORGNE(3)=1
C    pion and nucleon :from initialn nucleon
      IORGNE(4)=2
      IORGNE(5)=2

C
      IFLGNE(1)=-1
      IFLGNE(2)=-1
      IFLGNE(3)=0
      IFLGNE(4)=0
      IFLGNE(5)=0

      MODENE = JMOD
C ---------------------------------------------------------------
C
C
C
C -- SET KIND OF PARTICLE
C
      IF(IPAR.LT.0)GO TO 220

      IF(JMOD.EQ.23)GO TO 201
      IF(JMOD.EQ.44)GO TO 202
      IF(JMOD.EQ.45)GO TO 203
  220 IF(JMOD.EQ.-23)GO TO 204
      IF(JMOD.EQ.-44)GO TO 205
      IF(JMOD.EQ.-45)GO TO 206


      GO TO 999

  201 IPNE(2)=2112
      IPNE(3)=IPAR-1
      IPNE(4)=3122
      IPNE(5)=321
      IFLAG=1
      GO TO 300
  202 IPNE(2)=2112
      IPNE(3)=IPAR
      IPNE(4)=3122
      IPNE(5)=311
      IFLAG=3
      IFLGNE(3)=2
      GO TO 300
  203 IPNE(2)=2212
      IPNE(3)=IPAR
      IPNE(4)=3122
      IPNE(5)=321
      IFLAG=2
      IFLGNE(3)=2
      GO TO 300
  204 IPNE(2)=2212
      IPNE(3)=IPAR+1
      IPNE(4)=3122
      IPNE(5)=311
      IFLAG=11
      GO TO 300
  205 IPNE(2)=2112
      IPNE(3)=IPAR
      IPNE(4)=3122
      IPNE(5)=311
      IFLAG=13
      IFLGNE(3)=2
      GO TO 300
  206 IPNE(2)=2212
      IPNE(3)=IPAR
      IPNE(4)=3122
      IPNE(5)=321
      IFLAG=12
      IFLGNE(3)=2
      GO TO 300

C
C
C
C
 300  CONTINUE

CC  translate K0 -> K (LONG) or K (SHORT)
CC  K0 : 311 KL : 130 KS : 310

CC     IF (IPNE(5).eq.311) then
CC        IPNE(5) = 130
CC        RNTYP = RNDM(DUM)
CC        IF (RNTYP.ge.0.5) IPNE(5) = 310
CC     ENDIF 

cc
      DO 310 I=1,3
         PNCEPZ(I)=PNE(I,2)
 310  CONTINUE

      IF (IPNE(3).ne.IPAR) ICRNNE(3) = 1
      ICRNNE(4) = 1
      ICRNNE(5) = 1

C-ORiG    CALL RSDISTK(IFLAG,IPNE(1),E,PNCEPZ(1),RTPLEP,RTPK,RTPNC)
      
      DZ=DNEUT(3)
      DX=DNEUT(1)
      DY=DNEUT(2)

      COSTH=DZ
      if (abs(dz).eq.1.) then
         sinth = 0.
      else
         SINTH=SQRT(1.-DZ**2)
      endif
      IF (ABS(DX).GT.0.000001) THEN
        COSPH=COS(ATAN2(DY,DX))
        SINPH=SIN(ATAN2(DY,DX)) 
      ELSE
        COSPH=0.
        IF (ABS(DY).GT.0.00001) THEN
          SINPH=1.*(ABS(DY)/DY)
        ELSE
          SINPH=1.
        ENDIF
      ENDIF 
 
C-----------------------------------------------------------------
C      write(*,*) "NEUTRINO=",PNE(1,1),",",PNE(2,1),",",PNE(3,1)
C      CALL RSROTVEC(PNE(1,1),PNCEPZ(1),
C     $     1.,0.,COSPH,-1.*SINPH)
C      CALL RSROTVEC(PNCEPZ(1),PNCEPZ(1),
C     $     COSTH,-1.*SINTH,1.,0.)
C      write(*,*) "NEUTRINO=",PNCEPZ(1),",",PNCEPZ(2)
C     $     ,",",PNCEPZ(3)
C      write(*,*) "NUCLEON=",PNE(1,2),",",PNE(2,2),",",PNE(3,2)
C-----------------------------------------------------------------

      CALL RSROTVEC(PNE(1,2),PNCEPZ(1),
     $     1.,0.,COSPH,-1.*SINPH)
      CALL RSROTVEC(PNCEPZ(1),PNCEPZ(1),
     $     COSTH,-1.*SINTH,1.,0.)

C-----------------------------------------------------------------
C      write(*,*) "NUCLEON=",PNCEPZ(1),",",PNCEPZ(2)
C     $     ,",",PNCEPZ(3)
C-----------------------------------------------------------------
      CALL RSDISTK(IFLAG,IPNE(1),E,PNCEPZ(1),RTPLEP,RTPK,RTPNC,IERR)

      CALL RSROTVEC(RTPLEP,PNE(1,3),COSTH,SINTH,COSPH,SINPH)
      CALL RSROTVEC(RTPNC ,PNE(1,4),COSTH,SINTH,COSPH,SINPH)
      CALL RSROTVEC(RTPK ,PNE(1,5),COSTH,SINTH,COSPH,SINPH)

C------------------------------------------------------------------
C      write(*,*) "Conservation Laws"
C      write(*,*) "INCOMING E=",
C     $     (E+sqrt(PNE(1,2)**2+PNE(2,2)**2+PNE(3,2)**2+XMN2))
C      write(*,*) "OUTGOING E=",
C     $     ( sqrt(PNE(1,4)**2+PNE(2,4)**2+PNE(3,4)**2+XMN2)
C     $      +sqrt(PNE(1,5)**2+PNE(2,5)**2+PNE(3,5)**2+XMPI2)
C     $      +sqrt(PNE(1,3)**2+PNE(2,3)**2+PNE(3,3)**2+XMMU**2))
C      write(*,*) "INCOMING s=",
C     $     (E+sqrt(PNE(1,2)**2+PNE(2,2)**2+PNE(3,2)**2+XMN2))**2
C     $    -( (PNE(1,1)+PNE(1,2))**2+(PNE(2,1)+PNE(2,2))**2
C     $      +(PNE(3,1)+PNE(3,2))**2)
C      write(*,*) "OUTGOING s=",
C     $     ( sqrt(PNE(1,4)**2+PNE(2,4)**2+PNE(3,4)**2+XMN2)
C     $      +sqrt(PNE(1,5)**2+PNE(2,5)**2+PNE(3,5)**2+XMPI2)
C     $      +sqrt(PNE(1,3)**2+PNE(2,3)**2+PNE(3,3)**2+XMMU**2))**2
C     $    -( (PNE(1,4)+PNE(1,5)+PNE(1,3))**2
C     $      +(PNE(2,4)+PNE(2,5)+PNE(2,3))**2
C     $      +(PNE(3,4)+PNE(3,5)+PNE(3,3))**2)
C------------------------------------------------------------------

      RETURN

 999  WRITE(*,*)
     $ "RSSETVCT:*** Invalid JMOD(Interaction Code):JMOD=",JMOD," ***"
      STOP

      END
