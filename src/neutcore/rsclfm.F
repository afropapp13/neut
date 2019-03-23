***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE RSCLFM(IFLAG,IBLOCK,PAR,F)
*     ---------------------------------------------------
*
*     ( purpose )
*       calculate helicity amplitude
*
*     ( input )
*       iflag       : interaction flag
*       iblock      : resonance-code
*       par(14)     : coefficients(FKR modep parameter)
*
*     ( output )
*       f(6)        : helicity amplitude
*
*     ( creation date and author )
*       1993.Dec. ; for Super-Kamioka by Y.Hayato
*                   (Based on Rein and Sehgal's code)
*
*     ( comment )
*       Reference:Rein and Sehgal,Ann. of Phys. 133(1981),79-153
*                 Rein           ,Z.Phys.C 35(1987),43-64
*                 Feynman et al. ,Phys.Rev.D3(1971),2706-
***********************************************************************

      IMPLICIT NONE

      INTEGER IFLAG
      INTEGER IBLOCK
      REAL PAR(14)
      REAL F(6)

C     CONSTANTS

#include<rscons.h>

C      COMMON /RSCONS/Z,PI,MN,MN2,MPI,MPI2,MV,MV2,MA,MA2,OMEG,XW,EPSI
C      REAL  Z,PI,MN,MN2,MPI,MPI2,MV,MV2,MA,MA2,OMEG,XW,EPSI

c      REAL XW

      REAL XLAM,XLAM2
      REAL R,RA,RV,RM,RP,T,TA,TV,TP,TM,S,B,C
      REAL CO1,CO2,CO3,CO4,CO5,CO6,CO7,CO8
         
      RV  = PAR( 1)
      RA  = PAR( 2)
      R   = PAR( 3)
      RP  = PAR( 4)
      RM  = PAR( 5)
      TV  = PAR( 6)
      TA  = PAR( 7)
      T   = PAR( 8)
      TP  = PAR( 9)
      TM  = PAR(10)
      S   = PAR(11)
      B   = PAR(12)
      C   = PAR(13)
      XLAM= PAR(14)

      XLAM2= XLAM**2

      IF(IBLOCK.EQ. 0)  GOTO 500
      IF(IBLOCK.EQ. 1)  GOTO 510
      IF(IBLOCK.EQ. 2)  GOTO 520
      IF(IBLOCK.EQ. 3)  GOTO 530
      IF(IBLOCK.EQ. 4)  GOTO 540
      IF(IBLOCK.EQ. 5)  GOTO 550
      IF(IBLOCK.EQ. 6)  GOTO 560
      IF(IBLOCK.EQ. 7)  GOTO 570
      IF(IBLOCK.EQ. 8)  GOTO 580
      IF(IBLOCK.EQ. 9)  GOTO 590
      IF(IBLOCK.EQ.10)  GOTO 600
      IF(IBLOCK.EQ.11)  GOTO 610
      IF(IBLOCK.EQ.12)  GOTO 620
      IF(IBLOCK.EQ.13)  GOTO 630
      IF(IBLOCK.EQ.14)  GOTO 640
      IF(IBLOCK.EQ.15)  GOTO 650
      IF(IBLOCK.EQ.16)  GOTO 660
      IF(IBLOCK.EQ.17)  GOTO 670
      IF(IBLOCK.EQ.18)  GOTO 680
      IF(IBLOCK.EQ.19)  GOTO 690
      IF(IBLOCK.EQ.20)  GOTO 700
      IF(IBLOCK.EQ.21)  GOTO 710
      IF(IBLOCK.EQ.22)  GOTO 720
      IF(IBLOCK.EQ.23)  GOTO 730
      IF(IBLOCK.EQ.24)  GOTO 740
      IF(IBLOCK.EQ.25)  GOTO 750
      IF(IBLOCK.EQ.26)  GOTO 760
      IF(IBLOCK.EQ.27)  GOTO 770
      IF(IBLOCK.EQ.28)  GOTO 780
      IF(IBLOCK.EQ.29)  GOTO 790
      IF(IBLOCK.EQ.30)  GOTO 800
      
 500  CONTINUE
      
C     **********************************
C     IBLOCK = 0      P33(1234)P,N
C     **********************************
      IF(IFLAG.LE.3)  GOTO 501
      F(1) =-SQRT(6.)*(RM+2.*XW*R)
      F(2)=-SQRT(2.)*(RM+2.*XW*R)
      F(4)=+SQRT(2.)*(RP+2.*XW*R)
      F(3)=+SQRT(6.)*(RP+2.*XW*R)
      F(5)= 2.*SQRT(2.)*C
      F(6)= F(5)
      GOTO 2000
 501  CONTINUE
      F(1) = SQRT(6.)*RM
      F(2) = SQRT(2.)*RM
      F(4) = -SQRT(2.)*RP
      F(3) = -SQRT(6.)*RP
      F(5) = -2.*SQRT(2.)*C
      F(6) = F(5)
      GOTO 2000
 510  CONTINUE
C     **********************************
C     IBLOCK = 1       D13(1525)P
C     **********************************
      IF(IFLAG.LE.3)  GOTO 511
      F(1) =-SQRT(9./2.)*(TM+2.*XW*T)
      F(2)=-SQRT(3./2.)*(TM+2.*XW*T)+SQRT(4./3.)*XLAM*(RM
     $     +3.*XW*R)
      F(4)=-SQRT(3./2.)*(TP+2.*XW*T)+SQRT(4./3.)*XLAM*(RP
     $     +3.*XW*R)
      F(3)=-SQRT(9./2.)*(TP+2.*XW*T)
      F(5)= +SQRT(3.)*XLAM*S*(1.-2.*XW)+SQRT(4./3.)*XLAM*C
      F(6)= +SQRT(3.)*XLAM*S*(1.-2.*XW)-SQRT(4./3.)*XLAM*C
      GOTO 2000
 511  CONTINUE
      F(1) = 0.
      F(2) = 0.
      F(4) = 0.
      F(3) = 0.
      F(5) = 0.
      F(6) = 0.
      GOTO 2000
 520  CONTINUE
C     **********************************
C     IBLOCK = 2       D13(1525)N
C     **********************************
      IF(IFLAG.LE.3)  GOTO 521
      F(1) = +SQRT(9./2.)*(TM+2.*XW*T)
      F(2)= +SQRT(3./2.)*(TM+2.*XW*T)-SQRT(4./3.)*XLAM*(RM
     $     +1.*XW*R)
      F(4)= +SQRT(3./2.)*(TP+2.*XW*T)-SQRT(4./3.)*XLAM*(RP
     $     +1.*XW*R)
      F(3)= +SQRT(9./2.)*(TP+2.*XW*T)
      F(5)=-SQRT(3.)*XLAM*S*(1.-2.*XW)-SQRT(4./3.)*XLAM*C
      F(6)=-SQRT(3.)*XLAM*S*(1.-2.*XW)+SQRT(4./3.)*XLAM*C
      GOTO 2000
 521  CONTINUE
      F(1) =-2.*SQRT(9./2.)*TM
      F(2) =-SQRT(6.)*TM+SQRT(16./3.)*XLAM*RM
      F(4) =-SQRT(6.)*TP+SQRT(16./3.)*XLAM*RP
      F(3) =-2.*SQRT(9./2.)*TP
      F(5) = +SQRT(12.)*XLAM*S+SQRT(16./3.)*XLAM*C
      F(6) = +SQRT(12.)*XLAM*S-SQRT(16./3.)*XLAM*C
      GOTO 2000
 530  CONTINUE
C     **********************************
C     IBLOCK = 3       S11(1540)P
C     **********************************
      IF(IFLAG.LE.3)  GOTO 531
      F(1) = 0.
      F(2)=-SQRT(3.)*(TM+2.*XW*T)-SQRT(2./3.)*XLAM*(RM+3.
     $     *XW*R)
      F(4)=+SQRT(3.)*(TP+2.*XW*T)+SQRT(2./3.)*XLAM*(RP+3.
     $     *XW*R)
      F(3)= 0.
      F(5)=-SQRT(3./2.)*XLAM*S*(1.-2.*XW)-SQRT(2./3.)*(XLAM
     $     *C-3.*B)
      F(6)=+SQRT(3./2.)*XLAM*S*(1.-2.*XW)-SQRT(2./3.)*(XLAM
     $     *C-3.*B)
      GOTO 2000
 531  CONTINUE
      F(1) = 0.
      F(2) = 0.
      F(4) = 0.
      F(3) = 0.
      F(5) = 0.
      F(6) = 0.
      GOTO 2000
 540  CONTINUE
C     **********************************
C     IBLOCK = 4       S11(1540)N
C     **********************************
      IF(IFLAG.LE.3)  GOTO 541
      F(1) = 0.
      F(2)=+SQRT(3.)*(TM+2.*XW*T)+SQRT(2./3.)*XLAM*(RM+1.
     $     *XW*R)
      F(4)=-SQRT(3.)*(TP+2.*XW*T)-SQRT(2./3.)*XLAM*(RP+1.
     $     *XW*R)
      F(3)= 0.
      F(5)=+SQRT(3./2.)*XLAM*S*(1.-2.*XW)+SQRT(2./3.)*(XLAM
     $     *C-3.*B)
      F(6)=-SQRT(3./2.)*XLAM*S*(1.-2.*XW)+SQRT(2./3.)*(XLAM
     $     *C-3.*B)
      GOTO 2000
 541  CONTINUE
      F(1) = 0.
      F(2) =-SQRT(12.)*TM-SQRT(16./6.)*XLAM*RM
      F(4) = +SQRT(12.)*TP+SQRT(16./6.)*XLAM*RP
      F(5) =-SQRT(6.)*XLAM*S-SQRT(8./3.)*(XLAM*C-3.*B)
      F(6) = +SQRT(6.)*XLAM*S-SQRT(8./3.)*(XLAM*C-3.*B)
      F(3) = 0.
      GOTO 2000
 550  CONTINUE
C     **********************************
C     IBLOCK = 5       S31(1620)P,N        P=SQRT(3)*N
C     **********************************
      IF(IFLAG.LE.3)  GOTO 551
      F(1) = 0.
      F(2)= SQRT(3.)*(TM+2.*XW*T)-SQRT(1./6.)*XLAM*(RM+2.
     $     *XW*R)
      F(4)=-SQRT(3.)*(TP+2.*XW*T)+SQRT(1./6.)*XLAM*(RP+2.
     $     *XW*R)
      F(3)= 0.
      F(5)= SQRT(3./2.)*XLAM*S*(1.-2.*XW)-SQRT(1./6.)*(XLAM
     $     *C-3.*B)
      F(6)=-SQRT(3./2.)*XLAM*S*(1.-2.*XW)-SQRT(1./6.)*(XLAM
     $     *C-3.*B)
      GOTO 2000
 551  CONTINUE
      F(1) = 0.
      F(2) = -SQRT(3.)*TM+SQRT(1./6.)*XLAM*RM
      F(4) = SQRT(3.)*TP-SQRT(1./6.)*XLAM*RP
      F(3) = 0.
      F(5) = -SQRT(3./2.)*XLAM*S+SQRT(1./6.)*(XLAM*C-3.*B)
      F(6) = SQRT(3./2.)*XLAM*S+SQRT(1./6.)*(XLAM*C-3.*B)
      GOTO 2000
 560  CONTINUE
C     **********************************
C     IBLOCK = 6       S11(1640)P
C     **********************************
      IF(IFLAG.LE.3)  GOTO 561
      F(1) = 0.
      F(2)= SQRT(1./24.)*XLAM*RM
      F(4)=-SQRT(1./24.)*XLAM*RP
      F(3)= 0.
      F(5)= -SQRT(1./6.)*(XLAM*C-3.*B)
      F(6)= -SQRT(1./6.)*(XLAM*C-3.*B)
      GOTO 2000
 561  CONTINUE
      F(1) = 0.
      F(2) = 0.
      F(4) = 0.
      F(3) = 0.
      F(5) = 0.
      F(6) = 0.
      GOTO 2000
 570  CONTINUE
C     **********************************
C     IBLOCK = 7       S11(1640)N
C     **********************************
      IF(IFLAG.LE.3)  GOTO 571
      F(1) = 0.
      F(2)= -SQRT(1./24.)*XLAM*(RM+4.*XW*R)
      F(4)= SQRT(1./24.)*XLAM*(RP+4.*XW*R)
      F(3)= 0.
      F(5)= SQRT(1./6.)*(XLAM*C-3.*B)
      F(6)= SQRT(1./6.)*(XLAM*C-3.*B)
      GOTO 2000
 571  CONTINUE
      F(1) = 0.
      F(2) = SQRT(1./6.)*XLAM*RM
      F(4) = -SQRT(1./6.)*XLAM*RP
      F(3) = 0.
      F(5) = -SQRT(2./3.)*(XLAM*C-3.*B)
      F(6) = -SQRT(2./3.)*(XLAM*C-3.*B)
      GOTO 2000
 580  CONTINUE
C     **********************************
C     IBLOCK = 8       D13(1670)P
C     **********************************
      IF(IFLAG.LE.3)  GOTO 581
      F(1) =-SQRT(9./40.)*XLAM*RM
      F(2)=-SQRT(1./120.)*XLAM*RM
      F(4)=-SQRT(1./120.)*XLAM*RP
      F(3)=-SQRT(9./40.)*XLAM*RP
      F(5)= +SQRT(1./30.)*XLAM*C
      F(6)=-SQRT(1./30.)*XLAM*C
      GOTO 2000
 581  CONTINUE
      F(1) = 0.
      F(2) = 0.
      F(4) = 0.
      F(3) = 0.
      F(5) = 0.
      F(6) = 0.
      GOTO 2000
 590  CONTINUE
C     **********************************
C     IBLOCK = 9       D13(1670)N
C     **********************************
      IF(IFLAG.LE.3)  GOTO 591
      F(1) = +SQRT(9./40.)*XLAM*(RM+4.*XW*R)
      F(2)=+SQRT(1./120.)*XLAM*(RM+4.*XW*R)
      F(4)=+SQRT(1./120.)*XLAM*(RP+4.*XW*R)
      F(3)= +SQRT(9./40.)*XLAM*(RP+4.*XW*R)
      F(5)=-SQRT(1./30.)*XLAM*C
      F(6)=+SQRT(1./30.)*XLAM*C
      GOTO 2000
 591  CONTINUE
      F(1) =-SQRT(9./10.)*XLAM*RM
      F(2) =-SQRT(1./30.)*XLAM*RM
      F(4) =-SQRT(1./30.)*XLAM*RP
      F(3) =-SQRT(9./10.)*XLAM*RP
      F(5) = +SQRT(2./15.)*XLAM*C
      F(6) =-SQRT(2./15.)*XLAM*C
      GOTO 2000
 600  CONTINUE
C     **********************************
C     IBLOCK = 10      D15(1680)P
C     **********************************
      IF(IFLAG.LE.3)  GOTO 601
      F(1) = -SQRT(3./20.)*XLAM*RM
      F(2)= -SQRT(3./40.)*XLAM*RM
      F(4)= SQRT(3./40.)*XLAM*RP
      F(3)= SQRT(3./20.)*XLAM*RP
      F(5)= SQRT(3./10.)*XLAM*C
      F(6)= SQRT(3./10.)*XLAM*C
      GOTO 2000
 601  CONTINUE
      F(1) = 0.
      F(2) = 0.
      F(4) = 0.
      F(3) = 0.
      F(5) = 0.
      F(6) = 0.
      GOTO 2000
 610  CONTINUE
C     **********************************
C     IBLOCK = 11      D15(1680)N
C     **********************************
      IF(IFLAG.LE.3)  GOTO 611
      F(1) = SQRT(3./20.)*XLAM*(RM+4.*XW*R)
      F(2)= SQRT(3./40.)*XLAM*(RM+4.*XW*R)
      F(4)= -SQRT(3./40.)*XLAM*(RP+4.*XW*R)
      F(3)= -SQRT(3./20.)*XLAM*(RP+4.*XW*R)
      F(5)= -SQRT(3./10.)*XLAM*C
      F(6)= -SQRT(3./10.)*XLAM*C
      GOTO 2000
 611  CONTINUE
      F(1) = -SQRT(3./5.)*XLAM*RM
      F(2) = -SQRT(3./10.)*XLAM*RM
      F(4) = SQRT(3./10.)*XLAM*RP
      F(3) = SQRT(3./5.)*XLAM*RP
      F(5) = SQRT(6./5.)*XLAM*C
      F(6) = SQRT(6./5.)*XLAM*C
      GOTO 2000
 620  CONTINUE
C     **********************************
C     IBLOCK = 12      D33(1730)P,N        P=SQRT(3)*N
C     **********************************
      IF(IFLAG.LE.3)  GOTO 621
      F(1) = SQRT(9./2.)*(TM+2.*XW*T)
      F(2)= SQRT(3./2.)*(TM+2.*XW*T)+SQRT(1./3.)*XLAM*(RM
     $     +2.*XW*R)
      F(4)= SQRT(3./2.)*(TP+2.*XW*T)+SQRT(1./3.)*XLAM*(RP
     $     +2.*XW*R)
      F(3)= SQRT(9./2.)*(TP+2.*XW*T)
      F(5)= -SQRT(3.)*XLAM*S*(1.-2.*XW)+SQRT(1./3.)*XLAM*C
      F(6)= -SQRT(3.)*XLAM*S*(1.-2.*XW)-SQRT(1./3.)*XLAM*C
      GOTO 2000
 621  CONTINUE
      F(1) = -SQRT(9./2.)*TM
      F(2) = -SQRT(3./2.)*TM-SQRT(1./3.)*XLAM*RM
      F(4) = -SQRT(3./2.)*TP-SQRT(1./3.)*XLAM*RP
      F(3) = -SQRT(9./2.)*TP
      F(5) = SQRT(3.)*XLAM*S-SQRT(1./3.)*XLAM*C
      F(6) = SQRT(3.)*XLAM*S+SQRT(1./3.)*XLAM*C
      GOTO 2000
 630  CONTINUE
C     **********************************
C     IBLOCK = 13     P11(1450)P
C     **********************************
      IF(IFLAG.LE.3)  GOTO 631
      F(1) = 0.
      F(2)= -(5./12.)*SQRT(3.)*XLAM2*(RM+2.*XW*(6./5.)*R)
      F(4)= -(5./12.)*SQRT(3.)*XLAM2*(RP+2.*XW*(6./5.)*R)
      F(3)= 0.
      F(5)= -SQRT(3./16.)*XLAM2*S*(1.-4.*XW)-(5./12.)
     $     *SQRT(3.)*XLAM*(XLAM*C-2.*B)
      F(6)= -SQRT(3./16.)*XLAM2*S*(1.-4.*XW)+(5./12.)
     $     *SQRT(3.)*XLAM*(XLAM*C-2.*B)
      GOTO 2000
 631  CONTINUE
      F(1) = 0.
      F(2) = 0.
      F(4) = 0.
      F(3) = 0.
      F(5) = 0.
      F(6) = 0.
      GOTO 2000
 640  CONTINUE
C     **********************************
C     IBLOCK = 14     P11(1450)N
C     **********************************
      IF(IFLAG.LE.3)  GOTO 641
      F(1) = 0.
      F(2)= (5./12.)*SQRT(3.)*XLAM2*(RM+2.*XW*(4./5.)*R)
      F(4)= (5./12.)*SQRT(3.)*XLAM2*(RP+2.*XW*(4./5.)*R)
      F(3)= 0.
      F(5)= +SQRT(3./16.)*XLAM2*S+(5./12.)*SQRT(3.)*XLAM
     $     *(XLAM*C-2.*B)
      F(6)= +SQRT(3./16.)*XLAM2*S-(5./12.)*SQRT(3.)*XLAM
     $     *(XLAM*C-2.*B)
      GOTO 2000
 641  CONTINUE
      F(1) = 0.
      F(2) = -(5./6.)*SQRT(3.)*XLAM2*RM
      F(4) = -(5./6.)*SQRT(3.)*XLAM2*RP
      F(3) = 0.
      F(5) = -SQRT(3./4.)*XLAM2*S-(5./6.)*SQRT(3.)*XLAM
     $     *(XLAM*C-2.*B)
      F(6) = -SQRT(3./4.)*XLAM2*S+(5./6.)*SQRT(3.)*XLAM
     $     *(XLAM*C-2.*B)
      GOTO 2000
 650  CONTINUE
C     **********************************
C     IBLOCK = 15     P33(1640)P,N
C     **********************************
      IF(IFLAG.LE.3)  GOTO 651
      F(1) = SQRT(1./2.)*XLAM2*(RM+2.*XW*R)
      F(2)= F(1)/SQRT(3.)
      F(4)= -SQRT(1./6.)*XLAM2*(RP+2.*XW*R)
      F(3)= SQRT(3.)*F(4)
      F(5)= -SQRT(2./3.)*XLAM*(XLAM*C-2.*B)
      F(6)= F(5)
      GOTO 2000
 651  CONTINUE
      F(1) = -SQRT(1./2.)*XLAM2*RM
      F(2) = -SQRT(1./6.)*XLAM2*RM
      F(4) = SQRT(1./6.)*XLAM2*RP
      F(3) = SQRT(1./2.)*XLAM2*RP
      F(5) = SQRT(2./3.)*XLAM*(XLAM*C-2.*B)
      F(6) = F(5)
      GOTO 2000
 660  CONTINUE
C     **********************************
C     IBLOCK = 16     F15(1680)P
C     **********************************
      IF(IFLAG.LE.3)  GOTO 661
      F(1) = -SQRT(9./10.)*XLAM*(TM+4.*XW*T)
      F(2)= -SQRT(9./20.)*XLAM*(TM+4.*XW*T)
     1     +SQRT(5./8.)*XLAM2*(RM+2.*XW*(6./5.)*R)
      F(4)= -SQRT(9./20.)*XLAM*(TP+4.*XW*T)
     1     +SQRT(5./8.)*XLAM2*(RP+2.*XW*(6./5.)*R)
      F(3)= -SQRT(9./10.)*XLAM*(TP+4.*XW*T)
      F(5)= SQRT(9./40.)*XLAM2*S*(1.-4.*XW)+SQRT(5./8.)
     $     *XLAM2*C
      F(6)= SQRT(9./40.)*XLAM2*S*(1.-4.*XW)-SQRT(5./8.)
     $     *XLAM2*C
      GOTO 2000
 661  CONTINUE
      F(1) = 0.
      F(2) = 0.
      F(4) = 0.
      F(3) = 0.
      F(5) = 0.
      F(6) = 0.
      GOTO 2000
 670  CONTINUE
C     **********************************
C     IBLOCK = 17     F15(1680)N
C     **********************************
      IF(IFLAG.LE.3)  GOTO 671
      F(1) = SQRT(9./10.)*XLAM*TM
      F(2)= SQRT(9./20.)*XLAM*TM-SQRT(5./8.)*XLAM2*(RM+2.
     $     *XW*(4./5.)*R)
      F(4)= SQRT(9./20.)*XLAM*TP-SQRT(5./8.)*XLAM2*(RP+2.
     $     *XW*(4./5.)*R)
      F(3)= SQRT(9./10.)*XLAM*TP
      F(5)= -SQRT(9./40.)*XLAM2*S-SQRT(5./8.)*XLAM2*C
      F(6)= -SQRT(9./40.)*XLAM2*S+SQRT(5./8.)*XLAM2*C
      GOTO 2000
 671  CONTINUE
      F(1) = -SQRT(18./5.)*XLAM*TM
      F(2) = -SQRT(9./5.)*XLAM*TM+SQRT(5./2.)*XLAM2*RM
      F(4) = -SQRT(9./5.)*XLAM*TP+SQRT(5./2.)*XLAM2*RP
      F(3) = -SQRT(18./5.)*XLAM*TP
      F(5) = SQRT(9./10.)*XLAM2*S+SQRT(5./2.)*XLAM2*C
      F(6) = SQRT(9./10.)*XLAM2*S-SQRT(5./2.)*XLAM2*C
      GOTO 2000
 680  CONTINUE
C     **********************************
C     IBLOCK = 18     P11(1710)P
C     **********************************
      IF(IFLAG.LE.3)  GOTO 681
      F(1) = 0.
      F(2)=-SQRT(1./6.)*XLAM2*(RM+3.*XW*R)
      F(4)=-SQRT(1./6.)*XLAM2*(RP+3.*XW*R)
      F(3)= 0.
      F(5)=-SQRT(3./8.)*XLAM2*S*(1.-2.*XW)-SQRT(1./6.)*XLAM
     $     *(XLAM*C-2.*B)
      F(6)=-SQRT(3./8.)*XLAM2*S*(1.-2.*XW)+SQRT(1./6.)*XLAM
     $     *(XLAM*C-2.*B)
      GOTO 2000
 681  CONTINUE
      F(1) = 0.
      F(2) = 0.
      F(4) = 0.
      F(3) = 0.
      F(5) = 0.
      F(6) = 0.
      GOTO 2000
 690  CONTINUE
C     **********************************
C     IBLOCK = 19     P11(1710)N
C     **********************************
      IF(IFLAG.LE.3)  GOTO 691
      F(1) = 0.
      F(2)= +SQRT(1./6.)*XLAM2*(RM+1.*XW*R)
      F(4)= +SQRT(1./6.)*XLAM2*(RP+1.*XW*R)
      F(3)= 0.
      F(5)= +SQRT(3./8.)*XLAM2*S*(1.-2.*XW)+SQRT(1./6.)
     $     *XLAM*(XLAM*C-2.*B)
      F(6)= +SQRT(3./8.)*XLAM2*S*(1.-2.*XW)-SQRT(1./6.)
     $     *XLAM*(XLAM*C-2.*B)
      GOTO 2000
 691  CONTINUE
      F(1) = 0.
      F(2) =-SQRT(2./3.)*XLAM2*RM
      F(4) =-SQRT(2./3.)*XLAM2*RP
      F(3) = 0.
      F(5) =-SQRT(3./2.)*XLAM2*S-SQRT(2./3.)*XLAM*(XLAM*C-2.
     $     *B)
      F(6) =-SQRT(3./2.)*XLAM2*S+SQRT(2./3.)*XLAM*(XLAM*C-2.
     $     *B)
      GOTO 2000
 700  CONTINUE
C     **********************************
C     IBLOCK = 20     P13(1740)P
C     **********************************
      IF(IFLAG.LE.3)  GOTO 701
      F(1) = SQRT(9./40.)*XLAM*(TM+4.*XW*T)
      F(2)= -SQRT(27./40.)*XLAM*(TM+4.*XW*T)-SQRT(5./12.)
     $     *XLAM2*(RM+2.*XW*(6./5.)*R)
      F(4)= +SQRT(27./40.)*XLAM*(TP+4.*XW*T)+SQRT(5./12.)
     $     *XLAM2*(RP+2.*XW*(6./5.)*R)
      F(3)= -SQRT(9./40.)*XLAM*(TP+4.*XW*T)
      F(5)=-SQRT(3./20.)*XLAM2*S*(1.-4.*XW)-SQRT(5./12.)
     $     *XLAM*(XLAM*C-5.*B)
      F(6)= SQRT(3./20.)*XLAM2*S*(1.-4.*XW)-SQRT(5./12.)
     $     *XLAM*(XLAM*C-5.*B)
      GOTO 2000
 701  CONTINUE
      F(1) = 0.
      F(2) = 0.
      F(4) = 0.
      F(3) = 0.
      F(5) = 0.
      F(6) = 0.
      GOTO 2000
 710  CONTINUE
C     **********************************
C     IBLOCK = 21     P13(1740)N
C     **********************************
      IF(IFLAG.LE.3)  GOTO 711
      F(1) = -SQRT(9./40.)*XLAM*TM
      F(2)= SQRT(27./40.)*XLAM*TM+SQRT(5./12.)*XLAM2*(RP+2.
     $     *XW*(4./5.)*R)
      F(4)= -SQRT(27./40.)*XLAM*TP-SQRT(5./12.)*XLAM2*(RP+2
     $     .*XW*(4./5.)*R)
      F(3)= SQRT(9./40.)*XLAM*TP
      F(5)= SQRT(3./20.)*XLAM2*S+SQRT(5./12.)*XLAM*(XLAM*C-5
     $     .*B)
      F(6)=-SQRT(3./20.)*XLAM2*S+SQRT(5./12.)*XLAM*(XLAM*C-5
     $     .*B)
      GOTO 2000
 711  CONTINUE
      F(1) = SQRT(9./10.)*XLAM*TM
      F(2) = -SQRT(27./10.)*XLAM*TM-SQRT(5./3.)*XLAM2*RM
      F(4) = SQRT(27./10.)*XLAM*TP+SQRT(5./3.)*XLAM2*RP
      F(3) = -SQRT(9./10.)*XLAM*TP
      F(5) = -SQRT(3./5.)*XLAM2*S-SQRT(5./3.)*XLAM*(XLAM*C-5
     $     .*B)
      F(6) = SQRT(3./5.)*XLAM2*S-SQRT(5./3.)*XLAM*(XLAM*C-5.
     $     *B)
      GOTO 2000
 720  CONTINUE
C     **********************************
C     IBLOCK = 22     P31(1920)P,N
C     **********************************
      IF(IFLAG.LE.3)  GOTO 721
      F(1) = 0.
      F(2)= +SQRT(1./15.)*XLAM2*(RM+2.*XW*R)
      F(4)= +SQRT(1./15.)*XLAM2*(RP+2.*XW*R)
      F(3)= 0.
      F(5)=-SQRT(4./15.)*XLAM*(XLAM*C-5.*B)
      F(6)= +SQRT(4./15.)*XLAM*(XLAM*C-5.*B)
      GOTO 2000
 721  CONTINUE
      F(1) = 0.
      F(2) =-SQRT(1./15.)*XLAM2*RM
      F(4) =-SQRT(1./15.)*XLAM2*RP
      F(3) = 0.
      F(5) = +SQRT(4./15.)*XLAM*(XLAM*C-5.*B)
      F(6) = -F(5)
      GOTO 2000
 730  CONTINUE
C     **********************************
C     IBLOCK = 23     F35(1920)P,N
C     **********************************
      IF(IFLAG.LE.3)  GOTO 731
      F(1) =-SQRT(18./35.)*XLAM2*(RM+2.*XW*R)
      F(2)=-SQRT(1./35.)*XLAM2*(RM+2.*XW*R)
      F(4)=-SQRT(1./35.)*XLAM2*(RP+2.*XW*R)
      F(3)=-SQRT(18./35.)*XLAM2*(RP+2.*XW*R)
      F(5)= +SQRT(4./35.)*XLAM2*C
      F(6)= -F(5)
      GOTO 2000
 731  CONTINUE
      F(1) = +SQRT(18./35.)*XLAM2*RM
      F(2) = +SQRT(1./35.)*XLAM2*RM
      F(4) = +SQRT(1./35.)*XLAM2*RP
      F(3)= +SQRT(18./35.)*XLAM2*RP
      F(5) =-SQRT(4./35.)*XLAM2*C
      F(6) = -F(5)
      GOTO 2000
 740  CONTINUE
C     **********************************
C     IBLOCK = 24     F37(1950)P,N
C     **********************************
      IF(IFLAG.LE.3)  GOTO 741
      F(1) = -SQRT(2./7.)*XLAM2*(RM+2.*XW*R)
      F(2)= -SQRT(6./35.)*XLAM2*(RM+2.*XW*R)
      F(4)= SQRT(6./35.)*XLAM2*(RP+2.*XW*R)
      F(3)= SQRT(2./7.)*XLAM2*(RP+2.*XW*R)
      F(5)= SQRT(24./35.)*XLAM2*C
      F(6)= F(5)
      GOTO 2000
 741  CONTINUE
      F(1) = SQRT(2./7.)*XLAM2*RM
      F(2) = SQRT(6./35.)*XLAM2*RM
      F(4) = -SQRT(6./35.)*XLAM2*RP
      F(3) = -SQRT(2./7.)*XLAM2*RP
      F(5) = -SQRT(24./35.)*XLAM2*C
      F(6) = F(5)
      GOTO 2000
 750  CONTINUE
C     **********************************
C     IBLOCK = 25     P33(1960)P,N
C     **********************************
      IF(IFLAG.LE.3)  GOTO 751
      F(1) = -SQRT(1./5.)*XLAM2*(RM+2.*XW*R)
      F(2)= SQRT(1./15.)*XLAM2*(RM+2.*XW*R)
      F(4)= -SQRT(1./15.)*XLAM2*(RP+2.*XW*R)
      F(3)= SQRT(1./5.)*XLAM2*(RP+2.*XW*R)
      F(5)= -SQRT(4./15.)*XLAM*(XLAM*C-5.*B)
      F(6)= +F(5)
      GOTO 2000
 751  CONTINUE
      F(1) = SQRT(1./5.)*XLAM2*RM
      F(2) = -SQRT(1./15.)*XLAM2*RM
      F(4) = SQRT(1./15.)*XLAM2*RP
      F(3) = -SQRT(1./5.)*XLAM2*RP
      F(5) = SQRT(4./15.)*XLAM*(XLAM*C-5.*B)
      F(6) = F(5)
      GOTO 2000
 760  CONTINUE
C     **********************************
C     IBLOCK = 26     F17(1970)P
C     **********************************
      IF(IFLAG.LE.3)  GOTO 761
      F(1) = -SQRT(1./28.)*XLAM2*RM
      F(2)= -SQRT(3./140.)*XLAM2*RM
      F(4)= SQRT(3./140.)*XLAM2*RP
      F(3)= SQRT(1./28.)*XLAM2*RP
      F(5)= SQRT(3./70.)*XLAM2*C
      F(6)= F(5)
      GOTO 2000
 761  CONTINUE
      F(1) = 0.
      F(2) = 0.
      F(4) = 0.
      F(3) = 0.
      F(5) = 0.
      F(6) = 0.
      GOTO 2000
 770  CONTINUE
C     **********************************
C     IBLOCK = 27     F17(1970)N
C     **********************************
      IF(IFLAG.LE.3)  GOTO 771
      F(1) = SQRT(1./28.)*XLAM2*(RM+4.*XW*R)
      F(2)= SQRT(3./140.)*XLAM2*(RM+4.*XW*R)
      F(4)= -SQRT(3./140.)*XLAM2*(RP+4.*XW*R)
      F(3)= -SQRT(1./28.)*XLAM2*(RP+4.*XW*R)
      F(5)= -SQRT(3./70.)*XLAM2*C
      F(6)= F(5)
      GOTO 2000
 771  CONTINUE
      F(1) = -SQRT(1./7.)*XLAM2*RM
      F(2) = -SQRT(3./35.)*XLAM2*RM
      F(4) = SQRT(3./35.)*XLAM2*RP
      F(3) = SQRT(1./7.)*XLAM2*RP
      F(5) = SQRT(6./35.)*XLAM2*C
      F(6) = F(5)
      GOTO 2000
 780  CONTINUE
C     **********************************
C     IBLOCK = 28   BACKGROUND   PROTON
C     **********************************
      IF(IFLAG.LE.3)  GOTO 781
      F(1) = 0.
      F(2)= (5.*(RM+2.*XW*(6./5.)*R)/2.)
      F(4)= (5.*(RP+2.*XW*(6./5.)*R)/2.)
      F(3)= 0.
      F(5)= (3.*S*(1.-4.*XW)+5.*C)/2.
      F(6)= (3.*S*(1.-4.*XW)-5.*C)/2.
      GOTO 2000
 781  CONTINUE
      F(1) = 0.
      F(2) = 0.
      F(4) = 0.
      F(3) = 0.
      F(5) = 0.
      F(6) = 0.
      GOTO 2000
 790  CONTINUE
C     **********************************
C     IBLOCK = 29   BACKGROUND   NEUTRON
C     **********************************
      IF(IFLAG.LE.3)  GOTO 791
      F(1) = 0.
      F(2)= -5.*(RM+2.*XW*(4./5.)*R)/2.
      F(4)= -5.*(RP+2.*XW*(4./5.)*R)/2.
      F(3)= 0.
      F(5)= -(3.*S+5.*C)/2.
      F(6)= -(3.*S-5.*C)/2.
      GOTO 2000
 791  CONTINUE
      F(1) = 0.
      F(2) = 5.*RM
      F(4) = 5.*RP
      F(3) = 0.
      F(5) = (3.*S+5.*C)
      F(6) = (3.*S-5.*C)
      GOTO 2000
 800  CONTINUE
C     **********************************
C     IBLOCK = 30  BACKGROUND  J=3/2
C     **********************************
      CO1 = 0.
      CO2 = 0.
      CO3 = 0.
      CO4 = 0.
      CO5 = 0.
      CO6 = 0.
      CO7 = 0.
      CO8 = 0.
C     J = 3/2 BACKGROUND NOW INOPERATIVE
      IF(IFLAG.LE.3) GOTO 801
      F(1) = CO1*RM
      F(2)= CO2*RM
      F(4)= CO3*RP
      F(3)= CO4*RP
      F(5)= 0.
      F(6)= 0.
      GOTO 2000
 801  CONTINUE
      F(1) = CO5*RM
      F(2) = CO6*RM
      F(4) = CO7*RP
      F(3) = CO8*RP
      F(5) = 0.
      F(6) = 0.
      
 2000 CONTINUE
      

      END
      
