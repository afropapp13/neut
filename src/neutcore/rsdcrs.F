***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE RSDCRS(IMODE,OriginalF,Lambda,XMLEP,E,Q2,W,DIFCRS,
     $           DFCRS3)
*     ---------------------------------------------------
*     ( purpose )
*       calculate Differential Cross-section d^2(sig)/dwdq2.
*
*     ( input )
*       imode : Interaction mode
*               imode > 10  mans anti-neutrino
*       e     : Incident neutrino energy in lab. frame.
*               (Lab. frame is nucleon rest frame.)
*       q2    : sqaured momentum transfer.(minus sign)
*       w     : Rest mass of products(pion and nucreon).
*     ( output )
*       difcrs: differential cross-section.
*               If error occurs in this subroutine,difcrs
*               returns -9999.999
*       dfcrs3: differential cross-section of DELTA++(1232).
*               If error occurs in this subroutine,difcrs
*               returns -9999.999
*     ( creation date and author )
*       1993.Dec. ; for Super-Kamioka by Y.Hayato
*                   (Based on Rein and Sehgal's code)
*
*       2007.Dec. ; add lepton mass effects by G.Mitsuka
*       2013.10    ; P.Rodrigues & A.Bercellieadd
*                              new form facor for Rein&Sehgal 1pi
*
*     ( comment )
*       Reference:Rein and Sehgal,Ann. of Phys. 133(1981),79-153
*                 Rein           ,Z.Phys.C 35(1987),43-64
*                 Feynman et al. ,Phys.Rev.D3(1971),2706-
*                 K.S.Kuzmin et al., hep-ph/0312107
*                 Berger and Sehgal, hep-ph/0709.4378
C                 K.M. Graczyk and J.T. Sobczyk, Phys.Rev.D 77,053001 (2008)
C                    Equations are correction in erratum
C                 J.A. Nowak, arXiv:0909.3659
C                    Vector Form Factors
C                    CV1, CV2 are for simplicity, not true form factors
C                    Added the proportional factors of the vector form factors
C                    GV to CV3,CV4,CV5
C
***********************************************************************

      IMPLICIT NONE

      INTEGER IMODE
      REAL E,Q2,W,DIFCRS,DFCRS3

#include"rscons.h"
#include <neutmodel.h>

      REAL*4 XMASPI2, XMVSPI2

      INTEGER IFLAG
      REAL  XNU,XNU2,XNUS,XLAM,XLAM2,XMR
      REAL  AM3(50),AM1(50),AP1(50),AP3(50),A0P(50),A0M(50)
      REAL  AP1C(50),AP3C(50),AM1C(50),AM3C(50),A0PC(50),A0MC(50)
      REAL  FFBW(50),ABWR(50),ABWI(50)


      INTEGER IPP(31),NRR(31)

      INTEGER IB,NR,IP,IBLOCK,NCTR
      REAL WMIN,QMAX,SMAX
      REAL W2,X,QV,QV2,QVS,QVS2,A,ES
      REAL U,V,U2,V2,AK
      REAL FAKT,FRP1,FRP3,FRM1,FRM3,FR0P,FR0M,FIP1,FIP3,FIM1
C      REAL FIM3,FI0M,FI0P,BR,FBW,FBWNO,ETA,XE,RV,FFV,R,TA,FFA,T
      REAL FIM3,FI0M,FI0P,BR,FBW,FBWNO,XE,RV,FFV,R,TA,FFA,T
      REAL TV,RA,RP,RM,TP,TM,S,B,C,FM3C,FP3C,FM1C,FP1C,F0PC
      REAL F0MC,FM3,FP3,FM1,FP1,F0P,F0M,XXE,RAD,PQ,PQR,RADR,BRA
      REAL CO,SL,SR,SS,SL1,SR1,SS1,SL0,SR0,SS0,DSI,DSI1,DSI0

      REAL GV3,GV1,NFF,CV1,CV2,CV3A,CV3B,CV4,CV5
      REAL CA5

      COMPLEX ETABW(31)

      COMPLEX ASUM(2,6,14)
      REAL AM(6,31)

      COMPLEX S31M3N,S11M3N,P31M3N,P11M3N,P33M3N,P13M3N,D33M3N,D13M3N
      COMPLEX D35M3N,D15M3N,F35M3N,F15M3N,F37M3N,F17M3N

      COMPLEX S31P3N,S11P3N,P31P3N,P11P3N,P33P3N,P13P3N,D33P3N,D13P3N
      COMPLEX D35P3N,D15P3N,F35P3N,F15P3N,F37P3N,F17P3N

      COMPLEX S31M1N,S11M1N,P31M1N,P11M1N,P33M1N,P13M1N,D33M1N,D13M1N
      COMPLEX D35M1N,D15M1N,F35M1N,F15M1N,F37M1N,F17M1N

      COMPLEX S31P1N,S11P1N,P31P1N,P11P1N,P33P1N,P13P1N,D33P1N,D13P1N
      COMPLEX D35P1N,D15P1N,F35P1N,F15P1N,F37P1N,F17P1N

      COMPLEX S310MN,S110MN,P310MN,P110MN,P330MN,P130MN,D330MN,D130MN
      COMPLEX D350MN,D150MN,F350MN,F150MN,F370MN,F170MN

      COMPLEX S310PN,S110PN,P310PN,P110PN,P330PN,P130PN,D330PN,D130PN
      COMPLEX D350PN,D150PN,F350PN,F150PN,F370PN,F170PN



      COMPLEX S31M3P,S11M3P,P31M3P,P11M3P,P33M3P,P13M3P,D33M3P,D13M3P
      COMPLEX D35M3P,D15M3P,F35M3P,F15M3P,F37M3P,F17M3P

      COMPLEX S31P3P,S11P3P,P31P3P,P11P3P,P33P3P,P13P3P,D33P3P,D13P3P
      COMPLEX D35P3P,D15P3P,F35P3P,F15P3P,F37P3P,F17P3P

      COMPLEX S31M1P,S11M1P,P31M1P,P11M1P,P33M1P,P13M1P,D33M1P,D13M1P
      COMPLEX D35M1P,D15M1P,F35M1P,F15M1P,F37M1P,F17M1P

      COMPLEX S31P1P,S11P1P,P31P1P,P11P1P,P33P1P,P13P1P,D33P1P,D13P1P
      COMPLEX D35P1P,D15P1P,F35P1P,F15P1P,F37P1P,F17P1P

      COMPLEX S310MP,S110MP,P310MP,P110MP,P330MP,P130MP,D330MP,D130MP
      COMPLEX D350MP,D150MP,F350MP,F150MP,F370MP,F170MP

      COMPLEX S310PP,S110PP,P310PP,P110PP,P330PP,P130PP,D330PP,D130PP
      COMPLEX D350PP,D150PP,F350PP,F150PP,F370PP,F170PP

      EQUIVALENCE ( ASUM( 1 , 1 , 1 ) ,  S31M3P )
      EQUIVALENCE ( ASUM( 1 , 1 , 2 ) ,  S11M3P )
      EQUIVALENCE ( ASUM( 1 , 1 , 3 ) ,  P31M3P )
      EQUIVALENCE ( ASUM( 1 , 1 , 4 ) ,  P11M3P )
      EQUIVALENCE ( ASUM( 1 , 1 , 5 ) ,  P33M3P )
      EQUIVALENCE ( ASUM( 1 , 1 , 6 ) ,  P13M3P )
      EQUIVALENCE ( ASUM( 1 , 1 , 7 ) ,  D33M3P )
      EQUIVALENCE ( ASUM( 1 , 1 , 8 ) ,  D13M3P )
      EQUIVALENCE ( ASUM( 1 , 1 , 9 ) ,  D35M3P )
      EQUIVALENCE ( ASUM( 1 , 1 , 10 ) ,  D15M3P )
      EQUIVALENCE ( ASUM( 1 , 1 , 11 ) ,  F35M3P )
      EQUIVALENCE ( ASUM( 1 , 1 , 12 ) ,  F15M3P )
      EQUIVALENCE ( ASUM( 1 , 1 , 13 ) ,  F37M3P )
      EQUIVALENCE ( ASUM( 1 , 1 , 14 ) ,  F17M3P )

      EQUIVALENCE ( ASUM( 1 , 2 , 1 ) ,  S31M1P )
      EQUIVALENCE ( ASUM( 1 , 2 , 2 ) ,  S11M1P )
      EQUIVALENCE ( ASUM( 1 , 2 , 3 ) ,  P31M1P )
      EQUIVALENCE ( ASUM( 1 , 2 , 4 ) ,  P11M1P )
      EQUIVALENCE ( ASUM( 1 , 2 , 5 ) ,  P33M1P )
      EQUIVALENCE ( ASUM( 1 , 2 , 6 ) ,  P13M1P )
      EQUIVALENCE ( ASUM( 1 , 2 , 7 ) ,  D33M1P )
      EQUIVALENCE ( ASUM( 1 , 2 , 8 ) ,  D13M1P )
      EQUIVALENCE ( ASUM( 1 , 2 , 9 ) ,  D35M1P )
      EQUIVALENCE ( ASUM( 1 , 2 , 10 ) ,  D15M1P )
      EQUIVALENCE ( ASUM( 1 , 2 , 11 ) ,  F35M1P )
      EQUIVALENCE ( ASUM( 1 , 2 , 12 ) ,  F15M1P )
      EQUIVALENCE ( ASUM( 1 , 2 , 13 ) ,  F37M1P )
      EQUIVALENCE ( ASUM( 1 , 2 , 14 ) ,  F17M1P )

      EQUIVALENCE ( ASUM( 1 , 3 , 1 ) ,  S31P3P )
      EQUIVALENCE ( ASUM( 1 , 3 , 2 ) ,  S11P3P )
      EQUIVALENCE ( ASUM( 1 , 3 , 3 ) ,  P31P3P )
      EQUIVALENCE ( ASUM( 1 , 3 , 4 ) ,  P11P3P )
      EQUIVALENCE ( ASUM( 1 , 3 , 5 ) ,  P33P3P )
      EQUIVALENCE ( ASUM( 1 , 3 , 6 ) ,  P13P3P )
      EQUIVALENCE ( ASUM( 1 , 3 , 7 ) ,  D33P3P )
      EQUIVALENCE ( ASUM( 1 , 3 , 8 ) ,  D13P3P )
      EQUIVALENCE ( ASUM( 1 , 3 , 9 ) ,  D35P3P )
      EQUIVALENCE ( ASUM( 1 , 3 , 10 ) ,  D15P3P )
      EQUIVALENCE ( ASUM( 1 , 3 , 11 ) ,  F35P3P )
      EQUIVALENCE ( ASUM( 1 , 3 , 12 ) ,  F15P3P )
      EQUIVALENCE ( ASUM( 1 , 3 , 13 ) ,  F37P3P )
      EQUIVALENCE ( ASUM( 1 , 3 , 14 ) ,  F17P3P )

      EQUIVALENCE ( ASUM( 1 , 4 , 1 ) ,  S31P1P )
      EQUIVALENCE ( ASUM( 1 , 4 , 2 ) ,  S11P1P )
      EQUIVALENCE ( ASUM( 1 , 4 , 3 ) ,  P31P1P )
      EQUIVALENCE ( ASUM( 1 , 4 , 4 ) ,  P11P1P )
      EQUIVALENCE ( ASUM( 1 , 4 , 5 ) ,  P33P1P )
      EQUIVALENCE ( ASUM( 1 , 4 , 6 ) ,  P13P1P )
      EQUIVALENCE ( ASUM( 1 , 4 , 7 ) ,  D33P1P )
      EQUIVALENCE ( ASUM( 1 , 4 , 8 ) ,  D13P1P )
      EQUIVALENCE ( ASUM( 1 , 4 , 9 ) ,  D35P1P )
      EQUIVALENCE ( ASUM( 1 , 4 , 10 ) ,  D15P1P )
      EQUIVALENCE ( ASUM( 1 , 4 , 11 ) ,  F35P1P )
      EQUIVALENCE ( ASUM( 1 , 4 , 12 ) ,  F15P1P )
      EQUIVALENCE ( ASUM( 1 , 4 , 13 ) ,  F37P1P )
      EQUIVALENCE ( ASUM( 1 , 4 , 14 ) ,  F17P1P )

      EQUIVALENCE ( ASUM( 1 , 5 , 1 ) ,  S310PP )
      EQUIVALENCE ( ASUM( 1 , 5 , 2 ) ,  S110PP )
      EQUIVALENCE ( ASUM( 1 , 5 , 3 ) ,  P310PP )
      EQUIVALENCE ( ASUM( 1 , 5 , 4 ) ,  P110PP )
      EQUIVALENCE ( ASUM( 1 , 5 , 5 ) ,  P330PP )
      EQUIVALENCE ( ASUM( 1 , 5 , 6 ) ,  P130PP )
      EQUIVALENCE ( ASUM( 1 , 5 , 7 ) ,  D330PP )
      EQUIVALENCE ( ASUM( 1 , 5 , 8 ) ,  D130PP )
      EQUIVALENCE ( ASUM( 1 , 5 , 9 ) ,  D350PP )
      EQUIVALENCE ( ASUM( 1 , 5 , 10 ) ,  D150PP )
      EQUIVALENCE ( ASUM( 1 , 5 , 11 ) ,  F350PP )
      EQUIVALENCE ( ASUM( 1 , 5 , 12 ) ,  F150PP )
      EQUIVALENCE ( ASUM( 1 , 5 , 13 ) ,  F370PP )
      EQUIVALENCE ( ASUM( 1 , 5 , 14 ) ,  F170PP )

      EQUIVALENCE ( ASUM( 1 , 6 , 1 ) ,  S310MP )
      EQUIVALENCE ( ASUM( 1 , 6 , 2 ) ,  S110MP )
      EQUIVALENCE ( ASUM( 1 , 6 , 3 ) ,  P310MP )
      EQUIVALENCE ( ASUM( 1 , 6 , 4 ) ,  P110MP )
      EQUIVALENCE ( ASUM( 1 , 6 , 5 ) ,  P330MP )
      EQUIVALENCE ( ASUM( 1 , 6 , 6 ) ,  P130MP )
      EQUIVALENCE ( ASUM( 1 , 6 , 7 ) ,  D330MP )
      EQUIVALENCE ( ASUM( 1 , 6 , 8 ) ,  D130MP )
      EQUIVALENCE ( ASUM( 1 , 6 , 9 ) ,  D350MP )
      EQUIVALENCE ( ASUM( 1 , 6 , 10 ) ,  D150MP )
      EQUIVALENCE ( ASUM( 1 , 6 , 11 ) ,  F350MP )
      EQUIVALENCE ( ASUM( 1 , 6 , 12 ) ,  F150MP )
      EQUIVALENCE ( ASUM( 1 , 6 , 13 ) ,  F370MP )
      EQUIVALENCE ( ASUM( 1 , 6 , 14 ) ,  F170MP )

      EQUIVALENCE ( ASUM( 2 , 1 , 1 ) ,  S31M3N )
      EQUIVALENCE ( ASUM( 2 , 1 , 2 ) ,  S11M3N )
      EQUIVALENCE ( ASUM( 2 , 1 , 3 ) ,  P31M3N )
      EQUIVALENCE ( ASUM( 2 , 1 , 4 ) ,  P11M3N )
      EQUIVALENCE ( ASUM( 2 , 1 , 5 ) ,  P33M3N )
      EQUIVALENCE ( ASUM( 2 , 1 , 6 ) ,  P13M3N )
      EQUIVALENCE ( ASUM( 2 , 1 , 7 ) ,  D33M3N )
      EQUIVALENCE ( ASUM( 2 , 1 , 8 ) ,  D13M3N )
      EQUIVALENCE ( ASUM( 2 , 1 , 9 ) ,  D35M3N )
      EQUIVALENCE ( ASUM( 2 , 1 , 10 ) ,  D15M3N )
      EQUIVALENCE ( ASUM( 2 , 1 , 11 ) ,  F35M3N )
      EQUIVALENCE ( ASUM( 2 , 1 , 12 ) ,  F15M3N )
      EQUIVALENCE ( ASUM( 2 , 1 , 13 ) ,  F37M3N )
      EQUIVALENCE ( ASUM( 2 , 1 , 14 ) ,  F17M3N )

      EQUIVALENCE ( ASUM( 2 , 2 , 1 ) ,  S31M1N )
      EQUIVALENCE ( ASUM( 2 , 2 , 2 ) ,  S11M1N )
      EQUIVALENCE ( ASUM( 2 , 2 , 3 ) ,  P31M1N )
      EQUIVALENCE ( ASUM( 2 , 2 , 4 ) ,  P11M1N )
      EQUIVALENCE ( ASUM( 2 , 2 , 5 ) ,  P33M1N )
      EQUIVALENCE ( ASUM( 2 , 2 , 6 ) ,  P13M1N )
      EQUIVALENCE ( ASUM( 2 , 2 , 7 ) ,  D33M1N )
      EQUIVALENCE ( ASUM( 2 , 2 , 8 ) ,  D13M1N )
      EQUIVALENCE ( ASUM( 2 , 2 , 9 ) ,  D35M1N )
      EQUIVALENCE ( ASUM( 2 , 2 , 10 ) ,  D15M1N )
      EQUIVALENCE ( ASUM( 2 , 2 , 11 ) ,  F35M1N )
      EQUIVALENCE ( ASUM( 2 , 2 , 12 ) ,  F15M1N )
      EQUIVALENCE ( ASUM( 2 , 2 , 13 ) ,  F37M1N )
      EQUIVALENCE ( ASUM( 2 , 2 , 14 ) ,  F17M1N )

      EQUIVALENCE ( ASUM( 2 , 3 , 1 ) ,  S31P3N )
      EQUIVALENCE ( ASUM( 2 , 3 , 2 ) ,  S11P3N )
      EQUIVALENCE ( ASUM( 2 , 3 , 3 ) ,  P31P3N )
      EQUIVALENCE ( ASUM( 2 , 3 , 4 ) ,  P11P3N )
      EQUIVALENCE ( ASUM( 2 , 3 , 5 ) ,  P33P3N )
      EQUIVALENCE ( ASUM( 2 , 3 , 6 ) ,  P13P3N )
      EQUIVALENCE ( ASUM( 2 , 3 , 7 ) ,  D33P3N )
      EQUIVALENCE ( ASUM( 2 , 3 , 8 ) ,  D13P3N )
      EQUIVALENCE ( ASUM( 2 , 3 , 9 ) ,  D35P3N )
      EQUIVALENCE ( ASUM( 2 , 3 , 10 ) ,  D15P3N )
      EQUIVALENCE ( ASUM( 2 , 3 , 11 ) ,  F35P3N )
      EQUIVALENCE ( ASUM( 2 , 3 , 12 ) ,  F15P3N )
      EQUIVALENCE ( ASUM( 2 , 3 , 13 ) ,  F37P3N )
      EQUIVALENCE ( ASUM( 2 , 3 , 14 ) ,  F17P3N )

      EQUIVALENCE ( ASUM( 2 , 4 , 1 ) ,  S31P1N )
      EQUIVALENCE ( ASUM( 2 , 4 , 2 ) ,  S11P1N )
      EQUIVALENCE ( ASUM( 2 , 4 , 3 ) ,  P31P1N )
      EQUIVALENCE ( ASUM( 2 , 4 , 4 ) ,  P11P1N )
      EQUIVALENCE ( ASUM( 2 , 4 , 5 ) ,  P33P1N )
      EQUIVALENCE ( ASUM( 2 , 4 , 6 ) ,  P13P1N )
      EQUIVALENCE ( ASUM( 2 , 4 , 7 ) ,  D33P1N )
      EQUIVALENCE ( ASUM( 2 , 4 , 8 ) ,  D13P1N )
      EQUIVALENCE ( ASUM( 2 , 4 , 9 ) ,  D35P1N )
      EQUIVALENCE ( ASUM( 2 , 4 , 10 ) ,  D15P1N )
      EQUIVALENCE ( ASUM( 2 , 4 , 11 ) ,  F35P1N )
      EQUIVALENCE ( ASUM( 2 , 4 , 12 ) ,  F15P1N )
      EQUIVALENCE ( ASUM( 2 , 4 , 13 ) ,  F37P1N )
      EQUIVALENCE ( ASUM( 2 , 4 , 14 ) ,  F17P1N )

      EQUIVALENCE ( ASUM( 2 , 5 , 1 ) ,  S310PN )
      EQUIVALENCE ( ASUM( 2 , 5 , 2 ) ,  S110PN )
      EQUIVALENCE ( ASUM( 2 , 5 , 3 ) ,  P310PN )
      EQUIVALENCE ( ASUM( 2 , 5 , 4 ) ,  P110PN )
      EQUIVALENCE ( ASUM( 2 , 5 , 5 ) ,  P330PN )
      EQUIVALENCE ( ASUM( 2 , 5 , 6 ) ,  P130PN )
      EQUIVALENCE ( ASUM( 2 , 5 , 7 ) ,  D330PN )
      EQUIVALENCE ( ASUM( 2 , 5 , 8 ) ,  D130PN )
      EQUIVALENCE ( ASUM( 2 , 5 , 9 ) ,  D350PN )
      EQUIVALENCE ( ASUM( 2 , 5 , 10 ) ,  D150PN )
      EQUIVALENCE ( ASUM( 2 , 5 , 11 ) ,  F350PN )
      EQUIVALENCE ( ASUM( 2 , 5 , 12 ) ,  F150PN )
      EQUIVALENCE ( ASUM( 2 , 5 , 13 ) ,  F370PN )
      EQUIVALENCE ( ASUM( 2 , 5 , 14 ) ,  F170PN )

      EQUIVALENCE ( ASUM( 2 , 6 , 1 ) ,  S310MN )
      EQUIVALENCE ( ASUM( 2 , 6 , 2 ) ,  S110MN )
      EQUIVALENCE ( ASUM( 2 , 6 , 3 ) ,  P310MN )
      EQUIVALENCE ( ASUM( 2 , 6 , 4 ) ,  P110MN )
      EQUIVALENCE ( ASUM( 2 , 6 , 5 ) ,  P330MN )
      EQUIVALENCE ( ASUM( 2 , 6 , 6 ) ,  P130MN )
      EQUIVALENCE ( ASUM( 2 , 6 , 7 ) ,  D330MN )
      EQUIVALENCE ( ASUM( 2 , 6 , 8 ) ,  D130MN )
      EQUIVALENCE ( ASUM( 2 , 6 , 9 ) ,  D350MN )
      EQUIVALENCE ( ASUM( 2 , 6 , 10 ) ,  D150MN )
      EQUIVALENCE ( ASUM( 2 , 6 , 11 ) ,  F350MN )
      EQUIVALENCE ( ASUM( 2 , 6 , 12 ) ,  F150MN )
      EQUIVALENCE ( ASUM( 2 , 6 , 13 ) ,  F370MN )
      EQUIVALENCE ( ASUM( 2 , 6 , 14 ) ,  F170MN )

      REAL RESSGN(31)
      DATA RESSGN/1, -1, -1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, 1, 1,
     $     1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1, 1, 1, 1, 1. /

      REAL PAR(14)
      REAL RETF(6)

      REAL ASL(3),ASR(3),ASS(3)


      DATA IPP/3,5,5,1,1,1,1,1,5,5,5,5,5,3,3,3,7,7,3,3,3,3,3,7,7,3,7,7,
     1     1,1,3/
      DATA NRR/0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
     1     0,0,0/


      Integer OriginalF          ! 0:Lepton mass effects, 1:Original R-S model
      Integer Original
      Integer Lambda
      Real*4 XMLEP, XMLEP2
      Real*4 El, Plep, theta
      Real*4 Q_lab, Q2_lab, Q, nu_lab
      Real*4 jl(4), jl_n1, jl_n2, Norm, C_L, C_R, C_S, Qs, nus, ng
      Real*4 Ei, Ef, nu
      Real*4 B_BRS,C_BRS,Ctmp

C-------------------------------------------------------------------------

      XMASPI2= XMASPI * XMASPI
      XMVSPI2= XMVSPI * XMVSPI

      Original = OriginalF

      IF (IMODE.GT.10) THEN
         IFLAG = IMODE - 10
      ELSE
         IFLAG = IMODE
      ENDIF

      IF (IFLAG.GT.7) GOTO 9999

      WMIN = XMN+XMPI
      SMAX = (2 * XMN * E) + XMN2
      QMAX = (SMAX - XMN2)

      W2 = W**2

      IF ( ((2 * XMN * E) * (2 * XMN * E + XMN2 - W2) /
     $     (2 * XMN * E + XMN2)).LT.(-Q2) ) THEN
C        WRITE(*,*) "RSDCRS: Q2 Kinematically disallowed: ", ((2 * XMN * E) *
C     $  (2 * XMN * E + XMN2 - W2) / (2 * XMN * E + XMN2)), ".LT.", -Q2
        GOTO 999
      ENDIF

C     -------------------------------------------------------
C     NOW  W  AND  E  ARE FIXED
C     Q2MAX = 2.*XMN*E*(2.*XMN*E+XMN2-W2)/(2.*XMN*E+XMN2)
C     -------------------------------------------------------

C     KINEMATICAL QUANTITIES

      XNU = (W2 - XMN2 - Q2) / (2 * XMN)
      X = -Q2 / (2 * XMN * XNU)

      IF(XNU.GE.E) THEN
C        WRITE(*,*) "RSDCRS: Energy Kinematically disallowed: ", XNU, ".GE.", E
        GOTO 999
      ENDIF

      XNU2 = XNU**2
      ES = E-XNU
      QV = SQRT(XNU2-Q2)
      QV2 = QV**2
      QVS = QV*XMN/W
      QVS2 = QVS**2
      XNUS = ABS(W2-XMN2+Q2)/(2.*W)
      XLAM = QVS*SQRT(2./OMEG)
      XLAM2 = XLAM**2
      A = ((W+XMN)**2-Q2)/(2.*W)

      U = (E+ES+QV)/(2.*E)
      V = (E+ES-QV)/(2.*E)
      IF(V.LT.0.) V = 0.
      U2 = U**2
      V2 = V**2
      AK = (W2-XMN2)/(2.*XMN)
      FAKT = W*PI/(W2-XMN2)

C     ---------------------------
C     --- Lepton mass effects ---
      XMLEP2 = XMLEP**2

C     Lepton energy
      El = E+XMN-(XMN2+W**2-q2)/(2*XMN)
C
      Ei = (2.*XMN*E+Q2-XMLEP2)/(2.*W)
      Ef = (2.*XMN*El-Q2-XMLEP2)/(2.*W)
      nu = Ei - Ef
      Q = QVS

C     Lab frame
      if (El**2-XMLEP2.le.0.) then
C         Plep = 0.000001
         Plep = 1e-7
      else
         Plep = SQRT(El**2-XMLEP2)
      endif

      theta = (Q2-XMLEP2+2*E*El)/(2*E*Plep)
      if (theta.ge.0.999999 ) then
         theta = 0.999999
         goto 999
      endif
      if (theta.le.-0.999999) then
         theta =-0.999999
         Original = 1
      endif

      nu_lab = E - El
      Q_lab = SQRT(nu_lab**2-Q2)
      Q2_lab = Q_lab**2

      C_L = u
      C_R = v
      C_S = sqrt(2.*u*v)

      Qs = 1.
      nus = 1.

      if (Original.ne.1) then
C     Leptonic current in RRF
C     jl_0, jl_x, i*jl_y, jl_z
         jl(1) = (XMN-El-Lambda*Plep)/W
         jl(2) = (Plep-Lambda*E)/Q_lab
         jl(3) = -Lambda
         jl(4) = ((E+Lambda*Plep)*(XMN-El)+Plep
     $          *(Lambda*E+2.*E*theta-Plep))/(W*Q_lab)

         jl_n1 = SQRT(1.-Lambda*theta)*SQRT(E*(El-Lambda*Plep))
         jl_n2 = SQRT(1.+Lambda*theta)*SQRT(E*(El-Lambda*Plep))

         jl(1) = jl_n1*jl(1)
         jl(2) = jl_n2*jl(2)
         jl(3) = jl_n2*jl(3)
         jl(4) = jl_n1*jl(4)

C     Coefficients C_L, C_R, C_S and Qs, nus
         Norm = SQRT(2.)*E/Q_lab*SQRT(-Q2)
         C_L = (jl(2)+jl(3))/(SQRT(2.)*Norm)
         C_R = (jl(3)-jl(2))/(SQRT(2.)*Norm)
         C_S = SQRT(ABS(jl(1)**2-jl(4)**2))/Norm

         if (C_S.eq.0) then
            Original = 1
            goto 200
         endif

         Qs  = jl(1)/(Norm)*SQRT(-Q2)/C_S
         nus = jl(4)/(Norm)*SQRT(-Q2)/C_S
      endif


C     SUMMATION OVER RESONANCES

 200  FRP1 = 0.
      FRP3 = 0.
      FRM1 = 0.
      FRM3 = 0.
      FR0P = 0.
      FR0M = 0.
      FIP1 = 0.
      FIP3 = 0.
      FIM1 = 0.
      FIM3 = 0.
      FI0P = 0.
      FI0M = 0.

      DO 60  IB = 1,31

         IBLOCK = IB-1
         NR = NRR(IB)
         XMR = XMRR(IB)
         BR = BRR(IB)
C        SETS UP NORMALIZATION CONSTANTS
C        G&S don't use the BW normalizing factors for their model
C        NOTE: This should only work for G&S form factors only

         if(IFFSPI.ne.0.and.IB.eq.1) then
            FBWNO = 1
         else
            FBWNO = FBWNOO(IB)
         endif

C        FBWNOO is set in rscons.h

C         ETA = ETAA(IB)
         XE = XEE(IB)
         IP = IPP(IB)

C     **********************************

C     FORMFACTORS AND RELATED QUANTITIES

         if (IFFSPI.ne.0) then

C       K.M. Graczyk and J.T. Sobczyk, Phys.Rev.D 77,053001 (2008)
C       Equations are correction in erratum
C       CA5(0) and XMA are set by user

            CA5=   RCA5ISPI/((1-Q2/XMASPI2)**2)

C     J.A. Nowak, arXiv:0909.3659
C     Vector Form Factors
C     CV1, CV2 are for simplicity, not true form factors
C     Added the proportional factors of the vector form factors GV
C     to CV3,CV4,CV5

C     (PAR) Yargh, the variable Q2 here has a minus sign relative to
C     Q2 in the Nowak/Graczyk & Sobczyk papers

            CV1=  1./((1.-Q2/(4.*XMVSPI2))
     $           * ((1.-Q2/XMVSPI2)**2))
            CV2=  1./((1.-Q2/(4.*XMVSPI2))
     $           * ((1.-Q2/(0.766*XMVSPI2))**2))
            CV3A= 2.13*CV1*((W+XMN)/XMN)
            CV3B= 2.13*CV1*((XMN2-Q2+XMN*W)/(XMN*W))
            CV4= -1.51*CV1*((W2+Q2-XMN2)/(2.*XMN2))
            CV5=  0.48*CV2*((W2-Q2-XMN2)/(2.*XMN2))

            GV3=  (.5/SQRT(3.))*(CV4+CV5+CV3A)
            GV1=  (-.5/SQRT(3.))*(CV4+CV5-CV3B)

            NFF=  SQRT((3.*(GV3)**2)+GV1**2)

C     Note, J.A. Nowak has different form factors that those quoted
C     in G&S (eq. 58). G&S seems to match previous MC better

            FFV = 0.5 * (1 - Q2/(XMN+W)**2)**0.5
     &           * (1 - Q2/(4*XMN2))**(-1*NR)
     &           * NFF

C     Axial. Divide by Z, the arbitrary Rein-Sehgal factor, because
C     we multiply by it later

            FFA = 0.5 * SQRT(3.0) * (1 - Q2/(XMN+W)**2)**0.5
     &           * (1 - Q2/(4*XMN2))**(-1*NR)
     &           * (1 - (W2 + Q2 - XMN2)/(8*XMN2))
     &           * CA5 / Z

C     Rein-Seghal form factors
C     Difference between Rein-Seghal form factors are shown in
C     Phys Rev D 77 (053001), equation 12, 13
         else
            if(NRTYPESPI.eq.0) then
               FFV = ((1.-Q2/(4.*XMN2))**(.5-NR/2.))
     $              /(1.-Q2/XMVSPI2)**2
               FFA = ((1.-Q2/(4.*XMN2))**(.5-NR/2.))
     $              /(1.-Q2/XMASPI2)**2
            else
               FFV = ((1.-Q2/(4.*XMN2))**(.5-NR))
     $              /(1.-Q2/XMVSPI2)**2
               FFA = ((1.-Q2/(4.*XMN2))**(.5-NR))
     $              /(1.-Q2/XMASPI2)**2
            endif
         endif

C     ATTENTION: NR CHANGED INTO NR/2.

         RV = SQRT(2.)*FFV*QVS*(W+XMN)/(2.*W*A)
         RA = SQRT(2.)*FFA*Z*(W+XMN+NR*OMEG/A)/(6.*W)
         R  = RV
         RP = -(RV+RA)
         RM = -(RV-RA)
         TV = FFV*SQRT(OMEG/2.)/(3.*W)
         TA = FFA*SQRT(OMEG/2.)*Z*QVS/(3.*W*A)
         T  = TV
         TP = -(TV+TA)
         TM = -(TV-TA)

C     ---------------------------
C     --- Lepton mass effects ---
         ng = (W+XMN)**2-Q2

C     New definitions of S, B, C
         S = (Qs*Q-nus*nu)*(3.*W*XMN+Q2-XMN2)/(6.*Q2_lab*XMN2)*FFV
         B = 2.*SQRT(OMEG/2.)*(Qs+2.*W*nus*Q/ng)/(6.*W*Q)*Z*FFA
         C = ((Qs*Q-nus*nu)*(1./3.+2.*W*nu/ng)
     &        +nus*(2.*W/3.+(NR*OMEG/3.+Q2)*2.*W/ng))/(2.*W*Q)*Z*FFA

C     Pion-pole contribution (according Berger and Sehgal ArXiV:0709.4387)
C     Correction is added to B_BRS(2M_Na should be M_Na in the reference)
         B_BRS = (Qs*nu-nus*Q)/(3.*W*Q*(XMPI**2-Q2))
     &        *sqrt(OMEG/2.)*(nu+2.*Q**2*W/ng)*Z*FFA
         C_BRS = (Qs*nu-nus*Q)/(XMPI**2-Q2)
     &        *(1./3.+Q2/ng+NR*OMEG/(3.*ng))*Z*FFA

         B = B + B_BRS
         C = C + C_BRS
C     ---------------------------
C     ---------------------------

*     Original dynamical form factors
         if (Original.eq.1) then
            S = (FFV*(3.*W*XMN-XMN2+Q2)/(6.*W2))*(-Q2/QVS2)
            B = FFA*Z*SQRT(OMEG/2.)*(1.+XNUS/A)/(3.*W)
            C = FFA*Z*(W2-XMN2+NR*OMEG*XNUS/A)/(6.*W*QVS)
         endif

         PAR( 1) = RV
         PAR( 2) = RA
         PAR( 3) = R
         PAR( 4) = RP
         PAR( 5) = RM
         PAR( 6) = TV
         PAR( 7) = TA
         PAR( 8) = T
         PAR( 9) = TP
         PAR(10) = TM
         PAR(11) = S
         PAR(12) = B
         PAR(13) = C
         PAR(14) = XLAM

         CALL RSCLFM(IFLAG,IBLOCK,PAR,RETF)

         IF (IFLAG.LE.3) THEN
            FM3C = RETF(1)
            FM1C = RETF(2)
            FP3C = RETF(3)
            FP1C = RETF(4)
            F0PC = RETF(5)
            F0MC = RETF(6)
         ELSE
            FM3 = RETF(1)
            FM1 = RETF(2)
            FP3 = RETF(3)
            FP1 = RETF(4)
            F0P = RETF(5)
            F0M = RETF(6)
         ENDIF


C     **********************************

C     STORE RESONANCE AMPLITUDES

         XXE = SQRT(XE)

         RAD =((W2-XMN2-XMPI2)**2-4.*XMN2*XMPI2)/(4.*W2)
         IF(RAD.LE.0.) RAD = 0.0001
         PQ = SQRT(RAD)
         RADR = ((XMR**2-XMN2-XMPI2)**2-4.*XMN2*XMPI2)/(4.*XMR**2)
         IF(RADR.LE.0.) RADR = 0.0001
         PQR = SQRT(RADR)
         BRA = BR*(PQ/PQR)**IP
         FBW = (BRA/(2.*PI))/((W-XMR)**2+BRA**2/4.)

         ETABW(IB) = CMPLX(W-XMR,-1.*BRA/2.)
         ETABW(IB) = ETABW(IB) / ((W-XMR)**2+BRA**2/4.)
         ETABW(IB) = ETABW(IB)*SQRT(BRA / (2.*PI) / FBWNO)

         FFBW(IB) = FBW/FBWNO

         ABWR(IB) = FFBW(IB)*(W-XMR)/SQRT(BRA/(2.*PI))
         ABWR(IB) = ABWR(IB)*SQRT(FBWNO)
         ABWI(IB) = FFBW(IB)*SQRT(BRA*PI/2.)
         ABWI(IB) = ABWI(IB)*SQRT(FBWNO)
         IF(IFLAG.LE.3)  GOTO 991
         AM3(IB) = FM3*XXE
         AM1(IB) = FM1*XXE
         AP1(IB) = FP1*XXE
         AP3(IB) = FP3*XXE
         A0P(IB) = F0P*XXE
         A0M(IB) = F0M*XXE
         AM1C(IB) = AM1(IB)
         AM3C(IB) = AM3(IB)
         AP1C(IB) = AP1(IB)
         AP3C(IB) = AP3(IB)
         A0PC(IB) = A0P(IB)
         A0MC(IB) = A0M(IB)
         GOTO 992
 991     CONTINUE
         AM3C(IB) = FM3C*XXE
         AM1C(IB) = FM1C*XXE
         AP1C(IB) = FP1C*XXE
         AP3C(IB) = FP3C*XXE
         A0PC(IB) = F0PC*XXE
         A0MC(IB) = F0MC*XXE
 992     CONTINUE

C     WE CHOOSE I=1/2, J=1/2 NONRES. INCOHERENT BACKGROUND

         CO = ((W2-WMIN**2)/(4.*W2))**.25
         IF(IB.GE.29)  FFBW(IB) = 1.*(CO**2)
         IF(IB.eq.29)  ETABW(IB) = RBGSCLSPI*CMPLX(SQRT(FFBW(IB)))
         IF(IB.eq.30)  ETABW(IB) = RBGSCLSPI*CMPLX(SQRT(FFBW(IB)))

C     **********                 ****CHANGE WITH COH/INCOH BACKGROUND

 60   CONTINUE

      DO 61 NCTR=1,31
         AM(1,NCTR) = AM3C(NCTR)
         AM(2,NCTR) = AM1C(NCTR)
         AM(3,NCTR) = AP3C(NCTR)
         AM(4,NCTR) = AP1C(NCTR)
         AM(5,NCTR) = A0PC(NCTR)
         AM(6,NCTR) = A0MC(NCTR)
 61   CONTINUE


      CALL RSASUM(AM,ETABW,RESSGN,ASUM)


      CO = ((W2-(XMN+XMPI)**2)/(4.*W2))**.25

C     **********                 ****CHANGE WITH COH/INCOH BACKGROUND

C     SELECT APPROPRIATE NUCLEON - PION CHANNEL

      IF(IFLAG.EQ.1)  GOTO 2010
      IF(IFLAG.EQ.2)  GOTO 2020
      IF(IFLAG.EQ.3)  GOTO 2030
      IF(IFLAG.EQ.4)  GOTO 2040
      IF(IFLAG.EQ.5)  GOTO 2050
      IF(IFLAG.EQ.6)  GOTO 2060
      IF(IFLAG.EQ.7)  GOTO 2070


C     START SELECTIVE SUMMATION IN FIRST CHANNEL

 2010 CONTINUE

      CALL RSIFL1(ASUM,ETABW,RESSGN,AM,ASL,ASR,ASS)

      GOTO 69

 2020 CONTINUE

      CALL RSIFL2(ASUM,ETABW,RESSGN,AM,ASL,ASR,ASS)

      GOTO 69

 2030 CONTINUE

      CALL RSIFL3(ASUM,ETABW,RESSGN,AM,ASL,ASR,ASS)

      GOTO 69

 2040 CONTINUE

      CALL RSIFL4(ASUM,ETABW,RESSGN,AM,ASL,ASR,ASS)

      GOTO 69

 2050 CONTINUE

      CALL RSIFL5(ASUM,ETABW,RESSGN,AM,ASL,ASR,ASS)

      GOTO 69

 2060 CONTINUE

      CALL RSIFL6(ASUM,ETABW,RESSGN,AM,ASL,ASR,ASS)

      GOTO 69

 2070 CONTINUE

      CALL RSIFL7(ASUM,ETABW,RESSGN,AM,ASL,ASR,ASS)

      GOTO 69

 69   CONTINUE

      SL = ASL(1)
      SR = ASR(1)
      SS = ASS(1) * (-QVS2 / Q2)

      SL1 = ASL(2)
      SR1 = ASR(2)
      SS1 = ASS(2) * (-QVS2 / Q2)

      SL0 = ASL(3)
      SR0 = ASR(3)
      SS0 = ASS(3) * (-QVS2/ Q2)



      IF (IMODE.GT.10) THEN     ! for anti-neutrino
         Ctmp = C_L
         C_L = -C_R
         C_R = -Ctmp
      ENDIF

      DSI  = AK*(W/XMN)*(-Q2/QV2)*FAKT*((C_L**2)*SL  + (C_R**2)*SR
     $       + (C_S**2)*SS)
      DSI1 = AK*(W/XMN)*(-Q2/QV2)*FAKT*((C_L**2)*SL1 + (C_R**2)*SR1
     $       + (C_S**2)*SS1)
      DSI0 = AK*(W/XMN)*(-Q2/QV2)*FAKT*((C_L**2)*SL0 + (C_R**2)*SR0
     $       + (C_S**2)*SS0)

C     dsigma/(dQ2*dW)
      DSI  = DSI * 1.321
      DSI1 = DSI1 * 1.321
      DSI0 = DSI0 * 1.321

      if (Original.eq.1) then
         DSI  = DSI/2.
         DSI1 = DSI1/2.
         DSI0 = DSI0/2.
      endif

      GOTO 100

 999  CONTINUE
c      write(*,*) "Wrong arguments."
      DSI = -0.0000000E1
      DSI0 = -0.0000000E1
      GOTO 100

 9999 DSI = -9999.999
      DSI0 = -9999.999

 100  DIFCRS = DSI
      DFCRS3 = DSI0

      END
