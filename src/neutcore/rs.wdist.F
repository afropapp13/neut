***********************************************************************
*     ---------------------------------------------------
*      SUBROUTINE RSCTCS
*     ---------------------------------------------------
*
*     ( purpose )
*       calculate total cross-section of single-pion production
*     ( input )
*       none
*
*     ( output )
*       none
*       (creates file named "rstcrs.tbl".)
*
*     ( creation date and author )
*       1993.Dec. ; for Super-Kamioka by Y.Hayato
*                   (Based on Rein and Sehgal's code)
*
*     ( comment )
*       Reference:Rein and Sehgal,Ann. of Phys. 133(1981),79-153
*                 Rein           ,Z.Phys.C 35(1987),43-64
*                 Feynman et al. ,Phys.Rev.D3(1971),2706-
*
*       ----------- IFLAG means Interaction mode ----------
*         
*       IFLAG < 10 -> neutrino
*       IFLAG > 10 -> anti neutrino
*
*       IFLAG = 1:'CHARGED CURRENT, P - PI+'
*       IFLAG = 2:'CHARGED CURRENT, P - PI0'
*       IFLAG = 3:'CHARGED CURRENT, N - PI+'
*       IFLAG = 4:'NEUTRAL CURRENT, P - PI0'
*       IFLAG = 5:'NEUTRAL CURRENT, N - PI+'
*       IFLAG = 6:'NEUTRAL CURRENT, N - PI0'
*       IFLAG = 7:'NEUTRAL CURRENT, P - PI-'
*     
*       IFLAG =11 ... -> charge conjugate
*
***********************************************************************

      IMPLICIT NONE

#include<rscons.h>

      INTEGER IE,IFLAG,NQ,IQ,IW,NW,IMODE
      REAL WMIN,DW,DQ2,SIGMA,SIGMA3,QMAX,SMAX
      REAL SIQ,SIQ3,Q2,E,W
      REAL DSI,DSI3

C
      
      INTEGER IEMAX / 20/
      REAL ETBL(20)
      
      DATA ETBL/0.37,0.4,0.5,0.625,0.75,0.875,1.00,1.125,1.25,
     $     1.375,1.50,2.00,2.50,3.00,3.50,4.75,6.00,10.0,50.0,
     $     100.0/

      CHARACTER*30 SMODE(14)

      DATA SMODE/'NU P -> P PI+ L-','NU N -> P PI0 L-',
     $           'NU N -> N PI+ L-','NU P -> P PI0 NU',
     $           'NU P -> N PI+ NU','NU N -> N PI0 NU',
     $           'NU N -> P PI- NU',
     $           'NU-BAR N -> N PI- L+','NU-BAR P -> N PI0 L+',
     $           'NU-BAR P -> P PI- L+','NU-BAR P -> P PI0 NU-BAR',
     $           'NU-BAR P -> N PI+ NU-BAR','NU-BAR N -> N PI0 NU-BAR',
     $           'NU-BAR N -> P PI- NU-BAR'/
C
C
C Program begins from here.      
C
C

      WRITE(*,*) "#"
      WRITE(*,*) "#Neutrino Cross-section table (Rein Sehgal)"
      WRITE(*,*) "#"

C     set COMMON /RSCONS/

      WMIN = XMN+XMPI
      
      
      DW = .02
      DQ2 = 0.05
            
C     read constants
      
      WRITE(*,*) "# Parameters"
      WRITE(*,*) "#    W  :",WMAX
      WRITE(*,*) "#    DW :",DW
      WRITE(*,*) "#    DQ2:",DQ2
      WRITE(*,*) "#"
      
      DO 105  IMODE = 2,2

         IFLAG = IMODE
         IF (IMODE.GT.7) IFLAG=IMODE-7+10

         WRITE(*,*) "#"
         WRITE(*,*) "#MODE:",SMODE(IMODE)
         WRITE(*,*) "Now Calculating! MODE:",SMODE(IMODE)
         WRITE(*,*) IEMAX

         DO 95 IE = 1,1
            
            E = 4.
            
C     NOW  E IS FIXED
            
C     X-DEPENDENCE
            
            SIGMA = 0.
            SIGMA3 = 0.
            
            SMAX = (2 * XMN * E) + XMN2
            
            QMAX = (SMAX - XMN2)
            
            NQ = INT(QMAX / DQ2)
            
            NW = INT((WMAX - WMIN) / DW)
            
            DO 90  IW = 1,NW
                  
               W = WMIN+DW*(FLOAT(IW)-.5)
               
               SIQ = 0.
               SIQ3 = 0.

               DO 70 IQ = 1,NQ
                  Q2 = -DQ2 * (FLOAT(IQ) - 0.5)
                  CALL RSDCRS(IFLAG,E,Q2,W,DSI,DSI3)
                  IF (DSI.LT.0.0) GOTO 999
                  
                  SIQ = SIQ+DQ2*DSI
                  SIQ3 = SIQ3+DQ2*DSI3

 70            CONTINUE
 71            CONTINUE

               WRITE(20,*) E,W,SIQ,SIQ3
               WRITE(*,*) E,W,SIQ,SIQ3
               
               SIGMA = SIGMA + DW * SIQ
               SIGMA3 = SIGMA3 + DW * SIQ3
               
 90         CONTINUE
            
            write(*,*) "SIGMA=",sigma,"  : SIGMA3=",sigma3      
 95      CONTINUE
 105  CONTINUE
      GOTO 1000
      
 999  WRITE(*,*) "Error occured in rsdcrs."
      WRITE(*,*) "iflag=",IFLAG
      WRITE(*,*) "e    =",E
      WRITE(*,*) "q2   =",Q2
      WRITE(*,*) "w    =",W
      WRITE(*,*) "dsi  =",DSI
      
 1000 CONTINUE
      GOTO 5000
      
 9999 WRITE(*,*) "Can not write output file!"

 5000 CONTINUE
      END
      
