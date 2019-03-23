***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE RSTCRSG(IMODE,E,SIGMA,SIGMA3)
*     ---------------------------------------------------
*
*     ( purpose )
*       return cross-section for specified mode and energy.
*              for electron neutrino
*
*     ( input )
*       IMODE       : Interaction mode
*       E           : Energy of incident neutrino
*
*     ( output )
*       sigma       : cross-section
*       sigma3      : cross-section(through P33(1232)-resonance)
*
*     ( creation date and author )
*       1993.Dec. ; for Super-Kamioka by Y.Hayato
*                   (Based on Rein and Sehgal's code)
*
*       1998.Mar. : change W<1.4 -> W<2.0  by J.K.
*       2007.08.23 ; G.Mitsuka add delta->gamma decay
*       2007.12.05 ; G.Mitsuka add cross section of MA=1.21
*
*     ( comment )
*       Reference:Rein and Sehgal,Ann. of Phys. 133(1981),79-153
*                 Rein           ,Z.Phys.C 35(1987),43-64
*                 Feynman et al. ,Phys.Rev.D3(1971),2706-
***********************************************************************

      IMPLICIT NONE

#include "rscons.h"
#include "neutfilepath.h"
#include "neutmodel.h"

      INTEGER IMODE
      REAL E,SIGMA,SIGMA3

      INTEGER SFLAG
      SAVE SFLAG
      DATA SFLAG/0/

      INTEGER IEMAX
      PARAMETER (IEMAX=24)

      INTEGER LUNI
      PARAMETER (LUNI=88)

      INTEGER IP
      PARAMETER (IP=12)

      INTEGER I

      INTEGER ISPIMA

      SAVE XNUCRS
      
      REAL ETBL(IEMAX),XNUCRS(17,2,IEMAX)

      DATA ETBL/0.37,0.4,0.5,0.625,0.75,0.875,1.00,1.125,1.25,
     $     1.375,1.50,2.00,2.50,3.00,3.50,4.75,6.00,10.0,50.0,
     $     100.0,200.0,400.0,650.0,1000.0/

      CHARACTER*256 DUMSTR
      character*1024 fname
      integer*4      LENSTR
      external       lenstr
      integer*4      lenpath

      lenpath = lenstr(crstblpath)
      ISPIMA = INT(XMARES*100+0.5)
      
      IF (SFLAG.eq.0) then
         IF (ISPIMA.eq.121) then
            fname = 'sgam_nue_xsec_ma1.21.dat'
            fname = crstblpath(1:lenpath)//fname            
            open(LUNI,file=fname,
     $           form='formatted',status='OLD')
         else if (ISPIMA.eq.111) then
            fname = 'sgam_nue_xsec_ma1.11.dat'
            fname = crstblpath(1:lenpath)//fname            
            open(LUNI,file=fname,
     $           form='formatted',status='OLD')
         else
            write(*,*) 'Invalid MA value was specified: ',ISPIMA
            stop
         endif
         
         read(LUNI,*) DUMSTR

         call rstcrsgamrd(LUNI,12,IEMAX,ETBL,XNUCRS)
         
         close(LUNI)

         SFLAG=1
      endif
            
      IF (E.GE.ETBL(IEMAX)) GOTO 9999

      DO 110 I=1,IEMAX
         IF (E.LT.ETBL(I)) THEN
            IF (I.NE.1) THEN
               SIGMA = (E-ETBL(I-1))/(ETBL(I)-ETBL(I-1))*
     $              (XNUCRS(IMODE,1,I)-XNUCRS(IMODE,1,I-1))
     $              +XNUCRS(IMODE,1,I-1)
               SIGMA3 = (E-ETBL(I-1))/(ETBL(I)-ETBL(I-1))*
     $              (XNUCRS(IMODE,2,I)-XNUCRS(IMODE,2,I-1))
     $              +XNUCRS(IMODE,2,I-1)
               GOTO 120
            ELSE
               SIGMA = 0.000000E0
               SIGMA3= 0.000000E0
               GOTO 120
            ENDIF
         ENDIF
 110  CONTINUE

 120  CONTINUE
      SIGMA = SIGMA / 10.
      SIGMA3= SIGMA3/ 10.

      RETURN

 9999 SIGMA=XNUCRS(IMODE,1,IEMAX)/10.
      SIGMA3=XNUCRS(IMODE,2,IEMAX)/10.
      RETURN

      END
