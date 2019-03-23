***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE RSTCRST(IMODE,E,SIGMA,SIGMA3)
*     ---------------------------------------------------
*
*     ( purpose )
*       return cross-section for specified mode and energy.
*              for tau neutrino
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
*       2007.Dec. ; add lepton mass effects by G.Mitsuka
*       2009.Feb. ; change Z=0.74 -> Z=0.76 by Y.Nakajima
*       2010.Oct. ; use file for the cross-section table by Y.Hayato
*
*     ( comment )
*       Reference:Rein and Sehgal,Ann. of Phys. 133(1981),79-153
*                 Rein           ,Z.Phys.C 35(1987),43-64
*                 Feynman et al. ,Phys.Rev.D3(1971),2706-
*                 K.S.Kuzmin et al., hep-ph/0312107
*                 Berger and Sehgal, hep-ph/0709.4378
***********************************************************************

      IMPLICIT NONE

#include "neutmodel.h"
#include "rscons.h"
#include "neutfilepath.h"

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
      PARAMETER (IP=16)
      
      INTEGER I

      INTEGER ISPIMA

      SAVE XNUCRS

      REAL ETBL(IEMAX),XNUCRS(17,2,IEMAX)

      character*15 fname_ref(3)
      character*1024 fname

      data fname_ref/
     $     'spi_nue_xsec_ma',
     $     'spi_num_xsec_ma',
     $     'spi_nut_xsec_ma'/

      DATA ETBL/0.37,0.4,0.5,0.625,0.75,0.875,1.00,1.125,1.25,
     $     1.375,1.50,2.00,2.50,3.00,3.50,4.75,6.00,10.0,50.0,
     $     100.0,200.0,400.0,650.0,1000.0/

      CHARACTER*256 DUMSTR
      integer*4      LENSTR
      external       lenstr
      integer*4      lenpath

      lenpath = lenstr(crstblpath)
      ISPIMA = INT(XMASPI*100+0.5)

      IF (SFLAG.eq.0) then

         if (IFFSPI.ne.0) then

C- new form factor
            write(fname,'(A,F4.2,A,F4.2,A,F4.2,A)') fname_ref(3), 
     $           xmaspi, '_ca5I', RCA5ISPI,
     $           '_bgscl', RBGSCLSPI,'.dat'

         else

C- original form factor
            write(fname,'(A,F4.2,A)') fname_ref(3), xmaspi, '.dat'

         endif

         fname = crstblpath(1:lenpath)//fname         
         write(6,*) 'Opening file ', fname(1:lenstr(fname))

         open(LUNI,file=fname,
     $        form='formatted',status='OLD')

        read(LUNI,*) DUMSTR

        call rstcrsrd(LUNI,IP,IEMAX,ETBL,XNUCRS)

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

      IF((E.LT.3.92).and.((IMODE.lt.4).or.((IMODE.gt.10).and.(IMODE.
     $     lt.14)))) THEN
*        write(*,*) "RSTCTST: E<3.92GeV",E,SIGMA,SIGMA3
         SIGMA=0.
         SIGMA3=0.
      ENDIF

      RETURN

 9999 SIGMA=XNUCRS(IMODE,1,IEMAX)/10.
      SIGMA3=XNUCRS(IMODE,2,IEMAX)/10.
      RETURN

      END
