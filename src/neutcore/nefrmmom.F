************************************************************************
*     ------------------------
      SUBROUTINE NEFRMMOM(PFMEV,IWAVE)
*     ------------------------
*
*     (Purpose)
*       return fermi momentum according to the distribution
*       given by RNSWAV,RNPWAV (same as FVFRMGEV)
*
*     (Input)
*       none
*     
*     (Output)
*       PFMEV : Fermi momentum in MeV/c <- different from FVFRMGEV
*       IWAVE : S-WAVE : 0
*               P-WAVE : 1      
*
*     (Creation Date and Author)
*       2000.10.13 ; Y.Hayato
*     
************************************************************************
      IMPLICIT NONE
      REAL*4     PFMEV
      INTEGER*4  IWAVE

      REAL*4     RLU
      EXTERNAL   RLU

      INTEGER*4     DUMMY

#include "neutparams.h"

      real*4     XHIT
      INTEGER*4  IINIT
      SAVE       IINIT
      DATA       IINIT/0/

C--- constant case
      IF (IINIT.EQ.0) THEN
         write(*,*) 'Use Flat distribution for PF'
         IINIT = 1
      ENDIF
 100  continue
        xhit=pfmax**2*rlu(dummy)*10.
        PFMEV = (PFMAX*1000.)*RLU(DUMMY)
        if (xhit.gt.(pfmev/1000.)**2) goto 100
 200  continue

/* IWAVE does not have a meaning */      
      IWAVE = 1

      RETURN
      END

