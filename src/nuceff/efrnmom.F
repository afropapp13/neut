************************************************************************
*     ------------------------
      SUBROUTINE EFRNMOM(RADIUS,PFMEV)
*     ------------------------
*
*     (Purpose)
*       return fermi momentum ( local Fermi Gas in MeV )
*
*     (Input)
*       RADIUS: in fm
*
*     (Output)
*       PFMEV : Fermi momentum in MeV/c <- different from FVFRMGEV
*
*     (Creation Date and Author)
*       2017.05.27 ; Y.Hayato
*
************************************************************************
      IMPLICIT NONE
      REAL*4     RADIUS,PFMEV

      REAL*4     RLU,EFFRMGAS
      external   RLU,EFFRMGAS

      INTEGER*4  IDUM

      real*4     XHIT, PFMAX,RDUM1,RDUM2

      PFMAX = EFFRMGAS(RDUM1,RDUM2,RADIUS)

 100  continue
        xhit  = pfmax**2*rlu(IDUM)*10.
        PFMEV = PFMAX*RLU(IDUM)
        if (xhit.gt.(pfmev**2)) goto 100
 200  continue

      RETURN
      END

