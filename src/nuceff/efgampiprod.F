*     -----------------------------------------------------
      SUBROUTINE EFGAMPIPROD(EGAM,IDNUCI,IDPION,IDNUCF,COSPI,
     &                       IMODE,IERR)
*     -----------------------------------------------------
*
*     ( purpose )
*       Determine pion/eta type, direction and
C       weight after gam+N interaction
*       And outgoing nucleon type
*
*     ( input )
C         EGAM   : PHOTON ENERGY (MeV) (in nucleon rest frame)
C         IDNUCI : NUCLEON PID (PDG)
*
*     ( output )
c         IDPION : outgoing pion/eta type (PDG)
C         IDNUCF : outgoing nucleon type
C         COSPI  : COS(C.M. Angle of outgoing pion)
C         WEIGHT : Event weight = gam+N total cross section
C                  (must adjust for multi-nucleon nuclei in
C                   calling program)
C         IMODE  : Reaction mode label
*     (Creation Date and Author)
*         2010.01.18 ; P.de Perio
*
      IMPLICIT NONE

      REAL*4    RLU
      external  RLU

      REAL EFGAMPIMP
      EXTERNAL EFGAMPIMP
      REAL EFGAMPIMPT
      EXTERNAL EFGAMPIMPT

      REAL EFGAMPIPN
      EXTERNAL EFGAMPIPN
      REAL EFGAMPIPNT
      EXTERNAL EFGAMPIPNT

      REAL EFGAMPI0P
      EXTERNAL EFGAMPI0P
      REAL EFGAMPI0PT
      EXTERNAL EFGAMPI0PT

      REAL EFGAMPI0N
      EXTERNAL EFGAMPI0N
      REAL EFGAMPI0NT
      EXTERNAL EFGAMPI0NT

      integer IDNUCI,IDPION,IDNUCF,IDUM,I,IMODE,IERR
      real EGAM,COSPI,WEIGHT,dum

      real piDcrs,pi0Dcrs,etaDcrs
      real sumDcrs,maxDcrs

      integer nSteps
      real stepInt

      IERR = 0

C     Find maximum of the diff. xsec. at this energy.
C     For best results, make sure this section
C     corresponds to the steps in input data.

      cospi = -1.
      maxDcrs = 0.
      nSteps = 41
      stepInt = 2./(nSteps-1.)

      DO 100 I=1,nSteps

         if (cospi.GT.1.)  cospi=1.

         if (IDNUCI.EQ.2212) then
            sumDcrs = efgampipn(EGAM,cospi)+efgampi0p(EGAM,cospi)
         else if (IDNUCI.EQ.2112) then
            sumDcrs = efgampimp(EGAM,cospi)+efgampi0n(EGAM,cospi)
         else
            write (*,*) "Error gampi: Must be a nucleon"
            return
         end if

         if (sumDcrs.GE.maxDcrs) maxDcrs = sumDcrs

         cospi = cospi + stepInt

 100  CONTINUE




C     Below energy threshold or above input data range
      if (maxDcrs.LE.0.) then
         IERR = 1
         RETURN
      endif


C     Determine direction of pion/eta
 102  CONTINUE

      COSPI=2*RLU(IDUM)-1

      if (IDNUCI.EQ.2212) then
         piDcrs = efgampipn(EGAM,cospi)
         pi0Dcrs = efgampi0p(EGAM,cospi)
      else if (IDNUCI.EQ.2112) then
         piDcrs = efgampimp(EGAM,cospi)
         pi0Dcrs = efgampi0n(EGAM,cospi)
      end if

      sumDcrs = piDcrs + pi0Dcrs

C     Below energy threshold or above input data range
      if (sumDcrs.LE.0. .or. maxDcrs.LE.0) then
         IERR = 1
         RETURN
      endif


c     Accept/reject weighting for angle distribution
      IF (sumDcrs/maxDcrs.LE.RLU(IDUM)) GOTO 102


c     Determine pion type
      dum = RLU(IDUM)
      if (dum.LT.pi0Dcrs/sumDcrs) then
         IDPION = 111
         IDNUCF = IDNUCI
         if (IDNUCI.EQ.2212) then
            IMODE = 1
         else
            IMODE = -1
         endif
      else
         if (IDNUCI.EQ.2212) then
            IDPION = 211
            IDNUCF = 2112
            IMODE = 2
         else if (IDNUCI.EQ.2112) then
            IDPION = -211
            IDNUCF = 2212
            IMODE = -2
         endif
      end if

      RETURN
      END
