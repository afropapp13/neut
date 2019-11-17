*     -----------------------------------------------------
      DOUBLE PRECISION FUNCTION EVPIPROB()
*     -----------------------------------------------------
*
*     ( purpose )
C        Calculate probability of pion final state interaction event.
C        MUST BE KEPT SYNCHRONIZED WITH eftrace.F
*
*     ( input )
c        FSIHIST common block
c        NEFFPR common block
*
*     ( output )
C        Probability of the event
C
C     ( Updating record )
C         2010.09.28 ; P.de Perio
C                        - Creation
C         2010.10.13 ; Add condition to escape infinite loop
c
C

C#define DEBUG

      implicit none

      DOUBLE PRECISION PROB

      DOUBLE PRECISION PPI(3)

      INTEGER I,J, iNuc

      DOUBLE PRECISION POSDIFF

#include "vcwork.h"
#include "posinnuc.h"
#include "neutparams.h"
      REAL X(3),DIR(3), P(3)
      INTEGER ISEED
      COMMON/RNSEED/ ISEED
#include<efpion.h>
      REAL AMPI,PCUT
      DATA AMPI/139./
      DATA PCUT/500./

      REAL DL,ABSP,ABSPBST,DXXY,RMAX, dt
      INTEGER NQN,NQPAI,NECHARGE,IP
      REAL PESC,PINEL,PCOH,PABS,PQE,PCX
      REAL PESCP,PESCN,PINELP,PINELN,PCOHP,PCOHN
      REAL PQEP,PQEN,PCXP,PCXN
      REAL PESCL,PESCH,PESCHP,PESCHN
      REAL MODFRAC,NUCFRAC
      REAL DXXYFACTOR ! For fine tuning the number of steps between vertices
      DATA DXXYFACTOR/0.99/

      INTEGER VERBOSE

#include "necard.h"
#include "fsihist.h"
#include "efcoulforce.h"

      VERBOSE = 0

      EVPIPROB = 1
      FSIPROB  = 1

C     No intermediate particle
      IF (NVERT.le.1 .or. NVCVERT.lt.1) return

C     Not bound target
      IF (IBOUND.eq.0) return

C     Not tuned FSI model
      IF (NENEFMODL.EQ.0 .or. NENEFMODH.EQ.0) THEN
         write(*,*) "EVPIPROB Warning: Attemping to call with NENEFMODL!=0 and NENEFMODH!=0"
         return
      END IF

c     Weight for nucleon densities for high energy calculation
C     TODO: Need to test behaviour for non-isoscalar nuclei

C     Actual simulation uses:
C               NQN=0
C               IF(RLU(IDUM).LT.(float(NUMBNDP)/float(NUMATOM))) NQN=1
      NUCFRAC = float(NUMBNDP)/float(NUMATOM)


      if (FEFCOUL.eq.1) then
C     -- Load integrated density table for Coulomb force
         call loadIntDensTable

C     Get index in nucleusZ array for Coulomb force
         do iNuc = 1, nNuclei
            if (nucZ(iNuc).eq.NUMBNDP) exit
         end do

         if (iNuc.eq.nNuclei) then
            write(*,*) "EFTRACE Error: Unimplemented Z value = ",NUMBNDP
            stop
         end if
      end if

C     Array bounds check
      if (NVCVERT.gt.MAXVCVERT) then
         write(*,*) "EVPIPROB Warning: NVCVERT > ", MAXVCVERT
         NVCVERT=MAXVCVERT
      end if
      if (NVERT.gt.MAXVERT) then
         write(*,*) "EVPIPROB Warning: NVERT > ", MAXVERT
      end if

C     Loop over all intermediate particles
      DO J=1,NVCVERT

         IP = IPVERT(J)

C     Check particle type
         IF (ABS(IP).NE.211.AND.IP.NE.111) then
            write(*,*) "EVPIPROB Error: Particle #",J," is not a pion (",IP,")"
            cycle
         END IF

C     SKDETSIM implemention of NEUT cascade uses same FSI common block to
C     store elastic interactions on free hydrogen and nucleus so skip these
C     vertices here
         IF (IFLGVERT(IVERTI(J)).ge.100) then
C            write(*,*) "EVPIPROB: Skipping non-eftrace interaction flag = ",IFLGVERT(IVERTI(J))
            cycle
         end if

C     Initialize variables from FSIHIST
         DO I=1,3
            X(I)=POSVERT(I,IVERTI(J))
            DIR(I)=DIRVERT(I,J)
         end do
         ABSP=ABSPVERT(J)
         ABSPBST=ABSTPVERT(J)

         RMAX=CC2
         DXXY=RMAX/100.

C     For high-energy start with neutron target and alternate each step
         NQN=1

C     -- NOW TRACE
C
 100     X(1)=X(1)+DIR(1)*DXXY
         X(2)=X(2)+DIR(2)*DXXY
         X(3)=X(3)+DIR(3)*DXXY
         DL=SQRT(X(1)**2+X(2)**2+X(3)**2)

C         write(*,*) "EVPIPROB(): ",X(1), X(2), X(3), DL

C        This should not happen, unless particle direction doesn't correspond
C        between two vertices (bug). Stop routine to prevent infinite loop.
         IF(DL.GE.RMAX*1.1) then
            write(*,*) "EVPIPROB WARNING: Final vertex not found for current
     &                  particle direction. Returning value = 1"
            EVPIPROB = 1
            RETURN
         END IF


C        This alternates between proton/neutron for each step
C        if not using NUCFRAC weight above.
C         IF (NQN.eq.1) then
C            NQN=0
C         ELSE IF (NQN.eq.0) then
C            NQN=1
C         END IF

C        Calculate current distance from final vertex
         POSDIFF=0;
         DO I=1,3
            POSDIFF = POSDIFF + (POSVERT(I,IVERTF(J)) - X(I))**2
         END DO



C++   LOW MOMENTUM ( P < PCUT (500MEV/C) )

C     EFTRACE uses the following condition for smooth model transition:
C
         MODFRAC = (ABSP - 400)/100.
C     IF( (NENEFMODL.EQ.0 .and. ABSP.LE.PCUT) .or.
C     $    (NENEFMODL.EQ.1 .and. RLU(IDUM).GE.MODFRAC) ) THEN

         IF (MODFRAC.LT.0) then
            MODFRAC=0
         else if (MODFRAC.GT.1) then
            MODFRAC=1
         endif

C     So let's calculate escape probability for high and low energy
C     models and use the weighted (by MODFRAC) sum
         IF (NENEFMODL.EQ.1) THEN

C           Intermediate step
            IF (POSDIFF .ge. (DXXY*DXXYFACTOR)**2) THEN
               CALL EFPQEAB(ABSP,DL,DXXY,PESCL,PQE,PABS,PCOH)

               CALL EFPQEABH(ABSPBST,DL,DXXY,NQPAI,0,PESCHN,
     &              PQE,PINEL,PCOH,PCX)
               CALL EFPQEABH(ABSPBST,DL,DXXY,NQPAI,1,PESCHP,
     &              PQE,PINEL,PCOH,PCX)

               PESCH = NUCFRAC*PESCHP + (1-NUCFRAC)*PESCHN

               PESC = (1-MODFRAC)*PESCL + MODFRAC*PESCH

               EVPIPROB = EVPIPROB * PESC

#ifdef DEBUG
               if (EVPIPROB.le.0) write(*,*) "EVPIPROB: 1. ", J, IP, ABSP
#endif

               if (VERBOSE.eq.3) write (*,*) "Esc. step (blend) = ",EVPIPROB

               GO TO 4000
C           Vertex step
            ELSE
               IF (IFLGVERT(IVERTF(J)).ge.10) then
                  goto 300
               ELSE IF (IFLGVERT(IVERTF(J)).lt.10) then
                  goto 200
               END IF

            END IF


         ELSE IF (NENEFMODL.EQ.0 .and. ABSP.LE.PCUT) then

C     However, when comparing probabilities, let's pick one model for
C     the current particle (low energy if the final vertex was a low
C     energy interaction). ERROR: Does not conserve unitarity in
C     transition region:
c           (NENEFMODL.EQ.1 .and. ABSP.LE.PCUT .and.
c      $     IFLGVERT(IVERTF(J)).LT.10) then

C     Or the following which does conserve unitarity, but we lose smoothness
C     $        (NENEFMODL.EQ.1 .and. ABSP.LE.450) ) then

 200        CALL EFPQEAB(ABSP,DL,DXXY,PESC,PQE,PABS,PCOH)

C           Intermediate step
            IF (POSDIFF .ge. (DXXY*DXXYFACTOR)**2) THEN
               EVPIPROB = EVPIPROB * PESC

#ifdef DEBUG
               if (EVPIPROB.le.0) write(*,*) "EVPIPROB: 2. ", J, IP, ABSP
#endif

               if (VERBOSE.eq.3) write (*,*) "Esc. step (lowe) = ",EVPIPROB

               GO TO 4000
C           Vertex step
            else
               IF (IFLGVERT(IVERTF(J)).eq.-1) then
                  EVPIPROB = EVPIPROB * PESC
               ELSE IF (IFLGVERT(IVERTF(J)).eq.3) then
                  EVPIPROB = EVPIPROB * PABS
               ELSE IF (IFLGVERT(IVERTF(J)).eq.4 .or.
     &                  IFLGVERT(IVERTF(J)).eq.8) then

                  EVPIPROB = EVPIPROB * PQE

C           Charge exchange probability determined within QE channel
                  PPI(1) = ABSP
                  PPI(2) = 0
                  PPI(3) = 0

#ifdef DEBUG
                  write(*,*) "EVPIPROB: Before EVSETPF, EVPIPROB =  ",EVPIPROB
#endif

                  CALL EVSETPF(IP,PPI,DL,IFLGVERT(IVERTF(J)),PROB)

                  EVPIPROB = EVPIPROB * PROB

#ifdef DEBUG
                  write(*,*) "EVPIPROB: After EVSETPF ", PROB, EVPIPROB
#endif


               ELSE IF (IFLGVERT(IVERTF(J)).eq.9) then
                  EVPIPROB = EVPIPROB * PCOH

C              Should be a high energy calculation
               ELSE IF (IFLGVERT(IVERTF(J)).ge.10) then
                  goto 300

               ELSE
                  write(*,*) "EVPIPROB: Unknown low interaction mode =",
     &                 IFLGVERT(IVERTF(J))
               end if

#ifdef DEBUG
               if (EVPIPROB.le.0) write(*,*) "EVPIPROB: 3. ", J, IP, ABSP, IFLGVERT(IVERTF(J)), PPI(1), PPI(2), PPI(3), DL, EVPIPROB, PROB
#endif
               if (VERBOSE.eq.3) write (*,*) "Int. step (lowe) = ",IFLGVERT(IVERTF(J)),EVPIPROB


               GO TO 10
            end if



C
C++   HIGH MOMENTUM ( P > PCUT (500MEV/C) )
         ELSE

C           Charge of pion
 300        NQPAI=NECHARGE(IP)
C
C---  In Nucleus (not free protons i.e. hydrogens)
C

C---  Select nucleon depending on p/n density

C           Weight probabilities by nucleon density
            IF (POSDIFF.ge.(DXXY*DXXYFACTOR)**2) THEN

 310           CALL EFPQEABH(ABSPBST,DL,DXXY,NQPAI,1,PESCP,PQEP,PINELP,
     $                       PCOHP,PCXP)
               CALL EFPQEABH(ABSPBST,DL,DXXY,NQPAI,1,PESCN,PQEN,PINELN,
     $                       PCOHN,PCXN)

               PESC = NUCFRAC*PESCP + (1-NUCFRAC)*PESCN
               PQE = NUCFRAC*PQEP + (1-NUCFRAC)*PQEN
               PINEL = NUCFRAC*PINELP + (1-NUCFRAC)*PINELN
               PCOH = NUCFRAC*PCOHP + (1-NUCFRAC)*PCOHN
               PCX = NUCFRAC*PCXP + (1-NUCFRAC)*PCXN

               goto 320

C           Some interaction
            ELSE

c              Charge exchange vertex step (select appropriate nucleon charge
C              for calculation of cross section)
               if ( IFLGVERT(IVERTF(J)).eq.40 ) then
                  IF(NQPAI.EQ.1) THEN
                     NQN = 0
                  ELSE IF(NQPAI.EQ.-1) THEN
                     NQN = 1
C                 For pi0, cross section is same on proton or neutron
                  ELSE IF(NQPAI.EQ.0) THEN
                     NQN = 1
                  END IF

                  CALL EFPQEABH(ABSPBST,DL,DXXY,NQPAI,NQN,PESC,PQE,PINEL,PCOH,PCX)

C              non-CX interaction
               else
                  goto 310

               endif

            END IF

C           Intermediate step
 320        IF (POSDIFF .ge. (DXXY*DXXYFACTOR)**2) THEN
               EVPIPROB = EVPIPROB * PESC

#ifdef DEBUG
               if (EVPIPROB.le.0) write(*,*) "EVPIPROB: 4. ", J, IP, ABSP
#endif
               if (VERBOSE.eq.3) write (*,*) "Esc. step (hie) = ",EVPIPROB

               GO TO 4000
C     Vertex step
            else
               IF (IFLGVERT(IVERTF(J)).eq.-1) then
                  EVPIPROB = EVPIPROB * PESC
C               ELSE IF (IFLGVERT(IVERTF(J)).eq.30) then
C                  EVPIPROB = EVPIPROB * PABS
               ELSE IF (IFLGVERT(IVERTF(J)).eq.40) then
                  EVPIPROB = EVPIPROB * PCX
               ELSE IF (IFLGVERT(IVERTF(J)).eq.70) then
                  EVPIPROB = EVPIPROB * PINEL
               ELSE IF (IFLGVERT(IVERTF(J)).eq.80) then
                  EVPIPROB = EVPIPROB * PQE
               ELSE IF (IFLGVERT(IVERTF(J)).eq.90) then
                  EVPIPROB = EVPIPROB * PCOH

C              Should be a low energy calculation
               ELSE IF (IFLGVERT(IVERTF(J)).lt.10) then
                  goto 200

               ELSE
                  write(*,*) "EVPIPROB: Unknown hi interaction mode =",
     &                 IFLGVERT(IVERTF(J))
               end if

#ifdef DEBUG
               if (EVPIPROB.le.0) write(*,*) "EVPIPROB: 5. ", J, IP, ABSP, IFLGVERT(IVERTF(J))
#endif
               if (VERBOSE.eq.3) write (*,*) "Int. step (hie) = ",
     &                           IFLGVERT(IVERTF(J)),EVPIPROB

               GO TO 10
            end if

         end if

 10   end do

      goto 5000

C
C     -- No interaction, apply Coulomb force
C
C     Calculate time of step in units of (seconds * speed of light * fm)
C     See neututils/calc_nuc_dens/coulomb_factor.xls
 4000 if (FEFCOUL.eq.1) then
         dt = DXXY / sqrt(1 - 1/((AMPI**2+ABSP**2)/AMPI))

         do i = 1, 3
            P(i) = DIR(i)*ABSP
         end do

         CALL EFCOULFORCE(P,IP,X,DL,dt,iNuc)

         ABSP=SQRT(P(1)**2+P(2)**2+P(3)**2)
         IF(ABSP.GT.0.001) then
            DIR(1)=P(1)/ABSP
            DIR(2)=P(2)/ABSP
            DIR(3)=P(3)/ABSP
         else
            ABSP=  0.001
            DIR(1)=0.001
            DIR(2)=0.0
            DIR(3)=0.0
         end if
      end if

      goto 100

 5000 if (EVPIPROB .le. 0) then
         write(*,*) "Warning: EVPIPROB<=0 (",EVPIPROB, "), returning 1"
         EVPIPROB=1
      else if (EVPIPROB .ne. EVPIPROB) then
         write(*,*) "Warning: EVPIPROB=NaN (",EVPIPROB, "), returning 1"
         EVPIPROB=1
      end if

      FSIPROB = EVPIPROB

      RETURN
      END
