*     -----------------------------------------------------
      SUBROUTINE EFTRACE(APN,AP,KIN,XIN,PIN,NABS,NEL,NINEL,KFI,PFI)
*     -----------------------------------------------------
*     
*     ( purpose )
C        TRACE PION ZIGZAG PATH TO OUTSIDE OF NUCLEI OR TO ABSORPTION
*
*     ( input ) 
C         KIN  : PARTICLE KIND OF INITIAL PION
C                KIN= 211   PI+
C                KIN= 111   PI0
C                KIN=-211   PI-
C         XIN  : STARTING POINT OF PION  (FERMI)
C         PIN  : INITIAL MOMENTUM OF PION  (MEV)
C         APN  : # OF NUCLEONS  (DUMMY)
C         AP   : # OF PROTONS   (DUMMY)
*
*     ( output ) 
C          KFI : PARTICLE KIND OF FINAL PION
C                KFI=+211  PI+
C                KFI= 111  PI0
C                KFI=-211  PI-
C         NABS : ABSORPTION FLAG
C                NABS=0  PION ABSORBED
C                NABS=-1  PION ESCAPE OUTSIDE NUCLEI
C         NEL  : NUMBER OF ELASTIC SCATTERING BEFORE ESCAPE OR ABSORB
C         NINEL: NUMBER OF INELASTIC SCATTERING BEFORE ESCAPE OR ABSORB
C         PFI  : FINAL MOMENTUM OF PION  (MEV)
C
C     ( Updating record )
C         2005.03.16 ; M.Ishitsuka
C                        Change to lower the boundary energy of scattering 
C                        850 -> 500 MeV
C         2010.06    ; P.de Perio
C                        - Save absorbed pi momentum for nucleon ejection
C                        - Add inelastic scattering for p>500 MeV/c
C                        - Separate NZIG scatter counter into NEL elastic 
C                          and NINEL inelastic counters
C                        - Choose nucleon type depending on nuclear 
C                          density for that nucleon (p>500)
C                        - Add model for high momentum (p>500): 
C                            Distinguish between p/n cross sections and 
C                            separate charge exchange
C                        - Use charge of current pion at high momentum
C                          instead of initial pion:
C                              NQPAI=NECHARGE(KIN) -> NQPAI=NECHARGE(IP)
C                        - Add tuned qe/abs paths for p<500, including
C                          absorption at p>500 MeV/c; smooth transition
C                          in model select between 400-500 MeV/c
c
C        2010.09     ; P. de Perio
C                        - Store intermediate vertex and particle info
C
C        2010.10     ;   - Fix bug: decrement FSIHIST vertex and part. info
C                                   when EFHADPRO fails
C                        - Move scaling parameters into probbility calculating 
C                          functions
C
c

C      implicit none

C#define DEBUG

#include "neutparams.h"
      REAL XIN(3),PIN(3),PFI(3),X(3),P(3),DIR(3),INITDIR(3),PSCAT(3)
      REAL RLU
      INTEGER IERROR,IDUM,ISEED
      COMMON/RNSEED/ ISEED
#include<efpion.h>
      REAL AMPI,AMNUCL,PCUT
      DATA AMPI/139./
      DATA AMNUCL/939./
      DATA PCUT/500./
      
      INTEGER I, iNuc
      REAL APN,AP,HIT
      REAL RADIUS,DL,INITABSP,ABSP,DXXY,RMAX,dt
      INTEGER IPF,NQN,NQPAI,NECHARGE,IP,KIN,KFI,NABS,NEL,NINEL
      REAL PESC,PINEL,PCOH,PABS,PQE,PCX

      COMMON/YHDEBUG/LTINEL,LTABS,LTCOH
      integer LTINEL,LTABS,LTCOH

C     Error flag from EFINELS
      integer ierr
      COMMON/EFLAG/IERR

      REAL MOMFERM,NUCDIR(3),DUM,PNI(3)
      REAL TOTENUCL,ABSBETANUCL,BETANUCL(3),GAMMANUCL,PBST(3)
      REAL INVNUCDIR(3),ABSPBST

      integer ireaction

C     FSI History vars
      INTEGER IVERT,CURRVERT,VERTOFFSET
      parameter (epsilon=0.001)

      real MODFRAC
C function
      REAL GAMFCT
      EXTERNAL GAMFCT

#include "necard.h"
#include "fsihist.h"
#include "efcoulforce.h"

C
C --  PARTICLE CHECK   1995.03.23 K.K.
C
C-- Always copy PIN to PFI
      CALL UCOPY(PIN,PFI,3)

      IF ( (ABS(KIN).NE.211.AND.KIN.NE.111) 
     $     .or. NENEFO16.NE.0 ) THEN
C         KIN=KFI
         KFI=KIN
         NABS=-1
         NEL=0
         NINEL=0
C         CALL UCOPY(PIN,PFI,3)
         RETURN
      ENDIF
C
C
      IP=KIN
      DO I=1,3
         X(I)=XIN(I)
         P(I)=PIN(I)
      end do
      NEL=0
      NINEL=0
C     -- SET CALCULATING LIMIT OF RADIUS
      RMAX=CC2
      DXXY=RMAX/100.

      PCX=0


C     Initialize FSI history
      CURRVERT=0
      if (NVERT.eq.0) then
         CALL FSIFILLVERT(X,0)
         CURRVERT=NVERT

      else
         DO IVERT=1,NVERT
            if (ABS(X(1)-POSVERT(1,IVERT)) .gt. epsilon) cycle
            if (ABS(X(2)-POSVERT(2,IVERT)) .gt. epsilon) cycle
            if (ABS(X(3)-POSVERT(3,IVERT)) .gt. epsilon) cycle
            CURRVERT=IVERT

#ifdef DEBUG
            write(*,*) "Matched vertex: ", CURRVERT, IFLGVERT(CURRVERT), POSVERT(1,IVERT), POSVERT(2,IVERT), POSVERT(3,IVERT)
#endif

            exit
         end do
      end if

      if (CURRVERT.eq.0) then

#ifdef DEBUG
         write (*,*) "epsilon = ", epsilon
         
         DO IVERT=1,NVERT
            write (*,*) "Vertex ", IVERT, ": ", 
     &           ABS(X(1)-POSVERT(1,IVERT)), 
     &           ABS(X(2)-POSVERT(2,IVERT)), 
     &           ABS(X(3)-POSVERT(3,IVERT)),
     &           IFLGVERT(IVERT)
         end do

         write(*,*) "eftrace.F: Warning no parent vertex for 
     $ particle at ", X(1), X(2), X(3)
#endif

         CALL FSIFILLVERT(X,0)
         CURRVERT=NVERT

      end if

      if (FEFCOUL.eq.1) then
C     -- Load integrated density table for Coulomb force
         call loadIntDensTable

C     Get index in nucleusZ array for Coulomb force 
         do iNuc = 1, nNuclei
            if (nucZ(iNuc).eq.NUMBNDP) exit
         end do
         
         if (iNuc.gt.nNuclei) then
            write(*,*) "EFTRACE Error: Unimplemented Z value = ",NUMBNDP
            stop
         end if
      end if


C     -- NOW TRACE
C     
   10 ABSP=SQRT(P(1)**2+P(2)**2+P(3)**2)
      IF(ABSP.GT.0.001) GO TO 14
      ABSP=  0.001
      DIR(1)=0.001
      DIR(2)=0.0
      DIR(3)=0.0
C     T.KAJITA 12-NOV-85   FOR AVOIDING ZERO DIVIDE
C     
   14 DIR(1)=P(1)/ABSP
      DIR(2)=P(2)/ABSP
      DIR(3)=P(3)/ABSP

C     Store initial particle direction and momentum prior to Coulomb corrections
      do i = 1, 3
         INITDIR(i) = DIR(i)
      end do
      INITABSP = ABSP

  100 X(1)=X(1)+DIR(1)*DXXY
      X(2)=X(2)+DIR(2)*DXXY
      X(3)=X(3)+DIR(3)*DXXY
      DL=SQRT(X(1)**2+X(2)**2+X(3)**2)
C
C-- not necessary to consider fermi motion for lower energy:
C-- momentum is "incoming pion relative to the Oxygen.
C
C M.Ishitsuka 2005.03.16 ; Consider Fermi motion
C
      radius = DL
      IF (NEFKINVER.eq.0) THEN
         call rnferm(momferm,dum)
      else IF (NEFKINVER.eq.1) THEN
         CALL EFRNMOM(radius,momferm)
      else
         write(*,*) 'Invalid NEFKINVER =',NEFKINVER
         stop
      endif
      CALL RNDIR(NUCDIR)
      TOTENUCL=SQRT(AMNUCL**2+MOMFERM**2)
      DO I=1,3
         PNI(I)=MOMFERM*NUCDIR(I)
         BETANUCL(I)=PNI(I)/TOTENUCL
         INVNUCDIR(I)=-1.0*NUCDIR(I)
      ENDDO
      ABSBETANUCL=SQRT(BETANUCL(1)**2+BETANUCL(2)**2+BETANUCL(3)**2)
      GAMMANUCL=GAMFCT(ABSBETANUCL)
      DO I=1,3
         PBST(I)=P(I)
      ENDDO
      CALL MCVECBST(PBST,AMPI,INVNUCDIR,GAMMANUCL)
      ABSPBST=SQRT(PBST(1)**2+PBST(2)**2+PBST(3)**2)

C  100 X(1)=X(1)+DIR(1)*DXXY
C      X(2)=X(2)+DIR(2)*DXXY
C      X(3)=X(3)+DIR(3)*DXXY
C      DL=SQRT(X(1)**2+X(2)**2+X(3)**2)

C      write(*,*) "EFTRACE(): ",X(1), X(2), X(3), DL

      IF(DL.GE.RMAX)GO TO 1000
   
C     For smooth transition between energy regions
      MODFRAC = (ABSP - 400)/100.

C++   LOW MOMENTUM ( P < PCUT (500MEV/C) )
      IF( (NENEFMODL.EQ.0 .and. ABSP.LE.PCUT) .or.
     $    (NENEFMODL.EQ.1 .and. RLU(IDUM).GE.MODFRAC) ) THEN
C     $    (NENEFMODL.EQ.1 .and. ABSP.LE.450) ) THEN
        
         CALL EFPQEAB(ABSP,DL,DXXY,PESC,PQE,PABS,PCOH)

         HIT=RLU(IDUM)

C        No interaction in this step
         IF(HIT.LT.PESC)GO TO 4000

#ifdef DEBUG
         write(*,*) "After EFPQEAB: CURRVERT = ", CURRVERT
#endif

C        ABSORPTION
         IF(HIT.LT.(PESC+PABS)) then
            LTABS=LTABS+1
            
            CALL FSIFILLVERT(X,3)
            CALL FSIFILLPART(INITDIR,INITABSP,ABSPBST,IP,CURRVERT,NVERT)

            GO TO 2000
         endif

C        FORWARD (ELASTIC) SCATTER
         IF(HIT.LT.(PESC+PABS+PCOH)) then
            LTCOH=LTCOH+1

            CALL FSIFILLVERT(X,9)
            CALL FSIFILLPART(INITDIR,INITABSP,ABSPBST,IP,CURRVERT,NVERT)

            GO TO 3000
         endif

C        QUASI-ELASTIC SCATTER
         LTINEL=LTINEL+1

         GO TO 2010

C     
C++   HIGH MOMENTUM ( P > PCUT (500MEV/C) )
      ELSE
         
C        This is a bug in previous official versions
         IF (NENEFMODH.EQ.0) THEN
            NQPAI=NECHARGE(KIN)
         ELSE
            NQPAI=NECHARGE(IP)
         END IF

         NQN=0
         ireaction=1
C
C--- In Nucleus (not free protons i.e. hydrogens)
C

C---  Select nucleon depending on p/n density
c        Only iso-scalar nuclei considered for previous official versions
         IF (NENEFMODH.EQ.0) THEN
            IF(RLU(IDUM).LT.0.5) NQN=1
         ELSE
            IF(RLU(IDUM).LT.(float(NUMBNDP)/float(NUMATOM))) NQN=1
         END IF

C
C--   Here use ABSPBST (for consistency)

         CALL EFPQEABH(ABSPBST,DL,DXXY,NQPAI,NQN,PESC,PQE,PINEL,PCOH,PCX)

         HIT=RLU(IDUM)

C        No interaction in this step
         IF(HIT.LT.PESC) GO TO 4000
         
C        Quasi-elastic scatter
         IF(HIT.LT.(PESC+PQE)) then
            CALL FSIFILLVERT(X,80)
            CALL FSIFILLPART(INITDIR,INITABSP,ABSPBST,IP,CURRVERT,NVERT)

            GO TO 2030
         end if

C        Forward (elastic) scatter
         IF(HIT.LT.(PESC+PQE+PCOH)) then
            CALL FSIFILLVERT(X,90)
            CALL FSIFILLPART(INITDIR,INITABSP,ABSPBST,IP,CURRVERT,NVERT)

            GO TO 3000
         end if

C        Charge exchange
         IF(HIT.LT.(PESC+PQE+PCOH+PCX)) then
            CALL FSIFILLVERT(X,40)
            CALL FSIFILLPART(INITDIR,INITABSP,ABSPBST,IP,CURRVERT,NVERT)

            GO TO 2020
         end if
C     
C     -- INTERACTION WHICH MULTIPLIES THE NUMBER OF PARTICLE
C     
         CALL EFHADPRO(PBST,X,NQPAI,NQN,NUCDIR,GAMMANUCL,IERROR)

C        Error: Revert FSIHIST common block
         IF(IERROR.NE.0) then
            GO TO 100
         END IF

         CALL FSIFILLVERT(X,70)
         CALL FSIFILLPART(INITDIR,INITABSP,ABSPBST,IP,CURRVERT,NVERT)

         GO TO 2000
      END IF
C     
C     
C     
C     -- INELASTIC SCTTERING (p < 500)
 2010 radius = DL 
      CALL EFINELS(IP,P,IPF,PSCAT,radius)

C     Error in EFINELS = no interaction
      IF(IERR.NE.0) then

#ifdef DEBUG
         write(*,*) "EFINELS Error: CURRVERT = ", CURRVERT
#endif         
         go to 4000

      else

C     Charge-exchange occured 
         IF (IP.NE.IPF) then
            CALL FSIFILLVERT(X,4)

         else
            CALL FSIFILLVERT(X,8)
            
         end if

         CALL FSIFILLPART(INITDIR,INITABSP,ABSPBST,IP,CURRVERT,NVERT)

         NINEL=NINEL+1

         IP=IPF
         P(1)=PSCAT(1)
         P(2)=PSCAT(2)
         P(3)=PSCAT(3)

         GO TO 10

      end if      


C     -- CHARGE EXCHANGE (p > 500)
 2020 ireaction=3
C     -- INELASTIC SCTTERING (p > 500)
 2030 radius = DL

      IF (NENEFKINH.eq.1) THEN
         CALL EFHINELS_PS(IP,P,IPF,PSCAT,NQN,radius,ireaction)
      ELSE
         CALL EFHINELS_ISO(IP,P,IPF,PSCAT,NQN,radius,ireaction)
      ENDIF
      IP=IPF
      P(1)=PSCAT(1)
      P(2)=PSCAT(2)
      P(3)=PSCAT(3)
      NINEL=NINEL+1
      GO TO 10
C     
C     -- PION ESCAPED
C     
 1000 CALL FSIFILLVERT(X,-1)

C     Store escaped particle info
#ifdef DEBUG
      write(*,*) "Escape: CURRVERT = ", CURRVERT, ", NVERT = ",NVERT
#endif
      CALL FSIFILLPART(INITDIR,INITABSP,ABSPBST,IP,CURRVERT,NVERT)

      KFI=IP
      PFI(1)=P(1)
      PFI(2)=P(2)
      PFI(3)=P(3)
      NABS=-1
      RETURN
C     
C     -- PION ABSORBED
C     
 2000 KFI=IP
      PFI(1)=P(1)
      PFI(2)=P(2)
      PFI(3)=P(3)
      NABS=0
      RETURN
C     
C     -- COHERENT SCATTERING
C     
 3000 CALL EFCOHSCT(PBST,PSCAT)
      CALL MCVECBST(PSCAT,AMPI,NUCDIR,GAMMANUCL)
      P(1)=PSCAT(1)
      P(2)=PSCAT(2)
      P(3)=PSCAT(3)
      NEL=NEL+1
      GO TO 10

C
C     -- No interaction, apply Coulomb force
C
C     Calculate time of step in units of (seconds * speed of light * fm)
C     See neututils/calc_nuc_dens/coulomb_factor.xls
 4000 if (FEFCOUL.eq.1) then
         dt = DXXY / sqrt(1 - 1/((AMPI**2+ABSP**2)/AMPI))

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
      
      go to 100
      
      END
      
