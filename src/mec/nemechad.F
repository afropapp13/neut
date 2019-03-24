************************************************************************
*     --------------------------------------------
      SUBROUTINE NEMECHAD(Q0,Q3,PIN1,PIN2,IPMODE,POUT1,POUT2,IERR)
*     --------------------------------------------
*
*     (Purpose)
*       Set initial and final momenta of nucleons for 2p2h
*
*     (Input)
*       Q0     : Energy transferred from lepton to nucleon system
*       Q3     : 3-momentum transferred from lepton to nucleon system
*       IPMODE     : PID of all interacting particles
*
*     (Output)
*       PIN1(3)  : FERMI MOMENTUM for nucleon 1 ( GEV )
*       PIN2(3)  : FERMI MOMENTUM for nucleon 2 ( GEV )
*       POUT1(3) : Outgoing nucleon momentum(nucleon 1) ( GEV )
*       POUT2(3) : Outgoing nucleon momentum(nucleon 2) ( GEV )
*       IERR     : Error flag
*
*     (Creation Date and Author)
*       2013.04.17 ; P.Sinclair created losely based on neelsspff.F
*     
************************************************************************
      IMPLICIT NONE

#include "neutparams.h"
#include "vcwork.h"
#include "posinnuc.h"
#include "necard.h" 
     
      real PIN1(3),PIN2(3), POUT1(3), POUT2(3)
      REAL*4 PFABS1, PFABS2
      REAL*4 PFABSBST
      REAL*4 POUTLAB1, POUTLAB2
      REAL*4 HIT, MAXPROB
      REAL PI
      REAL costheta, phi
      real*4 INMASS1, INMASS2, OUTMASS1, OUTMASS2
      REAL Q0, Q3(3), W0, W3(3), W32, Q3MAG,INVMASSMAX,W3MIN
      REAL EBIND1, EBIND2
      REAL OUTMASS, INVMASS
      REAL BOOSTDIR(3)
      REAL*4 dum, fvfrmgev, rlu
      EXTERNAL fvfrmgev, rlu
      INTEGER IPMODE(6)
      REAL effrmgas
      REAL DUM1, DUM2, DUM3(3)
      REAL INTERACT_POS(3), INTERACT_RAD, SURFACE_MOM
      INTEGER POS_ITER, ACCEPT_POS, POS_ITER_LIM
      INTEGER ACCEPT_INITMOM
      INTEGER MOM_ITER, ACCEPT_MOM, MOM_ITER_LIM
      REAL PFMIN, PFMIN1, PFMIN2, EMIN, Q0MIN
      REAL SURFACE_MAX
      REAL Nnn,Nnp,Npp

      integer*4 ierr
      integer*4 i,j
      real*4 beta,gm
     

      IERR=0

      PI = 2.*ACOS(0.0)
      Q3MAG = SQRT(Q3(1)**2 + Q3(2)**2 + Q3(2)**2)

C Set limits for number of iterations over each variable until error is thrown
      POS_ITER_LIM=1000
      MOM_ITER_LIM=10000
      POS_ITER=0
      MOM_ITER=0


C Set isospin of nucleons. Choose np pair with probability based on combinatorics (i.e. number of nn/np/pp pairs)
C Fist count number of nn/np/pp pairs
      Nnn=0.5*(REAL(NUMBNDN)**2)
      Nnp=REAL(NUMBNDN*NUMBNDP)
      Npp=0.5*(REAL(NUMBNDP)**2)
C Then set isospin of nucleon which changes (defined by whether initial lepton is nu or nubar)
C Then set spectator isospin
      IF (IPMODE(1).gt.0) THEN
         IPMODE(2)=2112
         IPMODE(5)=2212
         HIT=RLU(DUM)
         IF(HIT>(Nnp/(Nnp+Nnn))) THEN
           IPMODE(3)=2112
           IPMODE(6)=2112
         ELSE
           IPMODE(3)=2212
           IPMODE(6)=2212
         ENDIF
      ELSE
         IPMODE(2)=2212
         IPMODE(5)=2112
         HIT=RLU(DUM)
         IF(HIT>(Nnp/(Nnp+Npp))) THEN
           IPMODE(3)=2212
           IPMODE(6)=2212
         ELSE
           IPMODE(3)=2112
           IPMODE(6)=2112
         ENDIF

      ENDIF
C Now count number of nn/np/pp pairs


C First set nucleon momentum to zero
      DO 101 i=1,3
        PIN1(i)=0
        PIN2(i)=0
        POUT1(i)=0
        POUT1(i)=0
 101  CONTINUE

      CALL MCMASS(IPMODE(2),INMASS1)
      CALL MCMASS(IPMODE(3),INMASS2)
      CALL MCMASS(IPMODE(5),OUTMASS1)
      CALL MCMASS(IPMODE(6),OUTMASS2)
      INMASS1=INMASS1*1E-3
      INMASS2=INMASS2*1E-3
      OUTMASS1=OUTMASS1*1E-3
      OUTMASS2=OUTMASS2*1E-3

C Set binding energy of two nucleons
      EBIND1=VNUINI
      EBIND2=VNUINI
C Find maximum nucleon momentum
      SURFACE_MAX = EFFRMGAS(DUM1,DUM2,0.0)*1E-3 
C Calculate maximum energy of hadronic system
      W0 = Q0 + SQRT(INMASS1**2+SURFACE_MAX**2) +
     $          SQRT(INMASS2**2+SURFACE_MAX**2) + EBIND1 + EBIND2
C Calculate minimum magnitude of hadronic state 3-momentum (corresponding to max invariant mass)
      W3MIN = SQRT(Q3(1)**2+Q3(2)**2+Q3(3)**2) - (2*SURFACE_MAX)

      IF(W3MIN<0) W3MIN=0
     
      IF(SQRT(W0**2-W3MIN**2).lt.(OUTMASS1+OUTMASS2)) GOTO 130

C Check that there's sufficient energy transfer to overcome binding energy and mass change
      Q0MIN = SQRT(OUTMASS1**2+SURFACE_MAX**2) 
     $     - SQRT(INMASS1**2+SURFACE_MAX**2)
     $     - EBIND1 - EBIND2
      IF(Q0.lt.Q0MIN) GOTO 150


      ACCEPT_POS=0
      POS_ITER=0
      DO WHILE((ACCEPT_POS.eq.0).and.(POS_ITER.lt.POS_ITER_LIM))

C First choose the interaction point. This defines the Pauli surface momentum.
        DUM1   =0.
        DUM3(1)=0.
        DUM3(2)=0.
        DUM3(3)=0.
C        CALL NERANBLL(DUM1,DUM3,DUM3,INTERACT_POS)
        CALL NERANBLL2(DUM1,DUM3,DUM3,INTERACT_POS)
        INTERACT_RAD = SQRT(INTERACT_POS(1)**2 
     $                    + INTERACT_POS(2)**2
     $                    + INTERACT_POS(3)**2)
        DUM1=0.
        DUM2=0.

        SURFACE_MOM = EFFRMGAS(DUM1,DUM2,INTERACT_RAD)*1E-3 

C Check that there can be sufficient final-state energy to allow 2 nucleons
C over fermi-momentum in the final state
        PFMIN =((SQRT(OUTMASS1**2+SURFACE_MOM**2)-Q0-EBIND1-EBIND2)**2)
     $            -(INMASS1**2)

        IF(PFMIN<0) PFMIN = 0
        PFMIN = SQRT(PFMIN)


C Check that this position is possible (suffucient surface momentum)
        W0 = Q0 + SQRT(INMASS1**2+SURFACE_MOM**2) +
     $            SQRT(INMASS2**2+SURFACE_MOM**2) + EBIND1 + EBIND2
C Calculate minimum magnitude of hadronic state 3-momentum (corresponding to max invariant mass)
        W3MIN = SQRT(Q3(1)**2+Q3(2)**2+Q3(3)**2) - (2*SURFACE_MOM)
        IF(W3MIN<0) W3MIN=0

        Q0MIN = SQRT(OUTMASS1**2+SURFACE_MOM**2) 
     $        - SQRT(INMASS1**2+SURFACE_MOM**2)
     $        - EBIND1 - EBIND2



        IF(PFMIN.gt.SURFACE_MOM) THEN
          ACCEPT_POS=0
        ELSE IF(W0**2-W3MIN**2.lt.((OUTMASS1+OUTMASS2)**2)) THEN
          ACCEPT_POS=0
        ELSE
          ACCEPT_POS = 1
        END IF 
        POS_ITER = POS_ITER+1

        IF(ACCEPT_POS.eq.1) THEN

C Pick random pF between 0 and PFMAX for both nucleons
C Pick according to flat probability for all states using FVFRMGEV
C Probability normalisation (find maximum in pdf-like distribution)
          MAXPROB = FVFRMGEV(SURFACE_MOM)

C To improve efficiency of pF throwing algorithm:
C Find minimum possible nucleon momentum with which event could be kinematically allowed



C Loop over various combinations of initial state nucleons and outgoing nucleons to find a viable combination.
C Failed combinations have viable initial state nucleons
C but produce one or mode final state nucleons under PFSURF.
C Continue loop while viable momenta not found and iteration less than limit
C or if too many initial nucleon tries.
          ACCEPT_MOM = 0
          MOM_ITER=0

          DO WHILE ((ACCEPT_MOM.eq.0).and.(MOM_ITER.lt.MOM_ITER_LIM))

            ACCEPT_INITMOM=0

C Loop over various pairs of nucleons to try to find a pair which, 
C when added to the momentum transfer from the hadronic system, will 
C allow the creation of two on-shell nucleons in the lab frame
C (total available energy > 2M + Ebind)

C Pick pF for 1st nucleon
            PFABS1 = 0
              HIT=1+MAXPROB
            DO WHILE ((HIT.gt.FVFRMGEV(PFABS1)))
              PFABS1 = (RLU(DUM)*(SURFACE_MOM-PFMIN)) + PFMIN
              HIT = RLU(DUM)*MAXPROB
            END DO

C Pick pF for 2nd nucleon
            PFABS2=0
            HIT=1
            DO WHILE ((HIT.gt.FVFRMGEV(PFABS2)))
              PFABS2 = (RLU(DUM)*(SURFACE_MOM-PFMIN)) + PFMIN
              HIT = RLU(DUM)*MAXPROB
            END DO

C Calculate energy of hadronic system
            W0 = Q0 + SQRT(INMASS1**2+PFABS1**2) + 
     $                SQRT(INMASS2**2+PFABS2**2) + EBIND1 + EBIND2

C See if it's possible to get sufficient invariant mass with this combination of momenta
            W3MIN = SQRT(Q3(1)**2+Q3(2)**2+Q3(3)**2) - PFABS1 - PFABS2
            IF(W3MIN.lt.0) W3MIN=0
              INVMASSMAX = SQRT(W0**2 - W3MIN**2)
  
            IF(INVMASSMAX.lt.(OUTMASS1+OUTMASS2)) THEN
              ACCEPT_INITMOM = 0

            ELSE IF(W0.lt.EMIN) THEN
              ACCEPT_INITMOM = 0

            ELSE
C Pick random direction for 1st nucleon
              phi=2.0*PI*RLU(DUM)
              costheta=(2*RLU(DUM))-1
              PIN1(1) = PFABS1*SIN(ACOS(costheta))*COS(phi)
              PIN1(2) = PFABS1*SIN(ACOS(costheta))*SIN(phi)
              PIN1(3) = PFABS1*costheta
  
C Pick random direction for 2nd nucleon
              phi=2.0*PI*RLU(DUM)
              costheta=(2*RLU(DUM))-1
              PIN2(1) = PFABS2*SIN(ACOS(costheta))*COS(phi)
              PIN2(2) = PFABS2*SIN(ACOS(costheta))*SIN(phi)
              PIN2(3) = PFABS2*costheta

C Calculate momentum of hadronic system
              W3(1) = Q3(1)+PIN1(1)+PIN2(1)
              W3(2) = Q3(2)+PIN1(2)+PIN2(2)
              W3(3) = Q3(3)+PIN1(3)+PIN2(3)
              W32   = W3(1)**2 + W3(2)**2 + W3(3)**2

              IF(((W0**2-W32)).lt.0) THEN 
                INVMASS=-1
              ELSE
                INVMASS=sqrt((W0*W0)-W32)
              END IF

              IF(INVMASS.gt.(OUTMASS1+OUTMASS2)) THEN
                ACCEPT_INITMOM=1
              ELSE
                ACCEPT_INITMOM=0
              END IF
            END IF

C If a valid pair has been picked, try to choose final momenta. If not, skip and re-choose interaction position
            IF(ACCEPT_INITMOM.eq.1) THEN

C Calculate magnitude of momentum in hadronic CoM frame
              PFABSBST = (1./(2.*(INVMASS)))*sqrt(
     $             (INVMASS**4)+(OUTMASS1**4)+(OUTMASS2**4)
     $            -( 2.*INVMASS**2  *OUTMASS1**2 )
     $            -( 2.*INVMASS**2  *OUTMASS2**2 )
     $            -( 2.*OUTMASS1**2 *OUTMASS2**2 ) )

C Pick random direction for nucleon axis
              phi=2.0*PI*RLU(DUM)
              costheta=(2.*RLU(DUM))-1
              POUT1(1) =  PFABSBST*SIN(ACOS(costheta))*COS(phi)
              POUT1(2) =  PFABSBST*SIN(ACOS(costheta))*SIN(phi)
              POUT1(3) =  PFABSBST*costheta
              POUT2(1) = -POUT1(1)
              POUT2(2) = -POUT1(2)
              POUT2(3) = -POUT1(3)

C Calculate boost to get outgoing hadron momenta in lab frame
              BETA=sqrt(W32)/W0
              GM=1./SQRT(1.-BETA**2)
              BOOSTDIR(1) = W3(1)/sqrt(W32) 
              BOOSTDIR(2) = W3(2)/sqrt(W32) 
              BOOSTDIR(3) = W3(3)/sqrt(W32) 

C Boost hadron momenta into lab frame
              CALL MCVECBST(POUT1,OUTMASS1,BOOSTDIR,GM)
              CALL MCVECBST(POUT2,OUTMASS2,BOOSTDIR,GM)
    
              POUTLAB1 = SQRT(POUT1(1)**2 + POUT1(2)**2 + POUT1(3)**2)
              POUTLAB2 = SQRT(POUT2(1)**2 + POUT2(2)**2 + POUT2(3)**2)

C Both nucleons must be above the Fermi-surface momentum if nucleon is to be Pauli-blocked

              IF ((POUTLAB1.gt.SURFACE_MOM)
     $            .and.(POUTLAB2.gt.SURFACE_MOM))
     $        THEN
                ACCEPT_MOM=1

              ELSE
                ACCEPT_MOM=0
              END IF

            ENDIF
          MOM_ITER = MOM_ITER+1
          END DO


C No good hadron kinematics found. Change position
          IF((ACCEPT_POS.eq.1).and.(ACCEPT_MOM.ne.1)) THEN
            ACCEPT_POS=0
          ELSE IF((ACCEPT_POS.eq.1).and.(ACCEPT_MOM.eq.1)) THEN

C Store interaction position
            DO 1030 i=1,4
              DO 1020 j=1,3
                POSNUC(j,i)=INTERACT_POS(j)
 1020         CONTINUE
 1030       CONTINUE

            GOTO 140
          END IF
        END IF
      END DO


C Generic 'no match found'
      IF(ACCEPT_POS.eq.0) THEN
 120    WRITE(6,121)
 121    FORMAT(' *** ERROR AT NEMECHAD: Lepton kinematics rejected'
     $         ' (too many vertices tried) ***')
        GOTO 111
      ENDIF

 130  WRITE(6,131)
 131  FORMAT(' *** ERROR AT NEMECHAD: No vertices possible'
     $       '(insufficient invariant mass of initial state) ***')
      GOTO 111

 150  WRITE(6,151)
 151  FORMAT(' *** ERROR AT NEMECHAD: No vertices possible '
     $       '(energy transfer too small) ***')
      GOTO 111

 111  DO 125 i=1,3
        PIN1(i)=0
        PIN2(i)=0
        POUT1(i)=0
        POUT2(i)=0
 125  CONTINUE
      IERR=1
      RETURN

C140  WRITE(6,141)
C141  FORMAT(' NEMECHAD: HADRON KINEMATICS SUCCESSFULY FOUND ')
 140  CONTINUE
      IERR=0
      RETURN

      END
