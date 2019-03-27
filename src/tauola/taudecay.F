***********************************************************************
*     ----------------------------
      SUBROUTINE TAUDECAY(IP,PTAU)
*     ----------------------------
*     
*     ( purpose )
*       tau decay simulation using tauola 2.4
*
*     ( input ) 
*       IP       :  15 - TAU-
*                  -15 - TAU+
*       PTAU(3)  : MOMENTUM OF TAU ( MeV/c )
*
*     ( output ) 
*       COMMON  /HEPEVT/
*
*     ( creation date and author )
*       1995.10.25 ; First version by K.Kaneyuki
*
*     ( comment )
*       test version.
*       tau have only Pz.  
*      
***********************************************************************
#include "nework.h"
#include "taucom.h"
#include "tacard.h"
#include "tahepevt.h"
#include "taumc_com.h"
      
      COMMON /PTAU/ PTAUGV
      REAL PTAU(3)
      REAL POL(4),DUMMY(4),AQF(4)
      DATA ICALL/0/
      
C     VARIABLES FOR TAUPOL     
      INTEGER          FLAG, TAU_DEBUG, INDEX
      DOUBLE PRECISION SIGN, ENU, ETAU
      DOUBLE PRECISION THETA, COSTHETA, SINTHETA, PX, PY
      DOUBLE PRECISION SIGMA, SX, SY, SZ 

      REAL NUCLEON_MASS, NUCLEON_MOMENTUM, NUCLEON_ENERGY
      REAL NUCLEON_FOURVECT(4), NU_FOURVECT(4), TAU_FOURVECT(4)
      REAL BOOSTED_TAU(4), BOOSTED_NU(4)
      REAL DELTAX, DELTAY, DELTATHETA, POLTHETA, ARROWLENGTH

C     CONTROL THE OUTPUT      
      TAU_DEBUG = 1

      IF (ABS(IP).NE.15) RETURN
C
C --  INITIALIZATION
C      
      IF (ICALL.EQ.0) THEN
         ICALL=1
         CALL TAINICOM
         CALL DEXAY(-1,DUMMY)
      ENDIF

C --  CALCULATE MOMENTUM OF TAU LEPTON

      PTAUGV = SQRT(PTAU(1)**2+PTAU(2)**2+PTAU(3)**2)*1.0E-3
C
C --  SET POLARITY
C      

C --  IF ITPOL = 1 USE VALUE SET IN CARD FILE ELSE CALCULATE

      IF (ITAPOL.EQ.1) THEN
         DO 1 I=1,3
            IF (IP.GE.0) THEN
               POL(I)=-TAPOL(I)
            ELSE
               POL(I)= TAPOL(I)
            ENDIF
 1       CONTINUE
      ELSE

C     CALCULATE USING TAUPOL
         
C --  DO A LORENTZ TRANSFORM INTO THE REST FRAME OF THE STRUCK NUCLEON

C     CALCULATE ENERGY AND MOMENTUM OF NUCLEON

         CALL MCMASS(IPNE(2), NUCLEON_MASS)

         NUCLEON_MOMENTUM = SQRT( PNE(1,2)**2.0 + PNE(2,2)**2.0
     $        + PNE(3,2)**2.0)*1000.0

         NUCLEON_ENERGY = SQRT(NUCLEON_MOMENTUM**2.0
     $        + NUCLEON_MASS**2.0)

         IF (TAU_DEBUG.EQ.1) THEN
            WRITE(*,*) "NUCLEON CODE: ", IPNE(2), " MASS ", NUCLEON_MASS
            WRITE(*,*) "NUCLEON MOMENTUM: ", PNE(1,2), PNE(2,2),
     $           PNE(3,2), NUCLEON_MOMENTUM
            WRITE(*,*) "NUCLEON ENERGY: ", NUCLEON_ENERGY
         ENDIF

C     CALCULATE ENERGY AND MOMENTUM OF NEUTRINO AND LEPTON

C     NEUTRINO
         IF (IP.GE.0) SIGN =  1.0D0
C     ANTI NEUTRINO
         IF (IP.LT.0) SIGN = -1.0D0

         ENU   = SQRT( PNE(1,1)**2.0 + PNE(2,1)**2.0 +PNE(3,1)**2.0 )
         ETAU  = SQRT ( PTAUGV**2.0 + AMTAU**2.0)

         COSTHETA = (PNE(1,1)*PTAU(1)+PNE(2,1)*PTAU(2)+PNE(3,1)*PTAU(3))
     $        /( ENU*(PTAUGV*1000.0) )

         IF ( (COSTHETA.GT.1.0) .AND. (COSTHETA-1.0).LT. .0001) THEN
            COSTHETA = 1.0
         ENDIF

         THETA = (180.0/3.1415927)*ACOS( COSTHETA )

         IF (TAU_DEBUG.EQ.1) THEN
            WRITE(*,*) "TAU SIGN: ", SIGN
            WRITE(*,*) "TAU ENERGY (NU, LEPTON): ", ENU, ETAU
            WRITE(*,*) "TAU THETA: ", THETA
         ENDIF

C     PREPARE NUCLEON 4 VECTOR
         
         NUCLEON_FOURVECT(1) = PNE(1,2)*1000.0
         NUCLEON_FOURVECT(2) = PNE(2,2)*1000.0
         NUCLEON_FOURVECT(3) = PNE(3,2)*1000.0
         NUCLEON_FOURVECT(4) = NUCLEON_ENERGY

C     PREPARE TAU LEPTON 4 VECTOR

         TAU_FOURVECT(1) = PTAU(1)
         TAU_FOURVECT(2) = PTAU(2)
         TAU_FOURVECT(3) = PTAU(3)
         TAU_FOURVECT(4) = ETAU*1000.0

C     PREPARE NEUTRINO 4 VECTOR

         NU_FOURVECT(1) = PNE(1,1)*1000.0
         NU_FOURVECT(2) = PNE(2,1)*1000.0
         NU_FOURVECT(3) = PNE(3,1)*1000.0
         NU_FOURVECT(4) = ENU*1000.0

C     BOOST NU AND LEPTON TO THE NUCLEON REST FRAME
         
         IF ( NUCLEON_FOURVECT(1).EQ.0.0 .AND.
     $        NUCLEON_FOURVECT(2).EQ.0.0 .AND.
     $        NUCLEON_FOURVECT(3).EQ.0.0 ) THEN

            DO INDEX=1,4
               BOOSTED_TAU(INDEX) = TAU_FOURVECT(INDEX)
               BOOSTED_NU(INDEX)  = NU_FOURVECT(INDEX)
            END DO

         ELSE
            CALL LOREN4( NUCLEON_FOURVECT, TAU_FOURVECT, BOOSTED_TAU )
            CALL LOREN4( NUCLEON_FOURVECT, NU_FOURVECT,  BOOSTED_NU  )
         ENDIF

         IF (TAU_DEBUG.EQ.1) THEN
            WRITE(*,*) "LORENTZ VECTORS"
            WRITE(*,*) "NUCLEON:   ", (NUCLEON_FOURVECT(INDEX),INDEX=1,4)
            WRITE(*,*) "TAU LAB:   ", (TAU_FOURVECT(INDEX),INDEX=1,4)
            WRITE(*,*) "NU  LAB:   ", (NU_FOURVECT(INDEX),INDEX=1,4)
            WRITE(*,*) "TAU BOOST: ", (BOOSTED_TAU(INDEX),INDEX=1,4)
            WRITE(*,*) "NU  BOOST: ", (BOOSTED_NU(INDEX),INDEX=1,4)
         ENDIF

C     NOW CALL TAUPOL

         ENU  = BOOSTED_NU(4)
         ETAU = BOOSTED_TAU(4) 

C     ADD POTENTIAL DIFFERENCE
         ETAU = ETAU + 27

         COSTHETA = ( BOOSTED_NU(1)*BOOSTED_TAU(1) +
     $        BOOSTED_NU(2)*BOOSTED_TAU(2) +
     $        BOOSTED_NU(3)*BOOSTED_TAU(3) )/
     $        (ENU*SQRT( BOOSTED_TAU(1)**2.0 + BOOSTED_TAU(2)**2.0 +
     $        BOOSTED_TAU(3)**2.0) )

         IF ( (COSTHETA.GT.1.0) .AND. (COSTHETA-1.0).LT. .0001) THEN
            COSTHETA = 1
         ENDIF

         THETA = (180.0/3.1415927)*ACOS(COSTHETA)

         IF (TAU_DEBUG.EQ.1) THEN
            WRITE(*,*) "BOOSTED ENU:   ", ENU
            WRITE(*,*) "BOOSTED ETAU:  ", ETAU
            WRITE(*,*) "BOOSTED THETA: ", THETA
         ENDIF

         CALL TAUPOL(SIGN, ENU/1000.0, ETAU/1000.0, THETA, SIGMA,
     $        SX, SY, SZ, FLAG, MODENE)

         POL(1) = 2.0*SX
         POL(2) = 2.0*SY
         POL(3) = 2.0*SZ

#define WRITE_PAW_ARROWS
#ifdef WRITE_PAW_ARROWS

C        WRITE OUT TABLES FOR MAKING PLOTS TO CHECK CODE

         SINTHETA = SQRT(1.0 - COSTHETA**2)

         DELTAX      = (POL(3)*COSTHETA - POL(1)/5.0*SINTHETA)
         DELTAY      = (POL(3)*SINTHETA + POL(1)/5.0*COSTHETA)
         DELTATHETA  = ATAN2(DELTAY,DELTAX)
         ARROWLENGTH = 0.75
         POLTHETA     = ATAN2(POL(1),POL(3))

         PX = VMOD(BOOSTED_TAU,3)/1000.0*COSTHETA
         PY = VMOD(BOOSTED_TAU,3)/1000.0*SINTHETA

         WRITE(*,*) "ARROW: ", PX, PY, POL(1), POL(3)
         WRITE(*,*) "ANGLE OF POL: ", POLTHETA*180./3.14
         WRITE(*,*) "ANGLE OF TAU: ", THETA
         WRITE(*,*) "ANGLE OF ROTATED POL: ", DELTATHETA*180./3.14

         IF ( ABS(MODENE).NE.1) THEN
            IF ( SIGN.EQ.1) THEN
               write(64,*) "ARROW ", 
     $              PX, PX + DELTAX*ARROWLENGTH,
     $              PY, PY + DELTAY*ARROWLENGTH, " 0.2"
            ELSEIF ( SIGN.EQ.-1) THEN
               write(65,*) "ARROW ",
     $              PX, PX + DELTAX*ARROWLENGTH,
     $              PY, PY + DELTAY*ARROWLENGTH, " 0.2"
            ENDIF
         ENDIF
#endif

      ENDIF
C
C --  CLEAR /HEPEVT/
C      
      CALL RESLUX
C
C --  DETERMINE MOMENTA
C
      DO 10 I=1,3
         AQF(I)=0.0
 10   CONTINUE
      AQF(4)=AMTAU
      
      IF (IP.GE.0) CALL TRALO4(2,AQF,AQF,AM)
      IF (IP.LT.0) CALL TRALO4(1,AQF,AQF,AM)
C
C --  FILL TAU TO /HEPEVT
C      
      CALL FILHEP(1,1,IP,1,2,0,0,AQF,AMTAU,.TRUE.)
C
C --  EXECUTE KOLAL-Z      
C      
      IF (IP.GE.0) CALL DEXAY(2,POL)
      IF (IP.LT.0) CALL DEXAY(1,POL)
C
C --  ROTATE DIRECTION
C

C
C --  SAVE INFORMATION ABOUT THE TAU DECAY TO THE SKTAU COMMON BLOCK
C --  MODE IS SET IN TAUOLA
C 
      TAUMOM    = PTAUGV
      TAUDIR(1) = PTAU(1)/(PTAUGV*1000.0)
      TAUDIR(2) = PTAU(2)/(PTAUGV*1000.0)
      TAUDIR(3) = PTAU(3)/(PTAUGV*1000.0)
      POLTAU(1) = POL(1)
      POLTAU(2) = POL(2)
      POLTAU(3) = POL(3)
      
      IF (TAU_DEBUG.EQ.1) THEN
         WRITE(*,*) "TAU MOMENTUM IS: ",   TAUMOM
         WRITE(*,*) "TAU DECAY MODE IS: ", TAUDCY
         WRITE(*,*) "TAU POL IS: ", POLTAU(1), " ", POLTAU(2),
     $        " ", POLTAU(3), VMOD(POLTAU,3)
      ENDIF
C     
C --  FILL /MCWORK/
C

C     
C --  FILL /LUDATA/ ( IF NEED )
C      
C     CALL LUHEPC(2)

C
C --  DUMP /LUDATA/ ( IF NEED )
C      
C     CALL LULIST(2)

      RETURN
      END
      
