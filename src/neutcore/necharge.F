**********************************************************************
*     ---------------------
      FUNCTION NECHARGE(IP)
*     ---------------------
*     
*     ( purpose )
*       CHARGE OF PARTICLE
*
*     ( input ) 
*       IP      : PARICLE CODE
*
*     ( output ) 
*       NECHARGE : CHARGE OF PARTICLE
*
*     ( creation date and author )
*       198?.??.?? ; First version by M.Nakahata
*       1995.03.03 ; K.Kaneyuki for S.K.
*                    ICHARG -> NECHARGE
*
**********************************************************************
      PARAMETER(MAXPAR=13)
      INTEGER NEUTRAL(MAXPAR)
C                  NEU-E NEU-M NEU-T GAMMA     N 
      DATA NEUTRAL/   12,   14,   16,   22, 2112,
C                    PI0    K0    KS    KL   ETA
     &               111,  311,  310,  130,  221,
C                   RHO0 OMEGA   K*0
     &               113,  223,  313/

     
      NECHARGE=0
      IPABS=IABS(IP)
      DO 10 I=1,MAXPAR
         IF (IPABS.EQ.NEUTRAL(I)) RETURN
 10   CONTINUE
C
C --  LEPTON
C
      IF(IPABS.EQ.11 .OR. IPABS.EQ.13 .OR. IPABS.EQ.15) THEN
         IF (IP.GT.0) THEN
            NECHARGE=-1
         ELSE
            NECHARGE= 1
         ENDIF
C
C --  OTHERS
C
      ELSE
         IF (IP.GT.0) THEN
            NECHARGE= 1
         ELSE
            NECHARGE=-1
         ENDIF
      ENDIF

      RETURN
      END
