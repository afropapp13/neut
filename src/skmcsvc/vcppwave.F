***********************************************************************
*     -------------------------------
      SUBROUTINE VCPPWAVE(PFERM,WGHT)
*     -------------------------------
*     
*     ( purpose )
*       calc. P distribution of P-wave
*
*     ( input ) 
*       PFERM : MOMENTUM ( MeV/c )
*
*     ( output ) 
*       WGHT  : weight
*
*     ( creation date and author )
*                  ; First version by M.Nakahata
*       1992. 9.18 ; Mod by K.Kaneyuki for S.K.
*                    PPWAVE -> VCPPWV
*       1995. 2. 7 ; change name
*                    VCPPWV -> VCPPWAVE      
*     ( comment )
*
***********************************************************************
      DIMENSION DIST(12)
      DATA DIST/0.8 , 1.6 , 3.4 , 5.9 , 5.5 , 4.9 , 3.8 , 2.5 , 1.5 ,
     &     0.8 , 0.4 , 0./
C     
      WGHT=0.
      IF(PFERM.GE.275.)RETURN
      I=IFIX(PFERM/25.)+1
      FRAC=(PFERM-25.*FLOAT(I-1))/25.
      WGHT=DIST(I)+(DIST(I+1)-DIST(I))*FRAC
      WGHT=WGHT*PFERM**2

      RETURN
      END
