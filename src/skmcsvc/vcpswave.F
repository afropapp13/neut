***********************************************************************
*     -------------------------------
      SUBROUTINE VCPSWAVE(PFERM,WGHT)
*     -------------------------------
*     
*     ( purpose )
*       calc. P distribution of S-wave
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
*                    PSWAVE -> VCPSWV
*       1995. 2. 7 ; change name
*                    VCPSWV -> VCPSWAVE
*
*     ( comment )
*
***********************************************************************
      DIMENSION DIST(10)
      DATA DIST/4. , 3.8 , 3.4 , 2.6 , 1.9 , 1.2 , 0.75 , 0.4 , 0.2,
     &     0./
C     
      WGHT=0.
      IF(PFERM.GE.225.)RETURN
      I=IFIX(PFERM/25.)+1
      FRAC=(PFERM-25.*FLOAT(I-1))/25.
      WGHT=DIST(I)+(DIST(I+1)-DIST(I))*FRAC
      WGHT=WGHT*PFERM**2

      RETURN
      END
