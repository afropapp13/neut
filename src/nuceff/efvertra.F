*     -----------------------------------------------------
      SUBROUTINE EFVERTRA(PMOM,TMX)
*     -----------------------------------------------------
*
*     ( purpose )
C        FROM ( 0, 0, 0 ) TO (PX,PY,PZ) ROTATION MATRIX
*
*     ( input )
*       PMOM    : input matrix(momentum of particle)
*
*     ( output )
*       TMX     : roration matrix
*
C
C     5-OCT-87
C     PRD. BY SATO,N
C
C******************************************************************
      DIMENSION PMOM(3),TMX(3,3)
      DATA PI/3.1415/
      HPI=PI/2.
      TPI=PI*2.
C------
      CALL VZERO(TMX,9)
      ABSP=PMOM(1)**2+PMOM(2)**2+PMOM(3)**2
      ABSP=SQRT(ABSP)
      IF(ABSP.EQ.0.) RETURN
C--------------------
C     THETA( RAD UNIT )
C--------------------
      THETA=ACOS(PMOM(3)/ABSP)
C--------------------
C     PHI( RAD UNIT )
C-------------------
      CD=SQRT(PMOM(1)**2+PMOM(2)**2)
      IF(CD.LT.0.0001) THEN
         ABACK=0.
      ELSE
         X=PMOM(1)/CD
         Y=PMOM(2)/CD
         IF(X.GT.0.0001) ABACK=ATAN(Y/X)
         IF(X.LT.-0.0001) ABACK=ATAN(Y/X)+PI
         IF(ABS(X).LE.0.0001) ABACK=SIGN(HPI,Y)
         IF(ABACK.LT.0.0) ABACK=ABACK+TPI
      END IF
      PHI=ABACK
C-------------
C     MATRIX
C------------
      COSP=COS(PHI)
      SINP=SIN(PHI)
      COST=COS(THETA)
      SINT=SIN(THETA)
C---
      TMX(1,1)=COSP*COST
      TMX(1,2)=-SINP
      TMX(1,3)=COSP*SINT
      TMX(2,1)=SINP*COST
      TMX(2,2)=COSP
      TMX(2,3)=SINP*SINT
      TMX(3,1)=-SINT
      TMX(3,2)=0.
      TMX(3,3)=COST
C-----
      RETURN
      END
