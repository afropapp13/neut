************************************************************************
*     -------------------------------------
      SUBROUTINE NECHADIR(COST,PHI,DIR,GDIR)
*     -------------------------------------
*
*     (Purpose)
*       CHANGE DIRECTION GIVEN BY COS(THETA) AND PHI
*
*     (Input)
*       COST   : COSINE OF THE ANGLE BETWEEN DIR & GDIR
*       PHI    : PHI    OF THE GDIR
*       DIR(3) : DIRECTION OF INITIAL VECTOR
*     
*     (Output)
*       GDIR   : DIRECTION OF FINAL   VECTOR
*
*     (Creation Date and Author)
*       1984.05.07 ; M.NAKAHATA
*       1995.02.08 ; K. KANEYUKI  FOR S.K.
*                    CHADIR -> NECHADIR
*     
************************************************************************
      REAL DIR(3),GDIR(3)

      if (abs(COST).ge.1.) then
         SINT=0.
      else
         SINT=SQRT(1.-COST**2)
      endif
      RDIR=DIR(1)**2+DIR(2)**2+DIR(3)**2
      RDIR=SQRT(RDIR)

      C=COS(PHI)
      S=SIN(PHI)
      X1=SINT*C*RDIR
      Y1=SINT*S*RDIR
      Z1=COST*RDIR

      IF(ABS(DIR(1)).LT.1.E-2.AND.ABS(DIR(2)).LT.1.E-2) THEN
         IF(DIR(3).GT.0.) THEN
            GDIR(1)=X1
            GDIR(2)=Y1
            GDIR(3)=Z1
         ELSE
            GDIR(1)=-X1
            GDIR(2)=-Y1
            GDIR(3)=-Z1
         ENDIF
         RETURN
      ENDIF

      COSA=DIR(3)/RDIR
      IF (abs(COSA).eq.1.) then
         SINA = 0.
      else
         SINA=SQRT(1.-COSA**2)
      endif

      IF(SINA.LT.1.E-7) THEN
         IF(DIR(3).GT.0.) THEN
            GDIR(1)=X1
            GDIR(2)=Y1
            GDIR(3)=Z1
         ELSE
            GDIR(1)=-X1
            GDIR(2)=-Y1
            GDIR(3)=-Z1
         ENDIF
         RETURN
      ENDIF
      
      COSB=DIR(1)/RDIR/SINA
      SINB=DIR(2)/RDIR/SINA

      X2=X1*COSA+Z1*SINA
      Y2=Y1
      Z2=-X1*SINA+Z1*COSA

      GDIR(1)=X2*COSB-Y2*SINB
      GDIR(2)=X2*SINB+Y2*COSB
      GDIR(3)=Z2

      RETURN
      END
