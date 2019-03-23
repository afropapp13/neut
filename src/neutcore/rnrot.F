************************************************************************
*     -------------------------------
      SUBROUTINE RNROT(COST,DIR,GDIR)
*     -------------------------------
*
*     (Purpose)
*       Give random direction with a given COSINE
*       Determin Cherenkov photon direction
*     
*     (Input)
*       COST    : COSINE OF THE ANGLE BETWEEN DIR & GDIR      
*       DIR(3)  : DIRECTION OF INPUT CHARGED PARTICLE
*     (Output)
*       GDIR(3) : DIRCTION OF OUTPUT GAMMA      
*
*     (Creation Date and Author)
*       1980.08.26 ; K.Arisaka
*       1995.02.03 ; K. KANEYUKI  FOR S.K.
*                    RANROT -> RNROT
*     
************************************************************************

      REAL DIR(3),GDIR(3)

      if (abs(cost).eq.1.) then
         sint=0.
      else
         SINT=SQRT(1.-COST**2)
      endif

      CALL RNAZI(C,S)
      X1=SINT*C
      Y1=SINT*S
      Z1=COST

      IF (ABS(DIR(1)).LT.1.E-2.AND.ABS(DIR(2)).LT.1.E-2) THEN
         IF (DIR(3).LT.0.) THEN
            GDIR(1)=-X1
            GDIR(2)=-Y1
            GDIR(3)=-Z1
         ELSE
            GDIR(1)=X1
            GDIR(2)=Y1
            GDIR(3)=Z1
         ENDIF
         RETURN
      ENDIF

      RDIR=DIR(1)**2+DIR(2)**2+DIR(3)**2
      RDIR=SQRT(RDIR)
      COSA=DIR(3)/RDIR
      if (abs(cosa).eq.1.) then
         sina = 0.
      else
         SINA=SQRT(1.-COSA**2)
      endif

      IF (SINA.LT.1.E-7) THEN
         IF (DIR(3).LT.0.) THEN
            GDIR(1)=-X1
            GDIR(2)=-Y1
            GDIR(3)=-Z1
         ELSE
            GDIR(1)=X1
            GDIR(2)=Y1
            GDIR(3)=Z1
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
