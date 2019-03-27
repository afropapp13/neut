          SUBROUTINE nrrODIS(RMIN,DR,NBIN,PROBM,ROT)
          implicit none
          real*4    rmin,dr
          integer*4 nbin
          REAL*4   ROT(NBIN),PROBM(NBIN)

          integer*4 i,j
          real*4    r

          real*4   nrroxy
          external nrroxy

          R=RMIN-DR
          DO I=1,NBIN
                 R=R+DR
                 ROT(I)=nrROXY(R)
                 PROBM(I)=4*3.1416*((R+DR)**3-R**3)/3
          END DO
          J=NBIN-1
          DO I=1,J
                 PROBM(I)=PROBM(I)*(ROT(I)+ROT(I+1))/2
          END DO
          PROBM(NBIN)=PROBM(NBIN)*ROT(NBIN)
          RETURN
          END
