C   18/10/85 510251024  MEMBER NAME  SIGETA   (FORT)        FORTRAN77
C
      SUBROUTINE SIGETA(RMASS,PCMS,SIG)
C              RMASS-- INV. MASS (MEV)
C              PCMS -- M0M. IN C.M.S. (MEV)
C              SIG -- CROSS SECTION (MB)
C
C   * CONSIDER 3 RESONANCES N(1440) GAMMA=200  BR(ETA,N)=10%
C                           N(1535) GAMMA=150  BR(ETA,N)=35%
C                           N(1710) GAMMA=110  BR(ETA,N)=25%
C                   REF.) PARTICLE DATA BOOKLET APRIL 1984
C
C   * THE TOTAL CROSS SECTION IS ESTIMATED BY THE FOLLOWING EQ.(B-W),
C              PAI    GAMMA(ETA,N)*GAMMA(TOTAL-ETA,N)
C      SIGM= ------*---------------------------------------
C             K**2  (M(RESONANCE)-RMASS)**2 + 1/4*GAMMA**2
C                  REF.) D.A.SPARROW, UNIV. OF PENN PREPRINT (1984)
C
C        T.KAJITA   18-OCT-85
C        T.KAJITA   22-OCT-85  MODIFICATION
C        J.KAMEDA   08-DEC-97  MODIFICATION
C          consider only N(1540)
C        J.KAMEDA   20-FEB-98  MODIFICATION
C          add N(1650) 
C        
      REAL K
C
      
      PAI=3.14
      K=PCMS/197.
      sig1 = 0.     
      sig2 = 0.     
      sig  = 0.

C
C (1) -- N(1535)     GAMMA (WIDTH) IS MODIFIED ACCORDING TO
C                REF) MIYAZAWA AND TAKEDA, ELEMENTARY PARTICLE PHYSICS
      A=1.54
      FACT=(A*K)**3/(1+(A*K)**2)
      GAMTOT=150.*FACT
      GAMETA=GAMTOT*0.35
      GAMPAI=GAMTOT*0.65
      AMAS=1535.

      SIG1 = PAI/K**2 *  GAMETA*GAMPAI   /
     & ( (AMAS-RMASS)**2 + 0.25*GAMTOT**2)  *10.

C
C (2) -- N(1650)     GAMMA (WIDTH) IS MODIFIED ACCORDING TO
C                REF) MIYAZAWA AND TAKEDA, ELEMENTARY PARTICLE PHYSICS
      A=1.54
      FACT=(A*K)**3/(1+(A*K)**2)
      GAMTOT=150.*FACT

      GAMETA=GAMTOT*0.07
      GAMPAI=GAMTOT*0.80
      AMAS=1650.
      SIG2 = PAI/K**2 *  GAMETA*GAMPAI   /
     & ( (AMAS-RMASS)**2 + 0.25*GAMTOT**2)  *10.

      sig = sig1 + sig2
c      sig = sig1
      RETURN
      END





