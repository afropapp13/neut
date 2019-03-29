C************************************************
C  $Id: lloren.F,v 1.1 2007-01-20 07:13:06 skrep Exp $
C
C  This is the kinematics used by IMB (west cost).
C
C  $Log: not supported by cvs2svn $
C  Revision 1.3  1997/11/22 09:53:59  bviren
C  Many changes:
C
C   - Fix proton.F to not give bogus kinematics for the 2nd produced pion.
C
C   - Added pdkin.c to use for (currently, only) p->e+pi0 kinematics
C  generation. This still uses partnuc for nuclear corrections.
C
C   - Added `fermib', which is just `fermid' but with out the determination
C  of nuclear/free protonness. Also values are in MeV and correct P(4)
C  is used (uses the effective mass). See fermit for simple test code.
C
C   - Added various C-wrappers and prototypes to call IMBWest routines.
C
C   - Added code to partnuc to track which interactions the pion underwent.
C
C   - Added program `picount' which will read in a track file and count
C  some statistics.
C
C   - Other things which I probably forgot.
C
C  Revision 1.2  1997/06/02 22:08:18  mcgrew
C  Changes to sync skrep with sksuny.
C
C  Revision 1.1  1996/10/06 01:11:36  fgoebel
C  This is basically the IMB (west) kinematics program. A few
C  things have been changed to make it compatible with the SuperK
C  analysis software.
C
C************************************************

C****************************************************
C   @(#)lloren.f        1.1 modified on 12/7/92
      SUBROUTINE lLOREN(PPRIME,EPRIME,P,E,BETA)
      DIMENSION PPRIME(3),P(3),BETA(3)
      BSQ=0.
      BDP=0.
      DO I=1,3
         BSQ=BSQ+BETA(I)**2
         BDP=BDP+BETA(I)*P(I)
      end do
      IF(BSQ) 30,30,15
 15   GAM=SQRT(1./(1.-BSQ))
      SHIFT=(GAM-1.)*BDP/BSQ -GAM*E
      DO I=1,3
         PPRIME(I)=P(I)+SHIFT*BETA(I)
      end do
      EPRIME=GAM*(E-BDP)
      RETURN
 30   DO I=1,3
         PPRIME(I)=P(I)
      end do
      EPRIME=E
      RETURN
      END
