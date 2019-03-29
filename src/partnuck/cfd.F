C************************************************
C  $Id: cfd.F,v 1.1 2007-01-20 07:13:06 skrep Exp $
C
C  This is the kinematics used by IMB (west cost).
C
C  $Log: not supported by cvs2svn $
C  Revision 1.2  1997/06/02 22:08:01  mcgrew
C  Changes to sync skrep with sksuny.
C
C  Revision 1.1  1996/10/06 01:11:46  fgoebel
C  This is basically the IMB (west) kinematics program. A few
C  things have been changed to make it compatible with the SuperK
C  analysis software.
C
C************************************************

#define VERSION "@(#)cfd.f      1.1 modified on 12/7/92"
C******************************
C  Convert the probablity distribution Y into a cumulative probablity
C     distribution YI.  Y does not need to be normalized.  YI will be
C     normalized correctly.
      SUBROUTINE CFD(YI,N,Y)
      INTEGER N
      REAL Y(N),YI(N+1)
      CF=0.0
      DO I=1,N
         CF=CF+Y(I)
         YI(I+1)=CF
      ENDDO
      YI(1)=0.0
      cf=1./cf
      DO I=2,n+1
         YI(I)=YI(I)*cf
      ENDDO
      RETURN
      END
