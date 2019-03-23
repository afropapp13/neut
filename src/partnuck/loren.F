C************************************************
C  $Id: loren.F,v 1.1 2007-01-20 07:13:07 skrep Exp $
C
C  This is the kinematics used by IMB (west cost).
C
C  $Log: not supported by cvs2svn $
C  Revision 1.4  1997/10/28 15:04:34  bviren
C  Removed superfluous do labels
C
C  Revision 1.3  1997/10/28 14:47:19  bviren
C  Added comments, indents, and the code (loren()) has been checked by
C  Gajewski and I for correctness.
C
C  Revision 1.2  1997/06/02 22:08:19  mcgrew
C  Changes to sync skrep with sksuny.
C
C  Revision 1.1  1996/10/06 01:11:39  fgoebel
C  This is basically the IMB (west) kinematics program. A few
C  things have been changed to make it compatible with the SuperK
C  analysis software.
C
C************************************************
C  Definitions of arguments
C  C = 4 momentum of decaying particle in LAB system
C  P = 4 momentum of decay products in decaying particle's rest frame
C  Q = 4 momentum of decay products in LAB system
C  N = number of decay products.
C  K = -1, (if K=+1 then def of P and Q are switched)
C
C***********************************************
C   @(#)loren.f	1.1 modified on 12/7/92
      SUBROUTINE LOREN(C,P,Q,N,K)
      DIMENSION C(4),P(4,N),Q(4,N)
      S=0.
      DO I=1,3
         S=S+C(I)**2
      end do

      IF(S.EQ.0.)GO TO	100
      S=C(4)**2-S
      W=SQRT(S)
      Z=W+C(4)
      DO J=1,N
         D=C(4)*P(4,J)
         DO I=1,3
            D=D-K*C(I)*P(I,J)
         end do
         E=D/W
         A=(P(4,J)+E)/Z
         Q(4,J)=E
         DO I=1,3
            Q(I,J)=P(I,J)-A*K*C(I)
         end do
      end do

      GO TO 1000
 100  DO J=1,N
         DO I=1,4
            Q(I,J)=P(I,J)
         end do
      end do
      GO TO 1000
 1000 RETURN
      END




