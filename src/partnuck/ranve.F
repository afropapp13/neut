C************************************************
C  $Id: ranve.F,v 1.1.1.1 2016/11/04 20:31:14 mdunkman Exp $
C
C  This is the kinematics used by IMB (west cost).
C
C  $Log: ranve.F,v $
C  Revision 1.1.1.1  2016/11/04 20:31:14  mdunkman
C  Import of NEUT 5.3.3 v1r27 for use with T2KReWeight v1r27 (Winter OA 2017)
C
C  Revision 1.2  1997/06/02 22:08:28  mcgrew
C  Changes to sync skrep with sksuny.
C
C  Revision 1.1  1996/10/06 01:11:46  fgoebel
C  This is basically the IMB (west) kinematics program. A few
C  things have been changed to make it compatible with the SuperK
C  analysis software.
C
C************************************************

C*****************************************************
C     @(#)ranve.f	1.1 modified on 12/7/92
       SUBROUTINE RANVE(VR,V,Q1,Q2)
C
C  ASSIGNS TO V(3) A RANDOM VECTOR OF LENGTH VR
	DIMENSION V(3)
       CT=Q1
       PHI=Q2
       IF(abs(CT).GT.1.0)CT=2.*ranf()-1.
       IF(PHI.LT.0..or.phi.gt.6.28318)PHI=6.28318*ranf()
       ST=SQRT(1.-CT**2)
       V(1)=VR*ST*COS(PHI)
       V(2)=VR*ST*SIN(PHI)
       V(3)=VR*CT
       RETURN
       END
