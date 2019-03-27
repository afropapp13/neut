C************************************************
C  $Id: cone.F,v 1.1 2007-01-20 07:13:06 skrep Exp $
C
C  This is the kinematics used by IMB (west cost).
C
C  $Log: not supported by cvs2svn $
C  Revision 1.2  1997/06/02 22:08:02  mcgrew
C  Changes to sync skrep with sksuny.
C
C  Revision 1.1  1996/10/06 01:11:48  fgoebel
C  This is basically the IMB (west) kinematics program. A few
C  things have been changed to make it compatible with the SuperK
C  analysis software.
C
C************************************************

C*******************************************
C  @(#)cone.f   1.1 modified on 12/7/92
       SUBROUTINE CONE(SAa,SB,CTH,PHI)

C  ASSIGNS TO SB(3) A VECTOR FORMING AN ANGLE ACOS(CTH) WITH SA(3)
C  IF PHI<0 THEN A RANDOM PHI IS TAKEN (FROM CONICAL SURFACE)
C
       DIMENSION SA(4),saa(3),SB(3),D(3,3),RK(3)
        r=0.
        do i=1,3
                r=r+saa(i)**2
        end do
        r=sqrt(r)
        do i=1,3
                sa(i)=saa(i)/r
        end do

       CALL RANVE(1.,RK,CTH,PHI)
       R=SQRT(SA(1)**2+SA(2)**2)

        IF(R.Eq.0.)then
                SB(1)=RK(1)
                sb(2)=rk(2)
                sb(3)=rk(3)*sa(3)
        else
                D(1,1)=SA(1)*SA(3)/R
                D(1,2)=SA(2)*SA(3)/R
                D(1,3)=-R
                D(2,1)=-SA(2)/R
                D(2,2)=SA(1)/R
                D(2,3)=0.
                D(3,1)=SA(1)
                D(3,2)=SA(2)
                D(3,3)=SA(3)
                DO 1 I=1,3
                        SB(I)=0.
                        DO 1 J=1,3
 1                              SB(I)=SB(I)+D(J,I)*RK(J)
        end if

 10      RETURN
       END
