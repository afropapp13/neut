C************************************************
C  $Id: prbin.F,v 1.1 2007-01-20 07:13:08 skrep Exp $
C
C  This is the kinematics used by IMB (west cost).
C
C  $Log: not supported by cvs2svn $
C  Revision 1.3  1997/06/02 22:08:26  mcgrew
C  Changes to sync skrep with sksuny.
C
C  Revision 1.2  1997/01/29 21:06:46  mcgrew
C  Improve the histograming and fix some overly crude binning effects.
C
C  Revision 1.1  1996/10/06 01:11:37  fgoebel
C  This is basically the IMB (west) kinematics program. A few
C  things have been changed to make it compatible with the SuperK
C  analysis software.
C
C************************************************

C*****************************************************
C  @(#)prbin.f	1.1 modified on 12/7/92
       FUNCTION	PRBIN(YI,N,XMIN,XSTEP,YR)
       DIMENSION YI(N+1)
       MIN=1
       MAX=N+1
       GO TO 3
    1  IF((MAX-MIN).EQ.1) GO TO	5
       IF(YI(M)-YR) 2,6,4
    2  MIN=M
    3  M=(MAX+MIN)/2
       GO TO 1
    4  MAX=M
       GO TO 3
    5  M=MIN
       FR=(YR-YI(MIN))/(YI(MAX)-YI(MIN))
       PRBIN=XMIN+XSTEP*(FLOAT(MIN)+FR-1.0)
       RETURN
    6  FR=0.0
       PRBIN=XMIN+XSTEP*(FLOAT(M)-1.0)
       RETURN
       END
