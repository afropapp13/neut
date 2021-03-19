***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE MKDCRS(IMODE, E, Q2, W, COSTHETA, PHI, XMLEP, DIFCRS)
*     ---------------------------------------------------
*     ( purpose )
*       calculate Differential Cross-section d^2(sig)/dwdQ2d(theta)d(phi) using Minoo Kabirnezhad's single pion model.
*
*       Code largely inspired by rsdcrs.F to make easy to follow
*       Calls different C functions for the different interaction modes
*
*     ( input )
*       imode : Interaction mode
*               imode > 10  mans anti-neutrino
*       e     : Incident neutrino energy in lab. frame.
*               (Lab. frame is nucleon rest frame.)
*       q2    : sqaured momentum transfer.(minus sign)
*       w     : Rest mass of products(pion and nucreon).
*       costheta: cos theta as defined in the Adler frame
*       phi   : as defined in the Adler frame
*       xmlep : mass of lepton
*
*     ( output )
*       difcrs: differential cross-section.
*               If error occurs in this subroutine,difcrs
*               returns -9999.999
*
*     ( creation date and author )
*       2016 Oct: Clarence Wret and Minoo Kabirnezhad
*                 Final implementation in about May 2017
*
***********************************************************************

      IMPLICIT NONE

C     Define the inputs
      INTEGER IMODE
      REAL E, Q2, W, COSTHETA, PHI
      REAL diff_xsec, difcrs

C     Use the neutmodel and Rein Sehgal constants
#include "rscons.h"
#include "neutmodel.h"

      REAL XMLEP

C     Checkers for valid kinematic
      REAL XNU 
      Real W2
      REAL K_0, Q_0, K_2, K_1, K_2_ISO

C-------------------------------------------------------------------------

      diff_xsec = 0.

C     Check the W range
      IF (W .LE. XMPI+XMN) THEN
        W = XMPI+XMN+1.E-5
      ENDIF
      W2 = W**2

c     KINEMATICAL QUANTITIES
c     Check that the Q2 is allowed
      IF ( ((2 * XMN * E) * (2 * XMN * E + XMN2 - W2) /
     $     (2 * XMN * E + XMN2) + Q2) .LT. -1.0E-5) THEN
c       write(*,*) "Q2 range"
c       write (*,*) "calcq2=", ( ((2 * XMN * E) * (2*XMN*E+XMN2-W2)
c    $  / (2 * XMN * E + XMN2)))
c       write (*,*) "q2=", q2
        GOTO 999
      ENDIF

      XNU = (W2 - XMN2 - Q2) / (2 * XMN)
      IF(XNU.GT.E) THEN
c       write(*,*) "XNU > E" 
        GOTO 999
      ENDIF

      Q_0 = (W2 - XMN2 + XMPI2)/(2.*W)
      if (Q_0*Q_0 - XMPI2 .LT. 0) THEN
c       write(*,*) "Q_02 < mpi2"
        GOTO 999
      endif

      K_0 = (W2 - XMN2 + Q2)/(2.*W)
      if (K_0*K_0 - Q2 .LT. 0) then
c       write(*,*) "K_02 < Q2"
        GOTO 999
      endif

      K_2 = (XMN2 + 2.*XMN*E - W2 -XMLEP*XMLEP)/(2.*W) ;
      K_1 =  k_2 + (W2 - XMN*XMN + Q2)/(2.*W) ;
      if (k_2*k_2-xmlep*xmlep .lt. 0) then
        goto 999
      endif
      K_2_ISO = sqrt(k_2*k_2 - XMLEP*XMLEP);
      if (abs(2.*k_1*k_2  + Q2 - XMLEP*XMLEP)/( 2.*k_1*k_2_iso) .gt. 1) then
c       write(*,*) "COSTHETA NEGATIVE"
        GOTO 999
      endif

C     SELECT Channels for neutrino and anti-neutrino
C     Neutrino
      IF(IMODE.EQ.1)  THEN 
        GOTO 1001  ! CC 'P - PI+'
      ELSE IF(IMODE.EQ.2) THEN  
        GOTO 1002  ! CC 'P - PI0'
      ELSE IF(IMODE.EQ.3) THEN  
        GOTO 1003  ! CC 'N - PI+'
      ELSE IF(IMODE.EQ.4) THEN  
        GOTO 1004  ! NC 'P - PI0'
      ELSE IF(IMODE.EQ.5) THEN  
        GOTO 1005  ! NC 'N - PI+'
      ELSE IF(IMODE.EQ.6) THEN  
        GOTO 1006  ! NC 'N - PI0'
      ELSE IF(IMODE.EQ.7) THEN  
        GOTO 1007  ! NC 'P - PI-'
C     anti-Neutrino
      ELSE IF(IMODE.EQ.11) THEN  
        GOTO 1011 ! CC 'N - PI-'
      ELSE IF(IMODE.EQ.12) THEN  
        GOTO 1012 ! CC 'N - PI0'
      ELSE IF(IMODE.EQ.13) THEN  
        GOTO 1013 ! CC 'P - PI-' 
      ELSE IF(IMODE.EQ.14) THEN  
        GOTO 1014 ! NC 'P - PI0'
      ELSE IF(IMODE.EQ.15) THEN  
        GOTO 1015 ! NC 'N - PI+' 
      ELSE IF(IMODE.EQ.16) THEN  
        GOTO 1016 ! NC 'N - PI0'
      ELSE IF(IMODE.EQ.17) THEN  
        GOTO 1017 ! NC 'P - PI-' 
      ELSE 
        write(*,*) "! WRONG IMODE FOR MK MODEL !"
        diff_xsec = -1.00
        GOTO 69
      ENDIF

 1001 CONTINUE
      CALL Xsec_IM1(diff_xsec, W, -Q2, E, COSTHETA, PHI, XMLEP)
      GOTO 69

 1002 CONTINUE
      CALL Xsec_IM2(diff_xsec, W, -Q2, E, COSTHETA, PHI, XMLEP)
      GOTO 69

 1003 CONTINUE
      CALL Xsec_IM3(diff_xsec, W, -Q2, E, COSTHETA, PHI, XMLEP)
      GOTO 69

C     NC interactions don't need lepton mass passed
 1004 CONTINUE
      CALL Xsec_IM4(diff_xsec, W, -Q2, E, COSTHETA)
      GOTO 69

 1005 CONTINUE
      CALL Xsec_IM5(diff_xsec, W, -Q2, E, COSTHETA)
      GOTO 69

 1006 CONTINUE
      CALL Xsec_IM6(diff_xsec, W, -Q2, E, COSTHETA)
      GOTO 69

 1007 CONTINUE
      CALL Xsec_IM7(diff_xsec, W, -Q2, E, COSTHETA)
      GOTO 69

 1011 CONTINUE
      CALL Xsec_IM11(diff_xsec, W, -Q2, E, COSTHETA, PHI, XMLEP)
      GOTO 69

 1012 CONTINUE
      CALL Xsec_IM12(diff_xsec, W, -Q2, E, COSTHETA, PHI, XMLEP)
      GOTO 69

 1013 CONTINUE
      CALL Xsec_IM13(diff_xsec, W, -Q2, E, COSTHETA, PHI, XMLEP)
      GOTO 69

C     NC interactions don't need lepton mass passed
 1014 CONTINUE
      CALL Xsec_IM14(diff_xsec, W, -Q2, E, COSTHETA)
      GOTO 69

 1015 CONTINUE
      CALL Xsec_IM15(diff_xsec, W, -Q2, E, COSTHETA)
      GOTO 69

 1016 CONTINUE
      CALL Xsec_IM16(diff_xsec, W, -Q2, E, COSTHETA)
      GOTO 69

 1017 CONTINUE
      CALL Xsec_IM17(diff_xsec, W, -Q2, E, COSTHETA)
      GOTO 69

 999  CONTINUE
C     write(*,*) "Found 999 flag in mkdcrs.F"
      diff_xsec = 0
      GOTO 69

C     Go here when we've calculated a cross-section
 69   CONTINUE
      difcrs = diff_xsec

      IF (ISNAN(difcrs)) THEN
        DIFCRS = 0
        WRITE(*,*) "FOUND NAN FOR CROSS-SECTION IN KABIRNEZHAD 1PI:"
        WRITE(*,*) "SETTING TO ZERO"
        WRITE(*,*) "MODE = ", IMODE
        WRITE(*,*) "E = ", E
        WRITE(*,*) "W = ", W
        WRITE(*,*) "Q2 = ", Q2
        WRITE(*,*) "COSTHETA = ", COSTHETA
        WRITE(*,*) "PHI = ", PHI
      ENDIF

      END
