*     -----------------------------------------------------
      SUBROUTINE EFTRANS(TMX,XX,XIN,XOUT,MODE)
*     -----------------------------------------------------
*
*     ( purpose )
C         TRANSFORMATION  ----- ROTATION + TRASLATION --------*
*
*     ( input )
C         TMX(3,3) : ROTATION MATRIX
C         XX(3)    : TRANSLATION
C         XIN(3)   : INPUT VECTOR
C         MODE     : 1: ROT + TRA
C                    0: ROT
*     ( output )
C         XOUT(3)  : OUTPUT VECTOR
C
C*********************************************************************
      DIMENSION TMX(3,3),XX(3),XIN(3),XOUT(3),XTR(3)
      IF(MODE.EQ.1) THEN
         XTR(1)=XX(1)
         XTR(2)=XX(2)
         XTR(3)=XX(3)
      ELSE
         IF(MODE.EQ.0) THEN
            XTR(1)=0.
            XTR(2)=0.
            XTR(3)=0.
         ELSE
            PRINT *,' +++++ MODE SELECT ERROR AT SUB-TRANS ++++'
         END IF
      END IF
      XOUT(1)=TMX(1,1)*XIN(1)+TMX(1,2)*XIN(2)+TMX(1,3)*XIN(3)
     &     +XTR(1)
      XOUT(2)=TMX(2,1)*XIN(1)+TMX(2,2)*XIN(2)+TMX(2,3)*XIN(3)
     &     +XTR(2)
      XOUT(3)=TMX(3,1)*XIN(1)+TMX(3,2)*XIN(2)+TMX(3,3)*XIN(3)
     &     +XTR(3)
      RETURN
      END
