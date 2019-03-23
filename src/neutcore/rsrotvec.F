***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE RSROTVEC(PINVEC,POUTVEC,COSP,SINP,COST,SINT)
*     ---------------------------------------------------
*
*     ( purpose )
*       rotate vector according to given value
*
*     ( input )
*       pinvec(3)   : input vector
*       cosp        : cos(phi)
*       sinp        : sin(phi)
*       cost        : cos(theta)
*       sint        : sin(theta)
*     
*     ( output )
*       poutvec(3)  : output vector
*       
*
*     ( creation date and author )
*       1993.Dec. ; for Super-Kamioka by Y.Hayato
*                   (Based on Rein and Sehgal's code)
*
*     ( comment )
*
***********************************************************************

      IMPLICIT NONE

      REAL PINVEC(3),POUTVEC(3)
      REAL COSP,SINP,COST,SINT

      REAL TMPVEC(3)

      TMPVEC(3)=COSP*PINVEC(3)-SINP*PINVEC(1)
      TMPVEC(1)=SINP*PINVEC(3)+COSP*PINVEC(1)
      TMPVEC(2)=PINVEC(2)

      POUTVEC(1)=COST*TMPVEC(1)-SINT*TMPVEC(2)
      POUTVEC(2)=SINT*TMPVEC(1)+COST*TMPVEC(2)
      POUTVEC(3)=TMPVEC(3)

      RETURN
      END
      
