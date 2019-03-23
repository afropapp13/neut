*     -----------------------------------------------------
      SUBROUTINE EFPARQ(IPART,Q)      
*     -----------------------------------------------------
*     
*     ( purpose )
C       CHARGE OF THE PARTICLE
*       
*     ( input ) 
C       IPART:  PARTICLE # OF GEANT ASIGNMENT
*
*     ( output ) 
C         Q  :  CHARGE OF THE PARTICLE( FOR EXAMPLE ELECTRON= -1
C     PROTON = +1     )
*
C     
C     5-OCT-87   PRD. BY SATO,N
C     
C**********************************************************************

      DIMENSION CHAR(48)
      DATA CHAR/  0.,  1., -1.,  0.,  1., -1.,  0.,  1., -1.,  0.,
     &     1., -1.,  0.,  1., -1.,  0.,  0.,  0.,  1.,  0.,
     &     -1.,  0., -1., -1.,  0.,  0., -1.,  0.,  1.,  0.,
     &     1.,  1.,  1.,  1.,  1., -1.,  0.,  0.,  1., -1.,
     &     1.,  1., -1.,  0.,  1.,  1.,  2.,  0./
C------
      Q=CHAR(IPART)
      RETURN
      END
