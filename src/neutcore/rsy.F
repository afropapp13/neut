C     A file to store all the Spherical Harmonics needed for the
C     Rein-Sehgal single pion code.
C
C     Thought I'd bump them together here instead of bloating the
C     directory with rsy**.F


************************************************************************
*     ---------------------------------
      REAL FUNCTION RSY00(THETA,PHI)
*     ---------------------------------
C
C     (Purpose)
C        Calculate spherical harmonics function Y00
C
C     (Input)
C        Theta : theta
C        phi   : phi
C     (Output)
C        Y00
C
C     (Creation Date and Author)
C        199?.??.?? ; Y.Hayato
C
************************************************************************      
      REAL*4 THETA,PHI

      PARAMETER(PI = 3.14159265)

      RSY00 = SQRT(1./(4.*PI))

      RETURN

      END

************************************************************************
*     ---------------------------------
      REAL FUNCTION RSY10(THETA,PHI)
*     ---------------------------------
C
C     (Purpose)
C        Calculate spherical harmonics function Y10
C
C     (Input)
C        Theta : theta
C        phi   : phi
C     (Output)
C        Y10
C
C     (Creation Date and Author)
C        2016 Oct. C Wret: For Rein-Sehgal pion ejection update
C
************************************************************************      

      REAL*4 THETA,PHI

      PARAMETER(PI = 3.14159265)

      RSY10 = 0.5 * SQRT(3./PI) 
     &      * COS(THETA)

      RETURN

      END

************************************************************************
*     ---------------------------------
      COMPLEX FUNCTION RSY11(THETA,PHI)
*     ---------------------------------
C
C     (Purpose)
C        Calculate spherical harmonics function Y11
C
C     (Input)
C        Theta : theta
C        phi   : phi
C     (Output)
C        Y21
C
C     (Creation Date and Author)
C       2016.Oct C. Wret: Add in the imaginary parts for Rein-Sehgal
C                         pion ejection update
C
************************************************************************

      REAL*4 THETA,PHI

      PARAMETER(PI = 3.14159265)

      RSY11 = (-1./2.) * SQRT(3./(2.*PI)) 
     &  * SIN(THETA)
     &  * CMPLX(COS(PHI), SIN(PHI))

      RETURN

      END

************************************************************************
*     ---------------------------------
      REAL FUNCTION RSY20(THETA,PHI)
*     ---------------------------------
C
C     (Purpose)
C        Calculate spherical harmonics function Y20
C
C     (Input)
C        Theta : theta
C        phi   : phi
C     (Output)
C        Y20
C
C     (Creation Date and Author)
C        199?.??.?? ; Y.Hayato
C
************************************************************************      

      REAL*4 THETA,PHI

      PARAMETER(PI = 3.14159265)

      RSY20 = (1./4.) * SQRT(5./PI) 
     &      * (3.*COS(THETA)**2 - 1)

      RETURN

      END

************************************************************************
*     ---------------------------------
      COMPLEX FUNCTION RSY21(THETA,PHI)
*     ---------------------------------
C
C     (Purpose)
C        Calculate spherical harmonics function Y21
C
C     (Input)
C        Theta : theta
C        phi   : phi
C     (Output)
C        Y21
C
C     (Creation Date and Author)
C       199?.??.?? ; Y.Hayato
C       2016.Oct C. Wret: Add in the imaginary parts for Rein-Sehgal
C                         pion ejection update
C
************************************************************************

      REAL*4 THETA,PHI

      PARAMETER(PI = 3.14159265)

      RSY21 = (-1./2.) * SQRT(15./(2.*PI)) 
     &  * COS(THETA)*SIN(THETA)
     &  * CMPLX(COS(PHI), SIN(PHI))

      RETURN

      END

************************************************************************
*     ---------------------------------
      COMPLEX FUNCTION RSY22(THETA,PHI)
*     ---------------------------------
C
C     (Purpose)
C        Calculate spherical harmonics function Y22
C
C     (Input)
C        Theta : theta
C        phi   : phi
C
C     (Output)
C        Y22
C
C     (Creation Date and Author)
C       199?.??.?? ; Y.Hayato
C       2016.Oct C. Wret: Add in the imaginary parts for pion ejection
C                         update to Rein-Sehgal
C
************************************************************************

      REAL*4 THETA,PHI

      PARAMETER(PI = 3.14159265)

      RSY22 = (1./4.) * SQRT(15./(2.*PI))
     &      * SIN(THETA)**2
     &      * CMPLX(COS(2.*PHI), SIN(2.*PHI))

      RETURN

      END

************************************************************************
*     ---------------------------------
      REAL FUNCTION RSY30(THETA,PHI)
*     ---------------------------------
C
C     (Purpose)
C        Calculate spherical harmonics function Y30
C
C     (Input)
C        Theta : theta      
C        phi   : phi
C     (Output)    
C        Y30
C
C     (Creation Date and Author)
C        2016 Oct. C Wret: For Rein-Sehgal pion ejection update
C
************************************************************************      

      REAL*4 THETA,PHI

      PARAMETER(PI = 3.14159265)

      RSY30 = (1./4.) * SQRT(7./PI)
     &      * (5*COS(THETA)**3 - 3*COS(THETA))

      RETURN

      END

************************************************************************
*     ---------------------------------
      COMPLEX FUNCTION RSY31(THETA,PHI)
*     ---------------------------------
C
C     (Purpose)
C        Calculate spherical harmonics function Y31
C
C     (Input)
C        Theta : theta      
C        phi   : phi
C     (Output)    
C        Y31
C
C     (Creation Date and Author)
C        2016 Oct C. Wret: Needed for pion ejection update Rein-Sehgal
C
************************************************************************      

      REAL*4 THETA,PHI

      PARAMETER(PI = 3.14159265)

      RSY31 = (-1./8.) * SQRT(21./PI)
     &      * SIN(THETA) * (5*COS(THETA)**2 - 1)
     &      * CMPLX(COS(PHI), SIN(PHI))

      RETURN

      END

************************************************************************
*     ---------------------------------
      COMPLEX FUNCTION RSY32(THETA,PHI)
*     ---------------------------------
C
C     (Purpose)
C        Calculate spherical harmonics function Y32
C
C     (Input)
C        Theta : theta      
C        phi   : phi
C     (Output)    
C        Y32
C
C     (Creation Date and Author)
C        2016 Oct C. Wret: needed for pion ejection update Rein-Sehgal
C
************************************************************************      

      REAL*4 THETA,PHI

      PARAMETER(PI = 3.14159265)

      RSY32 = (1./4.) * SQRT(105./(2.*PI))
     &      * SIN(THETA)**2 * COS(THETA)
     &      * CMPLX(COS(2.*PHI), SIN(2.*PHI))

      RETURN

      END

************************************************************************
*     ---------------------------------
      COMPLEX FUNCTION RSY33(THETA,PHI)
*     ---------------------------------
C
C     (Purpose)
C        Calculate spherical harmonics function Y33
C
C     (Input)
C        Theta : theta      
C        phi   : phi
C     (Output)    
C        Y33
C
C     (Creation Date and Author)
C        2016 Oct C. Wret: needed for pion ejection update Rein-Sehgal
C
************************************************************************      

      REAL*4 THETA,PHI

      PARAMETER(PI = 3.14159265)

      RSY33 = (-1./8.) * SQRT(35./PI)
     &      * SIN(THETA)**3
     &      * CMPLX(COS(3.*PHI), SIN(3.*PHI))

      RETURN

      END
