************************************************************************
*     ---------------------------------
      subroutine nrlibver(NUCCVER)
*     ---------------------------------
C
C     (Purpose)
C        Reterun the version of nuccorspl library
C
C     (Input)
C        NENUCCVER : DNENUCCVER in Imakefile
C
C     (Output)
C        NUCCVER   : version number
C
C     (Creation Date and Author)
C        2007.12.01 ; G.Mitsuka
C
************************************************************************

      integer NUCCVER, NUCCVERtmp
      data NUCCVERtmp/NENUCC_VERSION/

      NUCCVER = NUCCVERtmp

      return
      end
