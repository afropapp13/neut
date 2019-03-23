************************************************************************
*     ---------------------------------
      subroutine eflibver(NUCEVER)
*     ---------------------------------
C
C     (Purpose)
C        Reterun the version of nuceff library
C
C     (Input)
C        NENUCEVER : DNENUCEVER in Imakefile
C
C     (Output)    
C        NUCEVER   : version number
C
C     (Creation Date and Author)
C        2007.12.01 ; G.Mitsuka
C
************************************************************************      

      integer NUCEVER, NUCEVERtmp
      data NUCEVERtmp/NENUCEVER/

      NUCEVER = NUCEVERtmp

      return
      end
