************************************************************************
*     ---------------------------------
      subroutine nelibver(COREVER)
*     ---------------------------------
C
C     (Purpose)
C        Reterun the version of neutcore library
C
C     (Input)
C        NECOREVER : DNECOREVER in Imakefile
C
C     (Output)
C        COREVER   : version number
C
C     (Creation Date and Author)
C        2007.08.22 ; G.Mitsuka
C
************************************************************************

      integer COREVER, COREVERtmp
      data COREVERtmp/NECORE_VERSION/

      COREVER = COREVERtmp

      return
      end
