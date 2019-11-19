      SUBROUTINE NERDSEED
*
*        READ RANDOM-SEED from the file
*
*
      IMPLICIT NONE

      real*4     ranf
      external   ranf

      integer*4      LENSTR
      external       lenstr

      INTEGER*4  ISEED
      INTEGER*4  ISEED1,ISEED2,istau1,istau2,istau3
      CHARACTER*80 RANFILE

      real*4     DUMMY


      call getenv('RANFILE',ranfile)

*     Generate a new seed if the ranfile isn't set.
      IF (LENSTR(ranfile) .EQ. 0) THEN
        CALL NEGENSEED
        RETURN
      endif

      write(*,*) "Now opening random seed file(",ranfile,", len = ",LEN(ranfile),")"
      open(unit=99,file=ranfile,status='unknown')
      READ(99,*) ISEED1,ISEED2,istau1,istau2,istau3

      ISEED = ISEED1

C Now ranf is calling RLU (2008/11/22, hayato)
CC--- use same seed for RANLUX :
C      call seedranf(iseed)
C      dummy = ranf()

C--- Changed from RNDM to RANLUX
      CALL RLUXGO(3,ISEED1,0,0)

C--- Changed from RAMBO to RM48
      CALL RM48IN(ISEED2,0,0)

#ifdef NEUT_USE_TAUOLA
      CALL RMARIN(istau1,istau2,istau3)
#else
      write(*,*) "************************************************"
      write(*,*) "*****  RMARIN function was not called as   *****"
      write(*,*) "*****  tauola support was not built.       *****"
      write(*,*) "*****  If this is unexpected, rebuild with *****"
      write(*,*) "*****  tauola.                             *****"
      write(*,*) "************************************************"
#endif

      RETURN
      END
