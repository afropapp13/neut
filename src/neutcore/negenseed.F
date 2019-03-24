      SUBROUTINE NEGENSEED
*
*        Generate RANDOM-SEED from the file
*
*
      IMPLICIT NONE

      INTEGER*4  getpid,time
#ifndef gFortran
      external  getpid,time
#endif


      REAL*4     RANF
      EXTERNAL   RANF

      INTEGER*4  IPID,ITIME,IDIV
      INTEGER*4  ISEED
      INTEGER*4  ISEED1,ISEED2,istau1,istau2,istau3
      
      REAL*4     DUMMY
      
      
      IPID=getpid()
      ITIME=time()

C---
      IDIV=16384
 10   CONTINUE
         ISEED1 = ITIME+(IPID*IDIV)*(-1**MOD(IPID,2))
         IF (ISEED1.gt.0) goto 20
         IDIV = IDIV/2
      goto 10
            
 20   ISEED = ISEED1

C Now ranf is calling RLU (2008/11/22, hayato)
CC--- use same seed for RANLUX 
C      call seedranf(iseed)
C      dummy = ranf()
      CALL RLUXGO(3,ISEED1,0,0)

C---
      IDIV=32768
 110   CONTINUE
         ISEED2 = (ITIME/2)+(IPID*IDIV)*(-1**MOD(IPID,2))
         IF (ISEED2.gt.0) goto 120
         IDIV = IDIV/2
      goto 110
 120  CALL RM48IN(ISEED2,0,0)


      IDIV=32768
 210  CONTINUE
         ISTAU1 = ITIME+(IPID*IDIV)*(-1**MOD(IPID,2))*(-1)
         ISTAU1 = MOD(ISTAU1,65536*1024)
         IF (ISTAU1.gt.0) goto 220
         IDIV = IDIV/2
      goto 210

 220  IDIV=IDIV*2

 230  CONTINUE
      ISTAU2=0
C         ISTAU2 = (ITIME/2)+(IPID*IDIV)*(-1**MOD(IPID,2))*(-1)
C         ISTAU2 = MOD(ISTAU2,65536*1024)
C         IF (ISTAU2.gt.0) goto 240
C         IDIV = IDIV/2
C      goto 230

 240  IDIV=IDIV*4

 250  CONTINUE
      ISTAU3=0
C         ISTAU3 = (ITIME/4)+(IPID*IDIV)*(-1**MOD(IPID,2))*(-1)
C         ISTAU3 = MOD(ISTAU3,10)
C         IF (ISTAU3.gt.0) goto 300
C         IDIV = IDIV/2
C      goto 250


 300  CALL RMARIN(istau1,istau2,istau3)


      RETURN
      END
      
      
