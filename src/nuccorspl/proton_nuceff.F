C
C     TEST PROGRAM TO READ VECTOR
C
      IMPLICIT NONE

#include <f77/f77_floatingpoint.h>

#include "nework.h"
#include "vcwork.h"
#include "vcvrtx.h"
#include "skhead.h"

#include "skparm.h"
#include "sktq.h"


      CHARACTER*32 CNAME
      INTEGER IDATA(20)

      INTEGER*4 LUNI
      PARAMETER(LUNI=10)
      INTEGER*4 LUNO
      PARAMETER(LUNO=20)

      INTEGER*4  LRECL
      PARAMETER (LRECL=1024)

      INTEGER*4  ISIZE
      PARAMETER (ISIZE=5000000)

      REAL*4     H(ISIZE)
      COMMON/PAWC/H

      INTEGER*4  I,LBANK,IERR,IEVENT,NDATA,ISEED
      real*4 dummy
      real*4   ranf,rndm
      external ranf,rndm


      CALL ieee_handler('set','common',SIGFPE_ABORT)
      CALL ieee_handler('set','division',SIGFPE_ABORT)
      CALL ieee_handler('set','invalid',SIGFPE_ABORT)


C
C Initialize
C
      CALL NECARD
      CALL NECARDVC
C      CALL NETARGCD

      SK_GEOMETRY=1
      CALL KZINIT
      CALL HLIMIT(-1*ISIZE)
      CALL GEOSET
C
C Open Files
C
      CALL SKOPENF(LUNI,0,'Z',IERR)
      if (ierr.lt.0) then
         write(*,*) 'Failed to open input file.'
         if (ierr.eq.-1001) then
            write(*,*) 'Please specify the file in the shell script'
         endif
         stop
      endif

      CALL SKOPENF(LUNO,0,'Z',IERR)
      if (ierr.lt.0) then
         write(*,*) 'Failed to open output file.'
         if (ierr.eq.-1001) then
            write(*,*) 'Please specify the file in the shell script'
         endif
         stop
      endif

C
C --  READ VECTOR
C
      DO 10 I=1,200000
         CALL SKOPTN("31,30,29")
         CALL KZECLR
         CALL SKREAD(LUNI,1000,2000,10,10)
         CALL KZBLOC('MCVECT',LBANK)
         IF (LBANK.eq.0) GOTO 2000
         IF (IERR.EQ.2) GOTO 1000
         IF (IERR.EQ.1) GOTO 2000
         write(6,*)  ' *** ',i,' ***'

         IEVENT=I

C -- initialize random seed everytime
         IF (i.eq.1) then
             CNAME='MCHEAD'
             CALL KZBLOC(CNAME,LBANK)
             IF (LBANK.EQ.0) then
                write(*,*) '*** Can not locate MCHEAD ***'
                stop
             endif
             CALL KZGET1(CNAME,0,NDATA,IDATA)
             IF (NDATA.EQ.0 ) THEN
                write(*,*) '*** Can not read ISEED ***'
                stop
             endif
             IF (NDATA.EQ.8) THEN
                ISEED=IDATA(7)
             else
                write(*,*) '** correct seed **'
                call seedranf(idata(8))
                dummy = ranf()
                CALL RLUXGO(3,idata(8),idata(9),idata(10))
             endif
         endif

         CALL VCRDVCCM
         CALL KZBDEL('MCVECT')
C-- make bank for the backup ( MCVECTBF )
         call vcmkvcbf

         call nrsetxpos
         CALL nrintr

         CALL VCMKVC

         call KZWRIT(LUNO)

 10   CONTINUE
      GOTO 3000

 1000 WRITE(6,*) ' READ ERROR : ',LUNI
      GOTO 3000

 2000 WRITE(6,*) ' END OF FILE : ',LUNI

 3000 CONTINUE

      CALL SKCLOSEF(LUNI)
      CALL SKCLOSEF(LUNO)

      STOP
      END

