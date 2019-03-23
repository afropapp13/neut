***********************************************************************
*     -----------------
      SUBROUTINE nemknpos
*     -----------------
*
*     ( purpose )
*       Create NUCPOS bank
*
*     ( input ) 
*       none
*
*     ( output ) 
*       none
*
*     ( creation date and author )
*       1992.6.5 ; First version by K.Kaneyuki
*       1994.9.19; mod. for zbs V1.5
*       2006.6.12; Store all the positions
*
***********************************************************************
      IMPLICIT NONE

      INTEGER IDATA(20)
      REAL    RDATA(20)
      EQUIVALENCE(IDATA(1),RDATA(1))
      CHARACTER*32 CNAME

      integer*4 ISEG,NWMAST,LSEGM,NSEGM,IERR

#include "vcwork.h"
#include "nework.h"
#include "posinnuc.h"

C      COMMON /SAVXIN/INO16,XIN(3)

      IF (IBOUND.eq.0) then
         DO 100 ISEG=1,NVC
            POSNUC(1,ISEG)=9999.
            POSNUC(2,ISEG)=9999.
            POSNUC(3,ISEG)=9999.
 100     continue
      else

         if (ibound.eq.-1) then
            write(*,*) 'nemknpos: IBOUND flag is still -1..'
            write(*,*) 'nemknpos: mode=',modene
            stop
         endif
         DO 110 ISEG=1,NVC
            do 120 NWMAST=1,3
               if (POSNUC(NWMAST,ISEG).eq.-9999.) then
                  write(*,*) 'nemknpos: position is not set correctly.'
                  write(*,*) 'nemknpos: mode=',modene
                  write(*,*) 'nemknpos: Particle # ',ISEG, 
     $                 'PID=',IPVC(ISEG)
                  stop
               endif
 120        continue
 110     continue
         
      endif
C
C     -- MAKE EMB
C
C     LEM+1 : # OF VERTEX
C
C      IDATA(1)=1
      IDATA(1)=NVC

      CNAME='NUCPOS'

      NWMAST=1
      NSEGM=1
C      NSEGM=NVC
      LSEGM=3

      CALL KZBCR1(CNAME,NWMAST,'1I',NSEGM,LSEGM,'3F ',IERR)
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAME,0,NWMAST,IDATA)

      DO 200 ISEG=1,NVC
C     
C     -- SET LOCATION in the nucleus INFORMATION 
C     
C        IND+ 1 : X
C           + 2 : Y
C           + 3 : Z

         RDATA( 1)=POSNUC(1,ISEG)
         RDATA( 2)=POSNUC(2,ISEG)
         RDATA( 3)=POSNUC(3,ISEG)
C     
C     -- MAKE I-TH SEGMENT
C     
C         write(*,'(A,F12.5,A,F12.5,A,F12.5,A)')
C     $        'Position:(',posnuc(1,ISEG),
C     $                 ',',posnuc(2,ISEG),
C     $                 ',',posnuc(3,ISEG),')'
C
         CALL KZREP1(CNAME,ISEG,3,IDATA(1))
 200  continue

      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
  900 FORMAT(' *** ERROR *** ( BANK CREATION IN NUCPOS )')
      RETURN
      END
