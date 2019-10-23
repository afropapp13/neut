************************************************************************
*     -------------------
      SUBROUTINE NECARDBM
*     -------------------
*
*     (Purpose)
*       Read NUBM card and setup common blocks.
*
*     (Input)
*       card file
*
*     (Output)
*
*     (Creation Date and Author)
*       1995.01.30 ; Koshio ( SGCARD )
*       2007.11.18 ; Y.Hayato
*
************************************************************************
      IMPLICIT NONE

#include "necardbm.h"
#include "cardname.h"

      CHARACTER*40 TMPFIL, COMMAND
      INTEGER*4 IRET


C --  Initialize

C-- To support numu flux with nue interactions,
C     flavor and particle should be specified independently.
C
C-- flux flavor
      IDBMFLX = 0
C-- interacting particle * neutrino *
      IDBMPID = 0
C-- Number of events to be generated
      NEVBM    = 0
C-- Detector ID ( ND 1 to X, SK = 0 )
      IDBMDET  = -2
C-- Maximum Radius ( maybe limited by the flux ) in cm
      BMRADMX = 0.

C --- get card file name (1st argument of the program)

      WRITE(6,*)
      CALL RESOLVECARDNAME
      WRITE(6,*) '======= CONFIGURATION FOR ATMPD VECT ======='
      WRITE(6,*) 'NEUT(NEBM) card file you read is ',FCARD

C -- extract NEUT card

      CALL GETUNI(TMPFIL)
      OPEN(UNIT=90, FILE=FCARD, STATUS='OLD')
      OPEN(UNIT=91, FILE=TMPFIL,STATUS='UNKNOWN')
      CALL EXCARD('NEBM', 90, 91, IRET)
      CLOSE(UNIT=91)
      CLOSE(UNIT=90)

      OPEN(UNIT=91, FILE=TMPFIL,STATUS='OLD')
      CALL FFINIT(0)
      CALL FFSET('LINP', 91)
      CALL FFSET('SIZE', 10)
      CALL FFKEY('FLX'  ,IDBMFLX, 1,'INTEGER')
      CALL FFKEY('PID'  ,IDBMPID, 1,'INTEGER')
      CALL FFKEY('NEVNT',NEVBM,   1,'INTEGER')
      CALL FFKEY('DETID',IDBMDET, 1,'INTEGER')
      CALL FFKEY('RADMX',BMRADMX, 1,'REAL')
      CALL FFGO
      CLOSE(UNIT=91)
      COMMAND = '/bin/rm '//TMPFIL
      CALL SYSTEM(COMMAND)

      WRITE(6,*) 'NECARDBM : INPUTTED PARAMETER '
      WRITE(6,*) 'NECARDBM : # of events =',NEVBM

      if (IDBMFLX.ne.-1) then
         IF ((ABS(IDBMFLX).ne.12).and.(ABS(IDBMFLX).ne.14)) then
            WRITE(6,*) 'NECARDBM : INVALID FLUX (flavor) ',IDBMFLX,
     $           ' was specified.'
            STOP
         endif
      else
         WRITE(6,*) 'NECARDBM : FLUX FLAVOR =',IDBMFLX
      endif

      IF ((ABS(IDBMPID).ne.12).and.(ABS(IDBMPID).ne.14)) then
         WRITE(6,*) 'NECARDBM : INVALID particle ',IDBMPID,
     $        ' was specified.'
         STOP
      else
         WRITE(6,*) 'NECARDBM : PID         =',IDBMPID
      endif

      IF (IDBMDET.lt.-1) then
         WRITE(6,*) 'NECARDBM : INVALID DETECTOR ID ',IDBMDET,
     $        ' was specified.'
         STOP
      else
         WRITE(6,*) 'NECARDBM : DETECTOR ID =',IDBMDET
      endif

      IF (BMRADMX.lt.0.) then
         BMRADMX = 1.e6
      endif
      WRITE(6,*) 'NECARDBM : MAX RADIUS  =',BMRADMX, 'cm'

      RETURN
      END
