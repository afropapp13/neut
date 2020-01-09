************************************************************************
*     -------------------
      SUBROUTINE NECARDEV
*     -------------------
*
*     (Purpose)
*       Read EVCT card and setup common blocks.
*
*     (Input)
*       card file
*
*     (Output)
*
*     (Creation Date and Author)
*
************************************************************************
      IMPLICIT NONE

#include "necardev.h"
#include "cardname.h"

      CHARACTER*40 TMPFIL, COMMAND

      INTEGER*4     IRET

      INTEGER*4     IARGC
#ifndef gFortran
      EXTERNAL      IARGC
#endif

C --  Initialize

      NEVTEVCT   = 0
      IDPTEVCT   = 0

      MPOSEVCT   = 0
      CALL VZERO(3,POSEVCT)

      RADEVCT    = 0.0

      MDIREVCT   = 0
      CALL VZERO(3,DIREVCT)

      MPVEVCT    = 0
      CALL VZERO(2,PVEVCT)

C --- get card file name (1st argument of the program)

      WRITE(6,*)
      CALL RESOLVECARDNAME
      WRITE(6,*) 'NEUT(EVCT) card file you read is ', FCARD(1:LEN_TRIM(FCARD))

C -- extract NEUT card

      CALL GETUNI(TMPFIL)
      OPEN(UNIT=90, FILE=FCARD, STATUS='OLD')
      OPEN(UNIT=91, FILE=TMPFIL,STATUS='UNKNOWN')
      CALL EXCARD('EVCT', 90, 91, IRET)
      CLOSE(UNIT=91)
      CLOSE(UNIT=90)

      OPEN(UNIT=91, FILE=TMPFIL,STATUS='OLD')
      CALL FFINIT(0)
      CALL FFSET('LINP', 91)
      CALL FFSET('SIZE', 10)

      CALL FFKEY('NEVT',NEVTEVCT,1,'INTEGER')

      CALL FFKEY('IDPT',IDPTEVCT,1,'INTEGER')

      CALL FFKEY('MPOS',MPOSEVCT,1,'INTEGER')
      CALL FFKEY('POS' ,POSEVCT ,3,'REAL')
      CALL FFKEY('RAD' ,RADEVCT ,1,'REAL')

      CALL FFKEY('MDIR',MDIREVCT,1,'INTEGER')
      CALL FFKEY('DIR' ,DIREVCT ,3,'REAL')

      CALL FFKEY('MPV' ,MPVEVCT ,1,'INTEGER')
      CALL FFKEY('PV'  ,PVEVCT  ,2,'REAL')

      CALL FFKEY('FILENM',FILENMEVCT,80,'MIXED')
      CALL FFKEY('HISTNM',HISTNMEVCT,80,'MIXED')
      CALL FFKEY('INMEV' ,INMEVEVCT ,1,'INTEGER')
      CALL FFKEY('OUTPUTFORMAT' ,OUTPUTFORMAT ,1,'INTEGER')

      CALL FFGO
      CLOSE(UNIT=91)
      COMMAND = '/bin/rm '//TMPFIL
      CALL SYSTEM(COMMAND)

      WRITE(6,*) 'NECARDEV : INPUTTED PARAMETER '
      WRITE(6,*) 'NECARDEV : NEV  = ',NEVTEVCT

      IF (      (abs(IDPTEVCT).ne.12
     $     .and.(abs(IDPTEVCT).ne.14))) THEN
         WRITE(6,*) 'NECARDEV : IDPT in the card was set to ',
     $              IDPTEVCT,' but this should be +-12, +-14.'
         stop
      else
         WRITE(6,*) 'NECARDEV : IDPT = ',IDPTEVCT
      endif
      if ((mposevct.ne.1).and.(mposevct.ne.2)) then
         write(6,*) 'NECARDEV : MPOSEVCT in the card was set to ',
     $              MPOSEVCT,' but this should be 1 or 2.'
         stop
      else
         WRITE(6,*) 'NECARDEV : MPOS = ',MPOSEVCT
      endif
      IF (MPOSEVCT.EQ.1) WRITE(6,*) 'NECARDEV : POS = ',POSEVCT
      IF (MPOSEVCT.EQ.2) WRITE(6,*) 'NECARDEV : RAD = ',RADEVCT

      if ((mdirevct.ne.1).and.(mdirevct.ne.2)) then
         write(6,*) 'NECARDEV : MDIREVCT in the card was set to ',
     $              MDIREVCT,' but this should be 1 or 2.'
         stop
      else
         WRITE(6,*) 'NECARDEV : MDIR = ',MDIREVCT
      endif
      IF (MDIREVCT.EQ.1) WRITE(6,*) 'NECARDEV : DIR = ',DIREVCT


      if (MPVEVCT.eq.1) then
         WRITE(6,*) 'NECARDEV : MPV  = ',MPVEVCT
         WRITE(6,*) 'NECARDEV : PV  = ',PVEVCT(1)
      else if (MPVEVCT.eq.2) then
         WRITE(6,*) 'NECARDEV : MPV  = ',MPVEVCT
         WRITE(6,*) 'NECARDEV : PV  =[',PVEVCT(1),
     $        ',',PVEVCT(2),']'
      else if (MPVEVCT.eq.3) then
         WRITE(6,*) 'NECARDEV : MPV  = ',MPVEVCT
         WRITE(6,*) 'Random ENEU generation using flux histogram'
      else
         WRITE(6,*) 'NECARDEV : MPV in the card was set to ',
     $              MPVEVCT,' but this should be 1, 2 or 3.'
         stop
      endif

      if (OUTPUTFORMAT.eq.1) then
        write (*,*) 'Writing output vector in neutvect format'
#ifdef USE_HEPMC
      else if (OUTPUTFORMAT.eq.2) then
        write (*,*) 'Writing output vector in HepMC3 format'
#endif
      else
        write(*,*) 'Unknown output format', OUTPUTFORMAT, ' specified.'
        stop
      endif

      RETURN
      END
