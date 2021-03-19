************************************************************************
*     -----------------
      SUBROUTINE NRCARD
*     -----------------
*
*     (Purpose)
*       Read NUCLEON RESCATTERING card and setup common blocks.
*
*     (Input)
*       card file
*
*     (Output)
*
*     (Creation Date and Author)
*       1995.1.30  ; Koshio ( SGCARD )
*       2003.03.30 ; for nucleon rescattering
*
************************************************************************
#include "nrcard.h"
#include "cardname.h"

      CHARACTER*256 TMPFIL, COMMAND
      
      LOGICAL  CARDISTOML
      EXTERNAL CARDISTOML
      
C --  Initialize

      NUCRESCAT = 1
      XNUCELAFACT = 1
      XNUCSPIFACT = 1
      XNUCDPIFACT = 1
      XNUCFACT  = 1.
      NUCRESFLG = 2

C --- get card file name (1st argument of the program)

      CALL RESOLVECARDNAME
      WRITE(6,*) '======= CONFIGURATION FOR INTERACTIONS ======='
      IF ( CARDISTOML() ) THEN
         WRITE(6,*) 'Reading', FCARD(1:LEN_TRIM(FCARD)), ' as a toml file.'
         CALL NRCARDTOML
      ELSE
         WRITE(6,*) 'Reading', FCARD(1:LEN_TRIM(FCARD)), ' as a card file.'

C -- extract NUCRES card

        CALL GETUNI(TMPFIL)
        OPEN(UNIT=90, FILE=FCARD, STATUS='OLD')
        OPEN(UNIT=91, FILE=TMPFIL,STATUS='UNKNOWN')
        CALL EXCARD('NUCRES', 90, 91, IRET)
        CLOSE(UNIT=91)
        CLOSE(UNIT=90)

        OPEN(UNIT=91, FILE=TMPFIL,STATUS='OLD')
        CALL FFINIT(0)
        CALL FFSET('LINP', 91)
        CALL FFSET('SIZE', 9)
        CALL FFKEY('RESCAT',NUCRESCAT, 1,'INTEGER')
        CALL FFKEY('ELAFACT'  ,XNUCELAFACT, 1,'REAL')
        CALL FFKEY('SPIFACT'  ,XNUCSPIFACT, 1,'REAL')
        CALL FFKEY('DPIFACT'  ,XNUCDPIFACT, 1,'REAL')
        CALL FFKEY('FACT'  ,XNUCFACT, 1,'REAL')
        CALL FFKEY('FRMGAS',NUCRESFLG, 1,'INTEGER')

        CALL FFGO
        CLOSE(UNIT=91)
        COMMAND = '/bin/rm '//TMPFIL
        CALL SYSTEM(COMMAND)
      ENDIF

      WRITE(6,*) 'NUCRES : INPUTTED PARAMETER '
      WRITE(6,*) 'NUCRES : RESCAT   = ',NUCRESCAT, ' 0: off, 1: on'
      WRITE(6,*) 'NUCRES : XNUCELAFACT = ',XNUCELAFACT
      WRITE(6,*) 'NUCRES : XNUCSPIFACT = ',XNUCSPIFACT
      WRITE(6,*) 'NUCRES : XNUCDPIFACT = ',XNUCDPIFACT
      WRITE(6,*) 'NUCRES : XNUCFACT = ',XNUCFACT
      WRITE(6,*) 'NUCRES : NUCRESFLG= ',NUCRESFLG
      RETURN
      END
