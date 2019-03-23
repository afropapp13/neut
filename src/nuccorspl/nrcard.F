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

      CHARACTER*256 TMPFIL, COMMAND
      CHARACTER*256 FCARD

      INTEGER*4 IARGC
#ifndef gFortran
      EXTERNAL  IARGC
#endif
C --  Initialize

      NUCRESCAT = 1
      XNUCFACT  = 1.
      NUCRESFLG = 1

C --- get card file name (1st argument of the program)

      WRITE(6,*) 
      IF (IARGC() .GE. 1) THEN
         CALL GETARG(1, FCARD)
      ELSE
         FCARD = 'neut.card'
      END IF
      WRITE(6,*) 'NUCLEON RESCATTERING card file you read is ',FCARD

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
      CALL FFKEY('FACT'  ,XNUCFACT, 1,'REAL')
      CALL FFKEY('FRMGAS',NUCRESFLG, 1,'INTEGER')

      CALL FFGO
      CLOSE(UNIT=91)
      COMMAND = '/bin/rm '//TMPFIL
      CALL SYSTEM(COMMAND)

      WRITE(6,*) 'NUCRES : INPUTTED PARAMETER '
      WRITE(6,*) 'NUCRES : RESCAT   = ',NUCRESCAT, ' 0: off, 1: on'
      WRITE(6,*) 'NUCRES : XNUCFACT = ',XNUCFACT
      WRITE(6,*) 'NUCRES : NUCRESFLG= ',NUCRESFLG
      RETURN
      END
