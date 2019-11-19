************************************************************************
*     -----------------
      SUBROUTINE TACARD
*     -----------------
*
*     (Purpose)
*       Read TAU card and setup common blocks.
*
*     (Input)
*       card file
*
*     (Output)
*
*     (Creation Date and Author)
*       1995.1.30 ; Koshio ( SGCARD )
*       1996.1.11 ; K.Kaneyuki
*       2004.2.04 ; C.Saji
************************************************************************
#include "dsbnkf.h"
#include "tacard.h"
#include "cardname.h"

      CHARACTER*40 TMPFIL, COMMAND

C --  Initialize

      ITAMODE = 0
      CALL VZERO(22,TAUBRA)
      ITAPOL = 0
      CALL VZERO(3,TAPOL)
      CALL VZERO(3,ITAURNDM)

C --- get card file name (1st argument of the program)

      WRITE(6,*)
      CALL RESOLVECARDNAME
      WRITE(6,*) 'TAU card file you read is ', FCARD(1:LENSTR(FCARD))

C -- extract NEUT card

      CALL GETUNI(TMPFIL)
      write(6,'(A)') TMPFIL
      OPEN(UNIT=90, FILE=FCARD, STATUS='OLD')
      OPEN(UNIT=91, FILE=TMPFIL,STATUS='UNKNOWN')
      CALL EXCARD('TAU', 90, 91, IRET)
      CLOSE(UNIT=91)
      CLOSE(UNIT=90)

      OPEN(UNIT=91, FILE=TMPFIL,STATUS='OLD')
      CALL FFINIT(0)
      CALL FFSET('LINP', 91)
      CALL FFSET('SIZE', 5)
      CALL FFSET('LENG',100)
      CALL FFKEY('MODE',ITAMODE, 1,'INTE')
      CALL FFKEY('BRA' ,TAUBRA, 22,'REAL')
      CALL FFKEY('IPOL',ITAPOL,  1,'INTE')
      CALL FFKEY('POL' ,TAPOL,   3,'REAL')
      CALL FFKEY('RNDM',ITARNDM, 3,'INTE')
      CALL FFKEY('GENE',ITAGENE, 1,'INTE')
      CALL FFGO
      CLOSE(UNIT=91)
      COMMAND = '/bin/rm '//TMPFIL
      CALL SYSTEM(COMMAND)

      WRITE(6,*)
      WRITE(6,*) 'TACARD : INPUTTED PARAMETER '
      WRITE(6,*) 'TACARD : MODE = ',ITAMODE
      IF (ITAMODE.EQ.-1) WRITE(6,*) 'TACARD : BRA  = ',TAUBRA
      WRITE(6,*) 'TACARD : IPOL = ',ITPOL
      IF (ITAPOL.EQ.1) WRITE(6,*) 'TACARD : POL  = ',TAPOL
      IF (ITARNDM(1).NE.0) THEN
        WRITE(6,*) 'TACARD : RNDM  = ',ITARNDM
        CALL RMARIN(ITARNDM(1),ITARNDM(2),ITARNDM(3))
      ELSE
         WRITE(6,*) "!! SET TAU-RNDM VALUE IN CARD(tacard.F) !!"
         STOP
      ENDIF

      RETURN
      END
