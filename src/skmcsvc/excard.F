C========================================================================
      SUBROUTINE EXCARD(CCARD, LIN, LOUT, IRET)
C
C     FUNCTIONAL DESCRIPTION:
C       This routine extracts a card from a joined card.
C       The format of the joined card is
C
C         cardname-xxx
C         ...
C      
C       The characters after "-" are extracted.
C
C     INPUTS:
C       CARD : name of the card (e.g. "GEANT")
C       LIN  : logical unit number of the joined card
C       LOUT : logical unit number of the output card
C
C     OUTPUT:
C       IRET : return code 0=normal, -1=card does not exist
C========================================================================
C

      CHARACTER*132 LINE, LUPPER, CMATCH
      CHARACTER*(*) CCARD
      INTEGER IEND

      CMATCH = CCARD
      LM = INDEX(CMATCH,' ') - 1

      IEND = 0
      IRET = -1
      DO 100 WHILE (.TRUE.)
         READ(LIN,'(A)', END=1000) LINE

         NCH = 132
         DO 110 WHILE (NCH.GT.1 .AND. LINE(NCH:NCH).EQ.' ')
            NCH = NCH - 1
 110     CONTINUE

         NCT = 1
         DO 120 WHILE (NCT.NE.NCH .AND. LINE(NCT:NCT).NE.'-')
            NCT = NCT + 1
 120     CONTINUE
         IF(NCT.EQ.NCH) GOTO 100

         LUPPER = LINE
         CALL FFUPCA(LUPPER, 1, NCH)
         IF (LUPPER(1:NCT-1) .EQ. CMATCH(1:LM)) THEN
            IRET = 0
            WRITE(LOUT, '(A)') LINE(NCT+1:NCH)
            IF (LUPPER(NCT+1:NCH).EQ.'END') IEND = 1
         END IF
 100  CONTINUE
C
c end of input card
c    rewind input card
c    write a line 'END' if it does not exist in the output card
C
 1000 CONTINUE
      IF (IEND .EQ. 0) THEN
         WRITE(LOUT, '(A)') 'END'
      END IF
      REWIND LIN
      RETURN
      END
