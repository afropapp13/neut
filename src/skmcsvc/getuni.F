C====================================================================
C     GET AN UNIQUE NAME IN A SYSTEM
C     
C     (INPUT)
C        NONE
C     (OUTPUT)
C        STR : A STRING LIKE "/tmp/pidxxxx.tmp"
C====================================================================
      SUBROUTINE GETUNI(STR)

      CHARACTER*(*) STR
      INTEGER GETPID, PID

      PID = GETPID()

      WRITE(STR, '(A,I7.7,A)') '/tmp/pid', PID, '.tmp'

c unnecessary - I7.7 zero-pads the field
c      DO 100 I = 9, 14
c         IF (STR(I:I) .EQ. ' ') STR(I:I) = '0'
c 100  CONTINUE

      RETURN
      END
