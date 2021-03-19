      SUBROUTINE CCOPNF(FNAME,FSTAT,FACCES,FUNIT,IUNI)
      LOGICAL*1 TAPE,YES,NO
      CHARACTER*(*) FUNIT,FACCES,FSTAT,FNAME
      YES=.TRUE.
      NO=.FALSE.
      IF(FUNIT.EQ.'tape')THEN
        TAPE=YES
      ELSE
        TAPE=NO
      ENDIF
C
      IRONLY=0
      IWONLY=0
      IRDWR=0
      ICREAT=0
      IAPPND=0
      IEXCL=0
      IF(TAPE)THEN
        IF(FSTAT.EQ.'old')IRONLY=1
        IF(FSTAT.EQ.'new')IWONLY=1
      ELSE
        IF(FSTAT.EQ.'unknown' .OR. FSTAT.EQ.' ')THEN
          ICREAT=1
          IRDWR=1
        ENDIF
        IF(FACCES.EQ.'append')THEN
          IAPPND=1
          IF(FSTAT.EQ.'new')STOP
        ELSE
          IF(FSTAT.EQ.'old')IRDWR=1
          IF(FSTAT.EQ.'new')THEN
            ICREAT=1
            IEXCL=1
            IRDWR=1
          ENDIF
        ENDIF
      ENDIF
C
      CALL CCOPEN(FNAME,IRONLY,IWONLY,IRDWR,ICREAT,IAPPND,IEXCL,IUNI)
      IF(IUNI.LT.1)THEN
        IF(TAPE)STOP
        IF(FSTAT.EQ.'old')THEN
          IRONLY=1
          IRDWR=0
          CALL CCOPEN(FNAME,
     &                  IRONLY,IWONLY,IRDWR,ICREAT,IAPPND,IEXCL,IUNI)
          IF(IUNI.GE.1)GOTO 100
        ELSEIF(FSTAT.EQ.'new' .OR. FSTAT.EQ.'unknown'
     &                        .OR. FSTAT.EQ.' ')THEN
          IWONLY=1
          IRDWR=0
          CALL CCOPEN(FNAME,
     &                  IRONLY,IWONLY,IRDWR,ICREAT,IAPPND,IEXCL,IUNI)
          IF(IUNI.GE.1)GOTO 100
        ENDIF
        STOP
      ENDIF
C
  100 CONTINUE
      RETURN
      END
