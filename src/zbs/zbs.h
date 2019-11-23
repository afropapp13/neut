C----------------------------------------------------------------------
C             MACRO ZBS
C----------------------------------------------------------------------
      integer    NPDIV   ,NPOFFS
      PARAMETER (NPDIV=20,NPOFFS=8)
      integer    NPSTOR       ,NPDATA
C      PARAMETER (NPSTOR=250000,NPDATA=NPSTOR-NPOFFS)

C     NPSTOR 250000 -> 500000 96-Apr-24 K.Kaneyuki
C     NPSTOR 500000 -> 750000 96-Jun-01 Y.Itow
C     NPSTOR 750000 ->2000000 04-Mar-21 Y.Koshio

      PARAMETER (NPSTOR=2000000,NPDATA=NPSTOR-NPOFFS)

      integer       IXSTOR,IXDIV(NPDIV),IFENCE,LZBS(NPSTOR)
      COMMON /KZBS/ IXSTOR,IXDIV       ,IFENCE,LZBS
      SAVE /KZBS/

      INTEGER IZBS(NPDATA)
      REAL    RZBS(NPDATA)
      EQUIVALENCE (LZBS(NPOFFS+1),IZBS(1),RZBS(1))


      integer        IQUEST(100)
      COMMON /QUEST/ IQUEST
      SAVE /QUEST/

      integer         IFZEBI
      COMMON /KZFLAG/ IFZEBI
      SAVE /KZFLAG/
C------------------------ END OF MACRO ZBS -------------------------
