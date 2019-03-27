*     -----------------------------------------------------
      SUBROUTINE EFPNCRSH(PAIMOM,NQPAI,NQN,SIGTOT,RSIGEL,IERROR)
*     -----------------------------------------------------
*
*     ( purpose )
C        CROSS SECTION OF PAI-N INTERACTION AT HIGH ENERGIES (>800MEV).
C                    T.KAJITA  25-JUN-88
*
*     ( input )
C           PAIMOM : MOMENTUM OF PION IN MEV/C
C           NQPAI  : CHARGE OF PION  <=== NO MEANING AT PRESENT
C           NQN    : CHARGE OF NUCLEON<=== NO MEANING AT PRESENT
*
*     ( output )
C           SIGTOT : TOTAL CROSS SECTIONIN MB.
C           RSIGEL : (ELASTIC CROSS SECTION) / (TOTAL CROSS SECTION)
C           IERROR : ERROR CODE  =0 -- NO ERROR
C                                =1 -- ERROR ( PAIMOM < 800MEV )
C
C     (ASSUMPTIONS)
C  (1) INDIVIDUAL TOTAL CROSS SECTIONS ARE SIGMA(PAI(+/-)-DEUT) /2.0
C     FOR ALL COMBINATIONS.
C  (2) FRACTION OF ELASTIC CROSS SECTIONS
C           SIGMA(PAI(-)-N) IS SAME AS SIGMA(PAI(+)-P)
C           SIGMA(PAI(+)-N) AND SIGMA(PAI(0)-N) ARE SAME AS SIGMA(PAI(-)-P)
C
C   2008/03 Added cross sections for 340 & 450 MeV/c (nishino)
C
C
C
      IMPLICIT NONE
      REAL PAIMOM, SIGTOT, RSIGEL
      INTEGER NQPAI,NQN, IERROR
      INTEGER IE, IEL, IEH
      REAL DEL, DEH
      REAL STSUM, SESUM
      REAL PMEV,STPPMD, STPPP, STPMP, SEPPP, SEPMP
      DIMENSION PMEV(16), STPPMD(16), STPPP(16), STPMP(16),
     &     SEPPP(16), SEPMP(16)
C
C     -- PION MOMENTUM IN MEV/C
      DATA PMEV  /  340.,  450.,  500.,   600.,   700.,   800.,
     $     900.,  1000.,  1200., 1500.,  2000.,  3000.,  5000.,
     $     10000., 20000., 50000./
C
C     -- TOTAL CROSS SECTION PAI(+/-)-DEUT.
      DATA STPPMD/ 180.,  55.,   60.,    52.,   58.0,   56.0,
     $     62.0,   77.0, 68.0, 74.0,   65.0,   60.0,   56.0,
     $     50.0,   48.0,   47.0 /
C
C     -- TOTAL CROSS SECTION PAI(+)-P
      DATA STPPP / 134.,  48.,  35.,    20.,    17.0,   16.0,
     $     20.0,   25.0, 28.0,  40.0,   29.0,   29.0,   26.0,
     $     25.0,   24.0,   24.0  /
C
C     -- TOTAL CROSS SECTION PAI(-)-P
      DATA STPMP /   49., 26.,   27.,    30.,   49.0,   40.0,
     $     40.0,   57.0,  37.0, 36.0,   36.0,   32.0,   29.0,
     $     27.0,   25.0,   24.0  /
C
C     -- ELASTIC CROSS SECTION PAI(+)-P
      DATA SEPPP /   134., 46. , 35.,    20.,   13.0,    9.0,
     $     8.2,   12.0,   14.0, 18.0,   10.0,    7.8,    5.7,
     $     4.8,    4.0,    3.3  /
C
C     -- ELASTIC CROSS SECTION PAI(-)-P
      DATA SEPMP / 17.,  10.,   10.,    13.,   28.0,   16.0,
     $     16.0,   25.0, 15.0, 12.0,    9.5,    7.5,    6.0,
     $     4.6,    4.0,    3.4/
C
C
C
      RSIGEL=0.
      SIGTOT=0.
      IERROR=1
      IF(PAIMOM.LT.340.)RETURN
      IERROR=0
C
      DO 10 IE=2,16
         IF(PMEV(IE).GT.PAIMOM) THEN
            IEH=IE
            IEL=IE-1
            DEH=ABS(PMEV(IEH)-PAIMOM)
            DEL=ABS(PAIMOM-PMEV(IEL))
            GO TO 11
         END IF
   10 CONTINUE
      IEH=16
      IEL=15
      DEH=0.
      DEL= ABS (PMEV(IEH)-PMEV(IEL))
   11 CONTINUE
C
C
      SIGTOT = ( DEL*STPPMD(IEH) + DEH*STPPMD(IEL) ) / ( DEL+DEH ) / 2.
      IF(NQPAI.EQ.0) THEN
         STSUM = ( DEL*STPMP(IEH) + DEH*STPMP(IEL) ) / ( DEL+DEH )
         SESUM = ( DEL*SEPMP(IEH) + DEH*SEPMP(IEL) ) / ( DEL+DEH )
      ELSE IF((NQPAI+NQN).EQ.2 .OR. (NQPAI+NQN).EQ.-1) THEN
         STSUM = ( DEL*STPPP(IEH) + DEH*STPPP(IEL) ) / ( DEL+DEH )
         SESUM = ( DEL*SEPPP(IEH) + DEH*SEPPP(IEL) ) / ( DEL+DEH )
      ELSE IF((NQPAI+NQN).EQ.1 .OR. (NQPAI+NQN).EQ.0) THEN
         STSUM = ( DEL*STPMP(IEH) + DEH*STPMP(IEL) ) / ( DEL+DEH )
         SESUM = ( DEL*SEPMP(IEH) + DEH*SEPMP(IEL) ) / ( DEL+DEH )
      ELSE
         PRINT *, 'Error in total charge in nuclear effect'
         STOP
      ENDIF

CCCC  PRINT *,' PAIMOM SIGTOT STSUM SESUM',PAIMOM,SIGTOT,STSUM,SESUM
      RSIGEL=SESUM/STSUM
CCCC  PRINT *,' RSIGEL',RSIGEL
C
C
      RETURN
      END
