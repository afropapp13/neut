*     -----------------------------------------------------
      SUBROUTINE EFSCATT(IP,PPI,PNUC,MODE,IPF,PPIF,INUCF,PNUCF,radius)      
*     -----------------------------------------------------
*     
*     ( purpose )
*       SET SCATTERING PION MOMENTUM
*
*     ( input ) 
C         IP   : INITIAL PION KIND
C         PPI  : INITIAL PION MOMENTUM IN LAB. SYSTEM(MeV)
C         PNUC : INITIAL NUCLEON MOMENTUM(MeV)
C         MODE : INTERACTION MODE
*       
*     ( output ) 
C         IPF  : FINAL PION KIND
C         PPIF : FINAL PION MOMENTUM(MeV)
c         PNUCF: FINAL NUCLEON MOMENTUM
C         INUCF: Final Nucleon kind
*
C
C     ( update )
C         2009.12.09: Added ejected nucleon - P. de Perio

      DIMENSION PPI(3),PNUC(3),PPIF(3),PNUCF(3)
      DIMENSION SK(3),SKD(3),GK(3),EV(3),XMAP(11,11)
      DATA AMN/939./,AMPI/139./
C     
C     -- CALCULATE PI-N RELATIVE MOMENTUM IN C-M SYSTEM  I.E. SK(3)
      GK(1)=PNUC(1)+PPI(1)
      GK(2)=PNUC(2)+PPI(2)
      GK(3)=PNUC(3)+PPI(3)
      EPI=SQRT(AMPI**2+PPI(1)**2+PPI(2)**2+PPI(3)**2)
      ENU=SQRT(AMN**2+PNUC(1)**2+PNUC(2)**2+PNUC(3)**2)
      GKABS=SQRT(GK(1)**2+GK(2)**2+GK(3)**2)
      BETA=GKABS/(EPI+ENU)
      GM=1./SQRT(1.-BETA**2)
      IF(GKABS.EQ.0.)GO TO 140
      EV(1)=-GK(1)/GKABS
      EV(2)=-GK(2)/GKABS
      EV(3)=-GK(3)/GKABS
      GO TO 150
  140 EV(1)=1.
      EV(2)=0.
      EV(3)=0.
  150 SK(1)=PNUC(1)
      SK(2)=PNUC(2)
      SK(3)=PNUC(3)
      CALL MCVECBST(SK,AMN,EV,GM)
C     -- CALCULATE AMPLITUDE OF PHASE-SHIFT
      CALL EFCALAMP(SK,GK,PPI,radius)
C     -- MAKE COS(THETA)-PHI AMPLITUDE MAP
      CALL EFMAKMAP(PNUC,PPI,SK,GK,EV,GM,MODE,XMAP,radius)
C     
      XMAX=0.
      DO 300 I=1,11
         DO 310 J=1,11
            IF(XMAX.GT.XMAP(I,J))GO TO 310
            XMAX=XMAP(I,J)
  310    CONTINUE
  300 CONTINUE
C  320 COST=-1.+RNDM(DUM)*2.
  320 COST=-1.+RLU(IDUM)*2.
C      PHI=RNDM(DUM)*6.283185
      PHI=RLU(IDUM)*6.283185
C      HIT=XMAX*RNDM(DUM)
      HIT=XMAX*RLU(IDUM)
      COA=(COST+1.)/0.2
      PHA=PHI/0.6283185
      NC=IFIX(COA)+1
      NP=IFIX(PHA)+1
      IF(NC.GE.11 .OR. NP.GE.11) then
         WRITE(6,600)
         GO TO 320
      END IF
  600 FORMAT(' *** ERROR IN SCATT ... RETRY ***')
      DCO=COA+1.-FLOAT(NC)
      DPH=PHA+1.-FLOAT(NP)
      DAT1=(1.-DCO)*XMAP(NC,NP)+DCO*XMAP(NC+1,NP)
      DAT2=(1.-DCO)*XMAP(NC,NP+1)+DCO*XMAP(NC+1,NP+1)
      DISTRI=(1.-DPH)*DAT1+DPH*DAT2
      IF(HIT.GT.DISTRI)GO TO 320
C     
      CALL NECHADIR(COST,PHI,SK,PPIF)
      EV(1)=-EV(1)
      EV(2)=-EV(2)
      EV(3)=-EV(3)

C     -- SET NUCLEON FINAL MOMENTUM
      PNUCF(1)=PPIF(1)
      PNUCF(2)=PPIF(2)
      PNUCF(3)=PPIF(3)
      CALL MCVECBST(PNUCF,AMN,EV,GM)
C     -- SET FINAL NUCLEON KIND
      IF(MODE.EQ.1)INUCF=2212
      IF(MODE.EQ.2)INUCF=2112
      IF(MODE.EQ.3)INUCF=2212
      IF(MODE.EQ.4)INUCF=2212
      IF(MODE.EQ.5)INUCF=2112
      IF(MODE.EQ.6)INUCF=2112
      IF(MODE.EQ.7)INUCF=2212
      IF(MODE.EQ.8)INUCF=2212
      IF(MODE.EQ.9)INUCF=2112
      IF(MODE.EQ.10)INUCF=2112

C     -- SET PION FINAL MOMENTUM
      PPIF(1)=-PPIF(1)
      PPIF(2)=-PPIF(2)
      PPIF(3)=-PPIF(3)
      CALL MCVECBST(PPIF,AMPI,EV,GM)
C     -- SET FINAL PION KIND
      IF(MODE.EQ.1)IPF=211
      IF(MODE.EQ.2)IPF=211
      IF(MODE.EQ.3)IPF=111
      IF(MODE.EQ.4)IPF=111
      IF(MODE.EQ.5)IPF=111
      IF(MODE.EQ.6)IPF=211
      IF(MODE.EQ.7)IPF=-211
      IF(MODE.EQ.8)IPF=-211
      IF(MODE.EQ.9)IPF=-211
      IF(MODE.EQ.10)IPF=111
      RETURN
      END
