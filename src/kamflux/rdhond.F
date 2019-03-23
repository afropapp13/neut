C   01/07/91 312072024  MEMBER NAME  RDHOND   (FORT)        FORTRAN
C
      SUBROUTINE RDHOND
C
C ++ READ FLUX CALCULATED BY HONDA, AND FILL INTO COMMON /HONDAF/
C                            (AND INTO COMMON /HONDAH/ ; OCT 22, 1993)
C
C
C         T.KAJITA  JULY-01,1991
C         T.KAJITA  OCT -22,1993   HIGH ENERGY FLUX SUPPORTED
C
C
      character*80 fname,unit,lgname,access,stg,multi,comment,hname,f8
      data lin/51/
      PARAMETER (NE=41)
C                NUMBER OF HIGH-ENERGY FLUX'S ENERGY BINS
      DIMENSION FLINT(4),FLDIF(4)
      COMMON /HONDAF/ ENELOW(72), PCOSL(20),
     &       FEMIN(20,72), FEMAX(20,72), FEBMIN(20,72), FEBMAX(20,72),
     &       FMMIN(20,72), FMMAX(20,72), FMBMIN(20,72), FMBMAX(20,72)
      COMMON /HONDAH/ ENEHI(NE), PCOSH(20),
     &       FEHMN(20,NE), FEHMX(20,NE), FEBHMN(20,NE), FEBHMX(20,NE),
     &       FMHMN(20,NE), FMHMX(20,NE), FMBHMN(20,NE), FMBHMX(20,NE)
      DO 3 I=1,20
  3   PCOSL(I)= 1.05-0.10*FLOAT(I)
      PCOSL(1)= 1.001
      PCOSL(20)=-1.001
C
      DO 1 ICOS =1,20
      NFILE=ICOS +50
C
C ++ open files
C
      ipt = icos         
      call rgetv(lin,ipt,fname,unit,lgname,
     &        access,stg,multi,comment,hname,f8)
      if(fname.eq.' ')then
          write(0,*)'rdhond : cannot open table file ',lin
          stop
      endif
      open(lin,file=fname,form='formatted',status='old')
C
C
      DO 10 IE=1,72
COLD      READ(NFILE,500,END=100) ENEU,(FLINT(I),I=1,4),(FLDIF(J),J=1,4),
      READ(LIN,500,END=100) ENEU,(FLINT(I),I=1,4),(FLDIF(J),J=1,4),
     & FNENM, RNMNMB, RNENEB
 500  FORMAT(E10.4,4F8.2,4F9.2,3F6.3)
CC    WRITE(6,510) IE,PCOSL(ICOS),ENEU,(FLINT(I),I=1,4),
CC   &             (FLDIF(J),J=1,4), FNENM, RNMNMB, RNENEB
CC510 FORMAT(I3,F5.2,2X,E10.4,4F8.2,4F9.2,3F6.3)
C
      ENELOW(IE)=ENEU
      FEMIN(ICOS,IE)=FLDIF(3)
      FEBMIN(ICOS,IE)=FLDIF(4)
      FMMIN(ICOS,IE)=FLDIF(1)
      FMBMIN(ICOS,IE)=FLDIF(2)
      FEMAX(ICOS,IE)=FLDIF(3)
      FEBMAX(ICOS,IE)=FLDIF(4)
      FMMAX(ICOS,IE)=FLDIF(1)
      FMBMAX(ICOS,IE)=FLDIF(2)
c      write(0,*) (fldif(i),i=1,4)
C
  10  CONTINUE
C
C ++ CLOSE FILE
C
      CLOSE(LIN)
C
C
 100  CONTINUE
   1  CONTINUE
C
C
C READ HIGH-ENERGY ATMOSPHERIC-NEUTRINO FLUX   OCT 22, 1993
      CALL RDHONH
C
      RETURN
      END
