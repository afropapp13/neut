C   01/07/91 312072120  MEMBER NAME  RDHONH   (FORT)        FORTRAN
C
C
      SUBROUTINE RDHONH
C
C ++ READ HIGH-ENERGY NEUTRINO FLUX CALCULATED BY HONDA,
C ++ AND FILL INTO COMMON /HONDAH/.
C
C         T.KAJITA  OCT-14,1993
C
C
      character*80 fname,unit,lgname,access,stg,multi,comment,hname,f8
      data lin/52/
      PARAMETER (NE=41)
C                NUMBER OF ENERGY BINS
      DIMENSION FLDIF(4)
      COMMON /HONDAH/ ENEHI(NE), PCOSH(20),
     &       FEHMN(20,NE), FEHMX(20,NE), FEBHMN(20,NE), FEBHMX(20,NE),
     &       FMHMN(20,NE), FMHMX(20,NE), FMBHMN(20,NE), FMBHMX(20,NE)
c
c      do 200 i=81,90
c         write(fname,'(\'fort.\',i2)')i
c         open(i,file=fname,form='formatted',status='old')
c 200  continue
C
      DO 3 I=1,20
  3   PCOSH(I)= 1.05-0.10*FLOAT(I)
      PCOSH(1)= 1.001
      PCOSH(20)=-1.001
C
      DO 1 ICOS =1,10
      NFILE=ICOS +80
      ICDN=ICOS+10
      ICUP=-ICOS+11
C
C ++ open files
C
      ipt = icos
      call rgetv(lin,ipt,fname,unit,lgname,
     &        access,stg,multi,comment,hname,f8)
      if(fname.eq.' ')then
          write(0,*)'rdhond : cannot open table file ',ipt
          stop
      endif
      open(lin,file=fname,form='formatted',status='old')
C
C
      DO 10 IE=1,NE
COLD      READ(NFILE,500,END=100) ENEU,(FLDIF(I),I=1,4),
      READ(LIN,500,END=100) ENEU,(FLDIF(I),I=1,4),
     & FNEMMB, FNUEEB
 500  FORMAT(7(1X,E10.5))
C     WRITE(6,510) IE,PCOSH(ICOS),ENEU,
C    &             (FLDIF(J),J=1,4), FNEMMB, FNUEEB
C 510 FORMAT(I3,F5.2,2X,E11.4,6E11.4)
C
      ENEHI(IE)=ENEU
      FEHMN(ICDN,IE)=FLDIF(3)
      FEBHMN(ICDN,IE)=FLDIF(4)
      FMHMN(ICDN,IE)=FLDIF(1)
      FMBHMN(ICDN,IE)=FLDIF(2)
      FEHMX(ICDN,IE)=FLDIF(3)
      FEBHMX(ICDN,IE)=FLDIF(4)
      FMHMX(ICDN,IE)=FLDIF(1)
      FMBHMX(ICDN,IE)=FLDIF(2)
C
      FEHMN(ICUP,IE)=FLDIF(3)
      FEBHMN(ICUP,IE)=FLDIF(4)
      FMHMN(ICUP,IE)=FLDIF(1)
      FMBHMN(ICUP,IE)=FLDIF(2)
      FEHMX(ICUP,IE)=FLDIF(3)
      FEBHMX(ICUP,IE)=FLDIF(4)
      FMHMX(ICUP,IE)=FLDIF(1)
      FMBHMX(ICUP,IE)=FLDIF(2)
C
  10  CONTINUE
C
C ++ CLOSE FILE
C
      CLOSE(LIN)
C
 100  CONTINUE
   1  CONTINUE
C
c      do 201 i=81,90
c         close(i)
c 201  continue
c
      RETURN
      END
