************************************************************************
*     ----------------------------------
      SUBROUTINE BBBA05(Q2,GPE,GPM,GNE,GNM)
*
* Nucl. Phys. Proc. Suppl. 159;127-132,2006
*
*       (Porpose)
*         Calculate vector form factor
*
*       (input)
*         Q**2  ( Q**2 > 0 )
*
*       (output)
*         Vector form factor
*         GPE,GNE : charge vector form factor for proton and neutron
*         GPM,GNM : magnetic vector form factor for proton and neutron
*
*     ----------------------------------
      IMPLICIT NONE
  
      INTEGER I   
      REAL*8 PM,Q2,TAU

      REAL*8 GAPM,GAPE,GANM,GANE
      REAL*8 GBPM,GBPE,GBNM,GBNE
      REAL*8 GPE,GPM,GNE,GNM

      REAL*8 APM(3),APE(3),ANE(3),ANM(3)

      REAL*8 BPM(4),BPE(4),BNE(4),BNM(4)

      DATA PM/0.938/

      DATA APE/1., -0.0578, 0./
      DATA APM/1., 0.15,    0./
      DATA ANE/0., 1.25,    1.3/
      DATA ANM/1., 1.81,    0./

      DATA BPE/11.1,  13.6, 33.0,  0./
      DATA BPM/11.1,  19.6, 7.54,  0./
      DATA BNE/-9.86, 305., -758., 802./
      DATA BNM/14.1,  20.7, 68.7,  0./

C      write(*,*) "bbba05.F Q2 = ", Q2
      
      if (Q2.lt.0.) THEN
         write(*,*) 'Q2 should be larger than 0. (q2=',q2,')'
         stop
      endif

      TAU = 0.25*Q2/PM**2 
 
      GAPM = 0.
      GAPE = 0.
      GANM = 0.
      GANE = 0.
      DO I=0,2
       GAPM = GAPM + APM(I+1)*TAU**I
       GAPE = GAPE + APE(I+1)*TAU**I
       GANM = GANM + ANM(I+1)*TAU**I
       GANE = GANE + ANE(I+1)*TAU**I
      ENDDO

      GBPM = 0.
      GBPE = 0.
      GBNM = 0.
      GBNE = 0.
      DO I=1,4
       GBPM = GBPM + BPM(I)*TAU**I
       GBPE = GBPE + BPE(I)*TAU**I
       GBNM = GBNM + BNM(I)*TAU**I
       GBNE = GBNE + BNE(I)*TAU**I
      ENDDO

      GPM = 2.793*GAPM/(1.+GBPM)
      GPE =       GAPE/(1.+GBPE)
      GNM = -1.913*GANM/(1.+GBNM)
      GNE =        GANE/(1.+GBNE)

C      write(*,*) "bbba05.F GPE GPM GNE GNM = " , GPE   , GPM  ,  GNE , GNM

      RETURN
      END
