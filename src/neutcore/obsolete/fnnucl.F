************************************************************************
*     ------------------
      FUNCTION FNNUCL(P)
*     ------------------
*
*     (Purpose)
*       NUCLEON-NUCLEOUS POTENTIOAL AS A FUNCTION OF NUCLEON MOMENTUM.
*
*     (Input)
*       P      : NUCLEON MOMENTUM ( GEV/C )
*     (Output)
*       FNNUCL : POTENTIAL ( GEV/C )
*
*     (Creation Date and Author)
*       1988.09.17 ; T.KAJITA
*          FROM    F.A.BRIEVA AND A.DELLAFIORE  NUCL.PHYS.A292(1977)445
*                  (PAGE449, FIG.1)
*       1995.02.08 ; K. KANEYUKI  FOR S.K.
*                    VNUCL -> FNNUCL
*                    CHANGE UNIT ( MEV/C -> GEV/C )      
*       1995.02.08 ; K. KANEYUKI  CHANGE THIS ROUTINE TO DUMMY
*       2000.09.05 ; K. Kaneyuki CHANGE IT FOR GEV/C UNIT
*       2002.09.18 ; Y. Hayato  / Decided to use fixed potential!!!
*     
*       2003.03.17 ; Y. Hayato / Stop if this function was called.
************************************************************************

      DIMENSION V(11)

      integer*4  getpid
      external   getpid

      integer*4  ipid,iret

C
C     DATA ARE FROM FIG.1 OF THE ABOVE REF.
C
      DATA V/-0.0545, -0.0525, -0.0475, -0.0345, -0.022, -0.013,
     &       -0.008 , -0.0045, -0.002 , -0.001 ,    0. /
C      DATA V/-54.5, -52.5, -47.5, -34.5, -22., -13., -8., -4.5, -2.,
C     &        -1.,    0. /


      FNNUCL=0.0

      ipid = getpid()
C--- use SIGNAL 3 : SIGQUIT to make core-dump file
      iret = kill(ipid,3)
      stop
      RETURN

C      IF(P.GT.0.9999) RETURN
C
C      I=INT(P*10.)+1
C      DI=P*10.-FLOAT(I-1)
C      FNNUCL=(1.-DI)*V(I)+DI*V(I+1)
C     
C------------------      IF(P.GT.999.9) RETURN
C------------------
C------------------      I=INT(P/100.)+1
C------------------      DI=P/100.-FLOAT(I-1)
C------------------      FNNUCL=(1.-DI)*V(I)+DI*V(I+1)

      END
