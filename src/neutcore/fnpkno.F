**********************************************************************
*     ------------------------
      FUNCTION FNPKNO(N,XMEAN)
*     ------------------------
*     
*     ( purpose )
*        PROBABILITY OF PION MULTIPLICITY    M.NAKAHATA
*        (EXECUTED FROM KNO SCALING)
*
*     ( input ) 
*        N     ; PION MULTIPLICITY
*        XMEAN ; MEAN PION MULTIPLICITY
*
*     ( output )
*        FNPKNO : PROBABILIRY OF PION MULTIPLICITY
*
*     ( creation date and author )
*       198?.??.?? ; First version by M.Nakahata
*       1988.06.27 ; T.KAJITA      PKNO(N,W) ---> PKNO(N,XMEAN)
*       1995.03.03 ; K.Kaneyuki for S.K.
*                    PKNO -> FNPKNO
*
**********************************************************************
      DIMENSION PROB(16)
      DATA PROB/1.E-1 , 1.7E-1 , 5.5E-1 , 1.3 , 2. , 2.1 , 1.8 , 1.2 ,
     &     6.5E-1 , 2.6E-1 ,1.E-1 ,3.9E-2 , 1.5E-2 , 6.1E-3 ,
     &     2.5E-3 , 1.E-3/
C     
C     -- MEAN VALUE OF PION MULTIPLICITY
C     
      FNPKNO = 0.
      RATIO=FLOAT(N)/XMEAN
      IBIN=IFIX(RATIO/0.2)+1
      IF(IBIN.GE.16) RETURN
      XLOG=ALOG10(PROB(IBIN+1))
      XLOG1=ALOG10(PROB(IBIN))
      PLOG=(XLOG-XLOG1)*(RATIO-0.2*FLOAT(IBIN-1))/0.2+XLOG1
      FNPKNO=10**PLOG
      RETURN
      END
