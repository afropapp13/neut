**********************************************************************
*     --------------------
      FUNCTION NEMULTIP(W)
*     --------------------
*     
*     ( purpose )
*        GIVE PION MULTIPLICITY     M.NAKAHATA
*
*     ( input ) 
*        W : ( GEV )
*
*     ( output )
*        NEMULTIP : PION MULTIPLICITY
*
*     ( creation date and author )
*       198?.??.?? ; First version by M.Nakahata
*       1988.06.27 ; T.KAJITA      PKNO(N,W) ---> PKNO(N,XMEAN)
*       1995.03.03 ; K.Kaneyuki for S.K.
*                    MULTIP -> NEMULTIP
*
**********************************************************************
      DIMENSION BUF(10)

      XMEAN=1.83*ALOG(W**2)
C
C     MEAN PION MULTIPLICITY IN NEUTRINO INTERACTIONS
C     
      BUF(1)=FNPKNO(1,XMEAN)
      DO 100 I=2,10
         BUF(I)=BUF(I-1)+FNPKNO(I,XMEAN)
  100 CONTINUE
  110 RAN=RLU(DUM)*BUF(10)
      I=1
  150 IF(BUF(I).GE.RAN)GO TO 200
      I=I+1
      IF(I.EQ.11)GO TO 110
      GO TO 150
  200 NEMULTIP=I
      RETURN
      END
