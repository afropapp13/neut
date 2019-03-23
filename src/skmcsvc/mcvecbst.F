***********************************************************************
*     --------------------------------------
      SUBROUTINE MCVECBST(PMOM,XMASS,UV,GAM)
*     --------------------------------------
*     
*     ( purpose )
*       ROLENTZ BOOST
*
*     ( input ) 
*       PMOM(3) : INITIAL MOMENTUM
*       XMASS   : MASS OF PARTICLE
*       UV(3)   : DIRECTION OF BOOST ( UNIT VECTOR )
*       GAM     : GAMMA FACTOR
*
*     ( output ) 
*       PMOM(3) : FINAL MOMENTUM
*
*     ( creation date and author )
*       1979.10.11 ; First version by S.Odaka
*       1992. 9.18 ; Mod by K.Kaneyuki for S.K.
*                       VECBST -> VCBOST
*       1995. 2. 7 ; change name VCBOST -> MCVECBST
*     
*     ( comment )
*
***********************************************************************
      REAL PMOM(3),UV(3)
     
      PP=PMOM(1)*UV(1)+PMOM(2)*UV(2)+PMOM(3)*UV(3)
      ES=PMOM(1)**2+PMOM(2)**2+PMOM(3)**2+XMASS**2
      CE=(GAM-1.)*PP+SQRT((GAM*GAM-1.)*ES)
     
      DO 10 I=1,3
         PMOM(I)=CE*UV(I)+PMOM(I)
 10   CONTINUE
     
      RETURN
      END
