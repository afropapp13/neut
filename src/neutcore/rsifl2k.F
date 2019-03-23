***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE RSIFL2K(ASUM,ETABW,RESSGN,AM,ASL,ASR,ASS)
*     ---------------------------------------------------
*
*     ( purpose )
*       calculate helicity cross-sections
*
*     ( input )
*       asum(2,6,14): amplitude of each resonance state.
*       etabw(31)   : Breit-Wigner factor
*       ressgn(31)  : decay sign
*       am(6.31)    : amplitude of each resonance.
*     
*     ( output )
*       asl(3)      : helicity amplitude
*       asr(3)      : helicity amplitude
*       ass(3)      : helicity amplitude
*
*     ( creation date and author )
*       1993.Dec. ; for Super-Kamioka by Y.Hayato
*                   (Based on Rein and Sehgal's code)
*                   for K
*     ( comment )
*       Reference:Rein and Sehgal,Ann. of Phys. 133(1981),79-153
*                 Rein           ,Z.Phys.C 35(1987),43-64
*                 Feynman et al. ,Phys.Rev.D3(1971),2706-
***********************************************************************

      IMPLICIT NONE

      COMPLEX ASUM(2,6,14)
      COMPLEX ETABW(31)
      REAL  RESSGN(31)
      REAL AM(6,31)

      REAL ASL(3),ASR(3),ASS(3)

      COMPLEX S31M3N,S11M3N,P31M3N,P11M3N,P33M3N,P13M3N,D33M3N,D13M3N
      COMPLEX D35M3N,D15M3N,F35M3N,F15M3N,F37M3N,F17M3N
      COMPLEX S31P3N,S11P3N,P31P3N,P11P3N,P33P3N,P13P3N,D33P3N,D13P3N
      COMPLEX D35P3N,D15P3N,F35P3N,F15P3N,F37P3N,F17P3N
      COMPLEX S31M1N,S11M1N,P31M1N,P11M1N,P33M1N,P13M1N,D33M1N,D13M1N
      COMPLEX D35M1N,D15M1N,F35M1N,F15M1N,F37M1N,F17M1N
      COMPLEX S31P1N,S11P1N,P31P1N,P11P1N,P33P1N,P13P1N,D33P1N,D13P1N
      COMPLEX D35P1N,D15P1N,F35P1N,F15P1N,F37P1N,F17P1N
      COMPLEX S310MN,S110MN,P310MN,P110MN,P330MN,P130MN,D330MN,D130MN
      COMPLEX D350MN,D150MN,F350MN,F150MN,F370MN,F170MN
      COMPLEX S310PN,S110PN,P310PN,P110PN,P330PN,P130PN,D330PN,D130PN
      COMPLEX D350PN,D150PN,F350PN,F150PN,F370PN,F170PN

      COMPLEX S31M3P,S11M3P,P31M3P,P11M3P,P33M3P,P13M3P,D33M3P,D13M3P
      COMPLEX D35M3P,D15M3P,F35M3P,F15M3P,F37M3P,F17M3P
      COMPLEX S31P3P,S11P3P,P31P3P,P11P3P,P33P3P,P13P3P,D33P3P,D13P3P
      COMPLEX D35P3P,D15P3P,F35P3P,F15P3P,F37P3P,F17P3P
      COMPLEX S31M1P,S11M1P,P31M1P,P11M1P,P33M1P,P13M1P,D33M1P,D13M1P
      COMPLEX D35M1P,D15M1P,F35M1P,F15M1P,F37M1P,F17M1P
      COMPLEX S31P1P,S11P1P,P31P1P,P11P1P,P33P1P,P13P1P,D33P1P,D13P1P
      COMPLEX D35P1P,D15P1P,F35P1P,F15P1P,F37P1P,F17P1P
      COMPLEX S310MP,S110MP,P310MP,P110MP,P330MP,P130MP,D330MP,D130MP
      COMPLEX D350MP,D150MP,F350MP,F150MP,F370MP,F170MP
      COMPLEX S310PP,S110PP,P310PP,P110PP,P330PP,P130PP,D330PP,D130PP
      COMPLEX D350PP,D150PP,F350PP,F150PP,F370PP,F170PP
                  
      COMPLEX BKGM3P,BKGM1P,BKGP3P,BKGP1P,BKG0PP,BKG0MP

      S31M3P =ASUM( 1 , 1 , 1 )  
      S11M3P =ASUM( 1 , 1 , 2 )  
      P31M3P =ASUM( 1 , 1 , 3 )  
      P11M3P =ASUM( 1 , 1 , 4 )  
      P33M3P =ASUM( 1 , 1 , 5 )  
      P13M3P =ASUM( 1 , 1 , 6 )  
      D33M3P =ASUM( 1 , 1 , 7 )  
      D13M3P =ASUM( 1 , 1 , 8 )  
      D35M3P =ASUM( 1 , 1 , 9 )  
      D15M3P=ASUM( 1 , 1 , 10 ) 
      F35M3P=ASUM( 1 , 1 , 11 ) 
      F15M3P=ASUM( 1 , 1 , 12 ) 
      F37M3P=ASUM( 1 , 1 , 13 ) 
      F17M3P=ASUM( 1 , 1 , 14 ) 

      S31M1P =ASUM( 1 , 2 , 1 )  
      S11M1P =ASUM( 1 , 2 , 2 )  
      P31M1P =ASUM( 1 , 2 , 3 )  
      P11M1P =ASUM( 1 , 2 , 4 )  
      P33M1P =ASUM( 1 , 2 , 5 )  
      P13M1P =ASUM( 1 , 2 , 6 )  
      D33M1P =ASUM( 1 , 2 , 7 )  
      D13M1P =ASUM( 1 , 2 , 8 )  
      D35M1P =ASUM( 1 , 2 , 9 )  
      D15M1P=ASUM( 1 , 2 , 10 ) 
      F35M1P=ASUM( 1 , 2 , 11 ) 
      F15M1P=ASUM( 1 , 2 , 12 ) 
      F37M1P=ASUM( 1 , 2 , 13 ) 
      F17M1P=ASUM( 1 , 2 , 14 ) 

      S31P3P =ASUM( 1 , 3 , 1 )  
      S11P3P =ASUM( 1 , 3 , 2 )  
      P31P3P =ASUM( 1 , 3 , 3 )  
      P11P3P =ASUM( 1 , 3 , 4 )  
      P33P3P =ASUM( 1 , 3 , 5 )  
      P13P3P =ASUM( 1 , 3 , 6 )  
      D33P3P =ASUM( 1 , 3 , 7 )  
      D13P3P =ASUM( 1 , 3 , 8 )  
      D35P3P =ASUM( 1 , 3 , 9 )  
      D15P3P=ASUM( 1 , 3 , 10 ) 
      F35P3P=ASUM( 1 , 3 , 11 ) 
      F15P3P=ASUM( 1 , 3 , 12 ) 
      F37P3P=ASUM( 1 , 3 , 13 ) 
      F17P3P=ASUM( 1 , 3 , 14 ) 

      S31P1P =ASUM( 1 , 4 , 1 )  
      S11P1P =ASUM( 1 , 4 , 2 )  
      P31P1P =ASUM( 1 , 4 , 3 )  
      P11P1P =ASUM( 1 , 4 , 4 )  
      P33P1P =ASUM( 1 , 4 , 5 )  
      P13P1P =ASUM( 1 , 4 , 6 )  
      D33P1P =ASUM( 1 , 4 , 7 )  
      D13P1P =ASUM( 1 , 4 , 8 )  
      D35P1P =ASUM( 1 , 4 , 9 )  
      D15P1P=ASUM( 1 , 4 , 10 ) 
      F35P1P=ASUM( 1 , 4 , 11 ) 
      F15P1P=ASUM( 1 , 4 , 12 ) 
      F37P1P=ASUM( 1 , 4 , 13 ) 
      F17P1P=ASUM( 1 , 4 , 14 ) 

      S310PP =ASUM( 1 , 5 , 1 )  
      S110PP =ASUM( 1 , 5 , 2 )  
      P310PP =ASUM( 1 , 5 , 3 )  
      P110PP =ASUM( 1 , 5 , 4 )  
      P330PP =ASUM( 1 , 5 , 5 )  
      P130PP =ASUM( 1 , 5 , 6 )  
      D330PP =ASUM( 1 , 5 , 7 )  
      D130PP =ASUM( 1 , 5 , 8 )  
      D350PP =ASUM( 1 , 5 , 9 )  
      D150PP=ASUM( 1 , 5 , 10 ) 
      F350PP=ASUM( 1 , 5 , 11 ) 
      F150PP=ASUM( 1 , 5 , 12 ) 
      F370PP=ASUM( 1 , 5 , 13 ) 
      F170PP=ASUM( 1 , 5 , 14 ) 

      S310MP =ASUM( 1 , 6 , 1 )  
      S110MP =ASUM( 1 , 6 , 2 )  
      P310MP =ASUM( 1 , 6 , 3 )  
      P110MP =ASUM( 1 , 6 , 4 )  
      P330MP =ASUM( 1 , 6 , 5 )  
      P130MP =ASUM( 1 , 6 , 6 )  
      D330MP =ASUM( 1 , 6 , 7 )  
      D130MP =ASUM( 1 , 6 , 8 )  
      D350MP =ASUM( 1 , 6 , 9 )  
      D150MP=ASUM( 1 , 6 , 10 ) 
      F350MP=ASUM( 1 , 6 , 11 ) 
      F150MP=ASUM( 1 , 6 , 12 ) 
      F370MP=ASUM( 1 , 6 , 13 ) 
      F170MP=ASUM( 1 , 6 , 14 ) 

      S31M3N =ASUM( 2 , 1 , 1 )  
      S11M3N =ASUM( 2 , 1 , 2 )  
      P31M3N =ASUM( 2 , 1 , 3 )  
      P11M3N =ASUM( 2 , 1 , 4 )  
      P33M3N =ASUM( 2 , 1 , 5 )  
      P13M3N =ASUM( 2 , 1 , 6 )  
      D33M3N =ASUM( 2 , 1 , 7 )  
      D13M3N =ASUM( 2 , 1 , 8 )  
      D35M3N =ASUM( 2 , 1 , 9 )  
      D15M3N=ASUM( 2 , 1 , 10 ) 
      F35M3N=ASUM( 2 , 1 , 11 ) 
      F15M3N=ASUM( 2 , 1 , 12 ) 
      F37M3N=ASUM( 2 , 1 , 13 ) 
      F17M3N=ASUM( 2 , 1 , 14 ) 

      S31M1N =ASUM( 2 , 2 , 1 )  
      S11M1N =ASUM( 2 , 2 , 2 )  
      P31M1N =ASUM( 2 , 2 , 3 )  
      P11M1N =ASUM( 2 , 2 , 4 )  
      P33M1N =ASUM( 2 , 2 , 5 )  
      P13M1N =ASUM( 2 , 2 , 6 )  
      D33M1N =ASUM( 2 , 2 , 7 )  
      D13M1N =ASUM( 2 , 2 , 8 )  
      D35M1N =ASUM( 2 , 2 , 9 )  
      D15M1N=ASUM( 2 , 2 , 10 ) 
      F35M1N=ASUM( 2 , 2 , 11 ) 
      F15M1N=ASUM( 2 , 2 , 12 ) 
      F37M1N=ASUM( 2 , 2 , 13 ) 
      F17M1N=ASUM( 2 , 2 , 14 ) 

      S31P3N =ASUM( 2 , 3 , 1 )  
      S11P3N =ASUM( 2 , 3 , 2 )  
      P31P3N =ASUM( 2 , 3 , 3 )  
      P11P3N =ASUM( 2 , 3 , 4 )  
      P33P3N =ASUM( 2 , 3 , 5 )  
      P13P3N =ASUM( 2 , 3 , 6 )  
      D33P3N =ASUM( 2 , 3 , 7 )  
      D13P3N =ASUM( 2 , 3 , 8 )  
      D35P3N =ASUM( 2 , 3 , 9 )  
      D15P3N=ASUM( 2 , 3 , 10 ) 
      F35P3N=ASUM( 2 , 3 , 11 ) 
      F15P3N=ASUM( 2 , 3 , 12 ) 
      F37P3N=ASUM( 2 , 3 , 13 ) 
      F17P3N=ASUM( 2 , 3 , 14 ) 

      S31P1N =ASUM( 2 , 4 , 1 )  
      S11P1N =ASUM( 2 , 4 , 2 )  
      P31P1N =ASUM( 2 , 4 , 3 )  
      P11P1N =ASUM( 2 , 4 , 4 )  
      P33P1N =ASUM( 2 , 4 , 5 )  
      P13P1N =ASUM( 2 , 4 , 6 )  
      D33P1N =ASUM( 2 , 4 , 7 )  
      D13P1N =ASUM( 2 , 4 , 8 )  
      D35P1N =ASUM( 2 , 4 , 9 )  
      D15P1N=ASUM( 2 , 4 , 10 ) 
      F35P1N=ASUM( 2 , 4 , 11 ) 
      F15P1N=ASUM( 2 , 4 , 12 ) 
      F37P1N=ASUM( 2 , 4 , 13 ) 
      F17P1N=ASUM( 2 , 4 , 14 ) 

      S310PN =ASUM( 2 , 5 , 1 )  
      S110PN =ASUM( 2 , 5 , 2 )  
      P310PN =ASUM( 2 , 5 , 3 )  
      P110PN =ASUM( 2 , 5 , 4 )  
      P330PN =ASUM( 2 , 5 , 5 )  
      P130PN =ASUM( 2 , 5 , 6 )  
      D330PN =ASUM( 2 , 5 , 7 )  
      D130PN =ASUM( 2 , 5 , 8 )  
      D350PN =ASUM( 2 , 5 , 9 )  
      D150PN=ASUM( 2 , 5 , 10 ) 
      F350PN=ASUM( 2 , 5 , 11 ) 
      F150PN=ASUM( 2 , 5 , 12 ) 
      F370PN=ASUM( 2 , 5 , 13 ) 
      F170PN=ASUM( 2 , 5 , 14 ) 

      S310MN =ASUM( 2 , 6 , 1 )  
      S110MN =ASUM( 2 , 6 , 2 )  
      P310MN =ASUM( 2 , 6 , 3 )  
      P110MN =ASUM( 2 , 6 , 4 )  
      P330MN =ASUM( 2 , 6 , 5 )  
      P130MN =ASUM( 2 , 6 , 6 )  
      D330MN =ASUM( 2 , 6 , 7 )  
      D130MN =ASUM( 2 , 6 , 8 )  
      D350MN =ASUM( 2 , 6 , 9 )  
      D150MN=ASUM( 2 , 6 , 10 ) 
      F350MN=ASUM( 2 , 6 , 11 ) 
      F150MN=ASUM( 2 , 6 , 12 ) 
      F370MN=ASUM( 2 , 6 , 13 ) 
      F170MN=ASUM( 2 , 6 , 14 ) 


      ASL(2) = 3. *(  ABS(S11P3P)**2 + ABS(P11P3P)**2
     $              + ABS(P13P3P)**2 + ABS(D13P3P)**2
     $              + ABS(D15P3P)**2 + ABS(F15P3P)**2
     $              + ABS(F17P3P)**2 
     $              + ABS(S11P1P)**2 + ABS(P11P1P)**2
     $              + ABS(P13P1P)**2 + ABS(D13P1P)**2
     $              + ABS(D15P1P)**2 + ABS(F15P1P)**2
     $              + ABS(F17P1P)**2)

      ASR(2) = 3. *(  ABS(S11M3P)**2 + ABS(P11M3P)**2
     $              + ABS(P13M3P)**2 + ABS(D13M3P)**2
     $              + ABS(D15M3P)**2 + ABS(F15M3P)**2
     $              + ABS(F17M3P)**2 
     $              + ABS(S11M1P)**2 + ABS(P11M1P)**2
     $              + ABS(P13M1P)**2 + ABS(D13M1P)**2
     $              + ABS(D15M1P)**2 + ABS(F15M1P)**2
     $              + ABS(F17M1P)**2)
      
      ASS(2) = 3. *(  ABS(S110MP)**2 + ABS(P110MP)**2
     $              + ABS(P130MP)**2 + ABS(D130MP)**2
     $              + ABS(D150MP)**2 + ABS(F150MP)**2
     $              + ABS(F170MP)**2 
     $              + ABS(S110PP)**2 + ABS(P110PP)**2
     $              + ABS(P130PP)**2 + ABS(D130PP)**2
     $              + ABS(D150PP)**2 + ABS(F150PP)**2
     $              + ABS(F170PP)**2)
      
      
      BKGM3P = AM(1,29)*ETABW(29)*RESSGN(29)
      BKGM1P = AM(2,29)*ETABW(29)*RESSGN(29)
      BKGP3P = AM(3,29)*ETABW(29)*RESSGN(29)
      BKGP1P = AM(4,29)*ETABW(29)*RESSGN(29)
      BKG0PP = AM(5,29)*ETABW(29)*RESSGN(29)
      BKG0MP = AM(6,29)*ETABW(29)*RESSGN(29)
      

      ASL(1) = 3. *(  ABS(S11P3P)**2 + ABS(P11P3P)**2
     $              + ABS(P13P3P)**2 + ABS(D13P3P)**2
     $              + ABS(D15P3P)**2 + ABS(F15P3P)**2
     $              + ABS(F17P3P)**2 + ABS(BKGP3P)**2
     $              + ABS(S11P1P)**2 + ABS(P11P1P)**2
     $              + ABS(P13P1P)**2 + ABS(D13P1P)**2
     $              + ABS(D15P1P)**2 + ABS(F15P1P)**2
     $              + ABS(F17P1P)**2 + ABS(BKGP1P)**2)
      
      ASR(1) = 3. *(  ABS(S11M3P)**2 + ABS(P11M3P)**2
     $              + ABS(P13M3P)**2 + ABS(D13M3P)**2
     $              + ABS(D15M3P)**2 + ABS(F15M3P)**2
     $              + ABS(F17M3P)**2 + ABS(BKGM3P)**2
     $              + ABS(S11M1P)**2 + ABS(P11M1P)**2
     $              + ABS(P13M1P)**2 + ABS(D13M1P)**2
     $              + ABS(D15M1P)**2 + ABS(F15M1P)**2
     $              + ABS(F17M1P)**2 + ABS(BKGM1P)**2)
      
      ASS(1) = 3. *(  ABS(S110MP)**2 + ABS(P110MP)**2
     $              + ABS(P130MP)**2 + ABS(D130MP)**2
     $              + ABS(D150MP)**2 + ABS(F150MP)**2
     $              + ABS(F170MP)**2 + ABS(BKG0MP)**2
     $              + ABS(S110PP)**2 + ABS(P110PP)**2
     $              + ABS(P130PP)**2 + ABS(D130PP)**2
     $              + ABS(D150PP)**2 + ABS(F150PP)**2
     $              + ABS(F170PP)**2 + ABS(BKG0PP)**2)
      
      
      
      P33M3P = AM(1, 1)*ETABW( 1)*RESSGN( 1)
      P33M1P = AM(2, 1)*ETABW( 1)*RESSGN( 1)
      P33P3P = AM(3, 1)*ETABW( 1)*RESSGN( 1)
      P33P1P = AM(4, 1)*ETABW( 1)*RESSGN( 1)
      P330PP = AM(5, 1)*ETABW( 1)*RESSGN( 1)
      P330MP = AM(6, 1)*ETABW( 1)*RESSGN( 1)
      
      
      ASL(3) = 3. * (ABS(P33P3P)**2 + ABS(P33P1P)**2)
      ASR(3) = 3. * (ABS(P33M3P)**2 + ABS(P33M1P)**2)
      ASS(3) = 3. * (ABS(P330MP)**2 + ABS(P330PP)**2)                  

      END
