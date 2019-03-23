************************************************************************
*     --------
*     TACARD.H
*     --------
*
*     (Purpose)
*       COMMON for CARD on TAU
*
*     (Variables)
*       ITAMODE : TAU DECAY MODE
*                 0  : ALL DECAY MODE WITH NORMAL BR. ( default ) 
*                 -1 : ALL DECAY MODE WITH BR. IN TAUBRA(22)
*                 1-22 : ONLY SELLECTED MODE
*       TAUBRA(22) : BRANCHING RATION 
*
*       ITAPOL  : TAU+ POLARITY 
*                 0 : CONSIDER POLARITY ( default )
*                     NOW, TAPOL = ( 0., 0., 1. )
*                 1 : FIX POLARITY IN TAPOL
*       TAPOL(3) : POLARITY
*       ITARNDM(3) : RANDUM SEED FOR TAUOLA
*                     IF ITAURNDM(1)=0, RANDUM SEED IS NOT SET.
*     
*     (Creation Date and Author)
*       1996.1.11 ; K.Kaneyuki
*
************************************************************************

      COMMON/TAUCARD/ITAMODE,TAUBRA(22),ITAPOL,TAPOL(3),
     &               ITARNDM(3)


