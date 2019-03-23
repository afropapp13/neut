************************************************************************
*     ------------------
      FUNCTION FNNPOT(INOUT,P)
*     ------------------
*
*     (Purpose)
*       NUCLEON-NUCLEOUS POTENTIOAL AS A FUNCTION OF NUCLEON MOMENTUM.
*
*     (Input)
*       INOUT  : INITIAL : 1
*                FINAL   : 2  
*       P      : NUCLEON MOMENTUM ( GEV/C )
*     (Output)
*       FNNPOT : POTENTIAL ( GEV/C )
*
*     (Creation Date and Author)
*       2003.03.17 ; Y. Hayato / Taken fron fnnucl
************************************************************************

      IMPLICIT NONE

#include "neutparams.h"

      REAL*4    FNNPOT
      real*4    P
      INTEGER*4 INOUT
      
      FNNPOT=0.0E0

      if (INOUT.eq.1) THEN
C--initial
         FNNPOT=VNUINI
         return
      else if (INOUT.eq.2) THEN
C--final
         FNNPOT=VNUFIN
         return
      else
         write(*,*) 'FNNPOT: BAD argument. inout=',inout
         stop
      endif

      end
      
         
