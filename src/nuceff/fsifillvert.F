*     -----------------------------------------------------
      SUBROUTINE FSIFILLVERT(POS,IFLG)
*
*     ( purpose )
C        Fill FSI history common block with vertex info
*
*     ( input )
C         POS(3) : Position of vertex (fm)
C         IFLG   : Interaction type at vertex (*10 FOR HI-NRG model)
C                  -1 : ESCAPE
C                   0 : INITIAL
C                   1 :
C                   2 :
C                   3 : ABSORPTION
C                   4 : CHARGE EXCHANGE
C                   5 :
C                   6 :
C                   7 : HADRON PRODUCTION (hi-nrg only, i.e. 70)
C                   8 : QUASI-ELASTIC SCATTER
C                   9 : FORWARD (ELASTIC-LIKE) SCATTER
*
*     -----------------------------------------------------
      implicit none

      REAL*4 POS(3)
      INTEGER*4 IFLG, I

C#define DEBUG

#include "fsihist.h"

C     Do not add any more vertices if there is no more space for
C     the intermediate particle (or more vertices)
      if (NVERT.ge.MAXVERT.or.NVCVERT.ge.MAXVCVERT) then
        return
      end if

      NVERT=NVERT+1

#ifdef DEBUG
      write(*,*) "FSIFILLVERT(): NVERT = ",NVERT
#endif

      DO 1 I=1,3
         POSVERT(I,NVERT)=POS(I)
#ifdef DEBUG
         write (*,*) "POSVERT() = ", POSVERT(I,NVERT)
#endif

 1    continue

      IFLGVERT(NVERT)=IFLG

#ifdef DEBUG
      write (*,*) "IFLGVERT(", NVERT, ") = ",IFLGVERT(NVERT)
#endif
      END

