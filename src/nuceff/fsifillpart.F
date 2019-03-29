*     -----------------------------------------------------
      SUBROUTINE FSIFILLPART(DIR,ABSP,ABSPBST,IP,VERTI,VERTF)
*
*     ( purpose )
C        Fill FSI history common block with intermediate particle info
*
*     ( input )
C         DIR     : Direction of particle
C         ABSP    : Momentum in lab frame (MeV/c)
C         ABSPBST : Momentum in nucleon rest frame (MeV/c)
C         IP      : Particle ID (PDG)
C         VERTI   : Initial vertex index
C         VERTF   : Final vertex index
*
*     -----------------------------------------------------
      implicit none

      REAL*4 DIR(3),ABSP,ABSPBST
      INTEGER*4 VERTI,VERTF,IP, I

C#define DEBUG

#include "fsihist.h"

C     Do not add any more intermediate particles if there
C     was not enough space to add the final vertex
C     (and if there's no more space in the particle array)
      if (NVCVERT.ge.MAXVCVERT.or.NVERT.ge.MAXVERT) then
         return
      end if

      NVCVERT=NVCVERT+1

#ifdef DEBUG
      write(*,*) "FSIFILLPART(): NVCVERT = ",NVCVERT
#endif

      DO 1 I=1,3
         DIRVERT(I,NVCVERT)=DIR(I)

#ifdef DEBUG
         write (*,*) "DIRVERT() = ", DIRVERT(I,NVCVERT)
#endif

 1    continue

      ABSPVERT(NVCVERT)=ABSP
      ABSTPVERT(NVCVERT)=ABSPBST
      IPVERT(NVCVERT)=IP
      IVERTI(NVCVERT)=VERTI
      IVERTF(NVCVERT)=VERTF

      VERTI=VERTF

#ifdef DEBUG
      write(*,*) "ABSPVERT  = ",ABSPVERT(NVCVERT)
      write(*,*) "ABSTPVERT = ",ABSTPVERT(NVCVERT)
      write(*,*) "IPVERT    = ",IPVERT(NVCVERT)
      write(*,*) "IVERTI    = ",IVERTI(NVCVERT)
      write(*,*) "IVERTF    = ",IVERTF(NVCVERT)
#endif

      END

