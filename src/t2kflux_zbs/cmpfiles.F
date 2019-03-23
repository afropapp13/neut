C     
C     TEST PROGRAM TO READ VECTOR
C
      IMPLICIT NONE


#include "nework.h"
#include "vcwork.h"
#include "vcvrtx.h"

#include "skhead.h"
C#include "skparm.h"
C#include "sktq.h"

      INTEGER*4 LUNI1,LUNI2
      PARAMETER(LUNI1=10)
      PARAMETER(LUNI2=20)

      character*200 fname_in

      REAL POS(3)

      INTEGER*4  LRECL
      PARAMETER (LRECL=1024)

      INTEGER*4  ISIZE
      PARAMETER (ISIZE=5000000)

      INTEGER*4  IERR,I,N1,N2,ICYCLE,LBANK,IRET

      INTEGER*4  MODE1,MODE2,NPAR1,NPAR2
      REAL*4     ABSP1,ABSP2

      INTEGER*4  lenchr
      external   lenchr

C      integer*4 rlu
C      External  rlu
C      H(1)= rlu(i)

C
C --  INITIALIZATION
C
      CALL KZINIT

      call getarg(1, fname_in) 
      if (fname_in(1:5) .eq. '/disk') then
         call set_rflist(luni1,fname_in,'DISK',' ','RED',' ',' ',
     &        'recl=5670 status=old',' ',' ')
         print *, 'open as DISK: ', fname_in(1:lenchr(fname_in))
      else
         call set_rflist(luni1,fname_in,'LOCAL',' ','RED',' ',' ',
     &        'recl=5670 status=old',' ',' ')
         print *, 'open as LOCAL: ', fname_in(1:lenchr(fname_in))
      endif
      call skopenf(lunI1,1,'Z',iret)
      
      call getarg(2, fname_in) 
      if (fname_in(1:5) .eq. '/disk') then
         call set_rflist(luni2,fname_in,'DISK',' ','RED',' ',' ',
     &        'recl=5670 status=old',' ',' ')
         print *, 'open as DISK: ', fname_in(1:lenchr(fname_in))
      else
         call set_rflist(luni2,fname_in,'LOCAL',' ','RED',' ',' ',
     &        'recl=5670 status=old',' ',' ')
         print *, 'open as LOCAL: ', fname_in(1:lenchr(fname_in))
      endif
      call skopenf(luni2,1,'Z',iret)

C
C --  READ VECTOR
C

      DO 10 I=1,50000
 100     CALL KZECLR
         CALL KZREAD(LUNI1,IERR)
         CALL HEADSK(100)
         IF (IERR.eq.2) GOTO 1000
         IF (IERR.eq.1) GOTO 2000

         CALL KZBLOC('NEUT',LBANK)
         IF (LBANK.eq.0) GOTO 2000
         IF (MOD(I,1000).eq.0) then
             write(6,*)  ' *** ',i,' ***'
         endif

         CALL NERDNEBK(POS)

         MODE1 = MODENE         
         NPAR1 = NUMNE
         ABSP1=
     $        sqrt( PNE(1,1)**2 +
     $        PNE(2,1)**2 +
     $        PNE(3,1)**2)
            
 200     CALL KZECLR
         CALL KZREAD(LUNI2,IERR)
         CALL HEADSK(100)
         IF (IERR.eq.2) GOTO 1100
         IF (IERR.eq.1) GOTO 1200

         CALL KZBLOC('NEUT',LBANK)
         IF (LBANK.eq.0) GOTO 200
         IF (MOD(I,1000).eq.0) then
             write(6,*)  ' *** ',i,' ***'
         endif

         CALL NERDNEBK(POS)

         MODE2 = MODENE         
         NPAR2 = NUMNE
         ABSP2=
     $        sqrt( PNE(1,1)**2 +
     $        PNE(2,1)**2 +
     $        PNE(3,1)**2)
            
         if ( (mode1.ne.mode2).or.
     $        (npar1.ne.npar2).or.
     $        (absp1.ne.absp2)) then
            write(*,*) 'Different Event found at ',NEVSK
            goto 3000
         endif

 10   CONTINUE
      GOTO 3000

 1000 WRITE(6,*) ' READ ERROR : ',LUNI1
      GOTO 3000

 1100 WRITE(6,*) ' READ ERROR : ',LUNI2
      GOTO 3000

 1200 WRITE(6,*) ' END OF FILE : ',LUNI2
      GOTO 3000

 2000 WRITE(6,*) ' END OF FILE : ',LUNI1

 3000 CONTINUE

C      CALL SKCLOSEF(LUNI)

      STOP
      END
