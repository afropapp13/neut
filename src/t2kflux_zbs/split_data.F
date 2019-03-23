C     
C     TEST PROGRAM TO READ VECTOR
C
      IMPLICIT NONE


      INTEGER*4 LUNI1
      PARAMETER(LUNI1=10)
      INTEGER*4 LUNO2
      PARAMETER(LUNO2=15)

      INTEGER*4  LRECL
      PARAMETER (LRECL=1024)

      INTEGER*4  ISIZE
      PARAMETER (ISIZE=5000000)

      REAL*4     H(ISIZE)
      COMMON/PAWC/H

      INTEGER*4  I,LBANK,IERR,IEVENT,num,ifile

      character*100 char_num,in_fname,out_fname,suffix

      integer*4 atoi
      external  atoi

C
C Initialize 
C
      CALL KZINIT
C
C Open Files
C
      
      call getarg(1,in_fname)
      call getarg(2,char_num)
      num = atoi(char_num)
      write(*,*) 'Now opening input file',in_fname
      write(*,*) 'Split input file : Each of the new file has ',
     $           num, 'events.'
      call fort_fopen(LUNI1,in_fname,'RED','Z')

      ifile = 0
C
C --  READ VECTOR
C

      DO 10 I=1,500000

         if (mod(i,num).eq.1) then
            CALL SKCLOSEF(LUNO2)
            write(suffix,'(A1,I3.3)') '.',ifile
C            write(*,*) suffix
            call concat(in_fname,suffix,out_fname)
            write(*,*) 'Now opening output file',out_fname
           call fort_fopen(LUNO2,out_fname,'WRT','Z')
           ifile = ifile +1 
         endif

         CALL KZECLR

C-- read 1st file
         CALL KZREAD(LUNI1,IERR)
         if (IERR.ne.0) then 
            goto 1000
         endif

         call kzwrit(luno2)

 10   CONTINUE
      GOTO 3000

 1000 WRITE(6,*) ' READ ERROR : '
      GOTO 3000

 2000 WRITE(6,*) ' END OF FILE : '

 3000 CONTINUE

      CALL SKCLOSEF(LUNI1)
      CALL SKCLOSEF(LUNO2)

      STOP
      END
