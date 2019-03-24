      program dumpcrsmain

      implicit none
      
      integer*4 IPARTBL(4)
      data ipartbl/14,-14,12,-12/

      integer*4 ipar
      
      integer*4 I,J
      real*4    E

      real*4 cost,phi,tdir(3),tgdir(3)

      call necard

      do 100 I=1,4

         ipar = ipartbl(i)

         if (ipar.gt.0) then
            write(*,10) '# IP','E',
     $           '----','CCQE',
     $           'CCSPI1','CCSPI2','CCSPI3','CCCOH',
     $           'CCMPI1','CCMPI2','CCSETA','CCSK','CCSGAM',
     $           'NCELF','NCELP','NCELN',
     $           'NCSPI1','NCSPI2','NCSPI3','NCSPI4','NCCOH',
     $           'NCMPI1','NCMPI2',
     $           'NCSETA1','NCSETA2','NCSK1','NCSK2',
     $           'NCSGAM','NCSGAM2'
            
 10         format(A4,A8,A8,A8,A8,A8,A8,A8,A8,A12,A8,A8,A8,A8,
     $             A8,A8,A8,A8,A8,A8,A8,A8,A12,A8,A8,A8,A8,A8,A8)
         else
            write(*,20) '# IP','E',
     $           'CCQEF','CCQEP',
     $           'CCSPI1','CCSPI2','CCSPI3','CCCOH',
     $           'CCMPI1','CCMPI2','CCSETA','CCSK','CCSGAM',
     $           'NCELF','NCELP','NCELN',
     $           'NCSPI1','NCSPI2','NCSPI3','NCSPI4','NCCOH',
     $           'NCMPI1','NCMPI2',
     $           'NCSETA1','NCSETA2','NCSK1','NCSK2',
     $           'NCSGAM','NCSGAM2'
            
 20         format(A4,A8,A8,A8,A8,A8,A8,A8,A8,A12,A8,A8,A8,A8,
     $             A8,A8,A8,A8,A8,A8,A8,A8,A12,A8,A8,A8,A8,A8,A8)
         endif
            
         do 200 J=0,100
            E = 0.01*real(j)
            call dumpcrs(ipar,E)
 200     continue
         do 300 J=1,90
            E = 0.1*real(j)+1.
            call dumpcrs(ipar,E)
 300     continue
         do 400 J=1,989
            E = 1*real(j) + 10.
            call dumpcrs(ipar,E)
 400     continue
 100  continue

      call nechadir(cost,phi,tdir,tgdir)

      end
      
      subroutine rmarin(i,j,k)
      write(*,*) 'custom rmarin'
      
      return
      end
      
