      program dump_crs_for_osc_main

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

         write(*,10) '# IP','E',
     $        'TOTCRS','OSCCRS'
            
 10         format(A4,A7,A10,A10)

         do 200 J=1,600
            E = 0.05*real(j)-0.025
            call dump_crs_for_osc(ipar,E)
 200     continue

 100  continue

      call nechadir(cost,phi,tdir,tgdir)

      end
      