***********************************************************************
*
*  Subroutine reads xsec data from relevant file into arrays
*  Files are named nuPdg_targPdg_current.csv
*  i.e. numu, C12, ccqe is 14_1000060120_cc.csv
*
*  files are stored under qelSfData/maxDiff
*
***********************************************************************
*
*  input:
*    ipar - neutrino pdg code
*    targ - target pdg code
*    nc - true for nc, false for cc
*
*  output:
*    common block diffXsecArrays - contains:
*      energies - array of energies in table
*      diffXsecs - corresponding maximum differential xsec
*
***********************************************************************

      subroutine readEffSFXsecData(ipar, targ, nc, MDLQE, XMAQE)
      
        implicit none
        
#include "neutfilepath.h"

        integer ipar
        integer targ
        logical nc
        integer idum
        integer MDLQE
        real*4  XMAQE

        integer i

        integer*4      LENSTR
        external       lenstr
        integer*4      lenpath

c ----- setting up the filename to open
        character(len=1) sgn
        character(len=1024) fileName
        character(len=2) interaction

c ----- common block for arrays - filled by readData subroutine
c ----- first 2 indices nucleon(C,O,Fe,-,-,-,-,-,-) IPAR(nue,nuebar,numu,numuba,nutau,nutaubar)
c ----- energies doesn't need a nucleon index (nc don't need neutrino index)
c ----- nc has 2 neutrino types, (neutrino/anti)
c ----- nc array second index is proton(1)/neutron(2)
        double precision energiescc(6, 210), diffXsecscc(12, 6, 210)! nucleus, pdg
        double precision energiesnc(200), diffXsecsnc(2,2,12,200)! pdg, nucleon, nucleus
        common / diffXsecArrayscc / energiescc, diffXsecscc 
        common / diffXsecArraysnc / energiesnc, diffXsecsnc 

        lenpath = lenstr(CRSTBLPATH)

        if (nc .EQV. .TRUE.)then
          interaction = "nc"
        else
          interaction = "cc"
        endif

        if (ipar .LT. 0) then
          sgn = "-"
        else
          sgn = "+"
        endif

c ----- filename:
        if (mod(MDLQE,1000)/100 .EQ. 7) then
C           write (fileName, "(A38,A2,A4)") 
C     $     "qelSfData/temSFMaxDiff/ALL_1000000000_",interaction,".csv"
        write (fileName, "(A10,F4.2,A29,A2,A4)") 
     $          "qelSfData.",XMAQE,
     $          "/temSFMaxDiff/ALL_1000000000_",interaction,".csv"
        else
C           write (fileName, "(A38,A2,A4)") 
C     $     "qelSfData/effSFMaxDiff/ALL_1000000000_",interaction,".csv"
           write (fileName, "(A10,F4.2,A29,A2,A4)") 
     $          "qelSfData.",XMAQE,
     $          "/effSFMaxDiff/ALL_1000000000_",interaction,".csv"
        endif    

          
c ----- open file
        fileName = CRSTBLPATH(1:lenpath)//fileName

        open(99, file=fileName, status='old')
        write(*,*) "reading in effective sf xsec file ", fileName

c ----- read file into arrays
c        read(99,*) nBins
        if (nc.EQV..FALSE.) then
        do i = 1, 1090
          if (i.LE.210) then ! nue
            read(99,*) idum, energiescc(1,i), diffXsecscc(1,1,i),
     &            diffXsecscc(2,1,i), diffXsecscc(3,1,i),
     &            diffXsecscc(4,1,i), diffXsecscc(5,1,i),
     &            diffXsecscc(6,1,i), diffXsecscc(7,1,i),
     &            diffXsecscc(8,1,i), diffXsecscc(9,1,i),
     &            diffXsecscc(10,1,i), diffXsecscc(11,1,i),
     &            diffXsecscc(12,1,i)
          else if (i.GT.210 .AND. i.LE.420) then ! nue-bar
            read(99,*) idum, energiescc(2,i-210),
     &            diffXsecscc(1,2,i-210), diffXsecscc(2,2,i-210),
     &            diffXsecscc(3,2,i-210), diffXsecscc(4,2,i-210),
     &            diffXsecscc(5,2,i-210), diffXsecscc(6,2,i-210),
     &            diffXsecscc(7,2,i-210), diffXsecscc(8,2,i-210),
     &            diffXsecscc(9,2,i-210), diffXsecscc(10,2,i-210),
     &            diffXsecscc(11,2,i-210), diffXsecscc(12,2,i-210)

           else if (i.GT.420 .AND. i.LE.630) then ! numu
            read(99,*) idum, energiescc(3,i-420),
     &             diffXsecscc(1,3,i-420), diffXsecscc(2,3,i-420),
     &             diffXsecscc(3,3,i-420), diffXsecscc(4,3,i-420),
     &             diffXsecscc(5,3,i-420), diffXsecscc(6,3,i-420),
     &             diffXsecscc(7,3,i-420), diffXsecscc(8,3,i-420),
     &             diffXsecscc(9,3,i-420), diffXsecscc(10,3,i-420),
     &             diffXsecscc(11,3,i-420), diffXsecscc(12,3,i-420)

           else if (i.GT.630 .AND. i.LE.840) then ! numu-bar
            read(99,*) idum, energiescc(4,i-630),
     &             diffXsecscc(1,4,i-630), diffXsecscc(2,4,i-630),
     &             diffXsecscc(3,4,i-630), diffXsecscc(4,4,i-630),
     &             diffXsecscc(5,4,i-630), diffXsecscc(6,4,i-630),
     &             diffXsecscc(7,4,i-630), diffXsecscc(8,4,i-630),
     &             diffXsecscc(9,4,i-630), diffXsecscc(10,4,i-630),
     &             diffXsecscc(11,4,i-630), diffXsecscc(12,4,i-630)

           else if (i.GT.840 .AND. i.LE.965) then ! nutau
            read(99,*) idum, energiescc(5,i-840),
     &             diffXsecscc(1,5,i-840), diffXsecscc(2,5,i-840),
     &             diffXsecscc(3,5,i-840), diffXsecscc(4,5,i-840),
     &             diffXsecscc(5,5,i-840), diffXsecscc(6,5,i-840),
     &             diffXsecscc(7,5,i-840), diffXsecscc(8,5,i-840),
     &             diffXsecscc(9,5,i-840), diffXsecscc(10,5,i-840),
     &             diffXsecscc(11,5,i-840), diffXsecscc(12,5,i-840)

           else if (i.GT.965 .AND. i.LE.1090) then ! nutau
            read(99,*) idum, energiescc(6,i-965),
     &             diffXsecscc(1,6,i-965), diffXsecscc(2,6,i-965),
     &             diffXsecscc(3,6,i-965), diffXsecscc(4,6,i-965),
     &             diffXsecscc(5,6,i-965), diffXsecscc(6,6,i-965),
     &             diffXsecscc(7,6,i-965), diffXsecscc(8,6,i-965),
     &             diffXsecscc(9,6,i-965), diffXsecscc(10,6,i-965),
     &             diffXsecscc(11,6,i-965), diffXsecscc(12,6,i-965)

         endif
        enddo ! end loop through cc bits
        else if (nc.EQV..TRUE.) then
          do i = 1, 800
          if (i .LE. 200) then
            read(99,*) idum, energiesnc( i),
     &      diffXsecsnc(1, 1,1,i), ! neutrino, nucleon, nucleus, index(energy)
     &      diffXsecsnc(1, 1,2,i), 
     &      diffXsecsnc(1, 1,3,i), 
     &      diffXsecsnc(1, 1,4,i), 
     &      diffXsecsnc(1, 1,5,i), 
     &      diffXsecsnc(1, 1,6,i), 
     &      diffXsecsnc(1, 1,7,i), 
     &      diffXsecsnc(1, 1,8,i), 
     &      diffXsecsnc(1, 1,9,i), 
     &      diffXsecsnc(1, 1,10,i),
     &      diffXsecsnc(1, 1,11,i),
     &      diffXsecsnc(1, 1,12,i)

          else if (i.GT.200 .AND. i .LE. 400)  then
            read(99,*) idum, energiesnc( i-200),
     &      diffXsecsnc(1, 2,1,i-200),
     &      diffXsecsnc(1, 2,2,i-200),
     &      diffXsecsnc(1, 2,3,i-200),
     &      diffXsecsnc(1, 2,4,i-200),
     &      diffXsecsnc(1, 2,5,i-200),
     &      diffXsecsnc(1, 2,6,i-200),
     &      diffXsecsnc(1, 2,7,i-200),
     &      diffXsecsnc(1, 2,8,i-200),
     &      diffXsecsnc(1, 2,9,i-200),
     &      diffXsecsnc(1, 2,10,i-200),
     &      diffXsecsnc(1, 2,11,i-200),
     &      diffXsecsnc(1, 2,12,i-200)

          else if (i.GT.400 .AND. i .LE. 600)  then
            read(99,*) idum, energiesnc( i-400),
     &      diffXsecsnc(2, 1,1,i-400),
     &      diffXsecsnc(2, 1,2,i-400),
     &      diffXsecsnc(2, 1,3,i-400),
     &      diffXsecsnc(2, 1,4,i-400),
     &      diffXsecsnc(2, 1,5,i-400),
     &      diffXsecsnc(2, 1,6,i-400),
     &      diffXsecsnc(2, 1,7,i-400),
     &      diffXsecsnc(2, 1,8,i-400),
     &      diffXsecsnc(2, 1,9,i-400),
     &      diffXsecsnc(2, 1,10,i-400),
     &      diffXsecsnc(2, 1,11,i-400),
     &      diffXsecsnc(2, 1,12,i-400)

          else if (i.GT.600 .AND. i .LE. 800)  then
            read(99,*) idum, energiesnc( i-600),
     &      diffXsecsnc(2, 2,1,i-600),
     &      diffXsecsnc(2, 2,2,i-600),
     &      diffXsecsnc(2, 2,3,i-600),
     &      diffXsecsnc(2, 2,4,i-600),
     &      diffXsecsnc(2, 2,5,i-600),
     &      diffXsecsnc(2, 2,6,i-600),
     &      diffXsecsnc(2, 2,7,i-600),
     &      diffXsecsnc(2, 2,8,i-600),
     &      diffXsecsnc(2, 2,9,i-600),
     &      diffXsecsnc(2, 2,10,i-600),
     &      diffXsecsnc(2, 2,11,i-600),
     &      diffXsecsnc(2, 2,12,i-600)

          endif
          enddo
        endif ! ending if nc
        
        close(99)
      return
      end
