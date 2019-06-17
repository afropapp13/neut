******************************************************************
***** makeTables.exe builds the NEUT spectral function tables    *
***** tables are built individually for a combination of         *
***** neutrino, nucleus, and cc/nc.                              *
***** tables must then be combined into one for NEUT.            *
***** The value of MAQE must be hardcoded in this file currently *
******************************************************************

#include "effsfevent.F"

      program makeEffSFTables
        implicit none

#include "../neutcore/necard.h"
#include "../neutcore/neutmodel.h"
#include "effspecfunc.h"
#include "nework.h"

c ----- neutrino type, energy
        integer IPAR
        real E
        real DIRNEU(3)
        integer IERR, mode
        logical nc, single, first

        integer i

        double precision norm
        integer nBins
        double precision electronEnergy(210)
        double precision muonEnergy(210)
        double precision tauEnergy(125)
        
        double precision mean, maxValue, meanNC(2), maxValueNC(2)
        integer N_tries
        common /makeEffSFTable/ mean,maxValue,meanNC,maxValueNC, N_tries

        integer targ
        character(len=1) sgn
        character(len=2) interaction
        character(len=44) maxFile
        character(len=44) totFile

        character(10) argument
        character(10) argMDLQE


********************************************************************
c ----- initialise parameters in sfevent required for table building
        mean = 0
        maxValue = 0
        N_tries = 0
        meanNC(1) = 0
        meanNC(2) = 0
        maxValueNC(1) = 0
        maxValueNC(2) = 0
        DIRNEU(1) = 0
        DIRNEU(2) = 0
        DIRNEU(3) = 1
        single = .FALSE.
        first = .TRUE.
        ierr = 0

c       need to set FPQE and PF too - otherwise default is zero
        FPQE = 1
        PFSF = 209
*********************************************************************

c ----- set parameters for interaction - neutrino, target, interaction
c ----- these are the defaults - others read in on the command line
        IPAR = 14
        NUMBNDN = 6
        NUMBNDP = 6
        NUMATOM = 12

******* set value of MAQE to be used *******
        XMAQE = 1.21

******* set which form factors to use ******
        MDLQE = 602

******* set scc if desired *******
c        SCCFV = 1
c        SCCFA = 1

c ----- if commandline arguments given then read them in
c ----- run makeTables.exe IPAR(PDG) Nucleus(e.g C12) NC/CC
        if (iargc() .GT. 0) then
          call GETARG(1, argument)
c ------- neutrino pdg first          
          if (argument .EQ. "14") then
            IPAR = 14
          else if (argument .EQ. "-14") then
            IPAR = -14
          else if (argument .EQ. "12") then
            IPAR = 12
          else if (argument .EQ. "12") then
            IPAR = 12
          else if (argument .EQ. "-12") then
            IPAR = -12
          else if (argument .EQ. "16") then
            IPAR = 16
          else if (argument .EQ. "-16") then
            IPAR = -16
          else if (argument .EQ. "combine") then
c           option to combine all pre-made tables into one file for use.
            call GETARG(2, argument)
            if ((argument .EQ. "CC") .OR. (argument .EQ. "NC")) then
               call GETARG(3, argMDLQE)
               if (argMDLQE .EQ. "TEM") then
                  MDLQE = 702
               endif
              call combineEffSFTables(argument, MDLQE)
              stop
            else
              write (*,*) "When combining tables, a second argument
     &, CC or NC should be given"
            stop
            endif
          else
            write (*,*) "INVALID ARGUMENTS 1, MUST BE IPAR (pdg code)
     &, NUCLEUS (e.g. C12, Fe56), NC/CC, or combine, NC/CC"
            stop
          endif
c ------- now check nucleus
          call GETARG(2, argument)

          if (argument .EQ. "H2") then
            NUMBNDN = 1
            NUMBNDP = 1
            NUMATOM = 2
          else if (argument .EQ. "He3") then
            NUMBNDN = 1
            NUMBNDP = 2
            NUMATOM = 3
         else if (argument .EQ. "He4") then
            NUMBNDN = 2
            NUMBNDP = 2
            NUMATOM = 4
         else if (argument .EQ. "C12") then
            NUMBNDN = 6
            NUMBNDP = 6
            NUMATOM = 12
         else if (argument .EQ. "O16") then
            NUMBNDN = 8
            NUMBNDP = 8
            NUMATOM = 16
         else if (argument .EQ. "Ne20") then
            NUMBNDN = 10
            NUMBNDP = 10
            NUMATOM = 20
         else if (argument .EQ. "Al27") then
            NUMBNDN = 14
            NUMBNDP = 13
            NUMATOM = 27
         else if (argument .EQ. "Ar40") then
            NUMBNDN = 22
            NUMBNDP = 18
            NUMATOM = 40
         else if (argument .EQ. "Fe56") then
            NUMBNDN = 30
            NUMBNDP = 26
            NUMATOM = 56
         else if (argument .EQ. "Cu63") then
            NUMBNDN = 34
            NUMBNDP = 29
            NUMATOM = 63
         else if (argument .EQ. "Zn64") then
            NUMBNDN = 34
            NUMBNDP = 30
            NUMATOM = 64
         else if (argument .EQ. "Pb208") then
            NUMBNDN = 126
            NUMBNDP = 82
            NUMATOM = 208
          else
            write (*,*) "INVALID ARGUMENTS, MUST BE IPAR (pdg code)
     &, NUCLEUS (e.g. C12, Fe56), NC/CC, or combine, CC/NC"
            stop
         endif
c ------- finally, nc
          call GETARG(3, argument)
          if (argument .EQ. "NC") then
            nc = .TRUE.
          else if (argument .EQ. "CC") then
            nc = .FALSE.
          else
            write (*,*) "INVALID ARGUMENTS, MUST BE IPAR (pdg code)
     &, NUCLEUS (e.g. C12, Fe56), NC/CC, or combine, CC, NC"
            stop
          endif
          call GETARG(4, argument)
          if (argument .EQ. "TEM") then
             MDLQE = 702
          endif
        endif
          

        write (*,*) "creating effective SF maxDiff table"
        write (*,*) "creating effective SF total xsec table"

c ----- set true for nc tables
        if (nc .EQV..TRUE.) then
          mode = 51
        else
          mode = 1
        endif
**********************************************************************
********************************************************************************
c ----- set variables in filename - neutrino type, interaction type, target code
        if (IPAR .LT. 0) then
          sgn = "-"
        else
          sgn = "+"
        endif

        if (nc .EQV. .TRUE.) then
          interaction = "nc"
        else
          interaction = "cc"
        endif

        targ = 1000000000 + NUMBNDP*10000 + NUMATOM*10

c ----- create output fileNames
        if (mod(MDLQE,1000)/100 .EQ. 7) then
           write(maxFile, "(A23,A1, I2, A1, I10, A1, A2, A4)")
     &          "qelSfData/temSFMaxDiff/",sgn,abs(ipar),
     &          "_",targ,"_",interaction,".csv" 
           
           write (totFile, "(A23,A1, I2, A1, I10, A1, A2, A4)")
     &          "qelSfData/temSFTotXsec/",sgn,abs(ipar),
     &          "_",targ,"_",interaction,".csv"
        else
           write(maxFile, "(A23,A1, I2, A1, I10, A1, A2, A4)")
     &          "qelSfData/effSFMaxDiff/",sgn,abs(ipar),
     &          "_",targ,"_",interaction,".csv" 
           
           write (totFile, "(A23,A1, I2, A1, I10, A1, A2, A4)")
     &          "qelSfData/effSFTotXsec/",sgn,abs(ipar),
     &          "_",targ,"_",interaction,".csv"
        endif
        
        write (*,"(A1, I2, A1, I10, A1, A2)")
     &  sgn,abs(ipar),"_",targ,"_",interaction

c ----- open files     
        open(27, file=maxFile,status="new")
        open(35, file=totFile,status="new")

**************************************************************************
* define arrays of values energy for calculations to be performed
*
* electron, muon and tau can potentially use different binning
* currently electron and muon are the same, tau is different
**************************************************************************
c ----- energies given in MeV - sfevent takes them in GeV

      data electronEnergy /   25,  75, 125, 175, 225, 275,
     &                       325, 375, 425, 475, 525, 575,
     &                       625, 675, 725, 775, 825, 875,
     &                       925, 975,1025,1075,1125,1175,
     &                      1225,1275,1325,1375,1425,1475,
     &                      1525,1575,1625,1675,1725,1775,
     &                      1825,1875,1925,1975,2025,2075,
     &                      2125,2175,2225,2275,2325,2375,
     &                      2425,2475,2525,2575,2625,2675,
     &                      2725,2775,2825,2875,2925,2975,
     &                      3025,3075,3125,3175,3225,3275,
     &                      3325,3375,3425,3475,3525,3575,
     &                      3625,3675,3725,3775,3825,3875,
     &                      3925,3975,4025,4075,4125,4175,
     &                      4225,4275,4325,4375,4425,4475,
     &                      4525,4575,4625,4675,4725,4775,
     &                      4825,4875,4925,4975,5025,5075,
     &                      5125,5175,5225,5275,5325,5375,
     &                      5425,5475,5525,5575,5625,5675,
     &                      5725,5775,5825,5875,5925,5975,
     &                      6025,6075,6125,6175,6225,6275,
     &                      6325,6375,6425,6475,6525,6575,
     &                      6625,6675,6725,6775,6825,6875,
     &                      6925,6975,7025,7075,7125,7175,
     &                      7225,7275,7325,7375,7425,7475,
     &                      7525,7575,7625,7675,7725,7775,
     &                      7825,7875,7925,7975,8025,8075,
     &                      8125,8175,8225,8275,8325,8375,
     &                      8425,8475,8525,8575,8625,8675,
     &                      8725,8775,8825,8875,8925,8975,
     &                      9025,9075,9125,9175,9225,9275,
     &                      9325,9375,9425,9475,9525,9575,
     &                      9625,9675,9725,9775,9825,9875,
     &                      9925,9975,
     &                      11000,12000,13000,14000,15000,
     &                      16000,17000,18000,19000,20000/

      data muonEnergy /      25,  75, 125, 175, 225, 275,
     &                       325, 375, 425, 475, 525, 575,
     &                       625, 675, 725, 775, 825, 875,
     &                       925, 975,1025,1075,1125,1175,
     &                      1225,1275,1325,1375,1425,1475,
     &                      1525,1575,1625,1675,1725,1775,
     &                      1825,1875,1925,1975,2025,2075,
     &                      2125,2175,2225,2275,2325,2375,
     &                      2425,2475,2525,2575,2625,2675,
     &                      2725,2775,2825,2875,2925,2975,
     &                      3025,3075,3125,3175,3225,3275,
     &                      3325,3375,3425,3475,3525,3575,
     &                      3625,3675,3725,3775,3825,3875,
     &                      3925,3975,4025,4075,4125,4175,
     &                      4225,4275,4325,4375,4425,4475,
     &                      4525,4575,4625,4675,4725,4775,
     &                      4825,4875,4925,4975,5025,5075,
     &                      5125,5175,5225,5275,5325,5375,
     &                      5425,5475,5525,5575,5625,5675,
     &                      5725,5775,5825,5875,5925,5975,
     &                      6025,6075,6125,6175,6225,6275,
     &                      6325,6375,6425,6475,6525,6575,
     &                      6625,6675,6725,6775,6825,6875,
     &                      6925,6975,7025,7075,7125,7175,
     &                      7225,7275,7325,7375,7425,7475,
     &                      7525,7575,7625,7675,7725,7775,
     &                      7825,7875,7925,7975,8025,8075,
     &                      8125,8175,8225,8275,8325,8375,
     &                      8425,8475,8525,8575,8625,8675,
     &                      8725,8775,8825,8875,8925,8975,
     &                      9025,9075,9125,9175,9225,9275,
     &                      9325,9375,9425,9475,9525,9575,
     &                      9625,9675,9725,9775,9825,9875,
     &                      9925,9975,
     &                      11000,12000,13000,14000,15000,
     &                      16000,17000,18000,19000,20000/

      data tauEnergy /  2700, 2800, 2900, 3000, 3100, 3200,
     &                  3300, 3400, 3500, 3600, 3700, 3800,
     &                  3900, 4000, 4100, 4200, 4300, 4400,
     &                  4500, 4600, 4700, 4800, 4900, 5000,
     &                  5100, 5200, 5300, 5400, 5500, 5600,
     &                  5700, 5800, 5900, 6000, 6100, 6200,
     &                  6300, 6400, 6500, 6600, 6700, 6800,
     &                  6900, 7000, 7100, 7200, 7300, 7400,
     &                  7500, 7600, 7700, 7800, 7900, 8000,
     &                  8100, 8200, 8300, 8400, 8500, 8600,
     &                  8700, 8800, 8900, 9000, 9100, 9200,
     &                  9300, 9400, 9500, 9600, 9700, 9800,
     &                  9900,10000,10000,11000,12000,13000,
     &                 14000,15000,16000,17000,18000,19000,
     &                 20000,21000,22000,23000,24000,25000,
     &                 26000,27000,28000,29000,30000,31000,
     &                 32000,33000,34000,35000,36000,37000,
     &                 38000,39000,40000,41000,42000,43000,
     &                 44000,45000,46000,47000,48000,49000,
     &                 50000,51000,52000,53000,54000,55000,
     &                 56000,57000,58000,59000,60000/

***************************************************************************
* Calculate values and fill tables here
***************************************************************************

c ----- number of values used differs between neutrinos
        if (abs(IPAR) .EQ. 16) then
          nbins = 125
        elseif (abs(IPAR) .EQ. 14) then
          nbins = 210
        else
          nbins = 210
        endif

c ----- neutral current tables all have 200 entries
        if (nc.EQV..TRUE.)then
          nbins = 200
        endif

c ----- write number of bins
        write (27, "(I3)") nBins
        write (35, "(I3)") nBins

c ----- fill table with values for reading
        do i = 1, nBins

c ----- pick array for correct neutrino
          if (abs(IPAR) .EQ. 16) then
            E = tauEnergy(i)/1000
          elseif (abs(IPAR) .EQ. 14) then
            E = muonEnergy(i)/1000
          else
            E = electronEnergy(i)/1000
          endif
            
          if (nc .EQV..TRUE.) then
            E = i*0.05 - 0.025
          endif

c ----- norm - hold values in a nicer unit 10^-38 cm^2
          norm = 1E38

c ----- run calculation
          write(*,*) "calculating for energy", int(E*1000), "MeV"
          call effsfevent(ipar, e, dirneu, ierr, single, mode, first)
          if (nc .EQV. .FALSE.) then
            write (27, "(F10.4, F12.6)") dble(e), maxValue*norm
            write (35, "(F10.4, F12.6)") dble(e), mean*norm
          else
             write(*,*) dble(e), maxValueNC(1)*
     & norm, maxValueNC(2)*norm
             write (27, "(F10.4, F12.6, F12.6)") dble(e), maxValueNC(1)*
     & norm, maxValueNC(2)*norm
            write (35, "(F10.4, F12.6, F12.6)") dble(e), meanNC(1)*
     & norm, meanNC(2)*norm
          endif
c ----- reset values for next pass          
          mean = 0
          maxValue = 0
          meanNC(1) = 0
          meanNC(2) = 0
          maxValueNC(1) = 0
          maxValueNC(2) = 0
          N_tries = 0
        enddo
      
      stop
      end
