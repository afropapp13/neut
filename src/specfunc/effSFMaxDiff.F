***********************************************************************
*
*  function returns the maximum possible differential cross-section for
*  a given energy, flavour, and target.  event weights are then thrown
*  against the look up values to accept/reject them.
*  
*
***********************************************************************
*
*  input:
*    E - neutrino energy
*    IPAR - neutrino type (pdg)
*    targetType - target type (pdg)
*    nc - true if NC, false if CC
*
*  output:
*   maxDiff - maximum differential cross-section (10^-38 cm^2)
*
*
***********************************************************************
      double precision function effSFMaxDiff(E, ipar, targ, nc, MDLQE,
     $     XMAQE)
        implicit none
#include "nework.h"
        double precision E
        integer ipar
        integer targ
        integer itarg
        logical nc
        integer MDLQE
        real*4  XMAQE
        integer maxbin

c ----- looping index, i, and selected bin
        integer i
        integer bin

c ----- common block for arrays - filled by readData subroutine
c ----- 9 nuclei(C,O,Fe,-,-,-,-,-,-)
c ----- 6 neutrino flavours (all nc are equivalent so not used) (nue,nuebar,numu,numubar,nutau,nutaubar)
c ----- for nc, 2 nucleons (proton, neutron)
        double precision energiescc(6, 210), diffXsecscc(12,6,210)! nucleus, pdg
        double precision energiesnc(200), diffXsecsnc(2,2,12,200)! pdg, nucleon, nucleus

        common / diffXsecArrayscc / energiescc, diffXsecscc 
        common / diffXsecArraysnc / energiesnc, diffXsecsnc 

        integer nuc, neutrino

        integer i_loadDiffXseccc, i_loadDiffXsecnc
        save i_loadDiffXseccc, i_loadDiffXsecnc

        data i_loadDiffXseccc /0/
        data i_loadDiffXsecnc /0/

        if (nc.EQV..TRUE.)then
          if (i_loadDiffXsecnc .EQ. 0) then
            call readEffSFXsecData(ipar, targ, nc, MDLQE, XMAQE)
            i_loadDiffXsecnc = 1
          endif
        else
          if (i_loadDiffXseccc .EQ. 0) then
            call readEffSFXsecData(ipar, targ, nc, MDLQE, XMAQE)
            i_loadDiffXseccc = 1
          endif
        endif


        maxbin = 210

        if (ipar .EQ. 12) then
          neutrino = 1
        else if (ipar .EQ. -12) then
          neutrino = 2
        else if (ipar .EQ. 14) then
          neutrino = 3
        else if (ipar .EQ. -14) then
          neutrino = 4
        else if (ipar .EQ. 16) then
          neutrino = 5
          maxbin = 125
        else if (ipar .EQ. -16) then
          neutrino = 6
          maxbin = 125
        endif

        if (targ .EQ. 1000010020) then
           nuc = 1
        else if (targ .EQ. 1000020030) then
           nuc = 2
        else if (targ .EQ. 1000020040) then
           nuc = 3
        else if (targ .EQ. 1000060120) then
           nuc = 4
        else if (targ .EQ. 1000080160) then
           nuc = 5
        else if (targ .EQ. 1000100200) then
           nuc = 6
        else if (targ .EQ. 1000130270) then
           nuc = 7
        else if (targ .EQ. 1000180400) then
           nuc = 8
        else if (targ .EQ. 1000260560) then
           nuc = 9
        else if (targ .EQ. 1000290630) then
           nuc = 10
        else if (targ .EQ. 1000300640) then
           nuc = 11
        else if (targ .EQ. 1000822080) then
           nuc = 12
        endif

c ----- If using CC
        if (nc .EQV..FALSE.)then
c ----- check for under flow of calculated tables - return 0 as xsec is small and unknown
        if (E .LT. energiescc(neutrino, 1)) then
          write(*,*) "energy lower than minumum calculated - xsec = 0"
          effSFMaxDiff = 0
          return
        endif

c ----- loop through energy array to find correct point
C        do i = 0, 209
        do i = 0, maxbin - 1
          if (energiescc(neutrino, i+1) .GE. E) then
            bin = i
            goto 10
          endif
        enddo
c ----- if sat here - E is bigger than the max allowed - extrapolate linearly
c ----- extrapolation based on bins N and N-2
        bin = 210
c       write (*,*) "Selected energy greater than maximum calculated", E
        effSFMaxDiff = diffXsecscc(nuc, neutrino, bin-2)
     & + (diffXsecscc(nuc, neutrino, bin)-
     & diffXsecscc(nuc, neutrino, bin-2))
     & *( E - energiescc( neutrino, bin-2) ) /
     & ( energiescc(neutrino, bin) - 
     & energiescc(neutrino, bin-2) )
        return

10      continue

c ----- interpolate to correct point in effSFMaxDiff array
        effSFMaxDiff = diffXsecscc(nuc, neutrino, bin)
     & + (diffXsecscc(nuc, neutrino, bin+1)-
     & diffXsecscc(nuc, neutrino, bin))
     & *( E - energiescc( neutrino, bin) ) /
     & ( energiescc(neutrino, bin+1) - 
     & energiescc(neutrino, bin) )

        return

c ----- using NC
        else if (nc .EQV..TRUE.) then
        if (IPNE(2) .EQ. 2212) then
          itarg = 1
        else if (IPNE(2) .EQ. 2112) then
          itarg = 2
        endif

c ----- check for under flow of calculated tables - return 0 as xsec is small and unknown
        if (E .LT. energiesnc(1)) then
          write(*,*) "energy lower than minumum calculated - xsec = 0"
          effSFMaxDiff = 0
          return
        endif

c ----- loop through energy array to find correct point
        do i = 0, 199
          if (energiesnc(i+1) .GE. E) then
            bin = i
            goto 20
          endif
        enddo
c ----- if sat here - E is bigger than the max allowed - extrapolate linearly
c ----- extrapolation based on bins N and N-2
        bin = 200

c       write (*,*) "Selected energy greater than maximum calculated", E
        if ( ipar. gt. 0 ) then 
           neutrino = 1
        else
           neutrino = 2
        endif

        effSFMaxDiff = diffXsecsnc(neutrino,itarg, nuc,  bin-2)
     & + (diffXsecsnc(neutrino,itarg, nuc,  bin)-
     & diffXsecsnc(neutrino,itarg, nuc,  bin-2))
     & *( E - energiesnc( bin-2) ) /
     & ( energiesnc( bin) - 
     & energiesnc( bin-2) )
       return

20      continue


c ----- interpolate to correct point in maxDiff array
        if ( ipar. gt. 0 ) then 
           neutrino = 1
        else
           neutrino = 2
        endif

        effSFMaxDiff = diffXsecsnc(neutrino,itarg, nuc,  bin)
     & + (diffXsecsnc(neutrino,itarg, nuc,  bin+1)-
     & diffXsecsnc(neutrino,itarg, nuc,  bin))
     & *( E - energiesnc( bin) ) /
     & ( energiesnc( bin+1) - 
     & energiesnc( bin) )

        return

        endif

        end
