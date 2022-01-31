      program main

C     An adaptive integrator using the QUADPACK lib
C     Does the integral over phi, theta, q2, w using dqag
C     dqag uses a Gauss-Kronrod technique, see quadpack.f
C     for more info and
C     https://people.sc.fsu.edu/~jburkardt/f77_src
C                             /quadpack/quadpack.html
C
C     Author: Clarence Wret (cwret14@imperial.ac.uk,
C     clarence.wret@gmail.com)
C     N.B. could use CERNLIB integration routines too

      Implicit None

#include "neutmodel.h"
#include "rscons.h"
#include "necard.h"

C     Modified Rein-Seghal model by Wroclaw group
C     Lepton mass effects are taken into account
C     following the formalism by Ch.Berger and L.M.Sehgal
C     Phys. Rev. D 76, 113004 (2007)
C     Phys. Rev. D 77, 059901(E) (2008) !! Erratum !!

      integer imode, imode_com
      common imode_com
      integer iparlop

      integer iptbl(3)
      data iptbl /12,14,16/

      integer iq2
      integer nq2
      parameter(nq2=100)

      character*15 fname_ref(3)
      character*100 fname
      character*50 argv1
      character*50 argv2

      data fname_ref/
     $     'spi_nue_xsec_ma',
     $     'spi_num_xsec_ma',
     $     'spi_nut_xsec_ma'/

      character*25 typesn(7),typesa(7)
      data typesn/
     $     "nu p -> mu- p pi+",
     $     "nu n -> mu- p pi0",
     $     "nu n -> mu- n pi+",
     $     "nu p -> nu- p pi0",
     $     "nu p -> nu- n pi+",
     $     "nu n -> nu- n pi0",
     $     "nu n -> nu- p pi-"/
      data typesa/
     $     "nubar p -> mu+   n pi-",
     $     "nubar p -> mu+   n pi0",
     $     "nubar p -> mu+   p pi-",
     $     "nubar n -> nubar n pi0",
     $     "nubar n -> nubar p pi-",
     $     "nubar p -> nubar p pi0",
     $     "nubar p -> nubar n pi+"/
      character*30 ffn(2)
      data ffn/
     $     "Rein-Seghal default",
     $     "Graczyk & Sobczyk default"/

C     The iterator for the energy
      double precision dq2, q2_max, q2_min, q2
C     The neutrino energy, cross-section, delta cross-section,
C     user-specified maximum hadronic mass, minimum hadronic mass
      double precision E, sigma, wwmax
C     Mass of lepton (depends on user input)
      double precision xmlep, smax

C     The iterators for initially setting the arrays to zero
      integer i,j,k

C     We have neutrino and anti-neutrino
      integer nutypes
      parameter(nutypes=2)

C     And 7 modes for each (3 CC, 4 NC) -- see neutcore/rsdcrs.F for
C     details
      integer nmodes
      parameter(nmodes=7)

C     xsec(number of lepton generations, number of anti-part/part, # of
C     interaction modes, number of energies)
C     Xsec = x-sec, Xsec3 = x-sec from delta, Enu = neutrino E
c     real             Xsec(3,nutypes,nmodes,iemax),
c    &                 Xsec3(3,nutypes,nmodes,iemax), 
c    &                 Enu(3,nutypes,nmodes,iemax)

C     Have a "handy" common block to make sure all functions know what
C     Enu and W we're on (because we're integrating with a function)
      double precision enu_com, w_com, q2_com
      common /kine/ enu_com, w_com, q2_com

C     Another common block with the rscons.h but with double precision
      double precision xmn_con, xmlep_con, xmpi_con, w_max_con
      common /const/ xmn_con, xmlep_con, xmpi_con, w_max_con

C     Interaction flag, neutrino/anti-neutrino flag and the common block
C     which_com refers to if we're doing the total cross-section
C     integral or the Delta only integral. It's a bit of a hack now (do
C     two integration loops) but I can't bother re-writing quadpack...
      integer intflag, iorg, ipflag, iflag
      integer intflag_com, iorg_com, ipflag_com, which_com
      common /int/ intflag_com, iorg_com, ipflag_com, which_com

C     Include lepton mass-effects or not (ieffect = 0: ON, ieffect = 1:
C     OFF), gets written to iorg depending on lepton and energy
      integer ieffect
      parameter(ieffect=0)

C     The maximum value of W (hadronic mass)
      WWMAX=2.0
C     Set the necard.F settings here
C     Choose Minoo's model
      MDLSPI=1
C     Choose MA
      XMANFFRES=0.95
      XMARSRES = XMANFFRES
C     Choose MV
      XMVNFFRES=0.84
      XMVRSRES = XMVNFFRES
C     Choose CA5
      RNECA5I=1.01
C     Choose I12 background
      RNEBGSCL  = 1.30
C     Form factor models
C     Rein-Seghal: 0, Graczyk & Sobcyzk (default):1
      IFFSPI = 1

C     NRTYPESPI only affects RS model
C     Equations 12, 13 in Graczyk & Sobcyzk (Phys Rev D 77 053001)
C     Electroproduction: 0, Neutrino-production (default):1
      NRTYPESPI = 1

      call nefillmodel

C     Print the settings
      if (iargc().ne.1 .and. iargc().ne.4 .and. iargc().ne.5) then
         print *, 'Cross-section table calculator for'
         print *, 'Rein-Sehgal single pion model (adaptive)'
         print *, '*************************'
         print *, 'iargc is ', iargc()
         print *, 'args are either:'
         print *, 'nu type(1=e,2=m,3=t), MARES, CA50, I_12, WCUT'
         print *, 'nu type(1=e,2=m,3=t), MARES, CA50, I_12'
         print *, 'nu type(1=e,2=m,3=t)'
         print *, '*************************'
         print *, 'If unspecified,'
         print *, 'MARES =', XMASPI
         print *, 'MVRES =', XMVSPI
         print *, 'CA50  =', RCA5ISPI
         print *, 'I12   =', RBGSCLSPI
         print *, 'IFFSPI=', IFFSPI
         print *, 'NRTYPE=', NRTYPESPI
         print *, 'MASSEF=', IEFFECT
         print *, 'WCUT  =', WWMAX
         stop
      endif

C     Set M_A^RES from user input
      if (iargc().eq.4 .or. iargc().eq.5) then
        call getarg(2, argv2)
        read(argv2, '(F4.2)') XMASPI
      endif


C     Set CA5 from user input
      if (iargc().eq.4 .or. iargc().eq.5) then
        call getarg(3, argv2)
        read(argv2, '(F4.2)') RCA5ISPI
      endif

C     Set background scaling
      if (iargc().eq.4 .or. iargc().eq.5) then
        call getarg(4, argv2)
        read(argv2, '(F4.2)') RBGSCLSPI
      endif

C     WMAX input (W<2.0 GeV for "proper" R-S)
      if (iargc() .eq. 5) then
        call getarg(5, argv2)
        read(argv2, '(F4.2)') WWMAX
      endif

      print *, 'MARES =', XMASPI
      print *, 'MVRES =', XMVSPI
      print *, 'CA50  =', RCA5ISPI
      print *, 'I12   =', RBGSCLSPI
      print *, 'IFFSPI=', IFFSPI
      print *, 'NRTYPE=', NRTYPESPI
      print *, 'MASSEF=', IEFFECT
      print *, 'WCUT  =', WWMAX

C     Copy to the common block variable
      w_max_con = WWMAX

      print*, '==============================='
      if (ieffect.eq.0) then
         print*, '== LEPTON MASS EFFECT -> ON ==='
      elseif (ieffect.eq.1) then
         print*, '== LEPTON MASS EFFECT -> OFF =='
      else
         stop' please set LEPTONMASS'
      endif
      print*, '==============================='

C     Form factor variables
      XMVSPI = 0.84

      call getarg(1, argv1)
      read(argv1, '(I10)') IPARLOP
      print *, 'IPARLOP=', IPARLOP

      if ((IPARLOP.GT.3).or.(IPARLOP.lt.1)) then
        write(*,*) 'IPARLOP=', IPARLOP, ' invalid neutrino flavour'
        write(*,*) '1 = nue, 2 = numu, 3 = nutau'
        stop
      endif

      XMLEP = 0.


C     IPFLAG set neutrino or anti-neutrino
      DO 110  IPFLAG = 0,nutypes-1

C       Set the common block holder (needed for Lambda in xsec calc)
        ipflag_com = ipflag

C       IFLAG is looping what interaction we're doing (listed above in
C       data, e.g. CC1pi+1p)
         DO 105  IFLAG = 1,nmodes

            INTFLAG = IFLAG+IPFLAG*10
            intflag_com = intflag

            iOrg = ieffect ! Consider lepton mass effecs(0) or not(1)
            IF (IPARLOP.eq.1)  XMLEP=XME
            IF (IPARLOP.eq.2)  XMLEP=XMMU
            IF (IPARLOP.eq.3)  XMLEP=XMTAU
            IF (IFLAG.gt.3)    then
              iOrg = 1 ! Don't use lepton mass-effect for NC modes
              XMLEP=0.0E0
            ENDIF

C           Convert from REAL to double precision
C           Because NEUT insists on REAL rather than double...
C           whereas QUADPACK insists on double precision rather than
C           REAL, ha!
            xmn_con   = dble(XMN)
            xmlep_con = dble(XMLEP)
            xmpi_con  = dble(XMPI)

C           Write what mode we're doing to the user
            if (IPFLAG .EQ. 0) THEN
              write(*,*)  "#",TYPESn(IFLAG)
            else 
              write(*,*)  "#",TYPESa(IFLAG)
            endif

            E = 1.0
            enu_com = e
            q2_min = 0
            q2_max = 1
            dq2 = (q2_max-q2_min)/Nq2
            DO 1337 iq2=0,nq2
            q2 = q2_min+Iq2*Dq2
            q2_com = q2

C     Lepton mass effects are negligible in higher energy cases
            if (IPARLOP.eq.1.and.E.gt.1.5 ) iOrg = 1
            if (IPARLOP.eq.2.and.E.gt.100.) iOrg = 1
            if (IPARLOP.eq.3.and.E.gt.200.) iOrg = 1
            iorg_com = iorg

            call int_all(intflag, SIGMA)
            write(*,*) q2, sigma
 1337   CONTINUE

 105    CONTINUE ! Finish mode loop
 110  CONTINUE ! Finish lepton flavour loop

C ************************************************************
C     Now write the output
C ************************************************************


C      do I=iparlop,iparlop
C
C      if (WWMAX .ne. 2.0) then
C          write(fname,'(A,F4.2,A,F4.2,A,F4.2,A,F4.2,A)') 
C     $     fname_ref(I), xmaspi, '_ca5I', RCA5ISPI, '_bgscl', 
C     $        RBGSCLSPI, '_W', WWMAX,'.dat'
C          write(*,*) 'Opening file ', fname
C      else 
C          write(fname,'(A,F4.2,A,F4.2,A,F4.2,A)') 
C     $     fname_ref(I), xmaspi, '_ca5I', RCA5ISPI, '_bgscl', 
C     $        RBGSCLSPI, '.dat'
C          write(*,*) 'Opening file ', fname
C      endif
C
C          open(60,FILE=FNAME,ERR=9999)
C
C          write(60,'(A,A,A,F4.2,A,F4.2,A,F4.2,A,I1,A,F4.2,A,F4.2)') 
C     $        '# 1pi :', ffn(IFFSPI+1), ':MA=', XMAspi, 
C     $        ', CA5(0)=', RCA5ISPI, ', BGscale=', RBGSCLSPI, 
C     $        ', NRTYPESPI=', NRTYPESPI,
C     $        ', WMIN=', XMN+XMPI, ', WMAX=',WWMAX
C
C          do IE=1,IEMAX
C          write(60,'(I4,I2,F8.2,7E13.6)')
C     $           IPTBL(I),1,Enu(I,1,1,IE),
C     $           XSec(I,1,1,IE), Xsec(I,1,2,IE), Xsec(I,1,3,IE),
C     $           Xsec(I,1,4,IE), Xsec(I,1,5,IE), Xsec(I,1,6,IE),
C     $           Xsec(I,1,7,IE)
C          enddo
C          do IE=1,IEMAX
C          write(60,'(I4,I2,F8.2,7E13.6)')
C     $           IPTBL(I),2,Enu(I,1,1,IE),
C     $           Xsec3(I,1,1,IE), Xsec3(I,1,2,IE), Xsec3(I,1,3,IE),
C     $           Xsec3(I,1,4,IE), Xsec3(I,1,5,IE), Xsec3(I,1,6,IE),
C     $           Xsec3(I,1,7,IE)
C          enddo
C          do IE=1,IEMAX
C          write(60,'(I4,I2,F8.2,7E13.6)')
C     $        -1*IPTBL(I),1,Enu(I,1,1,IE),
C     $           Xsec(I,2,1,IE), Xsec(I,2,2,IE), Xsec(I,2,3,IE),
C     $           Xsec(I,2,4,IE), Xsec(I,2,5,IE), Xsec(I,2,6,IE),
C     $           Xsec(I,2,7,IE)
C          enddo
C          do IE=1,IEMAX
C          write(60,'(I4,I2,F8.2,7E13.6)')
C     $           -1*IPTBL(I),2,Enu(I,1,1,IE),
C     $           Xsec3(I,2,1,IE), Xsec3(I,2,2,IE), Xsec3(I,2,3,IE),
C     $           Xsec3(I,2,4,IE), Xsec3(I,2,5,IE), Xsec3(I,2,6,IE),
C     $           Xsec3(I,2,7,IE)
C
C          enddo
C          enddo
C          goto 1000

c9999     CONTINUE
c         write(*,*) "Error in opening file."
c         STOP

c1000     continue
c         write(*,*) "DONE!"
          END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
C     a subroutine to perform the integration over all the variables
      subroutine int_all(intflag, SIGMA)

      integer intflag
      double precision smax, sigma, w
      double precision w_min, w_max

      integer limit_w
      parameter (limit_w = 1000)
      integer lenw_w
      parameter (lenw_w = limit_w * 4)

      external f_w

      double precision epsabs_w, epsrel_w, abserr_w
      double precision result_w
      double precision work_w(lenw_w)

      integer neval_w, ier_w, iwork_w(limit_w), key_w, last_w

C     Have a "handy" common block to make sure all functions know what
C     Q2, W, Enu, theta, phi we're currently on
      double precision enu_com, w_com, q2_com
      common /kine/ enu_com, w_com, q2_com

      double precision xmn_con, xmlep_con, xmpi_con, w_max_con
      common /const/ xmn_con, xmlep_con, xmpi_con, w_max_con

      result_w = 0.
      epsabs_w = 0.001
      epsrel_w = 0.001
      key_w = 1
      
C     The maximum center of mass energy
      smax = 2.*enu_com*XMN_con + XMN_con**2

      w_min = XMN_con + XMPI_con
      w_max = MIN(SQRT(SMAX) - XMLEP_con, w_max_con)

C     This does the integral over W
      call dqag (f_w, w_min, w_max, epsabs_w, epsrel_w, 
     & key_w, result_w, abserr_w, neval_w, ier_w, limit_w, 
     & lenw_w, last_w, iwork_w, work_w)

      SIGMA = result_w

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function f_w(w)

C     Q2, returned dxsec sums
      double precision w, dxsec, dxsec3

      integer intflag_com, iorg_com, ipflag_com, which_com
      common /int/ intflag_com, iorg_com, ipflag_com, which_com

C     Another common block with the rscons.h stuff in
      double precision xmn_con, xmlep_con, xmpi_con
      common /const/ xmn_con, xmlep_con, xmpi_con

      double precision enu_com, w_com, q2_com
      common /kine/ enu_com, w_com, q2_com

C     The real counterparts
      real*4 enu_r, q2_r, w_r, xmlep_r, dxsec_r(2), dxsec3_r(2)
      real*4 smax_r, q2_min_r, q2_max_r

      integer*4 lambda(2)
      integer*4 ipflag_r, iorg_r, intflag_r

      dxsec = 0.
      dxsec3 = 0.
      dxsec_r(1) = 0.
      dxsec_r(2) = 0.
      dxsec3_r(1) = 0.
      dxsec3_r(2) = 0.

C     NEUT needs REAL as argument...
      enu_r = enu_com
      q2_r  = q2_com
      w_com = w
      w_r   = w_com
      xmlep_r = xmlep_con

C     Now that we have an enu and w we also have the q2 limits
C     The return is q2_min and q2_max
      smax_r = (2.*XMN_con*enu_com) + XMN_con**2
      call res_q2_limits(smax_r, w_r, xmlep_r,
     &                   q2_min_r, q2_max_r)

c     write(*,*) "q2min=",q2_min, "q2max=",q2_max
      if ((   q2_com .lt. -1.0*q2_min_r) 
     &  .or. (q2_com.gt.-1.0*q2_max_r)) then
        f_w = 0.0
        goto 123
      endif

      intflag_r = intflag_com
      iorg_r    = iorg_com
      ipflag_r  = ipflag_com

C     Need to loop over lambda sign
C     Need the neutrino type (anti-neutrino or not)

C     Neutrino type (neutrino or anti-neutrino)
      if (IPFLAG_r .eq. 0) then
        Lambda(1) =  1
        Lambda(2) = -1
      else
        Lambda(1) = -1
        Lambda(2) =  1
      endif

C     First call
      call rsdcrs(intflag_r, iorg_r, lambda(1), 
     &  xmlep_r, enu_r, -1.0*q2_r, w_r, 
     &  dxsec_r(1), dxsec3_r(1))

C     Second call
      call rsdcrs(intflag_r, iorg_r, lambda(2),
     &  xmlep_r, enu_r, -1.0*q2_r, w_r,
     &  dxsec_r(2), dxsec3_r(2))

C     Write the two single precisions to a double precision to return
      dxsec  = dxsec_r(1) + dxsec_r(2)
      dxsec3 = dxsec3_r(1) + dxsec3_r(2)

      f_w = dxsec

 123  return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
