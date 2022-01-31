      program crscalc_minoo

C     An adaptive integrator using the QUADPACK lib
C     Does the integral over phi, theta, q2, w using dqag
C     dqag uses a Gauss-Kronrod technique, see quadpack.f
C     for more info and
C     https://people.sc.fsu.edu/~jburkardt/f77_src
C                             /quadpack/quadpack.html
C
C     Author: Clarence Wret (cwret14@imperial.ac.uk,
C     clarence.wret@gmail.com)

      Implicit None

#include "neutmodel.h"
#include "rscons.h"
#include "necard.h"

      integer imode,ipflag, imode_com
      common imode_com

      integer iparlop

      integer*4 iptbl(3)
      data iptbl /12,14,16/

      integer i

      integer iemax
      parameter(IEMAX=18)
      double precision etbl(18)
c     parameter(IEMAX=24)
c     double precision etbl(24)

      DATA ETBL/
     $     0.370, 0.400, 0.500, 0.625, 0.750,
     $     0.875, 1.000, 1.125, 1.250, 1.375,
     $     1.500, 2.000, 2.500, 3.000, 3.500,
     $     4.750, 6.000, 10.00/
c    $     4.750, 6.000, 10.00, 50.00, 100.0,
C    $     200.0, 400.0, 650.0, 1000.0/

      character*18 fname_ref(3)
      character*100 fname
      character*50 argv1
      character*50 argv2

      data fname_ref/
     $     'spi_nue_xsec_mk_ma',
     $     'spi_num_xsec_mk_ma',
     $     'spi_nut_xsec_mk_ma'/

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
     $     "nubar n -> mu+   n pi-",
     $     "nubar p -> mu+   n pi0",
     $     "nubar p -> mu+   p pi-",
     $     "nubar n -> nubar n pi0",
     $     "nubar n -> nubar p pi-",
     $     "nubar p -> nubar p pi0",
     $     "nubar p -> nubar n pi+"/
      character*30 ffn(1)
      data ffn/
     $     "Minoo Kabirnezhad default"/

C     xsec(number of lepton generations, number of anti-part/part, # of
C     interaction modes, number of energies)
C     Xsec = x-sec, Xsec3 = x-sec from delta, Enu = neutrino E
      double precision Xsec(3,17,IEMAX), Enu(3,17,IEMAX)

      double precision E
      integer ie, iflag
      double precision sigma, w_max, w_min, wwmax, xmlep

C     Have a "handy" common block to make sure all functions know what
C     Q2, W, Enu, theta, phi we're currently on
      double precision enu_com, w_com, q2_com, th_com, phi_com
      common /kine/ enu_com, w_com, q2_com, th_com, phi_com

C     Another common block with the rscons.h stuff in
      double precision xmn_con, xmlep_con, xmpi_con, w_max_con
      common /const/ xmn_con, xmlep_con, xmpi_con, w_max_con

C     The maximum value of W (hadronic mass)
      WWMAX=2.0
C     Set the necard.F settings here
C     Choose Minoo's model
      MDLSPI=3
C     Choose MABKG
      XMABKGM=1.03
C     Choose MA
      XMANFFRES=1.06
C     Choose MV
      XMVNFFRES=0.84
C     Choose CA5
      RNECA5I=1.02
      call nefillmodel

      if (iargc().ne.1 .and. iargc().ne.3 .and. iargc().ne.4 .and.
     &    iargc().ne.5) then
         print *, 'Cross-section table calculator for'
         print *, 'Minoo Kabirnezhad single pion model'
         print *, '*************************'
         print *, 'iargc is ', iargc()
         print *, 'args are either:'
         print *, 'nu type(1=e,2=m,3=t), MARES, CA50, MABKG, WCUT'
         print *, 'nu type(1=e,2=m,3=t), MARES, CA50, WCUT'
         print *, 'nu type(1=e,2=m,3=t), MARES, CA50'
         print *, 'nu type(1=e,2=m,3=t)'
         print *, '*************************'
         print *, 'If unspecified, defaults are:'
         print *, 'MARES =', XMASPI
         print *, 'MVRES =', XMVSPI
         print *, 'CA50  =', RCA5ISPI
         print *, 'MABKG =', XMABKGM
         print *, 'WCUT  =', WWMAX
         print *, 'taken from neutcore/necard.F'
         stop
      endif

C     Set MARES
      if (iargc().eq.3 .or. iargc().eq.4 .or. iargc().eq.5) then
        call getarg(2, argv2)
        read(argv2, '(F4.2)') XMASPI
      endif

C     Set CA5(0)
      if (iargc().eq.3 .or. iargc().eq.4 .or. iargc().eq.5) then
        call getarg(3, argv2)
        read(argv2, '(F4.2)') RCA5ISPI
      endif

C     Set MABKG
      if (iargc().eq.5) then
        call getarg(4, argv2)
        read(argv2, '(F4.2)') XMABKGM
      endif

C     WMAX input (W < 2.0 GeV for single pion production)
      if (iargc() .eq. 4 .or. iargc().eq.5) then
        call getarg(iargc(), argv2)
        read(argv2, '(F4.2)') WWMAX
      endif

      print *, 'MARES =', XMASPI
      print *, 'MVRES =', XMVSPI
      print *, 'CA50  =', RCA5ISPI
      print *, 'MABKG =', XMABKGM
      print *, 'WCUT  =', WWMAX

C     Write the w max to the common block
      w_max_con = wwmax

      call getarg(1, argv1)
      read(argv1, '(I10)') IPARLOP
      write(*,'(A,I10)') 'IPARLOP =', IPARLOP

      if ((IPARLOP.GT.3).or.(IPARLOP.lt.1)) then
        write(*,*) 'IPARLOP=', IPARLOP, ' invalid neutrino flavour'
        write(*,*) '1 = nue, 2 = numu, 3 = nutau'
        goto 9999
      endif

C     Vector form factor mass
      XMVSPI = 0.84

      DO 110 IPFLAG = 0,1

C     IFLAG is looping what interaction we're doing (listed above in
C     data)
      DO 105  IFLAG = 1,7

C     Set the lepton mass depending on the neutrino flavour and
C     interaction mode
        XMLEP = 0.
        IMODE = IFLAG + 10*IPFLAG
        IF (IPARLOP.eq.1)  XMLEP=XME
        IF (IPARLOP.eq.2)  XMLEP=XMMU
        IF (IPARLOP.eq.3)  XMLEP=XMTAU
C       If NC interaction set lepton mass to zero
        IF (IFLAG.gt.3)    then
          XMLEP=0.0E0
        ENDIF

        imode_com = imode

C       Convert from REAL to double precision
C       Because NEUT insists on REAL rather than double...
        xmn_con   = dble(XMN)
        xmlep_con = dble(XMLEP)
        xmpi_con  = dble(XMPI)

        write(*,*) "Starting energy loop"

C       Loop over the Enu table
        DO 95 IE = 1,IEMAX
C       Read the energy from the table
          E = etbl(IE)

C         Set the cross-section for this Enu to 0
          SIGMA = 0.

C         Set common block Enu to Enu from the table
          enu_com = E

C         Integrate out everything to only leave sigma(Enu)
C         This essentially goes W->Q2->costheta->phi
          call int_all(imode, E, SIGMA)

C         Write the arrays with the cross-section and the energy
          Xsec(IPARLOP,IMODE,IE) = SIGMA
          Enu( IPARLOP,IMODE,IE) = E

          write(*,*) E, SIGMA

C      End here for the energy loop
 95    CONTINUE

C      Write stuff out
       if (IPFLAG .EQ. 0) THEN
         write(*,*)  "#",TYPESn(IFLAG)
       else 
         write(*,*)  "#",TYPESa(IFLAG)
       endif

 105  CONTINUE
 110  CONTINUE

C ************************************************************
C     Now write the output
C ************************************************************
      do I = iparlop,iparlop

      if (WWMAX .ne. 2.0) then
        write(fname,'(A,F4.2,A,F4.2,A,F4.2,A,F4.2,A)')
     $    fname_ref(I), XMASPI, '_ca5', RCA5ISPI,
     $    '_mabkg', XMABKGM, '_W', wwmax, '.dat'
        write(*,*) 'Opening file ', fname
      else 
        write(fname,'(A,F4.2,A,F4.2,A,F4.2,A)')
     $     fname_ref(I), XMASPI, '_ca5', RCA5ISPI,
     $     '_mabkg', XMABKGM, '.dat'
        write(*,*) 'Opening file ', fname
      endif

      open(60,FILE=FNAME,ERR=9999)

      write(60,'(A,A,A,F4.2,A,F4.2,A,F4.2,A,F4.2)')
     $        '# 1pi MK :', ffn(1), ':MA=', XMASPI, 
     $        ', CA50=;', RCA5ISPI,
     $        ', MABKG=;', XMABKGM,
     $        ', W_MAX=',WWMAX

      do IE=1,IEMAX
      write(60,'(I4,I2,F8.2,7E13.6)')
     $           IPTBL(I),1,Enu(I,1,IE), 
     $           XSec(I,1,IE), Xsec(I,2,IE), Xsec(I,3,IE),
     $           Xsec(I,4,IE), Xsec(I,5,IE), Xsec(I,6,IE),
     $           Xsec(I,7,IE)
      enddo

C     This is literally a copy of the above other than the "2"
C     This is done because the event generation code assumes the table
C     for Minoo's model has the same structure as Rein-Sehgal
C     (intentional to avoide code duplicates). The Rein-Sehgal code has
C     two cross-sections saved: 1) the total interaction mode cross-section
C     for all resonance, 2) the interaction mode cross-section coming
C     only from the Delta resonance (which is where the 2 comes from)
      do IE=1,IEMAX
      write(60,'(I4,I2,F8.2,7E13.6)')
     $           IPTBL(I),2,Enu(I,1,IE),
     $           Xsec(I,1,IE), Xsec(I,2,IE), Xsec(I,3,IE),
     $           Xsec(I,4,IE), Xsec(I,5,IE), Xsec(I,6,IE),
     $           Xsec(I,7,IE)
      enddo

      do IE=1,IEMAX
      write(60,'(I4,I2,F8.2,7E13.6)')
     $        -1*IPTBL(I),1,Enu(I,1,IE),
     $           Xsec(I,11,IE), Xsec(I,12,IE), Xsec(I,13,IE),
     $           Xsec(I,14,IE), Xsec(I,15,IE), Xsec(I,16,IE),
     $           Xsec(I,17,IE)
      enddo

C     Agan a copy of the above but with the 2 in place instead of 1
      do IE=1,IEMAX
      write(60,'(I4,I2,F8.2,7E13.6)')
     $        -1*IPTBL(I),2,Enu(I,1,IE),
     $           Xsec(I,11,IE), Xsec(I,12,IE), Xsec(I,13,IE),
     $           Xsec(I,14,IE), Xsec(I,15,IE), Xsec(I,16,IE),
     $           Xsec(I,17,IE)
      enddo

      enddo ! End the iparlop loop (which really does nothing...)
      goto 1000


 1000 continue
      STOP

 9999 write(*,*) "Error in opening file."
      END ! End the program crscalc


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
C     a subroutine to perform the integration over all the variables

      subroutine int_all(imode, E, SIGMA)

      integer limit_w
      parameter (limit_w = 1000)
      integer lenw_w
      parameter (lenw_w = limit_w * 4)

      double precision W_min ! The lower limit for W
      double precision W_max ! The upper limit for W
      external f_w

C     Variables for the integrator (dqag in quadpack)
      double precision epsabs_w, epsrel_w, abserr_w
      double precision result_w
      double precision work_w(lenw_w)

      integer neval_w, ier_w, iwork_w(limit_w), key_w, last_w

      double precision smax, E
      integer imode
      double precision sigma ! The return

C     Have a "handy" common block to make sure all functions know what
C     Q2, W, Enu, theta, phi we're currently on
      double precision enu_com, w_com, q2_com, th_com, phi_com
      common /kine/ enu_com, w_com, q2_com, th_com, phi_com

      double precision xmn_con, xmlep_con, xmpi_con, w_max_con
      common /const/ xmn_con, xmlep_con, xmpi_con, w_max_con

      enu_com = E
      imode_com = imode

      result_w = 0.
      
C     The maximum center of mass energy
      smax = 2.*enu_com*XMN_con + XMN_con**2

      w_min = XMN_con + XMPI_con
      w_max = MIN(SQRT(SMAX) - XMLEP_con, w_max_con)
      epsabs_w = 0.001
      epsrel_w = 0.001
      key_w = 1

C     This does the integral over W
      call dqag (f_w, w_min, w_max, epsabs_w, epsrel_w, key_w, result_w,
     &           abserr_w, neval_w, ier_w, limit_w, lenw_w, last_w,
     &           iwork_w, work_w)

C     write(*,*) "**********************"
C     write(*,*) "W INT"
C     write(*,*) "Result = ", result_w, "+-", abserr_w
C     write(*,*) "Evals = ", neval_w
C     write(*,*) "Error code = ", ier_w

      SIGMA = result_w

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Function which gives the integral over all the variables EXCEPT W
C     This is the function we integrate over when we want to integrate
C     over W
      double precision function f_w(w)

      double precision w,smax

      integer limit_q2
      parameter (limit_q2 = 1000)
      integer lenw_q2
      parameter (lenw_q2 = limit_q2 * 4)

      double precision q2_min ! The lower limit for q2
      double precision q2_max ! The upper limit for q2
      external f_q2

C     The DQAG variables
      double precision epsabs_q2, epsrel_q2, abserr_q2
      double precision result_q2
      double precision work_q2(lenw_q2)

      integer neval_q2, ier_q2, iwork_q2(limit_q2), key_q2, last_q2
      
C     Have a "handy" common block to make sure all functions know what
C     Q2, W, Enu, theta, phi we're currently on
      double precision enu_com, w_com, q2_com, th_com, phi_com
      common /kine/ enu_com, w_com, q2_com, th_com, phi_com

      double precision xmn_con, xmlep_con, xmpi_con
      common /const/ xmn_con, xmlep_con, xmpi_con

      real smax_r, w_com_r, xmlep_con_r, q2_min_r, q2_max_r

      result_q2 = 0.

      epsabs_q2 = 0.001
      epsrel_q2 = 0.001
      key_q2 = 1

C     Set the common block  W to the function w
      w_com = w
      w_com_r = w_com

      smax = (2.*XMN_con*enu_com) + XMN_con**2
      smax_r = smax

      xmlep_con_r = xmlep_con

C     Now that we have an enu and w we also have the q2 limits
C     The return is q2_min and q2_max
      call res_q2_limits(smax_r, w_com_r, xmlep_con_r,
     &                   q2_min_r, q2_max_r)

      q2_min = q2_min_r
      q2_max = q2_max_r

C     Integrate out Q2
C     N.B. reverse the q2_max and q2_min because both quantities are
C     negative
      call dqag (f_q2, q2_max, q2_min, epsabs_q2, epsrel_q2, key_q2,
     &           result_q2, abserr_q2, neval_q2, ier_q2, limit_q2,
     &           lenw_q2, last_q2, iwork_q2, work_q2)

C     write(*,*) "**********************"
C     write(*,*) "Q2 INT"
C     write(*,*) "Result = ", result_q2, "+-", abserr_q2
C     write(*,*) "Evals = ", neval_q2
C     write(*,*) "Error code = ", ier_q2

      f_w = result_q2
C     write(*,*) f_w

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Function which gives the integral over all the variables EXCEPT Q2
C     This is the function we integrate over when we want to integrate
C     over Q2
      double precision function f_q2(q2)

      double precision q2

      integer limit_th
      parameter (limit_th = 1000)
      integer lenw_th
      parameter (lenw_th = limit_th * 4)

C     Have a "handy" common block to make sure all functions know what
C     Q2, W, Enu, theta, phi we're currently on
      double precision thmin, thmax
      external f_th ! Function for theta

      double precision epsabs_th, epsrel_th, abserr_th
      double precision result_th
      double precision work_th(lenw_th)

      integer neval_th, ier_th, iwork_th(limit_th), key_th, last_th

      double precision enu_com, w_com, q2_com, th_com, phi_com
      common /kine/ enu_com, w_com, q2_com, th_com, phi_com

      epsabs_th = 0.001
      epsrel_th = 0.001
      key_th = 1

C     Set the common block q2 to the function input q2
      q2_com = q2

      result_th = 0.

C     Maximum and minimum of costheta
      thmin = -1.0
      thmax = 1.0

C     Integrate out the theta, phi
      call dqag (f_th, thmin, thmax, epsabs_th, epsrel_th, key_th,
     &           result_th, abserr_th, neval_th, ier_th, limit_th,
     &           lenw_th, last_th, iwork_th, work_th)

C     write(*,*) "**********************"
C     write(*,*) "th INT"
C     write(*,*) "Result = ", result_th, "+-", abserr_th
C     write(*,*) "Evals = ", neval_th
C     write(*,*) "Error code = ", ier_th

      f_q2 = result_th


      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function f_th(th)

      double precision th

      integer limit_phi
      parameter (limit_phi = 1000)
      integer lenw_phi
      parameter (lenw_phi = limit_phi * 4)

      double precision phimin, phimax
      external f_phi ! Function for phieta

      double precision epsabs_phi, epsrel_phi, abserr_phi
      double precision result_phi
      double precision work_phi(lenw_phi)

      integer neval_phi, ier_phi, iwork_phi(limit_phi)
      integer key_phi, last_phi

C     Have a "handy" common block to make sure all functions know what
C     Q2, W, Enu, theta, phi we're currently on
      double precision enu_com, w_com, q2_com, th_com, phi_com
      common /kine/ enu_com, w_com, q2_com, th_com, phi_com

      epsabs_phi = 0.001
      epsrel_phi = 0.001
      key_phi = 1

C     Set the common block  W to the function w
      th_com = th

C     The limits on phi is just 0 to 2pi
      phimin = 0.
      phimax = 2.*3.1415926

      result_phi = 0.

C     Integrate out the phi finally
      call dqag (f_phi, phimin, phimax, epsabs_phi, epsrel_phi, key_phi,
     &           result_phi, abserr_phi, neval_phi, ier_phi, limit_phi,
     &           lenw_phi, last_phi, iwork_phi, work_phi)

C     write(*,*) "**********************"
C     write(*,*) "phi INT"
C     write(*,*) "Result = ", result_phi, "+-", abserr_phi
C     write(*,*) "Evals = ", neval_phi
C     write(*,*) "Error code = ", ier_phi

      f_th = result_phi

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function f_phi(phi)

      double precision phi

C     Have a "handy" common block to make sure all functions know what
C     Q2, W, Enu, theta, phi we're currently on
      double precision enu_com, w_com, q2_com, th_com, phi_com
      common /kine/ enu_com, w_com, q2_com, th_com, phi_com
      double precision xmn_con, xmlep_con, xmpi_con
      common /const/ xmn_con, xmlep_con,xmpi_con

      integer imode_com
      common imode_com

      real enu_r, q2_r, w_r, th_r, phi_r, xmlep_r, dxsec_r
      double precision dxsec

C     Set the common block  W to the function w
      phi_com = phi

C     write(*,*) enu_com, w_com, q2_com, th_com, phi_com

C     Now we have a E, Q2, W, costheta and phi
C     Can finally calculate a cross-section!

C     NEUT needs REAL as argument so convert our double to real
      enu_r = enu_com
      q2_r = q2_com
      w_r = w_com
      th_r = th_com
      phi_r = phi
      xmlep_r = xmlep_con

      call mkdcrs(imode_com, enu_r, q2_r, w_r, 
     &            th_r, phi_r, xmlep_r, dxsec_r)

C     write(*,*) "**********************"
C     write(*,*) "One integration:"
C     write(*,*) "mode: ", imode_com
C     write(*,*) "enu: ", enu_r
C     write(*,*) "q2: ", q2_r
C     write(*,*) "had mass:", w_r
C     write(*,*) "theta: ", th_r
C     write(*,*) "phi: ", phi_r
C     write(*,*) "xmlep: ", xmlep_r
C     write(*,*) "xsec: ", dxsec_r
      
      dxsec = dxsec_r
      f_phi = dxsec
 
      return
      end
