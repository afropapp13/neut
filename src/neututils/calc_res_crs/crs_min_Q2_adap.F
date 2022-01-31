      program crscalc_minoo_q2

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

      integer iq2
      integer nq2
      parameter(nq2=100)
      real dq2

      character*21 fname_ref(3)
      character*100 fname
      character*50 argv1
      character*50 argv2

      data fname_ref/
     $     'spi_nue_xsec_mk_ma_Q2',
     $     'spi_num_xsec_mk_ma_Q2',
     $     'spi_nut_xsec_mk_ma_Q2'/

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
      double precision Xsec(3,17,nq2), Enu(3,17,nq2)

      double precision E
      integer ie, iflag
      double precision sigma, xmlep
      real wwmax
      double precision q2, q2_min,Q2_max

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
c     DO 110 IPFLAG = 0,0

C     IFLAG is looping what interaction we're doing (listed above in
C     data)
      DO 105  IFLAG = 1,7
c     DO 105  IFLAG = 4,7


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

       if (IPFLAG .EQ. 0) THEN
         write(*,*)  "#",TYPESn(IFLAG)
       else 
         write(*,*)  "#",TYPESa(IFLAG)
       endif

      imode_com = imode

C       Convert from REAL to double precision
C       Because NEUT insists on REAL rather than double...
      xmn_con   = dble(XMN)
      xmlep_con = dble(XMLEP)
      xmpi_con  = dble(XMPI)

C     Fixed Enu (GeV)
      E = 1
C     Set common block Enu to Enu from the table
      enu_com = E

C     Do Q2 between 0 and 1
      q2_min = 0
      q2_max = 1

      write(*,*) "Using fixed energy E=", E
      write(*,*) "Running between", q2_min, " < Q2 < ", q2_max

      dQ2 = (Q2_max-Q2_min)/NQ2
      DO 1337 iq2=0,nq2
      Q2 = Q2_min+dQ2*IQ2
      q2_com = q2
        
C       Set the cross-section for this Enu to 0
      SIGMA = 0.

      call int_all(imode, SIGMA)

C       Write the arrays with the cross-section and the energy
c       Xsec(IPARLOP,IMODE,IE) = SIGMA
c       W( IPARLOP,IMODE,IE) = W

      write(*,*) q2, SIGMA

C      End here for the energy loop
c95    CONTINUE

C      Write stuff out
c       if (IPFLAG .EQ. 0) THEN
c         write(*,*)  "#",TYPESn(IFLAG)
c       else 
c         write(*,*)  "#",TYPESa(IFLAG)
c       endif
 1337 CONTINUE
 105  CONTINUE
 110  CONTINUE

c ************************************************************
c     Now write the output
c ************************************************************
c       do I = iparlop,iparlop

c       if (WWMAX .ne. 2.0) then
c         write(fname,'(A,F4.2,A,F4.2,A,F4.2,A,F4.2,A)')
c    $    fname_ref(I), XMASPI, '_ca5', RCA5ISPI,
c    $    '_mabkg', XMABKGM, '_W', wwmax, '.dat'
c         write(*,*) 'Opening file ', fname
c       else 
c         write(fname,'(A,F4.2,A,F4.2,A,F4.2,A)')
c    $     fname_ref(I), XMASPI, '_ca5', RCA5ISPI,
c    $     '_mabkg', XMABKGM, '.dat'
c         write(*,*) 'Opening file ', fname
c       endif

c       open(60,FILE=FNAME,ERR=9999)

c       write(60,'(A,A,A,F4.2,A,F4.2,A,F4.2,A,F4.2)')
c    $        '# 1pi MK :', ffn(1), ':MA=', XMASPI, 
c    $        ', CA50=;', RCA5ISPI,
c    $        ', MABKG=;', XMABKGM,
c    $        ', W_MAX=',WWMAX

c       do IE=1,IEMAX
c       write(60,'(I4,I2,F8.2,7E13.6)')
c    $           IPTBL(I),1,Enu(I,1,IE), 
c    $           XSec(I,1,IE), Xsec(I,2,IE), Xsec(I,3,IE),
c    $           Xsec(I,4,IE), Xsec(I,5,IE), Xsec(I,6,IE),
c    $           Xsec(I,7,IE)
c       enddo

c     This is literally a copy of the above other than the "2"
c     This is done because the event generation code assumes the table
C     for Minoo's model has the same structure as Rein-Sehgal
C     (intentional to avoide code duplicates). The Rein-Sehgal code has
C     two cross-sections saved: 1) the total interaction mode cross-section
C     for all resonance, 2) the interaction mode cross-section coming
c     only from the Delta resonance (which is where the 2 comes from)
c       do IE=1,IEMAX
c       write(60,'(I4,I2,F8.2,7E13.6)')
c    $           IPTBL(I),2,Enu(I,1,IE),
c    $           Xsec(I,1,IE), Xsec(I,2,IE), Xsec(I,3,IE),
c    $           Xsec(I,4,IE), Xsec(I,5,IE), Xsec(I,6,IE),
c    $           Xsec(I,7,IE)
c       enddo

c       do IE=1,IEMAX
c       write(60,'(I4,I2,F8.2,7E13.6)')
c    $        -1*IPTBL(I),1,Enu(I,1,IE),
c    $           Xsec(I,11,IE), Xsec(I,12,IE), Xsec(I,13,IE),
c    $           Xsec(I,14,IE), Xsec(I,15,IE), Xsec(I,16,IE),
c    $           Xsec(I,17,IE)
c       enddo

c     Agan a copy of the above but with the 2 in place instead of 1
c       do IE=1,IEMAX
c       write(60,'(I4,I2,F8.2,7E13.6)')
c    $        -1*IPTBL(I),2,Enu(I,1,IE),
c    $           Xsec(I,11,IE), Xsec(I,12,IE), Xsec(I,13,IE),
c    $           Xsec(I,14,IE), Xsec(I,15,IE), Xsec(I,16,IE),
c    $           Xsec(I,17,IE)
c       enddo

c       enddo ! End the iparlop loop (which really does nothing...)
      goto 1000


 1000 continue
      STOP
 9999 write(*,*) "Error in opening file."
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Function which gives the integral over all the variables EXCEPT W
C     This is the function we integrate over when we want to integrate
C     over W
      subroutine int_all(imode, sigma)

      double precision w,smax
      double precision w_min, w_max

      integer limit_w
      parameter (limit_w = 1000)
      integer lenw_w
      parameter (lenw_w = limit_w * 4)

      double precision q2_min ! The lower limit for q2
      double precision q2_max ! The upper limit for q2
      external f_w

C     The DQAG variables
      double precision epsabs_w, epsrel_w, abserr_w
      double precision result_w
      double precision work_w(lenw_w)
      double precision sigma

      integer neval_w, ier_w, iwork_w(limit_w), key_w, last_w

C     Have a "handy" common block to make sure all functions know what
C     Q2, W, Enu, theta, phi we're currently on
      double precision enu_com, w_com, q2_com, th_com, phi_com
      common /kine/ enu_com, w_com, q2_com, th_com, phi_com

      double precision xmn_con, xmlep_con, xmpi_con, w_max_con
      common /const/ xmn_con, xmlep_con, xmpi_con, w_max_con

      result_w = 0.

      epsabs_w = 0.001
      epsrel_w = 0.001
      key_w = 1

      smax = (2.*XMN_con*enu_com) + XMN_con**2

      w_min = XMN_con + XMPI_con
      w_max = MIN(SQRT(SMAX) - XMLEP_con, w_max_con)

C     Integrate out Q2
C     N.B. reverse the q2_max and q2_min because both quantities are
C     negative
      call dqag (f_w, w_min, w_max, epsabs_w, epsrel_w, key_w,
     &           result_w, abserr_w, neval_w, ier_w, limit_w,
     &           lenw_w, last_w, iwork_w, work_w)

c     write(*,*) "**********************"
c     write(*,*) "W INT"
c     write(*,*) "Result = ", result_w, "+-", abserr_w
c     write(*,*) "Evals = ", neval_w
c     write(*,*) "Error code = ", ier_w

      SIGMA = result_w
C     write(*,*) f_w

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Function which gives the integral over all the variables EXCEPT Q2
C     This is the function we integrate over when we want to integrate
C     over Q2
      double precision function f_w(w)

      double precision w

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
      double precision xmn_con, xmlep_con, xmpi_con
      common /const/ xmn_con, xmlep_con, xmpi_con

      real smax_r, w_com_r, xmlep_con_r, q2_min_r, q2_max_r

      epsabs_th = 0.001
      epsrel_th = 0.001
      key_th = 1

C     Set the common block q2 to the function input q2
      w_com = w

      smax = (2.*XMN_con*enu_com) + XMN_con**2
      smax_r = smax

      xmlep_con_r = xmlep_con

C     Now that we have an enu and w we also have the q2 limits
C     The return is q2_min and q2_max
      call res_q2_limits(smax_r, w_com_r, xmlep_con_r,
     &                   q2_min_r, q2_max_r)

c     write(*,*) "q2min=",q2_min, "q2max=",q2_max
      if ((   q2_com .lt. -1.0*q2_min_r) 
     &  .or. (q2_com.gt.-1.0*q2_max_r)) then
        f_w = 0.0
        goto 123
      endif

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

      f_w = result_th


 123  return
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

      call mkdcrs(imode_com, enu_r, -1.0*q2_r, w_r, 
     &            th_r, phi_r, xmlep_r, dxsec_r)

c     write(*,*) "**********************"
c     write(*,*) "One integration:"
c     write(*,*) "mode: ", imode_com
c     write(*,*) "enu: ", enu_r
c     write(*,*) "q2: ", q2_r
c     write(*,*) "had mass:", w_r
c     write(*,*) "theta: ", th_r
c     write(*,*) "phi: ", phi_r
c     write(*,*) "xmlep: ", xmlep_r
c     write(*,*) "xsec: ", dxsec_r

      dxsec = dxsec_r
      f_phi = dxsec

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Get the wlimits to integrate out
c     subroutine wlimits

c       double E_1s(double E) {return M*E/(sqrt(M*M + 2.*M*E)) ;}

c       double sqrt_s(double E) {return sqrt(M*M + 2.*M*E) ;}

c       double lambda (double W, double E) {return (sqrt_s(E)*sqrt_s(E)
c       - m_l*m_l - W*W)*(sqrt_s(E)*sqrt_s(E) - m_l*m_l - W*W)
c       -4.*m_l*m_l*W*W ;}

c       double Wdown() {return M+m_pi;}

c        double Wup(double E) {return sqrt_s(E) -m_l; }

c         double Qmin(double W, double E) {
c            if (W<Wdown()|| W>Wup(E))
c                 return 0.;
c                    else
c                         return -m_l*m_l +
c                         (E_1s(E)/sqrt_s(E))*(sqrt_s(E)*sqrt_s(E) +
c                         m_l*m_l - W*W - sqrt(lambda(W,E)) ) ;      
c                         }
c                          double Qmax(double W, double E) {
c                              if (W<Wdown()|| W>Wup(E))
c                                   return 0.;
c                                      else
c                                           return -m_l*m_l +
c                                           (E_1s(E)/sqrt_s(E))*(sqrt_s(E)*sqrt_s(E)
c                                           + m_l*m_l - W*W +
c                                           sqrt(lambda(W,E)) ) ; }


