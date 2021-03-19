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

C     Input reader
      character*150 argv1

C     Differential cross-sectino
      real dxsec

      integer iparlop, imode
      real enu, w, q2, costh, phi
      real xmlep

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

      if (iargc().ne.7) then
         print *, 'Cross-section table calculator for'
         print *, 'Minoo Kabirnezhad single pion model'
         print *, '*************************'
         print *, 'iargc is ', iargc()
         print *, 'args are either:'
         print *, 'nu type(1=e,2=m,3=t), IMODE, Enu, W, Q2, costh, phi'
         print *, '*************************'
         stop
      endif

C     Set neutrino type
      call getarg(1, argv1)
      read(argv1, '(I2)') iparlop

C     Set interaction type
      call getarg(2, argv1)
      read(argv1, '(I2)') imode

C     Set Enu
      call getarg(3, argv1)
      read(argv1, '(F10.10)') enu

C     Set W
      call getarg(4, argv1)
      read(argv1, '(F10.10)') w

C     Set Q2
      call getarg(5, argv1)
      read(argv1, '(F10.10)') q2

C     Set costh
      call getarg(6, argv1)
      read(argv1, '(F10.10)') costh

C     Set phi
      call getarg(7, argv1)
      read(argv1, '(F10.10)') phi

      if ((IPARLOP.GT.3).or.(IPARLOP.lt.1)) then
        write(*,*) 'IPARLOP=', IPARLOP, ' invalid neutrino flavour'
        write(*,*) '1 = nue, 2 = numu, 3 = nutau'
        return
      endif

C     Vector form factor mass
      XMVSPI = 0.84

      XMLEP = 0.
      IF (IPARLOP.eq.1)  XMLEP=XME
      IF (IPARLOP.eq.2)  XMLEP=XMMU
      IF (IPARLOP.eq.3)  XMLEP=XMTAU
C       If NC interaction set lepton mass to zero
      IF (IMODE.gt.14 .or. (imode.gt.3 .and. imode.lt.8)) then
        XMLEP=0.0E0
      ENDIF

      call mkdcrs(imode, enu, -1.0*q2, w, 
     &            costh, phi, xmlep, dxsec)

      write(*,*) "**********************"
      write(*,*) "One integration:"
      write(*,*) "nu: ", iparlop
      write(*,*) "mode: ", imode
      write(*,*) "enu: ", enu
      write(*,*) "q2: ", q2
      write(*,*) "had mass:", w
      write(*,*) "theta: ", costh
      write(*,*) "phi: ", phi
      write(*,*) "xmlep: ", xmlep
      write(*,*) "xsec: ", dxsec

      return
      end
