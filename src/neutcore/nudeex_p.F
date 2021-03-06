      Subroutine nudeex_p(piabs,shell,nout,idout,pout)
C
C+  Return particles from de-excitation of Oxygen-16
C+
C+  Input:
C+            piabs       (logical)  .true. if a pion was absorbed
C+            shell       (integer)  0: determine shell randomly
C+                                   1: s(1/2)
C+                                   2: p(3/2)
C+                                   3: the other shell states
C+  Output:
C+            nout        (integer)  number of particles from de-excitation
C+            idout()     (integer)  PDG codes of emitted particles
C+            pout(4,)    (real)     momentum (in MeV) of emitted particles
C+
C+  Creation Date and Author:
C+       2011.02.04 ; First version by K.Ueno, pion absorption case by T.Mori
C+       2011.09.25 ; Update of spectroscopic factors by K.Ueno
C+       2011.10.27 ; Update of the excitation-energy region of the s-hole state by K.Ueno
C+       2014.06.24 : Change to Ankowski spectral factor 
C+       2014.08.11 : Revise Ankowski spectral factor 
C+       2017.11.22 ; move others from p1/2 to s1/2
C+       2018.11.10 ; small bug above change

      implicit none
      
      integer maxout
      parameter (maxout = 5)
      logical piabs
      integer shell, nout, idout(maxout)
      real pout(4,maxout)
      
      integer max_ngamma_neut
      parameter (max_ngamma_neut = maxout)
      
      real etotal_gamma, sint, cost, phi, ecum(max_ngamma_neut)
      real rmode, cum, pmom, el, eg1, eg2, ep, en, ex
      integer ip, initial, mode
      
      real rlu, xrndm
      integer dummy, lunout
      parameter (lunout = 6)
      integer proton_code, neutron_code, photon_code
      parameter (proton_code = 2212)
      parameter (neutron_code = 2112)
      parameter (photon_code = 22)
      real proton_mass, neutron_mass, photon_mass
      parameter (proton_mass = 938.27)
      parameter (neutron_mass = 939.56)
      parameter (photon_mass = 0.)
      real c14_mass, n14_mass
      parameter (c14_mass = 13043.94) ! Mass of C14 g.s. (MeV)
      parameter (n14_mass = 13043.78) ! Mass of N14 g.s. (MeV)
      real c14p_eth, n14n_eth
      parameter (c14p_eth = 10.21) ! Eth (MeV) of 2-body decay from N15 to C14+p
      parameter (n14n_eth = 10.83) ! Eth (MeV) of 2-body decay from N15 to N14+n
      real qq2pi
      parameter (qq2pi = 6.283185)
C     
C+  Data from Ejiri's paper & Benhar's paper
C+  Use 18.8% for s(1/2), 37.5% for p(3/2) 6.32 MeV state,
C+  and 6% for p(3/2) 9.93 + 10.7 MeV states.
C     
      integer max_shell
      parameter (max_shell = 3)
      real prob_shell(max_shell)
C      data prob_shell/0.1055,0.3515,0.543/ ! s(1/2), p(3/2), no gamma
C      data prob_shell/0.543,0.3515,0.1055/ ! s(1/2), p(3/2), no gamma
      data prob_shell/0.4905,0.3515,0.158/ ! s(1/2), p(3/2), no gamma
      
C+ Description of p-shell proton hole modes
      integer pp_hole
      parameter (pp_hole = 6)
      real br_ppmode(pp_hole)
      data br_ppmode/0.3053,0.013425,0.000848,0.000363,0.002664,0.0289/
      real norm_br_ppmode
      data norm_br_ppmode/0.3515/ ! Total branching ratio relative to those above
      real e_ppgamma(pp_hole,2)
      data e_ppgamma
     $     /6.32,9.93,6.32,7.30,5.3,0.,  ! 1st deex
     $     0.  ,0.  ,3.61,2.63,4.64,0./  ! 2nd deex
      real e_ppproton(pp_hole)
      data e_ppproton /0.,0.,0.,0.,0.,0.5/
            
C+ Description of s-shell proton hole modes (N15 --> C13+d, C12+t, N14+n, C14+p)
      integer s_hole
      parameter (s_hole = 21)
      real br_smode(s_hole)
      data br_smode
     $     /0.03,0.0417,0.0167,0.0288,0.058,
     $     0.0674,0.0504,0.0288,0.0162,0.0034,0.0012,0.0023,0.0196,0.0661,
     $     0.0115,0.0041,0.0279,0.0196,0.0095,0.319,0.1778/
      real norm_br_smode
      data norm_br_smode/1./    ! Total branching ratio relative to those above
      real e_sgamma(s_hole)
      data e_sgamma
     $     /3.09,3.68,3.68,3.85,4.44,
     $     0.,4.92,3.38,5.69,5.11,5.83,5.11,6.44,7.03,
     $     0.,6.73,6.09,6.73,7.34,0.,0./
      real e_slevel(s_hole)
      data e_slevel
     $     /3.09,3.68,3.85,3.85,4.44,
     $     0.,4.92,5.69,5.69,5.83,5.83,6.45,6.45,7.03,
     $     0.,6.73,7.34,7.34,7.34,0.,0./
      
C+ Description of pion absorprion
      integer max_gamma_piabs
      parameter (max_gamma_piabs = 11)
      real br_pigamma(max_gamma_piabs), e_pigamma(max_gamma_piabs,2)
      data br_pigamma
     $ /0.002,0.048,0.013,0.003,0.019,0.01,0.04,0.007,0.005,0.017,0.836/
      data e_pigamma 
     $     /2.79,1.63,0.72,2.31,3.68,3.85,4.44,5.11,5.27,6.13,0., ! 1st deex
     $      2.31,2.31,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0./ ! 2nd deex
      
C
C      include 'macro.inc'
C      include 'nextrn.inc'
C
C+  Entry point
      
      nout = 0
      call vzero(idout,maxout)
      call vzero(pout,4*maxout)
      eg1 = 0.
      eg2 = 0.
      ep  = 0.
      en  = 0.
      rmode = RLU(dummy)
      mode = 1
      if (.not.piabs) then
C     
C+  If no pion absorption, using the following papers:
C+  H.Ejiri, Phys. Rev. C48, 1442 (1993)
C+  O.Benhar et al., Phys. Rev. D72, 053005 (2005)
C+  A.M.Ankowski et al., nucl-th/1110.0679v1 (2011)
C+  K.Kobayashi et al., nucl-ex/0604006 (2006)
C+  M.Yosoi, PhD thesis, Univ. of Kyoto (2003)
C     
         if (shell.eq.0) then   ! Determine level of hole
            initial = 1
            xrndm = RLU(dummy)
            cum = prob_shell(initial)
            do while(xrndm.gt.cum.and.initial.le.2)
               initial = initial + 1
               cum = cum + prob_shell(initial)
            enddo
         else
            initial = shell
         endif   
C     
C+  Short-cut: neglect de-excitation from other than s(1/2) and p(3/2) states
C     
         if (initial.eq.3) return
C     
C+  Select final state based on flavor and level of hole
C     
         if(initial.eq.1) then  ! s(1/2) hole
C     
C+  To give the energy of the nucleon from 2-body decay of 15N or 15O,
C+  this code determines the excitation energy.
C+  Benhar's paper suggests the distribution is almost uniform in 10.65 < Ex < 50.15 MeV
C+  over the nucleon momentum range from 0 to 300 MeV/c.
C     
            ex = 39.5*RLU(dummy) + 10.65
            cum = br_smode(mode)/norm_br_smode
            do while (rmode.gt.cum.and.mode.le.s_hole)
               mode = mode + 1
               cum = cum + br_smode(mode)/norm_br_smode
            enddo
            if (mode.gt.s_hole) then
               write(lunout,*) 'nudeex: error locating s(1/2) hole'
               return
            endif
            eg1 = e_sgamma(mode)
            el  = e_slevel(mode)
            if(mode.ge.6.and.mode.le.14) then ! N15 --> N14+n
               en = max(0.,(ex-el-n14n_eth)
     $              *(n14_mass+el)/(neutron_mass+n14_mass+el))
            else if(mode.ge.15.and.mode.le.19) then ! N15 --> C14+p
               ep = max(0.,(ex-el-c14p_eth)
     $              *(c14_mass+el)/(proton_mass+c14_mass+el))
            else if(mode.eq.20) then ! neutron from 3-body decay of N15
               en = 5.*RLU(dummy)
            else if(mode.eq.21) then ! proton from 3-body decay of N15
               ep = 5.*RLU(dummy)
            endif
            
         else if(initial.eq.2) then ! p(3/2) hole
            
            cum = br_ppmode(mode)/norm_br_ppmode
            do while (rmode.gt.cum.and.mode.le.pp_hole)
               mode = mode + 1
               cum = cum + br_ppmode(mode)/norm_br_ppmode
            enddo
            if (mode.gt.pp_hole) then
               write(lunout,*) 'nudeex: error locating p(3/2) proton hole'
               return
            endif
            eg1 = e_ppgamma(mode,1)
            eg2 = e_ppgamma(mode,2)
            ep  = e_ppproton(mode)
         endif
      else
C     
C+  Come here if a pion was absorbed. Using the following paper:
C+  H.D.Engelhardt et al., Nucl. Phys. A258, 480 (1976)
C     
         cum = br_pigamma(mode)
         do while (rmode.gt.cum.and.mode.le.max_gamma_piabs)
            mode = mode + 1
            cum = cum + br_pigamma(mode)
         enddo
         eg1 = e_pigamma(mode,1)
         eg2 = e_pigamma(mode,2)
      endif
C     
C+  Fill output particles
C     
      if (eg1.gt.0.) then
         nout = nout + 1
         idout(nout) = photon_code
         pout(4,nout) = eg1      
      endif
      if (eg2.gt.0.) then
         nout = nout + 1
         idout(nout) = photon_code
         pout(4,nout) = eg2
      endif
      if (ep.gt.0.) then
         nout = nout + 1
         idout(nout) = proton_code
         pout(4,nout) = proton_mass + ep
      endif
      if (en.gt.0.) then
         nout = nout + 1
         idout(nout) = neutron_code
         pout(4,nout) = neutron_mass + en
      endif
      
      do ip = 1, nout
         cost = min(max(2.*(RLU(dummy)-0.5),-1.),1.)
         if (abs(cost).eq.1.) then
            sint = 0.
         else
            sint = sqrt(1.-cost**2)
         endif
         phi  = qq2pi*RLU(dummy)
         if (idout(ip).eq.proton_code) then
            pmom = sqrt(pout(4,ip)**2 - proton_mass**2)
         else if (idout(ip).eq.neutron_code) then
            pmom = sqrt(pout(4,ip)**2 - neutron_mass**2)
         else
            pmom = pout(4,ip)
         endif
         pout(1,ip) = pmom*cos(phi)*sint
         pout(2,ip) = pmom*sin(phi)*sint
         pout(3,ip) = pmom*cost
      enddo
      
      return
      end
      
