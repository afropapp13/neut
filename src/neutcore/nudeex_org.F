      Subroutine nudeex_org(piabs,idtarg,shell,nout,idout,pout)
C
C+  Return particles from de-excitation of Oxygen-16
C+
C+  Input:
C+            piabs       (logical)  .true. if a pion was absorbed
C+            idtarg      (integer)  PDG code of target (2112,2212)
C+            shell       (integer)  0: determine shell randomly
C+                                   1: s(1/2)
C+                                   2: p(3/2)
C+                                   3: p(1/2)
C+  Output:
C+            nout        (integer)  number of particles from de-excitation
C+            idout()     (integer)  PDG codes of emitted particles
C+            pout(4,)    (real)     momentum (in MeV) of emitted particles
C+
C+  Modified version (to be independent of Nuance routines)
C
      implicit none
      integer maxout
      parameter (maxout = 5)
      logical piabs
      integer idtarg, shell, nout, idout(maxout)
      real pout(4,maxout)
C
      integer max_ngamma_neut
      parameter (max_ngamma_neut = maxout)

      real etotal_gamma, sint, cost, phi, ecum(max_ngamma_neut)
      real rmode, cum, pmom, eg1, eg2, ep1, ep2, en
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
      real qq2pi
      parameter (qq2pi = 6.283185)
C
C+  Data from Ejiri's paper
C
      integer max_shell
      parameter (max_shell = 3)
      integer max_ngamma_ejiri
      parameter (max_ngamma_ejiri = 2)
      real cum_shell(max_shell)
      data cum_shell/0.25,0.75,1.0/  !25% s(1/2), 50% p(3/2), 25% p(1/2)
      real prob_shell(max_shell)
      data prob_shell/0.25,0.5,0.25/
      integer pshell_total, nshell_total
      parameter (pshell_total = 14)
      parameter (nshell_total = 10)
      real e_pgamma(2,pshell_total), e_ngamma(2,nshell_total)
      real e_pproton(2,pshell_total), e_nproton(2,nshell_total)
      real e_pneutron(pshell_total), e_nneutron(nshell_total)
C+ Description of proton hole modes
      data e_pgamma /0.,0.,6.32,0.,9.93,0.,6.32,3.61,0.,0.,0.,0.,
     1                7.03,0.,0.,0.,0.,0.,7.01,0.,0.,0.,3.5,0.,4.44,0.,
     1                15.1,0./
      data e_pproton/0.,0.,0.,0.,0.,0.,0.,0.,0.5,0.,0.,0.,0.,0.,1.6,0.,
     1              21.,0.,14.,0.,11.,0.,0.,0.,0.,0.,0.,0./
      data e_pneutron/5*0.,20.,13.,11.,0.,0.,2.,3*0./
C+ Description of neutron hole modes
      data e_ngamma/0.,0.,6.18,0.,0.,0.,7.03,0.,0.,0.,0.,0.,0.,0.,
     1               3.5,0.,4.44,0.,15.1,0./
      data e_nproton/0.,0.,0.,0.,24.,0.,17.,0.,14.5,1.5,0.,0.,2.,0.,
     1               0.,0.,0.,0.,0.,0./
      data e_nneutron/5*0.,18.,11.5,0.,0.,0./
      real br_pmode(pshell_total,3), br_nmode(nshell_total,3)
      data br_pmode/5*0.,0.02,0.02,0.01,0.02,0.02,0.03,0.104,0.013,0.013,
     1              0.,0.41,0.0234,0.0066,0.03,6*0.,0.024,0.003,0.003,
     1             0.25,13*0./
      data br_nmode/0.,0.,0.02,0.02,0.01,0.02,0.02,0.128,0.016,0.016,
     1             0.,0.44,5*0.,0.048,0.006,0.006,
     1              0.25,9*0./

C
C      include 'macro.inc'
C      include 'nextrn.inc'
C
C+  Entry point
C
      nout = 0
      call vzero(idout,maxout)
      call vzero(pout,4*maxout)
      if (idtarg .eq. 0) return
      if (.not. piabs) then
C
C+  If no pion absorption, use Ejiri's paper (Phys.Rev C48, 1442(1993))
C
          if (shell.eq.0) then  ! Determine level of hole
              initial = 1
              xrndm = RLU(dummy)
              do while(xrndm.gt.cum_shell(initial))
                  initial = initial + 1
              enddo
          else
              initial = shell
          endif
C
C+  Short-cut: no de-excitation from p(1/2) state
C
          if (initial.eq.3) return
C
C+  Select final state based on flavor and level of hole
C
          rmode = RLU(dummy)
          if (idtarg .eq. proton_code) then
              cum = 0.
              mode = 1
          do while ((rmode*prob_shell(initial).gt.
     1                   cum+br_pmode(mode,initial)).and.
     1                  (mode.le.pshell_total))
                  cum = cum + br_pmode(mode,initial)
                  mode = mode + 1
              enddo
              if (mode.gt.pshell_total) then
                  write(lunout,*) 'nudeex: error locating proton mode'
                  return
              endif
              eg1 = e_pgamma(1,mode)
              eg2 = e_pgamma(2,mode)
              ep1 = e_pproton(1,mode)
              ep2 = e_pproton(2,mode)
              en  = e_pneutron(mode)
          else if (idtarg .eq. neutron_code) then
              cum = 0.
              mode = 1
          do while((rmode*prob_shell(initial).gt.
     1                    cum+br_nmode(mode,initial)).and.
     1                (mode.le.nshell_total))
                  cum = cum + br_nmode(mode,initial)
                  mode = mode + 1
              enddo
              if (mode.gt.nshell_total) then
                  write(lunout,*) 'nudeex: error location neutron mode'
                  return
              endif
              eg1 = e_ngamma(1,mode)
              eg2 = e_ngamma(2,mode)
              ep1 = e_nproton(1,mode)
              ep2 = e_nproton(2,mode)
              en  = e_nneutron(mode)
          endif
C
C+  Fill output particles
C
          if (eg1.gt.0.) then
              nout = nout+1
              idout(nout) = photon_code
              pout(4,nout) = eg1      
          endif
          if (eg2.gt.0.) then
              nout = nout+1
              idout(nout) = photon_code
              pout(4,nout) = eg2
          endif
          if (ep1.gt.0.) then
              nout = nout + 1
              idout(nout) = proton_code
              pout(4,nout) = proton_mass + ep1
          endif
          if (ep2.gt.0.) then
              nout = nout + 1
              idout(nout) = proton_code
              pout(4,nout) = proton_mass + ep2
          endif
          if (en.gt.0.) then
              nout = nout + 1
              idout(nout) = neutron_code
              pout(4,nout) = neutron_mass + en
          endif
      else
C
C+  Come here if a pion was absorbed.
C+  This code copied from NEUT(T.Kajita/Y.Hayato)
C
         etotal_gamma = 5. + 15.*RLU(dummy)
         if (etotal_gamma.lt.10.) then
            nout = min(int(RLU(dummy)*2 + 1), 2)
         else
            nout = min(int(RLU(dummy)*3+2),4)
            do while (etotal_gamma.gt.16 .and. nout.eq.1)
               nout = min(int(RLU(dummy)*3+2),4)
            enddo
         endif
         if (nout.eq.1) then
            pout(4,1) = etotal_gamma
            idout(1) = photon_code
         else
            do ip = 1, nout - 1
               ecum(ip) = etotal_gamma * RLU(dummy)
            enddo
            ecum(nout) = etotal_gamma
            call sortrq(ecum,1,nout,1)
            pout(4,1) = ecum(1)
            idout(1) = photon_code
            do ip = 2, nout
               pout(4,ip) = ecum(ip) - ecum(ip-1)
               idout(ip) = photon_code
            enddo
         endif
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
