      Subroutine nudeex(piabs,idtarg,shell,nout,idout,pout)
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
C+
C+  2011/03/09 : Now use nudeex_p and nudeex_n. ( Y. Hayato )
C
      implicit none
      logical piabs
      integer maxout
      parameter (maxout = 5)
      integer idtarg, shell, nout, idout(maxout)
      real pout(4,maxout)

      if (idtarg.eq.2112) then
         call nudeex_n(piabs,shell,nout,idout,pout)
      else if (idtarg.eq.2212) then
         call nudeex_p(piabs,shell,nout,idout,pout)
      else
         write(*,*) 'NUDEEX: Invalid target particle ID : ',idtarg
         stop
      endif

      end
