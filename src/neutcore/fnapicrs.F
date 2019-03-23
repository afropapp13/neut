************************************************************************
*     ---------------------------------------------
      FUNCTION FNAPICRS(E, IPAR, NTYPE, RP, R1320)
*     ---------------------------------------------
*
*     (Purpose)
*       MULTI PI INTERACTION CROSS SECTION
*       ( W>2.0GeV with N(pi)>=1 & 1.3<W<2.0GeV with N(pi)>=2)
*                 invariant mass threshold was changed from 1.4 GeV
*                 to 2.0 GeV due to improvement of Rein-Sehgal model.
*                 nu_e, anti-nu_e, nu_mu, and anti-nu_mu are supported.
*                 Following cross section tables were calculated by 
*                 integrating dndxdy.F function.
*                 Mar.-2-1998   by M.Shiozawa
*
*     (Input)
*       E      : NEUTRINO ENERGY ( GEV )
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR            
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR            
*                  16 : NEU-TAU       (not supported yet, Mar.-2-1998)
*                 -16 : NEU-TAU-BAR   (not supported yet, Mar.-2-1998)
*       NTYPE  : INTERACTION TYPE
*         =0   : NEUTRAL CURRENT
*         =1   : CHARGED CURRENT
*       RP     : FRACTION OF PROTON
*         RP < 0. -> DEFAULT VALUE 0.556
*
*     (Output)
*       FNAPICRS : Total multi-pi production cross section ( 10**-38cm^2 )
*                  ( W>2.0GeV with N(pi)>=1 & 1.3<W<2.0GeV with N(pi)>=2)
*       R1320    : Ratio of cross section
*                  [1.3<W<2.0GeV with N(pi)>=2] / [total cross section]
*
*     (Creation Date and Author)
*       1983.??.?? ; M.NAKAHATA  
*       1987.08.?? ; N.SATO    FOR TAU
*       1988.08.31 ; T.KAJITA  DATA UPDATE
*       1988.09.06 ; T.KAJITA  R1314 IS ADDED
*       1988.09.19 ; T.KAJITA  DX/DY WAS CHANGED BY THAT OF N.SATO'S
*                              WHICH INCLUDE LEPTON MASS TERM AND SMALL
*                              TERMS
*       1988.10.05 ; T.KAJITA  SIGMA(NC)/SIGMA(CC) RATIOS AT HIGH ENERGIES
*                              J.E.KIM ET AL., REV.MOD.PHYS.53(1981)211
*       1989.07.21 ; K.KANEYUKI NEU-TAU C.C. CROSS SECTION WAS UPDATED SAME 
*                               AS NEU-E,NEU-MU
*                               NEU-TAU N.C. CROSS SECTION => 
*                               SAME AS NEU-E,NEU-MU
*       1998.03.02 ; M.Shiozawa  invariant mass threshold was changed due
*                                to new improved Rein-Sehgal model.
*       1998.09.?? ; J.Kameda   New Calculation based on New Structure 
*                               function GRV94 DIS
*                               Consider Nu_tau cross section
*       2006.08.04 ; G.Mitsuka  Cross section is culculated after loading
*                               cross section table
*       2007.11.05 ; G.Mitsuka  support target nucleus besides 16O
*       2007.11.10 ; T.Tanaka   add upmu mode
*       2007.12.05 ; G.Mitsuka  Maximum neutrno energy is extended to 100TeV
*                               even if not upmu mode
*       2008.11.17 ; R.Tacik    calculate inump and inumpn for each event
*
************************************************************************
*
*   -- MEANING OF THE DATA --
*    CRSNE      : TOTAL CROSS SECTION FOR NEUTRINO (W > 2.0 GeV)
*               : (N,1) --- FOR ELECTRON NEUTRINO
*               : (N,2) --- FOR MUON NEUTRINO
*               : (N,3) --- FOR TAU  NEUTRINO
*
*    CRSNB      : TOTAL CROSS SECTION FOR ANTI NEUTRINO (W > 2.0 GeV)
*               : (N,1) --- FOR ELECTRON ANTI NEUTRINO
*               : (N,2) --- FOR MUON ANTI NEUTRINO
*               : (N,3) --- FOR TAU  ANTI NEUTRINO
*
*    C13NE      : TOTAL CROSS SECTION FOR NEUTRINO (W = 1.3-2.0 GeV)
*               :                                  (N(PI)>=2)
*               : (N,1) --- FOR ELECTRON NEUTRINO
*               : (N,2) --- FOR MUON NEUTRINO
*               : (N,3) --- FOR TAU  NEUTRINO
*
*    C13NB      : TOTAL CROSS SECTION FOR ANTI NEUTRINO (W=1.3-2.0 GeV)
*               :                                       (N(PI)>=2)
*               : (N,1) --- FOR ELECTRON ANTI NEUTRINO
*               : (N,2) --- FOR MUON ANTI NEUTRINO
*               : (N,3) --- FOR TAU  ANTI NEUTRINO
************************************************************************

C #include "discrs.h"
#include "necard.h"
#include "neutfilepath.h"

      integer iloadcrs
      save    iloadcrs
      data    iloadcrs/0/

      integer ipnrat, dum, lunpcrs, lunpncrs
      parameter (ipnrat=0)
      real inump, inumpn

      real e1320p1(61),m1320p1(61),t1320p1(61)
      real c1320p1(61,3)
      equivalence(c1320p1(1,1)  ,e1320p1(1))
      equivalence(c1320p1(1,2)  ,m1320p1(1))
      equivalence(c1320p1(1,3)  ,t1320p1(1))

      real eb1320p1(61),mb1320p1(61),tb1320p1(61)
      real cb1320p1(61,3)
      equivalence(cb1320p1(1,1) ,eb1320p1(1))
      equivalence(cb1320p1(1,2) ,mb1320p1(1))
      equivalence(cb1320p1(1,3) ,tb1320p1(1))

      real e1320pn1(61),m1320pn1(61),t1320pn1(61)
      real c1320pn1(61,3)
      equivalence(c1320pn1(1,1) ,e1320pn1(1))
      equivalence(c1320pn1(1,2) ,m1320pn1(1))
      equivalence(c1320pn1(1,3) ,t1320pn1(1))

      real eb1320pn1(61),mb1320pn1(61),tb1320pn1(61)
      real cb1320pn1(61,3)
      equivalence(cb1320pn1(1,1) ,eb1320pn1(1))
      equivalence(cb1320pn1(1,2) ,mb1320pn1(1))
      equivalence(cb1320pn1(1,3) ,tb1320pn1(1))

      real e1320p2(38),m1320p2(38),t1320p2(38)
      real c1320p2(38,3)
      equivalence(c1320p2(1,1)  ,e1320p2(1))
      equivalence(c1320p2(1,2)  ,m1320p2(1))
      equivalence(c1320p2(1,3)  ,t1320p2(1))

      real eb1320p2(38),mb1320p2(38),tb1320p2(38)
      real cb1320p2(38,3)
      equivalence(cb1320p2(1,1) ,eb1320p2(1))
      equivalence(cb1320p2(1,2) ,mb1320p2(1))
      equivalence(cb1320p2(1,3) ,tb1320p2(1))

      real e1320pn2(38),m1320pn2(38),t1320pn2(38)
      real c1320pn2(38,3)
      equivalence(c1320pn2(1,1) ,e1320pn2(1))
      equivalence(c1320pn2(1,2) ,m1320pn2(1))
      equivalence(c1320pn2(1,3) ,t1320pn2(1))

      real eb1320pn2(38),mb1320pn2(38),tb1320pn2(38)
      real cb1320pn2(38,3)
      equivalence(cb1320pn2(1,1) ,eb1320pn2(1))
      equivalence(cb1320pn2(1,2) ,mb1320pn2(1))
      equivalence(cb1320pn2(1,3) ,tb1320pn2(1))
      
      real e20p1(61),m20p1(61),t20p1(61)
      real c20p1(61,3)
      equivalence(c20p1(1,1)  ,e20p1(1))
      equivalence(c20p1(1,2)  ,m20p1(1))
      equivalence(c20p1(1,3)  ,t20p1(1))

      real eb20p1(61),mb20p1(61),tb20p1(61)
      real cb20p1(61,3)
      equivalence(cb20p1(1,1) ,eb20p1(1))
      equivalence(cb20p1(1,2) ,mb20p1(1))
      equivalence(cb20p1(1,3) ,tb20p1(1))

      real e20pn1(61),m20pn1(61),t20pn1(61)
      real c20pn1(61,3)
      equivalence(c20pn1(1,1) ,e20pn1(1))
      equivalence(c20pn1(1,2) ,m20pn1(1))
      equivalence(c20pn1(1,3) ,t20pn1(1))

      real eb20pn1(61),mb20pn1(61),tb20pn1(61)
      real cb20pn1(61,3)
      equivalence(cb20pn1(1,1),eb20pn1(1))
      equivalence(cb20pn1(1,2),mb20pn1(1))
      equivalence(cb20pn1(1,3),tb20pn1(1))

      real e20p2(38),m20p2(38),t20p2(38)
      real c20p2(38,3)
      equivalence(c20p2(1,1)  ,e20p2(1))
      equivalence(c20p2(1,2)  ,m20p2(1))
      equivalence(c20p2(1,3)  ,t20p2(1))

      real eb20p2(38),mb20p2(38),tb20p2(38)
      real cb20p2(38,3)
      equivalence(cb20p2(1,1) ,eb20p2(1))
      equivalence(cb20p2(1,2) ,mb20p2(1))
      equivalence(cb20p2(1,3) ,tb20p2(1))

      real e20pn2(38),m20pn2(38),t20pn2(38)
      real c20pn2(38,3)
      equivalence(c20pn2(1,1) ,e20pn2(1))
      equivalence(c20pn2(1,2) ,m20pn2(1))
      equivalence(c20pn2(1,3) ,t20pn2(1))

      real eb20pn2(38),mb20pn2(38),tb20pn2(38)
      real cb20pn2(38,3)
      equivalence(cb20pn2(1,1),eb20pn2(1))
      equivalence(cb20pn2(1,2),mb20pn2(1))
      equivalence(cb20pn2(1,3),tb20pn2(1))

      real e20p3(26),m20p3(26),t20p3(26)
      real c20p3(26,3)
      equivalence(c20p3(1,1)  ,e20p3(1))
      equivalence(c20p3(1,2)  ,m20p3(1))
      equivalence(c20p3(1,3)  ,t20p3(1))

      real eb20p3(26),mb20p3(26),tb20p3(26)
      real cb20p3(26,3)
      equivalence(cb20p3(1,1) ,eb20p3(1))
      equivalence(cb20p3(1,2) ,mb20p3(1))
      equivalence(cb20p3(1,3) ,tb20p3(1))

      real e20pn3(26),m20pn3(26),t20pn3(26)
      real c20pn3(26,3)
      equivalence(c20pn3(1,1) ,e20pn3(1))
      equivalence(c20pn3(1,2) ,m20pn3(1))
      equivalence(c20pn3(1,3) ,t20pn3(1))

      real eb20pn3(26),mb20pn3(26),tb20pn3(26)
      real cb20pn3(26,3)
      equivalence(cb20pn3(1,1),eb20pn3(1))
      equivalence(cb20pn3(1,2),mb20pn3(1))
      equivalence(cb20pn3(1,3),tb20pn3(1))

      real e20p4(11),m20p4(11),t20p4(11)
      real c20p4(11,3)
      equivalence(c20p4(1,1)  ,e20p4(1))
      equivalence(c20p4(1,2)  ,m20p4(1))
      equivalence(c20p4(1,3)  ,t20p4(1))

      real eb20p4(11),mb20p4(11),tb20p4(11)
      real cb20p4(11,3)
      equivalence(cb20p4(1,1) ,eb20p4(1))
      equivalence(cb20p4(1,2) ,mb20p4(1))
      equivalence(cb20p4(1,3) ,tb20p4(1))

      real e20pn4(11),m20pn4(11),t20pn4(11)
      real c20pn4(11,3)
      equivalence(c20pn4(1,1) ,e20pn4(1))
      equivalence(c20pn4(1,2) ,m20pn4(1))
      equivalence(c20pn4(1,3) ,t20pn4(1))

      real eb20pn4(11),mb20pn4(11),tb20pn4(11)
      real cb20pn4(11,3)
      equivalence(cb20pn4(1,1),eb20pn4(1))
      equivalence(cb20pn4(1,2),mb20pn4(1))
      equivalence(cb20pn4(1,3),tb20pn4(1))

      real e20p5(25),m20p5(25),t20p5(25)
      real c20p5(25,3)
      equivalence(c20p5(1,1)  ,e20p5(1))
      equivalence(c20p5(1,2)  ,m20p5(1))
      equivalence(c20p5(1,3)  ,t20p5(1))

      real eb20p5(25),mb20p5(25),tb20p5(25)
      real cb20p5(25,3)
      equivalence(cb20p5(1,1) ,eb20p5(1))
      equivalence(cb20p5(1,2) ,mb20p5(1))
      equivalence(cb20p5(1,3) ,tb20p5(1))

      real e20pn5(25),m20pn5(25),t20pn5(25)
      real c20pn5(25,3)
      equivalence(c20pn5(1,1) ,e20pn5(1))
      equivalence(c20pn5(1,2) ,m20pn5(1))
      equivalence(c20pn5(1,3) ,t20pn5(1))

      real eb20pn5(25),mb20pn5(25),tb20pn5(25)
      real cb20pn5(25,3)
      equivalence(cb20pn5(1,1),eb20pn5(1))
      equivalence(cb20pn5(1,2),mb20pn5(1))
      equivalence(cb20pn5(1,3),tb20pn5(1))

      real e20p6(36),m20p6(36),t20p6(36)
      real c20p6(36,3)
      equivalence(c20p6(1,1)  ,e20p6(1))
      equivalence(c20p6(1,2)  ,m20p6(1))
      equivalence(c20p6(1,3)  ,t20p6(1))

      real eb20p6(36),mb20p6(36),tb20p6(36)
      real cb20p6(36,3)
      equivalence(cb20p6(1,1) ,eb20p6(1))
      equivalence(cb20p6(1,2) ,mb20p6(1))
      equivalence(cb20p6(1,3) ,tb20p6(1))

      real e20pn6(36),m20pn6(36),t20pn6(36)
      real c20pn6(36,3)
      equivalence(c20pn6(1,1) ,e20pn6(1))
      equivalence(c20pn6(1,2) ,m20pn6(1))
      equivalence(c20pn6(1,3) ,t20pn6(1))

      real eb20pn6(36),mb20pn6(36),tb20pn6(36)
      real cb20pn6(36,3)
      equivalence(cb20pn6(1,1),eb20pn6(1))
      equivalence(cb20pn6(1,2),mb20pn6(1))
      equivalence(cb20pn6(1,3),tb20pn6(1))

      real e20p7(31),m20p7(31),t20p7(31)
      real c20p7(31,3)
      equivalence(c20p7(1,1)  ,e20p7(1))
      equivalence(c20p7(1,2)  ,m20p7(1))
      equivalence(c20p7(1,3)  ,t20p7(1))

      real eb20p7(31),mb20p7(31),tb20p7(31)
      real cb20p7(31,3)
      equivalence(cb20p7(1,1) ,eb20p7(1))
      equivalence(cb20p7(1,2) ,mb20p7(1))
      equivalence(cb20p7(1,3) ,tb20p7(1))

      real e20pn7(31),m20pn7(31),t20pn7(31)
      real c20pn7(31,3)
      equivalence(c20pn7(1,1) ,e20pn7(1))
      equivalence(c20pn7(1,2) ,m20pn7(1))
      equivalence(c20pn7(1,3) ,t20pn7(1))

      real eb20pn7(31),mb20pn7(31),tb20pn7(31)
      real cb20pn7(31,3)
      equivalence(cb20pn7(1,1),eb20pn7(1))
      equivalence(cb20pn7(1,2),mb20pn7(1))
      equivalence(cb20pn7(1,3),tb20pn7(1))

      real e20p8(26),m20p8(26),t20p8(26)
      real c20p8(26,3)
      equivalence(c20p8(1,1)  ,e20p8(1))
      equivalence(c20p8(1,2)  ,m20p8(1))
      equivalence(c20p8(1,3)  ,t20p8(1))

      real eb20p8(26),mb20p8(26),tb20p8(26)
      real cb20p8(26,3)
      equivalence(cb20p8(1,1) ,eb20p8(1))
      equivalence(cb20p8(1,2) ,mb20p8(1))
      equivalence(cb20p8(1,3) ,tb20p8(1))

      real e20pn8(26),m20pn8(26),t20pn8(26)
      real c20pn8(26,3)
      equivalence(c20pn8(1,1) ,e20pn8(1))
      equivalence(c20pn8(1,2) ,m20pn8(1))
      equivalence(c20pn8(1,3) ,t20pn8(1))

      real eb20pn8(26),mb20pn8(26),tb20pn8(26)
      real cb20pn8(26,3)
      equivalence(cb20pn8(1,1),eb20pn8(1))
      equivalence(cb20pn8(1,2),mb20pn8(1))
      equivalence(cb20pn8(1,3),tb20pn8(1))
  
      real crsne(2), csneh1(2), csneh2(2), csneh3(2), csneh4(2)
      real           csneh5(2), csneh6(2), csneh7(2)

      real crsnb(2), csnbh1(2), csnbh2(2), csnbh3(2), csnbh4(2)
      real           csnbh5(2), csnbh6(2), csnbh7(2)

      real c13ne(2), c13neh1(2),c13neh2(2),c13neh3(2),c13neh4(2)
      real           c13neh5(2),c13neh6(2),c13neh7(2)

      real c13nb(2), c13nbh1(2),c13nbh2(2),c13nbh3(2),c13nbh4(2)
      real           c13nbh5(2),c13nbh6(2),c13nbh7(2)


      CHARACTER*1024 FNAME
      CHARACTER*80   DUMSTR

      integer*4      LENSTR
      external       lenstr

      integer*4      lenpath

      lenpath = lenstr(CRSTBLPATH)

      if (iloadcrs.eq.0) then
         print*, 'Loading Cross section table for Multi-Pi productin'
         if(nepdf.eq.7.and.nebodek.eq.0) then ! GRV94 original
            print*, '---> Loading GRV94 orginal version'
            fname = crstblpath(1:lenpath)//'94org_p.dat'
            open(81,file=fname,form='formatted',status='old')
            fname = crstblpath(1:lenpath)//'94org_pn.dat'
            open(82,file=fname,form='formatted',status='old')
            lunpcrs = 81
            lunpncrs = 82
         elseif(nepdf.eq.7.and.nebodek.eq.1) then ! GRV94 Bodek corr.
            print*, '---> Loading GRV94 modified version'
            fname = crstblpath(1:lenpath)//'94mod_p.dat'
            open(83,file=fname,form='formatted',status='old')
            fname = crstblpath(1:lenpath)//'94mod_pn.dat'
            open(84,file=fname,form='formatted',status='old')
            lunpcrs = 83
            lunpncrs = 84
         elseif(nepdf.eq.12.and.nebodek.eq.0) then ! GRV98 original
            print*, '---> Loading GRV98 original version'
            fname = crstblpath(1:lenpath)//'98org_p.dat'
            open(85,file=fname,form='formatted',status='old')
            fname = crstblpath(1:lenpath)//'98org_pn.dat'
            open(86,file=fname,form='formatted',status='old')
            lunpcrs = 85
            lunpncrs = 86
         elseif(nepdf.eq.12.and.nebodek.eq.1) then ! GRV98 Bodek corr.
            print*, '---> Loading GRV98 modified version'
            fname = crstblpath(1:lenpath)//'98mod_p.dat'
            open(87,file=fname,form='formatted',status='old')
            fname = crstblpath(1:lenpath)//'98mod_pn.dat'
            open(88,file=fname,form='formatted',status='old')
            lunpcrs = 87
            lunpncrs = 88
         endif

         read(lunpcrs,'(A)') DUMSTR
         read(lunpncrs,'(A)') DUMSTR
         do i = 1,61 ! 0.5 GeV to 6.5 GeV (61bin, 1.3 < W < 2.0)
            read(lunpcrs,*) dum,e1320p1(i),m1320p1(i),t1320p1(i),
     &           eb1320p1(i),mb1320p1(i),tb1320p1(i)
            read(lunpncrs,*) dum,e1320pn1(i),m1320pn1(i),t1320pn1(i),
     &           eb1320pn1(i),mb1320pn1(i),tb1320pn1(i)
         enddo
         read(lunpcrs,'(A)') DUMSTR
         read(lunpncrs,'(A)') DUMSTR
         do i = 1,38 ! 6.5 GeV to 25.0 GeV (38bin, 1.3 < W < 2.0)
            read(lunpcrs,*) dum,e1320p2(i),m1320p2(i),t1320p2(i),
     &           eb1320p2(i),mb1320p2(i),tb1320p2(i)
            read(lunpncrs,*) dum,e1320pn2(i),m1320pn2(i),t1320pn2(i),
     &           eb1320pn2(i),mb1320pn2(i),tb1320pn2(i)
         enddo
         read(lunpcrs,'(A)') DUMSTR
         read(lunpncrs,'(A)') DUMSTR
         do i = 1,61 ! 0.5 GeV to 6.5 GeV (61bin, 2.0 < W)
            read(lunpcrs,*) dum,e20p1(i),m20p1(i),t20p1(i),
     &           eb20p1(i),mb20p1(i),tb20p1(i)
            read(lunpncrs,*) dum,e20pn1(i),m20pn1(i),t20pn1(i),
     &           eb20pn1(i),mb20pn1(i),tb20pn1(i)
         enddo
         read(lunpcrs,'(A)') DUMSTR 
         read(lunpncrs,'(A)') DUMSTR
         do i = 1,38 ! 6.5 GeV to 25.0 GeV (38bin, 2.0 < W)
            read(lunpcrs,*) dum,e20p2(i),m20p2(i),t20p2(i),
     &           eb20p2(i),mb20p2(i),tb20p2(i)
            read(lunpncrs,*) dum,e20pn2(i),m20pn2(i),t20pn2(i),
     &           eb20pn2(i),mb20pn2(i),tb20pn2(i)
         enddo
         read(lunpcrs,'(A)') DUMSTR
         read(lunpncrs,'(A)') DUMSTR
         do i = 1,26 ! 25.0 GeV to 50.0 GeV (26bin, 2.0 < W)
            read(lunpcrs,*) dum,e20p3(i),m20p3(i),t20p3(i),
     &           eb20p3(i),mb20p3(i),tb20p3(i)
            read(lunpncrs,*) dum,e20pn3(i),m20pn3(i),t20pn3(i),
     &           eb20pn3(i),mb20pn3(i),tb20pn3(i)
         enddo
         read(lunpcrs,'(A)') DUMSTR
         read(lunpncrs,'(A)') DUMSTR
         do i = 1,11 ! 50.0 GeV to 100 GeV (11bin, 2.0 < W)
            read(lunpcrs,*) dum,e20p4(i),m20p4(i),t20p4(i),
     &           eb20p4(i),mb20p4(i),tb20p4(i)
            read(lunpncrs,*) dum,e20pn4(i),m20pn4(i),t20pn4(i),
     &           eb20pn4(i),mb20pn4(i),tb20pn4(i)
         enddo
         read(lunpcrs,'(A)') DUMSTR
         read(lunpncrs,'(A)') DUMSTR
         do i = 1,25 ! 100 GeV to 2500 GeV (25bin, 2.0 < W)
            read(lunpcrs,*) dum,e20p5(i),m20p5(i),t20p5(i),
     &           eb20p5(i),mb20p5(i),tb20p5(i)
            read(lunpncrs,*) dum,e20pn5(i),m20pn5(i),t20pn5(i),
     &           eb20pn5(i),mb20pn5(i),tb20pn5(i)
         enddo
         
         read(lunpcrs,'(A)') DUMSTR
         read(lunpncrs,'(A)') DUMSTR
         do i = 1,36 ! 2.5TeV to 20TeV (36bin, 2.0 < W)
            read(lunpcrs,*) dum,e20p6(i),m20p6(i),t20p6(i),
     &           eb20p6(i),mb20p6(i),tb20p6(i)
            read(lunpncrs,*) dum,e20pn6(i),m20pn6(i),t20pn6(i),
     &           eb20pn6(i),mb20pn6(i),tb20pn6(i)
         enddo
         read(lunpcrs,'(A)') DUMSTR
         read(lunpncrs,'(A)') DUMSTR
         do i = 1,31 ! 20 TeV to 50TeV (31bin, 2.0 < W)
            read(lunpcrs,*) dum,e20p7(i),m20p7(i),t20p7(i),
     &           eb20p7(i),mb20p7(i),tb20p7(i)
            read(lunpncrs,*) dum,e20pn7(i),m20pn7(i),t20pn7(i),
     &           eb20pn7(i),mb20pn7(i),tb20pn7(i)
         enddo
         read(lunpcrs,'(A)') DUMSTR
         read(lunpncrs,'(A)') DUMSTR
         do i = 1,26 ! 50TeV to 100TeV (26bin, 2.0 < W)
            read(lunpcrs,*) dum,e20p8(i),m20p8(i),t20p8(i),
     &           eb20p8(i),mb20p8(i),tb20p8(i)
            read(lunpncrs,*) dum,e20pn8(i),m20pn8(i),t20pn8(i),
     &           eb20pn8(i),mb20pn8(i),tb20pn8(i)
         enddo



c$$$         if(ipnrat.eq.0) then ! p/n=10:8 (i.e. H2O)
c$$$            print*, 'p/n ratio is 10:8'
c$$$            inump =  2./18.
c$$$            inumpn = 8./9.
c$$$         elseif(ipnrat.eq.1) then ! p/n=1:1 (i.e. O2)
c$$$            print*, 'p/n ratio is 1:1'
c$$$            inump =  0.
c$$$            inumpn = 1.            
c$$$         endif
C---  p/n ratio is calculated by NUMBNDP,NUMFREP, and so on (2007 Nov G.Mitsuka)
         
         iloadcrs = 1
         
      endif

*     Calculate the contrubution of each cross section table
         if (RP.eq.0.5) then
            inump = 0.
            inumpn = 1.
c            print*, 'p/n ratio is 1:1'
         else
            inump  = float(NUMBNDP+NUMFREP-NUMBNDN)/float(NUMBNDP+NUMFREP+NUMBNDN)
            inumpn = 2.*float(NUMBNDN)/float(NUMBNDP+NUMFREP+NUMBNDN)
c            print*, 'p/n ratio is ',NUMBNDP+NUMFREP,' :',NUMBNDN
c            print*, 'table contribution is ',inump,' :',inumpn
         endif

c         iloadcrs = 1
c
c      endif

C  DETERMINE FRACTION OF PROTON, AND CALCULATE CORRECTION FACTOR
C     (DEVIATION FROM ISOSCALAR) OF THE CROSS SECTION        31-AUG-88
      RPR=RP
C     RPR is given by the ratio of (total proton)/(total nucleon) (2007 Nov, G.Mitsuka)
      IF (RPR.LT.0.) RPR=float(NUMBNDP+NUMFREP)/float(NUMBNDP+NUMFREP+NUMBNDN)
      IF (RPR.GT.1.) RPR=1.0
      IF (IPAR.GT.0) CORREC=RPR*(2./3.)+(1.-RPR)*(4./3.)
      IF (IPAR.LT.0) CORREC=RPR*(4./3.)+(1.-RPR)*(2./3.)

C maybe.. ( Y.Hayato 2008/08/14 )
      correc = 1. 

      FNAPICRS=0.
      R1320=0.
C      If(ABS(IPAR).EQ.16 .AND. NTYPE.EQ.1) GO TO 7000

      IF(E.GE.6.5 .AND. E.LT.25. .AND. IPAR.GT.0)GO TO 1000
      IF(E.GE.6.5 .AND. E.LT.25. .AND. IPAR.LT.0)GO TO 1100

      IF(E.GE.25. .AND. E.LT.50. .AND. IPAR.GT.0)GO TO 2000
      IF(E.GE.25. .AND. E.LT.50. .AND. IPAR.LT.0)GO TO 2100

      IF(E.GE.50. .AND. E.LT.100..AND. IPAR.GT.0)GO TO 3000
      IF(E.GE.50. .AND. E.LT.100..AND. IPAR.LT.0)GO TO 3100

      IF(E.GE.100.  .AND. E.LT.2500. .AND. IPAR.GT.0)GO TO 4000
      IF(E.GE.100.  .AND. E.LT.2500. .AND. IPAR.LT.0)GO TO 4100

      IF(E.GE.2500. .AND. E.LT.20000. .AND. IPAR.GT.0)GO TO 5000
      IF(E.GE.2500. .AND. E.LT.20000. .AND. IPAR.LT.0)GO TO 5100

      IF(E.GE.20000..AND. E.LT.50000. .AND. IPAR.GT.0)GO TO 6000
      IF(E.GE.20000..AND. E.LT.50000. .AND. IPAR.LT.0)GO TO 6100

      IF(E.GE.50000..AND. IPAR.GT.0)GO TO 7000
      IF(E.GE.50000..AND. IPAR.LT.0)GO TO 7100

      IF(E.GE.100. .AND. IPAR.GT.0)GO TO 4000
      IF(E.GE.100. .AND. IPAR.LT.0)GO TO 4100

      IF(E.LE.0.5)RETURN
      IBIN=IFIX((E-0.5)/0.1)+1
      EBIN=0.5+FLOAT(IBIN-1)*0.1
      IF(IBIN.LT.61)GO TO 10
      IBIN=60

   10 IL=(IABS(IPAR)-10)/2
      IF(NTYPE.EQ.0)IL=1
      IF(IPAR.LT.0)GO TO 500
C     
C     -- NEUTRINO CASE (E < 6.5)
C     
      crsne(1) = inump*c20p1(ibin,il)   + inumpn*c20pn1(ibin,il)
      crsne(2) = inump*c20p1(ibin+1,il) + inumpn*c20pn1(ibin+1,il)

      FNAPICRS=CRSNE(1)+(CRSNE(2)-CRSNE(1))*
     &     (E-EBIN)/0.1
C     ADD  ( 1.3 < W < 2.0 GeV , with N(pi)>=2 )

      c13ne(1) = inump*c1320p1(ibin,il)   + inumpn*c1320pn1(ibin,il)
      c13ne(2) = inump*c1320p1(ibin+1,il) + inumpn*c1320pn1(ibin+1,il)

      ADDCRS=C13NE(1)+(C13NE(2)-C13NE(1))*
     &     (E-EBIN)/0.1
      FNAPICRS=FNAPICRS+ADDCRS
      IF (FNAPICRS.GT.1.0E-5) then
         R1320= ADDCRS/FNAPICRS
      else
         adcrs = 0.
         fnapicrs = 0.
      endif
      GO TO 330
C     
C     -- ANTI NEUTRINO CASE (E < 6.5)
C     
  500 continue
      CRSNB(1) = inump*cb20p1(ibin  ,il) + inumpn*cb20pn1(ibin  ,il)
      CRSNB(2) = inump*cb20p1(ibin+1,il) + inumpn*cb20pn1(ibin+1,il)
      
      FNAPICRS=CRSNB(1)+(CRSNB(2)-CRSNB(1))*
     &     (E-EBIN)/0.1
C     ADD  ( 1.3 < W < 2.0 GeV , with N(pi)>=2 )

      C13NB(1) = inump*cb1320p1(ibin  ,il)+inumpn*cb1320pn1(ibin  ,il)
      C13NB(2) = inump*cb1320p1(ibin+1,il)+inumpn*cb1320pn1(ibin+1,il)
      ADDCRS=C13NB(1)+(C13NB(2)-C13NB(1))*
     &     (E-EBIN)/0.1
      FNAPICRS=FNAPICRS+ADDCRS
      IF (FNAPICRS.GT.1.0E-5) R1320= ADDCRS/FNAPICRS
      GO TO 330

C     
C     -- NEUTRINO CASE FOR ( 6.5 < E < 25.)
C     
 1000 CONTINUE
 1100 IBIN=IFIX((E-6.5)/0.5)+1
      EBIN=6.5+FLOAT(IBIN-1)*0.5
      IF(IBIN.LT.38)GO TO 1200
      IBIN=38

 1200 IL=(IABS(IPAR)-10)/2
      IF(NTYPE.EQ.0)IL=1
      IF(IPAR.LT.0)GO TO 1500
C     
C     -- NEUTRINO CASE E
C     
      csneh1(1)=inump*c20p2(ibin  ,il) + inumpn*c20pn2(ibin  ,il)
      csneh1(2)=inump*c20p2(ibin+1,il) + inumpn*c20pn2(ibin+1,il)

      FNAPICRS=CSNEH1(1)+(CSNEH1(2)-CSNEH1(1))*
     &     (E-EBIN)/0.5
C     ADD  ( 1.3 < W < 2.0 GeV , with N(pi)>=2 )
C      ADDCRS=C13NE(60,IL)

      c13neh1(1)=inump*c1320p2(ibin  ,il)+inumpn*c1320pn2(ibin  ,il)
      c13neh1(2)=inump*c1320p2(ibin+1,il)+inumpn*c1320pn2(ibin+1,il)

      ADDCRS=C13NEH1(1)+(C13NEH1(2)-C13NEH1(1))*
     &     (E-EBIN)/0.5
      FNAPICRS=FNAPICRS+ADDCRS
      IF (FNAPICRS.GT.1.0E-5) R1320= ADDCRS/FNAPICRS
      GO TO 330
C     
C     -- ANTI NEUTRINO CASE
C     
 1500 continue
      csnbh1(1)=inump*cb20p2(ibin  ,il) + inumpn*cb20pn2(ibin  ,il)
      csnbh1(2)=inump*cb20p2(ibin+1,il) + inumpn*cb20pn2(ibin+1,il)

      FNAPICRS=CSNBH1(1)+(CSNBH1(2)-CSNBH1(1))*
     &     (E-EBIN)/0.5
C     ADD  ( 1.3 < W < 2.0 GeV , with N(pi)>=2 )
C      ADDCRS=C13NB(60,IL)

      c13nbh1(1)=inump*cb1320p2(ibin  ,il)+inumpn*cb1320pn2(ibin  ,il)
      c13nbh1(2)=inump*cb1320p2(ibin+1,il)+inumpn*cb1320pn2(ibin+1,il)

      ADDCRS=C13NBH1(1)+(C13NBH1(2)-C13NBH1(1))*
     &     (E-EBIN)/0.5
      FNAPICRS=FNAPICRS+ADDCRS
      IF (FNAPICRS.GT.1.0E-5) R1320= ADDCRS/FNAPICRS
      GO TO 330


C     
C     -- NEUTRINO CASE FOR ( 25. < E < 50.)
C     
 2000 CONTINUE
 2100 IBIN=IFIX((E-25)/1.0)+1
      EBIN=25.0+FLOAT(IBIN-1)*1.0
      IF(IBIN.LT.26)GO TO 2200
      IBIN=26

 2200 IL=(IABS(IPAR)-10)/2
      IF(NTYPE.EQ.0)IL=1
      IF(IPAR.LT.0)GO TO 2500
C     
C     -- NEUTRINO CASE E
C     
      CSNEH2(1) = inump*c20p3(ibin  ,il)+inumpn*c20pn3(ibin  ,il)
      CSNEH2(2) = inump*c20p3(ibin+1,il)+inumpn*c20pn3(ibin+1,il)

      FNAPICRS=CSNEH2(1)+(CSNEH2(2)-CSNEH2(1))*
     &     (E-EBIN)/1.0
C     ADD  ( 1.3 < W < 2.0 GeV , with N(pi)>=2 )
C      ADDCRS=C13NE(60,IL)

      C13NEH1(1) = inump*c1320p2(37,il) + inumpn*c1320pn2(37,il)

      ADDCRS=C13NEH1(1)
      FNAPICRS=FNAPICRS+ADDCRS
      IF (FNAPICRS.GT.1.0E-5) R1320= ADDCRS/FNAPICRS
      GO TO 330
C     
C     -- ANTI NEUTRINO CASE
C     
 2500 continue
      
      CSNBH2(1) = inump*cb20p3(ibin  ,il)+inumpn*cb20pn3(ibin  ,il)
      CSNBH2(2) = inump*cb20p3(ibin+1,il)+inumpn*cb20pn3(ibin+1,il)

      FNAPICRS=CSNBH2(1)+(CSNBH2(2)-CSNBH2(1))*
     &     (E-EBIN)/1.0
C     ADD  ( 1.3 < W < 2.0 GeV , with N(pi)>=2 )
C      ADDCRS=C13NB(60,IL)

      C13NBH1(1) = inump*cb1320p2(37,il)+inumpn*cb1320pn2(37,il)

      ADDCRS=C13NBH1(1)
      FNAPICRS=FNAPICRS+ADDCRS
      IF (FNAPICRS.GT.1.0E-5) R1320= ADDCRS/FNAPICRS
      GO TO 330

C     
C     -- NEUTRINO CASE FOR ( 50. < E < 100.)
C     
 3000 CONTINUE
 3100 IBIN=IFIX((E-50.)/5.0)+1
      EBIN=50.0+FLOAT(IBIN-1)*5.0
      IF(IBIN.LT.11)GO TO 3200
      IBIN=11

 3200 IL=(IABS(IPAR)-10)/2
      IF(NTYPE.EQ.0)IL=1
      IF(IPAR.LT.0)GO TO 3500
C     
C     -- NEUTRINO CASE E
C     
      CSNEH3(1) = inump*c20p4(ibin  ,il) + inumpn*c20pn4(ibin  ,il)
      CSNEH3(2) = inump*c20p4(ibin+1,il) + inumpn*c20pn4(ibin+1,il)

      FNAPICRS=CSNEH3(1)+(CSNEH3(2)-CSNEH3(1))*
     &     (E-EBIN)/5.0
C     ADD  ( 1.3 < W < 2.0 GeV , with N(pi)>=2 )
C      ADDCRS=C13NE(60,IL)

      C13NEH1(1) = inump*c1320p2(37,il) + inumpn*c1320pn2(37,il)

      ADDCRS=C13NEH1(1)
      FNAPICRS=FNAPICRS+ADDCRS
      IF (FNAPICRS.GT.1.0E-5) R1320= ADDCRS/FNAPICRS
      GO TO 330
C     
C     -- ANTI NEUTRINO CASE
C     
 3500 continue

      CSNBH3(1) = inump*cb20p4(ibin  ,il)+inumpn*cb20pn4(ibin  ,il)
      CSNBH3(2) = inump*cb20p4(ibin+1,il)+inumpn*cb20pn4(ibin+1,il)

      FNAPICRS=CSNBH3(1)+(CSNBH3(2)-CSNBH3(1))*
     &     (E-EBIN)/5.0
C     ADD  ( 1.3 < W < 2.0 GeV , with N(pi)>=2 )
C      ADDCRS=C13NB(60,IL)

      C13NBH1(1) = inump*cb1320p2(37,il)+inumpn*cb1320pn2(37,il)

      ADDCRS=C13NBH1(1)
      FNAPICRS=FNAPICRS+ADDCRS
      IF (FNAPICRS.GT.1.0E-5) R1320= ADDCRS/FNAPICRS
      GO TO 330
      

C     
C     -- NEUTRINO CASE FOR ( 100. < E < 2500.)
C     
 4000 CONTINUE
 4100 IBIN=IFIX((E-100.)/100.0)+1
      EBIN=100.0+FLOAT(IBIN-1)*100.0
      IF(IBIN.LT.25)GO TO 4200
      IBIN=25

 4200 IL=(IABS(IPAR)-10)/2
      IF(NTYPE.EQ.0)IL=1
      IF(IPAR.LT.0)GO TO 4500
C     
C     -- NEUTRINO CASE E
C     
      CSNEH4(1) = inump*c20p5(ibin  ,il)+inumpn*c20pn5(ibin  ,il)
      CSNEH4(2) = inump*c20p5(ibin+1,il)+inumpn*c20pn5(ibin+1,il)

      FNAPICRS=CSNEH4(1)+(CSNEH4(2)-CSNEH4(1))*
     &     (E-EBIN)/100.
C     ADD  ( 1.3 < W < 2.0 GeV , with N(pi)>=2 )
C      ADDCRS=C13NE(60,IL)

      C13NEH1(1) = inump*c1320p2(37,il)+inumpn*c1320pn2(37,il)

      ADDCRS=C13NEH1(1)
      FNAPICRS=FNAPICRS+ADDCRS
      IF (FNAPICRS.GT.1.0E-5) R1320= ADDCRS/FNAPICRS
      GO TO 330
C     
C     -- ANTI NEUTRINO CASE
C     
 4500 continue

      CSNBH4(1) = inump*cb20p5(ibin  ,il)+inumpn*cb20pn5(ibin  ,il)
      CSNBH4(2) = inump*cb20p5(ibin+1,il)+inumpn*cb20pn5(ibin+1,il)
      
      FNAPICRS=CSNBH4(1)+(CSNBH4(2)-CSNBH4(1))*
     &     (E-EBIN)/100.
C     ADD  ( 1.3 < W < 2.0 GeV , with N(pi)>=2 )
C      ADDCRS=C13NB(60,IL)

      C13NBH1(1) = inump*cb1320p2(37,il)+inumpn*cb1320pn2(37,il)

      ADDCRS=C13NBH1(1)
      FNAPICRS=FNAPICRS+ADDCRS
      IF (FNAPICRS.GT.1.0E-5) R1320= ADDCRS/FNAPICRS
      GO TO 330


C
C     -- NEUTRINO CASE FOR ( 2500. < E < 20000.)
C
 5000   CONTINUE
 5100      IBIN=IFIX((E-2500.)/500.0)+1
      EBIN=2500.0+FLOAT(IBIN-1)*500.0
      IF(IBIN.LT.36)GO TO 5200
      IBIN=36

 5200   IL=(IABS(IPAR)-10)/2
      IF(NTYPE.EQ.0)IL=1
      IF(IPAR.LT.0)GO TO 5500
C
C     -- NEUTRINO CASE E
C
      CSNEH5(1) = inump*c20p6(ibin  ,il)+inumpn*c20pn6(ibin  ,il)
      CSNEH5(2) = inump*c20p6(ibin+1,il)+inumpn*c20pn6(ibin+1,il)

      FNAPICRS=CSNEH5(1)+(CSNEH5(2)-CSNEH5(1))*
     &     (E-EBIN)/500.
C     ADD  ( W > 1.3 - 1.4 GEV , WITH N(PAI)=2 )  06-SEP-88
C     ADDCRS=C13NE(60,IL)
C      C13NEH1(1) = inump*c1320p1(60,il) + inumpn*c1320pn1(60,il)

      C13NEH1(1) = inump*c1320p2(37,il) + inumpn*c1320pn2(37,il)

      ADDCRS=C13NEH1(1)
      FNAPICRS=FNAPICRS+ADDCRS
      R1314= ADDCRS/FNAPICRS
      GO TO 330
C
C     -- ANTI NEUTRINO CASE
C
 5500 continue

      CSNBH5(1) = inump*cb20p6(ibin  ,il)+inumpn*cb20pn6(ibin  ,il)
      CSNBH5(2) = inump*cb20p6(ibin+1,il)+inumpn*cb20pn6(ibin+1,il)

      FNAPICRS=CSNBH5(1)+(CSNBH5(2)-CSNBH5(1))*
     &     (E-EBIN)/500.
C     ADD  ( W > 1.3 - 1.4 GEV , WITH N(PAI)=2 )  06-SEP-88
C     ADDCRS=C13NB(60,IL)
C      C13NBH1(1) = inump*cb1320p1(60,il)+inumpn*cb1320pn1(60,il)

      C13NBH1(1) = inump*cb1320p2(37,il)+inumpn*cb1320pn2(37,il)

      ADDCRS=C13NBH1(1)
      FNAPICRS=FNAPICRS+ADDCRS
      R1314= ADDCRS/FNAPICRS
      GO TO 330

C
C     -- NEUTRINO CASE FOR ( 20000. < E < 50000.)
C
 6000   CONTINUE
 6100      IBIN=IFIX((E-20000.)/1000.0)+1
      EBIN=20000.0+FLOAT(IBIN-1)*1000.0
      IF(IBIN.LT.31)GO TO 6200
      IBIN=31

 6200   IL=(IABS(IPAR)-10)/2
      IF(NTYPE.EQ.0)IL=1
      IF(IPAR.LT.0)GO TO 6500
C
C     -- NEUTRINO CASE E
C
      CSNEH6(1) = inump*c20p7(ibin  ,il)+inumpn*c20pn7(ibin  ,il)
      CSNEH6(2) = inump*c20p7(ibin+1,il)+inumpn*c20pn7(ibin+1,il)

      FNAPICRS=CSNEH6(1)+(CSNEH6(2)-CSNEH6(1))*
     &     (E-EBIN)/1000.
C     ADD  ( W > 1.3 - 1.4 GEV , WITH N(PAI)=2 )  06-SEP-88
C     ADDCRS=C13NE(60,IL)
C      C13NEH1(1) = inump*c1320p1(60,il) + inumpn*c1320pn1(60,il)

      C13NEH1(1) = inump*c1320p2(37,il) + inumpn*c1320pn2(37,il)

      ADDCRS=C13NEH1(1)
      FNAPICRS=FNAPICRS+ADDCRS
      R1314= ADDCRS/FNAPICRS
      GO TO 330
C
C     -- ANTI NEUTRINO CASE
C
 6500 continue

      CSNBH6(1) = inump*cb20p7(ibin  ,il)+inumpn*cb20pn7(ibin  ,il)
      CSNBH6(2) = inump*cb20p7(ibin+1,il)+inumpn*cb20pn7(ibin+1,il)

      FNAPICRS=CSNBH6(1)+(CSNBH6(2)-CSNBH6(1))*
     &     (E-EBIN)/1000.
C     ADD  ( W > 1.3 - 1.4 GEV , WITH N(PAI)=2 )  06-SEP-88
C     ADDCRS=C13NB(60,IL)
C      C13NBH1(1) = inump*cb1320p1(60,il)+inumpn*cb1320pn1(60,il)

      C13NBH1(1) = inump*cb1320p2(37,il) + inumpn*cb1320pn2(37,il)

      ADDCRS=C13NBH1(1)
      FNAPICRS=FNAPICRS+ADDCRS
      R1314= ADDCRS/FNAPICRS
      GO TO 330
C
C     -- NEUTRINO CASE FOR ( 50000. < E < 100000.)
C
 7000   CONTINUE
 7100      IBIN=IFIX((E-50000.)/2000.0)+1
      EBIN=50000.0+FLOAT(IBIN-1)*2000.0
      IF(IBIN.LT.25)GO TO 7200
      IBIN=25

 7200   IL=(IABS(IPAR)-10)/2
      IF(NTYPE.EQ.0)IL=1
      IF(IPAR.LT.0)GO TO 7500
C
C     -- NEUTRINO CASE E
C
      CSNEH7(1) = inump*c20p8(ibin  ,il)+inumpn*c20pn8(ibin  ,il)
      CSNEH7(2) = inump*c20p8(ibin+1,il)+inumpn*c20pn8(ibin+1,il)

      FNAPICRS=CSNEH7(1)+(CSNEH7(2)-CSNEH7(1))*
     &     (E-EBIN)/2000.
C     ADD  ( W > 1.3 - 1.4 GEV , WITH N(PAI)=2 )  06-SEP-88
C     ADDCRS=C13NE(60,IL)
C      C13NEH1(1) = inump*c1320p1(60,il) + inumpn*c1320pn1(60,il)      

      C13NEH1(1) = inump*c1320p2(37,il) + inumpn*c1320pn2(37,il)

      ADDCRS=C13NEH1(1)
      FNAPICRS=FNAPICRS+ADDCRS
      R1314= ADDCRS/FNAPICRS
      GO TO 330
C
C     -- ANTI NEUTRINO CASE
C
 7500 continue
      CSNBH7(1) = inump*cb20p8(ibin  ,il)+inumpn*cb20pn8(ibin  ,il)
      CSNBH7(2) = inump*cb20p8(ibin+1,il)+inumpn*cb20pn8(ibin+1,il)

      FNAPICRS=CSNBH7(1)+(CSNBH7(2)-CSNBH7(1))*
     &     (E-EBIN)/1000.
C     ADD  ( W > 1.3 - 1.4 GEV , WITH N(PAI)=2 )  06-SEP-88
C     ADDCRS=C13NB(60,IL)
C      C13NBH1(1) = inump*cb1320p1(60,il)+inumpn*cb1320pn1(60,il)

      C13NBH1(1) = inump*cb1320p2(37,il) + inumpn*cb1320pn2(37,il)

      ADDCRS=C13NBH1(1)
      FNAPICRS=FNAPICRS+ADDCRS
      R1314= ADDCRS/FNAPICRS
      GO TO 330


  330 IF(NTYPE.EQ.0)THEN
C     MODIFICATION  T.KAJITA  05-OCT-88
C     SIGMA(NC)/SIGMA(CC) RATIOS AT HIGH ENERGIES
C     J.E.KIM ET AL., REV.MOD.PHYS.53(1981)211
C     -- IF NEUTRAL CURRENT CASE MULTIPLY 0.26(NEUT, <3GEV)
C     --                                  0.30(NEUT, >6GEV)
C     --                                  0.39(NEUBAR, <3GEV)
C     --                                  0.37(NEUBAR, >6GEV)
         IF(IPAR.GT.0) THEN
            FACTOR=0.26
            IF(E.GT.6.)FACTOR=0.30
            IF(E.GT.3..AND.E.LE.6.)
     &           FACTOR= 0.26+(0.30-0.26)*(E-3.)/(6.-3.)
            FNAPICRS=FNAPICRS*FACTOR
         END IF
         IF(IPAR.LT.0) THEN
            FACTOR=0.39
            IF(E.GT.6.)FACTOR=0.37
            IF(E.GT.3..AND.E.LE.6.)
     &           FACTOR= 0.39+(0.37-0.39)*(E-3.)/(6.-3.)
            FNAPICRS=FNAPICRS*FACTOR
         END IF
      ELSE
         FNAPICRS=FNAPICRS*CORREC
      END IF
      RETURN
      END
