************************************************************************
*     --------
*     NECARD.H
*     --------
*
*     (Purpose)
*       COMMON for CARD on NEUT
*
*     (Variables)
*       NEFRMFLG : ( NEUT-FERM ) Flag on Fermi motion
*                     0 : on ( default ) 
*                     1 : off
*       NEPAUFLG : ( NEUT-PAUL ) Flag on Pauli brokking
*                     0 : on ( default ) 
*                     1 : off
*       NENEFO16 : ( NEUT-NEFF ) Flag on Nuclear effect in O16
*                     0 : on ( default ) 
*                     1 : off
*       NENEFMODL: ( NEUT-MODL ) Select model for low energy pion
*                     0 : default (Salcedo et. al)
*                     1 : Tuned to pion scattering data
*       NENEFMODH: ( NEUT-MODH ) Select model for high energy pion
*                     0 : default
*                     1 : proton/neutron separated
*       NENEFKINH: ( NEUT-KINH ) Select kinematical model for 
*                                   hi-nrg pion inelastic scattering
*                     0 : Isotropic resonance decay
*                     1 : SAID WI08 Phase Shifts
*       neabspiemit : ( NEUT-neabspiemit ) Flag on Nucleon Ejection
*                     0 : off 
*                     1 : on ( default ) 
*       NUSIM       : Neutrino Simulation Flag
*                     1 : Neutrino Simulation   (default)
*                     0 : Other (piscat, gampi)
*
*       NEFKINVER   : VERSION of pion nuclear effects kinematics
*                     0 : Original model ( default )
*                       :  some part uses measured momentum spectrum of proton
*                     1 : New model ( completely LFG for pion )
*
*       NEMODFLG : ( NEUT-MODE ) Flag on Interaction mode on neutrino int.
*                     0 : normal ( default )
*                    -1 : input cross section by CRSNEUT
*                     n : sellect one mode ( n > 0 )   See nemodsel.F
*                           n =  1 : charged current Q.E. 
*                           n = 11,12,13
*                                  : charged current Single pi production 
*                           n = 21 : charged current Multi pi production
*                           n = 31,32,33,34
*                                  : neutral current Single pi production 
*                           n = 41 : neutral current Multi pi production
*                           n = 51,52 : neutral current elastic
*                           n = 22,42,43 : single eta production 
*                           n = 23,44,45 : single  K  production 
* 
*       CRSNEUT(30)   : ( NEUT-CRS ) Multiplied factor to cross section
*                                    on each mode.   See nemodsel.F
*       CRSNEUTB(30)  : ( NEUT-CRSB ) Multiplied factor to cross section
*                                    on each mode.   See nemodsel.F
*
*       NECOHEPI      : ( NEUT_COHEPI ) Select Coherent pi model
*                            0 : Rein & Sehgal
*                            1 : Kartavtsev 
*                            2 : Berger & Sehgal
*
*       NEDIFPI      : ( NEUT_DIFPI ) Select Diffractive pi model
*                            0 : Rein
*
*       NEPDF         : ( NEUT-Select Parton distribution function
*                           n = 7 : GRV94
*                           n =12 : GRV98
*
*       NEBODEK       : ( NEUT- turn off/on Bodek-Yang correction
*                       0 :   off
*                       1 :   on - axial/vector unified, 2005 Bodek-Yang Modification
*       		2 :   on - with new A-C parameters in the 2019 Bodek-Yang paper, axial part seperated, scaling factor for xf3.
*                       3 :   on - removing the new H(x,Q2) scaling factor for xf3 from nebodek=2
*	          	4 :   on - further removing axial/vector separation from nebodek=3
*                       below are options that should not be used at present
*                       5 :   on - Add in a scaling factor K_LW for vector valence K(u/d)
*                       6 :   on - All updates and enable charm production
*
*       NEMULT        : (Hadron multiplicity model for multi-pi mode)
*                       0 :  old NEUT model
*			1 :  fit of deuterium bubble-chamber data (hep-ph:1607.06558)
*		        2 :  AGKY model (hep-ph:0904.4043 and Yang, T. et al. AIP Conf.Proc. 967 (2007) 269-275)
*
*       ITAUFLGCORE  : control Tau run mode(this is not set here)
*
*       NUMBNDN  : total number of neutron      (e.g. H2O => 8 , Fe => 30)
*       NUMBNDP  : total number of bound proton (e.g. H2O => 8 , Fe => 26)
*       NUMFREP  : total number of free proton  (e.g. H2O => 2 , Fe =>  0)
*       NUMATOM  : atomic number                (e.g.   O =>16 , Fe => 56)
C
C       IPILESSDCY: Pion less delta deacy             ( 1 : on / 0: off )
C       RPILESSDCY: Fraction of pion less delta deacy ( default 0.2 )
C
*       NEIFF:    Tells which form factors are used (see rsdcrs)
*       NENRTYPE: For R-S form factors, there is two separate
*
C       IRADCORR:  Radiative correction on/off ( 1 : on / 0: off )
*                 ( Currently off by default )
*
C       QUIET : Screen output verbosity
C         0 : Default (prints all initial state info)
C         1 : Print only neutrino energy
C         2 : Prints almost nothing (except PYTHIA output)
C
*
*     
*     (Creation Date and Author)
*       1995.11.17 ; K.Kaneyuki
*       1997.12.01 ; J.Kameda modify to consider eta
*       1998.03.01 ; J.Kameda modify to consider  K
*       2006.??.?? ; G.Mitsuka added NEPDF & NEBODEK 
*       2006.12.30 ; Y.Hayato move flux /detector dependent part to
*                              necardfx.h
*       2007.01.08 ; G.Mitsuka add NECOHEPI
*       2007.08.31 ; G.Mitsuka rename ITAUFLG to ITAUFLGCORE
*                    not to conflict with neutflux
*       2007.11.05 ; G.Mitsuka add target material information
*       2010.06    ; P.de Perio - add model select for high energy pions
*                               - add kinematical model select for high 
*                                 energy pion inelastic scattering
*       2012.12    ; P.Rodrigues added form factors
*       2019.10    ; J.Xia      - add more iterations for Bodek-Yang corrections
************************************************************************
      INTEGER NEFRMFLG,NEPAUFLG,NENEFO16,NENEFMODL,NENEFMODH,
     &        NENEFKINH,NEMODFLG,NESELMOD,ITAUFLGCORE,NUSIM, QUIET
      REAL   CRSNEUT,CRSNEUTB

      COMMON/NEUTCARD/NEFRMFLG,NEPAUFLG,NENEFO16,NENEFMODL,NENEFMODH,
     &                NENEFKINH,NEMODFLG,NESELMOD,CRSNEUT(30),
     &                CRSNEUTB(30),ITAUFLGCORE,NUSIM,QUIET

      INTEGER NEFKINVER
      COMMON/NUCEFFVER/NEFKINVER

      INTEGER NEPDF, NEBODEK, NEMULT
      COMMON/NEUTDIS/NEPDF,NEBODEK,NEMULT

      INTEGER NEIFF,   NENRTYPE
      REAL    RNECA5I, RNEBGSCL
      REAL    XMANFFRES,XMVNFFRES,XMARSRES,XMVRSRES
      COMMON/NEUT1PI/XMANFFRES,XMVNFFRES,XMARSRES,XMVRSRES,
     $               NEIFF,NENRTYPE,RNECA5I,RNEBGSCL

      INTEGER NEDIFPI
      COMMON/NEUTDIF/NEDIFPI

      INTEGER NECOHEPI
      COMMON/NEUTCOH/NECOHEPI

      INTEGER NUMBNDN,NUMBNDP,NUMFREP,NUMATOM
      COMMON/NEUTTARGET/NUMBNDN,NUMBNDP,NUMFREP,NUMATOM

      INTEGER NEABSPIEMIT
      COMMON/NEUTPIABS/NEABSPIEMIT

      INTEGER IPILESSDCY
      REAL*4  RPILESSDCY
      COMMON /NEUTPILESS/IPILESSDCY,RPILESSDCY

      INTEGER IRADCORR
      COMMON /NEUTRADCORR/IRADCORR
						  
