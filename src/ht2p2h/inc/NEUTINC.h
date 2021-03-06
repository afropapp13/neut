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
*       CRSNEUT(28)   : ( NEUT-CRS ) Multiplied factor to cross section
*                                    on each mode.   See nemodsel.F
*       CRSNEUTB(28)  : ( NEUT-CRSB ) Multiplied factor to cross section
*                                    on each mode.   See nemodsel.F
*
*       NECOHEPI      : ( NEUT_COHEPI ) Select Coherent pi model
*                            0 : Rein & Sehgal
*                            1 : Kartavtsev 
*
*       NEPDF         : ( NEUT-Select Parton distribution function
*                           n = 7 : GRV94
*                           n =12 : GRV98
*
*       NEBODEK       : ( NEUT- urn off/on Bodek-Yang correction
*                       0 :   off
*                       1 :   on 
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
*
************************************************************************
      INTEGER NEFRMFLG,NEPAUFLG,NENEFO16,NENEFMODL,NENEFMODH,
     &        NENEFKINH,NEMODFLG,NESELMOD,ITAUFLGCORE,NUSIM, QUIET
      REAL   CRSNEUT,CRSNEUTB

      COMMON/NEUTCARD/NEFRMFLG,NEPAUFLG,NENEFO16,NENEFMODL,NENEFMODH,
     &                NENEFKINH,NEMODFLG,NESELMOD,CRSNEUT(28),
     &                CRSNEUTB(28),ITAUFLGCORE,NUSIM,QUIET
      INTEGER NEPDF, NEBODEK
      COMMON/NEUTDIS/NEPDF,NEBODEK

      INTEGER NEIFF,   NENRTYPE
      REAL    RNECA5I, RNEBGSCL
      REAL    XMANFFRES,XMVNFFRES,XMARSRES,XMVRSRES
      COMMON/NEUT1PI/XMANFFRES,XMVNFFRES,XMARSRES,XMVRSRES,
     $               NEIFF,NENRTYPE,RNECA5I,RNEBGSCL

      INTEGER NECOHEPI
      COMMON/NEUTCOH/NECOHEPI

      INTEGER NUMBNDN,NUMBNDP,NUMFREP,NUMATOM
      COMMON/NEUTTARGET/NUMBNDN,NUMBNDP,NUMFREP,NUMATOM

      INTEGER NEABSPIEMIT
      COMMON/NEUTPIABS/NEABSPIEMIT

      INTEGER IPILESSDCY
      REAL*4  RPILESSDCY
      COMMON /NEUTPILESS/IPILESSDCY,RPILESSDCY
************************************************************************
*     neutparams.h
C
*     define NEUT parameter which is related to the
C     fermi momentum and nuclear potential
C
C---------------------------------------------------------
C
C    PFSURF : Fermi srface momentum
C               in (GeV)
C
C    PFMAX  : Maximum value of the Fermi momentum 
C                (usually, this value is same as PFSURF.)
C               in (GeV)
C
C    VNUINI : Nuclear potential (initial state) 
C               Negative VALUE (GeV)
C
C    VNUFIN : Nuclear potential (final state)
C               Negative VALUE (GeV)
C
C    IFORMLEN: FORMATION LENGTH effect ON/OFF
C              IFORMLEN=   1  : ALL ON (default)
C              IFORMLEN=   0  : ALL OFF
C              IFORMLEN= 110  : OFF for QE/Elastic
C              IFORMLEN= 100  : ON  for mPi/DIS only
C
C    FZMU2     FORMATION ZONE FREE PARAMETER (SCAT MODEL, DEFAULT 0.08) 
C---------------------------------------------------------
C
C    FEFQE  : Correction Factor to the probability (p<500)
C                              of quasielastic scattering ( single pi )
C    FEFQEH : Correction Factor to the probability (p>500)
C                              of quasielastic scattering ( single pi )
C
C    FEFINEL: Correction Factor to the probability (p>500)
C                              of hadron prod. (inel. scatter)
C
C    FEFABS : Correction Factor to the probability (p<500)
C                              of absorption
C
C    FEFCOH : Correction Factor to the probability (p>500)
C                              of coherent(forward)-scattering 
C
C    FEFCX  : Correction Factor to the charge exchange amplitude (p<500)
C
C    FEFCXH : Correction Factor to the charge exchange prob. (p>500)
C
C    FEFQEHF: Portion of QE scattering that has isotropic resonance decay 
C                              kinematics (p>500)
C
C    FEFCOHF: Amount of forward (elastic) scatter relative 
C                              to inelastic (p<500)
C
C    FEFCXHF: Portion of inel. scattering that includes true CX (p>500)
C
C    FEFCOUL: Coulomb correction to pion cascade trajectory (0: default, off)
C
C---------------------------------------------------------
C
      REAL*4 PFSURF,PFMAX,VNUINI,VNUFIN,FZMU2
      INTEGER*4 IFORMLEN

      COMMON/NENUPR/PFSURF,PFMAX,VNUINI,VNUFIN,IFORMLEN,FZMU2
C
C---------------------------------------------------------
C
      REAL*4 FEFQE,FEFQEH,FEFINEL,FEFABS,FEFCOH,FEFQEHF,FEFCOHF,FEFCX,
     $       FEFCXHF,FEFCXH,FEFCOUL

      COMMON/NEFFPR/FEFQE,FEFQEH,FEFINEL,FEFABS,FEFCOH,FEFQEHF,FEFCOHF,
     $              FEFCX,FEFCXHF,FEFCXH,FEFCOUL
C-------------------------------------------------------
C
C     INCLUDE FILE FOR NEUTRINO INTERACTION ( nework.h )
C
C     WORK AREA OF NEUTRINO INTERACTION
C  
C     CAUTION **** UNIT IS GEV ****
C
C-------------------------------------------------------
      INTEGER MAXNE, MODENE, NUMNE, IPNE, IORGNE, IFLGNE, ICRNNE
      REAL    PNE
      PARAMETER(MAXNE=100)
C
C     COMMON /NEWORK/
C  
C     MODENE      : MODE OF INTERACTION
C     NUMNE       : # OF PARTICLE
C     IPNE(I)     : PARTICLE CODE OF I-TH PARTICLE
C     PNE(3,I)    : MOMENTUM OF I-TH PARTICLE ( GEV/C )
C     IORGNE(I)   : ID OF ORIGIN PARTICLE
C     IFLGNE(I)   : FLAG OF FINAL STATE 
C                  -1 : initial particle
C                   0 : DETERMINED LATER PROCEDURE  
C                   1 : DECAY TO OTHER PARTICLE
C                   2 : ESCAPE FROM DETECTOR
C                   3 : ABSORPTION
C                   4 : CHARGE EXCHANGE
C                   5 : STOP AND NOT CONSIDER IN M.C. 
C                   6 : E.M. SHOWER
C     ICRNNE(I)  : FLAG OF TO CHASE OR NOT
C                   0 : DO NOT CHASE
C                   1 : CHASE
C
      COMMON /NEWORK/ MODENE,NUMNE,IPNE(MAXNE),
     $    PNE(3,MAXNE),IORGNE(MAXNE),IFLGNE(MAXNE),
     $    ICRNNE(MAXNE)

C-------------------------------------------------------------
#ifndef POSINNUC_INCLUDED

      INTEGER*4 IBOUND

#define MAXVC 100 

      REAL*4 POSNUC
      COMMON /POSINNUC/IBOUND,POSNUC(3,MAXVC)

#define POSINNUC_INCLUDED


***********************************************************************
*
*     neutfilepath.h
*  
*     ( purpose )
*       keep path name for the cross-section tables etc.
*  
*     ( creation date and author )
*       2013.Dec. ; Y.Hayato
*
***********************************************************************

      CHARACTER*1024 CRSTBLPATH

      COMMON /NEUTFILEPATH/CRSTBLPATH
			  
#endif
