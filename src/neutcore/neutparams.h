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
