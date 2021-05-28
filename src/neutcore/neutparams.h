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
C    SFEBSHIFT: Extra missing energy shift for SF events (GeV)
C                 Positive means more missing energy.
C    SFEBNEGBEH: Behavior if SFEb is below 0 (i.e. chosen nucleon is unbound)
C               0: Pin to EB = 0 and continue event generation
C               1: Redraw from SF (for large negative SFEBShift, will get stuck in infinite loops)
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
C    FEFALL:  Correction factor to the total mean free path (all momenta)
C
C---------------------------------------------------------
C
      REAL*4 PFSURF,PFMAX,VNUINI,VNUFIN,FZMU2,SFEBSHIFT
      INTEGER*4 IFORMLEN, SFEBNEGBEH

      COMMON/NENUPR/PFSURF,PFMAX,VNUINI,VNUFIN,IFORMLEN,FZMU2,
     $ SFEBSHIFT, SFEBNEGBEH
C
C---------------------------------------------------------
C
      REAL*4 FEFQE,FEFQEH,FEFINEL,FEFABS,FEFCOH,FEFQEHF,FEFCOHF,FEFCX,
     $       FEFCXHF,FEFCXH,FEFCOUL,FEFALL

      COMMON/NEFFPR/FEFQE,FEFQEH,FEFINEL,FEFABS,FEFCOH,FEFQEHF,FEFCOHF,
     $              FEFCX,FEFCXHF,FEFCXH,FEFCOUL,FEFALL
