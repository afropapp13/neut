C
C-- Model for QE
C  MDLQE   = XX1 ; CC : Smith-Moniz cross-section with classical kinematics
C            XX2 ; CC : Smith-Moniz cross-section with classical kinematics
C                      ( BBBA05 vector form factor )
C
C          = X0X ; NC : Original code
C          = X1X ; NC : Dipole
C          = X2X ; NC : BBBA05 calculation
C
C          = 0XX ; ALL; No correction on GM form factor
C          = 1XX ; ALL; Correction on GM form factor ( Bodek et al )
C
C          *** 4XX, 6XX and 7XX ignores the 1st and 2nd digits ***
C
C          = 4XX ; ALL; Use SF nuclear model for kinematics
C          = 6XX ; ALL; Use Effective SF nuclear model for kinematics
C          = 7XX ; ALL; Use TEM-effective SF nuclear model for kinematics
C
C          =1xxx ; apply RPA correction
C
C          *** 2XXX ignores 1st, 2nd and 3rd digits
C          =2xxx ; Nieves LFG
C
C
C-- Axial form factor for QE
C  MDLQEAF = 1 ; Simple dipole form factor
C  MDLQEAF = 0;  BBBA07
C
C-- MA (MV) value of QE
C  XMAQE (XMVQE)
C
C-- Kappa factor for QE (ala MiniBooNE)
C  KAPP
C
C-- MA value of NCEL
C  XMANCEL
C
C-- Nieves CCQE flags
C
C-- Use Bidning energy or not
C
C  NVQEBIND =1 ; Binding energy != 0
C  NVQEBIND =0 ; Binding energy == 0
C
C-- RFG or not
C
C  NVQERFG = 0 ; default : ( not used )
C
C-- Turn on / off RPA
C  NVQERPA =1 ; With RPA correction
C  NVQERPA =0 ; Without RPA correction
C
C-- RPA correction parameters
C  NVRPAFP0IN     = 0.33
C  NVRPAPF0EX     = 0.45
C  NVRPAF         = 1
C  NVRPAFSTAR     = 2.13
C  NVRPAPILAMBDA  = 1200.
C  NVRPACR0       = 2.0
C  NVRPARHOLAMBDA = 2500.
C  NVRPAGP        = 0.63
C  NVRPAXMPI      = 139.57
C  NVRPAXMRHO     = 770.0
C  NVRPAIREL      = 1.
C
C-- FermiMomentum used in SF model only (for pauli blocking)
C-- (if user doesn't define, uses default for nucleus)
C  PFSF
C
C-- Second class form factors for CCQE (both 0.0 as default)
C  SCCFV
C  SCCFA
C
C-- Error term for pseudoscalar form factor for CCQE
C  QEFP
C
C-- Model for 2p2h interaction
C  MDL2P2H = 1 ; Nieves model with table
C  MDL2P2H = 2 ; Nieves model with hadron tensor
C
C-- Model for Single pion production
C  MDLSPI = 1 ; Rein Sehgal
C  MDLSPI = 2 ; Rein Sehgal + Local Fermi gas ( yet to be implemented )
C  MDLSPI = 3 ; Minoo's model
C               Controlled by rca5ispi, xmaspi, xmvspi, xmabkgm
C
C-- Model for 1 pion ejection direction
C  MDLSPIEJ = 0 ; Isotropic
C  MDLSPIEJ = 1 ; Delta only
C  MDLSPIEJ = 2 ; All resonances 
C  MDLSPIEJ = 3 ; Delta + Isotropic (same as old NEUT)

C  SPIDELTA = 0 ; The event was generated with isotropic
C           = 1 ; The event was generated with Delta res
C           = -1; variable was not set (not applicable to event)
C
C-- MA (MV) for Rein-Sehgal Single meson production
C  XMARSRES (XMVRSRES)
C  Set in necard.h and copied over into xmaspi here
C
C-- MA for New Form factor Single meson production
C  XMANFFRES (XMVNFFRES)
C  Set in necard.h and copied over into xmaspi here
C
C-- Coherent pion
C  MDLCOH = 0 ; Rein&Sehgal w/ lepton mass corr.
C           1 ; Kartavtsev et al.
C           2 ; Berger&Sehgal
C
C  XMACOH  (Default = 1.0)
C  RAD0NU  (Default = 1 fm)
C  fA1COH  (Default = 0.0)
C  fb1COH  (Default = 0.0)
C
C-- DIS
C  MODELDIS = 70 ; GRV94 Original
C             71 ; GRV94 Bodek
C            120 ; GRV98 Original
C            121 ; GRV98 Bodek
C
C-- Diffractive pion
C  MDLDIF = 0 ; Rein&Sehgal w/ lepton mass corr.
C
C  XMADIF  (Default = 1.1)
C  NUCVOLDIF  (Default = 7 GeV^-2)
C

      INTEGER*4 MODELDIS,MODELCOH,MODELDIF
      INTEGER*4 MDLQE,MDL2P2H,MDLSPI,MDLDIS,MDLCOH,MDLDIF,MDLQEAF
      INTEGER*4 MDLSPIEJ
      REAL*4    XMAQE,XMASPI,XMARES,XMVQE,XMVSPI,XMVRES,
     $          KAPP,XMACOH,RAD0NU,fA1COH,fb1COH,XMABKGM
      INTEGER*4 IFFSPI,NRTYPESPI
      REAL*4    RCA5ISPI,RBGSCLSPI
      REAL*4    SCCFV, SCCFA, FPQE
      REAL*4    PFSF,XMADIF,NUCVOLDIF
      INTEGER*4 AXZEXPQ4, AXZEXPNT
      REAL*4    AXFFALPHA, AXFFGAMMA,
     $          AXFFTHETA, AXFFBETA
      REAL*4    AXZEXPT0, AXZEXPTC,
     $          AXZEXPA0,AXZEXPA1,AXZEXPA2,
     $          AXZEXPA3,AXZEXPA4,AXZEXPA5,
     $          AXZEXPA6,AXZEXPA7,AXZEXPA8,
     $          AXZEXPA9

      REAL*4    XMANCEL
      INTEGER*4 SPIDELTA

      COMMON /NEUTMODEL/MODELDIS,MODELCOH,MODELDIF
      COMMON /NEMDLS/MDLQE,MDLSPI,MDLDIS,MDLCOH,MDLDIF,MDLSPIEJ,
     $               MDLQEAF,XMAQE,XMASPI,XMVQE,XMVSPI,
     $               KAPP,XMACOH,RAD0NU,fA1COH,fb1COH,
     $               IFFSPI,NRTYPESPI,RCA5ISPI,RBGSCLSPI,
     $               XMARES,XMVRES,XMABKGM,
     $               SCCFV, SCCFA, FPQE,
     $               PFSF,XMADIF,NUCVOLDIF,
     $               AXFFALPHA, AXFFGAMMA,
     $               AXFFTHETA, AXFFBETA,
     $               AXZEXPQ4, AXZEXPNT,
     $               AXZEXPT0, AXZEXPTC,
     $               AXZEXPA0,AXZEXPA1,AXZEXPA2,
     $               AXZEXPA3,AXZEXPA4,AXZEXPA5,
     $               AXZEXPA6,AXZEXPA7,AXZEXPA8,
     $               AXZEXPA9,
     $               MDL2P2H,
     $               XMANCEL,
     $               SPIDELTA

#include "nieves1p1h.h"
#include "nieves2p2h.h"


