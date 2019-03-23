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
C          = 4XX ; ALL; Use SF nuclear model for kinematics
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
C-- Model for Single pion production
C  MDLSPI = 1 ; Rein Sehgal
C
C-- MA (MV) for Rein-Sehgal Single meson production
C  XMARSRES (XMVRSRES)
C
C-- MA for New Form factor Single meson production
C  XMANFFRES (XMVNFFRES)
C
C-- Coherent pion
C  MDLCOH = 0 ; Rein&Sehgal w/ lepton mass corr.
C             1 ; Kartavtsev et al.
C
C  XMACOH  (Default = 1.0)
C  RAD0NU  (Default = 1 fm)
C
C-- DIS
C  MODELDIS = 70 ; GRV94 Original
C             71 ; GRV94 Bodek
C            120 ; GRV98 Original
C            121 ; GRV98 Bodek

      INTEGER*4 MODELDIS,MODELCOH
      INTEGER*4 MDLQE,MDLSPI,MDLDIS,MDLCOH,MDLQEAF
      REAL*4    XMAQE,XMASPI,XMARES,XMVQE,XMVSPI,XMVRES,
     $          KAPP,XMACOH,RAD0NU
      INTEGER*4 IFFSPI,NRTYPESPI
      REAL*4    RCA5ISPI,RBGSCLSPI
      REAL*4    SCCFV, SCCFA, FPQE
      REAL*4    PFSF

      
	  COMMON /NEUTMODEL/MODELDIS,MODELCOH
      COMMON /NEMDLS/MDLQE,MDLSPI,MDLDIS,MDLCOH,
     $               MDLQEAF,XMAQE,XMASPI,XMVQE,XMVSPI,
     $               KAPP,XMACOH,RAD0NU,
     $               IFFSPI,NRTYPESPI,RCA5ISPI,RBGSCLSPI,
     $               XMARES,XMVRES,
     $               SCCFV, SCCFA, FPQE,
     $               PFSF
