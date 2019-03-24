************************************************************************
*     ---------------------------------
      subroutine nefillmodel
*     ---------------------------------
C
C     (Purpose)
C        Fill flags of reaction models to common block
C
C     (Input)
C        COMMON in necard
C
C     (Output)    
C        COMMON in neutmodel
C
C     (Creation Date and Author)
C        2007.08.22 ; G.Mitsuka
*        2013.10    ; P.Rodrigues & A.Bercellieadd 
*                              new form facor for Rein&Sehgal 1pi
C
************************************************************************      
      implicit none

#include "necard.h"
#include "neutmodel.h"
#include "rscons.h"

C-- DIS
C  MODELDIS = MDLDIS  
C  MDLDIS   = 70 ; GRV94 Original
C             71 ; GRV94 Bodek
C            120 ; GRV98 Original
C            121 ; GRV98 Bodek
      if (NEPDF.eq.7) then
         if (NEBODEK.eq.0) MDLDIS = 70
         if (NEBODEK.eq.1) MDLDIS = 71
      else if (NEPDF.eq.12) then
         if (NEBODEK.eq.0) MDLDIS = 120
         if (NEBODEK.eq.1) MDLDIS = 121
      end if
      MODELDIS = MDLDIS

C---------------------------------------------------------------------
C
C-- Coherent pion
C  MODELCOH = 0 ; Rein&Sehgal w/ lepton mass corr.
C  MDLCOH =   0 ; Rein&Sehgal w/ lepton mass corr.
C             1 ; Kartavtsev et al.
      MDLCOH   = NECOHEPI
      MODELCOH = MDLCOH

C
C-- Diffractive pion
C  MODELDIF = 0 ; Rein&Sehgal w/ lepton mass corr.
C  MDLDIF =   0 ; Rein&Sehgal w/ lepton mass corr.
C
      MDLDIF   = NEDIFPI
      MODELDIF = MDLDIF

C---------------------------------------------------------------------
C-- Model for QE
C
C   MDLQE now filled in necard.F
C
C-- MDLQE = 1 ; Smith-Moniz cross-section with classical kinematics
C  
C          = 0X ; NC : Original code
C          = 1X ; NC : Dipole
C          = 2X ; NC : BBBA05 calculation
C
C
C-- Vector Form factor for QE
C
C-- MDLQEAF = 1 ; Simple dipole form factor
C
C     MDLQEAF now filled in necard.F
C
C      MDLQEAF = 1;
C
C-- MA for QE      
C     XMAQE now filled in necard.F 
C
C     
C---------------------------------------------------------------------
C-- Model for Single pion production
C
C-- MDLSPI = 1 ; Rein Sehgal
C
      MDLSPI = 1;

C-- MA for Single pion production
C     XMASPI now filled in necard.F 


C-----------------------------------------------------------------------
C Copy to rscons.h

C-- for resonance productions 
      XMARES  = XMARSRES
      XMVRES  = XMVRSRES
      
C- Copy to neutmodel.h
      IFFSPI    = NEIFF
      NRTYPESPI = NENRTYPE
      RCA5ISPI  = RNECA5I
      RBGSCLSPI = RNEBGSCL

      if (NEIFF.ne.0) then
C-- new form factor
         XMASPI  = XMANFFRES
         XMVSPI  = XMVNFFRES
      else
C-- Rein-Sehgal
         XMASPI  = XMARSRES
         XMVSPI  = XMVRSRES
      endif

      return
      end
