************************************************************************
*     -------------------
C     INCLUDE FILE NECARDBM
*     -------------------
*
*     (Purpose)
*       To store NUBM card  ( NUBMGEN )
*
*     (Creation Date and Author)
*       2007.11.18 ; Y.Hayato
*
************************************************************************
C-- IDBMFLX  : flux flavor
C-- IDBMPID  : interacting particle * neutrino *
C-- NEVBM    : Number of events to be generated
C-- IDBMDET  : Detector ID ( SK : 0 , ND : 1 to X )
C-- BMRADMX  : Maximum Radius ( maybe limited by the flux ) in cm

      INTEGER*4 IDBMFLX,IDBMPID,NEVBM,IDBMDET
      REAL*4    BMRADMX

      COMMON /NUBMGEN/IDBMFLX, IDBMPID, NEVBM, IDBMDET, BMRADMX
