************************************************************************
*     -------------------
C     INCLUDE FILE NEUTCRS
*     -------------------
*
*     (Purpose)
*       To store NEUT cross-section
*
*     (Creation Date and Author)
*       2010.07.30 ; Y.Hayato
*
*     (Updates)
*       2011.05.13 ; P. de Perio
*                      - Added more kinematic variables
*
************************************************************************
C-- TOTCRSNE : Total cross-section
C-- DIFCRSNE : Differential cross-section

      REAL*4 CRSENERGY, CRSX, CRSY, CRSZ, CRSPHI, CRSQ2
      REAL*4 TOTCRSNE,DIFCRSNE(8)
      COMMON /NEUTCRSCOM/CRSENERGY,TOTCRSNE,DIFCRSNE,
     $                   CRSX,CRSY,CRSZ,CRSPHI,CRSQ2
