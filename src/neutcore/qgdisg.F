************************************************************************
*     --------------------------------------
      SUBROUTINE QGDISG(X0,Q2,G,U,D,AU,AD,S,C,B,T)
*
*     --------------------------------------
*
*     (Purpose)
*     ++ Gets PDF based on settings in card file
*
*     (Input)
*       X0      : Q**2/(2Mv) BJORKEN x
*       Q2      : Transferred 4-momentum
*     
*     (Output)
*       x*proba for gluon, up_total, down_total, up_sea, down_sea,
*       strange_sea, charm_sea, bottom_sea and top_sea
*
*     (Creation Date and Author)
*       ????.??.?? ; ?? first version
*       2006.08.?? ; G. Mitsuka
*                    Support for different PDFs and Bodek-Yang corrections
*       2016.07.22 ; C. Bronner
*                    Cleaning and commenting
*
************************************************************************

      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL verbose/.false./
#include "necard.h"

      COMMON /PARAM/EN,E0,IFQG,IFEL

C     Allowed Q2 range for the different PDFs
      ! ---> GRV94DI
      if(NEPDF.eq.7) then
         if(NEBODEK.eq.1) then ! use corrections by Bodek and Yang          
            QMIN2=0.0D0
            QMAX2=1.0D6
         else if(NEBODEK.eq.0) then ! original         
            QMIN2=0.23D0
            QMAX2=1.0D6
         endif
      ! ---> GRV98LO
      else if(NEPDF.eq.12) then
         if(NEBODEK.eq.1) then ! use corrections by Bodek and Yang 
            QMIN2=0.0D0
            QMAX2=1.0D6
         else if(NEBODEK.eq.0) then ! original
            QMIN2=0.8D0
            QMAX2=1.0D6
         endif
      else
         print *, 'You had better use GRV94 or GRV98 ...'
         MODE = 0.
         QMIN2=0.0D0
         QMAX2=1.0D6
      endif

C --- Force Q2 to be within allowed range for the PDF
      If(Q2.Gt.Qmax2) Q2=Qmax2
      If(Q2.Lt.Qmin2) Q2=Qmin2
      Q=Sqrt(Q2)

C --- Call PDF calculation routine
C --- GRV94 Leading Order
      If(NEPDF.Eq.5) Then
        CALL GRV94LO(X0,Q,UV,DV,AU,AD,AS,AC,AB,AT,G)
C --- GRV94 Next Leading Order (MS_bar)
      Else If(NEPDF.Eq.6) Then
        CALL GRV94HO(X0,Q,UV,DV,AU,AD,AS,AC,AB,AT,G)
C --- GRV94 Next Leading Order (Deep Inelastic Scattering)
      Else If(NEPDF.Eq.7) Then
        CALL GRV94DI(X0,Q,UV,DV,AU,AD,AS,AC,AB,AT,G)
      Else If(NEPDF.Eq.12) Then
C --- GRV98 Leading Order
        CALL GRV98_LO(X0,Q,UV,DV,AU,AD,AS,AC,AB,AT,G)
      Else
        Write(*,*) '%QGDISG-E: unknown PDF!',NEPDF
      End If

*     Probabilities to return
*     Return total probability (valence+sea) for up and down quarks
      U = UV+AU
      D = DV+AD
      S = AS
      C = AC
      B = AB
      T = AT
      


      Return
      End
