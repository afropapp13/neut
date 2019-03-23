      SUBROUTINE QGDISG(X0,Q2,G,U,D,AU,AD,S,C,B,T)
*
*     Quark/Gluon distribution function ["GRV94"]
*     Ref. M.Glueck, E.Reya, A.Vogt, DO-TH 94/24=DESY 94-206
*                                    (To Appear in Z.Phys.C)
*     Based on Mori_san's Subroutine "qgdim.f"
*
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL verbose/.false./
#include "necard.h"
*     *** QMIN CHANGED 0.4 -> 0.24(Saji Aug. 02) 
*     *** QMIN CHANGED 0.4 -> 0.  (Saji Dec. 03) 
*           QMIN is defined at grv94di.F
C#ifndef NO_GRV94CORR
C      PARAMETER (MODE=0,QMIN2=0.0D0,QMAX2=1.0D6)
C#else
C      PARAMETER (MODE=0,QMIN2=0.4D0,QMAX2=1.0D6)
C#endif
C      Parameter (Xmin=1.0D-5,Xmax=1.0D0)
      COMMON /PARAM/EN,E0,IFQG,IFEL

C     MODE,QMIN2,QMAX2 are separately defined for each PDFs
C     (Aug'06 G.Mitsuka)
      ! ---> GRV94DI
      if(NEPDF.eq.7) then
         if(NEBODEK.eq.1) then ! modified
            MODE = 0.
            QMIN2=0.0D0
            QMAX2=1.0D6
         else if(NEBODEK.eq.0) then ! original
            MODE = 0.
            QMIN2=0.23D0
            QMAX2=1.0D6
         endif
      ! ---> GRV98LO
      else if(NEPDF.eq.12) then
         if(NEBODEK.eq.1) then ! modified
            MODE = 0.
            QMIN2=0.0D0
            QMAX2=1.0D6
         else if(NEBODEK.eq.0) then ! original
            MODE = 0.
            QMIN2=0.8D0
            QMAX2=1.0D6
         endif
      else
         print *, 'You had better use GRV94 or GRV98 ...'
         MODE = 0.
         QMIN2=0.0D0
         QMAX2=1.0D6
      endif

C --- Q Setting
      If(Q2.Gt.Qmax2) Q2=Qmax2
      If(Q2.Lt.Qmin2) Q2=Qmin2
      Q=Sqrt(Q2)

C --- X Setting
      If(X.Lt.Xmin) X=Xmin

C --- Mode Setting
C     IFqg = 31 : Leading Order
C            32 : Next Leading Order (MS_bar)
C            33 : Next Leading Order (Deep Inelastic Scattering)
C      ifqg=33
      ifqg=26+NEPDF

C --- Leading Order
      If(IFqg.Eq.31) Then
        CALL GRV94LO(X0,Q,UV,DV,AU,AD,AS,AC,AB,AT,G)
C --- Next Leading Order (MS_bar)
      Else If(IFqg.Eq.32) Then
        CALL GRV94HO(X0,Q,UV,DV,AU,AD,AS,AC,AB,AT,G)
C --- Next Leading Order (Deep Inelastic Scatterin)
      Else If(IFqg.Eq.33) Then
        CALL GRV94DI(X0,Q,UV,DV,AU,AD,AS,AC,AB,AT,G)
      Else If(IFqg.Eq.38) Then
C --- Leading Order
        CALL GRV98_LO(X0,Q,UV,DV,AU,AD,AS,AC,AB,AT,G)
      Else
        Write(*,*) '%QGDISG-E: unknown IFQG!',IFQG
      End If

      U = UV+AU
      D = DV+AD
      S = AS
      C = AC
      B = AB
      T = AT
      
      Return
      End
