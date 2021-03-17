**********************************************************************
*     ------------------------
      SUBROUTINE BYKVEC(QS,KVU,KVD,KSU,KSD)
*     ------------------------
*     
*     ( purpose )
*      Computes vector K factors of Bodek-Yang model
*
*     (input)
*     QS: Q^2  
*
*     ( output )
*     KVU: valence up quark K factor
*     KVD: valence down quark K factor
*     KSU: sea up quark K factor
*     KSD: sea down quark K factor      
*
*     ( creation date and author )      
*       2020.11.13 ; C. Bronner
*                    First version
*
**********************************************************************
#include "necard.h"
      
      REAL*8 QS
      REAL KVU,KVD,KSU,KSD
      REAL BY_C_val1_d, BY_C_val1_u, BY_C_val2_d
      REAL BY_C_val2_u, BY_C_sea_d, BY_C_sea_u


C      PRINT *,"Q2",QS
      
      if(NEBODEK.eq.1) then
*     Parameters from arXiv:hep-ex/0301036 (default for NEUT 5.3.4)
C     BY_C_val1_d=0.544
C     BY_C_val1_u=0.544
C     BY_C_val2_d=0.431
C     BY_C_val2_u=0.431
C     BY_C_sea_d=0.380
C     BY_C_sea_u=0.380

*     Parameters from arXiv:hep-ex/0508007 (default for NEUT 5.4.0)
         BY_C_val1_d=0.202
         BY_C_val1_u=0.291
         BY_C_val2_d=0.255
         BY_C_val2_u=0.189
         BY_C_sea_d=0.621
         BY_C_sea_u=0.363         
         
      else if (NEBODEK.eq.2) then
*     Using new parameters from BY2019 paper in prep.       
         BY_C_val1_d=0.341
         BY_C_val2_d=0.323
         BY_C_val1_u=0.417
         BY_C_val2_u=0.264
         BY_C_sea_d=0.561
         BY_C_sea_u=0.369
      endif

      
*     Compute K factors from Q2 and C parameters
      if(NEBODEK.eq.1.or.NEBODEK.eq.2) then
         KVU=(1. - 1./(1. + QS/0.71)**4)*(QS +BY_C_val2_u )/
     &        (QS + BY_C_val1_u)
         KVD=(1. - 1./(1. + QS/0.71)**4)*(QS +BY_C_val2_d )/
     &        (QS + BY_C_val1_d)
         KSU=QS/(QS + BY_C_sea_u)
         KSD=QS/(QS + BY_C_sea_d)
      else
         KVU=1.0
         KVD=1.0
         KSU=1.0
         KSD=1.0
      endif

      
      RETURN 
      END
