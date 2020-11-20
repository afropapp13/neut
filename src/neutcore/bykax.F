**********************************************************************
*     ------------------------
      SUBROUTINE BYKAX(QS,KAV,KAS)
*     ------------------------
*     
*     ( purpose )
*      Computes vector K factors of Bodek-Yang model
*
*     (input)
*     QS: Q^2  
*
*     ( output )
*     KAV: valence quark K factor
*     KAS: sea quark K factor
*
*     ( creation date and author )      
*       2020.11.20 ; C. Bronner
*                    First version from J. Xia's work
*
**********************************************************************

      REAL*8 QS
      REAL KAV,KAS
      REAL C_axi_sea,P_axi_sea,P_axi_val


*     Parameters from BY2019 paper in prep.      
      C_axi_sea=0.75
      P_axi_sea=0.55
      P_axi_val=0.018


      KAS = (QS + P_axi_sea * C_axi_sea)/(QS + C_axi_sea)
      KAV = (QS + p_axi_val)/(QS + 0.18)

      RETURN 
      END
