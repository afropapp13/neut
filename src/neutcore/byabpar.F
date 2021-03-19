**********************************************************************
*     ------------------------
      SUBROUTINE BYABPAR(PARA,PARB)
*     ------------------------
*     
*     ( purpose )
*       Get scaling parameters for Bodek-Yang model
*
*   
*     ( output )
*     PARA: A parameter of BY model
*     PARB: B parameter of BY model
*
*     ( creation date and author )      
*       2020.11.13 ; C. Bronner
*                    First version
*
**********************************************************************
#include "necard.h"
      REAL PARA, PARB

      if(NEBODEK.eq.1) then
         
*     Parameters from arXiv:hep-ex/0301036 (default for NEUT 5.3.4)
C     PARA=0.419
C     PARB=0.223

*     Parameters from arXiv:hep-ex/0508007 (default for NEUT 5.4.0)        
         PARA=0.538
         PARB=0.305
         
      else if(NEBODEK.eq.2) then
*     Updated parameters from Bodek 2019 paper in prep:   
         PARA=0.631
         PARB=0.380
      else                      !-- shouldn't be called a priori
         PARA=0
         PARB=1
      endif

      RETURN
      END
