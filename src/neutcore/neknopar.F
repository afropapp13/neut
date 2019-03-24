**********************************************************************
*     --------------------
      FUNCTION NEKNOPAR(W,IPAR,ITARG,XMEAN,CKNO,ALPHAKNO)
*     --------------------
*     
*     ( purpose )
*        Give parameters to compute pion multiplicity
*
*     ( input ) 
*        W      : Invariant mass of hadronic system( GEV )
*        IPAR   : NEUTRINO SPECIES
*        ITARG  : TARGET NUCLEON
*
*     ( output )
*        XMEAN    : Average hadron multiplicity
*        CKNO     : Model parameter
*        ALPHAKNO : Model parameter
*
*     ( creation date and author )
*       2016.07.28 ; C.Bronner
*                    First version
*
**********************************************************************
#include "necard.h"
      real *8 AKNO, BKNO
      real *4 XMEAN, CKNO, PFRAC, ALPHAKNO

      PFRAC=0

*     Set parameters, based on model, neutrino type and target nucleon   
      IF (NEMULT.eq.0) then     ! --- old model from NEUT 5.3.4  
C     Seems to be from M. Derrick et al., PRD 17.1 (1978)
C CB         AKNO=0
         AKNO=0.06
         BKNO=1.22
         CKNO=0
         ALPHAKNO=0
      else if (nemult.eq.1) then ! --- deuterium fit, hep-ph:1607.06558
         if (ipar.gt.0) then    ! --- neutrino
            if (itarg.eq.2212) then ! --- interaction on proton
               AKNO=0.58
               BKNO=1.35
               CKNO=9.92
               ALPHAKNO=-0.58
               PFRAC=0.66
            else if (itarg.eq.2112) then ! --- interaction on neutron
               AKNO=0.35
               BKNO=1.24
               CKNO=8.02
               ALPHAKNO=-0.89
               PFRAC=0.5
            endif            
         else                   ! --- anti-neutrino
            if (itarg.eq.2212) then ! --- interaction on proton
               AKNO=0.41
               BKNO=1.18
               CKNO=4.10
               ALPHAKNO=0.68
               PFRAC=0.5
            else if (itarg.eq.2112) then ! --- interaction on neutron
               AKNO=0.80
               BKNO=0.94
               CKNO=3.02
               ALPHAKNO=0.38
               PFRAC=0.66
            endif
         endif
      else if (nemult.eq.2) then ! --- AGKY model, hep-ph:0904.4043
         if (ipar.gt.0) then    ! --- neutrino
            if (itarg.eq.2212) then ! --- interaction on proton
               AKNO=0.40
               BKNO=1.42
               CKNO=7.93
               ALPHAKNO=0
            else if (itarg.eq.2112) then ! --- interaction on neutron
               AKNO=-0.20
               BKNO=1.42
               CKNO=5.22
               ALPHAKNO=0
            endif            
         else                   ! --- anti-neutrino
            if (itarg.eq.2212) then ! --- interaction on proton
               AKNO=0.02
               BKNO=1.28
               CKNO=5.22
               ALPHAKNO=0
            else if (itarg.eq.2112) then ! --- interaction on neutron
               AKNO=0.80
               BKNO=0.95
               CKNO=7.93
               ALPHAKNO=0
            endif
         endif
      Else
        Write(*,*) '%NEMULTIP: unknown multiplicity model',NEMULT
      End If

*     Compute average hadron multiplicity
      if (nemult.eq.1) then
         XMEAN=1.5*(AKNO+BKNO*ALOG(W**2)-PFRAC)+1.
      else
         XMEAN=1.5*(AKNO+BKNO*ALOG(W**2))
      endif

      RETURN
      END
