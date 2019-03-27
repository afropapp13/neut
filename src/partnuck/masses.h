C***************************************************
C     @(#)masses.cdk    1.2 modified on 12/30/92
C     The particle masses
C     PROTON    - AMP
C     NEUTRON   - AMN
C     LAMBDA    - AML
C     SIGMA     - AMS
C     ELECTRON  - AMEL
C     MUON      - AMMU
C     PI 0      - AMPO
C     PI +/-    - AMPC
C     ETA       = AMETA
C     OMEGA     = AMOM
C     RHO       = AMRHO
C     K 0       = AMKO
C     K +-      = AMKC
C     K* 0      = AMKSO
C     K* +-     = AMKSC
C*** The leptons
      real AMEL, AMMU, AMTAU
C*** The baryons
      real AMP, AMN, AML, AMS
C*** The mesons
      real AMPO, AMPC, AMETA, AMOMG, AMRHO
C*** The strange mesons.
      real AMKO, AMKC, AMKSO, AMKSC
      COMMON/MASSES/AMEL,AMMU,
     $     AMP, AMN, AML, AMS,
     $     AMPO, AMPC, AMETA, AMRHO, AMOMG,
     $     AMKO, AMKC, AMKSO, AMKSC, AMTAU


