*     -----------------------------------------------------
      SUBROUTINE EFPFCOST(PF,PPI,SK,GK,EVP,GM,NMODE,MODE,AMP,radius)
*     -----------------------------------------------------
*
*     ( purpose )
C       ++ CALCULATE AMPLITUDE OF GIVEN PION AND NUCLEON MOMENTUM
C       INTEGRATE BY SCATTERING DIRECTION. AT THIS TIME PAULI PRINCIPLE
C       IS CONSIDERED.
*
*     ( input )
C         PF    : NUCLEON MOMENTUM IN LAB. SYSTEM(MeV)
C         PPI   : PION MOMENTUM IN LAB. SYSTEM(MeV)
C         SK    : PI-N RELATIVE MOMENTUM IN C-M SYSTEM(MeV)
C         GK    : PI-N TOTAL MOMENTUM IN LAB. SYSTEM(MeV)
C         EVP   : BOOST DIRECTION OF LAB. ---> C-M SYSTEM(MeV)
C         GM    : BOOST GAMMA
C         NMODE : NUMBER OF CALCULATING MODES
C         MODE  : INTERACTION MODE (SEE FUNCTION FSQUAR)
C         radius: radius(fm)
*     ( output )
C         AMP   : AMPLITUDE OF INELASTIC SCATTERING OF EACH MODES
C
      DIMENSION PF(3),PPI(3),SK(3),SKD(3),GK(3),EVP(3),EVM(3),PPF(3)
      DIMENSION MODE(5),AMP(5)
      DATA AMN/939.0/
      DATA MESH/10/
C
      DO 10 I=1,NMODE
   10 AMP(I)=0.
C
      FSURF=EFFRMGAS(16.,8.,radius)
      EVM(1)=-EVP(1)
      EVM(2)=-EVP(2)
      EVM(3)=-EVP(3)
C
C     -- NOW INTEGRATE BY SCATTERING DIRECTION
C
      DCOST=2./FLOAT(MESH)
      DPHI=6.283185/FLOAT(MESH)
      DO 100 I=1,MESH
         COST=-1.+DCOST*(FLOAT(I)-0.5)
         DO 200 J=1,MESH
            PHI=DPHI*(FLOAT(J)-0.5)
C     - SET SKD
            CALL NECHADIR(COST,PHI,SK,SKD)
C     - BOOST BACK SKD
            PPF(1)=SKD(1)
            PPF(2)=SKD(2)
            PPF(3)=SKD(3)
            CALL MCVECBST(PPF,AMN,EVM,GM)
            PPFABS=SQRT(PPF(1)**2+PPF(2)**2+PPF(3)**2)
            IF(PPFABS.LT.FSURF)GO TO 200
            DO 210 K=1,NMODE
               MODEE=MODE(K)
               AMP(K)=AMP(K)+EFFSQUAR(MODEE,SK,SKD,GK)
  210       CONTINUE
  200    CONTINUE
  100 CONTINUE
C
      RETURN
      END
