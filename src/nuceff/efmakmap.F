*     -----------------------------------------------------
      SUBROUTINE EFMAKMAP(PF,PPI,SK,GK,EVP,GM,MODE,XMAP,radius)      
*     -----------------------------------------------------
*     
*     ( purpose )
*       MAKE MAP OF SCATTERING ANGLE (2 DIMENSIONAL OF THETA AND PHI)
*
*     ( input ) 
C          PF    : NUCLEON MOMENTUM IN LAB. SYSTEM(MeV)
C          PPI   : PION MOMENTUM IN LAB. SYSTEM(MeV)
C          SK    : PI-N RELATIVE MOMENTUM IN C-M SYSTEM(MeV)
C          GK    : PI-N TOTAL MOMENTUM IN LAB. SYSTEM(MeV)
C          EVP   : BOOST DIRECTION OF LAB. ---> C-M SYSTEM(MeV)
C          GM    : BOOST GAMMA
C          MODE  : INTERACTION MODE (SEE FUNCTION FSQUAR)
*          radius: radius (fm)
*     ( output ) 
C          XMAP  : THETA PHI MAP OF SCATTERING

      DIMENSION PF(3),PPI(3),SK(3),SKD(3),GK(3),EVP(3),EVM(3),PPF(3)
      DIMENSION XMAP(11,11)
      DATA AMN/939.0/
      DATA MESH/10/
C     
      DO 10 I=1,11
         DO 10 J=1,11
   10 XMAP(I,J)=0.
C     
      FSURF=EFFRMGAS(16.,8.,radius)
      EVM(1)=-EVP(1)
      EVM(2)=-EVP(2)
      EVM(3)=-EVP(3)
C     
      DCOST=2./FLOAT(MESH)
      DPHI=6.283185/FLOAT(MESH)
      DO 100 I=1,MESH+1
         COST=-1.+DCOST*FLOAT(I-1)
         DO 200 J=1,MESH
            PHI=DPHI*FLOAT(J-1)
C     - SET SKD
            CALL NECHADIR(COST,PHI,SK,SKD)
C     - BOOST BACK SKD
            PPF(1)=SKD(1)
            PPF(2)=SKD(2)
            PPF(3)=SKD(3)
            CALL MCVECBST(PPF,AMN,EVM,GM)
            PPFABS=SQRT(PPF(1)**2+PPF(2)**2+PPF(3)**2)
            IF(PPFABS.LT.FSURF)GO TO 200
            XMAP(I,J)=EFFSQUAR(MODE,SK,SKD,GK)
  200    CONTINUE
  100 CONTINUE
C     
      DO 300 I=1,MESH+1
  300 XMAP(I,MESH+1)=XMAP(I,1)
C     
      RETURN
      END
