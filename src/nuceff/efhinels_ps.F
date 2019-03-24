*     -----------------------------------------------------
      SUBROUTINE EFHINELS_PS(IP,PIN,IPF,PFI,NQN,radius,ireaction)
*     -----------------------------------------------------
*     
*     ( purpose )
*       CALCULATE INELASTIC SCATTERING at high momentum
*
*     ( input ) 
C          IP   : INITIAL KIND OF PION
C          PIN  : INITIAL PION MOMENTUM(MeV)
C          NQN  : TARGET NUCLEON CHARGE
*
*     ( output ) 
C          IPF  : FINAL KIND OF PION
C          PFI  : FINAL PION MOMENTUM
*
c     ( creation & updates)
C          2010.05: P. de Perio
C
C
      implicit none

#include "vcwork.h"
#include "necard.h"
#include "neutparams.h"
C#include "efpiinth.h"

      REAL PIN(3),PFI(3),PNUC(3),PPI(3)
      REAL PNO(3),DIR(3), PBUF2(3)
      REAL PF,radius,PNOABS,PFIABS,PPIABS,COST,THETA,PHI,PINABS
      REAL ABSCM,BETA,GAMMA,TOTENUCL,TOTEPI
      INTEGER*4 IDUM,INUCF,NQPAI,NQN,IP,IPF,ireaction,INUC
      REAL tpi,thet_dsg,dsg,maxDcrs,PIPHI,thet_dsg_cm

      REAL BETAGAM(3),INVNUCDIR(3),PPIRST(3),BSTVEC(3)
      REAL ABSBETAGAM,GAMMANUCL,PPIRSTABS
      REAL GAMFCT
      EXTERNAL GAMFCT


      REAL s, NUMERATOR
      REAL PPICM(3)

      INTEGER IERR
      COMMON/EFLAG/IERR

      REAL AMN, AMPI
      DATA AMN/939./,AMPI/139./

      REAL FRAC

      REAL FSURF

      INTEGER I,J

      INTEGER NECHARGE
      EXTERNAL NECHARGE

      REAL EFFRMGAS, RLU
      EXTERNAL EFFRMGAS, RLU

C      logical first /.true./

C---  Initialize
      IPF=IP
      PFI(1)=PIN(1)
      PFI(2)=PIN(2)
      PFI(3)=PIN(3)

      if (NQN.eq.1) then
         INUC = 2212
      else 
         INUC = 2112
      end if
      INUCF=INUC

      IERR = 0

C---  Charge exchange
      IF (ireaction.eq.3) THEN
         
         CALL EFHINELSCX(IP,INUC,IPF,INUCF)
         
C---  Set reaction mode (see efdsg.F) if not charge exchange     
      ELSE 
         NQPAI=NECHARGE(IP)
         IF(NQPAI.EQ.0) THEN
            ireaction=4
         ELSE IF((NQPAI+NQN).EQ.2 .OR. (NQPAI+NQN).EQ.-1) THEN
            ireaction=1
         ELSE IF((NQPAI+NQN).EQ.1 .OR. (NQPAI+NQN).EQ.0) THEN
            ireaction=2
         ELSE
            PRINT *, 'Error in total charge in efhinels_ps'
            RETURN
         ENDIF
      END IF
      



C     Check that within valid energy range (see efdsgps.F)
      PINABS=SQRT(PIN(1)**2+PIN(2)**2+PIN(3)**2)
      TPI=SQRT(AMPI**2+PINABS**2)-AMPI

C      write(*,*) "EFHINELS Before:", PIN(1), PIN(2), PIN(3), sqrt(AMPI**2+PINABS**2)+AMN

C     TPI <= 5 MeV: no scatter
      IF (tpi.le.5) THEN
         RETURN
         
C     TPI >= 2000 MeV: Forward scattering with no nucleon ejection
      ELSE IF (tpi.ge.2000) THEN
 10      CALL EFCOHSCT(PIN,PFI)
         RETURN
         
C     400 < TPI < 2000 MeV:
C     Mix between PWA and Forward Scatter (since insufficient heavy
C     resonances included in efdsgps.F).
C     Fraction is quadratic from 0-100% between 400 and 2000 MeV.
      ELSE IF (400.lt.tpi .and. tpi.lt.2000) THEN
         
         FRAC = (-3.771790E-005*tpi**2 + 0.153023*tpi - 55.1743)/100.

c        Forward scattering with no nucleon ejection
         IF (RLU(IDUM).lt.FRAC) THEN
            CALL EFCOHSCT(PIN,PFI)
            RETURN
         END IF

c     ELSE continue to use SAID PWA below

      END IF



C     -- SET PION DIRECTION TO X-AXIS
      PPI(1)=PINABS
      PPI(2)=0.
      PPI(3)=0.
      COST=PIN(3)/PINABS
      THETA=ACOS(COST)
      IF(PIN(1).EQ.0. .AND. PIN(2).EQ.0.) THEN
         PHI=0.
      ELSE
         PHI=ATAN2(PIN(2),PIN(1))
      END IF


C     -- SET NUCLEON MOMENTUM AND DIRECTION in lab frame
      PF = 0.
      DIR(1) = 0.
      DIR(2) = 0.
      DIR(3) = 0.
      
      if(NUMATOM.gt.1) then
         IF (NEFKINVER.eq.0) THEN
 200        CALL RNFERM(PF,IDUM)
         
C     Should this depend on target nucleus size?
            FSURF=EFFRMGAS(16.,8.,radius)
            IF(PF.GT.FSURF)GO TO 200
         else IF (NEFKINVER.eq.1) THEN
            CALL EFRNMOM(radius,PF)
         else
            write(*,*) 'Invalid NEFKINVER =',NEFKINVER
            stop
         endif
         CALL RNDIR(DIR)
      end if

      PNUC(1) = DIR(1)*PF
      PNUC(2) = DIR(2)*PF
      PNUC(3) = DIR(3)*PF 

C     Boost pion to nucleon rest frame
      TOTENUCL=SQRT(AMN**2+PF**2)
      DO J=1,3
         BETAGAM(J)=PNUC(J)/TOTENUCL
         INVNUCDIR(J)=-1.0*DIR(J)
      ENDDO

      ABSBETAGAM=SQRT(BETAGAM(1)**2+BETAGAM(2)**2+BETAGAM(3)**2)
      GAMMANUCL=GAMFCT(ABSBETAGAM)

      DO J=1,3
         PPIRST(J)=PPI(J)
      ENDDO

      CALL MCVECBST(PPIRST,AMPI,INVNUCDIR,GAMMANUCL)

      PPIRSTABS=SQRT(PPIRST(1)**2+PPIRST(2)**2+PPIRST(3)**2)

      TOTEPI=SQRT(AMPI**2+PPIRSTABS**2)

      tpi=TOTEPI-AMPI

C     Debug
C      if (first) then
C         CALL HBOOK1(4,'INSIDE INPUT DIFF XSEC COS(THETA)$',
C     $        180,0,180,0.)  ! Debug
C         open (unit=50,file = "pip_efdsg_xsecs_250.txt")  ! Debug
C         first = .false.
C      end if


C     Find maximum of the diff. xsec. at this energy.
      thet_dsg = 0
      maxDcrs = 0.

      DO 100 I=1,180
         
         if (thet_dsg.GT.180)  thet_dsg=180         

         CALL EFDSG(ireaction,tpi,thet_dsg,dsg,ierr)
         if ( ierr .ne. 0 ) goto 10

C         CALL HF1(4,thet_dsg,dsg)  ! Debug

C     Cross section returned is normalized to solid angle: 
C     i.e. 2*pi*sin(theta)*dtheta
         dsg = sin(3.14159*thet_dsg/180.)*dsg

         if (dsg.GE.maxDcrs) maxDcrs = dsg

         thet_dsg = I
C         thet_dsg = thet_dsg + I  ! bug in < v5.1.4

 100  CONTINUE

C      write(*,*) "MAX dsg = ", maxDcrs
    

C     If something wrong
      if (maxDcrs.LE.0.) then
         WRITE(6,*)  "Error efhinels_ps: max DSG <= 0" 
         STOP
      endif

C     Calculate pion momentum in C.M. frame
      if (NUMATOM.eq.1) then
         s = AMPI**2 + AMN**2 
     $        + 2*TOTEPI*AMN

C     Consider binding energy in nucleus
      else
         s = AMPI**2 + (AMN+VNUINI*1000)**2 
     $        + 2*TOTEPI*(AMN+VNUINI*1000)
      end if

      NUMERATOR = s**2 - 2*s*(AMN**2+AMPI**2)
     &               + AMN**4 + AMPI**4 
     &               - 2*(AMN**2)*(AMPI**2)

C     Threshold for reaction
      if (NUMERATOR.LT.0) THEN
         WRITE(6,*)  "Error efhinels_ps: Below threshold for reaction" 
         STOP
      ENDIF

      PFIABS = SQRT( NUMERATOR / (4*s) )

C     Determine direction of pion/eta
 201  thet_dsg=180*RLU(IDUM)
            
C     Get differential cross section
      CALL EFDSG(ireaction,tpi,thet_dsg,dsg,ierr)
      if ( ierr .ne. 0 ) goto 10

      thet_dsg=3.14159*thet_dsg/180.

      dsg = sin(thet_dsg)*dsg

c     Accept/reject method for angle distribution
      IF (dsg/maxDcrs.LT.RLU(IDUM)) GOTO 201

C     Relative to X-axis
      PPICM(1)=PFIABS*COS(thet_dsg)
      PPICM(2)=0
      PPICM(3)=PFIABS*SIN(thet_dsg)

C     - RANDOM ROTATION AROUND X-AXIS
      PIPHI=RLU(IDUM)*2*3.14159
      PFI(1)=PPICM(1)
      PFI(2)=PPICM(2)*COS(PIPHI)-PPICM(3)*SIN(PIPHI)
      PFI(3)=PPICM(2)*SIN(PIPHI)+PPICM(3)*COS(PIPHI)
C      PFI(2)=PPICM(2)*COS(PIPHI) ! bug in < v5.1.4
C      PFI(3)=PPICM(3)*SIN(PIPHI) ! bug in < v5.1.4

C     Nucleon recoil in C.M. frame
      PNO(1)=-PFI(1)
      PNO(2)=-PFI(2)
      PNO(3)=-PFI(3)

C     Debug
C      PFIABS = sqrt(PFI(1)**2+PFI(2)**2+PFI(3)**2)
C      thet_dsg = acos((PIN(1)*PFI(1)+PIN(2)*PFI(2)
C     $        +PIN(3)*PFI(3))/PFIABS/PINABS)/3.14159*180.
C      write(*,*) "Before Boost: ", PFIABS, thet_dsg
C      CALL HF1(4,thet_dsg,1./sin(thet_dsg/180.*3.14159))  ! Debug
C      CALL HF1(4,thet_dsg,1.) ! Debug
C      thet_dsg_cm = thet_dsg

C     Boost from C.M. frame into (x-axis) lab frame 
      ABSCM = SQRT((PPI(1)+PNUC(1))**2+(PPI(2)+PNUC(2))**2 
     &     + (PPI(3)+PNUC(3))**2)

      PPIABS=SQRT(PPI(1)**2+PPI(2)**2+PPI(3)**2)
      TOTEPI=SQRT(AMPI**2+PPIABS**2)
      BETA = ABSCM/( TOTENUCL + TOTEPI )
      GAMMA = GAMFCT(BETA)
      BSTVEC(1) = (PPI(1)+PNUC(1))/ABSCM
      BSTVEC(2) = (PPI(2)+PNUC(2))/ABSCM
      BSTVEC(3) = (PPI(3)+PNUC(3))/ABSCM  

      CALL MCVECBST(PNO,AMN,BSTVEC,GAMMA)
      CALL MCVECBST(PFI,AMPI,BSTVEC,GAMMA)

C     Debug
C      PFIABS = sqrt(PFI(1)**2+PFI(2)**2+PFI(3)**2)
C      thet_dsg = acos((PIN(1)*PFI(1)+PIN(2)*PFI(2)
C     $        +PIN(3)*PFI(3))/PFIABS/PINABS)/3.14159*180.
C      write(*,*) " After Boost: ", PFIABS, thet_dsg
C      CALL HF1(1,thet_dsg,1./sin(thet_dsg/180.*3.14159))  ! Debug
C      CALL HF1(1,thet_dsg,1.) ! Debug
C      write(50,*) thet_dsg_cm, thet_dsg

C     Pauli-blocking
      if (NEPAUFLG.EQ.0 .and. NUMATOM.gt.1) then 
         PNOABS = SQRT(PNO(1)**2+PNO(2)**2+PNO(3)**2)
         IF(PNOABS.LT.FSURF)GO TO 200
      END IF


C     -- BACK TO INITIAL COORDINATE

C     - THETA ROTATION PION
      PBUF2(1)=SIN(THETA)*PFI(1)-COS(THETA)*PFI(3)
      PBUF2(2)=PFI(2)
      PBUF2(3)=COS(THETA)*PFI(1)+SIN(THETA)*PFI(3)
C     - PHI ROTATION
      PFI(1)=COS(PHI)*PBUF2(1)-SIN(PHI)*PBUF2(2)
      PFI(2)=SIN(PHI)*PBUF2(1)+COS(PHI)*PBUF2(2)
      PFI(3)=PBUF2(3)
      

C     - THETA ROTATION NUCLEON
      PBUF2(1)=SIN(THETA)*PNO(1)-COS(THETA)*PNO(3)
      PBUF2(2)=PNO(2)
      PBUF2(3)=COS(THETA)*PNO(1)+SIN(THETA)*PNO(3)
C     - PHI ROTATION
      PNO(1)=COS(PHI)*PBUF2(1)-SIN(PHI)*PBUF2(2)
      PNO(2)=SIN(PHI)*PBUF2(1)+COS(PHI)*PBUF2(2)
      PNO(3)=PBUF2(3)


C     Debug
C      PNOABS = SQRT(PNO(1)**2+PNO(2)**2+PNO(3)**2)
C      PPIABS=SQRT(PFI(1)**2+PFI(2)**2+PFI(3)**2)
C      TOTEPI=SQRT(AMPI**2+PPIABS**2)
C      TOTENUCL=SQRT(AMN**2+PNOABS**2)

C      write(*,*) "EFHINELS after: ", PFI(1)+PNO(1), PFI(2)+PNO(2), PFI(3)+PNO(3), TOTEPI+TOTENUCL
C      write(*,*) "EFHINELS after: ", PFI(1), PFI(2), PFI(3), TOTEPI+TOTENUCL


C     ADD NUCLEON TO VCWORK (non-functional, need better way
C                            to keep track of indexing)
C      if (neabspiemit.eq.0) then
C         return
C      endif
C
C     Do not store if stack is full (save one spot for pion)
C      if (nvc.GE.maxvc-1) then
C         return
C      end if
C      
C      ipvc(nvc+1)  = INUCF
C      icrnvc(nvc+1) = 1
C      iflgvc(nvc+1) = 0
C      ivtivc(nvc+1) = 1
C
C      iorgvc(nvc+1) = LOOPPI
C
C      pvc(1,nvc+1) = PNO(1)
C      pvc(2,nvc+1) = PNO(2)
C      pvc(3,nvc+1) = PNO(3)
C
C      nvc = nvc + 1

      RETURN
      END
