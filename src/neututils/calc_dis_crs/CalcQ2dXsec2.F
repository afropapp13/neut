      PROGRAM CALCQ2DXSEC2
C     
C     calculate DIS xsec as a function of x and Q2
C     
      IMPLICIT NONE
#include "../../neutcore/necard.h"

C      EXTERNAL STRUCTG
C      EXTERNAL QGDISG
C      EXTERNAL GRV98_LO
      REAL*8 DNDXDY
      EXTERNAL DNDXDY
      Real tpispir
      External tpispir
      Integer iargc
      INTEGER*4 IXX,IPLEP
      REAL*4 AMLEP,YMAX,RTT,QSMX,QSMN,XMAX
      REAL*4 PMASS,NMASS,NUCMASS
      REAL*4 ELEP,PLEP,Q2,COST
      REAL*4 W2,WMPI,WDIS
      REAL*4 XMPI,XDIS
      REAL*4 LMPI, LDIS
      REAL*4 dydis, dympi
      REAL*8 DCROSS
      INTEGER*4 CALCMPI,CALCDIS, NOXSEC

      PARAMETER (PMASS=0.938272)
      PARAMETER (NMASS=0.939565)

      INTEGER*4 IXXMAX
      REAL*4    XXMAX

C      REAL*8 DSIGMA,DSIGMA1320
      REAL*8 XX,YY,EE
      INTEGER*4 ITYPE,IPAR,IP,IT,IQ2
      INTEGER*4 IPTBL(6)
      DATA IPTBL/12,-12,14,-14,16,-16/

      Character*256 argc,output

      INTEGER*4 scale
      
      CALL NECARD
      CALL NECARDEV

*     Threshold for the 2 models
      WMPI=1.3
      WDIS=2.

      scale = 100
*     Number of points for numerical integration        
      IXXMAX=1000

      XXMAX=1000.
      
      IF(IARGC().NE.5) THEN
         WRITE(*,*) 'ARGUMENT FAILED!!'
         STOP
      ENDIF
      CALL GETARG(2,argc) ! target
      Read(argc,*) IT
      CALL GETARG(3,argc) ! neutrino type
      Read(argc,*) IP
      CALL GETARG(4,argc) ! Enu
      Read(argc,*) EE
      CALL GETARG(5,OUTPUT)

*     Compute only CC cross-section for now
      ITYPE = 1

      Open(99,file=output,form='formatted',status='unknown')
      IF (IT.EQ.2212) NUCMASS = PMASS
      IF (IT.EQ.2112) NUCMASS = NMASS
      print *,IT,ITYPE,IP, EE

*     Max Y given energy and lepton type
      IPLEP=IABS(IP)-ITYPE
      CALL MCMASS(IPLEP,AMLEP)
      AMLEP=AMLEP*1.0E-3
      YMAX=1.D0-AMLEP/EE


      NOXSEC=0
      IF(YMAX.LT.0.) THEN
         write(*,*) 'Ymax is negative'
         NOXSEC=1
      END IF

      
      
*     Q2 Loop
      DO 30 IQ2=1,1090
C     DO 30 IQ2=1,120
         IF (IQ2.le.100) THEN
            Q2=0.005D0+(0.01e0*DBLE(IQ2-1))
         ELSE
            Q2=1.05D0+(0.1e0*DBLE(IQ2-101))
         END IF

C         PRINT *, NOXSEC
*     Init xsec         
C         DSIGMA=0.D0
C         DSIGMA1320=0.D0    

         IF (NOXSEC.ne.1) THEN  !-- only calculate XSec if there's a possible value of Y
            
*     Max value of x for each mode
            XMPI=Q2/(WMPI**2-NUCMASS**2+Q2)
            XDIS=Q2/(WDIS**2-NUCMASS**2+Q2)

C     YY=Q2/(2.*NUCMASS*EE*XMPI)
C     W2=NUCMASS**2.+2.*NUCMASS*EE*YY*(1.0-XMPI)  
C     print *,'MPI',Q2,XMPI,W2

C     YY=Q2/(2.*NUCMASS*EE*XDIS)
C     W2=NUCMASS**2.+2.*NUCMASS*EE*YY*(1.0-XDIS)  
C     print *,'DIS',Q2,XDIS,W2

*     Check which xsec will be non zero
            CALCMPI=1
            CALCDIS=0
            if (XMPI.lt.0) THEN
               CALCMPI=0
C               CALCDIS=0
            else if (XMPI.gt.1) then
               XMPI = 1
            END IF
            
            if (XDIS.lt.0) THEN
               XDIS=0
               CALCDIS=0
            else if (XDIS.gt.1) THEN
               XDIS=1
               CALCMPI=0
            END IF
            
C     print *,CALCMPI, CALCDIS
            
*     Compute DIS xsec
            if (CALCDIS.eq.1) THEN
               LDIS=XDIS
               DO 40 IXX=1,IXXMAX
                  
                 XX=DBLE(IXX*(LDIS/XXMAX))
                 YY=Q2/(2.*NUCMASS*EE*XX)
C                 dydis=Q2/(2.*NUCMASS*EE)*(1/(XX-DBLE(0.5*(LDIS/XXMAX)))
C     $                 -1/(XX+DBLE(0.5*(LDIS/XXMAX))))
                  IF(YY.GT.YMAX) THEN
                     GOTO 40
                  END IF

*     Check lepton momentum is positive               
                  RTT=EE*EE*(1.E0-YY)**2-AMLEP**2
                  IF(RTT.LT.0.) THEN
C     write(*,*) 'RTT is negative'
                     GOTO 40
                  END IF

                  ELEP=EE-YY*EE
                  PLEP=ELEP**2-AMLEP**2
                  
                  PLEP=SQRT(PLEP)
                  COST=(2.*EE*ELEP-AMLEP*AMLEP-Q2)/2./EE/PLEP
                  
                  IF(ABS(COST).GT.1.) THEN
C     write(*,*) 'lepton angle is irregal'
                     GOTO 40
                  END IF

*     Check W is in the DIS mode range
                  W2=NUCMASS**2.+2.*NUCMASS*EE*YY*(1.0-XX)
                  IF (W2.LE.(WDIS**2-0.001)) THEN !- should not happen 
                     GOTO 40
                  ENDIF

                  
*     double differential cross section               
C                  print *, 'In X grid of integral: Bodek # is', nebodek 
                  DCROSS = DNDXDY(IP,IT,ITYPE,EE,XX,YY)
C     multiplied by the grid x and y width
C                  DCROSS=DNDXDY(IP,IT,ITYPE,EE,XX,YY)*dydis

C                  DSIGMA=DSIGMA+DCROSS
                  write(99,'(I4,F13.3,F13.3,F15.8,F15.8,I4)') 
     $                 IP,EE,Q2,XX,DCROSS,NEBODEK     
 40            CONTINUE
            ENDIF               !-- End of DIS xsec

*     Compute multi-pi xsec
            if (CALCMPI.eq.1) THEN
               LMPI=XMPI-XDIS
               DO 50 IXX=1,IXXMAX
                  
                 XX=XDIS+DBLE(IXX*(LMPI/XXMAX))
                 YY=Q2/(2.*NUCMASS*EE*XX)
C                 dympi=Q2/(2.*NUCMASS*EE)*(1/(XX-DBLE(0.5*(LMPI/XXMAX)))
C     $                  -1/(XX+DBLE(0.5*(LMPI/XXMAX))))

CIF (IQ2.eq.5) THEN                    
C     print *,LMPI,LMPI/XXMAX,XX,YY,YMAX
C     ENDIF
                  
                  IF(YY.GT.YMAX) THEN
                     GOTO 50
                  END IF

*     Check lepton momentum is positive               
                  RTT=EE*EE*(1.E0-YY)**2-AMLEP**2
                  IF(RTT.LT.0.) THEN
                     print *,'RTT is negative'
                     GOTO 50
                  END IF

                  PLEP=SQRT(RTT)
                  COST=(2.*EE*ELEP-AMLEP*AMLEP-Q2)/2./EE/PLEP
                  
                  IF(ABS(COST).GT.1.) THEN
C     print *,'lepton angle is irregal'
                     GOTO 50
                  END IF

*     Check W is in the multi-pi mode range
                  W2=NUCMASS**2.+2.*NUCMASS*EE*YY*(1.0-XX)
                  IF ((W2.GT.(WDIS**2+0.001)).or.
     $                 (W2.LT.(WMPI**2-0.001))) THEN !- should not happen
                     print *,'Wrong W2 -- should not see this',W2,Q2,XX
                     GOTO 50
                  ENDIF
                  
*     double differential cross section               
C                  DCROSS=DNDXDY(IP,IT,ITYPE,EE,XX,YY)*dympi
C                  print *, 'In X grid of integral: Bodek # is', nebodek 
                  DCROSS=DNDXDY(IP,IT,ITYPE,EE,XX,YY)*
     $                 tpispir(sqrt(W2),ip,it)
C                  DSIGMA1320=DSIGMA1320
C     $                 +DCROSS*tpispir(sqrt(W2),ip,it)
C     print *,DCROSS,DSIGMA1320
                  write(99,'(I4,F13.3,F13.3,F15.8,F15.8,I4)') 
     $                 IP,EE,Q2,XX,DCROSS,NEBODEK     
 50            CONTINUE
            ENDIF               !-- End of MPI xsec
         ENDIF
*     divide by size of intervals used  for the numerical integration

C         DSIGMA=DSIGMA*(LDIS/XXMAX)
C         DSIGMA1320=DSIGMA1320*(LMPI/XXMAX)
 
         
C         DSIGMA=DSIGMA*(LDIS/XXMAX)*scale
C         DSIGMA1320=DSIGMA1320*(LMPI/XXMAX)*scale
         
 
 30   CONTINUE
      END
