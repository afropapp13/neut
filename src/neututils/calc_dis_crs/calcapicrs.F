      PROGRAM CALCAPICRS
C     
C     CALCULATION multi pion Production cross-section
C     
      IMPLICIT NONE
      REAL*8 DNDXDY
      EXTERNAL DNDXDY
      Real tpispir
      External tpispir
      Integer iargc
      INTEGER*4 IXX,IYY,IPLEP
      REAL*4 AMLEP,YMAX,YMIN,RTT,QSMX,QSMN,XMAX,XMIN
      REAL*4 PMASS,NMASS,NUCMASS
      REAL*4 ELEP,PLEP,Q2,COST
      REAL*4 W2
      REAL*8 DCROSS
      real*8 u,d,f2,xf3

      PARAMETER (PMASS=0.938272)
      PARAMETER (NMASS=0.939565)

      INTEGER*4 IXXMAX,IYYMAX
      REAL*4    XXMAX,YYMAX
      
C      REAL*8 DSIGMA,DSIGMA2,DSIGMA1320,DSIGMA1314,DSIGMAAL
      REAL*8 DSIGMA,DSIGMA1320
      REAL*8 XX,YY,EE,XXX
      INTEGER*4 ITYPE,IPAR,IP,IT,IETBL
      INTEGER*4 IPTBL(6)
      DATA IPTBL/12,-12,14,-14,16,-16/

      Character*256 argc,output
      Integer iparmin, iparmax
      
      INTEGER*4 JMOD,JMDTBL(6)
      DATA JMDTBL/21,-21,21,-21,21,-21/

      CALL NECARD
      CALL NECARDEV

      IXXMAX=1000
      IYYMAX=1000
      XXMAX=1000.
      YYMAX=1000.
      IPARMIN = 3
      IPARMAX = 3

      print *, 'xmax, ymax=',ixxmax, iyymax
      print *, 'iparmin, iparmax=',iparmin, iparmax

C     read(5,*) IT
C     read(5,*) ITYPE
C     read(5,*) IP
C     read(5,*) output
*---- Argument check
      IF(IARGC().NE.4) THEN
         WRITE(*,*) 'ARGUMENT FAILED!!'
         STOP
      ENDIF
      CALL GETARG(2,argc)
      Read(argc,*) IT
      CALL GETARG(3,argc)
      Read(argc,*) IP
      CALL GETARG(4,OUTPUT)
C      IT = 2212
      ITYPE = 1
C      IP = 14
      Open(99,file=output,form='formatted',status='new')
      IF (IT.EQ.2212) NUCMASS = PMASS
      IF (IT.EQ.2112) NUCMASS = NMASS
      print *,IT,ITYPE,IP,NUCMASS
      
*---- energy loop
C      DO 30 IETBL=1,20
C      DO 30 IETBL=1,161
      DO 30 IETBL=1,254
C      DO 30 IETBL=162,254
C      DO 30 IETBL=254,254

C DATA from 0.5 GeV to 6.5GeV(0.1GeV step,61bins)
         IF (IETBL.le.61) THEN
            EE =0.5+(0.1e0*DBLE(IETBL-1))
         ENDIF
C DATA from 6.5 GeV to 25.0GeV(0.5GeV step,38bins)
         IF ((IETBL.ge.62).and.(IETBL.le.99))  THEN
            EE =6.5+(0.5*DBLE(IETBL-62))
         ENDIF
C DATA from 25.0GeV to 50.0GeV(1.0GeV step,26bins)
         IF ((IETBL.ge.100).and.(IETBL.le.125))  THEN
            EE =25.0+(1.0*DBLE(IETBL-100))
         ENDIF
C DATA from 50.0GeV to 100GeV(5.0GeV step,11bins)
         IF ((IETBL.ge.126).and.(IETBL.le.136))  THEN
            EE =50.0+(5.0*DBLE(IETBL-126))
         ENDIF
C DATA from 100GeV to 2500GeV(100GeV step,26bins)
         IF ((IETBL.ge.137).and.(IETBL.le.161))  THEN
            EE =100.0+(100.0*DBLE(IETBL-137))
         ENDIF
C-----> upmu
C     DATA from 2.5TeV to 20TeV(500GeV step,36bins)
         IF ((IETBL.ge.162).and.(IETBL.le.197))  THEN
            EE =2500.0+(500.0*DBLE(IETBL-162))
         ENDIF
C     DATA from 20TeV to 50TeV(1TeV step,31bins)
         IF ((IETBL.ge.198).and.(IETBL.le.228))  THEN
            EE =20000.0+(1000.0*DBLE(IETBL-198))
         ENDIF
C     DATA from 50TeV to 100TeV(2TeV step,26bins)
         IF (IETBL.ge.229)  THEN
            EE =50000.0+(2000.0*DBLE(IETBL-229))
         ENDIF

         DSIGMA=0.D0
C         DSIGMA2=0.D0
         DSIGMA1320=0.D0
C         DSIGMA1314=0.D0
C         DSIGMAAL=0.D0
         
         DO 40 IXX=1,IXXMAX
            XX=DBLE(IXX*(1./XXMAX))
            
            DO 50 IYY=1,IYYMAX
               YY=DBLE(IYY*(1./YYMAX))
               
               IPLEP=IABS(IP)-ITYPE
               CALL MCMASS(IPLEP,AMLEP)
               AMLEP=AMLEP*1.0E-3
               YMAX=1.-AMLEP/EE
               
               IF(YMAX.LT.0.) THEN
                  GOTO 50
               END IF
               
               IF(YY.GT.YMAX) THEN
                  GOTO 50
               END IF
               
               RTT=EE*EE*(1.E0-YY)**2-AMLEP**2
               IF(RTT.LT.0.) THEN
C     write(*,*) 'RTT is negative'
                  GOTO 50
               END IF
               
               QSMX=2.*EE*EE*(1.E0-YY)+2.*EE*SQRT(RTT)-AMLEP**2
               QSMN=2.*EE*EE*(1.E0-YY)-2.*EE*SQRT(RTT)-AMLEP**2
               IF(QSMX.LT.0.) THEN
C     write(*,*) 'Qmax is negative'
                  GOTO 50
               END IF
               
C     XMAX=QSMX/(PMASS**2+2.*PMASS*YY*EE)
               XMAX=QSMX/(2.*NUCMASS*YY*EE)
               XMIN=QSMN/(2.*NUCMASS*YY*EE)
               IF(XX.GE.XMAX) THEN
C     write(*,*) 'X is larger than xmax'
                  GOTO 50
               END IF
               
               Q2=XX*(2.*NUCMASS*YY*EE)
               
               ELEP=EE-YY*EE
               
               PLEP=ELEP**2-AMLEP**2
               
               IF(PLEP.LT.0.) THEN
C     write(*,*) 'plep is negative'
                  GOTO 50
               END IF
               
               PLEP=SQRT(PLEP)
               
               COST=(2.*EE*ELEP-AMLEP*AMLEP-Q2)/2./EE/PLEP
               
               IF(ABS(COST).GT.1.) THEN
C     write(*,*) 'lepton angle is irregal'
                  GOTO 50
               END IF
               
C     HAYATO-ORG     W2=.939**2+2*.939*EE*YY-XX*(2.*.939*YY*EE+.939**2)
               W2=NUCMASS**2.+2.*NUCMASS*EE*YY*(1.0-XX)
               
               IF (W2.LE.1.3**2) THEN
                  GOTO 50
               ENDIF
               
               call structg(ip,it,EE,XX,YY,f2,xf3,u,d)
               
CCC   CONVERT x' -> x' = x*(2.*NMASS*EE*YY + NMASS**2)/2/NMASS/EE/YY
C     XXX = XX*(2.*NMASS*EE*YY + NMASS**2)/2/NMASS/EE/YY
               
C     DCROSS = DNDXDY(IP,EE,XX,YY,ITYPE,JMOD)
               DCROSS = DNDXDY(IP,IT,ITYPE,EE,XX,YY)
               
C     IF ((mod(IXX,100).eq.0).and.(mod(IYY,100).eq.0)) then
C     write(*,*) 'Dndxdy',DCROSS,XX,YY,xmax,ymax
C     ENDIF
               IF (W2.gt.(1.3**2)) THEN
C                  DSIGMAAL=DSIGMAAL
C     $                 +DCROSS
                  IF (W2.le.(2.0**2)) THEN
                     DSIGMA1320=DSIGMA1320
     $                    +DCROSS*tpispir(sqrt(W2))
                  Else 
                     DSIGMA=DSIGMA
     $                    +DCROSS
                  Endif
                  
C                  IF (W2.le.(1.4**2)) THEN
C                     DSIGMA1314=DSIGMA1314
C     $                    +DCROSS*tpispir(sqrt(W2))
C                  Else 
C                     DSIGMA2=DSIGMA2
C     $                    +DCROSS
C                  ENDIF
                  
               Endif
 50         CONTINUE
 40      CONTINUE
C     write(*,*) 'Dndxdy',DSIGMA,XX,YY,XMAX,XMIN,YMAX,YMIN
         DSIGMA=DSIGMA*(1./XXMAX)*(1./YYMAX)
C         DSIGMA2=DSIGMA2*(1./XXMAX)*(1./YYMAX)
         DSIGMA1320=DSIGMA1320*(1./XXMAX)*(1./YYMAX)
C         DSIGMA1314=DSIGMA1314*(1./XXMAX)*(1./YYMAX)
C         DSIGMAAL=DSIGMAAL*(1./XXMAX)*(1./YYMAX)
C         write(*,'(I4,F13.3,F15.5,F15.5,F15.5,F15.5,F15.5)') 
C     $        IP,EE,DSIGMA2,DSIGMA1314,DSIGMA,DSIGMA1320,DSIGMAAL
C         write(99,'(I4,F13.3,F15.5,F15.5,F15.5,F15.5,F15.5)') 
C     $        IP,EE,DSIGMA2,DSIGMA1314,DSIGMA,DSIGMA1320,DSIGMAAL
         write(*,'(I4,F13.3,F15.5,F15.5)') 
     $        IP,EE,DSIGMA,DSIGMA1320
         write(99,'(I4,F13.3,F15.5,F15.5)') 
     $        IP,EE,DSIGMA,DSIGMA1320
 30   CONTINUE
      END

      SUBROUTINE RMARIN(istau1, istau2, istau3)

      write(*,*) "***********************************************"
      write(*,*) "*****  Dummy RMARIN function was called.  *****"
      write(*,*) "*****  This subroutine should not be used *****"
      write(*,*) "*****  if you are trying to use tauola.   *****"
      write(*,*) "***********************************************"

      return
      end
      
