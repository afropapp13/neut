************************************************************************
*     -------------------------------------------
      SUBROUTINE NEMECLPV(ENEU,ELEP,PLEP,IPMODE,IERR)
*     -------------------------------------------
*
*     (Purpose)
*       Set final lepton 4-momentum vector for MEC event,
*       as well as initial and final nucleon 4-momenta
*       and spectator nucleon isospin
*
*     (Input)
*       ENEU   : NEUTRINO ENERGY ( GEV )
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR
*                  16 : NEU-TAU
*                 -16 : NEU-TAU-BAR
*
*     (Output)
*       PLEP(3)   : FINAL LEPTON MOMENTUM ( GEV/C )
*       IPMODE(6) : PID OF PARTICLES IN INTERACTION
*       IERR    : ERROR CODE
*
*     (Creation Date and Author)
*       2013.02.26 ; P.Sinclair (based on NEELSLPV)
*       2013.08.23 ; P.Sinclair Added extension up to Enu=30GeV
*       2014.01.15 ; A.Redij Added extension of different target
************************************************************************
      IMPLICIT NONE
#include "necard.h"
#include "neutparams.h"
#include "neutfilepath.h"

      REAL    ENEU
      REAL TEMP0,TEMP1,TEMP2
      REAL*4 DUM, RLU
      INTEGER IPAR,IERR,I,J,IPMODE(6)
      REAL PLEP(3), ELEP
      CHARACTER*1024 FILENAME
      integer*4      LENSTR
      external       lenstr
      integer*4      lenpath
      INTEGER COUNTER
C  Lookup function return interger bin labels for a given Enu
      INTEGER ELOOKUP, PLOOKUP, TLOOKUP
      INTEGER EINT, PINT, TINT
      REAL XSECMAX
      REAL EXTPLT
C  Array to store xsecs at Enu bin edges
      REAL XSECBINS_INPUTL(2,31,31)
      REAL XSECBINS_INPUT(2,31,31)
      REAL XSECBINS(2,31,31)
C  Array to store interpolated xsec at precise Enu
      REAL XSEC(31,31)
      REAL XSECVAL
      REAL ETABLE(92),PTABLE(31),TTABLE(31)
      REAL P,T,HIT
C  Maximum and minimum possible values of Tmu and cos(thetaMu). Maximum Pmu somewhat arbitrarily defined here.
      REAL PMAX, PMIN, TMAX, TMIN
      INTEGER NEVALS, NPVALS, NTVALS
      REAL LEPMOM, LEPMASS, Q
C Some variables for the filename formatting
      INTEGER EINTLEN, NUINTLEN
      CHARACTER*80 FILEFORM
C Some variables to make bin interpolation more intelligible
      INTEGER PBININT_HIGH, PBININT_LOW
      REAL X1_LOWE, X1_HIGHE, X2_LOWE, X2_HIGHE
      REAL PEVAL, TEVAL
C Variable to input correct file for atom
      CHARACTER*4 ATOMSTRING
C Variables for high energy extension
      REAL Q3MAX_DATA
      REAL Q3MAX
      REAL PMAX_DATA,PMIN_DATA,TMAX_DATA,TMIN_DATA
      REAL Q3THROW
C Array of atomic mass number for which the cross-section tables are available.
      INTEGER*4 ATABLE(3)
      DATA ATABLE/12,16,40/
      INTEGER*4 IAM
      INTEGER*4 AMASS
C  Weight to scale the cross-section table for given value of Atomic mass
      REAL WEIGHT


C
C Table with hardcoded Enu values. These must correspond to points in the files in crsdat
C

      DATA  ETABLE / 0.10  , 0.1025, 0.1050, 0.1075, 0.11  ,
     $               0.1125, 0.1150, 0.1175, 0.12  , 0.1225,
     $               0.1250, 0.1275, 0.13  , 0.1325, 0.1350,
     $               0.1375, 0.14  , 0.1425, 0.145 , 0.1475,
     $               0.15  , 0.1525, 0.155 , 0.1575, 0.16  ,
     $               0.1625, 0.165 , 0.17  , 0.175 , 0.18  ,
     $               0.185 , 0.19  , 0.2   , 0.205 , 0.210 ,
     $               0.22  , 0.23  , 0.24  , 0.25  , 0.26  ,
     $               0.275 , 0.3   , 0.325 , 0.35  , 0.375 ,
     $               0.4   , 0.425 , 0.45  , 0.475 , 0.5   ,
     $               0.525 , 0.55  , 0.6   , 0.65  , 0.7   ,
     $               0.75  , 0.8   , 0.85  , 0.9   , 0.95  ,
     $               1     , 1.05  , 1.1   , 1.15  , 1.2   ,
     $               1.25  , 1.3   , 1.4   , 1.49  , 1.4999999,
     $            1.5000001, 1.502 , 1.5075, 1.513 , 1.515 ,
     $               1.52  , 1.525 , 1.55  , 1.575 , 1.6   ,
     $               1.65  , 1.8   , 1.9   , 2     , 2.5   ,
     $               3     , 4     , 6     , 7.5   , 10    ,
     $               20    , 30/

C NOTE: throughout this code: E = neutrino energy, P = Tmu (muon kinetic energy), T = cos(thetaMu) (neutrino-muon opening angle)

      NEVALS=92
      NPVALS=31
      NTVALS=31

C Q3MAX_DATA is the Q3 limit used to make the tables
C Q3MAX is the Q3 limit you want to actually use (must be <= Q3MAX_DATA)

      Q3MAX_DATA=1.5
      Q3MAX=1.2

      lenpath = lenstr(crstblpath)

      CALL MCMASS(IPMODE(4),LEPMASS)
      LEPMASS = LEPMASS*1.0E-3

C Interpolation code has problems if neutrino energy is exactly equal to the kinematic cut-off in the data because of a discontinuity in the bin definitions if Enu=q3max in data. So, if the neutrino energy is very near q3max_data, change neutrino energy by a tiny amount to make it OK. Note: change in neutrino energy is only by a maximum of 20eV!
      IF((ENEU.lt.1.5000001).and.(ENEU.gt.1.4999999)) THEN
        ENEU=1.5000001
      ENDIF

C AMASS is mass number of target element.  Check if the elements is C12, O16 or Ca40 for which tables are available set IAM to 0. If not, identify the elements on either side of the target element for interpolation. Set IAM =1 if atomic mass number of the element is between C12 and O16 and IAM to 2 if its greater than O16.

       AMASS=NUMBNDP+NUMBNDN

       IF((NUMBNDP.eq.6).and.(NUMBNDN.eq.6))THEN
       IAM=0
       GOTO 11
       ELSE IF((NUMBNDP.eq.8).and.(NUMBNDN.eq.8))THEN
       IAM=0
       GOTO 11
       ELSE IF((NUMBNDP.eq.20).and.(NUMBNDN.eq.20))THEN
       IAM=0
       GOTO 11
       ENDIF

       IF(AMASS.LT.ATABLE(1))THEN
       IAM=1
       GOTO 11
       ELSE IF(AMASS.GE.ATABLE(3))THEN
       IAM=3
       GOTO 11
       ENDIF

       DO 10 J=1,2
         IF((AMASS.GE.ATABLE(J))
     $        .and.
     $        (AMASS.LT.ATABLE(J+1))) THEN
            IAM = J
            GOTO 10
         ENDIF

 10    CONTINUE
 11    CONTINUE

      IF(IAM.EQ.1)WEIGHT=REAl(AMASS)/REAL(ATABLE(1))
      IF(IAM.EQ.2)WEIGHT=REAL(AMASS)/REAL(ATABLE(2))
      IF(IAM.EQ.3)WEIGHT=REAL(AMASS)/REAL(ATABLE(3))
      IF(IAM.EQ.0)WEIGHT=1.0


C Set Q value. This is an energy related to the difference between initial and final nuclear state in Nieves code. Card coded numbers come from values used by Nieves


      IF((IPMODE(1).gt.0).and.(((NUMBNDP.eq.8).and.(NUMBNDN.eq.8)).OR.IAM.eq.2)) THEN
C neutrino on o16
      Q=0.014906
      ELSE IF((IPMODE(1).gt.0).and.(((NUMBNDP.eq.6).and.(NUMBNDN.eq.6)).OR.IAM.eq.1)) THEN
C neutrino on c12
      Q=0.016827
      ELSE IF((IPMODE(1).lt.0).and.(((NUMBNDP.eq.6).and.(NUMBNDN.eq.6)).OR.IAM.eq.1 )) THEN
C antineutrino on c12
      Q=0.013880
      ELSE IF((IPMODE(1).lt.0).and.(((NUMBNDP.eq.8).and.(NUMBNDN.eq.8)).OR.IAM.eq.2 )) THEN
C antineutrino on o16
      Q=0.010931
      ELSE IF((IPMODE(1).gt.0).and.(((NUMBNDP.eq.20).and.(NUMBNDN.eq.20)).OR.IAM.eq.3 )) THEN
C neutrino on ca40
      Q=0.013809
      ELSE IF((IPMODE(1).lt.0).and.(((NUMBNDP.eq.20).and.(NUMBNDN.eq.20)).OR.IAM.eq.3)) THEN
C antineutrino on ca40
      Q=0.001822
      ELSE
        PRINT *,"ERROR IN nemeclpv.F: Unidentified neutrino or target."
        RETURN
      ENDIF


C Set limits for lepton kinematics (for run-time muon and muon in xsec data)

      TMAX_DATA = 1
      TMAX      = 1
      PMAX_DATA = ENEU-LEPMASS-Q
      PMAX      = ENEU-LEPMASS-Q

      IF(ENEU.lt.Q3MAX_DATA) THEN
        PMIN_DATA=0
        TMIN_DATA=-1
      ELSE
        PMIN_DATA=SQRT(LEPMASS**2+(ENEU-Q3MAX_DATA)**2)-LEPMASS
        TMIN_DATA=SQRT(1-(Q3MAX_DATA/ENEU)**2)
      ENDIF

      IF(ENEU.lt.Q3MAX) THEN
        PMIN=0
        TMIN=-1
      ELSE
        PMIN=SQRT(LEPMASS**2+(ENEU-Q3MAX)**2)-LEPMASS
        TMIN=SQRT(1-(Q3MAX/ENEU)**2)
      ENDIF

C
C Get cross-section data from files
C -------------------------------------------------------------
C

C First, get values at lower edge of Enu bin in which ENEU lies

      EINT=ELOOKUP(ENEU)


C    Some statements to unify formatting of xsec filenames
      IF(IAM.EQ.1 .OR. AMASS.EQ.12) ATOMSTRING="C12"
      IF(IAM.EQ.2 .OR. AMASS.EQ.16) ATOMSTRING="O16"
      IF(IAM.EQ.3 .OR. AMASS.EQ.40) ATOMSTRING="Ca40"

      IF(ABS(IPMODE(1)).LT.10) NUINTLEN=1
      IF(ABS(IPMODE(1)).GT.9) NUINTLEN=2
      IF(ABS(IPMODE(1)).GT.99) NUINTLEN=3
      IF(IPMODE(4).LT.0) NUINTLEN=NUINTLEN+1

      IF(EINT.LT.10) EINTLEN=1
      IF(EINT.GT.9) EINTLEN=2
      IF(EINT.GT.99) EINTLEN=3

C Now get xsec data from tables

      IF(IAM.EQ.3 .OR. AMASS.EQ.40)THEN
      WRITE(FILEFORM,'(A5,I1,A11,I1,A4)')
     & '(A7,I',NUINTLEN,',A1,A4,A1,I',EINTLEN,',A4)'
      ELSE
      WRITE(FILEFORM,'(A5,I1,A11,I1,A4)')
     & '(A7,I',NUINTLEN,',A1,A3,A1,I',EINTLEN,',A4)'
      ENDIF

      WRITE(FILENAME,FILEFORM)
     $      'Nieves_',IPMODE(1),'_',ATOMSTRING,'_',EINT, '.dat'

      FILENAME=CRSTBLPATH(1:lenpath)//FILENAME

C Get values from file
      OPEN(unit=21,file=FILENAME)
      DO 220 I=1,NTVALS
        DO 230 J=1,NPVALS
          READ(21,*)
     $    TEMP0,TEMP1,TEMP2,XSECBINS_INPUTL(1,I,J)

C   Scale the value for resp. target.
          XSECBINS_INPUT(1,I,J)=XSECBINS_INPUTL(1,I,J)*WEIGHT
 230    continue
 220  continue

      CLOSE(21)

C Now get values from file containing xsec at upper edge of Enu bin
C (or lower edge of bin above, if xsec is to be extrapolated outside range of Enu for which we have data)

      IF((EINT+1).LT.10) EINTLEN=1
      IF((EINT+1).GT.9) EINTLEN=2
      IF((EINT+1).GT.99) EINTLEN=3

      IF(IAM.EQ.3 .OR. AMASS.EQ.40)THEN
      WRITE(FILEFORM,'(A5,I1,A11,I1,A4)')
     & '(A7,I',NUINTLEN,',A1,A4,A1,I',EINTLEN,',A4)'
      ELSE
      WRITE(FILEFORM,'(A5,I1,A11,I1,A4)')
     & '(A7,I',NUINTLEN,',A1,A3,A1,I',EINTLEN,',A4)'
      ENDIF

      WRITE(FILENAME,FILEFORM)
     $      'Nieves_', IPMODE(1),'_',ATOMSTRING,'_',(EINT+1), '.dat'

      FILENAME=CRSTBLPATH(1:lenpath)//FILENAME

      OPEN(unit=31,file=FILENAME)
      DO 320 I=1,NTVALS
        DO 330 J=1,NPVALS
           READ(31,*)
     $      TEMP0,TEMP1,TEMP2,XSECBINS_INPUTL(2,I,J)

C   Scale the value for resp. target.
           XSECBINS_INPUT(2,I,J)=XSECBINS_INPUTL(2,I,J)*WEIGHT
 330    continue
 320  continue

      CLOSE(31)
C
C -------------------------------------------------------------
C

C Make list of Tmu, cosThetaMu values at this energy

      DO 335 I=1,NTVALS
        TTABLE(I) = TEVAL(I,ENEU,TMIN_DATA,TMAX_DATA,NTVALS)
 335  CONTINUE

      DO 336 I=1,NPVALS
        PTABLE(I) = PEVAL(I,ENEU,PMIN_DATA,PMAX_DATA,NPVALS)
 336  CONTINUE


C
C Linear interpolation of cross-section to precise value at ENEU
C ------------------------------------------------------------
C


C Linear inter/extra-polation of xsec at ENEU and find maximum value of xsec
C Formula: y = (((y1-y2)/(x1-x2)) * x ) + (((y2*x1)-(y1*x2)) / (x1-x2))
C where y = extra/inter-polated xsec, Enu = x
C y1, y2 = xsecs from tables at Enu = x1, x2.

C For this interpolation method, values are interpolated between points at the same fractional Tmu/cosTheta within allowed phase-space.
C In other words for grids at different Enu, cross-section is evaluated to different values of Tmu/cosTheta
C but the interpolation is always between, for example, points 'half way across allowed phase space in Tmu, a quarter of the way up allowd phase-space in cosTheta'
C At the same time, search for the maximum value of the cross-section (required to throw random number to select T-theta)

      XSECMAX=0
      DO 300 I=1,NTVALS
        DO 310 J=1,NPVALS
          XSEC(I,J)=EXTPLT(ENEU,ETABLE(EINT),ETABLE(EINT+1),
     $                     XSECBINS_INPUT(1,I,J),XSECBINS_INPUT(2,I,J))
          IF(XSEC(I,J).GT.XSECMAX) THEN
            XSECMAX=XSEC(I,J)
          ENDIF
 310    CONTINUE
 300  CONTINUE

C
C Pick a pair of lepton kinematics {Tmu, cos(thetaMu)} with probability proportional to the double-differential cross-section at ENEU
C ------------------------------------------------------------
C

      COUNTER=0
      HIT=1
      XSECVAL=0

      DO WHILE ((HIT.gt.XSECVAL).and.(COUNTER.lt.5000))
C Pick random values of Tmu and cos(thetaMu) within allowed phase-space
        P=(RLU(DUM)*(PMAX-PMIN))+PMIN
        T=(RLU(DUM)*(TMAX-TMIN))+TMIN
        LEPMOM = SQRT(P*(P+(2*LEPMASS)))
        Q3THROW=SQRT(LEPMOM**2+ENEU**2-(2.0*LEPMOM*ENEU*T))

        IF(Q3THROW.le.Q3MAX) THEN

C Inter/extra-polate xsec at these values of Tmu and cos(thetaMu)
          PINT = PLOOKUP(P,ENEU,NPVALS,PMIN_DATA,PMAX_DATA)
          TINT = TLOOKUP(T,ENEU,NTVALS,TMIN_DATA,TMAX_DATA)
C First, linear interpolation between Tmu points for the two different cos(thetaMu) values
          TEMP1 = EXTPLT(P,PTABLE(PINT),PTABLE(PINT+1),
     $                   XSEC(TINT,PINT),XSEC(TINT,PINT+1))
          TEMP2 = EXTPLT(P,PTABLE(PINT),PTABLE(PINT+1),
     $                   XSEC(TINT+1,PINT),XSEC(TINT+1,PINT+1))
C Now, linear interpolation between the two cos(thetaMu) points, for interpolated xsec at ENEU, P (muon kinetic energy), T (cos(thetaMu))
          XSECVAL = EXTPLT(T,TTABLE(TINT),TTABLE(TINT+1),TEMP1,TEMP2)
C Throw random to see if accepted
          HIT=RLU(DUM)*XSECMAX

        END IF

        COUNTER=COUNTER+1

      END DO

C Return an error if there have been too many tries. This should not really happen in practice.
      IF(COUNTER.ge.5000) GOTO 800

C
C ------------------------------------------------------------
C


C Calculate lepton momentum (Nieves predictions are of kinetic energy)
      LEPMOM = SQRT(P*(P+(2*LEPMASS)))

      PLEP(1)=LEPMOM*T
C At this point, it's OK to have all perpendiculat momentum in the +y direction
      PLEP(2)=LEPMOM*SIN(ACOS(T))
      PLEP(3)=0.
      ELEP = sqrt(LEPMASS**2 + PLEP(1)**2 + PLEP(2)**2 + PLEP(3)**2)

      IERR=0

      RETURN

 800  WRITE(6,801)
 801  FORMAT(' *** ERROR AT NEMECLPV (too many xsec tries) ***')
      IERR=1
      RETURN


      END

C This function returns an integer corresponding to the energy array entry directly below the input energy E.
C If E < lowest energy, return 1. If E > highest energy, return 15
C This is done to make it easier to access relevant xsec tables for xsec extrapolation later in the code
      INTEGER FUNCTION ELOOKUP(E)

      REAL E
      INTEGER I
      REAL ETABLE(92)

      DATA  ETABLE / 0.10  , 0.1025, 0.1050, 0.1075, 0.11  ,
     $               0.1125, 0.1150, 0.1175, 0.12  , 0.1225,
     $               0.1250, 0.1275, 0.13  , 0.1325, 0.1350,
     $               0.1375, 0.14  , 0.1425, 0.145 , 0.1475,
     $               0.15  , 0.1525, 0.155 , 0.1575, 0.16  ,
     $               0.1625, 0.165 , 0.17  , 0.175 , 0.18  ,
     $               0.185 , 0.19  , 0.2   , 0.205 , 0.210 ,
     $               0.22  , 0.23  , 0.24  , 0.25  , 0.26  ,
     $               0.275 , 0.3   , 0.325 , 0.35  , 0.375 ,
     $               0.4   , 0.425 , 0.45  , 0.475 , 0.5   ,
     $               0.525 , 0.55  , 0.6   , 0.65  , 0.7   ,
     $               0.75  , 0.8   , 0.85  , 0.9   , 0.95  ,
     $               1     , 1.05  , 1.1   , 1.15  , 1.2   ,
     $               1.25  , 1.3   , 1.4   , 1.49  , 1.4999999,
     $            1.5000001, 1.502 , 1.5075, 1.513 , 1.515 ,
     $               1.52  , 1.525 , 1.55  , 1.575 , 1.6   ,
     $               1.65  , 1.8   , 1.9   , 2     , 2.5   ,
     $               3     , 4     , 6     , 7.5   , 10    ,
     $               20    , 30/

      IF((E.LT.ETABLE(1))) THEN
        ELOOKUP=1
        RETURN
      ENDIF

      DO 20 I = 1,91
          IF((E.GE.ETABLE(I)).AND.(E.LT.ETABLE(I+1))) THEN
          ELOOKUP=I
          RETURN
        ENDIF
 20   CONTINUE

      ELOOKUP=91

      RETURN
      END


C Similar function to ELOOKUP, but for lepton kinetic energy.
      INTEGER FUNCTION PLOOKUP(P,E,NDATA,PMIN_DATA,PMAX_DATA)

      REAL    P, E
      INTEGER I, NDATA
      REAL    PMIN_DATA,PMAX_DATA

      IF(P.LT.PEVAL(1,E,PMIN_DATA,PMAX_DATA,NDATA)) THEN
        PLOOKUP=1
        RETURN
      ENDIF

      DO 21 I = 1,(NDATA-1)
          IF((P.GE.PEVAL(I,E,PMIN_DATA,PMAX_DATA,NDATA)).AND.
     $       (P.LT.PEVAL(I+1,E,PMIN_DATA,PMAX_DATA,NDATA))) THEN
          PLOOKUP=I
          RETURN
        ENDIF
 21   CONTINUE

      PLOOKUP=NDATA-1

      RETURN
      END

C Similar function to ELOOKUP, but for cos(thetaMu).
      INTEGER FUNCTION TLOOKUP(T,E,NDATA,TMIN_DATA,TMAX_DATA)

      REAL    T, E
      INTEGER I, NDATA
      REAL    TMIN_DATA,TMAX_DATA

      IF(T.LT.TEVAL(1,E,TMIN_DATA,TMAX_DATA,NDATA)) THEN
        TLOOKUP=0
        RETURN
      ENDIF

      DO 21 I = 1,(NDATA-1)
          IF((T.GE.TEVAL(I,E,TMIN_DATA,TMAX_DATA,NDATA)).AND.
     $       (T.LT.TEVAL(I+1,E,TMIN_DATA,TMAX_DATA,NDATA))) THEN
          TLOOKUP=I
          RETURN
        ENDIF
 21   CONTINUE

      TLOOKUP=NDATA-1

      RETURN
      END

C Linear extrapolation function:
C Extrapolate value of X at Y, given points (X1,Y1) and (X2,Y2)
      REAL FUNCTION EXTPLT(X,X1,X2,Y1,Y2)

      REAL X,X1,X2,Y1,Y2

      EXTPLT=(((Y1-Y2)/(X1-X2))*X)
     $       +(((Y2*X1)-(Y1*X2))/(X1-X2))

      RETURN
      END



C Evaluate Ith Tmu in data
      REAL FUNCTION PEVAL(I,ENU,PMIN_DATA,PMAX_DATA,NDATA)

C ek = ratio of point density at finish to point density at start of Tmu range
      REAL ek
      INTEGER I, NDATA
      REAL ENU, PMIN_DATA, PMAX_DATA
      ek=10

      PEVAL= (((PMAX_DATA-PMIN_DATA)/log(ek))*
     $       log(((REAL(I-1)*(ek-1.))/(REAL(NDATA)-1.))+1.))
     $       +PMIN_DATA

      RETURN
      END




C Evaluate Ith cosTheta in data
      REAL FUNCTION TEVAL(I,ENU,TMIN_DATA,TMAX_DATA,NDATA)

C ek = ratio of point density at finish to point density at start of Tmu range
      REAL ek
      INTEGER I, NDATA
      REAL ENU, TMIN_DATA, TMAX_DATA
      ek=25

      TEVAL= (((TMAX_DATA-TMIN_DATA)/log(ek))*
     $       log(((REAL(I-1)*(ek-1.))/(REAL(NDATA)-1.))+1.))
     $       +TMIN_DATA

      RETURN
      END




