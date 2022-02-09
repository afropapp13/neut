************************************************************************
*     -----------------
      SUBROUTINE NECARD
*     -----------------
*
*     (Purpose)
*       Read NEUT card and setup common blocks.
*
*     (Input)
*       card file
*
*     (Output)
*
*     (Creation Date and Author)
*       1995.1.30 ; Koshio ( SGCARD )
*       1995.11.17 ; K.Kaneyuki
*       1999.10.11 ; D.Casper - allow user to specify #years in card file
*       2003.03.30 ; add neutparams.h(PFSURF,PFMAX,VNUINI,VNUFIN)
*       2006.12.30 ; Y.Hayato remove flux and exposure related part
*       2007.08.23 ; G.Mitsuka add Delta->gamma decay mode
*       2007.11.05 ; G.Mitsuka add nesettarg.F for nuceff
*                                  nrsettarg.F for nuccorspl
*       2008.08.13 ; Y.Hayato add FEF*** to control mean free path of pis
*       2009.01.11 ; Y.Hayato Default potential is calculated
*                                         using Fermi surface momentum
*       2009.04.11 ; Y.Hayato Move the location to call default parameter
*                             initialization routine (nesetfgparams)
*       2010.06    ; P.de Perio add flag for nucleon ejection
*                               add mfp model select for low energy pions
*                               add mfp model select for high energy pions
*                               add kinematical model select for hinrg pi
*                               add more factors to FEF
*                               added flag to check if nu sim or other sim
*       2013.05    ; F.d.M. Blaszczyk add formation zone parameter FZMU2
*                                     default value for FZMU2 is 0.08
*       2013.10    ; P.Rodrigues & A.Bercellieadd
*                              new form facor for Rein&Sehgal 1pi
*       2016.07    ; C. Bronner
*
*                    Hadron multiplicity model for multi-pi mode
*       2016.10    ; C.Wret add Rein-Sehgal pion ejection methods (MDLSPIEJ)
*                    and Minoo's model. Added background MA for Minoo
*
*       2019.04    ; Y.Hayato add new variable for LFG parameters
*
************************************************************************
      IMPLICIT NONE

      REAL*4 SPACE
      COMMON /CFREAD/SPACE(1000)

#include "cardname.h"
#include "necard.h"
#include "neutparams.h"
#include "neutmodel.h"
#include "neutfilepath.h"

      CHARACTER*40 TMPFIL, COMMAND
      CHARACTER*1024 CRSPATH_ENV

      INTEGER*4     IRANDFLG

      INTEGER*4     I,IRET

      INTEGER*4     IARGC
#ifndef gFortran
      EXTERNAL      IARGC
#endif

      LOGICAL  CARDISTOML
      EXTERNAL CARDISTOML
      
      REAL*4        XMN

      REAL*4        DEF_MARES_GS, DEF_CA5_GS, DEF_MARES_MK, DEF_CA5_MK

      INTEGER IRSRESMA, INFFRESMA, IQEMA

      integer*4     LENCRSPATH

      integer*4     LENSTR
      external      LENSTR

C     Opening file error
      INTEGER       OPENERR

      logical HaveSpecFuncForNuclei
      external HaveSpecFuncForNuclei

C --  Initialize

      CALL MCMASSGV(2112, XMN)

      NEFRMFLG = 0
      NEPAUFLG = 0
      NENEFO16 = 0
      NUCDEXITE = 1
      NENEFMODH = 1
      NENEFMODL = 1
      NENEFKINH = 1
      NEMODFLG = 0
      CALL VZERO(CRSNEUT,30)
      CALL VZERO(CRSNEUTB,30)
      IRANDFLG = 0
      NECOHEPI = 2
      NEDIFPI  = 0
      NEPDF    = 12
      NEBODEK  = 1
      NEBYLW   = 0
      NEBYFORCET1 = 0
      NEMULT   = 0
      NEBYSEAQFAC = 1
      NEBYRLONGMOD = 0

      XMAQE   = 1.21
      XMVQE   = 0.84

C-NCEL
      XMANCEL  = 1.21

C-- MEC
      MDL2P2H = 2
C     NV2P2HQVAL 1 is original Valencia, 2 is FS Eb
      NV2P2HQVAL = 1

C      XMASPI  = 1.21
      XMARSRES = 1.21
      XMVRSRES = 0.84

C     Use G&S form factors ( 0: Rein-Sehgal, 1 : G&S )
      NEIFF     = 1
C     If not G&S form factors, use R&S with quark oscillator level (=0)
      NENRTYPE  = 0
C     Pion ejection in Rein-Sehgal model. 
C     0 use no resonances (isotropic)
C     1 use Delta only
C     2 use multiple resonances
C     3 use Delta+isotropic (old NEUT default)
      MDLSPIEJ = 3

      ROP3P3_MULT = 1
      ROM3M3_MULT = 1
      ROP1P1_MULT = 1
      ROM1M1_MULT = 1
      ROP3P1_MULT = 1
      ROM1M3_MULT = 1
      ROP3M1_MULT = 1
      ROP1M3_MULT = 1

C     Default pion model is Rein-Sehgal + Global Fermi-gas = 1
C     Default pion model is Rein-Sehgal + Local  Fermi-gas = 2 ( not yet )
C     Minoo's model                                        = 3
      MDLSPI   = 1
C     G&S Defaults
      XMVNFFRES = 0.84
      DEF_MARES_GS = 0.95
      DEF_CA5_GS = 1.01
      RNEBGSCL  = 1.30

C     For Minoo's model we also tune the MA of the background. Here use
C     the parameters saved in neutmodel.h directly
      DEF_MARES_MK = 1.02
      DEF_CA5_MK = 1.06
      XMABKGM = 1.03

C     Default values so we can check if they were set.
      RNECA5I   = -999
      XMANFFRES = -999

      KAPP    = 1.0
      XMACOH  = 1.0
      RAD0NU  = 1.0
      fA1COH  = 0.0
      fb1COH  = 0.0

C Default CCQE model and axial form factor parameterisation
      MDLQE     = 2002
      MDLQEAF   = 1

C 2/3 component, TN 315 p 25
C ANL/BNL 'fixed' fit
      AXFFALPHA = 0.23
      AXFFBETA = 1.22
      AXFFGAMMA = 0.44
      AXFFTHETA = -0.26
C ANL/BNL 'free' fit
C       AXFFALPHA = 0.11
C       AXFFBETA = 1.04
C       AXFFGAMMA = 0.41
C       AXFFTHETA = -0.26

C Z-expansion TN 315 p 21 <3GeV as recommended
C Use official z-expansion parameters, not Patrick's results
      AXZEXPQ4 = 0
      AXZEXPNT = 8
      AXZEXPT0 = -0.280
      AXZEXPTC = 0.1764

C This is automatically set
      AXZEXPA0 = 0.5
C Up to the user
      AXZEXPA1 = 2.28
      AXZEXPA2 = 0.26
      AXZEXPA3 = -5.20
      AXZEXPA4 = 2.60
C These are automatically set
      AXZEXPA5 = 0.5
      AXZEXPA6 = 0.5
      AXZEXPA7 = 0.5
      AXZEXPA8 = 0.5
      AXZEXPA9 = 0.5
C End of z-expansion parameters

      NVQERFG        = 0
      NVQEBIND       = 1

      NVQERPA        = 1
      XNVRPAFP0IN    = 0.33
      XNVRPAPF0EX    = 0.45
      XNVRPAFSTAR    = 2.13
      XNVRPAF        = 1.
      XNVRPAPILAMBDA = 1200.
      XNVRPACR0      = 2.0
      XNVRPARHOLAMBDA= 2500
      XNVRPAGP       = 0.63
      XNVRPAXMPI     = 139.57
      XNVRPAXMRHO    = 775.4
      XNVRPAIREL     = 1.
      NVBINDFERMICOR = 1

      NEABSPIEMIT = 1

      PFSURF =  0.225
      PFMAX  =  0.225
C      VNUINI = -0.027 ! defined below
      VNUFIN = -0.000
      IFORMLEN= 110 ! default ; skip CCQE/Elastic
      FZMU2 = 0.08

      SFEBSHIFT = 0
      SFEBNEGBEH = 0 ! by default set Eb to physical limit and continue
C       SFEBNEGBEH = 1 ! Draw from SF again if Eb unphysical

      NUSIM   = 1
      NEFKINVER= 0
      QUIET = 2
      ITAUFLGCORE = 0           ! dummy

      NUMBNDN = 8
      NUMBNDP = 8
      NUMFREP = 2
      NUMATOM = 16

C     Paramters of new model tuned to pi-scattering data
      FEFQE   = 1.069
      FEFQEH  = 1.824
      FEFINEL = 1.002
      FEFABS  = 1.404
      FEFCOH  = 1.
      FEFQEHF = 1.
      FEFCOHF = 0.
      FEFCX   = 0.697
      FEFCXH  = 1.800
      FEFCXHF = 0.
      FEFCOUL = 0
      FEFALL  = 1.

C     Pion less delta decay : with MEC, expected to be 0
C      IPILESSDCY = 1
C      RPILESSDCY = 0.20E0
      IPILESSDCY = 0
      RPILESSDCY = 0.00E0

C     Cross-section data table PATH
C     Check environemnt variable first.
      CALL GETENV('NEUT_CRSPATH',CRSPATH_ENV)
      WRITE(*,*) "CRSPATH_ENV: ", CRSPATH_ENV(1:LEN_TRIM(CRSPATH_ENV))
      IF ( LEN_TRIM(CRSPATH_ENV) .GT. 0) THEN
        CRSTBLPATH = CRSPATH_ENV
      ELSE
        CRSTBLPATH = './'
      ENDIF

C     FermiMomentum used for pauli blocking in SF model only
C     (NOTE: default is -1, gets set elsewhere correctly if not user-defined)
      PFSF = -1

C     relatve size of second class currents (vector and axial)
      SCCFV = 0.0
      SCCFA = 0.0

C     RT: adding error parameter for FP(Q2=0)
      FPQE = 1.0

C     Radiative correction on/off ( default off )
      IRADCORR = 0

C     RT: adding diffractive pion parameters
      XMADIF = 1.1
      NUCVOLDIF = 7.0

C --- get card file name (1st argument of the program)

      CALL RESOLVECARDNAME
      WRITE(6,*) '======= CONFIGURATION FOR INTERACTIONS ======='
      IF ( CARDISTOML() ) THEN
        WRITE(6,*) 'Reading', FCARD(1:LEN_TRIM(FCARD)), ' as a toml file.'
        CALL NECARDTOML
      ELSE
        WRITE(6,*) 'Reading', FCARD(1:LEN_TRIM(FCARD)), ' as a card file.'

C -- extract NEUT card
        CALL GETUNI(TMPFIL)
        OPEN(UNIT=90, FILE=FCARD, STATUS='OLD')
        OPEN(UNIT=91, FILE=TMPFIL,STATUS='UNKNOWN')
        CALL EXCARD('NEUT', 90, 91, IRET)
        CLOSE(UNIT=91)
        CLOSE(UNIT=90)

        OPEN(UNIT=91, FILE=TMPFIL,STATUS='OLD')
        CALL FFINIT(1000)
        CALL FFSET('LINP', 91)
        CALL FFSET('SIZE', 32)
        CALL FFSET('LENG', 100)
        CALL FFKEY('CRSPATH',CRSTBLPATH,80,'CHAR')
        CALL FFKEY('FERM',NEFRMFLG, 1,'INTEGER')
        CALL FFKEY('PAUL',NEPAUFLG, 1,'INTEGER')
        CALL FFKEY('NEFF',NENEFO16, 1,'INTEGER')
        CALL FFKEY('NUCDEXITE',NUCDEXITE, 1,'INTEGER')
        CALL FFKEY('MODH',NENEFMODH, 1,'INTEGER')
        CALL FFKEY('MODL',NENEFMODL, 1,'INTEGER')
        CALL FFKEY('KINH',NENEFKINH, 1,'INTEGER')
        CALL FFKEY('MODE',NEMODFLG, 1,'INTEGER')
        CALL FFKEY('COHEPI',NECOHEPI, 1,'INTEGER')
        CALL FFKEY('DIFPI',NEDIFPI, 1,'INTEGER')
        CALL FFKEY('PDF',NEPDF, 1,'INTEGER')
        CALL FFKEY('BODEK',NEBODEK, 1,'INTEGER')
        CALL FFKEY('MULT',NEMULT, 1,'INTEGER')
        CALL FFKEY('BYLW',NEBYLW, 1,'INTEGER')
        CALL FFKEY('BYFORCET1',NEBYFORCET1, 1,'INTEGER')
        CALL FFKEY('BYSEAQFAC',NEBYSEAQFAC, 1,'REAL')
        CALL FFKEY('BYRLONGMOD',NEBYRLONGMOD, 1,'REAL')

        CALL FFKEY('IFF',NEIFF, 1,'INTEGER')
        CALL FFKEY('NRTYPE',NENRTYPE, 1,'INTEGER')
        CALL FFKEY('MDLSPIEJ',MDLSPIEJ, 1,'INTEGER')
        CALL FFKEY('MDLSPI',MDLSPI, 1,'INTEGER')

        CALL FFKEY('MARSRES'  ,XMARSRES ,1,'REAL')
        CALL FFKEY('MARNFFRES',XMANFFRES,1,'REAL')
        CALL FFKEY('MVRSRES' ,XMVRSRES ,1,'REAL')
        CALL FFKEY('MVNFFRES' ,XMVNFFRES ,1,'REAL')
        CALL FFKEY('MABKGM',XMABKGM ,1,'REAL')
        CALL FFKEY('CA5', RNECA5I, 1,'REAL')
        CALL FFKEY('BGSCL',RNEBGSCL, 1,'REAL')

        CALL FFKEY('MAQE' ,XMAQE ,1,'REAL')
        CALL FFKEY('MVQE' ,XMVQE ,1,'REAL')
        CALL FFKEY('KAPP' ,KAPP ,1,'REAL')

        CALL FFKEY('MANCEL' ,XMANCEL ,1,'REAL')

        CALL FFKEY('MACOH' ,XMACOH ,1,'REAL')
        CALL FFKEY('R0COH' ,RAD0NU ,1,'REAL')
        CALL FFKEY('fA1COH',fA1COH ,1,'REAL')
        CALL FFKEY('fb1COH',fb1COH ,1,'REAL')

        CALL FFKEY('MDLQE',MDLQE,1,'INTEGER')
        CALL FFKEY('MDLQEAF',MDLQEAF,1,'INTEGER')

        CALL FFKEY('NVQERFG',NVQERFG,1,'INTEGER')
        CALL FFKEY('NVQEBIND',NVQEBIND,1,'INTEGER')

        CALL FFKEY('NVQERPA',NVQERPA,1,'INTEGER')

        CALL FFKEY('NVRPAFP0IN'    ,XNVRPAFP0IN,    1,'REAL')
        CALL FFKEY('NVRPAPF0EX'    ,XNVRPAPF0EX,    1,'REAL')
        CALL FFKEY('NVRPAFSTAR'    ,XNVRPAFSTAR,    1,'REAL')
        CALL FFKEY('NVRPAF'        ,XNVRPAF,        1,'REAL')
        CALL FFKEY('NVRPAPILAMBDA' ,XNVRPAPILAMBDA, 1,'REAL')
        CALL FFKEY('NVRPACR0'      ,XNVRPACR0,      1,'REAL')
        CALL FFKEY('NVRPARHOLAMBDA',XNVRPARHOLAMBDA,1,'REAL')
        CALL FFKEY('NVRPAGP'       ,XNVRPAGP,       1,'REAL')
        CALL FFKEY('NVRPAXMPI'     ,XNVRPAXMPI,     1,'REAL')
        CALL FFKEY('NVRPAXMRHO'    ,XNVRPAXMRHO,    1,'REAL')
        CALL FFKEY('NVRPAIREL'     ,XNVRPAIREL,     1,'REAL')

        CALL FFKEY('NVBINDFERMICOR',NVBINDFERMICOR, 1,'INTEGER')

        CALL FFKEY('MDL2P2H',MDL2P2H,1,'INTEGER')
        CALL FFKEY('NV2P2HQVAL',NV2P2HQVAL,1,'INTEGER')

        CALL FFKEY('CRS' ,CRSNEUT ,30,'REAL')
        CALL FFKEY('CRSB',CRSNEUTB,30,'REAL')
        CALL FFKEY('RAND',IRANDFLG,1,'INTEGER')

        CALL FFKEY('NUMBNDN',NUMBNDN,1,'INTEGER')
        CALL FFKEY('NUMBNDP',NUMBNDP,1,'INTEGER')
        CALL FFKEY('NUMFREP',NUMFREP,1,'INTEGER')
        CALL FFKEY('NUMATOM',NUMATOM,1,'INTEGER')

        CALL FFKEY('FEFQE'  ,FEFQE,1,'REAL')
        CALL FFKEY('FEFQEH' ,FEFQEH,1,'REAL')
        CALL FFKEY('FEFINEL',FEFINEL,1,'REAL')
        CALL FFKEY('FEFABS',FEFABS,1,'REAL')
        CALL FFKEY('FEFCOH',FEFCOH,1,'REAL')
        CALL FFKEY('FEFQEHF',FEFQEHF,1,'REAL')
        CALL FFKEY('FEFCOHF',FEFCOHF,1,'REAL')
        CALL FFKEY('FEFCX',FEFCX,1,'REAL')
        CALL FFKEY('FEFCXH',FEFCXH,1,'REAL')
        CALL FFKEY('FEFCXHF',FEFCXHF,1,'REAL')
        CALL FFKEY('FEFCOUL',FEFCOUL,1,'REAL')
        CALL FFKEY('FEFALL',FEFALL,1,'REAL')

        CALL FFKEY('XMADIF',XMADIF,1,'REAL')
        CALL FFKEY('NUCVOLDIF',NUCVOLDIF,1,'REAL')

        CALL FFKEY('IPILESSDCY',IPILESSDCY,1,'INTEGER')
        CALL FFKEY('RPILESSDCY',RPILESSDCY,1,'REAL')

        CALL FFKEY('PFSF',PFSF,1,'REAL')

        CALL FFKEY('SCCFV',SCCFV,1,'REAL')
        CALL FFKEY('SCCFA',SCCFA,1,'REAL')

        CALL FFKEY('FPQE',FPQE,1,'REAL')

        CALL FFGO
        CLOSE(UNIT=91)
        COMMAND = '/bin/rm '//TMPFIL
        CALL SYSTEM(COMMAND)

        call nesetfgparams

C----------------- Loop twice to override default parameters

        CALL GETUNI(TMPFIL)
        OPEN(UNIT=90, FILE=FCARD, STATUS='OLD')
        OPEN(UNIT=91, FILE=TMPFIL,STATUS='UNKNOWN')
        CALL EXCARD('NEUT', 90, 91, IRET)
        CLOSE(UNIT=91)
        CLOSE(UNIT=90)

        OPEN(UNIT=91, FILE=TMPFIL,STATUS='OLD')
        CALL FFINIT(0)
        CALL FFSET('LINP', 91)
        CALL FFSET('SIZE', 9)
        CALL FFSET('LENG', 100)
C-- Initialize parameters here to allow override

        CALL FFKEY('PFSURF',PFSURF,1,'REAL')
        CALL FFKEY('PFMAX',PFMAX,1,'REAL')

C Set VNUINI here ( when PFSURF is fixed)
        VNUINI = -1. * (sqrt(XMN * XMN + PFSURF * PFSURF) - XMN)

        CALL FFKEY('VNUINI',VNUINI,1,'REAL')
        CALL FFKEY('VNUFIN',VNUFIN,1,'REAL')

C       These could be above, but FFKEY complains about running out of
C       space, so read these on the second pass.
C       Alt AXFF params
        CALL FFKEY('AXFFALPHA',AXFFALPHA, 1, 'REAL')
        CALL FFKEY('AXFFBETA',AXFFBETA, 1, 'REAL')
        CALL FFKEY('AXFFGAMMA',AXFFGAMMA, 1, 'REAL')
        CALL FFKEY('AXFFTHETA',AXFFTHETA, 1, 'REAL')

        CALL FFKEY('AXZEXPQ4',AXZEXPQ4, 1, 'INTEGER')
        CALL FFKEY('AXZEXPNT',AXZEXPNT, 1, 'INTEGER')
        CALL FFKEY('AXZEXPT0',AXZEXPT0, 1, 'REAL')
        CALL FFKEY('AXZEXPTC',AXZEXPTC, 1, 'REAL')
        CALL FFKEY('AXZEXPA0',AXZEXPA0, 1, 'REAL')
        CALL FFKEY('AXZEXPA1',AXZEXPA1, 1, 'REAL')
        CALL FFKEY('AXZEXPA2',AXZEXPA2, 1, 'REAL')
        CALL FFKEY('AXZEXPA3',AXZEXPA3, 1, 'REAL')
        CALL FFKEY('AXZEXPA4',AXZEXPA4, 1, 'REAL')
        CALL FFKEY('AXZEXPA5',AXZEXPA5, 1, 'REAL')
        CALL FFKEY('AXZEXPA6',AXZEXPA6, 1, 'REAL')
        CALL FFKEY('AXZEXPA7',AXZEXPA7, 1, 'REAL')
        CALL FFKEY('AXZEXPA8',AXZEXPA8, 1, 'REAL')
        CALL FFKEY('AXZEXPA9',AXZEXPA9, 1, 'REAL')
        
        CALL FFKEY('IFORMLEN',IFORMLEN,1,'INTEGER')
        CALL FFKEY('FZMU2',FZMU2,1,'REAL')
        CALL FFKEY('ABSPIEMIT',NEABSPIEMIT,1,'INTEGER')
        CALL FFKEY('NUSIM',NUSIM,1,'INTEGER')
        CALL FFKEY('EFVER',NEFKINVER,1,'INTEGER')
        CALL FFKEY('IRADCORR',IRADCORR,1,'INTEGER')
        CALL FFKEY('QUIET',QUIET,1,'INTEGER')

        CALL FFKEY('SFEBSHIFT',SFEBSHIFT,1,'REAL')
        CALL FFKEY('SFEBNEGBEH',SFEBNEGBEH,1,'INTEGER')

        CALL FFGO
        CLOSE(UNIT=91)
        COMMAND = '/bin/rm '//TMPFIL
        CALL SYSTEM(COMMAND)
      ENDIF

C     If user has not specified CA5
      IF (RNECA5I .LT. 0) THEN
C       For Rein-Sehgal
        IF (mdlspi .EQ. 1) THEN
          RNECA5I = DEF_CA5_GS
C       For Minoo
        ELSE IF (mdlspi .EQ. 3) THEN
          RNECA5I = DEF_CA5_MK
        ENDIF
      ENDIF

C     If user has not specified MARES
      IF (XMANFFRES .LT. 0) THEN
C       For Rein-Sehgal
        IF (MDLSPI .EQ. 1) THEN
          XMANFFRES = DEF_MARES_GS
C       For Minoo
        ELSE IF (MDLSPI .EQ. 3) THEN
          XMANFFRES = DEF_MARES_MK
        ENDIF
      ENDIF

      IF (mod(MDLQE, 1000)/100 .EQ. 4) THEN
        IF (HaveSpecFuncForNuclei()) THEN
          write (6,*) "Using Spectral Function for elastic events"
          GOTO 200
        ELSE
          write (6,*) "Spectral Function not available"
          write (6,*) "- use RFG (NEUT-MDLQE 002)"
          STOP
        ENDIF
      ENDIF
200   CONTINUE

      LENCRSPATH = LENSTR(CRSTBLPATH)
      if (LENCRSPATH.eq.0) THEN
        write(*,*) 'Cross section table path undefined.'
        write(*,*) 'Set NEUT_CRSPATH environment variable'
        write (*,*) 'or NEUT-CRSPATH in card file.'
        STOP
      endif
      if (CRSTBLPATH(LENCRSPATH:LENCRSPATH).ne.'/') THEN
         CRSTBLPATH=CRSTBLPATH(1:LENCRSPATH)//'/'
         LENCRSPATH = LENCRSPATH + 1
      endif

      WRITE(6,*) 'NECARD : INPUTTED PARAMETER '
      WRITE(6,*) 'NECARD : CRSPATH= ',CRSTBLPATH(1:LENCRSPATH)
      WRITE(6,*) 'NECARD : FERM   = ',NEFRMFLG
      WRITE(6,*) 'NECARD : PAUL   = ',NEPAUFLG
      WRITE(6,*) 'NECARD : NEFF   = ',NENEFO16
      WRITE(6,*) 'NECARD : NUCDEXITE = ',NUCDEXITE
      WRITE(6,*) 'NECARD : MODL   = ',NENEFMODL
      WRITE(6,*) 'NECARD : MODH   = ',NENEFMODH
      WRITE(6,*) 'NECARD : KINH   = ',NENEFKINH
      WRITE(6,*) 'NECARD : MODE   = ',NEMODFLG
      WRITE(6,*) 'NECARD : COHEPI = ',NECOHEPI
      WRITE(6,*) 'NECARD : DIFPI  = ',NEDIFPI

      WRITE(6,*) '********* CCQE SETTINGS ***************'
      WRITE(6,*) 'NECARD : MDLQE   = ',MDLQE
      WRITE(6,*) 'NECARD : MDLQEAF = ',MDLQEAF
      WRITE(6,*) 'NECARD : MDL2P2H = ',MDL2P2H
      WRITE(6,*) 'NECARD : NV2P2HQVAL = ',NV2P2HQVAL
      WRITE(6,*) 'NECARD : MAQE    = ',XMAQE
      WRITE(6,*) 'NECARD : AXFFALPHA =',AXFFALPHA
      WRITE(6,*) 'NECARD : AXFFBETA=',AXFFBETA
      WRITE(6,*) 'NECARD : AXFFGAMMA=',AXFFGAMMA
      WRITE(6,*) 'NECARD : AXFFTHETA=',AXFFTHETA
C Clarence hack, z exp
      WRITE(6,*) 'NECARD : AXZEXPQ4=',AXZEXPQ4
      WRITE(6,*) 'NECARD : AXZEXPNT=',AXZEXPNT
      WRITE(6,*) 'NECARD : AXZEXPT0=',AXZEXPT0
      WRITE(6,*) 'NECARD : AXZEXPTC=',AXZEXPTC
      WRITE(6,*) 'NECARD : AXZEXPA0=',AXZEXPA0
      WRITE(6,*) 'NECARD : AXZEXPA1=',AXZEXPA1
      WRITE(6,*) 'NECARD : AXZEXPA2=',AXZEXPA2
      WRITE(6,*) 'NECARD : AXZEXPA3=',AXZEXPA3
      WRITE(6,*) 'NECARD : AXZEXPA4=',AXZEXPA4
      WRITE(6,*) 'NECARD : AXZEXPA5=',AXZEXPA5
      WRITE(6,*) 'NECARD : AXZEXPA6=',AXZEXPA6
      WRITE(6,*) 'NECARD : AXZEXPA7=',AXZEXPA7
      WRITE(6,*) 'NECARD : AXZEXPA8=',AXZEXPA8
      WRITE(6,*) 'NECARD : AXZEXPA9=',AXZEXPA9
      WRITE(6,*) 'NECARD : MVQE   = ',XMVQE
      WRITE(6,*) 'NECARD : KAPP   = ',KAPP
      WRITE(6,*) 'NECARD : SCCFV     = ',SCCFV
      WRITE(6,*) 'NECARD : SCCFA     = ',SCCFA
      WRITE(6,*) 'NECARD : FPQE      = ',FPQE
      WRITE(6,*) 'NECARD : PFSF      = ',PFSF
      WRITE(6,*) 'NECARD : IRADCORR  = ',IRADCORR
      WRITE(6,*) 'NECARD : IPILESSDCY= ',IPILESSDCY
      WRITE(6,*) 'NECARD : RPILESSDCY= ',RPILESSDCY

      WRITE(6,*) 'NECARD : MANCEL = ',XMANCEL

      WRITE(6,*) '********** SINGLE PION ***************'
      WRITE(6,*) 'NECARD : IFF       = ',NEIFF
      WRITE(6,*) 'NECARD : NRTYPE    = ',NENRTYPE
      WRITE(6,*) 'NECARD : MDLSPIEJ  = ',MDLSPIEJ
      WRITE(6,*) 'NECARD : MDLSPI    = ',MDLSPI
      WRITE(6,*) 'NECARD : MARSRES   = ',XMARSRES
      WRITE(6,*) 'NECARD : MANFFRES  = ',XMANFFRES
      WRITE(6,*) 'NECARD : MABKGM    = ',XMABKGM
      WRITE(6,*) 'NECARD : CA5       = ',RNECA5I
      WRITE(6,*) 'NECARD : BGSCL     = ',RNEBGSCL
      WRITE(6,*) 'NECARD : MVRSRES   = ',XMVRSRES
      WRITE(6,*) 'NECARD : MVNFFRES  = ',XMVNFFRES

      WRITE(6,*) '********* MULTIPI SETTINGS ***************'
      WRITE(6,*) 'NECARD : PDF    = ',NEPDF
      WRITE(6,*) 'NECARD : BODEK  = ',NEBODEK
      WRITE(6,*) 'NECARD : MULT   = ',NEMULT
      WRITE(6,*) 'NECARD : BYLW   = ',NEBYLW
      WRITE(6,*) 'NECARD : BYFORCET1   = ',NEBYFORCET1

      WRITE(6,*) '********** COHERENT *************'
      WRITE(6,*) 'NECARD : COHEPI   = ',NECOHEPI
      WRITE(6,*) 'NECARD : MACOH    = ',XMACOH
      WRITE(6,*) 'NECARD : R0COH    = ',RAD0NU
      WRITE(6,*) 'NECARD : A1COHerr = ',fA1COH
      WRITE(6,*) 'NECARD : b1COHerr = ',fb1COH
      WRITE(6,*) 'NECARD : XMADIF   = ',XMADIF
      WRITE(6,*) 'NECARD : NUCVOLDIF= ',NUCVOLDIF

      WRITE(6,*) '======= PION FSI ======='
      WRITE(6,*) 'NECARD : FEFQE     = ',FEFQE
      WRITE(6,*) 'NECARD : FEFQEH    = ',FEFQEH
      WRITE(6,*) 'NECARD : FEFINEL   = ',FEFINEL
      WRITE(6,*) 'NECARD : FEFABS    = ',FEFABS
      WRITE(6,*) 'NECARD : FEFCOH    = ',FEFCOH
      WRITE(6,*) 'NECARD : FEFCX     = ',FEFCX
      WRITE(6,*) 'NECARD : FEFCXH    = ',FEFCXH
      WRITE(6,*) 'NECARD : FEFCXHF   = ',FEFCXHF
      WRITE(6,*) 'NECARD : FEFQEHF   = ',FEFQEHF
      WRITE(6,*) 'NECARD : FEFCOHF   = ',FEFCOHF
      WRITE(6,*) 'NECARD : FEFCOUL   = ',FEFCOUL
      WRITE(6,*) 'NECARD : FEFALL    = ',FEFALL

      IF (NEMODFLG.EQ.-1) THEN
         WRITE(6,*) 'NECARD : CRS  = ',CRSNEUT
         WRITE(6,*) 'NECARD : CRSB = ',CRSNEUTB
      ELSE
         DO 10 I=1,30
            CRSNEUT(I) = 1.0
            CRSNEUTB(I)= 1.0
 10      CONTINUE
      ENDIF

      WRITE(6,*) 'NECARD : NUCEFF KINEMATICS VERSION  = ',NEFKINVER

      IF (IRANDFLG.eq.0) THEN
         write(6,*) "NECARD :Read RANDOM number from FILE"
         CALL NERDSEED
      ELSE
         write(6,*) "NECARD :Generating RANDOM SEED"
         CALL NEGENSEED
      ENDIF

      WRITE(6,*) '======= CONFIGURATION FOR TARGET MATERIAL ======='
      WRITE(6,*) 'NECARD : NUMBNDN = ',NUMBNDN
      WRITE(6,*) 'NECARD : NUMBNDP = ',NUMBNDP
      WRITE(6,*) 'NECARD : NUMFREP = ',NUMFREP
      WRITE(6,*) 'NECARD : NUMATOM = ',NUMATOM

      WRITE(6,*) 'NECARD : PFSURF = ',PFSURF
      WRITE(6,*) 'NECARD : PFMAX  = ',PFMAX
      WRITE(6,*) 'NECARD : VNUINI = ',VNUINI
      WRITE(6,*) 'NECARD : VNUFIN = ',VNUFIN
      WRITE(6,*) 'NECARD : SFEBSHIFT = ',SFEBSHIFT
      WRITE(6,*) 'NECARD : IFORMLEN = ',IFORMLEN
      WRITE(6,*) 'NECARD : FZMU2 = ',FZMU2
      WRITE(6,*) 'NECARD : ABSPIEMIT = ',NEABSPIEMIT
      WRITE(6,*) 'NECARD : NUSIM = ',NUSIM
      WRITE(6,*) 'NECARD : QUIET = ',QUIET

C     Parameter checks

      IQEMA = INT(XMAQE*100+0.5)
C      ISPIMA = INT(XMASPI*100+0.5)
      IRSRESMA  = INT(XMARSRES*100+0.5)
      INFFRESMA = INT(XMANFFRES*100+0.5)

      if (mod(MDLQE, 10000)/1000 .NE. 2) then
         IF ( IQEMA.ne.101 .and. IQEMA.ne.103 .and.
     $        IQEMA.ne.105 .and.
     $        IQEMA.ne.111 .and. IQEMA.ne.121 .and.
     &        IQEMA.ne.131 .and. IQEMA.ne.160 )
     $        THEN
            WRITE(6,*)
     &           'NECARD error: MAQE must be 1.01, 1.03, 1.05',
     $           ' 1.[1-3]1 or 1.60'
            STOP
         END IF
      endif

      IF (IRSRESMA.ne.95 .and. IRSRESMA.ne.101 .and. IRSRESMA.ne.103.
     $     .and. IRSRESMA.ne.111 .and. IRSRESMA.ne.121 ) THEN
         WRITE(6,*)
     $     'NECARD error: MARS must be 0.95, 1.01, 1.03, 1.11 or 1.21'
         STOP
      END IF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Check 1 pi model: can only generate with certain parameters
C     because we read a pre-generated cross-section table
C     For Rein-Sehgal original form factors
      if (MDLSPI .EQ. 1 .AND. NEIFF .EQ. 0) THEN
      IF ( INFFRESMA.ne.95 .and. INFFRESMA.ne.101 .and.
     $     INFFRESMA.ne.103. .and. INFFRESMA.ne.111 .and.
     $     INFFRESMA.ne.121 ) THEN
         WRITE(6,*)
     $    'NECARD error: MARNFF must be 0.95, 1.01, 1.03, 1.11 or 1.21'
         STOP
      END IF
      ENDIF

C     Check Rein-Sehgal model with Graczyk Sobczyk form factors
C     Can only generate with certain parameters
C     because we read a pre-generated cross-section table
      IF (MDLSPI .EQ. 1 .AND. NEIFF .EQ. 1) THEN
        IF ( XMANFFRES .ne. 0.95 ) then
          WRITE(6,*) 
     $    'NECARD: MARNFFRES must be 0.95 for generating Rein-Sehgal'
          STOP
        ELSE IF (RNECA5I .ne. 1.01) then
          WRITE(6,*) 
     $    'NECARD: CA5 must be 1.01 for generating Rein-Sehgal'
          STOP
        ELSE IF (RNEBGSCL .ne. 1.30) then
          WRITE(6,*) 
     $    'NECARD: BGSCL must be 1.30 for generating Rein-Sehgal'
          STOP
        ENDIF
      ENDIF

C     For Minoo model do the checks
      IF (MDLSPI .EQ. 3) THEN
        IF (XMABKGM .ne. 1.03 ) then
          WRITE(6,*) 
     $    'NECARD: MABKGM must be 1.03 for generating Minoo model'
          STOP
        ELSE IF (RNECA5I .ne. 1.02) then
          WRITE(6,*) 
     $    'NECARD: CA5 must be 1.02 for generating Minoo model'
          STOP
        ELSE IF (XMANFFRES .ne. 1.06) then
          WRITE(6,*)
     $    'NECARD: XMANFFRES must be 1.06 for generating Minoo mode'
          STOP
        ENDIF
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Check FSI rescattering probabilities

      IF (FEFQEHF.lt.0 .or. FEFQEHF.gt.1) THEN
         WRITE(6,*) 'NECARD error: FEFQEHF must be between 0 and 1!'
         STOP
      END IF

      IF (FEFCXHF.lt.0 .or. FEFCXHF.gt.1) THEN
         WRITE(6,*) 'NECARD error: FEFCXHF must be between 0 and 1!'
         STOP
      END IF

      IF ((MDLSPI .NE. 1) .AND. (MDLSPIEJ .NE. 1)) THEN
        WRITE(*,*) '******************************'
        WRITE(*,*) 'NECARD error: specified a Rein-Sehgal ejection'
        WRITE(*,*) 'method with a model that WAS NOT Rein-Sehgal!'
        WRITE(*,*) 'Model: ', MDLSPI, 'Method: ', MDLSPIEJ
        STOP
      ENDIF

C     Check pion ejection model
      IF (MDLSPIEJ.LT.0 .OR. MDLSPIEJ.GT.3) then
        WRITE(*,*) '******************************'
        WRITE(*,*) 'NECARD error: specified a Rein-Sehgal ejection that
     &              doesnt exist'
        write(*,*) 'MDLSPIEJ: ', MDLSPIEJ, ' but should be'
        write(*,*) 'MDLSPIEJ=0 Isotropic'
        write(*,*) 'MDLSPIEJ=1 Delta'
        write(*,*) 'MDLSPIEJ=2 Many'
        write(*,*) 'MDLSPIEJ=3 Delta+Isotrop (old NEUT)'
        STOP
      ENDIF

C-- Add nuceff/nuccorspl related card
C--      call nesettarg
C--      call nrsettarg
C-- ADD proton nuceff related card
      call nrcard

      call nefillmodel

      RETURN

 123  write(*,*) 'Card file was not found: ', fcard
      stop

      END
