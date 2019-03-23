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
*
************************************************************************
      IMPLICIT NONE
      
      REAL*4 SPACE
      COMMON /CFREAD/SPACE(1000)

#include "necard.h"
#include "neutparams.h"
#include "neutmodel.h"
#include "neutfilepath.h"

      CHARACTER*40 TMPFIL, COMMAND

      CHARACTER*128 FCARD

      INTEGER*4     IRANDFLG

      INTEGER*4     I,IRET

      INTEGER*4     IARGC
#ifndef gFortran 
      EXTERNAL      IARGC
#endif

      REAL*4        XMN

      INTEGER IRSRESMA, INFFRESMA, IQEMA

      integer*4     LENCRSPATH

      integer*4     LENSTR
      external      LENSTR

C --  Initialize

      CALL MCMASSGV(2112, XMN)

      NEFRMFLG = 0
      NEPAUFLG = 0
      NENEFO16 = 0
      NENEFMODH = 1
      NENEFMODL = 1
      NENEFKINH = 1
      NEMODFLG = 0
      CALL VZERO(28,CRSNEUT)
      CALL VZERO(28,CRSNEUTB)
      IRANDFLG = 0
      NECOHEPI = 0
      NEPDF   = 12
      NEBODEK = 1

      XMAQE   = 1.21
C      XMAQE   = 1.03
      XMVQE   = 0.84
C      XMASPI  = 1.21
      XMARSRES = 1.21
      XMVRSRES = 0.84

C     For form factors for Rein-Seighal update
      XMANFFRES = 0.95
      XMVNFFRES = 0.84
      NEIFF     = 1
      NENRTYPE  = 0
      RNECA5I   = 1.01
      RNEBGSCL  = 1.30

      KAPP    = 1.0
      XMACOH  = 1.0
      RAD0NU  = 1.0
      fA1COH  = 0.0
      fb1COH  = 0.0 

      MDLQE     = 22
      MDLQEAF   = 1

      NEABSPIEMIT = 1
      
      PFSURF =  0.225
      PFMAX  =  0.225
C      VNUINI = -0.027 ! defined below
      VNUFIN = -0.000
      IFORMLEN= 1
      FZMU2 = 0.08

      NUSIM   = 1
      QUIET = 0
      ITAUFLGCORE = 0           ! dummy

      NUMBNDN = 8
      NUMBNDP = 8
      NUMFREP = 2
      NUMATOM = 16

C     Paramters of new model tuned to pi-scattering data
      FEFQE   = 1.
      FEFQEH  = 1.8
      FEFINEL = 1.
      FEFABS  = 1.1
      FEFCOH  = 1.
      FEFQEHF = 1.
      FEFCOHF = 0.
      FEFCX   = 1.
      FEFCXH  = 1.8
      FEFCXHF = 0.
      FEFCOUL = 0
      FEFALL  = 1.

C     Pion less delta decay : with MEC, expected to be 0
C      IPILESSDCY = 1
C      RPILESSDCY = 0.20E0
      IPILESSDCY = 0
      RPILESSDCY = 0.00E0

C     Cross-section data table PATH
      CRSTBLPATH ='./'

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

C --- get card file name (1st argument of the program)

      WRITE(6,*) 
      IF (IARGC() .GE. 1) THEN
         CALL GETARG(1, FCARD)
      ELSE
         FCARD = 'neut.card'
      END IF
      WRITE(6,*) '======= CONFIGURATION FOR INTERACTIONS ======='
      WRITE(6,*) 'NEUT card file you read is ',FCARD

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
      CALL FFKEY('MODH',NENEFMODH, 1,'INTEGER')
      CALL FFKEY('MODL',NENEFMODL, 1,'INTEGER')
      CALL FFKEY('KINH',NENEFKINH, 1,'INTEGER')
      CALL FFKEY('MODE',NEMODFLG, 1,'INTEGER')
      CALL FFKEY('COHEPI',NECOHEPI, 1,'INTEGER')
      CALL FFKEY('PDF',NEPDF, 1,'INTEGER')
      CALL FFKEY('BODEK',NEBODEK, 1,'INTEGER')

      CALL FFKEY('IFF',NEIFF, 1,'INTEGER')
      CALL FFKEY('NRTYPE',NENRTYPE, 1,'INTEGER')
      CALL FFKEY('CA5',RNECA5I, 1,'REAL')
      CALL FFKEY('BGSCL',RNEBGSCL, 1,'REAL')

      CALL FFKEY('MAQE' ,XMAQE ,1,'REAL')
C      CALL FFKEY('MASPI' ,XMASPI ,1,'REAL')
      CALL FFKEY('MARSRES'  ,XMARSRES ,1,'REAL')
      CALL FFKEY('MARNFFRES',XMANFFRES ,1,'REAL')
      CALL FFKEY('MVQE' ,XMVQE ,1,'REAL')
C      CALL FFKEY('MVSPI' ,XMVSPI ,1,'REAL')
      CALL FFKEY('MVRSRES' ,XMVRSRES ,1,'REAL')
      CALL FFKEY('MVNFFRES' ,XMVNFFRES ,1,'REAL')
      CALL FFKEY('KAPP' ,KAPP ,1,'REAL')
      CALL FFKEY('MACOH' ,XMACOH ,1,'REAL')
      CALL FFKEY('R0COH' ,RAD0NU ,1,'REAL')
      CALL FFKEY('fA1COH',fA1COH ,1,'REAL')
      CALL FFKEY('fb1COH',fb1COH ,1,'REAL')
 
      CALL FFKEY('MDLQE',MDLQE,1,'INTEGER')
      CALL FFKEY('MDLQEAF',MDLQEAF,1,'INTEGER')

      CALL FFKEY('CRS' ,CRSNEUT ,28,'REAL')
      CALL FFKEY('CRSB',CRSNEUTB,28,'REAL')
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
      CALL FFKEY('IFORMLEN',IFORMLEN,1,'INTEGER')
      CALL FFKEY('FZMU2',FZMU2,1,'REAL')
      CALL FFKEY('ABSPIEMIT',NEABSPIEMIT,1,'INTEGER')
      CALL FFKEY('NUSIM',NUSIM,1,'INTEGER')
      CALL FFKEY('IRADCORR',IRADCORR,1,'INTEGER')
      CALL FFKEY('QUIET',QUIET,1,'INTEGER')

      CALL FFGO
      CLOSE(UNIT=91)
      COMMAND = '/bin/rm '//TMPFIL
      CALL SYSTEM(COMMAND)

      IF (mod(MDLQE, 1000)/100 .EQ. 4) THEN
        IF (NUMATOM .EQ. 12) THEN
          write (6,*) "Using Spectral Function for elastic events"
          GOTO 200
        ELSE IF (NUMATOM .EQ. 16) THEN
          write (6,*) "Using Spectral Function for elastic events"
          GOTO 200
        ELSE IF (NUMATOM .EQ. 56) THEN
          write (6,*) "Using Spectral Function for elastic events"
          GOTO 200
        ELSE
          write (6,*) "Spectral Function not available - use RFG"
          MDLQE = MDLQE - 400
        ENDIF
      ENDIF
200   CONTINUE

      LENCRSPATH = LENSTR(CRSTBLPATH)
      if (CRSTBLPATH(LENCRSPATH:LENCRSPATH).ne.'/') THEN
         CRSTBLPATH=CRSTBLPATH(1:LENCRSPATH)//'/'
         LENCRSPATH = LENCRSPATH + 1
      endif

      WRITE(6,*) 'NECARD : INPUTTED PARAMETER '
      WRITE(6,*) 'NECARD : CRSPATH= ',CRSTBLPATH(1:LENCRSPATH)
      WRITE(6,*) 'NECARD : FERM   = ',NEFRMFLG
      WRITE(6,*) 'NECARD : PAUL   = ',NEPAUFLG
      WRITE(6,*) 'NECARD : NEFF   = ',NENEFO16
      WRITE(6,*) 'NECARD : MODL   = ',NENEFMODL
      WRITE(6,*) 'NECARD : MODH   = ',NENEFMODH
      WRITE(6,*) 'NECARD : KINH   = ',NENEFKINH
      WRITE(6,*) 'NECARD : MODE   = ',NEMODFLG
      WRITE(6,*) 'NECARD : COHEPI = ',NECOHEPI
      WRITE(6,*) 'NECARD : PDF    = ',NEPDF
      WRITE(6,*) 'NECARD : BODEK  = ',NEBODEK
      WRITE(6,*) 'NECARD : MDLQE  = ',MDLQE
      WRITE(6,*) 'NECARD : MDLQEAF= ',MDLQEAF
      WRITE(6,*) 'NECARD : MAQE   = ',XMAQE
C      WRITE(6,*) 'NECARD : MASPI  = ',XMASPI
      WRITE(6,*) 'NECARD : MARSRES = ',XMARSRES
      WRITE(6,*) 'NECARD : MANFFRES= ',XMANFFRES

      WRITE(6,*) 'NECARD : MVQE   = ',XMVQE
C      WRITE(6,*) 'NECARD : MVSPI  = ',XMVSPI
      WRITE(6,*) 'NECARD : MVRSRES = ',XMVRSRES
      WRITE(6,*) 'NECARD : MVNFFRES= ',XMVNFFRES
      WRITE(6,*) 'NECARD : KAPP   = ',KAPP
      WRITE(6,*) 'NECARD : XMACOH   = ',XMACOH
      WRITE(6,*) 'NECARD : R0COH   = ',RAD0NU
      WRITE(6,*) 'NECARD : A1COHerr   = ',fA1COH
      WRITE(6,*) 'NECARD : b1COHerr   = ',fb1COH

      WRITE(6,*) 'NECARD : SCCFV     = ',SCCFV
      WRITE(6,*) 'NECARD : SCCFA     = ',SCCFA
      WRITE(6,*) 'NECARD : FPQE      = ',FPQE
      WRITE(6,*) 'NECARD : PFSF      = ',PFSF
 
      IF (NEMODFLG.EQ.-1) THEN
         WRITE(6,*) 'NECARD : CRS  = ',CRSNEUT
         WRITE(6,*) 'NECARD : CRSB = ',CRSNEUTB
      ELSE
         DO 10 I=1,28
            CRSNEUT(I) = 1.0
            CRSNEUTB(I)= 1.0
 10      CONTINUE
      ENDIF

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
      WRITE(6,*) 'NECARD : IFORMLEN = ',IFORMLEN
      WRITE(6,*) 'NECARD : FZMU2 = ',FZMU2      
      WRITE(6,*) 'NECARD : ABSPIEMIT = ',NEABSPIEMIT
      WRITE(6,*) 'NECARD : NUSIM = ',NUSIM
      WRITE(6,*) 'NECARD : QUIET = ',QUIET

      WRITE(6,*) 'NECARD : FEFQE   = ',FEFQE
      WRITE(6,*) 'NECARD : FEFQEH  = ',FEFQEH
      WRITE(6,*) 'NECARD : FEFINEL = ',FEFINEL
      WRITE(6,*) 'NECARD : FEFABS  = ',FEFABS
      WRITE(6,*) 'NECARD : FEFCOH  = ',FEFCOH
      WRITE(6,*) 'NECARD : FEFCX   = ',FEFCX
      WRITE(6,*) 'NECARD : FEFCXH  = ',FEFCXH
      WRITE(6,*) 'NECARD : FEFCXHF = ',FEFCXHF
      WRITE(6,*) 'NECARD : FEFQEHF = ',FEFQEHF
      WRITE(6,*) 'NECARD : FEFCOHF = ',FEFCOHF
      WRITE(6,*) 'NECARD : FEFCOUL = ',FEFCOUL
      WRITE(6,*) 'NECARD : FEFALL = ',FEFALL

      WRITE(6,*) 'NECARD : IPILESSDCY = ',IPILESSDCY
      WRITE(6,*) 'NECARD : RPILESSDCY = ',RPILESSDCY

      WRITE(6,*) 'NECARD : IFF =',NEIFF
      WRITE(6,*) 'NECARD : NRTYPE =',NENRTYPE
      WRITE(6,*) 'NECARD : CA5 =',RNECA5I
      WRITE(6,*) 'NECARD : BGSCL =',RNEBGSCL

      WRITE(6,*) 'NECARD : IRADCORR =',IRADCORR

C     Parameter checks

      IQEMA = INT(XMAQE*100+0.5)
C      ISPIMA = INT(XMASPI*100+0.5)
      IRSRESMA  = INT(XMARSRES*100+0.5)
      INFFRESMA = INT(XMANFFRES*100+0.5)

      IF (IQEMA.ne.101 .and. IQEMA.ne.103 .and. IQEMA.ne.111 .and. 
     &    IQEMA.ne.121 .and. IQEMA.ne.131 .and. IQEMA.ne.160 ) THEN
         WRITE(6,*)
     &    'NECARD error: MAQE must be 1.01, 1.03, 1.[1-3]1 or 1.60'
         STOP
      END IF

C      IF ( ISPIMA.ne.95 .and. ISPIMA.ne.101 .and. ISPIMA.ne.103. 
C     $     .and. ISPIMA.ne.111 .and. ISPIMA.ne.121 ) THEN
C         WRITE(6,*) 
C     $        'NECARD error: MASPI must be 0.95, 1.01, 1.03, 1.11 or 1.21'
C         STOP
C      END IF


      IF ( IRSRESMA.ne.95 .and. IRSRESMA.ne.101 .and. IRSRESMA.ne.103. 
     $     .and. IRSRESMA.ne.111 .and. IRSRESMA.ne.121 ) THEN
         WRITE(6,*) 
     $     'NECARD error: MARS must be 0.95, 1.01, 1.03, 1.11 or 1.21'
         STOP
      END IF

      IF ( INFFRESMA.ne.95 .and. INFFRESMA.ne.101 .and. 
     $     INFFRESMA.ne.103. .and. INFFRESMA.ne.111 .and. 
     $     INFFRESMA.ne.121 ) THEN
         WRITE(6,*) 
     $    'NECARD error: MARNFF must be 0.95, 1.01, 1.03, 1.11 or 1.21'
         STOP
      END IF

      IF (FEFQEHF.lt.0 .or. FEFQEHF.gt.1) THEN
         WRITE(6,*) 'NECARD error: FEFQEHF must be between 0 and 1!'
         STOP
      END IF
      
      IF (FEFCXHF.lt.0 .or. FEFCXHF.gt.1) THEN
         WRITE(6,*) 'NECARD error: FEFCXHF must be between 0 and 1!'
         STOP
      END IF

C-- Add nuceff/nuccorspl related card
C--      call nesettarg
C--      call nrsettarg
C-- ADD proton nuceff related card
      call nrcard

      call nefillmodel

      RETURN
      END
