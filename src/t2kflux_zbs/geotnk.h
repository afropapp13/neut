************************************************************************
*     -----------------
*     INCLUDE 'geotnk.h'
*     -----------------
*
*     (Purpose)
*       Set basic parameters for the Super-Kamioka water tank.
*
*     (Comments)
*       Unit in cm for entire programs
*
*     (Creation Date and Author)
*       1992.01.16 ; First version by Y.Suzuki
*       1998.03.19 ; ZMED modified from 45 to 55 cm
*       2001.05.04 ; Intermedidate version for true 1kt geom by Y.Itow
*                    STEP 70.7->70cm, ID hight 862.4->854
*                    CENTINTK 54->44, still RINTK is for 70.7 spacing
* 
************************************************************************
*
*     (Modification)
*       Basic geometry parameters are now set by PARAMETER statements.
*       And variable names are changed so as to include TK at their
*       end.
*
*     (Date and Author)
*       1992.05.12 ; J. Kanzaki
*
************************************************************************
*
*     (Modification)
*       Parameters are set to the value as of OCT '93
*        New parameters for the size of the wall frame
*     
*       By Y. Suzuki (10/23/93)
*
************************************************************************

      REAL DITKTK, HITKTK, DIINTK, HIINTK, HIWAL, TTATTK, TBATTK, 
     &   TWATTK, RTKTK, ZPTKTK, ZMTKTK, RINTK, ZPINTK, ZMINTK, ZPWAL,
     &   ZMWAL,ZCNTTK, RMED, ZMED
#ifdef ICHI_KILO
      REAL CENTINTK
C     Set basic tank parameters

C     DITKTK=1080.   ; Diameter of the water tank
      PARAMETER (DITKTK=1080.)
C     HITKTK=1085.   ; Height of the water tank ( to water surface )
      PARAMETER (HITKTK=1080.)
C     DIINTK=855.2   ; Diameter of the inner volume
      PARAMETER (DIINTK=860.0)
C     HIINTK=862.4   ; Height of the inner volume
      PARAMETER (HIINTK=854.)
C     HIWAL=848.4   ; Height of the barrel pmt frames
      PARAMETER (HIWAL=840.)
C     TTATTK=55.    ; Thickness of the top anti layer
      PARAMETER (TTATTK=66.)
C     TBATTK=52.6    ; Thickness of the bottom anti layer
      PARAMETER (TBATTK=64.)
C     TWATTK=67.4    ; Thickness of the wall anti layer
      PARAMETER (TWATTK=67.4)

C     CENTINTK = 54. ; Shift of center of inner volume
      parameter (CENTINTK = 44.)      

#else

C     Set basic tank parameters

C     DITKTK=3930.   ; Diameter of the water tank
      PARAMETER (DITKTK=3930.)
C     HITKTK=4140.   ; Height of the water tank ( to water surface )
      PARAMETER (HITKTK=4140.)
C     DIINTK=3380.   ; Diameter of the inner volume
      PARAMETER (DIINTK=3380.)
C     HIINTK=3620.   ; Height of the inner volume
      PARAMETER (HIINTK=3620.)
C     HIWAL=3605.7   ; Height of the barrel pmt frames
      PARAMETER (HIWAL=3605.7)
C     TTATTK=260.    ; Thickness of the top anti layer
      PARAMETER (TTATTK=260.)
C     TBATTK=260.    ; Thickness of the bottom anti layer
      PARAMETER (TBATTK=260.)
C     TWATTK=275.    ; Thickness of the wall anti layer
      PARAMETER (TWATTK=275.)

#endif

C     RTKTK=DITKTK/2. ; Radius of the water tank

      PARAMETER (RTKTK = DITKTK/2. )
C     ZPTKTK=HITKTK/2.  ; Height of the tank from the center
      PARAMETER (ZPTKTK=HITKTK/2.)
C     ZMTKTK=-HITKTK/2. ; -Height of the tank from the center
      PARAMETER (ZMTKTK=-HITKTK/2.)
C     RINTK = DIINTK/2. ; Radius of the inner volume
      PARAMETER (RINTK=DIINTK/2.)
C     ZPINTK=HIINTK/2. ; Height of the inner volume from the center
      PARAMETER (ZPINTK=HIINTK/2.)
C     ZMINTK=-HIINTK/2. ; -Height of the inner volume from the center
      PARAMETER (ZMINTK=-HIINTK/2.)
C     ZPWAL=HIWAL/2.   ; Height of the wall from the center
      PARAMETER (ZPWAL=HIWAL/2.)
C     ZMWAL=-HIWAL/2.   ; -Height of the wall from the center
      PARAMETER (ZMWAL=-HIWAL/2.)
C     ZCNTTK=HITKTK/2. ; Height of zero point (center) from the bottom
c                      ; of the tank
      PARAMETER (ZCNTTK=HITKTK/2.)

C     Thickness of the nonsensitive regeion
#ifdef ICHI_KILO
      PARAMETER ( RMED  = 58. )
      PARAMETER ( ZMED  = 56. )
#else
      PARAMETER ( RMED  = 55. )
      PARAMETER ( ZMED  = 55. )
#endif


