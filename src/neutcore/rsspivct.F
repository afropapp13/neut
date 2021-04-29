***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE RSSPIVCT(IPAR,JMOD,E,DNEUT,IERR)
*     ---------------------------------------------------
*
*     ( purpose )
*       VECTOR GENERATION FOR SINGLE PION PRODUCTION
*
*     ( input )
*       IPAR        : kind of incoming neutrino
*       JMOD        : interaction mode
*       E           : Energy of incident neutrino
*       dneut(3)    : direction of incoming neutrino
*     
*     ( output )
*       IERR        : ERROR CODE
*
*     ( creation date and author )
*       1993.Dec. ; for Super-Kamioka by Y.Hayato
*                   (Based on Rein and Sehgal's code)
*		1995.Sep. ; BUG fix(around calc. sinth,costh,sinphi,cosphi)
*       1995.Nov. ; BUG fix(set IORGNE,IFLGNE)
*       1996.Nov. ; BUG fix(Avoid Infinite LOOP,Check more carefully)
*       1997.Dec. ; BUG fix(add common INO16)
*       1998.09   ; ADD Pauli blocking by J.Kameda
*       2003.03.17; Change the scheme of the pauli blocking
*                   If the Pauli blocking occurs, that interaction
*                   is suppressed by setting 
*                      IFLGNE to 5 
*                   and  
*                      ICRNE  to 0 (Y.Hayato)
*                   Also use PFSURF instead of 226(constant)          
*       2003.03.30; Use nefrmmom instead of rnferm (Y.Hayato)
*       2006.06.12; Use IBOUND instead of INO16 and IFREE. (Y.H)
*                 
*
*     ( comment )
*       Reference:Rein and Sehgal,Ann. of Phys. 133(1981),79-153
*                 Rein           ,Z.Phys.C 35(1987),43-64
*                 Feynman et al. ,Phys.Rev.D3(1971),2706-
***********************************************************************


      IMPLICIT NONE

     
      REAL DNEUT(3),PNEU(3),DIRN(3)

C      INTEGER MAXNE
#include "necard.h"
#include "nework.h"     
#include "neutparams.h"

C      INTEGER MODENE,NUMNE,IPNE,IORGNE,IFLGNE,ICRNNE
C      REAL    PNE

      REAL*4  TMPEP,TMPEL,TMPEEE
      REAL*4  TMPSIG
      INTEGER*4 IMOD
      REAL*4  RSSPICRS
      EXTERNAL RSSPICRS

      REAL*4 evdifcrs, check
      external evdifcrs

C      INTEGER IFREE,IPAR,JMOD,IERR,IEC,IDECAY,I,INO16,IWAVE
      INTEGER IPAR,JMOD,IERR,IEC,IDECAY,I,IWAVE
      INTEGER IFLAG
      REAL    E,PFABS,RTPLEP(3),RTPPI(3),RTPNC(3),DX,DY,DZ
      REAL    PNCEPZ(3)
      REAL    COSTH,SINTH,COSPH,SINPH
      REAL    PTMP,POMES(3),PONUC(3),POUT
      INTEGER NTMP
      REAL    PROB,W,PTOT(3),E1,E2
      REAL    DUM,RLU
      EXTERNAL RLU

#include<rscons.h>

C-ADD definition --
C      REAL*4 XIN(3)
C      COMMON /SAVXIN/INO16,XIN
C      COMMON/FREPRO/IFREE

#include "vcwork.h"
#include "posinnuc.h"

C
      IERR=0
      IEC=0
      IDECAY=0
C
C -- SET IMOD AND LVECT
C
      NUMNE=5

C-CONVERT MODE
      IF (JMOD.eq.11)  IMOD=1
      IF (JMOD.eq.12)  IMOD=2
      IF (JMOD.eq.13)  IMOD=3
      IF (JMOD.eq.31)  IMOD=11
      IF (JMOD.eq.32)  IMOD=12
      IF (JMOD.eq.33)  IMOD=13
      IF (JMOD.eq.34)  IMOD=14
      IF (JMOD.eq.-11) IMOD=4
      IF (JMOD.eq.-12) IMOD=5
      IF (JMOD.eq.-13) IMOD=6
      IF (JMOD.eq.-31) IMOD=15
      IF (JMOD.eq.-32) IMOD=16
      IF (JMOD.eq.-33) IMOD=17
      IF (JMOD.eq.-34) IMOD=18
C
C -- SET NEUTRINO DIRECTION
C 
      DO 100 I=1,3
         PNE(I,1)=E*DNEUT(I)
         PNEU(I)=PNE(I,1)
  100 CONTINUE
      IPNE(1)=IPAR

      PFABS=0.

 1500 CONTINUE

C    If we have a bounded nucleon get the momentum of it
      IF(IBOUND.EQ.1 .AND. ABS(IPAR).NE.16 .AND. PFMAX.NE.0) THEN
        CALL nefrmmom(PFABS,IWAVE)
      endif
      PFABS = PFABS/1000.

      CALL RNDIR(DIRN)
      PNE(1,2)=DIRN(1)*PFABS
      PNE(2,2)=DIRN(2)*PFABS
      PNE(3,2)=DIRN(3)*PFABS

C     The energy of the nucleon
      TMPEP=sqrt(XMN2+(PNE(1,2)**2+PNE(2,2)**2+PNE(3,2)**2))
C     The energy of the neutrino
      TMPEL=E

      TMPEEE=( (TMPEP+TMPEL)**2
     $        -( (PNE(1,2)+PNEU(1))**2+(PNE(2,2)+PNEU(2))**2
     $          +(PNE(3,2)+PNEU(3))**2)-XMN2)/2./XMN

C-Check IF sigma < 1.e-42 then re-calculate fermi-momentum 
C     Get's the x-sec from free-nucleon
      TMPSIG=RSSPICRS(TMPEEE,IPAR,IMOD)
      IF (TMPSIG.le.1.e-4) GO TO 1500

C--- initialize/set some commons
C     Here set the parents
      IORGNE(1)=0
      IORGNE(2)=0
C    lepton :  from initial lepton
      IORGNE(3)=1
C    pion and nucleon :from initial nucleon
      IORGNE(4)=2
      IORGNE(5)=2

C     Set flags of initial and final state particles
      IFLGNE(1)=-1
      IFLGNE(2)=-1
      IFLGNE(3)=0
      IFLGNE(4)=0
      IFLGNE(5)=0

      MODENE = JMOD

C -- SET KIND OF PARTICLE
C
      IF(IPAR.LT.0)GO TO 220
      IF(JMOD.EQ.11)GO TO 201
      IF(JMOD.EQ.12)GO TO 202
      IF(JMOD.EQ.13)GO TO 203
      IF(JMOD.EQ.31)GO TO 204
      IF(JMOD.EQ.32)GO TO 205
      IF(JMOD.EQ.33)GO TO 206
      IF(JMOD.EQ.34)GO TO 207
  220 IF(JMOD.EQ.-11)GO TO 208
      IF(JMOD.EQ.-12)GO TO 209
      IF(JMOD.EQ.-13)GO TO 210
      IF(JMOD.EQ.-31)GO TO 211
      IF(JMOD.EQ.-32)GO TO 212
      IF(JMOD.EQ.-33)GO TO 213
      IF(JMOD.EQ.-34)GO TO 214

      GO TO 999

C --- SET THE PDG CODES FOR OUTGOING PARTICLES ---
C====================
C     NEUTRINO
C====================
C     CC1pi+ on proton
  201 IPNE(2)=2212
      IPNE(3)=IPAR-1
      IPNE(4)=2212
      IPNE(5)=211
      IFLAG=1
      GO TO 300
C     CC1pi0 on neutron
  202 IPNE(2)=2112
      IPNE(3)=IPAR-1
      IPNE(4)=2212
      IPNE(5)=111
      IFLAG=2
      GO TO 300
C     CC1pi+ on neutron
  203 IPNE(2)=2112
      IPNE(3)=IPAR-1
      IPNE(4)=2112
      IPNE(5)=211
      IFLAG=3
      GO TO 300
C     NC1pi0 on neutron
  204 IPNE(2)=2112
      IPNE(3)=IPAR
      IPNE(4)=2112
      IPNE(5)=111
      IFLAG=6
      IFLGNE(3)=2
      GO TO 300
C     NC1pi0 on neutron
  205 IPNE(2)=2212
      IPNE(3)=IPAR
      IPNE(4)=2212
      IPNE(5)=111
      IFLAG=4
      IFLGNE(3)=2
      GO TO 300
C     NC1pi- on neutron
  206 IPNE(2)=2112
      IPNE(3)=IPAR
      IPNE(4)=2212
      IPNE(5)=-211
      IFLAG=7
      IFLGNE(3)=2
      GO TO 300
C     NC1pi+ on proton
  207 IPNE(2)=2212
      IPNE(3)=IPAR
      IPNE(4)=2112
      IPNE(5)=211
      IFLAG=5
      IFLGNE(3)=2
      GO TO 300

C====================
C     ANTINEUTRINO
C====================
C     CC1pi- on n
  208 IPNE(2)=2112
      IPNE(3)=IPAR+1
      IPNE(4)=2112
      IPNE(5)=-211
      IFLAG=11
      GO TO 300
C     CC1pi0 on p
  209 IPNE(2)=2212
      IPNE(3)=IPAR+1
      IPNE(4)=2112
      IPNE(5)=111
      IFLAG=12
      GO TO 300
C     CC1pi- on p
  210 IPNE(2)=2212
      IPNE(3)=IPAR+1
      IPNE(4)=2212
      IPNE(5)=-211
      IFLAG=13
      GO TO 300
C     NC1pi0 on n
  211 IPNE(2)=2112
      IPNE(3)=IPAR
      IPNE(4)=2112
      IPNE(5)=111
      IFLAG=16
      IFLGNE(3)=2
      GO TO 300
C     NC1pi0 on p
  212 IPNE(2)=2212
      IPNE(3)=IPAR
      IPNE(4)=2212
      IPNE(5)=111
      IFLAG=14
      IFLGNE(3)=2
      GO TO 300
C     NC1pi- on n
  213 IPNE(2)=2112
      IPNE(3)=IPAR
      IPNE(4)=2212
      IPNE(5)=-211
      IFLAG=17
      IFLGNE(3)=2
      GO TO 300
C     NC1pi+ on p
  214 IPNE(2)=2212
      IPNE(3)=IPAR
      IPNE(4)=2112
      IPNE(5)=211
      IFLAG=15
      IFLGNE(3)=2

 300  CONTINUE

      DO 310 I=1,3
         PNCEPZ(I)=PNE(I,2)
 310  CONTINUE

      IF (IPNE(3).ne.IPAR) ICRNNE(3) = 1
      ICRNNE(4) = 1
      ICRNNE(5) = 1

      DX=DNEUT(1)
      DY=DNEUT(2)
      DZ=DNEUT(3)

      COSTH=DZ
      if (abs(dz).ge.1.) then
         sinth = 0.
      else
         SINTH=SQRT(1.-DZ**2)
      endif
      IF (ABS(DX).GT.0.000001) THEN
        COSPH=COS(ATAN2(DY,DX))
        SINPH=SIN(ATAN2(DY,DX)) 
      ELSE
        COSPH=0.
        IF (ABS(DY).GT.0.00001) THEN
          SINPH=1.*(ABS(DY)/DY)
        ELSE
          SINPH=1.
        ENDIF
      ENDIF 
 
      CALL RSROTVEC(PNE(1,2),PNCEPZ(1),
     $     1.,0.,COSPH,-1.*SINPH)
      CALL RSROTVEC(PNCEPZ(1),PNCEPZ(1),
     $     COSTH,-1.*SINTH,1.,0.)

 1000 continue

      CALL RSDIST(IFLAG,IPNE(1),E,PNCEPZ(1),RTPLEP,RTPPI,RTPNC,IERR)
      IF (ierr.ne.0) then
         return
      endif

C     Do Pauli blocking on the outgoing nucleon
C     CWRET add consistenly turn off pauli blocking for res
      IF (IBOUND.eq.1 .and. nepauflg.eq.0) THEN
         PTMP = sqrt(RTPNC(1)**2+RTPNC(2)**2+RTPNC(3)**2)
         if (PTMP.LE.PFSURF) THEN
            DO 800 I=3,5
               IFLGNE(I)=5
               ICRNNE(I) =0
 800        continue
         endif
      endif

      CALL RSROTVEC(RTPLEP,PNE(1,3),COSTH,SINTH,COSPH,SINPH)
      CALL RSROTVEC(RTPNC ,PNE(1,4),COSTH,SINTH,COSPH,SINPH)
      CALL RSROTVEC(RTPPI ,PNE(1,5),COSTH,SINTH,COSPH,SINPH)

CCC     Consider Non-pion Decay of Delta resonance
      IF ((IBOUND.eq.1).and.(IPILESSDCY.ne.0).and.ICRNNE(4).eq.1) THEN
         prob = RLU(dum)
         if (prob.le.RPILESSDCY) then
            if (QUIET.eq.0) then
               write(*,*) 'RSSPIVCT.F : DELTA NON-PI DECAY in nuclear'
            endif
            
C-- modify to put correct PID and # of particles for delta absorption
            NUMNE     = 4
C   Delta++                
             if (jmod.eq.11)then
                IPNE(4)   = 2224
C   Delta+                
             else if ((jmod.eq. 12).or.(jmod.eq. 13).or.
     $                (jmod.eq. 32).or.(jmod.eq. 34).or.
     $                (jmod.eq.-32).or.(jmod.eq.-34)) then
                IPNE(4)   = 2214
C   Delta0   
             else if ((jmod.eq. 31).or.(jmod.eq. 33).or.
     $                (jmod.eq.-12).or.(jmod.eq.-13).or.
     $                (jmod.eq.-31).or.(jmod.eq.-33)) then
                IPNE(4)   = 2114
C   Delta-  
             else if (jmod.eq.-11) then
                IPNE(4)   = 1114
             endif
             
             IFLGNE(4) = 3 
             ICRNNE(4) = 0
             IPNE(5)   = 0
             IFLGNE(5) = 0 
             ICRNNE(5) = 0
          endif          
       ENDIF

C       DONE WITH THE PIONLESS DELTA DECAY

C     check=evdifcrs()
C     write(*,*) "evdifcrs after particles done:", check
      RETURN

 999  WRITE(*,*)
     $ "RSSPIVCT:*** Invalid JMOD(Interaction Code):JMOD=",JMOD," ***"
      STOP

      END
