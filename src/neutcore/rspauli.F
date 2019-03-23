*************************************************************
*     -------------------------------------------------------
      SUBROUTINE RSPAULI(PINUC,PIMES,MNUC,MMES,PONUC,POMES)
*     -------------------------------------------------------
*
*    (Purpose)
*     Recaluculete Momentum of Hadron System
*       In such Reaction as
*         Resonance state  -> 2body decay
*
*     Decay is made isotropically at Hadron CMS 
* 
*    (Input)
*     PINUC  Input Momentum of Nucleon ( Lab frame , GeV/c)
*     PIMES  Input Momentum of Meson   ( Lab frame , GeV/c)
*     MNUC   Nucleon mass (GeV)
*     MMES   Meson mass   (GeV)
*  
*    (Output)
*     PONUC  Output Momentum of Nucleon ( Lab frame , GeV/c)
*     POMES  Output Momentum of Meson   ( Lab frame , GeV/c)
*
*     (Creation Date & Author )
*     1998.09.??  ; J.Kameda 
*
*************************************************************

       IMPLICIT NONE

       REAL PINUC(3),PIMES(3),PONUC(3),POMES(3)
       REAL MNUC,MMES
       
       REAL VABS,V(3),GAMMA,UV(3)
       REAL ERES,ENUC,EMES
       REAL PABS,PRES(3),MRES
       REAL PCMS,DIR(3)


       PRES(1) = PINUC(1) + PIMES(1)
       PRES(2) = PINUC(2) + PIMES(2)
       PRES(3) = PINUC(3) + PIMES(3)

       ENUC = SQRT(PINUC(1)**2+PINUC(2)**2+PINUC(3)**2+MNUC**2)
       EMES = SQRT(PIMES(1)**2+PIMES(2)**2+PIMES(3)**2+MMES**2)

       PABS = SQRT(PRES(1)**2+PRES(2)**2+PRES(3)**2)
       ERES = ENUC + EMES

CC --------  ReCalculate Momentum in Hadron CMS --------

       MRES = ERES**2 - PABS**2
         IF (MRES.lt.0) then
          WRITE(*,*) 'RSPAULI: ERROR,IRREGAL INPUT MRES<0'
          WRITE(*,*) ERES,PABS
          RETURN
         ENDIF
       MRES = SQRT(MRES)
       PCMS = ((MRES**2-(MNUC**2+MMES**2))**2 - (2*MNUC*MMES)**2)
     $     /4/MRES**2
         IF (PCMS.lt.0) then
          WRITE(*,*) 'RSPAULI: ERROR,IRREGAL INPUT PCMS<0'
          WRITE(*,*) 'MRES',MRES,'MNUC',MNUC,'MMES',MMES
          RETURN
         ENDIF

       PCMS = SQRT(PCMS)

       CALL RNDIR(DIR)

       PONUC(1) = PCMS*DIR(1)
       PONUC(2) = PCMS*DIR(2)
       PONUC(3) = PCMS*DIR(3)
       POMES(1) = -1.0*PCMS*DIR(1)
       POMES(2) = -1.0*PCMS*DIR(2)
       POMES(3) = -1.0*PCMS*DIR(3)

CC ---------  Lorentz Boost CMS -> Lab system


       V(1) = PRES(1)/ERES
       V(2) = PRES(2)/ERES
       V(3) = PRES(3)/ERES
       VABS = PABS/ERES
       GAMMA = 1.0/SQRT(1.0 - VABS**2)

      IF (VABS.gt.0.00001) then
       UV(1) = V(1)/VABS
       UV(2) = V(2)/VABS
       UV(3) = V(3)/VABS
      ELSE
       UV(1) = 1.0
       UV(2) = 0.0
       UV(3) = 0.0
      ENDIF

      CALL MCVECBST(PONUC,MNUC,UV,GAMMA)
      CALL MCVECBST(POMES,MMES,UV,GAMMA)

      RETURN

      END
      
