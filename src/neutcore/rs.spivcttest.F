c
c
c     $Log: rs.spivcttest.F,v $
c     Revision 1.1.1.1  2016/11/04 20:31:26  mdunkman
c     Import of NEUT 5.3.3 v1r27 for use with T2KReWeight v1r27 (Winter OA 2017)
c
c     Revision 1.2  1998/08/12 14:18:05  cmauger
c     Just change the name of the program from test to rsspiT.  Also, add
c     an rcs log.
c
c
c

      program rsspiT

#include<nework.h>
      real dir(3)
      real H
      COMMON/PAWC/H(100000)

      COMMON/FREPRO/IFREE
      IFREE=0

      CALL HLIMIT(100000)

      CALL HBOOK1(10,"Energy of incoming neutrino" ,500,0.,5.,0.)
      CALL HBOOK1(11,"Momentum of incident nucleon",100,0.,1.,0.)

      Call HBOOK1(20,"Momentum of outgoing lepton" ,1000,0.,10.,0.)
      Call HBOOK1(21,"Momentum of outgoing pion  " ,1000,0.,10.,0.)
      Call Hbook1(22,"Momentum of outgoing nucleon",1000,0.,10.,0.)

      Call HBOOK1(30,"Angle between neutrino and outgo-lepton",
     $            314,0.,3.14,0.)
      Call HBOOK1(31,"Angle between neutrino and outgo-pion",
     $            314,0.,3.14,0.)
      Call HBOOK1(32,"Angle between neutrino and outgo-nucleon",
     $            314,0.,3.14,0.)
      Call HBOOK1(33,"Angle between outgo-lepton and outgo-pion",
     $            314,0.,3.14,0.)
      
      Call HBOOK1(40,"Cos between neutrino and outgo-lepton",
     $            200,-1.,1.,0.)
      Call HBOOK1(41,"Cos between neutrino and outgo-pion",
     $            200,-1.,1.,0.)
      Call HBOOK1(42,"Cos between neutrino and outgo-nucleon",
     $            200,-1.,1.,0.)
      Call HBOOK1(43,"Cos between outgo-lepton and outgo-pion",
     $            200,-1.,1.,0.)

      Call HBOOK1(50,"Invariant Mass(incoming)",
     $            1000,0.,10.,0.)
      Call HBOOK1(51,"Invariant Mass(pion-nucleon)",
     $            1000,0.,10.,0.)
      Call HBOOK1(52,"Invariant Mass(outgoing)",
     $            1000,0.,10.,0.)
      Call HBOOK1(53,"Invariant Mass difference(incoming-outgoing)",
     $            1000,-.5,+.5,0.)

      JMOD = 12
      IPAR  = 14
      EMAX  = 4.0
      
      do 20 I=1,300
         write(*,*) "Now processing Event #",I 

 10      dir(1) = rand(0)
         dir(2) = rand(0)
         if ((dir(1)**2 + dir(2)**2).GT.1.) goto 10
         if ((dir(1)**2 + dir(2)**2).eq.1.) then
            dir(3) = 0.
         else
            dir(3) = sqrt(1. - (dir(1)**2 + dir(2)**2))
         endif
         
         EIN = EMAX
         CALL HF1(10,EIN,1.)

         
         call rsspivct(IPAR,JMOD,EIN,dir,IERR)
         
         INUM = 2
         XMASS = 939./1000.
         CALL GENP(XMASS,PNE(1,INUM),PPR)
         CALL HF1(11,PPR,1.)

         DO 15 INUM = 3,5
            IF (IPNE(INUM).eq.IPAR) THEN
               XMASS=0
               CALL GENP(XMASS,PNE(1,INUM),PLP)
               CALL HF1(20,PLP,1.)
               GOTO 15
            endif
            IF (ABS(IPNE(INUM)).eq.11) THEN
               XMASS=0.511/1000.
               CALL GENP(XMASS,PNE(1,INUM),PLP)
               CALL HF1(20,PLP,1.)
               GOTO 15
            endif
            IF (ABS(IPNE(INUM)).eq.13) THEN
               XMASS=105.6/1000.
               CALL GENP(XMASS,PNE(1,INUM),PLP)
               CALL HF1(20,PLP,1.)
               GOTO 15
            endif
            IF (ABS(IPNE(INUM)).eq.15) THEN
               XMASS=1.78
               CALL GENP(XMASS,PNE(1,INUM),PLP)
               CALL HF1(20,PLP,1.)
               GOTO 15
            endif

            if (ABS(IPNE(INUM)).eq.211) THEN
               XMASS = 139./1000.
               CALL GENP(XMASS,PNE(1,INUM),PPI)
               CALL HF1(21,PPI,1.)
               GOTO 15
            ENDIF
            if (IPNE(INUM).eq.111) THEN
               XMASS=135./1000.
               CALL GENP(XMASS,PNE(1,INUM),PPI)
               CALL HF1(21,PPI,1.)
               GOTO 15
            ENDIF
            IF ((IPNE(INUM).eq.2112).or.(IPNE(INUM).eq.2212)) THEN
               XMASS=939./1000.
               CALL GENP(XMASS,PNE(1,INUM),PPR)
               CALL HF1(22,PPR,1.)
               GOTO 15
            ENDIF
            write(*,*) "What??"
 15      continue

         cosine=dir(1)*PNE(1,3)+dir(2)*PNE(2,3)+dir(3)*PNE(3,3)
         cosine=cosine/sqrt(PNE(1,3)**2+PNE(2,3)**2+PNE(3,3)**2)
         CALL HF1(40,cosine,1.)
         CALL HF1(30,acos(cosine),1.)

         cosine=dir(1)*PNE(1,4)+dir(2)*PNE(2,4)+dir(3)*PNE(3,4)
         cosine=cosine/sqrt(PNE(1,4)**2+PNE(2,4)**2+PNE(3,4)**2)
         CALL HF1(42,cosine,1.)
         CALL HF1(32,acos(cosine),1.)

         cosine=dir(1)*PNE(1,5)+dir(2)*PNE(2,5)+dir(3)*PNE(3,5)
         cosine=cosine/sqrt(PNE(1,5)**2+PNE(2,5)**2+PNE(3,5)**2)
         CALL HF1(41,cosine,1.)
         CALL HF1(31,acos(cosine),1.)

         cosine= PNE(1,3)*PNE(1,5)+PNE(2,3)*PNE(2,5)
     $          +PNE(3,3)*PNE(3,5)
         cosine=cosine/sqrt(PNE(1,5)**2+PNE(2,5)**2+PNE(3,5)**2)
     $                /sqrt(PNE(1,3)**2+PNE(2,3)**2+PNE(3,3)**2)
         CALL HF1(43,cosine,1.)
         CALL HF1(33,acos(cosine),1.)


         XMASS = 0.
         CALL GENE(XMASS,PNE(1,1),ENE1)
         XMASS = 939./1000.
         CALL GENE(XMASS,PNE(1,2),ENE2)
         W2=  (ENE1+ENE2)**2
     $      - (PNE(1,1)+PNE(1,2))**2-(PNE(2,1)+PNE(2,2))**2
     $      - (PNE(3,1)+PNE(3,2))**2
         WINC =  sqrt(W2)
         CALL HF1(50,WINC,1.)
        
         XMASS = 939./1000.
         CALL GENE(XMASS,PNE(1,4),ENENC)
         if (ABS(IPNE(5)).eq.211) THEN
            XMASS = 139./1000.
            CALL GENE(XMASS,PNE(1,5),ENEPI)
         ENDIF
         if (IPNE(5).eq.111) THEN
            XMASS=135./1000.
            CALL GENE(XMASS,PNE(1,5),ENEPI)
         ENDIF
         W2=  (ENENC+ENEPI)**2 
     $      - (PNE(1,5)+PNE(1,4))**2-(PNE(2,5)+PNE(2,4))**2
     $      - (PNE(3,5)+PNE(3,4))**2
         W = sqrt(W2)
         CALL HF1(51,W,1.)


         IF (IPNE(3).eq.IPAR) THEN
            XMASS=0
            CALL GENE(XMASS,PNE(1,3),ENELEP)
         endif
         IF (ABS(IPNE(3)).eq.11) THEN
            XMASS=0.511/1000.
            CALL GENE(XMASS,PNE(1,3),ENELEP)
         endif
         IF (ABS(IPNE(3)).eq.13) THEN
            XMASS=105.6/1000.
            CALL GENE(XMASS,PNE(1,3),ENELEP)
         endif
         IF (ABS(IPNE(3)).eq.15) THEN
            XMASS=1.78/1000.
            CALL GENE(XMASS,PNE(1,3),ENELEP)
         endif

         W2=  (ENELEP+ENENC+ENEPI)**2 
     $      - (PNE(1,5)+PNE(1,4)+PNE(1,3))**2
     $      - (PNE(2,5)+PNE(2,4)+PNE(2,3))**2
     $      - (PNE(3,5)+PNE(3,4)+PNE(3,3))**2
         WOUT = sqrt(W2)
         CALL HF1(52,WOUT,1.)
         CALL HF1(53,WINC-WOUT,1.)


 20   continue

      CALL HRPUT(0,'single-pi.hbk','NT')

      STOP
      END
      

      
      SUBROUTINE GENP(XMASS,XMOM,ABSP)

      REAL XMASS,XMOM(3),ABSP

      ABSP = -9999.
      PTOT2 = 0.

      DO 100 I=1,3
         PTOT2 = PTOT2 + XMOM(I)**2
 100  CONTINUE
      ABSP = sqrt(PTOT2)

      return
      end
      

      SUBROUTINE GENE(XMASS,XMOM,ENERGY)

      REAL XMASS,XMOM(3),ENERGY

      ENERGY = -9999.
      PTOT2 = 0.

      DO 100 I=1,3
         PTOT2 = PTOT2 + XMOM(I)**2
 100  CONTINUE
      ENERGY = sqrt((XMASS**2)+PTOT2)

      return
      end
      

