C    ---------------------------------------------------------------------
      subroutine eftrcomega(xin,pin,kin,nflag)
C    ---------------------------------------------------------------------
C     [input]
C         xin(3)     : decay position in nucleus       [fm]
C         pin(3)     : initial omega momentum          [MeV/c^2]
C
C     [output]
C         common     : OMEGAOUT                            
C --------------------------------------------------------------
C  99/5/7 for proton decay MC of omega mode (Y.H)
C  03/4/4 modified for ATMnu MC (M.Ishitsuka)
C --------------------------------------------------------------
      implicit none
#include <efpion.h>

C -- function
      real effrmgas

C -- argument(input)      
      integer kin
      real xin(3),pin(3)

C -- argument(output)      
      integer fin_id(3),nump,log
      real xout(3),pfi(3,3)

C -- for tracking
      real x(3),p(3),absp,dir(3),step,distance,pferm,rmax
      real pf(3,3),pn(3)
      integer log_a,log_b,log_c,log_d,type,fpid(3)
*      data pferm/260/   ! use function effrmgas insted of this
      real dum
      integer i,j,ii
      INTEGER IPP,IP0,NL,NH
      INTEGER NFLAG
      INTEGER NOMEGATL, LOOPOMEGA
      INTEGER IPINTOMEGA(50), IPIORIOMEGA(50)
      REAL PINTOMEGA(3,50), XINTOMEGA(3,50)
      COMMON /OMEGAOUT/ NOMEGATL,IPINTOMEGA,PINTOMEGA
     &     ,XINTOMEGA,IPIORIOMEGA,LOOPOMEGA

C ---------------- start program
C ===============
C  initialize
C =============== 

      ii=1  ! pauli blocking counter
      type=0
      log_a=0
      log_b=0
      log_c=0
      log_d=0

C ===============
C  program start
C =============== 
      do i=1,3
         x(i)=xin(i)
         p(i)=pin(i)
      enddo

C -- determine the tracing limit from radius of nucleus
      rmax= cc2                    ! -- 2.5 * (radius)
*      step= rmax/100.              ! -- consider max 100 steps.
      step= rmax/1000.             ! -- consider max 1000 steps.
      
C ========================================================= 
C              start point of loop for tracing 
C ========================================================= 
 100  continue
      absp=sqrt(p(1)**2+p(2)**2+p(3)**2)
      if(absp .le. 0.001)then
         absp= 0.001
         dir(1)= 0.001
         dir(2)= 0.0
         dir(3)= 0.0
      else
         dir(1)=p(1)/absp
         dir(2)=p(2)/absp
         dir(3)=p(3)/absp
      endif
      
      x(1)=x(1)+dir(1)*step
      x(2)=x(2)+dir(2)*step
      x(3)=x(3)+dir(3)*step
      distance=sqrt(x(1)**2+x(2)**2+x(3)**2)

      if(distance .gt. rmax) then        ! escape from nucleus
         distance=rmax+0.001
c         do i=1,3,1
c            x(i) = -9999.
c         enddo
         log_b=3
         type=7
         goto 900                        ! and decay
      endif

C==========================================================
C  subroutine efprobom determine what will happen.
C----------------------------------------------------------
C event
C      1: omN - piN
C      2: omN - rhoN(pipiN)
C      3: omN - rhopiN(pipipiN)
C      4: omN - pipiN
C      5: omN - omN
C      6: omN - sigmaN(pipiN)
C      7: decay
C      8: nothing
C==========================================================
      call efprobom(absp,distance,step,type)
      if(type.eq.8)  goto 100

C===================================================================
C  determine kinematics for each process
C===================================================================
C efkinom
C   [input]
C      type: kinematics type (see above)
C      p(3):  momentum vector of omega
C   [output]
C      type: decay branch (only for decay)
C         11 pi+pi-pi0
C         12 pi0gamma
C         13 pi+pi-
C      pf(3,3):   final state momentum(pi,.. etc. not nucleon)
C      pn(3)  :   final state nucleon momentum(to check pauli blocking)
C      fpid(3):   final state particle ID(as PDG code)
C      nump: number of particles at final state
C-------------------------------------------------------------------
 900  continue
      call efkinom(type,p,pf,fpid,pn,nump)
      
*      print *,'efkinom: ',type,fpid,pf

C processing interaction 

*  get fermi surface momentum at given radii.      
      pferm=effrmgas(dum,dum,distance)

*  consider pauli blocking  -----------------

      if((type.ge.1).and.(type.le.6).and.
     &     (sqrt(pn(1)**2+pn(1)**2+pn(1)**2).lt.pferm)) then 

         ii=ii+1
         type=0

         goto 100

*  log varible ------------------------------
      else if((type.ge.1).and.(type.le.6)) then
         log_a=log_a+1
         log_d=log_d*10+type

         if(type.ne.5)then      ! finish at this interaction
            log_b=2
            log_c=4
         else                   ! need futher process
            type=0
            do i=1,3
               p(i)=pf(i,1)
            enddo
            goto 100
         endif
      endif
      
C processing decay
      if(type.gt.10)then
         log_c=type-10
         if(log_b.eq.0) log_b=1
      endif

C== finish (prepare output argument)
      log=log_d*10000+log_c*1000+log_b*100+log_a*10

      do j=1,nump
         fin_id(j)=fpid(j)
         do i=1,3
            pfi(i,j)=pf(i,j)
         enddo
      enddo

      do i=1,3
         xout(i)=x(i)
      enddo

      IF(NOMEGATL.GE.50)NOMEGATL=49
      NL=NOMEGATL+1
      NH=NOMEGATL+nump
      IF(NH.GT.50) THEN
         WRITE(6,600)
 600     FORMAT('  ********** ERROR  MULTIPLICITY > 50 IN OMEGA-DECAY ******')
         NH=50
      END IF
C     
      DO IPP=NL,NH
         IP0=IPP-NOMEGATL
         IPIORIOMEGA(IPP)=LOOPOMEGA
         PINTOMEGA(1,IPP)=pfi(1,IP0)
         PINTOMEGA(2,IPP)=pfi(2,IP0)
         PINTOMEGA(3,IPP)=pfi(3,IP0)
         XINTOMEGA(1,IPP)=xout(1)
         XINTOMEGA(2,IPP)=xout(2)
         XINTOMEGA(3,IPP)=xout(3)
         IPINTOMEGA(IPP)=fin_id(IP0)
      ENDDO
      NOMEGATL=NH
      NFLAG=1

      return 
      end

