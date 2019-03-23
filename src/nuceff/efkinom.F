* ==========================================================
      subroutine efkinom(type,p,pf,fpid,pn,nump)
* ----------------------------------------------------------
*
* calcurate kinematics for each process
*
*  [input]
*     type: kinematics type
*        1: omN - piN
*        2: omN - rhoN(pipiN)
*        3: omN - rhopiN(pipipiN)
*        4: omN - pipiN
*        5: omN - omN
*        6: omN - sigmaN(pipiN)
*        7: decay
*     p(3): momentum vector of omega
*
*  [output]
*     type: decay branch(only for decay)
*        11: pi+pi-pi0
*        12: pi0gamma
*        13: pi+pi-
*     pf(3,3): final momentum(pi+,pi-,....,not for nucleon)
*     pn(3): final state nucleon momentum
*     fpid(3): final state particle ID(as PDG code)
*     nump: number of particle at final state
*
* ==========================================================
      implicit none

* --- function
      real rlu,dummy,random
      
* ---- argument (input)
      integer type
      real p(3)

* ---- argument (output)
      integer nump,fpid(3)
      real pn(3),pf(3,3)

* ---- misc
      integer i,j
      real momferm,nucdir(3),dum,tmpdir(3)
      real invmass,pcms(3),pni(3)

      real rate_otoh
      data rate_otoh/0.5/
      
      real nucmass,ommass,pimass,pizmass,rhomass,sigmass
      data nucmass,ommass,pimass,pizmass,rhomass,sigmass
     &     /938.,782.,140.,135.,770.,780./
      
* ---- for N-body kinematics(vcphsvct)
      integer nbody,ierr
      real winv,amass(10),pmom(3,10),beta(3),absbeta,tote,gamma
      real ubeta(3),fngamfct

* ---- for omega decay
      integer imode,nbodyom(3),id_decay(10)
      data nbodyom/3,2,2/

      real decaymom(3,10),bratio(3)
      data bratio/0.89,0.09,0.02/

* -- decayid(i,j) i: branch,j: particle index
      integer decayid(3,3)
      data decayid/211,-211,111,
     &             111,22,0,
     &             211,-211,0/

*========== program start
      momferm=0
      do i=1,3
         fpid(i)=0
         nucdir(i)=0
         pni(i)=0
         pn(i)=0
         do j=1,3
            pf(i,j)=0
         enddo
      enddo

* ----------------------------------------
* get fermi momentum of nucleon 
* and get invariant mass & cms momentum
* ----------------------------------------
      if(type.lt.7) then
         call rnferm(momferm,dum)
         call rndir(nucdir)
         tote=sqrt(nucmass**2+momferm**2)+
     &        sqrt(ommass**2+p(1)**2+p(2)**2+p(3)**2)
         do i=1,3
            pni(i)=momferm*nucdir(i)
            beta(i)=(p(i)+pni(i))/tote
         enddo
         absbeta=sqrt(beta(1)**2+beta(2)**2+beta(3)**2)
         gamma=fngamfct(absbeta)
            
         do i=1,3
            ubeta(i)=beta(i)/absbeta
         enddo
         call ainvms(nucmass,ommass,pni,p,invmass,pcms)
      endif

      if(type.eq.1)goto 11 ! omN --> piN
      if(type.eq.2)goto 12 ! omN --> rhoN
      if(type.eq.3)goto 13 ! omN --> rhopiN
      if(type.eq.4)goto 14 ! omN --> pipiN
      if(type.eq.5)goto 15 ! omN --> omN
      if(type.eq.6)goto 16 ! omN --> sigmaN
      if(type.eq.7)goto 17 ! om --> decay
      
      nump=0
      type=-10
      
      goto 999

* *****  omN --> piN
 11   continue
      nump=1

      winv=invmass
      nbody=2
      amass(1)=pimass
      amass(2)=nucmass
      call vcphsvct(winv,amass,nbody,pmom,ierr)

      if(rlu(dummy).le.rate_otoh)then
         fpid(1)=211
      else
         fpid(1)=-211
      endif
      
      do i=1,3
         pf(i,1)=pmom(i,1)
         pn(i)=pmom(i,2)
      enddo
      call mcvecbst(pf(1,1),pimass,ubeta,gamma)
      call mcvecbst(pn,nucmass,ubeta,gamma)
      
      goto 999
      
* *****  omN --> rhoN
 12   continue

      nump=2

*----- interaction part(omN --> rhoN)
      winv=invmass
      nbody=2
      amass(1)=rhomass
      amass(2)=nucmass
      call vcphsvct(winv,amass,nbody,pmom,ierr)
      do i=1,3
         pf(i,1)=pmom(i,1)
         pn(i)=pmom(i,2)
      enddo
      call mcvecbst(pf(1,1),rhomass,ubeta,gamma)
      call mcvecbst(pn,nucmass,ubeta,gamma)

*----- rho decay part(rho --> pipi)
      do i=1,3
         beta(i)=pf(i,1)/
     &        sqrt(rhomass**2+pf(1,1)**2+pf(2,1)**2+pf(3,1)**2)
      enddo
      absbeta=sqrt(beta(1)**2+beta(2)**2+beta(3)**2)
      gamma=fngamfct(absbeta)
      do i=1,3
         ubeta(i)=beta(i)/absbeta
      enddo
      nbody=2
      amass(1)=pimass
      amass(2)=pizmass
      call vcphsvct(rhomass,amass,nbody,pmom,ierr)
      do i=1,3
         pf(i,1)=pmom(i,1)
         pf(i,2)=pmom(i,2)
      enddo
      call mcvecbst(pf(1,1),pimass,ubeta,gamma)
      call mcvecbst(pf(1,2),pizmass,ubeta,gamma)

      if(rlu(dummy).le.rate_otoh)then
         fpid(1)=211
         fpid(2)=111
      else
         fpid(1)=-211
         fpid(2)=111
      endif

      goto 999

* *****  omN --> rhopiN  %%%%%%%%%%%%
 13   continue

      nump=3

*----- interaction part(omN --> rhopiN) 
      winv=invmass
      nbody=3
      amass(1)=rhomass
      amass(2)=pizmass
      amass(3)=nucmass
      call vcphsvct(winv,amass,nbody,pmom,ierr)
      do i=1,3
         pf(i,1)=pmom(i,1)
         pf(i,2)=pmom(i,2)
         pn(i)=pmom(i,3)
      enddo
      call mcvecbst(pf(1,1),rhomass,ubeta,gamma)
      call mcvecbst(pf(1,2),pimass,ubeta,gamma)
      call mcvecbst(pn,nucmass,ubeta,gamma)


      if(rlu(dummy).le.rate_otoh)then
         fpid(1)=211
         fpid(2)=111
         fpid(3)=111
      else
         fpid(1)=-211
         fpid(2)=111
         fpid(3)=111
      endif

*----- rho decay part(rho --> pipi)
      do i=1,3
         beta(i)=pf(i,1)/
     &        sqrt(rhomass**2+pf(1,1)**2+pf(2,1)**2+pf(3,1)**2)
      enddo
      absbeta=sqrt(beta(1)**2+beta(2)**2+beta(3)**2)
      gamma=fngamfct(absbeta)
      do i=1,3
         ubeta(i)=beta(i)/absbeta
      enddo
      nbody=2
      amass(1)=pimass
      amass(2)=pizmass
      call vcphsvct(rhomass,amass,nbody,pmom,ierr)
      call rndir(tmpdir)
      do i=1,3
         pf(i,3)=pf(i,2)
         pf(i,1)=pmom(i,1)
         pf(i,2)=pmom(i,2)
      enddo
      call mcvecbst(pf(1,1),pimass,ubeta,gamma)
      call mcvecbst(pf(1,2),pizmass,ubeta,gamma)

      if(rlu(dummy).le.rate_otoh)then
         fpid(1)=211
         fpid(2)=111
         fpid(3)=111
      else
         fpid(1)=-211
         fpid(2)=111
         fpid(3)=111
      endif

      goto 999

* *****  omN --> pipiN
 14   continue

      nump=2

      winv=invmass
      nbody=3
      amass(1)=pimass
      amass(2)=pimass
      amass(3)=nucmass
      call vcphsvct(winv,amass,nbody,pmom,ierr)

      do i=1,3
         pf(i,1)=pmom(i,1)
         pf(i,2)=pmom(i,2)
         pn(i)=pmom(i,3)
      enddo
      call mcvecbst(pf(1,1),pimass,ubeta,gamma)
      call mcvecbst(pf(1,2),pimass,ubeta,gamma)
      call mcvecbst(pn,nucmass,ubeta,gamma)
      
      fpid(1)=211
      fpid(2)=-211

      goto 999

* *****  omN --> omN
 15   continue

      nump=1

      winv=invmass
      nbody=2
      amass(1)=ommass
      amass(2)=nucmass

      call rndir(tmpdir)

      do i=1,3
         pf(i,1)=sqrt(p(1)**2+p(2)**2+p(3)**2)*tmpdir(i)
         pn(i)=p(i)+pni(i)-pf(i,1)
      enddo
      
      fpid(1)=223

      goto 999

* *****  omN --> sigmaN
 16   continue

      nump=2
*----- interaction part(omN --> sigmaN)
      winv=invmass
      nbody=2
      amass(1)=sigmass
      amass(2)=nucmass
      call vcphsvct(winv,amass,nbody,pmom,ierr)
      do i=1,3
         pf(i,1)=pmom(i,1)
         pn(i)=pmom(i,2)
      enddo
      call mcvecbst(pf(1,1),sigmass,ubeta,gamma)
      call mcvecbst(pn,nucmass,ubeta,gamma)

*----- rho decay part(sigma --> pipi)
      do i=1,3
         beta(i)=pf(i,1)/
     &        sqrt(sigmass**2+pf(1,1)**2+pf(2,1)**2+pf(3,1)**2)
      enddo
      absbeta=sqrt(beta(1)**2+beta(2)**2+beta(3)**2)
      gamma=fngamfct(absbeta)
      do i=1,3
         ubeta(i)=beta(i)/absbeta
      enddo
      nbody=2
      amass(1)=pimass
      amass(2)=pimass
      call vcphsvct(sigmass,amass,nbody,pmom,ierr)
      do i=1,3
         pf(i,1)=pmom(i,1)
         pf(i,2)=pmom(i,2)
      enddo
      call mcvecbst(pf(1,1),pimass,ubeta,gamma)
      call mcvecbst(pf(1,2),pimass,ubeta,gamma)

      fpid(1)=-211
      fpid(2)=211

      goto 999

* *****  om --> decay
 17   continue

* deribe lorentz factor 
      do i=1,3
         beta(i)=p(i)/
     &        sqrt(ommass**2+p(1)**2+p(2)**2+p(3)**2)
      enddo
      absbeta=sqrt(beta(1)**2+beta(2)**2+beta(3)**2)
      gamma=fngamfct(absbeta)
      do i=1,3
         ubeta(i)=beta(i)/absbeta
      enddo

* determine decay mode
      random = rlu(dummy)
      if(random .lt. bratio(1)) then
         imode=1
      else if(random .lt. bratio(1)+bratio(2)) then
         imode=2
      else 
         imode=3
      endif

* calcurate kinematics
      do i=1,nbodyom(imode)
         id_decay(i)=decayid(i,imode)
         call mcmass(id_decay(i),amass(i))
      enddo
      
      call vcphsvct(ommass,amass,nbodyom(imode),decaymom,ierr)

      do i=1,nbodyom(imode)
         call mcvecbst(decaymom(1,i),amass(i),ubeta,gamma)
         fpid(i)=decayid(i,imode)
         do j=1,3
            pf(j,i)=decaymom(j,i)
         enddo
      enddo    

      nump=nbodyom(imode)      
      type=10+imode

*--------------------------
 999  continue
      return
      end

