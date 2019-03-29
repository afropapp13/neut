C==================================================================
      subroutine efprobom(absp,distance,step,type)
C==================================================================
C
C     determine what process will happen.
C
C     input:
C         absp - absolute momentum                          [MeV/c^2]
C         distance - distance from center of nucleus        [fm]
C         step - trace step                                 [fm]
C     output:
C         type - process ID which will happen hereafter.
C
C--------------------------------------------------------------------
C    Note:
C       this routine uses aproximation as 1-exp(-Ax)=Ax for each
C       process, so value "step" supposed to be very small.
C--------------------------------------------------------------------
      implicit none
C - function
      real efabrho,rlu,efxsctom

C - argument input
      real absp,distance,step

C - argument output
      integer type

C - misc variables
      integer i,j
      real prob(8)
      real dens ! nucleon density
      real dens_0 ! nucleon density at center of nucleus
      real omlife  ! omega lifetime at rest(fm)
      real ommass
      real random,dum,sumrand
      data ommass/782./
      data omlife/23.4/
      data dens_0/0.15/

*==============program start
      dens = dens_0*efabrho(distance)

* ------- interaction omN - xN probability
      do i=1,6
         prob(i)=efxsctom(absp,i)*dens*0.1*step
      enddo

* ------- decay probability
      prob(7)=ommass*step/(omlife*absp)

* %%%%%%% none probability
      prob(8)=1
      do i=1,7
         prob(8)=prob(8)-prob(i)
      enddo

*===== play a dise
      random=rlu(dum)

      do i=1,7
         sumrand=0
         do j=1,i
            sumrand=sumrand+prob(j)
         enddo
         if(random.le.sumrand) then
            type=i
            goto 999
         endif
      enddo

***** nothing
      type=8

 999  continue
***********************************************
      return
      end








