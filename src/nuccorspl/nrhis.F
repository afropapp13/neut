	  subroutine nrhis(x,rbin,num,ilo,ihi,erem)
	  implicit none

	  real*4    x,rbin,erem
	  integer*4 num,ilo,ihi

	  ilo=(x/rbin+0.01)
	  erem=(x-float(ilo)*rbin)/rbin
	  ilo=ilo+1
	  ihi=ilo+1
	  if(ilo.ge.num)then
		 ilo=num
		 ihi=num
		 erem=1.0000000
	  end if
	  if(ilo.lt.1)then
		 ilo=1
		 ihi=1
		 erem=1.0000000
	  end if
	  return
	  end
