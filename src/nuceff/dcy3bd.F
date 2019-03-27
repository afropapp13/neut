        subroutine dcy3bd(am0,am1,am2,am3,amomd1,amomd2,amomd3)
C
C   GET MOMENTUM OF DECAYED PARTICLE IN C.M.S. ONLY 3 BODY DECAY IS
C  CONSIDERED.
C      INPUT  : AM0 - MASS OF DECAYED PARTICLE
C               AM1 -    "    DECAY PRODUCT 1
C               AM2 -    "          "       2
C               AM3 -    "          "       3
C      OUTPUT : AMOMD1(3) - MOMENTUM OF DECAY PRODUCT 1
C               AMOMD2(3) -            "              2
C               AMOMD3(3) -            "              3
C
C    M.ISHITSUKA   02-AUG-02
C
        implicit none
        real am0,am1,am2,am3
        real amomd1(3),amomd2(3),amomd3(3)
        real p1,p2,p3             !total momenta of daughter 1
        real E1,E2,E3          !energy of daughters in COM frame
        real E1max,E2max,E3max
        real Etotal
        real cang           !COM frame angle of daughter1
        integer i
        real rlu

        E1max=(am0**2-(am2+am3)**2+am1**2)/(am0*2.)  !E1 max
        E2max=(am0**2-(am1+am3)**2+am2**2)/(am0*2.)  !E2 max
        E3max=(am0**2-(am1+am2)**2+am3**2)/(am0*2.)  !E2 max
5       E1=rlu()*(E1max-am1)+am1
        E2=rlu()*(E2max-am2)+am2
        E3=rlu()*(E3max-am3)+am3
        Etotal=E1+E2+E3
        if(Etotal.gt.am0) goto 5
        E1=E1*(am0/Etotal)
        if(E1.ge.E1max) goto 5
        E2=E2*(am0/Etotal)
        if(E2.ge.E2max) goto 5
        E3=E3*(am0/Etotal)
        if(E3.ge.E3max) goto 5
        p1=E1**2-am1**2
        p2=E2**2-am2**2
        p3=E3**2-am3**2

        if(p1.gt.0.) then
           p1=sqrt(p1)
        else
          goto 5
        end if
        if(p2.gt.0.) then
           p2=sqrt(p2)
        else
          goto 5
        end if
        if(p3.gt.0.) then
           p3=sqrt(p3)
        else
          goto 5
        end if
        call ranvector3body(amomd1,1.)
        cang=(p3**2-p1**2-p2**2)/2/p1/p2     !cos ang 1-2 com
        if(abs(cang).gt.1.) then
          goto 5
        end if
        call scatvec3body(amomd1,cang,amomd2)

        do i=1,3
           amomd1(i)=amomd1(i)*p1
           amomd2(i)=amomd2(i)*p2
           amomd3(i)=-amomd1(i)-amomd2(i)
        end do

        return
        end

        subroutine scatvec3body(x,c,y)
*************************************************
*       returns vector y random in phi to x but
*       at an angle cos=c
*************************************************
        implicit none
        real x(3),y(3),z(3),w(3)
        real c,s
        integer k
        c=min(c,1.)
        c=max(-1.,c)
        if (abs(c).eq.1.) then
           s=0.
        else
           s=sqrt(1.00000001-c**2)
        endif
        call ranvector3body(z,1.)
        w(1)=x(2)*z(3)-x(3)*z(2)
        w(2)=x(3)*z(1)-x(1)*z(3)
        w(3)=x(1)*z(2)-x(2)*z(1)
        call normalize3body(w)
        do k=1,3
                y(k)=c*x(k)+s*w(k)
        end do
        return
        end

        subroutine normalize3body(x)
        implicit none
        real x(3),size
        integer i
        size=0.
        do i=1,3
          size=size+x(i)**2
        end do
        size=1./(sqrt(size)+.0000001)
        do i=1,3
          x(i)=x(i)*size
        end do
        return
        end

        subroutine ranvector3body(x,r)
        implicit none
        real x(3),r
        real rmod,rlu
        integer k
1       continue
        do k=1,3
          x(k)=2.*(rlu()-.5)
        end do
        rmod=x(1)**2+x(2)**2+x(3)**2
        if(rmod.gt.1.) goto 1
        call normalize3body(x)
        do k=1,3
          x(k)=x(k)*r
        end do
        return
        end
