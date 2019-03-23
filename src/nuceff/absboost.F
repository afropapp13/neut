c ---------------------------------------------------------------
c     ------------------------------
      subroutine absboost(beta,pin,pout)
c     ------------------------------
c
c     ( purpose )
c       Boost 4-vectors between reference frames
c
c     ( input )
c       beta(3) relative velocity of two frames
c       pin(4)  four vector (E,px,py,pz) initial frame
c
c     ( output )
c       pout(4)  four vector (E,px,py,pz) finaal frame
c
      implicit none
      
      real*4 beta(3),pin(4),pout(4)

      real*4 bx,by,bz,b2,g,gm1

      bx  = beta(1)
      by  = beta(2)
      bz  = beta(3)
      b2  = bx**2 + by**2 + bz**2
      g   = sqrt(1./(1.-b2))
      gm1 = g - 1.

      pout(1) = g*pin(1) + g*bx*pin(2) + g*by*pin(3) + g*bz*pin(4)
      pout(2) = g*bx*pin(1) + (1.+gm1*bx**2/b2)*pin(2) +
     &          (gm1*bx*by/b2)*pin(3) + (gm1*bx*bz/b2)*pin(4)
      pout(3) = g*by*pin(1) + (gm1*bx*by/b2)*pin(2) +
     &          (1.+gm1*by**2/b2)*pin(3) + (gm1*by*bz/b2)*pin(4)
      pout(4) = g*bz*pin(1) + (gm1*bx*bz/b2)*pin(2) +
     &          (gm1*by*bz/b2)*pin(3) + (1.+gm1*bz**2/b2)*pin(4)

      end
