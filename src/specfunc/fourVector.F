c fourVector notation is used - fourVectors in arrays of (E, px, py, pz)
c arrays indexed from 1 such that time component is element 1
c space components x, y, z are elements 2, 3, 4


c small function that returns the 4-vector dot product
      double precision function dotProduct( a, b )
        double precision a(4), b(4)

        dotProduct = (a(1) * b(1)) - (a(2) * b(2)) - (a(3) * b(3)) -
     +  (a(4) * b(4))

       return
       end

c boost subroutine - takes a velocity vector as the space components
c of a 4-vector
      subroutine boost(vec, v)
        implicit none
        double precision vec(4), temp_vec(4)
        double precision v(4)
        double precision gam, v2, modMomentum
        integer i

        gam = 1. / sqrt( 1. - modMomentum(v)**2)
c        write(*,*) gam
        v2 = modMomentum(v)**2
c        write(*,*) v2

        temp_vec(1) = gam*(vec(1) - v(2)*vec(2) - v(3)*vec(3)
     +                - v(4)*vec(4) )
        temp_vec(2) = ( - vec(1)*v(2)*gam
     +           + vec(2)*(1 + (gam - 1)*(v(2)*v(2))/v2)
     +           + vec(3)*(gam - 1)*v(2)*v(3)/v2
     +           + vec(4)*(gam - 1)*v(2)*v(4)/v2 )
        temp_vec(3) = ( - vec(1)*v(3)*gam
     +           + vec(2)*(gam - 1)*(v(2)*v(3))/v2
     +           + vec(3)*(1 + (gam - 1)*v(3)*v(3)/v2)
     +           + vec(4)*(gam - 1)*v(3)*v(4)/v2 )
        temp_vec(4) = ( - vec(1)*v(4)*gam
     +           + vec(2)*(gam - 1)*(v(2)*v(4))/v2
     +           + vec(3)*(gam - 1)*v(3)*v(4)/v2
     +           + vec(4)*(1 + (gam - 1)*v(4)*v(4)/v2 ))

        do i = 1, 4
          vec(i) = temp_vec(i)
        enddo

        return
        end

c function calculates the magnitude of threeMomentum from a FourVector
      double precision function modMomentum(fourVector)
        implicit none
        double precision fourVector(*)
        double precision squareMomentum
        integer i
        squareMomentum = 0
        do i = 2, 4
          squareMomentum = squareMomentum + fourVector(i)**2
        enddo
        modMomentum = sqrt(squareMomentum)

        return
        end


c function returns the invariant mass of a four-vector
      double precision function getInvMass(vector)
        implicit none
        double precision vector(4)
        double precision invMass, modMomentum
        external modMomentum

        invMass = vector(1)**2 - modMomentum(vector)**2

        getInvMass = invMass

        return
        end

c subroutine takes a mass and a 4-vector with energy component, and sets the
c momentum components to a random phi direction, with known theta and the right
c magnitude
      subroutine pickDirection_fixtheta(vector, mass, theta)
        double precision vector(4)
        double precision theta, phi, mag, mass, dummy
        real RLU
        external RLU

        phi = dble(RLU(dummy)) * 2 * 3.14159265358
        mag = sqrt(vector(1)**2 - mass**2)

        vector(2) = mag*sin(theta)*cos(phi)
        vector(3) = mag*sin(theta)*sin(phi)
        vector(4) = mag*cos(theta)

        return
        end

c subroutine takes a mass and a 4-vector with energy component, and sets the
c momentum components to a random direction with the right magnitude
      subroutine pickDirection(vector, mass)
        double precision vector(4)
        double precision theta, mass
        real RLU
        external RLU

        theta = acos(dble(RLU(dummy))*2 - 1)
        call pickDirection_fixtheta(vector, mass, theta)
        end


c this subroutine rotates a vector about the z-axis
      subroutine setPhi(fourVector, phi)
        implicit none
        double precision fourVector(4)
        double precision phi
        double precision theta, length
        integer i
        theta = 0
        length = 0
        do i = 2, 4
          length = sqrt( length**2 + fourVector(i)**2 )
        enddo

        theta = acos( fourVector(4)/ length)
        fourVector(2) = length*sin(theta)*cos(phi)
        fourVector(3) = length*sin(theta)*sin(phi)
        fourVector(4) = length*cos(theta)
        return
        end

c this function returns the value of theta
      double precision function getTheta(vector)
        implicit none
        double precision vector(4)
        double precision mag, theta
        double precision modMomentum

        mag = modMomentum(vector)
        theta = acos(vector(4) / mag)
        getTheta = theta

        return
        end
