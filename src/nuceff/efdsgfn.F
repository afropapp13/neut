
      function sind(x)
      real x, pi
      pi = acos(-1.)
      sind = sin(x*pi/180.)
      end

      function cosd(x)
      real x, pi
      pi = acos(-1.)
      cosd = cos(x*pi/180.)
      end

      function acosd(x)
      real x, pi
      pi = acos(-1.)
      acosd = 180.*acos(x)/pi
      end

      function atand(x)
      real x, pi
      pi = acos(-1.)
      atand = 180.*atan(x)/pi
      end
