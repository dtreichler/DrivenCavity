c     Derrick Treichler
c     30 Oct 2007
c     Functions exau exav exap provide the exact solution to the Taylor-Green
c     vortex at supplied time t

      function exau(x,y,t)

      exau = -exp(-2*t)*cos(x)*sin(y)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function exav(x,y,t)

      exav = exp(-2*t)*sin(x)*cos(y)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function exap(x,y,t)

      exap = -0.25*exp(-4*t)*(cos(2*x)+cos(2*y))

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
