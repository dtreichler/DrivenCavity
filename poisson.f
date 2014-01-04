c     Derrick Treichler
c     30 Oct 2007
c     Function poisson solves Poisson's equation given an rhs and a guess
c     Uses Homogeneous BCs and BiCGSTAB with LU5 Preconditioning or SOR

      subroutine poisson(u,rhs,ti)
      include 'dimension.h'
      
      common /grid/ x, y, xc, yc, dx, dy
      common /pois/ ae, aw, an, as, ap
      common /LU5/ ke, kw, kn, ks, kp
      common /iter/ eps, pitmax, nout, BOS
      common /res/ nrmres, piter

      real ae(nx+2,ny+2), aw(nx+2,ny+2), an(nx+2,ny+2)
      real as(nx+2,ny+2), ap(nx+2,ny+2), ke(nx+2,ny+2)
      real kw(nx+2,ny+2), kn(nx+2,ny+2), ks(nx+2,ny+2)
      real kp(nx+2,ny+2), rhs(nx+2,ny+2)
      real x(nx+2), y(ny+2)
      real xc(nx+2), yc(ny+2)
      real u(nx+2,ny+2)

      real res(nx+2,ny+2)
      real nrmres, nrmres0
      real nrmdif, time
      real time0, dif(nx+2,ny+2), uold(nx+2,ny+2)

      real omega, ut, rhojac
      real BOS
      integer pitmax, piter

      real rhoold, alpha, w
      real rbar0(nx+2,ny+2), p(nx+2,ny+2), v(ny+2,ny+2)
      real rho, beta, psi(nx+2,ny+2)
      real eta(nx+2,ny+2), chi(nx+2,ny+2), zeta(nx+2,ny+2)
      real rbarv0, ktkt, ktks
      real kit(nx+2,ny+2), kis(nx+2,ny+2), t(nx+2,ny+2)
      real s(nx+2,ny+2), ts, tt
      real nrmreslim


990   format (i6,1p3e17.8)

      pi = 4.*atan(1.)
      nrmreslim = eps
      open(file='iter.out', unit=99, access='append')

c      write(6,*) 

      nrmres = 0.
      do i=2,nx+1
         do j=2,ny+1
            res(i,j) = rhs(i,j)-
     $           ap(i,j)*u(i,j)-an(i,j)*u(i,j+1)
     $           -as(i,j)*u(i,j-1)-ae(i,j)*u(i+1,j)
     $           -aw(i,j)*u(i-1,j)
            nrmres = nrmres + res(i,j)**2
         enddo
      enddo
      nrmres0 = sqrt(nrmres/ny/ny)

      call cpu_time(time0)

c      write(6,*) '  iter       CPU         ',
c     $      ' |res|/|res0|        |diff|'
c      write(6,990) 0, 0., nrmres0, 0.

      write(99,*) '  iter       CPU         ',
     $      ' |res|/|res0|        |diff|'
      write(99,990) 0, 0., nrmres0, 0.

c     BiCGSTAB
      if (BOS .eq. 0) then
         rho = 1
         alpha = 1
         w = 1
         do i=2,nx+1
            do j=2,ny+1
               rbar0(i,j)=res(i,j)
               p(i,j)=0.
               v(i,j)=0.
            enddo
         enddo
         
         do n=1,pitmax
            do i=2,nx+1
               do j=2,ny+1
                  uold(i,j) = u(i,j)
               enddo
            enddo
            rhoold = rho
            rho = 0.
            do i=2,nx+1
               do j=2,ny+1
                  rho=rho+rbar0(i,j)*res(i,j)
               enddo
            enddo

            beta=(rho/rhoold)*(alpha/w)

            do i=2,nx+1
               do j=2,ny+1
                  p(i,j)=res(i,j)+beta*(p(i,j)-w*v(i,j))
               enddo
            enddo
            
            
c     Solve K*eta=p --> LU*eta=p --> L*psi=p --> U*eta=psi
            do i=2,nx+1
               do j=2,ny+1
                  psi(i,j)=(1./kp(i,j))*(p(i,j)-ks(i,j)*psi(i-1,j)
     $                 -kw(i,j)*psi(i,j-1))
               enddo
            enddo
          
            do i=nx+1,2,-1
               do j=ny+1,2,-1
                  eta(i,j)=(psi(i,j)-ke(i,j)*eta(i,j+1)
     $                 -kn(i,j)*eta(i+1,j))
               enddo
            enddo
            
c     v=A*eta, calcluate inner product needed for alpha
            rbarv=0.
            do i=2,nx+1
               do j=2,ny+1
                  v(i,j)=ap(i,j)*eta(i,j)+ae(i,j)*eta(i+1,j)
     $                 +an(i,j)*eta(i,j+1)+as(i,j)*eta(i,j-1)
     $                 +aw(i,j)*eta(i-1,j)
                  rbarv=rbarv+rbar0(i,j)*v(i,j)
               enddo
            enddo
            alpha=rho/rbarv
                          
            do i=2,nx+1
               do j=2,ny+1
                  s(i,j)=res(i,j)-alpha*v(i,j)
               enddo
            enddo
            
c     Solve K*zeta=s --> LU*zeta=s --> L*chi=s --> U*zeta=chi
            do i=2,nx+1
               do j=2,ny+1
                  chi(i,j)=(1./kp(i,j))*(s(i,j)-as(i,j)*chi(i-1,j)
     $                 -aw(i,j)*chi(i,j-1))
               enddo
            enddo
            
            do i=nx+1,2,-1
               do j=ny+1,2,-1
                  zeta(i,j)=(chi(i,j)-ke(i,j)*zeta(i,j+1)
     $                 -kn(i,j)*zeta(i+1,j))
               enddo
            enddo
            
c     t=A*zeta, calculate inner products needed for new w
            ktks=0
            ktkt=0
            do i=2,nx+1
               do j=2,ny+1
                  t(i,j)=ap(i,j)*zeta(i,j)+ae(i,j)*zeta(i+1,j)
     $                 +an(i,j)*zeta(i,j+1)+as(i,j)*zeta(i,j-1)
     $                 +aw(i,j)*zeta(i-1,j)
                  kit(i,j)=(1./kp(i,j))*(t(i,j)-as(i,j)*kit(i-1,j)
     $                 -aw(i,j)*kit(i,j-1))
                  kis(i,j)=(1./kp(i,j))*(s(i,j)-as(i,j)*kis(i-1,j)
     $                 -aw(i,j)*kis(i,j-1))
                  ktks=ktks+kit(i,j)*kis(i,j)
                  ktkt=ktkt+kit(i,j)*kit(i,j)
               enddo
            enddo
            w=ktks/ktkt
            
            do i=2,nx+1
               do j=2,ny+1
                  u(i,j)=uold(i,j)+alpha*eta(i,j)+w*zeta(i,j)
               enddo
            enddo
            
            if (MOD(n,nout) .eq. 0) then
               nrmres = 0.
               nrmdif = 0.
               do i=2,nx+1
                  do j=2,ny+1
                     dif(i,j) = uold(i,j) - u(i,j)
                     nrmdif = nrmdif + dif(i,j)**2
                     res(i,j) = rhs(i,j)-
     $                    ap(i,j)*u(i,j)-an(i,j)*u(i,j+1)
     $                    -as(i,j)*u(i,j-1)-ae(i,j)*u(i+1,j)
     $                    -aw(i,j)*u(i-1,j)
                     nrmres = nrmres + res(i,j)**2
                  enddo
               enddo
               call cpu_time(time)
               nrmdif = sqrt(nrmdif/nx/ny)
               nrmres = sqrt(nrmres/nx/ny)
c               write(6,990) n, time-time0, nrmres/nrmres0, 
c     $              nrmdif
               write(99,990) n, time-time0, nrmres/nrmres0, 
     $              nrmdif
               
               if (nrmres .lt. eps) then
                  go to 100
               endif
               if ((nrmres/nrmres0).lt.nrmreslim) then
                  goto 100
               endif
            else
               do i=2,nx+1
                  do j=2,ny+1
                     res(i,j) = s(i,j)-w*t(i,j)
                  enddo
               enddo                     
            endif
         enddo
      elseif (BOS .eq. 1) then
         rhojac = 1. - pi**2/nx**2
         omega = 2./(1.+sqrt(1.-rhojac**2))

         do n=1,pitmax
            do i=1,nx+2
               do j=1,ny+2
                  uold(i,j)=u(i,j)
               enddo
            enddo

            do i=2,nx+1
               do j=2,ny+1
                  ut = (rhs(i,j)-(as(i,j)*u(i,j-1)+an(i,j)*u(i,j+1)
     $                +ae(i,j)*u(i+1,j)+aw(i,j)*u(i-1,j)))/ap(i,j)
                  u(i,j) = omega*ut + (1-omega)*u(i,j)
               enddo
            enddo
            if (MOD(n,nout) .eq. 0) then
               nrmres = 0.
               nrmdif = 0.
               do i=2,nx+1
                  do j=2,ny+1
                     dif(i,j) = uold(i,j) - u(i,j)
                     nrmdif = nrmdif + dif(i,j)**2
                     res(i,j) = rhs(i,j)-
     $                 ap(i,j)*u(i,j)-an(i,j)*u(i,j+1)
     $                 -as(i,j)*u(i,j-1)-ae(i,j)*u(i+1,j)
     $                 -aw(i,j)*u(i-1,j)
                     nrmres = nrmres + res(i,j)**2
                  enddo
               enddo
               call cpu_time(time)
               nrmdif = sqrt(nrmdif/nx/ny)
               nrmres = sqrt(nrmres/nx/ny)
               write(6,990) n, time-time0, nrmres/nrmres0, 
     $              nrmdif
               write(99,990) n, time-time0, nrmres/nrmres0, 
     $              nrmdif
               if (nrmres.lt.eps) then
                  goto 100
               endif
c              if ((nrmres/nrmres0).lt.nrmreslim) then
c                 goto 100
c              endif
            endif
         enddo
      endif


100   continue

      piter = n

      do i=1,nx+2
         u(i,1) = u(i,2)
         u(i,ny+2) = u(i,ny+1)
      enddo
      do j=1,ny+2
         u(1,j) = u(2,j)
         u(nx+2,j) = u(nx+1,j)
      enddo
         
      close(99)
      
      return
      end
