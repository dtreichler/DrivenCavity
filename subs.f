cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gridset

      include 'dimension.h'

      real x(nx+2), xc(nx+2), dx
      real xg(nx+2,ny+2), xcg(nx+2,ny+2)
      real y(ny+2), yc(ny+2), dy
      real yg(nx+2,ny+2), ycg(nx+2,ny+2)

      common /grid/ x, y, xc, yc, dx, dy, rat

910   format (e21.14)

      pi = 4.*atan(1.)

      dx = 1.0/nx
      dy = rat/ny

      do i=1,nx+2
         x(i) = (i-1)*dx
         xc(i) = x(i)-0.5*dx
      enddo
      do j=1,ny+2
         y(j) = (j-1)*dy
         yc(j) = y(j)-0.5*dy
      enddo
      
!       do i=1,nx+2
!          do j=1,ny+2
!             xg (i,j) = x (i)
!             xcg(i,j) = xc(i)
!             yg (i,j) = y (j)
!             ycg(i,j) = yc(j)
!          enddo
!       enddo

      open(file='xc.out', unit=220)
      write(220,910) xc
      close(220)
      open(file='x.out', unit=221)
      write(221,910) x
      close(221)
      open(file='yc.out', unit=222)
      write(222,910) yc
      close(222)
      open(file='y.out', unit=223)
      write(223,910) y
      close(223)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE readold

      include 'dimension.h'

      real t, p(nx+2,ny+2), u(nx+2,ny+2)
      real v(nx+2,ny+2)
      common /soln/ u, v, p, t
      
      open(file='tf.out', unit=881)
      read(881,*) t
      close(881)
      open(file='p.out', unit=882)
      read(882,*) p
      close(882)
      open(file='u.out', unit=883)
      read(883,*) u
      close(883)
      open(file='v.out', unit=884)
      read(884,*) v
      close(884)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE initpoisson

      include 'dimension.h'

      common /grid/ x, y, xc, yc, dx, dy
      common /pois/ ae, aw, an, as, ap
      common /LU5/ ke, kw, kn, ks, kp
      parameter (lusr=0.)

      real ae(nx+2,ny+2), aw(nx+2,ny+2), an(nx+2,ny+2)
      real as(nx+2,ny+2), ap(nx+2,ny+2), ke(nx+2,ny+2)
      real kw(nx+2,ny+2), kn(nx+2,ny+2), ks(nx+2,ny+2)
      real kp(nx+2,ny+2), in(nx+2,ny+2), rhs(nx+2,ny+2)
      real x(nx+2), y(ny+2)
      real xc(nx+2), yc(ny+2)

c     Generate A matrix
      do i=2,nx+1
         do j=2,ny+1
            as(i,j) = 1./((y(j+1)-y(j))*(yc(j)-yc(j-1)))
            aw(i,j) = 1./((x(i+1)-x(i))*(xc(i)-xc(i-1)))
            ae(i,j) = 1./((x(i+1)-x(i))*(xc(i+1)-xc(i)))
            an(i,j) = 1./((y(j+1)-y(j))*(yc(j+1)-yc(j)))
            ap(i,j) = -(as(i,j)+aw(i,j)+ae(i,j)+an(i,j))
         enddo
      enddo
      
c     Homogeneous Neumann BC
      do i=1,nx+2
         j=2
         ap(i,j) = ap(i,j) + as(i,j)
         as(i,j) = 0.
         j=ny+1
         ap(i,j) = ap(i,j) + an(i,j)
         an(i,j) = 0.
      enddo
      do j=1,ny+2
         i=2
         ap(i,j) = ap(i,j) + aw(i,j)
         aw(i,j) = 0.
         i=ny+1
         ap(i,j) = ap(i,j) + ae(i,j)
         ae(i,j) = 0.
      enddo

c     Generate LU5 Coefficients

      do i=2,nx+1
         do j=2,ny+1
            ks(i,j)=as(i,j)
            kw(i,j)=aw(i,j)
            kp(i,j)=ap(i,j)-kw(i,j)*kn(i,j-1)-ks(i,j)*ke(i-1,j)
     $           +lusr*(kw(i,j)*ke(i,j-1)+ks(i,j)*kn(i-1,j))
            kn(i,j)=an(i,j)/kp(i,j)
            ke(i,j)=ae(i,j)/kp(i,j)
         enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE timestepcalc

      include 'dimension.h'

      real x(nx+2), xc(nx+2), dx
      real xg(nx+2,ny+2), xcg(nx+2,ny+2)
      real y(ny+2), yc(ny+2), dy
      real yg(nx+2,ny+2), ycg(nx+2,ny+2)

c     Solution Variables
      real u(nx+2,ny+2), v(nx+2,ny+2), p(nx+2,ny+2)
      real ustar(nx+2,ny+2), vstar(nx+2,ny+2)
      real umax, vmax
      real dt, dtc, dtv

      common /grid/ x, y, xc, yc, dx, dy
      common /soln/ u, v, p, t, dt, ustar, vstar

      real BOS, new, gng
      real CFL, sigma, nu
      real eps, tmax
      integer pitmax, nitmax, nout
      common /iter/ eps, pitmax, nout, BOS, CFL, sigma, nu

      umax = 0.
      vmax = 0.
      do i=2,nx+1
         do j=2,ny+1
            umax = max(umax,abs(u(i,j)))
            vmax = max(vmax,abs(v(i,j)))
         enddo
      enddo
      dtc = CFL/max(umax/dx,vmax/dy)
      dtv = sigma*min(dx**2,dy**2)/nu
      dt = min(dtc,dtv)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE bound(uarg,varg)

      include 'dimension.h'

      common /wall/ uw
      real uw
      real uarg(nx+2,ny+2), varg(nx+2,ny+2)

      do j=1,ny+2
         i=1
         uarg(i,j) = 0.
         varg(i,j) = -1.*varg(i+1,j)
         i=nx+2
         varg(i,j) = -1.*varg(i-1,j)
         uarg(nx+1,j) = 0.
      enddo

      do i=1,nx+2
         j=1
         varg(i,j) = 0.
         uarg(i,j) = -1.*uarg(i,j+1)
         j=ny+2
         uarg(i,j) = 2.*uw-uarg(i,j-1)
         varg(i,ny+1) = 0.
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine fluxes(rk)

      include 'dimension.h'

      integer rk    !<------ BASTARD
      real x(nx+2), xc(nx+2), dx
      real y(ny+2), yc(ny+2), dy
      real u(nx+2,ny+2), v(nx+2,ny+2), p(nx+2,ny+2)
      real ustar(nx+2,ny+2), vstar(nx+2,ny+2)
      real dt, dtc, dtv
      common /grid/ x, y, xc, yc, dx, dy
      common /soln/ u, v, p, t, dt, ustar, vstar

      real Fcx(nx+2,ny+2), Fcy(nx+2,ny+2), Gx(nx+2,ny+2)
      real Fvx(nx+2,ny+2), Fvy(nx+2,ny+2), Gy(nx+2,ny+2)
      real Fpx(nx+2,ny+2), Fpy(nx+2,ny+2), rhs(nx+2,ny+2)
      real Gxold(nx+2,ny+2), Gyold(nx+2,ny+2)
      common /flux/ Fcx, Fcy, Fvx, Fvy, Fpx, Fpy, Gx, Gy, rhs, Gxold,
     $     Gyold

      real gc(3), dtg(3), dtrk(3)
      common /rk3/ gc, dtg, dtrk

      do i=2,nx+1
         do j=2,ny+1
            Fcx(i,j)=-0.25*(((u(i+1,j)+u(i,j))**2-
     $           (u(i,j)+u(i-1,j))**2)/dx + 
     $           ((u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j))-
     $           (u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)))/dy)
            Fvx(i,j)=nu*(u(i+1,j)-2*u(i,j)+u(i-1,j))/dx**2+
     $           (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy**2
            Gx(i,j) = Fcx(i,j)+Fvx(i,j)-gc(rk)*Gx(i,j)
            rhs(i,j) = 0.
            Fcy(i,j)=-0.25*(((v(i,j)+v(i,j+1))**2-
     $           (v(i,j)+v(i,j-1))**2)/dy + 
     $           ((v(i+1,j)+v(i,j))*(u(i,j+1)+u(i,j))-
     $           (v(i,j)+v(i-1,j))*(u(i-1,j+1)+u(i-1,j)))/dx)
            Fvy(i,j)=nu*(v(i+1,j)-2*v(i,j)+v(i-1,j))/dx**2+
     $           (v(i,j+1)-2*v(i,j)+v(i,j-1))/dy**2
            Gy(i,j) = Fcy(i,j)+Fvy(i,j)-gc(rk)*Gy(i,j)
            
         enddo
      enddo

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE timeadv

      include 'dimension.h'

      integer rk

c     Grid
      real x(nx+2), xc(nx+2), dx
      real xg(nx+2,ny+2), xcg(nx+2,ny+2)
      real y(ny+2), yc(ny+2), dy
      real yg(nx+2,ny+2), ycg(nx+2,ny+2)

c     Solution Variables
      real u(nx+2,ny+2), v(nx+2,ny+2), p(nx+2,ny+2)
      real ustar(nx+2,ny+2), vstar(nx+2,ny+2), uold(nx+2,ny+2)
      real umax, vmax, udiffnorm
      real dt, dtc, dtv
      real dtrhs

      common /grid/ x, y, xc, yc, dx, dy
      common /soln/ u, v, p, t, dt, ustar, vstar

c     Fluxes
      real Fcx(nx+2,ny+2), Fcy(nx+2,ny+2), Gx(nx+2,ny+2)
      real Fvx(nx+2,ny+2), Fvy(nx+2,ny+2), Gy(nx+2,ny+2)
      real Fpx(nx+2,ny+2), Fpy(nx+2,ny+2), rhs(nx+2,ny+2)
      real Gxold(nx+2,ny+2), Gyold(nx+2,ny+2)
      common /flux/ Fcx, Fcy, Fvx, Fvy, Fpx, Fpy, Gx, Gy, rhs, Gxold,
     $     Gyold

c     Stepping Coefficients
      real gc(3), dtg(3), dtrk(3)
      common /rk3/ gc, dtg, dtrk, udiffnorm

      do i=1,nx+2
         do j=1,ny+2
            uold(i,j) = u(i,j)
         enddo
      enddo

      do rk=1,3
         call fluxes(rk)

c         do i=1,nx+2
c            do j=1,ny+2
c               ustar(i,j) = 0.
c               vstar(i,j) = 0.
c               rhs(i,j) = 0.
c            enddo
c         enddo

         do i=2,nx
            do j=2,ny+1
               ustar(i,j) = u(i,j)+dtg(rk)*dt*Gx(i,j)
            enddo
         enddo

         do i=2,nx+1
            do j=2,ny
               vstar(i,j) = v(i,j)+dtg(rk)*dt*Gy(i,j)
            enddo
         enddo

         call bound(ustar,vstar)

         dtrhs = 1.D0/real(dtrk(rk)*dt)


         do i=2,nx+1
            do j=2,ny+1
               rhs(i,j) = ((ustar(i,j)-ustar(i-1,j))/dx
     $              +(vstar(i,j)-vstar(i,j-1))/dy)*dtrhs
            enddo
         enddo
         
         t = t+dtrk(rk)*dt

         call poisson(p,rhs,t)

         do i=2,nx+1
            do j=2,ny+1
               Fpx(i,j) = -(p(i+1,j)-p(i,j))*(1.D0/dx)
               Fpy(i,j) = -(p(i,j+1)-p(i,j))*(1.D0/dy)
            enddo
         enddo

         do i=2,nx
            do j=2,ny+1
               u(i,j) = ustar(i,j) + dtrk(rk)*dt*Fpx(i,j)
            enddo
         enddo

         do i=2,nx+1
            do j=2,ny
               v(i,j) = vstar(i,j) + dtrk(rk)*dt*Fpy(i,j)
            enddo
         enddo

         call bound(u,v)

      enddo

      udiffnorm = 0.
      do i=2,nx+1
         do j=2,ny+1
            udiffnorm = udiffnorm + (u(i,j)-uold(i,j))**2
         enddo
      enddo
      udiffnorm = sqrt(udiffnorm/nx/ny)
      
      return
      end
