c Derrick Treichler
c 8 Nov 07
c ENME646

c Navier-Stokes Solver, Take 2

      PROGRAM SOLVER2

      include 'dimension.h'

c     Grid
      real x(nx+2), xc(nx+2), dx
      real xg(nx+2,ny+2), xcg(nx+2,ny+2)
      real y(ny+2), yc(ny+2), dy
      real yg(nx+2,ny+2), ycg(nx+2,ny+2)

c     Solution Variables
      real u(nx+2,ny+2), v(nx+2,ny+2), p(nx+2,ny+2)
      real ustar(nx+2,ny+2), vstar(nx+2,ny+2), vort(nx+2,ny+2)
      real umax, vmax
      real dt, dtc, dtv
      real stime0, stime, rat

      common /grid/ x, y, xc, yc, dx, dy, rat
      common /soln/ u, v, p, t, dt, ustar, vstar
      common /res/ nrmres, piter
      integer piter
      real nrmres, udiffnormuw

c     Poisson
      real ap(nx+2,ny+2), ae(nx+2,ny+2), aw(nx+2,ny+2)
      real as(nx+2,ny+2), an(nx+2,ny+2), kp(nx+2,ny+2)
      real kn(nx+2,ny+2), ks(nx+2,ny+2), ke(nx+2,ny+2)
      real kw(nx+2,ny+2)

c     Stepping Coefficients
      real gc(3), dtg(3), dtrk(3)

c     Difference norm for each RK3 step
      real udiffnorm

      common /rk3/ gc, dtg, dtrk, udiffnorm

c     Fluxes
      real Fcx(nx+2,ny+2), Fcy(nx+2,ny+2), Gx(nx+2,ny+2)
      real Fvx(nx+2,ny+2), Fvy(nx+2,ny+2), Gy(nx+2,ny+2)
      real Fpx(nx+2,ny+2), Fpy(nx+2,ny+2), rhs(nx+2,ny+2)
      real Gxold(nx+2,ny+2), Gyold(nx+2,ny+2)
      common /flux/ Fcx, Fcy, Fvx, Fvy, Fpx, Fpy, Gx, Gy, rhs, Gxold,
     $     Gyold

c     Parameters
      real BOS, new, gng
      real CFL, sigma, nu
      real eps, tmax, uw
      real diffmin
      integer pitmax, nitmax, nout
      common /iter/ eps, pitmax, nout, BOS, CFL, sigma, nu
      common /wall/ uw
      real uexamax, vexamax, pexamax
      common /err/ uexamax, vexamax, pexamax
      

 910  format (e46.36)
 810  format (/,1x,A,i5)
 811  format (1x,2(A,g25.15))
 813  format (i6,1p5(e17.8),i7)


c     Read input
      open(file='var.inp', unit=11)
      read(11,*) BOS
      read(11,*) new
      read(11,*) dtcc
      read(11,*) dt
      read(11,*) CFL
      read(11,*) sigma
      read(11,*) nu
      read(11,*) nitmax
      read(11,*) pitmax
      read(11,*) nout
      read(11,*) eps
      read(11,*) tmax
      read(11,*) diffmin
      read(11,*) uw
      read(11,*) rat
      
      close(11)

      call gridset
     
      if (new .eq. 1) then
         t = 0.
         do i=2,nx+1
            do j=2,ny+1
               p(i,j) = 0.
               u(i,j) = 0.
               v(i,j) = 0.
            enddo
         enddo
         call bound(u,v)
         open(file='iter.out', unit=99, status='replace')
         close(99)
      else
         call readold
      endif

      open(file='gng.inp', unit=59)
      write(59,'(i)') 1
      close(59)

      do i=1,nx+2
         do j=1,ny+2
            Gx(i,j) = 0.
            Gy(i,j) = 0.
            Gxold(i,j) = 0.
            Gyold(i,j) = 0.
         enddo
      enddo

      call initpoisson

c     RK3 Coefficients
      gc  (1) = 0.
      gc  (2) = 5./9.
      gc  (3) = 153./128.
      dtg (1) = 1./3.
      dtg (2) = 15./16.
      dtg (3) = 8./15.
      dtrk(1) = 1./3.
      dtrk(2) = 5./12.
      dtrk(3) = 1./4.


      call cpu_time(stime0)
      do n = 1,nitmax

         open(file='gng.inp', unit=59)
         read(59,*) gng
         close(59)
         if (gng .eq. 0) then
            write(6,'(/A)') 'gng = 0, Writing variables and quitting.'
            goto 1000
         endif

c     Calculate time step
         if (dtcc .eq. 0) then
            call timestepcalc
         endif

         open(file='iter.out', unit=99, status='replace')
         write(99,810) 'n = ', n
         write(99,811) 't = ', t, ' dt = ', dt
         close(99)

         call timeadv

         call cpu_time(stime)

         udiffnormuw = udiffnorm/uw
         write(6,813) n, stime-stime0, t, dt, nrmres, udiffnormuw, piter

         if ((t.gt.tmax)) then!.or.(udiffnorm.lt.diffmin)) then
            if (t.gt.tmax) then
               write(6,*)
               write(6,*) 't > tmax, Writing variables and quitting.'
            else
               write(6,*)
               write(6,*) '|diff| < min, Writing',
     &             ' variables and quitting.'
            endif
            goto 1000
         endif
      enddo

 1000 continue

      do i=2,nx+1
         do j=2,ny+1
            vort(i,j) = (v(i+1,j)-v(i,j))/dx + (u(i,j+1)-u(i,j))/dy
         enddo
      enddo

      open(file='tf.out', unit=881)
      write(881,910) t
      close(881)
      open(file='p.out', unit=882)
      write(882,910) p
      close(882)
      open(file='u.out', unit=883)
      write(883,910) u
      close(883)
      open(file='v.out', unit=884)
      write(884,910) v
      close(884)
      open(file='rhs.out', unit=885)
      write(885,910) rhs  
      close(885)
      open(file='ustar.out', unit=886)
      write(886,910) ustar  
      close(886)
      open(file='vstar.out', unit=887)
      write(887,910) vstar
      close(887)
      open(file='vort.out', unit=888)
      write(888,910) vort
      close(888)
      open(file='Fcx.out', unit=885)
      write(885,910) Fcx  
      close(885)
      open(file='Fcy.out', unit=886)
      write(886,910) Fcy  
      close(886)
      open(file='Fvx.out', unit=887)
      write(887,910) Fvx
      close(887)
      open(file='Fvy.out', unit=888)
      write(888,910) Fvy
      close(888)

 601  format (A,e26.16)
 602  format (A,i5)

      write(6,*)
      write(6,*) '############# CLEAN STOP ##############'
      write(6,602) '         nx  =', nx
      write(6,602) '         ny  =', ny
      write(6,601) '       CPU  =', stime-stime0
      write(6,601) '        tf  =', t
      write(6,601) '       dtf  =', dt
      write(6,601) '       CFL  =', CFL
      write(6,601) '     sigma  =', sigma
      write(6,601) '     |res|  =', nrmres
      write(6,*) '#######################################'
      stop
      end
