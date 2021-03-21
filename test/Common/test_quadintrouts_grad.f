      implicit real *8 (a-h,o-z)
      real *8 xyztarg(3,1000),xyztarg0(3,1000)
      complex *16, allocatable :: slp_grad(:,:,:)
      complex *16, allocatable :: slp_grad_ex(:,:,:)
      integer ipars(10)
      integer igxyz(3,3)
      integer, allocatable :: iprint(:,:)

      real *8 dpars(5)
      complex * 16 zpars(3),ima

      external h3d_slp,h3d_sgradx,h3d_sgrady,h3d_sgradz

      character type

      ima = dcmplx(0.0d0,1.0d0)
      
      done = 1
      pi = atan(done)*4


      call prini(6,13)
      call prinf('enter n*',n,0)
      read *, n


      norder = 8
      npols = (norder+1)*norder/2


      
      ntarg = 3

c
cc       note the exact integrals are hardcoded for the
c        the following geometry
c
c        update documentation needed
c
c        see mathematica notebook quadintegrals.nb
c        
      

      zk = 1.1d0

      zpars(1) = zk
      xyztarg0(1,1) = 0
      xyztarg0(2,1) = 0.1d0
      xyztarg0(3,1) = 1.019d0

      xyztarg0(1,2) = 1.2d0
      xyztarg0(2,2) = 0.1d0
      xyztarg0(3,2) = -2.0d0/7.0d0

      xyztarg0(1,3) = 0.3d0
      xyztarg0(2,3) = 1.4d0
      xyztarg0(3,3) = -3.0d0/7.0d0

      xyztarg(1,1) = xyztarg0(1,1) 
      xyztarg(2,1) = xyztarg0(2,1)
      xyztarg(3,1) = xyztarg0(3,1) - 1.0d0
      igxyz(1,1) = 1
      igxyz(2,1) = 2
      igxyz(3,1) = 3

      xyztarg(1,2) = xyztarg0(2,2) 
      xyztarg(2,2) = xyztarg0(3,2)
      xyztarg(3,2) = xyztarg0(1,2) - 1.0d0
      igxyz(1,2) = 3
      igxyz(2,2) = 1
      igxyz(3,2) = 2

      xyztarg(1,3) = xyztarg0(1,3) 
      xyztarg(2,3) = xyztarg0(3,3)
      xyztarg(3,3) = xyztarg0(2,3) - 1.0d0

      igxyz(1,3) = 1
      igxyz(2,3) = 3
      igxyz(3,3) = 2

      allocate(slp_grad(npols,ntarg,3),slp_grad_ex(npols,ntarg,3))

      do l=1,3
        do i=1,ntarg
          do j=1,npols
            slp_grad(j,i,l) = 0
            slp_grad_ex(j,i,l) = 0
          enddo
        enddo
      enddo

      eps = 1.0d-6
      nqorder = 20
      intype = 2
      nquadmax = 5000

      type = 't'
      call cquadints_adap(eps,intype,norder,type,npols,ntarg,xyztarg,
     1 nquadmax,h3d_sgradx,dpars,zpars,ipars,nqorder,slp_grad(1,1,1))

      call cquadints_adap(eps,intype,norder,type,npols,ntarg,xyztarg,
     1 nquadmax,h3d_sgrady,dpars,zpars,ipars,nqorder,slp_grad(1,1,2))


      call cquadints_adap(eps,intype,norder,type,npols,ntarg,xyztarg,
     1 nquadmax,h3d_sgradz,dpars,zpars,ipars,nqorder,slp_grad(1,1,3))
      nn = norder + 1

      slp_grad_ex(1,1,1) = 0.0d0 
      slp_grad_ex(2,1,1) = 4.108907986135193d0+0.527199071666214d0*ima
      slp_grad_ex(nn,1,1) = 0.0d0 

      slp_grad_ex(1,1,2)=-0.4205878887206168d0-0.1503250920833929d0*ima
      slp_grad_ex(2,1,2)= 0
      slp_grad_ex(nn,1,2)=4.072371328901900d0+0.525911095874265d0*ima

      slp_grad_ex(1,1,3)=-6.244445243530355d0-0.031068048605975d0*ima 
      slp_grad_ex(2,1,3)= 0
      slp_grad_ex(nn,1,3)=-0.6164533544673438d0-0.0002506281110130d0*ima

      slp_grad_ex(1,2,1)=-5.691117319341798d0 - 0.322333025299182d0*ima
      slp_grad_ex(2,2,1)=-0.4911891850156081d0-0.0026109845472635d0*ima
      slp_grad_ex(nn,2,1)=1.3928919779596355d0+0.0074591230001994d0*ima 

      slp_grad_ex(1,2,2)=-0.3896127345928584d0-0.1481115899132737d0*ima
      slp_grad_ex(2,2,2)=2.989587349251013d0+0.518276183191058d0*ima
      slp_grad_ex(nn,2,2)=0.0510418033821068d0+0.00351414240705731d0*ima

      slp_grad_ex(1,2,3)=1.1657077092615330d0+0.4231801354264063d0*ima
      slp_grad_ex(2,2,3)=0.05629465013887204d0+0.00351455883161975d0*ima
      slp_grad_ex(nn,2,3)=2.744742426132103d0+0.509146891660142d0*ima

      slp_grad_ex(1,3,1)=-0.9956327486653673d0-0.4280336885881407d0*ima
      slp_grad_ex(2,3,1)=1.802182693818707d0+0.489594064704992d0*ima
      slp_grad_ex(nn,3,1)=0.1900417554895804d0+0.0153997531717577d0*ima

      slp_grad_ex(1,3,2)=-4.712302914772059d0-0.621584380952863d0*ima
      slp_grad_ex(2,3,2)=-1.0154377749654707d0-0.0152618388506026d0*ima
      slp_grad_ex(nn,3,2)=1.424266708898861d0+0.021799413779375d0*ima

      slp_grad_ex(1,3,3)=1.488229207865768d0+0.611484730858201d0*ima
      slp_grad_ex(2,3,3)=0.2098107684212962d0+0.0154021630343004d0*ima
      slp_grad_ex(nn,3,3)=1.553385382842890d0+0.477967230602385d0*ima

      do l=1,3
        do i=3,3
          print *, "itarg =", i
          print *, ""
          print *, ""

          ll = igxyz(l,i)

          err1s = abs(slp_grad(1,i,ll)-slp_grad_ex(1,i,l))
          err2s = abs(slp_grad(2,i,ll)-slp_grad_ex(2,i,l))
          err3s = abs(slp_grad(nn,i,ll)-slp_grad_ex(nn,i,l))

          print *, "0,0 int=",err1s
          print *, "1,0 int=",err2s
          print *, "0,1 int=",err3s
          print *, ""
          print *, ""
        enddo
      enddo

      stop
      end
c
c
c
c
c


      subroutine h3d_slp(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)

      f = exp(ima*zpars(1)*rr)/rr

      return
      end

c      
c      
c
c
c
c
c

      subroutine h3d_dlp(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,z

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)
      z = ima*zpars(1)*rr

      f = exp(z)*y(3)*(z-1.0d0)/rr**3

      return
      end

c
c
c
c     
c
c
      subroutine h3d_sgradx(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,z

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)
      z = ima*zpars(1)*rr

      f = exp(z)*(y(1)-x(1))*(z-1.0d0)/rr**3

      return
      end

c
c
c
c     
c
c
      subroutine h3d_sgrady(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,z

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)
      z = ima*zpars(1)*rr

      f = exp(z)*(y(2)-x(2))*(z-1.0d0)/rr**3

      return
      end

c
c
c
c
c     
c
c
      subroutine h3d_sgradz(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,z

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)
      z = ima*zpars(1)*rr

      f = exp(z)*y(3)*(z-1.0d0)/rr**3

      return
      end

c
c     
c
c
      subroutine h3d_dgradx(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,z

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)
      z = ima*zpars(1)*rr

      f = -exp(z)*(y(1)-x(1))*y(3)*(-3.0d0 + 3*z - z**2)/rr**5

      return
      end

c
c
c
c     
c
c
      subroutine h3d_dgrady(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,z

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)
      z = ima*zpars(1)*rr

      f = -exp(z)*(y(2)-x(2))*y(3)*(-3.0d0 + 3*z - z**2)/rr**5

      return
      end

c
c
c
c
c     
c
c
      subroutine h3d_dgradz(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,z

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)
      z = ima*zpars(1)*rr
      f=-exp(z)*((1.0d0-z)*rr**2+y(3)*y(3)*(-3.0d0 + 3*z - z**2))/rr**5
      


      return
      end

c
c
c
c
c
c
c

      subroutine hslp(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f
      
      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)

      f = exp(ima*zpars(1)*rr)/rr

      return
      end

c      
c      
c
c
c
c
c

      subroutine hdlp(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,z
      
      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)
      z = ima*zpars(1)*rr

      f = exp(z)*y(3)*(1-z)/rr**3

      return
      end

