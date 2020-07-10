      implicit real *8 (a-h,o-z)
      real *8 xyztarg(3,1000),xyztarg0(3,1000)
      complex *16, allocatable :: slp(:,:)
      complex *16, allocatable :: dlp(:,:)
      complex *16, allocatable :: slp_ex(:,:)
      complex *16, allocatable :: dlp_ex(:,:)
      integer ipars(10)
      integer, allocatable :: iprint(:,:)

      real *8 dpars(5)
      complex * 16 zpars(3),ima

      external hslp,hdlp

      character type

      ima = dcmplx(0.0d0,1.0d0)
      
      done = 1
      pi = atan(done)*4


      call prini(6,13)
      call prinf('enter n*',n,0)
      read *, n


      norder = 8
      npols = norder*norder
      
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

      xyztarg(1,2) = xyztarg0(2,2) 
      xyztarg(2,2) = xyztarg0(3,2)
      xyztarg(3,2) = xyztarg0(1,2) - 1.0d0

      xyztarg(1,3) = xyztarg0(1,3) 
      xyztarg(2,3) = xyztarg0(3,3)
      xyztarg(3,3) = xyztarg0(2,3) - 1.0d0

      allocate(slp(npols,ntarg),slp_ex(npols,ntarg))
      allocate(dlp(npols,ntarg),dlp_ex(npols,ntarg))

      do i=1,ntarg
        do j=1,npols
          slp(j,i) = 0
          slp_ex(j,i) = 0
          dlp(j,i) = 0
          dlp_ex(j,i) = 0
        enddo
      enddo

      eps = 1.0d-6
      nqorder = 20
      intype = 2
      nquadmax = 5000

      type = 't'
      call cquadints_adap(eps,intype,norder,type,npols,ntarg,xyztarg,
     1       nquadmax,hslp,dpars,zpars,ipars,nqorder,slp)

      call cquadints_adap(eps,intype,norder,type,npols,ntarg,xyztarg,
     1       nquadmax,hdlp,dpars,zpars,ipars,nqorder,dlp)
      nn = norder + 1

      slp_ex(1,1) =  5.2060662105466030d0 + ima*3.8329924409148406d0
      slp_ex(2,1) = 0 
      slp_ex(nn,1) = 0.4103351431258097d0 + ima*0.0527184106855233d0

      dlp_ex(1,1) = 6.244445243530355d0 + ima*0.031068048605975d0
      dlp_ex(2,1) = 0
      dlp_ex(nn,1) = 0.6164533544673438d0 + ima*0.0002506281110130d0

      slp_ex(1,2) = 3.951622201719054d0+3.739912706597636d0*ima
      slp_ex(2,2) = 0.3017426968132834d0+0.0519535961159559d0*ima 
      slp_ex(nn,2) = -0.8521748769111288d0-0.1484084689967149d0*ima

      dlp_ex(1,2) = 5.691117319341798d0+0.322333025299182d0*ima 
      dlp_ex(2,2) = 0.4911891850156081d0+0.0026109845472635d0*ima 
      dlp_ex(nn,2) = -1.3928919779596355d0-0.0074591230001994d0*ima

      slp_ex(1,3) = 2.568571723511053d0+3.512622598195679d0*ima 
      slp_ex(2,3) = 0.5990712670972985d0+0.1501911845426824d0*ima
      slp_ex(nn,3) = -0.8427156061556782d0-0.2145004609068398d0*ima

      dlp_ex(1,3) = 4.712302914772059d0+0.621584380952863d0*ima
      dlp_ex(2,3) = 1.0154377749654707d0+0.0152618388506026d0*ima
      dlp_ex(nn,3) = -1.424266708898861d0-0.021799413779375d0*ima


      do i=1,ntarg
        print *, "itarg =", i
        print *, ""
        print *, ""

        err1s = abs(slp(1,i)-slp_ex(1,i))
        err1d = abs(dlp(1,i)-dlp_ex(1,i))

        err2s = abs(slp(2,i)-slp_ex(2,i))
        err2d = abs(dlp(2,i)-dlp_ex(2,i))

        err3s = abs(slp(nn,i)-slp_ex(nn,i))
        err3d = abs(dlp(nn,i)-dlp_ex(nn,i))

        print *, "0,0 int=",err1s,err1d
        print *, "1,0 int=",err2s,err2d
        print *, "0,1 int=",err3s,err3d
        print *, ""
        print *, ""
      enddo

      stop
      end
c
c
c
c
c

      subroutine lslp(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(3),y(3),dpars(*)
      complex *16 zpars(*)
      integer ipars(*)
      complex *16 f
      
      rr = (x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2
      f = 1/sqrt(rr)

      return
      end


      
      


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

