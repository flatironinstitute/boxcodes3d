      implicit real *8 (a-h,o-z)
      real *8 xyztarg(3,1000),errs(1000)
      integer, allocatable :: iind2p(:,:)
      complex *16, allocatable :: pot(:,:)
      complex *16, allocatable :: pot_ex(:,:)
      complex * 16 zk,ima
      real *8 x,y,z
      integer ifail

      external fker

      character type

      ima = dcmplx(0.0d0,1.0d0)
      
      done = 1
      pi = atan(done)*4


      call prini(6,13)
      call prinf('enter n*',n,0)
      read *, n


      norder = 3
      
      ntarg = 2


      zk = 1.1d0+ima*0.01d0

      xyztarg(1,1) = 0
      xyztarg(2,1) = 0.1d0
      xyztarg(3,1) = 1.019d0

      xyztarg(1,2) = 0.2d0
      xyztarg(2,2) = 0.1d0
      xyztarg(3,2) = -2.0d0/7.0d0

      type = 't'
      ndeg = norder-1
      call legetens_npol_3d(ndeg,type,npols)

      allocate(iind2p(3,npols))
      call legetens_ind2pow_3d(ndeg,type,iind2p)

      allocate(pot(npols,ntarg),pot_ex(npols,ntarg))
      do i=1,ntarg
        do j=1,npols
          pot(j,i) = 0
          pot_ex(j,i) = 0
        enddo
      enddo

      call cpu_time(t1)
C$      t1 = omp_get_wtime      

      do i=1,ntarg
        x = xyztarg(1,i)
        y = xyztarg(2,i)
        z = xyztarg(3,i)
        print *, "itarg=",i
        do j=1,npols
           ix = iind2p(1,j)+1
           iy = iind2p(2,j)+1
           iz = iind2p(3,j)+1
           print *, "jpol=",j
           call mksurhelm3dp(x,y,z,ix,iy,iz,zk,norder,pot_ex(j,i),ifail)
        enddo
      enddo
      call cpu_time(t2)
C$       t2 = omp_get_wtime()     
      call prin2('old adap quad time=*',t2-t1,1)

 1111 continue

      eps = 1.0d-10
      nqorder = 11
      ncubemax = 5000

      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call ccubeints_adap(eps,norder,type,npols,ntarg,xyztarg,
     1       ncubemax,fker,dpars,zk,ipars,nqorder,pot)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()   

      call prin2('time in ccubeints_adap=*',t2-t1,1)

      do i=1,ntarg
        print *, "itarg =", i
        print *, ""
        print *, ""
        do j=1,npols
          errs(j) = abs(pot(j,i)-pot_ex(j,i)) 
        enddo
        call prin2('errors=*',errs,npols)
      enddo
      print *, ""
      print *, ""
      print *, "Starting split8int adap" 

      do i=1,ntarg
        do j=1,npols
          pot(j,i) = 0
        enddo
      enddo
      call ccubeints_split8int_adap(eps,norder,type,npols,
     1   xyztarg(1,2),ntarg,xyztarg,
     1   ncubemax,fker,dpars,zk,ipars,nqorder,pot)

      do i=1,ntarg
        print *, "itarg =", i
        print *, ""
        print *, ""
        do j=1,npols
          errs(j) = abs(pot(j,i)-pot_ex(j,i)) 
        enddo
        call prin2('errors=*',errs,npols)
      enddo
      

      stop
      end
c
c
c
c
c


      subroutine fker(x,y,dpars,zk,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(3),y(3)
      complex *16 zk,ima,f
      data ima/(0.0d0,1.0d0)/
      
      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + (x(3)-y(3))**2)

      f = exp(ima*zk*rr)/rr

      return
      end

