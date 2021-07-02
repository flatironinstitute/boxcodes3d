c
c  In this file we test the adaptive integration
c  routines in src/Common/cubeintrouts.f
c
c  There are two subroutines in cubeintrouts which
c  are being tested:
c
c  ccubeints_adap: adaptive integration [-1,1]^3
c  ccubeints_split8int_adap: adaptive integration on [-1,1]^3,
c     performed by doing adaptive integration on 8 sub cubes
c     split at xyzsplit \in [-1,1]^3 specified by the user
c
c  The integrals have been precomputed for the helmholtz
c  volume potential with some dissipation at two points
c  outside but close to [-1,1]^3 and contained in the mathematica
c  notebook 3dadapinttest.nb in this folder
c
c  This test can be run by either running make test in the 
c  main directory of this repo or by
c  running make -f test_cubeintrouts.make in the current
c  directory
c
c  
c  
c
c


      implicit real *8 (a-h,o-z)
      real *8 xyztarg(3,1000),errs(1000),xyzsplit(3)
      integer, allocatable :: iind2p(:,:)
      complex *16, allocatable :: pot(:,:)
      complex *16, allocatable :: pot_ex(:,:)
      complex * 16 zk,ima
      real *8 x,y,z

      external h3d_vslp

      character type

      ima = dcmplx(0.0d0,1.0d0)
      
      done = 1
      pi = atan(done)*4


      call prini(6,13)


      norder = 3
      
      ntarg = 2


      zk = 1.1d0+ima*0.01d0

      xyztarg(1,1) = 0
      xyztarg(2,1) = 0.1d0
      xyztarg(3,1) = 1.019d0

      xyztarg(1,2) = 0.2d0
      xyztarg(2,2) = 0.1d0
      xyztarg(3,2) = -1.0d0-2.0d0/7.0d0

      xyzsplit(1:3) = xyztarg(1:3,2)
      xyzsplit(3) = xyzsplit(3) + 1.0d0

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


      pot_ex(1,1) = 1.962567844230420d0 + 5.666932640107789d0*ima
      pot_ex(2,1) = 0.0d0
      pot_ex(4,1) = 0.2393609226238699d0 + 0.08941244906536621d0*ima

      pot_ex(1,2) = 0.3704600411814699d0 + 4.872164206040963d0*ima 
      pot_ex(2,2) = 0.2969329492428068d0 + 0.1639961756200836d0*ima 
      pot_ex(4,2) = 0.1490652437621909d0 + 0.08200892377305815d0*ima 

      eps = 1.0d-10
      nqorder = 11
      ncubemax = 5000

      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call ccubeints_adap(eps,norder,type,npols,ntarg,xyztarg,
     1       ncubemax,h3d_vslp,dpars,zk,ipars,nqorder,pot)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()   

      call prin2('time in ccubeints_adap=*',t2-t1,1)

      do i=1,ntarg
        print *, "itarg =", i
        print *, ""
        print *, ""
        err1 = abs(pot(1,i)-pot_ex(1,i))
        err2 = abs(pot(2,i)-pot_ex(2,i))
        err3 = abs(pot(4,i)-pot_ex(4,i))

        print *, "errors=",err1,err2,err3
        print *, ""
        print *, ""
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
     1   xyzsplit,ntarg,xyztarg,
     1   ncubemax,h3d_vslp,dpars,zk,ipars,nqorder,pot)

      do i=1,ntarg
        print *, "itarg =", i
        print *, ""
        print *, ""
        err1 = abs(pot(1,i)-pot_ex(1,i))
        err2 = abs(pot(2,i)-pot_ex(2,i))
        err3 = abs(pot(4,i)-pot_ex(4,i))

        print *, "errors=",err1,err2,err3
        print *, ""
        print *, ""
      enddo

      

      return
      end
c
c
c
c
c
