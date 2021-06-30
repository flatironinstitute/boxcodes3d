c
c  This file contains the following user callable subroutine:
c 
c     h3dtabp_ref_brute - compute reference quadrature 
c       tables using adaptive integration. The self
c       quadrature is handled by splitting the cube into
c       8 cubes while the rest are handled using 
c       standard adaptive integration
c
c
c
c
c
      subroutine h3dtab_ref_brute(ndeg,tol,grid,npt,tab,ldtab,fker,
     1   dpars,zpars,ipars)

c
c
c
c     generate the quadrature table at the reference
c     points using adaptive integration where the kernel of integration
c     is given by fker
c
c     Calling sequence arguments for fker should be
c     
c     subroutine fker(x,y,dpars,zpars,ipars,f)
c
c     Note that this subroutine does not compute the integrals
c     for the self cube and doesn't do any intelligent split
c     between near and far
c
c
c     input
c
c     ndeg - integer, highest degree of the basis polynomials
c                    (measured in total degree, e.g. x^0y^1z^2
c                     has total degree 3)
c     ldtab - integer, leading dimension of the output table
c     fker - function handle for evaluating kernel
c     dpars - real parameters to be used by fker
c     zpars - complex parameters to be used by fker
c     ipars - integer parameters to be used by fker
c
c     output
c      
c     tab - complex *16 array (ldtab,*), tab(i,j) is the integral
c     of the j-th tensor polynomial (in the ordering specified
c     in legetens.f) against given kernel
c     at the i-th reference target point (see tensrefpts3d)
c     generate the Helmholtz potential table at the reference
c     points using adaptive integration
     
      implicit real *8 (a-h,o-z)
      integer ndeg,ldtab
      real *8 grid(3,npt)
      complex *16 tab(ldtab,*), zpars(*)
      real *8 dpars(*)
      integer ipars(*)
      

c       local      
      complex *16, allocatable :: tab_t(:,:),tab_tmp(:,:)
      real *8 xyzc(3),bs,xq(100),w(100)
      real *8, allocatable :: xyztarg(:,:)
      real *8, allocatable :: xyztmp(:,:)
      integer, allocatable :: iindtmp(:)
      character ptype

      external fker



      n = ndeg + 1
      norder = n

      call legetens_npol_3d(ndeg,'t',npols)

      bs = 2.0d0
      xyzc(1) = -1
      xyzc(2) = -1
      xyzc(3) = -1

      ntarg = 10*npt
      allocate(xyztarg(3,ntarg))
      allocate(xyztmp(3,ntarg),iindtmp(ntarg))

      istart1 = 4*npt+1
      istart2 = 7*npt+1
      call tensrefpts3d_grid(npt,grid,bs,xyzc,xyztarg,
     1  xyztarg(1,istart1),xyztarg(1,istart2))

cc      call prin2('xyztarg=*',xyztarg,3*npt)
      
      allocate(tab_t(npols,ntarg))
c
c       vector initialization
c
      tab_t = 0


      eps = tol 
      ncmax = 30000
      nqorder = 11

      call prin2('zk=*',zpars,2)

      ntt = 1
      
      call cpu_time(t1)
C$     t1 = omp_get_wtime()
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(DYNAMIC)
      do i=1,npt
        call ccubeints_split8int_adap(eps,norder,'t',npols,
     1        xyztarg(1,i),ntt,xyztarg(1,i),
     1        ncmax,fker,dpars,zpars,ipars,nqorder,tab_t(1,i))
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
C$OMP$SCHEDULE(DYNAMIC)
      do i=npt+1,ntarg
        call ccubeints_adap(eps,norder,'t',npols,ntt,xyztarg(1,i),
     1        ncmax,fker,dpars,zpars,ipars,nqorder,tab_t(1,i))
      enddo
C$OMP END PARALLEL DO

      call cpu_time(t2)
C$      t2 = omp_get_wtime()     
      
c
c    now transpose tab_t
c      
      do i=1,ntarg
        do j=1,npols
          tab(i,j) = tab_t(j,i) 
        enddo
      enddo

      return
      end



c
c
c


c      
c      
c
c
c
      subroutine h3dtabbox_ref_brute(ndeg,tol,x,npt,tab,npols,fker,
     1   dpars,zpars,ipars)

c
c
c
c     generate the quadrature table at the reference
c     points using adaptive integration where the kernel of integration
c     is given by fker
c
c     Calling sequence arguments for fker should be
c     
c     subroutine fker(x,y,dpars,zpars,ipars,f)
c
c     Note that this subroutine does not compute the integrals
c     for the self cube and doesn't do any intelligent split
c     between near and far
c
c
c     input
c
c     ndeg - integer, highest degree of the basis polynomials
c                    (measured in total degree, e.g. x^0y^1z^2
c                     has total degree 3)
c     ldtab - integer, leading dimension of the output table
c     fker - function handle for evaluating kernel
c     dpars - real parameters to be used by fker
c     zpars - complex parameters to be used by fker
c     ipars - integer parameters to be used by fker
c
c     output
c      
c     tab - complex *16 array (ldtab,*), tab(i,j) is the integral
c     of the j-th tensor polynomial (in the ordering specified
c     in legetens.f) against given kernel
c     at the i-th reference target point (see tensrefpts3d)
c     generate the Helmholtz potential table at the reference
c     points using adaptive integration
     
      implicit real *8 (a-h,o-z)
      integer ndeg,ldtab
      real *8 x(3,npt)
      complex *16 tab(npt,*), zpars(*)
      real *8 dpars(*)
      integer ipars(*)
      

c       local      
      complex *16, allocatable :: tab_t(:,:),tab_tmp(:,:)
      real *8 xyzc(3),bs,xq(100),w(100)
      character ptype

      external fker



      n = ndeg + 1
      norder = n
      
      allocate(tab_t(npols,npt))
c
c       vector initialization
c
      tab_t = 0


      eps = tol 
      ncmax = 30000
      nqorder = 11

      call prin2('zk=*',zpars,2)

      ntt = 1
      

      if(ifself.eq.1) then
        call cpu_time(t1)
C$       t1 = omp_get_wtime()
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(DYNAMIC)
        do i=1,npt
          call ccubeints_split8int_adap(eps,norder,'t',npols,
     1          x(1,i),ntt,x(1,i),
     1          ncmax,fker,dpars,zpars,ipars,nqorder,tab_t(1,i))
        enddo
C$OMP END PARALLEL DO

      else


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
C$OMP$SCHEDULE(DYNAMIC)
        do i=1,npt
          call ccubeints_adap(eps,norder,'t',npols,ntt,x(1,i),
     1          ncmax,fker,dpars,zpars,ipars,nqorder,tab_t(1,i))
        enddo
C$OMP END PARALLEL DO

        call cpu_time(t2)
C$        t2 = omp_get_wtime()     
      endif
c
c    now transpose tab_t
c      
      do i=1,npt
        do j=1,npols
          tab(i,j) = tab_t(j,i) 
        enddo
      enddo

      return
      end



c
c
c

c
c
c
c
      subroutine h3d_vslp(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(3),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + (x(3)-y(3))**2)

      f = exp(ima*zpars(1)*rr)/rr

      return
      end

c
c
      subroutine h3d_vslp_gradx(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(3),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,zk

      zk = zpars(1)

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + (x(3)-y(3))**2)

      f = -(y(1)-x(1))*(1.0d0-ima*zk*rr)*exp(ima*zpars(1)*rr)/rr**3

      return
      end

c      
c      
c
c
c
c
      subroutine h3d_vslp_grady(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(3),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,zk

      zk = zpars(1)

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + (x(3)-y(3))**2)

      f = -(y(2)-x(2))*(1.0d0-ima*zk*rr)*exp(ima*zpars(1)*rr)/rr**3

      return
      end
c
c
c
c
      subroutine h3d_vslp_gradz(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(3),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,zk

      zk = zpars(1)

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + (x(3)-y(3))**2)

      f = -(y(3)-x(3))*(1.0d0-ima*zk*rr)*exp(ima*zpars(1)*rr)/rr**3

      return
      end

c      
c      
c
c
c
c
