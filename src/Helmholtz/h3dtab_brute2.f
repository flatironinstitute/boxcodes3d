      subroutine h3dtabp_ref_brute(ndeg,zk,tab,ldtab)

c
c
c     generate the Helmholtz potential table at the reference
c     points using adaptive integration
c
c
c     input
c
c     ndeg - integer, highest degree of the basis polynomials
c                    (measured in total degree, e.g. x^0y^1z^2
c                     has total degree 3)
c     zk - complex*16, helmholtz parameter
c     ldtab - integer, leading dimension of the output table
c
c     output
c      
c     tab - complex *16 array (ldtab,*), tab(i,j) is the integral
c     of the j-th tensor polynomial (in the ordering specified
c     in legetens.f) against the scaled green's function
c     exp(ikr)/r at the i-th reference target point (see tensrefpts3d)
     
      implicit real *8 (a-h,o-z)
      integer ndeg,ldtab
      complex *16 tab(ldtab,*), zk

c       local      
      complex *16, allocatable :: tab_t(:,:),tab_tmp(:,:)
      real *8 xyzc(3),bs,xq(100),w(100)
      real *8, allocatable :: xyztarg(:,:)
      real *8, allocatable :: xyztmp(:,:)
      integer, allocatable :: iindtmp(:)
      character ptype

      external h3d_vslp


      n = ndeg + 1
      norder = n

      call legetens_npol_3d(ndeg,'t',npols)

      bs = 2.0d0
      xyzc(1) = -1
      xyzc(2) = -1
      xyzc(3) = -1

      itype = 0
      call legeexps(itype,n,xq,u,v,w)

      do i=1,norder
        xq(i) = xq(i) + 1
      enddo


      npt = n**3
      ntarg = 10*npt
      allocate(xyztarg(3,ntarg))
      allocate(xyztmp(3,ntarg),iindtmp(ntarg))

      istart1 = 4*npt+1
      istart2 = 7*npt+1
      call tensrefpts3d(xq,norder,bs,xyzc,xyztarg,xyztarg(1,istart1),
     1  xyztarg(1,istart2))

      call prin2('xyztarg=*',xyztarg,3*npt)
      
      allocate(tab_t(npols,ntarg),tab_tmp(npols,ntarg))


      ntmp = npt
      do i=1,npt
        xyztmp(1,i) = xyztarg(1,i)
        xyztmp(2,i) = xyztarg(2,i)
        xyztmp(3,i) = xyztarg(3,i)
      enddo


      nbatches = 10

      call prinf('nbatches=*',nbatches,1)

      nttpcore = ceiling((ntmp+0.0d0)/nbatches)

      eps = 1.0d-11
      ncmax = 30000
      nqorder = 11

      call prin2('zk=*',zk,2)
      
      call cpu_time(t1)
C$     t1 = omp_get_wtime()
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,istart,iend,ntt)
C$OMP$SCHEDULE(DYNAMIC)
      do i=1,nbatches
        istart = (i-1)*nttpcore+1
        iend = min(i*nttpcore,ntmp)
        ntt = iend-istart+1
        print *,i,ntt,istart,iend
        call ccubeints_adap(eps,norder,'t',npols,ntt,xyztmp(1,istart),
     1        ncmax,h3d_vslp,dpars,zk,ipars,nqorder,tab_t(1,istart))
        print *, "Finished batch# ", i
      enddo
C$OMP END PARALLEL DO

      call cpu_time(t2)
C$      t2 = omp_get_wtime()     
      
      call prin2('self time=*',t2-t1,1)
      call prinf('finished self*',i,0)
      
c
c   gather near targets
c     

      rfac = 1.25d0
      ntmp = 0
      do i=npt+1,ntarg
        iinc = 0
        if(abs(xyztarg(1,i)).ge.1.25d0.or.
     1      abs(xyztarg(2,i)).ge.1.25d0.or.
     2      abs(xyztarg(3,i)).ge.1.25d0) iinc = 1
        
        if(iinc.eq.0) then
          ntmp = ntmp + 1
          xyztmp(1,ntmp) = xyztarg(1,i)
          xyztmp(2,ntmp) = xyztarg(2,i)
          xyztmp(3,ntmp) = xyztarg(3,i)
          iindtmp(ntmp) = i
        endif
      enddo
      
      call prinf('nnear=*',ntmp,1)

      nttpcore = ceiling((ntmp+0.0d0)/nbatches)

      call cpu_time(t1)
C$     t1 = omp_get_wtime()
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,istart,iend,ntt)
C$OMP$SCHEDULE(DYNAMIC)
      do i=1,nbatches
        istart = (i-1)*nttpcore+1
        iend = min(i*nttpcore,ntmp)
        ntt = iend-istart+1
        print *, i,istart,iend,ntt
        call ccubeints_adap(eps,norder,'t',npols,ntt,xyztmp(1,istart),
     1        ncmax,h3d_vslp,dpars,zk,ipars,nqorder,tab_tmp(1,istart))
        print *, "finished batch #",i
      enddo
C$OMP END PARALLEL DO

      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      call prin2('near time=*',t2-t1,1)
      call prinf('finished near*',i,0)
      
c
c   resort near
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,itarg,j)
      do i=1,ntmp
        itarg = iindtmp(i)
        do j=1,npols
          tab_t(j,itarg) = tab_tmp(j,i)
        enddo
      enddo
C$OMP END PARALLEL DO      

c
c   gather near targets
c     

      rfac = 1.25d0
      ntmp = 0
      do i=npt+1,ntarg
        iinc = 0
        if(abs(xyztarg(1,i)).ge.1.25d0.or.
     1      abs(xyztarg(2,i)).ge.1.25d0.or.
     2      abs(xyztarg(3,i)).ge.1.25d0) iinc = 1
        
        if(iinc.eq.1) then
          ntmp = ntmp + 1
          xyztmp(1,ntmp) = xyztarg(1,i)
          xyztmp(2,ntmp) = xyztarg(2,i)
          xyztmp(3,ntmp) = xyztarg(3,i)
          iindtmp(ntmp) = i
        endif
      enddo
      
      call prinf('nfar=*',ntmp,1)

      nttpcore = ceiling((ntmp+0.0d0)/nbatches)

      
      call cpu_time(t1)
C$     t1 = omp_get_wtime()
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,istart,iend,ntt)
C$OMP$SCHEDULE(DYNAMIC)
      do i=1,nbatches
        istart = (i-1)*nttpcore+1
        iend = min(i*nttpcore,ntmp)
        ntt = iend-istart+1
        print *, i,istart,iend,ntt
        call ccubeints_adap(eps,norder,'t',npols,ntt,xyztmp(1,istart),
     1        ncmax,h3d_vslp,dpars,zk,ipars,nqorder,tab_tmp(1,istart))
      enddo
C$OMP END PARALLEL DO

      call cpu_time(t2)
C$      t2 = omp_get_wtime()  
      call prin2('far time=*',t2-t1,1)
      call prinf('finished far*',i,0)
      
c
c   resort near
c     
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,itarg)
      do i=1,ntmp
        itarg = iindtmp(i)
        do j=1,npols
          tab_t(j,itarg) = tab_tmp(j,i)
        enddo
      enddo
C$OMP END PARALLEL DO      


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
      subroutine h3dtabp_ref_brute_new(ndeg,zk,tol,tab,ldtab)

c
c
c     generate the Helmholtz potential table at the reference
c     points using adaptive integration
c
c
c     input
c
c     ndeg - integer, highest degree of the basis polynomials
c                    (measured in total degree, e.g. x^0y^1z^2
c                     has total degree 3)
c     zk - complex*16, helmholtz parameter
c     tol - input tolerance
c     ldtab - integer, leading dimension of the output table
c
c     output
c      
c     tab - complex *16 array (ldtab,*), tab(i,j) is the integral
c     of the j-th tensor polynomial (in the ordering specified
c     in legetens.f) against the scaled green's function
c     exp(ikr)/r at the i-th reference target point (see tensrefpts3d)
     
      implicit real *8 (a-h,o-z)
      integer ndeg,ldtab
      complex *16 tab(ldtab,*), zk

c       local      
      complex *16, allocatable :: tab_t(:,:),tab_tmp(:,:)
      real *8 xyzc(3),bs,xq(100),w(100)
      real *8, allocatable :: xyztarg(:,:)
      real *8, allocatable :: xyztmp(:,:)
      integer, allocatable :: iindtmp(:)
      character ptype

      external h3d_vslp


      n = ndeg + 1
      norder = n

      call legetens_npol_3d(ndeg,'t',npols)

      bs = 2.0d0
      xyzc(1) = -1
      xyzc(2) = -1
      xyzc(3) = -1

      itype = 0
      call legeexps(itype,n,xq,u,v,w)

      do i=1,norder
        xq(i) = xq(i) + 1
      enddo


      npt = n**3
      ntarg = 10*npt
      allocate(xyztarg(3,ntarg))
      allocate(xyztmp(3,ntarg),iindtmp(ntarg))

      istart1 = 4*npt+1
      istart2 = 7*npt+1
      call tensrefpts3d(xq,norder,bs,xyzc,xyztarg,xyztarg(1,istart1),
     1  xyztarg(1,istart2))

cc      call prin2('xyztarg=*',xyztarg,3*npt)
      
      allocate(tab_t(npols,ntarg))
c
c       vector initialization
c
      tab_t = 0


      eps = tol
      ncmax = 30000
      nqorder = 11

      call prin2('zk=*',zk,2)

      ntt = 1
      
      call cpu_time(t1)
C$     t1 = omp_get_wtime()
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
C$OMP$SCHEDULE(DYNAMIC)
      do i=npt+1,ntarg
        call ccubeints_adap(eps,norder,'t',npols,ntt,xyztarg(1,i),
     1        ncmax,h3d_vslp,dpars,zk,ipars,nqorder,tab_t(1,i))
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
c
c
c
c
