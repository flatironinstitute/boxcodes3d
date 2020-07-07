      implicit real *8 (a-h,o-z)
      real *8 dpars(100)
      integer iptr(9), ipars(1)
      integer, allocatable :: itree(:)
      real *8, allocatable :: fvals(:,:,:),centers(:,:),boxsize(:)
      real *8, allocatable :: umat(:,:),vmat(:,:),xref(:,:),wts(:)

      complex *16, allocatable :: qval(:,:),rhsval(:,:),uval(:,:)
      real *8 xyztmp(3),rintl(0:200)
      real *8 timeinfo(6),tprecomp(3)
      complex *16 zk,zpars(4)

      complex *16, allocatable :: pot(:,:),soln(:,:)
      complex *16 ima,zz,ztmp, zbeam(3)

      real *8 alpha,beta
      real *8 errs(1000)

      character *1, ttype
      data ima/(0.0d0,1.0d0)/

      external ftree_eaton_and_beam_ls
      logical flag

      call prini(6,13)

      done = 1
      pi = atan(done)*4

c
c      initialize problem parameters
c


c     size and wavelength

      boxlen = 1.0d0
      nwave = 6

      zk = 2*pi*nwave/boxlen

c     beam parameters

c     coordinate direction
      ibeam = 1
      zbeam(1) = -boxlen/2 - 0.01d0 + ima*0.5d0
c
      zbeam(2) = 0.27d0*boxlen
      zbeam(3) = 0.01d0

c     load into zpars, ipars

      ipars(1) = ibeam
      zpars(1) = zbeam(1)
      zpars(2) = zbeam(2)
      zpars(3) = zbeam(3)

      zpars(4) = zk


      
c     eaton lens parameters

c     support
      alpha = 0.45d0

c     q cut-off
      beta = 2.0d0

c     precompute eaton lens q evaluator structure
      ldpars = 100
      nlege = 40
      call eaton_pre(alpha,beta,dpars,ldpars,nlege)


      
c     discretization parameters
      
      norder = 4
      iptype = 1
      eta = 2


      zkeff = zk*boxlen

      npbox = norder*norder*norder

      eps = 1.0d-3
      call cpu_time(t1)
C$      t1 = omp_get_wtime()

      nd = 3
      call vol_tree_mem(eps,zk,boxlen,norder,iptype,eta,
     1     ftree_eaton_and_beam_ls,nd,dpars,zpars,ipars,nlevels,
     2     nboxes,ltree,rintl)

      call prinf('nboxes=*',nboxes,1)
      call prinf('nlevels=*',nlevels,1)

      allocate(fvals(nd,npbox,nboxes),centers(3,nboxes))
      allocate(boxsize(0:nlevels),itree(ltree))

      call prin2('dpars=*',dpars,2)
      call prin2('zk=*',zk,2)
      call prinf('nd=*',nd,1)

      call prin2('eps=*',eps,1)
      call prin2('eta=*',eta,1)
      call prinf('iptype=*',iptype,1)
      call prin2('rintl=*',rintl,nlevels+1)

      call vol_tree_build(eps,zk,boxlen,norder,iptype,eta,
     1     ftree_eaton_and_beam_ls,nd,dpars,zpars,ipars,nlevels,
     2     nboxes,ltree,rintl,itree,iptr,fvals,centers,boxsize)
      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      

      print *, "done building tree"
      allocate(qval(npbox,nboxes),uval(npbox,nboxes),
     1    rhsval(npbox,nboxes))
      do i=1,nboxes
        do j=1,npbox
          qval(j,i) = fvals(1,j,i)
          rhsval(j,i) = fvals(2,j,i) + ima*fvals(3,j,i)
        enddo
      enddo


      call prin2('time taken to build tree=*',t2-t1,1)
      call prin2('speed in points per sec=*',
     1     (nboxes*norder**3+0.0d0)/(t2-t1),1)


      stop




c
c
c       convert values to coefs
c
      
      npols = norder*(norder+1)*(norder+2)/6


      allocate(soln(npbox,nboxes),pot(npbox,nboxes))
      numit = 200
      niter = 0
      irep = 1
      eps_gmres = 1.0d-7

      ttype = 'T'
      call ls_solver_guru(eps,zk,nboxes,nlevels,ltree,itree,iptr,
     1   norder,npols,ttype,qval,centers,boxsize,npbox,
     2   rhsval,irep,eps_gmres,numit,niter,errs,rres,soln)

      call prinf('niter=*',niter,1)
      call prin2('errs=*',errs,niter)

      if(irep.eq.1) then
      call cpu_time(t1) 
C$     t1 = omp_get_wtime()      
        call helmholtz_volume_fmm(eps,zk,nboxes,nlevels,ltree,itree,
     1    iptr,norder,npols,ttype,soln,centers,boxsize,npbox,
     2    pot,timeinfo,tprecomp)
      call cpu_time(t2) 
C$     t2 = omp_get_wtime()      
      call prin2('time taken in fmm=*',t2-t1,1)
      endif

      if(irep.eq.2) then
        do ibox=1,nboxes
          do j=1,npbox
            pot(j,ibox) = soln(j,ibox)
          enddo
        enddo
      endif
      
      erra = 0
      ra = 0  
      do ibox=1,nboxes
        if(itree(iptr(4)+ibox-1).eq.0) then
          do j=1,npbox
            ra = ra + abs(uval(j,ibox))**2
            erra = erra + abs(uval(j,ibox)-pot(j,ibox))**2
          enddo
        endif
      enddo

      print *, "erra=*",erra
      print *, "ra=",ra
      erra = sqrt(erra/ra)

      call prin2('erra = *',erra,1)

      stop
      end
c
c
c
c 
      subroutine ftree_eaton_and_beam_ls(nd,xyz,dpars,zpars,ipars,f)
c
c     compute the value of q for the eaton lens (f(1)), the real and
c     imaginary parts of the incoming gaussian beam (f(2),f(3))
c
c     eaton lens parameters are in dpars
c        dpars(*) fast eaton evaluator structure
c
c     gaussian beam parameters are in zpars(1:3)
c        zpars(1:3) - complex center for beam
c        zpars(4) - zk, helmholtz wave number
c        ipars(1) - beam direction
c      
      implicit none
      integer nd,ipars(*)
      complex *16 zpars(*)
      integer ii
      real *8 dpars(*),f(*),xyz(3),r2,r,yi
      complex *16 zdiff(3), z2, z, zk, ima, zf
      data ima / (0.0d0,1.0d0) /

      r2 = xyz(1)**2 + xyz(2)**2 + xyz(3)**2
      r = sqrt(r2)

c      call eaton_post(r,dpars,f(1))

      f(1) = 1.0d0

      zdiff(1) = xyz(1) - zpars(1)
      zdiff(2) = xyz(2) - zpars(2) 
      zdiff(3) = xyz(3) - zpars(3)

      z2 = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
      z = sqrt(z2)
      
      zk = zpars(4)

      ii = ipars(1)
      yi = -real(ima*zpars(ii))
      
      zf = exp(ima*zk*z)/z*exp(-zk*yi)

      f(2) = real(zf)
      f(3) = -real(ima*zf)

      return
      end

c
c
c
c 


