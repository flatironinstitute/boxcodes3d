      implicit real *8 (a-h,o-z)
      real *8 dpars(1000)
      integer iptr(9)
      integer, allocatable :: itree(:)
      real *8, allocatable :: fvals(:,:,:),centers(:,:),boxsize(:)
      real *8, allocatable :: umat(:,:),vmat(:,:),xref(:,:),wts(:)

      complex *16, allocatable :: qval(:,:),rhsval(:,:),uval(:,:)
      real *8 xyztmp(3),rintl(0:200)
      real *8 timeinfo(6),tprecomp(3)
      complex *16 zk,zpars

      complex *16, allocatable :: pot(:,:),soln(:,:),potcoefs(:,:)
      complex *16 ima,zz,ztmp

      real *8 alpha,beta
      real *8 errs(1000)

      character *1, ttype
      data ima/(0.0d0,1.0d0)/

      external fgauss3_ls
      logical flag

      call prini(6,13)

      done = 1
      pi = atan(done)*4

c
c      initialize function parameters
c
      boxlen = 1.0d0

      nd = 3
      dpars(1) = 0.01d0
      dpars(2) = 0.0d0
      dpars(3) = 0.0d0
      do i=1,3
        dpars(3+i) = dpars(i) -0.07d0*hkrand(0)
      enddo

      rsig = 1.0d0/13.0d0
cc      rsig = 0.005d0

      dpars(7) = rsig 
      dpars(8) = rsig

     

      zk = 2.0d0
      norder = 8
      iptype = 1
      eta = 3


      zkeff = zk*boxlen

      npbox = norder*norder*norder

      eps = 1.0d-5
      call cpu_time(t1)
C$      t1 = omp_get_wtime()


      call vol_tree_mem(eps,zk,boxlen,norder,iptype,eta,
     1   fgauss3_ls,nd,dpars,zk,ipars,nlevels,nboxes,ltree,rintl)

      call prinf('nboxes=*',nboxes,1)
      call prinf('nlevels=*',nlevels,1)


      allocate(fvals(nd,npbox,nboxes),centers(3,nboxes))
      allocate(boxsize(0:nlevels),itree(ltree))

      call prin2('dpars=*',dpars,8)
      call prin2('zk=*',zk,2)
      call prinf('nd=*',nd,1)

      call prin2('eps=*',eps,1)
      call prin2('eta=*',eta,1)
      call prinf('iptype=*',iptype,1)
      call prin2('rintl=*',rintl,nlevels+1)

      call vol_tree_build(eps,zk,boxlen,norder,iptype,eta,fgauss3_ls,
     1  nd,dpars,zk,ipars,nlevels,nboxes,ltree,rintl,itree,iptr,
     2  fvals,centers,boxsize)
      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      

      print *, "done building tree"
      allocate(qval(npbox,nboxes),uval(npbox,nboxes),
     1    rhsval(npbox,nboxes))
      do i=1,nboxes
        do j=1,npbox
          qval(j,i) = fvals(1,j,i)
          uval(j,i) = fvals(2,j,i)
          rhsval(j,i) = fvals(3,j,i)
        enddo
      enddo


      call prin2('time taken to build tree=*',t2-t1,1)
      call prin2('speed in points per sec=*',
     1   (nboxes*norder**3+0.0d0)/(t2-t1),1)




c
c
c       convert values to coefs
c
      
      npols = norder*(norder+1)*(norder+2)/6


      allocate(soln(npbox,nboxes),pot(npbox,nboxes))
      allocate(potcoefs(npols,nboxes))
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
     2    pot,potcoefs,timeinfo,tprecomp)
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
      subroutine fgauss3_ls(nd,xyz,dpars,zpars,ipars,f)
c
c       compute three gaussians, their
c       centers are given in dpars(1:6), and their 
c       variances in dpars(7:8), the third function
c       is (\Delta + k^2(1+q)) applied to the second
c        function where u is the first function
c
      implicit real *8 (a-h,o-z)
      integer nd,ipars
      complex *16 zpars
      real *8 dpars(*),f(3),xyz(3),f2lap,rr(2),sigma(2)


      do i=1,2
        idp = (i-1)*3 
        rr(i) = (xyz(1) - dpars(idp+1))**2 + 
     1     (xyz(2) - dpars(idp+2))**2 + 
     1     (xyz(3) - dpars(idp+3))**2

        sigma(i) = (dpars(6+i)**2)*2
        f(i) = exp(-rr(i)/sigma(i))
      enddo

      f2lap = -6.0d0*f(2)/sigma(2) + 4.0d0*rr(2)*f(2)/sigma(2)**2 
      f(3) = real(zpars)**2*(1.0d0+f(1))*f(2) + f2lap  
      
      

      return
      end

c
c
c
c 


