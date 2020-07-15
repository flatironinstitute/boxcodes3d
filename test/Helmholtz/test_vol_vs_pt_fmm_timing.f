      implicit real *8 (a-h,o-z)
      complex *16 zk

      zk = 1.2d0
      iprec = 2
      norder = 4
      
      
      call get_fmm_timing(zk,iprec,icase,norder,tvtree,tvpre,tvol,
     1   tfmm,errvfmm,errpfmm)

      

      stop
      end
      
      
      
      subroutine get_fmm_timing(zk,iprec,icase,norder,tvtree,tvpre,tvol,
     1   tfmm,errvfmm,errpfmm)
      implicit real *8 (a-h,o-z)
      real *8 dpars(1000)
      integer iptr(9)
      integer, allocatable :: itree(:)
      real *8, allocatable :: fvals(:,:,:),centers(:,:),boxsize(:)
      real *8, allocatable :: umat(:,:),vmat(:,:),xref(:,:),wts(:)
      real *8 xyztmp(3),rintl(0:200)
      real *8 timeinfo(6),tprecomp(3)

c
c       point fmm variables
c
      real *, allocatable :: sources(:,:)
      complex *16, allocatable :: charges(:)
      complex *16, allocatable :: potfmm(:),potfmmex(:)
      complex *16 zk,zpars

      complex *16, allocatable :: pot(:,:),potex(:,:)
      complex *16 ima,zz,ztmp

      real *8 alpha,beta

      procedure (), pointer :: fker

      character *1, type
      data ima/(0.0d0,1.0d0)/

      external fgaussn,fgauss1
      logical flag

      call prini(6,13)

      done = 1
      pi = atan(done)*4

      boxlen = 1.0d0

      if(icase.eq.1) then

        nd = 2
        dpars(1) = 0.01d0
        dpars(2) = 0.0d0
        dpars(3) = 0.0d0
        do i=1,3
          dpars(3+i) = dpars(i)
        enddo

        rsig = 1.0d0/13.0d0

        dpars(7) = rsig 
        dpars(8) = rsig

        fker => fgaussn
      endif

      iptype = 0
      eta = 3

      zkeff = zk*boxlen
      npbox = norder*norder*norder

      if(iprec.eq.0) eps = 0.51d-2
      if(iprec.eq.1) eps = 0.51d-3
      if(iprec.eq.2) eps = 0.51d-6
      if(iprec.eq.3) eps = 0.51d-9
      if(iprec.eq.4) eps = 0.51d-12

      call cpu_time(t1)
C$      t1 = omp_get_wtime()

      call vol_tree_mem(eps,zk,boxlen,norder,iptype,eta,
     1   fgaussn,nd,dpars,zpars,ipars,nlevels,nboxes,ltree,rintl)

      call prinf('nboxes=*',nboxes,1)
      call prinf('nlevels=*',nlevels,1)


      allocate(fvals(nd,npbox,nboxes),centers(3,nboxes))
      allocate(boxsize(0:nlevels),itree(ltree))

      call vol_tree_build(eps,zk,boxlen,norder,iptype,eta,fgaussn,nd,
     1  dpars,zpars,ipars,nlevels,nboxes,ltree,rintl,itree,iptr,fvals,
     2  centers,boxsize)
      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      

      tvtree = t2-t1
      
      if(icase.eq.1) then
        do i=1,nboxes
          do j=1,npbox
            fvals(2,j,i) = 0
          enddo
        enddo
      endif
      call prin2('time taken to build tree=*',t2-t1,1)
      call prin2('speed in points per sec=*',
     1   (nboxes*norder**3+0.0d0)/(t2-t1),1)
c
c
c       convert values to coefs
c
      
      npols = norder*(norder+1)*(norder+2)/6


      allocate(pot(npbox,nboxes))

      do i=1,nboxes
        do j=1,npbox
          pot(j,i) = 0
        enddo
      enddo

      type = 'T'
      
      call cpu_time(t1) 
C$     t1 = omp_get_wtime()      
      call helmholtz_volume_fmm(eps,zk,nboxes,nlevels,ltree,itree,
     1   iptr,norder,npols,type,fvals,centers,boxsize,npbox,
     2   pot,timeinfo,tprecomp)
      call cpu_time(t2) 
C$     t2 = omp_get_wtime()      
      call prin2('time taken in fmm=*',t2-t1,1)
      
      tvpre = sum(tprecomp)
      tvfmm = sum(timeinfo)

      call prin2('precomp time=*',tprecomp,3)
      call prin2('timeinfo=*',timeinfo,6)

      nlfbox = 0
      do ilevel=1,nlevels
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) nlfbox = nlfbox+1
        enddo
      enddo
      call prinf('nlfbox=*',nlfbox,1)
      call prin2('speed in pps=*',(npbox*nlfbox+0.0d0)/(tvfmm),1)

      allocate(xref(3,npbox),umat(npols,npbox),vmat(npbox,npols),
     1    wts(npbox))
      call legetens_exps_3d(itype,norder,'t',xref,umat,1,vmat,1,wts)
      nsrc = npbox*nlfbox

      allocate(sources(3,nsrc),charges(nsrc),potfmm(nsrc))
      ipt = 1
      do ilevels=1,nlevels
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
            do j=1,npbox
              x = centers(1,ibox) + xref(1,j)*boxsize(ilevel)/2.0d0
              y = centers(2,ibox) + xref(2,j)*boxsize(ilevel)/2.0d0
              z = centers(3,ibox) + xref(3,j)*boxsize(ilevel)/2.0d0
              sources(1,ipt) = x
              sources(2,ipt) = y
              sources(3,ipt) = z
              charges(ipt) = fval(1,j,ibox) + ima*fval(2,j,ibox)
              potfmm(ipt) = 0
              ipt = ipt + 1
            enddo
          endif
        enddo
      enddo

      call cpu_time(t1)
C$       t1 = omp_get_wtime()      

      call hfmm3d_s_c_p(eps,zk,nsrc,sources,charges,potfmm)
      call cpu_time(t2)
C$       t2 = omp_get_wtime()     
      tfmm = t2-t1


    

c
c
c   write potential to file, error will be compared in matlab/python
c   since error function of complex argument is requried
c

      if(icase.eq.1) then
        erra = 0.0d0
        ra = 0.0d0

        allocate(potex(npbox,nboxes))

        ss = dpars(7)
        c = exp(-ss**2*zk**2/2.0d0)*ss**3*(2.0d0*pi)**1.5d0
        itype = 0
        do ilevel=1,nlevels
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
              do j=1,npbox
                x = centers(1,ibox) + xref(1,j)*boxsize(ilevel)/2.0d0
                y = centers(2,ibox) + xref(2,j)*boxsize(ilevel)/2.0d0
                z = centers(3,ibox) + xref(3,j)*boxsize(ilevel)/2.0d0

                dx = x - dpars(1)
                dy = y - dpars(2)
                dz = z - dpars(3)

                r = sqrt(dx**2 + dy**2 + dz**2)

                zz = (ss*ima*zk - r/ss)/sqrt(2.0d0)
                nn = 6*abs(real(zz)*imag(zz))+20
                call zerrf(zz,ztmp,nn)

                zz = exp(-ima*zk*r)*ztmp
                potex(j,ibox) = c/r*(-real(zz)+ima*sin(zk*r))

                erra = erra + abs(pot(j,ibox)-potex(j,ibox))**2
                ra = ra + abs(potex(j,ibox))**2
              enddo
            endif
          enddo
        enddo

        erra = sqrt(erra/ra)
        call prin2('erra vol fmm=*',erra,1)
        
        errvfmm = erra
       

c
c         test accuracy at 10 targets 
c
        ntest = 10
        allocate(potfmmex(ntest))
        call h3ddirectcp(1,zk,sources,charges,ns,sources,ntest,
     1     potfmmex)

        erra = 0
        ra = 0
        do i=1,ntest
          ra = ra + abs(potfmmex(i))**2
          erra = erra + abs(potfmm(i)-potfmmex(i))**2
        enddo
        erra = sqrt(erra/ra)
        call prin2('erra pt fmm=*',erra,1)
        errpfmm = erra


      endif

      return
      end
c
c
c
c 
      subroutine fgaussn(nd,xyz,dpars,zpars,ipars,f)
c
c       compute three gaussians, their
c       centers are given in dpars(1:3*nd), and their 
c       variances in dpars(3*nd+1:4*nd)
c
      implicit real *8 (a-h,o-z)
      integer nd,ipars
      complex *16 zpars
      real *8 dpars(*),f(nd),xyz(3)


      do i=1,nd
        idp = (i-1)*3 
        rr = (xyz(1) - dpars(idp+1))**2 + 
     1     (xyz(2) - dpars(idp+2))**2 + 
     1     (xyz(3) - dpars(idp+3))**2

        sigma = (dpars(nd*3+i)**2)*2
        f(i) = exp(-rr/sigma)
      enddo
      f(2) = 0
      
      

      return
      end

c
c
c
c 
      subroutine fgauss1(nd,xyz,dpars,zpars,ipars,f)
c
c       compute a single gaussian, their
c       centers are given in dpars(1:3), and their 
c       variances in dpars(4)
c
      implicit real *8 (a-h,o-z)
      integer nd,ipars
      complex *16 zpars
      real *8 dpars(4),f,xyz(3)

      rr = (xyz(1) - dpars(1))**2 + 
     1     (xyz(2) - dpars(2))**2 + 
     1     (xyz(3) - dpars(3))**2

      sigma = (dpars(4)**2)*2
      f = exp(-rr/sigma)
      
      

      return
      end


