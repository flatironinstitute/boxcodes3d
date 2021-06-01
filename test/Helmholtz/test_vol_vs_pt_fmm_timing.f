      implicit real *8 (a-h,o-z)
      character *2 arg_comm
      complex *16 zk

      call prini(6,13)

      zk = 6.28d0
      call get_command_argument(1,arg_comm)
      read(arg_comm,*) iprec

      call get_command_argument(2,arg_comm)
      read(arg_comm,*) norder

      call get_command_argument(3,arg_comm)
      read(arg_comm,*) icase

      call prinf('starting computation*',i,0)

      if(iprec.eq.3.and.norder.eq.4) goto 1111
      if(iprec.eq.4.and.norder.eq.4) goto 1111
      if(iprec.eq.4.and.norder.eq.6) goto 1111
      if(iprec.eq.0.and.norder.eq.8) goto 1111
      if(iprec.eq.1.and.norder.eq.8) goto 1111
      if(iprec.eq.0.and.norder.eq.12) goto 1111
      if(iprec.eq.1.and.norder.eq.12) goto 1111

      errvfmm = 0
      errpfmm = 0
      
      
      call get_fmm_timing(zk,iprec,icase,norder,n,tvtree,tvpre,tvol,
     1   tfmm,errvfmm,errpfmm)

      call prinsum(zk,iprec,icase,norder,n,tvtree,tvpre,tvol,tfmm,
     1  errvfmm,errpfmm)

 1111 continue     
      

      stop
      end
      
      
      
      subroutine get_fmm_timing(zk,iprec,icase,norder,nnn,tvtree,tvpre,
     1   tvol,tfmm,errvfmm,errpfmm)
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
      real *8, allocatable :: sources(:,:)
      complex *16, allocatable :: charges(:)
      complex *16, allocatable :: potfmm(:),potfmmex(:)
      complex *16 zk,zpars

      complex *16, allocatable :: pot(:,:),potex(:,:)
      complex *16 ima,zz,ztmp

      real *8 alpha,beta

      procedure (), pointer :: fker

      character *1, type
      data ima/(0.0d0,1.0d0)/

      external fgaussn,fgauss1,fgaussmulti
      logical flag


      done = 1
      pi = atan(done)*4

      boxlen = 1.0d0

      if(icase.eq.1.or.icase.eq.2) then

        nd = 2
        dpars(1) = 0.01d0
        dpars(2) = 0.0d0
        dpars(3) = 0.0d0
        do i=1,3
          dpars(3+i) = dpars(i)
        enddo

        if(icase.eq.1) rsig = 1.0d0/13.0d0/2.0d0
        if(icase.eq.2) rsig = 1.0d0/13.0d0/4.0d0

        dpars(7) = rsig 
        dpars(8) = rsig

        fker => fgaussn
      endif

      if(icase.eq.3) then
        ngau = 25
        ipars = ngau
        do i=1,ngau
          rtmp = 0.125d0 + 0.25d0*hkrand(0)
          phi = hkrand(0)*2*pi
          thet = hkrand(0)*pi
          dpars(4*i-3) = rtmp*sin(thet)*cos(phi)
          dpars(4*i-2) = rtmp*sin(thet)*sin(phi)
          dpars(4*i-1) = rtmp*cos(thet)
          dpars(4*i) = 0.01d0 + hkrand(0)*0.04d0
          dpars(4*i) = 1.0d0/13.0d0
        enddo
        fker => fgaussmulti
        nd = 2
      endif

      iptype = 0
      eta = 2

      zkeff = zk*boxlen
      npbox = norder*norder*norder

      if(iprec.eq.0) eps = 0.51d-2
      if(iprec.eq.1) eps = 0.51d-3
      if(iprec.eq.2) eps = 0.51d-6
      if(iprec.eq.3) eps = 0.51d-9
      if(iprec.eq.4) eps = 0.51d-12

      call prin2('eps=*',eps,1)
      call prin2('zk=*',zk,2)
      call prinf('norder=*',norder,1)
      call prinf('nd=*',nd,1)

      

      call cpu_time(t1)
C$      t1 = omp_get_wtime()

      call vol_tree_mem(eps,zk,boxlen,norder,iptype,eta,
     1   fker,nd,dpars,zpars,ipars,nlevels,nboxes,ltree,rintl)

      call prinf('nboxes=*',nboxes,1)
      call prinf('nlevels=*',nlevels,1)


      allocate(fvals(nd,npbox,nboxes),centers(3,nboxes))
      allocate(boxsize(0:nlevels),itree(ltree))

      call vol_tree_build(eps,zk,boxlen,norder,iptype,eta,fker,nd,
     1  dpars,zpars,ipars,nlevels,nboxes,ltree,rintl,itree,iptr,fvals,
     2  centers,boxsize)
      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      

      tvtree = t2-t1
      
      if(icase.eq.1.or.icase.eq.2) then
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
      tvol = sum(timeinfo)

      call prin2('precomp time=*',tprecomp,3)
      call prin2('timeinfo=*',timeinfo,6)

      nlfbox = 0
      do ilevel=1,nlevels
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) nlfbox = nlfbox+1
        enddo
      enddo
      call prinf('nlfbox=*',nlfbox,1)
      call prin2('speed in pps=*',(npbox*nlfbox+0.0d0)/(tvol),1)

      allocate(xref(3,npbox),umat(npols,npbox),vmat(npbox,npols),
     1    wts(npbox))
      itype = 1
      call legetens_exps_3d(itype,norder,'t',xref,umat,1,vmat,1,wts)
      nsrc = npbox*nlfbox
      nnn = nsrc

      allocate(sources(3,nsrc),charges(nsrc),potfmm(nsrc))
      ipt = 1
      call prinf('npbox=*',npbox,1)
      call prinf('nlfbox=*',nlfbox,1)

      do ilevel=1,nlevels
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) then
            do j=1,npbox
              x = centers(1,ibox) + xref(1,j)*boxsize(ilevel)/2.0d0
              y = centers(2,ibox) + xref(2,j)*boxsize(ilevel)/2.0d0
              z = centers(3,ibox) + xref(3,j)*boxsize(ilevel)/2.0d0
              sources(1,ipt) = x
              sources(2,ipt) = y
              sources(3,ipt) = z
              charges(ipt) = fvals(1,j,ibox) + ima*fvals(2,j,ibox)
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

      if(icase.eq.1.or.icase.eq.2) then
        erra = 0.0d0
        ra = 0.0d0

        nboxtest = 10

        allocate(potex(npbox,nboxes))

        ss = dpars(7)
        c = exp(-ss**2*zk**2/2.0d0)*ss**3*(2.0d0*pi)**1.5d0
        itype = 0
        iboxtest = 0
        do ilevel=1,nlevels
          rsc = (boxsize(ilevel)/2.0d0)**3
          do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
            if(itree(iptr(4)+ibox-1).eq.0) then
              iboxtest = iboxtest+1
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

                erra = erra + abs(pot(j,ibox)-potex(j,ibox))**2*
     1             wts(j)*rsc
                ra = ra + abs(potex(j,ibox))**2*wts(j)*rsc
              enddo
              if(iboxtest.ge.nboxtest) goto 1111
            endif
          enddo
        enddo
 1111 continue        

        erra = sqrt(erra/ra)
        call prin2('erra vol fmm=*',erra,1)
        
        errvfmm = erra
       

c
c         test accuracy at 10 targets 
c
        ntest = 10
        allocate(potfmmex(ntest))
        potfmmex = 0
        thresh = boxlen*2.0d0**(-51)
        call h3ddirectcp(1,zk,sources,charges,nsrc,sources,ntest,
     1     potfmmex,thresh)

        erra = 0
        ra = 0
        do i=1,ntest
          ra = ra + abs(potfmmex(i))**2
          erra = erra + abs(potfmm(i)-potfmmex(i))**2
        enddo

        erra = sqrt(erra/ra)
        call prin2('erra pt fmm=*',erra,1)

        call prin2('potfmm=*',potfmm,20)
        call prin2('potfmmex=*',potfmmex,20)
        
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
      subroutine fgaussmulti(nd,xyz,dpars,zpars,ipars,f)
c
c       compute three gaussians, their
c       centers are given in dpars(1:3*nd), and their 
c       variances in dpars(3*nd+1:4*nd)
c
      implicit real *8 (a-h,o-z)
      integer nd,ipars
      complex *16 zpars
      real *8 dpars(*),f(2),xyz(3)


      ngau = ipars
      f(1) = 0
      do i=1,ngau
        idp = (i-1)*4 
        rr = (xyz(1) - dpars(idp+1))**2 + 
     1     (xyz(2) - dpars(idp+2))**2 + 
     1     (xyz(3) - dpars(idp+3))**2

        sigma = (dpars(idp+4)**2)*2
        f(1) = f(1)+exp(-rr/sigma)
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


      subroutine prinsum(zk,iprec,icase,norder,n,tvtree,tvpre,tvol,tfmm,
     1  errvfmm,errpfmm)
      implicit real *8 (a-h,o-z)
      complex *16 zk

      svtree = (n+0.0d0)/tvtree
      svpre = (n+0.0d0)/tvpre
      svol = (n+0.0d0)/tvol
      sfmm = (n+0.0d0)/tfmm
      call prin2(" *",i,0)
      call prin2(" *",i,0)
      call prin2("==============*",i,0)
      call prin2('zk=*',zk,2)
      call prinf('iprec=*',iprec,1)
      call prinf('norder=*',norder,1)
      call prinf('N = ',n,1)
      call prin2('volume tree generation speed=*',svtree,1)
      call prin2('volume precomputation speed=*',svpre,1)
      call prin2('volume fmm speed=*',svol,1)
      call prin2('point fmm speed=*',sfmm,1)

      if(icase.eq.1.or.icase.eq.2) then
        call prin2('volume fmm error=*',errvfmm,1)
        call prin2('pt fmm error=*',errpfmm,1)
      endif

      
      return
      end
