      implicit real *8 (a-h,o-z)
      real *8 dpars(1000)
      integer iptr(9)
      integer, allocatable :: itree(:)
      real *8, allocatable :: fvals(:,:,:),centers(:,:),boxsize(:)
      real *8, allocatable :: umat(:,:),vmat(:,:),xref(:,:),wts(:)
      real *8 xyztmp(3)
      complex *16 zk,zpars

      complex *16, allocatable :: pot(:,:),potex(:,:)
      complex *16 ima,zz

      real *8 alpha,beta

      character *1, type
      real *8, allocatable :: fcoefs(:,:,:)
      data ima/(0.0d0,1.0d0)/

      external fgaussn,fgauss1
      logical flag

      call prini(6,13)

      done = 1
      pi = atan(done)*4

c
c      initialize function parameters
c
      boxlen = 1.0d0

      nd = 2
      dpars(1) = 0.01d0
      dpars(2) = 0.0d0
      dpars(3) = 0.0d0
      do i=1,3
        dpars(3+i) = dpars(i)
      enddo

      dpars(7) = 1.0d0/13.0d0
      dpars(8) = 1.0d0/13.0d0

      zk = 2.1d0
      norder = 8
      iptype = 1
      eta = 2.0d0


      zkeff = zk*boxlen

      npbox = norder*norder*norder

      eps = 1.0d-6
      call cpu_time(t1)
C$      t1 = omp_get_wtime()


      call vol_tree_mem(eps,zk,boxlen,norder,iptype,eta,
     1   fgaussn,nd,dpars,zpars,ipars,nlevels,nboxes,ltree)

      call prinf('nboxes=*',nboxes,1)
      call prinf('nlevels=*',nlevels,1)


      allocate(fvals(nd,npbox,nboxes),centers(3,nboxes))
      allocate(boxsize(0:nlevels),itree(ltree))

      call vol_tree_build(eps,zk,boxlen,norder,iptype,eta,fgaussn,nd,
     1  dpars,zpars,ipars,nlevels,nboxes,ltree0,itree,iptr,fvals,
     2  centers,boxsize)
      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      

      do i=1,nboxes
        do j=1,npbox
          fvals(2,j,i) = 0
        enddo
      enddo
      call prin2('time taken to build tree=*',t2-t1,1)
      call prin2('speed in points per sec=*',
     1   (nboxes*norder**3+0.0d0)/(t2-t1),1)

      eps = 1.0d-3


c
c
c       convert values to coefs
c
      
      npols = norder*(norder+1)*(norder+2)/6

      allocate(xref(3,npbox),umat(npols,npbox),vmat(npbox,npols))
      allocate(wts(npbox))

      allocate(fcoefs(nd,npols,nboxes))

      type = 'T'
      itype = 2
      call legetens_exps_3d(itype,norder,type,xref,umat,npols,vmat,
     1   npbox,wts)
      call prinf('norder=*',norder,1)
      call prinf('npbox=*',npbox,1)



cc      print *,npols,npbox,nd

      alpha = 1.0d0
      beta = 0
      do ibox=1,nboxes
        call dgemm('n','t',nd,npols,npbox,alpha,fvals(1,1,ibox),
     1     nd,umat,npols,beta,fcoefs(1,1,ibox),nd)

c
c       test fcoefs values at a few random points in the box
c  
        xyztmp(1) = centers(1,ibox) + (hkrand(0)-0.5d0)*0.125d0
        xyztmp(2) = centers(1,ibox) + (hkrand(0)-0.5d0)*0.125d0
        xyztmp(3) = centers(1,ibox) + (hkrand(0)-0.5d0)*0.125d0

        x = (xyztmp(1) - centers(1,ibox))*1.0d0 
      enddo


      


      allocate(pot(npbox,nboxes))

      do i=1,nboxes
        do j=1,npbox
          pot(j,i) = 0
        enddo
      enddo


      call cpu_time(t1) 
      call helmholtz_volume_fmm(eps,zk,nboxes,nlevels,ltree,itree,
     1   iptr,norder,npols,type,fcoefs,centers,boxsize,npbox,
     2   pot)
      call cpu_time(t2) 
      call prin2('time taken in fmm=*',t2-t1,1)

      nlfbox = itree(2*nlevels+2)-itree(2*nlevels+1)+1

      print *, npbox,nlfbox
      call prin2('speed in pps=*',(npbox*nlfbox+0.0d0)/(t2-t1),1)
 1000 continue

    

c
c
c   write potential to file, error will be compared in matlab/python
c   since error function of complex argument is requried
c

 2623 format(6(2x,e11.5))
 2625 format(2(2x,e11.5))

      erra = 0.0d0
      ra = 0.0d0

      allocate(potex(npbox,nboxes))

      ss = dpars(7)
      c = exp(-ss**2*zk**2/2.0d0)*ss**3*(2.0d0*pi)**1.5d0

      do ibox=itree(2*nlevels+1),itree(2*nlevels+2)
        do j=1,npbox
          x = centers(1,ibox) + xref(1,j)*boxsize(nlevels)/2.0d0
          y = centers(2,ibox) + xref(2,j)*boxsize(nlevels)/2.0d0
          z = centers(3,ibox) + xref(3,j)*boxsize(nlevels)/2.0d0

          dx = x - dpars(1)
          dy = y - dpars(2)
          dz = z - dpars(3)

          r = sqrt(dx**2 + dy**2 + dz**2)

          zz = ima*(ss*ima*zk/sqrt(2.0d0)-r/sqrt(2.0d0)/ss)
          uu = 0
          vv = 0
          xx = real(zz)
          yy = imag(zz)
          call wofz(xx,yy,uu,vv,flag)


          zz = (1.0d0-(uu + ima*vv)*exp(zz**2))*exp(-ima*zk*r)
          potex(j,ibox) = c/r*(-real(zz)+ima*sin(zk*r))

          rr1 = imag(potex(j,ibox))
          rr2 = imag(pot(j,ibox))


          erra = erra + abs(pot(j,ibox)-potex(j,ibox))**2
cc          erra = erra + abs(rr1-rr2)**2
          ra = ra + abs(potex(j,ibox))**2
          write(33,2623) x,y,z,rr1,rr2,rr1/rr2
          write(34,2625) real(pot(j,ibox)),imag(pot(j,ibox))
          write(35,2625) real(potex(j,ibox)),imag(potex(j,ibox))
        enddo
      enddo

      erra = sqrt(erra/ra)
      call prin2('erra=*',erra,1)

      stop
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


