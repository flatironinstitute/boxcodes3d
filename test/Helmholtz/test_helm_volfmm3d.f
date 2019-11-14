      implicit real *8 (a-h,o-z)
      real *8 dpars(1000)
      integer iptr(9)
      integer, allocatable :: itree(:)
      real *8, allocatable :: fvals(:,:,:),centers(:,:),boxsize(:)
      real *8, allocatable :: umat(:,:),vmat(:,:),xref(:,:),wts(:)
      complex *16 zk,zpars

      complex *16, allocatable :: pot(:,:),potex(:,:)

      real *8 alpha,beta

      character *1, type
      real *8, allocatable :: fcoefs(:,:,:)

      external fgaussn,fgauss1

      call prini(6,13)

c
c      initialize function parameters
c
      boxlen = 1.0d0

      nd = 2
      dpars(1) = 0.03d0
      dpars(2) = -0.17d0
      dpars(3) = 0.23d0
      do i=1,3
        dpars(3+i) = dpars(i)
      enddo

      dpars(7) = 3.0d0
      dpars(8) = 3.0d0

      zk = 2.1d0
      norder = 4
      iptype = 1
      eta = 1


      zkeff = zk*boxlen

      npbox = norder*norder*norder

      eps = 1.0d-5
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

      alpha = 1
      beta = 0
      do ibox=1,nboxes
        call dgemm('n','t',nd,npols,npbox,alpha,fvals(1,1,ibox),
     1     nd,umat,npols,beta,fcoefs(1,1,ibox),nd)

      enddo


      allocate(pot(npbox,nboxes))

      do i=1,nboxes
        do j=1,npbox
          pot(j,i) = 0
        enddo
      enddo



      
      call helmholtz_volume_fmm(eps,zk,nboxes,nlevels,ltree,itree,
     1   iptr,norder,npols,type,fcoefs,centers,boxsize,npbox,
     2   pot)
 1000 continue

c
c
c   write potential to file, error will be compared in matlab/python
c   since error function of complex argument is requried
c

 2623 format(4(2x,d22.16))
      do ibox=itree(2*nlevels+1),itree(2*nlevels+2)
        print *, ibox,npbox,boxsize(nlevels)
        print *, centers(1:3,ibox)
        do j=1,npbox
          if(ibox.eq.9) print *, xref(1,j),xref(2,j),xref(3,j)
          x = centers(1,ibox) + xref(1,j)*boxsize(nlevels)/2.0d0
          y = centers(2,ibox) + xref(2,j)*boxsize(nlevels)/2.0d0
          z = centers(3,ibox) + xref(3,j)*boxsize(nlevels)/2.0d0
          write(33,2623) x,y,z,real(pot(j,ibox))
        enddo
      enddo

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


