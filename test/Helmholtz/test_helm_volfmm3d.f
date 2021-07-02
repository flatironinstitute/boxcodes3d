      implicit real *8 (a-h,o-z)
      real *8 dpars(1000)
      integer iptr(9)
      integer, allocatable :: itree(:)
      real *8, allocatable :: fvals(:,:,:),centers(:,:),boxsize(:)
      real *8, allocatable :: umat(:,:),vmat(:,:),xref(:,:),wts(:)
      real *8 xyztmp(3),rintl(0:200)
      real *8 timeinfo(6),tprecomp(3)
      complex *16 zk,zpars

      complex *16, allocatable :: pot(:,:),potex(:,:)
      complex *16, allocatable :: potcoefs(:,:)
      complex *16 ima,zz,ztmp

      real *8 alpha,beta

      character *1, type
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

      rsig = 1.0d0/13.0d0
cc      rsig = 0.005d0

      dpars(7) = rsig 
      dpars(8) = rsig

     

      zk = 2.0d0
      norder = 8
      iptype = 0
      eta = 2


      zkeff = zk*boxlen

      npbox = norder*norder*norder

      eps = 1.0d-3
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


      allocate(pot(npbox,nboxes))
      allocate(potcoefs(npols,nboxes))

      do i=1,nboxes
        do j=1,npbox
          pot(j,i) = 0
        enddo

        do j=1,npols
          potcoefs(j,i) = 0

        enddo
      enddo

      type = 'T'
      
      call cpu_time(t1) 
C$     t1 = omp_get_wtime()      
      call helmholtz_volume_fmm(eps,zk,nboxes,nlevels,ltree,itree,
     1   iptr,norder,npols,type,fvals,centers,boxsize,npbox,
     2   pot,potcoefs,timeinfo,tprecomp)
      call cpu_time(t2) 
C$     t2 = omp_get_wtime()      
      call prin2('time taken in fmm=*',t2-t1,1)
      nlfbox = 0
      do ilevel=1,nlevels
        do ibox=itree(2*ilevel+1),itree(2*ilevel+2)
          if(itree(iptr(4)+ibox-1).eq.0) nlfbox = nlfbox+1
        enddo
      enddo
      call prinf('nlfbox=*',nlfbox,1)
      call prin2('speed in pps=*',(npbox*nlfbox+0.0d0)/(t2-t1),1)

    

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
      print *, "ss=",ss


      itype = 0
      allocate(xref(3,npbox),umat(npols,npbox),vmat(npbox,npols),
     1    wts(npbox))
      call legetens_exps_3d(itype,norder,'t',xref,umat,1,vmat,1,wts)
      do ilevel=1,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,j,x,y,z,dx,dy,dz,r)
C$OMP$PRIVATE(zz,nn,ztmp,rr1,rr2) REDUCTION(+:erra,ra)
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
              

              rr1 = imag(potex(j,ibox))
              rr2 = imag(pot(j,ibox))

              erra = erra + abs(pot(j,ibox)-potex(j,ibox))**2
              ra = ra + abs(potex(j,ibox))**2
            enddo
          endif
        enddo
C$OMP END PARALLEL DO        
      enddo


      erra = sqrt(erra/ra)
      call prin2('erra=*',erra,1)
      call prin2('ra=*',ra,1)

      i1 = 0
      ntests = 1
      if(erra.lt.eps) i1 = 1


      nsuccess = i1

      open(unit=33,file='../../print_testres.txt',access='append')
      write(33,'(a,i1,a,i1,a)') 'Successfully completed ',nsuccess,
     1  ' out of ',ntests,' in h3d volume fmm testing suite'
      write(*,'(a,i1,a,i1,a)') 'Successfully completed ',nsuccess,
     1  ' out of ',ntests,' in h3d volume fmm testing suite'
      close(33)
      

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


