      implicit real *8 (a-h,o-z)
      real *8 dpars(12)
      integer iptr(9)
      integer, allocatable :: itree(:)
      real *8, allocatable :: fvals(:,:,:),centers(:,:),boxsize(:)
      real *8 rintl(0:200)
      complex *16 zk,zpars

      external fgauss3,fgauss1

      call prini(6,13)

c
c      initialize function parameters
c
      boxlen = 1.0d0
      zk = 2.0d0

      nd = 1
      do i=1,3
        dpars(i) = (hkrand(0)-0.5d0)*boxlen*0.01
      enddo
 1100 format(2x,e11.5,3(2x,i1),2(2x,i6),2x,i9,6(2x,e11.5))              

      dpars(4) = 0.05d0
cc      dpars(4) = 0.75d0
      call prin2('dpars=*',dpars,4)

      norder = 8 
      iptype = 0 
      npbox = norder*norder*norder
      iprec = 2
              
      print *, ""
      print *, ""
      print *, "========================="

      print *, norder,iptype,iprec

      eta = 0.0d0

      eps = 10.0d0**(-iprec*3)
      call cpu_time(t1)
C$    t1 = omp_get_wtime()      

      call vol_tree_mem(eps,zk,boxlen,norder,iptype,eta,
     1      fgauss1,nd,dpars,zpars,ipars,nlevels,nboxes,ltree,rintl)

      call cpu_time(t2)
C$    t2 = omp_get_wtime()     
      
      tmem = t2-t1
      call prin2('time taken in memory routine=*',t2-t1,1)



      call prinf('nboxes=*',nboxes,1)
      call prinf('nlevels=*',nlevels,1)
      call prin2('rintl=*',rintl,nlevels+1)

      allocate(fvals(nd,npbox,nboxes),centers(3,nboxes))
      allocate(boxsize(0:nlevels),itree(ltree))

      call cpu_time(t1)
C$    t1 = omp_get_wtime()
      
      call vol_tree_build(eps,zk,boxlen,norder,iptype,eta,
     1         fgauss1,nd,dpars,zpars,ipars,nlevels,nboxes,ltree0,
     2         rintl,itree,iptr,fvals,centers,boxsize)
      call cpu_time(t2)
C$    t2 = omp_get_wtime()      
              
      tgen  =t2-t1
      ttot = tgen+tmem
      call prin2('time taken to build tree=*',t2-t1,1)
      call prin2('speed in points per sec=*',
     1        (nboxes*norder**3+0.0d0)/(t2-t1+tmem),1)
      n = nboxes*npbox
      smem = (n+0.0d0)/tmem
      sgen = (n+0.0d0)/tgen
      stot = (n+0.0d0)/ttot
      write(*,1100) dpars(4),norder,iptype,iprec,nboxes,npbox,
     1         n,tmem,tgen,ttot,smem,sgen,stot

      

      stop
      end
c
c
c
c 
      subroutine fgauss3(nd,xyz,dpars,zpars,ipars,f)
c
c       compute three gaussians, their
c       centers are given in dpars(1:9), and their 
c       variances in dpars(10:12)
c
      implicit real *8 (a-h,o-z)
      integer nd,ipars
      complex *16 zpars
      real *8 dpars(12),f(3),xyz(3)

      if(nd.ne.3) then
        call prinf('someone changed nd - crashing*',i,0)
        stop
      endif

      do i=1,3
        idp = (i-1)*3 
        rr = (xyz(1) - dpars(idp+1))**2 + 
     1     (xyz(2) - dpars(idp+2))**2 + 
     1     (xyz(3) - dpars(idp+3))**2

        sigma = (dpars(9+i)**2)*2
        f(i) = exp(-rr/sigma)
      enddo
      
      

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

