      implicit real *8 (a-h,o-z)
      complex *16 zk, im, zero, one
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /            
      
      call prini(6,13)
      
      zk = 2.5d0 + im*0.05d0
      zk = 2.5d0 + im*0.1d0
      zk = 1.6d0

      iseed = 123
c      call prinf('iseed*',iseed,1)
      a = hkrand(iseed)
      
      tol = 1.0d-9
      
      ndeg = 7
      n = ndeg + 1
      npt = n**3
c     point to test (should be in range 1:10*n**3)
      itest = 0*npt + 3

c     basis polynomial to test
      ix = (hkrand(0)*1d5)
      iy = (hkrand(0)*1d5)
      iz = (hkrand(0)*1d5)
      ndeg1 = ndeg
      ix = ix-ndeg1*(ix/ndeg1)
      ndeg1 = ndeg1-ix
      if (ndeg1 .gt. 0) then
         iy = iy-ndeg1*(iy/ndeg1)
         ndeg1 = ndeg1-iy
         if (ndeg1 .gt. 0) then
            iz = iz-ndeg1*(iz/ndeg1)
         else
            iz = 0
         endif
      else
         iy = 0
         iz = 0
      endif
      ix = ix + 1
      iy = iy + 1
      iz = iz + 1

      call prinf('ix is *',ix,1)
      call prinf('iy is *',iy,1)
      call prinf('iz is *',iz,1)
      
      call test1(zk,ndeg,itest,ix,iy,iz,tol)
      
      stop
      end
c
c      

      subroutine test1(zk,ndeg,itest,ix,iy,iz,tol)
      implicit real *8 (a-h,o-z)
      integer, allocatable :: ip2ind(:,:,:)
      real *8 xq(200), u, v, w, xyzc(3), xtest(3), tol
      real *8, allocatable :: xyz(:,:)
      complex *16, allocatable :: tab(:,:)
      complex *16 zk, im, zero, one, cv, ccv
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /            
      character type

      type = 't'
      call legetens_npol_3d(ndeg,type,npol3)

      call prin2('zk is *',zk,2)
      call prinf('ndeg is *',ndeg,1)      

      n = ndeg+1

      allocate(ip2ind(n,n,n))
      call legetens_pow2ind_3d(ndeg,type,ip2ind)
      
c
c     get reference points
c

      itype = 0
      call legeexps(itype,n,xq,u,v,w)

      do i = 1,n
         xq(i) = xq(i)+1.0d0
      enddo

      
      ntarg0 = 10*n**3
      allocate(xyz(3,ntarg0))

      xyzc(1) = -1.0d0
      xyzc(2) = -1.0d0
      xyzc(3) = -1.0d0
      bs = 2.0d0

      istart1 = 4*n**3+1
      istart2 = 7*n**3+1
      
      call tensrefpts3d(xq,n,bs,xyzc,xyz,xyz(1,istart1),
     1        xyz(1,istart2))

c
c     call table generator
c      

      allocate(tab(ntarg0,npol3))
      call cpu_time(t1)
c$    t1 = omp_get_wtime()      
      call h3dtabp_ref(ndeg,zk,tol,tab,ntarg0)
      call cpu_time(t2)
c$    t2 = omp_get_wtime()            

      call prin2('time to generate table *',t2-t1,1)
      
c
c     brute force one entry
c     

      xtest(1) = xyz(1,itest)
      xtest(2) = xyz(2,itest)
      xtest(3) = xyz(3,itest)

      ipol = ip2ind(ix,iy,iz)

      call mksurhelm3dp(xtest(1),xtest(2),xtest(3),ix,iy,iz,
     1     zk,n,cv,ifail)

      write(*,*) itest, ipol

      ccv = tab(itest,ipol)

c
c     compare
c

      err1 = abs(cv-ccv)

      call prinf('itest is *',itest,1)
      call prinf('ix is *',ix,1)
      call prinf('iy is *',iy,1)
      call prinf('iz is *',iz,1)
      call prinf('ifail is *',ifail,1)            
      call prin2_long('error is *',err1,1)

      return
      end
