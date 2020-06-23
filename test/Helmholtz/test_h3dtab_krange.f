      implicit real *8 (a-h,o-z)
      character *24 fname
      integer npol_test(3,100)
      complex *16 zk, im, zero, one
      complex *16 zkvals(10)

      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /            
      
      call prini(6,13)

      nzkvals = 9
      do i=1,nzkvals
        zkvals(i) = 0.5d0 + (i-1)*0.5d0
      enddo

      iseed = 123
      a = hkrand(iseed)
      
      tol = 1.0d-3
      toltest = 1.0d-2
      
      ndeg = 2
      ntest = ndeg+1
      do i=0,ndeg
        npol_test(1,i+1) = irand_int(i)
        npol_test(2,i+1) = irand_int(i-npol_test(1,i+1))
        npol_test(3,i+1) = irand_int(i-npol_test(1,i+1)-
     1     npol_test(2,i+1))
      enddo

      do i=0,ndeg
        npol_test(1,i+1) = npol_test(1,i+1)+1
        npol_test(2,i+1) = npol_test(2,i+1)+1
        npol_test(3,i+1) = npol_test(3,i+1)+1
      enddo

      ntest = 2
      ntest = 1
      call cpu_time(t1)

      do izk=1,2
        write(fname,'(a,i1,a)') 'zkrange-res/zk',izk,'.txt'
        call test1(zkvals(izk),ndeg,ntest,npol_test,tol,fname)
      enddo
      call cpu_time(t2)
      
      stop
      end
c
c      

      subroutine test1(zk,ndeg,ntest,npol_test,tol,fname)
      implicit real *8 (a-h,o-z)
      integer npol_test(3,ntest)
      integer, allocatable :: ip2ind(:,:,:)
      integer, allocatable :: ifail(:,:)
      complex *16, allocatable :: cv(:,:),ccv(:,:),ccv2(:,:)
      complex *16, allocatable :: absest(:,:)
      real *8 xq(200), u, v, w, xyzc(3), xtest(3), tol
      real *8, allocatable :: xyz(:,:)
      complex *16, allocatable :: tab(:,:),tab2(:,:)
      complex *16 zk, im, zero, one
      logical ifsphere
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /           
      character (len=*) fname
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

      allocate(tab(ntarg0,npol3),tab2(ntarg0,npol3))
      call cpu_time(t1)
c$    t1 = omp_get_wtime()      

      ifsphere = .true.
      call h3dtabp_ref0(ndeg,zk,tol,tab,ntarg0,ifsphere)

      ifsphere = .false.
      call h3dtabp_ref0(ndeg,zk,tol,tab2,ntarg0,ifsphere)
      call cpu_time(t2)
c$    t2 = omp_get_wtime()            

      call prin2('time to generate table *',t2-t1,1)
      
c
c     brute force all entries
c     


      allocate(cv(ntarg0,ntest),ifail(ntarg0,ntest),
     1   absest(ntarg0,ntest),ccv(ntarg0,ntest),ccv2(ntarg0,ntest))

      do ii=1,ntest

        ix = npol_test(1,ii)
        iy = npol_test(2,ii)
        iz = npol_test(3,ii)

        ipol = ip2ind(ix,iy,iz)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ixyz)        
        do ixyz = 1,ntarg0
          print *, ixyz,ix,iy,iz
          xtest(1) = xyz(1,ixyz)
          xtest(2) = xyz(2,ixyz)
          xtest(3) = xyz(3,ixyz)

          call mksurhelm3dp0(xtest(1),xtest(2),xtest(3),ix,iy,iz,
     1     zk,n,cv(ixyz,ii),ifail(ixyz,ii),absest(ixyz,ii))

          ccv(ixyz,ii) = tab(ixyz,ipol)
          ccv2(ixyz,ii) = tab2(ixyz,ipol)


        enddo
C$OMP END PARALLEL DO        
      enddo

 1111 format(5(2x,i4),10(2x,e11.5))
      open(unit=33,file=trim(fname),access='append')
      do ii=1,ntest
        ix = npol_test(1,ii)
        iy = npol_test(2,ii)
        iz = npol_test(3,ii)

        ipol = ip2ind(ix,iy,iz)

        do ixyz = 1,ntarg0
          err1 = abs(cv(ixyz,ii)-ccv(ixyz,ii))
          err2 = abs(cv(ixyz,ii)-ccv2(ixyz,ii))
          err3 = abs(ccv(ixyz,ii)-ccv2(ixyz,ii))

          write(33,1111) ix,iy,iz,ixyz,ifail(ixyz,ii),
     1       real(cv(ixyz,ii)),imag(cv(ixyz,ii)),real(ccv(ixyz,ii)),
     2       imag(ccv(ixyz,ii)),real(ccv2(ixyz,ii)),
     3       imag(ccv2(ixyz,ii)),abs(absest(ixyz,ii)),err1,err2,err3
        enddo
      enddo

      return
      end




      integer function irand_int(nmax)
      implicit real *8 (a-h,o-z)
      
      irand_int = ceiling(nmax*hkrand(0))
      irand_int = min0(irand_int,nmax)
      irant_int = max0(irand_int,1)
      
      return
      end
c
c
c
c
c
      subroutine mksurhelm3dp0(x0,y0,z0,ix,iy,iz,zk0,norder0,
     1     value,ifail,absest)
c
c       This subroutine is the wrapper for computing
c       the integrals  
c     \int_{[-1,1]^3} e^{ikr}/r *
c          P_{ix-1}(x)*P_{iy-1}(y)*P_{iz-1}(z) dx dy dz \, ,
c        
c       where P_{n}(x) are legendre polynomials
c       and r = \sqrt{(x-x_{0})^2 + (y-y_{0})^2 + (z-z_{0})^2} 
c
c       Be careful, there are a few global variables
c
c       The kernel can be changed by appropriately changing
c       the function "fhelmgreen3d"
c
c       Input parameters:
c          x0,y0,z0 - coordinates of target location
c          ix,iy,iz - polynomial order in x,y, and z variables
c          norder - order of box code generation. 
c                   Must be greater than max(ix,iy,iz)
c       
c       Output parameters:
c          value - value of integral
c          ifail - ifail = 0 for successful computation of integral
c                     check other routines for error code otherwise
c          
c
c
c
c

      implicit none
      external fhelm3dp
      integer key, n, nf, ndim, mincls, maxcls, ifail, neval, nw
      parameter (ndim = 3, nw = 4000000, nf = 2)
      real *8 a(ndim), b(ndim)
      real *8, allocatable :: wrkstr(:)
      real *8 absest(nf), finest(nf), absreq, relreq,absest0
      real *8 xtarg,ytarg,ztarg
      real *8 x0, y0, z0
      complex *16 value, zk, zk0
      complex *16 im
      data im / (0.0d0,1.0d0) /
      integer ix,iy,iz,ixpol,iypol,izpol,norder,norder0
      common /cbh3dtab_brute/ xtarg,ytarg,ztarg,ixpol,iypol,izpol,
     1     norder,zk
c$omp threadprivate(/cbh3dtab_brute/)      

      xtarg = x0
      ytarg = y0
      ztarg = z0

      allocate(wrkstr(nw))

      ixpol = ix
      iypol = iy
      izpol = iz

      norder = norder0
      zk = zk0


      do 10 n = 1,ndim
         a(n) = -1.0d0
         b(n) =  1.0d0
   10 continue
      mincls = 0
      maxcls = 4000000
      key = 0
      absreq = 1d-12
      relreq = 1d-12
      ifail = 0

      call dcuhre(ndim, nf, a, b, mincls, maxcls, fhelm3dp, 
     1      absreq, relreq, key, nw, 0, finest, absest, neval,
     2      ifail, wrkstr)


      value = finest(1) + im*finest(2)

      return
      end



      


      subroutine h3dtabp_ref0(ndeg,zk,tol,tab,ldtab,ifsphere)
c
c     generate the Helmholtz potential table at the reference
c     points
c
c
c     input
c
c     ndeg - integer, highest degree of the basis polynomials
c                    (measured in total degree, e.g. x^0y^1z^2
c                     has total degree 3)
c     zk - complex*16, helmholtz parameter
c     tol - tolerance for error in table entries
c     ldtab - integer, leading dimension of the output table
c     ifsphere - flag for using spherical polynomials for
c        constructing tables
c
c     output
c      
c     tab - complex *16 array (ldtab,*), tab(i,j) is the integral
c     of the j-th tensor polynomial (in the ordering specified
c     in legetens.f) against the scaled green's function
c     exp(ikr)/r at the i-th reference target point (see tensrefpts3d)
      
c      
      implicit none
      integer ndeg, ldtab
      complex *16 tab(ldtab,*), zk
      real *8 tol
c     local
      integer idims(6), ndeg2, ii, jj, iface, npol2, npol3, ifdiff
      integer ntarg0, ntarg, n, idim, istart, npt, itype
      real *8 slicevals(6), flipd(6), flips(6), abszk, abszktol
      real *8 rcond, val, derscale, tol2, r, theta, phi
      parameter (abszktol = 2.5d0)

      real *8, allocatable :: x(:,:), w(:), pols(:), v(:,:)
      real *8 u, pi4
      integer ldu, ldv
      
      complex *16 zero, im, one
      complex *16, allocatable :: tabtemp(:,:), ahc(:,:), zv(:,:)
      complex *16, allocatable :: ahderc(:,:), ahcleg(:,:)
      complex *16, allocatable :: ahdercleg(:,:), ahelm(:), ahelms(:,:)
      complex *16, allocatable :: ahcleg3(:,:), leg2sph(:,:)
      complex *16, allocatable :: slp_pots(:,:), dlp_pots(:,:)

      
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /      
      data im / (0.0d0,1.0d0) /
      data slicevals / -1.0d0, 1.0d0, -1.0d0, 1.0d0, -1.0d0, 1.0d0 /
      data flipd / -1.0d0, 1.0d0, -1.0d0, 1.0d0, -1.0d0, 1.0d0 /
      data flips / 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0 /
      data idims / 1, 1, 2, 2, 3, 3 /
      
      logical ifsphere
      character type

      pi4 = 16*atan(1.0d0)
      
      n = ndeg+1
      
      type = 'T'
      
      abszk = abs(zk)

      if(ifsphere) ndeg2 = 20
      if(.not.ifsphere) ndeg2 = ndeg

      print *, ifsphere, ndeg2

c     memory for coeffs, etc.
      
      call legetens_npol_3d(ndeg,type,npol3)
      call legetens_npol_2d(ndeg2,type,npol2)

      allocate(ahc(npol2,npol3),ahderc(npol2,npol3))
      allocate(ahcleg3(npol3,npol3),leg2sph(npol3,npol3))
      allocate(ahcleg(npol2,npol3),ahdercleg(npol2,npol3))

      npt = n**3
      ldu = 1
      ldv = npt
      itype = 4
      allocate(x(3,npt),w(npt),v(npt,npol3))
      allocate(zv(npt,npol3))
      allocate(ahelm(npol3),ahelms(npol3,npt))
      call legetens_exps_3d(itype,n,type,x,u,ldu,v,ldv,w)

      
      if (.not. ifsphere) then
         tol2 = 1.0d-12
         call h3danti_legetens_form(ndeg,type,zk,tol2,ahcleg3,
     1        npol3,rcond)
         do ii = 1,npol3
            do jj = 1,npt
               zv(jj,ii) = -pi4*v(jj,ii)
            enddo
         enddo
         
         call zgemm('N','N',npt,npol3,npol3,one,zv,npt,
     1        ahcleg3,npol3,zero,tab,ldtab)
         
      else
         call legetens_spherepol(n,npol3,leg2sph)
         
         do jj = 1,npt
            call sphcart2polar(x(1,jj),r,theta,phi)
            call h3danti_sphere(ndeg,zk,r,theta,phi,ahelm)
            do ii = 1,npol3
               ahelms(ii,jj) = -pi4*ahelm(ii)
            enddo
         enddo

         call zgemm('T','N',npt,npol3,npol3,one,ahelms,npol3,
     1        leg2sph,npol3,zero,tab,ldtab)
         
      endif

c     memory depending on ntarg

      ntarg0 = 10*npt
      ntarg = 6*ntarg0
      allocate(slp_pots(npol2,ntarg),dlp_pots(npol2,ntarg))
      allocate(tabtemp(ntarg0,npol3))

      call h3d_facelayerpot_eval(tol,zk,ndeg,ndeg2,type,slp_pots,
     1  dlp_pots,npol2)

      do ii = 1,npol3
         do jj = npt+1,ntarg0
            tab(jj,ii) = zero
         enddo
      enddo

      do iface = 1,6

         idim = idims(iface)
         val = slicevals(iface)
         derscale = val

c     get poly coeffs of anti-helmholtzian of each polynomial
c     on face
         if (ifsphere) then
            call h3danti_sphere_slicecoeffs(ndeg2,type,ndeg,zk,
     1           idim,val,derscale,ahc,ahderc,npol2)
            call zgemm('N','N',npol2,npol3,npol3,one,ahc,npol2,
     1           leg2sph,npol3,zero,ahcleg,npol2)
            call zgemm('N','N',npol2,npol3,npol3,one,ahderc,npol2,
     1           leg2sph,npol3,zero,ahdercleg,npol2)
            
         else
            ifdiff = 1
            call legetens_slicezcoeffs_3d(ndeg,type,idim,val,
     1           ahcleg3,npol3,npol3,ahcleg,npol2,ifdiff,ahdercleg,
     2           npol2)
            do ii = 1,npol3
               do jj = 1,npol2
                  ahdercleg(jj,ii) = derscale*ahdercleg(jj,ii)
               enddo
            enddo
         endif
         

         istart = (iface-1)*ntarg0 + 1
         call zgemm('T','N',ntarg0,npol3,npol2,one,
     1        slp_pots(1,istart),
     1        npol2,ahdercleg,npol2,zero,tabtemp,ntarg0)

         do ii = 1,npol3
            do jj = 1,ntarg0
               tab(jj,ii) = tab(jj,ii) +
     1              flips(iface)*tabtemp(jj,ii)
            enddo
         enddo

         istart = (iface-1)*ntarg0 + 1 
         call zgemm('T','N',ntarg0,npol3,npol2,one,
     1        dlp_pots(1,istart),
     1        npol2,ahcleg,npol2,zero,tabtemp,ntarg0)

         do ii = 1,npol3
            do jj = 1,ntarg0
               tab(jj,ii) = tab(jj,ii) +
     1              flipd(iface)*tabtemp(jj,ii)
            enddo
         enddo
         
      enddo         
      
      return
      end

