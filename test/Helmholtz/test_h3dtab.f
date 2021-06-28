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
      a = hkrand(iseed)
      
      tol = 1.0d-9
      toltest = 1.0d-8
      
      ndeg = 3

c     number of points/polys to test per box
      npt_box = 2
      npol_box = 2

      
      call testall(zk,ndeg,npt_box,npol_box,tol,toltest)
      
      stop
      end
c
c      

      subroutine testall(zk,ndeg,npt_box,npol_box,tol,tol_test)
      implicit real *8 (a-h,o-z)
      integer, allocatable :: ip2ind(:,:,:), ind2p(:,:)
      integer iref(100), idimp(3,100), iflip(3,100)
      real *8 xq(200), u, v, w, xyzc(3), xtest(3), tol
      real *8, allocatable :: xyz(:,:), xcoll(:,:,:), xbtos(:,:,:)
      real *8, allocatable :: xstob(:,:,:)      
      complex *16, allocatable :: zints(:)
      complex *16, allocatable :: tab(:,:)
      complex *16, allocatable :: tabtemp(:,:)      
      complex *16, allocatable :: tabcoll(:,:,:)
      complex *16, allocatable :: tabbtos(:,:,:)
      complex *16, allocatable :: tabstob(:,:,:)      
      complex *16 zk, im, zero, one, cv, ccv
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /            
      character type
      logical all_pass, all_pass_coll, all_pass_btos, all_pass_stob

      all_pass = .true.
      all_pass_coll = .true.
      all_pass_btos = .true.
      all_pass_stob = .true.      
      
c     generate reference table

      n = ndeg + 1
      npt = n**3
      ntarg0 = 10*npt
      type = 't'
      call legetens_npol_3d(ndeg,type,npol)

      allocate(zints(npol))
      allocate(tab(ntarg0,npol))
      call cpu_time(t1)
c$    t1 = omp_get_wtime()      
      call h3dtabp_ref(ndeg,zk,tol,tab,ntarg0)
      call cpu_time(t2)
c$    t2 = omp_get_wtime()            

      call prin2('time to generate table *',t2-t1,1)
      
      call prin2('zk is *',zk,2)
      call prinf('ndeg is *',ndeg,1)      

      n = ndeg+1

      allocate(ip2ind(n,n,n))
      call legetens_pow2ind_3d(ndeg,type,ip2ind)

      allocate(ind2p(3,npol))
      call legetens_ind2pow_3d(ndeg,type,ind2p)
      
c
c     get all target points
c

      itype = 0
      call legeexps(itype,n,xq,u,v,w)

      do i = 1,n
         xq(i) = xq(i)+1.0d0
      enddo

      xyzc(1) = -1.0d0
      xyzc(2) = -1.0d0
      xyzc(3) = -1.0d0
      bs = 2.0d0

      nbtos = 56
      nstob = 56      
      ncoll = 27

      allocate(xcoll(3,npt,ncoll),xbtos(3,npt,nbtos),xstob(3,npt,nstob))
      
      call alltargs3d(xq,n,bs,xyzc,xcoll,xbtos,xstob)


      allocate(tabcoll(npt,npol,4),tabbtos(npt,npol,3),
     1     tabstob(npt,npol,3),tabtemp(npt,npol))
      ldtab = ntarg0
      call splitreftab3d(tab,ldtab,tabcoll,tabbtos,tabstob,npt,npol)


c     load symmetries and test

      call loadsymsc(iref,idimp,iflip)

      do i = 1,ncoll
         call buildtabfromsyms3d(ndeg,type,iref(i),idimp(1,i),
     1        iflip(1,i),tabcoll,tabtemp,npt,npol)

         do jj = 1,npol_box
            ipol = irand_int(npol)
            ix = ind2p(1,ipol)+1
            iy = ind2p(2,ipol)+1
            iz = ind2p(3,ipol)+1          
            do kk = 1,npt_box
               ipt = irand_int(npt)
               xtest(1) = xcoll(1,ipt,i)
               xtest(2) = xcoll(2,ipt,i)
               xtest(3) = xcoll(3,ipt,i)

c     brute force
               call mksurhelm3dp(xtest(1),xtest(2),xtest(3),ix,iy,iz,
     1              zk,n,cv,ifail)

c     table entry
               ccv = tabtemp(ipt,ipol)

               err1 = abs(cv-ccv)

               all_pass_coll =(all_pass_coll.and.(err1 .lt. tol_test))

               call prinf('colleague, box is *',i,1)               
               call prinf('ipt is *',ipt,1)
               call prinf('ipol is *',ipol,1)               
               call prinf('ix is *',ix,1)
               call prinf('iy is *',iy,1)
               call prinf('iz is *',iz,1)
               call prinf('ifail is *',ifail,1)
               call prin2_long('cv is *',cv,2)
               call prin2_long('ccv is *',ccv,2)               
               call prin2_long('error is *',err1,1)

            enddo
         enddo
      enddo

      call loadsymsbtos(iref,idimp,iflip)

      do i = 1,nbtos
         call buildtabfromsyms3d(ndeg,type,iref(i),idimp(1,i),
     1        iflip(1,i),tabbtos,tabtemp,npt,npol)

         do jj = 1,npol_box
            ipol = irand_int(npol)
            ix = ind2p(1,ipol)+1
            iy = ind2p(2,ipol)+1
            iz = ind2p(3,ipol)+1          
            do kk = 1,npt_box
               ipt = irand_int(npt)
               xtest(1) = xbtos(1,ipt,i)
               xtest(2) = xbtos(2,ipt,i)
               xtest(3) = xbtos(3,ipt,i)

c     brute force
               call mksurhelm3dp(xtest(1),xtest(2),xtest(3),ix,iy,iz,
     1              zk,n,cv,ifail)

c     table entry
               ccv = tabtemp(ipt,ipol)

               err1 = abs(cv-ccv)

               all_pass_btos = (all_pass_btos.and.(err1 .lt. tol_test))

               call prinf('big to small, box is *',i,1)               
               call prinf('ipt is *',ipt,1)
               call prinf('ipol is *',ipol,1)               
               call prinf('ix is *',ix,1)
               call prinf('iy is *',iy,1)
               call prinf('iz is *',iz,1)
               call prinf('ifail is *',ifail,1)
               call prin2_long('cv is *',cv,2)
               call prin2_long('ccv is *',ccv,2)               
               call prin2_long('error is *',err1,1)

            enddo
         enddo
      enddo


      call loadsymsstob(iref,idimp,iflip)

      do i = 1,nstob
         call buildtabfromsyms3d(ndeg,type,iref(i),idimp(1,i),
     1        iflip(1,i),tabstob,tabtemp,npt,npol)

         do jj = 1,npol_box
            ipol = irand_int(npol)
            ix = ind2p(1,ipol)+1
            iy = ind2p(2,ipol)+1
            iz = ind2p(3,ipol)+1          
            do kk = 1,npt_box
               ipt = irand_int(npt)
               xtest(1) = xstob(1,ipt,i)
               xtest(2) = xstob(2,ipt,i)
               xtest(3) = xstob(3,ipt,i)

c     brute force
               call mksurhelm3dp(xtest(1),xtest(2),xtest(3),ix,iy,iz,
     1              zk,n,cv,ifail)

c     table entry
               ccv = tabtemp(ipt,ipol)

               err1 = abs(cv-ccv)

               all_pass_stob = (all_pass_stob .and.(err1 .lt. tol_test))

               call prinf('small to big, box is *',i,1)               
               call prinf('ipt is *',ipt,1)
               call prinf('ipol is *',ipol,1)               
               call prinf('ix is *',ix,1)
               call prinf('iy is *',iy,1)
               call prinf('iz is *',iz,1)
               call prinf('ifail is *',ifail,1)
               call prin2_long('cv is *',cv,2)
               call prin2_long('ccv is *',ccv,2)               
               call prin2_long('error is *',err1,1)

            enddo
         enddo
      enddo

      all_pass = (all_pass_coll .and. all_pass_btos .and. all_pass_stob)

      
      write(*,*) 'all passed, coll?  ', all_pass_coll
      write(*,*) 'all passed, btos?  ', all_pass_btos
      write(*,*) 'all passed, stob?  ', all_pass_stob
      write(*,*) 'all passed?  ', all_pass      
      
      return
      end



      integer function irand_int(nmax)
      implicit real *8 (a-h,o-z)
      
      irand_int = ceiling(nmax*hkrand(0))
      irand_int = min0(irand_int,nmax)
      irant_int = max0(irand_int,1)
      
      return
      end
      
