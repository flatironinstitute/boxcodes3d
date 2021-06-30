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
      
      tol = 1.0d-5
      toltest = 1.0d-4
      
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
      real *8, allocatable :: x(:,:)
      integer, allocatable :: ipt2depth(:,:),idepth2pt(:,:,:)
      real *8, allocatable :: xyz(:,:), xcoll(:,:,:), xbtos(:,:,:)
      real *8, allocatable :: xstob(:,:,:) 
      real *8, allocatable :: pmat(:,:),pmat_qr(:,:),pmat_tau(:),pols(:)
      integer, allocatable :: pmat_jpvt(:)
      complex *16, allocatable :: zints(:,:),zintsp(:)
      complex *16, allocatable :: tab(:,:)
      complex *16, allocatable :: tabtemp(:,:)      
      complex *16, allocatable :: tabcoll(:,:,:)
      complex *16, allocatable :: tabbtos(:,:,:)
      complex *16, allocatable :: tabstob(:,:,:)   
      complex *16, allocatable :: tab_ref(:,:)
      complex *16, allocatable :: tabcoll_ref(:,:,:)
      complex *16, allocatable :: tabbtos_ref(:,:,:)
      complex *16, allocatable :: tabstob_ref(:,:,:)
      complex *16, allocatable :: tabtemp_ref(:,:)
      complex *16, allocatable :: tabtemp_ref_coefs(:,:)
      complex *16 zk, im, zero, one, cv, ccv
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /            
      character type
      logical all_pass, all_pass_coll, all_pass_btos, all_pass_stob
      external h3d_vslp

      all_pass = .true.
      all_pass_coll = .true.
      all_pass_btos = .true.
      all_pass_stob = .true.      
      
c     generate reference table

      n = ndeg + 1
      nptmax = n**3
      allocate(x(3,nptmax))

      
      allocate(ipt2depth(3,nptmax),idepth2pt(n,n,n))
      type = 't'
      call legetens_npol_3d(ndeg,type,npol)

      ifloor = 1.75d0*npol
      print *, "n=",n
      print *, "ifloor=",ifloor
      call fakepolya3d(n,ifloor,npt,x,ipt2depth,idepth2pt)
      itype = 0
      call prinf('npt=*',npt,1)
      call prin2('rat=*',(npt+0.0d0)/(nptmax+0.0d0),1)
cc      call prin2('x=*',x,3*npt)



      ntarg0 = 10*npt
      npol0 = 10*npol
c
c  get coefs2vals matrix and store it in pmat
c
c

      allocate(pmat(npt,npol),pmat_qr(npt,npol),pols(npol))
      allocate(pmat_tau(npol),pmat_jpvt(npol))

      do i=1,npt
        call legetens_pols_3d(x(1,i),ndeg,type,pols)
        do j=1,npol
          pmat(i,j) = pols(j)
        enddo
      enddo

      call get_qrdecomp(npt,npol,pmat,pmat_qr,pmat_jpvt,pmat_tau)

      allocate(zints(npt,npol),zintsp(npol))
      allocate(tab(npol0,npol))
      allocate(tab_ref(ntarg0,npol))


      nup = 8
      iflg = 1
      tab = 0
      call cpu_time(t1)
c$    t1 = omp_get_wtime()      
      call h3dtabp_ref(ndeg,zk,tol,npt,x,tab,npol0,nup,iflg,
     1   pmat_qr,npol,pmat_jpvt,pmat_tau,tt1,tt2,tt3)
      call cpu_time(t2)
c$    t2 = omp_get_wtime()            

      call prin2('time to generate table *',t2-t1,1)
      call prin2('tab=*',tab,12)

      call prin2('zk is *',zk,2)
      call prinf('ndeg is *',ndeg,1)      

      n = ndeg+1
c
c     get all target points
c

      xyzc(1) = -1.0d0
      xyzc(2) = -1.0d0
      xyzc(3) = -1.0d0
      bs = 2.0d0

      nbtos = 56
      nstob = 56      
      ncoll = 27

      allocate(xcoll(3,npt,ncoll),xbtos(3,npt,nbtos),xstob(3,npt,nstob))
      
      call alltargs3d_grid(x,npt,bs,xyzc,xcoll,xbtos,xstob)


      allocate(tabcoll(npol,npol,4),tabbtos(npol,npol,3),
     1     tabstob(npol,npol,3),tabtemp(npol,npol))
      allocate(tabtemp_ref(npt,npol),tabtemp_ref_coefs(npol,npol))
      ldtab = npol0
      call splitreftab3dcc(tab,ldtab,tabcoll,tabbtos,tabstob,npol)

c     load symmetries and test

      call loadsymsc(iref,idimp,iflip)
      ifprint = 0
      print *, "npt=",npt
      do i = 1,ncoll
         call buildtabfromsyms3dcc(ndeg,type,iref(i),idimp(1,i),
     1        iflip(1,i),tabcoll,tabtemp,npol)
         ifself = 0
         if(i.eq.13) ifself = 1
         call h3dtabbox_ref_brute(ndeg,tol,xcoll(1,1,i),npt,
     1     tabtemp_ref,npol,h3d_vslp,dpars,zk,ipars)

         call zqrsolv(npt,npol,pmat_qr,pmat_jpvt,pmat_tau,npol,
     1     tabtemp_ref,tabtemp_ref_coefs)
         do j=1,npol
            do l=1,npol
              err1 = abs(tabtemp_ref_coefs(l,j) - tabtemp(l,j))
              all_pass_coll =(all_pass_coll.and.(err1 .lt. tol_test))
              if(err1.gt.tol_test) then
                print *, i,l,j,err1
                print *, real(tabtemp_ref_coefs(l,j)),real(tabtemp(l,j))
                print *, imag(tabtemp_ref_coefs(l,j)),imag(tabtemp(l,j))
                print *, tabtemp_ref_coefs(l,j)/tabtemp(l,j)
              endif
            enddo
         enddo

         
      enddo

      
      call loadsymsbtos(iref,idimp,iflip)
      do i = 1,nbtos
         call buildtabfromsyms3dcc(ndeg,type,iref(i),idimp(1,i),
     1        iflip(1,i),tabbtos,tabtemp,npol)
         ifself = 0
         call h3dtabbox_ref_brute(ndeg,tol,xbtos(1,1,i),npt,
     1     tabtemp_ref,npol,h3d_vslp,dpars,zk,ipars)

         call zqrsolv(npt,npol,pmat_qr,pmat_jpvt,pmat_tau,npol,
     1     tabtemp_ref,tabtemp_ref_coefs)
         do j=1,npol
            do l=1,npol
              err1 = abs(tabtemp_ref_coefs(l,j) - tabtemp(l,j))
              all_pass_btos =(all_pass_btos.and.(err1 .lt. tol_test))
              if(err1.gt.tol_test) then
                print *, i,j,l,err1
                print *, real(tabtemp_ref_coefs(l,j)),real(tabtemp(l,j))
                print *, imag(tabtemp_ref_coefs(l,j)),imag(tabtemp(l,j))
              endif
            enddo
         enddo
      enddo



      call loadsymsstob(iref,idimp,iflip)

      do i = 1,nstob
         call buildtabfromsyms3dcc(ndeg,type,iref(i),idimp(1,i),
     1        iflip(1,i),tabstob,tabtemp,npol)
         ifself = 0
         call h3dtabbox_ref_brute(ndeg,tol,xstob(1,1,i),npt,
     1     tabtemp_ref,npol,h3d_vslp,dpars,zk,ipars)

         call zqrsolv(npt,npol,pmat_qr,pmat_jpvt,pmat_tau,npol,
     1     tabtemp_ref,tabtemp_ref_coefs)
         do j=1,npol
            do l=1,npol
              err1 = abs(tabtemp_ref_coefs(l,j) - tabtemp(l,j))
              all_pass_stob =(all_pass_stob.and.(err1 .lt. tol_test))
              if(err1.gt.tol_test) then
                print *, i,j,l,err1
                print *, real(tabtemp_ref_coefs(l,j)),real(tabtemp(l,j))
                print *, imag(tabtemp_ref_coefs(l,j)),imag(tabtemp(l,j))
              endif
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
      
