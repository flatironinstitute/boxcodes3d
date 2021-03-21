      implicit real *8 (a-h,o-z)
      real *8 tts(3,100,2),errm_rel1(100,2),errm_rel2(100,2)
      real *8 errm_abs(100,2)
      real *8 tt3d(3),err_print(3)
      complex *16 zk, im, zero, one
      complex *16, allocatable :: tab_ref(:,:,:),tab(:,:,:),tab2(:,:,:)
      logical ifsphere
      character type
      character *2 arg_comm

      character *100 fname
      procedure (), pointer :: fker
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /           
      external h3d_vslp,h3d_vslp_gradx,h3d_vslp_grady,h3d_vslp_gradz
      
      call prini(6,13)
      

      call get_command_argument(1,arg_comm)
      read(arg_comm,*) iprec

      call get_command_argument(2,arg_comm)
      read(arg_comm,*) ndeg

      n = ndeg + 1

      type = 't'
      call legetens_npol_3d(ndeg,type,npol3)

      if(iprec.eq.0) tol = 1.0d-2
      if(iprec.eq.1) tol = 1.0d-3
      if(iprec.eq.2) tol = 1.0d-6
      if(iprec.eq.3) tol = 1.0d-9
      if(iprec.eq.4) tol = 1.0d-12

      npt = n**3

      ntarg0 = 10*npt



      allocate(tab_ref(ntarg0,npol3,3),tab(ntarg0,npol3,3))
      allocate(tab2(ntarg0,npol3,3))

c
c
c
c
      ifgen = 0
      izk = 4
      if(ifgen.eq.1) then
         if(izk.eq.1) zk = 1.6d0
         if(izk.eq.2) zk = 0.16d0
         if(izk.eq.3) zk = im*1.6d0
         if(izk.eq.4) zk = 5.0d0
         print *, "here"
         print *, "zk=",zk
         do ii=1,3
           if(ii.eq.1) fker=> h3d_vslp_gradx
           if(ii.eq.2) fker=> h3d_vslp_grady
           if(ii.eq.3) fker=> h3d_vslp_gradz
           call h3dtabp_ref_brute_fker(ndeg,tab_ref(1,1,ii),ntarg0,fker,
     1        dpars,zk,ipars)
         enddo
          write(fname,'(a,i2.2,a,i1,a)') 'tabref_grad_',n,'_izk_',izk,
     1       'new.dat'
          open(unit=33,file=trim(fname),action='readwrite',
     1       form='unformatted',access='stream')
         write(unit=33) tab_ref
         close(33)
         stop
      endif

      if(ifgen.ne.1) then

        do izk=4,4

          if(izk.eq.1) zk = 1.6d0
          if(izk.eq.2) zk = 0.16d0
          if(izk.eq.3) zk = im*1.6d0
          if(izk.eq.4) zk = 5.0d0
          write(fname,'(a,i2.2,a,i1,a)') 'tabref_grad_',n,'_izk_',izk,
     1       'new.dat'
          open(unit=33,file=trim(fname),action='readwrite',
     1       form='unformatted',access='stream')
          read(unit=33) tab_ref

          print *, "done reading"

          close(33)

          ifsphere = .false.
          if(ifsphere) ndeg2 = 20
          if(.not.ifsphere) ndeg2 = ndeg
          
          iflg = 1

          call h3dtabp_ref2_grad(ndeg,zk,tol,tab,ntarg0,npol3,
     1      ifsphere,ndeg2,iflg,tts(1,izk,iflg),tts(2,izk,iflg),
     2      tts(3,izk,iflg))
          
          
          iflg = 2

          call h3dtabp_ref2_grad(ndeg,zk,tol,tab,ntarg0,npol3,
     1      ifsphere,ndeg2,iflg,tts(1,izk,iflg),tts(2,izk,iflg),
     2      tts(3,izk,iflg))
          print *, "Done with 2d adap quad"
        
          call prin2('tts iflg1=*',tts(1,izk,1),3)
          call prin2('tts iflg2=*',tts(1,izk,2),3)

          call cpu_time(t1)
C$           t1 = omp_get_wtime()          
          do idir=1,3
             if(idir.eq.1) fker=> h3d_vslp_gradx
             if(idir.eq.2) fker=> h3d_vslp_grady
             if(idir.eq.3) fker=> h3d_vslp_gradz
            call h3dtabp_ref_brute_new_fker(ndeg,tol,tab2(1,1,idir),
     1         ntarg0,fker,dpars,zk,ipars)
          enddo
          call cpu_time(t2)
C$           t2 = omp_get_wtime()          
cc     
          print *, "done with 3d adap quad"
          errmax1 = 0
          errmax2 = 0
          errmaxa = 0
          rmax = 0
          do i=1,npol3
            do j=npt+1,ntarg0
             ra = max(abs(tab_ref(j,i,1)),1.0d0)
             erra = abs(tab(j,i,1)-tab_ref(j,i,1))/ra
             if(erra.gt.errmax1) errmax1 = erra

             
             erra = abs(tab(j,i,1)-tab_ref(j,i,1))/abs(tab_ref(j,i,1))
             if(erra.gt.errmax2) errmax2 = erra
             if(i.lt.3.and.j-npt.lt.5) then
               call prin2('tab=*',tab(j,i,1),6)
               call prin2('tab_ref=*',tab_ref(j,i,1),6)
             endif


             erra = abs(tab(j,i,1)-tab_ref(j,i,1))
             if(erra.gt.errmaxa) errmaxa = erra
             ra = abs(tab_ref(j,i,1))
             if(ra.gt.rmax) rmax = ra
            enddo
          enddo

          errm_rel1(izk,1) = errmax1
          errm_rel2(izk,1) = errmax2
          errm_abs(izk,1) = errmaxa/rmax
          err_print(1) = errmax1
          err_print(2) = errmax2
          err_print(3) = errmaxa/rmax

          call prin2(' *',i,0)
          call prin2(' *',i,0)
          call prin2('2d anti helmholtzian*',i,0)
          call prin2(' *',i,0)
          call prin2(' *',i,0)
          call prin2('max errors=*',err_print,3)
          call prin2('time taken=*',tts(:,izk,:),6)


          errmax1 = 0
          errmax2 = 0
          errmaxa = 0
          rmax = 0
          do i=1,npol3
            do j=npt+1,ntarg0
             ra = max(abs(tab_ref(j,i,1)),1.0d0)
             erra = abs(tab2(j,i,1)-tab_ref(j,i,1))/ra
             if(erra.gt.errmax1) errmax1 = erra

             
             erra = abs(tab2(j,i,1)-tab_ref(j,i,1))/abs(tab_ref(j,i,1))

             if(erra.gt.errmax2) errmax2 = erra


             erra = abs(tab2(j,i,1)-tab_ref(j,i,1))
             if(erra.gt.errmaxa) errmaxa = erra
             ra = abs(tab_ref(j,i,1))
             if(ra.gt.rmax) rmax = ra
            enddo
          enddo

          call prin2(' *',i,0)
          call prin2(' *',i,0)
          call prin2('3d adap integration comp*',i,0)
          call prin2(' *',i,0)
          call prin2(' *',i,0)
          errm_rel1(izk,2) = errmax1
          errm_rel2(izk,2) = errmax2
          errm_abs(izk,2) = errmaxa/rmax
          err_print(1) = errmax1
          err_print(2) = errmax2
          err_print(3) = errmaxa/rmax
          call prin2('max error=*',err_print,3)
          tt3d(izk) = t2-t1
          call prin2('time taken=*',tt3d(izk),1)

 1111     continue          
        enddo
        call prin2(' *',i,0)
        call prin2(' *',i,0)
        call prin2('results summary:*',i,0)
        call prin2(' *',i,0)
        call prin2(' *',i,0)

        call prin2('max errors relative 1=*',errm_rel1(4,1:2),2)
        call prin2('max errors relative 2=*',errm_rel2(4,1:2),2)
        call prin2('max errors absolute=*',errm_abs(4,1:2),2)
        call prin2('time taken 2d=*',tts(:,4,:),6)
        call prin2('time taken 3d adap=*',tt3d(4),1)
        stop
      endif



      
      stop
      end
c
c      
      
      subroutine h3dtabp_ref2_grad(ndeg,zk,tol,tab,ldtab1,ldtab2,
     1   ifsphere,ndeg2,iflg,tpat,tlpadap,tlpspr)
c
c     generate the Helmholtz gradient table at the reference
c     points.
c
c
c     input
c
c     ndeg - integer, highest degree of the basis polynomials
c                    (measured in total degree, e.g. x^0y^1z^2
c                     has total degree 3)
c     zk - complex*16, helmholtz parameter
c     tol - tolerance for error in table entries
c     ldtab1 - integer, leading dimension of the output table
c     ldtab2 - integer, second leading dimension of the output table
c     ifsphere - whether to use spherical polynomials          
c     ndeg2 - order used in adaptive integration
c     iflg - flag for determining which version of adaptive
c              integration to use
c            if iflg = 1, one with precomputation is used
c            if iflg = 2, one without precomputation is used
          
          
c
c     output
c      
c     tab - complex *16 array (ldtab,*,3), tab(i,j,1) is the integral
c     of the j-th tensor polynomial (in the ordering specified
c     in legetens.f) against the scaled green's function
c     d/dx (exp(ikr)/r) at the i-th reference target point (see tensrefpts3d)
c     tab(i,j,2:3) are the corresponding values of d/dy and d/dz
c
c     tpat - time taken in getting particular solution
c     tlpadap - time taken in layer potential computation in adaptive
c                integration
c     tlpspr - time taken to spread layer potential to all targets
c
c      
      implicit none
      integer ndeg, ldtab1,ldtab2
      complex *16 tab(ldtab1,ldtab2,3), zk
      real *8 tol
c     local
      integer idims(6), ndeg2, ii, jj, iface, npol2, npol3, ifdiff
      integer idimgx(6),idimgy(6),idimgz(6),idimgcomb(3,6)
      integer ntarg0, ntarg, n, idim, istart, npt, itype
      real *8 slicevals(6), flipd(6), flips(6), abszk, abszktol
      real *8 rcond, val, derscale, tol2, r, theta, phi
      real *8 tpat,tlpadap,tlpspr,t1,t2,omp_get_wtime
      integer i,j,idir
      parameter (abszktol = 2.5d0)

      real *8, allocatable :: x(:,:), w(:), pols(:), v(:,:)
      real *8 u, pi4
      integer ldu, ldv
      integer iflg
      
      complex *16 zero, im, one
      complex *16, allocatable :: tabtemp(:,:), ahc(:,:), zv(:,:)
      complex *16, allocatable :: ahderc(:,:), ahcleg(:,:)
      complex *16, allocatable :: ahdercleg(:,:), ahelm(:), ahelms(:,:)
      complex *16, allocatable :: ahcleg3(:,:), leg2sph(:,:)
      complex *16, allocatable :: slp_pots(:,:), slp_grad(:,:,:)
      complex *16, allocatable :: dlp_pots(:,:), dlp_grad(:,:,:)
      complex *16, allocatable :: ahcleg3grad(:,:,:)
      real *8, allocatable :: dmat(:,:),vtmp1(:),vtmp2(:)

      
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /      
      data im / (0.0d0,1.0d0) /
      data slicevals / -1.0d0, 1.0d0, -1.0d0, 1.0d0, -1.0d0, 1.0d0 /
      data flipd / -1.0d0, 1.0d0, -1.0d0, 1.0d0, -1.0d0, 1.0d0 /
      data flips / 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0 /
      data idims / 1, 1, 2, 2, 3, 3 /
      data idimgx /3,3,1,1,1,1/
      data idimgy /1,1,3,3,2,2/
      data idimgz /2,2,2,2,3,3/
      
      logical ifsphere
      character type

      pi4 = 16*atan(1.0d0)

      do i=1,6
        idimgcomb(1,i) = idimgx(i)
        idimgcomb(2,i) = idimgy(i)
        idimgcomb(3,i) = idimgz(i)
      enddo

      call prinf('idimgcomb=*',idimgcomb,18)
      
      n = ndeg+1
      
      type = 'T'
      
      abszk = abs(zk)

c     determine anti-helmholtzian technique
      


      call cpu_time(t1)
C$      t1 = omp_get_wtime()
c     memory for coeffs, etc.
      
      call legetens_npol_3d(ndeg,type,npol3)
      call legetens_npol_2d(ndeg2,type,npol2)

      allocate(ahc(npol2,npol3),ahderc(npol2,npol3))
      allocate(ahcleg3(npol3,npol3),leg2sph(npol3,npol3))
      allocate(ahcleg(npol2,npol3),ahdercleg(npol2,npol3))
      allocate(ahcleg3grad(npol3,npol3,3),dmat(n,n))
      allocate(vtmp1(npol3),vtmp2(npol3))

      dmat = 0
      call legecoeff_dmat(ndeg,dmat,ndeg+1)

      npt = n**3
      ldu = 1
      ldv = npt
      itype = 4
      allocate(x(3,npt),w(npt),v(npt,npol3))
      allocate(zv(npt,npol3))
      allocate(ahelm(npol3),ahelms(npol3,npt))
      call legetens_exps_3d(itype,n,type,x,u,ldu,v,ldv,w)

      
      tol2 = 1.0d-12
      call h3danti_legetens_form(ndeg,type,zk,tol2,ahcleg3,
     1        npol3,rcond)
      do ii = 1,npol3
         do jj = 1,npt
            zv(jj,ii) = -pi4*v(jj,ii)
         enddo
      enddo

      do idir = 1,3
         do i=1,npol3
           do j=1,npol3
             vtmp1(j) = real(ahcleg3(j,i))
             vtmp2(j) = 0
           enddo
           call legediff_3d(ndeg,type,vtmp1,idir,dmat,vtmp2)
           do j=1,npol3
             ahcleg3grad(j,i,idir) = vtmp2(j)
           enddo

           do j=1,npol3
             vtmp1(j) = imag(ahcleg3(j,i))
             vtmp2(j) = 0
           enddo
           call legediff_3d(ndeg,type,vtmp1,idir,dmat,vtmp2)
           do j=1,npol3
             ahcleg3grad(j,i,idir) = ahcleg3grad(j,i,idir)+im*vtmp2(j)
           enddo
         enddo
      enddo
       
c
c
c   update this part to handle the gradient
c
      do idir=1,3
        call zgemm('N','N',npt,npol3,npol3,one,zv,npt,
     1          ahcleg3grad(1,1,idir),npol3,zero,tab(1,1,idir),ldtab1)
      enddo
         
      call cpu_time(t2)
C$       t2 = omp_get_wtime()      

      tpat = t2-t1
c     memory depending on ntarg

      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      ntarg0 = 10*npt
      ntarg = 6*ntarg0
      allocate(slp_pots(npol2,ntarg),slp_grad(npol2,ntarg,3))
      allocate(dlp_pots(npol2,ntarg),dlp_grad(npol2,ntarg,3))
      allocate(tabtemp(ntarg0,npol3))

      call h3d_facelayergrad_eval_new(tol,zk,ndeg,ndeg2,type,iflg,
     1     slp_pots,slp_grad(1,1,1),slp_grad(1,1,2),
     2     slp_grad(1,1,3),dlp_pots,dlp_grad(1,1,1),dlp_grad(1,1,2),
     3     dlp_grad(1,1,3),npol2)

cc      do idir=1,3
cc         do ii = 1,npol3
cc            do jj = npt+1,ntarg0
cc               tab(jj,ii,idir) = zero
cc            enddo
cc         enddo
cc      enddo
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      tlpadap = t2-t1


      call cpu_time(t1)
C$      t1 = omp_get_wtime()      

      do iface = 1,6
         idim = idims(iface)
         val = slicevals(iface)
         derscale = val

         ifdiff = 1
         call legetens_slicezcoeffs_3d(ndeg,type,idim,val,
     1     ahcleg3,npol3,npol3,ahcleg,npol2,ifdiff,
     2     ahdercleg,npol2)
         do ii = 1,npol3
            do jj = 1,npol2
               ahdercleg(jj,ii) = derscale*ahdercleg(jj,ii)
            enddo
         enddo

         istart = (iface-1)*ntarg0 + 1
         do idir=1,3
            call zgemm('T','N',ntarg0,npol3,npol2,one,
     1         slp_grad(1,istart,idimgcomb(idir,iface)),
     1         npol2,ahdercleg,npol2,zero,tabtemp,ntarg0)

            do ii = 1,npol3
               do jj = 1,ntarg0
                  tab(jj,ii,idir) = tab(jj,ii,idir) +
     1                 flips(iface)*tabtemp(jj,ii)
               enddo
            enddo

            istart = (iface-1)*ntarg0 + 1 
            call zgemm('T','N',ntarg0,npol3,npol2,one,
     1           dlp_grad(1,istart,idimgcomb(idir,iface)),
     1           npol2,ahcleg,npol2,zero,tabtemp,ntarg0)

            do ii = 1,npol3
               do jj = 1,ntarg0
                  tab(jj,ii,idir) = tab(jj,ii,idir) +
     1              flipd(iface)*tabtemp(jj,ii)
               enddo
            enddo
         enddo
      enddo         
      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      tlpspr = t2-t1

      return
      end
c      
c
      
      subroutine h3dtabp_ref2_grad_nw(ndeg,zk,tol,tab,ldtab1,ldtab2,
     1   ifsphere,ndeg2,iflg,tpat,tlpadap,tlpspr)
c
c     generate the Helmholtz gradient table at the reference
c     points.
c
c
c     Not functioning, due to incorrect jump in normal
c     component on the boundary
c
c
c     input
c
c     ndeg - integer, highest degree of the basis polynomials
c                    (measured in total degree, e.g. x^0y^1z^2
c                     has total degree 3)
c     zk - complex*16, helmholtz parameter
c     tol - tolerance for error in table entries
c     ldtab1 - integer, leading dimension of the output table
c     ldtab2 - integer, second leading dimension of the output table
c     ifsphere - whether to use spherical polynomials          
c     ndeg2 - order used in adaptive integration
c     iflg - flag for determining which version of adaptive
c              integration to use
c            if iflg = 1, one with precomputation is used
c            if iflg = 2, one without precomputation is used
          
          
c
c     output
c      
c     tab - complex *16 array (ldtab,*,3), tab(i,j,1) is the integral
c     of the j-th tensor polynomial (in the ordering specified
c     in legetens.f) against the scaled green's function
c     d/dx (exp(ikr)/r) at the i-th reference target point (see tensrefpts3d)
c     tab(i,j,2:3) are the corresponding values of d/dy and d/dz
c
c     tpat - time taken in getting particular solution
c     tlpadap - time taken in layer potential computation in adaptive
c                integration
c     tlpspr - time taken to spread layer potential to all targets
c
c      
      implicit none
      integer ndeg, ldtab1,ldtab2
      complex *16 tab(ldtab1,ldtab2,3), zk
      real *8 tol
c     local
      integer idims(6), ndeg2, ii, jj, iface, npol2, npol3, ifdiff
      integer ntarg0, ntarg, n, idim, istart, npt, itype
      real *8 slicevals(6), flipd(6), flips(6), abszk, abszktol
      real *8 rcond, val, derscale, tol2, r, theta, phi
      real *8 tpat,tlpadap,tlpspr,t1,t2,omp_get_wtime
      integer i,j,idir
      parameter (abszktol = 2.5d0)

      real *8, allocatable :: x(:,:), w(:), pols(:), v(:,:)
      real *8 u, pi4
      integer ldu, ldv
      integer iflg
      
      complex *16 zero, im, one
      complex *16, allocatable :: tabtemp(:,:), ahc(:,:), zv(:,:)
      complex *16, allocatable :: ahderc(:,:), ahcleg(:,:)
      complex *16, allocatable :: ahdercleg(:,:), ahelm(:), ahelms(:,:)
      complex *16, allocatable :: ahcleg3(:,:), leg2sph(:,:)
      complex *16, allocatable :: slp_pots(:,:), dlp_pots(:,:)
      complex *16, allocatable :: ahcleg3grad(:,:,:)
      real *8, allocatable :: dmat(:,:),vtmp1(:),vtmp2(:)

      
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

c     determine anti-helmholtzian technique
      


      call cpu_time(t1)
C$      t1 = omp_get_wtime()
c     memory for coeffs, etc.
      
      call legetens_npol_3d(ndeg,type,npol3)
      call legetens_npol_2d(ndeg2,type,npol2)

      allocate(ahc(npol2,npol3),ahderc(npol2,npol3))
      allocate(ahcleg3(npol3,npol3),leg2sph(npol3,npol3))
      allocate(ahcleg(npol2,npol3),ahdercleg(npol2,npol3))
      allocate(ahcleg3grad(npol3,npol3,3),dmat(n,n))
      allocate(vtmp1(npol3),vtmp2(npol3))

      dmat = 0
      call legecoeff_dmat(ndeg,dmat,ndeg+1)

      npt = n**3
      ldu = 1
      ldv = npt
      itype = 4
      allocate(x(3,npt),w(npt),v(npt,npol3))
      allocate(zv(npt,npol3))
      allocate(ahelm(npol3),ahelms(npol3,npt))
      call legetens_exps_3d(itype,n,type,x,u,ldu,v,ldv,w)

      
      tol2 = 1.0d-12
      call h3danti_legetens_form(ndeg,type,zk,tol2,ahcleg3,
     1        npol3,rcond)
      do ii = 1,npol3
         do jj = 1,npt
            zv(jj,ii) = -pi4*v(jj,ii)
         enddo
      enddo

      do idir = 1,3
         do i=1,npol3
           do j=1,npol3
             vtmp1(j) = real(ahcleg3(j,i))
             vtmp2(j) = 0
           enddo
           call legediff_3d(ndeg,type,vtmp1,idir,dmat,vtmp2)
           do j=1,npol3
             ahcleg3grad(j,i,idir) = vtmp2(j)
           enddo

           do j=1,npol3
             vtmp1(j) = imag(ahcleg3(j,i))
             vtmp2(j) = 0
           enddo
           call legediff_3d(ndeg,type,vtmp1,idir,dmat,vtmp2)
           do j=1,npol3
             ahcleg3grad(j,i,idir) = ahcleg3grad(j,i,idir)+im*vtmp2(j)
           enddo
         enddo
      enddo
       
c
c
c   update this part to handle the gradient
c
      do idir=1,3
        call zgemm('N','N',npt,npol3,npol3,one,zv,npt,
     1          ahcleg3grad(1,1,idir),npol3,zero,tab(1,1,idir),ldtab1)
      enddo
         
      call cpu_time(t2)
C$       t2 = omp_get_wtime()      

      tpat = t2-t1
c     memory depending on ntarg

      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      ntarg0 = 10*npt
      ntarg = 6*ntarg0
      allocate(slp_pots(npol2,ntarg),dlp_pots(npol2,ntarg))
      allocate(tabtemp(ntarg0,npol3))

      call h3d_facelayerpot_eval_new(tol,zk,ndeg,ndeg2,type,iflg,
     1     slp_pots,dlp_pots,npol2)

      do idir=1,3
         do ii = 1,npol3
            do jj = npt+1,ntarg0
               tab(jj,ii,idir) = zero
            enddo
         enddo
      enddo
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      tlpadap = t2-t1


      call cpu_time(t1)
C$      t1 = omp_get_wtime()      

      do idir=1,3
         do iface = 1,6

            idim = idims(iface)
            val = slicevals(iface)
            derscale = val

            ifdiff = 1
            call legetens_slicezcoeffs_3d(ndeg,type,idim,val,
     1        ahcleg3grad(1,1,idir),npol3,npol3,ahcleg,npol2,ifdiff,
     2        ahdercleg,npol2)
            do ii = 1,npol3
               do jj = 1,npol2
                  ahdercleg(jj,ii) = derscale*ahdercleg(jj,ii)
               enddo
           enddo
         

            istart = (iface-1)*ntarg0 + 1
            call zgemm('T','N',ntarg0,npol3,npol2,one,
     1           slp_pots(1,istart),
     1           npol2,ahdercleg,npol2,zero,tabtemp,ntarg0)

            do ii = 1,npol3
               do jj = 1,ntarg0
                  tab(jj,ii,idir) = tab(jj,ii,idir) +
     1                 flips(iface)*tabtemp(jj,ii)
               enddo
            enddo

            istart = (iface-1)*ntarg0 + 1 
            call zgemm('T','N',ntarg0,npol3,npol2,one,
     1           dlp_pots(1,istart),
     1           npol2,ahcleg,npol2,zero,tabtemp,ntarg0)

            do ii = 1,npol3
               do jj = 1,ntarg0
                  tab(jj,ii,idir) = tab(jj,ii,idir) +
     1              flipd(iface)*tabtemp(jj,ii)
               enddo
            enddo
         enddo
      enddo         
      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      tlpspr = t2-t1

      return
      end
