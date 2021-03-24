      implicit real *8 (a-h,o-z)
      real *8 tts(3,100,2),errm_rel1(100,2),errm_rel2(100,2)
      real *8 errm_abs(100,2)
      real *8 tt3d(3),err_print(3)
      complex *16 zk, im, zero, one
      complex *16, allocatable :: tab_ref(:,:),tab(:,:),tab2(:,:)
      logical ifsphere
      character type
      character *2 arg_comm

      character *100 fname
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /           
      external h3d_vslp
      
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



      allocate(tab_ref(ntarg0,npol3),tab(ntarg0,npol3))
      allocate(tab2(ntarg0,npol3))

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
         call h3dtabp_ref_brute_fker(ndeg,tab_ref,ntarg0,h3d_vslp,
     1      dpars,zk,ipars)
          write(fname,'(a,i2.2,a,i1,a)') 'tabref_',n,'_izk_',izk,
     1       'new.dat'
          open(unit=33,file=trim(fname),action='readwrite',
     1       form='unformatted',access='stream')
         write(unit=33) tab_ref
         close(33)
         stop
      endif

      if(igen.ne.1) then

        do izk=4,4

          if(izk.eq.1) zk = 1.6d0
          if(izk.eq.2) zk = 0.16d0
          if(izk.eq.3) zk = im*1.6d0
          if(izk.eq.4) zk = 5.0d0
          write(fname,'(a,i2.2,a,i1,a)') 'tabref_',n,'_izk_',izk,
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

          call h3dtabp_ref2(ndeg,zk,tol,tab,ntarg0,ifsphere,ndeg2,
     1      iflg,tts(1,izk,iflg),tts(2,izk,iflg),tts(3,izk,iflg))
          
          
          iflg = 2

          call h3dtabp_ref2(ndeg,zk,tol,tab,ntarg0,ifsphere,ndeg2,
     1      iflg,tts(1,izk,iflg),tts(2,izk,iflg),tts(3,izk,iflg))
          print *, "Done with 2d adap quad"
        
          call prin2('tts iflg1=*',tts(1,izk,1),3)
          call prin2('tts iflg2=*',tts(1,izk,2),3)

          call cpu_time(t1)
C$           t1 = omp_get_wtime()          
          call h3dtabp_ref_brute_new_fker(ndeg,tol,tab2,ntarg0,
     1        h3d_vslp,dpars,zk,ipars)
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
             ra = max(abs(tab_ref(j,i)),1.0d0)
             erra = abs(tab(j,i)-tab_ref(j,i))/ra
             if(erra.gt.errmax1) errmax1 = erra

             
             erra = abs(tab(j,i)-tab_ref(j,i))/abs(tab_ref(j,i))
             if(erra.gt.errmax2) errmax2 = erra


             erra = abs(tab(j,i)-tab_ref(j,i))
             if(erra.gt.errmaxa) errmaxa = erra
             ra = abs(tab_ref(j,i))
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
             ra = max(abs(tab_ref(j,i)),1.0d0)
             erra = abs(tab2(j,i)-tab_ref(j,i))/ra
             if(erra.gt.errmax1) errmax1 = erra

             
             erra = abs(tab2(j,i)-tab_ref(j,i))/abs(tab_ref(j,i))
             if(erra.gt.errmax2) errmax2 = erra


             erra = abs(tab2(j,i)-tab_ref(j,i))
             if(erra.gt.errmaxa) errmaxa = erra
             ra = abs(tab_ref(j,i))
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
c      
c
      
      subroutine h3dtabp_ref2(ndeg,zk,tol,tab,ldtab,ifsphere,ndeg2,
     1   iflg,tpat,tlpadap,tlpspr)
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
c     ifsphere - whether to use spherical polynomials          
c     ndeg2 - order used in adaptive integration
c     iflg - flag for determining which version of adaptive
c              integration to use
c            if iflg = 1, one with precomputation is used
c            if iflg = 2, one without precomputation is used
          
          
c
c     output
c      
c     tab - complex *16 array (ldtab,*), tab(i,j) is the integral
c     of the j-th tensor polynomial (in the ordering specified
c     in legetens.f) against the scaled green's function
c     exp(ikr)/r at the i-th reference target point (see tensrefpts3d)
c
c     tpat - time taken in getting particular solution
c     tlpadap - time taken in layer potential computation in adaptive
c                integration
c     tlpspr - time taken to spread layer potential to all targets
c
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
      real *8 tpat,tlpadap,tlpspr,t1,t2,omp_get_wtime
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

      do ii = 1,npol3
         do jj = npt+1,ntarg0
            tab(jj,ii) = zero
         enddo
      enddo
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      tlpadap = t2-t1


      call cpu_time(t1)
C$      t1 = omp_get_wtime()      

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
      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      tlpspr = t2-t1

      return
      end
