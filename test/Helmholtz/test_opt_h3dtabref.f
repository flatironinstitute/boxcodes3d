      implicit real *8 (a-h,o-z)
      real *8 tts(3,3),errm_rel1(3),errm_rel2(3),errm_abs(3)
      complex *16 zk, im, zero, one
      complex *16, allocatable :: tab_ref(:,:),tab(:,:)
      logical ifsphere
      character type

      character *100 fname
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /            
      
      call prini(6,13)
      

      izk = 1




      ndeg = 11
      n = ndeg + 1

      type = 't'
      call legetens_npol_3d(ndeg,type,npol3)

      npt = n**3

      ntarg0 = 10*npt

      print *, ndeg,ntarg0,npol3

      tol = 1.0d-6


      allocate(tab_ref(ntarg0,npol3),tab(ntarg0,npol3))

c
c
c
c
      ifgen = 0
      if(ifgen.eq.1) then
         call h3dtabp_ref_brute(ndeg,zk,tab_ref,ntarg0)
         do i=1,npol3
           do j=1,ntarg0
              write(33,*) real(tab_ref(j,i)),imag(tab_ref(j,i))
           enddo
         enddo
         close(33)
         stop
      endif

      if(igen.ne.1) then

        do izk=1,3

          if(izk.eq.1) zk = 1.6d0
          if(izk.eq.2) zk = 0.16d0
          if(izk.eq.3) zk = im*1.6d0
          write(fname,'(a,i2.2,a,i1,a)') 'tabref_',n,'_izk_',izk,'.dat'
          open(unit=33,file=trim(fname))


          do i=1,npol3
            do j=1,ntarg0
              read(33,*) tmp1,tmp2
              tab_ref(j,i) = tmp1+im*tmp2
            enddo
          enddo

          ifsphere = .true.
          if(ifsphere) ndeg2 = 20
          if(.not.ifsphere) ndeg2 = ndeg

          call h3dtabp_ref2(ndeg,zk,tol,tab,ntarg0,ifsphere,ndeg2,
     1      tts(1,izk),tts(2,izk),tts(3,izk))
        
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

          call prin2('max error=*',errmax1,1)
          errm_rel1(izk) = errmax1
          errm_rel2(izk) = errmax2
          errm_abs(izk) = errmaxa/rmax
          call prin2('time taken=*',tts,3)
          close(33)
        enddo

        call prin2('max errors relative 1=*',errm_rel1,3)
        call prin2('max errors relative 2=*',errm_rel2,3)
        call prin2('max errors absolute=*',errm_abs,3)
        call prin2('time taken=*',tts,9)
        stop
      endif



      
      stop
      end
c
c      
c      
c
      
      subroutine h3dtabp_ref2(ndeg,zk,tol,tab,ldtab,ifsphere,ndeg2,
     1   tpat,tlpadap,tlpspr)
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

      call h3d_facelayerpot_eval(tol,zk,ndeg,ndeg2,type,slp_pots,
     1  dlp_pots,npol2)

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
