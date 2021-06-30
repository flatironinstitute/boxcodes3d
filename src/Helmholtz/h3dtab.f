cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c     TABLE GENERATION UTILITIES
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This file contains the routines for evaluating 
c     the near-field volume potential using an 
c     anti-Helmholtzian
c      
c
c
c
c
c      
c      
c
      
      subroutine h3dtabp_ref(ndeg,zk,tol,npt,x,tab,ldtab,nup,
     1   iflg,pmat_qr,npol3,pmat_jpvt,pmat_tau,tpat,tlpadap,tlpspr)
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
c     nup - anti helmholtzian upward recurrence max
c     iflg - flag for determining which version of adaptive
c              integration to use
c            if iflg = 1, one with precomputation is used
c            if iflg = 2, one without precomputation is used
c     pmat_qr - qr decomposition of values to coefs matrix
c         computed using dgeqp3
c     npol3 - number of polynomials of total degree <= ndeg
c       and also the tailing dimension of pmat_qr
c     pmat_jpvt - permuation matrix in pivoted qr deocomposition
c     pmat_tau - scaling factors in householde decompositions of
c       qr decomposition of the matrix
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
      integer ntarg0, ntarg, n, idim, istart, npt, itype,nup
      integer iend,istart2,iend2
      real *8 slicevals(6), flipd(6), flips(6), abszk, abszktol
      real *8 rcond, val, derscale, tol2, r, theta, phi
      real *8 tpat,tlpadap,tlpspr,t1,t2,omp_get_wtime
      real *8 x(3,npt),pmat_qr(npt,npol3),pmat_tau(npol3)
      integer pmat_jpvt(npol3)
      parameter (abszktol = 2.5d0)

      real *8, allocatable :: w(:), pols(:), v(:,:)
      real *8, allocatable :: grid(:,:)
      real *8, allocatable :: vmat1d(:,:),x1d(:)
      real *8 v1d,w1d,u1d
      real *8 u, pi4
      integer ldu, ldv
      integer iflg
      
      
      complex *16 zero, im, one
      complex *16, allocatable :: tabtemp(:,:), ahc(:,:), zv(:,:)
      complex *16, allocatable :: tabtemp2(:,:)
      complex *16, allocatable :: ahcleg(:,:)
      complex *16, allocatable :: ahdercleg(:,:)
      complex *16, allocatable :: ahcleg3(:,:)
      complex *16, allocatable :: slp_pots(:,:), dlp_pots(:,:)

      real *8, allocatable :: errsup(:),errsdown(:)

      integer i,npolout,j,ipt
      real *8 derrmax

      
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
      
      ndeg2 = ndeg + 2*nup
      call legetens_npol_3d(ndeg2,type,npolout)
      call legetens_npol_2d(ndeg2,type,npol2)

      allocate(pols(npolout))


      allocate(ahcleg3(npolout,npol3))
      allocate(ahcleg(npol2,npol3),ahdercleg(npol2,npol3))

      allocate(errsup(npol3),errsdown(npol3))

      ntarg0 = 10*npt
      allocate(tabtemp2(ntarg0,npol3))


      call h3danti_form(ndeg,nup,type,zk,ahcleg3,npolout,derrmax,
     1   errsup,errsdown)
      call prin2('derrmax=*',derrmax,1)
      call prin2('errsup=*',errsup,12)
      call prin2('errsdown=*',errsdown,12)
c
c  get 1d interpolation matrix
c
c
      itype = 0
      allocate(x1d(n))
      call legeexps(itype,n,x1d,u1d,u1d,w1d)

      allocate(vmat1d(ndeg2+1,n))

      do i=1,n
        call legepols(x1d(i),ndeg2,vmat1d(1,i))
      enddo

c
c  note vectorized initialization
c

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ipt,j,pols)
      do ipt=1,npt
        call legetens_pols_3d(x(1,ipt),ndeg2,type,pols)
        do i=1,npol3
          tabtemp2(ipt,i) = 0
          do j=1,npolout
            tabtemp2(ipt,i) = tabtemp2(ipt,i) + ahcleg3(j,i)*pols(j)
          enddo
          tabtemp2(ipt,i) = -tabtemp2(ipt,i)*pi4
        enddo
      enddo
C$OMP END PARALLEL DO      
      

      call cpu_time(t2)
C$       t2 = omp_get_wtime()      

      tpat = t2-t1
c     memory depending on ntarg

      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      ntarg = 6*ntarg0

      allocate(slp_pots(npol2,ntarg),dlp_pots(npol2,ntarg))
      allocate(tabtemp(ntarg0,npol3))

      call h3d_facelayerpot_eval(tol,zk,ndeg,ndeg2,type,iflg,
     1     npt,x,slp_pots,dlp_pots,npol2)

      do ii = 1,npol3
         do jj = npt+1,ntarg0
            tabtemp2(jj,ii) = zero
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
          ifdiff = 1
          call legetens_slicezcoeffs_3d(ndeg2,type,idim,val,
     1         ahcleg3,npolout,npol3,ahcleg,npol2,ifdiff,ahdercleg,
     2         npol2)
          do ii = 1,npol3
             do jj = 1,npol2
                ahdercleg(jj,ii) = derscale*ahdercleg(jj,ii)
             enddo
          enddo
         

         istart = (iface-1)*ntarg0 + 1
         call zgemm('T','N',ntarg0,npol3,npol2,one,
     1        slp_pots(1,istart),
     1        npol2,ahdercleg,npol2,zero,tabtemp,ntarg0)

         do ii = 1,npol3
            do jj = 1,ntarg0
               tabtemp2(jj,ii) = tabtemp2(jj,ii) +
     1              flips(iface)*tabtemp(jj,ii)
            enddo
         enddo

         istart = (iface-1)*ntarg0 + 1 
         call zgemm('T','N',ntarg0,npol3,npol2,one,
     1        dlp_pots(1,istart),
     1        npol2,ahcleg,npol2,zero,tabtemp,ntarg0)

         do ii = 1,npol3
            do jj = 1,ntarg0
               tabtemp2(jj,ii) = tabtemp2(jj,ii) +
     1              flipd(iface)*tabtemp(jj,ii)
            enddo
         enddo
         
      enddo         
      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      tlpspr = t2-t1

c
c  tabtemp2 now stores the volume potential at all the polya points. Convert 
c  the volume potential to coefficients now
c
c
      do i=1,10
        istart = (i-1)*npt + 1
        iend = i*npt

        istart2 = (i-1)*npol3 + 1
        iend2 = i*npol3
        call zqrsolv(npt,npol3,pmat_qr,pmat_jpvt,pmat_tau,npol3,
     1    tabtemp2(istart:iend,1:npol3),tab(istart2:iend2,1:npol3))
      enddo

      return
      end
c
c
c     
c
c
      subroutine h3d_facelayerpot_eval(tol,zk,ndeg,ndegp,type,
     1  iflg,nppbox,grid,slp_pots,dlp_pots,lda)
      implicit real *8 (a-h,o-z)
      real *8 tol
      complex *16 zk
      integer ndeg,ndegp,lda,ntarg0
      real *8 grid(3,nppbox)
      complex *16 slp_pots(lda,*),dlp_pots(lda,*)
      complex *16, allocatable :: slp_near(:,:),dlp_near(:,:)
      complex *16, allocatable :: slp_far(:,:),dlp_far(:,:)
      character type

      real *8, allocatable :: xyztarg(:,:),xyztarg_near(:,:)
      real *8, allocatable :: xyztarg_far(:,:)
      integer, allocatable :: ifar_ind(:),inear_ind(:)
      real *8, allocatable :: xyztmp(:,:),qnodes(:,:),qwts(:)
      real *8 xq(ndeg+1),xyzc(3),u,v,w
      integer ipars(10)
      real *8 dpars(5)
      complex *16 zpars(3)

      external h3d_slp,h3d_dlp

      integer norder, norder_p, ncores

      done = 1
      pi = atan(done)*4

      norder = ndeg+1

      norder_p = ndegp+1

      call legetens_npol_2d(ndegp,type,npols)
cc      call prinf('npols=*',npols,1)
c     npols = norder_p*norder_p

      bs = 2.0d0
      xyzc(1) = -1
      xyzc(2) = -1
      xyzc(3) = -1 


      itype = 0
      call legeexps(itype,norder,xq,u,v,w)
      do i=1,norder
        xq(i) = xq(i) + 1
      enddo
c
      ntarg0 = 10*nppbox
      allocate(xyztmp(3,ntarg0))

      istart1 = 4*nppbox+1
      istart2 = 7*nppbox+1
      
      call tensrefpts3d_grid(nppbox,grid,bs,xyzc,xyztmp,
     1   xyztmp(1,istart1),xyztmp(1,istart2))

      ntarg = 6*ntarg0
      allocate(xyztarg(3,ntarg),xyztarg_near(3,ntarg))
      allocate(xyztarg_far(3,ntarg),inear_ind(ntarg))
      allocate(ifar_ind(ntarg))

      do i=1,ntarg0
        xyztarg(1,i+0*ntarg0) = xyztmp(2,i)
        xyztarg(2,i+0*ntarg0) = xyztmp(3,i)
        xyztarg(3,i+0*ntarg0) = xyztmp(1,i)+1

        xyztarg(1,i+1*ntarg0) = xyztmp(2,i)
        xyztarg(2,i+1*ntarg0) = xyztmp(3,i)
        xyztarg(3,i+1*ntarg0) = xyztmp(1,i)-1

        xyztarg(1,i+2*ntarg0) = xyztmp(1,i)
        xyztarg(2,i+2*ntarg0) = xyztmp(3,i)
        xyztarg(3,i+2*ntarg0) = xyztmp(2,i)+1
        
        xyztarg(1,i+3*ntarg0) = xyztmp(1,i)
        xyztarg(2,i+3*ntarg0) = xyztmp(3,i)
        xyztarg(3,i+3*ntarg0) = xyztmp(2,i)-1

        xyztarg(1,i+4*ntarg0) = xyztmp(1,i)
        xyztarg(2,i+4*ntarg0) = xyztmp(2,i)
        xyztarg(3,i+4*ntarg0) = xyztmp(3,i)+1

        xyztarg(1,i+5*ntarg0) = xyztmp(1,i)
        xyztarg(2,i+5*ntarg0) = xyztmp(2,i)
        xyztarg(3,i+5*ntarg0) = xyztmp(3,i)-1
      enddo
      

      nquadmax = 8000
      nqorder = 20
      eps = tol
      call h3d_get_eps_nqorder_nqmax(tol,norder,eps,nqorder,nquadmax,
     1        nqorderf)
      
      intype = 2
cc      call prinf("Starting adap quad for near*",i,0)


      zpars(1) = zk

      if(iflg.eq.1) then
      
        ntarg_f = 0
        ntarg_n = 0
      
        znear = 0.6d0
        xynear = 1.6d0

        do i=1,ntarg
          x=xyztarg(1,i)
          y=xyztarg(2,i)
          z=xyztarg(3,i)
          if((abs(z).le.znear.and.abs(x).le.xynear.
     1        and.abs(y).le.xynear)) then
            ntarg_n = ntarg_n + 1
            xyztarg_near(1,ntarg_n) = xyztarg(1,i)
            xyztarg_near(2,ntarg_n) = xyztarg(2,i)
            xyztarg_near(3,ntarg_n) = xyztarg(3,i)
            inear_ind(ntarg_n) = i
          else
            ntarg_f = ntarg_f + 1
            xyztarg_far(1,ntarg_f) = xyztarg(1,i)
            xyztarg_far(2,ntarg_f) = xyztarg(2,i)
            xyztarg_far(3,ntarg_f) = xyztarg(3,i)
            ifar_ind(ntarg_f) = i
          endif
        enddo

        allocate(slp_near(npols,ntarg_n),dlp_near(npols,ntarg_n))
        allocate(slp_far(npols,ntarg_f),dlp_far(npols,ntarg_f))

        do i=1,ntarg_n
          do j=1,npols
            slp_near(j,i) = 0
            dlp_near(j,i) = 0
          enddo
        enddo

        do i=1,ntarg_f
          do j=1,npols
            slp_far(j,i) = 0
            dlp_far(j,i) = 0
          enddo
        enddo

        nquadmax = 8000
        nqorder = 20
        eps = tol
        call h3d_get_eps_nqorder_nqmax(tol,norder,eps,nqorder,nquadmax,
     1        nqorderf)
      
        intype = 2



        zpars(1) = zk

        nbatches = 48
        nttpcore = ceiling((ntarg_n+0.0d0)/nbatches)
        ntt = 1

        t1 = second()
C$       t1 = omp_get_wtime()  
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,istart,iend,ntt)
c$OMP& SCHEDULE(DYNAMIC)      
        do i=1,nbatches
          istart = (i-1)*nttpcore+1
          iend = min(i*nttpcore,ntarg_n)
          ntt = iend-istart+1

          call cquadints_adap(eps,intype,norder_p,type,npols,ntt,
     1       xyztarg_near(1,istart),nquadmax,h3d_slp,dpars,zpars,ipars,
     2       nqorder,slp_near(1,istart))

          call cquadints_adap(eps,intype,norder_p,type,npols,ntt,
     1       xyztarg_near(1,istart),nquadmax,h3d_dlp,dpars,zpars,ipars,
     2       nqorder,dlp_near(1,istart))
        enddo
C$OMP END PARALLEL DO      
        t2 = second()
C$       t2 = omp_get_wtime()      

cc      call prin2('time taken in evaluating near=*',t2-t1,1)


        call squarearbq_pts(nqorderf,nnodes)

        nu = 3
        nqpts = nnodes*nu*nu
        allocate(qnodes(2,nqpts),qwts(nqpts))


        call gen_xg_uniftree_nodes(nqorderf,nnodes,nu,nqpts,qnodes,qwts)

        call cquadints_wnodes(norder_p,type,npols,ntarg_f,xyztarg_far,
     1       h3d_slp,dpars,zpars,ipars,nqpts,qnodes,qwts,slp_far)

        call cquadints_wnodes(norder_p,type,npols,ntarg_f,xyztarg_far,
     1        h3d_dlp,dpars,zpars,ipars,nqpts,qnodes,qwts,dlp_far)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
        do i=1,ntarg_n
          do j=1,npols
            slp_pots(j,inear_ind(i)) = slp_near(j,i)
            dlp_pots(j,inear_ind(i)) = dlp_near(j,i)
          enddo
        enddo
C$OMP END PARALLEL DO


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
        do i=1,ntarg_f
          do j=1,npols
            slp_pots(j,ifar_ind(i)) = slp_far(j,i)
            dlp_pots(j,ifar_ind(i)) = dlp_far(j,i)
          enddo
        enddo
C$OMP END PARALLEL DO
      endif


      if(iflg.eq.2) then
        ntt = 1

        t1 = second()
C$       t1 = omp_get_wtime()  
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
c$OMP& SCHEDULE(DYNAMIC)      
        do i=1,ntarg

          call cquadints_adap2(eps,intype,norder_p,type,npols,ntt,
     1     xyztarg(1,i),nquadmax,h3d_slp,dpars,zpars,ipars,
     2     nqorder,slp_pots(1,i))

          call cquadints_adap2(eps,intype,norder_p,type,npols,ntt,
     1     xyztarg(1,i),nquadmax,h3d_dlp,dpars,zpars,ipars,
     2     nqorder,dlp_pots(1,i))
        enddo
C$OMP END PARALLEL DO      
        t2 = second()
C$       t2 = omp_get_wtime()      
      endif


      return
      end
c
c
c
c     
c
c
      subroutine h3d_facelayergrad_eval_new(tol,zk,ndeg,ndegp,type,
     1  iflg,slp_pots,slp_gradx,slp_grady,slp_gradz,dlp_pots,
     2  dlp_gradx,dlp_grady,dlp_gradz,lda)
      implicit real *8 (a-h,o-z)
      real *8 tol
      complex *16 zk
      integer ndeg,ndegp,lda
      complex *16 slp_pots(lda,*),slp_gradx(lda,*)
      complex *16 slp_grady(lda,*),slp_gradz(lda,*)
      complex *16 dlp_pots(lda,*),dlp_gradx(lda,*)
      complex *16 dlp_grady(lda,*),dlp_gradz(lda,*)
      complex *16, allocatable :: slp_near(:,:),slp_far(:,:)
      complex *16, allocatable :: slp_gradx_near(:,:),slp_gradx_far(:,:)
      complex *16, allocatable :: slp_grady_near(:,:),slp_grady_far(:,:)
      complex *16, allocatable :: slp_gradz_near(:,:),slp_gradz_far(:,:)
      complex *16, allocatable :: dlp_near(:,:),dlp_far(:,:)
      complex *16, allocatable :: dlp_gradx_near(:,:),dlp_gradx_far(:,:)
      complex *16, allocatable :: dlp_grady_near(:,:),dlp_grady_far(:,:)
      complex *16, allocatable :: dlp_gradz_near(:,:),dlp_gradz_far(:,:)
      character type

      real *8, allocatable :: xyztarg(:,:),xyztarg_near(:,:)
      real *8, allocatable :: xyztarg_far(:,:)
      integer, allocatable :: ifar_ind(:),inear_ind(:)
      real *8, allocatable :: xyztmp(:,:),qnodes(:,:),qwts(:)
      real *8 xq(ndeg+1),xyzc(3),u,v,w
      integer ipars(10)
      real *8 dpars(5)
      complex *16 zpars(3)

      external h3d_slp,h3d_sgradx,h3d_sgrady,h3d_sgradz
      external h3d_dlp,h3d_dgradx,h3d_dgrady,h3d_dgradz

      integer norder, norder_p, ncores

      done = 1
      pi = atan(done)*4

      norder = ndeg+1

      norder_p = ndegp+1

      call legetens_npol_2d(ndegp,type,npols)
cc      call prinf('npols=*',npols,1)
c     npols = norder_p*norder_p

      bs = 2.0d0
      xyzc(1) = -1
      xyzc(2) = -1
      xyzc(3) = -1 


      itype = 0
      call legeexps(itype,norder,xq,u,v,w)
      do i=1,norder
        xq(i) = xq(i) + 1
      enddo
c
      nppbox = norder*norder*norder
      ntarg0 = 10*nppbox
      allocate(xyztmp(3,ntarg0))

      istart1 = 4*nppbox+1
      istart2 = 7*nppbox+1
      
      call tensrefpts3d(xq,norder,bs,xyzc,xyztmp,xyztmp(1,istart1),
     1        xyztmp(1,istart2))

      ntarg = 6*ntarg0
      allocate(xyztarg(3,ntarg),xyztarg_near(3,ntarg))
      allocate(xyztarg_far(3,ntarg),inear_ind(ntarg))
      allocate(ifar_ind(ntarg))

      do i=1,ntarg0
        xyztarg(1,i+0*ntarg0) = xyztmp(2,i)
        xyztarg(2,i+0*ntarg0) = xyztmp(3,i)
        xyztarg(3,i+0*ntarg0) = xyztmp(1,i)+1

        xyztarg(1,i+1*ntarg0) = xyztmp(2,i)
        xyztarg(2,i+1*ntarg0) = xyztmp(3,i)
        xyztarg(3,i+1*ntarg0) = xyztmp(1,i)-1

        xyztarg(1,i+2*ntarg0) = xyztmp(1,i)
        xyztarg(2,i+2*ntarg0) = xyztmp(3,i)
        xyztarg(3,i+2*ntarg0) = xyztmp(2,i)+1
        
        xyztarg(1,i+3*ntarg0) = xyztmp(1,i)
        xyztarg(2,i+3*ntarg0) = xyztmp(3,i)
        xyztarg(3,i+3*ntarg0) = xyztmp(2,i)-1

        xyztarg(1,i+4*ntarg0) = xyztmp(1,i)
        xyztarg(2,i+4*ntarg0) = xyztmp(2,i)
        xyztarg(3,i+4*ntarg0) = xyztmp(3,i)+1

        xyztarg(1,i+5*ntarg0) = xyztmp(1,i)
        xyztarg(2,i+5*ntarg0) = xyztmp(2,i)
        xyztarg(3,i+5*ntarg0) = xyztmp(3,i)-1
      enddo
      

      nquadmax = 8000
      nqorder = 20
      eps = tol
      call h3d_get_eps_nqorder_nqmax(tol,norder,eps,nqorder,nquadmax,
     1        nqorderf)
      
      intype = 2
cc      call prinf("Starting adap quad for near*",i,0)


      zpars(1) = zk

      if(iflg.eq.1) then
      
        ntarg_f = 0
        ntarg_n = 0
      
        znear = 0.6d0
        xynear = 1.6d0

        do i=1,ntarg
          x=xyztarg(1,i)
          y=xyztarg(2,i)
          z=xyztarg(3,i)
          if((abs(z).le.znear.and.abs(x).le.xynear.
     1        and.abs(y).le.xynear)) then
            ntarg_n = ntarg_n + 1
            xyztarg_near(1,ntarg_n) = xyztarg(1,i)
            xyztarg_near(2,ntarg_n) = xyztarg(2,i)
            xyztarg_near(3,ntarg_n) = xyztarg(3,i)
            inear_ind(ntarg_n) = i
          else
            ntarg_f = ntarg_f + 1
            xyztarg_far(1,ntarg_f) = xyztarg(1,i)
            xyztarg_far(2,ntarg_f) = xyztarg(2,i)
            xyztarg_far(3,ntarg_f) = xyztarg(3,i)
            ifar_ind(ntarg_f) = i
          endif
        enddo

        allocate(slp_near(npols,ntarg_n),slp_gradx_near(npols,ntarg_n))
        allocate(slp_far(npols,ntarg_f),slp_gradx_far(npols,ntarg_f))
        allocate(slp_grady_near(npols,ntarg_n))
        allocate(slp_grady_far(npols,ntarg_f))
        allocate(slp_gradz_near(npols,ntarg_n))
        allocate(slp_gradz_far(npols,ntarg_f))

        allocate(dlp_near(npols,ntarg_n),dlp_gradx_near(npols,ntarg_n))
        allocate(dlp_far(npols,ntarg_f),dlp_gradx_far(npols,ntarg_f))
        allocate(dlp_grady_near(npols,ntarg_n))
        allocate(dlp_grady_far(npols,ntarg_f))
        allocate(dlp_gradz_near(npols,ntarg_n))
        allocate(dlp_gradz_far(npols,ntarg_f))

        do i=1,ntarg_n
          do j=1,npols
            slp_near(j,i) = 0
            slp_gradx_near(j,i) = 0
            slp_grady_near(j,i) = 0
            slp_gradz_near(j,i) = 0

            dlp_near(j,i) = 0
            dlp_gradx_near(j,i) = 0
            dlp_grady_near(j,i) = 0
            dlp_gradz_near(j,i) = 0
          enddo
        enddo

        do i=1,ntarg_f
          do j=1,npols
            slp_far(j,i) = 0
            slp_gradx_far(j,i) = 0
            slp_grady_far(j,i) = 0
            slp_gradz_far(j,i) = 0

            dlp_far(j,i) = 0
            dlp_gradx_far(j,i) = 0
            dlp_grady_far(j,i) = 0
            dlp_gradz_far(j,i) = 0
          enddo
        enddo

        nquadmax = 8000
        nqorder = 20
        eps = tol
        call h3d_get_eps_nqorder_nqmax(tol,norder,eps,nqorder,nquadmax,
     1        nqorderf)
      
        intype = 2



        zpars(1) = zk

        nbatches = 48
        nttpcore = ceiling((ntarg_n+0.0d0)/nbatches)
        ntt = 1

        t1 = second()
C$       t1 = omp_get_wtime()  
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,istart,iend,ntt)
c$OMP& SCHEDULE(DYNAMIC)      
        do i=1,nbatches
          istart = (i-1)*nttpcore+1
          iend = min(i*nttpcore,ntarg_n)
          ntt = iend-istart+1

          call cquadints_adap(eps,intype,norder_p,type,npols,ntt,
     1       xyztarg_near(1,istart),nquadmax,h3d_slp,dpars,zpars,ipars,
     2       nqorder,slp_near(1,istart))

          call cquadints_adap(eps,intype,norder_p,type,npols,ntt,
     1       xyztarg_near(1,istart),nquadmax,h3d_sgradx,dpars,zpars,
     2       ipars,nqorder,slp_gradx_near(1,istart))

          call cquadints_adap(eps,intype,norder_p,type,npols,ntt,
     1       xyztarg_near(1,istart),nquadmax,h3d_sgrady,dpars,zpars,
     2       ipars,nqorder,slp_grady_near(1,istart))

          call cquadints_adap(eps,intype,norder_p,type,npols,ntt,
     1       xyztarg_near(1,istart),nquadmax,h3d_sgradz,dpars,zpars,
     2       ipars,nqorder,slp_gradz_near(1,istart))

          call cquadints_adap(eps,intype,norder_p,type,npols,ntt,
     1       xyztarg_near(1,istart),nquadmax,h3d_dlp,dpars,zpars,ipars,
     2       nqorder,dlp_near(1,istart))

          call cquadints_adap(eps,intype,norder_p,type,npols,ntt,
     1       xyztarg_near(1,istart),nquadmax,h3d_dgradx,dpars,zpars,
     2       ipars,nqorder,dlp_gradx_near(1,istart))

          call cquadints_adap(eps,intype,norder_p,type,npols,ntt,
     1       xyztarg_near(1,istart),nquadmax,h3d_dgrady,dpars,zpars,
     2       ipars,nqorder,dlp_grady_near(1,istart))

          call cquadints_adap(eps,intype,norder_p,type,npols,ntt,
     1       xyztarg_near(1,istart),nquadmax,h3d_dgradz,dpars,zpars,
     2       ipars,nqorder,dlp_gradz_near(1,istart))

        enddo
C$OMP END PARALLEL DO      
        t2 = second()
C$       t2 = omp_get_wtime()      

cc      call prin2('time taken in evaluating near=*',t2-t1,1)


        call squarearbq_pts(nqorderf,nnodes)

        nu = 3
        nqpts = nnodes*nu*nu
        allocate(qnodes(2,nqpts),qwts(nqpts))


        call gen_xg_uniftree_nodes(nqorderf,nnodes,nu,nqpts,qnodes,qwts)

        call cquadints_wnodes(norder_p,type,npols,ntarg_f,xyztarg_far,
     1      h3d_slp,dpars,zpars,ipars,nqpts,qnodes,qwts,slp_far)

        call cquadints_wnodes(norder_p,type,npols,ntarg_f,xyztarg_far,
     1   h3d_sgradx,dpars,zpars,ipars,nqpts,qnodes,qwts,slp_gradx_far)

        call cquadints_wnodes(norder_p,type,npols,ntarg_f,xyztarg_far,
     1   h3d_sgrady,dpars,zpars,ipars,nqpts,qnodes,qwts,slp_grady_far)

        call cquadints_wnodes(norder_p,type,npols,ntarg_f,xyztarg_far,
     1   h3d_sgradz,dpars,zpars,ipars,nqpts,qnodes,qwts,slp_gradz_far)

        call cquadints_wnodes(norder_p,type,npols,ntarg_f,xyztarg_far,
     1      h3d_dlp,dpars,zpars,ipars,nqpts,qnodes,qwts,dlp_far)

        call cquadints_wnodes(norder_p,type,npols,ntarg_f,xyztarg_far,
     1   h3d_dgradx,dpars,zpars,ipars,nqpts,qnodes,qwts,dlp_gradx_far)

        call cquadints_wnodes(norder_p,type,npols,ntarg_f,xyztarg_far,
     1   h3d_dgrady,dpars,zpars,ipars,nqpts,qnodes,qwts,dlp_grady_far)

        call cquadints_wnodes(norder_p,type,npols,ntarg_f,xyztarg_far,
     1   h3d_dgradz,dpars,zpars,ipars,nqpts,qnodes,qwts,dlp_gradz_far)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
        do i=1,ntarg_n
          do j=1,npols
            slp_pots(j,inear_ind(i)) = slp_near(j,i)
            slp_gradx(j,inear_ind(i)) = slp_gradx_near(j,i)
            slp_grady(j,inear_ind(i)) = slp_grady_near(j,i)
            slp_gradz(j,inear_ind(i)) = slp_gradz_near(j,i)
            dlp_pots(j,inear_ind(i)) = dlp_near(j,i)
            dlp_gradx(j,inear_ind(i)) = dlp_gradx_near(j,i)
            dlp_grady(j,inear_ind(i)) = dlp_grady_near(j,i)
            dlp_gradz(j,inear_ind(i)) = dlp_gradz_near(j,i)
          enddo
        enddo
C$OMP END PARALLEL DO


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
        do i=1,ntarg_f
          do j=1,npols
            slp_pots(j,ifar_ind(i)) = slp_far(j,i)
            slp_gradx(j,ifar_ind(i)) = slp_gradx_far(j,i)
            slp_grady(j,ifar_ind(i)) = slp_grady_far(j,i)
            slp_gradz(j,ifar_ind(i)) = slp_gradz_far(j,i)
            dlp_pots(j,ifar_ind(i)) = dlp_far(j,i)
            dlp_gradx(j,ifar_ind(i)) = dlp_gradx_far(j,i)
            dlp_grady(j,ifar_ind(i)) = dlp_grady_far(j,i)
            dlp_gradz(j,ifar_ind(i)) = dlp_gradz_far(j,i)
          enddo
        enddo
C$OMP END PARALLEL DO
      endif


      if(iflg.eq.2) then
        ntt = 1

        t1 = second()
C$       t1 = omp_get_wtime()  
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
c$OMP& SCHEDULE(DYNAMIC)      
        do i=1,ntarg

          call cquadints_adap2(eps,intype,norder_p,type,npols,ntt,
     1     xyztarg(1,i),nquadmax,h3d_slp,dpars,zpars,ipars,
     2     nqorder,slp_pots(1,i))

          call cquadints_adap2(eps,intype,norder_p,type,npols,ntt,
     1     xyztarg(1,i),nquadmax,h3d_sgradx,dpars,zpars,ipars,
     2     nqorder,slp_gradx(1,i))

          call cquadints_adap2(eps,intype,norder_p,type,npols,ntt,
     1     xyztarg(1,i),nquadmax,h3d_sgrady,dpars,zpars,ipars,
     2     nqorder,slp_grady(1,i))

          call cquadints_adap2(eps,intype,norder_p,type,npols,ntt,
     1     xyztarg(1,i),nquadmax,h3d_sgradz,dpars,zpars,ipars,
     2     nqorder,slp_gradz(1,i))

          call cquadints_adap2(eps,intype,norder_p,type,npols,ntt,
     1     xyztarg(1,i),nquadmax,h3d_dlp,dpars,zpars,ipars,
     2     nqorder,dlp_pots(1,i))

          call cquadints_adap2(eps,intype,norder_p,type,npols,ntt,
     1     xyztarg(1,i),nquadmax,h3d_dgradx,dpars,zpars,ipars,
     2     nqorder,dlp_gradx(1,i))

          call cquadints_adap2(eps,intype,norder_p,type,npols,ntt,
     1     xyztarg(1,i),nquadmax,h3d_dgrady,dpars,zpars,ipars,
     2     nqorder,dlp_grady(1,i))

          call cquadints_adap2(eps,intype,norder_p,type,npols,ntt,
     1     xyztarg(1,i),nquadmax,h3d_dgradz,dpars,zpars,ipars,
     2     nqorder,dlp_gradz(1,i))

        enddo
C$OMP END PARALLEL DO      
        t2 = second()
C$       t2 = omp_get_wtime()      
      endif


      return
      end
c
c
c
c
c
c
      subroutine h3d_get_eps_nqorder_nqmax(tol,norder,eps,nqorder,
     1   nqmax,nqorderf)
      implicit none
      real *8 tol,eps
      integer norder,nqorder,nqmax,iprec,nqorderf
c
c
c        fix this routine to optimize performance
c

      iprec = 0
      eps = tol
      nqorder = 10

      if(tol.lt.0.49d-2) iprec = 1
      if(tol.lt.0.49d-3) iprec = 2
      if(tol.lt.0.49d-6) iprec = 3
      if(tol.lt.0.49d-9) iprec = 4

cc      print *, "iprec=",iprec

      if(iprec.eq.0) then
c
c   norder, max(err/max(true,1)), max(err/true),max(err)/max(true)
c       4, 0.1e-2, 0.7e-2, 0.1e-2 
c       6, 0.1e-2, 0.8, 0.1e-2
c       8, 0.5e-3, 0.8e+2, 0.3e-3
c      12, 0.1e-2, 0.1e+5, 0.6e-3 
c
         nqmax = 1500
         nqorderf = 10
         eps = 0.5d-2
         if(norder.le.4) then
           nqorder = 7
         else if(norder.gt.4.and.norder.le.6) then
           nqmax = 2000
           eps = 0.5d-2
           nqorder = 12
           nqorderf = 12
         else if(norder.gt.6.and.norder.le.8) then
           nqorder = 14
           eps = 0.1d-2
           nqmax = 5000
           nqorderf = 12
         else if(norder.gt.8) then
           nqorder = 20
           eps = 0.5d-3
           nqmax = 12000
           nqorderf = 20
         endif

      endif

      if(iprec.eq.1) then
c
c   norder, max(err/max(true,1)), max(err/true),max(err)/max(true)
c       4, 0.2e-4, 0.7e-2, 0.1e-4 
c       6, 0.8e-4, 0.2, 0.6e-4
c       8, 0.1e-4, 0.8e+2, 0.8e-4
c      12, 0.2e-3, 0.1e+5, 0.1e-3 
c
         nqmax = 1500
         nqorderf = 10
         eps = 0.3d-2
         if(norder.le.4) then
           nqorder = 7
         else if(norder.gt.4.and.norder.le.6) then
           nqmax = 2000
           eps = 0.3d-2
           nqorder = 12
           nqorderf = 12
         else if(norder.gt.6.and.norder.le.8) then
           nqorder = 14
           eps = 0.6d-3
           nqmax = 5000
           nqorderf = 12
         else if(norder.gt.8) then
           nqorder = 20
           eps = 0.1d-3
           nqmax = 12000
           nqorderf = 20
         endif

      endif

      if(iprec.eq.2) then
c
c   norder, max(err/max(true,1)), max(err/true),max(err)/max(true)
c       4, 0.3e-7, 0.2e-5, 0.3e-7 
c       6, 0.3e-6, 0.3e-3, 0.2e-6
c       8, 0.5e-6, 0.8e-1, 0.3e-6
c      12, 0.5e-6, 6.6e3, 0.3e-6 
c
         nqmax = 5000
         nqorderf = 16
         eps = 0.3d-4
         if(norder.le.4) then
           nqorder = 12
         else if(norder.gt.4.and.norder.le.6) then
           nqorder = 14
         else if(norder.gt.6.and.norder.le.8) then
           nqorder = 14
           eps = 0.5d-5
           nqmax = 7000
         else if(norder.gt.8) then
           nqorder = 20
           eps = 0.1d-5
           nqmax = 20000
           nqorderf = 20
         endif

      endif

      if(iprec.eq.3) then
c
c   norder, max(err/max(true,1)), max(err/true),max(err)/max(true)
c       4, 0.4e-10, 0.3e-7, 0.3e-10 
c       6, 0.2e-9, 0.1e-5, 0.1e-9
c       8, 0.6e-10, 0.8e-3, 0.4e-10
c      12, 0.4e-10, 0.1e5, 0.3e-10 
c
         nqorder = 20
         eps = 3.0d-7
         nqmax = 4000
         nqorderf = 20
         if(norder.le.4) then
         else if(norder.gt.4.and.norder.le.6) then
           nqorderf = 22
         else if(norder.gt.6.and.norder.le.8) then
           nqmax = 7000
           eps = 3.0d-8
           nqorderf = 24
         else if(norder.gt.8) then
           nqorder = 24
           eps = 3.0d-9
           nqmax = 25000
           nqorderf = 30
         endif
      endif

      if(iprec.eq.4) then
c
c   NOT TESTED for accuracy
c
c
c   norder, max(err/max(true,1)), max(err/true),max(err)/max(true)
c
         nqorder = 24
         eps = 3.0d-10
         nqmax = 6000
         nqorderf = 24
         if(norder.le.4) then
         else if(norder.gt.4.and.norder.le.6) then
           nqorderf = 26
         else if(norder.gt.6.and.norder.le.8) then
           nqmax = 11000
           eps = 3.0d-11
           nqorderf = 28
         else if(norder.gt.8) then
           nqorder = 28
           eps = 3.0d-12
           nqmax = 30000
           nqorderf = 30
         endif
        
      endif



      return
      end

c
c
c
c
c
c
c
      subroutine h3d_get_nqorder_far(tol,nqorder)
      implicit none
      real *8 tol,eps
      integer norder,nqorder
c
c
c        fix this routine to optimize performance
c
      nqorder = 12
      if(tol.le.0.5d-3) nqorder = 16
      if(tol.le.0.5d-6) nqorder = 20
      if(tol.le.0.5d-9) nqorder = 25
      if(tol.le.0.5d-12) nqorder = 30
      return
      end

c
c
c
c
c

      subroutine h3d_slp(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)

      f = exp(ima*zpars(1)*rr)/rr

      return
      end

c      
c      
c
c
c
c
c

      subroutine h3d_dlp(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,z

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)
      z = ima*zpars(1)*rr

      f = exp(z)*y(3)*(z-1.0d0)/rr**3

      return
      end

c
c
c
c     
c
c
      subroutine h3d_sgradx(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,z

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)
      z = ima*zpars(1)*rr

      f = exp(z)*(y(1)-x(1))*(z-1.0d0)/rr**3

      return
      end

c
c
c
c     
c
c
      subroutine h3d_sgrady(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,z

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)
      z = ima*zpars(1)*rr

      f = exp(z)*(y(2)-x(2))*(z-1.0d0)/rr**3

      return
      end

c
c
c
c
c     
c
c
      subroutine h3d_sgradz(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,z

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)
      z = ima*zpars(1)*rr

      f = exp(z)*y(3)*(z-1.0d0)/rr**3

      return
      end

c
c     
c
c
      subroutine h3d_dgradx(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,z

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)
      z = ima*zpars(1)*rr

      f = -exp(z)*(y(1)-x(1))*y(3)*(-3.0d0 + 3*z - z**2)/rr**5

      return
      end

c
c
c
c     
c
c
      subroutine h3d_dgrady(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,z

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)
      z = ima*zpars(1)*rr

      f = -exp(z)*(y(2)-x(2))*y(3)*(-3.0d0 + 3*z - z**2)/rr**5

      return
      end

c
c
c
c
c     
c
c
      subroutine h3d_dgradz(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,z

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)
      z = ima*zpars(1)*rr
      f=-exp(z)*((1.0d0-z)*rr**2+y(3)*y(3)*(-3.0d0 + 3*z - z**2))/rr**5
      


      return
      end

c
c
c
