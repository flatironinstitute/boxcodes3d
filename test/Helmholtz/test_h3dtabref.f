      implicit real *8 (a-h,o-z)
      real *8 tts(3,100,2),errm_rel1(100,2),errm_rel2(100,2)
      real *8 errm_abs(100,2)
      real *8 tt3d(3),err_print(3)
      real *8, allocatable :: x(:,:)
      complex *16 zk, im, zero, one
      complex *16, allocatable :: tab_ref(:,:),tab(:,:),tab2(:,:)
      integer, allocatable :: ipt2depth(:,:),idepth2pt(:,:,:)
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

      nup = 8

      type = 't'
      call legetens_npol_3d(ndeg,type,npol3)

      if(iprec.eq.0) tol = 1.0d-2
      if(iprec.eq.1) tol = 1.0d-3
      if(iprec.eq.2) tol = 1.0d-6
      if(iprec.eq.3) tol = 1.0d-9
      if(iprec.eq.4) tol = 1.0d-12

      nptmax = n**3
      allocate(x(3,nptmax))
      
      allocate(ipt2depth(3,nptmax),idepth2pt(n,n,n))

      ifloor = 1.75d0*npol3
      call fakepolya3d(n,ifloor,npt,x,ipt2depth,idepth2pt)
      itype = 0
      call prinf('npt=*',npt,1)
      call prin2('rat=*',(npt+0.0d0)/(nptmax+0.0d0),1)
      call prin2('x=*',x,3*npt)

      
cc      call legetens_exps_3d(itype,n,'t',x,u,1,v,1,w)
      

      ntarg0 = 10*npt



      allocate(tab_ref(ntarg0,npol3),tab(ntarg0,npol3))
      allocate(tab2(ntarg0,npol3))

c
c
c
c
      ifgen = 1
      izk = 2
      eps_exact = 1.0d-8
      ifwrite = 0
      ifread = 0
      if(ifgen.eq.1) then
         if(izk.eq.1) zk = 1.6d0
         if(izk.eq.2) zk = 0.16d0
         if(izk.eq.3) zk = im*1.6d0
         if(izk.eq.4) zk = 5.0d0
         print *, "here"
         print *, "zk=",zk
         call h3dtab_ref_brute2(ndeg,eps_exact,x,npt,tab_ref,ntarg0,
     1      h3d_vslp,
     1      dpars,zk,ipars)
         if(ifwrite.eq.1) then
          write(fname,'(a,i2.2,a,i1,a)') 'tabref_data/tabref_',
     1       n,'_izk_',izk,'new.dat'
          open(unit=33,file=trim(fname),action='readwrite',
     1       form='unformatted',access='stream')
         write(unit=33) tab_ref
         close(33)
         stop
        endif
      endif


      

      if(igen.ne.1) then

        do izk=2,2

          if(izk.eq.1) zk = 1.6d0
          if(izk.eq.2) zk = 0.16d0
          if(izk.eq.3) zk = im*1.6d0
          if(izk.eq.4) zk = 5.0d0
          if(ifread.eq.1) then
            write(fname,'(a,i2.2,a,i1,a)') 'tabref_data/tabref_',n,
     1         '_izk_',izk,'new.dat'
            open(unit=33,file=trim(fname),action='readwrite',
     1         form='unformatted',access='stream')
            read(unit=33) tab_ref

            print *, "done reading"

            close(33)
          endif

          ndeg2 = ndeg + 2*nup 
          iflg = 1

          call h3dtabp_ref2(ndeg,zk,tol,npt,x,tab,ntarg0,nup,
     1      iflg,tts(1,izk,iflg),tts(2,izk,iflg),tts(3,izk,iflg))
          
          
          iflg = 2

          call h3dtabp_ref2(ndeg,zk,tol,npt,x,tab,ntarg0,nup,
     1      iflg,tts(1,izk,iflg),tts(2,izk,iflg),tts(3,izk,iflg))
          print *, "Done with 2d adap quad"
        
          call prin2('tts iflg1=*',tts(1,izk,1),3)
          call prin2('tts iflg2=*',tts(1,izk,2),3)

          call cpu_time(t1)
C$           t1 = omp_get_wtime()          

          call h3dtab_ref_brute2(ndeg,tol,x,npt,tab2,ntarg0,
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
            do j=1,ntarg0
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
            do j=1,ntarg0
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

        call prin2('max errors relative 1=*',errm_rel1(2,1:2),2)
        call prin2('max errors relative 2=*',errm_rel2(2,1:2),2)
        call prin2('max errors absolute=*',errm_abs(2,1:2),2)
        call prin2('time taken 2d=*',tts(:,2,:),6)
        call prin2('time taken 3d adap=*',tt3d(2),1)
        stop
      endif



      
      stop
      end
c
c      
c      
c
      
      subroutine h3dtabp_ref2(ndeg,zk,tol,npt,x,tab,ldtab,nup,
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
c     nup - anti helmholtzian upward recurrence max
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
      integer ntarg0, ntarg, n, idim, istart, npt, itype,nup
      real *8 slicevals(6), flipd(6), flips(6), abszk, abszktol
      real *8 rcond, val, derscale, tol2, r, theta, phi
      real *8 tpat,tlpadap,tlpspr,t1,t2,omp_get_wtime
      real *8 x(3,npt)
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
      call legetens_npol_3d(ndeg,type,npol3)
      call legetens_npol_3d(ndeg2,type,npolout)
      call legetens_npol_2d(ndeg2,type,npol2)

      allocate(pols(npolout))


      allocate(ahcleg3(npolout,npol3))
      allocate(ahcleg(npol2,npol3),ahdercleg(npol2,npol3))

      allocate(errsup(npol3),errsdown(npol3))


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
          tab(ipt,i) = 0
          do j=1,npolout
            tab(ipt,i) = tab(ipt,i) + ahcleg3(j,i)*pols(j)
          enddo
          tab(ipt,i) = -tab(ipt,i)*pi4
        enddo
      enddo
C$OMP END PARALLEL DO      
      

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

      call h3d_facelayerpot_eval_new2(tol,zk,ndeg,ndeg2,type,iflg,
     1     npt,x,slp_pots,dlp_pots,npol2)

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
c
c
c     
c
c
      subroutine h3d_facelayerpot_eval_new2(tol,zk,ndeg,ndegp,type,
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

      subroutine tensrefpts3d_grid(ngrid,grid,bs,xyzc,wc,wbtos,wstob)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     generate reference points for the limited subset of
c     points in colleague, big-to-small, and small-to-big
c     interactions which can be used to obtain the other
c     interactions by symmetries, where points on unit box
c     are specified
c
      
      implicit real *8 (a-h,o-z)
      dimension wc(3,ngrid,4), wbtos(3,ngrid,3), wstob(3,ngrid,3)
      dimension xshift(3), yshift(3), zshift(3), xyzc(3)
      dimension grid(3,ngrid)

c     lowest corner of cube (each coordinate is smallest)
      
      xc = xyzc(1)
      yc = xyzc(2)
      zc = xyzc(3)      

c     get corresponding meshgrid

      do i = 1,ngrid
         wc(1,i,1) = xc + grid(1,i)+1
         wc(2,i,1) = yc + grid(2,i)+1
         wc(3,i,1) = zc + grid(3,i)+1
      enddo


      xshift(1) = -1
      xshift(2) = 0
      xshift(3) = 0
      yshift(1) = -1
      yshift(2) = -1
      yshift(3) = 0
      zshift(1) = -1
      zshift(2) = -1
      zshift(3) = -1

      bsh = bs/2
      
      do ii = 1,3
         do i = 1,ngrid
            wc(1,i,ii+1) = xc + xshift(ii)*bs + grid(1,i)+1
            wc(2,i,ii+1) = yc + yshift(ii)*bs + grid(2,i)+1
            wc(3,i,ii+1) = zc + zshift(ii)*bs + grid(3,i)+1          

            wbtos(1,i,ii) = xc + xshift(ii)*bsh + (grid(1,i)+1)/2
            wbtos(2,i,ii) = yc + yshift(ii)*bsh + (grid(2,i)+1)/2
            wbtos(3,i,ii) = zc + zshift(ii)*bsh + (grid(3,i)+1)/2

            wstob(1,i,ii) = xc + (xshift(ii)-1)*bs + (grid(1,i)+1)*2
            wstob(2,i,ii) = yc + (yshift(ii)-1)*bs + (grid(2,i)+1)*2
            wstob(3,i,ii) = zc + (zshift(ii)-1)*bs + (grid(3,i)+1)*2

         enddo
      enddo
         


      return
      end


      
c
c  This file contains the following user callable subroutine:
c 
c     h3dtabp_ref_brute - compute reference quadrature 
c       tables using adaptive integration. The self
c       quadrature is handled by splitting the cube into
c       8 cubes while the rest are handled using 
c       standard adaptive integration
c
c
c
      subroutine h3dtab_ref_brute2(ndeg,tol,grid,npt,tab,ldtab,fker,
     1   dpars,zpars,ipars)

c
c
c
c     generate the quadrature table at the reference
c     points using adaptive integration where the kernel of integration
c     is given by fker
c
c     Calling sequence arguments for fker should be
c     
c     subroutine fker(x,y,dpars,zpars,ipars,f)
c
c     Note that this subroutine does not compute the integrals
c     for the self cube and doesn't do any intelligent split
c     between near and far
c
c
c     input
c
c     ndeg - integer, highest degree of the basis polynomials
c                    (measured in total degree, e.g. x^0y^1z^2
c                     has total degree 3)
c     ldtab - integer, leading dimension of the output table
c     fker - function handle for evaluating kernel
c     dpars - real parameters to be used by fker
c     zpars - complex parameters to be used by fker
c     ipars - integer parameters to be used by fker
c
c     output
c      
c     tab - complex *16 array (ldtab,*), tab(i,j) is the integral
c     of the j-th tensor polynomial (in the ordering specified
c     in legetens.f) against given kernel
c     at the i-th reference target point (see tensrefpts3d)
c     generate the Helmholtz potential table at the reference
c     points using adaptive integration
     
      implicit real *8 (a-h,o-z)
      integer ndeg,ldtab
      real *8 grid(3,npt)
      complex *16 tab(ldtab,*), zpars(*)
      real *8 dpars(*)
      integer ipars(*)
      

c       local      
      complex *16, allocatable :: tab_t(:,:),tab_tmp(:,:)
      real *8 xyzc(3),bs,xq(100),w(100)
      real *8, allocatable :: xyztarg(:,:)
      real *8, allocatable :: xyztmp(:,:)
      integer, allocatable :: iindtmp(:)
      character ptype

      external fker



      n = ndeg + 1
      norder = n

      call legetens_npol_3d(ndeg,'t',npols)

      bs = 2.0d0
      xyzc(1) = -1
      xyzc(2) = -1
      xyzc(3) = -1

      ntarg = 10*npt
      allocate(xyztarg(3,ntarg))
      allocate(xyztmp(3,ntarg),iindtmp(ntarg))

      istart1 = 4*npt+1
      istart2 = 7*npt+1
      call tensrefpts3d_grid(npt,grid,bs,xyzc,xyztarg,
     1  xyztarg(1,istart1),xyztarg(1,istart2))

cc      call prin2('xyztarg=*',xyztarg,3*npt)
      
      allocate(tab_t(npols,ntarg))
c
c       vector initialization
c
      tab_t = 0


      eps = tol 
      ncmax = 30000
      nqorder = 11

      call prin2('zk=*',zpars,2)

      ntt = 1
      
      call cpu_time(t1)
C$     t1 = omp_get_wtime()
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(DYNAMIC)
      do i=1,npt
        call ccubeints_split8int_adap(eps,norder,'t',npols,
     1        xyztarg(1,i),ntt,xyztarg(1,i),
     1        ncmax,fker,dpars,zpars,ipars,nqorder,tab_t(1,i))
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
C$OMP$SCHEDULE(DYNAMIC)
      do i=npt+1,ntarg
        call ccubeints_adap(eps,norder,'t',npols,ntt,xyztarg(1,i),
     1        ncmax,fker,dpars,zpars,ipars,nqorder,tab_t(1,i))
      enddo
C$OMP END PARALLEL DO

      call cpu_time(t2)
C$      t2 = omp_get_wtime()     
      
c
c    now transpose tab_t
c      
      do i=1,ntarg
        do j=1,npols
          tab(i,j) = tab_t(j,i) 
        enddo
      enddo

      return
      end



c
c
c
