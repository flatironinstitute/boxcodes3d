      implicit real *8 (a-h,o-z)
      real *8 dpars(100)
      integer iptr(9), ipars(1)
      integer, allocatable :: itree(:)
      real *8, allocatable :: fvals(:,:,:),centers(:,:),boxsize(:)
      real *8, allocatable :: umat(:,:),vmat(:,:),xref(:,:),wts(:)
      real *8, allocatable :: pts(:,:), errstest(:,:)
      real *8, allocatable :: pts1(:,:), pts2(:,:), pts3(:,:)      

      complex *16, allocatable :: qval(:,:),rhsval(:,:),uval(:,:)
      real *8 xyztmp(3),rintl(0:200)
      real *8 timeinfo(6),tprecomp(3)
      complex *16 zk,zpars(4)

      complex *16, allocatable :: pot(:,:),soln(:,:)
      complex *16, allocatable :: potcoef(:,:), pottest(:)
      complex *16, allocatable :: pottestref(:), pot1(:), pot2(:),
     1     pot3(:)
      complex *16, allocatable :: rhscoef(:,:), rhs1(:), rhs2(:),
     1     rhs3(:)
      complex *16, allocatable :: qcoef(:,:), q1(:), q2(:),
     1     q3(:)
      
      complex *16 ima,zz,ztmp, zbeam(3)

      real *8 alpha,beta
      real *8 errs(1000), eps_step(20), f3(3)

      character *1, ttype

      character *30 :: fextra
      character *80 :: fbase, fout1, fout2, fout3, fout4, fout5, fpre
      character *80 :: fout6, fout7
      
      character *3 :: snwave
      
      data ima/(0.0d0,1.0d0)/

      external ftree_eaton_and_beam_ls
      logical flag

      integer icontinue


      call prini(6,13)

      done = 1
      pi = atan(done)*4

      do i = 1,10
         eps_step(i) = 1d-1**(i+1)
      enddo
      

c     set up output file names
      
      fbase = "ls_eaton_"

      write(*,*) 'ENTER PREFIX: '
      read(*,'(A)') fpre
      write(*,*) 'YOU ENTERED: ', fpre

      write(*,*) 'ENTER FILENAME MODIFIER (E.G. DATETIME) MAX 30 CHAR: '
      read(*,'(A)') fextra

      write(*,*) 'ENTER NUMBER OF WAVE LENGTHS ON A SIDE (NWAVE): '
      read(*,*) nwave


      write(*,*) 'ENTER NUMBER OF CONVERGENCE STEPS (STARTS 1e-2): '
      read(*,*) nstep
      
      epsmin = 1d1
      do i = 1,nstep
         epsmin = min(epsmin,eps_step(i))
      enddo
      

      write(snwave,20) nwave
 20   format(I3.3)

      fbase =
     1   trim(fpre)//trim(fbase)//trim(fextra)//'_nwave_'//trim(snwave)
      
      write(*,*) 'FILENAME BASE: ', fbase

c
c      initialize problem parameters
c


c     size and wavelength

      boxlen = 1.0d0
      zk = 2*pi*nwave/boxlen

c     beam parameters

c     coordinate direction
      ibeam = 1
      zbeam(1) = -boxlen/2 - 0.01d0 + ima*0.5d0
c
      zbeam(2) = 0.27d0*boxlen
      zbeam(3) = 0.01d0

c     load into zpars, ipars

      ipars(1) = ibeam
      zpars(1) = zbeam(1)
      zpars(2) = zbeam(2)
      zpars(3) = zbeam(3)

      zpars(4) = zk


      
c     eaton lens parameters

c     support
      alpha = 0.45d0

c     q cut-off
      beta = 2.0d0

c     precompute eaton lens q evaluator structure
      ldpars = 100
      nlege = 40
      call eaton_pre(alpha,beta,dpars,ldpars,nlege)



c     number of random points to test

      ntest = 200
      iseed = 12345
      aaa = hkrand(iseed)
      allocate(pts(3,ntest),errstest(ntest,nstep),pottest(ntest),
     1     pottestref(ntest))
      do i = 1,200
         pts(1,i) = -boxlen/2.0d0 + boxlen*hkrand(0)
         pts(2,i) = -boxlen/2.0d0 + boxlen*hkrand(0)
         pts(3,i) = -boxlen/2.0d0 + boxlen*hkrand(0)
      enddo


c     plotting parameters

      nplot = 500
      nplot2 = nplot**2
      allocate(pts1(3,nplot2),pts2(3,nplot2),pts3(3,nplot2))
      allocate(pot1(nplot2),pot2(nplot2),pot3(nplot2))
      allocate(rhs1(nplot2),rhs2(nplot2),rhs3(nplot2))
      allocate(q1(nplot2),q2(nplot2),q3(nplot2))

      h = boxlen/nplot
      do i = 1,nplot
         do j = 1,nplot
            ind = (i-1)*nplot + j

c     xy plane slice
            pts1(1,ind) = -boxlen/2 + h/2 + (i-1)*h
            pts1(2,ind) = -boxlen/2 + h/2 + (j-1)*h
            pts1(3,ind) = 0
      
c     xz plane slice
            pts2(1,ind) = -boxlen/2 + h/2 + (i-1)*h
            pts2(3,ind) = -boxlen/2 + h/2 + (j-1)*h
            pts2(2,ind) = 0
      
c     yz plane slice
            pts3(2,ind) = -boxlen/2 + h/2 + (i-1)*h
            pts3(3,ind) = -boxlen/2 + h/2 + (j-1)*h
            pts3(1,ind) = 0

         enddo
      enddo
      
      
c     discretization parameters
      
      norder = 4
      iptype = 1
      eta = 2
      ttype = 'T'


c     solver parameters

      numit = max(500,2*nwave**2)
      irep = 1
      eps_gmres = 1.0d-7
      


      zkeff = zk*boxlen

      npbox = norder*norder*norder


      write(*,*) 'ENTER REFERENCE EPS: (SMALLER THAN ',
     1     epsmin,')'
      read(*,*) eps
      write(*,*) 'you entered ', eps
      if (eps .gt. epsmin) then
         write(*,*) 'error, ref should be higher accuracy '
         write(*,*) 'BAILING ...'
         stop
      endif
      
      call cpu_time(t1)
C$      t1 = omp_get_wtime()

      nd = 3
      call vol_tree_mem(eps,zk,boxlen,norder,iptype,eta,
     1     ftree_eaton_and_beam_ls,nd,dpars,zpars,ipars,nlevels,
     2     nboxes,ltree,rintl)

      
      write(*,*) 'nboxes= ',nboxes
      write(*,*) 'nlevels= ',nlevels

      allocate(fvals(nd,npbox,nboxes),centers(3,nboxes))
      allocate(boxsize(0:nlevels),itree(ltree))

      write(*,*) 'dpars= ',dpars(1:2)
      write(*,*) 'zk= ',zk
      write(*,*) 'nd= ',nd

      write(*,*) 'eps= ',eps
      write(*,*) 'eta= ',eta
      write(*,*) 'iptype= ',iptype
      write(*,*) 'rintl= ',rintl(0:nlevels)

      call vol_tree_build(eps,zk,boxlen,norder,iptype,eta,
     1     ftree_eaton_and_beam_ls,nd,dpars,zpars,ipars,nlevels,
     2     nboxes,ltree,rintl,itree,iptr,fvals,centers,boxsize)
      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      

      print *, "done building tree"
      allocate(qval(npbox,nboxes),uval(npbox,nboxes),
     1    rhsval(npbox,nboxes))
      do i=1,nboxes
        do j=1,npbox
          qval(j,i) = fvals(1,j,i)
          rhsval(j,i) = -zk**2*qval(j,i)*
     1         (fvals(2,j,i) + ima*fvals(3,j,i))
        enddo
      enddo

      write(*,*) 'time taken to build tree= ',t2-t1
      write(*,*) 'speed in points per sec= ',
     1     (nboxes*norder**3+0.0d0)/(t2-t1)

      write(*,*) 'iptr *',iptr(1:8)
      write(*,*) 'ltree *',ltree
      write(*,*) 'nboxes *', nboxes
      write(*,*) 'norder *',norder


      
      write(*,*) 'FOR REFERENCE SOLUTION ------- '
      write(*,*) 'NBOXES: ', nboxes
      write(*,*) 'NPTS:   ', nboxes*npbox
      write(*,*) 'CONTINUE? (ENTER 1 for yes, OTHER NUMBER for no)'
      read(*,*) icontinue

      if (icontinue .ne. 1) then
         write(*,*) 'STOPPING... NO FILES CREATED.'
         stop
      endif



      write(*,*) 'OK, HERE WE GO.'


      write(*,*) 'OUTPUT WILL BE WRITTEN TO THE FILES: '

      fout1 = trim(fbase)//'_pars.txt'
      fout2 = trim(fbase)//'_errs.txt'
      fout3 = trim(fbase)//'_prin.txt'
      fout4 = trim(fbase)//'_pltin.txt'
      fout5 = trim(fbase)//'_pltsc.txt'
      fout6 = trim(fbase)//'_pltq.txt'
      fout7 = trim(fbase)//'_errv.txt'            

      write(*,*) '    ', fout1
      write(*,*) '    ', fout2
      write(*,*) '    ', fout3
      write(*,*) '    ', fout4
      write(*,*) '    ', fout5
      write(*,*) '    ', fout6
      write(*,*) '    ', fout7
      
      open(unit=101,file=fout1)
      open(unit=102,file=fout2)
      open(unit=13,file=fout3)
      open(unit=104,file=fout4)
      open(unit=105,file=fout5)
      open(unit=106,file=fout6)
      open(unit=107,file=fout7)      

c     print incident wave to file for plotting
      
      ndeg = norder-1
      call legetens_npol_3d(ndeg,ttype,ncbox)
      write(*,*) 'NCBOX ',ncbox
      allocate(rhscoef(ncbox,nboxes),qcoef(ncbox,nboxes))

      write(*,*) 'printing incident wave plot and q plot to file ...'
      
      ndz=2
      call vol_tree_coef(nboxes,norder,ttype,rhsval,ndz,npbox,rhscoef,
     1     ncbox)
      call vol_tree_coef(nboxes,norder,ttype,qval,ndz,npbox,qcoef,
     1     ncbox)
      write(*,*) 'after coefs'
      call vol_tree_eval(norder,ttype,ndz,ncbox,rhscoef,nlevels,nboxes,
     1     itree,iptr,centers,boxsize,pts1,nplot2,rhs1)

      call vol_tree_eval(norder,ttype,ndz,ncbox,rhscoef,nlevels,nboxes,
     1     itree,iptr,centers,boxsize,pts2,nplot2,rhs2)
      
      call vol_tree_eval(norder,ttype,ndz,ncbox,rhscoef,nlevels,nboxes,
     1     itree,iptr,centers,boxsize,pts3,nplot2,rhs3)

      call vol_tree_eval(norder,ttype,ndz,ncbox,qcoef,nlevels,nboxes,
     1     itree,iptr,centers,boxsize,pts1,nplot2,q1)

      call vol_tree_eval(norder,ttype,ndz,ncbox,qcoef,nlevels,nboxes,
     1     itree,iptr,centers,boxsize,pts2,nplot2,q2)
      
      call vol_tree_eval(norder,ttype,ndz,ncbox,qcoef,nlevels,nboxes,
     1     itree,iptr,centers,boxsize,pts3,nplot2,q3)

      
      iunit = 104
      call write3dsliceplotz(iunit,pts1,rhs1,pts2,rhs2,
     1     pts3,rhs3,nplot2)
      close(104)

      iunit = 106
      call write3dsliceplotz(iunit,pts1,q1,pts2,q2,
     1     pts3,q3,nplot2)
      close(106)

c     print params to file
      
      call prini(6,101)
      write(*,*) 'writing parameters to file '

      call prinf('domain and wavelength params -------- *',ier,0)
      call prinf('nwave *',nwave,1)
      call prin2('boxlen *',boxlen,1)
      call prin2('zk *',zk,2)
      
      call prinf('beam params ------------------------- *',ier,0)
      call prin2('zbeam *',zbeam,6)
      call prinf('ibeam *',ibeam,1)

      call prinf('lens params ------------------------- *',ier,0)
      call prin2('alpha *',alpha,1)
      call prin2('beta *',beta,1)
      call prinf('nlege *',nlege,1)
      call prin2('dpars *',dpars,ldpars)

      call prinf('discretization params --------------- *',ier,0)
      call prinf('norder *',norder,1)
      call prinf('iptype *',iptype,1)
      call prin2('eta *',eta,1)
      call prina('ttype *',ttype,1)

      call prinf('solver params ----------------------- *',ier,0)
      call prinf('numit *',numit,1)
      call prinf('irep *',irep,1)
      call prin2('eps_gmres *',eps_gmres,1)

      call prin2('convergence params ------------------ *',ier,0)
      call prinf('nstep *',nstep,1)
      call prin2('eps per step *',eps_step,nstep)
      call prinf('iseed *',iseed,1)
      call prinf('ntest *',ntest,1)
      call prin2('points tested *',pts,3*ntest)
      
      call prin2('plotting params --------------------- *',ier,0)
      call prinf('nplot *',nplot,1)
      
      call prini(6,13)

      call prini(6,102)
      call prin2('reference solution ------------------ *',ier,0)
      call prin2('eps for reference solution *',eps,1)
      call prin2( 'time taken to build tree= *',t2-t1,1)
      call prin2( 'speed in points per sec= *',
     1        (nboxes*norder**3+0.0d0)/(t2-t1),1)
      call prinf('nboxes *',nboxes,1)
      call prinf('nlevels *',nlevels,1)
      call prinf('ltree *',nlevels,1)      
      call prinf('iptr *',iptr,8)      
      call prinf('npbox *',npbox,1)
      call prinf('num points *',nboxes*npbox,1)
      call prini(6,13)
      
c
c
c     Compute reference solution
c

      call prinf('COMPUTING REFERENCE SOLUTION ----*',ier,0)
      
      npols = norder*(norder+1)*(norder+2)/6


      allocate(soln(npbox,nboxes),pot(npbox,nboxes))


      call cpu_time(t1) 
C$    t1 = omp_get_wtime()      

      niter=0
      eps_bicgs = eps_gmres
      restart_tol = 1d-6
      call ls_solver_guru_bicgs(eps,zk,nboxes,nlevels,ltree,
     1     itree,iptr,norder,npols,ttype,qval,centers,boxsize,npbox,
     2     rhsval,irep,eps_bicgs,restart_tol,numit,niter,errs,rres,soln)


c      call ls_solver_guru(eps,zk,nboxes,nlevels,ltree,itree,iptr,
c     1   norder,npols,ttype,qval,centers,boxsize,npbox,
c     2   rhsval,irep,eps_gmres,numit,niter,errs,rres,soln)
      
      call cpu_time(t2) 
C$    t2 = omp_get_wtime()      


      call prinf('...DONE COMPUTING REFERENCE SOLUTION ----*',ier,0)
      
      call prinf('niter=*',niter,1)
      call prin2('errs=*',errs,niter)
      call prin2('time for solve *',t2-t1,1)
      call prin2('time per fmm *',(t2-t1)/niter,1)

      call prini(6,102)
      call prinf('niter=*',niter,1)
      call prin2('errs=*',errs,niter)
      call prin2('rres=*',rres,1)      
      call prin2('time for solve *',t2-t1,1)
      call prin2('time per fmm *',(t2-t1)/niter,1)
      call prini(6,13)

      if(irep.eq.1) then
         call cpu_time(t1) 
C$       t1 = omp_get_wtime()      
         call helmholtz_volume_fmm(eps,zk,nboxes,nlevels,ltree,itree,
     1        iptr,norder,npols,ttype,soln,centers,boxsize,npbox,
     2        pot,timeinfo,tprecomp)
         call cpu_time(t2) 
C$       t2 = omp_get_wtime()

         call prini(6,102)
         call prin2('time taken in fmm (eval)=*',t2-t1,1)
         call prin2('timeinfo fmm (eval)=*',timeinfo,6)
         call prin2('tprecomp fmm (eval)=*',tprecomp,3)
         
         call prini(6,13)
      endif

      if(irep.eq.2) then
        do ibox=1,nboxes
          do j=1,npbox
            pot(j,ibox) = soln(j,ibox)
          enddo
        enddo
      endif


c     print scattered wave to file for plotting
      
      ndeg = norder-1
      call legetens_npol_3d(ndeg,ttype,ncbox)
      write(*,*) 'NCBOX ',ncbox
      allocate(potcoef(ncbox,nboxes))

      write(*,*) 'printing incident wave plot to file ...'
      
      ndz=2
      call vol_tree_coef(nboxes,norder,ttype,pot,ndz,npbox,potcoef,
     1     ncbox)
      write(*,*) 'after coefs'
      call vol_tree_eval(norder,ttype,ndz,ncbox,potcoef,nlevels,nboxes,
     1     itree,iptr,centers,boxsize,pts1,nplot2,pot1)

      call vol_tree_eval(norder,ttype,ndz,ncbox,potcoef,nlevels,nboxes,
     1     itree,iptr,centers,boxsize,pts2,nplot2,pot2)
      
      call vol_tree_eval(norder,ttype,ndz,ncbox,potcoef,nlevels,nboxes,
     1     itree,iptr,centers,boxsize,pts3,nplot2,pot3)

      
      iunit = 105
      call write3dsliceplotz(iunit,pts1,pot1,pts2,pot2,
     1     pts3,pot3,nplot2)

      close(105)

c     evaluate at reference points

      call vol_tree_eval(norder,ttype,ndz,ncbox,potcoef,nlevels,nboxes,
     1     itree,iptr,centers,boxsize,pts,ntest,pottestref)
      

      do iii = 1,nstep

         call prini(6,102)
         call prinf('testing convergence, step = --------- *',iii,1)
         call prin2('eps_step = *',eps_step(iii),1)
         call prini(6,13)
         
         call prinf('building tree ...*',ier,0)

         call cpu_time(t1)
C$       t1 = omp_get_wtime()

         nd = 3
         eps1 = eps_step(iii)
         call vol_tree_mem(eps1,zk,boxlen,norder,iptype,eta,
     1        ftree_eaton_and_beam_ls,nd,dpars,zpars,ipars,nlevels1,
     2        nboxes1,ltree1,rintl)

      
         call prinf( 'nboxes1=*',nboxes1,1)
         call prinf( 'nlevels1=*',nlevels1,1)
         call prinf( 'ltree1=*',ltree1,1)
         
         call prin2( 'eps1=*',eps1,1)
         call prin2( 'eta=*',eta,1)
         call prinf( 'iptype=*',iptype,1)
         call prin2( 'rintl=*',rintl,nlevels1+1)

         if (nboxes1 .gt. nboxes) then
            write(*,*) 'more boxes for easier problem? bomb'
            stop
         endif
         if (nlevels1 .gt. nlevels) then
            write(*,*) 'more levels for easier problem? bomb'
            stop
         endif
         if (ltree1 .gt. ltree) then
            write(*,*) 'bigger tree for easier problem? bomb'
            stop
         endif
         
         call vol_tree_build(eps1,zk,boxlen,norder,iptype,eta,
     1        ftree_eaton_and_beam_ls,nd,dpars,zpars,ipars,nlevels1,
     2        nboxes1,ltree1,rintl,itree,iptr,fvals,centers,boxsize)
         
         call cpu_time(t2)
C$       t2 = omp_get_wtime()      

         call prinf('done building tree *',ier,0)
         do i=1,nboxes1
            do j=1,npbox
               qval(j,i) = fvals(1,j,i)
               rhsval(j,i) = -zk**2*qval(j,i)*
     1              (fvals(2,j,i) + ima*fvals(3,j,i))
            enddo
         enddo


         call prini(6,102)
         call prin2( 'time taken to build tree= *',t2-t1,1)
         call prin2( 'speed in points per sec= *',
     1        (nboxes1*norder**3+0.0d0)/(t2-t1),1)
         call prinf('nboxes1 *', nboxes1,1)
         call prinf('nlevels1 *', nlevels1,1)
         call prinf('ltree1 *',ltree1,1)
         call prinf('iptr *',iptr,8)
         call prinf('npbox *',npbox,1)
         call prinf('num points *',nboxes1*npbox,1)
         call prini(6,13)


         call prinf('COMPUTING SOLUTION, STEP ----*',iii,1)
         
         npols = norder*(norder+1)*(norder+2)/6
         
         niter = 0
         call cpu_time(t1) 
C$       t1 = omp_get_wtime()      

         call ls_solver_guru_bicgs(eps1,zk,nboxes1,nlevels1,ltree1,
     1        itree,iptr,norder,npols,ttype,qval,centers,boxsize,npbox,
     2        rhsval,irep,eps_bicgs,restart_tol,numit,niter,errs,rres,
     3        soln)
         
c         call ls_solver_guru(eps1,zk,nboxes1,nlevels1,ltree1,itree,iptr,
c     1        norder,npols,ttype,qval,centers,boxsize,npbox,
c     2        rhsval,irep,eps_gmres,numit,niter,errs,rres,soln)

         call cpu_time(t2) 
C$       t2 = omp_get_wtime()      
         
         call prinf('...DONE COMPUTING SOLUTION, STEP ----*',iii,1)
         
         call prinf('niter=*',niter,1)
         call prin2('errs=*',errs,niter)
         call prin2('time for solve *',t2-t1,1)
         call prin2('time per fmm *',(t2-t1)/niter,1)

         call prini(6,102)
         call prinf('niter=*',niter,1)
         call prin2('errs=*',errs,niter)
         call prin2('rres=*',rres,1)         
         call prin2('time for solve *',t2-t1,1)
         call prin2('time per fmm *',(t2-t1)/niter,1)
         call prini(6,13)

         do ibox=1,nboxes1
            do j=1,npbox
               pot(j,ibox) = 0.0d0
            enddo
         enddo

         if(irep.eq.1) then
            call cpu_time(t1) 
C$          t1 = omp_get_wtime()      
            call helmholtz_volume_fmm(eps1,zk,nboxes1,nlevels1,ltree1,
     1           itree,iptr,norder,npols,ttype,soln,centers,boxsize,
     2           npbox,pot,timeinfo,tprecomp)
            call cpu_time(t2) 
C$          t2 = omp_get_wtime()

            call prini(6,102)
            call prin2('time taken in fmm (eval)=*',t2-t1,1)
            call prin2('timeinfo fmm (eval)=*',timeinfo,6)
            call prin2('tprecomp fmm (eval)=*',tprecomp,3)
            call prini(6,13)
         endif


         
         if(irep.eq.2) then
            do ibox=1,nboxes1
               do j=1,npbox
                  pot(j,ibox) = soln(j,ibox)
               enddo
            enddo
         endif

c     evaluate at reference points

         ndz = 2
         call vol_tree_coef(nboxes1,norder,ttype,pot,ndz,npbox,potcoef,
     1        ncbox)
         call vol_tree_eval(norder,ttype,ndz,ncbox,potcoef,nlevels1,
     1        nboxes1,itree,iptr,centers,boxsize,pts,ntest,pottest)

         s1 = 0
         s2 = 0
         do jjj = 1,ntest
            errstest(jjj,iii) = abs(pottest(jjj)-pottestref(jjj))
            s1 = s1 + errstest(jjj,iii)
            s2 = s2 + abs(pottestref(jjj))

            write(107,30) errstest(jjj,iii), real(pottest(jjj)),
     1           -real(ima*pottest(jjj)), real(pottestref(jjj)),
     2           -real(ima*pottestref(jjj))
 30         format(' ',5E12.4)
         enddo
         
         call prini(6,102)
         call prin2('pottestref *',pottestref,2*ntest)
         call prin2('pottest *',pottest,2*ntest)
         call prin2('errstest *',errstest(1,iii),ntest)
         call prin2('relerr sum *',s1/s2,1)
         call prini(6,13)
            

      enddo      


      write(*,*) '... FINISHED. OUTPUT WRITTEN TO THE FILES: '
      
      write(*,*) '    ', fout1
      write(*,*) '    ', fout2
      write(*,*) '    ', fout3
      write(*,*) '    ', fout4
      write(*,*) '    ', fout5
      write(*,*) '    ', fout6
      write(*,*) '    ', fout7      


            
      
      close(101)
      close(102)
      close(13)
      close(107)
      
      stop
      end
c
c
c
c 
      subroutine ftree_eaton_and_beam_ls(nd,xyz,dpars,zpars,ipars,f)
c
c     compute the value of q for the eaton lens (f(1)), the real and
c     imaginary parts of the incoming gaussian beam (f(2),f(3))
c
c     eaton lens parameters are in dpars
c        dpars(*) fast eaton evaluator structure
c
c     gaussian beam parameters are in zpars(1:3)
c        zpars(1:3) - complex center for beam
c        zpars(4) - zk, helmholtz wave number
c        ipars(1) - beam direction
c      
      implicit none
      integer nd,ipars(*)
      complex *16 zpars(*)
      integer ii
      real *8 dpars(*),f(*),xyz(3),r2,r,yi
      complex *16 zdiff(3), z2, z, zk, ima, zf
      data ima / (0.0d0,1.0d0) /

      r2 = xyz(1)**2 + xyz(2)**2 + xyz(3)**2
      r = sqrt(r2)

      call eaton_post(r,dpars,f(1))

c      f(1) = 1.0d0

      zdiff(1) = xyz(1) - zpars(1)
      zdiff(2) = xyz(2) - zpars(2) 
      zdiff(3) = xyz(3) - zpars(3)

      z2 = zdiff(1)**2 + zdiff(2)**2 + zdiff(3)**2
      z = sqrt(z2)
      
      zk = zpars(4)

      ii = ipars(1)
      yi = -real(ima*zpars(ii))
      
      zf = exp(ima*zk*z)/z*exp(-zk*yi)

      f(2) = real(zf)
      f(3) = -real(ima*zf)

      return
      end

c
c
c
c 

      subroutine write3dsliceplotz(iunit,pts1,vals1,pts2,
     1     vals2,pts3,vals3,npts)
      implicit real *8 (a-h,o-z)
      real *8 pts1(3,*), pts2(3,*), pts3(3,*)
      real *8 vals1(2,*), vals2(2,*), vals3(2,*)

      do i = 1,npts
         write(iunit,1001) pts1(1,i), pts1(2,i), pts1(3,i), vals1(1,i),
     1        vals1(2,i)
 1001    format(' ',5E12.4)
      enddo

      do i = 1,npts
         write(iunit,1001) pts2(1,i), pts2(2,i), pts2(3,i), vals2(1,i),
     1        vals2(2,i)
 1002    format(' ',5E12.4)
      enddo
      do i = 1,npts
         write(iunit,1001) pts3(1,i), pts3(2,i), pts3(3,i), vals3(1,i),
     1        vals3(2,i)
 1003    format(' ',5E12.4)
      enddo
      
      return
      end
      

      

