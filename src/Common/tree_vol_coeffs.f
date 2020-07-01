c
c
c    generate level restricted oct tree based on resolving a 
c    function to desired precision
c
c
c    The function handle is of the form
c    call fun(nd,xyz,dpars,zpars,ipars,f)
c
c      where xyz is the location in (-L/2,L/2)^3
c
c    dpars are a collection of real parameters, zpars are complex
c    parameters and ipars are integer parameters, and f is a 
c    real array of size nd
c     
c
c    For the Helmholtz/Maxwell tree, the boxes are refined until 
c     Re(zk)*boxsize<5, beyond which the function resolution criterion
c    kicks in. 
c
c    A function is said to be resolved if it's interpolant at the 8
c    children nodes, agrees with the function values at those nodes
c    upto the user specified tolerance. 
c    The error is scaled by h**(eta)
c    where eta is user specified and h is the boxsize. If there is 
c    any confusion, the user should seet \eta to 0
c    Let \tilde{f} denote the interpolant of f, then
c    the refinement criterion is 
c      \int_{B_{j} |\tilde{f}-f|^{p} *h^{\eta} < 
c        \varepsilon V_{j}^{1/p}/V_{0}^{1/p}/(\int_{B_{0}}|f|^{p})^{1/p}
c    This implies that
c       \int_{B_{0}} |\tilde{f}-f|^{p} =
c          \sum_{j} \int_{B_{j}} |\tilde{f}-f|^{p}
c          \leq \sum_{j} \eps^{p}*h^{\eta p} 
c                  V_{j}/V_{0}/(\int_{B_{0}} |f|^p)
c      If \eta = 0,
c          \leq \eps^{p}/(\int_{B_{0}} |f|^{p})
c
c    i.e., this strategy guarantees that the interpolated function
c      approximates the function with relative lp accuracy of \eps
c      
c    This code has 2 main user callable routines
c      make_vol_tree_mem -> Returns the memory requirements, 
c          tree length, number of boxes, number of levels
c      make_vol_tree -> Makes the actual tree, returns centers of boxes,
c          colleague info, function values on leaf boxes
c       
c          
c          iptr(1) - laddr
c          iptr(2) - ilevel
c          iptr(3) - iparent
c          iptr(4) - nchild
c          iptr(5) - ichild
c          iptr(6) - ncoll
c          iptr(7) - coll
c          iptr(8) - ltree
c




      subroutine vol_tree_mem(eps,zk,boxlen,norder,iptype,eta,fun,nd,
     1  dpars,zpars,ipars,nlevels,nboxes,ltree,rintl)
c
c      get memory requirements for the tree
c
c      input parameters:
c        eps - double precision
c           precision requested
c        zk - double complex
c           Helmholtz parameter
c        boxlen - double precision
c           length of box in which volume is contained, 
c           if boxlen = L, then volume is [-L/2,L/2]^3
c        norder - integer
c           order of discretization
c        eta - double precision
c           scaling parameter for error
c        fun - function handle
c           function to evalute it everywhere in the volume
c        nd - integer
c           number of real functions returned by fun
c        dpars - double precision
c           real parameters for function evaluation
c        zpars - double complex
c           complex parameters for function evaluation
c        ipars - integer
c           integer parameters for function evaluation
c 
c        output:
c           nlevels - integer
c             number of levels
c           nboxes - integer
c             number of boxes
c           nlboxes - integer
c             number of leaf boxes
c           ltree - integer
c             length of tree
c           rintl(0:nlevels) - real *8
c             lp norm to scale the functions by
c             (on input rintl should be of size(0:200)
c
c      
c

      implicit none
      real *8 eps,boxlen,eta,dpars(*)
      real *8, allocatable :: fvals(:,:,:)
      complex *16 zk,zpars(*)
      integer nd,ipars(*),iptype
      integer nlevels,nboxes,ltree,norder

      external fun

      integer, allocatable :: laddr(:,:),ilevel(:),iparent(:),nchild(:)
      integer, allocatable :: ichild(:,:),ncoll(:),icoll(:,:)
      real *8, allocatable :: centers(:,:)
      integer, allocatable :: nbors(:,:),nnbors(:)

      integer, allocatable :: ilevel2(:),iparent2(:),nchild2(:),
     1    ichild2(:,:)
      real *8, allocatable :: centers2(:,:),fvals2(:,:,:)

      integer nbmax,nlmax,npbox,npc
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8 xq(norder),wts(norder),umat,vmat,xyz(3)
      real *8, allocatable :: wts3(:),xref3(:,:),umat3(:,:)
      real *8 rintl(0:200)
      real *8 rint
      real *8, allocatable :: rintbs(:),rintbs2(:)
      integer i,itype,j,npols

      real *8, allocatable :: fval1(:,:,:,:),centerstmp(:,:,:)
      real *8, allocatable :: boxsize(:)
      integer, allocatable :: irefinebox(:)

      real *8 rsc,ra
      integer nbloc,nbctr,nbadd,irefine,ilev,ifirstbox,ilastbox
      integer nbtot,iii,idim
      

      nbmax = 100000
      nlmax = 200

      allocate(boxsize(0:nlmax))

      
      allocate(laddr(2,0:nlmax),ilevel(nbmax),iparent(nbmax))
      allocate(nchild(nbmax),ichild(8,nbmax))

      allocate(fvals(nd,norder**3,nbmax),centers(3,nbmax))

      allocate(rintbs(nbmax))

c
c      set tree info for level 0
c
      laddr(1,0) = 1
      laddr(2,0) = 1
      ilevel(1) = 0
      iparent(1) = -1
      nchild(1) = 0
      do i=1,8
        ichild(i,1) = -1
      enddo

      centers(1,1) = 0
      centers(2,1) = 0
      centers(3,1) = 0

c
c
      npbox = norder**3
      allocate(grid(3,npbox))
c
c     Generate a grid on the box [-1/2,1/2]^3
c
c
      itype = 1
      call legeexps(itype,norder,xq,umat,vmat,wts)
      do i=1,norder
        xq(i) = xq(i)/2
      enddo

      call mesh3d(xq,norder,xq,norder,xq,norder,grid)

      npols = norder*(norder+1)*(norder+2)/6
      
      allocate(wts3(npbox),xref3(3,npbox),umat3(npols,npbox))
      itype = 3
      call legetens_exps_3d(itype,norder,'t',xref3,umat3,npols,
     1   vmat,1,wts3)


c
c       compute fvals at the grid
c

      boxsize(0) =boxlen

      rint = 0
      rintbs(1) = 0

c
c   note extra factor of 8 sincee wts3 are on [-1,1]^3 
c   as opposed to [-1/2,1/2]^3
c
      rsc = boxlen**3/8

      do i=1,npbox
        xyz(1) = grid(1,i)*boxlen
        xyz(2) = grid(2,i)*boxlen
        xyz(3) = grid(3,i)*boxlen
        call fun(nd,xyz,dpars,zpars,ipars,fvals(1,i,1))
        if(iptype.eq.0) then
          do idim=1,nd
            if(abs(fvals(idim,i,1)).gt.rintbs(1)) rintbs(1) = 
     1          fvals(idim,i,1)
          enddo
        endif

        if(iptype.eq.1) then
          do idim=1,nd
            rintbs(1) = rintbs(1) + abs(fvals(idim,i,1))*wts3(i)*rsc
          enddo
        endif

        if(iptype.eq.2) then
          do idim=1,nd
            rintbs(1) = rintbs(1) + fvals(idim,i,1)**2*wts3(i)*rsc
          enddo
        endif
      enddo

      if(iptype.eq.0.or.iptype.eq.1) rint = rintbs(1)
      if(iptype.eq.2) rint = sqrt(rintbs(1))

      rintl(0) = rint

c
c       compute the interpolation matrix
c
      npc = 8*npbox
      allocate(ximat(npc,npbox),qwts(npc))
      call get_children_interp_mat(norder,npbox,npc,ximat)
      call get_children_qwts(norder,npc,wts,qwts)

      nbctr = 1
     

      do ilev=0,nlmax-1
        irefine = 0

        ifirstbox = laddr(1,ilev) 
        ilastbox = laddr(2,ilev)

        nbloc = ilastbox-ifirstbox+1

        allocate(fval1(nd,norder**3,8,nbloc))
        allocate(centerstmp(3,8,nbloc))
        allocate(irefinebox(nbloc))

        
        if(iptype.eq.2) rsc = sqrt(1.0d0/boxsize(0)**3)
        if(iptype.eq.1) rsc = (1.0d0/boxsize(0)**3)
        if(iptype.eq.0) rsc = 1.0d0

        rsc = rsc*rint

        print *, ilev, rint,rsc

        call vol_tree_find_box_refine(nd,iptype,eta,eps,zk,norder,npbox,
     1       fvals,npols,umat3,boxsize(ilev),nbmax,ifirstbox,
     2       nbloc,rsc,irefinebox,irefine)
     

c
c
c          figure out if current set of boxes is sufficient
c

        nbadd = 0 
        do i=1,nbloc
          if(irefinebox(i).eq.1) nbadd = nbadd+8
        enddo

        nbtot = nbctr+nbadd

c
c         if current memory is not sufficient reallocate
c
        if(nbtot.gt.nbmax) then
          print *, "Reallocating"
          allocate(centers2(3,nbmax),ilevel2(nbmax),iparent2(nbmax))
          allocate(nchild2(nbmax),ichild2(8,nbmax))
          allocate(fvals2(nd,npbox,nbmax),rintbs2(nbmax))

          call tree_copy(nd,nbctr,npbox,centers,ilevel,iparent,nchild,
     1            ichild,fvals,centers2,ilevel2,iparent2,nchild2,
     2            ichild2,fvals2)
          call dcopy(nbctr,rintbs,1,rintbs2,1)

          deallocate(centers,ilevel,iparent,nchild,ichild,fvals,rintbs)

          nbmax = nbtot
          allocate(centers(3,nbmax),ilevel(nbmax),iparent(nbmax))
          allocate(nchild(nbmax),ichild(8,nbmax),fvals(nd,npbox,nbmax))
          allocate(rintbs(nbmax))

          call tree_copy(nd,nbctr,npbox,centers2,ilevel2,iparent2,
     1            nchild2,ichild2,fvals2,centers,ilevel,iparent,nchild,
     2            ichild,fvals)
          call dcopy(nbctr,rintbs2,1,rintbs,1)

          deallocate(centers2,ilevel2,iparent2,nchild2,ichild2,fvals2)
          deallocate(rintbs2)
        endif


        if(irefine.eq.1) then
          boxsize(ilev+1) = boxsize(ilev)/2
          laddr(1,ilev+1) = nbctr+1

          call vol_tree_refine_boxes(irefinebox,nd,npbox,fvals,
     1      fun,dpars,zpars,ipars,grid,nbmax,ifirstbox,nbloc,centers,
     2      boxsize(ilev+1),nbctr,ilev+1,ilevel,iparent,nchild,ichild)


          rsc = boxsize(ilev+1)**3/8
          call update_rints(nd,npbox,nbmax,fvals,ifirstbox,nbloc,
     1       iptype,nchild,ichild,wts3,rsc,rintbs,rint)
          
          rintl(ilev+1) = rint

          
          laddr(2,ilev+1) = nbctr
        else
          exit
        endif

        deallocate(fval1,irefinebox,centerstmp)
      enddo

      nboxes = nbctr
      nlevels = ilev

      if(nlevels.ge.2) then

        nbtot = 16*nboxes
        if(nbtot.gt.nbmax) then
          print *, "Reallocating"
          allocate(centers2(3,nbmax),ilevel2(nbmax),iparent2(nbmax))
          allocate(nchild2(nbmax),ichild2(8,nbmax))
          allocate(fvals2(nd,npbox,nbmax),rintbs2(nbmax))

          call tree_copy(nd,nboxes,npbox,centers,ilevel,iparent,nchild,
     1            ichild,fvals,centers2,ilevel2,iparent2,nchild2,
     2            ichild2,fvals2)
          call dcopy(nboxes,rintbs,1,rintbs2,1)

          deallocate(centers,ilevel,iparent,nchild,ichild,fvals,rintbs)

          nbmax = nbtot
          allocate(centers(3,nbmax),ilevel(nbmax),iparent(nbmax))
          allocate(nchild(nbmax),ichild(8,nbmax),fvals(nd,npbox,nbmax))
          allocate(rintbs(nbmax))

          call tree_copy(nd,nboxes,npbox,centers2,ilevel2,iparent2,
     1          nchild2,ichild2,fvals2,centers,ilevel,iparent,nchild,
     2          ichild,fvals)
          call dcopy(nboxes,rintbs2,1,rintbs,1)

          deallocate(centers2,ilevel2,iparent2,nchild2,ichild2,fvals2)
          deallocate(rintbs2)
        endif

        allocate(nnbors(nbmax))
        allocate(nbors(27,nbmax))

        do i=1,nboxes
          nnbors(i) = 0
          do j=1,27
            nbors(j,i) = -1
          enddo
        enddo

        call computecoll(nlevels,nboxes,laddr,boxsize,centers,
     1        iparent,nchild,ichild,nnbors,nbors)

        if(nlevels.ge.2) then
          call vol_tree_fix_lr(fun,nd,dpars,zpars,ipars,norder,npbox,
     1         fvals,grid,centers,nlevels,nboxes,boxsize,nbmax,nlmax,
     2         laddr,ilevel,iparent,nchild,ichild,nnbors,nbors)
        endif
      endif

      ltree = 39*nboxes + 2*(nlevels+1)

      return
      end
c
c
c
c
c

      subroutine vol_tree_build(eps,zk,boxlen,norder,iptype,eta,
     1  fun,nd,dpars,zpars,ipars,nlevels,nboxes,ltree,rintl,itree,iptr,
     2  fvals,centers,boxsize)
c
c      compute the tree
c
c      input parameters:
c        eps - double precision
c           precision requested
c        zk - double complex
c           Helmholtz parameter
c        boxlen - double precision
c           length of box in which volume is contained, 
c           if boxlen = L, then volume is [-L/2,L/2]^3
c        norder - integer
c           order of discretization
c        iptype - integer
c           error norm
c           iptype = 0 - linf
c           iptype = 1 - l1
c           iptype = 2 - l2
c        eta - double precision
c           scaling parameter for error
c        fun - function handle
c           function to evalute it everywhere in the volume
c        nd - integer
c           number of real functions returned by fun
c        dpars - double precision
c           real parameters for function evaluation
c        zpars - double complex
c           complex parameters for function evaluation
c        ipars - integer
c           integer parameters for function evaluation
c        nlevels - integer
c          number of levels
c        nboxes - integer
c          number of boxes
c        ltree - integer
c          length of tree = 2*(nlevels+1)+39*nboxes
c        rintl - real *8 (0:nlevels)
c          estimate of lp norm for scaling the errors
c          at various levels. 
c          We require the estimate at each level to make sure
c          that the memory estimate code is consitent
c          with the build code else there could be potential
c          memory issues 
c         
c
c      output:
c        itree - integer(ltree)
c          tree info
c        iptr - integer(8)
c          iptr(1) - laddr
c          iptr(2) - ilevel
c          iptr(3) - iparent
c          iptr(4) - nchild
c          iptr(5) - ichild
c          iptr(6) - ncoll
c          iptr(7) - coll
c          iptr(8) - ltree
c        fvals - double precision (nd,norder**3,nboxes)
c          function values at discretization nodes
c        centers - double precision (3,nboxes)
c          xyz coordinates of box centers in the oct tree
c        boxsize - double precision (0:nlevels)
c          size of box at each of the levels
c

      implicit none
      real *8 eps,boxlen,eta,dpars(*)
      complex *16 zk,zpars(*)
      integer nd,ipars(*),iptype
      integer nlevels,nboxes,ltree,norder
      integer iptr(8)
      integer itree(ltree),ier
      real *8 fvals(nd,norder**3,nboxes),centers(3,nboxes)
      real *8, allocatable :: fval1(:,:,:,:),centerstmp(:,:,:)
      integer, allocatable :: irefinebox(:)
      real *8 boxsize(0:nlevels)
      real *8 xq(norder),wts(norder),umat(norder,norder)
      real *8 vmat(norder,norder)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8 rintl(0:nlevels)
      real *8 xyz(3)

      integer i,ilev,irefine,itype,nbmax,nlmax,npbox,npc,ii
      integer ifirstbox,ilastbox,nbctr,nbloc
      real *8 rsc
      real *8, allocatable :: wts3(:),xref3(:,:),umat3(:,:)

      real *8 ra
      integer j,nboxes0,npols

      external fun
c
      iptr(1) = 1
      iptr(2) = 2*(nlevels+1)+1
      iptr(3) = iptr(2) + nboxes
      iptr(4) = iptr(3) + nboxes
      iptr(5) = iptr(4) + nboxes
      iptr(6) = iptr(5) + 8*nboxes
      iptr(7) = iptr(6) + nboxes
      iptr(8) = iptr(7) + 27*nboxes



      boxsize(0) = boxlen

      centers(1,1) = 0
      centers(2,1) = 0
      centers(3,1) = 0

c
c      set tree info for level 0
c
      itree(1) = 1
      itree(2) = 1
      itree(iptr(2)) = 0
      itree(iptr(3)) = -1
      itree(iptr(4)) = 0
      do i=1,8
        itree(iptr(5)+i-1) = -1
      enddo
c
c
      npbox = norder**3
      allocate(grid(3,npbox))
c
c     Generate a grid on the box [-1/2,1/2]^3
c
c
      itype = 1
      call legeexps(itype,norder,xq,umat,vmat,wts)
      do i=1,norder
        xq(i) = xq(i)/2
      enddo

      call mesh3d(xq,norder,xq,norder,xq,norder,grid)
      npols = norder*(norder+1)*(norder+2)/6

      allocate(wts3(npbox),xref3(3,npbox),umat3(npols,npbox))
      itype = 3
      call legetens_exps_3d(itype,norder,'t',xref3,umat3,npols,vmat,
     1  1,wts3)




      do i=1,npbox
        xyz(1) = grid(1,i)*boxlen
        xyz(2) = grid(2,i)*boxlen
        xyz(3) = grid(3,i)*boxlen
        call fun(nd,xyz,dpars,zpars,ipars,fvals(1,i,1))
      enddo


c
c       compute the interpolation matrix
c
      npc = 8*npbox
      allocate(ximat(npc,npbox),qwts(npc))
      call get_children_interp_mat(norder,npbox,npc,ximat)
      call get_children_qwts(norder,npc,wts,qwts)


c
c       Reset nlevels, nboxes
c
      nbctr = 1

      do ilev=0,nlevels-1
        irefine = 0

        ifirstbox = itree(2*ilev+1) 
        ilastbox = itree(2*ilev+2)

        nbloc = ilastbox-ifirstbox+1

        allocate(fval1(nd,norder**3,8,nbloc))
        allocate(centerstmp(3,8,nbloc))
        allocate(irefinebox(nbloc))

        
        if(iptype.eq.2) rsc = sqrt(1.0d0/boxsize(0)**3)
        if(iptype.eq.1) rsc = (1.0d0/boxsize(0)**3)
        if(iptype.eq.0) rsc = 1.0d0
        rsc = rsc*rintl(ilev)
        call vol_tree_find_box_refine(nd,iptype,eta,eps,zk,norder,npbox,
     1       fvals,npols,umat3,boxsize(ilev),nboxes,ifirstbox,
     2       nbloc,rsc,irefinebox,irefine)
        

        if(irefine.eq.1) then
          boxsize(ilev+1) = boxsize(ilev)/2
          itree(2*ilev+3) = nbctr+1

          call vol_tree_refine_boxes(irefinebox,nd,npbox,fvals,
     1      fun,dpars,zpars,ipars,grid,nboxes,ifirstbox,nbloc,centers,
     2      boxsize(ilev+1),nbctr,ilev+1,itree(iptr(2)),itree(iptr(3)),
     3      itree(iptr(4)),itree(iptr(5)))
          
          itree(2*ilev+4) = nbctr
        else
          exit
        endif

        deallocate(fval1,irefinebox,centerstmp)
      enddo

      nboxes0 = nbctr
      nlevels = ilev


      do i=1,nboxes0
        itree(iptr(6)+i-1) = 0
        do j=1,27
          itree(iptr(7)+27*(i-1)+j-1) = -1
        enddo
      enddo

      call computecoll(nlevels,nboxes0,itree(iptr(1)),boxsize,centers,
     1        itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),
     2        itree(iptr(6)),itree(iptr(7)))

      if(nlevels.ge.2) then
         call vol_tree_fix_lr(fun,nd,dpars,zpars,ipars,norder,npbox,
     1       fvals,grid,centers,nlevels,nboxes0,boxsize,nboxes,nlevels,
     2       itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),
     3       itree(iptr(4)),itree(iptr(5)),itree(iptr(6)),
     4       itree(iptr(7)))
      endif

      call prinf('nboxes0=*',nboxes0,1)
      call prinf('nlevels=*',nlevels,1)
      return
      end
c
c
c      
c
c
      subroutine vol_tree_find_box_refine(nd,iptype,eta,eps,zk,norder,
     1  npbox,fvals,npols,umat,
     2  boxsize,nboxes,ifirstbox,nbloc,rsc,irefinebox,irefine)
      implicit none
      integer nd,npc,npbox,norder,iptype,npols
      integer nboxes,nbloc
      real *8 eta,eps,fvals(nd,npbox,nboxes)
      real *8 umat(npols,npbox)
      real *8 rsc
      real *8, allocatable :: fcoefs(:,:,:),rmask(:)
      integer, allocatable :: iind2p(:,:)
      real *8 alpha,beta,boxsize,rsum
      integer irefinebox(nbloc),xind(8),yind(8),zind(8)
      complex *16 zk
      integer ifirstbox

      integer irefine

      integer i,j,k,l,ibox,ifunif,i1
      real *8 rscale2,err,bs,bs2
      data xind/-1,1,-1,1,-1,1,-1,1/ 
      data yind/-1,-1,1,1,-1,-1,1,1/ 
      data zind/-1,-1,-1,-1,1,1,1,1/ 

      character *1 transa,transb

      external fun

      ifunif = 0

      transa = 'n'
      transb = 't'


      allocate(fcoefs(nd,npols,nbloc))
      allocate(iind2p(3,npols))
      call legetens_ind2pow_3d(norder-1,'t',iind2p)

      
      allocate(rmask(npols))

      rsum = 0
      do i=1,npols
        rmask(i) = 0.0d0
        i1=iind2p(1,i)+iind2p(2,i)+iind2p(3,i)
        if(i1.eq.norder-1) then
          rmask(i) = 1.0d0
          rsum = rsum + 1
        endif
      enddo

      if(iptype.eq.2) rsum = sqrt(rsum)
      if(iptype.eq.0) rsum = 1
      



      alpha = 1
      beta = 0

      bs = boxsize/4.0d0
      bs2 = 2*bs
      rscale2 = bs2**eta

      if(real(zk)*boxsize.gt.5) then
        do i=1,nbloc
          irefinebox(i) = 1
        enddo
        goto 1000
      endif


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,err)
      do i=1,nbloc
        irefinebox(i) = 0

        ibox = ifirstbox + i-1

        call dgemm(transa,transb,nd,npols,npbox,alpha,fvals(1,1,ibox),
     1    nd,umat,npols,beta,fcoefs(1,1,i),nd)

        call fun_err(nd,npols,fcoefs(1,1,i),rmask,
     1     iptype,rscale2,err)
     
        err = err/rsum

        
        if(err.gt.eps*rsc) then
          irefinebox(i) = 1
        endif
      enddo
C$OMP END PARALLEL DO     

 1000 continue
      irefine = maxval(irefinebox(1:nbloc))


c
c       make tree uniform
c

      if(ifunif.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nbloc
          irefinebox(i) = irefine
        enddo
C$OMP END PARALLEL DO      
      endif

      return
      end
c
c
c
c
c


      subroutine vol_tree_refine_boxes(irefinebox,nd,npbox,fvals,
     1  fun,dpars,zpars,ipars,grid,nboxes,ifirstbox,nbloc,centers,
     2  bs,nbctr,nlctr,ilevel,iparent,nchild,ichild)
      implicit none
      integer nd,npbox
      integer ipars(*)
      real *8 dpars(*)
      complex *16 zpars(*)
      real *8 fvals(nd,npbox,nboxes)
      integer nboxes,nbloc,nbctr,nlctr
      real *8 grid(3,npbox),bs,xyz(3)
      real *8 centers(3,nboxes)
      integer ilevel(nboxes),iparent(nboxes)
      integer ichild(8,nboxes),nchild(nboxes)
      integer irefinebox(nbloc)
      integer ifirstbox
      integer, allocatable :: isum(:)

      integer i,ibox,nel0,j,l,jbox,nel1,nbl
      integer xind(8),yind(8),zind(8)

      real *8 bsh
      data xind/-1,1,-1,1,-1,1,-1,1/ 
      data yind/-1,-1,1,1,-1,-1,1,1/ 
      data zind/-1,-1,-1,-1,1,1,1,1/


      external fun

      allocate(isum(nbloc))
      call cumsum(nbloc,irefinebox,isum)

      bsh = bs/2

      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,nbl,j,jbox,l,xyz)
      do i = 1,nbloc
        ibox = ifirstbox + i-1
        if(irefinebox(i).eq.1) then
          nbl = nbctr + (isum(i)-1)*8
          nchild(ibox) = 8
          do j=1,8
            jbox = nbl+j
            centers(1,jbox) = centers(1,ibox)+xind(j)*bsh
            centers(2,jbox) = centers(2,ibox)+yind(j)*bsh
            centers(3,jbox) = centers(3,ibox)+zind(j)*bsh

            do l=1,npbox
              xyz(1) = centers(1,jbox) + grid(1,l)*bs
              xyz(2) = centers(2,jbox) + grid(2,l)*bs
              xyz(3) = centers(3,jbox) + grid(3,l)*bs
              call fun(nd,xyz,dpars,zpars,ipars,fvals(1,l,jbox))
            enddo
            iparent(jbox) = ibox
            nchild(jbox) = 0
            do l=1,8
              ichild(l,jbox) = -1
            enddo
            ichild(j,ibox) = jbox
            ilevel(jbox) = nlctr 
          enddo
        endif
      enddo
C$OMP END PARALLEL DO      

      nbctr = nbctr + isum(nbloc)*8


      return
      end
c
c
c
c
c
      subroutine fun_err(nd,n,fcoefs,rmask,iptype,rscale,err)
c       this subroutine estimates the error based on the expansion
c       coefficients in a given basis
c       
c       input
c        nd - integer
c          number of functions
c        n -  integer
c           number of points at which function is tabulated
c        fcoefs: double precision(nd,n) 
c          tensor product legendre coeffs of func
c        rmask: double precision(n)
c           coefficients which will be accounted for in the error
c        iptype: integer
c           type of error to be computed
c           iptype = 0 - linf error
c           iptype = 1 - l1 error
c           iptype = 2 - l2 error
c        rscale: double precision
c          scaling factor
c
c       output
c         err: double precision
c           max scaled error in the functions
c
        implicit none
        integer n,i,iptype,idim,nd
        real *8 rscale,err
        real *8 fcoefs(nd,n),rmask(n)
        real *8, allocatable :: errtmp(:),ftmp(:,:)
        real *8 alpha, beta

        allocate(errtmp(nd),ftmp(nd,n))

        alpha = 1.0d0
        beta = 0.0d0
   
        err = 0
        do idim=1,nd
          errtmp(idim) = 0
        enddo
        if(iptype.eq.0) then
          do i=1,n
            if(rmask(i).gt.0.5d0) then
              do idim = 1,nd
                if(errtmp(idim).lt.abs(fcoefs(idim,i))) 
     1             errtmp(idim)=abs(fcoefs(idim,i))
              enddo
            endif
          enddo
        endif

        if(iptype.eq.1) then
          do i=1,n 
            do idim=1,nd
              ftmp(idim,i) = abs(fcoefs(idim,i))
            enddo
          enddo
          call dgemv('n',nd,n,alpha,ftmp,nd,rmask,1,beta,errtmp,1)
        endif


        if(iptype.eq.2) then
          do i=1,n 
            do idim=1,nd
              ftmp(idim,i) = fcoefs(idim,i)**2
            enddo
          enddo
          call dgemv('n',nd,n,alpha,ftmp,nd,rmask,1,beta,errtmp,1)

          do idim=1,nd
            errtmp(idim) = sqrt(errtmp(idim))
          enddo
        endif

        err = 0
        do idim=1,nd
          if(errtmp(idim).gt.err) err = errtmp(idim)
        enddo

        err = err*rscale

        return
        end

c
c
c
c       
c
      subroutine update_rints(nd,npbox,nbmax,fvals,ifirstbox,nbloc,
     1       iptype,nchild,ichild,wts,rsc,rintbs,rint)
c
c------------------------
c  This subroutine updates the integrals of the function to
c  be resolved on the computational domain. It subtracts the
c  integral of the boxes which have been refined and adds
c  in the integrals corresponding to the function values 
c  tabulated at the children
c
c  Input arguments:
c  
c    - nd: integer
c        number of functions
c    - npbox: integer
c        number of points per box where the function is tabulated
c    - nbmax: integer
c        max number of boxes
c    - fvals: real *8 (nd,npbox,nbmax)
c        tabulated function values
c    - ifirstbox: integer
c        first box in the list of boxes to be processed
c    - nbloc: integer
c        number of boxes to be processed
c    - iptype: integer
c        Lp version of the scheme
c        * iptype = 0, linf
c        * iptype = 1, l1
c        * iptype = 2, l2
c    - nchild: integer(nbmax)
c        number of children 
c    - ichild: integer(8,nbmax)
c        list of children
c    - wts: real *8 (npbox)
c        quadrature weights for intgegrating functions 
c    - rsc: real *8
c        scaling parameter for computing integrals
c  
c  Inout arguemnts:
c
c     - rintbs: real *8(nbmax)
c         the integral for the new boxes cretated will be updated
c     - rint: real *8
c         the total integral will be updated
c    
c  
c      
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: nd,npbox,nbmax
      real *8, intent(in) :: fvals(nd,npbox,nbmax)
      integer, intent(in) :: ifirstbox,nbloc,iptype
      integer, intent(in) :: nchild(nbmax),ichild(8,nbmax)
      real *8, intent(in) :: wts(npbox),rsc
      real *8, intent(inout) :: rintbs(nbmax),rint

c
c
c      compute the integrals for the newly formed boxes
c   and update the overall integral
c
      if(iptype.eq.0) then
        do i=1,nbloc
          ibox = ifirstbox+i-1
          if(nchild(ibox).gt.0) then
            do j=1,8
              jbox = ichild(j,ibox)
              rintbs(jbox) = maxval(fvals(1:nd,1:npbox,jbox))
              if(rintbs(jbox).gt.rint) rint = rintbs(jbox)
            enddo
          endif
        enddo
      endif

      if(iptype.eq.1) then
        do i=1,nbloc
          ibox=ifirstbox+i-1
          if(nchild(ibox).gt.0) then
c     subtract contribution of ibox from rint
            rint = rint - rintbs(ibox) 
          endif
        enddo

c
c     add back contribution of children
c
        do i=1,nbloc
          ibox = ifirstbox+i-1
          if(nchild(ibox).gt.0) then
            do j=1,8
              jbox = ichild(j,ibox)
              rintbs(jbox) = 0
              do l=1,npbox
                do idim=1,nd
                  rintbs(jbox) = rintbs(jbox) + 
     1               abs(fvals(idim,l,jbox))*wts(l)*rsc
                enddo
              enddo
              rint = rint + rintbs(jbox)
            enddo
          endif
        enddo
      endif

      if(iptype.eq.2) then
        rintsq = rint**2
        do i=1,nbloc
          ibox=ifirstbox+i-1
          if(nchild(ibox).gt.0) then
c
c    note that if iptype = 2, then rintbs stores squares
c    of the integral on the box
c
             rintsq = rintsq - rintbs(ibox)
          endif
        enddo

        do i=1,nbloc
          ibox = ifirstbox+i-1
          if(nchild(ibox).gt.0) then
            do j=1,8
              jbox = ichild(j,ibox)
              rintbs(jbox) = 0
              do l=1,npbox
                do idim=1,nd
                  rintbs(jbox) = rintbs(jbox) + 
     1               fvals(idim,l,jbox)**2*wts(l)*rsc
                enddo
              enddo
              rintsq = rintsq + rintbs(jbox)
            enddo
          endif
        enddo
        rint = sqrt(rintsq)
      endif
          

      return
      end
c
c
c
c
c


      subroutine get_children_interp_mat(norder,npbox,npc,ximat)
c
c
c       construct interpolation matrix from nodes on parent box to 
c       nodes on 8 children boxes
c
c        input
c        norder: integer
c           order of discretization
c        npbox: integer
c           number of points per box = norder**3
c        npc: integer
c           number of points in chilren boxe = 8*npbox
c
c        output
c          ximat - double precision(npc,npbox)
c             interpolation matrix
c        
c
      implicit real *8 (a-h,o-z)
      integer norder,npbox,npc
      real *8 ximat(npc,npbox),cref(3,8),xyz(3)
      real *8, allocatable :: umat(:,:),vmat(:,:),pmat(:,:)
      real *8, allocatable :: xref(:,:),wts(:)
      character *1 type,transa,transb

      npols = norder*(norder+1)*(norder+2)/6

      allocate(xref(3,npbox),umat(npols,npbox),vmat(npbox,npols))
      allocate(pmat(npols,npc),wts(npbox))

      type = 'T'
      itype = 2
      call legetens_exps_3d(itype,norder,type,xref,umat,npols,vmat,
     1   npbox,wts)
      
      
      do ic=1,8
        ii = 2
        jj = 2
        if(ic.eq.1.or.ic.eq.2.or.ic.eq.5.or.ic.eq.6) ii=1
        if(ic.lt.5) jj = 1
        cref(1,ic) = (-1)**ic*0.5d0
        cref(2,ic) = (-1)**ii*0.5d0
        cref(3,ic) = (-1)**jj*0.5d0
        do j=1,npbox
          ipt = (ic-1)*npbox+j
          do l=1,3
            xyz(l) = cref(l,ic) + xref(l,j)*0.5d0
          enddo
          call legetens_pols_3d(xyz,norder-1,type,pmat(1,ipt))
        enddo
      enddo

      transa = 't'
      transb = 'n'

      alpha = 1
      beta = 0

      call dgemm(transa,transb,npc,npbox,npols,alpha,pmat,npols,
     1   umat,npols,beta,ximat,npc)
      

      return
      end
c
c
c
c
c
      subroutine get_children_qwts(norder,npc,wts,qwts)
      implicit real *8 (a-h,o-z)
      real *8 wts(norder),qwts(npc)
      
      ipt = 1
      do i=1,norder
        do j=1,norder
          do k=1,norder
            qwts(ipt) = wts(i)*wts(j)*wts(k)/8
            ipt = ipt+1
          enddo
        enddo
      enddo

      npbox = norder**3
      do i=1,7
        do j=1,npbox
          qwts(i*npbox+j) = qwts(j)
        enddo
      enddo


      return
      end
c
c
c
c
c
c
       subroutine tree_copy(nd,nb,npb,centers,ilevel,iparent,
     1            nchild,ichild,fvals,centers2,ilevel2,iparent2,nchild2,
     2            ichild2,fvals2)

       implicit none
       integer nd,nb,npb
       real *8 centers(3,nb),centers2(3,nb)
       integer ilevel(nb),ilevel2(nb)
       integer iparent(nb),iparent2(nb)
       integer nchild(nb),nchild2(nb)
       integer ichild(8,nb),ichild2(8,nb)
       real *8 fvals(nd,npb,nb),fvals2(nd,npb,nb)

       integer i,j,nel

       nel = nd*npb*nb
       call dcopy(nel,fvals,1,fvals2,1)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
       do i=1,nb
         centers2(1,i) = centers(1,i)
         centers2(2,i) = centers(2,i)
         centers2(3,i) = centers(3,i)
         ilevel2(i) = ilevel(i)
         iparent2(i) = iparent(i)
         nchild2(i) = nchild(i)
         do j=1,8
           ichild2(j,i) = ichild(j,i)
         enddo
       enddo
C$OMP END PARALLEL DO       
       

       return
       end
c
c
c
c
c

      subroutine computecoll(nlevels,nboxes,laddr,boxsize,
     1                       centers,iparent,nchild,ichild,
     2                       nnbors,nbors)

c     This subroutine computes the colleagues for an adaptive
c     pruned tree. box j is a colleague of box i, if they share a
c     vertex or an edge and the two boxes are at the same
c     level in the tree
c
c     INPUT arguments
c     nlevels     in: integer
c                 Number of levels
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     laddr       in: integer(2,0:nlevels)
c                 indexing array providing access to boxes at
c                 each level. 
c                 the first box on level i is laddr(1,i)
c                 the last box on level i is laddr(2,i)
c
c     boxsize     in: double precision(0:;nlevels)
c                 Array of boxsizes
c 
c     centers     in: double precision(3,nboxes)
c                 array of centers of boxes
c   
c     iparent     in: integer(nboxes)
c                 iparent(i) is the box number of the parent of
c                 box i
c
c     nchild      in: integer(nboxes)
c                 nchild(i) is the number of children of box i
c
c     ichild      in: integer(8,nboxes)
c                 ichild(j,i) is the box id of the jth child of
c                 box i
c
c----------------------------------------------------------------
c     OUTPUT
c     nnbors      out: integer(nboxes)
c                 nnbors(i) is the number of colleague boxes of
c                 box i
c
c     nbors       out: integer(27,nboxes)
c                 nbors(j,i) is the box id of the jth colleague
c                 box of box i
c---------------------------------------------------------------
      implicit none
      integer nlevels,nboxes
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(3,nboxes)
      integer iparent(nboxes), nchild(nboxes), ichild(8,nboxes)
      integer nnbors(nboxes)
      integer nbors(27,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,dad
      integer i,j,ifirstbox,ilastbox


c     Setting parameters for level = 0
      nnbors(1) = 1
      nbors(1,1) = 1
      do ilev = 1,nlevels
c        Find the first and the last box at level ilev      
         ifirstbox = laddr(1,ilev)
         ilastbox = laddr(2,ilev)
c        Loop over all boxes to evaluate neighbors, list1 and updating
c        hunglists of targets

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox)
         do ibox = ifirstbox,ilastbox
c           Find the parent of the current box         
            dad = iparent(ibox)
c           Loop over the neighbors of the parent box
c           to find out list 1 and list 2
            do i=1,nnbors(dad)
                jbox = nbors(i,dad)
                do j=1,8
c               ichild(j,jbox) is one of the children of the
c               neighbors of the parent of the current
c               box
                   kbox = ichild(j,jbox)
                   if(kbox.gt.0) then
c               Check if kbox is a nearest neighbor or in list 2
                      if((abs(centers(1,kbox)-centers(1,ibox)).le.
     1                   1.05*boxsize(ilev)).and.
     2                   (abs(centers(2,kbox)-centers(2,ibox)).le.
     3                   1.05*boxsize(ilev)).and.
     4                   (abs(centers(3,kbox)-centers(3,ibox)).le.
     5                   1.05*boxsize(ilev))) then
                     
                         nnbors(ibox) = nnbors(ibox)+1
                         nbors(nnbors(ibox),ibox) = kbox
                      endif
                   endif
                enddo
            enddo
c           End of computing colleagues of box i
         enddo
C$OMP END PARALLEL DO         
      enddo

      return
      end
c-------------------------------------------------------------      
      subroutine vol_tree_fix_lr(fun,nd,dpars,zpars,ipars,norder,
     1       npbox,fvals,grid,centers,nlevels,nboxes,boxsize,
     2       nbmax,nlmax,laddr,ilevel,iparent,nchild,ichild,nnbors,
     3       nbors)
c
c
c       convert an adaptive tree into a level restricted tree
c
      implicit none
      integer nd,ipars(*),norder,npbox,nlevels,nboxes,nlmax
      integer nbmax
      real *8 dpars(*),fvals(nd,npbox,nbmax),grid(3,npbox)
      real *8 centers(3,nbmax),boxsize(0:nlmax)
      complex *16 zpars(*)
      integer laddr(2,0:nlmax),ilevel(nbmax),iparent(nbmax)
      integer nchild(nbmax),ichild(8,nbmax),nnbors(nbmax)
      integer nbors(27,nbmax)
      integer laddrtail(2,0:nlmax),isum
      integer, allocatable :: iflag(:)

      integer i,j,k,l,ibox,jbox,kbox,ilev,idad,igranddad
      integer nbloc,ict
      real *8 xdis,ydis,zdis,distest

      external fun

      allocate(iflag(nbmax))

c     Initialize flag array
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nboxes
         iflag(i) = 0
      enddo
C$OMP END PARALLEL DO     



c     Flag boxes that violate level restriction by "1"
c     Violatioin refers to any box that is directly touching
c     a box that is more than one level finer
c
c     Method:
c     1) Carry out upward pass. For each box B, look at
c     the colleagues of B's grandparent
c     2) See if any of those colleagues are childless and in
c     contact with B.
c
c     Note that we only need to get up to level two, as
c     we will not find a violation at level 0 and level 1
c
c     For such boxes, we set iflag(i) = 1
c
      do ilev=nlevels,2,-1
c        This is the distance to test if two boxes separated
c        by two levels are touching
         distest = 1.05d0*(boxsize(ilev-1) + boxsize(ilev-2))/2.0d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,idad,igranddad,i,jbox)         
C$OMP$ PRIVATE(ict,xdis,ydis,zdis)
         do ibox = laddr(1,ilev),laddr(2,ilev) 
            idad = iparent(ibox)
            igranddad = iparent(idad)
            
c           Loop over colleagues of granddad            
            do i=1,nnbors(igranddad)
               jbox = nbors(i,igranddad)
c              Check if the colleague of grandad
c              is a leaf node. This automatically
c              eliminates the granddad
               if(nchild(jbox).eq.0.and.iflag(jbox).eq.0) then
                   xdis = centers(1,jbox) - centers(1,idad)
                   ydis = centers(2,jbox) - centers(2,idad)
                   zdis = centers(3,jbox) - centers(3,idad)
                   ict = 0
                   if(abs(xdis).le.distest) ict = ict + 1
                   if(abs(ydis).le.distest) ict = ict + 1
                   if(abs(zdis).le.distest) ict = ict + 1
                   if(ict.eq.3) then
                      iflag(jbox) = 1
                   endif
               endif
c              End of checking criteria for the colleague of
c              granddad
            enddo
c           End of looping over colleagues of
c           granddad
         enddo
c        End of looping over boxes at ilev         
C$OMP END PARALLEL DO
      enddo
c     End of looping over levels and flagging boxes


c     Find all boxes that need to be given a flag+
c     A flag+ box will be denoted by setting iflag(box) = 2
c     This refers to any box that is not already flagged and
c     is bigger than and is contacting a flagged box
c     or another box that has already been given a flag +.
c     It is found by performing an upward pass and looking
c     at the flagged box's parents colleagues and a flag+
c     box's parents colleagues and seeing if they are
c     childless and present the case where a bigger box 
c     is contacting a flagged or flag+ box.

      do ilev = nlevels,1,-1
c        This is the distance to test if two boxes separated
c        by one level are touching
         distest = 1.05d0*(boxsize(ilev) + boxsize(ilev-1))/2.0d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,idad,i,jbox,xdis,ydis)
C$OMP$PRIVATE(zdis,ict)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            if(iflag(ibox).eq.1.or.iflag(ibox).eq.2) then
               idad = iparent(ibox)
c              Loop over dad's colleagues               
               do i=1,nnbors(idad)
                  jbox = nbors(i,idad)
c                 Check if the colleague of dad
c                 is a leaf node. This automatically
c                 eliminates the dad
                  if(nchild(jbox).eq.0.and.iflag(jbox).eq.0) then
                     xdis = centers(1,jbox) - centers(1,ibox)
                     ydis = centers(2,jbox) - centers(2,ibox)
                     zdis = centers(3,jbox) - centers(3,ibox)
                     ict = 0
                     if(abs(xdis).le.distest) ict = ict + 1
                     if(abs(ydis).le.distest) ict = ict + 1
                     if(abs(zdis).le.distest) ict = ict + 1
                     if(ict.eq.3) then
                        iflag(jbox) = 2
                     endif
                  endif
c                 End of checking criteria for the colleague of
c                dad
               enddo
c              End of looping over dad's colleagues               
            endif
c           End of checking if current box is relevant for
c           flagging flag+ boxes
         enddo
c        End of looping over boxes at ilev        
C$OMP END PARALLEL DO 
      enddo
c     End of looping over levels

c     Subdivide all flag and flag+ boxes. Flag all the children
c     of flagged boxes as flag++. Flag++ boxes are denoted
c     by setting iflag(box) = 3. The flag++ boxes need 
c     to be checked later to see which of them need further
c     refinement. While creating new boxes, we will
c     need to update all the tree structures as well.
c     Note that all the flagged boxes live between
c     levels 1 and nlevels - 2. We process the boxes via a
c     downward pass. We first determine the number of boxes
c     that are going to be subdivided at each level and 
c     everything else accordingly
      do ilev = 0,nlevels
         laddrtail(1,ilev) = 0
         laddrtail(2,ilev) = -1
      enddo
 
      do ilev = 1,nlevels-2
c        First subdivide all the flag and flag+
c        boxes with boxno nboxes+1, nboxes+ 2
c        and so on. In the second step, we reorganize
c        all the structures again to bring it back
c        in the standard format

         laddrtail(1,ilev+1) = nboxes+1

         nbloc = laddr(2,ilev)-laddr(1,ilev)+1
         call vol_tree_refine_boxes_flag(iflag,nd,npbox,fvals,
     1    fun,dpars,zpars,ipars,grid,nbmax,laddr(1,ilev),nbloc,
     2    centers,boxsize(ilev+1),nboxes,ilev,ilevel,iparent,nchild,
     3    ichild)

         laddrtail(2,ilev+1) = nboxes
      enddo
c     Reorganize the tree to get it back in the standard format

      call vol_tree_reorg(nboxes,nd,npbox,centers,nlevels,laddr,
     1          laddrtail,ilevel,iparent,nchild,ichild,fvals,iflag)

c     Compute colleague information again      

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
      do i=1,nboxes
         nnbors(i) = 0
         do j=1,27
            nbors(j,i) = -1
         enddo
      enddo
C$OMP END PARALLEL DO      
      call computecoll(nlevels,nboxes,laddr, boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,nnbors,nbors)

c     Processing of flag and flag+ boxes is done
c     Start processing flag++ boxes. We will use a similar
c     strategy as before. We keep checking the flag++
c     boxes that require subdivision if they still
c     violate the level restriction criterion, create
c     the new boxes, append them to the end of the list to begin
c     with and in the end reorganize the tree structure.
c     We shall accomplish this via a downward pass
c     as new boxes that get added in the downward pass
c     will also be processed simultaneously.
c     We shall additionally also need to keep on updating
c     the colleague information as we proceed in the 
c     downward pass

c     Reset the flags array to remove all the flag and flag+
c     cases. This is to ensure reusability of the subdivide
c     _flag routine to handle the flag++ case

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox)
      do ibox=1,nboxes
         if(iflag(ibox).ne.3) iflag(ibox) = 0
      enddo
C$OMP END PARALLEL DO      
 
      do ilev = 0,nlevels
         laddrtail(1,ilev) = 0
         laddrtail(2,ilev) = -1
      enddo

      do ilev = 2,nlevels-2

c     Step 1: Determine which of the flag++ boxes need
c     further division. In the even a flag++ box needs
c     further subdivision then flag the box with iflag(box) = 1
c     This will again ensure that the subdivide_flag routine
c     will take care of handling the flag++ case
         call updateflags(ilev,nboxes,nlevels,laddr,nchild,ichild,
     1                    nnbors,nbors,centers,boxsize,iflag)

         call updateflags(ilev,nboxes,nlevels,laddrtail,nchild,ichild,
     1                    nnbors,nbors,centers,boxsize,iflag)
         
c      Step 2: Subdivide all the boxes that need subdivision
c      in the laddr set and the laddrtail set as well
         laddrtail(1,ilev+1) = nboxes + 1

         nbloc = laddr(2,ilev)-laddr(1,ilev)+1
         call vol_tree_refine_boxes_flag(iflag,nd,npbox,fvals,
     1    fun,dpars,zpars,ipars,grid,nbmax,laddr(1,ilev),nbloc,
     2    centers,boxsize(ilev+1),nboxes,ilev,ilevel,iparent,nchild,
     3    ichild)

         nbloc = laddrtail(2,ilev)-laddrtail(1,ilev)+1
         call vol_tree_refine_boxes_flag(iflag,nd,npbox,fvals,
     1    fun,dpars,zpars,ipars,grid,nbmax,laddrtail(1,ilev),nbloc,
     2    centers,boxsize(ilev+1),nboxes,ilev,ilevel,iparent,nchild,
     3    ichild)

          laddrtail(2,ilev+1) = nboxes         
c      Step 3: Update the colleague information for the newly
c      created boxes

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,i,idad,jbox,j,kbox)
          do ibox = laddrtail(1,ilev+1),laddrtail(2,ilev+1)
            nnbors(ibox) = 0
c           Find the parent of the current box         
            idad = iparent(ibox)
c           Loop over the neighbors of the parent box
c           to find out colleagues
            do i=1,nnbors(idad)
                jbox = nbors(i,idad)
                do j=1,8
c               ichild(j,jbox) is one of the children of the
c               neighbors of the parent of the current
c               box
                   kbox = ichild(j,jbox)
                   if(kbox.gt.0) then
c               Check if kbox is a nearest neighbor or in list 2
                      if((abs(centers(1,kbox)-centers(1,ibox)).le.
     1                   1.05*boxsize(ilev+1)).and.
     2                   (abs(centers(2,kbox)-centers(2,ibox)).le.
     3                   1.05*boxsize(ilev+1)).and.
     4                   abs(centers(3,kbox)-centers(3,ibox)).le.
     5                   1.05*boxsize(ilev+1)) then
                     
                         nnbors(ibox) = nnbors(ibox)+1
                         nbors(nnbors(ibox),ibox) = kbox
                      endif
                   endif
                enddo
            enddo
c           End of computing colleagues of box i
         enddo
C$OMP END PARALLEL DO         
      enddo

c     Reorganize tree once again and we are all done      
      call vol_tree_reorg(nboxes,nd,npbox,centers,nlevels,laddr,
     1          laddrtail,ilevel,iparent,nchild,ichild,fvals,iflag)

c     Compute colleague information again      

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
      do i=1,nboxes
         nnbors(i) = 0
         do j=1,27
            nbors(j,i) = -1
         enddo
      enddo
C$OMP END PARALLEL DO    

      call computecoll(nlevels,nboxes,laddr, boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,nnbors,nbors)
      

      return
      end
      

c-------------------------------------------------------------      
      subroutine vol_tree_reorg(nboxes,nd,npbox,centers,nlevels,laddr,
     1     laddrtail,ilevel,iparent,nchild,ichild,fvals,iflag)

c    This subroutine reorganizes the current data in all the tree
c    arrays to rearrange them in the standard format.
c    The boxes on input are assumed to be arranged in the following
c    format
c    boxes on level i are the boxes from laddr(1,i) to 
c    laddr(2,i) and also from laddrtail(1,i) to laddrtail(2,i)
c
c    At the end of the sorting, the boxes on level i
c    are arranged from laddr(1,i) to laddr(2,i)  
c
c    INPUT/OUTPUT arguments
c    nboxes         in: integer
c                   number of boxes
c
c    nd             in: integer
c                   number of real value functions
c
c    npbox          in: integer
c                   number of grid points per function
c
c    centers        in/out: double precision(3,nboxes)
c                   x and y coordinates of the center of boxes
c
c    nlevels        in: integer
c                   Number of levels in the tree
c
c    laddr          in/out: integer(2,0:nlevels)
c                   boxes at level i are numbered between
c                   laddr(1,i) to laddr(2,i)
c
c    laddrtail      in: integer(2,0:nlevels)
c                   new boxes to be added to the tree
c                   structure are numbered from
c                   laddrtail(1,i) to laddrtail(2,i)
c
c    ilevel      in/out: integer(nboxes)
c                ilevel(i) is the level of box i
c
c    iparent     in/out: integer(nboxes)
c                 iparent(i) is the parent of box i
c
c    nchild      in/out: integer(nboxes)
c                nchild(i) is the number of children 
c                of box i
c
c    ichild       in/out: integer(8,nboxes)
c                 ichild(j,i) is the jth child of box i
c
c    iflag        in/out: integer(nboxes)
c                 iflag(i) is a flag for box i required to generate
c                 level restricted tree from adaptive tree

      implicit none
c     Calling sequence variables and temporary variables
      integer nboxes,nlevels,npbox,nd
      double precision centers(3,nboxes)
      integer laddr(2,0:nlevels), tladdr(2,0:nlevels)
      integer laddrtail(2,0:nlevels)
      integer ilevel(nboxes)
      integer iparent(nboxes)
      integer nchild(nboxes)
      integer ichild(8,nboxes)
      integer iflag(nboxes)
      double precision fvals(nd,npbox,nboxes)
      
      integer, allocatable :: tilevel(:),tiparent(:),tnchild(:)
      integer, allocatable :: tichild(:,:),tiflag(:)
      integer, allocatable :: iboxtocurbox(:),ilevptr(:),ilevptr2(:)

      double precision, allocatable :: tfvals(:,:,:),tcenters(:,:)



c     Temporary variables
      integer i,j,k,l
      integer ibox,ilev, curbox,idim,nblev

      allocate(tilevel(nboxes),tiparent(nboxes),tnchild(nboxes))
      allocate(tichild(8,nboxes),tiflag(nboxes),iboxtocurbox(nboxes))
      allocate(tfvals(nd,npbox,nboxes),tcenters(3,nboxes))

      do ilev = 0,nlevels
         tladdr(1,ilev) = laddr(1,ilev)
         tladdr(2,ilev) = laddr(2,ilev)
      enddo
      call tree_copy(nd,nboxes,npbox,centers,ilevel,iparent,nchild,
     1            ichild,fvals,tcenters,tilevel,tiparent,tnchild,
     2            tichild,tfvals)

      do ibox=1,nboxes
         tiflag(ibox) = iflag(ibox)
      enddo
     
c     Rearrange old arrays now

      do ilev = 0,1
         do ibox = laddr(1,ilev),laddr(2,ilev)
           iboxtocurbox(ibox) = ibox
         enddo
      enddo

      allocate(ilevptr(nlevels+1),ilevptr2(nlevels))

      ilevptr(2) = laddr(1,2)


      do ilev=2,nlevels
        nblev = laddr(2,ilev)-laddr(1,ilev)+1
        ilevptr2(ilev) = ilevptr(ilev) + nblev
        nblev = laddrtail(2,ilev)-laddrtail(1,ilev)+1
        ilevptr(ilev+1) = ilevptr2(ilev) + nblev
      enddo

      curbox = laddr(1,2)
      do ilev=2,nlevels
         laddr(1,ilev) = curbox
         do ibox = tladdr(1,ilev),tladdr(2,ilev)
            ilevel(curbox) = tilevel(ibox)
            nchild(curbox) = tnchild(ibox)
            centers(1,curbox) = tcenters(1,ibox)
            centers(2,curbox) = tcenters(2,ibox)
            centers(3,curbox) = tcenters(3,ibox)
            do i=1,npbox
              do idim=1,nd
                fvals(idim,i,ibox) = tfvals(idim,i,ibox)
              enddo
            enddo
            iflag(curbox) = tiflag(ibox)
            iboxtocurbox(ibox) = curbox

            curbox = curbox + 1
         enddo
         do ibox = laddrtail(1,ilev),laddrtail(2,ilev)
            ilevel(curbox) = tilevel(ibox)
            centers(1,curbox) = tcenters(1,ibox)
            centers(2,curbox) = tcenters(2,ibox)
            centers(3,curbox) = tcenters(3,ibox)
            nchild(curbox) = tnchild(ibox)
            do i=1,npbox
              do idim=1,nd
                fvals(idim,i,ibox) = tfvals(idim,i,ibox)
              enddo
            enddo
            iflag(curbox) = tiflag(ibox)
            iboxtocurbox(ibox) = curbox

            curbox = curbox + 1
         enddo
         laddr(2,ilev) = curbox-1
      enddo

c     Handle the parent children part of the tree 
c     using the mapping iboxtocurbox

      do ibox=1,nboxes
         if(tiparent(ibox).eq.-1) iparent(iboxtocurbox(ibox)) = -1
         if(tiparent(ibox).gt.0) 
     1    iparent(iboxtocurbox(ibox)) = iboxtocurbox(tiparent(ibox))
         do i=1,8
            if(tichild(i,ibox).eq.-1) ichild(i,iboxtocurbox(ibox)) = -1
            if(tichild(i,ibox).gt.0) 
     1      ichild(i,iboxtocurbox(ibox)) = iboxtocurbox(tichild(i,ibox))
         enddo
      enddo

      return
      end
c--------------------------------------------------------------------      
      subroutine updateflags(curlev,nboxes,nlevels,laddr,nchild,ichild,
     1                    nnbors,nbors,centers,boxsize,iflag)

c      This subroutine is to check the boxes flagged as flag++
c      and determine which of the boxes need refinement. The flag
c      of the box which need refinement is updated to iflag(box)=1
c      and that of the boxes which do not need refinement is
c      updated to iflag(box) = 0
c
c      INPUT arguments
c      curlev         in: integer
c                     the level for which boxes need to be processed
c
c      nboxes         in: integer
c                     total number of boxes
c
c      nlevels        in: integer
c                     total number of levels
c
c      laddr          in: integer(2,0:nlevels)
c                     boxes from laddr(1,ilev) to laddr(2,ilev)
c                     are at level ilev
c
c      nchild         in: integer(nboxes)
c                     nchild(ibox) is the number of children
c                     of box ibox
c
c      ichild         in: integer(4,nboxes)
c                     ichild(j,ibox) is the box id of the jth
c                     child of box ibox
c
c      nnbors         in: integer(nboxes)
c                     nnbors(ibox) is the number of colleagues
c                     of box ibox
c
c      nbors          in: integer(27,nboxes)
c                     nbors(j,ibox) is the jth colleague of box
c                     ibox
c
c      centers        in: double precision(3,nboxes)
c                     x and y coordinates of the box centers
c
c      boxsize        in: double precision(0:nlevels)
c                     boxsize(i) is the size of the box at level i
c
c      iflag          in/out: integer(nboxes)
c                     iflag(ibox)=3 if it is flag++. iflag(ibox) =1
c                     or 0 at the end of routine depending on
c                     whether box needs to be subdivided or not
c
      implicit none
c     Calling sequence variables
      integer curlev, nboxes, nlevels
      integer laddr(2,0:nlevels),nchild(nboxes),ichild(8,nboxes)
      integer nnbors(nboxes), nbors(27,nboxes)
      integer iflag(nboxes)
      double precision centers(3,nboxes),boxsize(0:nlevels)

c     Temporary variables
      integer i,j,k,l,ibox,jbox,kbox,lbox, ict
      double precision distest,xdis,ydis,zdis

      distest = 1.05d0*(boxsize(curlev) + boxsize(curlev+1))/2.0d0
c     Loop over all boxes at the current level     

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,i,jbox,j,kbox,xdis,ydis)
C$OMP$PRIVATE(zdis,ict)
      do ibox = laddr(1,curlev),laddr(2,curlev)
         if(iflag(ibox).eq.3) then
            iflag(ibox) = 0
c           Loop over colleagues of the current box      
            do i=1,nnbors(ibox)
c              Loop over colleagues of flag++ box        
               jbox = nbors(i,ibox)
              
c              Loop over the children of the colleague box
c              Note we do not need to exclude self from
c              the list of colleagues as a self box which
c              is flag++ does not have any children 
c              and will not enter the next loop
               do j=1,8
                  kbox = ichild(j,jbox)
                  if(kbox.gt.0) then
                     if(nchild(kbox).gt.0) then
                        xdis = centers(1,kbox) - centers(1,ibox)
                        ydis = centers(2,kbox) - centers(2,ibox)
                        zdis = centers(3,kbox) - centers(3,ibox)
                        ict = 0
                        if(abs(xdis).le.distest) ict = ict + 1
                        if(abs(ydis).le.distest) ict = ict + 1
                        if(abs(zdis).le.distest) ict = ict + 1
                        if(ict.eq.3) then
                           iflag(ibox) = 1
                           goto 1111
                        endif
                     endif
                  endif
c                 End of looping over the children of the child
c                 of the colleague box
               enddo
c              End of looping over the children of the colleague box       
            enddo
c           End of looping over colleagues            
 1111       continue        
         endif
c        End of testing if the current box needs to checked for         
      enddo
c     End of looping over boxes at the current level      
C$OMP END PARALLEL DO      

      return
      end
c
c
c
c
c


      subroutine vol_tree_refine_boxes_flag(iflag,nd,npbox,fvals,
     1  fun,dpars,zpars,ipars,grid,nboxes,ifirstbox,nbloc,centers,
     2  bs,nbctr,nlctr,ilevel,iparent,nchild,ichild)
      implicit none
      integer nd,npbox,nboxes
      real *8 fvals(nd,npbox,nboxes)
      integer nbloc,nbctr,nlctr
      real *8 centers(3,nboxes),bs,grid(3,npbox),xyz(3)
      integer ilevel(nboxes),iparent(nboxes)
      integer ichild(8,nboxes),nchild(nboxes)
      integer iflag(nboxes)
      integer ifirstbox,ilastbox
      integer, allocatable :: isum(:)
      real *8 dpars(*)
      integer ipars(*)
      complex *16 zpars(*)

      integer i,ibox,nel,j,l,jbox,nbl,ii
      integer xind(8),yind(8),zind(8)

      real *8 bsh
      data xind/-1,1,-1,1,-1,1,-1,1/ 
      data yind/-1,-1,1,1,-1,-1,1,1/ 
      data zind/-1,-1,-1,-1,1,1,1,1/

      external fun


      ilastbox = ifirstbox+nbloc-1


      bsh = bs/2

      allocate(isum(nbloc))

      call cumsum_nz(nbloc,iflag(ifirstbox),isum)

      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,j,jbox,nbl,l)
C$OMP$PRIVATE(xyz)
      do ibox = ifirstbox,ilastbox
        if(iflag(ibox).gt.0) then
          nchild(ibox) = 8
          nbl = nbctr + (isum(ibox-ifirstbox+1)-1)*8
          do j=1,8
            jbox = nbl+j
            centers(1,jbox) = centers(1,ibox)+xind(j)*bsh
            centers(2,jbox) = centers(2,ibox)+yind(j)*bsh
            centers(3,jbox) = centers(3,ibox)+zind(j)*bsh
            do l=1,npbox
              xyz(1) = centers(1,jbox) + grid(1,l)*bs
              xyz(2) = centers(2,jbox) + grid(2,l)*bs
              xyz(3) = centers(3,jbox) + grid(3,l)*bs
              call fun(nd,xyz,dpars,zpars,ipars,fvals(1,l,jbox))
            enddo
            iparent(jbox) = ibox
            nchild(jbox) = 0
            do l=1,8
              ichild(l,jbox) = -1
            enddo
            ichild(j,ibox) = jbox
            ilevel(jbox) = nlctr+1 
            if(iflag(ibox).eq.1) iflag(jbox) = 3
            if(iflag(ibox).eq.2) iflag(jbox) = 0
          enddo
        endif
      enddo
C$OMP END PARALLEL DO
      
      if(nbloc.gt.0) nbctr = nbctr + isum(nbloc)*8


      return
      end
c
c
c
c
c
      subroutine cumsum(n,a,b)
c
c        this subroutine computes the cumulative sum of an array
c
c
c       TODO: need to make this routine openmp
c
c       input:
c         n - number of elements
c         a - integer(n)
c              input array
c
c       output:
c         b - integer(n)
c            b(i) = sum_{j=1}^{i} a(j)
c
      implicit none
      integer a(n),b(n),n,i,isum
      isum = 0


      do i=1,n
        isum = isum + a(i)
        b(i) = isum
      enddo
      
      return
      end
c
c
c
c
c
      subroutine cumsum_nz(n,a,b)
c
c        this subroutine computes the cumulative sum of postive
c        entries of an array
c
c
c       TODO: need to make this routine openmp
c
c       input:
c         n - number of elements
c         a - integer(n)
c              input array
c
c       output:
c         b - integer(n)
c            b(i) = sum_{j=1}^{i} I_{a(j)>0}
c
      implicit none
      integer a(n),b(n),n,i,isum

      isum = 0
      do i=1,n
        if(a(i).gt.0) isum = isum+1
        b(i) = isum
      enddo
      
      return
      end
c
c
c
c
c
c
      subroutine get_list1(nboxes,nlevels,itree,ltree,iptr,
     1     centers,boxsize,nlist1,list1)
      implicit real *8 (a-h,o-z)
c
c       this subroutine returns the list1 of a given a tree
c
c       input:
c         nboxes - integer
c            total number of boxes
c         nlevels - integer
c            number of levels
c         itree - integer(ltree)
c            the tree structure
c         ltree - integer
c            length of tree
c         iptr - integer(8)
c           iptr(1) - laddr
c           iptr(2) - ilevel
c           iptr(3) - iparent
c           iptr(4) - nchild
c           iptr(5) - ichild
c           iptr(6) - ncoll
c           iptr(7) - coll
c           iptr(8) - ltree
c         centers - real *8 (3,nboxes)
c           centers of boxes in the tree structure
c         boxsize - real *8 (0:nlevels)
c           size of the boxes at each level
c         
c       output:
c         nlist1 - integer(nboxes)
c           number of boxes in list1
c         list1 - integer(139,nboxes)
c           list of boxes in list1
c
      integer nboxes,nlevels,itree(ltree),iptr(8),ltree
      integer list1(139,nboxes),nlist1(nboxes)
      real *8 centers(3,nboxes),boxsize(0:nlevels)
      integer firstbox,lastbox,dad

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,j)
      do ibox=1,nboxes
        nlist1(ibox) = 0
        do j=1,139
          list1(j,ibox) = -1
        enddo
      enddo
C$OMP END PARALLEL DO      
    
      if(itree(iptr(4)).eq.0) then
        nlist1(1) = 1
        list1(14,1) = 1
        return
      endif


      do ilev=1,nlevels
        firstbox = itree(iptr(1)+2*ilev)
        lastbox = itree(iptr(1)+2*ilev+1)
        distest0 = 1.05d0*(boxsize(ilev)+boxsize(ilev-1))/2.0d0
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nchild,nnbors,i,jbox,nchildj,ix,iy,iz,iind)
C$OMP$PRIVATE(j,kbox,xdis,ydis,zdis,dad)
        do ibox=firstbox,lastbox
          nc = itree(iptr(4)+ibox-1)
          if(nc.eq.0) then
            nnbors = itree(iptr(6)+ibox-1)
            do i=1,nnbors
              jbox = itree(iptr(7)+27*(ibox-1)+i-1)
              nchildj = itree(iptr(4)+jbox-1)
              if(nchildj.eq.0) then
                nlist1(ibox) = nlist1(ibox)+1
                
                ix = (centers(1,jbox)-centers(1,ibox))/boxsize(ilev)
                iy = (centers(2,jbox)-centers(2,ibox))/boxsize(ilev)
                iz = (centers(3,jbox)-centers(3,ibox))/boxsize(ilev)

                iind = (iz+1)*3*3 + (iy+1)*3 + ix+2
                list1(iind,ibox) = jbox
              else
                do j=1,8
                  distest = 1.05d0*(boxsize(ilev)+boxsize(ilev+1))/2.0d0
                  kbox = itree(iptr(5)+(jbox-1)*8+j-1)
                  if(kbox.gt.0) then
                    xdis = centers(1,kbox)-centers(1,ibox)
                    ydis = centers(2,kbox)-centers(2,ibox)
                    zdis = centers(3,kbox)-centers(3,ibox)
                    if(abs(xdis).lt.distest.and.abs(ydis).lt.distest.
     1                  and.abs(zdis).lt.distest) then
                      nlist1(ibox) = nlist1(ibox)+1
                      ix = (xdis + boxsize(ilev)/4)/boxsize(ilev)*2.0d0
                      iy = (ydis + boxsize(ilev)/4)/boxsize(ilev)*2.0d0
                      iz = (zdis + boxsize(ilev)/4)/boxsize(ilev)*2.0d0

                      ix = ix + 1
                      iy = iy + 1
                      iz = iz + 1

                      call get_iind_list1(ix,iy,iz,iind)
                      iind = iind + 27

                      list1(iind,ibox) = kbox
                    endif
                  endif
                enddo
              endif
            enddo
c
c
c                compute list1 at level ilev-1
c
            dad = itree(iptr(3)+ibox-1)
            nnbors = itree(iptr(6)+dad-1)
            do i=1,nnbors
              jbox = itree(iptr(7)+27*(dad-1)+i-1)
              nchildj = itree(iptr(4)+jbox-1)
              if(nchildj.eq.0) then
                xdis = centers(1,jbox)-centers(1,ibox)
                ydis = centers(2,jbox)-centers(2,ibox)
                zdis = centers(3,jbox)-centers(3,ibox)

                if(abs(xdis).lt.distest0.and.abs(ydis).lt.distest0.and.
     1               abs(zdis).lt.distest0) then
                  ix = (xdis - boxsize(ilev)/2.0d0)/boxsize(ilev)
                  iy = (ydis - boxsize(ilev)/2.0d0)/boxsize(ilev)
                  iz = (zdis - boxsize(ilev)/2.0d0)/boxsize(ilev)
                  ix = ix + 2
                  iy = iy + 2
                  iz = iz + 2
                  call get_iind_list1(ix,iy,iz,iind)
                  iind = iind + 27 + 56
                  nlist1(ibox) = nlist1(ibox)+1
                  list1(iind,ibox) = jbox
                endif
              endif
            enddo
          endif
        enddo
C$OMP END PARALLEL DO         
      enddo

      return
      end
c
c
c
c
      subroutine get_iind_list1(ix,iy,iz,iind)
      implicit none
      integer ix,iy,iz,iind

      iind = iz*16 + iy*4 + ix + 1
      if(iz.ge.2) iind = iind-(iz-1)*4
      if((iz.eq.1.or.iz.eq.2).and.iy.ge.2) iind = iind-(iy-1)*2
      if((iz.eq.2.or.iz.eq.1).and.(iy.eq.2.or.iy.eq.1).and.ix.eq.3) 
     1    iind = iind - 2
       
      return
      end

c
c
c
c
c 
      subroutine get_list1boxes_type(itype,istart,iend,nboxes,
     1   nlist1,list1,ijlist,n)
      implicit real *8 (a-h,o-z)
      integer itype,istart,iend,nlist1(nboxes),list1(139,nboxes),
     1   ijlist(2,nboxes)

      n = 0
      do i=istart,iend
        if(list1(itype,i).gt.0) then
          n = n+1
          ijlist(1,n) = i
          ijlist(2,n) = list1(itype,i)
        endif
      enddo

      return
      end

      
      
c-----------------------------------------------------------------
      subroutine computemnlists(nlevels,nboxes,laddr,boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,isep,nnbors,mnbors,nbors,mnlist1,
     3                   mnlist2,mnlist3,mnlist4)
c     Compute max nuber of boxes in list1,list2,list3,list4
      implicit none
      integer nlevels,nboxes
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(3,nboxes)
      integer iparent(nboxes),nchild(nboxes),ichild(8,nboxes)
      integer mnbors,isep
      integer nnbors(nboxes),nbors(mnbors,nboxes)
      integer mnlist1,mnlist2,mnlist3,mnlist4
      integer nlist1(nboxes),nlist2(nboxes),nlist3(nboxes)
      integer nlist4(nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,i,j,k,l
      integer firstbox,lastbox,dad
      double precision xdis,ydis,zdis,distest


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nboxes
         nlist1(i) = 0
         nlist2(i) = 0
         nlist3(i) = 0
         nlist4(i) = 0
      enddo
C$OMP END PARALLEL DO
      
      nlist1(1) = 1
      nlist2(1) = 0
      nlist3(1) = 0
      nlist4(1) = 0

      do ilev = 1,nlevels
         firstbox = laddr(1,ilev)
         lastbox = laddr(2,ilev)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox,distest,xdis,ydis,zdis)
         do ibox = firstbox,lastbox
            dad = iparent(ibox)
            do i=1,nnbors(dad)
               jbox = nbors(i,dad)
               do j=1,8
                  kbox = ichild(j,jbox)
                  if(kbox.gt.0) then
                     if((abs(centers(1,kbox)-centers(1,ibox)).ge.
     1                  1.05d0*isep*boxsize(ilev)).or.
     2                  (abs(centers(2,kbox)-centers(2,ibox)).ge.
     3                  1.05d0*isep*boxsize(ilev)).or.
     4                  (abs(centers(3,kbox)-centers(3,ibox)).ge.
     5                  1.05d0*isep*boxsize(ilev))) then

                        nlist2(ibox) = nlist2(ibox) + 1
                     endif
                  endif
               enddo
            enddo
c           Compute list1 and list3 of ibox if it is childless
            if(nchild(ibox).eq.0) then
               do i=1,nnbors(ibox)
                  jbox = nbors(i,ibox)

c
cc                     check for list1 at the same level
c
                  if(nchild(jbox).eq.0) then
                     nlist1(ibox) = nlist1(ibox) + 1
                  endif
c
cc                     check for list1 and list3 at one ilev+1
                  if(nchild(jbox).gt.0) then
                     distest = 1.05d0*(boxsize(ilev)+boxsize(ilev+1))/
     1                         2.0d0*isep
                     do j=1,8
                        kbox = ichild(j,jbox)
                        if(kbox.gt.0) then
                           xdis = dabs(centers(1,kbox)-centers(1,ibox))
                           ydis = dabs(centers(2,kbox)-centers(2,ibox))
                           zdis = dabs(centers(3,kbox)-centers(3,ibox))

                           if(xdis.lt.distest.and.ydis.lt.distest.and.
     1                        zdis.lt.distest) then
                              nlist1(ibox) = nlist1(ibox)+1
                           else
                              nlist3(ibox) = nlist3(ibox)+1
                           endif
                        endif
                     enddo
                  endif
               enddo
c
cc               compute list1 and list4 for boxes at level ilev-1 
               do i=1,nnbors(dad)
                   jbox = nbors(i,dad)
                   if(nchild(jbox).eq.0) then
                      distest = 1.05d0*(boxsize(ilev)+boxsize(ilev-1))/
     1                         2.0d0*isep
                      xdis = dabs(centers(1,jbox)-centers(1,ibox))
                      ydis = dabs(centers(2,jbox)-centers(2,ibox))
                      zdis = dabs(centers(3,jbox)-centers(3,ibox))
                      if(xdis.lt.distest.and.ydis.lt.distest.and.
     1                  zdis.lt.distest) then
                         nlist1(ibox) = nlist1(ibox)+1
                      endif
                   endif
               enddo
            endif
c
cc           compute list 4 at level ilev-1
c
            do i=1,nnbors(dad)
               jbox = nbors(i,dad)
               if(nchild(jbox).eq.0) then
                   distest = 1.05d0*(boxsize(ilev)+boxsize(ilev-1))/
     1                 2.0d0*isep
                    xdis = dabs(centers(1,jbox)-centers(1,ibox))
                    ydis = dabs(centers(2,jbox)-centers(2,ibox))
                    zdis = dabs(centers(3,jbox)-centers(3,ibox))
                    if(xdis.gt.distest.or.ydis.gt.distest.or.
     1                 zdis.gt.distest) then
                       nlist4(ibox) = nlist4(ibox)+1
                    endif
               endif
            enddo
         enddo
C$OMP END PARALLEL DO         
      enddo

      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) 
C$OMP$REDUCTION(max:mnlist1,mnlist2,mnlist3,mnlist4)      
      do i=1,nboxes
         if(nlist1(i).gt.mnlist1) mnlist1 = nlist1(i)
         if(nlist2(i).gt.mnlist2) mnlist2 = nlist2(i)
         if(nlist3(i).gt.mnlist3) mnlist3 = nlist3(i)
         if(nlist4(i).gt.mnlist4) mnlist4 = nlist4(i)
      enddo
C$OMP END PARALLEL DO      

      return
      end
c---------------------------------------------------------------      
      subroutine computelists(nlevels,nboxes,laddr,boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,isep,nnbors,mnbors,nbors,nlist1,
     3                   mnlist1,list1,nlist2,mnlist2,list2,
     4                   nlist3,mnlist3,list3,nlist4,mnlist4,list4)
c     Compute max nuber of boxes in list1,list2,list3,list4
      implicit none
      integer nlevels,nboxes
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(3,nboxes)
      integer iparent(nboxes),nchild(nboxes),ichild(8,nboxes)
      integer mnbors,isep
      integer nnbors(nboxes),nbors(mnbors,nboxes)
      integer mnlist1,mnlist2,mnlist3,mnlist4
      integer nlist1(nboxes),nlist2(nboxes),nlist3(nboxes)
      integer nlist4(nboxes)
      integer list1(mnlist1,nboxes),list2(mnlist2,nboxes)
      integer list3(mnlist3,nboxes),list4(mnlist4,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,i,j,k,l
      integer firstbox,lastbox,dad
      double precision xdis,ydis,zdis,distest

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,nboxes
         nlist1(i) = 0
         nlist2(i) = 0
         nlist3(i) = 0
         nlist4(i) = 0
      enddo
C$OMP END PARALLEL DO      
      if(nchild(1).eq.0) then
         nlist1(1) = 1
         list1(1,1) = 1
      else
         nlist1(1) = 0
      endif
      nlist2(1) = 0
      nlist3(1) = 0
      nlist4(1) = 0

      do ilev = 1,nlevels
         firstbox = laddr(1,ilev)
         lastbox = laddr(2,ilev)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox,xdis,ydis,zdis,distest)
         do ibox = firstbox,lastbox
            dad = iparent(ibox)
            do i=1,nnbors(dad)
               jbox = nbors(i,dad)
               do j=1,8
                  kbox = ichild(j,jbox)
                  if(kbox.gt.0) then
                     if((abs(centers(1,kbox)-centers(1,ibox)).ge.
     1                  1.05d0*isep*boxsize(ilev)).or.
     2                  (abs(centers(2,kbox)-centers(2,ibox)).ge.
     3                  1.05d0*isep*boxsize(ilev)).or.
     4                  (abs(centers(3,kbox)-centers(3,ibox)).ge.
     5                  1.05d0*isep*boxsize(ilev))) then

                        nlist2(ibox) = nlist2(ibox) + 1
                        list2(nlist2(ibox),ibox) = kbox
                     endif
                  endif
               enddo
            enddo
c           Compute list1 and list3 of ibox if it is childless
            if(nchild(ibox).eq.0) then
               do i=1,nnbors(ibox)
                  jbox = nbors(i,ibox)

c
cc                boxes in list 1 at the same level
c

                  if(nchild(jbox).eq.0) then
                     nlist1(ibox) = nlist1(ibox) + 1
                     list1(nlist1(ibox),ibox) = jbox
                  else
c
cc                     boxes in list1 and list3 at level ilev+1
c
                     distest = 1.05d0*(boxsize(ilev)+boxsize(ilev+1))/
     1                         2.0d0*isep
                     do j=1,8
                        kbox = ichild(j,jbox)
                        if(kbox.gt.0) then
                           xdis = dabs(centers(1,kbox)-centers(1,ibox))
                           ydis = dabs(centers(2,kbox)-centers(2,ibox))
                           zdis = dabs(centers(3,kbox)-centers(3,ibox))

                           if(xdis.lt.distest.and.ydis.lt.distest.and.
     1                        zdis.lt.distest) then
                              nlist1(ibox) = nlist1(ibox)+1
                              list1(nlist1(ibox),ibox) = kbox
                           else
                              nlist3(ibox) = nlist3(ibox)+1
                              list3(nlist3(ibox),ibox) = kbox
                           endif
                        endif
                     enddo
                  endif
               enddo
c
cc               compute list1 at level ilev-1 
               do i=1,nnbors(dad)
                   jbox = nbors(i,dad)
                   if(nchild(jbox).eq.0) then
                      distest = 1.05d0*(boxsize(ilev)+boxsize(ilev-1))/
     1                         2.0d0*isep
                      xdis = dabs(centers(1,jbox)-centers(1,ibox))
                      ydis = dabs(centers(2,jbox)-centers(2,ibox))
                      zdis = dabs(centers(3,jbox)-centers(3,ibox))
                      if(xdis.lt.distest.and.ydis.lt.distest.and.
     1                  zdis.lt.distest) then
                         nlist1(ibox) = nlist1(ibox)+1
                         list1(nlist1(ibox),ibox) = jbox
                      endif
                   endif
                enddo
            endif
c
cc           compute list 4 at level ilev-1
c
            do i=1,nnbors(dad)
               jbox = nbors(i,dad)
               if(nchild(jbox).eq.0) then
                   distest = 1.05d0*(boxsize(ilev)+boxsize(ilev-1))/
     1                 2.0d0*isep
                    xdis = dabs(centers(1,jbox)-centers(1,ibox))
                    ydis = dabs(centers(2,jbox)-centers(2,ibox))
                    zdis = dabs(centers(3,jbox)-centers(3,ibox))
                    if(xdis.gt.distest.or.ydis.gt.distest.or.
     1                 zdis.gt.distest) then
                       nlist4(ibox) = nlist4(ibox)+1
                       list4(nlist4(ibox),ibox)=jbox
                    endif
               endif
            enddo
         enddo
C$OMP END PARALLEL DO         
      enddo

      return
      end
c----------------------------------------------------------------

      subroutine getpwlistall(ibox,bs,nboxes,nnbors,nbors,
     1           nchild,ichild,centers,isep,nuall,uall,ndall,dall,nnall,
     2           nall,nsall,sall,neall,eall,nwall,wall,nu1234,u1234,
     3           nd5678,d5678,nn1256,n1256,ns3478,s3478,ne1357,e1357,
     4           nw2468,w2468,nn12,n12,nn56,n56,ns34,s34,ns78,s78,ne13,
     5           e13,ne57,e57,nw24,w24,nw68,w68,ne1,e1,ne3,e3,ne5,e5,
     6           ne7,e7,nw2,w2,nw4,w4,nw6,w6,nw8,w8)
c-------------------------------------------------------------------
      implicit none
      integer ibox
      double precision boxsize,bs
      integer nboxes,nnbors,nbors(nnbors)
      integer nchild, ichild(8,nboxes)
      double precision centers(3,nboxes)
      integer isep
      integer nuall,ndall,nnall,nsall,neall,nwall,nu1234
      integer nd5678,nn1256,ns3478,ne1357,nw2468
      integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
      integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8
      integer uall(1),dall(1),nall(1),sall(1),eall(1),wall(1)
      integer u1234(1),d5678(1),n1256(1),s3478(1),e1357(1),w2468(1)
      integer n12(1),n56(1),s34(1),s78(1),e13(1),e57(1),w24(1),w68(1)
      integer e1(1),e3(1),e5(1),e7(1)
      integer w2(1),w4(1),w6(1),w8(1)

      integer jbox,kbox
      integer c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c16
      integer i,j

      nuall = 0
      ndall = 0
      nnall = 0
      nsall = 0
      neall = 0
      nwall = 0
      nu1234 = 0
      nd5678 = 0
      nn1256 = 0
      ns3478 = 0
      ne1357 = 0
      nw2468 = 0
      nn12 = 0
      nn56 = 0
      ns34 = 0
      ns78 = 0
      ne13 = 0
      ne57 = 0
      nw24 = 0
      nw68 = 0
      ne1 = 0
      ne3 = 0
      ne5 = 0
      ne7 = 0
      nw2 = 0
      nw4 = 0
      nw6 = 0
      nw8 = 0
      do i=1,nnbors
         jbox = nbors(i)
         do j=1,8
            kbox = ichild(j,jbox)
            if(kbox.gt.0) then
               c1 = 0
               c2 = 0
               c3 = 0
               c4 = 0
               c5 = 0
               c6 = 0
               c7 = 0
               c8 = 0
               c9 = 0
               c10 = 0
               c11 = 0
               c12 = 0
               if((centers(3,kbox)-centers(3,ibox)).ge.
     1              1.01d0*isep*bs+bs/2.0d0) c1 = 1
               if((centers(3,kbox)-centers(3,ibox)).le.
     1              -1.01d0*isep*bs-bs/2.0d0) c2 = 1
               if((centers(2,kbox)-centers(2,ibox)).ge.
     1              1.01d0*isep*bs+bs/2.0d0) c3 = 1
               if((centers(2,kbox)-centers(2,ibox)).le.
     1              -1.01d0*isep*bs-bs/2.0d0) c4 = 1
               if((centers(1,kbox)-centers(1,ibox)).ge.
     1              1.01d0*isep*bs+bs/2.0d0) c5 = 1
               if((centers(1,kbox)-centers(1,ibox)).le.
     1              -1.01d0*isep*bs-bs/2.0d0) c6 = 1
               if((centers(3,kbox)-centers(3,ibox)).ge.
     1              1.01d0*isep*bs-bs/2.0d0) c7 = 1
               if((centers(3,kbox)-centers(3,ibox)).le.
     1              -1.01d0*isep*bs+bs/2.0d0) c8 = 1
               if((centers(2,kbox)-centers(2,ibox)).ge.
     1              1.01d0*isep*bs-bs/2.0d0) c9 = 1
               if((centers(2,kbox)-centers(2,ibox)).le.
     1              -1.01d0*isep*bs+bs/2.0d0) c10 = 1
               if((centers(1,kbox)-centers(1,ibox)).ge.
     1              1.01d0*isep*bs-bs/2.0d0) c11 = 1
               if((centers(1,kbox)-centers(1,ibox)).le.
     1              -1.01d0*isep*bs+bs/2.0d0) c12 = 1
               if(c1.eq.1) then
                  nuall = nuall + 1
                  uall(nuall) = kbox
               endif

               if(c2.eq.1) then
                  ndall = ndall + 1
                  dall(ndall) = kbox
               endif

               if(c3.eq.1.and.c1.eq.0.and.c2.eq.0) then
                  nnall = nnall + 1
                  nall(nnall) = kbox
               endif

               if(c4.eq.1.and.c1.eq.0.and.c2.eq.0) then   
                  nsall = nsall + 1
                  sall(nsall) = kbox
               endif

               if(c5.eq.1.and.c1.eq.0.and.c2.eq.0.and.c3.eq.0.and.
     1             c4.eq.0) then
                  neall = neall + 1
                  eall(neall) = kbox
               endif

               if(c6.eq.1.and.c1.eq.0.and.c2.eq.0.and.c3.eq.0.and.
     1            c4.eq.0) then
                  nwall = nwall + 1
                  wall(nwall) = kbox
               endif

               c16 = c1 + c2 + c3 + c4 + c5 +c6
               if(c16.eq.0.and.c7.eq.1) then
                  nu1234 = nu1234 + 1
                  u1234(nu1234) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1) then
                  nd5678 = nd5678 + 1
                  d5678(nd5678) = kbox
               endif

               if(c16.eq.0.and.c9.eq.1.and.c7.eq.0.and.c8.eq.0) then
                  nn1256 = nn1256 + 1
                  n1256(nn1256) = kbox
               endif

               if(c16.eq.0.and.c10.eq.1.and.c7.eq.0.and.c8.eq.0) then
                  ns3478 = ns3478 + 1
                  s3478(ns3478) = kbox
               endif

               if(c16.eq.0.and.c11.eq.1.and.(c7+c8+c9+c10).eq.0) then
                  ne1357 = ne1357 + 1
                  e1357(ne1357) = kbox
               endif

               if(c16.eq.0.and.c12.eq.1.and.(c7+c8+c9+c10).eq.0) then
                  nw2468 = nw2468 + 1
                  w2468(nw2468) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1.and.c9.eq.1) then
                  nn12 = nn12 + 1
                  n12(nn12) = kbox
               endif

               if(c16.eq.0.and.c7.eq.1.and.c9.eq.1) then
                  nn56 = nn56 + 1
                  n56(nn56) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1.and.c10.eq.1) then
                  ns34 = ns34 + 1
                  s34(ns34) = kbox
               endif

               if(c16.eq.0.and.c7.eq.1.and.c10.eq.1) then
                  ns78 = ns78 + 1
                  s78(ns78) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1.and.c11.eq.1
     1           .and.c9.eq.0.and.c10.eq.0) then
                  ne13 = ne13 + 1
                  e13(ne13) = kbox
               endif

               if(c16.eq.0.and.c7.eq.1.and.c11.eq.1
     1           .and.c9.eq.0.and.c10.eq.0) then
                  ne57 = ne57 + 1
                  e57(ne57) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1.and.c12.eq.1
     1           .and.c9.eq.0.and.c10.eq.0) then
                  nw24 = nw24 + 1
                  w24(nw24) = kbox
               endif

               if(c16.eq.0.and.c7.eq.1.and.c12.eq.1
     1           .and.c9.eq.0.and.c10.eq.0) then
                  nw68 = nw68 + 1
                  w68(nw68) = kbox
               endif

               if(c16.eq.0.and.c7.eq.0.and.c10.eq.1.and.c11.eq.1) then
                  ne1 = ne1 + 1
                  e1(ne1) = kbox
               endif
               if(c16.eq.0.and.c7.eq.0.and.c9.eq.1.and.c11.eq.1) then
                  ne3 = ne3 + 1
                  e3(ne3) = kbox
               endif
               if(c16.eq.0.and.c8.eq.0.and.c10.eq.1.and.c11.eq.1) then
                  ne5 = ne5 + 1
                  e5(ne5) = kbox
               endif
               if(c16.eq.0.and.c8.eq.0.and.c9.eq.1.and.c11.eq.1) then
                  ne7 = ne7 + 1
                  e7(ne7) = kbox
               endif
               if(c16.eq.0.and.c7.eq.0.and.c10.eq.1.and.c12.eq.1) then
                  nw2 = nw2 + 1
                  w2(nw2) = kbox
               endif
               if(c16.eq.0.and.c7.eq.0.and.c9.eq.1.and.c12.eq.1) then
                  nw4 = nw4 + 1
                  w4(nw4) = kbox
               endif
               if(c16.eq.0.and.c8.eq.0.and.c10.eq.1.and.c12.eq.1) then
                  nw6 = nw6 + 1
                  w6(nw6) = kbox
               endif
               if(c16.eq.0.and.c8.eq.0.and.c9.eq.1.and.c12.eq.1) then
                  nw8 = nw8 + 1
                  w8(nw8) = kbox
               endif
            endif
         enddo
      enddo
      

      return
      end
c--------------------------------------------------------------------     
c
c
c
c--------------------------------------------------------------------      

      subroutine getlist3pwlistall(ibox,bs,nboxes,nlist3,list3,
     1           isep,centers,nuall,uall,ndall,dall,nnall,
     2           nall,nsall,sall,neall,eall,nwall,wall)
c-------------------------------------------------------------------
      implicit none
      integer ibox
      integer isep
      integer nboxes,nlist3,list3(nlist3)
      double precision centers(3,nboxes)
      double precision sepdist,bs
      integer nuall,ndall,nnall,nsall,neall,nwall
      integer uall(1),dall(1),nall(1),sall(1),eall(1),wall(1)

      integer jbox
      integer c1,c2,c3,c4,c5,c6
      integer j

      nuall = 0
      ndall = 0
      nnall = 0
      nsall = 0
      neall = 0
      nwall = 0
      sepdist = 1.01d0*isep*bs+bs/2.0d0
      do j=1,nlist3
         jbox = list3(j)
C         print *,"jbox: ",jbox
         if(jbox.gt.0) then
            c1 = 0
            c2 = 0
            c3 = 0
            c4 = 0
            c5 = 0
            c6 = 0
            if((centers(3,jbox)-centers(3,ibox)).ge.sepdist) c1 = 1
            if((centers(3,jbox)-centers(3,ibox)).le.-sepdist) c2 = 1
            if((centers(2,jbox)-centers(2,ibox)).ge.sepdist) c3 = 1
            if((centers(2,jbox)-centers(2,ibox)).le.-sepdist) c4 = 1
            if((centers(1,jbox)-centers(1,ibox)).ge.sepdist) c5 = 1
            if((centers(1,jbox)-centers(1,ibox)).le.-sepdist) c6 = 1

C            print *,c1,c2,c3,c4,c5,c6

            if(c1.eq.1) then
               nuall = nuall + 1
               uall(nuall) = jbox
            endif

            if(c2.eq.1) then
               ndall = ndall + 1
               dall(ndall) = jbox
            endif

            if(c3.eq.1.and.c1.eq.0.and.c2.eq.0) then
               nnall = nnall + 1
               nall(nnall) = jbox
            endif

            if(c4.eq.1.and.c1.eq.0.and.c2.eq.0) then   
               nsall = nsall + 1
               sall(nsall) = jbox
            endif

            if(c5.eq.1.and.c1.eq.0.and.c2.eq.0.and.c3.eq.0.and.
     1         c4.eq.0) then
               neall = neall + 1
               eall(neall) = jbox
            endif

            if(c6.eq.1.and.c1.eq.0.and.c2.eq.0.and.c3.eq.0.and.
     1         c4.eq.0) then
               nwall = nwall + 1
               wall(nwall) = jbox
            endif
         endif
      enddo
      

      return
      end
c--------------------------------------------------------------------      
c
c
c
c
c--------------------------------------------------------------
c
      subroutine subdividebox(pos,npts,center,boxsize,
     1           isorted,iboxfl,subcenters)
      implicit none
      double precision pos(3,npts)
      double precision center(3)
      double precision subcenters(3,8)
      double precision boxsize
      integer npts
      integer isorted(*)
      integer iboxfl(2,8)

c     Temporary variables
      integer isortedtmp(npts)
      integer i,j,i12,i34,istart,jstart,kstart,ii,iii,nss,nee
      integer jj,irefinebox,ntt
      integer i56, i78, i1234, i5678
      integer ibox,ifirstbox,ilastbox,nbfirst
      integer is,it,ie
      integer nc(8)

      i1234 = 0
      i5678 = 0
      do i = 1,npts
         if(pos(3,i)-center(3).lt.0) then
            i1234 = i1234 + 1
            isorted(i1234) = i
         else
            i5678 = i5678 + 1
            isortedtmp(i5678) = i
         endif
      enddo
      do i=1,i5678
         isorted(i1234+i) = isortedtmp(i)
      enddo

c     Sort i1234 into i12 and i34         
      i12 = 0
      i34 = 0
      do i = 1,i1234
         if(pos(2,isorted(i))-center(2).lt.0) then
            i12 = i12 + 1
            isorted(i12) = isorted(i)
         else
            i34 = i34 + 1
            isortedtmp(i34) = isorted(i)
         endif
      enddo
c     Note at the end of the loop, i12 is where the particles
c     in part 12 of the box end
c
c     Reorder sources to include 34 in the array
      do i=1,i34
         isorted(i12+i) = isortedtmp(i)
      enddo

c     sort i5678 into i56 and i78
      i56 = i1234
      i78 = 0
      do i=i1234+1,npts
         if(pos(2,isorted(i))-center(2).lt.0) then
            i56 = i56 + 1
            isorted(i56) = isorted(i)
         else
            i78 = i78 + 1
            isortedtmp(i78) = isorted(i)
         endif
      enddo
c     Reorder sources to include 78 in the array
      do i=1,i78
         isorted(i56+i) = isortedtmp(i)
      enddo
c     End of reordering i5678         

      nc = 0
c     Sort into boxes 1 and 2
      do i = 1,i12
         if(pos(1,isorted(i))-center(1).lt.0) then
            nc(1) = nc(1) + 1
            isorted(nc(1)) = isorted(i)
         else
            nc(2) = nc(2) + 1
            isortedtmp(nc(2)) = isorted(i)
         endif
      enddo
c     Reorder sources so that sources in 2 are at the
c     end of this part of the array
      do i=1,nc(2)
         isorted(nc(1)+i) = isortedtmp(i)
      enddo

c     Sort into boxes 3 and 4
      do i = i12+1, i1234
         if(pos(1,isorted(i))-center(1).lt.0) then
            nc(3) = nc(3) + 1
            isorted(i12+nc(3)) = isorted(i)
         else
            nc(4) = nc(4)+1
            isortedtmp(nc(4)) = isorted(i)
         endif
      enddo
c     Reorder sources so that sources in 4 are at the
c     end of this part of the array
      do i=1,nc(4)
         isorted(i12+nc(3)+i) = isortedtmp(i)
      enddo

c     Sort into boxes 5 and 6
      do i = i1234+1,i56
         if(pos(1,isorted(i))-center(1).lt.0) then
            nc(5) = nc(5) + 1
            isorted(i1234+nc(5)) = isorted(i)
         else
            nc(6) = nc(6) + 1
            isortedtmp(nc(6)) = isorted(i)
         endif
      enddo
c     Reorder sources so that sources in 6 are at the
c     end of this part of the array
      do i=1,nc(6)
         isorted(i1234+nc(5)+i) = isortedtmp(i)
      enddo
c     End of sorting sources into boxes 5 and 6

c     Sort into boxes 7 and 8
      do i=i56+1,npts
         if(pos(1,isorted(i))-center(1).lt.0) then
            nc(7) = nc(7) + 1
            isorted(i56+nc(7)) = isorted(i)
         else
            nc(8) = nc(8) + 1
            isortedtmp(nc(8)) = isorted(i)
         endif
      enddo
c     Reorder sources so that sources in 8 are at the
c     end of the array
      do i=1,nc(8)
         isorted(i56+nc(7)+i) = isortedtmp(i)
      enddo
      
      istart=1
      iboxfl=0
      subcenters=0.0d0
      do i=1,8
         ii = 2
         jj = 2
         if(i.eq.1.or.i.eq.2.or.i.eq.5.or.i.eq.6) ii = 1
         if(i.lt.5) jj = 1
         if(nc(i).gt.0) then
            subcenters(1,i) = center(1)+(-1)**i*boxsize/2.0d0
            subcenters(2,i) = center(2)+(-1)**ii*boxsize/2.0d0
            subcenters(3,i) = center(3)+(-1)**jj*boxsize/2.0d0

            iboxfl(1,i) = istart
            iboxfl(2,i) = istart+nc(i)-1

            istart = istart+nc(i)
          endif
      enddo

      return
      end
c------------------------------------------------------------------      
c--------------------------------------------------------------------     
c
c
c
c--------------------------------------------------------------------      

      subroutine getlist4pwdir(dir,censrc,centrg,boxsize)
c-------------------------------------------------------------------
      implicit none
ccc   input/output variables
      integer dir
      double precision censrc(3)
      double precision centrg(3)
      double precision boxsize
ccc   scopded function variables
      double precision sepdist
      double precision ctmp(3)
      ctmp(1) = censrc(1)-0*boxsize/2.0d0
      ctmp(2) = censrc(2)-0*boxsize/2.0d0
      ctmp(3) = censrc(3)-0*boxsize/2.0d0

      sepdist=1.51d0*boxsize

      if(ctmp(3)-centrg(3).ge.sepdist) then
        dir=1
      else if(ctmp(3)-centrg(3).le.-sepdist) then
        dir=2
      else if(ctmp(2)-centrg(2).ge.sepdist) then
        dir=3
      else if(ctmp(2)-centrg(2).le.-sepdist) then
        dir=4
      else if(ctmp(1)-centrg(1).ge.sepdist) then
        dir=5
      else if(ctmp(1)-centrg(1).le.-sepdist) then
        dir=6
      else
        dir=0
      end if
C      dir=0

      return
      end
      
      subroutine getlist4pwdirtest(dir,censrc,centrg,boxsize)
c-------------------------------------------------------------------
      implicit none
ccc   input/output variables
      integer dir
      double precision censrc(3)
      double precision centrg(3)
      double precision boxsize
ccc   scopded function variables
      double precision sepdist
      double precision ctmp(3)
      ctmp(1) = censrc(1)-0*boxsize/2.0d0
      ctmp(2) = censrc(2)-0*boxsize/2.0d0
      ctmp(3) = censrc(3)-0*boxsize/2.0d0

      sepdist=1.51d0*boxsize

      if((ctmp(3)-centrg(3)).ge.sepdist) then
        dir=1
      else if((ctmp(3)-centrg(3)).le.-sepdist) then
        dir=2
      else if((ctmp(2)-centrg(2)).ge.sepdist) then
        dir=3
      else if((ctmp(2)-centrg(2)).le.-sepdist) then
        dir=4
      else if((ctmp(1)-centrg(1)).ge.sepdist) then
        dir=5
      else if((ctmp(1)-centrg(1)).le.-sepdist) then
        dir=6
      else
        dir=0
        print *,"dir:",dir
      end if
C      dir=0

      return
      end
