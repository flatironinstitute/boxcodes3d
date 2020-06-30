
c
c    this file contains subroutines for a lippman schwinger
c     solver for the Helmholtz equation 
c
c

      subroutine ls_solver_guru(eps,zk,nboxes,nlevels,ltree,itree,iptr,
     1   norder,ncbox,ttype,qvals,centers,boxsize,npbox,
     2   rhs,irep,eps_gmres,numit,niter,errs,rres,sigma)
c
c      this subroutine is a solver for the lippman schwinger solver 
c
c      input arguments:
c         eps - real *8 
c           precision requested
c         zk - double complex
c            Helmholtz wave number
c         nboxes - integer
c            number of boxes
c         nlevels - integer
c            number of levels
c         ltree - integer
c            length of array containing the tree structure
c         itree - integer(ltree)
c            array containing the tree structure
c         iptr - integer(8)
c            pointer to various parts of the tree structure
c           iptr(1) - laddr
c           iptr(2) - ilevel
c           iptr(3) - iparent
c           iptr(4) - nchild
c           iptr(5) - ichild
c           iptr(6) - ncoll
c           iptr(7) - coll
c           iptr(8) - ltree
c         norder - integer
c           order of expansions for input coefficients array
c         ncbox - integer
c           number of coefficients of expansions of functions
c           in each of the boxes
c         ttype - character *1
c            type of coefs provided, total order ('t') or full order('f')
c         qvals - double complex (npbox,nboxes)
c            material property tabulated on the tree
c         rhs - double complex (npbox,nboxes)
c            right hand side on the tree
c         centers - double precision (3,nboxes)
c           xyz coordintes of boxes in the tree structure
c         boxsize - double precision (0:nlevels)
c           size of boxes at each of the levels
c         npbox - integer
c           number of points per box where potential is to be dumped = (norder**3)
c         irep - integer
c           representation of lippman schwinger
c             if irep = 1, u = V_{k} \sigma
c               use \sigma + k^2 q V_{k} \sigma = f
c             if irep = 2,
c                u + V_{k} [k^2 q u] = V_{k}[f] 
c         eps_gmres - real *8
c           gmres tolerance requested
c         numit - integer
c           max number of gmres iterations
c 
c     output:
c        niter - integer
c          number of gmres iterations
c        errs(1:iter) - relative residual as a function of iteration
c          number
c        rres - real *8
c          relative residual as a function of iteration number
c        sigma - complex *16 (npbox,nboxes)
c          solution of the integral equation
c
c
      implicit real *8 (a-h,o-z)
      real *8 eps
      complex *16 zk
      integer nboxes,nlevels,ltree,itree(ltree),iptr(8)
      integer norder,ncbox,npbox
      real *8 errs(numit+1)
      character *1 ttype
      complex *16 qvals(npbox,nboxes),rhs(npbox,nboxes)
      real *8 centers(3,nboxes),boxsize(0:nlevels)
      real *8 tprecomp(3)
      integer irep
      complex *16 sigma(npbox,nboxes)
      complex *16 temp,ztmp
c
c       fmm precomp variables
c
c
      complex *16, allocatable :: rhsuse(:,:)
      complex *16, allocatable :: mpcoefsmat(:),tamat(:),tab(:)
      integer impcoefsmat(0:nlevels+1),itamat(0:nlevels+1)
      integer itab(0:nlevels+1)
      integer lmpcoefsmat,ltamat,ltab

c
c        gmres variables
c
      complex *16, allocatable :: vmat(:,:,:),hmat(:,:)
      complex *16, allocatable :: cs(:),sn(:)
      complex *16, allocatable :: svec(:),yvec(:),wtmp(:,:)
      complex *16, allocatable :: wtmp2(:,:)


      done =  1
      pi = atan(done)*4
      

      allocate(vmat(npbox,nboxes,numit+1),hmat(numit,numit))
      allocate(cs(numit),sn(numit))
      allocate(wtmp(npbox,nboxes),svec(numit+1),yvec(numit+1))
      allocate(wtmp2(npbox,nboxes))

c
c       get precomputation arrays for current tree structure
c
c
      call helmholtz_volume_fmm_mem(eps,zk,nboxes,nlevels,ltree,
     1   itree,iptr,boxsize,centers,npbox,ncbox,impcoefsmat,
     2   lmpcoefsmat,itamat,ltamat,itab,ltab)
      
      allocate(mpcoefsmat(lmpcoefsmat),tamat(ltamat),tab(ltab))
      
      call helmholtz_volume_fmm_init(eps,zk,nboxes,nlevels,boxsize,
     1  norder,npbox,ncbox,impcoefsmat,lmpcoefsmat,itamat,ltamat,
     2  itab,ltab,ttype,mpcoefsmat,tamat,tab,tprecomp)

      call prinf('finished precomp*',i,0)

      allocate(rhsuse(npbox,nboxes))

      if(irep.eq.1) then
        do ibox=1,nboxes
          do j=1,npbox
            rhsuse(j,ibox) = rhs(j,ibox)
          enddo
        enddo
      endif
      if(irep.eq.2) then
        call helmholtz_volume_fmm_wprecomp(eps,zk,nboxes,nlevels,
     1  ltree,itree,iptr,norder,ncbox,ttype,rhs,centers,
     2  boxsize,mpcoefsmat,impcoefsmat,lmpcoefsmat,tamat,itamat,
     3  ltamat,tab,itab,ltab,npbox,rhsuse,timeinfo)
       endif

c
c
c     start gmres code here
c
c     NOTE: matrix equation should be of the form (z*I + K)x = y
c       the identity scaling (z) is defined via zid below,
c       and K represents the action of the principal value 
c       part of the matvec
c
      zid = -4*pi
      print *, zid


      niter=0

c
c      compute norm of right hand side and initialize v
c 
      rb = 0

      do i=1,numit
        cs(i) = 0
        sn(i) = 0
      enddo
c
      do ibox=1,nboxes
        if(itree(iptr(4)+ibox-1).eq.0) then
          do j=1,npbox
            rb = rb + abs(rhsuse(j,ibox))**2
          enddo
        endif
      enddo
      rb = sqrt(rb)

      do ibox=1,nboxes
        do j=1,npbox
          vmat(j,ibox,1) = rhsuse(j,ibox)/rb
        enddo
      enddo

      svec(1) = rb

      do it=1,numit
        it1 = it + 1

        print *, "iter number=",it

c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c

        if(irep.eq.2) then
          do ibox=1,nboxes
            do j=1,npbox
              wtmp2(j,ibox) = vmat(j,ibox,it)*zk**2*qvals(j,ibox)
            enddo
          enddo
        endif

        if(irep.eq.1) then
          do ibox=1,nboxes
            do j=1,npbox
              wtmp2(j,ibox) = vmat(j,ibox,it)
            enddo
          enddo
        endif

        do ibox=1,nboxes
          do j=1,npbox
            wtmp(j,ibox) = 0
          enddo
        enddo


        call helmholtz_volume_fmm_wprecomp(eps,zk,nboxes,nlevels,
     1  ltree,itree,iptr,norder,ncbox,ttype,wtmp2,centers,
     2  boxsize,mpcoefsmat,impcoefsmat,lmpcoefsmat,tamat,itamat,
     3  ltamat,tab,itab,ltab,npbox,wtmp,timeinfo)


        if(irep.eq.1) then
          do ibox=1,nboxes
            do j=1,npbox
              wtmp(j,ibox) = wtmp(j,ibox)*zk**2*qvals(j,ibox)
            enddo
          enddo
        endif

        do k=1,it
          hmat(k,it) = 0
          do ibox=1,nboxes
            if(itree(iptr(4)+ibox-1).eq.0) then
              do j=1,npbox
                hmat(k,it) = hmat(k,it) +
     1            wtmp(j,ibox)*conjg(vmat(j,ibox,k))
              enddo
            endif
          enddo

          do ibox=1,nboxes
            if(itree(iptr(4)+ibox-1).eq.0) then
              do j=1,npbox
                wtmp(j,ibox) = wtmp(j,ibox)-hmat(k,it)*vmat(j,ibox,k)
              enddo
            endif
          enddo
        enddo
          
        hmat(it,it) = hmat(it,it)+zid
        wnrm2 = 0
        do ibox=1,nboxes
          if(itree(iptr(4)+ibox-1).eq.0) then
            do j=1,npbox
              wnrm2 = wnrm2 + abs(wtmp(j,ibox))**2
            enddo
          endif
        enddo
        wnrm2 = sqrt(wnrm2)
        
        do ibox=1,nboxes
          if(itree(iptr(4)+ibox-1).eq.0) then
            do j=1,npbox
              vmat(j,ibox,it1) = wtmp(j,ibox)/wnrm2
            enddo
          else
            do j=1,npbox
              vmat(j,ibox,it1) = 0
            enddo
          endif
        enddo

        do k=1,it-1
          temp = cs(k)*hmat(k,it)+sn(k)*hmat(k+1,it)
          hmat(k+1,it) = -sn(k)*hmat(k,it)+cs(k)*hmat(k+1,it)
          hmat(k,it) = temp
        enddo

        ztmp = wnrm2

        call zrotmat_gmres(hmat(it,it),ztmp,cs(it),sn(it))
          
        hmat(it,it) = cs(it)*hmat(it,it)+sn(it)*wnrm2
        svec(it1) = -sn(it)*svec(it)
        svec(it) = cs(it)*svec(it)
        rmyerr = abs(svec(it1))/rb
        errs(it) = rmyerr
        print *, "iter=",it,errs(it)

        if(rmyerr.le.eps_gmres.or.it.eq.numit) then

c
c            solve the linear system corresponding to
c            upper triangular part of hmat to obtain yvec
c
c            y = triu(H(1:it,1:it))\s(1:it);
c
          do j=1,it
            iind = it-j+1
            yvec(iind) = svec(iind)
            do l=iind+1,it
              yvec(iind) = yvec(iind) - hmat(iind,l)*yvec(l)
            enddo
            yvec(iind) = yvec(iind)/hmat(iind,iind)
          enddo



c
c          estimate x
c
          do ibox=1,nboxes
            if(itree(iptr(4)+ibox-1).eq.0) then
              do j=1,npbox
                sigma(j,ibox) = 0
                do i=1,it
                  sigma(j,ibox) = sigma(j,ibox) + yvec(i)*vmat(j,ibox,i)
                enddo
              enddo
            else
              do j=1,npbox
                sigma(j,ibox) = 0
              enddo
            endif
          enddo


          rres = 0
          do ibox=1,nboxes
            do j=1,npbox
              wtmp(j,ibox) = 0
            enddo
          enddo
c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c


          if(irep.eq.2) then
            do ibox=1,nboxes
              do j=1,npbox
                wtmp2(j,ibox) = sigma(j,ibox)*zk**2*qvals(j,ibox)
              enddo
            enddo
          endif

          if(irep.eq.1) then
            do ibox=1,nboxes
              do j=1,npbox
                wtmp2(j,ibox) = sigma(j,ibox)
              enddo
            enddo
          endif
          call helmholtz_volume_fmm_wprecomp(eps,zk,nboxes,nlevels,
     1     ltree,itree,iptr,norder,ncbox,ttype,wtmp2,centers,
     2     boxsize,mpcoefsmat,impcoefsmat,lmpcoefsmat,tamat,itamat,
     3     ltamat,tab,itab,ltab,npbox,wtmp,timeinfo)
       
          if(irep.eq.1) then
            do ibox=1,nboxes
              do j=1,npbox
                wtmp(j,ibox) = wtmp(j,ibox)*zk**2*qvals(j,ibox)
              enddo
            enddo
          endif

          do ibox=1,nboxes
            if(itree(iptr(4)+ibox-1).eq.0) then
              do j=1,npbox
                rres = rres + abs(zid*sigma(j,ibox) + 
     1            wtmp(j,ibox)-rhsuse(j,ibox))**2
              enddo
            endif
          enddo

          rres = sqrt(rres)/rb
          niter = it
          return
        endif
      enddo
c

      return
      end

