c
c
c        routines in this file
c         
c           ccubeints_adap - completely adaptive integration
c                       for the quad patch
c
c
c         We integrate against legendre polynomials on the standard
c          patch [-1,1]^3
c                       
c
c
c
c
c
c
c
c
c

      subroutine ccubeints_adap(eps,
     1     norder,ttype,npols,ntarg,xyztarg,ncmax,
     3     fker,dpars,zpars,ipars,nqorder,cintvals)

c
c       this subroutine computes the integrals
c
c       \int_{[-1,1]^3} K(x_{i},y) P_{n}(y_{1}) P_{m}(y_{2})P_{m}(y_{3}) 
c          dy \, ,
c
c        P_{n}(y) are Legendre polynomials on [-1,1]
c
c        using adaptive integration, no precomputation
c        and storage requirements
c
c
c        input arguments:
c        eps:     requested precision 
c
c        norder: order of polynomials on the patch 
c        npols = norder*norder*norder if (type=f)
c                norder*(norder+1)*(norder+2)/6 if(type=t)
c
c        ttype = total degree ('t') or full degree ('f') 
c          polynomials to be integrated
c
c        ntarg - total number of target points
c        xyztarg(3,ntarg) - location of target points 
c        ncmax - max number of cubes in adaptive integration 
c        fker - function handle for evaluating the kernel k
c 
c               expected calling sequence
c               fker(x,y,dpars,zpars,ipars,f)
c               x \in \mathbb{R}^{3}, y \in \mathbb{R}^{3}
c               the output is complex *16 

c         dpars(*) - real parameters for the fker routine
c         zpars(*) - complex parameters for the fker routine
c         ipars(*) - integer parameters for the fker routine
c         nqorder - order of quadrature nodes on each subquad
c                   to be used
c
c         output:
c
c         cintvals(npols,ntarg) - integrals at all targets
c                                  for all tensor product
c                                  legendre polynomials
c
c
c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer norder,npols
      
      integer ntarg
      real *8 xyztarg(3,ntarg)
      
      external fker
      real *8 dpars(*)
      complex *16 zpars(*)
      integer ipars(*)

      integer nqorder

      integer ncmax

      complex *16 cintvals(npols,ntarg)

      character ttype

c
c       tree variables
c
      integer nlmax,ltree
      real *8, allocatable :: tvs(:,:,:),da(:)
      integer, allocatable :: ichild_start(:)

      integer ncube,nlev,icube,istart,i,j,k
      integer ier,itarg,jj,jstart,npts
      integer iqquad,ii


      integer npmax

      real *8, allocatable :: uvsq(:,:),wts(:),uvtmp(:,:)
      real *8, allocatable :: umattmp(:,:),vmattmp(:,:)
      integer nqpols
      real *8, allocatable :: sigvals(:,:)
      real *8, allocatable :: uvvals(:,:),qwts(:)
      integer itmp

      character *1 transa,transb
      real *8 alpha,beta,ra
      integer lda,ldb,ldc
      real *8 u, v
      integer ldu, ldv, itype,ndeg
      
c
c       for each quad, we just store three pieces
c       of info
c         quad vertices
c         area of quad
c         ichild_start = index for first child, 
c         
c
      allocate(ichild_start(ncmax),tvs(3,4,ncmax))
      allocate(da(ncmax))

      do i=1,ncmax
        ichild_start(i) = -1
        da(i) = 0
        do j=1,3
          do k=1,2
            tvs(k,j,i) = 0
          enddo
        enddo
      enddo

      da(1) = 1.0d0
c
c
c        quad vertices nomenclature
c
c       v3
c        ________ 
c        |       |
c        |       |
c        |       |
c        ---------
c        v1       v2
c
c   and v4 is the +len in z direction from v1
c
c
      tvs(1,1,1) = -1
      tvs(2,1,1) = -1
      tvs(3,1,1) = -1

      tvs(1,2,1) = 1
      tvs(2,2,1) = -1
      tvs(3,2,1) = -1

      tvs(1,3,1) = -1
      tvs(2,3,1) = 1
      tvs(3,3,1) = -1

      tvs(1,4,1) = -1
      tvs(2,4,1) = -1
      tvs(3,4,1) = 1



c
c       get quadrature nodes and weights on the base quad
c       based on quadrature type
c

      nqpols = nqorder*nqorder*nqorder
      allocate(uvsq(3,nqpols),wts(nqpols))

      ldu = 1
      ldv = 1
      itype = 1
      call legetens_exps_3d(itype,nqorder,ttype,uvsq,u,ldu,v,ldv,wts) 


      allocate(uvtmp(3,nqpols))
     
      npmax = ncmax*nqpols
      allocate(sigvals(npols,nqpols))
      allocate(uvvals(3,nqpols),qwts(nqpols))

cs
c      current number of cubes in the adaptive structure
c


c
cc       for the current patch compute geometry info for base quad 
c

       
      nlmax = 20 
      do itarg=1,ntarg
        ncube = 1
c
c        intialize sigvals for root quad
c

        ndeg = norder-1
        call mapuv_cube(tvs(1,1,1),nqpols,uvsq,uvvals)
        do i=1,nqpols
          call legetens_pols_3d(uvvals(1,i),ndeg,ttype,sigvals(1,i))
          qwts(i) = wts(i)
        enddo
        
        call cubeadap(eps,nqorder,nqpols,nlmax,ncmax,ncube,
     1    ichild_start,tvs,da,uvsq,wts, 
     1    norder,ttype,npols,npmax,uvvals,qwts,sigvals,xyztarg(1,itarg),
     3    fker,dpars,zpars,ipars,cintvals(1,itarg))
      enddo


      return
      end

c
c
c
c
c
      subroutine cubeadap(eps,m,kpols,nlmax,ncmax,ncube,
     1             ichild_start,tvs,da,uvsq,wts,
     1             norder,ttype,npols,npmax,uvvals,qwts,
     2             sigvals,xt,fker,dpars,zpars,
     3             ipars,cintall)

c
c       this subroutine adaptively computes the integral
c        of the functions
c   
c        \int_{(-1,1)^3} 1/|xt- y(u,v,w)| P_{n}(u)P_{m}(v)P_{\ell}(w) 
c            du dv dw
c        n,m,\ell = 1,2\ldots npols
c
c        P_{n}(u) are the Legendre polynomials on (-1,1) 
c
c
c        relevant parameters on a cube refined
c        grid along the way
c
c        IN:
c        eps - precision requested
c        m - quadrature order
c        kpols - number of quarature nodes (kpols = m*m*m) 
c        nlmax - max level of refinement for geometry
c        ncmax - max number of quads
c        ncube - current number of quads in adaptive structure
c        ichild_start(i) - first child of quad i
c        tvs(3,4,ncmax) - vertices of hierarchy of quads
c        da(ncmax) - area of quads
c        uvsq(kpols) - integration nodes on standard quad
c        wts(kpols) - integration weights on standard quad
c        npols - total number of koornwinder polynomials to be integrated
c        norder - order of polynomials to be integrated 
c        ttype - ttype (full or total) 
c        npols - norder*norder*norder if ttype =f, 
c                norder*(norder+1)*(norder+2)/6 if ttype =t
c        npmax - max number of points = ncmax*kpols
c        uvvals(3,kpols) - geometry info on heirarchy of meshes
c        qwts(kpols) - quadrature weights 
c        sigvals(npols,kpols) - 
c                   tensor product GL polynomials computed along the adaptive grid
c        
c        OUT:
c        cintall(npols) - computed integral 
c
c         

      implicit real *8 (a-h,o-z)
      integer, allocatable :: istack(:)
      integer ichild_start(ncmax)
      real *8 da(ncmax)
      real *8 tvs(3,4,ncmax), uvsq(3,kpols),wts(kpols)
      integer nproclist0, nproclist
      integer idone
      real *8 sigvals(npols,kpols)
      real *8 uvvals(3,kpols),qwts(kpols)
      complex *16, allocatable :: xkernvals(:)
      real *8 xt(3),xs(3)
      complex *16 cintall(npols),fval,ctmp(npols)
      complex *16, allocatable :: cvals(:,:)

      real *8 dpars(*)
      complex *16 zpars(*)
      integer ipars(*)

      character ttype
      
      external fker

c
c         for historic reasons
c
      ksigpols = npols
      allocate(istack(ncmax))
      nproclist0 = 1
      istack(1) = 1
      allocate(cvals(ksigpols,ncmax))

      do i=1,ksigpols
         cvals(i,1) = 0
      enddo

      allocate(xkernvals(kpols))

c
cc      compute integral at level 0
c
      do i=1,kpols
         call fker(uvvals(1,i),xt,dpars,zpars,ipars,fval)
         xkernvals(i) = fval*qwts(i)
         do j=1,ksigpols
            cvals(j,1) = cvals(j,1)+xkernvals(i)*sigvals(j,i)
         enddo
      enddo

      
      do i=1,ksigpols
        cintall(i) = cvals(i,1)
      enddo



      call cubeadap_main(eps,kpols,nlmax,ncmax,ncube,ichild_start,
     1      tvs,da,uvsq,wts,norder,ttype,npols,
     2      npmax,uvvals,qwts,sigvals,xt,fker,dpars,
     3      zpars,ipars,cvals,istack,nproclist0,
     4      xkernvals,cintall)

      
      return
      end
c
c
c
c
c
       
      subroutine cubeadap_main(eps,kpols,nlmax,ncmax,ncube,
     1    ichild_start,tvs,da,uvsq,wts,norder,ttype,npols,
     2    npmax,uvvals,qwts,sigvals,xt,fker,dpars,
     3    zpars,ipars,cvals,istack,nproclist0,xkernvals,
     4    cintall)
      

      implicit real *8 (a-h,o-z)
      integer istack(*),nproclist0
      integer ichild_start(ncmax)
      real *8 da(ncmax)
      real *8 tvs(3,4,ncmax), uvsq(3,kpols),wts(kpols)
      real *8 uvvals(3,kpols),qwts(kpols)
      integer  nproclist
      integer idone
      real *8 sigvals(npols,kpols)
      complex *16 xkernvals(kpols)
      real *8 xt(3)
      complex *16 cintall(npols),fval,ctmp(npols)
      complex *16 cvals(npols,ncmax)

      character ttype

      real *8 dpars(*)
      complex *16 zpars(*)
      integer ipars(*)

      character *1 transa,transb
      integer lda,ldb,ldc
      external fker
      


c
c         for historic reasons
c
      ksigpols = npols
      kfine = 8*kpols


      do ilev=0,nlmax
        idone = 1
        nproclist = 0
        
        do iproc = 1,nproclist0
          icube = istack(iproc)

          if(ncube+8.gt.ncmax) then
            print *, "Too many quads in cquadadap"
            print *, "Exiting without computing anything"

            return
          endif
            
          ichild_start(icube) = ncube+1
          call getcubechildren(tvs(1,1,icube),tvs(1,1,ncube+1),
     1           tvs(1,1,ncube+2),tvs(1,1,ncube+3),tvs(1,1,ncube+4),
     2           tvs(1,1,ncube+5),tvs(1,1,ncube+6),tvs(1,1,ncube+7),
     3           tvs(1,1,ncube+8))
            

          ndeg = norder-1
c
cc         subtract contribution of current quad
c          
          do isig=1,ksigpols
            cintall(isig) = cintall(isig)-cvals(isig,icube)
            ctmp(isig) = 0
          enddo


          rr = 0.125d0*da(icube)
          do j=ncube+1,ncube+8
            da(j) = rr
            call mapuv_cube(tvs(1,1,j),kpols,uvsq,uvvals)
            do i=1,kpols
              call legetens_pols_3d(uvvals(1,i),ndeg,ttype,
     1               sigvals(1,i))
              call fker(uvvals(1,i),xt,dpars,
     1         zpars,ipars,fval)
              xkernvals(i) = fval*wts(i)*rr
            enddo
            do isig=1,ksigpols
              cvals(isig,j) = 0
            enddo

            do k=1,kpols
              do isig=1,ksigpols
                cvals(isig,j) = cvals(isig,j) + xkernvals(k)*
     1             sigvals(isig,k)
              enddo
            enddo
            do isig=1,ksigpols
              cintall(isig) = cintall(isig) + cvals(isig,j)
              ctmp(isig) = ctmp(isig)+cvals(isig,j)
            enddo
          enddo


          errmax = 0
          do isig=1,ksigpols
            if(abs(ctmp(isig)-cvals(isig,icube)).gt.errmax) 
     1          errmax = abs(ctmp(isig)-cvals(isig,icube))
          enddo

          if(errmax.gt.eps) then
            idone = 0
            do j=1,8
              istack(nproclist0+nproclist+j) = ncube+j
            enddo
            nproclist = nproclist+8
          endif
          ncube = ncube+8
cc        end of looping over all quads at current stage
        enddo
cc         if idone is still 1, that means that no more refinement
c          is needed
         if(idone.eq.1) goto 1111
         
        do i=1,nproclist
          istack(i) = istack(nproclist0+i)
        enddo
        nproclist0 = nproclist
      enddo
 1111 continue



      return
      end
c
c
c
c
c
c
c
c
c--------------------------------------------------------------------------------
        
      subroutine mapuv_cube(verts,kpols,uvs,uvout)
      implicit real *8 (a-h,o-z)
      integer kpols
      real *8 verts(3,4),uvs(3,kpols),uvout(3,kpols)

      dx = verts(1,2)-verts(1,1)
      dy = verts(2,3)-verts(2,1) 
      dz = verts(3,4)-verts(3,1)

      do i=1,kpols
        uvout(1,i) = verts(1,1) + dx*(uvs(1,i)+1)/2
        uvout(2,i) = verts(2,1) + dy*(uvs(2,i)+1)/2
        uvout(3,i) = verts(3,1) + dz*(uvs(3,i)+1)/2
      enddo

      return
      end
c-----------------------------------------      
      
c
c
c

      subroutine getcubechildren(v0,v1,v2,v3,v4,v5,v6,v7,v8)
c  
cc       given the four vertices of a cube v0,
c        this subroutine returns the vertices of 8 
c        smaller cubes constructed using the
c        midpoints of the quads
c 
c        input:
c        v0 - real *8 (3,4)  
c              vertices of parent cube
c
c        output:
c        v1,v2,v3,v4,v5,v6,v7,v8 - real *8 (3,4)
c                 vertices of children quads
      
      implicit real *8 (a-h,o-z)
      real *8 v0(3,4),v1(3,4),v2(3,4),v3(3,4),v4(3,4)
      real *8 v5(3,4),v6(3,4),v7(3,4),v8(3,4)

      dx = v0(1,2)-v0(1,1)
      dy = v0(2,3)-v0(2,1)
      dz = v0(3,4)-v0(3,1)


c
cc     first cube
c
      v1(1,1) = v0(1,1) 
      v1(2,1) = v0(2,1)
      v1(3,1) = v0(3,1)

      call cube_vert_to_dxyz(v1,dx,dy,dz)
     
c
cc      second cube
c
      v2(1,1) = v0(1,1)+dx/2 
      v2(2,1) = v0(2,1)
      v2(3,1) = v0(3,1)

      call cube_vert_to_dxyz(v2,dx,dy,dz)

c
cc      third cube
c
      v3(1,1) = v0(1,1) 
      v3(2,1) = v0(2,1)+dy/2
      v3(3,1) = v0(3,1)

      call cube_vert_to_dxyz(v3,dx,dy,dz)

c
cc      fourth cube
c
      v4(1,1) = v0(1,1)+dx/2 
      v4(2,1) = v0(2,1)+dy/2
      v4(3,1) = v0(3,1)

      call cube_vert_to_dxyz(v4,dx,dy,dz)


c
cc     fifth cube
c
      v5(1,1) = v0(1,1) 
      v5(2,1) = v0(2,1)
      v5(3,1) = v0(3,1)+dz/2

      call cube_vert_to_dxyz(v5,dx,dy,dz)
     
c
cc      sixth cube
c
      v6(1,1) = v0(1,1)+dx/2 
      v6(2,1) = v0(2,1)
      v6(3,1) = v0(3,1)+dz/2

      call cube_vert_to_dxyz(v6,dx,dy,dz)

c
cc      seventh cube
c
      v7(1,1) = v0(1,1) 
      v7(2,1) = v0(2,1)+dy/2
      v7(3,1) = v0(3,1)+dz/2

      call cube_vert_to_dxyz(v7,dx,dy,dz)

c
cc      eighth cube
c
      v8(1,1) = v0(1,1)+dx/2 
      v8(2,1) = v0(2,1)+dy/2
      v8(3,1) = v0(3,1)+dz/2

      call cube_vert_to_dxyz(v8,dx,dy,dz)



      return
      end

c----------------------------------


      subroutine cube_vert_to_dxyz(v1,dx,dy,dz)
      implicit none
      real *8 v1(3,4),dx,dy,dz

      v1(1,2) = v1(1,1) + dx/2
      v1(2,2) = v1(2,1)
      v1(3,2) = v1(3,1)

      v1(1,3) = v1(1,1)
      v1(2,3) = v1(2,1) + dy/2
      v1(3,3) = v1(3,1)

      v1(1,4) = v1(1,1)
      v1(2,4) = v1(2,1)
      v1(3,4) = v1(3,1) + dz/2


      return
      end
