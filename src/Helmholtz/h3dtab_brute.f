      subroutine mksurhelm3dp(x0,y0,z0,ix,iy,iz,zk0,norder0,
     1     value,ifail)
c
c       This subroutine is the wrapper for computing
c       the integrals  
c     \int_{[-1,1]^3} e^{ikr}/r *
c          P_{ix-1}(x)*P_{iy-1}(y)*P_{iz-1}(z) dx dy dz \, ,
c        
c       where P_{n}(x) are legendre polynomials
c       and r = \sqrt{(x-x_{0})^2 + (y-y_{0})^2 + (z-z_{0})^2} 
c
c       Be careful, there are a few global variables
c
c       The kernel can be changed by appropriately changing
c       the function "fhelmgreen3d"
c
c       Input parameters:
c          x0,y0,z0 - coordinates of target location
c          ix,iy,iz - polynomial order in x,y, and z variables
c          norder - order of box code generation. 
c                   Must be greater than max(ix,iy,iz)
c       
c       Output parameters:
c          value - value of integral
c          ifail - ifail = 0 for successful computation of integral
c                     check other routines for error code otherwise
c          
c
c
c
c

      implicit none
      external fhelm3dp
      integer key, n, nf, ndim, mincls, maxcls, ifail, neval, nw
      parameter (ndim = 3, nw = 4000000, nf = 2)
      real *8 a(ndim), b(ndim)
      real *8, allocatable :: wrkstr(:)
      real *8 absest(nf), finest(nf), absreq, relreq
      real *8 xtarg,ytarg,ztarg
      real *8 x0, y0, z0
      complex *16 value, zk, zk0
      complex *16 im
      data im / (0.0d0,1.0d0) /
      integer ix,iy,iz,ixpol,iypol,izpol,norder,norder0
      common /cbh3dtab_brute/ xtarg,ytarg,ztarg,ixpol,iypol,izpol,
     1     norder,zk
c$omp threadprivate(/cbh3dtab_brute/)      

      xtarg = x0
      ytarg = y0
      ztarg = z0

      allocate(wrkstr(nw))

      ixpol = ix
      iypol = iy
      izpol = iz

      norder = norder0
      zk = zk0


      do 10 n = 1,ndim
         a(n) = -1.0d0
         b(n) =  1.0d0
   10 continue
      mincls = 0
      maxcls = 4000000
      key = 0
      absreq = 1d-12
      relreq = 1d-12
      ifail = 0

      call dcuhre(ndim, nf, a, b, mincls, maxcls, fhelm3dp, 
     1      absreq, relreq, key, nw, 0, finest, absest, neval,
     2      ifail, wrkstr)


      value = finest(1) + im*finest(2)

      return
      end


      subroutine fhelm3dp(ndim, z, nfun, f)
      implicit none
      integer ndim, nfun
      real *8 z(ndim), f(nfun), rx, ry, rz
      real *8 rr,reps,dgreen
      real *8 xtarg,ytarg,ztarg
      integer ixpol,iypol,izpol,norder
      real *8 polsx(norder),polsy(norder),polsz(norder)
      real *8 h2
      complex *16 im, zf, zk
      data im / (0.0d0,1.0d0) /
      common /cbh3dtab_brute/ xtarg,ytarg,ztarg,ixpol,iypol,izpol,
     1     norder,zk
c$omp threadprivate(/cbh3dtab_brute/)      
      
      
      rx = z(1) - xtarg
      ry = z(2) - ytarg
      rz = z(3) - ztarg

      call legepols(z(1),norder-1,polsx)
      call legepols(z(2),norder-1,polsy)
      call legepols(z(3),norder-1,polsz)



      rr = dsqrt(rx*rx + ry*ry + rz*rz)

      zf = exp(im*zk*rr)*polsx(ixpol)*polsy(iypol)*polsz(izpol)/rr
      f(1) = dreal(zf)
      f(2) = dimag(zf)

      return
      end

      
