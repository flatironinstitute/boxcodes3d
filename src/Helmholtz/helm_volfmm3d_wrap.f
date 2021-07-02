      subroutine hvfmm_v_v_p(eps,zk,nboxes,nlevels,ltree,itree,
     1    iptr,norder,ncbox,ttype,fvals,centers,boxsize,npbox,
     2    pot)
c--------------------------     
c  This code applies the Helmholtz volume layer potential
c     \int_{B} exp(ik|x-y|)/|x-y| f(y) where
c  f is defined as samples at tensor product legendre nodes 
c  on a level restricted tree and the potential is returned
c  as samples on the same grid.
c 
c  Input arugments
c    - eps: double precision
c       tolerance requested
c    - zk: double complex
c        Helmholtz wave number
c    - nboxes: integer
c        number of boxes
c    - nlevels: integer
c        number of levels
c    - ltree: integer
c        length of array containing the tree structure
c    - itree: integer(ltree)
c        array containing the tree structure
c    - iptr: integer(8)
c        pointer to various parts of the tree structure
c        * iptr(1) - laddr
c        * iptr(2) - ilevel
c        * iptr(3) - iparent
c        * iptr(4) - nchild
c        * iptr(5) - ichild
c        * iptr(6) - ncoll
c        * iptr(7) - coll
c        * iptr(8) - ltree
c    - norder: integer
c        order of expansions for input coefficients array
c    - ncbox: integer
c        number of coefficients of expansions of functions
c        in each of the boxes
c    - ttype: character *1
c        type of coefs provided, total order ('t') or full order('f')
c    - fvals: double complex (npbox,nboxes)
c        function tabulated on a grid
c    - centers: double precision (3,nboxes)
c        xyz coordintes of boxes in the tree structure
c    - boxsize: double precision (0:nlevels)
c        size of boxes at each of the levels
c    - npbox: integer
c        number of points per box where potential is to 
c        be computed = (norder**3)
c
c  Output arguments:
c    - pot: double complex (npbox,nboxes)
c        volume potential on the tree structure (note that 
c        the potential is non-zero only in the leaf boxes
c
      implicit real *8 (a-h,o-z)
      real *8 eps
      complex *16 zk
      integer nboxes,nlevels,ltree,itree(ltree),iptr(8)
      integer norder,ncbox
      character *1 ttype
      complex *16 fvals(npbox,nboxes)
      real *8 centers(3,nboxes),boxsize(0:nlevels)
      integer npbox
      complex *16 pot(npbox,nboxes)
      complex *16, allocatable :: potcoefs(:,:)
      real *8 timeinfo(6),tprecomp(3)
      complex *16, allocatable :: fcoefs(:,:)
      real *8, allocatable :: xref(:,:),umat(:,:),vmat(:,:),wts0(:)

      print *, "npbox=",npbox
      print *, "ncbox=",ncbox
      print *, "norder=",norder
      allocate(fcoefs(ncbox,nboxes),xref(3,npbox),umat(ncbox,npbox),
     1   vmat(npbox,ncbox),wts0(npbox))
     
      itype = 2 
      call legetens_exps_3d(itype,norder,ttype,xref,umat,ncbox,vmat,
     1  npbox,wts0)
      
      ra = 1.0d0
      rb = 0.0d0
      do ibox=1,nboxes
        call dgemm('n','t',2,ncbox,npbox,ra,fvals(1,ibox),
     1    2,umat,ncbox,rb,fcoefs(1,ibox),2)
      enddo
      
      allocate(potcoefs(ncbox,nboxes))
      call helmholtz_volume_fmm(eps,zk,nboxes,nlevels,ltree,itree,
     1  iptr,norder,ncbox,ttype,fcoefs,centers,boxsize,npbox,pot,
     2  potcoefs,timeinfo,tprecomp)

      return
      end
c
c
c
c
      subroutine hvfmm_v_c_p(eps,zk,nboxes,nlevels,ltree,itree,
     1    iptr,norder,ncbox,ttype,fvals,centers,boxsize,npbox,
     2    potcoefs)
c--------------------------     
c  This code applies the Helmholtz volume layer potential
c     \int_{B} exp(ik|x-y|)/|x-y| f(y) where
c  f is defined as samples at tensor product legendre nodes 
c  on a level restricted tree and the potential is returned
c  as polynomial expansion on the same grid.
c 
c  Input arugments
c    - eps: double precision
c       tolerance requested
c    - zk: double complex
c        Helmholtz wave number
c    - nboxes: integer
c        number of boxes
c    - nlevels: integer
c        number of levels
c    - ltree: integer
c        length of array containing the tree structure
c    - itree: integer(ltree)
c        array containing the tree structure
c    - iptr: integer(8)
c        pointer to various parts of the tree structure
c        * iptr(1) - laddr
c        * iptr(2) - ilevel
c        * iptr(3) - iparent
c        * iptr(4) - nchild
c        * iptr(5) - ichild
c        * iptr(6) - ncoll
c        * iptr(7) - coll
c        * iptr(8) - ltree
c    - norder: integer
c        order of expansions for input coefficients array
c    - ncbox: integer
c        number of coefficients of expansions of functions
c        in each of the boxes
c    - ttype: character *1
c        type of coefs provided, total order ('t') or full order('f')
c    - fvals: double complex (npbox,nboxes)
c        function tabulated on a grid
c    - centers: double precision (3,nboxes)
c        xyz coordintes of boxes in the tree structure
c    - boxsize: double precision (0:nlevels)
c        size of boxes at each of the levels
c    - npbox: integer
c        number of points per box where potential is to 
c        be computed = (norder**3)
c
c  Output arguments:
c    - potcoefs: double complex (ncbox,nboxes)
c        volume potential coefficients on the tree structure (note that 
c        the potential is non-zero only in the leaf boxes
c
      implicit real *8 (a-h,o-z)
      real *8 eps
      complex *16 zk
      integer nboxes,nlevels,ltree,itree(ltree),iptr(8)
      integer norder,ncbox
      character *1 ttype
      complex *16 fvals(npbox,nboxes)
      real *8 centers(3,nboxes),boxsize(0:nlevels)
      integer npbox
      complex *16 potcoefs(ncbox,nboxes)
      complex *16, allocatable :: pot(:,:)
      real *8 timeinfo(6),tprecomp(3)

      complex *16, allocatable :: fcoefs(:,:)
      real *8, allocatable :: xref(:,:),umat(:,:),vmat(:,:),wts0(:)

      allocate(fcoefs(ncbox,nboxes),xref(3,npbox),umat(ncbox,npbox),
     1   vmat(npbox,ncbox),wts0(npbox))
     
      itype = 2 
      call legetens_exps_3d(itype,norder,ttype,xref,umat,ncbox,vmat,
     1  npbox,wts0)
      
      ra = 1.0d0
      rb = 0.0d0
      do ibox=1,nboxes
        call dgemm('n','t',2,ncbox,npbox,ra,fvals(1,ibox),
     1    2,umat,ncbox,rb,fcoefs(1,ibox),2)
      enddo
      
      allocate(pot(npbox,nboxes))
      call helmholtz_volume_fmm(eps,zk,nboxes,nlevels,ltree,itree,
     1  iptr,norder,ncbox,ttype,fcoefs,centers,boxsize,npbox,pot,
     2  potcoefs,timeinfo,tprecomp)

      return
      end
c
c
c
c
      subroutine hvfmm_v_vc_p(eps,zk,nboxes,nlevels,ltree,itree,
     1    iptr,norder,ncbox,ttype,fvals,centers,boxsize,npbox,
     2    pot,potcoefs)
c--------------------------     
c  This code applies the Helmholtz volume layer potential
c     \int_{B} exp(ik|x-y|)/|x-y| f(y) where
c  f is defined as samples at tensor product legendre nodes 
c  on a level restricted tree and the potential is returned
c  as polynomial expansion and samples on the same grid.
c 
c  Input arugments
c    - eps: double precision
c       tolerance requested
c    - zk: double complex
c        Helmholtz wave number
c    - nboxes: integer
c        number of boxes
c    - nlevels: integer
c        number of levels
c    - ltree: integer
c        length of array containing the tree structure
c    - itree: integer(ltree)
c        array containing the tree structure
c    - iptr: integer(8)
c        pointer to various parts of the tree structure
c        * iptr(1) - laddr
c        * iptr(2) - ilevel
c        * iptr(3) - iparent
c        * iptr(4) - nchild
c        * iptr(5) - ichild
c        * iptr(6) - ncoll
c        * iptr(7) - coll
c        * iptr(8) - ltree
c    - norder: integer
c        order of expansions for input coefficients array
c    - ncbox: integer
c        number of coefficients of expansions of functions
c        in each of the boxes
c    - ttype: character *1
c        type of coefs provided, total order ('t') or full order('f')
c    - fvals: double complex (npbox,nboxes)
c        function tabulated on a grid
c    - centers: double precision (3,nboxes)
c        xyz coordintes of boxes in the tree structure
c    - boxsize: double precision (0:nlevels)
c        size of boxes at each of the levels
c    - npbox: integer
c        number of points per box where potential is to 
c        be computed = (norder**3)
c
c  Output arguments:
c    - pot: double complex (npbox,nboxes)
c        volume potential on the tree structure (note that 
c        the potential is non-zero only in the leaf boxes
c    - potcoefs: double complex (ncbox,nboxes)
c        volume potential coefficients on the tree structure (note that 
c        the potential is non-zero only in the leaf boxes
c
      implicit real *8 (a-h,o-z)
      real *8 eps
      complex *16 zk
      integer nboxes,nlevels,ltree,itree(ltree),iptr(8)
      integer norder,ncbox
      character *1 ttype
      complex *16 fvals(npbox,nboxes)
      real *8 centers(3,nboxes),boxsize(0:nlevels)
      integer npbox
      complex *16 potcoefs(ncbox,nboxes),pot(npbox,nboxes)
      real *8 timeinfo(6),tprecomp(3)
      

      complex *16, allocatable :: fcoefs(:,:)
      real *8, allocatable :: xref(:,:),umat(:,:),vmat(:,:),wts0(:)

      allocate(fcoefs(ncbox,nboxes),xref(3,npbox),umat(ncbox,npbox),
     1   vmat(npbox,ncbox),wts0(npbox))
     
      itype = 2 
      call legetens_exps_3d(itype,norder,ttype,xref,umat,ncbox,vmat,
     1  npbox,wts0)
      
      ra = 1.0d0
      rb = 0.0d0
      do ibox=1,nboxes
        call dgemm('n','t',2,ncbox,npbox,ra,fvals(1,ibox),
     1    2,umat,ncbox,rb,fcoefs(1,ibox),2)
      enddo
      
      call helmholtz_volume_fmm(eps,zk,nboxes,nlevels,ltree,itree,
     1  iptr,norder,ncbox,ttype,fcoefs,centers,boxsize,npbox,pot,
     2  potcoefs,timeinfo,tprecomp)

      return
      end
c
c
c
c
      subroutine hvfmm_c_v_p(eps,zk,nboxes,nlevels,ltree,itree,
     1    iptr,norder,ncbox,ttype,fcoefs,centers,boxsize,npbox,
     2    pot)
c--------------------------     
c  This code applies the Helmholtz volume layer potential
c     \int_{B} exp(ik|x-y|)/|x-y| f(y) where
c  f is defined as samples at tensor product legendre nodes 
c  on a level restricted tree and the potential is returned
c  as samples on the same grid.
c 
c  Input arugments
c    - eps: double precision
c       tolerance requested
c    - zk: double complex
c        Helmholtz wave number
c    - nboxes: integer
c        number of boxes
c    - nlevels: integer
c        number of levels
c    - ltree: integer
c        length of array containing the tree structure
c    - itree: integer(ltree)
c        array containing the tree structure
c    - iptr: integer(8)
c        pointer to various parts of the tree structure
c        * iptr(1) - laddr
c        * iptr(2) - ilevel
c        * iptr(3) - iparent
c        * iptr(4) - nchild
c        * iptr(5) - ichild
c        * iptr(6) - ncoll
c        * iptr(7) - coll
c        * iptr(8) - ltree
c    - norder: integer
c        order of expansions for input coefficients array
c    - ncbox: integer
c        number of coefficients of expansions of functions
c        in each of the boxes
c    - ttype: character *1
c        type of coefs provided, total order ('t') or full order('f')
c    - fcoefs: double complex (npbox,nboxes)
c        function coefs tabulated on a grid
c    - centers: double precision (3,nboxes)
c        xyz coordintes of boxes in the tree structure
c    - boxsize: double precision (0:nlevels)
c        size of boxes at each of the levels
c    - npbox: integer
c        number of points per box where potential is to 
c        be computed = (norder**3)
c
c  Output arguments:
c    - pot: double complex (npbox,nboxes)
c        volume potential on the tree structure (note that 
c        the potential is non-zero only in the leaf boxes
c
      implicit real *8 (a-h,o-z)
      real *8 eps
      complex *16 zk
      integer nboxes,nlevels,ltree,itree(ltree),iptr(8)
      integer norder,ncbox
      character *1 ttype
      complex *16 fcoefs(ncbox,nboxes)
      real *8 centers(3,nboxes),boxsize(0:nlevels)
      integer npbox
      complex *16 pot(npbox,nboxes)
      complex *16, allocatable :: potcoefs(:,:)
      real *8 timeinfo(6),tprecomp(3)
      
      allocate(potcoefs(ncbox,nboxes))
      call helmholtz_volume_fmm(eps,zk,nboxes,nlevels,ltree,itree,
     1  iptr,norder,ncbox,ttype,fcoefs,centers,boxsize,npbox,pot,
     2  potcoefs,timeinfo,tprecomp)

      return
      end
c
c
c
c
      subroutine hvfmm_c_c_p(eps,zk,nboxes,nlevels,ltree,itree,
     1    iptr,norder,ncbox,ttype,fcoefs,centers,boxsize,npbox,
     2    potcoefs)
c--------------------------     
c  This code applies the Helmholtz volume layer potential
c     \int_{B} exp(ik|x-y|)/|x-y| f(y) where
c  f is defined as samples at tensor product legendre nodes 
c  on a level restricted tree and the potential is returned
c  as polynomial expansion on the same grid.
c 
c  Input arugments
c    - eps: double precision
c       tolerance requested
c    - zk: double complex
c        Helmholtz wave number
c    - nboxes: integer
c        number of boxes
c    - nlevels: integer
c        number of levels
c    - ltree: integer
c        length of array containing the tree structure
c    - itree: integer(ltree)
c        array containing the tree structure
c    - iptr: integer(8)
c        pointer to various parts of the tree structure
c        * iptr(1) - laddr
c        * iptr(2) - ilevel
c        * iptr(3) - iparent
c        * iptr(4) - nchild
c        * iptr(5) - ichild
c        * iptr(6) - ncoll
c        * iptr(7) - coll
c        * iptr(8) - ltree
c    - norder: integer
c        order of expansions for input coefficients array
c    - ncbox: integer
c        number of coefficients of expansions of functions
c        in each of the boxes
c    - ttype: character *1
c        type of coefs provided, total order ('t') or full order('f')
c    - fcoefs: double complex (ncbox,nboxes)
c        function coeffs tabulated on a grid
c    - centers: double precision (3,nboxes)
c        xyz coordintes of boxes in the tree structure
c    - boxsize: double precision (0:nlevels)
c        size of boxes at each of the levels
c    - npbox: integer
c        number of points per box where potential is to 
c        be computed = (norder**3)
c
c  Output arguments:
c    - potcoefs: double complex (ncbox,nboxes)
c        volume potential coefficients on the tree structure (note that 
c        the potential is non-zero only in the leaf boxes
c
      implicit real *8 (a-h,o-z)
      real *8 eps
      complex *16 zk
      integer nboxes,nlevels,ltree,itree(ltree),iptr(8)
      integer norder,ncbox
      character *1 ttype
      complex *16 fcoefs(ncbox,nboxes)
      real *8 centers(3,nboxes),boxsize(0:nlevels)
      integer npbox
      complex *16 potcoefs(ncbox,nboxes)
      complex *16, allocatable :: pot(:,:)
      real *8 timeinfo(6),tprecomp(3)

      
      allocate(pot(npbox,nboxes))
      call helmholtz_volume_fmm(eps,zk,nboxes,nlevels,ltree,itree,
     1  iptr,norder,ncbox,ttype,fcoefs,centers,boxsize,npbox,pot,
     2  potcoefs,timeinfo,tprecomp)

      return
      end
c
c
c
c
c
      subroutine hvfmm_c_vc_p(eps,zk,nboxes,nlevels,ltree,itree,
     1    iptr,norder,ncbox,ttype,fcoefs,centers,boxsize,npbox,
     2    pot,potcoefs)
c--------------------------     
c  This code applies the Helmholtz volume layer potential
c     \int_{B} exp(ik|x-y|)/|x-y| f(y) where
c  f is defined as samples at tensor product legendre nodes 
c  on a level restricted tree and the potential is returned
c  as polynomial expansion and samples on the same grid.
c 
c  Input arugments
c    - eps: double precision
c       tolerance requested
c    - zk: double complex
c        Helmholtz wave number
c    - nboxes: integer
c        number of boxes
c    - nlevels: integer
c        number of levels
c    - ltree: integer
c        length of array containing the tree structure
c    - itree: integer(ltree)
c        array containing the tree structure
c    - iptr: integer(8)
c        pointer to various parts of the tree structure
c        * iptr(1) - laddr
c        * iptr(2) - ilevel
c        * iptr(3) - iparent
c        * iptr(4) - nchild
c        * iptr(5) - ichild
c        * iptr(6) - ncoll
c        * iptr(7) - coll
c        * iptr(8) - ltree
c    - norder: integer
c        order of expansions for input coefficients array
c    - ncbox: integer
c        number of coefficients of expansions of functions
c        in each of the boxes
c    - ttype: character *1
c        type of coefs provided, total order ('t') or full order('f')
c    - fcoefs: double complex (ncbox,nboxes)
c        function coeffs tabulated on a grid
c    - centers: double precision (3,nboxes)
c        xyz coordintes of boxes in the tree structure
c    - boxsize: double precision (0:nlevels)
c        size of boxes at each of the levels
c    - npbox: integer
c        number of points per box where potential is to 
c        be computed = (norder**3)
c
c  Output arguments:
c    - pot: double complex (npbox,nboxes)
c        volume potential on the tree structure (note that 
c        the potential is non-zero only in the leaf boxes
c    - potcoefs: double complex (ncbox,nboxes)
c        volume potential coefficients on the tree structure (note that 
c        the potential is non-zero only in the leaf boxes
c
      implicit real *8 (a-h,o-z)
      real *8 eps
      complex *16 zk
      integer nboxes,nlevels,ltree,itree(ltree),iptr(8)
      integer norder,ncbox
      character *1 ttype
      complex *16 fcoefs(ncbox,nboxes)
      real *8 centers(3,nboxes),boxsize(0:nlevels)
      integer npbox
      complex *16 potcoefs(ncbox,nboxes),pot(npbox,nboxes)
      real *8 timeinfo(6),tprecomp(3)
      

      call helmholtz_volume_fmm(eps,zk,nboxes,nlevels,ltree,itree,
     1  iptr,norder,ncbox,ttype,fcoefs,centers,boxsize,npbox,pot,
     2  potcoefs,timeinfo,tprecomp)

      return
      end
c
c
c
c
