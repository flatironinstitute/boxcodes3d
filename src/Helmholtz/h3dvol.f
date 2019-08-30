c
c
c     this file contains subroutines relevant to volume
c     FMM calculations for the Helmholtz kernel in 3D
c
c     DEPENDS: legetens.f
c
c     h3ddensmpmat - form the matrix taking coefficients
c                    to a multipole expansion
c     h3dtaevalgridmatp - form the matrix mapping a local
c                         expansion on a box to values of the
c                         corresponding potential at a tensor
c                         Legendre grid
c
c     TODO: faster h3dtaevalgridmatp?


      subroutine h3ddensmpmat(zk,rscale,nterms,bs,type,n,nq,wlege,nlege,
     1     mat,ldmat)
c
c     Form the matrix mapping coefficients on a box to a
c     multipole expansion
c
c----------------------------------------------------------------------
c     INPUT:
c
c     zk          : Helmholtz parameter
c     rscale      : scaling parameter of multipole expansions
c     nterms      : number of terms in multipole expansions
c     bs          : box size
c     type        : character, tensor Legendre polynomial type
c                   't' gives total degree, 'f' gives full
c     n           : use polynomials up to degree n-1
c     nq          : generate table using nq x nq x nq grid of
c                   Legendre nodes
c     wlege       : precomputed array, see ylgndrfwini
c     nlege       : parameter of wlege, see ylgndrfwini
c
c     OUTPUT:
c
c     mat         : matrix mapping coefficients of a density to
c                   a multipole expansion of length nterms for the
c                   given zk, box size, etc
c
      implicit none
      integer nterms, n, nq, nlege, ldmat
      complex *16 zk,mat(ldmat,*)
      real *8 rscale, wlege(*), bs
      character type
c     local
      integer ndeg, npol, nq3, nmp, nstemp, nd, i, j, ldu, itype
      real *8, allocatable :: rs(:,:), ws(:), pols(:)
      real *8 center(3), utemp, vtemp
      complex *16, allocatable :: mpmat(:,:), fvalmatt(:,:)
      complex *16 zero, one, ctemp, bsh, bsh3
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /
      data center / 3*0.0d0 /

      ndeg = n-1
      nq3 = nq**3
      nmp = (nterms+1)*(2*nterms+1)

      bsh = bs/2.0d0
      bsh3 = bsh**3

      call legetens_npol_3d(ndeg,type,npol)

      allocate(rs(3,nq3),ws(nq3))
      allocate(mpmat(nmp,nq3),fvalmatt(npol,nq3),pols(npol))

      itype = 1
      ldu = 1
      call legetens_exps_3d(itype,nq,type,rs,utemp,ldu,
     1     vtemp,ldu,ws)

      do i = 1,nq3
         ws(i) = ws(i)*bsh3
         rs(1,i) = rs(1,i)*bsh
         rs(2,i) = rs(2,i)*bsh
         rs(3,i) = rs(3,i)*bsh
      enddo
      
      nstemp = 1
      nd = 1
      do i = 1,nq3
         do j = 1,nmp
            mpmat(j,i) = zero
         enddo
         ctemp = ws(i)
         call h3dformmpc(nd,zk,rscale,rs(1,i),ctemp,
     1        nstemp,center,nterms,mpmat(1,i),wlege,nlege)
      enddo
      
      do i = 1,nq3
         call legetens_pols_3d(rs(1,i),ndeg,type,pols)
         do j = 1,npol
            fvalmatt(j,i) = pols(j)
         enddo
      enddo

      call zgemm('N','T',nmp,npol,nq3,one,mpmat,nmp,
     1     fvalmatt,npol,zero,mat,ldmat)

      
      return
      end
      

      subroutine h3dtaevalgridmatp(zk,rscale,nterms,bs,
     1     n,wlege,nlege,mat,ldmat)
c
c     Form the matrix mapping a local expansion on a box to values
c     of the corresponding potential at a tensor Legendre grid
c
c----------------------------------------------------------------------
c     INPUT:
c
c     zk          : Helmholtz parameter
c     rscale      : scaling parameter of multipole expansions
c     nterms      : number of terms in multipole expansions
c     bs          : box size
c     n           : eval at n x n x n tensor Legendre grid
c     wlege       : precomputed array, see ylgndrfwini
c     nlege       : parameter of wlege, see ylgndrfwini
c     ldmat       : leading dimension of mat
c
c     OUTPUT:
c
c     mat         : matrix mapping local expansion coefficients
c     in a box of specified size to values of the
c                   corresponding potential at tensor
c                   Legendre points of order n
      implicit none
      complex *16 zk, mat(ldmat,*)
      integer nterms, n, nlege, ldmat
      real *8 wlege(*), rscale, bs
c     local
      real *8 x(3,n**3), w, u, v, bsh, center(3)
      integer ldu, itype, i, n3, nd, i, j
      character type
      complex *16, allocatable :: pot(:,:), locexp(:,:)
      complex *16 zero, one
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /
      data center / 3*0.0d0 /

      n3 = n**3
      bsh = bs/2.0d0

c     grab points and scale them

      itype = 0
      ldu = 1
      type = 't'
      call legetens_exps_3d(itype,n,type,x,u,ldu,
     1     v,ldu,w)

      do i = 1,n3
         x(1,i) = x(1,i)*bsh
         x(2,i) = x(2,i)*bsh
         x(3,i) = x(3,i)*bsh
      enddo

c     use vectorized routine to fill in (do by hand later if too slow)

      nd = (nterms+1)*(2*nterms+1)
      allocate(pot(nd,n3),locexp(nd,nd))

      do i = 1,nd
         do j = i,nd
            locexp(j,i) = zero
            if (j .eq. i) locexp(j,i) = one
         enddo
      enddo

      call h3dtaevalp(nd,zk,rscale,center,locexp,nterms,x,
     1            n3,pot,wlege,nlege)
      
      do i = 1,nd
         do j = 1,n3
            mat(j,i) = pot(i,j)
         enddo
      enddo

      return
      end

