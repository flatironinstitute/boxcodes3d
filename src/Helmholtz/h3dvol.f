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
      real *8, allocatable :: rs(:,:), ws(:), pols(:),rs2(:,:)
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


      allocate(rs(3,nq3),ws(nq3),rs2(3,nq3))
      allocate(mpmat(nmp,nq3),fvalmatt(npol,nq3),pols(npol))

      itype = 1
      ldu = 1
      call legetens_exps_3d(itype,nq,type,rs,utemp,ldu,
     1     vtemp,ldu,ws)

      do i = 1,nq3
         ws(i) = ws(i)*bsh3
         rs2(1,i) = rs(1,i)*bsh
         rs2(2,i) = rs(2,i)*bsh
         rs2(3,i) = rs(3,i)*bsh
      enddo
      
      nstemp = 1
      nd = 1
      do i = 1,nq3
         do j = 1,nmp
            mpmat(j,i) = zero
         enddo
         ctemp = ws(i)
         call h3dformmpc(nd,zk,rscale,rs2(1,i),ctemp,
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
      real *8  w, u, v, bsh, center(3)
      real *8, allocatable :: x(:,:)
      integer ldu, itype, i, n3, nd, j,l,ii
      character type
      complex *16, allocatable :: pot(:,:), locexp(:,:,:)
      complex *16 zero, one
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /


      n3 = n**3
      allocate(x(3,n3))
      bsh = bs/2.0d0

      center(1) = 0
      center(2) = 0
      center(3) = 0



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
      allocate(pot(nd,n3),locexp(nd,0:nterms,-nterms:nterms))

      do i = 1,nd
         do j = 0,nterms
           do l=-nterms,nterms
             locexp(i,j,l) = 0
           enddo
         enddo
      enddo

      ii = 0
      do j=-nterms,nterms
        do l=0,nterms
          ii = ii+1
          if(abs(j).le.l) then
            locexp(ii,l,j) = one
          endif
        enddo
      enddo


      do i=1,n3
        do j=1,nd
          pot(j,i) = 0
        enddo
      enddo

      call h3dtaevalp(nd,zk,rscale,center,locexp,nterms,x,
     1            n3,pot,wlege,nlege)
      
      do i = 1,nd
         do j = 1,n3
            mat(j,i) = pot(i,j)
         enddo
      enddo

cc      call prin2('mat=*',mat,2*nd*n3)

      return
      end

c
c
c
c
c
c----------------------------------------------------------------------
      subroutine h3dlist4pw_vol(ilev,nd,nexptotp,nexptot,nterms,nn,
     1           nlams,nlevels,ilevlist4,itree,
     2           nfourier,nphysical,ncbox,nmax,
     3           rdminus,rdplus,rlsc,
     4           xshift,yshift,zshift,
     5           fexp,
     6           mexpf1,mexpf2,tmp,tmp2,rsc,pgboxwexp,
     7           cntlist4,fcoefs,fimat,mpcoefsmatall)
      implicit none
cccccc input/output variables
      integer ilev
      integer nd
      integer nexptotp,nexptot,ncbox,nmax
      integer nterms,nn,nlams,nlevels,cntlist4
      integer ilevlist4(*),itree(*)
      integer nfourier(*)
      integer nphysical(*)
      double precision rdminus(0:nn,0:nn,-nn:nn)
      double precision rdplus(0:nn,0:nn,-nn:nn)
      double precision rsc(*)
      double precision fimat(ncbox,*)
      double complex rlsc(0:nterms,0:nterms,nlams)
      double complex xshift(-5:5,nexptotp)
      double complex yshift(-5:5,nexptotp)
      double complex zshift(5,nexptotp)
      double complex fexp(*)
      double complex mexpf1(nd,nexptot)
      double complex mexpf2(nd,nexptot)
      double complex tmp(nd,0:nterms,-nterms:nterms)
      double complex tmp2(nd,0:nterms,-nterms:nterms)
      double complex pgboxwexp(nd,nexptotp,cntlist4,6)
      double complex mpcoefsmatall((nmax+1)*(2*nmax+1),ncbox,2:nlevels)
      double complex fcoefs(ncbox,*)
cccccc scoped function variables
      integer ibox,jbox,i,idim,nlist3,nmp,ncc,fdim
      integer istart,iend,npts,jstart,jend,npts0
      double precision time1,time2,omp_get_wtime
      double precision alpha,beta
      double complex, allocatable :: gboxfcoefs(:,:)
      double complex, allocatable :: gboxmexp(:,:)
      double complex, allocatable ::  gboxwexp(:,:,:,:)
      double complex ac,bc

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      pgboxwexp=0d0
      nmp=(nterms+1)*(2*nterms+1)
      ncc=ncbox*8
      fdim=2
      ac=1.0d0
      bc=0.0d0
      alpha=1
      beta=0
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,jbox,jstart,jend,npts,npts0,i)
C$OMP$PRIVATE(gboxwexp,gboxmexp,gboxfcoefs)
C$OMP$PRIVATE(mexpf1,mexpf2,tmp,tmp2)
      do ibox=itree(2*ilev+1),itree(2*ilev+2)
        if(ilevlist4(ibox).gt.0) then
          allocate(gboxfcoefs(ncbox,8))
          allocate(gboxmexp(nd*nmp,8))
          allocate(gboxwexp(nd,nexptotp,6,8))
          call dgemm('n','n',fdim,ncc,ncbox,alpha,fcoefs(1,ibox),fdim,
     1         fimat,ncbox,beta,gboxfcoefs(1,1),fdim)
cccccc bad code, note gboxmexp is an array not scalar
          gboxmexp=0
          jbox=ilevlist4(ibox)
          do i=1,8
            call zgemv('n',nmp,ncbox,ac,mpcoefsmatall(1,1,ilev+1),
     1           nmp,gboxfcoefs(1,i),1,bc,gboxmexp(1,i),1)
ccc    convert to plane wave
            call mpscale(nd,nterms,gboxmexp(1,i),
     1           rsc,tmp)
c
cc              process up down for current box
c
            call hmpoletoexp(nd,tmp,nterms,nlams,nfourier,
     1           nexptot,mexpf1,mexpf2,rlsc)

            call hftophys(nd,mexpf1,nlams,nfourier,nphysical,
     1           gboxwexp(1,1,1,i),fexp)

            call hftophys(nd,mexpf2,nlams,nfourier,nphysical,
     1           gboxwexp(1,1,2,i),fexp)

            call hprocessgboxudexp(nd,gboxwexp(1,1,1,i),
     1           gboxwexp(1,1,2,i),i,nexptotp,
     2           pgboxwexp(1,1,jbox,1),
     3           pgboxwexp(1,1,jbox,2),
     4           xshift,yshift,zshift)
c
cc              process north-south for current box
c
            call rotztoy(nd,nterms,tmp,tmp2,rdminus)
            call hmpoletoexp(nd,tmp2,nterms,nlams,nfourier,
     1           nexptot,mexpf1,mexpf2,rlsc)

            call hftophys(nd,mexpf1,nlams,nfourier,nphysical,
     1           gboxwexp(1,1,3,i),fexp)

            call hftophys(nd,mexpf2,nlams,nfourier,nphysical,
     1           gboxwexp(1,1,4,i),fexp)

            call hprocessgboxnsexp(nd,gboxwexp(1,1,3,i),
     1           gboxwexp(1,1,4,i),i,nexptotp,
     2           pgboxwexp(1,1,jbox,3),
     3           pgboxwexp(1,1,jbox,4),
     4           xshift,yshift,zshift)

c
cc              process east-west for current box

            call rotztox(nd,nterms,tmp,tmp2,rdplus)
            call hmpoletoexp(nd,tmp2,nterms,nlams,nfourier,
     1           nexptot,mexpf1,mexpf2,rlsc)

            call hftophys(nd,mexpf1,nlams,nfourier,nphysical,
     1           gboxwexp(1,1,5,i),fexp)

            call hftophys(nd,mexpf2,nlams,nfourier,nphysical,
     1           gboxwexp(1,1,6,i),fexp)

            call hprocessgboxewexp(nd,gboxwexp(1,1,5,i),
     1           gboxwexp(1,1,6,i),i,nexptotp,
     2           pgboxwexp(1,1,jbox,5),
     3           pgboxwexp(1,1,jbox,6),
     4           xshift,yshift,zshift)
          enddo
          deallocate(gboxfcoefs)
          deallocate(gboxmexp)
          deallocate(gboxwexp)
        endif
      enddo
C$OMP END PARALLEL DO
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      return
      end

c
c
c
c
c
c----------------------------------------------------------------------
      subroutine h3dsplitboxp_vol(n,n3,iboxsrcind,iboxfl,
     1           iboxsubcenters,iboxsrc)
      implicit none
      integer n,n3
      real *8  w, u, v, center(3)
      real *8 x(3,n3)
      integer ldu, itype
      character type
      integer iboxsrcind(n3)
      integer iboxfl(2,8)
      double precision iboxsubcenters(3,8)
      double precision iboxsrc(3,n3)
      double precision bs

      center(1) = 0
      center(2) = 0
      center(3) = 0
      bs = 1.0

      itype = 0
      ldu = 1
      type = 't'
      call legetens_exps_3d(itype,n,type,x,u,ldu,
     1     v,ldu,w)

      call subdividebox(x,n3,center,bs,iboxsrcind,iboxfl,iboxsubcenters)
      call dreorderf(3,n3,x,iboxsrc,iboxsrcind)

      return
      end

c
c
c
c
c---------------------------------------------------------------------
      subroutine scale_points(xin,xout,npts,rs)
      implicit none
      integer npts,i
      double precision rs
      double precision xin(3,npts)
      double precision xout(3,npts)

      do i = 1,npts
         xout(1,i) = xin(1,i)*rs
         xout(2,i) = xin(2,i)*rs
         xout(3,i) = xin(3,i)*rs
      enddo

      return
      end
