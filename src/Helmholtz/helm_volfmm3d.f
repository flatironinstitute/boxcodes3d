      subroutine helmholtz_volume_fmm(eps,zk,nboxes,nlevels,ltree,
     1   itree,iptr,norder,ncbox,type,fcoefs,centers,boxsize,npbox,
     2   xgrid,pot)
c
c       This code applies the Helmholtz volume layer potential
c       to a collection of right hand sides
c 
c       input
c         eps - double precision
c            tolerance requested
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
c         fcoefs - double complex (ncbox,nboxes)
c           tensor product legendre expansions of the right hand side
c         centers - double precision (3,nboxes)
c           xyz coordintes of boxes in the tree structure
c         boxsize - double precision (0:nlevels)
c           size of boxes at each of the levels
c         npbox - integer
c           number of points per box where potential is to be dumped
c         xgrid - double precision (3,npbox)
c           location of the targets on the standard box [-0.5,0.5]^3
c
c     output:
c         pot - double complex (npbox,nboxes)
c            volume potential on the tree structure (note that 
c            the potential is non-zero only in the leaf boxes
c

      implicit real *8 (a-h,o-z)
      real *8 eps
      complex *16 zk
      integer nboxes,nlevels,ltree
      integer itree(ltree),iptr(8),ncbox,npbox
      complex *16 fcoefs(ncbox,nboxes)
      real *8 xgrid(3,npbox)
      complex *16 pot(npbox,nboxes)

      double precision, allocatable :: scales(:)
      intger, allocatable :: nterms(:)
      double precision, allocatable :: rmlexp(:)
      integer *8, allocatable :: iaddr(:,:)
      integer lmptemp
      integer *8 lmptot
      double precision, allocatable :: mptemp(:),mptemp2(:)

      double precision, allocatable :: wlege(:)
      character *1 type

c
cc      pw stuff
c
      integer nuall,ndall,nnall,nsall,neall,nwall
      integer nu1234,nd5678,nn1256,ns3478,ne1357,nw2468
      integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
      integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8

      integer uall(200),dall(200),nall(120),sall(120),eall(72),wall(72)
      integer u1234(36),d5678(36),n1256(24),s3478(24)
      integer e1357(16),w2468(16),n12(20),n56(20),s34(20),s78(20)
      integer e13(20),e57(20),w24(20),w68(20)
      integer e1(20),e3(5),e5(5),e7(5),w2(5),w4(5),w6(5),w8(5)

      integer ntmax, nexpmax, nlams, nmax, nthmax, nphmax
      double precision, allocatable :: carray(:,:), dc(:,:)
      double precision, allocatable :: rdplus(:,:,:)
      double precision, allocatable :: rdminus(:,:,:), rdsq3(:,:,:)
      double precision, allocatable :: rdmsq3(:,:,:)
      double complex, allocatable :: rdminus2(:,:,:),zeyep(:)
      double complex, allocatable :: rdplus2(:,:,:)
      double precision, allocatable :: zmone(:)
      integer nn,nnn
  
      double complex, allocatable :: rlams(:),whts(:)

      double complex, allocatable :: rlsc(:,:,:)
      integer, allocatable :: nfourier(:), nphysical(:)
      integer nexptot, nexptotp
      double complex, allocatable :: xshift(:,:),yshift(:,:),zshift(:,:)

      double complex, allocatable :: fexp(:),fexpback(:)

      double complex, allocatable :: mexp(:,:,:,:)
      double complex, allocatable :: tmp(:,:,:),tmp2(:,:,:)
      double complex, allocatable :: mexpf1(:,:),mexpf2(:,:)
      double complex, allocatable :: mexpp1(:,:),mexpp2(:,:),
     1    mexppall(:,:,:)

      double precision, allocatable :: rsc(:)i
      integer, allocatable :: ilevrel(:)
      complex *16, allocatable :: mpcoeffsmat(:,:)
      complex *16 ac,bc



      ifprint = 1

c
c       initialize potential
c 
      do i=1,nboxes
        do j=1,npbox
          pot(j,i) = 0 
        enddo
      enddo
c
c       find scales and number of terms required at each of
c       the levels
c

      allocate(scales(0:nlevels),nterms(0:nlevels))
 
      nmax = 0
      do ilev = 0,nlevels
        scales(ilev) = boxsize(ilev)*abs(zk)i
        call h3dterms(boxsize(ilev),zk,eps,nterms(ilev))
        if(nterms(ilev).gt.nmax) nmax = nterms(ilev)
      enddo

c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(2,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c     imptemp is pointer for single expansion (dimensioned by nmax)
c
c       ... allocate iaddr and temporary arrays
c

      allocate(iaddr(2,nboxes))
      lmptemp = (nmax+1)*(2*nmax+1)*2
      allocate(mptemp(lmptemp),mptemp2(lmptemp))

      nd = 1

      call mpalloc(nd,itree(iptr(1)),iaddr,nlevels,lmptot,nterms)
      if(ifprint.ge. 1) print *, "lmptot =",lmptot/1.0d9

      iert = 0

      allocate(rmlexp(lmptot),stat=iert)
      if(iert.ne.0) then
         print *, "Cannot allocate mpole expansion workspace"
         print *, "lmptot=", lmptot
         stop
      endif
c
c       initialize wlege
c
      nlege = 100
      lw7 = (nlege+1)**2*4
      allocate(wlege(lw7))
      call ylgndrfwini(nlege,wlege,lw7,lused7)
c
c       ... set all multipole and local expansions to zero
c
      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
            call mpzero(nd,rmlexp(iaddr(1,ibox)),nterms(ilev))
            call mpzero(nd,rmlexp(iaddr(2,ibox)),nterms(ilev))
         enddo
C$OMP END PARALLEL DO          
       enddo


      

c
c
c        step 1: convert coeffs to multipole expansions
c
    
      if(ifprint.ge.1) 
     $   call prinf("=== STEP 1 (coefs -> mp) ====*",i,0) 

      allocate(ilevrel(0:nlevels))
      ilevrel(0) = 0
      ilevrel(1) = 0
      do ilev = 2,nlevels
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.gt.0) ilevrel(ilev) = 1
        enddo
      enddo
     
      ac = 1.0d0
      bc = 1.0d0
      do ilev=2,nlevels
        nmp  = (nterms(ilev)+1)**2
        if(ilevrel(ilev).eq.1) then
          nq = 10
          allocate(mpcoeffsmat(nmp,ncbox))
          call h3ddensmpmat(zk,scales(ilev),nterms(ilev),
     1     boxsize(ilev),type,norder,nq,wlege,nlege,mpcoefsmat,
     2     nmp)
        endif
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4)+ibox-1)
          if(nchild.gt.0) then
            call zgemv('n',nmp,ncbox,ac,mpcoefsmat,nmp,fcoefs(1,ibox),
     1        1,bc,rmlexp(iaddr(1,ibox)),1)  
          endif
        enddo
      enddo

      return
      end
