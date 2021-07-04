      subroutine helmholtz_volume_fmm(eps,zk,nboxes,nlevels,ltree,
     1   itree,iptr,norder,ncbox,ttype,fcoefs,centers,boxsize,npbox,
     2   pot,potcoefs,timeinfo,tprecomp)

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
c         ttype - character *1
c            type of coefs provided, total order ('t') or full order('f')
c         fcoefs - double complex (ncbox,nboxes)
c           legendre expansion of function tabulated on a grid
c         centers - double precision (3,nboxes)
c           xyz coordintes of boxes in the tree structure
c         boxsize - double precision (0:nlevels)
c           size of boxes at each of the levels
c         npbox - integer
c           number of points per box where potential is to be dumped = (norder**3)
c
c     output:
c         pot - double complex (npbox,nboxes)
c            volume potential on the tree structure (note that 
c            the potential is non-zero only in the leaf boxes
c
c         pot_coefs - double complex (ncbox,nboxes)
c            volume potential coefficients on the tree structure 
c            (note that 
c            the potential is non-zero only in the leaf boxes
c
      implicit real *8 (a-h,o-z)
      real *8 eps
      complex * 16 zk
      integer nboxes,nlevels,ltree
      integer itree(ltree),iptr(8),norder,ncbox,npbox
      character *1 ttype
      complex *16 fcoefs(ncbox,nboxes)
      complex *16 pot(npbox,nboxes)
      complex *16 potcoefs(ncbox,nboxes)
      real *8 centers(3,nboxes)
      real *8 boxsize(0:nlevels)
      real *8 timeinfo(6),tprecomp(3)

      integer impcoefsmat(0:nlevels+1),itamat(0:nlevels+1)
      integer itab(0:nlevels+1)
      integer lmpcoefsmat,ltamat,ltab
      complex *16, allocatable :: mpcoefsmat(:),tamat(:),tab(:)


      call helmholtz_volume_fmm_mem(eps,zk,nboxes,nlevels,ltree,
     1   itree,iptr,boxsize,centers,npbox,ncbox,impcoefsmat,lmpcoefsmat,
     2   itamat,ltamat,itab,ltab)
      call prinf('impcoefsmat=*',impcoefsmat,nlevels+2)
      call prinf('itamat=*',itamat,nlevels+2)
      call prinf('itab=*',itab,nlevels+2)
      print *, lmpcoefsmat,ltamat,ltab

      
      allocate(mpcoefsmat(lmpcoefsmat),tamat(ltamat),tab(ltab))

      call helmholtz_volume_fmm_init(eps,zk,nboxes,nlevels,
     1   boxsize,norder,npbox,ncbox,impcoefsmat,lmpcoefsmat,
     2   itamat,ltamat,itab,ltab,ttype,mpcoefsmat,tamat,tab,tprecomp)

      call helmholtz_volume_fmm_wprecomp(eps,zk,nboxes,nlevels,
     1   ltree,itree,iptr,norder,ncbox,ttype,fcoefs,centers,boxsize,
     2   mpcoefsmat,impcoefsmat,lmpcoefsmat,tamat,itamat,
     3   ltamat,tab,itab,ltab,npbox,pot,potcoefs,timeinfo)
      

      return
      end
c
c
c
c
c
      subroutine helmholtz_volume_fmm_mem(eps,zk,nboxes,nlevels,ltree,
     1   itree,iptr,boxsize,centers,npbox,ncbox,impcoefsmat,
     2   lmpcoefsmat,itamat,ltamat,itab,ltab)
c    This subroutine estimates the memory requirements for the volume
c    fmm precmoputation
c
c    input arguments:
c         eps - real *8
c           precision
c         zk - complex *16
c           wave number for the problem
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
c      boxsize - real *8(0:nlevels)
c         size of box at each level
c
c      npbox = integer
c        number of points per box
c 
c      ncbox - integrer
c        number of polynomials per box
c
c   output arguemnts:
c      impcoefsmat - integer(0:nlevels+1)
c        array estimating the size of mpcoefsmat at each level
c      lmpcoefsmat - integer
c        size of mpcoefsmat
c      itamat - integer(0:nlevels+1)
c        array estimating the size of tamat at each level
c      ltamat - integer
c        size of tamat
c      itab - integer(0:nlevels+1)
c        array estimating the size of local quadratue table
c         at each level
c      ltab - integer
c        length of local correction tables
c   

      implicit real *8 (a-h,o-z)
      integer nboxes,nlevels,ltree
      integer itree(ltree),iptr(8),impcoefsmat(0:nlevels+1)
      integer itamat(0:nlevels+1),itab(0:nlevels+1)
      integer lmpcoefsmat,ltamat,ltab
      integer npbox,ncbox
      complex *16 zk
      real *8 boxsize(0:nlevels),eps,centers(3,nboxes)

      integer nterms(0:nlevels),ilevrel(0:nlevels)
      real *8 rscales(0:nlevels)

      integer, allocatable :: nlist1(:),list1(:,:)
      integer, allocatable :: nlist2(:),list2(:,:)
      integer, allocatable :: nlist3(:),list3(:,:)
      integer, allocatable :: nlist4(:),list4(:,:)
c
c
c
c      get number of terms
c

      nmax = 0
      do ilev = 0,nlevels
        rscales(ilev) = boxsize(ilev)*abs(zk)

        if(rscales(ilev).gt.1) rscales(ilev) = 1.0d0
        call h3dterms(boxsize(ilev),zk,eps,nterms(ilev))
        if(nterms(ilev).gt.nmax) nmax = nterms(ilev)
      enddo


c
c   generate mpcoefsmat
c
c
      lmpcoefsmat = 0
      impcoefsmat(0) = 1
      do ilev=0,nlevels
        ilevrel(ilev) = 0
        do ibox=itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4)+ibox-1)
          if(nchild.eq.0) ilevrel(ilev) = 1
        enddo
        llev = 0
        if(ilevrel(ilev).eq.1) then
           nmp = (nterms(ilev)+1)*(2*nterms(ilev)+1)
           llev = nmp*ncbox
        endif
        impcoefsmat(ilev+1) = impcoefsmat(ilev)+llev
        lmpcoefsmat = lmpcoefsmat + llev
      enddo

      print *, "lmpcoefsmat=",lmpcoefsmat


c
c
c      generate tamat
c
      itamat(0) = 1
      itamat(1) = 1
      itamat(2) = 1
      ltamat = 0
      do ilev=2,nlevels
         llev = 0
         if(ilevrel(ilev).eq.1) then
           nmp = (nterms(ilev)+1)*(2*nterms(ilev)+1)
           llev = nmp*npbox
         endif
         itamat(ilev+1) = itamat(ilev) + llev
         ltamat = ltamat + llev
      enddo
c
c  end of getting tamat
c
c
c   compute list info
c
c
      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0

      isep = 1
      mnbors = 27
      iper = 0
      call computemnlists(nlevels,nboxes,itree(iptr(1)),
     1       boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2       itree(iptr(5)),isep,itree(iptr(6)),mnbors,itree(iptr(7)),
     3       iper,mnlist1,mnlist2,mnlist3,mnlist4)

      allocate(list1(mnlist1,nboxes),list2(mnlist2,nboxes))
      allocate(list3(mnlist3,nboxes),list4(mnlist4,nboxes))
      allocate(nlist1(nboxes))
      allocate(nlist2(nboxes))
      allocate(nlist3(nboxes))
      allocate(nlist4(nboxes))

      call computelists(nlevels,nboxes,itree(iptr(1)),
     1   boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2   itree(iptr(5)),isep,itree(iptr(6)),mnbors,itree(iptr(7)),
     3   iper,nlist1,mnlist1,list1,nlist2,mnlist2,list2,
     4   nlist3,mnlist3,list3,nlist4,mnlist4,list4)


c
c   get near tables
c
c

      itab(0) = 1
      ltab = 0
      llev0 = 10*ncbox*ncbox
      do ilev=0,nlevels
         ilevrel(ilev) = 0
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
            if(nlist1(ibox).gt.0) ilevrel(ilev) = 1
         enddo

         llev = 0
         if(ilevrel(ilev).eq.1) llev = llev0
         
         itab(ilev+1) = itab(ilev) + llev
         ltab = ltab + llev
      enddo

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
c
c
      subroutine helmholtz_volume_fmm_init(eps,zk,nboxes,nlevels,
     1   boxsize,norder,npbox,ncbox,impcoefsmat,lmpcoefsmat,
     2   itamat,ltamat,itab,ltab,ttype,mpcoefsmat,tamat,tab,tprecomp)
c    This subroutine estimates the memory requirements for the volume
c    fmm precmoputation
c
c    input arguments:
c         eps - real *8
c           precision
c         zk - complex *16
c           wave number for the problem
c         nboxes - integer
c            number of boxes
c         nlevels - integer
c            number of levels
c      boxsize - real *8(0:nlevels)
c         size of box at each level
c      norder - integer
c        order of discretization
c      npbox = integer
c        number of points per box
c      ncbox - integrer
c        number of polynomials per box
c      impcoefsmat - integer(0:nlevels+1)
c        array estimating the size of mpcoefsmat at each level
c      lmpcoefsmat - integer
c        size of mpcoefsmat
c      itamat - integer(0:nlevels+1)
c        array estimating the size of tamat at each level
c      ltamat - integer
c        size of tamat
c      itab - integer(0:nlevels+1)
c        array estimating the size of local quadratue table
c         at each level
c      ltab - integer
c        length of local correction tables
c      ttype - character
c        type of expansions
c
c   output arguments:
c     mpcoefsmat - complex *16 (lmpcoefsmat)
c        coefs -> mp mat at relevant levels
c     tamat - complex *16 (ltamat)
c        ta -> pot mat at relevant levels
c     tab - complex *16(ltab)
c        near quadrature correction at each level
c     tprecomp - real *8 (3)
c        time taken in various precomputation steps
c   

      implicit real *8 (a-h,o-z)
      integer nboxes,nlevels,ltree
      integer impcoefsmat(0:nlevels+1)
      integer itamat(0:nlevels+1),itab(0:nlevels+1)
      integer lmpcoefsmat,ltamat,ltab
      integer npbox,ncbox
      complex *16 zk,zk2
      real *8 boxsize(0:nlevels),eps

      complex *16 mpcoefsmat(lmpcoefsmat),tab(ltab)
      complex *16 tamat(ltamat)

      character *1 ttype

      integer nterms(0:nlevels),ilevrel(0:nlevels)
      real *8 rscales(0:nlevels),tprecomp(3)

      real *8, allocatable :: x(:,:),pmat(:,:),pmat_qr(:,:)
      real *8, allocatable :: pmat_tau(:),pols(:)
      integer, allocatable :: pmat_jpvt(:),ipt2depth(:,:)
      integer, allocatable :: idepth2pt(:,:,:)

      real *8, allocatable :: wlege(:)

c
c
c
c      get number of terms
c

      nmax = 0
      do ilev = 0,nlevels
        rscales(ilev) = boxsize(ilev)*abs(zk)

        if(rscales(ilev).gt.1) rscales(ilev) = 1.0d0
        call h3dterms(boxsize(ilev),zk,eps,nterms(ilev))
        if(nterms(ilev).gt.nmax) nmax = nterms(ilev)
      enddo

c
c       initialize wlege
c
      nlege = nmax+10
      lw7 = (nlege+1)**2*4
      allocate(wlege(lw7))
      call ylgndrfwini(nlege,wlege,lw7,lused7)

      call prinf('Finished initializing wlege*',i,0)


c
c
c   generate mpcoefsmat
c
c
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      

      do ilev=0,nlevels

        nnn = impcoefsmat(ilev+1)-impcoefsmat(ilev)
        print *, ilev,nnn
        if(nnn.gt.0) then
          nq = 20
          nmp = (nterms(ilev)+1)*(2*nterms(ilev)+1)
          call h3ddensmpmat(zk,rscales(ilev),nterms(ilev),
     1      boxsize(ilev),ttype,norder,nq,wlege,nlege,
     2      mpcoefsmat(impcoefsmat(ilev)),nmp)
        endif
      enddo
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      tprecomp(1) = t2-t1
c
c   end of getting mpcoefsmat
c
      print *, "finished generating mpcoefsmat"

c
c
c      generate tamat
c
      call cpu_time(t1)
C$      t1 = omp_get_wtime()     

      do ilev=2,nlevels
        nnn = itamat(ilev+1)-itamat(ilev)
        if(nnn.gt.1) then
          nmp = (nterms(ilev)+1)*(2*nterms(ilev)+1)
          call h3dtaevalgridmatp_fast(zk,rscales(ilev),nterms(ilev),
     1      boxsize(ilev),norder,wlege,nlege,tamat(itamat(ilev)),npbox)
        endif
      enddo

      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      tprecomp(2) = t2-t1
c
c  end of getting tamat
c
c   get near tables
c
c

      call cpu_time(t1)
C$      t1 = omp_get_wtime()      


      ndeg = norder-1
      nptmax = npbox
      allocate(x(3,nptmax))
      allocate(ipt2depth(3,nptmax),idepth2pt(norder,norder,norder))

      ifloor = 1.75d0*ncbox
      call fakepolya3d(norder,ifloor,npt,x,ipt2depth,idepth2pt)
      itype = 0
      call prinf('npt=*',npt,1)
      call prin2('rat=*',(npt+0.0d0)/(nptmax+0.0d0),1)


      allocate(pmat(npt,ncbox),pmat_qr(npt,ncbox),pols(ncbox))
      allocate(pmat_tau(ncbox),pmat_jpvt(ncbox))

      do i=1,npt
        call legetens_pols_3d(x(1,i),ndeg,ttype,pols)
        do j=1,ncbox
          pmat(i,j) = pols(j)
        enddo
      enddo

      call get_qrdecomp(npt,ncbox,pmat,pmat_qr,pmat_jpvt,pmat_tau)




      ldtab = 10*ncbox
c
c  This needs fixing..
c
      nup = max(ndeg,5)
      iflg = 1
      do ilev=0,nlevels
        nnn = itab(ilev+1)-itab(ilev)
        if(nnn.gt.0) then
          zk2 = zk*boxsize(ilev)/2.0d0 
          call h3dtabp_ref(ndeg,zk2,eps,npt,x,tab(itab(ilev)),ldtab,
     1       nup,iflg,pmat_qr,ncbox,pmat_jpvt,pmat_tau,tt1,tt2,tt3)
        endif
      enddo

      call cpu_time(t2)
C$      t2 = omp_get_wtime()      

      tprecomp(3) = t2-t1
      print *, "here2"

      return
      end
c
c
c
c
c



      subroutine helmholtz_volume_fmm_wprecomp(eps,zk,nboxes,nlevels,
     1   ltree,itree,iptr,norder,ncbox,ttype,fcoefs,centers,boxsize,
     2   mpcoefsmat,impcoefsmat,lmpcoefsmat,tamat,itamat,
     3   ltamat,tab,itab,ltab,npbox,pot,potcoefs,timeinfo)

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
c         ttype - character *1
c            type of coefs provided, total order ('t') or full order('f')
c         fcoefs - double complex (ncbox,nboxes)
c            function coeffs tabulated on the tree
c         centers - double precision (3,nboxes)
c           xyz coordintes of boxes in the tree structure
c         boxsize - double precision (0:nlevels)
c           size of boxes at each of the levels
c         mpcoefsmat - double complex(lmpcoefsmat)
c           precomputed coefs -> mp matrices
c         impcoefsmat - integer(0:nlevels+1)
c           impcoefsmat(ilev) indicates location in mpcoefsmat where
c           coefs -> mp matrix for level ilev begins
c         lmpcoefsmat -> length of mpcoefsmat
c         tamat - double complex(ltamat)
c           precomputed ta -> pot matrix for various levels
c         itamat -> integer(0:nlevels+1)
c           itamat(ilev) is the location in tamat array where
c           ta -> pot matrix for level ilev begins
c         ltamat - integer
c           length of tamat array
c         tab - double complex(ltab)
c           precomputed local quadrature for all levels
c         itab - integer(0:nlevels+1)
c           itab(ilev) is the location in tab array where
c           near field tables for level ilev
c         npbox - integer
c           number of points per box where potential is to be dumped = (norder**3)
c
c     output:
c         pot - double complex (npbox,nboxes)
c            volume potential on the tree structure (note that 
c            the potential is non-zero only in the leaf boxes
c

      implicit real *8 (a-h,o-z)
      integer nd
      real *8 eps
      complex *16 zk,zk2
      integer nboxes,nlevels,ltree
      integer itree(ltree),iptr(8),ncbox,npbox,ncc
      complex *16 fcoefs(ncbox,nboxes)
      complex *16 pot(npbox,nboxes),potcoefs(ncbox,nboxes)
      double precision boxsize(0:nlevels),centers(3,nboxes)

      integer ltamat,ltab,lmpcoefsmat
      integer itamat(0:nlevels+1),itab(0:nlevels+1)
      integer impcoefsmat(0:nlevels+1)

      complex *16 tamat(ltamat),tab(ltab),mpcoefsmat(lmpcoefsmat)

      double precision, allocatable :: rscales(:)
      integer, allocatable :: nterms(:)
      double precision, allocatable :: rmlexp(:)
      integer *8, allocatable :: iaddr(:,:)
      integer lmptemp
      integer *8 lmptot
      double precision, allocatable :: mptemp(:),mptemp2(:)

      double precision, allocatable :: wlege(:)

      double precision xtargtmp(3)
      complex *16 pottmp,pottmpex,pottmp2
      character *1 ttype
      double precision, allocatable :: xnodes(:),wts(:)

      real *8, allocatable :: xref(:,:),umat(:,:),vmat(:,:),wts0(:)
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
      integer nn,nnn,ii
  
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

      double precision, allocatable :: rsc(:)
      integer, allocatable :: ilevrel(:)
      complex *16, allocatable :: tabcoll(:,:,:),tabbtos(:,:,:),
     1   tabstob(:,:,:)
      complex *16, allocatable :: tabtmp(:,:)
      complex *16, allocatable :: rhs(:,:),vals(:,:),coefs(:,:)
      complex *16 ac,bc,ima

      real *8 ra,rb

      integer nquad2,ifinit2

c
cc        temporary list info
c

      integer, allocatable :: nlist1(:),list1(:,:)
      integer, allocatable :: nlist1_detailed(:),list1_detailed(:,:)
      integer, allocatable :: nlist2(:),list2(:,:)
      integer, allocatable :: nlist3(:),list3(:,:)
      integer, allocatable :: nlist4(:),list4(:,:)

      integer, allocatable :: ijboxlist(:,:)
      double precision timeinfo(6)

      integer iref(100),idimp(3,100),iflip(3,100)
      integer irefbtos(100),idimpbtos(3,100),iflipbtos(3,100)
      integer irefstob(100),idimpstob(3,100),iflipstob(3,100)

      integer cntlist4
      integer, allocatable :: ilevlist4(:)
      double complex, allocatable :: pgboxwexp(:,:,:,:)
      double precision, allocatable :: fimat(:,:)

      double complex, allocatable :: iboxlexp(:,:)
      double precision iboxsubcenters(3,8)
      double precision subcenters(3,8)
      double precision iboxsrc(3,npbox)
      double precision subpts(3,npbox)
      double complex iboxpot(1,npbox)
      integer iboxsrcind(npbox)
      integer iboxfl(2,8)

      data ima/(0.0d0,1.0d0)/



      ifprint = 1

      allocate(ilevlist4(nboxes))
      do i=1,nboxes
        ilevlist4(i)=0
      enddo
      cntlist4 = 0

      done = 1
      pi = atan(done)*4

      ntmax = 1000
      allocate(rlams(ntmax),whts(ntmax),nfourier(ntmax),
     1   nphysical(ntmax))


      max_nodes = 10000
      allocate(xnodes(max_nodes))
      allocate(wts(max_nodes))


c
c       compute coefs
c
      
      allocate(xref(3,npbox),umat(ncbox,npbox),
     1   vmat(npbox,ncbox),wts0(npbox))
     
      itype = 2 
      call legetens_exps_3d(itype,norder,ttype,xref,umat,ncbox,vmat,
     1  npbox,wts0)
      


c
c      temporary measure for code consistency
c
      nd = 1

c
c       initialize potential
c 
      do i=1,nboxes
        do j=1,npbox
          pot(j,i) = 0 
        enddo
      enddo

c
c
c       compute list info
c
      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0

      isep = 1
      mnbors = 27
      iper = 0
      call computemnlists(nlevels,nboxes,itree(iptr(1)),
     1       boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2       itree(iptr(5)),isep,itree(iptr(6)),mnbors,itree(iptr(7)),
     3       iper,mnlist1,mnlist2,mnlist3,mnlist4)

      allocate(list1(mnlist1,nboxes),list2(mnlist2,nboxes))
      allocate(list3(mnlist3,nboxes),list4(mnlist4,nboxes))
      allocate(nlist1(nboxes))
      allocate(nlist2(nboxes))
      allocate(nlist3(nboxes))
      allocate(nlist4(nboxes))

      call computelists(nlevels,nboxes,itree(iptr(1)),
     1   boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2   itree(iptr(5)),isep,itree(iptr(6)),mnbors,itree(iptr(7)),
     3   iper,nlist1,mnlist1,list1,nlist2,mnlist2,list2,
     4   nlist3,mnlist3,list3,nlist4,mnlist4,list4)

      call prinf('mnlist4=*',mnlist4,1)


      allocate(ijboxlist(2,nboxes))

c
c       find scales and number of terms required at each of
c       the levels
c

      allocate(rscales(0:nlevels),nterms(0:nlevels))


 
      nmax = 0
      do ilev = 0,nlevels
        rscales(ilev) = boxsize(ilev)*abs(zk)

        if(rscales(ilev).gt.1) rscales(ilev) = 1.0d0
        call h3dterms(boxsize(ilev),zk,eps,nterms(ilev))
        if(nterms(ilev).gt.nmax) nmax = nterms(ilev)
      enddo
      call prinf('nterms=*',nterms,nlevels+1)
      call prin2('rscales=*',rscales,nlevels+1)

      allocate(rsc(0:nmax))


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

      call prinf('Finished initializing wlege*',i,0)

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

       call prinf('finished initializing rmlexp*',i,0)


      

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      ncc = 8*ncbox
      allocate(fimat(ncbox,ncc))
      call get_children_fcoef_interp_mat(norder,ncbox,ncc,fimat)
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      print *, "coefs interp mat time: ", time2-time1

      call h3dsplitboxp_vol(norder,npbox,iboxsrcind,iboxfl,
     1     iboxsubcenters,iboxsrc)

c
c
c        step 1: convert coeffs to multipole expansions
c
    
      if(ifprint.ge.1) 
     $   call prinf("=== STEP 1 (coefs -> mp) ====*",i,0)


      allocate(ilevrel(0:nlevels))
      ilevrel(0) = 0
      ilevrel(1) = 0
      do ilev = 0,nlevels
        ilevrel(ilev) = 0
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0) ilevrel(ilev) = 1
        enddo
      enddo

      ac = 1.0d0
      bc = 0.0d0
      do ilev=1,nlevels
        nmp  = (nterms(ilev)+1)*(2*nterms(ilev)+1)
        istart = impcoefsmat(ilev)
        print *, ilev,istart
        if(ilevrel(ilev).eq.1) then
          do ibox = itree(2*ilev+1),itree(2*ilev+2)
            nchild = itree(iptr(4)+ibox-1)
            if(nchild.eq.0) then
              call zgemv('n',nmp,ncbox,ac,mpcoefsmat(istart),nmp,
     1         fcoefs(1,ibox),1,bc,rmlexp(iaddr(1,ibox)),1) 
            endif
          enddo
        endif
      enddo

      call cpu_time(time2)
C$       time2 = omp_get_wtime()
      timeinfo(1) = time2-time1

      call cpu_time(time1)
C$      time1 = omp_get_wtime()     
      do ilev=nlevels-1,0,-1
         nquad2 = nterms(ilev)*2.5
         nquad2 = max(6,nquad2)
         ifinit2 = 1
         call legewhts(nquad2,xnodes,wts,ifinit2)
         neval = 0
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
           nchild = itree(iptr(4)+ibox-1)
           if(nchild.gt.0) then
             neval = neval + 1
           endif
         enddo

         print *, ilev,neval
         radius = boxsize(ilev)/2*sqrt(3.0d0)


C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,i,jbox,nchild) SCHEDULE(DYNAMIC)
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
            nchild = itree(iptr(4)+ibox-1)
            if(nchild.gt.0) then
               do i=1,8
                 jbox = itree(iptr(5) + 8*(ibox-1)+i-1)
                 call h3dmpmp(nd,zk,rscales(ilev+1),
     1             centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2             nterms(ilev+1),rscales(ilev),centers(1,ibox),
     3             rmlexp(iaddr(1,ibox)),nterms(ilev),
     4             radius,xnodes,wts,nquad2)
                enddo
            endif
         enddo
C$OMP END PARALLEL DO          
      enddo
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(2)=time2-time1

      call prin2('mpmp time=*',timeinfo(2),1)


      if(ifprint.ge.1)
     $    call prinf('=== Step 3 (mp to loc) + formta+mpeval ===*',i,0)
c      ... step 3, convert multipole expansions into local
c       expansions

      call cpu_time(time1)
C$        time1=omp_get_wtime()
      do ilev = 2,nlevels

c
cc       load the necessary quadrature for plane waves
c
      
         zk2 = zk*boxsize(ilev)
         if(real(zk2).le.16*pi.and.imag(zk2).le.12*pi) then
            ier = 0
c
c             get new pw quadrature
c
            
            call hwts3e(ier,eps,zk2,rlams,whts,nlams)
            call hnumfour(eps,zk2,nlams,nfourier)
            call hnumphys(eps,zk2,nlams,nphysical)


            
            nphmax = 0
            nthmax = 0
            nexptotp = 0
            nexptot = 0
            nn = 0
            do i=1,nlams
               nexptotp = nexptotp + nphysical(i)
               nexptot = nexptot + 2*nfourier(i)+1
               nn = nn + nfourier(i)*nphysical(i)
               if(nfourier(i).gt.nthmax) nthmax = nfourier(i)
               if(nphysical(i).gt.nphmax) nphmax = nphysical(i)
            enddo
            allocate(fexp(nn),fexpback(nn))

            allocate(xshift(-5:5,nexptotp))
            allocate(yshift(-5:5,nexptotp))
            allocate(zshift(5,nexptotp))
            allocate(rlsc(0:nterms(ilev),0:nterms(ilev),nlams))
            allocate(tmp(nd,0:nterms(ilev),-nterms(ilev):nterms(ilev)))
            allocate(tmp2(nd,0:nterms(ilev),-nterms(ilev):nterms(ilev)))
 
            allocate(mexpf1(nd,nexptot),mexpf2(nd,nexptot),
     1          mexpp1(nd,nexptotp))
            allocate(mexpp2(nd,nexptotp),mexppall(nd,nexptotp,16))


c
cc      NOTE: there can be some memory savings here
c
            bigint = 0
            bigint = nboxes
            bigint = bigint*6
            bigint = bigint*nexptotp*nd

            if(ifprint.ge.1) print *, "mexp memory=",bigint/1.0d9


            allocate(mexp(nd,nexptotp,nboxes,6),stat=iert)
            if(iert.ne.0) then
              print *, "Cannot allocate pw expansion workspace"
              print *, "bigint=", bigint
              stop
            endif


            nn = nterms(ilev)
            allocate(carray(4*nn+1,4*nn+1))
            allocate(dc(0:4*nn,0:4*nn))
            allocate(rdplus(0:nn,0:nn,-nn:nn))
            allocate(rdminus(0:nn,0:nn,-nn:nn))
            allocate(rdsq3(0:nn,0:nn,-nn:nn))
            allocate(rdmsq3(0:nn,0:nn,-nn:nn))

c     generate rotation matrices and carray
            call getpwrotmat(nn,carray,rdplus,rdminus,rdsq3,rdmsq3,dc)


            call hrlscini(rlsc,nlams,rlams,rscales(ilev),zk2,
     1         nterms(ilev))
            call hmkexps(rlams,nlams,nphysical,nexptotp,zk2,xshift,
     1           yshift,zshift)
            
            call hmkfexp(nlams,nfourier,nphysical,fexp,fexpback)


c
cc      zero out mexp
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(idim,i,j,k)
            do k=1,6
               do i=1,nboxes
                  do j=1,nexptotp
                     do idim=1,nd
                        mexp(idim,j,i,k) = 0.0d0
                     enddo
                  enddo
               enddo
            enddo
C$OMP END PARALLEL DO    


c
cc         compute powers of scaling parameter
c          for rescaling the multipole expansions
c
c          note: the scaling for helmholtz has been eliminated
c         since it is taken care in the scaling of the legendre
c         functions
c
          
           r1 = 1.0d0
           rsc(0) = 1.0d0
           do i=1,nterms(ilev)
             rsc(i) = rsc(i-1)*r1
           enddo

           cntlist4=0
           do ibox=itree(2*(ilev-1)+1),itree(2*(ilev-1)+2)
             if(nlist3(ibox).gt.0) then
               cntlist4=cntlist4+1
               ilevlist4(ibox)=cntlist4
             endif
           enddo
           allocate(pgboxwexp(nd,nexptotp,cntlist4,6))
           print *,"cnlist4:",cntlist4,"ilev",ilev
           call h3dlist4pw_vol(ilev-1,nd,nexptotp,nexptot,nterms(ilev),
     1          nn,nlams,nlevels,ilevlist4,itree,nfourier,nphysical,
     2          ncbox,nmax,rdminus,rdplus,rlsc,xshift,yshift,zshift,
     3          fexp,mexpf1,mexpf2,tmp,tmp2,rsc,pgboxwexp,cntlist4,
     4          fcoefs,fimat,mpcoefsmat(impcoefsmat(ilev)))

           call prinf('before starting mp to pw*',i,0)

c
cc         create multipole to plane wave expansion for
c          all boxes at this level
c
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,tmp,mexpf1,mexpf2,tmp2)
            do ibox = itree(2*ilev+1),itree(2*ilev+2)
c           rescale multipole expansion
               call mpscale(nd,nterms(ilev),rmlexp(iaddr(1,ibox)),
     1               rsc,tmp)
                
               call hmpoletoexp(nd,tmp,nterms(ilev),
     1                  nlams,nfourier,nexptot,mexpf1,mexpf2,rlsc) 

               call hftophys(nd,mexpf1,nlams,nfourier,nphysical,
     1                 mexp(1,1,ibox,1),fexp)           

               call hftophys(nd,mexpf2,nlams,nfourier,nphysical,
     1                 mexp(1,1,ibox,2),fexp)


c             form mexpnorth, mexpsouth for current box

c             Rotate mpole for computing mexpnorth and
c             mexpsouth
               call rotztoy(nd,nterms(ilev),tmp,
     1                           tmp2,rdminus)

               call hmpoletoexp(nd,tmp2,nterms(ilev),nlams,
     1                  nfourier,nexptot,mexpf1,mexpf2,rlsc)

               call hftophys(nd,mexpf1,nlams,nfourier,
     1                 nphysical,mexp(1,1,ibox,3),fexp)           

               call hftophys(nd,mexpf2,nlams,nfourier,
     1                 nphysical,mexp(1,1,ibox,4),fexp)   


c             Rotate mpole for computing mexpeast, mexpwest
               call rotztox(nd,nterms(ilev),tmp,
     1                              tmp2,rdplus)
               call hmpoletoexp(nd,tmp2,nterms(ilev),nlams,
     1                  nfourier,nexptot,mexpf1,mexpf2,rlsc)

               call hftophys(nd,mexpf1,nlams,nfourier,
     1                 nphysical,mexp(1,1,ibox,5),fexp)

               call hftophys(nd,mexpf2,nlams,nfourier,
     1                 nphysical,mexp(1,1,ibox,6),fexp)           

            enddo
C$OMP END PARALLEL DO      

           call prinf('finished converting mp to pw*',i,0)
           


c
cc         loop over parent boxes and ship plane wave
c          expansions to the first child of parent 
c          boxes. 
c          The codes are now written from a gathering perspective
c
c          so the first child of the parent is the one
c          recieving all the local expansions
c          coming from all the lists
c
c          
c

C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,nchild)
C$OMP$PRIVATE(mexpf1,mexpf2,mexpp1,mexpp2,mexppall)
C$OMP$PRIVATE(nuall,uall,ndall,dall,nnall,nall,nsall,sall)
C$OMP$PRIVATE(neall,eall,nwall,wall,nu1234,u1234,nd5678,d5678)
C$OMP$PRIVATE(nn1256,n1256,ns3478,s3478,ne1357,e1357,nw2468,w2468)
C$OMP$PRIVATE(nn12,n12,nn56,n56,ns34,s34,ns78,s78,ne13,e13,ne57,e57)
C$OMP$PRIVATE(nw24,w24,nw68,w68,ne1,e1,ne3,e3,ne5,e5,ne7,e7)
C$OMP$PRIVATE(nw2,w2,nw4,w4,nw6,w6,nw8,w8)
C$OMP$PRIVATE(jstart,jend,i)
C$OMP$PRIVATE(iboxlexp,subcenters,subpts,iboxpot) SCHEDULE(DYNAMIC)
            do ibox = itree(2*ilev-1),itree(2*ilev)
           
               nchild = itree(iptr(4)+ibox-1)
               if(nchild.gt.0) then

              
                  call getpwlistall(ibox,boxsize(ilev),nboxes,
     1            itree(iptr(6)+ibox-1),itree(iptr(7)+
     2            27*(ibox-1)),nchild,itree(iptr(5)),centers,
     3            isep,nuall,uall,ndall,dall,nnall,nall,nsall,sall,
     4            neall,eall,nwall,wall,nu1234,u1234,nd5678,d5678,
     5            nn1256,n1256,ns3478,s3478,ne1357,e1357,nw2468,w2468,
     6            nn12,n12,nn56,n56,ns34,s34,ns78,s78,ne13,e13,ne57,
     7            e57,nw24,w24,nw68,w68,ne1,e1,ne3,e3,ne5,e5,ne7,e7,
     8            nw2,w2,nw4,w4,nw6,w6,nw8,w8)
                 

                  call hprocessudexp(nd,zk2,ibox,ilev,nboxes,centers,
     1            itree(iptr(5)),rscales(ilev),boxsize(ilev),
     2            nterms(ilev),
     2            iaddr,rmlexp,rlams,whts,
     3            nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4            nuall,uall,nu1234,u1234,ndall,dall,nd5678,d5678,
     5            mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1,1),
     6            mexppall(1,1,2),mexppall(1,1,3),mexppall(1,1,4),
     7            xshift,yshift,zshift,fexpback,rlsc,pgboxwexp,
     8            cntlist4,ilevlist4,nlist4,list4,mnlist4)
                  
                  call hprocessnsexp(nd,zk2,ibox,ilev,nboxes,centers,
     1            itree(iptr(5)),rscales(ilev),boxsize(ilev),
     2            nterms(ilev),
     2            iaddr,rmlexp,rlams,whts,
     3            nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4            nnall,nall,nn1256,n1256,nn12,n12,nn56,n56,nsall,sall,
     5            ns3478,s3478,ns34,s34,ns78,s78,
     6            mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1,1),
     7            mexppall(1,1,2),mexppall(1,1,3),mexppall(1,1,4),
     8            mexppall(1,1,5),mexppall(1,1,6),mexppall(1,1,7),
     9            mexppall(1,1,8),rdplus,xshift,yshift,zshift,
     9            fexpback,rlsc,pgboxwexp,cntlist4,ilevlist4,nlist4,
     9            list4,mnlist4)


                  call hprocessewexp(nd,zk2,ibox,ilev,nboxes,centers,
     1            itree(iptr(5)),rscales(ilev),boxsize(ilev),
     2            nterms(ilev),
     2            iaddr,rmlexp,rlams,whts,
     3            nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexp,
     4            neall,eall,ne1357,e1357,ne13,e13,ne57,e57,ne1,e1,
     5            ne3,e3,ne5,e5,ne7,e7,nwall,wall,
     5            nw2468,w2468,nw24,w24,nw68,w68,
     5            nw2,w2,nw4,w4,nw6,w6,nw8,w8,
     6            mexpf1,mexpf2,mexpp1,mexpp2,mexppall(1,1,1),
     7            mexppall(1,1,2),mexppall(1,1,3),mexppall(1,1,4),
     8            mexppall(1,1,5),mexppall(1,1,6),
     8            mexppall(1,1,7),mexppall(1,1,8),mexppall(1,1,9),
     9            mexppall(1,1,10),mexppall(1,1,11),mexppall(1,1,12),
     9            mexppall(1,1,13),mexppall(1,1,14),mexppall(1,1,15),
     9            mexppall(1,1,16),rdminus,xshift,yshift,zshift,
     9            fexpback,rlsc,pgboxwexp,cntlist4,ilevlist4,nlist4,
     9            list4,mnlist4)
               endif

c
c handle mp eval for list3
c
               if(nlist3(ibox).gt.0) then
                 call getlist3pwlistall(ibox,boxsize(ilev),nboxes,
     1                nlist3(ibox),list3(1,ibox),isep,centers,
     2                nuall,uall,ndall,dall,nnall,nall,
     3                nsall,sall,neall,eall,nwall,wall)

                 allocate(iboxlexp(nd*(nterms(ilev)+1)*
     1                   (2*nterms(ilev)+1),8))
                 iboxlexp=0
                 call hprocesslist3udexplong(nd,zk2,ibox,nboxes,centers,
     1                boxsize(ilev),nterms(ilev),iboxlexp,rlams,whts,
     2                nlams,nfourier,nphysical,nthmax,nexptot,
     3                nexptotp,mexp,nuall,uall,ndall,dall,
     4                mexpf1,mexpf2,mexpp1,mexpp2,
     5                mexppall(1,1,1),mexppall(1,1,2),
     6                xshift,yshift,zshift,fexpback,rlsc)

                 call hprocesslist3nsexplong(nd,zk2,ibox,nboxes,centers,
     1                boxsize(ilev),nterms(ilev),iboxlexp,rlams,whts,
     2                nlams,nfourier,nphysical,nthmax,nexptot,
     3                nexptotp,mexp,nnall,nall,nsall,sall,
     4                mexpf1,mexpf2,mexpp1,mexpp2,
     5                mexppall(1,1,1),mexppall(1,1,2),rdplus,
     6                xshift,yshift,zshift,fexpback,rlsc)

                 call hprocesslist3ewexplong(nd,zk2,ibox,nboxes,centers,
     1                boxsize(ilev),nterms(ilev),iboxlexp,rlams,whts,
     2                nlams,nfourier,nphysical,nthmax,nexptot,
     3                nexptotp,mexp,neall,eall,nwall,wall,
     4                mexpf1,mexpf2,mexpp1,mexpp2,
     5                mexppall(1,1,1),mexppall(1,1,2),rdminus,
     6                xshift,yshift,zshift,fexpback,rlsc)

ccccccc          evaluate iboxlexp
                 call scale_points(iboxsubcenters,subcenters,
     1                8,boxsize(ilev))
                 call scale_points(iboxsrc,subpts,
     1                npbox,boxsize(ilev))
                 call dreorderf(2*nd,npbox,pot(1,ibox),
     1                iboxpot,iboxsrcind)
                 do i=1,8
                   jstart=iboxfl(1,i)
                   jend=iboxfl(2,i)
                   npts=jend-jstart+1
                   call h3dtaevalp(nd,zk,rscales(ilev),
     1                  subcenters(1,i),iboxlexp(1,i),
     2                  nterms(ilev),subpts(1,jstart),npts,
     3                  iboxpot(1,jstart),wlege,nlege)
                 enddo
                 call dreorderi(2*nd,npbox,iboxpot,pot(1,ibox),
     1                iboxsrcind)
                 deallocate(iboxlexp)
               endif
c and of handel mp eval for list3
            enddo
C$OMP END PARALLEL DO       

            deallocate(xshift,yshift,zshift,rlsc,tmp,tmp2)
            deallocate(carray,dc,rdplus,rdminus,rdsq3,rdmsq3)

            deallocate(mexpf1,mexpf2,mexpp1,mexpp2,mexppall,mexp)
            deallocate(fexp,fexpback)

            deallocate(pgboxwexp)
         else
            nquad2 = nterms(ilev)*2.2
            nquad2 = max(6,nquad2)
            ifinit2 = 1
            ier = 0

            call legewhts(nquad2,xnodes,wts,ifinit2)

            radius = boxsize(ilev)/2*sqrt(3.0d0)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,nl2,i,jbox) SCHEDULE(DYNAMIC)
            do ibox = itree(2*ilev+1),itree(2*ilev+2)

               nl2 = nlist2(ibox) 
               do i =1,nl2
                 jbox = list2(i,ibox) 
                   call h3dmploc(nd,zk,rscales(ilev),
     1               centers(1,jbox),
     1               rmlexp(iaddr(1,jbox)),nterms(ilev),
     2               rscales(ilev),centers(1,ibox),
     2               rmlexp(iaddr(2,ibox)),nterms(ilev),
     3               radius,xnodes,wts,nquad2)
               enddo
           enddo
C$OMP END PARALLEL DO        
         endif
      enddo
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      timeinfo(3) = time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== Step 4 (split loc) ===*',i,0)

      call cpu_time(time1)
C$        time1=omp_get_wtime()
      do ilev = 2,nlevels-1

        nquad2 = nterms(ilev)*2
        nquad2 = max(6,nquad2)
        ifinit2 = 1
        call legewhts(nquad2,xnodes,wts,ifinit2)
        radius = boxsize(ilev+1)/2*sqrt(3.0d0)

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,i,jbox,istart,iend,npts)
         do ibox = itree(2*ilev+1),itree(2*ilev+2) 
           do i=1,8
             jbox = itree(iptr(5)+8*(ibox-1)+i-1)
             if(jbox.gt.0) then
               call h3dlocloc(nd,zk,rscales(ilev),
     1           centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2           nterms(ilev),rscales(ilev+1),centers(1,jbox),
     3           rmlexp(iaddr(2,jbox)),nterms(ilev+1),
     4           radius,xnodes,wts,nquad2)
              endif
            enddo
         enddo
C$OMP END PARALLEL DO         
      enddo
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      timeinfo(4) = time2-time1




c
c
c       step  evaluate local expansions
c
      call cpu_time(time1)
C$       time1 = omp_get_wtime()      
      if(ifprint.ge.1)
     $    call prinf('=== Step 5 (loc eval) ===*',i,0)
      do ilev=2,nlevels
        print *, "ilev=",ilev,ilev
        neval = 0
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4)+ibox-1)
          if(nchild.eq.0) then
            neval = neval + 1
            ijboxlist(1,neval) = ibox
            ijboxlist(2,neval) = ibox
          endif
        enddo

        if(neval.gt.0) then
          nmp = (nterms(ilev)+1)*(2*nterms(ilev)+1) 
          allocate(rhs(nmp,neval),vals(npbox,neval))
          
          imp = 2
          call gather_mploc_vals(neval,ijboxlist,rmlexp,iaddr,imp,
     1      itree(iptr(2)),nboxes,nterms,nmp,rhs)

          ac = 1
          bc = 0
          call zgemm('n','n',npbox,neval,nmp,ac,tamat(itamat(ilev)),
     1       npbox,rhs,nmp,bc,vals,npbox)

          call scatter_vals(neval,ijboxlist,pot,npbox,nboxes,vals)

          deallocate(rhs,vals)
        endif
      enddo
      call cpu_time(time2)
C$      time2 = omp_get_wtime()     
      
      timeinfo(5) = time2-time1


c
c
c       step 6, handle list 1 procesing
c 
      call cpu_time(time1)
C$       time1 = omp_get_wtime()      
      if(ifprint.ge.1)
     $    call prinf('=== Step 6 (list1) ===*',i,0)

      allocate(nlist1_detailed(nboxes),list1_detailed(139,nboxes))

      call get_list1(nboxes,nlevels,itree,ltree,iptr,
     1   centers,boxsize,nlist1_detailed,list1_detailed)


      ldtab = 10*ncbox
      potcoefs = 0
      allocate(tabcoll(ncbox,ncbox,4))
      allocate(tabbtos(ncbox,ncbox,3),tabstob(ncbox,ncbox,3))
      allocate(tabtmp(ncbox,ncbox))

c
c      load table symmetries
c
      call loadsymsc(iref,idimp,iflip)
      call loadsymsbtos(irefbtos,idimpbtos,iflipbtos)
      call loadsymsstob(irefstob,idimpstob,iflipstob)

      ndeg = norder - 1
      do ilev=0,nlevels

c
c         check how many boxes in list 1 at this level
c
        nlist1lev = 0
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
           if(nlist1(ibox).gt.0) nlist1lev = nlist1lev + 1
        enddo

c
c        if number of boxes in list1 > 0 at this level,
c          then compute near field quadrature

        if(nlist1lev.gt.0) then
          print *, "ilev= ",ilev
          zk2 = zk*boxsize(ilev)/2.0d0
          ac = boxsize(ilev)**2/4.0d0
          bc = 0

          call splitreftab3dcc(tab(itab(ilev)),ldtab,tabcoll,tabbtos,
     1        tabstob,ncbox)
          
          call prin2('done splitting table*',i,0)

c
c           extract subtype of boxes in list1
c
          iboxstart = itree(2*ilev+1)
          iboxend = itree(2*ilev+2) 

c
c           handle colleagues
c
          do ibtype=1,27
            ntype = 0
            call get_list1boxes_type(ibtype,iboxstart,iboxend,
     1            nboxes,nlist1_detailed,
     1            list1_detailed,ijboxlist,ntype)

            if(ntype.gt.0) then
               allocate(rhs(ncbox,ntype),coefs(ncbox,ntype))
cc               print *,iref(ibtype),idimp(1:3,ibtype),iflip(1:3,ibtype)
cc               call prinf('iref=*',iref(ibtype),1)

               call buildtabfromsyms3dcc(ndeg,ttype,iref(ibtype),
     1           idimp(1,ibtype),iflip(1,ibtype),tabcoll,tabtmp,
     2           ncbox)
cc               call prinf('ibtype=*',ibtype,1)
cc               call prin2('tabtmp=*',tabtmp,npbox*ncbox*2)
               
               call gather_vals(ntype,ijboxlist,fcoefs,ncbox,nboxes,rhs)

cc               call prin2('rhs=*',rhs,ncbox*ntype*2)
              
               call zgemm('n','n',ncbox,ntype,ncbox,ac,tabtmp,ncbox,
     1             rhs,ncbox,bc,coefs,ncbox)
               
cc               call prin2('vals=*',vals,2*npbox*ntype) 

               call scatter_vals(ntype,ijboxlist,potcoefs,ncbox,
     1            nboxes,coefs)

               deallocate(rhs,coefs)
               
            endif
          enddo

          print *, "Done with neighbors"
c
c           handle big to small
c
          do ibtype=1,56
            ntype = 0
            ii = ibtype + 27
            call get_list1boxes_type(ii,iboxstart,iboxend,
     1            nboxes,nlist1_detailed,
     1            list1_detailed,ijboxlist,ntype)

            if(ntype.gt.0) then
               allocate(rhs(ncbox,ntype),coefs(ncbox,ntype))

               call buildtabfromsyms3dcc(ndeg,ttype,irefbtos(ibtype),
     1           idimpbtos(1,ibtype),iflipbtos(1,ibtype),tabbtos,tabtmp,
     2           ncbox)
cc               call prinf('ibtype=*',ibtype,1)
cc               call prin2('tabtmp=*',tabtmp,npbox*ncbox*2)
               
               call gather_vals(ntype,ijboxlist,fcoefs,ncbox,nboxes,rhs)

cc               call prin2('rhs=*',rhs,ncbox*ntype*2)
              
               call zgemm('n','n',ncbox,ntype,ncbox,ac,tabtmp,ncbox,
     1             rhs,ncbox,bc,coefs,ncbox)
               
cc               call prin2('vals=*',vals,2*npbox*ntype) 

               call scatter_vals(ntype,ijboxlist,potcoefs,ncbox,
     1            nboxes,coefs)

               deallocate(rhs,coefs)
               
            endif
          enddo
c
c           handle small to big
c
          do ibtype=1,56
            ntype = 0
            ii = ibtype + 27 + 56
            call get_list1boxes_type(ii,iboxstart,iboxend,
     1            nboxes,nlist1_detailed,
     1            list1_detailed,ijboxlist,ntype)

            if(ntype.gt.0) then
               allocate(rhs(ncbox,ntype),coefs(ncbox,ntype))

               call buildtabfromsyms3dcc(ndeg,ttype,irefstob(ibtype),
     1           idimpstob(1,ibtype),iflipstob(1,ibtype),tabstob,tabtmp,
     2           ncbox)
cc               call prinf('ibtype=*',ibtype,1)
cc               call prin2('tabtmp=*',tabtmp,npbox*ncbox*2)
               
               call gather_vals(ntype,ijboxlist,fcoefs,ncbox,nboxes,rhs)

cc               call prin2('rhs=*',rhs,ncbox*ntype*2)
              
               call zgemm('n','n',ncbox,ntype,ncbox,ac,tabtmp,ncbox,
     1             rhs,ncbox,bc,coefs,ncbox)
               
cc               call prin2('vals=*',vals,2*npbox*ntype) 

               call scatter_vals(ntype,ijboxlist,potcoefs,ncbox,
     1             nboxes,coefs)

               deallocate(rhs,coefs)
               
            endif
          enddo
        endif
      enddo



cc      call prin2('pot=*',pot,2*npbox*nboxes)

      deallocate(fimat)


      call prin2('done with fmm*',i,0)
c
c
c  convert coefs back to vals and add to potential
c
      ra = 1.0d0
      rb = 1.0d0
      do ibox=1,nboxes
        call dgemm('n','t',2,npbox,ncbox,ra,potcoefs(1,ibox),
     1    2,vmat,npbox,rb,pot(1,ibox),2)
      enddo
      
c
c   recompute coefs correctly now
c
c
      potcoefs = 0

      ra = 1.0d0
      rb = 0.0d0
      do ibox=1,nboxes
        call dgemm('n','t',2,ncbox,npbox,ra,pot(1,ibox),
     1    2,umat,ncbox,rb,potcoefs(1,ibox),2)
      enddo
      
      call cpu_time(time2)
C$      time2 = omp_get_wtime()      

      timeinfo(6) = time2-time1
      call prin2('fmm timeinfo=*',timeinfo,6)


      return
      end

c
c
c
c
c

      subroutine scatter_vals(n,ijlist,pot,npbox,nboxes,vals)
      implicit real *8 (a-h,o-z)

      integer ijlist(2,n),ncbox,nboxes
      complex *16 pot(npbox,nboxes),vals(npbox,n)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,j)
      do i=1,n
        ibox = ijlist(2,i)
        do j=1,npbox
          pot(j,ibox) = pot(j,ibox) + vals(j,i)
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

      subroutine gather_vals(n,ijlist,fcoefs,ncbox,nboxes,rhs)
      implicit real *8 (a-h,o-z)

      integer ijlist(2,n),ncbox,nboxes
      complex *16 rhs(ncbox,n),fcoefs(ncbox,nboxes)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,j)
      do i=1,n
        ibox = ijlist(1,i)
        do j=1,ncbox
          rhs(j,i) = fcoefs(j,ibox)
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

      subroutine gather_mploc_vals(n,ijlist,rmlexp,iaddr,imp,
     1   ilevel,nboxes,nterms,nmp,mpvals)
      implicit real *8 (a-h,o-z)
      integer n,ijlist(2,n),imp,ilevel(nboxes)
      integer *8 iaddr(2,nboxes)
      integer nboxes,nterms(*)
      real *8 rmlexp(*)
      complex *16 mpvals(nmp,n),ima

      data ima/(0.0d0,1.0d0)/

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,istart,j)
      do i=1,n
        ibox = ijlist(1,i)
        istart = iaddr(imp,ibox)
        do j=1,nmp
          mpvals(j,i) = rmlexp(istart+2*j-2) + ima*rmlexp(istart+2*j-1)
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

      subroutine scatter_mploc_vals(n,ijlist,rmlexp,iaddr,imp,
     1   ilevel,nboxes,nterms,nmp,mpvals)
      implicit real *8 (a-h,o-z)
      integer n,ijlist(2,n),imp,ilevel(nboxes)
      integer *8 iaddr(2,nboxes)
      integer nboxes,nterms(*)
      real *8 rmlexp(*)
      complex *16 mpvals(nmp,n),ima

      data ima/(0.0d0,1.0d0)/

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,istart,j)
      do i=1,n
        ibox = ijlist(2,i)
        istart = iaddr(imp,ibox)
        do j=1,nmp
          rmlexp(istart+2*j-2) = rmlexp(istart+2*j-2)+real(mpvals(j,i)) 
          rmlexp(istart+2*j-1) = rmlexp(istart+2*j-1)+imag(mpvals(j,i))
        enddo
      enddo
C$OMP END PARALLEL DO      

      return
      end

c
c

      subroutine gather_mploc_vals_t(n,ijlist,rmlexp,iaddr,imp,
     1   ilevel,nboxes,nterms,nmp,mpvals)
      implicit real *8 (a-h,o-z)
      integer n,ijlist(2,n),imp,ilevel(nboxes)
      integer *8 iaddr(2,nboxes)
      integer nboxes,nterms(*)
      real *8 rmlexp(*)
      complex *16 mpvals(n,nmp),ima

      data ima/(0.0d0,1.0d0)/

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,istart,j)
      do i=1,n
        ibox = ijlist(1,i)
        istart = iaddr(imp,ibox)
        do j=1,nmp
          mpvals(i,j) = rmlexp(istart+2*j-2) + ima*rmlexp(istart+2*j-1)
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

      subroutine scatter_mploc_vals_t(n,ijlist,rmlexp,iaddr,imp,
     1   ilevel,nboxes,nterms,nmp,mpvals)
      implicit real *8 (a-h,o-z)
      integer n,ijlist(2,n),imp,ilevel(nboxes)
      integer *8 iaddr(2,nboxes)
      integer nboxes,nterms(*)
      real *8 rmlexp(*)
      complex *16 mpvals(n,nmp),ima

      data ima/(0.0d0,1.0d0)/

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,istart,j)
      do i=1,n
        ibox = ijlist(2,i)
        istart = iaddr(imp,ibox)
        do j=1,nmp
          rmlexp(istart+2*j-2) = rmlexp(istart+2*j-2)+real(mpvals(i,j)) 
          rmlexp(istart+2*j-1) = rmlexp(istart+2*j-1)+imag(mpvals(i,j))
        enddo
      enddo
C$OMP END PARALLEL DO      

      return
      end

