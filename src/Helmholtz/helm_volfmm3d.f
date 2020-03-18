      subroutine helmholtz_volume_fmm(eps,zk,nboxes,nlevels,ltree,
     1   itree,iptr,norder,ncbox,type,fcoefs,centers,boxsize,npbox,
     2   pot)
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
c         type - character *1
c            type of coefs provided, total order ('t') or full order('f')
c         fcoefs - double complex (ncbox,nboxes)
c           tensor product legendre expansions of the right hand side
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

      implicit real *8 (a-h,o-z)
      real *8 eps
      complex *16 zk,zk2
      integer nboxes,nlevels,ltree
      integer itree(ltree),iptr(8),ncbox,npbox
      complex *16 fcoefs(ncbox,nboxes)
      complex *16 pot(npbox,nboxes)
      double precision boxsize(0:nlevels),centers(3,nboxes)

      double precision, allocatable :: rscales(:)
      integer, allocatable :: nterms(:)
      double precision, allocatable :: rmlexp(:)
      integer *8, allocatable :: iaddr(:,:)
      integer lmptemp
      integer *8 lmptot
      double precision, allocatable :: mptemp(:),mptemp2(:)

      double precision, allocatable :: wlege(:)

      double precision xtargtmp(3)
      complex *16 pottmp,pottmpex
      character *1 type
      double precision, allocatable :: xnodes(:),wts(:)

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

      double precision, allocatable :: rsc(:)
      integer, allocatable :: ilevrel(:)
      complex *16, allocatable :: mpcoefsmat(:,:),tab(:,:)
      complex *16, allocatable :: tabcoll(:,:,:),tabbtos(:,:,:),
     1   tabstob(:,:,:)
      complex *16, allocatable :: tabtmp(:,:),tamat(:,:)
      complex *16, allocatable :: rhs(:,:),vals(:,:)
      complex *16 ac,bc,ima

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
      double precision timeinfo(8)

      integer iref(100),idimp(3,100),iflip(3,100)

      data ima/(0.0d0,1.0d0)/



      ifprint = 1

      done = 1
      pi = atan(done)*4

      ntmax = 1000
      allocate(rlams(ntmax),whts(ntmax),nfourier(ntmax),
     1   nphysical(ntmax))


      max_nodes = 10000
      allocate(xnodes(max_nodes))
      allocate(wts(max_nodes))

      do i=1,8
        timeinfo(i) = 0
      enddo


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
      call computemnlists(nlevels,nboxes,itree(iptr(1)),
     1       boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2       itree(iptr(5)),isep,itree(iptr(6)),mnbors,itree(iptr(7)),
     3       mnlist1,mnlist2,mnlist3,mnlist4)

      allocate(list1(mnlist1,nboxes),list2(mnlist2,nboxes))
      allocate(list3(mnlist3,nboxes),list4(mnlist4,nboxes))
      allocate(nlist1(nboxes))
      allocate(nlist2(nboxes))
      allocate(nlist3(nboxes))
      allocate(nlist4(nboxes))

      call computelists(nlevels,nboxes,itree(iptr(1)),
     1   boxsize,centers,itree(iptr(3)),itree(iptr(4)),
     2   itree(iptr(5)),isep,itree(iptr(6)),mnbors,itree(iptr(7)),
     3   nlist1,mnlist1,list1,nlist2,mnlist2,list2,
     4   nlist3,mnlist3,list3,nlist4,mnlist4,list4)


      allocate(ijboxlist(2,nboxes))

c
c       find scales and number of terms required at each of
c       the levels
c

      allocate(rscales(0:nlevels),nterms(0:nlevels))
cc      call prinf('nboxes=*',nboxes,1)
cc      call prinf('nlevels=*',nlevels,1)
cc      call prin2('zk=*',zk,2)


 
      nmax = 0
      do ilev = 0,nlevels
        rscales(ilev) = boxsize(ilev)*abs(zk)

        if(rscales(ilev).gt.1) rscales(ilev) = 1.0d0
        call h3dterms(boxsize(ilev),zk,eps,nterms(ilev))
        if(nterms(ilev).gt.nmax) nmax = nterms(ilev)
      enddo
cc      call prinf('nterms=*',nterms,nlevels+1)
cc      call prin2('rscales=*',rscales,nlevels+1)

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


      

c
c
c        step 1: convert coeffs to multipole expansions
c
      call cpu_time(time1)
C$     time1 = omp_get_wtime()      
    
      if(ifprint.ge.1) 
     $   call prinf("=== STEP 1 (coefs -> mp) ====*",i,0)

cc      call prinf('ltree=*',ltree,1)
cc      call prinf('iptr=*',iptr,8)

      allocate(ilevrel(0:nlevels))
      ilevrel(0) = 0
      ilevrel(1) = 0
      do ilev = 2,nlevels
        ilevrel(ilev) = 0
        do ibox = itree(2*ilev+1),itree(2*ilev+2)
          nchild = itree(iptr(4) + ibox-1)
          if(nchild.eq.0) ilevrel(ilev) = 1
        enddo
      enddo

cc      call prinf('ilevrel=*',ilevrel,nlevels+1)
     
      do ilev=2,nlevels
        nmp  = (nterms(ilev)+1)*(2*nterms(ilev)+1)
        if(ilevrel(ilev).eq.1) then
          nq = 8
          allocate(mpcoefsmat(nmp,ncbox))

cc          call prinf('ilev=*',ilev,1)
cc          call prin2('zk=*',zk,2)
cc          call prin2('rscales=*',rscales(ilev),1)
cc          call prinf('nterms=*',nterms(ilev),1)
cc          call prin2('boxsize=*',boxsize(ilev),1)
cc          call prina('type=*',type,1)
cc          call prinf('nq=*',nq,1)
cc          call prinf('ncbox=*',ncbox,1)
cc          call prinf('nlege=*',nlege,1)
cc          print *, size(mpcoefsmat)


          call h3ddensmpmat(zk,rscales(ilev),nterms(ilev),
     1     boxsize(ilev),type,norder,nq,wlege,nlege,mpcoefsmat,
     2     nmp)
          do ibox = itree(2*ilev+1),itree(2*ilev+2)
            nchild = itree(iptr(4)+ibox-1)
            if(nchild.eq.0) then
               ac = 1.0d0
               bc = 0.0d0
              call zgemv('n',nmp,ncbox,ac,mpcoefsmat,nmp,
     1         fcoefs(1,ibox),1,bc,rmlexp(iaddr(1,ibox)),1) 
            endif
          enddo
        endif
      enddo

      call cpu_time(time2)
C$     time2 = omp_get_wtime()   

      timeinfo(1) = time2-time1


c       
      if(ifprint .ge. 1)
     $      call prinf('=== STEP 3 (merge mp) ====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c
c
c
c       note: faster multipole to multipole operator
c        possible by storing matrix from children
c        to parents
c

      do ilev=nlevels-1,0,-1
         nquad2 = nterms(ilev)*2.5
         nquad2 = max(6,nquad2)
         ifinit2 = 1
         call legewhts(nquad2,xnodes,wts,ifinit2)
cc         call prinf('nquad2=*',nquad2,1)
cc         call prin2('xnodes=*',xnodes,nquad2)
cc         call prin2('wts=*',wts,nquad2)
cc         call prinf('nd=*',nd,1)
         radius = boxsize(ilev)/2*sqrt(3.0d0)

cc         call prin2('radius=*',radius,1)


C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,i,jbox,nchild)
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
      timeinfo(3)=time2-time1


      xtargtmp(1) = 50.1d0
      xtargtmp(2) = 3.2d0
      xtargtmp(3) = 2.1d0

      x = xtargtmp(1)-0.01d0
      y = xtargtmp(2)
      z = xtargtmp(3)
      
      r = sqrt(x**2 + y**2 + z**2)

      ss = 1.0d0/13.0d0

      
      c = exp(-zk**2/2.0d0*ss**2)*(2.0d0*pi)**1.5d0*ss**3
      pottmpex = ima*c/r*sin(zk*r)

      print *, c
      print *, ima
      print *, r
      print *, zk

      thresh = 1.0d-16

      pottmp = 0.0d0
      
      call h3dmpevalp(nd,zk,rscales(0),centers,rmlexp(iaddr(1,1)),
     1   nterms(0),xtargtmp,1,pottmp,wlege,nlege,thresh)

      call prin2('pottmp=*',pottmp,2)
      call prin2('pottmpex=*',pottmpex,2)
      erra = abs(imag(pottmpex-pottmp))/abs(pottmpex)
      call prin2('error=*',erra,1)

      if(ifprint.ge.1)
     $    call prinf('=== Step 4 (mp to loc) ===*',i,0)
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
cc         if(1.eq.0) then
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
          
cc           r1 = rscales(ilev)
           r1 = 1.0d0
           rsc(0) = 1.0d0
           do i=1,nterms(ilev)
             rsc(i) = rsc(i-1)*r1
           enddo

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
     7            xshift,yshift,zshift,fexpback,rlsc)
                  
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
     9            fexpback,rlsc)

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
     9            fexpback,rlsc)
               endif
            enddo
C$OMP END PARALLEL DO       

            deallocate(xshift,yshift,zshift,rlsc,tmp,tmp2)
            deallocate(carray,dc,rdplus,rdminus,rdsq3,rdmsq3)

            deallocate(mexpf1,mexpf2,mexpp1,mexpp2,mexppall,mexp)
            deallocate(fexp,fexpback)

         else
            nquad2 = nterms(ilev)*2.2
            nquad2 = max(6,nquad2)
            ifinit2 = 1
            ier = 0

            call legewhts(nquad2,xnodes,wts,ifinit2)

            radius = boxsize(ilev)/2*sqrt(3.0d0)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,nl2,i,jbox)
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
      timeinfo(4) = time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== Step 5 (split loc) ===*',i,0)

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
      timeinfo(5) = time2-time1





c
c
c       step 7 evaluate local expansions
c

      call cpu_time(time1)
C$      time1 = omp_get_wtime()      
      if(ifprint.ge.1)
     $    call prinf('=== Step 7 (loc eval) ===*',i,0)
      do ilev=0,nlevels
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
          allocate(tamat(npbox,nmp),rhs(nmp,neval),vals(npbox,neval))
          call h3dtaevalgridmatp(zk,rscales(ilev),nterms(ilev),
     1      boxsize(ilev),norder,wlege,nlege,tamat,npbox)
          
cc          call prin2('tamat=*',tamat(1,+1),npbox*2)
          
cc          call prinf('neval=*',neval,1)
cc          call prinf('nboxes=*',nboxes,1)
cc          call prinf('nmp=*',nmp,1)
          imp = 2
          call gather_mploc_vals(neval,ijboxlist,rmlexp,iaddr,imp,
     1      itree(iptr(2)),nboxes,nterms,nmp,rhs)

          ac = 1
          bc = 0
          call zgemm('n','n',npbox,neval,nmp,ac,tamat,npbox,
     1       rhs,nmp,bc,vals,npbox)

cc          call prin2('vals=*',vals,2*neval*npbox)
          call scatter_vals(neval,ijboxlist,pot,npbox,nboxes,vals)

          deallocate(rhs,vals)
        endif
      enddo
      call cpu_time(time2)
C$      time2 = omp_get_wtime() 

      timeinfo(7) = time2-time1



c
c
c       step 8, handle list 1 procesing
c 

      call cpu_time(time1)
C$      time1 = omp_get_wtime()      
      allocate(nlist1_detailed(nboxes),list1_detailed(139,nboxes))

      call get_list1(nboxes,nlevels,itree,ltree,iptr,
     1   centers,boxsize,nlist1_detailed,list1_detailed)
      
cc      call prinf('nlist1_detailed=*',nlist1_detailed,nboxes)


      


      ntarg0 = 10*npbox
      allocate(tab(ntarg0,ncbox),tabcoll(npbox,ncbox,4))
      allocate(tabbtos(npbox,ncbox,3),tabstob(npbox,ncbox,3))
      allocate(tabtmp(npbox,ncbox))

c
c      load table symmetries
c
      call loadsymsc(iref,idimp,iflip)

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
          call h3dtabp_ref(ndeg,zk2,eps,tab,ntarg0)
          call splitreftab3d(tab,ntarg0,tabcoll,tabbtos,tabstob,
     1        npbox,ncbox)
          
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
               allocate(rhs(ncbox,ntype),vals(npbox,ntype))
cc               print *,iref(ibtype),idimp(1:3,ibtype),iflip(1:3,ibtype)
cc               call prinf('iref=*',iref(ibtype),1)

               call buildtabfromsyms3d(ndeg,type,iref(ibtype),
     1           idimp(1,ibtype),iflip(1,ibtype),tabcoll,tabtmp,
     2           npbox,ncbox)
cc               call prinf('ibtype=*',ibtype,1)
cc               call prin2('tabtmp=*',tabtmp,npbox*ncbox*2)
               
               call gather_vals(ntype,ijboxlist,fcoefs,ncbox,nboxes,rhs)

cc               call prin2('rhs=*',rhs,ncbox*ntype*2)
              
               call zgemm('n','n',npbox,ntype,ncbox,ac,tabtmp,npbox,
     1             rhs,ncbox,bc,vals,npbox)
               
cc               call prin2('vals=*',vals,2*npbox*ntype) 

               call scatter_vals(ntype,ijboxlist,pot,npbox,nboxes,vals)


               deallocate(rhs,vals)
               
            endif
          enddo
        endif
      enddo
      call cpu_time(time2)
C$      time2 = omp_get_wtime()    

      timeinfo(8) = time2-time1

      call prin2('time = *',timeinfo,8)
      d = 0
      do i=1,8
        d = d + timeinfo(i)
      enddo
      call prin2('total time=*',d,1)



cc      call prin2('pot=*',pot,2*npbox*nboxes)

      call prin2('done with fmm*',i,0)

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

      do i=1,n
        ibox = ijlist(2,i)
        do j=1,npbox
          pot(j,ibox) = pot(j,ibox) + vals(j,i)
        enddo
      enddo

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

      do i=1,n
        ibox = ijlist(1,i)
        do j=1,ncbox
          rhs(j,i) = fcoefs(j,ibox)
        enddo
      enddo

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


      do i=1,n
        ibox = ijlist(1,i)
        istart = iaddr(imp,ibox)
        do j=1,nmp
          mpvals(j,i) = rmlexp(istart+2*j-2) + ima*rmlexp(istart+2*j-1)
        enddo

      enddo

      return
      end

