cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c     TABLE GENERATION UTILITIES
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This file contains the routines for evaluating 
c     the near-field volume potential using an 
c     anti-Helmholtzian
c      
c     None of these routines is really designed to be user-callable
c      
c     
c     



      subroutine h3danti_form(ndeg,nup,type,zk,ahmat,ldahmat,derrmax,
     1     errsup,errsdown)
c
c     This routine forms the coefficients of an anti-helmholtzian
c     for each polynomial in the specified basis by computing both
c     upward (with nup steps) and downward recurrences and choosing
c     the option with the best backward error.
c      
c      
c     input:
c
c     ndeg - integer
c        the degree of the basis
c     nup - the number of upward iterations to run 
c     type - character
c        the type of the basis (for now, total degree
c        is the only implemented option)
c        type = 't', total degree
c        type = 'f', full tensor
c     zk - complex double
c        Helmholtz parameter
c     ldahmat - integer
c         leading dimension of ahmat, must be at least 
c         the number of tensor polynomials of degree ndeg+2*nup
c         see legetens_npol_3d
c      
c     output:
c    
c     ahmat - complex array (ldahmat,*)
c        column i gives the coefficients of the anti-Helmholtzian of
c        the ith function in the polynomial basis determined by ndeg
c        and type. The coefficients are in the poylnomial basis
c        determined by ndegout=ndeg+2*nup and type.
c     derrmax - real *8 worst case root mean square (backward) error
c        for the output coefficients (i.e. the difference between
c        spectral Helmholtz operator applied to the jth column of
c     coefficients and the standard basis vector e_j)
c     errsup - real *8 array (number of polynomials for ndeg and type)
c        error (including stability fudge factor) for upward
c        recurrence
c     errsdown - real *8 array (number of polynomials for ndeg and type)
c        error (including stability fudge factor) for downward
c        recurrence
      
      

      implicit real *8 (a-h,o-z)
      complex *16 ahmat(ldahmat,*)
      real *8 errsup(*), errsdown(*)
      complex *16, allocatable :: ahmatup(:,:), ahmatdown(:,:)
      complex *16, allocatable :: resup(:,:), resdown(:,:)      
      integer, allocatable :: ind2p(:,:), ind2pout(:,:)
      integer, allocatable :: ip2ind(:,:,:), ip2indout(:,:,:)      
      
      real *8, allocatable :: sumsup(:), sumsdown(:), stabup(:),
     1     stabdown(:)

      complex *16 im, zero, one, zk
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /

      character type
      
      ndegout = ndeg + nup*2
      
      call legetens_npol_3d(ndeg,type,npol)
      call legetens_npol_3d(ndegout,type,npolout)

      allocate(ind2p(3,npol),ind2pout(3,npolout))
      allocate(ip2ind(0:ndeg,0:ndeg,0:ndeg),
     1     ip2indout(0:ndegout,0:ndegout,0:ndegout))

      call legetens_ind2pow_3d(ndeg,type,ind2p)
      call legetens_ind2pow_3d(ndegout,type,ind2pout)      
      call legetens_pow2ind_3d(ndeg,type,ip2ind)
      call legetens_pow2ind_3d(ndegout,type,ip2indout)      
      allocate(ahmatup(npolout,npol),ahmatdown(npol,npol))

      call h3danti_legeup(ndeg,nup,type,zk,ahmatup,npolout)
      call h3danti_legedownfast(ndeg,type,zk,ahmatdown,npol)

c     
      allocate(resup(npolout,npol),resdown(npol,npol))

c$omp parallel do private(i,ndc,j)
      do i = 1,npol
         ndc = 2
         call legetens_lape_3d(ndc,ndeg,type,ahmatdown(1,i),
     1        resdown(1,i))
         call legetens_lape_3d(ndc,ndegout,type,ahmatup(1,i),
     1        resup(1,i))

         do j = 1,npol
            resdown(j,i) = resdown(j,i)+zk**2*ahmatdown(j,i)
         enddo
         do j = 1,npolout
            resup(j,i) = resup(j,i)+zk**2*ahmatup(j,i)
         enddo
      enddo
c$omp end parallel do      


      
c     determine backward error (l2)

      allocate(sumsup(npol),sumsdown(npol),stabup(npol),stabdown(npol))

c$omp parallel do private(i,iin,jin,kin,iout,jout,kout,j)      
      do i = 1,npol
         sumsup(i) = 0
         sumsdown(i) = 0
         stabup(i) = 0
         stabdown(i) = 0
         iin = ind2p(1,i)
         jin = ind2p(2,i)
         kin = ind2p(3,i)
         do j = 1,npolout
            iout = ind2pout(1,j)
            jout = ind2pout(2,j)
            kout = ind2pout(3,j)
            if (iin .eq. iout .and. jin .eq. jout
     1           .and. kin .eq. kout) then
               sumsup(i) = sumsup(i) + abs(one-resup(j,i))**2
            else
               sumsup(i) = sumsup(i) + abs(resup(j,i))**2
            endif
            stabup(i) = stabup(i)+abs(ahmatup(j,i))**2
         enddo

         do j = 1,npol
            if (j .eq. i) then
               sumsdown(i) = sumsdown(i) + abs(one-resdown(j,i))**2
            else
               sumsdown(i) = sumsdown(i) + abs(resdown(j,i))**2
            endif
            stabdown(i) = stabdown(i)+abs(ahmatdown(j,i))**2
         enddo

         sumsdown(i) = sqrt(sumsdown(i)/npol)
         sumsup(i) = sqrt(sumsup(i)/npol)
         stabup(i) = sqrt(stabup(i)/npol)
         stabdown(i) = sqrt(stabdown(i)/npol)         
      enddo
c$end parallel do


c     choose best for each input poly and get max error

      derrmax = 0.0d0
      epsfac = 2d-16
      
c$omp parallel do private(i,j,i1,i2,i3,jout) reduction(max:derrmax)      
      do i = 1,npol
         do j = 1,npolout
            ahmat(j,i) = 0
         enddo

         errdown = sumsdown(i) + stabdown(i)*epsfac
         errup = stabup(i)*epsfac+sumsup(i)
         errsdown(i) = errdown
         errsup(i) = errup         
         if (errdown .lt. errup) then
            do j = 1,npol
               i1 = ind2p(1,j)
               i2 = ind2p(2,j)
               i3 = ind2p(3,j)
               jout = ip2indout(i1,i2,i3)
               ahmat(jout,i) = ahmatdown(j,i)
            enddo
            derrmax = max(derrmax,sumsdown(i)+epsfac*stabdown(i))            
         else
            do j = 1,npolout
               ahmat(j,i) = ahmatup(j,i)
            enddo
            derrmax = max(derrmax,sumsup(i)+epsfac*stabup(i)) 
         endif

      enddo
c$omp end parallel do      

      return
      end
      
      subroutine h3danti_legedown(ndeg,type,zk,ahmat,ldahmat)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Apply the downward recurrence formula to obtain the
c     anti-Helmholtzian on Legendre polynomials
c
c     input:
c
c     ndeg - integer
c        the degree of the basis
c     type - character
c        the type of the basis (for now, total degree
c        is the only implemented option)
c        type = 't', total degree
c        type = 'f', full tensor
c     zk - complex double
c        Helmholtz parameter
c     ldahmat - integer
c        leading dimension of ahmat
c     output:
c     
c     ahmat - complex array (ldahmat,*)
c        column i gives the coefficients of the anti-Helmholtzian of
c        the ith function in the polynomial basis determined by ndeg
c        and type. The coefficients are in that same basis for the
c        downward recurrence, i.e. this function gives an anti-
c        Helmholtzian that's of the same degree.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer ndeg, ldahmat
      complex *16 ahmat(ldahmat,*), zk
      character type
c     local
      complex *16, allocatable :: lapmat(:,:), temp1(:,:), temp2(:,:)
      complex *16 one, zero, zk2, zk2inv, alpha, beta
      real *8, allocatable :: dlap(:,:)
      data one / (1.0d0,0.0d0) /
      data zero / (0.0d0,0.0d0) /
      
      integer npol, i, j, iii, niter

      call legetens_npol_3d(ndeg,type,npol)
      
      allocate(lapmat(npol,npol),dlap(npol,npol),temp1(npol,npol),
     1     temp2(npol,npol))

      call legetens_lapmat_3d(ndeg,ndeg,type,dlap,npol)

      zk2 = zk**2
      zk2inv = one/zk2
      
      do i = 1,npol
         do j = 1,npol
            lapmat(j,i) = dlap(j,i)*zk2inv
            temp1(j,i) = lapmat(j,i)*zk2inv
            ahmat(j,i) = zero
            if (j .eq. i) ahmat(j,i) = zk2inv
         enddo
      enddo

      niter = ndeg/2

      beta = one

      if (ndeg .lt. 2) return

c     first one manually

      do i = 1,npol
         do j = 1,npol
            ahmat(j,i) = ahmat(j,i) - temp1(j,i)
         enddo
      enddo

      alpha = -one
      beta = zero
      do iii = 2,niter
         
c     power up scaled laplacian matrix
         call zgemm('N','N',npol,npol,npol,alpha,lapmat,npol,
     1        temp1,npol,beta,temp2,npol)
         do i = 1,npol
            do j = 1,npol
               temp1(j,i) = temp2(j,i)
               ahmat(j,i) = ahmat(j,i) - temp1(j,i)
            enddo
         enddo
         
      enddo
      
      return
      end
c
c
c
c     


      subroutine h3danti_legedownfast(ndeg,type,zk,ahmat,ldahmat)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Apply the downward recurrence formula to obtain the
c     anti-Helmholtzian on Legendre polynomials
c
c     input:
c
c     ndeg - integer
c        the degree of the basis
c     type - character
c        the type of the basis (for now, total degree
c        is the only implemented option)
c        type = 't', total degree
c        type = 'f', full tensor
c     zk - complex double
c        Helmholtz parameter
c     ldahmat - integer
c        leading dimension of ahmat
c     output:
c     
c     ahmat - complex array (ldahmat,*)
c        column i gives the coefficients of the anti-Helmholtzian of
c        the ith function in the polynomial basis determined by ndeg
c        and type. The coefficients are in that same basis for the
c        downward recurrence, i.e. this function gives an anti-
c        Helmholtzian that's of the same degree.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer ndeg, ldahmat
      complex *16 ahmat(ldahmat,*), zk
      character type
c     local
      complex *16, allocatable :: lapmats1d(:,:,:), temp1(:,:)
      integer, allocatable :: ip2ind(:,:,:), ind2p(:,:)
      real *8, allocatable :: choose(:,:), rtemp1(:,:)
      complex *16 one, zero, zk2, zk2inv, alpha, beta, z2i,z2j,z2k
      real *8, allocatable :: dlap(:,:)
      data one / (1.0d0,0.0d0) /
      data zero / (0.0d0,0.0d0) /

      real *8 chl
      
      integer npol, i, j, iii, niter, n
      integer niterp, id,jd,kd,il,jl,kl,ilap,ip,kp,jp,indd

      call legetens_npol_3d(ndeg,type,npol)

      n = ndeg + 1
      
      allocate(ip2ind(0:ndeg,0:ndeg,0:ndeg),ind2p(3,npol))
      
      call legetens_pow2ind_3d(ndeg,type,ip2ind)
      call legetens_ind2pow_3d(ndeg,type,ind2p)

      niter = ndeg/2

c

      zk2 = zk**2
      zk2inv = one/zk2
      
      do i = 1,npol
         do j = 1,npol
            ahmat(j,i) = zero
            if (j .eq. i) ahmat(j,i) = one
         enddo
      enddo

      if (ndeg .lt. 2) goto 1000


c     power up to get necessary mutliples of d_xx matrix
      
      allocate(lapmats1d(0:ndeg,0:ndeg,0:niter),rtemp1(0:ndeg,0:ndeg))
      allocate(temp1(0:ndeg,0:ndeg))

      do i = 0,ndeg
         do j = 0,ndeg
            rtemp1(j,i) = 0.0d0
            lapmats1d(j,i,0) = zero
            if (j .eq. i) lapmats1d(j,i,0) = one
         enddo
      enddo
      
      call legecoeff_d2mat(ndeg,rtemp1,n)
      
      do i = 0,ndeg
         do j = 0,ndeg
            temp1(j,i) = -rtemp1(j,i)*zk2inv
            lapmats1d(j,i,1) = temp1(j,i)
         enddo
      enddo

      alpha = one
      beta = zero
      do iii = 2,niter
         call zgemm('N','N',n,n,n,alpha,temp1,n,
     1        lapmats1d(0,0,iii-1),n,beta,lapmats1d(0,0,iii),n)
      enddo
      
      allocate(choose(niter+1,niter+1))

C     compute choose(i,j) to be [(i-1) choose (j-1)].
C     (trinomial expansion of -lap/zk^2)
      
      choose(1,1) = 1.0d0
      do i = 2,niter+1
         choose(1,i) = 0.0d0
         choose(i,1) = 1.0d0
      enddo
      do i = 2, niter + 1
         do j = 2, niter + 1
            choose(i,j) = choose(i-1,j-1) + choose(i-1,j)
         enddo
      enddo


c     grab appropriate entries

      do iii = 1,npol
         ip = ind2p(1,iii)
         jp = ind2p(2,iii)
         kp = ind2p(3,iii)

         do ilap = 1,niter
            do il = 0,ilap
               if (il .gt. ip/2) cycle
               
               do jl = 0,ilap-il
                  
                  if (jl .gt. jp/2) cycle
                  kl = ilap-il-jl
                  if (kl .gt. kp/2) cycle
                  
                  chl = choose(ilap+1,il+1)*choose(ilap-il+1,jl+1)
                  
                  do kd = 0,kp-2*kl
                     z2k = lapmats1d(kd,kp,kl)
                     do jd = 0,jp-2*jl
                        z2j = lapmats1d(jd,jp,jl)
                        do id = 0,ip-2*il
                           z2i = lapmats1d(id,ip,il)
                           
                           indd = ip2ind(id,jd,kd)
                           ahmat(indd,iii) = ahmat(indd,iii) +
     1                          chl*z2k*z2j*z2i

                        enddo
                     enddo
                  enddo
                  
               enddo
            enddo
         enddo
      enddo

 1000 continue
      
      do i = 1,npol
         do j = 1,npol
            ahmat(j,i) = ahmat(j,i)*zk2inv
         enddo
      enddo
      
      return
      end
c
c
c
c     
      
      subroutine h3danti_legeup(ndeg,nup,type,zk,ahmat,ldahmat)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Apply the downward recurrence formula to obtain the
c     anti-Helmholtzian on Legendre polynomials
c
c     input:
c
c     ndeg - integer
c        the degree of the basis
c     nup - the number of upward iterations to run 
c     type - character
c        the type of the basis (for now, total degree
c        is the only implemented option)
c        type = 't', total degree
c        type = 'f', full tensor
c     zk - complex double
c        Helmholtz parameter
c     ldahmat - integer
c        leading dimension of ahmat
c     output:
c     
c     ahmat - complex array (ldahmat,*)
c        column i gives the coefficients of the anti-Helmholtzian of
c        the ith function in the polynomial basis determined by ndeg
c        and type. The coefficients are in that same basis for the
c        downward recurrence, i.e. this function gives an anti-
c        Helmholtzian that's of the same degree.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer ndeg, ldahmat, nup
      complex *16 ahmat(ldahmat,*), zk
      character type
c     local
      real *8, allocatable :: lapmats2d(:,:,:), temp1(:,:), dints(:,:,:)
      integer, allocatable :: ip2ind(:,:,:), ind2p(:,:)
      integer, allocatable :: ip2indout(:,:,:)
      integer, allocatable :: ip2ind2d(:,:), ind2p2d(:,:)
      real *8, allocatable :: choose(:,:), tempv(:)
      complex *16 one, zero, zk2, zk2inv, zkpow
      real *8, allocatable :: dlap(:,:)
      data one / (1.0d0,0.0d0) /
      data zero / (0.0d0,0.0d0) /

      real *8 chp, alpha, beta, sign, lapval, di, dj, dk
      
      integer npol, i, j, iii, niter, n, npol2d, ndegout, npolout
      integer niterp, id,jd,kd,il,jl,kl,ilap,ip,kp,jp,indd
      integer nord, io, jo, ko, niterloc, jjj, jint, jdeg
      integer ind2d, i2d, iint, kint, indout

      ndegout = ndeg + 2*nup
      
      call legetens_npol_3d(ndeg,type,npol)
      call legetens_npol_3d(ndegout,type,npolout)
      call legetens_npol_2d(ndeg,type,npol2d)

      n = ndeg + 1
      
      allocate(ip2ind(0:ndeg,0:ndeg,0:ndeg),ind2p(3,npol))
      allocate(ip2indout(0:ndegout,0:ndegout,0:ndegout))
      allocate(ip2ind2d(0:ndeg,0:ndeg),ind2p2d(2,npol2d))
      
      call legetens_pow2ind_3d(ndeg,type,ip2ind)
      call legetens_ind2pow_3d(ndeg,type,ind2p)
      call legetens_pow2ind_3d(ndegout,type,ip2indout)
      call legetens_pow2ind_2d(ndeg,type,ip2ind2d)
      call legetens_ind2pow_2d(ndeg,type,ind2p2d)

      niter = ndeg/2

c

      zk2 = zk**2
      zk2inv = one/zk2
      


c     power up to get necessary mutliples of lap2d and Ixx matrix
      
      allocate(lapmats2d(npol2d,npol2d,0:niter))
      allocate(temp1(npol2d,npol2d))

      do i = 1,npol2d
         do j = 1,npol2d
            temp1(j,i) = 0.0d0
         enddo
      enddo
            
      
      call legetens_lapmat_2d(ndeg,ndeg,type,temp1,npol2d)
      
      do i = 1,npol2d
         do j = 1,npol2d
            lapmats2d(j,i,0) = 0.0d0
            if (j .eq. i) lapmats2d(j,i,0) = 1.0d0
            if (niter .gt. 0) lapmats2d(j,i,1) = temp1(j,i)
         enddo
      enddo

      alpha = 1.0d0
      beta = 0.0d0
      do iii = 2,niter
         call dgemm('N','N',npol2d,npol2d,npol2d,alpha,temp1,npol2d,
     1       lapmats2d(1,1,iii-1),npol2d,beta,lapmats2d(1,1,iii),npol2d)
      enddo
      
      allocate(dints(0:ndegout+2*niter,0:ndeg,0:(niter+nup)),
     1     tempv(0:ndegout+2*niter))

      do i = 0,ndeg
         do j = 0,ndegout+2*niter
            dints(j,i,0) = 0.0d0
            if (j .eq. i) dints(j,i,0) = 1.0d0
         enddo
      enddo
      
      do i = 1,(nup+niter)
         do j = 0,ndeg
            jdeg = j + 2*(i-1)
            call legeinte_rect0(dints(0,j,i-1),jdeg,tempv)
            jdeg = jdeg + 1
            call legeinte_rect0(tempv,jdeg,dints(0,j,i))
         enddo
      enddo

C     compute choose(i,j) to be [(i-1) choose (j-1)].

      nord = nup + niter
      allocate(choose(nord+1,nord+1))

      choose(1,1) = 1.0d0
      do i = 2,nord+1
         choose(1,i) = 0.0d0
         choose(i,1) = 1.0d0
      enddo
      do i = 2, nord + 1
         do j = 2, nord + 1
            choose(i,j) = choose(i-1,j-1) + choose(i-1,j)
         enddo
      enddo
      
c     combine terms using expanded expressions for
c     lap^(-1)


      do i = 1,npol

         do j = 1,npolout
            ahmat(j,i) = zero
         enddo
         
         ip = ind2p(1,i)
         jp = ind2p(2,i)
         kp = ind2p(3,i)


         if (ip .ge. jp .and. ip .ge. kp) then
c     x pow is greatest
            
            niterloc = max(jp,kp)/2
            ind2d = ip2ind2d(jp,kp)
            
            zkpow = one
            do iii = 1,nup
               
               sign = 1.0d0
               
               do jjj = 0,niter

                  chp = choose(iii+jjj-1+1,iii-1+1)
                  do i2d = 1,npol2d
                     lapval = lapmats2d(i2d,ind2d,jjj)
                     jo = ind2p2d(1,i2d)
                     ko = ind2p2d(2,i2d)
                     if (jo + ko .gt. jp + kp - 2*jjj) cycle
                     do iint = 0,ip+2*iii+2*jjj
                        di = dints(iint,ip,iii+jjj)
                        indout = ip2indout(iint,jo,ko)
                        ahmat(indout,i) = ahmat(indout,i) +
     1                       di*lapval*zkpow*sign*chp
                     enddo
                  enddo
                  sign = -sign
               enddo

               zkpow = -zkpow*zk2

            enddo

         else if (jp .ge. ip .and. jp .ge. kp) then
c     y pow is greatest

            niterloc = max(ip,kp)/2
            ind2d = ip2ind2d(ip,kp)
            
            zkpow = one
            do iii = 1,nup
               
               sign = 1.0d0
               
               do jjj = 0,niter

                  chp = choose(iii+jjj-1+1,iii-1+1)
                  do i2d = 1,npol2d
                     
                     lapval = lapmats2d(i2d,ind2d,jjj)
                     io = ind2p2d(1,i2d)
                     ko = ind2p2d(2,i2d)
                     if (io + ko .gt. ip + kp - 2*jjj) cycle
                     
                     do jint = 0,jp+2*iii+2*jjj
                        dj = dints(jint,jp,iii+jjj)
                        indout = ip2indout(io,jint,ko)
                        
                        ahmat(indout,i) = ahmat(indout,i) +
     1                       dj*lapval*zkpow*sign*chp
                     enddo
                  enddo
                  sign = -sign
               enddo

               zkpow = -zkpow*zk2

            enddo
            
            
         else
c     z pow is greatest

            niterloc = max(ip,jp)/2
            ind2d = ip2ind2d(ip,jp)
            zkpow = one
            do iii = 1,nup
               
               sign = 1.0d0
               
               do jjj = 0,niter

                  chp = choose(iii+jjj-1+1,iii-1+1)
                  do i2d = 1,npol2d
                     lapval = lapmats2d(i2d,ind2d,jjj)
                     io = ind2p2d(1,i2d)
                     jo = ind2p2d(2,i2d)
                     if (io + jo .gt. ip + jp - 2*jjj) then
                        cycle
                     endif
                     do kint = 0,kp+2*iii+2*jjj
                        dk = dints(kint,kp,iii+jjj)
                        indout = ip2indout(io,jo,kint)
                        
                        ahmat(indout,i) = ahmat(indout,i) +
     1                       dk*lapval*zkpow*sign*chp
                     enddo
                  enddo
                  sign = -sign
               enddo

               zkpow = -zkpow*zk2

            enddo
            
            
         endif

      enddo
            
            
      return
      end
c
c
c
c     
      

      
      subroutine h3danti_legeupfast(ndeg,nup,type,zk,ahmat,ldahmat)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Apply the downward recurrence formula to obtain the
c     anti-Helmholtzian on Legendre polynomials
c
c     input:
c
c     ndeg - integer
c        the degree of the basis
c     nup - the number of upward iterations to run 
c     type - character
c        the type of the basis (for now, total degree
c        is the only implemented option)
c        type = 't', total degree
c        type = 'f', full tensor
c     zk - complex double
c        Helmholtz parameter
c     ldahmat - integer
c        leading dimension of ahmat
c     output:
c     
c     ahmat - complex array (ldahmat,*)
c        column i gives the coefficients of the anti-Helmholtzian of
c        the ith function in the polynomial basis determined by ndeg
c        and type. The coefficients are in that same basis for the
c        downward recurrence, i.e. this function gives an anti-
c        Helmholtzian that's of the same degree.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer ndeg, ldahmat, nup
      complex *16 ahmat(ldahmat,*), zk
      character type
c     local
      real *8, allocatable :: dxxs(:,:,:), temp1(:,:), dints(:,:,:)
      integer, allocatable :: ip2ind(:,:,:), ind2p(:,:)
      integer, allocatable :: ip2indout(:,:,:)
      real *8, allocatable :: choose(:,:), tempv(:)
      complex *16 one, zero, zk2, zk2inv, zkpow
      real *8, allocatable :: dlap(:,:)
      data one / (1.0d0,0.0d0) /
      data zero / (0.0d0,0.0d0) /

      real *8 chp, alpha, beta, sign, lapval, di, dj
      real *8 dk, dxxk, dxxj, dxxi, chdxx
      
      integer npol, i, j, iii, niter, n, npol1d, ndegout, npolout
      integer niterp, id,jd,kd,il,jl,kl,ilap,ip,kp,jp,indd
      integer nord, io, jo, ko, niterloc, jjj, jint, jdeg
      integer iint, kint, indout, ni, nj, nk

      ndegout = ndeg + 2*nup
      
      call legetens_npol_3d(ndeg,type,npol)
      call legetens_npol_3d(ndegout,type,npolout)

      npol1d = ndeg + 1

      n = ndeg + 1
      
      allocate(ip2ind(0:ndeg,0:ndeg,0:ndeg),ind2p(3,npol))
      allocate(ip2indout(0:ndegout,0:ndegout,0:ndegout))
      
      call legetens_pow2ind_3d(ndeg,type,ip2ind)
      call legetens_ind2pow_3d(ndeg,type,ind2p)
      call legetens_pow2ind_3d(ndegout,type,ip2indout)

      niter = ndeg/2

c

      zk2 = zk**2
      zk2inv = one/zk2
      


c     power up to get necessary mutliples of dxx and Ixx matrix
      
      allocate(dxxs(0:ndeg,0:ndeg,0:niter))
      allocate(temp1(0:ndeg,0:ndeg))

      do i = 0,ndeg
         do j = 0,ndeg
            temp1(j,i) = 0.0d0
         enddo
      enddo
            
      
      call legecoeff_d2mat(ndeg,temp1,npol1d)
      
      do i = 0,ndeg
         do j = 0,ndeg
            dxxs(j,i,0) = 0.0d0
            if (j .eq. i) dxxs(j,i,0) = 1.0d0
            if (niter .gt. 0) dxxs(j,i,1) = temp1(j,i)
         enddo
      enddo

      alpha = 1.0d0
      beta = 0.0d0
      do iii = 2,niter
         call dgemm('N','N',npol1d,npol1d,npol1d,alpha,temp1,npol1d,
     1       dxxs(0,0,iii-1),npol1d,beta,dxxs(0,0,iii),npol1d)
      enddo
      
      allocate(dints(0:ndegout+2*niter,0:ndeg,0:(niter+nup)),
     1     tempv(0:ndegout+2*niter))

      do i = 0,ndeg
         do j = 0,ndegout+2*niter
            dints(j,i,0) = 0.0d0
            if (j .eq. i) dints(j,i,0) = 1.0d0
         enddo
      enddo
      
      do i = 1,(nup+niter)
         do j = 0,ndeg
            jdeg = j + 2*(i-1)
            call legeinte_rect0(dints(0,j,i-1),jdeg,tempv)
            jdeg = jdeg + 1
            call legeinte_rect0(tempv,jdeg,dints(0,j,i))
         enddo
      enddo

C     compute choose(i,j) to be [(i-1) choose (j-1)].

      nord = nup + niter
      allocate(choose(nord+1,nord+1))

      choose(1,1) = 1.0d0
      do i = 2,nord+1
         choose(1,i) = 0.0d0
         choose(i,1) = 1.0d0
      enddo
      do i = 2, nord + 1
         do j = 2, nord + 1
            choose(i,j) = choose(i-1,j-1) + choose(i-1,j)
         enddo
      enddo
      
c     combine terms using expanded expressions for
c     lap^(-1)


      do i = 1,npol

         do j = 1,npolout
            ahmat(j,i) = zero
         enddo
         
         ip = ind2p(1,i)
         jp = ind2p(2,i)
         kp = ind2p(3,i)


         if (ip .ge. jp .and. ip .ge. kp) then
c     x pow is greatest
            
            zkpow = one
            do iii = 1,nup
               
               sign = 1.0d0
               
               do jjj = 0,niter

                  chp = choose(iii+jjj-1+1,iii-1+1)

                  do nj = 0,jjj
                     nk = jjj-nj
                     chdxx = choose(jjj+1,nk+1)
c                     if (nj .gt. jp/2 .or. nk .gt. kp/2) cycle
                     do jo = 0,jp-2*nj
                        dxxj = dxxs(jo,jp,nj)
                        do ko = 0,kp-2*nk
                           dxxk = dxxs(ko,kp,nk)
                           
                           do iint = 0,ip+2*iii+2*jjj
                              di = dints(iint,ip,iii+jjj)
                              indout = ip2indout(iint,jo,ko)
                              
                              ahmat(indout,i) = ahmat(indout,i) +
     1                             di*dxxj*dxxk*zkpow*sign*chp*chdxx
                           enddo
                        enddo
                     enddo
                  enddo
                  sign = -sign
               enddo
               
               zkpow = -zkpow*zk2

            enddo

         else if (jp .ge. ip .and. jp .ge. kp) then
c     y pow is greatest

            zkpow = one
            do iii = 1,nup
               
               sign = 1.0d0
               
               do jjj = 0,niter

                  chp = choose(iii+jjj-1+1,iii-1+1)
                  do ni = 0,jjj
                     nk = jjj-ni
                     if (ni .gt. ip/2 .or. nk .gt. kp/2) cycle
                     chdxx = choose(jjj+1,nk+1)
                     do io = 0,ip-2*ni
                        dxxi = dxxs(io,ip,ni)
                        do ko = 0,kp-2*nk
                           dxxk = dxxs(ko,kp,nk)
                           
                           do jint = 0,jp+2*iii+2*jjj
                              dj = dints(jint,jp,iii+jjj)
                              indout = ip2indout(io,jint,ko)
                              ahmat(indout,i) = ahmat(indout,i) +
     1                             dj*dxxi*dxxk*zkpow*sign*chp*chdxx
                           enddo
                        enddo
                     enddo
                  enddo
                  sign = -sign
               enddo

               zkpow = -zkpow*zk2

            enddo
            
            
         else
c     z pow is greatest

            zkpow = one
            do iii = 1,nup
               
               sign = 1.0d0
               
               do jjj = 0,niter

                  chp = choose(iii+jjj-1+1,iii-1+1)
                  do ni = 0,jjj
                     nj = jjj-ni
                     if (ni .gt. ip/2 .or. nj .gt. jp/2) cycle
                     chdxx = choose(jjj+1,nj+1)
                     do io = 0,ip-2*ni
                        dxxi = dxxs(io,ip,ni)
                        do jo = 0,jp-2*nj
                           dxxj = dxxs(jo,jp,nj)
                           
                           do kint = 0,kp+2*iii+2*jjj
                              dk = dints(kint,kp,iii+jjj)
                              indout = ip2indout(io,jo,kint)
                              ahmat(indout,i) = ahmat(indout,i) +
     1                             dk*dxxi*dxxj*zkpow*sign*chp*chdxx
                           enddo
                        enddo
                     enddo
                  enddo
                  sign = -sign
               enddo

               zkpow = -zkpow*zk2

            enddo
            
            
         endif

      enddo
            
            
      return
      end
c
c
c
c     
      

      
