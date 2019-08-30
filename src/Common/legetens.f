c
c
c
c     this file contains routines for working with tensor
c     product gauss-legendre nodes and polynomials in 2 and
c     3 dimensions; also provides some related spectral
c     differentiation operators
c
c--------------------------------------------------------------
c
c     CONVENTIONS
c
c     ndeg refers to polynomial degree
c      
c     n refers to order of approximation and number of nodes
c     used in each direction (n = ndeg+1)
c
c     type character specifies the set of tensor
c     polynomials to use. current options
c
c     type = 'F', full degree polynomials (T_i(x)T_j(y)T_k(z) where
c                 each of i, j, and k goes from 0 to ndeg)
c     type = 'T', total degree polynomials (T_i(x)T_j(y)T_k(z) where
c                 sum of i,j,k>=0 is less than or equal to ndeg)
c
c     TODO: implement "Euclidean degree"
c
c     a tensor grid of points is traversed with x on the inner
c        loop, y on the loop above that, and z above that.
c        e.g. for a 2x2x2 grid, we have the order:
c            (x1,y1,z1), (x2,y1,z1), (x1,y2,z1), (x2,y2,z1)
c     (x1,y1,z2), (x2,y1,z2), (x1,y2,z2), (x2,y2,z2)
c      
c     tensor product polynomials are numbered analagously
c        e.g. for full degree at ndeg = 3, we would have the order
c        1, T_1(x), T_2(x), T_1(y), T_1(x)T_1(y), T_2(x)T_1(y),
c        T_2(y), T_1(x)T_2(y), T_2(x)T_2(y), T_1(z), T_1(x)T_1(z),
c        T_2(x) T_1(z), T_1(y)T_1(z), T_1(x)T_1(y)T_1(z), T_2(x)
c        T_1(y)T_1(z), T_2(y) T_1(z), T_1(x)T_2(y) T_1(z),
c        T_2(x) T_2(y) T_1(z), T_2(z), T_1(x)T_2(z), T_2(x) T_2(z),
c        T_1(y)T_2(z), T_1(x)T_1(y)T_2(z), T_2(x) T_1(y) T_2(z),
c        T_2(y) T_2(z), T_1(x)T_2(y) T_2(z), T_2(x) T_2(y) T_2(z)
c 
c        e.g. for total degree, we would have the order
c        1, T_1(x), T_2(x), T_1(y), T_1(x)T_1(y), T_2(y),
c        T_1(z), T_1(x)T_1(z), T_1(y)T_1(z), T_2(z)
c
c-------------------------------------------------------------     
c
c     ROUTINES
c      
c     legetens_npol_*d - total number of polynomials up
c                        to a given degree of specified
c                        type in 2 or 3 dimensions
c     legetens_exps_*d - get nodes, weights, and, optionally,
c                        projection and evaluation matrices      
c                        in 2 or 3 dimensions
c     legecoeff_dmat - rectangular 1d spectral diff mat
c     legecoeff_d2mat - rectangular 1d spectral 2nd deriv mat
c     legetens_pols_*d - evaluate the tensor legendre polys at
c                        a given point      
c     legetens_lapmat_*d - rectangular spectral laplacian in
c                          2 and 3 dims
c     legetens_eyemat_*d - rectangular "identity" matrix between
c                          different polynomial orders in
c                          2 and 3 dims 

      subroutine legetens_npol_2d(ndeg,type,npol)
c
c     return the number of polynomials of a given type
c     up to degree ndeg
c      
      integer n, npol, ndeg
      character type

      n = ndeg + 1
      
      if (type .eq. 'f' .or. type .eq. 'F') then
         npol = n**2
      else if (type .eq. 't' .or. type .eq. 'T') then
         npol = n*(n+1)/2
      endif

      return
      end

c
      subroutine legetens_npol_3d(ndeg,type,npol)
c
c     return the number of polynomials of a given type
c     up to degree ndeg
c      
      integer n, npol, ndeg
      character type

      n = ndeg + 1
      
      if (type .eq. 'f' .or. type .eq. 'F') then
         npol = n**3
      else if (type .eq. 't' .or. type .eq. 'T') then
         npol = n*(n+1)*(n+2)/6
      endif

      return
      end


      subroutine legetens_pow2ind_3d(ndeg,type,ip2ind)
      integer ndeg, ip2ind(ndeg+1,ndeg+1,ndeg+1)
      character type

      integer i, j, k, ipol

      do i = 1,ndeg+1
         do j = 1,ndeg+1
            do k = 1,ndeg+1
               ip2ind(k,j,i) = -1
            enddo
         enddo
      enddo

      if (type .eq. 'f' .or. type .eq. 'F') then
         ipol = 0
         do i = 1,ndeg+1
            do j = 1,ndeg+1
               do k = 1,ndeg+1
                  ipol = ipol+1
                  ip2ind(k,j,i) = ipol
               enddo
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         ipol = 0
         do i = 1,ndeg+1
            do j = 1,ndeg+1+1-i
               do k = 1,ndeg+1+2-i-j
                  ipol = ipol+1
                  ip2ind(k,j,i) = ipol
               enddo
            enddo
         enddo
      endif

      return
      end
c
      subroutine legetens_ind2pow_3d(ndeg,type,iind2p)
      integer ndeg, iind2p(3,*)
      character type

      integer i, j, k, ipol


      if (type .eq. 'f' .or. type .eq. 'F') then
         ipol = 0
         do i = 1,ndeg+1
            do j = 1,ndeg+1
               do k = 1,ndeg+1
                  ipol = ipol+1
                  iind2p(1,ipol) = k-1
                  iind2p(2,ipol) = j-1
                  iind2p(3,ipol) = i-1                
               enddo
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         ipol = 0
         do i = 1,ndeg+1
            do j = 1,ndeg+1+1-i
               do k = 1,ndeg+1+2-i-j
                  ipol = ipol+1
                  iind2p(1,ipol) = k-1
                  iind2p(2,ipol) = j-1
                  iind2p(3,ipol) = i-1
               enddo
            enddo
         enddo
      endif

      return
      end
c
c
      subroutine legetens_exps_2d(itype,n,type,x,u,ldu,v,ldv,w)
c                 input parameters:
c
c  itype - the type of the calculation to be performed
c          itype=0 means that only the gaussian nodes are 
c                  to be constructed. 
c          itype=1 means that only the nodes and the weights 
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c          itype=3 only construct u
c          itype=4 only construct v
      
      implicit none
      integer itype, n, ldu, ldv
      character type
      real *8 x(2,*),w(*)
      real *8 u(ldu,*), v(ldv,*)
      real *8 x1d(n), w1d(n), u1d(n,n), v1d(n,n)
      integer i,j,ipt,itype1d,io, jo,ipol
      
      itype1d = 0
      if (itype .ge. 1) then
         itype1d = 1
      endif
      if (itype .ge. 2) then
         itype1d = 2
      endif
      
      call legeexps(itype1d,n,x1d,u1d,v1d,w1d)

      ipt = 0
      do i=1,n
         do j=1,n
            ipt = ipt + 1
            x(1,ipt) = x1d(j)
            x(2,ipt) = x1d(i)               
        enddo
      enddo

      if (itype .ge. 1) then
         ipt = 0
         do i=1,n
            do j=1,n
               ipt = ipt + 1
               w(ipt) = w1d(i)*w1d(j)
            enddo
         enddo
      endif


      if (itype .eq. 2 .or. itype .eq. 3) then
c     construct u from 1d u
         if (type .eq. 'f' .or. type .eq. 'F') then         
            ipt = 0
            do io = 1,n
            do jo = 1,n
               ipt = ipt + 1
               ipol = 0
               do i=1,n
               do j=1,n
                  ipol = ipol + 1
                  u(ipol,ipt) =  u1d(i,io)*u1d(j,jo)
               enddo
               enddo
            enddo
            enddo

         else if (type .eq. 't' .or. type .eq. 'T') then

            ipt = 0
            do io = 1,n
            do jo = 1,n
               ipt = ipt + 1
               ipol = 0
               do i=1,n
               do j=1,n+1-i
                  ipol = ipol + 1
                  u(ipol,ipt) = u1d(i,io)*u1d(j,jo)
               enddo
               enddo
            enddo
            enddo
         endif
      endif
      
      if (itype .eq. 2 .or. itype .eq. 4) then
c     construct v from 1d v
         if (type .eq. 'f' .or. type .eq. 'F') then         
            ipol = 0
            do io = 1,n
            do jo = 1,n
               ipol = ipol + 1
               ipt = 0
               do i=1,n
               do j=1,n
                  ipt = ipt + 1
                  v(ipt,ipol) =  v1d(i,io)*v1d(j,jo)
               enddo
               enddo
            enddo
            enddo

         else if (type .eq. 't' .or. type .eq. 'T') then

            ipol = 0
            do io = 1,n
            do jo = 1,n+1-io
               ipol = ipol + 1
               ipt = 0
               do i=1,n
               do j=1,n
                  ipt = ipt + 1
                  v(ipt,ipol) = v1d(i,io)*v1d(j,jo)
               enddo
               enddo
            enddo
            enddo
         endif
      endif         

      return
      end
c
c
      subroutine legetens_exps_3d(itype,n,type,x,u,ldu,v,ldv,w)
c                 input parameters:
c
c  itype - the type of the calculation to be performed
c          itype=0 means that only the gaussian nodes are 
c                  to be constructed. 
c          itype=1 means that only the nodes and the weights 
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c          itype=3 only construct x,w, and u
c          itype=4 only construct x,w, and v
      
      implicit none
      integer itype, n, ldu, ldv
      character type
      real *8 x(3,*),w(*)
      real *8 u(ldu,*), v(ldv,*)
      real *8 x1d(n), w1d(n), u1d(n,n), v1d(n,n)
      integer i,j,ipt,itype1d, k, io, jo, ko, ipol
      
      itype1d = 0
      if (itype .ge. 1) then
         itype1d = 1
      endif
      if (itype .ge. 2) then
         itype1d = 2
      endif
      
      call legeexps(itype1d,n,x1d,u1d,v1d,w1d)

      ipt = 0
      do i=1,n
         do j=1,n
            do k = 1,n
               ipt = ipt + 1
               x(1,ipt) = x1d(k)
               x(2,ipt) = x1d(j)               
               x(3,ipt) = x1d(i)
            enddo
        enddo
      enddo

      if (itype .ge. 1) then
         ipt = 0
         do i=1,n
            do j=1,n
               do k = 1,n
                  ipt = ipt + 1
                  w(ipt) = w1d(i)*w1d(j)*w1d(k) 
               enddo
            enddo
         enddo
      endif


      if (itype .eq. 2 .or. itype .eq. 3) then
c     construct u from 1d u
         if (type .eq. 'f' .or. type .eq. 'F') then         
            ipt = 0
            do io = 1,n
            do jo = 1,n
            do ko = 1,n
               ipt = ipt + 1
               ipol = 0
               do i=1,n
               do j=1,n
               do k = 1,n
                  ipol = ipol + 1
                  u(ipol,ipt) =  u1d(i,io)*u1d(j,jo)*u1d(k,ko)
               enddo
               enddo
               enddo
            enddo
            enddo
            enddo

         else if (type .eq. 't' .or. type .eq. 'T') then

            ipt = 0
            do io = 1,n
            do jo = 1,n
            do ko = 1,n
               ipt = ipt + 1
               ipol = 0
               do i=1,n
               do j=1,n+1-i
               do k = 1,n+2-i-j
                  ipol = ipol + 1
                  u(ipol,ipt) = u1d(i,io)*u1d(j,jo)*u1d(k,ko)
               enddo
               enddo
               enddo
            enddo
            enddo
            enddo
         endif
      endif
      
      if (itype .eq. 2 .or. itype .eq. 4) then
c     construct v from 1d v
         if (type .eq. 'f' .or. type .eq. 'F') then         
            ipol = 0
            do io = 1,n
            do jo = 1,n
            do ko = 1,n
               ipol = ipol + 1
               ipt = 0
               do i=1,n
               do j=1,n
               do k=1,n
                  ipt = ipt + 1
                  v(ipt,ipol) =  v1d(i,io)*v1d(j,jo)*v1d(k,ko)
               enddo
               enddo
               enddo
            enddo
            enddo
            enddo

         else if (type .eq. 't' .or. type .eq. 'T') then

            ipol = 0
            do io = 1,n
            do jo = 1,n+1-io
            do ko = 1,n+2-io-jo
               ipol = ipol + 1
               ipt = 0
               do i=1,n
               do j=1,n
               do k=1,n
                  ipt = ipt + 1
                  v(ipt,ipol) = v1d(i,io)*v1d(j,jo)*v1d(k,ko)
               enddo
               enddo
               enddo
            enddo
            enddo
            enddo
         endif
      endif         

      return
      end
c
c

      subroutine legetens_pols_2d(x,ndeg,type,pols)
c
c     type = 'F' full degree polynomials ((ndeg+1)**3)
c     type = 'T' total degree polynomials ((ndeg+1)*(ndeg+2)*(ndeg+3)/6)
c      
      implicit none
      integer ndeg,npols
      real *8 x(2),pols(*),px(ndeg+1),py(ndeg+1)
      integer i,j,ipol,n
      character type

      n = ndeg + 1
      
      call legepols(x(1),ndeg,px)
      call legepols(x(2),ndeg,py)

      if (type .eq. 'f' .or. type .eq. 'F') then
         ipol = 0
         do i=1,n
            do j =1,n
               ipol = ipol + 1
               pols(ipol) = px(j)*py(i)
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         ipol = 0
         do i=1,n
            do j=1,n+1-i
               ipol = ipol + 1
               pols(ipol) = px(j)*py(i)
            enddo
         enddo
      endif
      
      return
      end

      subroutine legetens_pols_3d(x,ndeg,type,pols)
c
c     type = 'F' full degree polynomials ((ndeg+1)**3)
c     type = 'T' total degree polynomials ((ndeg+1)*(ndeg+2)*(ndeg+3)/6)
c      
      implicit none
      integer ndeg,npols
      real *8 x(3),pols(*),px(ndeg+1),py(ndeg+1),pz(ndeg+1)
      integer i,j,ipol, k, n
      character type

      n = ndeg + 1
      
      call legepols(x(1),ndeg,px)
      call legepols(x(2),ndeg,py)
      call legepols(x(3),ndeg,pz)      

      if (type .eq. 'f' .or. type .eq. 'F') then
         ipol = 0
         do i=1,n
            do j=1,n
               do k = 1,n
                  ipol = ipol + 1
                  pols(ipol) = px(k)*py(j)*pz(i)
               enddo
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         ipol = 0
         do i=1,n
            do j=1,n+1-i
               do k = 1,n+2-i-j
                  ipol = ipol + 1
                  pols(ipol) = px(k)*py(j)*pz(i)
               enddo
            enddo
         enddo
      endif
      
      return
      end

      
      subroutine legetens_ders_3d(x,ndeg,type,idim,ders)
c
c     type = 'F' full degree polynomials ((ndeg+1)**3)
c     type = 'T' total degree polynomials ((ndeg+1)*(ndeg+2)*(ndeg+3)/6)
c      
      implicit none
      integer ndeg,npols,idim
      real *8 x(3),ders(*),px(ndeg+1),py(ndeg+1),pz(ndeg+1),tmp(ndeg+1)
      integer i,j,ipol, k, n
      character type

      n = ndeg + 1
      
      call legepols(x(1),ndeg,px)
      call legepols(x(2),ndeg,py)
      call legepols(x(3),ndeg,pz)

      if (idim .eq. 1) then
         call legepolders(x(1),tmp,px,ndeg)
      endif

      if (idim .eq. 2) then
         call legepolders(x(2),tmp,py,ndeg)
      endif

      if (idim .eq. 3) then
         call legepolders(x(3),tmp,pz,ndeg)
      endif

      if (type .eq. 'f' .or. type .eq. 'F') then
         ipol = 0
         do i=1,n
            do j=1,n
               do k = 1,n
                  ipol = ipol + 1
                  ders(ipol) = px(k)*py(j)*pz(i)
               enddo
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         ipol = 0
         do i=1,n
            do j=1,n+1-i
               do k = 1,n+2-i-j
                  ipol = ipol + 1
                  ders(ipol) = px(k)*py(j)*pz(i)
               enddo
            enddo
         enddo
      endif
      
      return
      end

      
      subroutine legecoeff_dmat(ndeg,dmat,lddmat)
c
c     
c      

      implicit none
      integer ndeg, lddmat
      real *8 dmat(lddmat,ndeg+1)
c     local
      real *8 polin(ndeg+1), polout(ndeg+1)
      integer i, j

      do i = 1,ndeg+1
         do j = 1,ndeg+1
            polin(j) = 0.0d0
            if (j .eq. i) polin(j) = 1.0d0
         enddo
         call legediff(polin,ndeg,polout)
         do j = 1,ndeg
            dmat(j,i) = polout(j)
         enddo
      enddo

      return
      end
      

      subroutine legecoeff_d2mat(ndeg,dmat,lddmat)
      implicit none
      integer ndeg, lddmat
      real *8 dmat(lddmat,ndeg+1)
c     local
      real *8 polin(ndeg+1), polout(ndeg+1)
      integer i, j

      do i = 1,ndeg+1
         do j = 1,ndeg+1
            polin(j) = 0.0d0
            if (j .eq. i) polin(j) = 1.0d0
         enddo
         call legediff(polin,ndeg,polout)
         call legediff(polout,ndeg-1,polin)         
         do j = 1,ndeg-1
            dmat(j,i) = polin(j)
         enddo
      enddo

      return
      end

      subroutine legetens_lapmat_2d(ndegin,ndegout,type,
     1     lapmat,ldlapmat)
      implicit none
      integer ndegin, ndegout, ldlapmat
      character type
      real *8 lapmat(ldlapmat,*)
c     local
      real *8 lap1d(ndegout+1,ndegin+1)
      integer ldlap1d, io, ii, jo, ji, iopol, iipol

      do ii = 1,ndegin+1
         do io = 1,ndegout+1
            lap1d(io,ii) = 0.0d0
         enddo
      enddo
      
      
      ldlap1d = ndegout+1

      call legecoeff_d2mat(ndegin,lap1d,ldlap1d)

      if (type .eq. 'f' .or. type .eq. 'F') then
         iipol = 0
         do ii = 1,ndegin+1
            do ji = 1,ndegin+1
               iipol = iipol + 1
               iopol = 0
               do io = 1,ndegout+1
                  do jo = 1,ndegout+1
                     iopol = iopol + 1
                     lapmat(iopol,iipol) = 0.0d0
                     if (io .eq. ii) then
                        lapmat(iopol,iipol) =
     1                       lapmat(iopol,iipol) + lap1d(jo,ji)
                     endif
                     if (jo .eq. ji) then
                        lapmat(iopol,iipol) =
     1                       lapmat(iopol,iipol) + lap1d(io,ii)
                     endif
                  enddo
               enddo
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         iipol = 0
         do ii = 1,ndegin+1
            do ji = 1,ndegin+1+1-ii
               iipol = iipol + 1
               iopol = 0
               do io = 1,ndegout+1
                  do jo = 1,ndegout+1+1-io
                     iopol = iopol + 1
                     lapmat(iopol,iipol) = 0.0d0
                     if (io .eq. ii) then
                        lapmat(iopol,iipol) =
     1                       lapmat(iopol,iipol) + lap1d(jo,ji)
                     endif
                     if (jo .eq. ji) then
                        lapmat(iopol,iipol) =
     1                       lapmat(iopol,iipol) + lap1d(io,ii)
                     endif
                     
                  enddo
               enddo
            enddo
         enddo
      endif
      
      return
      end

      subroutine legetens_lapmat_3d(ndegin,ndegout,type,
     1     lapmat,ldlapmat)
      implicit none
      integer ndegin, ndegout, ldlapmat
      character type
      real *8 lapmat(ldlapmat,*)
c     local
      real *8 lap1d(ndegout+1,ndegin+1)
      integer ldlap1d, io, ii, jo, ji, ko, ki, iopol, iipol

      ldlap1d = ndegout+1

      do ii = 1,ndegin+1
         do io = 1,ndegout+1
            lap1d(io,ii) = 0.0d0
         enddo
      enddo
      
      call legecoeff_d2mat(ndegin,lap1d,ldlap1d)

      if (type .eq. 'f' .or. type .eq. 'F') then
         iipol = 0
         do ii = 1,ndegin+1
            do ji = 1,ndegin+1
               do ki = 1,ndegin+1
                  iipol = iipol + 1
                  iopol = 0
                  do io = 1,ndegout+1
                     do jo = 1,ndegout+1
                        do ko = 1,ndegout+1
                           iopol = iopol + 1
                           lapmat(iopol,iipol) = 0.0d0
                           if (jo .eq. ji .and. io .eq. ii) then
                              lapmat(iopol,iipol) =
     1                             lapmat(iopol,iipol) + lap1d(ko,ki)
                           endif
                           if (ko .eq. ki .and. io .eq. ii) then
                              lapmat(iopol,iipol) =
     1                             lapmat(iopol,iipol) + lap1d(jo,ji)
                           endif
                           if (jo .eq. ji .and. ko .eq. ki) then
                              lapmat(iopol,iipol) =
     1                             lapmat(iopol,iipol) + lap1d(io,ii)
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         iipol = 0
         do ii = 1,ndegin+1
            do ji = 1,ndegin+1+1-ii
               do ki = 1,ndegin+1+2-ii-ji
                  iipol = iipol + 1
                  iopol = 0
                  do io = 1,ndegout+1
                     do jo = 1,ndegout+1+1-io
                        do ko = 1,ndegout+1+2-io-jo
                           iopol = iopol + 1
                           lapmat(iopol,iipol) = 0.0d0
                           if (jo .eq. ji .and. io .eq. ii) then
                              lapmat(iopol,iipol) =
     1                             lapmat(iopol,iipol) + lap1d(ko,ki)
                           endif
                           if (ko .eq. ki .and. io .eq. ii) then
                              lapmat(iopol,iipol) =
     1                             lapmat(iopol,iipol) + lap1d(jo,ji)
                           endif
                           if (jo .eq. ji .and. ko .eq. ki) then
                              lapmat(iopol,iipol) =
     1                             lapmat(iopol,iipol) + lap1d(io,ii)
                           endif

                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
      
      return
      end


      subroutine legetens_eyemat_2d(ndegin,ndegout,type,
     1     eyemat,ldeyemat)
      implicit none
      integer ndegin, ndegout, ldeyemat
      character type
      real *8 eyemat(ldeyemat,*)
c     local
      integer io, ii, jo, ji, iopol, iipol

      if (type .eq. 'f' .or. type .eq. 'F') then
         iipol = 0
         do ii = 1,ndegin+1
            do ji = 1,ndegin+1
               iipol = iipol + 1
               iopol = 0
               do io = 1,ndegout+1
                  do jo = 1,ndegout+1
                     iopol = iopol + 1
                     eyemat(iopol,iipol) = 0.0d0
                     if (jo .eq. ji .and. io .eq. ii) then
                        eyemat(iopol,iipol) = 1.0d0
                     endif
                  enddo
               enddo
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         iipol = 0
         do ii = 1,ndegin+1
            do ji = 1,ndegin+1+1-ii
               iipol = iipol + 1
               iopol = 0
               do io = 1,ndegout+1
                  do jo = 1,ndegout+1+1-io
                     iopol = iopol + 1
                     eyemat(iopol,iipol) = 0.0d0
                     if (jo .eq. ji .and. io .eq. ii) then
                        eyemat(iopol,iipol) = 1.0d0
                     endif
                  enddo
               enddo
            enddo
         enddo
      endif

      return
      end
      
      subroutine legetens_eyemat_3d(ndegin,ndegout,type,
     1     eyemat,ldeyemat)
      implicit none
      integer ndegin, ndegout, ldeyemat
      character type
      real *8 eyemat(ldeyemat,*)
c     local
      integer io, ii, jo, ji, ko, ki, iopol, iipol

      if (type .eq. 'f' .or. type .eq. 'F') then
         iipol = 0
         do ii = 1,ndegin+1
            do ji = 1,ndegin+1
               do ki = 1,ndegin+1
                  iipol = iipol + 1
                  iopol = 0
                  do io = 1,ndegout+1
                     do jo = 1,ndegout+1
                        do ko = 1,ndegout+1
                           iopol = iopol + 1
                           eyemat(iopol,iipol) = 0.0d0
                           if (jo .eq. ji .and. io .eq. ii .and.
     1                          ko .eq. ki) then
                              eyemat(iopol,iipol) = 1.0d0
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         iipol = 0
         do ii = 1,ndegin+1
            do ji = 1,ndegin+1+1-ii
               do ki = 1,ndegin+1+2-ii-ji
                  iipol = iipol + 1
                  iopol = 0
                  do io = 1,ndegout+1
                     do jo = 1,ndegout+1+1-io
                        do ko = 1,ndegout+1+2-io-jo
                           iopol = iopol + 1
                           eyemat(iopol,iipol) = 0.0d0
                           if (jo .eq. ji .and. io .eq. ii .and.
     1                          ko .eq. ki) then
                              eyemat(iopol,iipol) = 1.0d0
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
      
      return
      end
      

      subroutine legetens_slicezcoeffs_3d(ndeg,type,idim,val,cin,
     1     ldcin,nfun,cout,ldcout,ifdiff,dcout,lddcout)
c
c     complex valued coeffs, for real see legetens_slicedcoeffs...
c      
c     given a set of polynomials defined by the coefficients in
c     cin, get the 2d coefficients of the slice of that polynomial
c     where x_idim = val
c     
      implicit none
      integer ndeg, idim, ldcin, nfun, ldcout, ifdiff, lddcout
      character type
      real *8 val
      complex *16 cin(ldcin,*), cout(ldcout,*), dcout(lddcout,*)
c     local
      integer iipol, iopol, ii, ji, ki, io, jo, ko, n, npol2, jpol
      integer ilist(ndeg+1,ndeg+1), ifun
      real *8 pols(ndeg+1), ders(ndeg+1)
      complex *16 cv, zero, dv
      data zero / (0.0d0,0.0d0) /

      call legepolders(val,pols,ders,ndeg)

      call legetens_npol_2d(ndeg,type,npol2)

      n = ndeg+1

      if (type .eq. 'f' .or. type .eq. 'F') then
         
         iopol = 0
         do io = 1,n
            do jo = 1,n
               iopol = iopol+1
               ilist(jo,io) = iopol
            enddo
         enddo

         do ifun = 1,nfun

            do jpol = 1,npol2
               cout(jpol,ifun) = zero
            enddo

            if (ifdiff .eq. 1) then
               do jpol = 1,npol2
                  dcout(jpol,ifun) = zero
               enddo
            endif
            
            iipol = 0
            do ii = 1,n
               do ji = 1,n
                  do ki = 1,n
                     iipol = iipol+1
                     if (idim .eq. 1) then
                        io = ii
                        jo = ji
                        cv = pols(ki)*cin(iipol,ifun)
                        dv = ders(ki)*cin(iipol,ifun)
                     else if (idim .eq. 2) then
                        io = ii
                        jo = ki
                        cv = pols(ji)*cin(iipol,ifun)
                        dv = ders(ji)*cin(iipol,ifun)                        
                     else if (idim .eq. 3) then
                        io = ji
                        jo = ki
                        cv = pols(ii)*cin(iipol,ifun)
                        dv = ders(ii)*cin(iipol,ifun)
                     endif
                     iopol = ilist(jo,io)
                     cout(iopol,ifun) = cout(iopol,ifun)+cv
                     if (ifdiff .eq. 1) then
                        dcout(iopol,ifun) = dcout(iopol,ifun)+dv
                     endif
                  enddo
               enddo
            enddo
         enddo
         
      endif
      
      if (type .eq. 't' .or. type .eq. 'T') then
         
         iopol = 0
         do io = 1,n
            do jo = 1,n+1-io
               iopol = iopol+1
               ilist(jo,io) = iopol
            enddo
         enddo

         do ifun = 1,nfun

            do jpol = 1,npol2
               cout(jpol,ifun) = zero
            enddo

            if (ifdiff .eq. 1) then
               do jpol = 1,npol2
                  dcout(jpol,ifun) = zero
               enddo
            endif
            
            iipol = 0
            do ii = 1,n
               do ji = 1,n+1-ii
                  do ki = 1,n+2-ii-ji
                     iipol = iipol+1
                     if (idim .eq. 1) then
                        io = ii
                        jo = ji
                        cv = pols(ki)*cin(iipol,ifun)
                        dv = ders(ki)*cin(iipol,ifun)                        
                     else if (idim .eq. 2) then
                        io = ii
                        jo = ki
                        cv = pols(ji)*cin(iipol,ifun)
                        dv = ders(ji)*cin(iipol,ifun)                        
                     else if (idim .eq. 3) then
                        io = ji
                        jo = ki
                        cv = pols(ii)*cin(iipol,ifun)
                        dv = ders(ii)*cin(iipol,ifun)                        
                     endif
                     iopol = ilist(jo,io)
                     cout(iopol,ifun) = cout(iopol,ifun)+cv
                     if (ifdiff .eq. 1) then
                        dcout(iopol,ifun) = dcout(iopol,ifun)+dv
                     endif
                  enddo
               enddo
            enddo
         enddo
         
      endif
      
      return
      end

      
