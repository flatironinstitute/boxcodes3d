
      subroutine eaton(r,a,b,val)
      implicit none
      real *8 r, a, b, val
c     local
      complex *16 cmat(4,4), w(27*4+50), one, zero, roots(4), ima
      real *8 rsmall, rri, rii
      data rsmall / 1.0d-15 /
      data one / (1.0d0,0.0d0) /
      data zero / (0.0d0,0.0d0) /
      data ima / (0.0d0,1.0d0) /      
      integer i, ier, n4, nroots, j
      
      val = b
      if (r .lt. rsmall) return

      val = 0.0d0
      if (r .gt. a) return
      
      do i = 1,4
         do j = 1,4
            cmat(j,i) = zero
         enddo
      enddo
      
      do i = 1,3
         cmat(i+1,i) = one
      enddo

c     companion matrix for y^4 - 2*a*y/r + 1
      
      cmat(1,4) = -one
      cmat(2,4) = 2.0d0*a/r

      n4 = 4
      ier = 0

c      call prin2('cmat *',cmat,32)
      
      call nonsym_eigvals(ier,cmat,n4,roots,nroots,w)

      if (ier .ne. 0) then
         write(*,*) 'ier not zero '
      endif
c     all prinf('ier *',ier,1)
c      call prinf('ier nonsym_eigvals *',ier,1)
c     call prin2('roots *',roots,2*n4)

      val = 0.0d0
      do i = 1,n4
         rri = real(roots(i))
         rii = abs(real(ima*roots(i)))
         if (rri .gt. 1 .and. rii .lt. 1d-7) then
            val = rri**2 - 1.0d0
         endif
      enddo

      if (val .gt. b) val = b

      return
      end
      
      

      subroutine eatonprol_form(ier,a,b,c,r0,rend,w,lenw,lused,clege,
     1     nlege)
      implicit real *8 (a-h,o-z)
      real *8 a,b,c,w(*),r0,r1,clege(*)
      integer lenw
c
c     form the smoothed eaton lens evaluator in the prolate
c     basis
c
c     NOTE: we recommend a small buffer so that
c     [0,rmax] is contained in [r0,rend] (gets rid of some
c     gibbsing near endpoints)
c      
c
c     input
c
c     a - support of eaton lens
c     b - maximum of eaton lens (cut off)
c     c - bandlimit of prolate basis
c     r0 - lower bound of fitting region (not necessarily lower
c                 bound of where it will be used)
c     rend - upper bound of fitting region (not necessarily
c                 upper bound of where it will eventually
c                 be evaluated)
c     lenw - length of array w
c     w - work array storing everything needed to evaluate
c       prolate approximation of lens in [r0,rend]
c     lused - number of entries in w that are used
c
c     output
c
c     clege - legendre expansion of smoothed function
c     nlege - order of legendre expansion
c     ier - if ier = 4 lenw insufficient
c     

      real *8, allocatable :: wlst(:),amat(:,:),rhs(:),rnorms(:)
      
      ier = 0
      if (lenw .lt. 30) then
         ier = 4
         return
      endif

      w(1) = a
      w(2) = b
      w(3) = c
      w(4) = r0
      w(5) = rend

      nlam = 200 + c*2
      idkh = 21
      ilam1 = idkh + nlam
      ilam2 = ilam1 + nlam*2
      iw1 = ilam2 + nlam
      lenw1 = lenw - iw1 + 1

      w(6) = idkh + 0.1d0
      w(7) = ilam1 + 0.1d0
      w(8) = ilam2 + 0.1d0
      w(9) = iw1 + 0.1d0

      call prolcrea(ier1,c,w(iw1),lenw1,nvects,nhigh,w(idkh),
     1     w(ilam1),w(ilam2),keep1,lused1)

      if (lused1 .gt. lenw1) then
         lused = lused1 + iw1
         ier = 4
         return
      endif

      call prinf('eatonprol_form, nvects *',nvects,1)
      call prinf('eatonprol_form, nhigh *',nhigh,1)
      
      nfun = nvects + 1

      w(10) = nvects + 0.1d0
      w(11) = nfun + 0.1d0
      
      icoefs = iw1 + keep1 + 10

      keep = icoefs + nfun

      lused = max(lused,keep)
      
      w(12) = icoefs + 0.1d0

      if (keep .gt. lenw .or. lused .gt. lenw) then
         ier = 4
         return
      endif

      nsamp = 4*nfun + 400
      
      allocate(rhs(nsamp),amat(nsamp,nfun),
     1     wlst(nsamp*nfun*8 + 500),
     2     rnorms(nsamp+nfun))

      pa = -1.0d0
      pb = 1.0d0
      ph = (pb-pa)/(nsamp-1)

      rh = (rend-r0)/(nsamp-1)
      
      do i = 1,nsamp
         r = r0 + rh*(i-1)
         t = pa + ph*(i-1)

         call eaton(r,a,b,rhs(i))
         
         do j = 0,nvects
            call proleval(j,t,w(iw1),val,der)
            amat(i,j+1) = val
         enddo
      enddo

      eps = 1d-15
      call nrleastsq(amat,nsamp,nfun,eps,ncolsout,rnorms,wlst)
      call nrleasts2(wlst,rhs,w(icoefs))

      call protoleg(w(icoefs),nfun,w(iw1),clege,nlege)

      call prinf('eatonprol_form, nlege *',nlege,1)

      
      return
      end

      subroutine eatonprol_eval(rs,nr,clege,nlege,r0,rend,vals)
      implicit real *8 (a-h,o-z)
      real *8 rs(*), vals(*), clege(*)
c
c     must follow a call to eatonprol_form
c
c     evaluate the smoothed representation of the
c     eaton lens in the prolate basis at the specified
c     points
c
c     input:
c
c     rs - array of r values
c     nr - number of r values
c     clege - legendre coefficients
c     nlege - order of legendre coefs
c     r0, rend - endpoints used in call to eatonprol_form
c
c     output:
c
c     vals - array of values of smoothed function at given
c         points
c     


      do i = 1,nr
         r = rs(i)
         t = (r-r0)*2.0d0/(rend-r0)-1.0d0
         call legeexev(t,vals(i),clege,nlege)
      enddo

      return
      end

      subroutine eatonprol_eval0(rs,nr,w1,r0,rend,nvects,coefs,vals)
      implicit real *8 (a-h,o-z)
      real *8 rs(*), w1(*), coefs(*), vals(*)

      do i = 1,nr
         r = rs(i)
         t = (r-r0)*2.0d0/(rend-r0)-1.0d0

         vals(i) = 0.0d0
         do j = 0,nvects
            call proleval(j,t,w1,val,der)
            vals(i) = vals(i) + val*coefs(j+1)
         enddo
      enddo

      return
      end

      
