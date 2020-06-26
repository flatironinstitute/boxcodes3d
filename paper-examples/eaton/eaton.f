
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
      
      

      subroutine eatonprol_form(ier,a,b,c,r0,rend,w,lenw,lused,keep)
      implicit real *8 (a-h,o-z)
      real *8 a,b,c,w(*),r0,r1,rb
      integer lenw
c
c     form the smoothed eaton lens evaluator in the prolate
c     basis
c
c     NOTE: we recommend a small buffer so that
c     [0,rmax] is contained in [r0,rend], with rmax the
c     largest value where it will be evaluated
c     (gets rid of some gibbsing near endpoints)
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
c     lused - number of entries in w that are used
c
c     output
c
c     w - work array storing everything needed to evaluate
c       prolate approximation of lens in [r0,rend], lege
c       interpolant, and other info used by eatonprol_eval
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

      if (lused .gt. lenw) then
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

      ilegc = keep
      w(13) = ilegc + 0.1d0

      keep = ilegc + nhigh
      
      call protoleg(w(icoefs),nfun,w(iw1),w(ilegc),nlege)

      w(14) = nlege + 0.1d0

      call prinf('eatonprol_form, nlege *',nlege,1)


c     find point where it transitions

      nbis = 50

      r1 = r0
      r2 = a

      do i = 1,nbis
         rb = (r1+r2)/2
         call eaton(rb,a,b,val0)
         if (val0 .ge. b) then
            r1 = rb
         else
            r2 = rb
         endif
      enddo

      call prin2('b-val at rb *',b-val0,1)

      w(15) = rb
      
      
      return
      end

      subroutine eatonprol_eval(ifbell,rsupp,rs,nr,w,vals)
      implicit real *8 (a-h,o-z)
      real *8 rs(*), vals(*), w(*)
c
c     must follow a call to eatonprol_form
c
c     evaluate the smoothed representation of the
c     eaton lens in the prolate basis at the specified
c     points
c
c     input:
c
c     ifbell - ifbell = 1, then apply a bell that smoothly
c       transitions from b to the prolate interpolator on
c       [0,rb] and from the interpolator to 0 on [a,rsupp].
c       Note, rb is stored in w and is the max r such that 
c       eaton(r) = b.
c     rsupp - if ifbell = 1, then this is the true numerical
c       support of the output. rsupp > a in call to
c       eatonprol_form (a is in w(1)).
c     rs - array of r values
c     nr - number of r values
c     w - array output of eatonprol_form. 
c
c     output:
c
c     vals - array of values of smoothed function at given
c         points
c     

      a = w(1)
      b = w(2)
      r0 = w(4)
      rend = w(5)

      ilegc = w(13)
      nlege = w(14)
      rb = w(15)

      if (ifbell .eq. 1 .and. rsupp .lt. a) then
         write(*,*) 'bad value of rsupp! turning off bell'
         ifbell = 1
      endif
      
      call eatonprol_eval0(ifbell,rsupp,rb,r0,rend,a,b,rs,nr,
     1     w(ilegc),nlege,vals)
         
      return
      end

      subroutine eatonprol_eval0(ifbell,rsupp,rb,r0,rend,a,b,rs,nr,
     1     clege,nlege,vals)
      implicit real *8 (a-h,o-z)
      real *8 rs(*), vals(*), clege(*)

      if (ifbell .eq. 1) then
 
         do i = 1,nr
            r = rs(i)
            t = (r-r0)*2.0d0/(rend-r0)-1.0d0
            call legeexev(t,val0,clege,nlege)

            vals(i) = val0
            
            if (r .lt. rb) then
               x = 12.0d0*(r-rb/2.0d0)/rb
               call qerrfun(x,val1)
               eta = (1.0d0 + val1)/2.0d0
               vals(i) = val0*eta + (1.0d0-eta)*b
            endif

            if (r .gt. a) then
               x = 12.0d0*(r-(a+rsupp)/2.0d0)/(rsupp-a)
               call qerrfun(x,val1)
               eta = (1.0d0 + val1)/2.0d0
               vals(i) = (1.0d0-eta)*val0
            endif
         enddo
      else

         do i = 1,nr
            r = rs(i)
            t = (r-r0)*2.0d0/(rend-r0)-1.0d0
            call legeexev(t,val0,clege,nlege)

            vals(i) = val0
         enddo
      endif

      return
      end
