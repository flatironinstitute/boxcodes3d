
      subroutine zerrf(z,derfz,n)
c                
C     Double precision cmplx error function ...
c     z = arg
c     n = number of quad nodes (approx 4*abs(x*y))
C
      implicit real *8 (a-h,o-z)
      integer n
      complex *16 z,derfz,delta,a,b,eye
      real *8, allocatable ::  xs(:),ws(:)
C
c     erf(x+iy) = erf(x) + 
c      i \frac{2}{sqrt{pi}} e^{-x^2} int_0^y e^{u^2} cos 2xu du +
c      \frac{2}{sqrt{pi}} e^{-x^2} int_0^y e^{u^2} sin 2xu du.
      
c
c---- set order of quadrature
c
      eye = dcmplx(0,1)
      pi = 4.0d0*datan(1.0d0)
      allocate(xs(n))
      allocate(ws(n))
      xx = dreal(z)
      yy = dimag(z)
      er = erf(xx)
      fac = 2.0d0/dsqrt(pi)*dexp(-xx*xx)
      rint1 = 0.0d0
      rint2 = 0.0d0
c
c---- compute N Chebyshev nodes :
c
      ifwhts = 1
      call legewhts(n,xs,ws,ifwhts)
      a = 0.0d0
      b = yy
      do i=1,n
        xs(i) = (b-a)*xs(i)/2.0d0 + (b+a)/2.0d0
        ws(i) = ws(i)*(b-a)/2.0d0
      enddo
C
      do 200 i=1,n
         u = xs(i)
         rint1 = rint1 + dexp(u*u)*dcos(2*xx*u)*ws(i)
         rint2 = rint2 + dexp(u*u)*dsin(2*xx*u)*ws(i)
200   continue
ccc      write(6,*) 'rint1 = ',rint1
ccc      write(6,*) 'rint2 = ',rint2
      derfz = er + rint2 + eye*rint1
c      i \frac{2}{sqrt{pi}} e^{-x^2} int_0^y e^{u^2} cos 2xu du +
c      \frac{2}{sqrt{pi}} e^{-x^2} int_0^y e^{u^2} sin 2xu du.
      return
      end
