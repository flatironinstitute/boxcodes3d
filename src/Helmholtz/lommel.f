c
c
c     subroutines for evaluating Lommel functions
c     
c


      subroutine zloml1ps_dbg(mu,nu,z,tol,val,stab,ip)
c
c     debugging routine
c      
c     evaluate the 1st kind Lommel function s_{mu,nu}(z) where
c     mu, nu, and z are double precision complex numbers
c     using the standard power series. we must have
c     that mu \pm nu \neq  -(2m-1) for m a positive
c     integer (so that -(2m-1) is a negative odd int)
c
c     def
c
c     val = z^{mu+1} sum_{j=0}^p (-1)^j z^{2j}/a_{k+1}(mu,nu)
c
c     where p is chosen so that the last included
c     term is less than tol and
c
c     a_{k}(mu,nu) = \prod_{m=1}^k ( (mu+2*m-1)^2 - nu^2 )
c
c     a maximum number of terms is hardcoded (200).
c     This routine is only recommended for z with a
c     reasonably small modulus (say, less than 4)
c
c     input
c
c     mu - complex *16 parameter
c     nu - complex *16 parameter
c     z - complex *16 argument
c     tol - real *8 requested tolerance
c
c     output
c
c     val - complex *16 lommel function as defined above
c     stab - largest term in sum (gives a sense of stability
c     of calculation)
c     ip - integer, highest term (p) in the formula above
c     
      implicit none
      integer maxterms, ip 
      parameter (maxterms = 200)
      complex *16 mu, nu, z, val
      real *8 stab, tol
c     local
      complex *16 z2, ak, zterm, nu2, mup2m1, one
      real *8 szterm
      data one /(1.0d0, 0.0d0)/

      nu2 = nu*nu
      z2 = z*z
      mup2m1 = mu+1
      ak = mup2m1*mup2m1-nu2
      zterm = one/ak
      
      val = zterm
      stab = abs(val)

      do ip = 1,maxterms
         mup2m1 = mup2m1+2
         ak = mup2m1*mup2m1-nu2
         zterm = -zterm*z2/ak
         val = val + zterm
         szterm = abs(zterm)         
         stab = max(stab,szterm)
         if (szterm .lt. tol) exit
      enddo

      val = val*z**(mu+1)
      stab = sqrt(stab)

      return
      end
         
      
      subroutine zloml1ps_forpnzpow(mu,nu,z,ip,val)
c
c     fixed order, real parameter version. omits power of
c     z^{mu+1} (can often be computed outside to reduce cost
c     of multiple calls with same mu and z)
c      
c     evaluate a scaled version of the 1st kind Lommel
c     function s_{mu,nu}(z), where
c     z is a double precision complex number and mu and
c     nu are double precision numbers,
c     using the standard power series but omitting the leading
c     factor of z^{mu+1}. we must have
c     that mu \pm nu \neq  -(2m-1) for m a positive
c     integer (so that -(2m-1) is a negative odd int)
c
c     def
c
c     val = sum_{j=0}^p (-1)^j z^{2j}/a_{k+1}(mu,nu)
c
c     where p is specified by the user and
c
c     a_{k}(mu,nu) = \prod_{m=1}^k ( (mu+2*m-1)^2 - nu^2 )
c
c     This routine is only recommended for z with a
c     reasonably small modulus (say, less than 4)
c
c     input
c
c     mu - real *8 parameter
c     nu - real *8 parameter
c     z - complex *16 argument
c     ip - number of terms
c
c     output
c
c     val - complex *16 lommel function as defined above
c     
      implicit none
      integer ip 
      real *8 mu, nu
      complex *16 z, val
c     local
      integer i
      complex *16 z2, zterm, one
      real *8 ak, mup2m1, nu2
      data one /(1.0d0, 0.0d0)/

      nu2 = nu*nu
      z2 = z*z
      mup2m1 = mu+1
      ak = mup2m1*mup2m1-nu2
      zterm = one/ak
      
      val = zterm

      do i = 1,ip
         mup2m1 = mup2m1+2
         ak = mup2m1*mup2m1-nu2
         zterm = -zterm*z2/ak
         val = val + zterm
      enddo

      return
      end
         
      
      subroutine zloml1ps_forpnzpow_der(mu,nu,z,ip,val,der)
c
c     fixed order, real parameter version. omits power of
c     z^{mu+1} (can often be computed outside to reduce cost
c     of multiple calls with same mu and z)
c      
c     evaluate a scaled version of the 1st kind Lommel
c     function s_{mu,nu}(z), where
c     z is a double precision complex number and mu and
c     nu are double precision numbers,
c     using the standard power series but omitting the leading
c     factor of z^{mu+1}. we must have
c     that mu \pm nu \neq  -(2m-1) for m a positive
c     integer (so that -(2m-1) is a negative odd int)
c
c     def
c
c     val = sum_{j=0}^p (-1)^j z^{2j}/a_{k+1}(mu,nu)
c
c     where p is specified by the user and
c
c     a_{k}(mu,nu) = \prod_{m=1}^k ( (mu+2*m-1)^2 - nu^2 )
c
c     This routine is only recommended for z with a
c     reasonably small modulus (say, less than 4)
c
c     input
c
c     mu - real *8 parameter
c     nu - real *8 parameter
c     z - complex *16 argument
c     ip - number of terms
c
c     output
c
c     val - complex *16 lommel function as defined above
c     der - derivative with respect to z
c     
      implicit none
      integer ip 
      real *8 mu, nu
      complex *16 z, val, der
c     local
      integer i
      complex *16 z2, zterm, one, zero
      real *8 ak, mup2m1, nu2
      data one /(1.0d0, 0.0d0)/
      data zero / (0.0d0,0.0d0) /

      nu2 = nu*nu
      z2 = z*z
      mup2m1 = mu+1
      ak = mup2m1*mup2m1-nu2
      zterm = one/ak
      
      val = zterm
      der = zero

      do i = 1,ip
         mup2m1 = mup2m1+2
         ak = mup2m1*mup2m1-nu2
         zterm = -zterm*z/ak
         der = der + zterm*2*i
         zterm = zterm*z
         val = val + zterm
      enddo

      return
      end
         
      
