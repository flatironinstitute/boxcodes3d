c
c
c    This file contains the routines for reprensting polynomials
c    in spherical coordinates
c
c     x^{n}y^{m}z^{\ell} can be expressed in terms of functions
c      of the form Y_{n',m'}(\theta,\phi) r^{n'+2l'}
c
c    This file contains the following user callable routines
c       cart2polar - convert from cartesian to spherical coordinates
c       eval_spherepol_basis - evaluate the basis 
c          r^{n'+2l'} Y_{n',m')(\theta,\phi) at a given point
c       
c          The basis functions are ordered as:
c             Y_{0,0} r^{0}, Y_{0,0} r^{2} ... Y_{0,0} r^{p}
c             Y_{1,-1} r^{1}, Y_{1,-1} r^{3} ... Y_{1,-1} r^{p-1}
c             Y_{1,0} r^{1}, Y_{1,0} r^{3} ... Y_{1,0} r^{p-1}
c             Y_{1,1} r^{1}, Y_{1,1} r^{3} ... Y_{1,1} r^{p-1}
c                         .
c                         .
c                         .
c                         .
c             Y_{p,-p} r^{0} 
c             Y_{p,-p+1} r^{0} 
c                         .
c                         .
c                         .
c             Y_{p,p-1} r^{0} 
c             Y_{p,p} r^{0} 
c
c      legetens_spherepol - stores the map from tensor product legendre
c        polynomials to coeffs in the spherical polynomial basis
c   


c**********************************************************************
      subroutine cart2polar(zat,r,theta,phi)
c**********************************************************************
c
c     Convert from Cartesian to spherical coordinates.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c	zat   :  Cartesian vector
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c	r     :  |zat|
c	theta : angle subtended with respect to z-axis
c	phi   : angle of (zat(1),zat(2)) subtended with 
c               respect to x-axis
c
c-----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 zat(3)
      complex *16 ephi,eye
      data eye/(0.0d0,1.0d0)/
c
c 
      r= sqrt(zat(1)**2+zat(2)**2+zat(3)**2)
      proj = sqrt(zat(1)**2+zat(2)**2)
c
      theta = datan2(proj,zat(3))
      if( abs(zat(1)) .eq. 0 .and. abs(zat(2)) .eq. 0 ) then
      phi = 0
      else
      phi = datan2(zat(2),zat(1))
      endif
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

      subroutine eval_spherepol_basis(n,npols,r,theta,phi,pols)
c
c        computes the spherical polynomials at one point
c        see top of document for polynomial order
c
c
      implicit real *8 (a-h,o-z)
      real *8 ynm(0:n,0:n)
      complex *16 pols(*),ima

      ima = dcmplx(0.0d0,1.0d0)
      ct = cos(theta)

      call ylgndr(n,ct,ynm)

      ipol = 1

      do j=0,n-1
        do m=-j,j
          do l=0,n/2
            if(j+2*l.ge.n) goto 1000
            pols(ipol)=r**(j+2*l)*ynm(j,abs(m))*exp(ima*m*phi)
            ipol = ipol+1
 1000       continue
          enddo
        enddo
      enddo


      return
      end
c
c
c
c
c
c

      subroutine eval_spherepol_basis_wder(n,npols,r,theta,phi,pols,der)
c
c        computes the spherical polynomials at one point
c        see top of document for polynomial order
c
c
      implicit real *8 (a-h,o-z)
      real *8, allocatable :: ynm(:,:),ynmd(:,:)
      complex *16 pols(*),ima,der(3,*),dpdr,dpdt,dpdphi

      allocate(ynm(0:n,0:n),ynmd(0:n,0:n))

      ima = dcmplx(0.0d0,1.0d0)

      ctheta = cos(theta)
      cphi = cos(phi)
      stheta = sin(theta)
      sphi = sin(phi)


      call ylgndr2s(n,ctheta,ynm,ynmd)



cc      call prin2('ynmd=*',ynmd,(n+1)*(n+1))
c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
      rx = stheta*cphi
      thetax = ctheta*cphi
      phix = -sphi
      ry = stheta*sphi
      thetay = ctheta*sphi
      phiy = cphi
      rz = ctheta
      thetaz = -stheta
      phiz = 0.0d0


      ipol = 1

      der(1,1) = 0
      der(2,1) = 0
      der(3,1) = 0

      do j=0,n-1
        do m=-j,j
          do l=0,n/2
            if(j+2*l.ge.n) goto 1000
            if(m.eq.0) then
              pols(ipol)=r**(j+2*l)*ynm(j,0)

              if(j+2*l.gt.0) then
                dpdr = (j+2*l)*r**(j+2*l-1)*ynm(j,0)
                dpdt = -r**(j+2*l-1)*ynmd(j,0)*stheta
                der(1,ipol) = dpdr*rx + dpdt*thetax
                der(2,ipol) = dpdr*ry + dpdt*thetay
                der(3,ipol) = dpdr*rz + dpdt*thetaz
              endif
            else
              pols(ipol)=r**(j+2*l)*ynm(j,abs(m))*exp(ima*m*phi)*stheta

              if(j+2*l.gt.0) then
                dpdr = (j+2*l)*r**(j+2*l-1)*ynm(j,abs(m))*
     1              exp(ima*m*phi)*stheta
                dpdt = -r**(j+2*l-1)*ynmd(j,abs(m))*exp(ima*m*phi)
                dpdphi = r**(j+2*l-1)*ynm(j,abs(m))*exp(ima*m*phi)*
     1             ima*m 
                der(1,ipol) = dpdr*rx + dpdt*thetax + dpdphi*phix
                der(2,ipol) = dpdr*ry + dpdt*thetay + dpdphi*phiy
                der(3,ipol) = dpdr*rz + dpdt*thetaz
              endif
            endif
            
            ipol = ipol+1
 1000       continue
          enddo
        enddo
      enddo

cc      call prin2('der=*',der,6*npols)


      return
      end
c
c
c
c
c
c
c
      subroutine legetens_spherepol(n,npols,xmat)
c
c      stores the map from tensor product legendre
c      polynomials to coeffs in the spherical polynomial basis
c
c     both bases are assumed to be 'total degree' (see below)
c      
c      input parameters
c        n - max total order of tensor product legendre polynomials
c            all polynomials of total degree < n are considered
c        npols - total number of polynomials of order less than
c           total order n. npols = n(n+1)(n+2)/6
c
c      output parameters
c         xmat - complex *16 (npols,npols)
c           tensor product legendre polynomials are ordered as in
c           legetens.f
c
      implicit real *8 (a-h,o-z)
      integer n,npols
      complex *16 xmat(npols,npols)
      complex *16, allocatable :: xlsq(:,:),frhs(:),pols(:)
      complex *16, allocatable :: w(:)
      real *8, allocatable :: xx(:),xr(:),xt(:),xphi(:),xpolg(:,:)
      real *8, allocatable :: rnorms(:)
      complex *16 ima
      integer ix,iy,iz,i,j,k,l,m,itype,nn,nw,ipt,ncols,ipolno
      real *8 done,pi,u,v,w0,r,theta,phi,eps

      done = 1
      pi = atan(done)*4
      allocate(xx(n),xr(n),xt(n),xphi(2*n+1))
      
      itype = 0
      call legeexps(itype,n,xx,u,v,w0)

      do i=1,n
        xr(i) = (xx(i)+1)/2*sqrt(3.0d0)
        xt(i) = (xx(i)+1)/2*pi
      enddo

      do i=1,2*n+1
        xphi(i) = (i-1.0d0)/(2*n+1)*2*pi
      enddo
      
      nn = n*n*(2*n+1)
      allocate(xpolg(3,nn),frhs(nn))
c
c      generate the grid
c

      ipt = 1
      do i=1,n
        do j=1,n
          do l=1,2*n+1
            xpolg(1,ipt) = xr(i)*sin(xt(j))*cos(xphi(l))
            xpolg(2,ipt) = xr(i)*sin(xt(j))*sin(xphi(l))
            xpolg(3,ipt) = xr(i)*cos(xt(j))
            ipt = ipt + 1
          enddo
        enddo
      enddo


      nw = nn*npols*10 + 1000
      allocate(xlsq(nn,npols),w(nw),pols(npols))



      do ipt=1,nn
        call cart2polar(xpolg(1,ipt),r,theta,phi)
        call eval_spherepol_basis(n,npols,r,theta,phi,pols)


        do i=1,npols
          xlsq(ipt,i) = pols(i)
        enddo
      enddo

      ncols = 0
      allocate(rnorms(nn+10))
      eps = 1.0d-14
      call ncleastsq(xlsq,nn,npols,eps,ncols,rnorms,w)

      ipolno = 1
      do iz=0,n-1
        do iy=0,n-1-iz
          do ix=0,n-1-iz-iy
c
c             generate the right hand side for this
c             particular polynomial
c
            do ipt=1,nn
              call legepol(xpolg(1,ipt),ix,polx,der)
              call legepol(xpolg(2,ipt),iy,poly,der)
              call legepol(xpolg(3,ipt),iz,polz,der)
              frhs(ipt) = polx*poly*polz
            enddo

c
c         initialize the solution
c
            do i=1,npols
              xmat(i,ipolno) = 0
            enddo

            call ncleasts2(w,frhs,xmat(1,ipolno))

            ipolno = ipolno + 1
          enddo
        enddo
      enddo



      return
      end


