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
c     This file currently contains the following routines
c     legecoeff_zantihelm3d_form - Form the anti-Helmholtzian
c     
c     facelayerpot_eval - evaluate the potential
c     due to single and double layer potential
c     with arbitrary polynomial densities, to all
c     points in the near field
c     


      subroutine h3dtabp_ref(ndeg,zk,tol,tab,ldtab)
c
c     generate the Helmholtz potential table at the reference
c     points
c
c
c     input
c
c     ndeg - integer, highest degree of the basis polynomials
c                    (measured in total degree, e.g. x^0y^1z^2
c                     has total degree 3)
c     zk - complex*16, helmholtz parameter
c     tol - tolerance for error in table entries
c     ldtab - integer, leading dimension of the output table
c
c     output
c      
c     tab - complex *16 array (ldtab,*), tab(i,j) is the integral
c     of the j-th tensor polynomial (in the ordering specified
c     in legetens.f) against the scaled green's function
c     exp(ikr)/r at the i-th reference target point (see tensrefpts3d)
      
c      
      implicit none
      integer ndeg, ldtab
      complex *16 tab(ldtab,*), zk
      real *8 tol
c     local
      integer idims(6), ndeg2, ii, jj, iface, npol2, npol3, ifdiff
      integer ntarg0, ntarg, n, idim, istart, npt, itype
      real *8 slicevals(6), flipd(6), flips(6), abszk, abszktol
      real *8 rcond, val, derscale, tol2, r, theta, phi
      parameter (abszktol = 2.5d0)

      real *8, allocatable :: x(:,:), w(:), pols(:), v(:,:)
      real *8 u, pi4
      integer ldu, ldv
      
      complex *16 zero, im, one
      complex *16, allocatable :: tabtemp(:,:), ahc(:,:), zv(:,:)
      complex *16, allocatable :: ahderc(:,:), ahcleg(:,:)
      complex *16, allocatable :: ahdercleg(:,:), ahelm(:), ahelms(:,:)
      complex *16, allocatable :: ahcleg3(:,:), leg2sph(:,:)
      complex *16, allocatable :: slp_pots(:,:), dlp_pots(:,:)

      
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /      
      data im / (0.0d0,1.0d0) /
      data slicevals / -1.0d0, 1.0d0, -1.0d0, 1.0d0, -1.0d0, 1.0d0 /
      data flipd / -1.0d0, 1.0d0, -1.0d0, 1.0d0, -1.0d0, 1.0d0 /
      data flips / 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0 /
      data idims / 1, 1, 2, 2, 3, 3 /
      
      logical ifsphere
      character type

      pi4 = 16*atan(1.0d0)
      
      n = ndeg+1
      
      type = 'T'
      
      abszk = abs(zk)

c     determine anti-helmholtzian technique
      
      if (abszk .gt. abszktol) then
c     compute on coefficients
         ifsphere = .false.
         ndeg2 = ndeg
      else
c     compute using analytic solution on spherical polys
         ifsphere = .true.
         ndeg2 = 20
      endif

c     memory for coeffs, etc.
      
      call legetens_npol_3d(ndeg,type,npol3)
      call legetens_npol_2d(ndeg2,type,npol2)

      allocate(ahc(npol2,npol3),ahderc(npol2,npol3))
      allocate(ahcleg3(npol3,npol3),leg2sph(npol3,npol3))
      allocate(ahcleg(npol2,npol3),ahdercleg(npol2,npol3))

      npt = n**3
      ldu = 1
      ldv = npt
      itype = 4
      allocate(x(3,npt),w(npt),v(npt,npol3))
      allocate(zv(npt,npol3))
      allocate(ahelm(npol3),ahelms(npol3,npt))
      call legetens_exps_3d(itype,n,type,x,u,ldu,v,ldv,w)

      
      if (.not. ifsphere) then
         tol2 = 1.0d-12
         call h3danti_legetens_form(ndeg,type,zk,tol2,ahcleg3,
     1        npol3,rcond)
         do ii = 1,npol3
            do jj = 1,npt
               zv(jj,ii) = -pi4*v(jj,ii)
            enddo
         enddo
         
         call zgemm('N','N',npt,npol3,npol3,one,zv,npt,
     1        ahcleg3,npol3,zero,tab,ldtab)
         
      else
         call legetens_spherepol(n,npol3,leg2sph)
         
         do jj = 1,npt
            call sphcart2polar(x(1,jj),r,theta,phi)
            call h3danti_sphere(ndeg,zk,r,theta,phi,ahelm)
            do ii = 1,npol3
               ahelms(ii,jj) = -pi4*ahelm(ii)
            enddo
         enddo

         call zgemm('T','N',npt,npol3,npol3,one,ahelms,npol3,
     1        leg2sph,npol3,zero,tab,ldtab)
         
      endif

c     memory depending on ntarg

      ntarg0 = 10*npt
      ntarg = 6*ntarg0
      allocate(slp_pots(npol2,ntarg),dlp_pots(npol2,ntarg))
      allocate(tabtemp(ntarg0,npol3))

      call h3d_facelayerpot_eval(tol,zk,ndeg,ndeg2,type,slp_pots,
     1  dlp_pots,npol2)

      do ii = 1,npol3
         do jj = npt+1,ntarg0
            tab(jj,ii) = zero
         enddo
      enddo

      do iface = 1,6

         idim = idims(iface)
         val = slicevals(iface)
         derscale = val

c     get poly coeffs of anti-helmholtzian of each polynomial
c     on face
         if (ifsphere) then
            call h3danti_sphere_slicecoeffs(ndeg2,type,ndeg,zk,
     1           idim,val,derscale,ahc,ahderc,npol2)
            call zgemm('N','N',npol2,npol3,npol3,one,ahc,npol2,
     1           leg2sph,npol3,zero,ahcleg,npol2)
            call zgemm('N','N',npol2,npol3,npol3,one,ahderc,npol2,
     1           leg2sph,npol3,zero,ahdercleg,npol2)
            
         else
            ifdiff = 1
            call legetens_slicezcoeffs_3d(ndeg,type,idim,val,
     1           ahcleg3,npol3,npol3,ahcleg,npol2,ifdiff,ahdercleg,
     2           npol2)
            do ii = 1,npol3
               do jj = 1,npol2
                  ahdercleg(jj,ii) = derscale*ahdercleg(jj,ii)
               enddo
            enddo
         endif
         

         istart = (iface-1)*ntarg0 + 1
         call zgemm('T','N',ntarg0,npol3,npol2,one,
     1        slp_pots(1,istart),
     1        npol2,ahdercleg,npol2,zero,tabtemp,ntarg0)

         do ii = 1,npol3
            do jj = 1,ntarg0
               tab(jj,ii) = tab(jj,ii) +
     1              flips(iface)*tabtemp(jj,ii)
            enddo
         enddo

         istart = (iface-1)*ntarg0 + 1 
         call zgemm('T','N',ntarg0,npol3,npol2,one,
     1        dlp_pots(1,istart),
     1        npol2,ahcleg,npol2,zero,tabtemp,ntarg0)

         do ii = 1,npol3
            do jj = 1,ntarg0
               tab(jj,ii) = tab(jj,ii) +
     1              flipd(iface)*tabtemp(jj,ii)
            enddo
         enddo
         
      enddo         
      
      return
      end

      subroutine h3danti_sphere(ndeg,zk,r,theta,phi,ahelm)
c
c     this subroutine returns an anti-Helmholtzian of
c     each of the spherical polynomials (n = ndeg+1)
c
c     p_jml = Y_jm(cos(theta)) exp(i m phi) r^(j+2l), where
c
c     j = 0, ..., n-1
c     m = -j,...,j
c     l = 0,...,ceil((n-j)/2)-1
c
c     these can be thought of as the total degree polynomials
c     on the sphere, up to degree n-1
c
c     The order traversed has l on the innermost loop, m
c     on the loop outside that, and j on the outermost loop
c      
      implicit real *8 (a-h,o-z)
      real *8 ynm(0:ndeg+1,0:ndeg+1)
      complex *16 ahelm(*),ima,vals(0:ndeg+1),zk,z,zph1,zphm,zphj
      real *8 mu,nu,rj,r2,rl
      integer ip, k, n
      parameter (ip = 14)

      n = ndeg+1
      
      ima = dcmplx(0.0d0,1.0d0)
      ct = cos(theta)

      z=zk*r

      call ylgndr(n,ct,ynm)

c      call prin2('r=*',r,1)
c      call prin2('theta=*',theta,1)
c      call prin2('phi=*',phi,1)
      ipol = 1

      zph1 = exp(ima*phi)
      zphj = 1
      zphm = 1

      rj = 1
      r2 = r*r
      do j=0,n-1
         nu = j+0.5d0
         do k = 0,(n-j)/2
            mu = j+2*k+1.5d0
            call zloml1ps_forpnzpow(mu,nu,z,ip,vals(k))
         enddo
         zphm = zphj
         do m=-j,j
            rl = rj
            do l=0,n/2
               rl = rl*r2
               if(j+2*l.ge.n) goto 2000
c     pols(ipol)=r**(j+2*l)*ynm(j,abs(m))*exp(ima*m*phi)
c     mu = j+2*l+1.5d0
c     call zloml1ps_forpnzpow(mu,nu,z,ip,val)
               ahelm(ipol)=rl*vals(l)*ynm(j,abs(m))*
     1              zphm
               ipol = ipol+1
 2000          continue
            enddo
            zphm = zphm*zph1
         enddo
         zphj = zphj/zph1
         rj = rj*r
      enddo
c      call prinf('ipol=*',ipol,1)

      return
      end
      

      subroutine h3danti_sphere_wder(ndeg,zk,r,theta,phi,ahelm,der)
c
c     this subroutine returns an anti-Helmholtzian of
c     each of the spherical polynomials (n = ndeg+1)
c
c     p_jml = Y_jm(cos(theta)) exp(i m phi) r^(j+2l), where
c
c     j = 0, ..., n-1
c     m = -j,...,j
c     l = 0,...,ceil((n-j)/2)-1
c
c     these can be thought of as the total degree polynomials
c     on the sphere, up to degree n-1
c
c     The order traversed has l on the innermost loop, m
c     on the loop outside that, and j on the outermost loop
c      
      implicit real *8 (a-h,o-z)
      real *8 ynm(0:(ndeg+1),0:(ndeg+1)), ynmd(0:(ndeg+1),0:(ndeg+1))
      complex *16 ahelm(*),ima,vals(0:ndeg+1),zk,z,zph1,zphm,zphj
      complex *16 ders(0:ndeg+1)
      complex *16 der(3,*),dpdr,dpdt,dpdphi
      real *8 mu,nu,rj,r2,rl,rlm1,rjm1
      integer ip, k, n
      parameter (ip = 14)

      n = ndeg+1
      
      ima = dcmplx(0.0d0,1.0d0)
      ctheta = cos(theta)
      cphi = cos(phi)
      stheta = sin(theta)
      sphi = sin(phi)

      z=zk*r

      call ylgndr2s(n,ctheta,ynm,ynmd)

      ipol = 1

      zph1 = exp(ima*phi)
      zphj = 1
      zphm = 1

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

      der(1,1) = 0
      der(2,1) = 0
      der(3,1) = 0

      rj = 1
      rjm1 = r
      r2 = r*r
      do j=0,n-1
         nu = j+0.5d0
         do k = 0,(n-j)/2
            mu = j+2*k+1.5d0
            call zloml1ps_forpnzpow_der(mu,nu,z,ip,vals(k),ders(k))
            ders(k) = ders(k)*zk
         enddo
         zphm = zphj
         do m=-j,j
            rl = rj
            rlm1 = rjm1
            do l=0,n/2
               rl = rl*r2
               if(j+2*l.ge.n) goto 2000
               if(m.eq.0) then
                  ahelm(ipol)=rl*vals(l)*ynm(j,0)

                  dpdr = ((j+2*l+2)*rlm1*vals(l)+rl*ders(l))
     1                 *ynm(j,0)
                  dpdt = -rlm1*vals(l)*ynmd(j,0)*stheta
                  der(1,ipol) = dpdr*rx + dpdt*thetax
                  der(2,ipol) = dpdr*ry + dpdt*thetay
                  der(3,ipol) = dpdr*rz + dpdt*thetaz

               else
                  ahelm(ipol)=rl*vals(l)*ynm(j,abs(m))*
     1                 zphm*stheta

                  dpdr = ((j+2*l+2)*rlm1*vals(l)+rl*ders(l))
     1                 *ynm(j,abs(m))*zphm*stheta
                  dpdt = -rlm1*vals(l)*ynmd(j,abs(m))*zphm
                  dpdphi = rlm1*vals(l)*ynm(j,abs(m))*zphm*ima*m 
                  der(1,ipol) = dpdr*rx + dpdt*thetax + dpdphi*phix
                  der(2,ipol) = dpdr*ry + dpdt*thetay + dpdphi*phiy
                  der(3,ipol) = dpdr*rz + dpdt*thetaz

                  
               endif
               ipol = ipol+1
               rlm1 = rlm1*r2
               
 2000          continue
            enddo
            zphm = zphm*zph1
         enddo
         zphj = zphj/zph1
         rj = rj*r
         rjm1 = rjm1*r
      enddo
c      call prinf('ipol=*',ipol,1)

      return
      end



      subroutine h3danti_sphere_slicecoeffs(ndegslice,type,ndegahelm,
     1     zk,idim,val,derscale,ahc,ahderc,ldahc)
c
c     get the coefficients of the antihelmholtzian as a function
c     on the face defined by idim and val
c
c     idim - denotes which coordinate is held fixed (e.g. idim = 1
c            gives x = val, idim = 2 gives y = val)
c     val - the fixed value of the coordinate specified by idim
c           
c
c      
      implicit none
      integer ndegahelm, idim, ndegslice, ldahc
      complex *16 zk, ahc(ldahc,*), ahderc(ldahc,*)
      real *8 val, derscale
      character type
c     local
      complex *16, allocatable :: der3(:,:)
      complex *16, allocatable :: ahpolsall(:,:), ahderall(:,:)
      real *8 r, theta, phi, x(3,(ndegslice+1)*(ndegslice+1))
      real *8 u((ndegslice+1)**4), v((ndegslice+1)**4)
      real *8 w((ndegslice+1)**2)
      complex *16 zu( (ndegslice+1)**4 ), one, zero
      real *8 x2(2,(ndegslice+1)*(ndegslice+1))
      integer itype, n, ldu, ldv, npol2, npt2, npol
      integer isub(2,3), isub1, isub2, j, i
      data isub / 2, 3, 1, 3, 1, 2 /
      data one / (1.0d0,0.0d0) /
      data zero / (0.0d0,0.0d0) /      

c     figure out where to load 2d data in 3d
      
      isub1 = isub(1,idim)
      isub2 = isub(2,idim)

c     get points, transforms for a face
      
      call legetens_npol_2d(ndegslice,type,npol2)

      n = ndegslice+1
      npt2 = n**2
      itype = 2
      ldu = npol2
      ldv = npt2
      call legetens_exps_2d(itype,n,type,x2,u,ldu,v,ldv,w)

      do i = 1,npt2
         x(idim,i) = val
         x(isub1,i) = x2(1,i)
         x(isub2,i) = x2(2,i)
      enddo

c     grab antihelmholtzian values

      call legetens_npol_3d(ndegahelm,'t',npol)

      allocate(ahpolsall(npol,npt2),der3(3,npol),ahderall(npol,npt2))
      
      do i = 1,npt2
         call sphcart2polar(x(1,i),r,theta,phi)
         call h3danti_sphere_wder(ndegahelm,zk,r,theta,phi,
     1        ahpolsall(1,i),der3)
         do j = 1,npol
            ahderall(j,i) = derscale*der3(idim,j)
         enddo
      enddo
         
c     tranform these values

      do i = 1,npt2*npol2
         zu(i) = u(i)
      enddo
      
      call zgemm('N','T',npol2,npol,npt2,one,zu,ldu,ahpolsall,npol,
     1     zero,ahc,ldahc)
      call zgemm('N','T',npol2,npol,npt2,one,zu,ldu,ahderall,npol,
     1     zero,ahderc,ldahc)
      

      return
      end
      
      
      
      
      subroutine h3danti_legetens_form(ndeg,type,zk,tol,ahmat,
     1     ldahmat,rcond)
      implicit none
      integer ndeg, ldahmat
      complex *16 ahmat(ldahmat,*), zk
      real *8 rcond, tol
      character type
c     local
      complex *16, allocatable :: work(:), helm(:,:), u(:,:), vt(:,:)
      complex *16 one, zero
      real *8, allocatable :: rwork(:), dlap(:,:), eye(:,:)
      real *8, allocatable :: sval(:)
      real *8 s1, sc, tols1
      integer lwork, j, i, npol, info, ndub, incx
      data one / (1.0d0,0.0d0) /
      data zero / (0.0d0,0.0d0) /

      call legetens_npol_3d(ndeg,type,npol)
      npol = npol
      
      lwork = 10*npol
      allocate(work(lwork),rwork(lwork))
      allocate(u(npol,npol),vt(npol,npol),sval(npol))
      allocate(dlap(npol,npol),helm(npol,npol),eye(npol,npol))

      call legetens_lapmat_3d(ndeg,ndeg,type,dlap,npol)
      call legetens_eyemat_3d(ndeg,ndeg,type,eye,npol)

      do i = 1,npol
         do j = 1,npol
            helm(j,i) = dlap(j,i) + zk*zk*eye(j,i)
         enddo
      enddo

      call zgesvd('S','S',npol,npol,helm,npol,sval,u,npol,
     1     vt,npol,work,lwork,rwork,info)
      rcond = sval(1)/sval(npol)

cc      call prin2('zk=*',zk,2)
cc      call prin2('condition number is*',rcond,1)
cc      call prinf('npol=*',npol,1)
cc      call prinf('info=*',info,1)
cc      call prin2('sval=*',sval(1),1)
cc      call prin2('sval=*',sval(npol),1)
cc      call prin2('svals all=*',sval,npol)
      

c      call prinf('after svd, info = *',info,1)

c     form u*diag(s^{-1}) (reduced)

      s1 = sval(1)
      ndub = 2*npol
      incx = 1
      tols1 = tol*s1
      do i = 1,npol
         sc = 0.0d0
         if (sval(i) .ge. tols1) then
            sc = 1.0d0/sval(i)
            rcond = sc*s1
         endif
c     kinda hacky, use real valued routine
         call dscal(ndub,sc,u(1,i),incx)
      enddo

c      call prin2('rcond *',rcond,1)

c     form (vt)^T*(u * diag(s^{-1}))^T
      call zgemm('C','C',npol,npol,npol,one,vt,npol,
     1     u,npol,zero,ahmat,ldahmat)
      
      return
      end
c
c
c
c     
c
c
      subroutine h3d_facelayerpot_eval_new(tol,zk,ndeg,ndegp,type,
     1  slp_pots,dlp_pots,lda)
      implicit real *8 (a-h,o-z)
      real *8 tol
      complex *16 zk
      integer ndeg,ndegp,lda
      complex *16 slp_pots(lda,*),dlp_pots(lda,*)
      complex *16, allocatable :: slp_near(:,:),dlp_near(:,:)
      complex *16, allocatable :: slp_far(:,:),dlp_far(:,:)
      character type

      real *8, allocatable :: xyztarg(:,:),xyztarg_near(:,:)
      real *8, allocatable :: xyztarg_far(:,:)
      integer, allocatable :: ifar_ind(:),inear_ind(:)
      real *8, allocatable :: xyztmp(:,:),qnodes(:,:),qwts(:)
      real *8 xq(ndeg+1),xyzc(3),u,v,w
      integer ipars(10)
      real *8 dpars(5)
      complex *16 zpars(3)

      external h3d_slp,h3d_dlp

      integer norder, norder_p, ncores

      done = 1
      pi = atan(done)*4

      norder = ndeg+1

      norder_p = ndegp+1

      call legetens_npol_2d(ndegp,type,npols)
cc      call prinf('npols=*',npols,1)
c     npols = norder_p*norder_p

      bs = 2.0d0
      xyzc(1) = -1
      xyzc(2) = -1
      xyzc(3) = -1 


      itype = 0
      call legeexps(itype,norder,xq,u,v,w)
      do i=1,norder
        xq(i) = xq(i) + 1
      enddo
c
      nppbox = norder*norder*norder
      ntarg0 = 10*nppbox
      allocate(xyztmp(3,ntarg0))

      istart1 = 4*nppbox+1
      istart2 = 7*nppbox+1
      
      call tensrefpts3d(xq,norder,bs,xyzc,xyztmp,xyztmp(1,istart1),
     1        xyztmp(1,istart2))

      ntarg = 6*ntarg0
      allocate(xyztarg(3,ntarg))

      do i=1,ntarg0
        xyztarg(1,i+0*ntarg0) = xyztmp(2,i)
        xyztarg(2,i+0*ntarg0) = xyztmp(3,i)
        xyztarg(3,i+0*ntarg0) = xyztmp(1,i)+1

        xyztarg(1,i+1*ntarg0) = xyztmp(2,i)
        xyztarg(2,i+1*ntarg0) = xyztmp(3,i)
        xyztarg(3,i+1*ntarg0) = xyztmp(1,i)-1

        xyztarg(1,i+2*ntarg0) = xyztmp(1,i)
        xyztarg(2,i+2*ntarg0) = xyztmp(3,i)
        xyztarg(3,i+2*ntarg0) = xyztmp(2,i)+1
        
        xyztarg(1,i+3*ntarg0) = xyztmp(1,i)
        xyztarg(2,i+3*ntarg0) = xyztmp(3,i)
        xyztarg(3,i+3*ntarg0) = xyztmp(2,i)-1

        xyztarg(1,i+4*ntarg0) = xyztmp(1,i)
        xyztarg(2,i+4*ntarg0) = xyztmp(2,i)
        xyztarg(3,i+4*ntarg0) = xyztmp(3,i)+1

        xyztarg(1,i+5*ntarg0) = xyztmp(1,i)
        xyztarg(2,i+5*ntarg0) = xyztmp(2,i)
        xyztarg(3,i+5*ntarg0) = xyztmp(3,i)-1
      enddo
      

      nquadmax = 8000
      nqorder = 20
      eps = tol
      call h3d_get_eps_nqorder_nqmax(tol,norder,eps,nqorder,nquadmax,
     1        nqorderf)
      
      intype = 2
cc      call prinf("Starting adap quad for near*",i,0)


      zpars(1) = zk
      ntt = 1


      t1 = second()
C$       t1 = omp_get_wtime()  
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
c$OMP& SCHEDULE(DYNAMIC)      
      do i=1,ntarg

        call cquadints_adap(eps,intype,norder_p,type,npols,ntt,
     1     xyztarg(1,i),nquadmax,h3d_slp,dpars,zpars,ipars,
     2     nqorder,slp_pots(1,i))

        call cquadints_adap(eps,intype,norder_p,type,npols,ntt,
     1     xyztarg(1,i),nquadmax,h3d_dlp,dpars,zpars,ipars,
     2     nqorder,dlp_pots(1,i))
      enddo
C$OMP END PARALLEL DO      
      t2 = second()
C$       t2 = omp_get_wtime()      

      return
      end
c
c
c
c
c
      subroutine h3d_get_eps_nqorder_nqmax(tol,norder,eps,nqorder,
     1   nqmax,nqorderf)
      implicit none
      real *8 tol,eps
      integer norder,nqorder,nqmax,iprec,nqorderf
c
c
c        fix this routine to optimize performance
c

      iprec = 0
      eps = tol
      nqorder = 10

      if(tol.lt.0.49d-2) iprec = 1
      if(tol.lt.0.49d-3) iprec = 2
      if(tol.lt.0.49d-6) iprec = 3
      if(tol.lt.0.49d-9) iprec = 4

cc      print *, "iprec=",iprec

      if(iprec.eq.0) then
c
c   norder, max(err/max(true,1)), max(err/true),max(err)/max(true)
c       4, 0.1e-2, 0.7e-2, 0.1e-2 
c       6, 0.1e-2, 0.8, 0.1e-2
c       8, 0.5e-3, 0.8e+2, 0.3e-3
c      12, 0.1e-2, 0.1e+5, 0.6e-3 
c
         nqmax = 1500
         nqorderf = 10
         eps = 0.5d-2
         if(norder.le.4) then
           nqorder = 7
         else if(norder.gt.4.and.norder.le.6) then
           nqmax = 2000
           eps = 0.5d-2
           nqorder = 12
           nqorderf = 12
         else if(norder.gt.6.and.norder.le.8) then
           nqorder = 14
           eps = 0.1d-2
           nqmax = 5000
           nqorderf = 12
         else if(norder.gt.8) then
           nqorder = 20
           eps = 0.5d-3
           nqmax = 12000
           nqorderf = 20
         endif

      endif

      if(iprec.eq.1) then
c
c   norder, max(err/max(true,1)), max(err/true),max(err)/max(true)
c       4, 0.2e-4, 0.7e-2, 0.1e-4 
c       6, 0.8e-4, 0.2, 0.6e-4
c       8, 0.1e-4, 0.8e+2, 0.8e-4
c      12, 0.2e-3, 0.1e+5, 0.1e-3 
c
         nqmax = 1500
         nqorderf = 10
         eps = 0.3d-2
         if(norder.le.4) then
           nqorder = 7
         else if(norder.gt.4.and.norder.le.6) then
           nqmax = 2000
           eps = 0.3d-2
           nqorder = 12
           nqorderf = 12
         else if(norder.gt.6.and.norder.le.8) then
           nqorder = 14
           eps = 0.6d-3
           nqmax = 5000
           nqorderf = 12
         else if(norder.gt.8) then
           nqorder = 20
           eps = 0.1d-3
           nqmax = 12000
           nqorderf = 20
         endif

      endif

      if(iprec.eq.2) then
c
c   norder, max(err/max(true,1)), max(err/true),max(err)/max(true)
c       4, 0.3e-7, 0.2e-5, 0.3e-7 
c       6, 0.3e-6, 0.3e-3, 0.2e-6
c       8, 0.5e-6, 0.8e-1, 0.3e-6
c      12, 0.5e-6, 6.6e3, 0.3e-6 
c
         nqmax = 5000
         nqorderf = 16
         eps = 0.3d-4
         if(norder.le.4) then
           nqorder = 12
         else if(norder.gt.4.and.norder.le.6) then
           nqorder = 14
         else if(norder.gt.6.and.norder.le.8) then
           nqorder = 14
           eps = 0.5d-5
           nqmax = 7000
         else if(norder.gt.8) then
           nqorder = 20
           eps = 0.1d-5
           nqmax = 20000
           nqorderf = 20
         endif

      endif

      if(iprec.eq.3) then
c
c   norder, max(err/max(true,1)), max(err/true),max(err)/max(true)
c       4, 0.4e-10, 0.3e-7, 0.3e-10 
c       6, 0.2e-9, 0.1e-5, 0.1e-9
c       8, 0.6e-10, 0.8e-3, 0.4e-10
c      12, 0.4e-10, 0.1e5, 0.3e-10 
c
         nqorder = 20
         eps = 3.0d-7
         nqmax = 4000
         nqorderf = 20
         if(norder.le.4) then
         else if(norder.gt.4.and.norder.le.6) then
           nqorderf = 22
         else if(norder.gt.6.and.norder.le.8) then
           nqmax = 7000
           eps = 3.0d-8
           nqorderf = 24
         else if(norder.gt.8) then
           nqorder = 24
           eps = 3.0d-9
           nqmax = 25000
           nqorderf = 30
         endif
      endif

      if(iprec.eq.4) then
c
c   NOT TESTED for accuracy
c
c
c   norder, max(err/max(true,1)), max(err/true),max(err)/max(true)
c
         nqorder = 24
         eps = 3.0d-10
         nqmax = 6000
         nqorderf = 24
         if(norder.le.4) then
         else if(norder.gt.4.and.norder.le.6) then
           nqorderf = 26
         else if(norder.gt.6.and.norder.le.8) then
           nqmax = 11000
           eps = 3.0d-11
           nqorderf = 28
         else if(norder.gt.8) then
           nqorder = 28
           eps = 3.0d-12
           nqmax = 30000
           nqorderf = 30
         endif
        
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
      subroutine h3d_get_nqorder_far(tol,nqorder)
      implicit none
      real *8 tol,eps
      integer norder,nqorder
c
c
c        fix this routine to optimize performance
c
      nqorder = 12
      if(tol.le.0.5d-3) nqorder = 16
      if(tol.le.0.5d-6) nqorder = 20
      if(tol.le.0.5d-9) nqorder = 25
      if(tol.le.0.5d-12) nqorder = 30
      return
      end

c
c
c
c
c

      subroutine h3d_slp(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)

      f = exp(ima*zpars(1)*rr)/rr

      return
      end

c      
c      
c
c
c
c
c

      subroutine h3d_dlp(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(2),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,z

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + y(3)**2)
      z = ima*zpars(1)*rr

      f = exp(z)*y(3)*(z-1.0d0)/rr**3

      return
      end

c
c
c
c     
c
c
      subroutine h3d_facelayerpot_eval(tol,zk,ndeg,ndegp,type,
     1  slp_pots,dlp_pots,lda)
      implicit real *8 (a-h,o-z)
      real *8 tol
      complex *16 zk
      integer ndeg,ndegp,lda
      complex *16 slp_pots(lda,*),dlp_pots(lda,*)
      complex *16, allocatable :: slp_near(:,:),dlp_near(:,:)
      complex *16, allocatable :: slp_far(:,:),dlp_far(:,:)
      character type

      real *8, allocatable :: xyztarg(:,:),xyztarg_near(:,:)
      real *8, allocatable :: xyztarg_far(:,:)
      integer, allocatable :: ifar_ind(:),inear_ind(:)
      real *8, allocatable :: xyztmp(:,:),qnodes(:,:),qwts(:)
      real *8 xq(ndeg+1),xyzc(3),u,v,w
      integer ipars(10)
      real *8 dpars(5)
      complex *16 zpars(3)

      external h3d_slp,h3d_dlp

      integer norder, norder_p, ncores

      done = 1
      pi = atan(done)*4

      norder = ndeg+1

      norder_p = ndegp+1

      call legetens_npol_2d(ndegp,type,npols)
cc      call prinf('npols=*',npols,1)
c     npols = norder_p*norder_p

      bs = 2.0d0
      xyzc(1) = -1
      xyzc(2) = -1
      xyzc(3) = -1 


      itype = 0
      call legeexps(itype,norder,xq,u,v,w)
      do i=1,norder
        xq(i) = xq(i) + 1
      enddo
c
      nppbox = norder*norder*norder
      ntarg0 = 10*nppbox
      allocate(xyztmp(3,ntarg0))

      istart1 = 4*nppbox+1
      istart2 = 7*nppbox+1
      
      call tensrefpts3d(xq,norder,bs,xyzc,xyztmp,xyztmp(1,istart1),
     1        xyztmp(1,istart2))

      ntarg = 6*ntarg0
      allocate(xyztarg(3,ntarg),xyztarg_near(3,ntarg),
     1   xyztarg_far(3,ntarg),inear_ind(ntarg),ifar_ind(ntarg))

      do i=1,ntarg0
        xyztarg(1,i+0*ntarg0) = xyztmp(2,i)
        xyztarg(2,i+0*ntarg0) = xyztmp(3,i)
        xyztarg(3,i+0*ntarg0) = xyztmp(1,i)+1

        xyztarg(1,i+1*ntarg0) = xyztmp(2,i)
        xyztarg(2,i+1*ntarg0) = xyztmp(3,i)
        xyztarg(3,i+1*ntarg0) = xyztmp(1,i)-1

        xyztarg(1,i+2*ntarg0) = xyztmp(1,i)
        xyztarg(2,i+2*ntarg0) = xyztmp(3,i)
        xyztarg(3,i+2*ntarg0) = xyztmp(2,i)+1
        
        xyztarg(1,i+3*ntarg0) = xyztmp(1,i)
        xyztarg(2,i+3*ntarg0) = xyztmp(3,i)
        xyztarg(3,i+3*ntarg0) = xyztmp(2,i)-1

        xyztarg(1,i+4*ntarg0) = xyztmp(1,i)
        xyztarg(2,i+4*ntarg0) = xyztmp(2,i)
        xyztarg(3,i+4*ntarg0) = xyztmp(3,i)+1

        xyztarg(1,i+5*ntarg0) = xyztmp(1,i)
        xyztarg(2,i+5*ntarg0) = xyztmp(2,i)
        xyztarg(3,i+5*ntarg0) = xyztmp(3,i)-1
      enddo
      
      ntarg_f = 0
      ntarg_n = 0
      
      znear = 0.6d0
      xynear = 1.6d0

      do i=1,ntarg
        x=xyztarg(1,i)
        y=xyztarg(2,i)
        z=xyztarg(3,i)
        if((abs(z).le.znear.and.abs(x).le.xynear.
     1      and.abs(y).le.xynear)) then
          ntarg_n = ntarg_n + 1
          xyztarg_near(1,ntarg_n) = xyztarg(1,i)
          xyztarg_near(2,ntarg_n) = xyztarg(2,i)
          xyztarg_near(3,ntarg_n) = xyztarg(3,i)
          inear_ind(ntarg_n) = i
        else
          ntarg_f = ntarg_f + 1
          xyztarg_far(1,ntarg_f) = xyztarg(1,i)
          xyztarg_far(2,ntarg_f) = xyztarg(2,i)
          xyztarg_far(3,ntarg_f) = xyztarg(3,i)
          ifar_ind(ntarg_f) = i
        endif
      enddo

      allocate(slp_near(npols,ntarg_n),dlp_near(npols,ntarg_n))
      allocate(slp_far(npols,ntarg_f),dlp_far(npols,ntarg_f))

      do i=1,ntarg_n
        do j=1,npols
          slp_near(j,i) = 0
          dlp_near(j,i) = 0
        enddo
      enddo

      do i=1,ntarg_f
        do j=1,npols
          slp_far(j,i) = 0
          dlp_far(j,i) = 0
        enddo
      enddo

      nquadmax = 8000
      nqorder = 20
      eps = tol
      call h3d_get_eps_nqorder_nqmax(tol,norder,eps,nqorder,nquadmax,
     1        nqorderf)
      
      intype = 2
cc      call prinf("Starting adap quad for near*",i,0)


      zpars(1) = zk

      nbatches = 24
      nttpcore = ceiling((ntarg_n+0.0d0)/nbatches)
      ntt = 1

      t1 = second()
C$       t1 = omp_get_wtime()  
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
c$OMP& SCHEDULE(DYNAMIC)      
      do i=1,ntarg_n

        call cquadints_adap(eps,intype,norder_p,type,npols,ntt,
     1     xyztarg_near(1,i),nquadmax,h3d_slp,dpars,zpars,ipars,
     2     nqorder,slp_near(1,i))

        call cquadints_adap(eps,intype,norder_p,type,npols,ntt,
     1     xyztarg_near(1,i),nquadmax,h3d_dlp,dpars,zpars,ipars,
     2     nqorder,dlp_near(1,i))
      enddo
C$OMP END PARALLEL DO      
      t2 = second()
C$       t2 = omp_get_wtime()      

cc      call prin2('time taken in evaluating near=*',t2-t1,1)


      call squarearbq_pts(nqorderf,nnodes)

      nu = 3
      nqpts = nnodes*nu*nu
      allocate(qnodes(2,nqpts),qwts(nqpts))

cc      call prinf('ntarg_f=*',ntarg_f,1)
cc      call prinf('ntarg_n=*',ntarg_n,1)

      call gen_xg_uniftree_nodes(nqorderf,nnodes,nu,nqpts,qnodes,qwts)

      call cquadints_wnodes(norder_p,type,npols,ntarg_f,xyztarg_far,
     1       h3d_slp,dpars,zpars,ipars,nqpts,qnodes,qwts,slp_far)

      call cquadints_wnodes(norder_p,type,npols,ntarg_f,xyztarg_far,
     1       h3d_dlp,dpars,zpars,ipars,nqpts,qnodes,qwts,dlp_far)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
      do i=1,ntarg_n
        do j=1,npols
          slp_pots(j,inear_ind(i)) = slp_near(j,i)
          dlp_pots(j,inear_ind(i)) = dlp_near(j,i)
        enddo
      enddo
C$OMP END PARALLEL DO


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
      do i=1,ntarg_f
        do j=1,npols
          slp_pots(j,ifar_ind(i)) = slp_far(j,i)
          dlp_pots(j,ifar_ind(i)) = dlp_far(j,i)
        enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
