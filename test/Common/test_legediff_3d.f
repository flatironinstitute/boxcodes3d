      implicit real *8 (a-h,o-z)
      real *8, allocatable :: x(:,:),w(:),u(:,:),v(:,:),dmat(:,:)
      real *8, allocatable :: f(:),fcoef(:),fderxcoef(:)
      real *8, allocatable :: fderycoef(:),fderzcoef(:),pol(:)
      real *8 xpt(3)
      character ttype

      call prini(6,13)

      itype = 2
      ndeg = 21
      ttype = 't'
      npt = (ndeg+1)**3


      call legetens_npol_3d(ndeg,ttype,npol3)
      allocate(x(3,npt),w(npt),v(npt,npol3),u(npol3,npt))
      allocate(dmat(ndeg+1,ndeg+1))
      call legetens_exps_3d(itype,ndeg+1,ttype,x,u,npol3,v,npt,
     1  w)
      call legecoeff_dmat(ndeg,dmat,ndeg+1)
      allocate(f(npt),fcoef(npol3),fderxcoef(npol3))
      allocate(fderycoef(npol3),fderzcoef(npol3),pol(npol3))
      do i=1,npt
        call fval(x(1,i),x(2,i),x(3,i),f(i))
      enddo

      fcoef = 0
      do i=1,npol3
        fcoef(i) = 0
        do j=1,npt
          fcoef(i) = fcoef(i) + u(i,j)*f(j)
        enddo
      enddo

      call prin2('fcoef=*',fcoef,24)

      xpt(1) = hkrand(0)*2-1.0d0
      xpt(2) = hkrand(0)*2-1.0d0
      xpt(3) = hkrand(0)*2-1.0d0
      call legetens_pols_3d(xpt,ndeg,ttype,pol)

      call fval(xpt(1),xpt(2),xpt(3),fex)
      ftest = 0
      do i=1,npol3 
        ftest = ftest + pol(i)*fcoef(i)
      enddo
      call prin2('ftest=*',ftest,1)
      call prin2('fex=*',fex,1)
      call prin2('error =*', abs(ftest-fex),1)
c
c   nwo test derivative
c
c
      idir = 1
      call legediff_3d(ndeg,ttype,fcoef,idir,dmat,fderxcoef)

      
      call fder(xpt(1),xpt(2),xpt(3),idir,fex)
      ftest = 0
      do i=1,npol3 
        ftest = ftest + pol(i)*fderxcoef(i)
      enddo
      call prin2('ftest=*',ftest,1)
      call prin2('fex=*',fex,1)
      call prin2('error dx=*', abs(ftest-fex),1)
c
c
c
      idir = 2
      call legediff_3d(ndeg,ttype,fcoef,idir,dmat,fderycoef)

      
      call fder(xpt(1),xpt(2),xpt(3),idir,fex)
      ftest = 0
      do i=1,npol3 
        ftest = ftest + pol(i)*fderycoef(i)
      enddo
      call prin2('ftest=*',ftest,1)
      call prin2('fex=*',fex,1)
      call prin2('error dy=*', abs(ftest-fex),1)
c
c
c
c
      idir = 3
      call legediff_3d(ndeg,ttype,fcoef,idir,dmat,fderzcoef)

      
      call fder(xpt(1),xpt(2),xpt(3),idir,fex)
      ftest = 0
      do i=1,npol3 
        ftest = ftest + pol(i)*fderzcoef(i)
      enddo
      call prin2('ftest=*',ftest,1)
      call prin2('fex=*',fex,1)
      call prin2('error dz=*', abs(ftest-fex),1)

      


      return
      end
c
c
c
c
c
      subroutine fval(x,y,z,f)
      implicit real *8 (a-h,o-z)

      f = sin(x*y + 2.1d0*y*z + 3.1d0*z*x)

      return
      end
c
c
c
c
c
      subroutine fder(x,y,z,idir,f)
      implicit real *8 (a-h,o-z)

      if(idir.eq.1) then
        f = cos(x*y + 2.1d0*y*z + 3.1d0*z*x)*(y + 3.1d0*z)
      else if(idir.eq.2) then
        f = cos(x*y + 2.1d0*y*z + 3.1d0*z*x)*(x + 2.1d0*z)
      else
        f = cos(x*y + 2.1d0*y*z + 3.1d0*z*x)*(2.1d0*y + 3.1d0*x)
      endif

      return
      end
c
c
c
c
