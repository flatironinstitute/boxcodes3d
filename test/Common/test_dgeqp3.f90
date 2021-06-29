      implicit real *8 (a-h,o-z)
      real *8, allocatable :: x(:,:)
      real *8 done
      integer, allocatable :: ipt2depth(:,:),idepth2pt(:,:,:)
      real *8, allocatable :: fcoefs(:,:),fcoefsex(:,:)
      real *8, allocatable :: frhs(:,:)
      real *8, allocatable :: amat(:,:),pols(:),amat_copy(:,:)
      real *8, allocatable :: tau(:),work(:),work2(:)
      complex *16, allocatable :: zcoefs(:,:),zrhs(:,:),zcoefsex(:,:)
      integer, allocatable :: jpvt(:)
      complex *16 ima
      character *1 type

      data ima/(0.0d0,1.0d0)/

      call prini(6,13)
      done = 1.0d0

      ndeg = 7
      n = ndeg+1
      type = 't'
      call legetens_npol_3d(ndeg,type,npol3)
      

      nptmax = n**3
      allocate(x(3,nptmax))
      
      allocate(ipt2depth(3,nptmax),idepth2pt(n,n,n))

      ifloor = 1.75d0*npol3
      call fakepolya3d(n,ifloor,npt,x,ipt2depth,idepth2pt)
      itype = 0
      call prinf('npt=*',npt,1)
      call prin2('rat=*',(npt+0.0d0)/(nptmax+0.0d0),1)
      print *, "npt=",npt

      

      allocate(amat(npt,npol3),pols(npol3))
      allocate(amat_copy(npt,npol3))
      nrhs = 3
      allocate(fcoefs(npol3,nrhs),frhs(npt,nrhs))
      allocate(fcoefsex(npol3,nrhs))

      allocate(zcoefs(npol3,nrhs),zrhs(npt,nrhs))
      allocate(zcoefsex(npol3,nrhs))

      do j=1,nrhs
        do i=1,npol3
          fcoefsex(i,j) = hkrand(0) 
          zcoefsex(i,j) = hkrand(0) + ima*hkrand(0)
        enddo
        do i=1,npt
          frhs(i,j) = 0
          zrhs(i,j) = 0
        enddo
      enddo


      do i=1,npt
        call legetens_pols_3d(x(1,i),ndeg,type,pols)
        do j=1,npol3
          amat(i,j) = pols(j)
          do l=1,nrhs
            frhs(i,l) = frhs(i,l) + fcoefsex(j,l)*pols(j)
            zrhs(i,l) = zrhs(i,l) + zcoefsex(j,l)*pols(j)
          enddo
        enddo
      enddo


      allocate(tau(npol3),jpvt(npol3))
      call get_qrdecomp(npt,npol3,amat,amat_copy,jpvt,tau)

      call dqrsolv(npt,npol3,amat_copy,jpvt,tau,nrhs,frhs,fcoefs)
      call zqrsolv(npt,npol3,amat_copy,jpvt,tau,nrhs,zrhs,zcoefs)

      erra = 0.0d0
      ra = 0.0d0

      errz = 0.0d0
      rz = 0.0d0
      do j=1,nrhs
        do i=1,npol3
          erra = erra + (fcoefsex(i,j) - fcoefs(i,j))**2
          ra = ra + fcoefsex(i,j)**2
          errz = errz + abs(zcoefsex(i,j) - zcoefs(i,j))**2
          rz = rz + abs(zcoefsex(i,j))**2
        enddo
      enddo

      erra = sqrt(erra/ra)
      errz = sqrt(errz/rz)
      call prin2('error in coefs=*',erra,1)
      call prin2('error in z coefs=*',errz,1)



      

      stop
      end