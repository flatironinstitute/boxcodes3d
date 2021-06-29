      subroutine get_qrdecomp(m,n,a,aqr,jpvt,tau)
      implicit real *8 (a-h,o-z)
      real *8 a(m,n),aqr(m,n),tau(n)
      integer jpvt(n)
      real *8, allocatable :: work(:)
      lwork = 100*(m+n)
      allocate(work(lwork))

      aqr = a

      call dgeqp3(m,n,aqr,m,jpvt,tau,work,lwork,info)

      return
      end
!
!
!
!
!
!
      subroutine dqrsolv(m,n,aqr,jpvt,tau,nrhs,frhs,fsol)
      implicit real *8 (a-h,o-z)
      real *8 aqr(m,n),frhs(m,nrhs),fsol(n,nrhs)
      real *8 tau(n),done
      integer jpvt(n)
      real *8, allocatable :: ftmp(:,:),work(:)
      done = 1.0d0

      lwork = 10*(m+n)
      allocate(ftmp(m,nrhs),work(lwork))
      ftmp = frhs
      info = 0
      call dormqr('l','t',m,nrhs,n,aqr,m,tau,ftmp,m, &
        work,lwork,info)
      
      call dtrsm('l','u','n','n',n,nrhs,done,aqr,m,ftmp,m)
      ftmp(n+1:m,1:nrhs) = 0.0d0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(work,j)      
      do j=1,nrhs
        work(jpvt(1:n)) = ftmp(1:n,j)
        fsol(1:n,j) = work(1:n)
      enddo
!$OMP END PARALLEL DO      
      
      return
      end
!
!
!
!
!

      subroutine zqrsolv(m,n,aqr,jpvt,tau,nrhs,zrhs,zsol)
      implicit real *8 (a-h,o-z)
      real *8 aqr(m,n)
      complex *16 zrhs(m,nrhs),zsol(n,nrhs)
      real *8 tau(n),done
      integer jpvt(n)
      real *8, allocatable :: ftmp(:,:),work(:)
      complex *16 ima
      data ima/(0.0d0,1.0d0)/
      done = 1.0d0

      lwork = 10*(m+n)

      nrhsuse = 2*nrhs
      allocate(ftmp(m,nrhsuse),work(lwork))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
      do i=1,nrhs
        ftmp(1:m,i) = real(zrhs(1:m,i))
      enddo
!$OMP END PARALLEL DO     

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
      do i=1,nrhs
        ftmp(1:m,i+nrhs) = imag(zrhs(1:m,i))
      enddo
!$OMP END PARALLEL DO     

      info = 0
      call dormqr('l','t',m,nrhsuse,n,aqr,m,tau,ftmp,m, &
        work,lwork,info)
      
      call dtrsm('l','u','n','n',n,nrhsuse,done,aqr,m,ftmp,m)
      ftmp(n+1:m,1:nrhsuse) = 0.0d0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(work,j)      
      do j=1,nrhs
        work(jpvt(1:n)) = ftmp(1:n,j)
        zsol(1:n,j) = work(1:n)
      enddo
!$OMP END PARALLEL DO      
      
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(work,j)      
      do j=1,nrhs
        work(jpvt(1:n)) = ftmp(1:n,j+nrhs)
        zsol(1:n,j) = zsol(1:n,j) + ima*work(1:n)
      enddo
!$OMP END PARALLEL DO      
      
      return
      end
     

