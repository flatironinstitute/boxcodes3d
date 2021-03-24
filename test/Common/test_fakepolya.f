      implicit real *8 (a-h,o-z)

      integer, allocatable :: ipt2depth(:,:), idepth2pt(:,:,:)
      real *8, allocatable :: pts(:,:), amat(:,:)
      
      character type
      character(len=100) fname,fform
      call prini(6,13)

c
      write(*,*) 'enter polynomial order '
      read(*,*) nord
      write(*,*) 'enter thinning ratio (npts used/ npts optimal)'
      read(*,*) dfac

      fform = '(A16,I0.3,A4)'
      write(fname,fform) "output/fakepolya", nord, ".txt"
      write(*,*) trim(fname)

      open(unit=101,file=trim(fname),access='append')
      
      ng = nord
      ndeg = nord-1
      type = 't'
      nptmax = ng**3
      
      call legetens_npol_3d(ndeg,type,npol)
      
      allocate(pts(3,nptmax),ipt2depth(3,nptmax),
     1     idepth2pt(ng,ng,ng))

      ifloor = dfac*npol+0.1d0
      call fakepolya3d(ng,ifloor,npt,pts,ipt2depth,idepth2pt)
      call ichecksymms(idepth2pt,ipt2depth,ng,npt,ipass1,ipass2,
     1     ipass3)      

      write(*,*) 'npt found, npt requested ', npt, ifloor
      write(*,*) 'symm check (want all 1s)', ipass1, ipass2, ipass3      

      mm = npol
      nn = npt
      allocate(amat(mm,nn))
      do i = 1,npt
         call legetens_pols_3d(pts(1,i),ndeg,type,amat(1,i))
      enddo
      
      call dcomputecond(dkappa,amat,mm,nn)

      write(*,*) 'condition number (and scaled)', dkappa,
     1     dkappa/sqrt(npol*1d0)

      write(101,*) nord, dfac, npt, npol, dkappa, dkappa/sqrt(npol*1d0)

      close(101)
      
      stop
      end
c
c
c


      subroutine dcomputecond(dkappa,amat,m,n)
      implicit real *8 (a-h,o-z)
      real *8 amat(m,n), dkappa

      real *8 u,vt
      real *8, allocatable :: s(:), work(:)
      character jobu,jobvt

      lda = m

      mn = min(m,n)
      mnbig = max(m,n)
      lwork = 5*mn + mnbig
      ldu = 1
      ldvt = 1
      allocate(s(mn),work(lwork))

      jobu='N'
      jobvt='N'
      call dgesvd(jobu,jobvt,m,n,amat,lda,s,u,ldu,vt,ldvt,work,lwork,
     1     info)

      write(*,*) 'info ', info

      dkappa = s(1)/s(mn)

      return
      end
      
