      implicit real *8 (a-h,o-z)
      integer p,q
      real *8, allocatable :: fvalsp(:),fvalsq(:)
      real *8, allocatable :: fvalsqex(:),vmat(:,:)
      real *8, allocatable :: fcoefst(:),fcoefsf(:)
      real *8 pols(1000)
      real *8, allocatable :: xp(:,:),umatpf(:,:),vmatpf(:,:),wp(:)
      real *8, allocatable :: umatpt(:,:),vmatpt(:,:)
      real *8, allocatable :: xq(:,:),umatqf(:,:),vmatqf(:,:),wq(:)
      real *8, allocatable :: umatqt(:,:),vmatqt(:,:)
      real *8, allocatable :: xq1(:),wq1,uq1,vq1
      character *1 ttype

      call prini(6,13)

      ttype = 't'

      p = 15
      q = 18

      call legetens_npol_3d(p-1,ttype,npt)
      call legetens_npol_3d(q-1,ttype,nqt)

      ttype = 'f'


      call legetens_npol_3d(p-1,ttype,npf)
      call legetens_npol_3d(q-1,ttype,nqf)
      call prinf('npf=*',npf,1)
      call prinf('nqf=*',nqf,1)
      call prinf('npt=*',npt,1)
      call prinf('nqt=*',nqt,1)
      call prinf('p=*',p,1)
      call prinf('q=*',q,1)


      allocate(fvalsp(npf),fvalsq(nqf))
      allocate(fvalsqex(nqf),vmat(p,q))
      allocate(fcoefst(npt),fcoefsf(npf))

      allocate(xp(3,npf),wp(npf),umatpf(npf,npf),vmatpf(npf,npf)) 
      allocate(umatpt(npt,npf),vmatpt(npf,npt))

      allocate(xq(3,nqf),wq(nqf),umatqf(nqf,nqf),vmatqf(nqf,nqf))
      allocate(umatqt(nqt,nqf),vmatqt(nqf,nqt))

      itype = 2
      ttype = 't'
      call legetens_exps_3d(itype,p,ttype,xp,umatpt,npt,
     1  vmatpt,npf,wp)

      ttype = 'f'
      call legetens_exps_3d(itype,p,ttype,xp,umatpf,npf,
     1  vmatpf,npf,wp)

      

      ttype = 't'
      call legetens_exps_3d(itype,q,ttype,xq,umatqt,nqt,
     1  vmatqt,nqf,wq)

      ttype = 'f'
      call legetens_exps_3d(itype,q,ttype,xq,umatqf,nqf,
     1  vmatqf,nqf,wq)


      do i=1,npf
        call fval(xp(1,i),fvalsp(i))
      enddo

      do i=1,npf
        fcoefsf(i) = 0
        do j=1,npf
          fcoefsf(i) = fcoefsf(i) + umatpf(i,j)*fvalsp(j)
        enddo
      enddo


      do i=1,npt
        fcoefst(i) = 0
        do j=1,npf
          fcoefst(i) = fcoefst(i) + umatpt(i,j)*fvalsp(j)
        enddo
      enddo


      call prin2('fvals=*',fvalsp,24)
      call prin2('fcoefsf=*',fcoefsf,24)
      call prin2('fcoefst=*',fcoefst,24)

      
      do i=1,nqf
        call fval(xq(1,i),fvalsqex(i))
      enddo

      call prin2('fvalsqex=*',fvalsqex,24)

c
c  compute vmat
c
      itype = 0
      allocate(xq1(q))
      call prinf('q=*',q,1)
      call legeexps(itype,q,xq1,uq1,vq1,wq1)
      call prin2('xq1=*',xq1,q)

      do i=1,q
        call legepols(xq1(i),p-1,vmat(1,i))
      enddo


      ttype = 'f'
      call legetens_eval_3d_pq(p,ttype,q,fcoefsf,npf,fvalsq,nqf,vmat)

      call prin2('fvalsqex=*',fvalsqex,24)
      call prin2('fvalsq=*',fvalsq,24)

      erra = 0
      ra = 0
      do i=1,nqf
        erra = erra + abs(fvalsq(i)-fvalsqex(i))**2*wq(i)
        ra = ra + fvalsqex(i)**2*wq(i)
      enddo

      erra = sqrt(erra/ra)
      call prin2('error f=*',erra,1)


      ttype = 't'
      call legetens_eval_3d_pq(p,ttype,q,fcoefst,npt,fvalsq,nqf,vmat)

      call prin2('fvalsqex=*',fvalsqex,24)
      call prin2('fvalsq=*',fvalsq,24)

      erra = 0
      ra = 0
      do i=1,nqf
        erra = erra + abs(fvalsq(i)-fvalsqex(i))**2*wq(i)
        ra = ra + fvalsqex(i)**2*wq(i)
      enddo

      erra = sqrt(erra/ra)
      call prin2('error t=*',erra,1)

      stop
      end




      subroutine fval(xyz,f)
      implicit real *8 (a-h,o-z)
      real *8 xyz(3)

      x = xyz(1)
      y = xyz(2) 
      z = xyz(3) 
      f = sin(x*y + y*z/2 + z*x*1.1d0) 

      return
      end
