      implicit real *8 (a-h,o-z)
      real *8 tts(3,100,2),errm_rel1(100,2),errm_rel2(100,2)
      real *8 errm_abs(100,2)
      real *8 tt3d(3),err_print(3)
      real *8, allocatable :: x(:,:)
      real *8, allocatable :: pmat(:,:),pmat_qr(:,:)
      real *8, allocatable :: pmat_tau(:),pols(:)
      integer, allocatable :: pmat_jpvt(:)
      complex *16 zk, im, zero, one
      complex *16, allocatable :: tab_ref(:,:),tab(:,:),tab2(:,:)
      complex *16, allocatable :: tab_ref_coefs(:,:),tab2_coefs(:,:)
      integer, allocatable :: ipt2depth(:,:),idepth2pt(:,:,:)
      logical ifsphere
      character type
      character *2 arg_comm

      character *100 fname
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /           
      external h3d_vslp
      
      call prini(6,13)
      

      call get_command_argument(1,arg_comm)
      read(arg_comm,*) iprec

      call get_command_argument(2,arg_comm)
      read(arg_comm,*) ndeg

      n = ndeg + 1

      nup = 8

      type = 't'
      call legetens_npol_3d(ndeg,type,npol3)

      if(iprec.eq.0) tol = 1.0d-2
      if(iprec.eq.1) tol = 1.0d-3
      if(iprec.eq.2) tol = 1.0d-6
      if(iprec.eq.3) tol = 1.0d-9
      if(iprec.eq.4) tol = 1.0d-12

      nptmax = n**3
      allocate(x(3,nptmax))
      
      allocate(ipt2depth(3,nptmax),idepth2pt(n,n,n))

      ifloor = 1.75d0*npol3
      call fakepolya3d(n,ifloor,npt,x,ipt2depth,idepth2pt)
      itype = 0
      call prinf('npt=*',npt,1)
      call prin2('rat=*',(npt+0.0d0)/(nptmax+0.0d0),1)
      call prin2('x=*',x,3*npt)


      allocate(pmat(npt,npol3),pmat_qr(npt,npol3),pols(npol3))
      allocate(pmat_tau(npol3),pmat_jpvt(npol3))

      do i=1,npt
        call legetens_pols_3d(x(1,i),ndeg,type,pols)
        do j=1,npol3
          pmat(i,j) = pols(j)
        enddo
      enddo

      call get_qrdecomp(npt,npol3,pmat,pmat_qr,pmat_jpvt,pmat_tau)

      
cc      call legetens_exps_3d(itype,n,'t',x,u,1,v,1,w)
      

      ntarg0 = 10*npt

      npol0 = 10*npol3

      allocate(tab_ref(ntarg0,npol3),tab(npol0,npol3))
      allocate(tab2(ntarg0,npol3))
      allocate(tab_ref_coefs(npol0,npol3),tab2_coefs(npol0,npol3))

c
c
c
c
      ifgen = 0
      izk = 2
      eps_exact = 1.0d-11
      ifwrite = 0
      ifread = 1
      if(ifgen.eq.1) then
         if(izk.eq.1) zk = 1.6d0
         if(izk.eq.2) zk = 0.16d0
         if(izk.eq.3) zk = im*1.6d0
         if(izk.eq.4) zk = 5.0d0
         print *, "here"
         print *, "zk=",zk
         call h3dtab_ref_brute(ndeg,eps_exact,x,npt,tab_ref,ntarg0,
     1      h3d_vslp,dpars,zk,ipars)
         if(ifwrite.eq.1) then
          write(fname,'(a,i2.2,a,i1,a)') 'tabref_data/tabref_',
     1       n,'_izk_',izk,'new.dat'
          open(unit=33,file=trim(fname),action='readwrite',
     1       form='unformatted',access='stream')
         write(unit=33) tab_ref
         close(33)
        endif
      endif


      

      if(igen.ne.1) then

        do izk=2,2

          if(izk.eq.1) zk = 1.6d0
          if(izk.eq.2) zk = 0.16d0
          if(izk.eq.3) zk = im*1.6d0
          if(izk.eq.4) zk = 5.0d0
          if(ifread.eq.1) then
            write(fname,'(a,i2.2,a,i1,a)') 'tabref_data/tabref_',n,
     1         '_izk_',izk,'new.dat'
            open(unit=33,file=trim(fname),action='readwrite',
     1         form='unformatted',access='stream')
            read(unit=33) tab_ref

            print *, "done reading"

            close(33)
          endif
c
c
c  convert table reference values to coefs
c
          do i=1,10
            istart = (i-1)*npt+1
            iend = i*npt

            istart2 = (i-1)*npol3+1
            iend2 = i*npol3
            print *,i,istart,iend,istart2,iend2
            call zqrsolv(npt,npol3,pmat_qr,pmat_jpvt,pmat_tau,npol3,
     1        tab_ref(istart:iend,1:npol3),
     2        tab_ref_coefs(istart2:iend2,1:npol3))
          enddo


          ndeg2 = ndeg + 2*nup 
          iflg = 1

          call h3dtabp_ref(ndeg,zk,tol,npt,x,tab,npol0,nup,
     1      iflg,pmat_qr,npol3,pmat_jpvt,pmat_tau,
     2      tts(1,izk,iflg),tts(2,izk,iflg),tts(3,izk,iflg))
          
          
          iflg = 2

          call h3dtabp_ref(ndeg,zk,tol,npt,x,tab,npol0,nup,
     1      iflg,pmat_qr,npol3,pmat_jpvt,pmat_tau,
     2      tts(1,izk,iflg),tts(2,izk,iflg),tts(3,izk,iflg))
          print *, "Done with 2d adap quad"
        
          call prin2('tts iflg1=*',tts(1,izk,1),3)
          call prin2('tts iflg2=*',tts(1,izk,2),3)

          call cpu_time(t1)
C$           t1 = omp_get_wtime()          

          call h3dtab_ref_brute(ndeg,tol,x,npt,tab2,ntarg0,
     1        h3d_vslp,dpars,zk,ipars)
          do i=1,10
            istart = (i-1)*npt+1
            iend = i*npt

            istart2 = (i-1)*npol3+1
            iend2 = i*npol3
            call zqrsolv(npt,npol3,pmat_qr,pmat_jpvt,pmat_tau,npol3,
     1        tab2(istart:iend,1:npol3),
     2        tab2_coefs(istart2:iend2,1:npol3))
          enddo
          call cpu_time(t2)
C$           t2 = omp_get_wtime()          
cc    


          print *, "done with 3d adap quad"
          errmax1 = 0
          errmax2 = 0
          errmaxa = 0
          rmax = 0
          do i=1,npol3
            do j=1,npol0
             ra = max(abs(tab_ref_coefs(j,i)),1.0d0)
             erra = abs(tab(j,i)-tab_ref_coefs(j,i))/ra
             if(erra.gt.errmax1) errmax1 = erra

             
             erra = abs(tab(j,i)-tab_ref_coefs(j,i))/
     1          abs(tab_ref_coefs(j,i))
             if(erra.gt.errmax2) errmax2 = erra


             erra = abs(tab(j,i)-tab_ref_coefs(j,i))
             if(erra.gt.errmaxa) errmaxa = erra
             ra = abs(tab_ref_coefs(j,i))
             if(ra.gt.rmax) rmax = ra
            enddo
          enddo

          errm_rel1(izk,1) = errmax1
          errm_rel2(izk,1) = errmax2
          errm_abs(izk,1) = errmaxa/rmax
          err_print(1) = errmax1
          err_print(2) = errmax2
          err_print(3) = errmaxa/rmax

          call prin2(' *',i,0)
          call prin2(' *',i,0)
          call prin2('2d anti helmholtzian*',i,0)
          call prin2(' *',i,0)
          call prin2(' *',i,0)
          call prin2('max errors=*',err_print,3)
          call prin2('time taken=*',tts(:,izk,:),6)


          errmax1 = 0
          errmax2 = 0
          errmaxa = 0
          rmax = 0
          do i=1,npol3
            do j=1,npol0
             ra = max(abs(tab_ref_coefs(j,i)),1.0d0)
             erra = abs(tab2_coefs(j,i)-tab_ref_coefs(j,i))/ra
             if(erra.gt.errmax1) errmax1 = erra

             
             erra = abs(tab2_coefs(j,i)-tab_ref_coefs(j,i))/
     1          abs(tab_ref_coefs(j,i))
             if(erra.gt.errmax2) errmax2 = erra


             erra = abs(tab2_coefs(j,i)-tab_ref_coefs(j,i))
             if(erra.gt.errmaxa) errmaxa = erra
             ra = abs(tab_ref_coefs(j,i))
             if(ra.gt.rmax) rmax = ra
            enddo
          enddo

          call prin2(' *',i,0)
          call prin2(' *',i,0)
          call prin2('3d adap integration comp*',i,0)
          call prin2(' *',i,0)
          call prin2(' *',i,0)
          errm_rel1(izk,2) = errmax1
          errm_rel2(izk,2) = errmax2
          errm_abs(izk,2) = errmaxa/rmax
          err_print(1) = errmax1
          err_print(2) = errmax2
          err_print(3) = errmaxa/rmax
          call prin2('max error=*',err_print,3)
          tt3d(izk) = t2-t1
          call prin2('time taken=*',tt3d(izk),1)

 1111     continue          
        enddo
        call prin2(' *',i,0)
        call prin2(' *',i,0)
        call prin2('results summary:*',i,0)
        call prin2(' *',i,0)
        call prin2(' *',i,0)

        call prin2('max errors relative 1=*',errm_rel1(2,1:2),2)
        call prin2('max errors relative 2=*',errm_rel2(2,1:2),2)
        call prin2('max errors absolute=*',errm_abs(2,1:2),2)
        call prin2('time taken 2d=*',tts(:,2,:),6)
        call prin2('time taken 3d adap=*',tt3d(2),1)
        stop
      endif



      
      stop
      end
c
c
c
c
c     
c
