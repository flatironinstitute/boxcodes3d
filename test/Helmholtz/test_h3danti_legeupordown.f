      implicit real *8 (a-h,o-z)
      complex *16 zk, im, zero, one, alpha, beta

      real *8, allocatable :: dlap(:,:), dlapbig(:,:)
      
      complex *16, allocatable :: ahmatup(:,:), ahmatdown(:,:)
      complex *16, allocatable :: hmat(:,:), hmatbig(:,:)
      complex *16, allocatable :: resup(:,:), resdown(:,:)


      integer, allocatable :: ind2p(:,:), ind2pbig(:,:)
      
      real *8, allocatable :: sumsup(:), sumsdown(:)
      
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /

      character type
      
      call prini(6,13)


      ndeg = 8
      nup = 10

      ndegbig = ndeg + nup*2
       
      tol = 1d-15

      type = 't'

      r = 1.21d0
      zk = r/sqrt(2.0d0) + im*r/sqrt(2.0d0)
c      zk = r
c      zk = r*im
      
      
      call legetens_npol_3d(ndeg,type,npol)
      call legetens_npol_3d(ndegbig,type,npolbig)

      call prinf('npol *',npol,1)
      call prinf('npolbig *',npolbig,1)

      allocate(ind2p(3,npol),ind2pbig(3,npolbig))

      call legetens_ind2pow_3d(ndeg,type,ind2p)
      call legetens_ind2pow_3d(ndegbig,type,ind2pbig)      
      
      allocate(ahmatup(npolbig,npol),ahmatdown(npol,npol))

      call cpu_time(t1)
      call h3danti_legeup(ndeg,nup,type,zk,ahmatup,npolbig)
      call cpu_time(t2)

      write(*,*) 'time upward recurrence ', t2-t1

      call cpu_time(t1)
      call h3danti_legedown(ndeg,type,zk,ahmatdown,npol)
      call cpu_time(t2)

      write(*,*) 'time downward recurrence ', t2-t1
      
      allocate(dlapbig(npolbig,npolbig),hmatbig(npolbig,npolbig),
     1     hmat(npol,npol),dlap(npol,npol),resup(npolbig,npol),
     2     resdown(npol,npol))
      call legetens_lapmat_3d(ndegbig,ndegbig,type,dlapbig,npolbig)      
      call legetens_lapmat_3d(ndeg,ndeg,type,dlap,npol)
      
      do i= 1,npol
         do j = 1,npol
            hmat(j,i) = dlap(j,i)
            if (j .eq. i) hmat(j,i) = hmat(j,i) + zk**2
         enddo
      enddo
      
      do i= 1,npolbig
         do j = 1,npolbig
            hmatbig(j,i) = dlapbig(j,i)
            if (j .eq. i) hmatbig(j,i) = hmatbig(j,i) + zk**2
         enddo
      enddo
      
      alpha = one
      beta = zero
      call zgemm('N','N',npolbig,npol,npolbig,alpha,hmatbig,npolbig,
     1     ahmatup,npolbig,beta,resup,npolbig)
      call zgemm('N','N',npol,npol,npol,alpha,hmat,npol,
     1     ahmatdown,npol,beta,resdown,npol)
      
c      call prin2('res 1 *',res1,npol*npol*2)
c      call prin2('res 2 *',res2,npol*npol*2)
      
      sum1 = 0
      
      sumrel = 0

      allocate(sumsup(npol),sumsdown(npol))
      
      do i = 1,npol
         sumsup(i) = 0
         sumsdown(i) = 0
         iin = ind2p(1,i)
         jin = ind2p(2,i)
         kin = ind2p(3,i)
         do j = 1,npolbig
            ibig = ind2pbig(1,j)
            jbig = ind2pbig(2,j)
            kbig = ind2pbig(3,j)
            if (iin .eq. ibig .and. jin .eq. jbig
     1           .and. kin .eq. kbig) then
               sumsup(i) = sumsup(i) + abs(one-resup(j,i))**2
            else
               sumsup(i) = sumsup(i) + abs(resup(j,i))**2
            endif
         enddo

         do j = 1,npol
            if (j .eq. i) then
               sumsdown(i) = sumsdown(i) + abs(one-resdown(j,i))**2
            else
               sumsdown(i) = sumsdown(i) + abs(resdown(j,i))**2
            endif
         enddo

         sumsdown(i) = sqrt(sumsdown(i))
         sumsup(i) = sqrt(sumsup(i))
      enddo

      call prin2('zk *',zk,2)

c      call prin2('residual recurrence *',sqrt(sum2/npol),1)


c      call prin2('res up rec by col *',sumsup,npol)
c      call prin2('res down rec by col *',sumsdown,npol)

      summax = 0.0d0
      do i = 1,npol
         write(*,*) sumsup(i), sumsdown(i), min(sumsup(i),sumsdown(i))
         summax = max(summax,min(sumsup(i),sumsdown(i)))
      enddo

      call prin2('summax *',summax,1)
      call prinf('nup *',nup,1)
      
      stop
      end
