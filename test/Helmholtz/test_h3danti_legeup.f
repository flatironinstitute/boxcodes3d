c$    use omp_lib
      implicit real *8 (a-h,o-z)
      complex *16 zk, im, zero, one, alpha, beta

      real *8, allocatable :: dlap(:,:)
      
      complex *16, allocatable :: ahmat1(:,:), ahmat2(:,:), ahmat3(:,:)
      complex *16, allocatable :: hmat(:,:)
      complex *16, allocatable :: res1(:,:), res2(:,:), res3(:,:)


      integer, allocatable :: ind2p(:,:), ind2pout(:,:)
      
      real *8, allocatable :: sums1(:), sums2(:), sums3(:)
      
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /

      character type
      
      call prini(6,13)


      ndeg = 8
      nup = 12

      ndegout = ndeg + nup*2
      
      tol = 1d-15

      type = 't'

      zk = sqrt(1.5d0) + im*sqrt(1.5d0)
      zk = 1.5d0
      
      call legetens_npol_3d(ndeg,type,npol)
      call legetens_npol_3d(ndegout,type,npolout)

      call prinf('npol *',npol,1)
      call prinf('npolout *',npolout,1)

      allocate(ind2p(3,npol),ind2pout(3,npolout))

      call legetens_ind2pow_3d(ndeg,type,ind2p)
      call legetens_ind2pow_3d(ndegout,type,ind2pout)      
      
      allocate(ahmat1(npolout,npol),ahmat2(npolout,npol))

      call cpu_time(t1)
c$    t1 = omp_get_wtime()      
      call h3danti_legeup(ndeg,nup,type,zk,ahmat1,npolout)
      call cpu_time(t2)
c$    t2 = omp_get_wtime()            

      write(*,*) 'time upward recurrence ', t2-t1
      
      call cpu_time(t1)
c$    t1 = omp_get_wtime()      
      call h3danti_legeupfast(ndeg,nup,type,zk,ahmat2,npolout)
      call cpu_time(t2)
c$    t2 = omp_get_wtime()            

      write(*,*) 'time upward recurrence (fast?)', t2-t1

      allocate(dlap(npolout,npolout),hmat(npolout,npolout),
     1     res1(npolout,npol),res2(npolout,npol))
      call legetens_lapmat_3d(ndegout,ndegout,type,dlap,npolout)      
      
      do i= 1,npolout
         do j = 1,npolout
            hmat(j,i) = dlap(j,i)
            if (j .eq. i) hmat(j,i) = hmat(j,i) + zk**2
         enddo
      enddo
      
      alpha = one
      beta = zero
      call zgemm('N','N',npolout,npol,npolout,alpha,hmat,npolout,
     1     ahmat1,npolout,beta,res1,npolout)
      call zgemm('N','N',npolout,npol,npolout,alpha,hmat,npolout,
     1     ahmat2,npolout,beta,res2,npolout)
      
c      call prin2('res 1 *',res1,npol*npol*2)
c      call prin2('res 2 *',res2,npol*npol*2)
      
      sum1 = 0
      
      sumrel = 0

      allocate(sums1(npol),sums2(npol))
      
      do i = 1,npol
         sums1(i) = 0
         sums2(i) = 0
         iin = ind2p(1,i)
         jin = ind2p(2,i)
         kin = ind2p(3,i)
         do j = 1,npolout
            iout = ind2pout(1,j)
            jout = ind2pout(2,j)
            kout = ind2pout(3,j)
            if (iin .eq. iout .and. jin .eq. jout
     1           .and. kin .eq. kout) then
               sums1(i) = sums1(i) + abs(one-res1(j,i))**2
               sums2(i) = sums2(i) + abs(one-res2(j,i))**2
            else
               sums1(i) = sums1(i) + abs(res1(j,i))**2
               sums2(i) = sums2(i) + abs(res2(j,i))**2               
            endif
         enddo

         sums1(i) = sqrt(sums1(i))
         sums2(i) = sqrt(sums2(i))
      enddo

      call prin2('zk *',zk,2)

c      call prin2('residual recurrence *',sqrt(sum2/npol),1)


c      call prin2('res rec by col *',sums1,npol)
c      call prin2('res rec by col *',sums2,npol)

      smax1 = 0
      smax2 = 0
      do i = 1,npol
         smax1 = max(smax1,sums1(i))
         smax2 = max(smax2,sums2(i))
      enddo


      call prin2('max err, legeup *',smax1,1)
      call prin2('max err, legeupfast *',smax2,1)      
      

      stop
      end
