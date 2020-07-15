c$    use omp_lib
      implicit real *8 (a-h,o-z)
      complex *16 zk, im, zero, one, alpha, beta

      real *8, allocatable :: dlap(:,:)
      
      complex *16, allocatable :: ahmat1(:,:), ahmat2(:,:), ahmat3(:,:)
      complex *16, allocatable :: hmat(:,:), resnew2(:,:)
      complex *16, allocatable :: res1(:,:), res2(:,:), res3(:,:)


      integer, allocatable :: ind2p(:,:), ind2pout(:,:),ip2indout(:,:,:)
      
      real *8, allocatable :: sums1(:), sums2(:), sums3(:)
      
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /

      character type
      
      call prini(6,13)


      ndeg = 16
      nup = 12

      ndegout = ndeg + nup*2
      
      type = 't'

      rr = 2.5d0
      zk = rr*(1/sqrt(2.0d0) + im/sqrt(2.0d0))
c      zk = rr
      
      call legetens_npol_3d(ndeg,type,npol)
      call legetens_npol_3d(ndegout,type,npolout)

      call prinf('npol *',npol,1)
      call prinf('npolout *',npolout,1)

      allocate(ind2p(3,npol),ind2pout(3,npolout),
     1     ip2indout(0:ndegout,0:ndegout,0:ndegout))

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

      allocate(res1(npolout,npol),res2(npolout,npol))

      call cpu_time(t1)
c$    t1 = omp_get_wtime()      
      do j = 1,npol
         nd = 2
         call legetens_lape_3d(nd,ndegout,type,ahmat1(1,j),res1(1,j))
         call legetens_lape_3d(nd,ndegout,type,ahmat2(1,j),res2(1,j))
         
         do i = 1,npolout
            res1(i,j) = res1(i,j) + ahmat1(i,j)*zk**2
            res2(i,j) = res2(i,j) + ahmat2(i,j)*zk**2
         enddo
      enddo
      call cpu_time(t2)
c$    t2 = omp_get_wtime()            
      
      write(*,*) 'time to test residuals', t2-t1

      
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
               sums1(i) = sums1(i) + abs(one-res1(j,i))
               sums2(i) = sums2(i) + abs(one-res2(j,i))
            else
               sums1(i) = sums1(i) + abs(res1(j,i))
               sums2(i) = sums2(i) + abs(res2(j,i))               
            endif
         enddo

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
