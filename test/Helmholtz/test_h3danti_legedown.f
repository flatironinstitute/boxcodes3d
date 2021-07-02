      implicit real *8 (a-h,o-z)
      complex *16 zk, im, zero, one, alpha, beta

      real *8, allocatable :: dlap(:,:)
      
      complex *16, allocatable :: ahmat2(:,:), ahmat3(:,:)
      complex *16, allocatable :: hmat(:,:)
      complex *16, allocatable :: res2(:,:), res3(:,:)


      integer, allocatable :: ind2p(:,:)
      
      real *8, allocatable :: sums2(:), sums3(:)
      
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /

      character type
      
      call prini(6,13)


      ndeg = 8

      tol = 1d-15
      toltest = 1.0d-11
      ntests = 2
      i1 = 0
      i2 = 0

      type = 't'

      zk = 1.5d0 + im*0.0d0
      
      call legetens_npol_3d(ndeg,type,npol)

      call prinf('npol *',npol,1)

      allocate(ahmat2(npol,npol),ahmat3(npol,npol))
      call cpu_time(t1)
      call h3danti_legedown(ndeg,type,zk,ahmat2,npol)
      call cpu_time(t2)

      write(*,*) 'time downward recurrence ', t2-t1

      call cpu_time(t1)
      call h3danti_legedownfast(ndeg,type,zk,ahmat3,npol)
      call cpu_time(t2)

      write(*,*) 'time downward recurrence (fast?)', t2-t1


      allocate(dlap(npol,npol),hmat(npol,npol),
     1     res2(npol,npol),res3(npol,npol))
      call legetens_lapmat_3d(ndeg,ndeg,type,dlap,npol)      

      do i= 1,npol
         do j = 1,npol
            hmat(j,i) = dlap(j,i)
            if (j .eq. i) hmat(j,i) = hmat(j,i) + zk**2
            
         enddo
      enddo
      
      alpha = one
      beta = zero
      call zgemm('N','N',npol,npol,npol,alpha,hmat,npol,
     1     ahmat2,npol,beta,res2,npol)
      call zgemm('N','N',npol,npol,npol,alpha,hmat,npol,
     1     ahmat3,npol,beta,res3,npol)
      

      
      sum = 0
      sum1 = 0
      sum2 = 0
      sum3 = 0
      
      sumrel = 0

      allocate(sums2(npol),sums3(npol))
      
      do i = 1,npol
         sums2(i) = 0
         sums3(i) = 0
         do j = 1,npol
            if (j .eq. i) then
               sum2 = sum2 + abs(one-res2(j,i))**2
               sums2(i) = sums2(i) + abs(one-res2(j,i))**2
               
               sum3 = sum3 + abs(one-res3(j,i))**2
               sums3(i) = sums3(i) + abs(one-res3(j,i))**2
               
            else
               sum2 = sum2 + abs(res2(j,i))**2
               sums2(i) = sums2(i) + abs(res2(j,i))**2
               
               sum3 = sum3 + abs(res3(j,i))**2
               sums3(i) = sums3(i) + abs(res3(j,i))**2
            endif

         enddo

         sums2(i) = sqrt(sums2(i))
         sums3(i) = sqrt(sums3(i))
      enddo

      call prin2('zk *',zk,2)

      call prin2('residual recurrence *',sqrt(sum2/npol),1)
      call prin2('residual recurrence (fast?) *',sqrt(sum3/npol),1)

      if(sqrt(sum2/npol).lt.toltest) i1 = 1
      if(sqrt(sum3/npol).lt.toltest) i2 = 1

      nsuccess = i1 + i2

      open(unit=33,file='../../print_testres.txt',access='append')
      write(33,'(a,i1,a,i1,a)') 'Successfully completed ',nsuccess,
     1  ' out of ',ntests,' in h3d anti legedown testing suite'
      write(*,'(a,i1,a,i1,a)') 'Successfully completed ',nsuccess,
     1  ' out of ',ntests,' in h3d anti legedown testing suite'
      close(33)
      

      
      stop
      end
