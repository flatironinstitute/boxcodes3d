      implicit real *8 (a-h,o-z)
      complex *16 zk, im, zero, one, alpha, beta

      complex *16, allocatable :: ahmat(:,:)
      complex *16, allocatable :: res(:,:)

      integer, allocatable :: ind2p(:,:), ind2pout(:,:)
      integer, allocatable :: ip2ind(:,:,:), ip2indout(:,:,:)      
      
      real *8, allocatable :: sums(:), errsup(:), errsdown(:)
      
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /

      character type
      
      call prini(6,13)


      ndeg = 7
      nup = 8

      ndegout = ndeg + nup*2
       
      type = 't'

      r = 2.421d0
c      r = 1.1d0
c      r = 9d0
c      r = 0.00001d0

      zk = r/sqrt(2.0d0) + im*r/sqrt(2.0d0)
      zk = r
c      zk = r*im
      
      
      call legetens_npol_3d(ndeg,type,npol)
      call legetens_npol_3d(ndegout,type,npolout)

      call prinf('npol *',npol,1)
      call prinf('npolout *',npolout,1)

      allocate(ind2p(3,npol),ind2pout(3,npolout))
      allocate(ip2ind(ndeg+1,ndeg+1,ndeg+1),
     1     ip2indout(ndegout+1,ndegout+1,ndegout+1))

      call legetens_ind2pow_3d(ndeg,type,ind2p)
      call legetens_ind2pow_3d(ndegout,type,ind2pout)      
      call legetens_pow2ind_3d(ndeg,type,ip2ind)
      call legetens_pow2ind_3d(ndegout,type,ip2indout)      

      allocate(ahmat(npolout,npol),errsup(npol),errsdown(npol))

      t1 = timewrap()
      call h3danti_form(ndeg,nup,type,zk,ahmat,npolout,derrmax,
     1     errsup,errsdown)
      t2 = timewrap()

      write(*,*) 'time to form anti-helmholtzian ', t2-t1
      write(*,*) 'computed derrmax ', derrmax

c     
      allocate(res(npolout,npol))


      t1 = timewrap()
c$omp parallel do private(i,ndc,j)
      do i = 1,npol
         ndc = 2
         call legetens_lape_3d(ndc,ndegout,type,ahmat(1,i),
     1        res(1,i))
         
         do j = 1,npolout
            res(j,i) = res(j,i)+zk**2*ahmat(j,i)
         enddo
      enddo
c$omp end parallel do      
      
      t2 = timewrap()

      write(*,*) 'time to test anti-helmholtzian ', t2-t1
      
      
      allocate(sums(npol))
      
      do i = 1,npol
         sums(i) = 0
         iin = ind2p(1,i)
         jin = ind2p(2,i)
         kin = ind2p(3,i)
         do j = 1,npolout
            iout = ind2pout(1,j)
            jout = ind2pout(2,j)
            kout = ind2pout(3,j)
            if (iin .eq. iout .and. jin .eq. jout
     1           .and. kin .eq. kout) then
               sums(i) = sums(i) + abs(one-res(j,i))**2
            else
               sums(i) = sums(i) + abs(res(j,i))**2
            endif
         enddo

         sums(i) = sqrt(sums(i))
      enddo

      call prin2('zk *',zk,2)

      summax = 0.0d0
      do i = 1,npol
         summax = max(summax,sums(i))
      enddo

      call prin2('summax *',summax,1)
      call prin2('summax/sqrt(npol) *',summax/sqrt(npol*1d0),1)      
      call prinf('nup *',nup,1)
      
      stop
      end

      real *8 function timewrap()
      implicit none
c$    real *8 omp_get_wtime      
      timewrap = 0d0
      call cpu_time(timewrap)
c$    timewrap = omp_get_wtime()
      return
      end
