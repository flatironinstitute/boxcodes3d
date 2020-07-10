      implicit real *8 (a-h,o-z)
      complex *16 zk, im, zero, one

      real *8, allocatable :: allerrsup(:,:,:), allerrsdown(:,:,:)
      
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /

      character type

      complex *16 zks(20)
      integer :: nups(20)
      integer, allocatable :: ind2p(:,:)
      
      call prini(6,13)

      nnup = 17
      do i = 1,16
         nups(i) = i
      enddo
      nups(17) = 40

      nzk = 10
      do i = 1,nzk
         zks(i) = (1.0d0+0.5d0*(i-1))*(1/sqrt(2.0d0) + im/sqrt(2.0d0))
         write(*,*) zks(i), abs(zks(i))
      enddo

      type = 't'
      ndeg = 16
      call legetens_npol_3d(ndeg,type,npol)

      call prinf('npol *',npol,1)

      allocate(ind2p(3,npol))

      call legetens_ind2pow_3d(ndeg,type,ind2p)
      
      open(unit=101,file='h3danti_compare_dims.txt')
      open(unit=102,file='h3danti_compare_zk.txt')
      open(unit=103,file='h3danti_compare_nup.txt')
      open(unit=104,file='h3danti_compare_allup.txt')
      open(unit=105,file='h3danti_compare_alldown.txt')
      open(unit=106,file='h3danti_compare_ind2p.txt')

      write(101,*) 'anti-Helmholtzian comparison test.... '
      write(101,*) 'main files are ...allup and ...alldown'
      write(101,*) '      3D arrays of size npol,nzk,nnup'
      write(101,*) 'nnup - number of nups used', nnup
      write(101,*) 'nzk - number of zks used', nzk
      write(101,*) 'npol - number of polynomials',npol
      write(101,*) 'ndeg - degree of polynomials',ndeg

      close(101)

      do i = 1,nnup
         write(103,*) nups(i)
      enddo

      close(103)

      do i = 1,nzk
         write(102,*) real(zks(i)), -real(im*zks(i))
      enddo

      close(102)

      do i = 1,npol
         do j = 1,3
            write(106,*) ind2p(j,i)
         enddo
      enddo

      close(106)

      allocate(allerrsup(npol,nzk,nnup),allerrsdown(npol,nzk,nnup))

      do ii = 1,nnup
         nup = nups(ii)
         do jj = 1,nzk
            zk = zks(jj)

            write(*,*) '-----------------------------------------------'
            write(*,*) 'Running nup ', nup
            write(*,*) '        which is ', ii, 'out of ', nnup
            write(*,*) 'Running zk ', zk
            write(*,*) '        which is ', jj, 'out of ', nzk
            call runit(nup,zk,ndeg,allerrsup(1,jj,ii),
     1           allerrsdown(1,jj,ii))

            call prin2('allerrsup *',allerrsup(1,jj,ii),20)
            call prin2('allerrsdown *',allerrsdown(1,jj,ii),20)            

            do kk = 1,npol
               write(104,*) allerrsup(kk,jj,ii)
               write(105,*) allerrsdown(kk,jj,ii)
            enddo
            
         enddo
      enddo
      
           
      close(104)
      close(105)
      
      stop
      end



      subroutine runit(nup,zk,ndeg,errsup,errsdown)
      implicit real *8 (a-h,o-z)
      real *8 errsup(*), errsdown(*)
      complex *16 zk, im, zero, one

      complex *16, allocatable :: ahmatup(:,:), ahmatdown(:,:)
      complex *16, allocatable :: resup(:,:), resdown(:,:)


      integer, allocatable :: ind2p(:,:), ind2pout(:,:),ip2indout(:,:,:)

      data zero /(0.0d0,0.0d0)/
      data one /(1.0d0,0.0d0)/
      data im /(0.0d0,1.0d0)/            

      character type
      
      ndegout = ndeg + nup*2

      call prin2('zk *',zk,2)
      call prin2('abs(zk) *',abs(zk),1)
      call prinf('ndeg *',ndeg,1)
      call prinf('nup *',nup,1)      
      call prinf('ndegout ',ndegout,1)
      
      type = 't'

      call legetens_npol_3d(ndeg,type,npol)
      call legetens_npol_3d(ndegout,type,npolout)

      call prinf('npol *',npol,1)
      call prinf('npolout *',npolout,1)

      allocate(ind2p(3,npol),ind2pout(3,npolout),
     1     ip2indout(0:ndegout,0:ndegout,0:ndegout))

      call legetens_ind2pow_3d(ndeg,type,ind2p)
      call legetens_ind2pow_3d(ndegout,type,ind2pout)      
      
      allocate(ahmatdown(npol,npol),ahmatup(npolout,npol))

      call cpu_time(t1)
      call h3danti_legeupfast(ndeg,nup,type,zk,ahmatup,npolout)
      call cpu_time(t2)

      write(*,*) 'time upward recurrence', t2-t1

      call cpu_time(t1)
      call h3danti_legedownfast(ndeg,type,zk,ahmatdown,npol)
      call cpu_time(t2)

      write(*,*) 'time downward recurrence', t2-t1

      allocate(resup(npolout,npol),resdown(npol,npol))

      call cpu_time(t1)
      do j = 1,npol
         nd = 2
         call legetens_lape_3d(nd,ndegout,type,ahmatup(1,j),resup(1,j))
         call legetens_lape_3d(nd,ndeg,type,ahmatdown(1,j),resdown(1,j))
         
         do i = 1,npolout
            resup(i,j) = resup(i,j) + ahmatup(i,j)*zk**2
         enddo

         do i = 1,npol
            resdown(i,j) = resdown(i,j) + ahmatdown(i,j)*zk**2
         enddo
      enddo
      call cpu_time(t2)
      
      write(*,*) 'time to test residuals', t2-t1

      
      do i = 1,npol
         errsup(i) = 0
         errsdown(i) = 0
         iin = ind2p(1,i)
         jin = ind2p(2,i)
         kin = ind2p(3,i)
         do j = 1,npolout
            iout = ind2pout(1,j)
            jout = ind2pout(2,j)
            kout = ind2pout(3,j)
            if (iin .eq. iout .and. jin .eq. jout
     1           .and. kin .eq. kout) then
               errsup(i) = errsup(i) + abs(one-resup(j,i))
            else
               errsup(i) = errsup(i) + abs(resup(j,i))
            endif
         enddo
         do j = 1,npol
            if (i .eq. j) then
               errsdown(i) = errsdown(i) + abs(one-resdown(j,i))
            else
               errsdown(i) = errsdown(i) + abs(resdown(j,i))
            endif
         enddo
      enddo

      smaxup = 0
      smaxdown = 0
      do i = 1,npol
         smaxup = max(smaxup,errsup(i))
         smaxdown = max(smaxdown,errsdown(i))
      enddo


      call prin2('max err, legeup *',smaxup,1)
      call prin2('max err, legedown *',smaxdown,1)      
      
      return
      end
