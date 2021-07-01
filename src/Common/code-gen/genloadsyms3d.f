cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     determine symmetries relating all of the table types to
c     the reference tables (for colleague, stob, and btos)
c
c     this is done by brute force
c
c     subroutines are printed to the file loadsyms3d.f
      
      implicit real *8 (a-h,o-z)

      iw = 20
      nq = 3

      open(unit = iw, file='loadsyms3d.f')
      call prini(6,13)

      call testit(iw,nq)

      stop
      end




      subroutine testit(iw,nq)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      implicit real *8 (a-h,o-z)
      dimension wc(3,nq**3,27), wbtos(3,nq**3,56), wstob(3,nq**3,56)
      dimension wcref(3,nq**3,4), wbtosref(3,nq**3,3)
      dimension wstobref(3,nq**3,3)
      dimension xt(4*nq**3), yt(4*nq**3), zt(4*nq**3)
      dimension xts(nq**3), yts(nq**3), zts(nq**3)
      dimension xq(nq), xyzc(3)
      dimension dlens(nq**3), dlensref(nq**3)
      dimension alldist(nq**3,nq**3), alldistref(nq**3,nq**3)
      dimension ipperm(nq**3), icperm(nq**3), icsign(nq**3)
      dimension ippermold(nq**3)
      dimension ibtosdimperm(56), ibtosflipopt(56), ibtosref(56)
      dimension istobdimperm(56), istobflipopt(56), istobref(56)
      dimension icdimperm(56), icflipopt(56), icref(56)
      dimension idimperms(3,6)
      data idimperms / 1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1 /
      data ndimperms / 6 /
      dimension iflipopts(3,8)      
      data iflipopts / 1,1,1, -1,1,1, 1,-1,1, -1,-1,1, 1,1,-1, -1,1,-1,
     1     1,-1,-1, -1,-1,-1 /
      
      data nflipopts / 8 /

      character *32 title
      character cp
      logical iwork

      cp = 't'
      

      nq3 = nq**3

c     tolerance for a match

      tol = 1.0d-14

c     define box size, box corner for self

      bs = 1.0d0
      xyzc(1) = -0.5d0
      xyzc(2) = -0.5d0
      xyzc(3) = -0.5d0

c     create equispaced grid
      
      h = bs/nq
      hh = h/2
      do i = 1,nq
         xq(i) = -hh + i*h
      enddo

c     grab reference points
      
      call tensrefpts3d(xq,nq,bs,xyzc,wcref,wbtosref,wstobref)

c     grab all targets 

      call alltargs3d(xq,nq,bs,xyzc,wc,wbtos,wstob)
      
c     find permutations which match distances for btos
      
      do i = 1,56
         iworkever = 0
         call dists3d(wbtos(1,1,i),nq3,wcref(1,1,1),nq3,alldist)
         do j = 1,3
            call dists3d(wbtosref(1,1,j),nq3,wcref(1,1,1),nq3,
     1           alldistref)
            do i1 = 1,ndimperms
               do i2 = 1,nflipopts
c     for given permutation of dimensions and flipping of
c     index orders, check if you've got it
                  call buildperm3d(idimperms(1,i1),
     1                 iflipopts(1,i2),nq,cp,ipperm,icperm,icsign)
c                  call prinf('ipperm *',ipperm,nq**3)
                  iwork = .true.
                  do jj = 1,nq3
                     do ii = 1,nq3
                        iwork = iwork .and.
     1                       (abs(alldist(ipperm(ii),ipperm(jj))
     1                       -alldistref(ii,jj)) .lt. tol)
                     enddo
                  enddo

c     save if new
                  if (iwork .and. iworkever .eq. 0) then
                     ibtosref(i) = j
                     ibtosdimperm(i) = i1
                     ibtosflipopt(i) = i2
                     do iii = 1,nq3
                        ippermold(iii) = ipperm(iii)
                     enddo
                  endif
                  
                  if (iwork) iworkever = iworkever + 1
               enddo
            enddo
         enddo
         if (iworkever .eq. 0) then
            write(*,*) 'fail', i
         else
            write(*,*) 'success', i, ibtosref(i), iworkever
         endif
         
      enddo
      
c     find permutations which match distances for stob
      
      do i = 1,56
         iworkever = 0
         call dists3d(wstob(1,1,i),nq3,wcref(1,1,1),nq3,alldist)
         do j = 1,3
            call dists3d(wstobref(1,1,j),nq3,wcref(1,1,1),nq3,
     1           alldistref)
            do i1 = 1,ndimperms
               do i2 = 1,nflipopts
c     for given permutation of dimensions and flipping of
c     index orders, check if you've got it
                  call buildperm3d(idimperms(1,i1),
     1                 iflipopts(1,i2),nq,cp,ipperm,icperm,icsign)
c                  call prinf('ipperm *',ipperm,nq**3)
                  iwork = .true.
                  do jj = 1,nq3
                     do ii = 1,nq3
                        iwork = iwork .and.
     1                       (abs(alldist(ipperm(ii),ipperm(jj))
     1                       -alldistref(ii,jj)) .lt. tol)
                     enddo
                  enddo

c     save if new
                  if (iwork .and. iworkever .eq. 0) then
                     istobref(i) = j
                     istobdimperm(i) = i1
                     istobflipopt(i) = i2
                     do iii = 1,nq3
                        ippermold(iii) = ipperm(iii)
                     enddo
                  endif
                  
                  if (iwork) iworkever = iworkever + 1
               enddo
            enddo
         enddo
         if (iworkever .eq. 0) then
            write(*,*) 'fail', i
         else
            write(*,*) 'success', i, istobref(i), iworkever
         endif
         
      enddo
      
c     find permutations which match distances for colleague
      
      do i = 1,27
         iworkever = 0
         call dists3d(wc(1,1,i),nq3,wcref(1,1,1),nq3,alldist)
         do j = 1,4
            call dists3d(wcref(1,1,j),nq3,wcref(1,1,1),nq3,
     1           alldistref)
            do i1 = 1,ndimperms
               do i2 = 1,nflipopts
c     for given permutation of dimensions and flipping of
c     index orders, check if you've got it
                  call buildperm3d(idimperms(1,i1),
     1                 iflipopts(1,i2),nq,cp,ipperm,icperm,icsign)
c                  call prinf('ipperm *',ipperm,nq**3)
                  iwork = .true.
                  do jj = 1,nq3
                     do ii = 1,nq3
                        iwork = iwork .and.
     1                       (abs(alldist(ipperm(ii),ipperm(jj))
     1                       -alldistref(ii,jj)) .lt. tol)
                     enddo
                  enddo

c     save if new
                  if (iwork .and. iworkever .eq. 0) then
                     icref(i) = j
                     icdimperm(i) = i1
                     icflipopt(i) = i2
                     do iii = 1,nq3
                        ippermold(iii) = ipperm(iii)
                     enddo
                  endif
                  
                  if (iwork) iworkever = iworkever + 1
               enddo
            enddo
         enddo
         if (iworkever .eq. 0) then
            write(*,*) 'fail', i
         else
            write(*,*) 'success', i, icref(i), iworkever
         endif
         
      enddo

c     print to file

      write(iw,'(a)') 'c'
      write(iw,'(a)') 'c'            
      write(iw,'(a)') 'c       this file was generated automatically'
      write(iw,'(a)') 'c       it contains subroutines which load'
      write(iw,'(a)') 'c       symmetry info for volume code tables'
      write(iw,'(a)') 'c'      
      write(iw,'(a)') 'c'      
      write(iw,'(a)') ''      
      write(iw,'(a)') ''      
      write(iw,'(a)') ''      
      write(iw,'(a)') ''      
      write(iw,*) '      subroutine loadsymsc(iref,idimp,iflip)'
      write(iw,*) '      implicit real *8 (a-h,o-z)'
      write(iw,*) ''
      write(iw,*) '      dimension iref(*), idimp(3,*), iflip(3,*)'
      write(iw,*) ''
      do i = 1,27
         i1 = icdimperm(i)
         i2 = icflipopt(i)         
         write(iw,'(a,I2,a,I2)') '      iref(', i, ')  = ', icref(i)
         write(iw,'(a,I2,a,I2)') '      idimp(1,', i, ') = ',
     1        idimperms(1,i1)
         write(iw,'(a,I2,a,I2)') '      idimp(2,', i, ') = ',
     1        idimperms(2,i1)
         write(iw,'(a,I2,a,I2)') '      idimp(3,', i, ') = ',
     1        idimperms(3,i1)
         write(iw,'(a,I2,a,I2)') '      iflip(1,', i, ') = ',
     1        iflipopts(1,i2)
         write(iw,'(a,I2,a,I2)') '      iflip(2,', i, ') = ',
     1        iflipopts(2,i2)
         write(iw,'(a,I2,a,I2)') '      iflip(3,', i, ') = ',
     1        iflipopts(3,i2)
      enddo
      write(iw,*) ''
      write(iw,*) '      return'
      write(iw,*) '      end'
      write(iw,*) ''      
      write(iw,*) ''
      
      write(iw,*) '      subroutine loadsymsbtos(iref,idimp,iflip)'
      write(iw,*) '      implicit real *8 (a-h,o-z)'
      write(iw,*) ''
      write(iw,*) '      dimension iref(*), idimp(3,*), iflip(3,*)'
      write(iw,*) ''
      do i = 1,56
         i1 = ibtosdimperm(i)
         i2 = ibtosflipopt(i)         
         write(iw,'(a,I2,a,I2)') '      iref(', i, ')  = ', ibtosref(i)
         write(iw,'(a,I2,a,I2)') '      idimp(1,', i, ') = ',
     1        idimperms(1,i1)
         write(iw,'(a,I2,a,I2)') '      idimp(2,', i, ') = ',
     1        idimperms(2,i1)
         write(iw,'(a,I2,a,I2)') '      idimp(3,', i, ') = ',
     1        idimperms(3,i1)
         write(iw,'(a,I2,a,I2)') '      iflip(1,', i, ') = ',
     1        iflipopts(1,i2)
         write(iw,'(a,I2,a,I2)') '      iflip(2,', i, ') = ',
     1        iflipopts(2,i2)
         write(iw,'(a,I2,a,I2)') '      iflip(3,', i, ') = ',
     1        iflipopts(3,i2)
      enddo
      write(iw,*) ''
      write(iw,*) '      return'
      write(iw,*) '      end'
      write(iw,*) ''      
      write(iw,*) ''
      
      write(iw,*) '      subroutine loadsymsstob(iref,idimp,iflip)'
      write(iw,*) '      implicit real *8 (a-h,o-z)'
      write(iw,*) ''
      write(iw,*) '      dimension iref(*), idimp(3,*), iflip(3,*)'
      write(iw,*) ''
      do i = 1,56
         i1 = istobdimperm(i)
         i2 = istobflipopt(i)         
         write(iw,'(a,I2,a,I2)') '      iref(', i, ')  = ', istobref(i)
         write(iw,'(a,I2,a,I2)') '      idimp(1,', i, ') = ',
     1        idimperms(1,i1)
         write(iw,'(a,I2,a,I2)') '      idimp(2,', i, ') = ',
     1        idimperms(2,i1)
         write(iw,'(a,I2,a,I2)') '      idimp(3,', i, ') = ',
     1        idimperms(3,i1)
         write(iw,'(a,I2,a,I2)') '      iflip(1,', i, ') = ',
     1        iflipopts(1,i2)
         write(iw,'(a,I2,a,I2)') '      iflip(2,', i, ') = ',
     1        iflipopts(2,i2)
         write(iw,'(a,I2,a,I2)') '      iflip(3,', i, ') = ',
     1        iflipopts(3,i2)
      enddo
      write(iw,*) ''
      write(iw,*) '      return'
      write(iw,*) '      end'
      write(iw,*) ''      
      write(iw,*) ''
      
      return
      end

