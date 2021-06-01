      subroutine fakepolya3d(ng,ifloor,npt,pts,ipt2depth,idepth2pt)
c
c     This routine builds a thinned out tensor product grid of
c     Legendre nodes of order ng while maintaining symmetries
c     with respect to flipping coordinate directions and
c     permuting the order of x,y,z coordinates. The nodes are
c     removed uniformly (in terms of indices) at random.
c
c     input:
c
c     ng - integer, order of starting legendre tensor grid
c     ifloor - integer, return a set of nodes with at least
c        ifloor nodes
c
c     output:
c      
c     npt - integer, number of points in thinned out rule
c     pts - real*8 (3,npts), point coordinates
c     ipt2depth - integer (3,npts). ipt2depth(j,i) is the
c         index of point i in the jth coordinate direction
c         in the original tensor grid 
c     idepth2pt - integer (ng,ng,ng). let (i,j,k) be the
c         indices for a point in the original tensor array
c         if this point was kept, idepth2pt(i,j,k) is the
c         index in the output pts array corresponding to
c         (i,j,k) in the original array. If (i,j,k) was 
c         removed, then idepth2pt(i,j,k) = -1
c
      implicit real *8 (a-h,o-z)
      integer npt, ipt2depth(3,*),idepth2pt(ng,ng,ng)
      real *8 pts(3,*)

      real *8 u,v,w
      real *8, allocatable :: t(:)
      integer nperm, iperms(3,6), itemp(3)
      parameter (nperm = 6)

      data iperms / 1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1 /

      iseed = 1234

      dd = hkrand(iseed)

      itype = 0
      allocate(t(ng))
      call legeexps(itype,ng,t,u,v,w)

      npt = ng**3

c     for now use idepth2pt to store in/out

c     start all in

      do iz = 1,ng
         do iy = 1,ng
            do ix = 1,ng
               idepth2pt(ix,iy,iz) = 1
            enddo
         enddo
      enddo
      
      ntry2 = 1000
      ntry = 100

      k = ng/2

      write(*,*) 'k ng ', k, ng
      
c     remove x,y,z coords at random while you can
      do iii = 1,ntry2
         
         do i = 1,ntry
            if (mod(ng,2) .eq. 0) then
               ix = 1+hkrand(0)*k
               iy = 1+hkrand(0)*k
               iz = 1+hkrand(0)*k

               if (npt .lt. ifloor + 48 .and. npt .ge. ifloor + 24) then
                  iy = ix
               endif
               if (npt .lt. ifloor + 24 .and. npt .ge. ifloor + 8) then
                  iy = ix
                  iz = ix
               endif

               if (npt .lt. ifloor + 8) exit

               if ( (ix .lt. 1) .or. (iy.lt. 1) .or. (iz .lt. 1) .or.
     1              (ix .gt. k) .or. (iy.gt. k) .or. (iz .gt. k)) then
                  write(*,*) 'uh oh '
                  cycle
               endif
               
            else
               ix = hkrand(0)*(k+1)
               iy = hkrand(0)*(k+1)
               iz = hkrand(0)*(k+1)

               if (npt .lt. ifloor + 48 .and. npt .ge. ifloor + 24) then
                  iy = ix
               endif
               if (npt .lt. ifloor + 24 .and. npt .ge. ifloor + 12) then
                  iy = ix
                  iz = 0
               endif
               if (npt .lt. ifloor + 12 .and. npt .ge. ifloor + 8) then
                  iy = ix
                  iz = ix
               endif
               if (npt .lt. ifloor + 8 .and. npt .ge. ifloor + 6) then
                  iy = 0
                  ix = 0
               endif
               if (npt .lt. ifloor + 6 .and. npt .ge. ifloor + 1) then
                  ix = 0
                  iy = 0
                  iz = 0
               endif

               if (npt .lt. ifloor + 6 .and.
     1              idepth2pt(k+1,k+1,k+1) .eq. -1) exit
               
               
               if ( (ix .lt. 0) .or. (iy.lt. 0) .or. (iz .lt. 0) .or.
     1              (ix .gt. k) .or. (iy.gt. k) .or. (iz .gt. k)) then
                  cycle
               endif
            endif
            

            if (mod(ng,2) .eq. 0) then            
               ix2 = ix+k
               iy2 = iy+k
               iz2 = iz+k
            else
               ix2 = ix+k+1
               iy2 = iy+k+1
               iz2 = iz+k+1
            endif
            
            if (idepth2pt(ix2,iy2,iz2) .eq. -1) then
c     already gone
c               write(*,*) 'already gone'
               cycle
            endif

c     figure out how many are removed
            
            nrem = 48
            if (ix .eq. 0) nrem = nrem/2
            if (iy .eq. 0) nrem = nrem/2
            if (iz .eq. 0) nrem = nrem/2
            if ( (ix .eq. iy) .and. (ix .eq. iz)) then
c               write(*,*) 'all same'
               nrem = nrem/6
            else
               if ( (ix .eq. iy) .or. (ix .eq. iz) .or.
     1              (iy .eq. iz)) then
c                  write(*,*) 'two same'                  
                  nrem = nrem/2
               endif
            endif
c     write(*,*) ix,iy,iz,nrem

            if (npt - nrem .lt. ifloor) write(*,*) 'uho h', npt, nrem
            npt = npt-nrem

c     remove them

            do ip = 1,nperm
               itemp(iperms(1,ip)) = ix
               itemp(iperms(2,ip)) = iy
               itemp(iperms(3,ip)) = iz
               do ixs = 1,2
                  do iys = 1,2
                     do izs = 1,2
                        if (mod(ng,2) .eq. 0) then
                           ix2 = itemp(1) + k
                           iy2 = itemp(2) + k
                           iz2 = itemp(3) + k
                           if (ixs .eq. 2) ix2 = k-itemp(1)+1
                           if (iys .eq. 2) iy2 = k-itemp(2)+1
                           if (izs .eq. 2) iz2 = k-itemp(3)+1
                        else
                           ix2 = itemp(1) + k + 1
                           iy2 = itemp(2) + k + 1
                           iz2 = itemp(3) + k + 1
                           if (ixs .eq. 2) ix2 = k + 1 -itemp(1)
                           if (iys .eq. 2) iy2 = k + 1 -itemp(2)
                           if (izs .eq. 2) iz2 = k + 1 -itemp(3)
                        endif
                        idepth2pt(ix2,iy2,iz2) = -1
                     enddo
                  enddo
               enddo
            enddo

            exit
            
         enddo

         if (npt .lt. ifloor) exit
      enddo



c     get point ordering and such

      npt = 0
      do iz = 1,ng
         do iy = 1,ng
            do ix = 1,ng
               if (idepth2pt(ix,iy,iz) .eq. 1) then
                  npt = npt+1
                  ipt2depth(1,npt) = ix
                  ipt2depth(2,npt) = iy
                  ipt2depth(3,npt) = iz
                  pts(1,npt) = t(ix)
                  pts(2,npt) = t(iy)
                  pts(3,npt) = t(iz)
                  idepth2pt(ix,iy,iz) = npt
               endif
            enddo
         enddo
      enddo

               
      
      return
      end
      


      subroutine ichecksymms(idepth2pt,ipt2depth,ng,npt,ipass1,ipass2,
     1     ipass3)
c
c     correctness checking utility for a symmetry preserving thinning
c     of a tensor product grid
c
c     input:
c     
c     idepth2pt - ng x ng x ng array where idepth2pt(i,j,k) = -1 indicates
c                 that (i,j,k) has been removed and idepth2pt(i,j,k) is
c                 a unique label between 1 and npt otherwise.
c     ng - dimension of idepth2pt
c     npt - number of points that should be included in idepth2pt
c
c     output:
c
c     ipass1 - 1 if symmetries are satisfied for every included point,
c              -(i+ng*(j-1)+ng**2*(k-1)) if (i,j,k) failed symmetry check
c     ipass2 - 1 if all points accounted for, 0 if double label, -ipt if
c              point ipt is missing
c     ipass3 - 1 if ipt2depth is consistent with idepth2pt, -ipt if
c              point ipt is inconsistent
c
      implicit real *8 (a-h,o-z)
      integer idepth2pt(ng,ng,ng), ipt2depth(3,npt)

      parameter (nperm=6)
      integer itemp(3), iperms(3,nperm)
      integer, allocatable :: ihit(:)
      data iperms / 1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1 /

      allocate(ihit(npt))

      do i = 1,npt
         ihit(i) = 0
      enddo

      ipass1 = 1
      
      do iz = 1,ng
         do iy = 1,ng
            do ix = 1,ng
               ipt = idepth2pt(ix,iy,iz)
               if (ipt .le. 0) cycle

               ihit(ipt) = ihit(ipt) + 1
               
               do ip = 1,nperm
                  itemp(iperms(1,ip)) = ix
                  itemp(iperms(2,ip)) = iy
                  itemp(iperms(3,ip)) = iz
                  do ixs = 1,2
                     do iys = 1,2
                        do izs = 1,2
                           ix2 = itemp(1)
                           iy2 = itemp(2)
                           iz2 = itemp(3)
                           if (ixs .eq. 2) ix2 = ng-itemp(1)+1
                           if (iys .eq. 2) iy2 = ng-itemp(2)+1
                           if (izs .eq. 2) iz2 = ng-itemp(3)+1
                           ipt2 = idepth2pt(ix2,iy2,iz2)

                           if (ipt2 .le. 0 .and. ipass1 .eq. 1) then
                              ipass1 = -(ix+ng*(iy-1)+ng**2*(iz-1))
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      ipass3 = 1
      do i = 1,npt
         ix = ipt2depth(1,i)
         iy = ipt2depth(2,i)
         iz = ipt2depth(3,i)
         if (ix .lt. 0 .or. iy .lt. 0 .or. iz .lt. 0) then
            write(*,*) 'here fail'
            exit
         endif
         if (ix .gt. ng .or. iy .gt. ng .or. iz .gt. ng) then
            write(*,*) 'here fail 2', ix,iy,iz,i
            exit
         endif
         if (i .ne. idepth2pt(ix,iy,iz) .and. ipass3 .eq. 1) ipass3 = -i
      enddo
         
         

      ipass2 = 1
      do i = 1,npt
         if (ihit(i) .eq. 0 .and. ipass2 .eq. 1) ipass2=0
         if (ihit(i) .gt. 1 .and. ipass2 .eq. 1) ipass2=-i
      enddo
      
      return
      end


      subroutine fakepolya3dredblack(k,npt,pts,ipt2depth,idepth2pt)
c
c     This routine builds a thinned out tensor product grid of
c     Legendre nodes of order 2*k+1 while maintaining symmetries
c     with respect to flipping coordinate directions and
c     permuting the order of x,y,z coordinates.
c
c     Let (i,j,k) be the indices of a standard tensor grid
c     point. This point is kept if i+j+k is odd.
c      
c      
      implicit real *8 (a-h,o-z)
      integer k, npt, ipt2depth(3,*),idepth2pt(2*k+1,2*k+1,2*k+1)
      real *8 pts(3,*)

      real *8 u,v,w
      real *8, allocatable :: t(:)

      itype = 0
      n = 2*k+1
      allocate(t(n))
      call legeexps(itype,n,t,u,v,w)

      npt = 0
      do iz = 1,2*k+1
         do iy = 1,2*k+1
            do ix = 1,2*k+1
               idepth2pt(ix,iy,iz) = -1
               if (mod(ix+iy+iz,2) .eq. 1) then
                  npt = npt + 1
                  pts(1,npt) = t(ix)
                  pts(2,npt) = t(iy)
                  pts(3,npt) = t(iz)
                  ipt2depth(1,npt) = ix
                  ipt2depth(2,npt) = iy
                  ipt2depth(3,npt) = iz
                  idepth2pt(ix,iy,iz) = npt
               endif
            enddo
         enddo
      enddo
      
      return
      end
      
      subroutine fakepolya3dv2(k,npt,pts,ipt2depth,idepth2pt)
c
c     This routine builds a thinned out tensor product grid of
c     Legendre nodes of order 2*k+1 while maintaining symmetries
c     with respect to flipping coordinate directions and
c     permuting the order of x,y,z coordinates.
c
c     This routine begins with a red-black kind of removal.
c     Let (i,j,k) be the indices of a standard tensor grid
c     point. This point is kept if i+j+k is odd, removed if
c     i+j+k is even.
c     
c      
c      
      implicit real *8 (a-h,o-z)
      integer k, npt, ipt2depth(3,*),idepth2pt(2*k+1,2*k+1,2*k+1)
      real *8 pts(3,*)

      real *8 u,v,w
      real *8, allocatable :: t(:)
      integer nperm, iperms(3,6), itemp(3)
      parameter (nperm = 6)

      data iperms / 1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1 /


      ifloor = 2*(2*k)*(2*k+1)*(2*k+2)/6

      write(*,*) 'ifloor', ifloor
      
      iseed = 1234

      dd = hkrand(iseed)

      itype = 0
      n = 2*k+1
      allocate(t(n))
      call legeexps(itype,n,t,u,v,w)

      npt = (2*k+1)**3

c      write(*,*) 'ifloor ',  ifloor
c      write(*,*) 'npt start ', npt

c     for now use idepth2pt to store in/out

c     start all in

      do iz = 1,2*k+1
         do iy = 1,2*k+1
            do ix = 1,2*k+1
               idepth2pt(ix,iy,iz) = 1
            enddo
         enddo
      enddo
      
c     do the red-black thinning first
      
      do iz = 1,2*k+1
         do iy = 1,2*k+1
            do ix = 1,2*k+1
               if (mod(ix+iy+iz,2) .eq. 0) then
                  npt = npt-1
                  idepth2pt(ix,iy,iz) = -1
               endif
            enddo
         enddo
      enddo

c      write(*,*) 'npt RB', npt


      ntry = 100
      
      
c     remove x,y,z coords at random while you can
      do iii = 1,ntry
         
         do i = 1,ntry
            ix = hkrand(0)*(k+1)
            iy = hkrand(0)*(k+1)
            iz = hkrand(0)*(k+1)

            if ( (ix .lt. 0) .or. (iy.lt. 0) .or. (iz .lt. 0) .or.
     1           (ix .gt. k) .or. (iy.gt. k) .or. (iz .gt. k)) then
c               write(*,*) 'uh oh '
               cycle
            endif

            ix2 = ix+k+1
            iy2 = iy+k+1
            iz2 = iz+k+1

            if (idepth2pt(ix2,iy2,iz2) .eq. -1) then
c     already gone
c               write(*,*) 'already gone'
               cycle
            endif

c     figure out how many are removed
            
            nrem = 48
            if (ix .eq. 0) nrem = nrem/2
            if (iy .eq. 0) nrem = nrem/2
            if (iz .eq. 0) nrem = nrem/2
            if ( (ix .eq. iy) .and. (ix .eq. iz)) then
c               write(*,*) 'all same'
               nrem = nrem/6
            else
               if ( (ix .eq. iy) .or. (ix .eq. iz) .or.
     1              (iy .eq. iz)) then
c                  write(*,*) 'two same'                  
                  nrem = nrem/2
               endif
            endif
c            write(*,*) ix,iy,iz,nrem
            npt = npt-nrem

c     remove them

            do ip = 1,nperm
               itemp(iperms(1,ip)) = ix
               itemp(iperms(2,ip)) = iy
               itemp(iperms(3,ip)) = iz
               do ixs = 1,2
                  do iys = 1,2
                     do izs = 1,2
                        ix2 = itemp(1) + k + 1
                        iy2 = itemp(2) + k + 1
                        iz2 = itemp(3) + k + 1
                        if (ixs .eq. 2) ix2 = k + 1 -itemp(1)
                        if (iys .eq. 2) iy2 = k + 1 -itemp(2)
                        if (izs .eq. 2) iz2 = k + 1 -itemp(3)
                        idepth2pt(ix2,iy2,iz2) = -1
                     enddo
                  enddo
               enddo
            enddo

            exit
            
         enddo
c         write(*,*) 'npt ', iii, npt
         if (npt .lt. ifloor + 48) exit
      enddo


c     get point ordering and such

      npt = 0
      do iz = 1,2*k+1
         do iy = 1,2*k+1
            do ix = 1,2*k+1
               if (idepth2pt(ix,iy,iz) .eq. 1) then
                  npt = npt+1
                  ipt2depth(1,npt) = ix
                  ipt2depth(2,npt) = iy
                  ipt2depth(3,npt) = iz
                  pts(1,npt) = t(ix)
                  pts(2,npt) = t(iy)
                  pts(3,npt) = t(iz)
                  idepth2pt(ix,iy,iz) = npt
               endif
            enddo
         enddo
      enddo

c      write(*,*) 'npt end', npt
               
      
      
      
      return
      end
      
