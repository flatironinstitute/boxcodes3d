c
c
c     routines for working with tensor products of points and
c     tensor product polynomial bases
c
c     this file takes the following conventions
c
c     a tensor grid of points is traversed with x on the inner
c        loop, y on the loop above that, and z above that.
c        e.g. for a 2x2x2 grid, we have the order:
c            (x1,y1,z1), (x2,y1,z1), (x1,y2,z1), (x2,y2,z1)
c     (x1,y1,z2), (x2,y1,z2), (x1,y2,z2), (x2,y2,z2)
c      
c     tensor product polynomials are numbered analagously
c        e.g. for full degree at order 3, we would have the order
c           1, x, x^2, y, xy, x^2y, y^2, xy^2, x^2 y^2, 
c           z, xz, x^2 z, yz, xyz, x^2 yz, y^2 z, xy^2 z, x^2 y^2 z,
c           z^2, xz^2, x^2 z^2, yz^2, xyz^2, x^2 y z^2,
c           y^2 z^2, xy^2 z^2, x^2 y^2 z^2
c        e.g. for total degree, we would have the order
c           1, x, x^2, y, xy, y^2, z, xz, yz, z^2
c
c      
c     mesh3d - like matlab's mesh generation function in 3d
c     split3d - separate components of an array of points in 3d
c     lens3d - compute the magnitude of each point in an array of
c                3d points
c     dists3d - compute pairwise distances between the 3d points
c     in two arrays

      subroutine mesh3d(x,nx,y,ny,z,nz,xyz)
      implicit real *8 (a-h,o-z)
      dimension x(*), y(*), z(*), xyz(3,*)

      ind = 0
      do iz = 1,nz
         do iy = 1,ny
            do ix = 1,nx
               ind = ind+1
               xyz(1,ind) = x(ix)
               xyz(2,ind) = y(iy)
               xyz(3,ind) = z(iz)
            enddo
         enddo
      enddo

      return
      end

      subroutine split3d(xyz,n,x,y,z)
      implicit real *8 (a-h,o-z)
      dimension x(*), y(*), z(*), xyz(3,*)

      do i = 1,n
         x(i) = xyz(1,i)
         y(i) = xyz(2,i)
         z(i) = xyz(3,i)
      enddo

      return
      end

      subroutine lens3d(xyz,n,dlens)
      implicit real*8 (a-h,o-z)
      dimension xyz(3,*), dlens(*)

      do i = 1,n
         dlens(i) = 0.0d0
         do j = 1,3
            dlens(i) = dlens(i) + xyz(j,i)**2
         enddo
         dlens(i) = sqrt(dlens(i))
      enddo

      return
      end

      subroutine dists3d(a,n,b,m,dists)
      implicit real*8 (a-h,o-z)
      dimension a(3,*), b(3,*), dists(n,m)

      do i = 1,m
         do j = 1,n
            dists(j,i) = 0
            do k = 1,3
               dists(j,i) = dists(j,i) + (a(k,j)-b(k,i))**2
            enddo
            dists(j,i) = sqrt(dists(j,i))
         enddo
      enddo

      return
      end

      subroutine buildperm3d(idimp,iflip,nq,cp,ipperm,icperm,
     1     icsign)
c
c
c     combines dimension permutation and flip information into
c     index permutation (for target points) and coefficient
c     permutation (for source densities) information. also returns
c     overall sign flips per coefficient of a source density
c
c     TODO: coefficient permutations and sign flips
c
c            
      implicit real*8 (a-h,o-z)
      dimension idimp(3), iflip(3), ipperm(*), icperm(*)
      dimension icsign(*)
      dimension indmap(nq,3), iskip(3), ilist(nq,nq,nq)
      dimension itemp(3)
      character cp

      iskip(1) = 1
      iskip(2) = nq
      iskip(3) = nq**2

      do ii = 1,3
         do jj = 1,nq
            if (iflip(ii) .eq. -1) then
               indmap(jj,ii) = nq-jj+1
            else
               indmap(jj,ii) = jj
            endif
         enddo
      enddo

      do iz = 1,nq
         jz = indmap(iz,idimp(3))
         do iy = 1,nq
            jy = indmap(iy,idimp(2))
            do ix = 1,nq
               jx = indmap(ix,idimp(1))
               ind1 = (iz-1)*nq**2 + (iy-1)*nq + ix
               ind2 = (jz-1)*iskip(idimp(3)) +
     1              (jy-1)*iskip(idimp(2)) +
     1              (jx-1)*iskip(idimp(1)) + 1
               ipperm(ind1) = ind2
            enddo
         enddo
      enddo

      if (cp .eq. 't' .or. cp .eq. 'T') then
c     total degree order
         ii = 0
         do iz = 1,nq
            do iy = 1,nq+1-iz
               do ix = 1,nq+1-iy+1-iz
                  ii = ii+1
                  ilist(ix,iy,iz) = ii
               enddo
            enddo
         enddo

         ii = 0
         do iz = 1,nq
            do iy = 1,nq+1-iz
               do ix = 1,nq+1-iy+1-iz
                  ii = ii+1
                  itemp(1) = ix
                  itemp(2) = iy
                  itemp(3) = iz
                  jx = itemp(idimp(1))
                  jy = itemp(idimp(2))
                  jz = itemp(idimp(3))

                  jj = ilist(jx,jy,jz)
                  
                  icperm(ii) = jj
                  icsign(ii) = 
     1                 (iflip(idimp(1)))**(jx-1)*
     2                 (iflip(idimp(2)))**(jy-1)*
     2                 (iflip(idimp(3)))**(jz-1)
               enddo
            enddo
         enddo
      endif
         
         
      return
      end
      
      
      subroutine alltargs3d(xq,nq,bs,xyzc,wc,wbtos,wstob)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     generate reference points for all targets in colleague,
c     big to small, and small to big interactions
c
      
      implicit real *8 (a-h,o-z)
      dimension wc(3,nq**3,27), wbtos(3,nq**3,56), wstob(3,nq**3,56)
      dimension xshift(3), yshift(3), zshift(3), xyzc(3)
      dimension xshiftbtos1(4), xshiftbtos2(2)
      dimension xshiftstob1(4), xshiftstob2(2)      
      dimension grid(3,nq**3), xq(*)

      nshift1 = 2
      nshift2 = 4

c     lowest corner of cube (each coordinate is smallest)
      
      xc = xyzc(1)
      yc = xyzc(2)
      zc = xyzc(3)      

c     get corresponding meshgrid

      call mesh3d(xq,nq,xq,nq,xq,nq,grid)

      ngrid = nq**3

c     colleagues are straightforward

      ind = 0
      do iz = -1,1
         do iy = -1,1
            do ix = -1,1
               ind = ind+1
               do i = 1,ngrid
                  wc(1,i,ind) = xc + ix*bs + grid(1,i)
                  wc(2,i,ind) = yc + iy*bs + grid(2,i)
                  wc(3,i,ind) = zc + iz*bs + grid(3,i)
               enddo
            enddo
         enddo
      enddo


c     stob and btos are analogous but a
c     little more complicated

      bsh = bs/2

      ind = 0
      do iz = -1,2
         do iy = -1,2
            if (iz .eq. -1 .or. iz .eq. 2 .or.
     1           iy .eq. -1 .or. iy .eq. 2) then
c     all x positions possible
               do ix = -1,2
                  ind = ind+1
                  do i = 1,ngrid
                     wbtos(1,i,ind) = xc + ix*bsh + grid(1,i)/2
                     wbtos(2,i,ind) = yc + iy*bsh + grid(2,i)/2
                     wbtos(3,i,ind) = zc + iz*bsh + grid(3,i)/2

                     wstob(1,i,ind) = xc + (ix-1)*bs + grid(1,i)*2
                     wstob(2,i,ind) = yc + (iy-1)*bs + grid(2,i)*2
                     wstob(3,i,ind) = zc + (iz-1)*bs + grid(3,i)*2
                  enddo
               enddo
            else
c     only outer x positions possible
               do ix = -1,2,3
                  ind = ind+1
                  do i = 1,ngrid
                     wbtos(1,i,ind) = xc + ix*bsh + grid(1,i)/2
                     wbtos(2,i,ind) = yc + iy*bsh + grid(2,i)/2
                     wbtos(3,i,ind) = zc + iz*bsh + grid(3,i)/2

                     wstob(1,i,ind) = xc + (ix-1)*bs + grid(1,i)*2
                     wstob(2,i,ind) = yc + (iy-1)*bs + grid(2,i)*2
                     wstob(3,i,ind) = zc + (iz-1)*bs + grid(3,i)*2
                  enddo
               enddo
            endif
         enddo
      enddo
         
      return
      end

      subroutine tensrefpts3d(xq,nq,bs,xyzc,wc,wbtos,wstob)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     generate reference points for the limited subset of
c     points in colleague, big-to-small, and small-to-big
c     interactions which can be used to obtain the other
c     interactions by symmetries
c
      
      implicit real *8 (a-h,o-z)
      dimension wc(3,nq**3,4), wbtos(3,nq**3,3), wstob(3,nq**3,3)
      dimension xshift(3), yshift(3), zshift(3), xyzc(3)
      dimension grid(3,nq**3), xq(*)

c     lowest corner of cube (each coordinate is smallest)
      
      xc = xyzc(1)
      yc = xyzc(2)
      zc = xyzc(3)      

c     get corresponding meshgrid

      call mesh3d(xq,nq,xq,nq,xq,nq,grid)

      ngrid = nq**3
      
      do i = 1,ngrid
         wc(1,i,1) = xc + grid(1,i)
         wc(2,i,1) = yc + grid(2,i)
         wc(3,i,1) = zc + grid(3,i)
      enddo


      xshift(1) = -1
      xshift(2) = 0
      xshift(3) = 0
      yshift(1) = -1
      yshift(2) = -1
      yshift(3) = 0
      zshift(1) = -1
      zshift(2) = -1
      zshift(3) = -1

      bsh = bs/2
      
      do ii = 1,3
         do i = 1,ngrid
            wc(1,i,ii+1) = xc + xshift(ii)*bs + grid(1,i)
            wc(2,i,ii+1) = yc + yshift(ii)*bs + grid(2,i)
            wc(3,i,ii+1) = zc + zshift(ii)*bs + grid(3,i)          

            wbtos(1,i,ii) = xc + xshift(ii)*bsh + grid(1,i)/2
            wbtos(2,i,ii) = yc + yshift(ii)*bsh + grid(2,i)/2
            wbtos(3,i,ii) = zc + zshift(ii)*bsh + grid(3,i)/2

            wstob(1,i,ii) = xc + (xshift(ii)-1)*bs + grid(1,i)*2
            wstob(2,i,ii) = yc + (yshift(ii)-1)*bs + grid(2,i)*2
            wstob(3,i,ii) = zc + (zshift(ii)-1)*bs + grid(3,i)*2

         enddo
      enddo
         


      return
      end

