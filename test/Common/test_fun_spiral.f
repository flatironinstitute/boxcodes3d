      implicit real *8 (a-h,o-z)
      real *8, allocatable :: xyzs(:,:),f(:)
      call prini(6,13)

      nlat = 300
      nsc = 1

 1211 format(4(2x,e11.5))
      ntarg = nlat*nlat*nsc
      allocate(xyzs(3,ntarg),f(ntarg))
      rsig = 0.04d0
      i = 0
      do ilat=1,nsc
        if(ilat.eq.1) z = 0.0d0
        if(ilat.eq.2) z = 0.3d0
        if(ilat.eq.3) z = 0.6d0
        do jlat=1,nlat
          do llat = 1,nlat
            i = i+1
            xyzs(1,i) = -3 + 6*(jlat-1.0d0)/(nlat-1.0d0)
            xyzs(2,i) = 3 - 6*(llat-1.0d0)/(nlat-1.0d0)
            xyzs(3,i) = z 
            call fspiral1(nd,xyzs(1,i),rsig,zpars,ipars,f(i))
            write(44+ilat,1211) xyzs(1,i),xyzs(2,i),xyzs(3,i),f(i)


          enddo
        enddo
      enddo


       

      stop
      end


      subroutine fspiral1(nd,xyz,rsig,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
c
c   expects box to be in [-3,3] (there is a buffer just for
c   safety
c
c   rsig is the standard deviation of the gaussian used
c   to mollify the curve.
c
c   f(x,y,z) = \int_{\Gamma} exp(-((x-x(t))**2 + (y-y(t))**2 + 
c     (z-z(t))**2)/rsig**2) dl(t) 
c
c   where $\Gamma$ is the union of two spirals
c   
c
 
      real *8 xyz(3),rsig
      complex *16 zpars
      real *8 f,rint1,rint2

      external fker_spr1

      done = 1
      pi = atan(done)*4


      
      nd = 1
      ier = 0
      a = 0
      b = 2*pi
      eps = 1.0d-12
      m = 26
      maxrec = 200
      numint = 100000
      rint1 = 0
      ier = 0
      call adapgaus_new(ier,a,b,fker_spr1,xyz,rsig,m,eps,rint1,
     1   maxrec,numint)
      
      rint2 = 0
      call adapgaus_new(ier,a,b,fker_spr1,xyz,rsig,m,eps,rint2,
     1   maxrec,numint)
      
      f = rint1 + rint2*0

      
      return
      end




      subroutine fker_spr1(t,xyz,rsig,f)
      implicit real *8 (a-h,o-z)
      real *8 xyz(3),f,rsig,t

      rdec = 0.3d0 
      nosc = 10
      rz = 0.001d0*0

      rexp = exp(-rdec*t)
      
      xt = rexp*cos(nosc*t)
      yt = rexp*sin(nosc*t)
      zt = rz*t


      dxdt = -rdec*xt - nosc*yt
      dydt = -rdec*yt + nosc*xt
      dzdt = rz

      dst = sqrt(dxdt**2 + dydt**2 + dzdt**2)

      dx = xt-xyz(1)
      dy = yt-xyz(2)
      dz = zt-xyz(3)

      rr = dx**2 + dy**2 + dz**2

      rr2= rr/rsig**2

      if(rr2.lt.50) f = exp(-rr2)*dst
      if(rr2.ge.50) f = 0

      return
      end





      
