
      implicit real *8 (a-h,o-z)
      real *8, allocatable :: w(:), rs(:)
      real *8, allocatable :: eplot(:),eplot2(:)
      call prini(6,13)

      cp = 40.0d0

c     parameters for Eaton lens. a is support, b is max
      
      a = 0.45d0
      b = 2.0d0

      lenw = 1000000
      lenv = 10000
      allocate(w(lenw))

c     recommend choosing r0 slightly smaller than 0 and
c     rend a little bigger than the largest you'll call it
c     with
      
      r0 = -0.05d0
      rend = sqrt(0.375d0) + 0.05d0
      call eatonprol_form(ier,a,b,cp,r0,rend,w,lenw,lused,keep)


      rpa = 0.0d0
      rpb = sqrt(0.375d0)

      nplot = 4000

      rh = (rpb-rpa)/(nplot-1)

      allocate(rs(nplot),eplot(nplot),eplot2(nplot))
      
      do i = 1,nplot
         r = rpa + rh*(i-1)
         call eaton(r,a,b,eplot(i))
         rs(i) = r
      enddo

      call eatonprol_eval(rs,nplot,w,eplot2)
      
      open(unit=33,file='eaton_interp.txt')

      do i = 1,nplot
         r = ra + rh*(i-1)
         write(33,*) r, eplot(i), eplot2(i)
      enddo

      close(33)
         
      stop
      end

      
