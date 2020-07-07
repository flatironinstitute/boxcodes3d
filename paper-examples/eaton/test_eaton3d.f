
      implicit real *8 (a-h,o-z)
      real *8, allocatable :: w(:), rs(:), w2(:)
      real *8, allocatable :: eplot(:),eplot2(:),eplot3(:)
      real *8 clege(1000), ab(2,3)
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

      allocate(rs(nplot),eplot(nplot),eplot2(nplot),eplot3(nplot))
      
      do i = 1,nplot
         r = rpa + rh*(i-1)
         call eaton(r,a,b,eplot(i))
         rs(i) = r
      enddo

      ifbell = 1
      rsupp = 0.5d0

c     belling performed at both ends. rsupp should generally
c     be set to something between a and the box half-width.
c     

      call cpu_time(t1)
      call eatonprol_eval(ifbell,rsupp,rs,nplot,w,eplot2)
      call cpu_time(t2)

      write(*,*) 'pts per second ', nplot/(t2-t1)


      
      nlege = 40
      lenw2 = nlege + 10
      allocate(w2(lenw2))
      
      call eaton_pre(a,b,w2,lenw2,nlege)
      call cpu_time(t1)
      do i = 1,nplot
         call eaton_post(rs(i),w2,eplot3(i))
      enddo
      call cpu_time(t2)
      
      write(*,*) 'pts per second ', nplot/(t2-t1)
      
      open(unit=33,file='eaton_interp.txt')

      do i = 1,nplot
         r = ra + rh*(i-1)
         write(33,*) r, eplot(i), eplot2(i), eplot3(i)
      enddo

      close(33)


         
      stop
      end

      
