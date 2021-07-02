c
c
c
c
      subroutine h3d_vslp(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(3),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + (x(3)-y(3))**2)

      f = exp(ima*zpars(1)*rr)/rr

      return
      end

c
c
      subroutine h3d_vslp_gradx(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(3),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,zk

      zk = zpars(1)

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + (x(3)-y(3))**2)

      f = -(y(1)-x(1))*(1.0d0-ima*zk*rr)*exp(ima*zpars(1)*rr)/rr**3

      return
      end

c      
c      
c
c
c
c
      subroutine h3d_vslp_grady(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(3),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,zk

      zk = zpars(1)

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + (x(3)-y(3))**2)

      f = -(y(2)-x(2))*(1.0d0-ima*zk*rr)*exp(ima*zpars(1)*rr)/rr**3

      return
      end
c
c
c
c
      subroutine h3d_vslp_gradz(x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(3),y(3),dpars(*)
      complex *16 zpars(*),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(*)
      complex *16 f,zk

      zk = zpars(1)

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + (x(3)-y(3))**2)

      f = -(y(3)-x(3))*(1.0d0-ima*zk*rr)*exp(ima*zpars(1)*rr)/rr**3

      return
      end

c      
c      
c
c
c
c
