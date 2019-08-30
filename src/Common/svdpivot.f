cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        This is the end of the debugging code, and the beginning
c        of the SVD code proper. This file contains two user-
c        callable subroutines: svdpivot and svdpivot2. Following 
c        is a brief description of these subroutines.
c 
c   svdpivot - constructs the singular value decomposition of the 
c        user-specified matrix a, preferabluy under the assumption 
c        that the rank of a is considerably smaller than its
c        dimensionalities n,m. 
c 
c   svdpivot2 - constructs the singular values of the user-specified 
c        matrix a, preferabluy under the assumption that the rank of 
c        a is considerably smaller than its dimensionalities n,m. 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine svdpivot2(ier,a,n,m,s,ncols,eps,
     1      w,lw,ltot)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),s(1),w(1)
c 
c        this subroutine constructs the singular values
c        of the user-specified matrix a, under the assumption
c        that the rank of a is considerably smaller than its
c        dimensionalities n,m. This subroutine is a much 
c        castrated version of the subroutine svdpivot, in the
c        sense that the latter produces the whole singular value
c        decomposition of a.
c 
c                   input parameters:
c 
c  a - the matrix to be svd-decomposed; note that it is NOT damaged
c       by this subroutine in any way
c  n,m - dimensionalities of a
c  eps - the accuracy with which the calculation will be performed;
c        should not be less than the machine precision; recommended
c        value: 1.0d-16
c  lw - the number of real *8 words provided by the user in the work
c        array w. must be at least m+n+10 (otherwise, the subroutine
c        can not even estimate the actual memory requirements).
c 
c                   output parameters:
c 
c  ier - error return code
c      ier=0 means successful execution
c      ier=4 means that the amount of space lw allocated supplied
c            in the work array w is insufficient
c      ier=1024 means that the amount of space lw allocated supplied
c            in the work array w is grossly insufficient
c  s - the array of singular values of a
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  ltot - the total amount of space (in terms of real *8 words)
c        needed by the subroutine in the array w
c 
c                   work arrays:
c 
c  w - must be sufficiently long. 
c 
c 
c        using gram-schmidt process with pivoting, decompose
c        the matrix v in the form
c 
c          v=w t^*,
c 
c        with  w an orthogonal matrix of minimum rank,
c        and t a triangular matrix of dimensionality ncols * ncols
c 
c 
c        . . . using gram-schmidt process with pivoting, decompose
c              the matrix a in the form
c 
c                    a=U  V^*,
c 
c              with  u an orthogonal matrix of minimum rank,
c              and v whatever it wants to be
c 
        if( (lw .gt. m) .and. (lw. gt. n)) goto 1200
        ier=1024
        return
 1200 continue
c 
c        allocate memory
c
        iv=1
        lv=n*m+2
c      
        iu=iv+lv
        lu=n*m+2
c
        irnorms=iu+lu  
c
        lrnorms=n
        if(m .gt. n) lrnorms=m
        lrnorms=lrnorms+10
c 
        ltot=irnorms+lrnorms
        if(ltot .le. lw) goto 1300
        ier=512
        return
 1300 continue
c 
        ifpivot=1
        call svdpiv0(a,n,m,w(iu),w(iv),ncols,w(irnorms),eps,ifpivot)
c 
c       allocate memory for the subsequent operations
c 
        lv=m*ncols+2
c
        iw=iv+lv
        liw=m*ncols+10
c 
        it=iw+liw
        lt=ncols**2+10

        npn=48*ncols+100
        if(lt .lt. npn) lt=npn
c 
        ida=it+lt
        lda=ncols+10 
c
        isuba=ida+lda
        lsuba=ncols+10       
c
        irnorms=isuba+lsuba
        lrnorms=n
        if(m .gt. n) lrnorms=m
        lrnorms=lrnorms+10
c 
        ltot=irnorms+lrnorms
        if(ltot .le. lw) goto 1400
        ier=4
        return
 1400 continue
c 
c       conclude the construction of the singular value decomposition
c 
        call svdpiv2(w(iv),w(iw),w(it),n,m,ncols,
     1    w(irnorms),eps,s,w(ida),w(isuba) )
c
        return
        end
c 
c 
c 
c 
c 
        subroutine svdpiv2(v,w,t,n,m,ncols,rnorms,
     1    eps,s,da,suba)
        implicit real *8 (a-h,o-z)
        save
        dimension v(m,1),w(m,1),t(ncols,1),rnorms(1),
     1      suba(1),da(1),s(1)
c 
c        this subroutine performs the second part of the 
c        construction of the singular values of the 
c        user-specified matrix a, under the assumption
c        that the rank of a is considerably smaller than its
c        dimensionalities n,m. 
c 
c                   input parameters:
c 
c  v - the matrix produced by a preceding call to svdpiv0 (see)
c  n,m - dimensionalities of a
c  eps - the accuracy with which the calculation will be performed;
c        should not be less than the machine precision; recommended
c        value: 1.0d-16
c 
c                   output parameters:
c 
c  ier - error return code
c      ier=0 means successful execution
c      ier=4 means that the amount of space lw allocated supplied
c            in the work array w is insufficient
c      ier=1024 means that the amount of space lw allocated supplied
c            in the work array w is grossly insufficient
c  s - the array of singular values of a
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  ltot - the total amount of space (in terms of real *8 words)
c        needed by the subroutine in the array w
c 
c                   work arrays:
c 
c  w - must be sufficiently long. 
c 
c 
c        using gram-schmidt process without pivoting, decompose
c        the matrix v in the form
c 
c          v=w t^*,
c 
c        with  w an orthogonal matrix of minimum rank,
c        and t a triangular matrix of dimensionality ncols * ncols
c 
        ifpivot=0
        eps2=0
        call svdpiv0(v,m,ncols,w,t,ncols2,rnorms,eps2,ifpivot)
c 
c        use jacobi rotations to sv-decompose the matrix t
c 
        call prinf('m=*',m,1)


        call svdpiv_info_out(ind)
        if( (ncols .le. 7) .or. (ind .eq. 1) ) then
            ifuv=0
            call d2mudv(t,ncols,u27,v27,eps,ifuv,numrows)
c
            do 1800 i=1,ncols
            s(i)=t(i,1)
 1800 continue
c
            return
        endif
c
c        bidiagonalize the triangular matrix T
c
cccc        call svd_rbidiag2(t,ncols,n,m)
        call svd_rbidiag2(t,ncols)
c
c
        do 2400 i=1,ncols-1
        suba(i)=t(i+1,i)
        da(i)=t(i,i)
 2400 continue
c 
        da(ncols)=t(ncols,ncols)
c 
c        SVD the real matrix
c 
cc        call invsing_bidiag_vects(ier,ncols,da,suba,
cc     1      s,w,v,t)

        call prinf('in svdpivot2 before _vects, ncols=*',ncols,1)

        call qleigen_bidiag_vects(ier,ncols,da,suba,
     1      s,w,v,t)

c     
        return
        end
c 
c 
c 
c 
c 
        subroutine svdpivot(ier,a,n,m,u,v,s,ncols,eps,
     1      w,lw,ltot)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),u(n,1),v(m,1),s(1),w(1)
c 
c        this subroutine constructs the singular value decomposition
c        of the user-specified matrix a, under the assumption
c        that the rank of a is considerably smaller than its
c        dimensionalities n,m. the output of this subroutine
c        consists of two matrices u,v, and a vector s such that
c 
c                 a = u d v^*                                   (1)
c 
c        with d the diagonal matrix with the vector s on the
c        diagonal, and u, v two matrices whose columns are
c        orthonormal. the reason for the existence of this
c        subroutine is the fact that on exit, the matrices u,v
c        have dimensionalities (n,ncols), (m,ncols), respectively,
c        with ncols the numerical rank of a.
c 
c                   input parameters:
c 
c  a - the matrix to be svd-decomposed; note that it is NOT damaged
c       by this subroutine in any way
c  n,m - dimensionalities of a
c  eps - the accuracy with which the calculation will be performed;
c        should not be less than the machine precision; recommended
c        value: 1.0d-16
c  lw - the number of real *8 words provided by the user in the work
c        array w. must be at least m+n+10 (otherwise, the subroutine
c        can not even estimate the actual memory requirements).
c 
c                   output parameters:
c 
c  ier - error return code
c      ier=0 means successful execution
c      ier=4 means that the amount of space lw allocated supplied
c            in the work array w is insufficient
c      ier=1024 means that the amount of space lw allocated supplied
c            in the work array w is grossly insufficient
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  v - the second matrix in (1). it is dimensioned v(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so v should be dimensioned v(m,n) by the user. however,
c        on exit only the first ncols columns of v are meaningful
c  s - the array of singular values of a
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  ltot - the total amount of space (in terms of real *8 words)
c        needed by the subroutine in the array w
c 
c                   work arrays:
c 
c  w - must be sufficiently long. 3*n*m+ 2*(min(m,n))**2 +2*n+2*m+100
c        are always sufficient. however, must of the time it is much
c        less than that, and depends on the rank of a. the actual
c        amount of space available in w is communicated to the
c        subroutine via the parameter lw (see above), and if this
c        amount is insufficient, the output error code ier is set
c        accordingly, and the memory requirements of the subroutine
c        are returned in the parameter ltot (see above).
c 
c        . . . using gram-schmidt process with pivoting, decompose
c              the matrix a in the form
c 
c                    a=U  V^*,
c 
c              with  u an orthogonal matrix of minimum rank,
c              and v whatever it wants to be
c 
        if( (lw .gt. m) .and. (lw. gt. n)) goto 1200
        ier=1024
        return
 1200 continue
c 
        ifpivot=1
        call svdpiv0(a,n,m,u,v,ncols,w,eps,ifpivot)
c 
c       allocate memory for the subsequent operations
c 
        iw=1
        liw=m*ncols+10
c 
        iu2=iw+liw
        lu2=ncols**2+10
c 
        iv2=iu2+lu2
        lv2=ncols**2+10
c 
        it=iv2+lv2
        lt=lv2

cccc        lt=lt*2+20*ncols
c
cccc        npn=24*ncols+100
        npn=48*ncols+100
        if(lt .lt. npn) lt=npn
c 
        irnorms=it+lt 
        lrnorms=n
        if(m .gt. n) lrnorms=m
        lrnorms=lrnorms+10
c 
        iww=irnorms+lrnorms 
        lww=n
        if(m .gt. n) lww=m
        lww=lww+10
c 
        ltot=iww+lww
        if(ltot .le. lw) goto 1400
        ier=4
        return
 1400 continue
c 
c       conclude the construction of the singular value decomposition
c 
        iffirst=0
        call svdpiv1(a,u,v,w(iw),w(it),n,m,ncols,
     1    w(irnorms),eps,w(iv2),w(iu2),w(iww),iffirst,s)
c 
c       copy the right singluar matrix into the array v
c 
        jj=0
        do 1800 i=1,ncols
        do 1600 j=1,m
        jj=jj+1
        v(j,i)=w(jj)
 1600 continue
 1800 continue
c 
c       now, copy the singular values into the array s
c 
        do 2000 i=1,ncols
        s(i)=w(it+i-1)
 2000 continue
        return
        end
c
c
c
c
c
        subroutine svdpiv1(a,u,v,w,t,n,m,ncols,rnorms,
     1    eps,v2,u2,ww,iffirst,s)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),u(n,1),v(m,1),w(m,1),t(1),rnorms(1),
     1      v2(1),u2(1),ww(1)
c 
c        using gram-schmidt process with pivoting, decompose
c        the matrix a in the form
c 
c          a=U  V^*,
c 
c        with  u an orthogonal matrix of minimum rank,
c        and v whatever it wants to be
c 
        if(iffirst .eq. 0) goto 1200
        ifpivot=1
        call svdpiv0(a,n,m,u,v,ncols,rnorms,eps,ifpivot)
 1200 continue
c 
c        using gram-schmidt process without pivoting, decompose
c        the matrix v in the form
c 
c          v=w t^*,
c 
c        with  w an orthogonal matrix of minimum rank,
c        and t a triangular matrix of dimensionality ncols * ncols
c 
        ifpivot=0
        eps2=0
        call svdpiv0(v,m,ncols,w,t,ncols2,rnorms,eps2,ifpivot)
c
        call svdpiv_info_out(ind)

cccc        ind=1
        if(ind .eq. 1) goto 1600
c
        if(ncols .lt. 7) goto 1600
c
c        bidiagonalize the triangular matrix T
c
        call svd_rbidiag(t,ncols,u,w,n,m)
c
c        diagonalize the bidiagonal matrix
c
        call svd_bidiag_real_qr(t,ncols,u2,v2,rnorms,ww,s)
c
        goto 1800
 1600 continue
c 
c        use jacobi rotations to sv-decompose the matrix t
c 
        ifuv=1
        call d2mudv(t,ncols,u2,v2,eps/10000/10000,ifuv,numrows)

 1800 continue
c 
c       multiply the matrix u from the right by u2, and the
c       matrix w from the right by v2^*
c 
        call svdfinm(u,u2,w,v2,n,m,ncols,ww)
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_bidiag_real_qr(t,ncols,u3,v3,suba,da,rlams)
        implicit real *8 (a-h,o-z)
        save
        dimension u3(ncols,1),v3(ncols,1),suba(1),da(1),rlams(1)
c 
        real *8 t(ncols,ncols)

cccc        dimension u7(100000),v7(100000),t7(1000 000)
c 
c       construct the diagonal and subdiagonal
c 
        do 1400 i=1,ncols-1
        suba(i)=t(i+1,i)
        da(i)=t(i,i)
 1400 continue
c 
        da(ncols)=t(ncols,ncols)
c 
c        SVD the real matrix
c 
cc        call invsing_bidiag_vects(ier,ncols,da,suba,
cc     1      rlams,u3,v3,t)




        call qleigen_bidiag_vects(ier,ncols,da,suba,
     1      rlams,u3,v3,t)



        call prin2('after invsing_bidiag_vects, rlams=*',rlams,ncols)
 
cccc        stop



        do 2200 j=1,ncols
c
        da(j)=rlams(j)
 2200 continue
c
        do 2600 i=1,ncols
        do 2400 j=1,ncols
c
        t(j,i)=v3(i,j)
 2400 continue
 2600 continue
c
        call qleigen_copy(t,v3,ncols**2)
c 
 2800 continue
c     
        do 3200 i=1,ncols
c
        t(i,1)=da(i)
 3200 continue
c
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_rbidiag2(a,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),b(3),ub(2,2)
c 
c        starting with the left lower corner, use Jacobi
c        rotations to bidiagonalize a
c 
        do 2000 i=n,3,-1
c 
        do 1800 j=1,i-2
c 
c       find the 2*2-jacobi matrix to be applied to the columns
c 
        b(1)=a(i,j+1)
        b(2)=a(i,j)
c  
        call svd_rotfnd(b,ub)
c 
c       transform the columns
c 
cccc        call svd_cols_rotate(ub,a,n,j+1,j)
        call svd_cols_rotate4(ub,a,n,j+1,j,j,i)
c 
c       find the 2*2-jacobi matrix to be applied to the rows
c 
        b(1)=a(j+1,j+1)
        b(2)=a(j,j+1)
c 
        call svd_rotfnd(b,ub)
c 
cccc        call svd_rows_rotate(ub,a,n,j+1,j)
        call svd_rows_rotate4(ub,a,n,j+1,j,1,j+1)
 1800 continue
c 
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_rbidiag(a,n,u,v,nn,mm)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),b(3),ub(2,2),u(nn,n),v(mm,n)
c 
c        starting with the left lower corner, use Jacobi
c        rotations to bidiagonalize a
c 
        do 2000 i=n,3,-1
c 
        do 1800 j=1,i-2
c 
c       find the 2*2-jacobi matrix to be applied to the columns
c 
        b(1)=a(i,j+1)
        b(2)=a(i,j)
c
        call svd_rotfnd(b,ub)
c 
c       transform the columns
c 
cc        call svd_cols_rotate(ub,a,n,j+1,j)
        call svd_cols_rotate4(ub,a,n,j+1,j,j,i)
c 
c        apply the inverse of the column transformation
c        to the rows of matrix v
c 
        call svd_cols_rotate(ub,v,mm,j+1,j)
c 
c       find the 2*2-jacobi matrix to be applied to the rows
c 
        b(1)=a(j+1,j+1)
        b(2)=a(j,j+1)
c 
        call svd_rotfnd(b,ub)
c 
cccc        call svd_rows_rotate(ub,a,n,j+1,j)
        call svd_rows_rotate4(ub,a,n,j+1,j,1,j+1)
c 
c        apply the inverse of the row transformation
c        to the columns of matrix u
c 
        call svd_cols_rotate(ub,u,nn,j+1,j)
c 
 1800 continue
c 
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_rotfnd(a,u)
        implicit real *8 (a-h,o-z)
        save
        dimension a(2),u(2,2)
c 
        u21=-a(2)
        u22=a(1)
c 
        d=sqrt(u22*u22+u21*u21)
c 
        if(d .eq. 0) then
c 
            u(2,2)=1
            u(1,2)=0
            u(1,1)=1
            u(2,1)=0
            return
        endif
c 
        u(2,2)=u22/d
        u(2,1)=u21/d
c 
        u(1,1)=-u(2,2)
        u(1,2)=u(2,1)
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_colmult(a,n,ii,coef)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n)
c 
        do 1200 i=1,n
c 
        a(i,ii)=a(i,ii)*coef
 1200 continue
c 
        return
c 
c 
c 
c 
        entry svd_rowmult(a,n,ii,coef)
c 
        do 2200 i=1,n
c 
        a(ii,i)=a(ii,i)*coef
 2200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_cols_rotate4(u,a,n,ii,jj,i1,i2)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),u(2,2)
c 
cccc        do 1200 i=1,n
        do 1200 i=i1,i2
c 
        d1=u(1,1)*a(i,ii)+u(1,2)*a(i,jj)
        d2=u(2,1)*a(i,ii)+u(2,2)*a(i,jj)
c 
        a(i,ii)=d1
        a(i,jj)=d2
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_cols_rotate(u,a,n,ii,jj)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),u(2,2)
c 
        do 1200 i=1,n
c 
        d1=u(1,1)*a(i,ii)+u(1,2)*a(i,jj)
        d2=u(2,1)*a(i,ii)+u(2,2)*a(i,jj)
c 
        a(i,ii)=d1
        a(i,jj)=d2
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_rows_rotate4(u,a,n,ii,jj,i1,i2)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),u(2,2)
c 
        do 1200 i=1,n
c 
        d1=u(1,1)*a(ii,i)+u(1,2)*a(jj,i)
        d2=u(2,1)*a(ii,i)+u(2,2)*a(jj,i)
c 
        a(ii,i)=d1
        a(jj,i)=d2
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine svd_rows_rotate(u,a,n,ii,jj)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),u(2,2)
c 
        do 1200 i=1,n
c 
        d1=u(1,1)*a(ii,i)+u(1,2)*a(jj,i)
        d2=u(2,1)*a(ii,i)+u(2,2)*a(jj,i)
c 
        a(ii,i)=d1
        a(jj,i)=d2
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine svdfinm(u,u2,w,v2,n,m,ncols,ww)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,ncols),w(m,ncols),u2(ncols,ncols),
     1      v2(ncols,ncols),ww(1)
c 
c       multiply u from the right by u2
c 
        do 1600 i=1,n
        do 1400 j=1,ncols
c 
        d=0
        do 1200 k=1,ncols
        d=d+u(i,k)*u2(k,j)
 1200 continue
        ww(j)=d
 1400 continue
        do 1500 j=1,ncols
        u(i,j)=ww(j)
 1500 continue
 1600 continue
c 
c         multiply w from the right by the adjoint of v2
c 
        do 2600 i=1,m
        do 2400 j=1,ncols
c 
        d=0
        do 2200 k=1,ncols
        d=d+w(i,k)*v2(j,k)
 2200 continue
        ww(j)=d
 2400 continue
        do 2500 j=1,ncols
        w(i,j)=ww(j)
 2500 continue
 2600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine svdpiv0(a,n,m,b,v,ncols,rnorms,eps,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),b(n,m),v(m,1),rnorms(1)
c 
c       this matrix applies the compressing gram-schmidt
c       process to the matrix a, obtaining its decomposition
c       in the form
c 
c             a=b v^T,                                        (1)
c 
c       with the matrices b, v having dimensionalities
c       b(n,ncols), v(m,ncols), respectively. the reason for the
c       existence of this subroutine is the hope that the
c       dimensionality ncols, determined by this subroutine,
c       is comsiderably lower than either n or m
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c 
c                     output parameters:
c 
c  b - the matrix in (1). note that the matrix b has to be dimensioned
c        b(n,m), but on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  v - the second matrix in (1). it is dimensioned v(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so v should be dimensioned v(m,n) by the user. however,
c        on exit only the first ncols columns of v are meaningful
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c        . . . copy the user-supplied matrix a into b
c 
        do 1400 i=1,m
        do 1200 j=1,n
        b(j,i)=a(j,i)
 1200 continue
 1400 continue
c 
c        apply the gram-schmidt proces (with pivoting) to b
c 
         call svdgrm(b,n,m,rnorms,eps,ncols,ifpivot)
c 
c        project the original matrix on the obtained orthogonal
c        columns
c 
        do 2400 i=1,ncols
        do 2200 j=1,m
c 
        call svdscapr(a(1,j),b(1,i),n,prod)
        v(j,i)=prod
 2200 continue
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine svdgrm(b,n,m,rnorms,eps,ncols,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        dimension b(n,m),rnorms(1)
c 
c       this subroutine applies a pivoted double gram-schmidt
c       procedure to the matrix b, resulting in a collection of
c       orthogonal vectors. the number of these vectors is
c       equal to the numerical rank of b, given the precision
c       of computations eps, which is supplied by the user.
c 
c                    input paramneters:
c 
c  b - the matrix to be gram-schmidt'ed. it is destroyed by this
c       subroutine
c  n,m - dimensionalities of the matrix b
c  eps - the machine accuracy
c 
c                     output parameters:
c 
c  b - the matrix of gram-schmidt vectors of the matrix a. note
c        that on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrix b on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c       . . . prepare the array of values of norms
c             of columns
c 
        done=1
        dtot=0
        do 1400 i=1,m
        d=0
        do 1200 j=1,n
        d=d+b(j,i)**2
 1200 continue
        rnorms(i)=sqrt(d)
        dtot=dtot+d
 1400 continue
c 
c       . . . conduct gram-schmidt iterations
c 
        thresh=dtot*eps**2
        do 4000 i=1,m
c 
c       find the pivot
c 
         if(ifpivot .eq. 0) goto 2700
        ipivot=i
        rn=rnorms(i)
c 
        do 2200 j=i+1,m
        if(rnorms(j) .le. rn) goto 2200
        rn=rnorms(j)
        ipivot=j
 2200 continue
 2400 continue
c 
c       put the column number ipivot in the i-th place
c 
        do 2600 j=1,n
        d=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=d
 2600 continue
c 
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
 2700 continue
c 
c       orthogonalize the i-th column to all preceeding ones
c 
        if(i .eq. 1) goto 2790
        do 2780 j=1,i-1
c 
        call svdscapr(b(1,i),b(1,j),n,d)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*d
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call svdscapr(b(1,i),b(1,i),n,d)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        if( (ifpivot .eq. 1) .and.
     1      (d .lt. thresh) ) return
        ncols=i
c 
        d=done/dsqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,m
c 
        if( (rnorms(j) .lt. thresh) .and. 
     1      (ifpivot .eq. 1) ) goto 3200
c 
        call svdscapr(b(1,i),b(1,j),n,d)
c 
         rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*d
        rrn=rrn+b(l,j)**2
 3000 continue
        rnorms(j)=dsqrt(rrn)
 3200 continue
 3400 continue
 4000 continue
  
 4200 continue
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine svdscapr(x,y,n,prod)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*y(i)
 1200 continue
        return
        end
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c            this is the end of the debuging code and the beginning
c            of the actual SVD routines
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
        subroutine d2mudv(a,n,u,v,eps,ifuv,numrows)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),u0(2,2),v0(2,2),u(n,n),v(n,n),
     1      b(2,2),rlam(2)
c 
c        this subroutine uses the classical jacobi iteration
c        to construct the singular value decomposition of a
c        general real matrix, of the form
c 
c          a = u s v,                                        (1)
c 
c        with u,v orthogonal matrices, and s diagonal, with
c        the diagonal elements monotonically decreasing.
c 
c                     input parameters:
c 
c  a - the n * n matrix to be sv-decomposed, is destroyed
c        by the subroutine
c  n - the dimensionality of a
c  eps - the accuracy to which the computations will be performed
c  ifuv - tells the subroutine whether the singular vectors are
c        to be computed, or only singular values
c 
c                     output parameters:
c 
c  a - the first column of a contains the (ordered) singular values
c  u,v - orthogonal matrices in (1)
c  numrows - the number of jacobi iterations (i.e. single rotations)
c        that the algorithm has taken
c 
c       . . . one sweep after another, conduct jacobi
c             iterations
c 
        maxsweep=10
        if(ifuv .eq. 0) goto 1150
        zero=0
        do 1100 i=1,n
        do 1050 j=1,n
        u(j,i)=zero
        v(j,i)=zero
 1050 continue
        u(i,i)=1
        v(i,i)=1
 1100 continue
 1150 continue
        numrows=0
        thresh=dsqrt(eps)
        thresh=dsqrt(thresh)
        thresh=dsqrt(thresh)
        do 3400 iijjkk=1,4
c 
        do 3000 ijk=1,maxsweep
        ifany=0
c 
        do 2000 i=1,n
        do 1800 j=1,n
         if(j .eq. i) goto 1800
         if( (dabs(a(i,j)) .lt. thresh) .and.
     1       (dabs(a(j,i)) .lt. thresh) ) goto 1800
        ifany=1
        numrows=numrows+1
c 
c       construct the svd of the (i,j)-th 2 * 2 submatrix
c 
        b(1,1)=a(i,i)
        b(1,2)=a(i,j)
        b(2,1)=a(j,i)
        b(2,2)=a(j,j)
        call d2m_oneuv(b,rlam,u0,v0)
c 
c       transform the big matrix
c 
c       . . . the rows
c 
        do 1200 l=1,n
        di=u0(1,1)*a(i,l)-u0(1,2)*a(j,l)
        dj=-u0(2,1)*a(i,l)+u0(2,2)*a(j,l)
        a(i,l)=di
        a(j,l)=dj
 1200 continue
c 
c       . . . and the columns
c 
        do 1400 l=1,n
        di=v0(1,1)*a(l,i)+v0(1,2)*a(l,j)
        dj=v0(2,1)*a(l,i)+v0(2,2)*a(l,j)
        a(l,i)=di
        a(l,j)=dj
 1400 continue
c 
c 
c       if the user so requested, adjust the matrices u,v
c 
        if(ifuv .eq. 0) goto 1800
c 
        do 1600 l=1,n
c 
c       the matrix u
c 
        di=u0(1,1)*u(l,i)-u0(1,2)*u(l,j)
        dj=-u0(2,1)*u(l,i)+u0(2,2)*u(l,j)
        u(l,i)=di
        u(l,j)=dj
c 
c       . . . and the columns
c 
        di=v0(1,1)*v(i,l)+v0(1,2)*v(j,l)
        dj=v0(2,1)*v(i,l)+v0(2,2)*v(j,l)
        v(i,l)=di
        v(j,l)=dj
 1600 continue
c 
 1800 continue
 2000 continue
c 
      if(ifany .eq. 0) goto 3200
 3000 continue
 3200 continue
        thresh=thresh**2
 3400 continue
c 
c        now, reorder the singular values of a to put them
c        in increasing order
c 
        do 4200 i=1,n
        a(i,1)=a(i,i)
 4200 continue
c 
       call d2msvsrt(a(1,1),n,u,v,a(1,2),ifuv)
        return
        end
c 
c 
c 
c 
c 
        subroutine d2msvsrt(s,n,u,v,iw,ifuv)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,n),v(n,n),s(1),iw(1)
c 
c        if some of the alleged singular values
c        are negative - change their sign, and that
c        of the corresponding left eigenvectors
c 
        do 1100 i=1,n
        if(s(i) .ge. 0) goto 1100
        s(i)=-s(i)
        if(ifuv .eq. 0) goto 1100
        do 1050 j=1,n
        v(i,j)=-v(i,j)
 1050 continue
 1100 continue
  
c        using the bubble sort, reorder the elements
c        of array s
c 
        do 1200 i=1,n
        iw(i)=i
 1200 continue
c 
        do 2000 i=1,n-1
c 
        ifmoved=0
        do 1600 j=1,n-i
        if(s(j). ge. s(j+1)) goto 1600
        jj=iw(j)
        iw(j)=iw(j+1)
        iw(j+1)=jj
c 
        d=s(j)
        s(j)=s(j+1)
        s(j+1)=d
 1600 continue
 2000 continue
        if(ifuv .eq. 0) return
c 
c       now, reorder the rows of u and columns of v
c 
        do 2400 i=1,n-1
        ii=iw(i)
c 
c       exchange the rows number i and ii in the matrix u
c 
        do 2200 j=1,n
        d=v(i,j)
        v(i,j)=v(ii,j)
        v(ii,j)=d
c 
c       exchange columns number i and ii in the matrix v
c 
        d=u(j,i)
        u(j,i)=u(j,ii)
        u(j,ii)=d
 2200 continue
c 
        do 2300 j=i+1,n
        if(iw(j) .ne. i) goto 2300
        iw(j)=ii
        goto 2350
 2300 continue
 2350 continue
c 
 2400 continue
        return
        end
  
c
c
c
c        
c
        subroutine d2m_oneuv(a,rlam,u,v)
        implicit real *8 (a-h,o-z)
        save
        dimension a(2,2),u(2,2),v(2,2),rlam(2),w1(2,2),w2(2,2)
c
        d1=atan2(a(2,1)+a(1,2),a(2,2)-a(1,1))
        d2=atan2(a(2,1)-a(1,2),a(2,2)+a(1,1))

        psi=(d2+d1)/2
        phi=(d2-d1)/2
c
        u(1,1)=cos(phi) 
        u(2,2)=cos(phi) 
        u(1,2)=sin(phi)
        u(2,1)=-sin(phi)
c
        v(1,1)=cos(psi) 
        v(2,2)=cos(psi) 
        v(1,2)=sin(psi)
        v(2,1)=-sin(psi)
c
        call d2m_prod2(u,a,w1)
        call d2m_prod2(w1,v,w2)
c
        rlam(1)=w2(1,1)
        rlam(2)=w2(2,2)
c
        u(1,2)=-u(1,2)
        u(2,1)=-u(2,1)
c       
        v(1,2)=-v(1,2)
        v(2,1)=-v(2,1)
c
        return
        end
c 
c 
c 
c 
c 
        subroutine d2m_oneuv_old(a,rlam,u,v)
        implicit real *8 (a-h,o-z)
        save
        dimension a(2,2),u(2,2),v(2,2),rlam(2),
     1      b(2,2),w(2,2),vstar(2,2),z(2,2),rlams(2,2)
        data eps/1.0d-10/,eps2/1.0d-1/,done/1.0d0/
c 
c       this subroutine produces the singular value decomposition
c       of a 2 * 2 real matrix a, so that on exit,
c 
c       a = u d v,
c 
c       with u and v orthogonal, and d a diagonal matrix with the
c       vector rlam=(rlam(1),rlam(2)) on the diagonal.
c 
c       . . . simmetrize a
c 
        den=a(2,2)+a(1,1)
        rnum=a(2,1)-a(1,2)
c 
c       if denominator is not too small, use exact formulae
c 
        if(dabs(den) .lt. eps*dabs(rnum)) goto 1200
        tgphi=rnum/den
c 
        if(dabs(tgphi) .lt. eps2) goto 1100
c 
        dd=done/(done+tgphi**2)
        alpha=dsqrt(dd)
        beta=dsqrt(done-dd)
        if(tgphi .lt. 0) beta=-beta
        goto 1400
 1100 continue
c 
        dd=tgphi**2/(1+tgphi**2)
c 
        beta=tgphi-tgphi**3/2
        alpha=dsqrt(done-beta**2)
c 
cccc        beta=dsqrt(dd)
cccc        alpha=dsqrt(done-dd)
cccc        if(tgphi .lt. 0) beta=-beta
c 
cccc        phi=datan(tgphi)
cccc        alpha=dcos(phi)
cccc        beta=dsin(phi)
        goto 1400
 1200 continue
c 
c       denominator is too small. use taylor series
c 
        alpha=den/rnum
        beta=1
 1400 continue
c 
c      calculate the simmetrizing matrix and the simmetric one
c 
        w(1,1)=alpha
        w(1,2)=beta
        w(2,1)=-beta
        w(2,2)=alpha
        call d2m_prod2(w,a,b)
c 
c       now, diagonalize the simmetrized matrix
c 
        den=b(2,2)-b(1,1)
        rnum=b(1,2)*2
c 
c       if denominator is not too small, use exact formulae
c 
        if(dabs(den) .lt. eps*dabs(rnum)) goto 1600
        tg2phi=-rnum/den
        if(dabs(tg2phi) .lt. 1.0d-3) goto 1500
c 
        dd=dsqrt(4/tg2phi**2+4)
        if(tg2phi .lt. 0) tgphi=(-2/tg2phi-dd)/2
        if(tg2phi .gt. 0) tgphi=(-2/tg2phi+dd)/2
c 
        phi=datan(tgphi)
        alpha=dcos(phi)
        beta=dsin(phi)
        goto 1800
 1500 continue
c 
        beta=tg2phi/2
        alpha=dsqrt(done-beta**2)
  
cccc        phi=datan(tg2phi)/2
cccc        alpha=dcos(phi)
cccc        beta=dsin(phi)
        goto 1800
 1600 continue
c 
c       denominator is too small. use taylor expansion
c 
        alpha=dsqrt((1-den/rnum)/2)
        beta=dsqrt((1+den/rnum)/2)
 1800 continue
c 
c       construct the diagonalizing matrix
c       and the resulting diagonal
c 
        v(1,1)=alpha
        v(1,2)=beta
        v(2,1)=-beta
        v(2,2)=alpha
c 
        vstar(1,1)=v(1,1)
        vstar(1,2)=v(2,1)
        vstar(2,1)=v(1,2)
        vstar(2,2)=v(2,2)
c 
        call d2m_prod2(v,w,u)
c 
c       finally, compute the resulting diagonal elements
c 
        call d2m_prod2(u,a,z)
        call d2m_prod2(z,vstar,rlams)
        rlam(1)=rlams(1,1)
        rlam(2)=rlams(2,2)
c 
        u(1,2)=-u(1,2)
        u(2,1)=-u(2,1)
        return
        end
c 
c 
c 
c 
c 
        subroutine d2m_prod2(a,b,c)
        implicit real *8 (a-h,o-z)
        save
        dimension a(2,2),b(2,2),c(2,2)
c 
        c(1,1)=a(1,1)*b(1,1)+a(1,2)*b(2,1)
        c(1,2)=a(1,1)*b(1,2)+a(1,2)*b(2,2)
        c(2,1)=a(2,1)*b(1,1)+a(2,2)*b(2,1)
        c(2,2)=a(2,1)*b(1,2)+a(2,2)*b(2,2)
        return
        end
c
c
c
c
c
        SUBROUTINE svdpiv_info_in(ind)
        save
        DATA ind7/0/
c
        ind7=ind
        RETURN
c
c
c
c
        entry svdpiv_info_out(ind)
c
        ind=ind7
        return
        end
