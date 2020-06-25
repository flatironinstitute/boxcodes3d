c
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        This is the end of the debugging code, and the beginning of the
c        diagonalizing code proper.
c        This file contains three user-callable subroutines: 
c        nonsym_eigvals, nonsym_eigvects, and nonsym_ratios_check; the 
c        last of the three is not properly a diagonalizing routine, but
c        is a TESTING routine, meant to improve the user's self-esteem.
c        Following is a brief description of the three subroutines.
c
c   nonsym_eigvals - finds the spectrum of a (generally, non-Hermitian) 
c        complex matrix a. It uses a standard reduction to the upper 
c        Hessenberg form, followed by the QR procedure.
c
c   nonsym_eigvects - finds BOTH the eigenvalues and eigenvectors of a 
c       (generally, non-Hermitian) complex matrix a. In order to find 
c       the spectrum of a, it uses a standard reduction to the upper 
c       Hessenberg form, followed by the QR procedure. After that, it 
c       finds the eigenvectors via a straightforward inverse power 
c       scheme.
c
c   nonsym_ratios_check - tests the "quality" of the putative 
c       eigenvectors u(1,i) and roots(i). Provided purely for the user's
c       edification.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine nonsym_ratios_check(a,u,n,roots,errs,z)
        implicit real *8 (a-h,o-z) 
        save
        complex *16 u(n,n),z(1),roots(1),a(1,1)
        dimension errs(1)
c
c        This subroutine tests the "quality" of the putative
c        eigenvectors u(1,i) and roots(i) of the matrix a by
c        calculating the norms of the differences
c
c        A (u_i) - roots_i * u_i;
c
c        These are returned in array errs.
c
c                   Input parameters:
c
c  a - the matrix whose eigendecomposition we are supposed to
c        be testing
c  u - the putative eigenvectors of a
c  n - the dimensionality of matrices a, u
c  roots - the putative eigenvalues of a
c
c                   Output parameters:
c
c  errors - the errors in the eigenvector/eigenvalue pair 
c        (in this subroutine's opinion)
c
c                   Work array:
c
c  z - must be at least n complex *16 elements long
c        
c       
c
        do 3000 i=1,n
c
c       check how good this eigenvector/eigenvalue pair is
c
        d=-1
c
        call nonsym_matvec(a,n,u(1,i),z)
c
        d=0
        do 1400 j=1,n
c
        d=d+(z(j)-u(j,i)*roots(i))*conjg(z(j)-u(j,i)*roots(i))        
 1400 continue
c
        errs(i)=sqrt(d)
 3000 continue
c
        return
        end
c
c
c
c
c
        subroutine nonsym_eigvals(ier,a,n,roots,nroots,w)
        implicit real *8 (a-h,o-z) 
        complex *16 a(1),roots(1)
        real *8 w(1)
c
c       This subroutine finds the spectrum of a (generally, 
c       non-Hermitian) complex matrix a. It uses a standard 
c       reduction to the upper Hessenberg form, followed by 
c       the QR procedure.
c
c               Input parameters:
c
c  a - the n \times n matrix whose spectrum is to be found;
c       PLEASE NOTE THAT THIS MATRIX IS DESTROYED BY THIS
c       SUBROUTINE - UTTERLY!
c  n - the dimensionality of a
c  
c               Output parameters:
c
c  ier - error return code;
c    ier=0 means normal conclusion
c    ier=16 means that the iterations have failed to converge
c  roots - the (complex) eigenvalues of the matrix a
c  nroots - the number of eigenvalues actually calculated and 
c       and returned in array roots; if ier=0 then nroots=n
c
c               Work array:
c
c  w - must be at least 27*n+50 real *8 locations
c
c  
c
c        . . . reduce the matrix to the upper Hessenberg form
c
        ifv=0
        call nonsym_tohessen(a,n,v,ifv)
c
c       find the spectrum
c
        call nonsym_all_recurse(ier,a,n,roots,nroots,w)
c
        return
        end     
c
c
c
c
c
        subroutine nonsym_eigvects(ier,a,n,roots,u,dmin,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(1),u(n,n),w(1),roots(1)
c
c       This subroutine finds the eigenvalues and eigenvectors
c       of a (generally, non-Hermitian) complex matrix a. In 
c       order to find the spectrum of a, it uses a standard 
c       reduction to the upper Hessenberg form, followed by 
c       the QR procedure. After that, it finds the eigenvectors
c       via a straightforward inverse power scheme.
c
c               Input parameters:
c
c  a - the n \times n matrix whose spectrum is to be found;
c       PLEASE NOTE THAT THIS MATRIX IS DESTROYED BY THIS
c       SUBROUTINE - UTTERLY!
c  n - the dimensionality of a
c  
c               Output parameters:
c
c  ier - error return code;
c    ier=0 means normal conclusion
c    ier=16 means that the subroutine found thye spectrum 
c       successfully, but gagged on the eigenvectors - 
c       can happen to the best of us
c    ier=64 means that the subroutine failed to find the 
c       spectrum. This indicates either a user bug, or
c       a failure of the QR process.
c  roots - the (complex) eigenvalues of the matrix a
c  nroots - the number of eigenvalues actually calculated and 
c       and returned in array roots; if ier=0 then nroots=n
c  u - the i-th column of this matrix is the i-th eigenvector 
c       of A

c               Work array:
c
c  w - must be at least 5/2 * n**2+6*n+500 complex *16 elements
c  
c
c       . . . reduce the input matrix to upper Hessenberg form
c
        iv=1
        lv=n*n +10
c
        ifv=1
        call nonsym_tohessen(a,n,w(iv),ifv)
c
c        determine the machine zero
c
        call nonsym_mach_zero(zero_mach)
c
c       determine eps
c
        d=0
        do 1200 i=1,n*n
c
        d=d+a(i)*conjg(a(i))
 1200 continue
c
        d=sqrt(d)
        delta=zero_mach*d  *10
c
c        allocate memory
c
        ib=iv+lv
        lb=n*n+10
c
        iw=ib+lb
c
c       find the eigenvalues
c
        call nonsym_copy(a,w(ib),n**2*2)
c
        call nonsym_all_recurse(jer,w(ib),n,roots,nroots,w(iw) )
c
        if(jer .ne. 0) then
            ier=64
            return
        endif
c
c        find the distances between roots
c
        dmin=1.0d40
        do 1800 i=2,nroots
        do 1600 j=1,i-1
c
        d=(roots(i)-roots(j))*conjg(roots(i)-roots(j))
        if(d .lt. dmin) dmin=d
 1600 continue
 1800 continue
c
        dmin=sqrt(dmin)
c
c        allocate memory
c
        ib=iv+lv
        lb=n*n+10
c
        iww=ib+lb 
        lww=n**2/2+3*n+300
c
        iy=iww+lww
        ly=n+2
c
        ixold=iy+ly
        lxold=n+2
c
        iz=ixold+lxols
        lz=n+2
c
c       find the eigenvectors
c
        valmax=0
        icount=0
c
        do 2400 i=1,n
c
        icount=icount+1
        if(icount .eq. 20) then
            call prinf('in nonsym_eigvects, i=*',i,1)
            icount=0
        endif
c
        call nonsym_copy(a,w(ib),n**2*2)
c
        call nonsym_vector_find(jer,a,n,w(ib),roots(i),val,u(1,i),
     1      delta,w(iy),w(ixold),w(iww),w(iz) )
c
        if(jer .ne. 0) then
            ier=16
            return
        endif
c
        if(val .gt. valmax) valmax=val
 2400 continue
c
        call nonsym_matmul4(w(iv),u,n,w(ib))
c
        call nonsym_copy(w(ib),u,n**2*2)
c
        return
        end
c
c        
c
c
c
        subroutine nonsym_tohessen(a,n,v,ifv)
        implicit real *8 (a-h,o-z) 
        save
        complex *16 a(n,n),aa(3),u(2,2),us(2,2),v(n,n)
c
c       set be to identity
c
        if(ifv .eq. 0) goto 1300
c
        do 1200 i=1,n
        do 1100 j=1,n
c
        v(j,i)=0
 1100 continue
        v(i,i)=1
 1200 continue
c
 1300 continue
c
c       eliminate rows
c
        do 2000 i=1,n-2
c
        do 1400 j=n,i+2,-1
c
        aa(1)=a(j-1,i)
        aa(2)=a(j,i)
c
        call nonsym_rotfnd(aa,us(1,1))
        call nonsym_rows_rotate(us(1,1),a,n,j-1,j)
        if(ifv .ne. 0) call nonsym_rows_rotate(us(1,1),v,n,j-1,j)
c
c       apply the same operator to the columns
c
        u(1,1)=conjg(us(1,1))
        u(2,2)=conjg(us(2,2))
        u(1,2)=conjg(us(1,2))
        u(2,1)=conjg(us(2,1))
c
        call nonsym_cols_rotate(u,a,n,j-1,j)
c
 1400 continue
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
        subroutine nonsym_matmul4(a,b,n,c)
        implicit real *8 (a-h,o-z) 
        complex *16 a(n,n),cd,b(n,n),c(n,n)
c
        do 1600 i=1,n
        do 1400 j=1,n
        cd=0
        do 1200 k=1,n
c
        cd=cd+conjg(a(k,i))*b(k,j)
 1200 continue
        c(i,j)=cd
 1400 continue
 1600 continue
c
        return
        end
c
c
c
c
c
        subroutine nonsym_vector_find(ier,a,n,b,clam,discrep,x,
     1      delta,y,xold,ww,z)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,n),clam,b(n,n),
     1      x(1),y(1),clam2,z(1),ww(1),xold(1)
c
c       prepare the matrix to be inverted
c
        ier=0
        call nonsym_copy(a,b,n**2*2)
c
        eps=delta*100*100
c
        do 1200 i=1,n
c
        b(i,i)=b(i,i)-clam +delta
 1200 continue
c
c       invert b
c
        call nonsym_hessen_invert(b,n,ww,nww)
c
c       use inverse power method to calculate the eigenvector
c
        call nonsym_corrand(N*2,y)
c
        ifout=0
        do 2000 i=1,100
c
        call nonsym_copy(y,x,n*2) 
        call nonsym_hessen_inv_apply(n,ww,x)
c
        call nonsym_copy(x,y,n*2) 
        call nonsym_hessen_inv_apply(n,ww,y)
c
        d=0
        do 1400 j=1,n
        d=d+x(j)*conjg(x(j))
 1400 continue
c
        ddd=0
        d=1/sqrt(d)
        do 1600 j=1,n
c
        x(j)=x(j)*d
 1600 continue
c
        call nonsym_matvec(a,n,x,z)
        call nonsym_coef_fnd(x,z,n,err,clam2)
c
        call nonsym_copy(x,xold,n*2)
        err=(err+abs(clam2-clam))/( abs(clam2)+abs(clam) )
c
        call prin2('err=*',err,1)

        if(err .lt. eps) ifout=ifout+1
c
        call nonsym_copy(x,y,n*2) 
c
        if(ifout .eq. 2) goto 2200
 2000 continue
c
 2200 continue
c
        discrep=abs(clam2-clam)
c
        return
        end
c
c        
c
c
c
        subroutine nonsym_coef_fnd(x,y,n,err,cc)
        implicit real *8 (a-h,o-z) 
        complex *16 cd,x(1),y(1),cc
c
c        eliminate the subdiagonal
c
        cd=0
        d=0
        do 1200 i=1,n
c
        d=d+x(i)*conjg(x(i))
        cd=cd+y(i)*conjg(x(i))
 1200 continue
c
        cc=cd/d
c
        err=0
        do 1400 i=1,n
c
        cd=y(i)-cc*x(i)
        err=err+cd*conjg(cd)
 1400 continue
c
        err=sqrt(err)
        return
        end
c
c        
c
c
c
        subroutine nonsym_hessen_inv_apply(n,ww,y)
        implicit real *8 (a-h,o-z) 
        complex *16 cd,ww(1),y(1)
c
c        eliminate the subdiagonal
c
        ii=0
        do 2000 i=1,n
c
        ii=ii+1
        cd=ww(ii)
c
        y(i)=y(i)*cd
c
        if(i .eq. n) goto 2000
c
        ii=ii+1
        cd=ww(ii)
        j=i+1
c
        y(j)=y(j)-cd*y(i)
c
 2000 continue
c
c        eliminate the upper triangle
c
        do 3000 i=n,2,-1
c
        ii=ii+1
        cd=ww(ii)
c
        y(i)=y(i)*cd
c
        do 2800 j=1,i-1
c
        ii=ii+1
        cd=ww(ii)
c
        y(j)=y(j)-cd*y(i)
 2800 continue
 3000 continue
c
        return
        end
c
c        
c
c
c
        subroutine nonsym_hessen_invert(a,n,ww,ii)
        implicit real *8 (a-h,o-z) 
        complex *16 a(n,n),cd,ww(1)  
c
        ii=0
c
c        eliminate the subdiagonal
c
        do 2000 i=1,n
c
        cd=1/a(i,i)
c
        ii=ii+1
        ww(ii)=cd
c
        do 1400 j=1,n
c
        a(i,j)=a(i,j)*cd
c
 1400 continue
c
        if(i .eq. n) goto 2000
c
        j=i+1
c
        cd=a(j,i)
c
        ii=ii+1
        ww(ii)=cd
c
        do 1600 k=i,n
c
        a(j,k)=a(j,k)-cd*a(i,k)
 1600 continue
c
 2000 continue
c
c        eliminate the upper triangle
c
        do 3000 i=n,2,-1
c
        cd=1/a(i,i)
c
        ii=ii+1
        ww(ii)=cd
c
        do 2400 j=1,n
c
        a(i,j)=a(i,j)*cd
 2400 continue
c
        do 2800 j=1,i-1
c
        cd=a(j,i)
c
        ii=ii+1
        ww(ii)=cd
c
        a(j,i)=a(j,i)-cd*a(i,i)
 2800 continue
 3000 continue
c
        return
        end
c
c        
c
c
c
        subroutine nonsym_matvec(a,n,x,y)
        implicit real *8 (a-h,o-z) 
        complex *16 a(n,n),cd,x(1),y(1)
c
        do 1400 i=1,n
        cd=0
        do 1200 j=1,n
c
        cd=cd+a(i,j)*x(j)
 1200 continue
        y(i)=cd
 1400 continue
c
        return
        end
c
c        
c
c
c
        subroutine nonsym_all_recurse(ier,a,n,roots,nroots,w)
        implicit real *8 (a-h,o-z) 
        complex *16 a(1)
        dimension w(1)
c
c        determine the machine zero
c
        call nonsym_mach_zero(zero_mach)
c
c       determine eps
c
        d=0
        do 1200 i=1,n*n
c
        d=d+a(i)*conjg(a(i))
 1200 continue
c
        d=sqrt(d)
        eps=zero_mach*d
c
c        allocate memory
c
        iistack=1
        listack=n*4*4
c
        ius=iistack+listack
        lus=8*n+10
c
        ierrs=ius+lus
        lerrs=3*n+10    
c
c       find the spectrum
c
        call nonsym_all_recurse0(ier,a,n,eps,roots,nroots,
     1      w(iistack),w(ius),w(ierrs) )
c
        return
        end
c
c        
c
c
c
        subroutine nonsym_all_recurse0(ier,a,n,eps,roots,nroots,
     1      istack,us,errs)
        implicit real *8 (a-h,o-z) 
        integer *4 istack(4,1)
        complex *16 a(1),us(1)
        dimension errs(3,1)
c
c   istack(1,i)=i
c   istack(2,i) - tells whether this block has been processed
c   istack(3,i) - this block's location in the big array
c   istack(4,i) - this block's size
c
c
c
c        initialize the recursion
c
        ier=0
        istack(1,1)=1
        istack(2,1)=0
        istack(3,1)=1
        istack(4,1)=n
c
c        one after another, recursively subdivide the matrix
c
        nroots=0
        nblocks=1
        icount=0
        do 4000 iijjkk=1,1000 000
c
        icount=icount+1
        if(icount .eq. 20) then 
            call prinf('in nonsym_all_recurse0, iijjkk=*',iijjkk,1)
            icount=0
        endif

        do 3000 i=1,nblocks 
c
        if(istack(2,i) .ne. 0) goto 3000
c
c       this block has not been processed; act accordingly 
c
        iaddr=istack(3,i)
        m=istack(4,i)
c
        call nonsym_deal_withit(jer,a,m,eps,roots,nroots,istack,
     1      i,iaddr,nblocks,a(iaddr),us,errs)
        if(jer .ne. 0) then
            ier=16
            return
        endif
c
        if(nroots .eq. n) goto 4200
c
 3000 continue
 4000 continue
 4200 continue
c
        return
        end
c
c        
c
c
c
        subroutine nonsym_deal_withit(ier,aoned,n,eps,roots,nroots,
     1      istack,instack,iaddr,nblocks,a2,us,errs)
        implicit real *8 (a-h,o-z) 
        integer *4 istack(4,1)
        complex *16 aa(2,2),clam1,clam2,roots(1),aoned(1),
     1      a2(n,n),us(1)
        dimension errs(3,1)
c
c       if the dimensionality of a is 1 or 2, act accordingly
c
        ier=0  
        if(n .eq. 1) then
            nroots=nroots+1
            roots(nroots)=a2(1,1)
            istack(2,instack)=1
            return
        endif
c
        if(n .eq. 2) then
            call nonsym_quadr_solve(aoned(iaddr),clam1,clam2)
            nroots=nroots+1
            roots(nroots)=clam1
            nroots=nroots+1
            roots(nroots)=clam2
            istack(2,instack)=1
            return
        endif
c
c        write(*,*) eps
        call nonsym_eigen0(jer,aoned(iaddr),n,us,eps,niter,izero,errs)
c
        if(jer .ne. 0) then
            ier=16
            return
        endif
c
c        if neither of the chunks is smaller than 3, store
c        both of them as new blocks
c
        m=izero
c
        if( (izero .gt. n-3) .or. (izero .lt. 3) ) goto 1400
c
        m=izero
        iaddr1=iaddr
        iaddr2=iaddr1+m**2
        call nonsym_reform(aoned(iaddr),aoned(iaddr1),aoned(iaddr2),n,m)
c
        nblocks=nblocks+1
        nb=nblocks
c
        istack(1,nb)=nb
        istack(2,nb)=0
        istack(3,nb)=iaddr1
        istack(4,nb)=m
c
        nblocks=nblocks+1
        nb=nblocks
c
        istack(1,nb)=nb
        istack(2,nb)=0
        istack(3,nb)=iaddr2
        istack(4,nb)=n-m
c
        istack(2,instack)=1
c
        return
c
 1400 continue
c
c       if the last chunk has size 1, act accordingly
c
        if(izero .ne. n-1) goto 1600
c
        m=n-1
        iaddr1=iaddr
c
        nblocks=nblocks+1
        nb=nblocks
c
        istack(1,nb)=nb
        istack(2,nb)=0
        istack(3,nb)=iaddr1
        istack(4,nb)=m
c
        istack(2,instack)=1
c
        nroots=nroots+1
        roots(nroots)=a2(n,n)
c
        iaddr1=iaddr
        call nonsym_reform3(aoned(iaddr),aoned(iaddr1),n,m)
c
        return
c
 1600 continue
c
c        if the last chunk has size 2, act accordingly
c
        if(m .ne. n-2) goto 1800
c
        nblocks=nblocks+1
        nb=nblocks

        iaddr1=iaddr
c
        istack(1,nb)=nb
        istack(2,nb)=0
        istack(3,nb)=iaddr1
        istack(4,nb)=m
c
        istack(2,instack)=1
c
        aa(1,1)=a2(n-1,n-1)
        aa(2,2)=a2(n,n)
        aa(1,2)=a2(n-1,n)
        aa(2,1)=a2(n,n-1)
c
        call nonsym_quadr_solve(aa,clam1,clam2)
c
        nroots=nroots+1
        roots(nroots)=clam1
        nroots=nroots+1
        roots(nroots)=clam2
c
        iaddr1=iaddr
        iaddr2=iaddr1+m**2
        call nonsym_reform3(aoned(iaddr),aoned(iaddr1),n,m)
c
        return
c
 1800 continue
c
c       if the first chunk has size 1, act accordingly
c
        if(izero .ne. 1) goto 2200
c
        m=1
        iaddr1=iaddr
c
        nblocks=nblocks+1
        nb=nblocks
c
        istack(1,nb)=nb
        istack(2,nb)=0
        istack(3,nb)=iaddr1
        istack(4,nb)=n-m
c
        istack(2,instack)=1
c
        nroots=nroots+1
        roots(nroots)=a2(1,1)
c
        iaddr1=iaddr
        iaddr2=iaddr1+m**2
        call nonsym_reform2(aoned(iaddr),aoned(iaddr),n,m)
c       
        return
c
 2200 continue
c
c       if the first chunk has size 2, act accordingly
c
        if(izero .ne. 2) goto 2400
c
        iaddr1=iaddr
        nblocks=nblocks+1
        nb=nblocks
c
        istack(1,nb)=nb
        istack(2,nb)=0
        istack(3,nb)=iaddr1
        istack(4,nb)=n-m
c
        aa(1,1)=a2(1,1)
        aa(2,2)=a2(2,2)
        aa(1,2)=a2(1,2)
        aa(2,1)=a2(2,1)
c
        call nonsym_quadr_solve(aa,clam1,clam2)
c
        nroots=nroots+1
        roots(nroots)=clam1
        nroots=nroots+1
        roots(nroots)=clam2
c
        m=2
        iaddr2=iaddr1+m**2

        istack(2,instack)=1
c
        call nonsym_reform2(aoned(iaddr),aoned(iaddr1),n,m)
c
        return
c
 2400 continue
c
        return
        end
c
c
c
c
c
        subroutine nonsym_reform2(a,b,n,m)
        implicit real *8 (a-h,o-z) 
        complex *16 a(n,n),b(n-m,n-m)
c
        do 1400 i=1,n-m
        do 1200 j=1,n-m
c
        b(j,i)=a(j+m,i+m)
 1200 continue 
 1400 continue
c
        return
        end
c
c
c
c
c
        subroutine nonsym_reform3(a,b,n,m)
        implicit real *8 (a-h,o-z) 
        complex *16 a(n,n),b(m,m)
c
        do 1400 i=1,m
        do 1200 j=1,m
c
        b(j,i)=a(j,i)
 1200 continue 
 1400 continue
c
        return
        end
c
c
c
c
c
        subroutine nonsym_reform(a,b,c,n,m)
        implicit real *8 (a-h,o-z) 
        complex *16 a(n,n),b(m,m),c(n-m,n-m)
c
        do 1400 i=1,m
        do 1200 j=1,m
c
        b(j,i)=a(j,i)
 1200 continue 
 1400 continue
c
c
        do 2400 i=1,n-m
        do 2200 j=1,n-m
c
        c(j,i)=a(j+m,i+m)
 2200 continue 
 2400 continue
c
        return
        end
c
c
c
c
c
        subroutine nonsym_eigen0(ier,a,n,us,eps,niter,izero,errs)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,n),aa(2,2),us(2,2,1),clam1,clam2
        dimension errs(3,1)
c
        ier=0
        do 1200 i=1,n
c
        errs(1,i)=1.0d40
        errs(2,i)=1.0d40
        errs(3,i)=1.0d40
 1200 continue
c
        ier=0
        ifout=0
        do 4000 ijk=1,1000
c
        niter=ijk
c
c       determine the shifts
c
        aa(1,1)=a(n-1,n-1)
        aa(2,2)=a(n,n)
        aa(1,2)=a(n-1,n)
        aa(2,1)=a(n,n-1)
c
        call nonsym_quadr_solve(aa,clam1,clam2)
c
c       . . . shift and eliminate
c
        do 2200 i=1,n
c
        a(i,i)=a(i,i)-clam1
 2200 continue
c
        call nonsym_eigen00(a,n,us)
c
        do 2400 i=1,n
c
        a(i,i)=a(i,i)+clam1-clam2
 2400 continue
c
        call nonsym_eigen00(a,n,us)
c
        do 2600 i=1,n
c
        a(i,i)=a(i,i)+clam2
 2600 continue    
c
c       check if we have converged
c
        errmin=1d1*eps
        ifout=0
        do 2800 i=1,n-1
c
        d=abs(a(i+1,i))
        if(d .lt. errmin) errmin=d
        errs(1,i)=errs(2,i)
        errs(2,i)=errs(3,i)
        errs(3,i)=d

cccc        call prin2('errs(1,i)=*',errs(1,i),3)
c
        if( (errs(1,i) .lt. eps) .and. (errs(2,i) .lt. eps)
     1      .and. (errs(3,i) .lt. eps) ) then
c
            ifout=1 
        endif
 2800 continue

        if(ifout .eq. 1) goto 4200
c       
 4000 continue
c
        ier=16
 4200 continue
c
c        look at the subdiagonals more carefully
c
        d=1d1*eps
        do 4400 i=1,n-1
c
        if(errs(3,i) .lt. d) then
            d=errs(3,i)
            izero=i
        endif
 4400 continue
c
        return
        end
c
c
c
c
c
        subroutine nonsym_quadr_solve(aa,clam1,clam2)
        implicit real *8 (a-h,o-z)
        save
        complex *16 aa(2,2),clam1,clam2,discr
        data four/4.0d0/,half/0.5d0/
c
        discr=sqrt( (aa(1,1)+aa(2,2))**2+
     1      four*( aa(1,2)*aa(2,1) -aa(1,1)*aa(2,2)) )

        clam1=aa(1,1)+aa(2,2)+discr
     1      
        clam1=clam1*half
c
        clam2=aa(1,1)+aa(2,2)-discr
        clam2=clam2*half
c
        return
        end
c
c
c
c
c
        subroutine nonsym_eigen00(a,n,us)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,n),aa(3),u(2,2),us(2,2,1)
c
c       eliminate the subdiagonal
c
        do 1600 i=1,n-1
c
c        determine the complex rotation
c
        aa(1)=a(i,i)
        aa(2)=a(i+1,i)
c
        call nonsym_rotfnd(aa,us(1,1,i))
        call nonsym_rows_rotate(us(1,1,i),a,n,i,i+1)
c
 1600 continue
c
c       restore the subdiagonal
c
        do 1800 i=1,n-1
c
        u(1,1)=conjg(us(1,1,i))
        u(2,2)=conjg(us(2,2,i))
        u(1,2)=conjg(us(1,2,i))
        u(2,1)=conjg(us(2,1,i))
c
        call nonsym_cols_rotate(u,a,n,i,i+1)
c
 1800 continue
c
        return
        end
c
c
c
c
c
        subroutine nonsym_copy(a,b,n)
        implicit real *8 (a-h,o-z)
        dimension a(1),b(1)
c
        do 1200 i=1,n
        b(i)=a(i)
 1200 continue
c
        return
        end
c 
c 
c 
c 
c 
        subroutine nonsym_rotfnd(a,u)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(2),u(2,2)
c 
        u21=-a(2)
        u22=a(1)
c 
        d=sqrt(u22*conjg(u22)+u21*conjg(u21))
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
        u(1,1)=-conjg(u(2,2))
        u(1,2)=conjg(u(2,1))
        return
        end
c 
c 
c 
c 
c 
        subroutine nonsym_cols_rotate(u,a,n,ii,jj)
        implicit complex *16 (a-h,o-z)
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
        subroutine nonsym_rows_rotate(u,a,n,ii,jj)
        implicit complex *16 (a-h,o-z)
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
        subroutine nonsym_mach_zero(zero_mach)
        implicit real *8 (a-h,o-z)
        save
c
        zero_mach=100       
c
        d1=1.1
        d3=1.1
        d=1.11
        do 1200 i=1,1000
c

        d=d/2
        d2=d1+d
        call nonsym_mach_zero0(d2,d3,d4)
c
        if(d4 .eq. 0) goto 1400
c
 1200 continue
 1400 continue
c
        zero_mach=d
        return
        end
c 
c 
c 
c 
c 
        subroutine nonsym_mach_zero0(a,b,c)
        implicit real *8 (a-h,o-z)
        save
c
        c=b-a

        return
        end
c
c
c
c
c
        SUBROUTINE nonsym_corrand(N,Y)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION Y(1)
        DATA IFCALL /0/
c
c        This subroutine returns to the user a collection of
c        reasonably random numbers uniformly distributed on
c        the interval [0,1]
c
C       . . . CONDUCT PRELIMINARY RANDOMIZATION
C 
        IF (IFCALL.EQ.1 ) GOTO  1300
        IFCALL =1
c
c        construct parameters for the first process
c
        LAMBDA=2**7 +9
        MU=1
        IP=2**20
        M=17
c
c        construct parameters for the second process
c
        LAMBDAq=2**5 +29
        MUq=3
        IPq=2**22
        Mq=19
c
c        construct parameters for the third process
c
        LAMBDAqq=2**3 +31
        MUqq=13
        IPqq=2**18
        Mqq=31
c
        DO 1200 I=1,100
        M1=M*LAMBDA  +MU
        J=M1/IP
        M=M1-J*IP
c
        M1=Mq*LAMBDAq  +MUq
        J=M1/IPq
        Mq=M1-J*IPq
c
        M1=Mqq*LAMBDAqq  +MUqq
        J=M1/IPqq
        Mqq=M1-J*IPqq
c
 1200 CONTINUE
c
 1300 CONTINUE
c
        D=1
        D=D/IP
        Dq=1
        Dq=Dq/IPq
        Dqq=1
        Dqq=Dqq/IPqq
C 
C       GENERATE PSEUDO-RANDOM NUMBERS
C 
        DO 1400 I=1,N
c
        M1=M*LAMBDA +MU
        J=M1/IP
        M=M1-J*IP
c
        M1=Mq*LAMBDAq +MUq
        J=M1/IPq
        Mq=M1-J*IPq
c
        M1=Mqq*LAMBDAqq +MUqq
        J=M1/IPqq
        Mqq=M1-J*IPqq
c
        Y(I)=M*D+Mq*Dq+mqq*dqq
c
        call corrand_comp(y(i),rint)
c
        y(i)=rint
 1400 CONTINUE
c
        RETURN
        END
