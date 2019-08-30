c
c       This is the end of the debugging code, and the beginning of the 
c       code for the diagonalization of a tridiagonal symmetric matrix,
c       and for the SVD of a bidiagonal real matrix
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains three user-callable subroutines: qleigen_trid
c       qleigen_bidiag_vects, and qleigen_bidiag_vals. Below is a brief 
c       description of these three subroutines.
c
c   qleigen_bidiag_vects - For a user-supplied bidiagonal real matrix, 
c       constructs the said matrice's Singular Value Decomposition
c       (SVD). It converts the input bidiagonal matrix into a 
c       tridiagonal matrix of twice the size, and diagonalizes the
c       latter via a (fairly standard) Q-L algorithm followed by
c       a fairly standard inverse power method.
c
c   qleigen_bidiag_vals - For a user-supplied bidiagonal real matrix, 
c       this subroutine constructs the said matrix's Singular Values. 
c       It converts the input bidiagonal matrix into a tridiagonal 
c       matrix of twice the size, and finds the eigenvalues of the 
c       latter via a (fairly standard) Q-L algorithm. The cost is of 
c       the order N^2 operations, with a moderate coefficient.
c
c   qleigen_trid - uses a version of Q-L algorithm to find
c        the spectrum of a user-supplied tridiagonal symmetric
c        (real) matrix. After that (if requested by the user), it 
c        uses an inverse power scheme to find the corresponding 
c        eigenvectors. The cost is of the order N^2 operations, 
c        with a smallish coefficient. 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine qleigen_bidiag_vals(ier,nn,da,suba,
     1      rlams,w)
        implicit real *8 (a-h,o-z)
        real *8 da(1),suba(1),rlams(1),w(1)
c
c       For a user-supplied bidiagonal real matrix, this subroutine 
c       constructs the said matrix's Singular Values. It converts the 
c       input bidiagonal matrix into a tridiagonal matrix of twice 
c       the size, and finds the eigenvalues of the latter via a 
c       (fairly standard) Q-L algorithm. The cost is of the order 
c       N^2 operations, with a moderate coefficient.
c
c                Input parameters:
c
c  nn - the dimensionality of the matrix
c  da - the diagonal of the matrix to be SWD'ed
c  suba - the subdiagonal of the matrix to be SWD'ed
c
c                Output parameters:
c
c  ier - error return code:
c      ier=0 means normal return
c      ier=16 means that the scheme failed to converge at some point.
c       This is a whorrible error; it has never been observed, and is 
c       not expected. If it has happened, it is most likely a result 
c       of a bug. The author encourages the user to assume that it is 
c       HIS bug, as opposed to the author's. An exceeded array would
c       be the most likely explanation.
c        
c  rlams - the n singular values of the matrix A
c
c                 Work array:
c
c  w - must be at least 18*nn+100 real *8 elements long
c
c
        n=nn*2
c
        isub7=1
        lsub7=n+2
c
        idiag7=isub7+lsub7
        ldiag7=n+2
c
        irlams7=idiag7+ldiag7
        lrlams7=n+2
c
        iw7=irlams7+lrlams7
        lw7=6*n+100
c
        ltot=iw7+lw7
c        
        call qleigen_bidiag_vals0(ier,nn,da,suba,rlams,
     1      w(isub7),w(idiag7),w(irlams7),w(iw7))
c
        return
        end
c
c
c
c
c
        subroutine qleigen_bidiag_vals0(ier,nn,da,suba,
     1      rlams,sub7,diag7,rlams7,w7)
        implicit real *8 (a-h,o-z)
        save
        real *8 da(1),suba(1),rlams(1),
     1      rlams7(1),w7(1),diag7(1),sub7(1)
c
c
c       . . . create the tridiagonal matrix
c
        ier=0
        n=nn*2
c
        ii=0
        do 1200 i=1,nn-1
c
        ii=ii+1
        sub7(ii)=da(i)
        ii=ii+1
        sub7(ii)=suba(i)
 1200 continue
c
        sub7(2*nn-1)=da(nn)
c
        do 1250 i=1,nn*2
c
        diag7(i)=0
 1250 continue
c
        nvects=0
c
        call qleigen_trid_svd(ier,nn*2,diag7,sub7,rlams7,
     1      nvects,w7,u,v)
c
        do 1400 i=1,nn
c
        rlams(i)=rlams7(i)
 1400 continue
c
        return
        end
c
c
c
c
c
        subroutine qleigen_bidiag_vects(ier,nn,da,suba,
     1      rlams,u,v,w)
        implicit real *8 (a-h,o-z)
        save
        real *8 da(1),suba(1),rlams(1),u(nn,1),v(nn,1),w(1)
c
c       For a user-supplied bidiagonal real matrix, this subroutine 
c       constructs the said matrix's Singular Value Decomposition
c       (SVD). It converts the input bidiagonal matrix into a 
c       tridiagonal matrix of twice the size, and diagonalizes the
c       latter via a (fairly standard) Q-L algorithm followed by
c       a fairly standard inverse power method.In other words, 
c       given an input matrix A, the subroutine produces the 
c       matrices U, D, V, such that
c
c               A = U D V^*;
c
c       Needless to say, only the diagonal elements of D 
c       are returned. The subroutine uses a simple and 
c       effective version of the inverse power method. In 
c       most cases, the cost is of the order N^2 operations, 
c       with a largish coefficient.
c
c                Input parameters:
c
c  nn - the dimensionality of the matrix
c  da - the diagonal of the matrix to be SWD'ed
c  suba - the subdiagonal of the matrix to be SWD'ed
c
c                Output parameters:
c
c  ier - error return code:
c      ier=0 means normal return
c      ier=16 means that the scheme failed to converge at some point.
c       This is a whorrible error; it has never been observed, and is 
c       not expected. If it has happened, it is most likely a result 
c       of a bug. The author encourages the user to assume that it is 
c       HIS bug, as opposed to the author's. An exceeded array would
c       be the most likely explanation.
c        
c  rlams - the n singular values of the matrix A
c  u -left singular vectors of the input matrix A
c  v - right singular vectors of the matrix A
c
c                 Work array:
c
c  w - must be at least 18*nn+100 real *8 elements long
c
        n=nn*2
c
        isub7=1
        lsub7=n+2
c
        idiag7=isub7+lsub7
        ldiag7=n+2
c
        irlams7=idiag7+ldiag7
        lrlams7=n+2
c
        iw7=irlams7+lrlams7
        lw7=6*n+100
c
        ltot=iw7+lw7
c        
        call qleigen_bidiag_vects0(ier,nn,da,suba,rlams,
     1      u,v,w(isub7),w(idiag7),w(irlams7),w(iw7))
c
        return
        end
c
c
c
c
c
        subroutine qleigen_trid(ier,n,diag,sub,rlams,vects,nvects,w)
        implicit real *8(a-h,o-z)
        save
        dimension sub(1),diag(1),rlams(1),w(1),vects(1)
c
c
c        This subroutine uses a version of Q-L algorithm to find
c        the spectrum of a user-supplied tridiagonal symmetric
c        (real) matrix. After that (if requested by the user), it 
c        uses an inverse power scheme to find the corresponding 
c        eigenvectors. The cost is of the order N^2 operations, 
c        with a smallish coefficient. 
c
c                Input parameters:
c
c  n - the dimensionality of the matrix
c  diag - the diagonal of the matrix to be diagonalized (n elements)
c  sub - the subdiagonal of the matrix to be diagonalized (n-1 elements)
c  nvects - integer parameter telling the subroutine how many 
c        eigenvectors are to be constructed (it always constructs 
c        all n eigenvalues)
c    nvects=k > 0 tells the subroutine to construct k eigenvectors
c    nvects=0 tells the subroutine to only calculate the spectrum
c
c                Output parameters:
c
c  ier - error return code:
c      ier=0 means normal return
c      ier=1000000 means that the Q-L scheme failed to converge 
c        after recovering i eigenvalues. This is a whorrible error; 
c        it has never been observed, and is not expected. If it has 
c        happened, it is most likely a result of a bug. The author 
c        encourages the user to assume that it is HIS bug, as 
c        opposed to the author's. An exceeded array would be the 
c        most likely explanation.
c        
c  rlams - the n eigenvalues of the matrix a
c  vects - the eigenvectors corresponding to the eigenvalues rlams
c                 Work array:
c
c  w - must be at least 6*n+100 real *8 elements long 
c
c
c       . . . allocate memory
c
        ier=0
c
        idiag=1
        ldiag=n+2
c
        isub1=idiag+ldiag
        lsub1=n+2
c
        isuper1=isub1+lsub1
        lsuper1=n+2
c
        isub2=isuper1+lsuper1
        lsub2=n+2
c
        iuu=isub2+lsub2
        luu=n*2+4
c
        ltot=iuu+luu
cccc        call prinf('and ltot=*',ltot,1)
c
        call qleigen_copy(diag,w(idiag),n)
        call qleigen_copy(sub,w(isub1),n-1)
        call qleigen_copy(sub,w(isuper1),n-1)
c
        call qleigen_trid0(ier,n,w(isub1),w(isub2),w(idiag),
     1      w(isuper1),w(iuu),rlams)
c  
        if(ier .ne. 0) return
c
        call qleigen_bubble(rlams,n)
c
cccc        call prin2('after bubbling, rlams=*',rlams,n)

cccc          stop
        if(nvects .eq. 0) return
c
c       construct the eigenvectors
c
        idiag1=1
        ldiag1=n+2
c
        isub1=idiag1+ldiag1
        lsub1=n+2
c
        irhs=isub1+lsub1
        lrhs=n+2
c       
        iww=irhs+lrhs
        lww=6*n+30 
c
        inums=iww+lww
        lnums=n+2
c
        call qleigen_vectors(diag,sub,rlams,n,vects,nvects,
     1      w(idiag1),w(isub1),w(irhs),w(iww),w(inums) )
c
        return
        end
c
c
c
c
c
        subroutine qleigen_bidiag_vects0(ier,nn,da,suba,
     1      rlams,u,v,sub7,diag7,rlams7,w7)
        implicit real *8 (a-h,o-z)
        save
        real *8 da(1),suba(1),rlams(1),u(nn,1),v(nn,1),
     1      rlams7(1),w7(1),diag7(1),sub7(1)
c
c
c       . . . create the tridiagonal matrix
c
        ier=0
        n=nn*2
c
        ii=0
        do 1200 i=1,nn-1
c
        ii=ii+1
        sub7(ii)=da(i)
        ii=ii+1
        sub7(ii)=suba(i)
 1200 continue
c
        sub7(2*nn-1)=da(nn)
c
        do 1250 i=1,nn*2
c
        diag7(i)=0
 1250 continue
c
        nvects=nn
c
        call qleigen_trid_svd(ier,nn*2,diag7,sub7,rlams7,
     1      nvects,w7,u,v)
c
        done=1
        sqrt2=sqrt(2*done)
        do 1400 i=1,nn
c
        jj=0 
        do 1300 j=1,nn
c
        u(j,i)=u(j,i)*sqrt2
        v(j,i)=v(j,i)*sqrt2
 1300 continue
c
        rlams(i)=rlams7(i)
 1400 continue
c
        return
        end
c
c
c
c
c
        subroutine qleigen_trid_svd(ier,n,diag,sub,rlams,nvects,w,
     1      u,v)
        implicit real *8(a-h,o-z)
        save
        dimension sub(1),diag(1),rlams(1),w(1),u(1),v(1)
c
c       . . . allocate memory
c
        ier=0
c
        idiag=1
        ldiag=n+2
c
        isub1=idiag+ldiag
        lsub1=n+2
c
        isuper1=isub1+lsub1
        lsuper1=n+2
c
        isub2=isuper1+lsuper1
        lsub2=n+2
c
        iuu=isub2+lsub2
        luu=n*2+4
c
        ltot=iuu+luu
        call prinf('and ltot=*',ltot,1)
c
        call qleigen_copy(diag,w(idiag),n)
        call qleigen_copy(sub,w(isub1),n-1)
        call qleigen_copy(sub,w(isuper1),n-1)
c
        call qleigen_trid0(ier,n,w(isub1),w(isub2),w(idiag),
     1      w(isuper1),w(iuu),rlams)
c  
        if(ier .ne. 0) return
c
        call qleigen_bubble(rlams,n)
c
cccc        call prin2('after bubbling, rlams=*',rlams,n)

cccc          stop
        if(nvects .eq. 0) return
c
c       construct the eigenvectors
c
        idiag1=1
        ldiag1=n+2
c
        isub1=idiag1+ldiag1
        lsub1=n+2
c
        irhs=isub1+lsub1
        lrhs=n+2
c       
        iww=irhs+lrhs
        lww=6*n+30 
c
        inums=iww+lww
        lnums=n+2
c
        call qleigen_vectors_svd(diag,sub,rlams,n,nvects,
     1      w(idiag1),w(isub1),w(irhs),w(iww),w(inums),u,v)
c
        return
        end
c
c 
c 
c 
c 
        subroutine qleigen_vectors_svd(diag,sub,rlams,n,nvects,
     1      diag1,sub1,rhs,ww,nums,u,v)
        implicit real *8 (a-h,o-z)
        save
        real *8 diag(1),sub(1),rhs(1),ww(1),
     1      diag1(1),rlams(1),sub1(1),u(n/2,1),v(n/2,1)
        integer *4 nums(1)
c
c       for each eigenvector, find the list of preceding ones to
c       which it will have to be orthogonalized
c
        thresh=rlams(1)/n 
c
        nums(1)=0
        do 1600 i=2,nvects
c
        ii=0
        num0=nums(i-1)
        do 1400 j=i-1,1,-1
c
        d=rlams(j)-rlams(i)
        if(d .gt. thresh) goto 1500
c
        ii=ii+1        
 1400 continue
c
 1500 continue
c
        nums(i)=ii
 1600 continue
c
        call prinf('nums=*',nums,nvects)
cccc        call prin2('thresh=*',thresh,1)
c
c       one eigenvector after another, find them things
c
        call mach_zero(delta)
        delta=delta*100
c
        do 5000 ijk=1,nvects
c
c       . . . shift
c
        call qleigen_copy(diag,diag1,n)
        call qleigen_copy(sub,sub1,n-1)
c
        d=rlams(ijk)
        do 2400 i=1,n
c
        diag1(i)=diag(i)-d
 2400 continue
c
c       . . .  construct the shifted inverse 
c
        call qleigen_fact(diag1,sub,n,ww,delta)
c
        call qleigen_rand(N,rhs)
        num=nums(ijk)
        if(num .ne. 0) call qleigen_orthog_em_svd(ijk,rhs,n,num,u,v)
c
        do 3000 ij=1,2
c
c       . . . apply the inverse
c
        call qleigen_comp_solv(rhs,n,ww)
c
c       if needed, orthogonalize the current (putative) 
c       eigenvector to the preceding nums(ijk) eigenvectors
c
        if(num .ne. 0) call qleigen_orthog_em_svd(ijk,rhs,n,num,u,v)
c
c       . . . normalize
c
        call qleigen_rscap(rhs,rhs,n,d)
c
        d=1/sqrt(d)
        do 2600 i=1,n
c
        rhs(i)=rhs(i)*d
 2600 continue
c
 3000 continue
c
        done=1
        sqrt2=sqrt(2*done)
        nn=n/2
c
        jj=0 
        do 4300 j=1,nn
c
        jj=jj+1
        u(j,ijk)=rhs(jj)
        jj=jj+1
        v(j,ijk)=rhs(jj)
 4300 continue
c

        call qleigen_rscap(u(1,ijk),u(1,ijk),nn,d1)
        call qleigen_rscap(v(1,ijk),v(1,ijk),nn,d2)

        d1=1/sqrt(d1) /sqrt2
        d2=1/sqrt(d2) /sqrt2
        do 4600 j=1,nn
c
        u(j,ijk)=u(j,ijk)*d1
        v(j,ijk)=v(j,ijk)*d2 
 4600 continue
c
 5000 continue
c
        return
        end
c
c
c
c
c
        subroutine qleigen_orthog_em_svd(ijk,rhs,n,num,u,v)
        implicit real *8(a-h,o-z)
        real *8 rhs(1),u(n/2,1),v(n/2,1)
c
        do 2000 j=1,num
c
        d=0
        d1=0
        d2=0

        ijj=ijk-j
        ii=0
c
        do 1200 jj=1,n/2
c
        ii=ii+1
        d1=d1+rhs(ii)*u(jj,ijj)
        ii=ii+1
        d2=d2+rhs(ii)*v(jj,ijj)
 1200 continue

        d1=d1*2
        d2=d2*2
c
        ii=0
        do 1600 jj=1,n/2
c
        ii=ii+1
        rhs(ii)=rhs(ii)-d1*u(jj,ijj)
        ii=ii+1
        rhs(ii)=rhs(ii)-d2*v(jj,ijj)
 1600 continue
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
        subroutine qleigen_vectors(diag,sub,rlams,n,vects,nvects,
     1      diag1,sub1,rhs,ww,nums)
        implicit real *8 (a-h,o-z)
        save
        real *8 diag(1),sub(1),vects(n,n),rhs(1),ww(1),
     1      diag1(1),rlams(1),sub1(1)
        integer *4 nums(1)
c
c       for each eigenvector, find the list of preceding ones to
c       which it will have to be orthogonalized
c
        thresh=rlams(1)/n 
c
        nums(1)=0
        do 1600 i=2,nvects
c
        ii=0
        num0=nums(i-1)
        do 1400 j=i-1,1,-1
c
        d=rlams(j)-rlams(i)
        if(d .gt. thresh) goto 1500
c
        ii=ii+1        
 1400 continue
c
 1500 continue
c
        nums(i)=ii
 1600 continue
c
cccc        call prinf('nums=*',nums,n)
cccc        call prin2('thresh=*',thresh,1)
c
c       one eigenvector after another, find them things
c
        call mach_zero(delta)
        delta=delta*100
c
        do 4000 ijk=1,nvects
c
c       . . . shift
c
        call qleigen_copy(diag,diag1,n)
        call qleigen_copy(sub,sub1,n-1)
c
        d=rlams(ijk)
        do 2400 i=1,n
c
        diag1(i)=diag(i)-d
 2400 continue
c
c       . . .  construct the shifted inverse 
c
        call qleigen_fact(diag1,sub,n,ww,delta)
c
        call qleigen_rand(N,rhs)
        num=nums(ijk)
        if(num .ne. 0) call qleigen_orthog_em(ijk,rhs,vects,n,num)
c
        do 3000 ij=1,2
c
c       . . . apply the inverse
c
        call qleigen_comp_solv(rhs,n,ww)
c
c       if needed, orthogonalize the current (putative) 
c       eigenvector to the preceding nums(ijk) eigenvectors
c
        if(num .ne. 0) call qleigen_orthog_em(ijk,rhs,vects,n,num)
c
c       . . . normalize
c
        call qleigen_rscap(rhs,rhs,n,d)
c
        d=1/sqrt(d)
        do 2600 i=1,n
c
        rhs(i)=rhs(i)*d
 2600 continue
c
 3000 continue

        do 3700 i=1,n
c
        vects(i,ijk)=rhs(i)
 3700 continue
c
 4000 continue

cccc        call prin2('and vects=*',vects,n*n)
c
        return
        end
c
c
c
c
c
        subroutine qleigen_orthog_em(ijk,rhs,vects,n,num)
        implicit real *8(a-h,o-z)
        real *8 rhs(1),vects(n,1)
c
        do 2600 j=1,num
c
        call qleigen_rscap(rhs,vects(1,ijk-j),n,d)
c
        do 2400 jj=1,n
c
        rhs(jj)=rhs(jj)-d*vects(jj,ijk-j)
 2400 continue
 2600 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine qleigen_rscap(x,y,n,cd)
        implicit real *8 (a-h,o-z)
        save
        real *8 x(1),y(1),cd
c
        cd=0
        do 1200 i=1,n
        cd=cd+x(i)*y(i)
 1200 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine qleigen_fact(da,suba,n,ww,delta)
        implicit real *8 (a-h,o-z)
        save
        real *8 ww(1),da(1),rhs(1)
        real *8 suba(1)
c 
c        allocate memory
c
        i0=0
c
        iu=1
        lu=n+2
c
        iv=iu+lu
        lv=n+2
c
        iw=iv+lv
        lw=n+2
c
        ia=iw+lw
        la=n+2
c
        call qleigen_copy(da,ww(ia),n*2)
c
        call qleigen_fact0(ww(ia),suba,suba(i0),n,
     1      ww(iu),ww(iv),ww(iw),delta)
c
        return
c
c
c
c
        entry qleigen_comp_solv(rhs,n,ww)
c
        iu=1
        lu=n+2
c
        iv=iu+lu
        lv=n+2
c
        iw=iv+lv
        lw=n+2
c
        call qleigen_comp_solv0(ww(iu),ww(iv),ww(iw),n,rhs)
        return
        end
c 
c 
c 
c 
c 
        subroutine qleigen_fact0(a,b,c,n,u,v,w,delta)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(1),u(1),v(1),w(1),rhs(1),b(1),c(1)
c 
c        eliminate down
c 
        done=1
        w(1)=done/(a(1)+sign(delta,a(1)) )
        do 1200 i=1,n-1
        d=c(i+1)*w(i)
        a(i+1)=a(i+1)-b(i)*d
c
        a(i+1)=a(i+1)+sign(delta,a(i+1))
        u(i)=d
        w(i+1)=done/a(i+1)
 1200 continue
c 
c        eliminate up
c 
        do 1400 i=n,2,-1
        v(i)=b(i-1)*w(i)
 1400 continue
c 
        return
c 
c 
c 
c 
        entry qleigen_comp_solv0(u,v,w,n,rhs)
c 
c        eliminate down
c 
        do 2400 i=1,n-1
        rhs(i+1)=rhs(i+1)-u(i)*rhs(i)
 2400 continue
c 
c        eliminate up
c 
        do 2600 i=n,2,-1
        rhs(i-1)=rhs(i-1)-rhs(i)*v(i) 
 2600 continue
c
c       scale
c 
        do 2800 i=1,n
        rhs(i)=rhs(i)*w(i)
 2800 continue
        return
        end
c
c
c
c
c
        SUBROUTINE qleigen_rand(N,Y)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION Y(1)
        DATA IFCALL /0/
c
        IF (IFCALL. EQ. 1) GOTO  1200
C
C       GENERATE PSEUDO-RANDOM NUMBERS
C
        lambda=13
        mu=1
        ip=2**20
        d=1.0d0/ip
        half=0.5d0
c
        ifcall=1
 1200 continue
c
        DO 1400 I=1,N
        M1=M*LAMBDA +MU
        J=M1/IP
        M=M1-J*IP
        Y(I)=M*D-half
 1400 CONTINUE
        RETURN
        end
c
c
c
c
c
        subroutine qleigen_trid0(ier,n,sub1,sub2,diag,super1,uu,rlams)
        implicit real *8(a-h,o-z)
        save
        dimension uu(2,2,1),sub1(1),sub2(1),diag(1),super1(1),
     1      rlams(1)
c
        data ijkmax/0/
c
c       find the machine epsilon
c
        ier=0
        d1=0
        d2=0
        do 1200 i=1,n-1
c
        sub2(i)=0
        d1=d1+sub1(i)**2
        d2=d2+diag(i)**2
 1200 continue
c
        d2=d2+diag(n)**2
c
        call mach_zero(eps)
c
        eps=eps*sqrt((2*d1+d2)/n)/10
c
c       find the first n-2 eigenvalues
c
        do 2000 i=1,n-2
c
cccc        call prinf('i=*',i,1)
c
        call qleigen_one_size(jer,n-i+1,sub1(i),sub2(i),
     1      diag(i),super1(i),uu,eps,rlamout,rlam,i,iijk)
c
cccc        if(iijk .eq. 10) stop

        if(jer .ne. 0) then
            ier=1000000 + i-1
            return
        endif
c
        if(ijkmax .lt. iijk) ijkmax=iijk
c
        rlams(i)=rlam
 2000 continue
c  
c       calculate the last two eigenvalues
c
        d1=( diag(n-1)+diag(n) + 
     1      sqrt ( (diag(n-1)-diag(n))**2+4*sub1(n-1)**2 ) )/2
c
        d2=( diag(n-1)+diag(n) -
     1      sqrt ( (diag(n-1)-diag(n))**2+4*sub1(n-1)**2 ) )/2
c
        rlams(n-1)=d1+rlamout
        rlams(n)=d2+rlamout
c
        call prinf('exiting qleigen_trid, ijkmax=*',ijkmax,1)
        ijkmax=0
c
        return
        end
c
c
c
c
c
        subroutine qleigen_one_size(ier,n,sub1,sub2,diag,super1,
     1      uu,eps,rlamout,rlam,iii,iijk)
        implicit real *8(a-h,o-z)
        save
        dimension uu(2,2,1),sub1(1),sub2(1),diag(1),super1(1)
c
        ifout=0
        ier=0
        rlamout=0
        if(iii .ne. 1) rlamout=rlam_old
        if(abs(sub1(1)) .lt. eps*100) ifout=ifout+1
        if(abs(sub1(1)) .lt. eps/1000) then
            rlam=rlamout+diag(1)
            return
        endif
c
        iijk=0
        do 2000 ijk=1,200
c
        iijk=ijk
c
        d1=( diag(1)+diag(2) + 
     1      sqrt ( (diag(1)-diag(2))**2+4*sub1(1)**2 ) )/2
c
        d2=( diag(1)+diag(2) - 
     1      sqrt ( (diag(1)-diag(2))**2+4*sub1(1)**2 ) )/2
c
        d=d1
        if(abs(d1-diag(1)) .gt. abs(d2-diag(1)) ) d=d2
c
cccc        if(abs(d-diag(1)) .lt. sqrt(eps) ) d=diag(1)
c

        if(abs(d) .lt. abs(diag(1)+diag(2))*1.0d-6) d=d1

        do 1400 i=1,n
c
        diag(i)=diag(i)-d
 1400 continue
c
        rlamout=rlamout+d
c
        call qleigen_one_oper(n,sub1,sub2,diag,super1,uu)
c
cccc        call prin2('eps=*',eps,1)

        if(abs(sub1(1)) .lt. eps*100) ifout=ifout+1
ccc        if(abs(sub1(1)) .lt. eps/n/1000 000/1000) goto 2200
ccccccc        if(abs(sub1(1)) .lt. eps/100) goto 2200
cccc        if(abs(sub1(1)) .le. 1.0d-60) goto 2200
        if(abs(sub1(1)) .le. eps**2) goto 2200
cccccc        if(ifout .eq. 6) goto 2200
 2000 continue
c
        ier=64
 2200 continue
c
        rlam_old=rlamout
        rlam=rlamout+diag(1)
c
cccc        call prinf('iijk=*',iijk,1)
        return
        end
c
c
c
c
c
        subroutine qleigen_one_oper(n,sub1,sub2,diag,super1,uu)
        implicit real *8(a-h,o-z)
        save
        dimension uu(2,1),sub1(1),sub2(1),diag(1),super1(1)
c
c       starting at the bottom, eliminate the superdiagonal
c
        ii=0
        do 2400 i=n,2,-1
c
        ii=ii+1
        alpha=super1(i-1)
        beta=diag(i)
c
c       construct an eliminating 2 \times 2 - matrix
c
        sinphi=-alpha/sqrt(alpha**2+beta**2)
        cosphi=beta/sqrt(alpha**2+beta**2)
c
        u11=cosphi
        u12=sinphi
        u21=-sinphi
        u22=cosphi
c
        uu(1,ii)=u11
        uu(2,ii)=u21
c
        diag(i)=u21*super1(i-1)+u22*diag(i)
c
        d=u11*diag(i-1)+u12*sub1(i-1)
        dd=u21*diag(i-1)+u22*sub1(i-1)
        diag(i-1)=d
        sub1(i-1)=dd       
c
        if(i .ne. 2) sub1(i-2)=u11*sub1(i-2)+u12*sub2(i-2)
 2400 continue
c
        ii=0
        do 3400 i=n,2,-1
c
        ii=ii+1
c
        u11=uu(1,ii)
        u12=-uu(2,ii)
        u21=-u12
        u22=u11
c
        dd=u21*sub1(i-1)+u22*diag(i)
        diag(i)=dd
c
        d=u11*diag(i-1)
        dd=u21*diag(i-1)
        diag(i-1)=d
        super1(i-1)=dd
 3400 continue
c
        call qleigen_copy(super1,sub1,n-1)
c
        return
        end
c
c
c
c
c
        subroutine qleigen_setzero(a,n)
        implicit real *8 (a-h,o-z)
        dimension a(1)
c
        do 1200 i=1,n
        a(i)=0
 1200 continue
c
        return
        end
c
c
c
c
c
        subroutine qleigen_copy(a,b,n)
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
        subroutine qleigen_bubble(roots,n)
        implicit real *8 (a-h,o-z)
        save
        real *8  roots(1),cd
c
        do 1400 i=1,n
c
        nj=n-i
        if(nj .lt. 3) nj=3
        do 1200 j=1,nj
c
        d1=roots(j)
        d2=roots(j+1)
        if(d1 .gt. d2) goto 1200
c
        cd=roots(j)
        roots(j)=roots(j+1)
        roots(j+1)=cd
c
 1200 continue
 1400 continue
c
        return
        end
