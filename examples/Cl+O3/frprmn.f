      SUBROUTINE frprmn(p,n,ftol,iter,fret,hd)
      implicit double precision (a-h,o-z)
      INTEGER iter,n,NMAX,ITMAX
      dimension p(n)
      external func
c      PARAMETER (NMAX=50,ITMAX=300,EPS=1.d-10)
      PARAMETER (NMAX=50,ITMAX=400,EPS=1.d-09)
CU    USES dfunc,func,linmin
      INTEGER its,j
      dimension g(NMAX),h(NMAX),xi(NMAX)
      fp=func(p)
      call dfunc(p,xi,n,hd)
      do 11 j=1,n
        g(j)=-xi(j)
        h(j)=g(j)
        xi(j)=h(j)
11    continue
      do 14 its=1,ITMAX
        iter=its
c        write(6,*) 'p',p(1),p(2),p(3)
        call linmin(p,xi,n,fret)
        if(2.*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+EPS))return
        fp=func(p)
        call dfunc(p,xi,n,hd)
        gg=0.
        dgg=0.
        do 12 j=1,n
          gg=gg+g(j)**2
C         dgg=dgg+xi(j)**2
          dgg=dgg+(xi(j)+g(j))*xi(j)
12      continue
        if(gg.eq.0.)return
        gam=dgg/gg
        do 13 j=1,n
          g(j)=-xi(j)
          h(j)=g(j)+gam*h(j)
          xi(j)=h(j)
13      continue
14    continue
c      pause 'frprmn maximum iterations exceeded'
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software <%8

