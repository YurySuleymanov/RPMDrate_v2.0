      double precision FUNCTION f1dim(x)
      INTEGER NMAX
      double precision x
      PARAMETER (NMAX=50)
CU    USES func
      INTEGER j,ncom
      double precision pcom(nmax),xicom(nmax),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
c        xt(j) = dabs(xt(j))
c        write(6,*) j,xt(j)
11    continue
      f1dim=func(xt)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software <%8.
