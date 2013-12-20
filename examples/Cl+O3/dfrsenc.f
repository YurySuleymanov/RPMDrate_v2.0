      double precision FUNCTION dfrsenc(ir,r,n,h,err)
      INTEGER NTAB,ir
      double precision err,h,CON,CON2,BIG,SAFE,func
      PARAMETER (CON=1.4d0,CON2=CON*CON,BIG=1.D30,NTAB=10,SAFE=2.d0)
      EXTERNAL func
CU    USES func
      INTEGER i,j
      double precision errt,fac,hh,a(NTAB,NTAB),r(n),ep,em
      if(h.eq.0.d0) pause 'h must be nonzero in dfridr'
      hh=h

         r(ir) = r(ir) + hh
         ep=func(r)
         r(ir) = r(ir) - hh
         r(ir) = r(ir) - hh
         em=func(r)
         r(ir) = r(ir) + hh
         a(1,1) = (ep-em)/(2.0d0*hh) 

c      a(1,1)=(func(x+hh)-func(x-hh))/(2.0*hh)

      err=BIG
      do 12 i=2,NTAB
        hh=hh/CON
         r(ir) = r(ir) + hh
         ep=func(r)
         r(ir) = r(ir) - hh
         r(ir) = r(ir) - hh
         em=func(r)
         r(ir) = r(ir) + hh
         a(1,i) = (ep-em)/(2.0d0*hh) 
c        a(1,i)=(func(x+hh)-func(x-hh))/(2.0*hh)
        fac=CON2
        do 11 j=2,i
          a(j,i)=(a(j-1,i)*fac-a(j-1,i-1))/(fac-1.d0)
          fac=CON2*fac
          errt=max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
          if (errt.le.err) then
            err=errt
            dfrsenc=a(j,i)
          endif
11      continue
        if(abs(a(i,i)-a(i-1,i-1)).ge.SAFE*err)return
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software <%8.
c
