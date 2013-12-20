      subroutine dfunc(p,xi,n,h)
      implicit double precision (a-h,o-z)
      parameter(nmax=50)
      dimension p(n),xi(nmax)
      do i = 1,n
          xi(i)=dfrsenc(i,p,n,h,err)
      enddo
      return 
      end
