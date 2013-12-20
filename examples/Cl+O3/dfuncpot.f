      subroutine dfunc(p,xi,n,h)
      implicit double precision (a-h,o-z)
      parameter(nmax=3)
       common/cart2/rr(6),xb,yb,cost3,alpha,cosalp,theta
      dimension p(n),xi(nmax)
      dimension derpot(6)
c optimization of
c optimization of
          rr(4) = p(1)
          rr(2) = p(2)
c          rr(5) = p(3)
       call enerder_system(rr,vocl,derpot)
          xi(1)=derpot(4)
          xi(2)=derpot(2)
c          xi(3)=derpot(6)
c          xi(3)=derpot(5)

c      do i = 1,n
c          xi(i)=dfrsend_i,p,n,h,err)
c      enddo
      return 
      end
