       double precision function func(x)
       implicit double precision (a-h,o-z) 
c      Calcula la funcion numerica  dado un punto de 6 dimensiones 
       PARAMETER(ND1=100,NDP=10)
      COMMON/PRLIST/T,VP,H,TIME,NTZ,NT,ISEED0(8),NC,NX,IRESTART,NRESTART
       parameter (nmax=4)
       dimension x(nmax)
       dimension derpot(6)
       common/cart2/rr(6),xb,yb,cost3,alpha,cosalp,theta,pi





c optimization of 
          rr(4) = abs(x(4))
c          rr(2) = abs(x(2))
c          rr(5) = abs(x(3))
          rr(6) = abs(x(1))
       
          if(rr(4).le.1.0d0) rr(4)=1.0d0
          if(rr(6).le.1.0d0) rr(6)=1.5d0

          RR(5)= DSQRT(RR(1)**2+RR(3)**2- 2.0d0*RR(1)*RR(3)*DCOS(THETA))
          RR(2)= DSQRT(RR(4)**2+RR(1)**2 - 2.0d0*RR(4)*RR(1)*COSALP)

       call enerder_system(rr,vocl,derpot)
c       NC=NC+1
       func=vocl
       return
       end





