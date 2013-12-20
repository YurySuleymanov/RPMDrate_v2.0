       double precision function func(x)
       implicit double precision (a-h,o-z) 
c      Calcula la funcion numerica  dado un punto de 3 dimensiones 
       parameter(ndis=6)
       COMMON/PRLIST/T,VP,H,TIME,NTZ,NT,ISEED0(8),NC,NX,IRESTART,NRESTART
       parameter (nmax=3)
       dimension x(nmax)
       dimension derpot(ndis)
       common/cart2/rr(ndis),xb,yb,cost3,alpha,cosalp,theta,pi,beta


c optimization of 
          rr(4) = abs(x(1))
          rr(6) = abs(x(2))

            RR(2)= DSQRT(RR(4)**2+RR(1)**2 - 2.0d0*RR(4)*RR(1)*COSALP)
c  this is for 80 degrees dihedral angle
            RR(6)= DSQRT(RR(3)**2+RR(2)**2-2.0d0*RR(2)*RR(3)*DCOS(BETA))
c  this is for coplanar reaction plots
c            RR(6)= DSQRT(RR(5)**2+RR(4)**2- 2.0d0*RR(5)*RR(4)*cost3)


c try beta optimization
c         cosalp= RR(2)**2-RR(4)**2-RR(1)**2
c         cosalp=-cosalp/(2.0d0*RR(4)*RR(1))
c         RR(2)= DSQRT(RR(4)**2+RR(1)**2 - 2.0d0*RR(4)*RR(1)*COSALP)

c try beta optimization
         cosbet= RR(6)**2-RR(3)**2-RR(2)**2
         cosbet=-cosbet/(2.0d0*RR(3)*RR(2))
         if (cosbet.lt.-1.0d0) cosbet=-1.0d0

c  this is for 80 degrees dihedral angle (i believe)
            RR(6)= DSQRT(RR(3)**2+RR(2)**2-2.0d0*RR(2)*RR(3)*cosbet)
c this is for coplanar reaction plots
c           RR(6)= DSQRT(RR(5)**2+RR(4)**2- 2.0d0*RR(5)*RR(4)*cost3)


       call enerder_system(rr,vocl,derpot)
       func=vocl
       return
       end





