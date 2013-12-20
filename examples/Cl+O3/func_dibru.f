       double precision function func(x)
       implicit double precision (a-h,o-z) 
c      Calcula la funcion numerica  dado un punto de 3 dimensiones 
       parameter(ndis=6)
       COMMON/PRLIST/T,VP,H,TIME,NTZ,NT,ISEED0(8),NC,NX,IRESTART,NRESTART
       parameter (nmax=3)
       dimension x(nmax)
       dimension derpot(ndis)
       common/cart2/rr(ndis),xb,yb,senalp,cosalp,theta,pi,phi


c optimization of 
          rr(4) = abs(x(1))
c          rr(2) = abs(x(2))
c          rr(6) = abs(x(3))
c            cosalp=dcos(alpha) 
c            alp=dacos(cosalp) 
c            senalp=dsin(alp)
          RR(2)= DSQRT(RR(4)**2+RR(1)**2 - 2.0d0*RR(4)*RR(1)*COSALP)
c  this is for 80 degrees dihedral angle
c         RR(5)=(RR(3)*DSIN(THETA)*DSIN(PHI))**2
c         RR(5)=RR(5)+(-RR(3)*DSIN(THETA)*DCOS(PHI))**2
c         RR(5)=RR(5)+(-RR(3)*DCOS(THETA)-RR(1))**2
c         RR(5)=DSQRT(RR(5))

c  alpha calculation
c         cosalp = RR(1)**2 + RR(4)**2 - RR(2)**2
c         cosalp = cosalp/(2.0d0*RR(1)*RR(4))
c         alp=dacos(cosalp) 
c         senalp=dsin(alp)


       RR(6)=(RR(3)*DSIN(THETA)*DSIN(PHI))**2
       RR(6)=RR(6)+(RR(3)*DSIN(THETA)*DCOS(PHI)-RR(4)*SENALP)**2
       RR(6)=RR(6)+(-RR(3)*DCOS(THETA)-RR(1)+RR(4)*COSALP)**2
       RR(6)=DSQRT(RR(6))

       NC=0
       call enerder_system(rr,vocl,derpot)
c       NC=NC+1
       func=vocl
       return
       end





