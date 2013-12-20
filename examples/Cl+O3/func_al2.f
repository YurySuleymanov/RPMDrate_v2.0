       double precision function func(x)
       implicit double precision (a-h,o-z) 
c      Calcula la funcion numerica  dado un punto de 6 dimensiones 
       PARAMETER(ND1=100,NDP=10)
      COMMON/PRLIST/T,VP,H,TIME,NTZ,NT,ISEED0(8),NC,NX,IRESTART,NRESTART
       parameter (nmax=3)
       dimension x(nmax)
       dimension derpot(6)
       common/cart2/rr(6),xb,yb,cost3,alpha,cosalp,theta,pi





c optimization of 
          rr(4) = abs(x(1))
          rr(2) = abs(x(2))
          rr(6) = abs(x(3))
c          rr(5) = abs(x(3))

            RR(2)= DSQRT(RR(4)**2+RR(1)**2 - 2.0d0*RR(4)*RR(1)*COSALP)
            RR(6)= DSQRT(RR(5)**2+RR(4)**2- 2.0d0*RR(5)*RR(4)*cost3)


         cosalp= RR(2)**2-RR(4)**2-RR(1)**2
         cosalp=-cosalp/(2.0d0*RR(4)*RR(1))
          if (cosalp.lt.-1.0d0) then
             alp=pi-dacos(cosalp)
             cosalp=dcos(alp)
          endif
c comment for no alpha optimization
            cosalp=dcos(alpha) 

            alp=dacos(cosalp) 
            
              

          cost2=RR(3)**2-RR(5)**2-RR(1)**2
          cost2=-cost2/(2.0d0*RR(5)*RR(1))
          if (cost2.lt.-1.0d0) then
              theta2=pi-dacos(cost2)
          else
          theta2 = dacos(cost2)
          endif

          cost3=dcos(alp-theta2)

             RR(2)= DSQRT(RR(4)**2+RR(1)**2 - 2.0d0*RR(4)*RR(1)*COSALP)
             RR(6)= DSQRT(RR(5)**2+RR(4)**2- 2.0d0*RR(5)*RR(4)*cost3)

c              write(6,*) rr(2),rr(3),rr(4)

       call enerder_system(rr,vocl,derpot)
c       NC=NC+1
       func=vocl
       return
       end





