       double precision function func(x)
       implicit double precision (a-h,o-z) 
c      Calcula la funcion numerica  dado un punto de 6 dimensiones 
       PARAMETER(ND1=100,NDP=10)
       parameter (nmax=3,n=4)
       dimension x(nmax)
       dimension r(6),rr(6),derpot(6),q(3*n)
       common/cart2/tr(nd1*(nd1-1)/2,3)
       parameter(xlim1=11.03386,xlim2=10.0,xlim3=8.876)

      save icall
      data icall/0/



c        parameter(xlim=9.0d0)
c      x(1) y x(2) y x(3) se introducen desde fuera.       
        func=0.0d0
      if (icall.eq.0) then
   
        
       q(1)= 0.687275233490091        
       q(2)=1.09564438564590
       q(3)=0.0
       q(4)=0.0
       q(5)=0.0
       q(6)=0.0
       q(7)= 0.687275233490091        
       q(8)=-1.09564438564590
       q(9)=0.0
       q(10)=0.0
       q(11)=0.0
       q(12)=8.0



       natom=n
       i=0
       j=0
       do ni=1,natom-1
       do m=ni+1,natom
       i=i+1
       rr(i)=0.d0
       do k=1,3
        ki=ni+2*(ni-1)+k-1
        ji=ni-1+j+2*m+k-1
        tr(i,k)=q(ki)-q(ji)
        rr(i)=rr(i)+tr(i,k)**2
       enddo
       rr(i)=dsqrt(rr(i))
       j=j+1
       enddo
       j=0
       enddo

           x(1)=rr(1)
           x(2)=rr(2)
           x(3)=rr(4)
           write(6,*) 'r',rr(1),rr(2),rr(3)
           write(6,*) 'r',rr(4),rr(5),rr(6)
           write(8,*)  n
           write(8,*) 
           write(8,*) q(1),q(2),q(3)
           write(8,*) q(4),q(5),q(6)
           write(8,*) q(7),q(8),q(9)
           write(8,*) q(10),q(11),q(12)
        icall=1
      endif


         rr(1) = x(1)
         rr(2) = x(2)
         rr(3) = xlim1
         rr(4) = x(3)
         rr(5) = xlim2
         rr(6) = xlim3
c  geometry fro o + oclo
c          rr(1) = xlim
c          rr(2) = x(1)
c          rr(3) = xlim
c          rr(4) = x(3)
c          rr(5) = x(2)
c          rr(6) = xlim
c           call  clo3sur(r1,r2,r3,r4,r5,r6,vocl)
           call enerder_system(rr,vocl,derpot)
c        call  clo3sur_mod(r1,r2,r3,r4,r5,r6,vocl)
c   original farantos
c        vocl = vocl+ 0.294600844
c   gross-billing leps
c        vocl = vocl + 0.293299909092498
c
c        vocl = vocl*27.2114d0*23.06054d0
c  valor de o + oclo -124.42976124768823 
c        vocl= vocl+124.42976124768823         
       func=vocl
       return
       end





