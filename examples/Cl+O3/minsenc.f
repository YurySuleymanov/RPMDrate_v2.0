       program minsenc
       implicit double precision (a-h,o-z)
       parameter (nmax=3,ftol=5.0d-8)
       dimension p(nmax)
       open (unit=1,file='datos')
       pi = dacos(-1.0d0)
       read (1,*) p(1),p(2),p(3),hd
       write (*,*) '(x,y,z) in=',p(1),p(2),p(3)
       write (*,*) 'f(x,y,z) in=',func(p)
       iter=0
       n=3
       call frprmn(p,n,ftol,iter,fret,hd)
       write (*,*) 'punto del min x,y,z'
       write (*,*) p(1),p(2),p(3)
       write (*,*) 'angulo hoh en el minimo'
       arg = p(1)**2 + p(3)**2 - p(2)**2
       arg = arg/(2.0d0*p(1)*p(3))
       ang = dacos(arg)
       x = p(1)*cos(ang*0.50d0)  
       y = p(1)*sin(ang*0.50d0)  
       write (*,*) x,y
       ang = dacos(arg)*180.d0/pi
       write (*,*) ang
       write (*,*) 'valor del min en r1,r2,r3'
       write (*,*) fret
       write (*,*) 'nº de iteraciones'
       write (*,*) iter
       write (*,*) 'senc'
       write (*,*) func(p)
       stop
       end

