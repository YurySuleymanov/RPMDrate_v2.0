      subroutine confidradtraj

      include 'traj.inc'

       dimension z(nintm),dei(nintm),dtdrm(nbondm),
     .           delrv(nbondm),dtdr(nbondm),asb(nbondm)
     .           ,bigsort(ndatam)

c      sigma = 0.0001d0
      sigma2=sigma*sigma

      open(unit=3,file='CONTRAJ',status='unknown')

      do j=1,ndata


         rcut=rads(j)+1.d-12

       sum=0.d0
       do i=1,nbond
       sumb(j,i)=0.d0
       enddo

      do m=1,ndata
      if(m.eq.j)go to 10

      do n=1,ngroup

       vt=0.d0

       do i=1,nbond
        t=rda(m,ip(i,n))-rda(j,i)
       vt=t*t+vt
       enddo
c  rda is a reciprocal bondlength

c  only consider the first nneigh neighbours



      if (vt.gt.rcut) go to 222


c the Taylor expansion

      do i=1,nint
         z(i)=0.d0
         do l=1,nbond
            z(i) = ut(j,i,ip(l,1))*rda(m,ip(l,n)) + z(i)
         enddo
          z(i)=(z(i)-zda(j,i))
      enddo

      taydif=v0(m)-v0(j)
      do i=1,nint
       taydif=taydif-z(i)*(v1(j,i)+0.5d0*v2(j,i)*z(i))
	 dei(i) = z(i)*v2(j,i) + v1(j,i)
      enddo

      sum=sum+taydif**2/vt**3


c  transform derivatives of taylor poly wrt normal local coords into 
c  derivs wrt inverse bondlengths

      do i=1,nbond
            dtdr(i)=0.d0
            dtdrm(i)=0.d0
            do l=1,nint

               dtdr(i)=dei(l)*ut(j,l,ip(i,1)) + dtdr(i)
               dtdrm(i)=v1(m,l)*ut(m,l,ip(i,n)) + dtdrm(i)

            enddo
      enddo

      do i=1,nbond
      delrv(i)=(rda(m,ip(i,n))-rda(j,i))*(dtdrm(i)-dtdr(i))
      sumb(j,i)=sumb(j,i)+delrv(i)**2/vt**3
      enddo


c end the ngroup loop
222   continue
      enddo

10    continue

c end the m ndata loop
      enddo

c     do i=1,nbond
c     if(sumb(j,i).lt.1.d-15)then
c     write(6,*)j,i,sumb(j,i)
c     stop
c     endif
c     enddo

      rcon=(nneigh*sigma2/sum)**(1.d0/3.d0)

      do i=1,nbond
      sumb(j,i)=(nneigh*sigma2/sumb(j,i))**(1.d0/3.d0)
      enddo

c  test to see if some sumb have been wiildly extrapolated

120    format(9e8.2)


c  replace rads by rcon
      if(ipart.eq.2)then
      rads(j)=rcon
      else
      rads(j)=1.d0
      endif

c end the j ndata loop
      enddo

c      stop

c  try to eliminate abnormally large radii and sumb

c  first evaluate averages

      ard=0.d0
      do m=1,ndata
      ard=ard+rads(m)
      enddo
      ard=ard/dfloat(ndata)

      do n=1,nbond
      asb(n)=0.d0
      do m=1,ndata
      asb(n)=asb(n)+sumb(m,n)
      enddo
      asb(n)=asb(n)/dfloat(ndata)
      enddo

c     do m=1,ndata
c     if(rads(m).gt.1.5d0*ard)rads(m)=1.5d0*ard
c     do n=1,nbond
c     if(sumb(m,n).gt.1.5d0*asb(n))sumb(m,n)=1.5d0*asb(n)
c     enddo
c     enddo

c      do n=1,nbond
c      write(6,*) n,asb(n)
c      enddo

      ndt2=ndata/2

      do n=1,nbond

      do k=1,ndata
       bigsort(k)=sumb(k,n)
      enddo
      call sort(ndata,bigsort)
      do k=1,ndata
       if(sumb(k,n).gt.2.5d0*bigsort(ndt2))sumb(k,n)=2.5d0*bigsort(ndt2)
      enddo

      write(3,*) n,bigsort(ndt2)
      enddo

c  write rads and sumb to CON

      do j=1,ndata
      write(3,*) j
      write(3,*)rads(j),(sumb(j,n),n=1,nbond)
      enddo

      close(unit=3)

234   format(15e9.2)
7763  continue 
      return
      end
