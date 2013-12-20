
      subroutine neighbour(j,rinv)

      include 'traj.inc'

      dimension vt(ndatam,ngroupm)
      dimension rinv(nbondm)

c    assume a weighting function of 1/r**(2*ipow)


      totsum=0.d0

c  !ocl vi(k)
      do m=1,ndata
      do n=1,ngroup

         vt(m,n)=0.d0

         do i=1,nbond
c            t=1.d0/r(j,i)-1.d0/rda(m,ip(i,n))
             t=rinv(i)-1.d0/rda(m,ip(i,n))
c             write(6,*) rinv(i),1.d0/rda(m,ip(i,n))
c            t=rinv(i)-rda(m,ip(i,n))
c           t=r(j,i)-rda(m,ip(i,n))
            vt(m,n)=t*t+vt(m,n)
         enddo
         
         vt(m,n)=1.d0/vt(m,n)**ipow

         totsum=totsum+vt(m,n)

      enddo
      enddo

c  we set up the massive neighbour list
c  wtol is passed in common/shepdata
c  nforc is the number of neighbours
c  mda identifies the data point in POT that neighbour k comes from
c  nda identifies which permutation neighbour k is

      nforc=0

      do m=1,ndata
      do n=1,ngroup

         wt1=vt(m,n)/totsum

         if (wt1.gt.wtol) then
             nforc=nforc+1
             mda(nforc)=m
             nda(nforc)=n
         endif

      enddo
      enddo

      if (nforc.gt.maxf) then
         write(10,*)' *** ABORT ***'
         write(10,*)' nforc.gt.maxf'
         stop
      endif

c  debug


c     write(6,*)
c     write(6,*) ' nforc =',nforc
c     write(6,*) ' mda(k) =',(mda(k),k=1,nforc)
c     write(6,*) ' nda(k) =',(nda(k),k=1,nforc)

      return
      end

