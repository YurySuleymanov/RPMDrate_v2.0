
      subroutine neighbour1(j,rinv)

c  the original neighbour subroutine which calculates totsum from scratch

      include 'traj.inc'

      dimension vt(ndatam,ngroupm),rinv(nbondm)

c    assume a weighting function of 1/r**(2*ipow)

      totsum=0.d0

c      do i=1,nbond
c         rinv(i)=1.d0/r(j,i)
c      enddo

      do m=1,ndata
      do n=1,ngroup

         vt(m,n)=0.d0

         do i=1,nbond
            i1=ip(i,n)
            t=rinv(i)-rda(m,i1)
            vt(m,n)=t*t/sumb(m,i1)+vt(m,n)
         enddo

c patch to match in with dvdir.f

c         vt(m,n)=vt(m,n)
         vt(m,n)=1.d0/(vt(m,n)**lowp+vt(m,n)**ipow)

c  end of patch

         totsum=totsum+vt(m,n)

      enddo
      enddo

c  we set up the massive neighbour list
c  wtol is passed in common/shepdata
c  nforc is the number of neighbours
c  mda identifies the data point in POT that neighbour k comes from
c  nda identifies which permutation neighbour k is

      nforco=0
      test=wtol*outer
55    continue

      do m=1,ndata
      do n=1,ngroup

         wt1=vt(m,n)/totsum

         if (wt1.gt.test) then
             nforco=nforco+1

      if (nforco.gt.maxfo) then
      nforco=0
      test=2.d0*test
      go to 55
      
c        write(10,*)' *** ABORT ***'
c        write(10,*)' nforco.gt.maxfo'
c        write(10,*) ' neigh1: nforco =',nforco
c        stop
      endif

             mdao(nforco)=m
             ndao(nforco)=n
         endif

      enddo
      enddo


c  debug

c     write(6,*) 
c     write(6,*) ' neigh1: nforco =',nforco
c     write(6,*) ' mda(k) =',(mda(k),k=1,nforc)
c     write(6,*) ' nda(k) =',(nda(k),k=1,nforc)

      return
      end

