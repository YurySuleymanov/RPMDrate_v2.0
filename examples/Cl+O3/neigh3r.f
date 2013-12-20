
      subroutine neighbour3(j,rinv)

c  modified 19-7-96 by kt to approximate totsum by totsum from previous dvdi
c  and to calculate an inner neighbour list from the outer one

      include 'traj.inc'

      dimension rinv(nbondm)

c      do i=1,nbond
c         rinv(i)=1.d0/r(j,i)
c      enddo

c  we set up the inner neighbour list
c  wtol is the cutoff
c  nforc is the number of neighbours
c  mda identifies the data point in POT that neighbour k comes from
c  nda identifies which permutation neighbour k is

      nforc=0
      test=wtol*totsum
      if(test.gt.1.d0)then
      npow=lowp
      else
      npow=ipow
      endif
      test=1.d0/test**(1.d0/npow)

      do k=1,nforco

         wt1=0.d0

         do i=1,nbond
            i1=mdao(k)
            i2=ip(i,ndao(k))
            t=rinv(i)-rda(i1,i2)
            wt1=t*t/sumb(i1,i2)+wt1
         enddo


         if (wt1.lt.test) then
             nforc=nforc+1

      if (nforc.gt.maxf) then
         write(10,*)' *** ABORT ***'
         write(10,*)' nforc.gt.maxf'
         write(10,*) ' neigh3: nforc =',nforc
         stop
      endif

             mda(nforc)=mdao(k)
             nda(nforc)=ndao(k)
         endif

      enddo

c  debug

c     write(6,*)
c     write(6,*) ' neigh3: nforc =',nforc
c     write(6,*) ' mda(k) =',(mda(k),k=1,nforc)
c     write(6,*) ' nda(k) =',(nda(k),k=1,nforc)

      return
      end

