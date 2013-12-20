
      subroutine neighbour2(j,rinv)

c  modified 19-7-96 by kt to approximate totsum by totsum from previous dvdi
c  and to calculate an outer neighbour list

      include 'traj.inc'

      dimension rinv(nbondm)

c      do i=1,nbond
c         rinv(i)=1.d0/r(j,i)
c      enddo

c  we set up the outer neighbour list
c  outer is the factor by which the outer wtol is less than the inner wtol
c  wtol is the cutoff for allowing points into the neibour list
c  nforco is the number of outer neighbours
c  mdao identifies the outer data point in POT that neighbour k comes from
c  ndao identifies which permutation outer neighbour k is

      nforco=0
      test=wtol*outer*totsum
      if(test.gt.1.d0)then
      npow=lowp
      else
      npow=ipow
      endif
      test=1.d0/test**(1.d0/npow)
55    continue
      do m=1,ndata
      do n=1,ngroup

         wt1=0.d0

         do i=1,nbond
            i1=ip(i,n)
            t=rinv(i)-rda(m,i1)
            wt1=t*t/sumb(m,i1)+wt1
         enddo


         if (wt1.lt.test) then
             nforco=nforco+1

      if (nforco.gt.maxfo) then
      nforco=0
      test=0.5d0*test
      go to 55
c        write(10,*)' *** ABORT ***'
c        write(10,*)' nforco.gt.maxfo'
c        write(10,*) ' neigh2: nforco =',nforco
c        stop
      endif

             mdao(nforco)=m
             ndao(nforco)=n
         endif

      enddo
      enddo


c  debug

c     write(6,*)
c     write(6,*) ' neigh2: nforco =',nforco
c     write(6,*) ' mda(k) =',(mda(k),k=1,nforc)
c     write(6,*) ' nda(k) =',(nda(k),k=1,nforc)

      return
      end

