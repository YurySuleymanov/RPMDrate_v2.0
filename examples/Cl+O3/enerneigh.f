      subroutine enerder_system(rr,vocl,derpot)
c       program enerder_system
c  this program runs a potential energy surface interpolation scheme
c  to evaluate the energy at one or many geometries supplied

      include 'traj.inc'
      PARAMETER(ND1=100)
      COMMON/PRLIST/T,VP,H,TIME,NTZ,NT,ISEED0(8),NC,NX,IRESTART,NRESTART
c  @@@@@@@@@@@@ NOTE ARRAY RSTORE @@@@@@@@@@

      dimension derpot(nbondm),rr(nbondm)
      dimension rinv(nbondm)
      parameter(anstobohr=0.529177D0)
c
	DATA COND/49.61402939d0/
      save icall
      data icall/0/
      save nout,nin,ncpre
c
      if (icall.eq.0) then

        nout=0
        nin=0
c  read in the parameters which define the system

      call readccm

c   set up the bond definitions

      call readz

c  include permutational symmetry

      call makebp

c   read in run control parameters - to get potmax and ipart, lowp, ipow

      call cntl

80    format(a80)

c   read in the initial table of Taylor series potential coefficients 

      call readt
c      write(6,*) 'readt'

      if(ipart.eq.2)then
      call radius
c      write(6,*) 'radius'
      call confidradtraj
c      write(6,*) 'confid'
      endif
c
 	close(11)
        icall=1
      endif


       do k=1,nbond
c       r(1,k)=rr(k)/anstobohr
	rinv(k)=anstobohr/rr(k)
c        write(6,*) 'r',k, rr(k)
       enddo


      if (nc.eq.0) then

      if(ipart.eq.1)then
       call neighbour(1,rinv)
      else
       call neighbour1(1,rinv)
       call neighbour2(1,rinv)
       call neighbour3(1,rinv)
      endif
      endif




c
c
c uptdating neighbourgs lists
c
        ne=nc/neighco
        nout=ne*neighco



c uptdating neighbourgs lists 
c
      if(ipart.eq.2)then
        ne=nc/neighco
        nout=ne*neighco

       if (nc.eq.nout) then
         call neighbour2(1,rinv)
       endif 

        ne=nc/neighci
        nin=ne*neighci

       if (nc.eq.nin) then
       call neighbour3(1,rinv)
       endif

      endif
c  calculate the energy at these geometries

      if(ipart.eq.2)then
       call dvdir(1,rinv)
      else
       call neighbour(1,rinv)
       call dvdi(1,rinv)
      endif

        do i=1,nbond
          derpot(i) = dvr(1,i)*cond
c          write(6,*) i,derpot(i),rr(i)
        enddo

c        write(6,*) 'en',en(1)

c  the energy kcal/mol
          vocl =en(1) + 684.553030850739d0 
         vocl  = vocl*627.50959d0
      end
