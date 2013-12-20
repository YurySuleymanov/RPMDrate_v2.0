      subroutine readt

      include 'traj.inc'

      open (unit=7,file='POT',status='unknown')
c      open (unit=7,file='/home/jfc/hh2oas/POT',status='unknown') 
      call lisym(11,'-')
      write(11,*) ' '
      write(11,*) ' *** subroutine readt:',
     .  ' read in potential energy surface parameters ***'
      write(11,*) ' '

c   read comment/header

      read(7,80)comlin
   80 format(a80)
      write(11,81)comlin
   81 format(1x,a79)

c  read in the bond lengths, the matrix defining the normal local coords, 
c  the potential energy, the gradients and diagonal force constants (wrt to 
c  the normal local coords) at all of the Ndata points.

c  temporarily store vmin which was read in cntl

      vminstore=vmin

      j=1
10    continue
         read(7,80,end=20)comlin

         read(7,*) (rda(j,k),k=1,nbond)
c         write(6,1971) (rda(j,k),k=1,nbond)

         do i=1,nint
            read(7,*) (ut(j,i,k),k=1,nbond)
         enddo

         read(7,*) v0(j)
c         write(6,*) rda(j,55)*0.529177,v0(j)

c  find the minimum potential in the original data set

         if (j.eq.1) vmin=v0(j)
         if (v0(j).lt.vmin) vmin=v0(j)

         read(7,*) (v1(j,k),k=1,nint)
    
         read(7,*) (v2(j,k),k=1,nint)

         if (v0(j).gt.potmax) goto 10
            
         j=j+1
c      if (j.eq.711) goto 20
         go to 10

20    continue
1971  format(15(1x,g10.5))
    
      ndata=j-1
      write(11,*)' Number of PES data points = ',ndata

c set vmin to be the minimum of the readin value in cntl and the
c minimum found in the POT file

      write(11,*)' lowest potential energy in the POT file   = ',vmin

      if(vminstore.lt.vmin)vmin=vminstore

      write(11,*)' lowest potential energy expected from IN_CNT  = '
     .             ,vminstore

      write(11,*)' lowest potential energy   = ',vmin

      close (unit=7)

c  construct the local coords for each data point
c  remembering that rda are read in as bond lengths NOT inverses 

      do n=1,ndata
      do j=1,nint
         zda(n,j)=0.d0
         do k=1,nbond
            zda(n,j)=ut(n,j,k)/rda(n,k) + zda(n,j)
         enddo
      enddo
      enddo

      if(ipart.eq.2)then

c  invert data bonds once and for all

      do n=1,ndata
      do i=1,nbond
         rda(n,i)=1.d0/rda(n,i)
      enddo
      enddo

      endif

c debug

c     do n=1,ndata
c        write(6,*) ' '
c        write(6,1812) (zda(n,i),i=1,6)
c        write(6,1812) (zda(n,i),i=7,12)
c     enddo

1812  format(2x,6g12.5)

c     stop

      return
      end

