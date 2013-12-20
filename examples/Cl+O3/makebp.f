      subroutine makebp

      include 'traj.inc'

      character*80 icomm

c  read in the atomic permutations into the nperm array

       open(unit=8, file='IN_ATOMPERMS',status='old')


cc read in a header

       read(8,80)icomm

80    format(a80)

c  read in a comment line

       read(8,80)icomm

c  read in the order of the group

       read(8,*)ngroup

       nstore=ngroup

c  read in the atomic perms

       do i=1,ngroup

c  read in a comment line or separator
       read(8,80)icomm
       do k=1,natom
       read(8,*)idummy,nperm(i,k)
       enddo
c  idummy will be equal to k

       enddo

       close(unit=8)

c  now set up the ip array of bond perms by looking at nperm

      do n=1,ngroup

      ic=0
      do i=1,natom-1
      do j=i+1,natom
      ic=ic+1
      do k=1,nbond
      if(mb(k).eq.nperm(n,i).and.nb(k).eq.nperm(n,j))ip(k,n)=ic
      if(nb(k).eq.nperm(n,i).and.mb(k).eq.nperm(n,j))ip(k,n)=ic
      enddo

      enddo
      enddo

      enddo

      end
