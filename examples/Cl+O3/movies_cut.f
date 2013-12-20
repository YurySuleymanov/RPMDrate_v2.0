      program movies
      implicit double precision (a-h,o-z)
      character*2 alabel(4)

      parameter (npo=191,npf=95)
      dimension o1c(npo,3)
      dimension o2c(npo,3)
      dimension o3c(npo,3)
      dimension clc(npo,3)


      open(unit=1,file='traj_2.xyz')
      open(unit=2,file='traj_2_mov.xyz')

      do i=1,npo
       read(1,*) natom
       read(1,*) itraj
       read(1,*) alabel(1),(o1c(i,j),j=1,3)
       read(1,*) alabel(2),(o2c(i,j),j=1,3)
       read(1,*) alabel(3),(o3c(i,j),j=1,3)
       read(1,*) alabel(4),(clc(i,j),j=1,3)
      enddo

      do i=1,npo,2
       write(2,*) natom
       write(2,*) itraj
       write(2,*) alabel(1),(o1c(i,j),j=1,3)
       write(2,*) alabel(2),(o2c(i,j),j=1,3)
       write(2,*) alabel(3), (o3c(i,j),j=1,3)
       write(2,*) alabel(4),(clc(i,j),j=1,3)
      enddo

       stop
       end
