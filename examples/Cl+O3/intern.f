
      subroutine intern(j)
 
c  calculates bondlengths and their derivatives - NOT inverses !!!!!!!!

      include 'traj.inc'

c  bond lengths

      do i=1,nbond
         r(j,i)=sqrt((c(j,mb(i),1)-c(j,nb(i),1))**2
     .              +(c(j,mb(i),2)-c(j,nb(i),2))**2
     .              +(c(j,mb(i),3)-c(j,nb(i),3))**2)
      enddo

c  bond length derivatives

      do i=1,nbond
         dr(j,i,mb(i),1)=(c(j,mb(i),1)-c(j,nb(i),1))/r(j,i)
         dr(j,i,mb(i),2)=(c(j,mb(i),2)-c(j,nb(i),2))/r(j,i)
         dr(j,i,mb(i),3)=(c(j,mb(i),3)-c(j,nb(i),3))/r(j,i)
         dr(j,i,nb(i),1)=-dr(j,i,mb(i),1)
         dr(j,i,nb(i),2)=-dr(j,i,mb(i),2)
         dr(j,i,nb(i),3)=-dr(j,i,mb(i),3)
      enddo

      return
      end
