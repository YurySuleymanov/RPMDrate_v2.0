
	subroutine readccm
	include 'traj.inc'

	character*80 icomm, ielements, ialpha

        open (unit=7,file='IN_SYSTEM',status='old')

c  read comment ; a header for the SYSTEM file

        read(7,80)icomm
80      format(a80)

c  read in the actual number of atoms

        read(7,80)icomm

        read(7,*) natom

        nbond=natom*(natom-1)/2

        n3=3*natom

        nint=n3-6
        if(natom.eq.2)nint=1

c  read in the element labels for the atoms

        read(7,80)icomm

        read(7,80) ielements

c extract the elemental labels from ielements

        ialpha=ielements

        do i=1,natom-1
        length=len(ialpha)
        ifin=index(ialpha,',')
         if(ifin.gt.3)then
         write(6,*)'an atomic label in SYSTEM has too many characters'
         stop
         endif
        lab(i)=ialpha(1:ifin-1)
        ialpha=ialpha(ifin+1:length)
        enddo
        ifin=index(ialpha,' ')
        lab(natom)=ialpha(1:ifin-1)

c        length=len(ialpha)
c         if(length.gt.2)then
c         write(6,*)'last atomic label in SYSTEM has too many characters'
c         stop
c         endif
c       lab(natom)=ialpha(1:length)

c  read in the atomic masses

        read(7,80)icomm
        read(7,*)(amas(i),i=1,natom)

c  read in the cutoff lengths rmin(i,j)

        read(7,80)icomm
        read(7,80)icomm
        do i=1,natom-1
        do j=i+1,natom
        read(7,*) rmin(i,j)
        enddo
        enddo


        close (unit=7)

        return
        end

