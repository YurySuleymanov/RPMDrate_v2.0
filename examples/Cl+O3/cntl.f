      subroutine cntl

c   this subroutine reads in run control parameters from file on unit 5

      include 'traj.inc'


      call lisym(11,'-')
      write(11,*) ' *** subroutine cntl:',
     .  ' read in run control parameters ***'

      open (unit=7,file ='IN_CNT',status='old')

c   read in title of run

      read(7,80) title
   80 format(a80)

c  read in the number 1...7 defining the units

c     read(7,80) comlin
c     read(7,*) nunits

      nunits=5
c---------------------------------------------------------------------
c
c   calculate and scale model parameters
c   units for inputing  masses,energies and lengths:
c   if (nunits) then (mass  energy   length     scale factor) are:
c         1           g/mol kcal/mol angstroms  2.045482828d13
c         2           g/mol kj/mol   angstroms  1.0d13
c         3           g/mol cm-1     angstroms  1.093739383d12
c         4           g/mol eV       angstroms  9.822693571d13
c         5           g/mol hartrees bohr       9.682886569d14
c         6           g/mol aJ           angstroms  9.124087591d13
c         7           g/mol 1.d+05 j/mol angstroms  1.0d15
c
c on output masses, energies, lengths and times are scaled by
c amsc, esc, rsc and tsc (in seconds), respectively and
c must be multiplied by these scale factors to be converted
c to the stated units.

      amsc=1.0d0
      esc=1.0d0
      rsc=1.0d0
      scalef=1.d0
      if(nunits.eq.1)then
        scalef=2.045482828d13
        msunit=' g/mol'
        enunit=' kcal/mol'
        lnunit=' angstroms'
        eukjpm=4.184d0
      elseif(nunits.eq.2)then
        scalef=1.0d13
        msunit=' g/mol'
        enunit=' kJ/mol'
        lnunit=' angstroms'
        eukjpm=1.00d0
      elseif(nunits.eq.3)then
        scalef=1.093739383d12
        msunit=' g/mol'
        enunit=' cm-1'
        lnunit=' angstroms'
        eukjpm=1.196265839d-02
      elseif(nunits.eq.4)then
        scalef=9.822693571d13
        msunit=' g/mol'
        enunit=' eV'
        lnunit=' angstroms'
        eukjpm=96.48530899d0
      elseif(nunits.eq.5)then
        scalef=9.682886569d14
        msunit=' g/mol'
        enunit=' hartree'
        lnunit=' bohr'
        eukjpm=2625.499964d0
      elseif(nunits.eq.6)then
        scalef=9.124087591d13
        msunit=' g/mol'
        enunit=' 0.1382 aJ'
        lnunit=' angstroms'
        eukjpm=83.24897441d0
      elseif(nunits.eq.7)then
        scalef=1.0d14
        msunit=' g/mol'
        enunit=' 1.0d+05 J/mol'
        lnunit=' angstroms'
        eukjpm=1.00d+02
      endif
      tsc=dsqrt(amsc/esc)*rsc/scalef
      tups=tsc*1.0d12

c   Avagadro's number

      avanum=6.0221367d23

c   speed of light in cm s-1 and cm tu-1

      ccmps=2.99792458d10
      ccmptu=ccmps*tsc

c   Planck's constant in J s, kJ/mol s and eu tu

      hjs=6.6260755d-34
      hkjpms=hjs*avanum*1.0d-03
      heutu=(hkjpms/tsc)/eukjpm

c---------------------------------------------------------------------
c   echo some parameters

      write(11,*) title

c echo info from readccm

      write(11,*) 'The number of atoms = ',natom
      write(11,*)
      write(11,*) 'The atomic labels and masses are:'
      do i=1,natom
      write(11,*) lab(i), amas(i)
      enddo
      write(11,*)
      write(11,*) 'The bonded atoms'
      do i=1,nbond
      write(11,*) mb(i), nb(i)
      enddo
      write(11,*)
      write(11,*) 'The number of group elements = ',ngroup
      write(11,*)

c      write(11,*) 'The permutations of atoms in the group are:'
c      write(11,*)
c      do i=1,ngroup
c       do k=1,natom
c       write(11,*)k,nperm(i,k)
c       enddo
c       write(11,*)
c      enddo
      write(11,*) 'The permutations of bonds in the group are:'
      write(11,*)
      do i=1,ngroup
       write(11,*)
       write(11,8769)' Reordered atoms ',(nperm(i,k),k=1,natom)
       do k=1,nbond
       write(11,*)k, ip(k,i)
       enddo
       write(11,*)
      enddo
8769  format(a18,20i3)

      write(11,*) 'The units are defined as atomic units in cntl.f'
      write(11,*)

      write(11,40)esc,enunit,rsc,lnunit,amsc,msunit,
     .  tsc,tups,scalef
   40 format(
     ./,'  unit scaling parameters',
     ./,'  esc=',g12.5,a10,'  rsc=',g12.5,a10,'  amsc=',g12.5,a10,
     ./,'  tsc=',g17.10,' s =',g17.10,' ps  scalef=',g17.10,' s-1')

      write(11,45)avanum,ccmps,ccmptu,hjs,hkjpms,heutu,eukjpm
   45 format(
     ./,'  Avagadro''s number   =',g17.10,
     ./,'  Speed of light       =',g17.10,' cm s-1  =',
     .g17.10,' cm tu-1',
     ./,'  Planck''s constant   =',g17.10,' J s     =',
     .g17.10,' kJ/mol s',
     ./,'                                           =',g17.10,' eu tu',
     ./,'  1 eu = 1 energy unit = ',g17.10,' kJ/mol')

      call lisym(11,'-')


c     call lisym(11,'-')
      write(11,81) title
   81 format(1x,a79)
c     call lisym(11,'-')

c   read in time step,no. of traj. step, nprint,no. of molecule etc.

      read(7,80) comlin
      write(11,81) comlin
      read(7,*) delta,nstep,nprint,nmol
      write(11,100) delta,nstep,nprint,nmol
  100 format(1x,g18.10,3i10)

      read(7,80) comlin
      write(11,81) comlin
      read(7,*) efraga,efragb,etrans
      write(11,140) efraga,efragb,etrans

      read(7,80) comlin
      write(11,81) comlin
      read(7,*)econ
      write(11,140)econ

  140 format(1x,3g18.10)

      read(7,80) comlin
      write(11,81) comlin
      read(7,*) bipmax,distab
      write(11,140) bipmax,distab
      if(bipmax.gt.distab)then
        write(11,*) ' WARNING: bipmax =',bipmax,'.gt.distab=',distab
        write(11,*) ' RESET distab = bipmax = ',bipmax
        distab=bipmax
        write(11,81) comlin
        write(11,140) bipmax,distab
      endif


c  determine one or twp part weight function is used, ipart=1,2

       read(7,80) comlin
       write(11,81)comlin
       read(7,*)ipart
       write(11,91)ipart
91     format(2x,'we will use the ',i1,' part weight function')

      if(ipart.eq.1.or.ipart.eq.2)go to 92
      write(11,*)' incorrect value for 1 or 2 part weight function'
      stop
92    continue

c  read in the parameters defining the weight function, the neighbour
c  list and the number of timesteps between calls to neighbour

       read(7,80) comlin
       write(11,81)comlin
       read(7,*)lowp,ipow
       write(11,*)lowp,ipow

       read(7,80) comlin
       write(11,81)comlin
       read(7,*)wtol,outer
       write(11,*)wtol,outer

       read(7,80) comlin
       write(11,81)comlin
       read(7,*)neighco,neighci
       write(11,*)neighco,neighci

c  to allow one part neighbour to work

      neighc=neighci


c  read in the sampling fraction used in outp to determine number
c  of trajectory points output in file TOUT, and number of
c  data points chosen

       read(7,80) comlin
       write(11,81)comlin
       read(7,*) sample
       write(11,*) sample

       read(7,80) comlin
       write(11,81)comlin
       read(7,*) nsel
       write(11,*) nsel

c  read in energy maximum and minimum with which to screen data

       read(7,80)comlin
       write(11,81)comlin
       read(7,*) potmax,vmin
       write(11,*) potmax,vmin

c  read in flag to determine whether low energy trajectories are stopped

       read(7,80)comlin
       write(11,81)comlin
       read(7,*) nlowstop
       write(11,*) nlowstop

c read in the multipicitive factor used in radius.f, confidrad.f

       read(7,80)comlin
       write(11,81)comlin
       read(7,*) nneigh
       write(11,*) nneigh

c  Enter the energy error tolerance for confidrad.f

       read(7,80)comlin
       write(11,81)comlin
       read(7,*) sigma
       write(11,*) sigma

c calculate pi, conversion factors

      pi=2.d0*dacos(0.d0)
      todeg=45.d0/atan(1.d0)
      toang=0.52917706d0
      sq2=dsqrt(2.d0)

      write(11,*) ' conversion factors: '
      write(11,*) ' todeg=',todeg,' toang=',toang
      write(11,*) ' pi   =',pi,' sq2  =',sq2


c  scale atomic masses

      do 150 i=1,natom
        amas(i)=amas(i)/amsc
  150 continue


      close(unit=7)

      return
      end

c  --------------------------------------------------------------

      subroutine readz

c  modified jan 7 1996 by kt 

      include 'traj.inc'
 
      call lisym(11,'-')
      write(11,*) ' *** subroutine readz: creat mb,nb read in fragment data'
      write(11,*) ' '

c  given that we need all the atom-atom bonds, we just assign the indirect
c  addresses of the bonds, rather than the traditional read from ZMA

        k=1
        do i=1,natom-1
        do j=i+1,natom

           mb(k)=i
           nb(k)=j
           k=k+1

        enddo
        enddo

c   read fragment atom numbers

      open (unit=7,file='IN_ZMA',status='old')

80    format(a80)
81    format(1x,a79)

      read(7,80)comlin
      write(11,81)comlin

      read(7,*)nfraga,nfragb
      write(11,*)nfraga,nfragb

      if(nfraga+nfragb.ne.natom)then
        write(11,*) 'ABORT: nfraga+nfragb.ne.natom'
        write(11,*) nfraga,'+',nfragb,'.ne.',natom
        stop
      endif

      read(7,80)comlin
      write(11,81)comlin

      if(nfraga.gt.0)then
      read(7,*)(ifraga(i),i=1,nfraga)
      write(11,*)(ifraga(i),i=1,nfraga)
      endif

      read(7,80)comlin
      write(11,81)comlin

      if(nfragb.gt.0)then
      read(7,*)(ifragb(i),i=1,nfragb)
      write(11,*)(ifragb(i),i=1,nfragb)
      endif
      if (nfraga .le. 1) go to 3100

      nbfraga=0
      do 3000 i=1,nbond
	do 2900 j=1,nfraga
	  do 2800 k=1,nfraga
            if ((mb(i).eq.ifraga(j)).and.(nb(i).eq.ifraga(k))) then
              nbfraga=nbfraga+1
              ibfraga(nbfraga)=i
            endif
 2800     continue 
 2900   continue
 3000 continue
	
      write (11,*) 'no. of bonds in fragment a is', nbfraga
      write (11,*) 'bonds in fragment a are'
      write (11,*) (ibfraga(i),i=1,nbfraga)

 3100 continue

      if (nfragb .le. 1) go to 4100

      nbfragb=0
      do 4000 i=1,nbond
	do 3900 j=1,nfragb
	  do 3800 k=1,nfragb
            if ((mb(i).eq.ifragb(j)).and.(nb(i).eq.ifragb(k))) then
              nbfragb=nbfragb+1
              ibfragb(nbfragb)=i
              endif
 3800     continue 
 3900   continue
 4000 continue

      write (11,*) 'no. of bonds in fragment b is', nbfragb
      write (11,*) 'bonds in fragment b are'
      write (11,*) (ibfragb(i),i=1,nbfragb)
      
 4100 continue

      close(unit=7)

      return
      end

