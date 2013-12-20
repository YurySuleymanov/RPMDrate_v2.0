      subroutine initialize_potential()
!        implicit double precision (a-h,o-z)
       include 'traj.inc'
	double precision :: T,vp,h,time
	integer :: NTZ,NT,Iseed0, NC,NX,IRESTART,NRESTART
        INTEGER :: init
        parameter(ndis=6)
        dimension derpot(ndis),rr(ndis)

        COMMON/PRLIST/T,VP,H,TIME,NTZ,NT,ISEED0(8),NC,NX,IRESTART,NRESTART

        data init /0/
        save init

* Energy in kcal/mol, distances in Angstrom

c     rr(1)  O_1 - O_2
c     rr(2)  O_1 - O_3
c     rr(3)  O_1 - Cl
c     rr(4)  O_2 - O_3
c     rr(5)  O_2 - Cl
c     rr(6)  O_3 - Cl


c structure of TS ab initio

      rr(1)=1.29484
      rr(2)=2.19044
      rr(3)=2.44837
      rr(4)=1.31476
      rr(5)=3.14225
      rr(6)=3.58405


      NC=0 
      call enerder_clo3(rr,vocl,derpot)
      NC=NC+1
  
       NC = 0 

	write(*,*) 'V0=',vocl 
* energy of Cl+O3 channel (with Cl at infinity with respect to O3)

        init = 1
      
      end subroutine initialize_potential


       subroutine get_potential(q,Natoms,Nbeads,V,dVdq,info)
        implicit double precision (a-h,o-z)
        INTEGER :: Natoms, Nbeads
        DOUBLE PRECISION :: q(3,Natoms,Nbeads)
        DOUBLE PRECISION :: dVdq(3,Natoms,Nbeads), V(Nbeads)
        integer :: i, j, k, i1, j1, k1
        double precision :: derpot(6),rr(6)
	double precision :: T,vp,h,time
	integer :: NTZ,NT,Iseed0, NC,NX,IRESTART,NRESTART
	double precision :: XO1O2,YO1O2,ZO1O2
	double precision :: XO1O3,YO1O3,ZO1O3
	double precision :: XO2O3,YO2O3,ZO2O3
	double precision :: XO1Cl,YO1Cl,ZO1Cl
	double precision :: XO2Cl,YO2Cl,ZO2Cl
	double precision :: XO3Cl,YO3Cl,ZO3Cl
	double precision :: RO1O2,RO1O3,RO2O3,RO1Cl,RO2Cl,RO3Cl, vocl
        COMMON/PRLIST/T,VP,H,TIME,NTZ,NT,ISEED0(8),NC,NX,IRESTART,NRESTART
	info = 0 
	ndis = 6 

	do k = 1, Nbeads

	XO1O2 = q(1,3,k) - q(1,2,k)  
	YO1O2 = q(2,3,k) - q(2,2,k)  
	ZO1O2 = q(3,3,k) - q(3,2,k)  
	RO1O2 = dsqrt(XO1O2**2+YO1O2**2+ZO1O2**2)

	XO1O3 = q(1,3,k) - q(1,1,k)  
	YO1O3 = q(2,3,k) - q(2,1,k)  
	ZO1O3 = q(3,3,k) - q(3,1,k)  
	RO1O3 = dsqrt(XO1O3**2+YO1O3**2+ZO1O3**2)

	XO1Cl = q(1,3,k) - q(1,4,k)  
	YO1Cl = q(2,3,k) - q(2,4,k)  
	ZO1Cl = q(3,3,k) - q(3,4,k)  
	RO1Cl = dsqrt(XO1Cl**2+YO1Cl**2+ZO1Cl**2)

	XO2O3 = q(1,2,k) - q(1,1,k)  
	YO2O3 = q(2,2,k) - q(2,1,k)  
	ZO2O3 = q(3,2,k) - q(3,1,k)  
	RO2O3 = dsqrt(XO2O3**2+YO2O3**2+ZO2O3**2)

	XO2Cl = q(1,2,k) - q(1,4,k)  
	YO2Cl = q(2,2,k) - q(2,4,k)  
	ZO2Cl = q(3,2,k) - q(3,4,k)  
	RO2Cl = dsqrt(XO2Cl**2+YO2Cl**2+ZO2Cl**2)

	XO3Cl = q(1,1,k) - q(1,4,k)  
	YO3Cl = q(2,1,k) - q(2,4,k)  
	ZO3Cl = q(3,1,k) - q(3,4,k)  
	RO3Cl = dsqrt(XO3Cl**2+YO3Cl**2+ZO3Cl**2)

	rr(1) = RO1O2*0.529177d0
	rr(2) = RO1O3*0.529177d0
	rr(3) = RO1Cl*0.529177d0
	rr(4) = RO2O3*0.529177d0
	rr(5) = RO2Cl*0.529177d0
	rr(6) = RO3Cl*0.529177d0

        NC=0 

        call enerder_clo3(rr,vocl,derpot)

        NC=NC+1

        V(k) = vocl/627.50959d0
!
!
	derpot(1) = derpot(1)/49.61402939d0 !*0.529177d0/627.50959d0
	derpot(2) = derpot(2)/49.61402939d0 !*0.529177d0/627.50959d0
	derpot(3) = derpot(3)/49.61402939d0 !*0.529177d0/627.50959d0
	derpot(4) = derpot(4)/49.61402939d0 !*0.529177d0/627.50959d0
	derpot(5) = derpot(5)/49.61402939d0 !*0.529177d0/627.50959d0
	derpot(6) = derpot(6)/49.61402939d0 !*0.529177d0/627.50959d0

	dVdq(1,4,k) = derpot(3)*(-XO1Cl/RO1Cl) + derpot(5)*(-XO2Cl/RO2Cl) + derpot(6)*(-XO3Cl/RO3Cl)
	dVdq(2,4,k) = derpot(3)*(-YO1Cl/RO1Cl) + derpot(5)*(-YO2Cl/RO2Cl) + derpot(6)*(-YO3Cl/RO3Cl)
	dVdq(3,4,k) = derpot(3)*(-ZO1Cl/RO1Cl) + derpot(5)*(-ZO2Cl/RO2Cl) + derpot(6)*(-ZO3Cl/RO3Cl)

	dVdq(1,3,k) = derpot(1)*(XO1O2/RO1O2) + derpot(2)*(XO1O3/RO1O3) + derpot(3)*(XO1Cl/RO1Cl)
	dVdq(2,3,k) = derpot(1)*(YO1O2/RO1O2) + derpot(2)*(YO1O3/RO1O3) + derpot(3)*(YO1Cl/RO1Cl)
	dVdq(3,3,k) = derpot(1)*(ZO1O2/RO1O2) + derpot(2)*(ZO1O3/RO1O3) + derpot(3)*(ZO1Cl/RO1Cl)

	dVdq(1,2,k) = derpot(1)*(-XO1O2/RO1O2) + derpot(4)*(XO2O3/RO2O3) + derpot(5)*(XO2Cl/RO2Cl)
	dVdq(2,2,k) = derpot(1)*(-YO1O2/RO1O2) + derpot(4)*(YO2O3/RO2O3) + derpot(5)*(YO2Cl/RO2Cl)
	dVdq(3,2,k) = derpot(1)*(-ZO1O2/RO1O2) + derpot(4)*(ZO2O3/RO2O3) + derpot(5)*(ZO2Cl/RO2Cl)

	dVdq(1,1,k) = derpot(2)*(-XO1O3/RO1O3) + derpot(4)*(-XO2O3/RO2O3) + derpot(6)*(XO3Cl/RO3Cl)
	dVdq(2,1,k) = derpot(2)*(-YO1O3/RO1O3) + derpot(4)*(-YO2O3/RO2O3) + derpot(6)*(YO3Cl/RO3Cl)
	dVdq(3,1,k) = derpot(2)*(-ZO1O3/RO1O3) + derpot(4)*(-ZO2O3/RO2O3) + derpot(6)*(ZO3Cl/RO3Cl)

	end do
!
!
      return
      end 
!



!
!	FINITE DIFFERENCE CHECK OF DERIVATIVES
!
!  test 
!x	double precision :: shift, dvdqr(3,Natoms), dvdql(3,Natoms)
!x        DOUBLE PRECISION :: dVdq_test(3,Natoms,Nbeads)
!x
!x
!x	shift = 1.d-5
!x	do i1 = 1,3
!x	do j1 = 1,Natoms
!x	q(i1,j1,k) = q(i1,j1,k) - shift		 
!x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!x	XO1O2 = q(1,3,k) - q(1,2,k)  
!x	YO1O2 = q(2,3,k) - q(2,2,k)  
!x	ZO1O2 = q(3,3,k) - q(3,2,k)  
!x	RO1O2 = dsqrt(XO1O2**2+YO1O2**2+ZO1O2**2)
!x
!x	XO1O3 = q(1,3,k) - q(1,1,k)  
!x	YO1O3 = q(2,3,k) - q(2,1,k)  
!x	ZO1O3 = q(3,3,k) - q(3,1,k)  
!x	RO1O3 = dsqrt(XO1O3**2+YO1O3**2+ZO1O3**2)
!x
!x	XO1Cl = q(1,3,k) - q(1,4,k)  
!x	YO1Cl = q(2,3,k) - q(2,4,k)  
!x	ZO1Cl = q(3,3,k) - q(3,4,k)  
!x	RO1Cl = dsqrt(XO1Cl**2+YO1Cl**2+ZO1Cl**2)
!x
!x	XO2O3 = q(1,2,k) - q(1,1,k)  
!x	YO2O3 = q(2,2,k) - q(2,1,k)  
!x	ZO2O3 = q(3,2,k) - q(3,1,k)  
!x	RO2O3 = dsqrt(XO2O3**2+YO2O3**2+ZO2O3**2)
!x
!x	XO2Cl = q(1,2,k) - q(1,4,k)  
!x	YO2Cl = q(2,2,k) - q(2,4,k)  
!x	ZO2Cl = q(3,2,k) - q(3,4,k)  
!x	RO2Cl = dsqrt(XO2Cl**2+YO2Cl**2+ZO2Cl**2)
!x
!x	XO3Cl = q(1,1,k) - q(1,4,k)  
!x	YO3Cl = q(2,1,k) - q(2,4,k)  
!x	ZO3Cl = q(3,1,k) - q(3,4,k)  
!x	RO3Cl = dsqrt(XO3Cl**2+YO3Cl**2+ZO3Cl**2)
!x
!x	rr(1) = RO1O2*0.529177d0
!x	rr(2) = RO1O3*0.529177d0
!x	rr(3) = RO1Cl*0.529177d0
!x	rr(4) = RO2O3*0.529177d0
!x	rr(5) = RO2Cl*0.529177d0
!x	rr(6) = RO3Cl*0.529177d0
!x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!x        call enerder_clo3(rr,vocl,derpot)         
!x	dvdql(i1,j1) = vocl/627.50959d0
!x
!x	q(i1,j1,k) = q(i1,j1,k) + 2.d0*shift		 
!x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!x	XO1O2 = q(1,3,k) - q(1,2,k)  
!x	YO1O2 = q(2,3,k) - q(2,2,k)  
!x	ZO1O2 = q(3,3,k) - q(3,2,k)  
!x	RO1O2 = dsqrt(XO1O2**2+YO1O2**2+ZO1O2**2)
!x
!x	XO1O3 = q(1,3,k) - q(1,1,k)  
!x	YO1O3 = q(2,3,k) - q(2,1,k)  
!x	ZO1O3 = q(3,3,k) - q(3,1,k)  
!x	RO1O3 = dsqrt(XO1O3**2+YO1O3**2+ZO1O3**2)
!x
!x	XO1Cl = q(1,3,k) - q(1,4,k)  
!x	YO1Cl = q(2,3,k) - q(2,4,k)  
!x	ZO1Cl = q(3,3,k) - q(3,4,k)  
!x	RO1Cl = dsqrt(XO1Cl**2+YO1Cl**2+ZO1Cl**2)
!x
!x	XO2O3 = q(1,2,k) - q(1,1,k)  
!x	YO2O3 = q(2,2,k) - q(2,1,k)  
!x	ZO2O3 = q(3,2,k) - q(3,1,k)  
!x	RO2O3 = dsqrt(XO2O3**2+YO2O3**2+ZO2O3**2)
!x
!x	XO2Cl = q(1,2,k) - q(1,4,k)  
!x	YO2Cl = q(2,2,k) - q(2,4,k)  
!x	ZO2Cl = q(3,2,k) - q(3,4,k)  
!x	RO2Cl = dsqrt(XO2Cl**2+YO2Cl**2+ZO2Cl**2)
!x
!x	XO3Cl = q(1,1,k) - q(1,4,k)  
!x	YO3Cl = q(2,1,k) - q(2,4,k)  
!x	ZO3Cl = q(3,1,k) - q(3,4,k)  
!x	RO3Cl = dsqrt(XO3Cl**2+YO3Cl**2+ZO3Cl**2)
!x
!x	rr(1) = RO1O2*0.529177d0
!x	rr(2) = RO1O3*0.529177d0
!x	rr(3) = RO1Cl*0.529177d0
!x	rr(4) = RO2O3*0.529177d0
!x	rr(5) = RO2Cl*0.529177d0
!x	rr(6) = RO3Cl*0.529177d0
!x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!x        call enerder_clo3(rr,vocl,derpot)         
!x	dvdqr(i1,j1) = vocl/627.50959d0
!!!!
!x	dvdq_test(i1,j1,k) = (dvdqr(i1,j1) - dvdql(i1,j1))/(2.d0*shift)
!x	q(i1,j1,k) = q(i1,j1,k) - shift		 
!x	end do	
!x	end do
!x	write(2345,*) '*****************************'	
!x
!x	Do i = 1, 3
!x	Do j = 1, Natoms
!x	Do k = 1, Nbeads	
!x	write(2345,*) dvdq(i,j,k), dvdq_test(i,j,k)	
!x	end do
!x	end do
!x	end do
