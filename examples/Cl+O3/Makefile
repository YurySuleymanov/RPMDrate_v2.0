################################################################################
#
#   This Makefile will compile the potential energy surface module to expose it
#   to Python using f2py.
#
################################################################################

F2PY=f2py

.PHONY: all build pyf clean

all: build

build:
	$(F2PY) -c --f77flags='-extend_source -O3 -r8  -static  -fltconsistency  -C -assume nobuffered_io ' cntl.f dvdirnew.f intern.f neigh1r.f neigh2r.f neigh3r.f radius.f makebp.f readt.f sort.f neigh.f dvdi.f readccm.f confidradtraj.f enerder_clo3.f  lisym.f ranims.f clo3.f clo3.pyf 

pyf:
	$(F2PY)  -m PES -h clo3.pyf clo3.f --overwrite-signature

clean:
	rm -f *.o *.so
