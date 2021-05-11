###
###		 DJ subroutine library gamma.a
###

FC = gfortran

SOURCES = atg.f \
	  depth-ns.f \
	  depth-scv.f \
	  derthe.f \
	  eosall.f \
	  alphabeta.f \
	  eos8d.f \
	  e-solve.f \
	  gamma-errors.f \
	  gamma-n.f \
	  goor.f \
	  gamma-qdr.f \
	  get-lunit.f \
	  goor-solve.f \
	  index.f \
	  neutral-surfaces.f \
	  ocean-test.f \
	  scv-solve.f \
	  sig-vals.f \
	  stp-interp.f \
	  svan.f \
	  theta.f \
	  read-nc.f


## for profiling, append -pg to FFLAGS and to example compilation.
## Run ./a.out after building example.
## Then run gprof a.out gmon.out > profile.txt

OBJECTS = $(SOURCES:.f=.o)

FFLAGS =  -O4 -fPIC -ffixed-line-length-none -g -fbacktrace -fcheck=all -Wall

.KEEP_STATE :
.PRECIOUS :


################################################################################


gamma.a: gamma.a($(OBJECTS))
	ranlib gamma.a

example: example.f gamma.a
	$(FC) example.f gamma.a -L/usr/local/lib -lnetcdff -lhdf5_hl

## The -L option above is for the benefit of OS X; see ../setup.py
## for an explanation.

