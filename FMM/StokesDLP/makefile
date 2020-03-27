MEX=/share/apps/matlab/r2011b/bin/mex
#MEX=/share/apps/matlab/r2013a/bin/mex
#location of mex compiler on ronaldo

FORTRAN = ifort
FFLAGS = -c -fpic -openmp -O3
#LIB = -L/share/apps/intel/Compiler/12.1/mkl/lib/intel64 -lmkl_blas95_lp64

MEXLIBDIR = /share/apps/intel/Compiler/12.1/composer_xe_2011_sp1.7.256/compiler/lib/intel64/
#MEXLIBDIR = /share/apps/intel/Compiler/13.1/composer_xe_2013/lib/intel64
#for Ronaldo
MEXLIB = -lifcore -lirc -lsvml -limf -lgfortran -liomp5 

OBJS=\
lfmm2dpart.o\
lfmm2drouts.o\
cfmm2dpart.o\
rfmm2dpart.o\
d2mtreeplot.o\
d2tstrcr_omp.o\
l2dterms.o\
laprouts2d.o\
prini.o\
stokesDLP.o

all: cfmm2dpart d2mtreeplot d2tstrcr_omp l2dterms laprouts2d lfmm2drouts lfmm2dpart prini rfmm2dpart stokesDLP mex

mex: stokesDLPfmm.F $(OBJS) 
	$(MEX) -O stokesDLPfmm.F $(OBJS) -L $(MEXLIBDIR) $(MEXLIB)
	#$(MEX) -O -v FC="ifort" stokesDLPfmm.F $(OBJS) -L $(MEXLIBDIR) $(MEXLIB)

cfmm2dpart: cfmm2dpart.f
	$(FORTRAN) $(FFLAGS) cfmm2dpart.f

d2mtreeplot: d2mtreeplot.f
	$(FORTRAN) $(FFLAGS) d2mtreeplot.f

d2tstrcr_omp: d2tstrcr_omp.f
	$(FORTRAN) $(FFLAGS) d2tstrcr_omp.f

l2dterms: l2dterms.f
	$(FORTRAN) $(FFLAGS) l2dterms.f

laprouts2d: laprouts2d.f
	$(FORTRAN) $(FFLAGS) laprouts2d.f

lfmm2dpart: lfmm2dpart.f
	$(FORTRAN) $(FFLAGS) lfmm2dpart.f

lfmm2drouts: lfmm2drouts.f
	$(FORTRAN) $(FFLAGS) lfmm2drouts.f

prini: prini.f
	$(FORTRAN) $(FFLAGS) prini.f

rfmm2dpart: rfmm2dpart.f
	$(FORTRAN) $(FFLAGS) rfmm2dpart.f

stokesDLP: stokesDLP.f
	$(FORTRAN) $(FFLAGS) stokesDLP.f 

clean: 
	rm *.o
	cp *.mex* ../src
	rm *.mex*
