MEX=/Applications/MATLAB_R2019a.app/bin/mex

FORTRAN = gfortran-9
FFLAGS = -c -fpic -fopenmp -O3
#FFLAGS = -c -fpic -O3

#MEXLIBDIR = /Applications/MATLAB_R2019a.app/bin/maci64/:/Applications/MATLAB_R2019a.app/sys/os/maci64/
#MEXLIB = -libifcore -libirc -libsvml -libimf -libgfortran -libiomp5 
MEXLIB = -lgomp

OBJS=\
lfmm2drouts.o\
cfmm2dpart.o\
rfmm2dpart.o\
d2mtreeplot.o\
d2tstrcr_omp.o\
l2dterms.o\
laprouts2d.o\
prini.o\
stokesSLP.o

all: cfmm2dpart d2mtreeplot d2tstrcr_omp l2dterms laprouts2d lfmm2drouts prini rfmm2dpart stokesSLP mex

mex: stokesSLPfmm.F $(OBJS) 
	$(MEX) -O stokesSLPfmm.F $(OBJS) $(MEXLIB)

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

lfmm2drouts: lfmm2drouts.f
	$(FORTRAN) $(FFLAGS) lfmm2drouts.f

prini: prini.f
	$(FORTRAN) $(FFLAGS) prini.f

rfmm2dpart: rfmm2dpart.f
	$(FORTRAN) $(FFLAGS) rfmm2dpart.f

stokesSLP: stokesSLP.f
	$(FORTRAN) $(FFLAGS) stokesSLP.f

clean: 
	rm *.o
	cp *.mex* ../src
	rm *.mex*
