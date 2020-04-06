# Doubly Periodic Spectral Ewald

This library provides tools to compute sums of infinitely doubly periodic Green functions. The current library contains the single- and double-layer potentials of the Stokes equations. 

## Compilation Instructions

In principle, compiling should be as simple as running the cmake script to create the make files, and then running make. You will have to change the directories of your Matlab installation in `CMakeLists.txt`. By running everything in the `build` directory all the cmake files will be created in one place. The mex files will be placed in a new `bin` directory. The cmake script has been tested on the following architectures:

### Ubuntu 16.04 LTS / MATLAB 2019a

	cd src/build
	cmake ..
	make

### MAC OS Mojave / MATLAB 2018a

Here you have to make sure that gcc and g++ are called instead of clang. I also had to manually add the include and library directories to help the compiler find the GNU scientific library.

	cd src/build
	export CPATH=/usr/local/include/
	export LIBRARY_PATH=/usr/local/lib/
	export LD_LIBRARY_PATH=/usr/local/lib/
	cmake -DCMAKE_C_COMPILER=/usr/local/bin/gcc-9 -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-9 .. 
	make

## Testing

In the `tests` directory, there are several tests that can be used to verify the compilation worked correctly:
* consistency_test.m: checks that changing the Ewald parameters and enlarging the periodic box by adding replicates of the reference cell don't change the results
* direct_sums_test.m: compares the spectral Ewald implementation to matlab direct sums of the real and Fourier parts. The Matlab direct sum does not truncate in real space, and in Fourier space it does not spread the data to a uniform grid and thus does not use FFTs
* timings_test.m: checks the timings of the code for increasing numbers of source and target points. The timing should scale as O(N log N), where N is total number of points
* stresslet_indentity_test.m: verifies the stresslet identity for points inside and outside a circle

## References

Ewald decomposition of two-dimensional Stokeslet and stresslet along with error estimates:

Pålsson and Torberg, 2019. *An integral equation method for closely interacting surfactant-covered droplets in wall-confined Stokes flow*. <arxiv.org/pdf/1909.12581.pdf>

Lindbo and Tornberg, 2010. *Spectrally accurate fast summation for periodic Stokes potentials*. J. Comp. Phys 229(23). <sciencedirect.com/science/article/pii/S0021999110004730>
