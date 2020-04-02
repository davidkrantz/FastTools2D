# Compilation Instructions

In principle, compiling should be as simple as running the cmake script to create the make files, and then running make. You will have to change the directories of your Matlab installation in `CMakeLists.txt`. By running everything in the `build` directory all the cmake files will be created in one place. The mex files will be placed in the a new `bin` directory. The cmake script has been tested on the following architectures:

## Ubuntu 16.04 LTS

	cd build
	cmake ..
	make

## MAC OS Mojave

Here you have to make sure that gcc and g++ are called instead of clang. I also had to manually add the include and library directories to help the compiler find the GNU scientific library.

	cd build
	export CPATH=/usr/local/include/
	export LIBRARY_PATH=/usr/local/lib/
	export LD_LIBRARY_PATH=/usr/local/lib/
	cmake -DCMAKE_C_COMPILER=/usr/local/bin/gcc-9 -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-9 .. 
	make
