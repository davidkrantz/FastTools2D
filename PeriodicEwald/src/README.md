# Compilation Instructions

In principle, compiling should be just running the cmake script to create the make files, and then running make. You will have to change the directories of your Matlab installation in `CMakeLists.txt`.

## Ubuntu

	cmake .
	make

## MAC OS Mojave

Here you have to make sure that gcc and g++ are called instead of clang. I also had to manually add the include and library directories to help the compiler find the GNU scientific library.

	export CPATH=/usr/local/include/
	export LIBRARY_PATH=/usr/local/lib/
	export LD_LIBRARY_PATH=/usr/local/lib/
	cmake -DCMAKE_C_COMPILER=/usr/local/bin/gcc-9 -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-9 . 
	make
