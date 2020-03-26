include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
)


## MEX functions
matlab_add_mex(
	NAME mex_K0_real
	SRC bessellib.cpp ewald_tools.cpp mex_K0_real.cpp
)
target_link_libraries(mex_K0_real gsl gslcblas)

## MEX functions
matlab_add_mex(
	NAME mex_se_K0_2p_spreading
	SRC mex_se_K0_2p_spreading.c
)

## MEX functions
matlab_add_mex(
	NAME mex_se_K0_2p_quadrature
	SRC mex_se_K0_2p_quadrature.c
)

## MEX functions
matlab_add_mex(
	NAME mex_compK0xy
	SRC bessellib.cpp mex_compK0xy.cpp
)
target_link_libraries(mex_compK0xy gsl gslcblas)

