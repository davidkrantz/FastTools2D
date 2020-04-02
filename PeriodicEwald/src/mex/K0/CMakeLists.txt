include_directories(
  ${CMAKE_SOURCE_DIR}/mex/common
)

## MEX functions
matlab_add_mex(
	NAME mex_K0_real
	SRC ${CMAKE_SOURCE_DIR}/mex/common/bessellib.cpp ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_K0_real.cpp
)
target_link_libraries(mex_K0_real gsl gslcblas)

matlab_add_mex(
	NAME mex_se_K0_2p_spreading
	SRC mex_se_K0_2p_spreading.c
)

matlab_add_mex(
	NAME mex_se_K0_2p_quadrature
	SRC mex_se_K0_2p_quadrature.c
)

matlab_add_mex(
	NAME mex_compK0xy
	SRC ${CMAKE_SOURCE_DIR}/mex/common/bessellib.cpp mex_compK0xy.cpp
)
target_link_libraries(mex_compK0xy gsl gslcblas)

