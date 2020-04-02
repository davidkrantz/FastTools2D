include_directories(
  ${CMAKE_SOURCE_DIR}/mex/common
)

## MEX functions
matlab_add_mex(
	NAME mex_stokes_dlp_real
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_dlp_real.cpp
)

matlab_add_mex(
	NAME mex_stokes_dlp_kspace
	SRC mex_stokes_dlp_kspace.cpp
)
