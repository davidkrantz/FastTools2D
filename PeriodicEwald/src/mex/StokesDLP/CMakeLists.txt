include_directories(
  ${CMAKE_SOURCE_DIR}/mex/common
)

## MEX functions
matlab_add_mex(
	NAME mex_stokes_dlp_real
	mex_stokes_dlp_real.cpp
)

matlab_add_mex(
	NAME mex_stokes_dlp_kspace
	mex_stokes_dlp_kspace.cpp
)
