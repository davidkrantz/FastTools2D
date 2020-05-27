include_directories(
  ${CMAKE_SOURCE_DIR}/mex/common
)

## MEX functions
matlab_add_mex(
	NAME mex_stokes_slp_real
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_slp_real.cpp
)

matlab_add_mex(
	NAME mex_stokes_slp_kspace
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_slp_kspace.cpp
)

matlab_add_mex(
	NAME mex_stokes_slp_pressure_real
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_slp_pressure_real.cpp
)

matlab_add_mex(
	NAME mex_stokes_slp_pressure_kspace
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_slp_pressure_kspace.cpp
)

matlab_add_mex(
	NAME mex_stokes_slp_gradient_kspace
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_slp_gradient_kspace.cpp
)

matlab_add_mex(
	NAME mex_stokes_slp_gradient_real
	SRC ${CMAKE_SOURCE_DIR}/mex/common/ewald_tools.cpp mex_stokes_slp_gradient_real.cpp
)

target_link_libraries(mex_stokes_slp_real gomp)
target_link_libraries(mex_stokes_slp_kspace gomp)
