# ##############################################################################
# Fujitsu Fortran (frtpx / tcsds) — CPU, no CUDA
# ##############################################################################

include(GEPSLibPaths)

add_compile_options(
  "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Kfast,ocl,autoobjstack>")
add_compile_options(
  "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-X9 -Free -CcdRR8 -Cpp -x- -fw -Knofp_relaxed>"
)
add_compile_options(
  "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Ec -Nlst=a,lst=d,lst=i,lst=p,lst=t,lst=x -Koptmsg=2 >"
)
# Fujitsu SSL2 BLAMP (compiler math); not OpenBLAS/MKL
message(STATUS "BLAS/LAPACK: Fujitsu SSL2BLAMP")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-SSL2BLAMP >")

set(CMAKE_C_FLAGS_RELEASE "")
add_compile_options("$<$<COMPILE_LANGUAGE:C>:SHELL:-Xg>")

if(${USE_OMP})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Kopenmp>")
endif()
if(${USE_PAR})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Kparallel>")
endif()

if(GEPS_CPL_ROOT)
  set(CPL_INCLUDE_DIR "${GEPS_CPL_ROOT}/include")
  if(NOT EXISTS "${CPL_INCLUDE_DIR}")
    set(CPL_INCLUDE_DIR "${GEPS_CPL_ROOT}")
  endif()
endif()
if(GEPS_MCT_ROOT)
  set(MCT_INCLUDE_DIR "${GEPS_MCT_ROOT}/include")
endif()

find_package(NetCDF REQUIRED COMPONENTS Fortran)
link_libraries(-lnetcdf -lnetcdff -SSL2BLAMP)
find_package(DMS REQUIRED)

if(${USE_NDMS} STREQUAL "OFF")
  if(GEPS_DMS38KEY_ROOT)
    link_directories(${GEPS_DMS38KEY_ROOT}/lib)
  endif()
  link_libraries(-lrdms -lgdbm -ltirpc)
endif()

if(${USE_MPMD})
  find_package(MPMD REQUIRED)
endif()

if(GEPS_FFTW_ROOT)
  link_directories(${GEPS_FFTW_ROOT}/lib)
endif()
link_libraries(-lfftw3_threads -lfftw3 -lfftw3f_threads -lfftw3f)

if(GEPS_CPL_ROOT)
  link_directories(${GEPS_CPL_ROOT}/lib)
endif()
link_libraries(-lcpl)

if(GEPS_MCT_ROOT)
  link_directories(${GEPS_MCT_ROOT}/lib)
endif()
link_libraries(-lmct -lmpeu)

if(CPL_INCLUDE_DIR)
  include_directories(${CPL_INCLUDE_DIR})
endif()
if(MCT_INCLUDE_DIR)
  include_directories(${MCT_INCLUDE_DIR})
endif()
