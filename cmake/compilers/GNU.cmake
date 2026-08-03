# ##############################################################################
# GNU gfortran + OpenMPI — CPU, no CUDA
# ##############################################################################

include(GEPSLibPaths)

# Ensure .f90 runs through cpp (#ifdef). gfortran uses traditional cpp:
# '#' must be in column 1 (sources normalized accordingly).
set(CMAKE_Fortran_PREPROCESS ON)

add_compile_options(
  "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-O3 -ffree-form -fdefault-real-8 -cpp>")
add_compile_options(
  "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-fconvert=big-endian -fno-range-check -ffree-line-length-none>"
)
# -fallow-argument-mismatch is GCC>=10; RHEL8 system gfortran is 8.x
if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
  add_compile_options(
    "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-fallow-argument-mismatch>")
else()
  add_compile_options(
    "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Wno-argument-mismatch>")
endif()
add_compile_options("$<$<COMPILE_LANGUAGE:C>:SHELL:-O3 -fPIC>")

if(${USE_OMP})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-fopenmp>")
  add_link_options(-fopenmp)
endif()

if(GEPS_NETCDF_ROOT)
  set(NetCDF_Fortran_INCLUDE_DIRS "${GEPS_NETCDF_ROOT}/include")
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

# Prefer GEPS_LIB OpenBLAS (GNU has no vendor BLAS; intel/nvidia/fujitsu use theirs)
if(GEPS_OPENBLAS_ROOT)
  find_library(GEPS_OPENBLAS NAMES openblas
    HINTS "${GEPS_OPENBLAS_ROOT}/lib" "${GEPS_OPENBLAS_ROOT}/lib64"
    NO_DEFAULT_PATH)
endif()
if(NOT GEPS_OPENBLAS)
  find_library(GEPS_OPENBLAS NAMES openblas openblasp PATHS /usr/lib64 /lib64)
endif()
if(GEPS_OPENBLAS)
  message(STATUS "BLAS/LAPACK: OpenBLAS ${GEPS_OPENBLAS}")
  link_libraries(${GEPS_OPENBLAS})
else()
  message(FATAL_ERROR
    "GNU build requires OpenBLAS (GEPS_LIB openblas/* or system libopenblas). "
    "Other compilers use vendor math (MKL / NVHPC / SSL2).")
endif()

if(${USE_NDMS} STREQUAL "OFF")
  if(GEPS_DMS38KEY_ROOT)
    link_directories(${GEPS_DMS38KEY_ROOT}/lib)
  endif()
  link_libraries(-lrdms -lgdbm)
endif()

if(GEPS_ZLIB_ROOT)
  link_directories(${GEPS_ZLIB_ROOT}/lib)
endif()
link_libraries(-lz)

if(GEPS_NETCDF_ROOT)
  link_directories(${GEPS_NETCDF_ROOT}/lib)
endif()
find_package(NetCDF REQUIRED COMPONENTS Fortran)
link_libraries(-lnetcdf -lnetcdff)

if(GEPS_FFTW_ROOT)
  link_directories(${GEPS_FFTW_ROOT}/lib)
endif()
link_libraries(-lfftw3_threads -lfftw3 -lfftw3f_threads -lfftw3f)

if(GEPS_OPERLIB_ROOT)
  link_directories(${GEPS_OPERLIB_ROOT}/lib)
endif()
link_libraries(-lnwp)

if(GEPS_CPL_ROOT)
  link_directories(${GEPS_CPL_ROOT}/lib)
endif()
link_libraries(-lcpl)

if(GEPS_MCT_ROOT)
  link_directories(${GEPS_MCT_ROOT}/lib)
endif()
link_libraries(-lmct -lmpeu)

link_libraries(-ltirpc -lm -lcurl -lhdf5_hl -lhdf5)

if(CPL_INCLUDE_DIR)
  include_directories(${CPL_INCLUDE_DIR})
endif()
if(MCT_INCLUDE_DIR)
  include_directories(${MCT_INCLUDE_DIR})
endif()
