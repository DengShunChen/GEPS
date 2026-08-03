# ##############################################################################
# Intel oneAPI (ifx/ifort + MPI) — CPU, no CUDA
# ##############################################################################

include(GEPSLibPaths)

add_compile_options(
  "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-O3 -free -r8 -fpp -convert big_endian>")
add_compile_options(
  "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-align all -fno-alias -fp-model precise>")
add_compile_options("$<$<COMPILE_LANGUAGE:C>:SHELL:-O3 -fPIC>")

if(${USE_OMP})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-qopenmp>")
  add_link_options(-qopenmp)
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

# Intel oneAPI Math Kernel Library (required for intel builds)
if(NOT DEFINED ENV{MKLROOT} OR "$ENV{MKLROOT}" STREQUAL "")
  message(FATAL_ERROR
    "Intel build requires MKL (MKLROOT unset). "
    "Load mkl/2023.1.0 (geps/1.0 should do this) or set MKLROOT.")
endif()
message(STATUS "MKLROOT: $ENV{MKLROOT}")
message(STATUS "BLAS/LAPACK: Intel MKL (-qmkl=sequential)")
# sequential MKL: app already uses OpenMP; avoid nested oversubscription
link_libraries(-qmkl=sequential)
if(DEFINED ENV{MKLROOT})
  include_directories("$ENV{MKLROOT}/include")
  link_directories("$ENV{MKLROOT}/lib/intel64" "$ENV{MKLROOT}/lib")
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
message(STATUS "Looking for NetCDF (triggers FindMPI)...")
find_package(NetCDF REQUIRED COMPONENTS Fortran)
message(STATUS "NetCDF found")
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
