# ##############################################################################
# NVIDIA HPC SDK (nvfortran) + OpenACC/CUDA
# ##############################################################################

include(GEPSLibPaths)

add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Mfree -Ofast>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-r8>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Mpreprocess>")
add_compile_options(
  "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Mbyteswapio -Minline>")

if(${USE_OMP})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-mp=multicore>")
endif()
if(${USE_PAR})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Mconcur>")
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

# Host BLAS/LAPACK: NVHPC-shipped libs (not system OpenBLAS)
set(_nvhpc_math_lib "")
if(DEFINED ENV{NVHPC_ROOT} AND NOT "$ENV{NVHPC_ROOT}" STREQUAL "")
  set(_nvhpc_math_lib "$ENV{NVHPC_ROOT}/compilers/lib")
elseif(DEFINED ENV{NVHPC} AND NOT "$ENV{NVHPC}" STREQUAL "")
  # $NVHPC may be SDK root; prefer active compilers/lib via nvfortran sibling
  get_filename_component(_nvfc "${CMAKE_Fortran_COMPILER}" DIRECTORY)
  get_filename_component(_nvhpc_math_lib "${_nvfc}/../lib" ABSOLUTE)
endif()
if(NOT _nvhpc_math_lib OR NOT EXISTS "${_nvhpc_math_lib}")
  get_filename_component(_nvfc "${CMAKE_Fortran_COMPILER}" DIRECTORY)
  get_filename_component(_nvhpc_math_lib "${_nvfc}/../lib" ABSOLUTE)
endif()

find_library(GEPS_NVHPC_BLAS NAMES blas
  HINTS "${_nvhpc_math_lib}" NO_DEFAULT_PATH)
find_library(GEPS_NVHPC_LAPACK NAMES lapack
  HINTS "${_nvhpc_math_lib}" NO_DEFAULT_PATH)
if(GEPS_NVHPC_BLAS AND GEPS_NVHPC_LAPACK)
  message(STATUS "BLAS/LAPACK: NVHPC ${_nvhpc_math_lib}")
  link_directories("${_nvhpc_math_lib}")
  link_libraries(${GEPS_NVHPC_BLAS} ${GEPS_NVHPC_LAPACK})
else()
  message(FATAL_ERROR
    "NVIDIA build requires NVHPC BLAS/LAPACK under compilers/lib "
    "(looked in '${_nvhpc_math_lib}'). Load nvhpc module / set NVHPC_ROOT.")
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

link_libraries(-ltirpc -lm -lcurl -lhdf5_hl -lhdf5 -lgfortran -lcusparse
               -lcudart)

add_compile_options(-DUSE_CUDA=1)
add_compile_options(-acc=gpu -gpu=cc${GPU_ARCHS},cuda${CUDA_RUNTIME_VERSION}
                    -Minfo=accel -cuda -cudalib=cublas,cufft,cusolver,nccl)
link_libraries(-acc=gpu -gpu=cc${GPU_ARCHS},cuda${CUDA_RUNTIME_VERSION} -cuda
               -cudalib=cublas,cufft,cusolver,nccl)
if(${USE_PCAST})
  add_compile_options(-gpu=redundant)
  add_compile_options(-DUSE_PCAST=1)
  link_libraries(-gpu=redundant)
endif()

if(CPL_INCLUDE_DIR)
  include_directories(${CPL_INCLUDE_DIR})
endif()
if(MCT_INCLUDE_DIR)
  include_directories(${MCT_INCLUDE_DIR})
endif()
