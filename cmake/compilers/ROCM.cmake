# ##############################################################################
# AMD ROCm (amdclang / amdflang) + OpenMPI
#   USE_HIP=OFF  — CPU OpenMP (original path)
#   USE_HIP=ON   — OpenMP offload + hipBLAS/hipFFT/hipSOLVER/RCCL (MI300X/gfx942)
# NVIDIA CUDA/OpenACC stays on compilers/NVIDIA.cmake (GEPS_COMPILER=nvidia).
# ##############################################################################

include(GEPSLibPaths)

# Ensure .f90 runs through cpp (#ifdef).
set(CMAKE_Fortran_PREPROCESS ON)

# gfortran-only flags that amdflang rejects:
#   -ffree-line-length-none  -fallow-argument-mismatch  -fno-range-check
add_compile_options(
  "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-O3 -ffree-form -fdefault-real-8 -fdefault-double-8 -cpp -fPIC>"
)
add_compile_options(
  "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-fconvert=big-endian>")
add_compile_options("$<$<COMPILE_LANGUAGE:C>:SHELL:-O3 -fPIC>")
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:SHELL:-O3 -fPIC>")

if(${USE_OMP})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-fopenmp>")
  add_compile_options("$<$<COMPILE_LANGUAGE:C>:SHELL:-fopenmp>")
  add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:SHELL:-fopenmp>")
  add_link_options(-fopenmp)
endif()

# lld (amdflang) errors on sbyte_ in both libw3.a and libg2.a; bfd ld
# typically still links. Same symbol lives in w3 sbyte.o and g2 gbytesc.o.
add_link_options(-Wl,--allow-multiple-definition)

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

# Prefer GEPS_LIB OpenBLAS (same as GNU; intel/nvidia/fujitsu use vendor BLAS)
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
    "ROCm build requires OpenBLAS (GEPS_LIB openblas/* or system libopenblas). "
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

# ---------------------------------------------------------------------------
# AMD GPU: OpenMP offload (amdflang OpenACC cannot target AMDGPU) + HIP libs
# ---------------------------------------------------------------------------
if(${USE_HIP})
  if(NOT GEPS_ROCM_ROOT)
    message(FATAL_ERROR
      "USE_HIP=ON requires ROCm (set ROCM_PATH or HIP_PATH, typically /opt/rocm).")
  endif()

  if(AMD_GPU_ARCHS STREQUAL "" OR AMD_GPU_ARCHS STREQUAL "80")
    set(AMD_GPU_ARCHS "gfx942")
  endif()

  message(STATUS "AMD GPU offload: arch=${AMD_GPU_ARCHS}  ROCm=${GEPS_ROCM_ROOT}")

  # Reuse existing `#ifdef USE_CUDA` GPU code paths; shims live in src/rocm/.
  add_compile_options(-DUSE_CUDA=1 -DUSE_HIP=1 -DUSE_GPU=1)
  include_directories("${GEPS_ROCM_ROOT}/include"
                      "${GEPS_ROCM_ROOT}/include/hip"
                      "${CMAKE_SOURCE_DIR}/src/rocm")
  # amdflang does not inject OpenMPI's module dir; mpi.mod lives in lib/.
  foreach(_geps_mpi_hint
      "$ENV{INSTALL_DIR}/openmpi-4.1.6"
      "$ENV{MPI_DIR}"
      "$ENV{OPENMPI_ROOT}")
    if(_geps_mpi_hint AND EXISTS "${_geps_mpi_hint}/lib/mpi.mod")
      include_directories("${_geps_mpi_hint}/include" "${_geps_mpi_hint}/lib")
      break()
    endif()
  endforeach()

  # Host compile stays `-fopenmp` (set above). `--offload-arch` on every
  # translation unit makes device LTO link hundreds of empty amdgcn objects
  # and can take hours. GPU files get the arch flag in src/CMakeLists.txt.
  add_compile_options("$<$<COMPILE_LANGUAGE:C>:SHELL:-D__HIP_PLATFORM_AMD__>")
  add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:SHELL:-D__HIP_PLATFORM_AMD__>")
  add_link_options(-fopenmp --offload-arch=${AMD_GPU_ARCHS})

  link_directories("${GEPS_ROCM_ROOT}/lib" "${GEPS_ROCM_ROOT}/lib64")
  foreach(_geps_hip_lib hipblas hipfft hipsolver hipsparse rccl amdhip64)
    find_library(GEPS_${_geps_hip_lib}_LIB NAMES ${_geps_hip_lib}
      HINTS "${GEPS_ROCM_ROOT}/lib" "${GEPS_ROCM_ROOT}/lib64"
      NO_DEFAULT_PATH)
    if(NOT GEPS_${_geps_hip_lib}_LIB)
      find_library(GEPS_${_geps_hip_lib}_LIB NAMES ${_geps_hip_lib})
    endif()
    if(NOT GEPS_${_geps_hip_lib}_LIB)
      message(FATAL_ERROR "ROCm HIP lib not found: ${_geps_hip_lib} under ${GEPS_ROCM_ROOT}")
    endif()
    message(STATUS "  HIP lib ${_geps_hip_lib} = ${GEPS_${_geps_hip_lib}_LIB}")
    link_libraries(${GEPS_${_geps_hip_lib}_LIB})
  endforeach()
endif()
