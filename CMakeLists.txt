# ##############################################################################
# Copyright (c) 2021, Fujitsu Deng-Shun Chen dschen@cwb.gov.tw Copyright (c)
# 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
#
# See LICENSE for license information.
# ##############################################################################

cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(TCO LANGUAGES Fortran CXX C)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# ##############################################################################
# Compiler Setting, TCO is built on top of MPI
# ##############################################################################

set(CMAKE_Fortran_COMPILER $ENV{COMP_MP})
set(CMAKE_CXX_COMPILER $ENV{CXX_COMP_MP})
set(CMAKE_C_COMPILER $ENV{C_COMP_MP})

# ##############################################################################
# Build Options
# ##############################################################################

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      "${default_build_type}"
      CACHE STRING "Choose the type of build." FORCE)
endif()
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

option(USE_ACC "enable OpenACC" OFF)
option(USE_CUDA "enable CUDA" OFF)
option(USE_NVTX "enable NVTX" OFF)
option(USE_OMP "enable OpenMP" ON)
option(USE_PAR "enable Auto Parallel" ON)
option(USE_RSM "enable output RSM initial/basefield data" OFF)
option(USE_MPMD "enable output RSM via MPMD " OFF)
option(USE_PCAST "enable OpenACC data correctness check" OFF)
option(USE_SP "enable partial computation with single precision" OFF)
option(TIMCOMCPL "enable CPL TIMCOM" ON)
option(CWBSUM "enable CWBSUM" ON)
option(USE_aeroclx "enable Readaeroclx" OFF)

set(GPU_ARCHS
    "80"
    CACHE STRING "GPU architecture(V100=70;A100=80)")
set(CUDA_RUNTIME_VERSION
    "12.6"
    CACHE STRING "CUDA runtime version")
# set(TCO_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/src" CACHE STRING "TCO destination")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(TCO_OUTPUT_PATH
    "${CMAKE_BINARY_DIR}/bin"
    CACHE STRING "TCO destination")

if(NOT ${USE_RSM} AND ${USE_MPMD})
  message(FATAL_ERROR "Invalid configuration. \
        `USE_RSM` must be ON before set `USE_MPMD` to ON")
endif()

message(STATUS "Build Summary")
message(STATUS "USE_ACC:      ${USE_ACC}")
message(STATUS "USE_CUDA:     ${USE_CUDA}")
message(STATUS "USE_NVTX:     ${USE_NVTX}")
message(STATUS "USE_OMP:      ${USE_OMP}")
message(STATUS "USE_PAR:      ${USE_PAR}")
message(STATUS "USE_RSM:      ${USE_RSM}")
message(STATUS "USE_MPMD:     ${USE_MPMD}")
message(STATUS "USE_PCAST:    ${USE_PCAST}")
message(STATUS "USE_SP:       ${USE_SP}")
message(STATUS "TIMCOMCPL:      ${TIMCOMCPL}")
message(STATUS "CWBSUM:      ${CWBSUM}")
message(STATUS "USE_aeroclx:  ${USE_aeroclx}")
message(STATUS "GPU_ARCHS:    ${GPU_ARCHS}")
message(STATUS "CUDA_RUNTIME_VERSION: ${CUDA_RUNTIME_VERSION}")
message(STATUS "Output path:  ${TCO_OUTPUT_PATH}")
message(STATUS "Install path: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "BUILD_TYPE : ${CMAKE_BUILD_TYPE}")

if(${USE_CUDA})
  include(OpenACC)
  if(${USE_NVTX})
    include(NVTX)
  endif()
else()
  include(NoCUDA)
endif()

if(NOT ${USE_CUDA})
  set(CMAKE_Fortran_MODDIR_FLAG "-M")
endif()
# ##############################################################################
# Compiler Options (-c -D<macro> -flag)
# ##############################################################################

add_definitions(-DTIMING -DUSE_FFTW -DW3TAG -Dcheck_info)
# add_definitions(-DTIMING -DUSE_FFTW -DW3TAG -DVERBOSE -DNO_OUT)
if(${USE_I38K})
  add_definitions(-DI38K)
endif()
if(${USE_O38K})
  add_definitions(-DO38K)
endif()
if(${USE_38K})
  add_definitions(-DI38K -DO38K)
endif()
if(${USE_RSM})
  add_definitions(-DRSM)
endif()
if(${USE_MPMD})
  add_definitions(-DCWB_MPM)
endif()
if(${USE_SP})
  add_definitions(-DSP)
endif()
if(${CWBSUM})
  add_definitions(-DCWBSUM)
endif()
if(${TIMCOMCPL})
  add_definitions(-DTIMCOMCPL)
endif()
if(${USE_aeroclx})
  add_definitions(-DReadaeroclx -DLUT_aero)
endif()

# ##############################################################################
# Include Directory (-I<path>)
# ##############################################################################

# NetCDF
include_directories(${NetCDF_Fortran_INCLUDE_DIRS})

# FFTW
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDE_DIR})

# G2
include_directories($ENV{G2_INC_PATH})


set(TCOGFS_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/src
    ${CMAKE_Fortran_MODULE_DIRECTORY})

# ##############################################################################
# Link Options (-L<path> -l<lib>), for dms, netcdf and Intel librares
# ##############################################################################

# W3
link_libraries($ENV{W3LIB_LIB})

# CWB operlib
link_libraries($ENV{OPERLIB_LIB})

# G2
link_directories($ENV{G2_LIB_PATH})
link_libraries(-lg2)

# PNG
link_directories($ENV{PNG_LIB_PATH})
link_libraries(-lpng)

# Jasper
link_directories($ENV{JAPSER_LIB_PATH})
link_libraries(-ljasper)

# ##############################################################################
# Build Target, TCO
# ##############################################################################
add_subdirectory(src)

# ##############################################################################
# Build Unit Tests
# ##############################################################################
include(CTest)
add_subdirectory(test)
