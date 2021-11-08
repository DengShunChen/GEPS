# 
# Copyright (c) 2021, NVIDIA CORPORATION.
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(TCO LANGUAGES Fortran CXX C)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

###########################################################################
# Compiler Setting, TCO is built on top of MPI
###########################################################################

set(CMAKE_Fortran_COMPILER  mpifrtpx)
set(CMAKE_CXX_COMPILER      mpiFCCpx)
set(CMAKE_C_COMPILER        mpifccpx)

###########################################################################
# Build Options
###########################################################################
set(BUILD_TYPE "Release" CACHE STRING "Build type")

#option(USE_ACC   "enable OpenACC"         OFF)
#option(USE_CUDA  "enable CUDA"            OFF)
#option(USE_NVTX  "enable NVTX"            OFF)
option(USE_OMP   "enable OpenMP"          ON)
option(USE_PAR   "enable Auto Parallel"   ON)
#set(GPU_ARCHS "70" CACHE STRING "GPU architecture(V100=70;A100=80)")
#set(TCO_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/src" CACHE STRING "TCO destination")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(TCO_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin" CACHE STRING "TCO destination")

#if(NOT ${USE_ACC} AND ${USE_CUDA})
#    message(FATAL_ERROR "Invalid configuration. \
#        `USE_ACC` must be ON before set `USE_CUDA` to ON")
#endif()

message(STATUS "Build Summary")
#message(STATUS "USE_ACC:      ${USE_ACC}")
#message(STATUS "USE_CUDA:     ${USE_CUDA}")
#message(STATUS "USE_NVTX:     ${USE_NVTX}")
message(STATUS "USE_OMP:      ${USE_OMP}")
message(STATUS "USE_PAR:      ${USE_PAR}")
#message(STATUS "GPU_ARCHS:    ${GPU_ARCHS}")
message(STATUS "Output path:  ${TCO_OUTPUT_PATH}")
message(STATUS "Install path: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "BUILD_TYPE : ${BUILD_TYPE}")

# include CMake modules in `${PROJECT_SOURCE_DIR}/cmake/Modules`

#if(${USE_ACC})
#    include(OpenACC)
#endif()

#if(NOT ${USE_CUDA})
    # Use FFTW and Intel MKL
    #include(NoCUDA)
#endif()

#if(${USE_NVTX})
#    include(NVTX)
#endif()

###########################################################################
# Compiler Options (-c -D<macro> -flag)
###########################################################################

add_definitions(-DTIMING -DUSE_FFTW -DW3TAG)


# Specific flags for Fortran only:
if(BUILD_TYPE STREQUAL "Release" )
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Kfast,ocl,autoobjstack>")
endif()
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-X9 -Free -CcdRR8 -Cpp -x- -fw>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Ec -Qa,d,i,p,t,x -Koptmsg=2 >")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-SSL2BLAMP >")
if (${USE_OMP})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Kopenmp>")
endif()
if (${USE_PAR})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Kparallel>")
endif()

# Specific flags for C only:
add_compile_options("$<$<COMPILE_LANGUAGE:C>:SHELL:-Xg>")

###########################################################################
# Include Directory (-I<path>)
###########################################################################
# NetCDF
find_package( NetCDF REQUIRED COMPONENTS Fortran )
include_directories(${NetCDF_Fortran_INCLUDE_DIRS})

# FFTW
find_package( FFTW REQUIRED )
include_directories(${FFTW_INCLUDE_DIR})

set ( TCOGFS_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src
                          ${CMAKE_CURRENT_BINARY_DIR}/src
                          ${CMAKE_Fortran_MODULE_DIRECTORY} )


###########################################################################
# Link Options (-L<path> -l<lib>), for dms, netcdf and Intel librares
###########################################################################

# DMS
find_package( DMS REQUIRED )

# W3
link_libraries($ENV{W3NCO_LIBd})

# CWB operlib
link_libraries($ENV{OPERLIB_LIB})

# Fujitsu scientific and math libraries.
link_libraries(-SSL2BLAMP)

# Auto-parallel and OpenMP
if (${USE_OMP})
  link_libraries(-Kopenmp)
endif()
if (${USE_PAR})
  link_libraries(-Kparallel)
endif()

###########################################################################
# Build Target, TCO
###########################################################################
add_subdirectory(src)
add_subdirectory(test)

