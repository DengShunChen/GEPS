# ##############################################################################
# Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
#
# See LICENSE for license information.
# ##############################################################################

# Specific flags for Fortran only
add_compile_options(
  "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Kfast,ocl,autoobjstack>")
add_compile_options(
  "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-X9 -Free -CcdRR8 -Cpp -x- -fw -Knofp_relaxed>"
)
add_compile_options(
  "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Ec -Nlst=a,lst=d,lst=i,lst=p,lst=t,lst=x -Koptmsg=2 >"
)
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-SSL2BLAMP >")

# Specific flags for C only
set(CMAKE_C_FLAGS_RELEASE "")
add_compile_options("$<$<COMPILE_LANGUAGE:C>:SHELL:-Xg>")

#set(NetCDF_Fortran_INCLUDE_DIRS "/package/fx1000/netcdf-4.9.0/include")
set(CPL_INCLUDE_DIR "/data/common/gfs/GEPSv3_lib/coupler/fx1000")
set(MCT_INCLUDE_DIR "/data/common/gfs/GEPSv3_lib/fx1000/MCT/include")

# Auto-parallel and OpenMP
if(${USE_OMP})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Kopenmp>")
endif()
if(${USE_PAR})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Kparallel>")
endif()

# Find and link libraries
set(NETCDF_INCLUDE_DIR "/package/fx1000/netcdf-4.1.3/include")
set(NETCDF_LIB_DIR "/package/fx1000/netcdf-4.1.3/lib")

# Include directories for NetCDF
#include_directories(${NETCDF_INCLUDE_DIR})

# Link directories for NetCDF
#link_directories(${NETCDF_LIB_DIR})

# Link libraries for NetCDF
#link_libraries(-lnetcdf -lnetcdff)

find_package(NetCDF REQUIRED COMPONENTS Fortran)
link_libraries(-lnetcdf -lnetcdff -SSL2BLAMP)
find_package(DMS REQUIRED)
link_directories(/users/xa09/pkg/fx1000/dms38key/lib)
link_libraries(-lrdms -lgdbm -ltirpc)
if(${USE_MPMD})
  find_package(MPMD REQUIRED)
endif()

# Link library NetCDF
#link_directories(/package/fx1000/netcdf-4.9.0/lib)
#link_libraries(-lnetcdf -lnetcdff)

# Link library FFTW
link_directories(/users/xa09/pkg/fx1000/fftw-3.3.10/lib)
link_libraries(-lfftw3_threads -lfftw3 -lfftw3f_threads -lfftw3f)

link_directories(/data/common/gfs/GEPSv3_lib/coupler)
link_libraries(-lcpl)

link_directories(/data/common/gfs/GEPSv3_lib/fx1000/MCT/lib)
link_libraries(-lmct -lmpeu)

# CPL
include_directories(${CPL_INCLUDE_DIR})
# MCT
include_directories(${MCT_INCLUDE_DIR})
# NetCDF
#include_directories(${NetCDF_Fortran_INCLUDE_DIRS})

