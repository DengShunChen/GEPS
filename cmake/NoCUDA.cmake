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

# Auto-parallel and OpenMP
if(${USE_OMP})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Kopenmp>")
endif()
if(${USE_PAR})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Kparallel>")
endif()

# Find and link libraries
find_package(NetCDF REQUIRED COMPONENTS Fortran)
link_libraries(-SSL2BLAMP)
find_package(DMS REQUIRED)
link_directories(/users/xa09/pkg/fx1000/dms38key/lib)
link_libraries(-lrdms -lgdbm -ltirpc)
if(${USE_MPMD})
  find_package(MPMD REQUIRED)
endif()
# Link library FFTW
link_directories(/users/xa09/pkg/fx1000/fftw-3.3.10/lib)
link_libraries(-lfftw3_threads -lfftw3 -lfftw3f_threads -lfftw3f)
