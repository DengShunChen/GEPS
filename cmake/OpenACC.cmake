# ##############################################################################
# Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
#
# See LICENSE for license information.
# ##############################################################################

# Specific flags for Fortran only
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Mfree -Ofast>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-r8>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Mpreprocess>")
add_compile_options(
  "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Mbyteswapio -Minline>")

# Auto-parallel and OpenMP
if(${USE_OMP})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-mp=multicore>")
endif()
if(${USE_PAR})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Mconcur>")
endif()

# Set variable for NetCDF and W3 libraries
set(NetCDF_Fortran_INCLUDE_DIRS "/package/x86_64/nvidia/netcdf-4.9.0/include")
set(W3_LIBRARIES "/package/x86_64/w3lib-2.0.2/lib/libw3.a")
set(CPL_INCLUDE_DIR "/data/common/gfs/GEPSv3_lib/coupler/a100")
set(MCT_INCLUDE_DIR "/package/x86_64/nvidia/mct-2.11.0/include")

# Link library for dgemm
link_directories(/home/xa09/pkg/openmpi-4.0.1/lib)
link_libraries(-lblas -llapack)

# Link library dms library
link_directories(/users/xa09/pkg/x86_64/dms38key/lib)
link_libraries(-lrdms -lgdbm)

# Link library zlib
link_directories(/usr/lib64)
link_libraries(-lz)

# Link library NetCDF
link_directories(/package/x86_64/nvidia/netcdf-4.9.0/lib)
link_libraries(-lnetcdf -lnetcdff)

# Link library FFTW
link_directories(/usr/lib64)
link_libraries(-lfftw3_threads -lfftw3 -lfftw3f_threads -lfftw3f)

# Link library operlib
link_directories(/package/x86_64/operlib/lib)
link_libraries(-lnwp)

link_directories(/package/x86_64/nvidia/geps_coupler/lib)
link_libraries(-lcpl)

link_directories(/package/x86_64/nvidia/mct-2.11.0/lib)
link_libraries(-lmct -lmpeu)

# Additional link
link_libraries(-ltirpc -lm -lcurl -lhdf5_hl -lhdf5 -lgfortran -lcusparse -lcudart)

# Add OpenACC options
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

# CPL
include_directories(${CPL_INCLUDE_DIR})
# MCT
include_directories(${MCT_INCLUDE_DIR})

