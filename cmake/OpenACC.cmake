# Specific flags for Fortran only
add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-Mfree>)
add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-r8>)
add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-Mpreprocess>)
add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-Mbyteswapio>)

# Auto-parallel and OpenMP
if (${USE_OMP})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-mp=multicore>")
endif()
if (${USE_PAR})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-stdpar>")
endif()

# Set variable for NetCDF and W3 libraries
set(NetCDF_Fortran_INCLUDE_DIRS "/package/x86_64/nvidia/netcdf-4.9.0/include")
set(W3_LIBRARIES "/package/x86_64/w3lib-2.0.2/lib/libw3.a")

# Link library for dgemm
link_directories(/home/xa09/pkg/openmpi-4.0.1/lib)
link_libraries(-lblas -llapack)

# Link library dms library
link_directories(/package/x86_64/dms/dms.v4/lib)
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

# Additional link
link_libraries(-ltirpc -lm -lcurl -lhdf5_hl -lhdf5 -lgfortran)
