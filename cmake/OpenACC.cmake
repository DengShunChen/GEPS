set(CMAKE_Fortran_COMPILER  mpif90)
set(CMAKE_CXX_COMPILER      mpic++)
set(CMAKE_C_COMPILER        mpicc)
add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-Mfree>)
add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-r8>)
add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-Mpreprocess>)
add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-Mbyteswapio>)

if (${USE_OMP})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-mp=multicore>")
endif()
if (${USE_PAR})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-stdpar>")
endif()

set(NetCDF_Fortran_INCLUDE_DIRS "/package/x86_64/nvidia/netcdf-4.9.0/include")
set(FFTW_INCLUDE_DIR "/usr/include")

link_directories(/home/xa09/pkg/openmpi-4.0.1/lib)
link_libraries(-lblas -llapack)
link_directories(/usr/lib64)
link_libraries(-lz)
link_directories(/package/x86_64/dms/dms.v4/lib)
link_libraries(-lrdms -lgdbm)

link_directories(/package/x86_64/nvidia/netcdf-4.9.0/lib)
link_libraries(-lnetcdf -lnetcdff)

link_directories(/usr/lib64)
link_libraries(-lfftw3_threads -lfftw3 -lfftw3f_threads -lfftw3f)

link_libraries(-ltirpc)

link_directories(/package/x86_64/operlib/lib)
link_libraries(-lnwp)

link_libraries(-lm -lcurl)

link_libraries(-lhdf5_hl -lhdf5)
link_libraries(-lgfortran)
set(W3_LIBRARIES "/package/x86_64/w3lib-2.0.2/lib/libw3.a")
