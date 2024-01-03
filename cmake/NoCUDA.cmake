# Specific flags for Fortran only
set(CMAKE_Fortran_FLAGS_RELEASE "-Kfast,ocl,autoobjstack,simd=2")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-X9 -Free -CcdRR8 -Cpp -Cfpp -x- -fw -Knofp_relaxed>")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Ec -Nlst=a,lst=d,lst=i,lst=p,lst=t,lst=x -Koptmsg=2 >")
add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-SSL2BLAMP >")

# Specific flags for C only
set(CMAKE_C_FLAGS_RELEASE "")
add_compile_options("$<$<COMPILE_LANGUAGE:C>:SHELL:-Xg>")

# Auto-parallel and OpenMP
if (${USE_OMP})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Kopenmp>")
endif()
if (${USE_PAR})
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-Kparallel>")
endif()

# Find and link libraries
find_package( NetCDF REQUIRED COMPONENTS Fortran )
link_libraries(-SSL2BLAMP)
find_package( DMS REQUIRED )
if (${USE_MPMD})
 find_package( MPMD REQUIRED )
endif()
link_libraries(-SSL2BLAMP)

# For special source files
set_source_files_properties(fftx.f90 ndslfv_pack.f90
	PROPERTIES COMPILE_FLAGS "-Knoparallel"
)
set_source_files_properties(pbl_noah.f90 rcloud.f90
	PROPERTIES COMPILE_FLAGS "-Knosimd"
)
