program t
  use mpi
  implicit none
  integer :: ierr, s, s2, s3, s4, s5
  call MPI_INIT(ierr)
  call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, s, ierr)
  call MPI_TYPE_SIZE(MPI_REAL8, s2, ierr)
  call MPI_TYPE_SIZE(MPI_REAL, s3, ierr)
  call MPI_TYPE_SIZE(MPI_REAL4, s4, ierr)
  call MPI_TYPE_SIZE(MPI_INTEGER, s5, ierr)
  print *, 'sizes: DOUBLE_PRECISION=', s, ' REAL8=', s2, ' REAL=', s3, ' REAL4=', s4, ' INTEGER=', s5
  print *, 'handles: DP=', MPI_DOUBLE_PRECISION, ' REAL8=', MPI_REAL8, ' REAL=', MPI_REAL, ' REAL4=', MPI_REAL4, ' INTEGER=', MPI_INTEGER
  print *, 'storage_size real(8)=', storage_size(1.0d0), ' default real=', storage_size(1.0)
  call MPI_FINALIZE(ierr)
end program
