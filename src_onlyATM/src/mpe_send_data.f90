      subroutine mpe_send_data(sbuf, n, j, ist)
!
#if defined(RSM) && defined(CWB_MPMD)
      use rank, only : root_io,MPI_COMM_gfs_all
#else
      use rank, only : root_io
#endif
      use mpi
      real*8  sbuf(n)

#if defined(RSM) && defined(CWB_MPMD)
      call MPI_SEND( SBUF, n, MPI_DOUBLE_PRECISION, root_io, J, &
                     MPI_COMM_gfs_all, ist )
#else
      call MPI_SEND( SBUF, n, MPI_DOUBLE_PRECISION, root_io, J, &
                     MPI_COMM_WORLD, ist )
#endif
 
      return
      end
