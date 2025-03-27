      subroutine mpe_send_data(sbuf, n, j, ist)
!
#if defined(RSM) && defined(CWB_MPMD)
      use rank, only : root_io,MPI_COMM_gfs_all
#else
      use rank, only : root_io,MPI_COMM_atm
#endif
      use mpi

#ifdef SP
      real*4  sbuf(n)
#else
      real*8  sbuf(n)
#endif

#if defined(RSM) && defined(CWB_MPMD)

#ifdef SP
      call MPI_SEND( SBUF, n, MPI_REAL, root_io, J, &
                     MPI_COMM_gfs_all, ist )
#else
      call MPI_SEND( SBUF, n, MPI_DOUBLE_PRECISION, root_io, J, &
                     MPI_COMM_gfs_all, ist )
#endif
!=======================
#else

#ifdef SP
      call MPI_SEND( SBUF, n, MPI_REAL, root_io, J, &
                     MPI_COMM_WORLD, ist )
#else
      call MPI_SEND( SBUF, n, MPI_DOUBLE_PRECISION, root_io, J, &
                     MPI_COMM_atm, ist )
#endif

#endif
 
      return
      end
