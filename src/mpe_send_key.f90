      subroutine mpe_send_key(key, j, ist)
!
#if defined(RSM) && defined(CWB_MPMD)
      use rank, only : root_io,MPI_COMM_gfs_all
#else
      use rank, only : root_io,MPI_COMM_atm
#endif
      use const,only : keyo,KLENO
      use mpi

#if defined(RSM) && defined(CWB_MPMD)
      call MPI_SEND( keyo, KLENO, MPI_CHARACTER, root_io, J, &
                     MPI_COMM_gfs_all, ist )
#else
      call MPI_SEND( keyo, KLENO, MPI_CHARACTER, root_io, J, &
                     MPI_COMM_WORLD, ist )
#endif
      return
      end
